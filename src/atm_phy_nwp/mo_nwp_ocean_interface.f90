!>
!! @brief Interface between NWP physics and the ocean, through a coupler
!!
!! @author Martin Koehler (DWD) and Rene Redler (MPI-M) based on ECHAM version by Marco Giorgetta (MPI-M)
!!
!! @par Revision History
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_nwp_ocean_interface

  USE mo_kind                ,ONLY: wp
  USE mo_model_domain        ,ONLY: t_patch
  USE mo_nonhydro_types      ,ONLY: t_nh_diag
  USE mo_nwp_phy_types       ,ONLY: t_nwp_phy_diag
  USE mo_nwp_lnd_types       ,ONLY: t_wtr_prog, t_lnd_diag
  USE mo_ext_data_types      ,ONLY: t_external_data
  USE mo_ccycle_config       ,ONLY: ccycle_config
                                
  USE mo_parallel_config     ,ONLY: nproma
  USE mo_impl_constants_grf  ,ONLY: grf_bdywidth_c
  USE mo_impl_constants      ,ONLY: min_rlcell, min_rlcell_int
  USE mo_loopindices         ,ONLY: get_indices_c

  USE mo_run_config          ,ONLY: ltimer, ico2, nlev
  USE mo_timer               ,ONLY: timer_start, timer_stop,                &
       &                            timer_coupling_put, timer_coupling_get, &
       &                            timer_coupling_1stget

  USE mo_lnd_nwp_config      ,ONLY: lseaice, isub_water, isub_seaice, isub_lake

  USE mo_sync                ,ONLY: sync_c, sync_patch_array

  USE mo_bc_greenhouse_gases ,ONLY: ghg_co2mmr

#ifdef YAC_coupling
#if !defined(__NO_JSBACH__) && !defined(__NO_JSBACH_HD__)
  USE mo_interface_hd_ocean  ,ONLY: jsb_fdef_hd_fields
#endif
#endif

  USE mo_parallel_config     ,ONLY: nproma

  USE mo_coupling_config     ,ONLY: is_coupled_run
#ifdef YAC_coupling
  USE mo_atmo_coupling_frame ,ONLY: lyac_very_1st_get, nbr_inner_cells,     &
    &                               mask_checksum, field_id
  USE mo_yac_finterface      ,ONLY: yac_fput, yac_fget,                     &
    &                               COUPLING, OUT_OF_BOUND
#endif

  USE mo_exception           ,ONLY: warning, message, finish
  USE mo_util_dbg_prnt       ,ONLY: dbg_print
  USE mo_dbg_nml             ,ONLY: idbg_mxmn, idbg_val
  USE mo_physical_constants  ,ONLY: amd, amco2
  USE mo_lnd_nwp_config      ,ONLY: hice_max                      ! maximum sea-ice thickness [m]

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: nwp_couple_ocean

  CHARACTER(len=*), PARAMETER :: str_module = 'mo_nwp_ocean_interface' ! Output of module for debug

CONTAINS


  !>
  !! SUBROUTINE nwp_couple_ocean -- the interface between
  !! NWP physics and the ocean, through a coupler
  !!
  !! This subroutine is called from nwp_nh_interface.

  SUBROUTINE nwp_couple_ocean( p_patch, pt_diag, lnd_diag, &
    &                          wtr_prog_now, wtr_prog_new, prm_diag, ext_data )

    ! Arguments

    TYPE(t_patch),   TARGET, INTENT(INOUT)  :: p_patch
    TYPE(t_nh_diag), TARGET, INTENT(INOUT)  :: pt_diag
    TYPE(t_wtr_prog),        INTENT(INOUT)  :: wtr_prog_now, wtr_prog_new
    TYPE(t_lnd_diag),        INTENT(INOUT)  :: lnd_diag
    TYPE(t_nwp_phy_diag),    INTENT(INOUT)  :: prm_diag
    TYPE(t_external_data),   INTENT(INOUT)  :: ext_data

    ! Local variables

    LOGICAL               :: write_coupler_restart
    INTEGER               :: nbr_hor_cells         ! inner points
    INTEGER               :: jg                    ! grid index
    INTEGER               :: nn                    ! block offset
    INTEGER               :: jb                    ! block loop count
    INTEGER               :: jc                    ! nproma loop count
    INTEGER               :: ic                    ! nproma loop count
    INTEGER               :: nlen                  ! nproma/npromz
    INTEGER               :: info, ierror          ! return values from cpl_put/get calls
    INTEGER               :: no_arr                ! no of arrays in bundle for put/get calls
    INTEGER               :: rl_start, rl_end
    INTEGER               :: i_startblk, i_endblk  ! blocks
    INTEGER               :: i_startidx, i_endidx  ! slices
    INTEGER               :: i_nchdom              ! domain index
    INTEGER               :: ncount                ! buffer counter

    REAL(wp)              :: scr(nproma,p_patch%alloc_cell_blocks)
    REAL(wp)              :: frac_oce(nproma,p_patch%alloc_cell_blocks)
    REAL (wp), PARAMETER  :: csmall = 1.0E-5_wp    ! small number (security constant)
    REAL(wp), ALLOCATABLE :: buffer(:,:)           ! buffer transferred to YAC coupler

    IF ( .NOT. is_coupled_run() ) RETURN

#ifndef YAC_coupling
    CALL finish('nwp_couple_ocean: unintentionally called. Check your source code and configure.')
#else

    ALLOCATE(buffer(nproma*p_patch%nblks_c,5))

    jg = p_patch%id
    i_nchdom  = MAX(1,p_patch%n_childdom)

    ! include boundary interpolation zone of nested domains and halo points
    rl_start = 1
    rl_end   = min_rlcell

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

    nbr_hor_cells = p_patch%n_patch_cells

    !-------------------------------------------------------------------------
    ! If running in atm-oce coupled mode, exchange information 
    !-------------------------------------------------------------------------
    !
    ! Possible fields that contain information to be sent to the ocean include
    !
    ! 1. prm_diag%umfl_s_t(:,:,:)                       zonal resolved surface stress and  [N/m2] 
    ! 2. prm_diag%vmfl_s_t(:,:,:)                       meridional resolved surface stress [N/m2]
    !
    ! 3. prm_diag%rain_con_rate(:,:)                    convective surface rain rate    [kg/m2/s]
    !    prm_diag%rain_gsp_rate(:,:)                    grid-scale surface rain rate    [kg/m2/s]
    !
    !    prm_diag%snow_con_rate(:,:)                    convective surface snow_rate    [kg/m2/s]
    !    prm_diag%snow_gsp_rate    (:,:)                grid_scale surface snow rate    [kg/m2/s]
    !    prm_diag%ice_gsp_rate     (:,:)                grid_scale surface ice rate     [kg/m2/s]
    !    prm_diag%graupel_gsp_rate (:,:)                grid_scale surface graupel rate [kg/m2/s]
    !    prm_diag%hail_gsp_rate    (:,:)                grid_scale surface hail rate    [kg/m2/s]
    !
    !    prm_diag%qhfl_s_t(:,:,isub_water/isub_seaice)  moisture flux (surface) aka evaporation rate at surface [Kg/m2/s]
    
    ! 4. prm_diag%swflxsfc_t (:,:,:)                    tile-based shortwave net flux at surface [W/m2]
    !    prm_diag%lwflxsfc_t (:,:,:)                    tile-based longwave net flux at surface  [W/m2]
    !    prm_diag%shfl_s_t   (:,:,:)                    tile-based sensible heat flux at surface [W/m2]
    !    prm_diag%lhfl_s_t   (:,:,:)                    tile-based latent   heat flux at surface [W/m2]
    !
    ! Possible fields to receive from the ocean include
    !
    ! 1. lnd_diag%t_seasfc     (:,:)   SST
    !    ... tbc
    !
    !
    !  Send fields to ocean:
    !   field_id(1)  represents "surface_downward_eastward_stress" bundle  - zonal wind stress component over ice and water
    !   field_id(2)  represents "surface_downward_northward_stress" bundle - meridional wind stress component over ice and water
    !   field_id(3)  represents "surface_fresh_water_flux" bundle          - liquid rain, snowfall, evaporation
    !   field_id(4)  represents "total heat flux" bundle                   - short wave, long wave, sensible, latent heat flux
    !   field_id(5)  represents "atmosphere_sea_ice_bundle"                - sea ice surface and bottom melt potentials
    !   field_id(10) represents "10m_wind_speed"                           - atmospheric wind speed
    !   field_id(11) represents "qtrc(nlev,co2)"                           - co2 mixing ratio
    !   field_id(13) represents "pres_msl"                                 - sea level pressure
    !
    !  Receive fields from ocean:
    !   field_id(6)  represents "sea_surface_temperature"                  - SST
    !   field_id(7)  represents "eastward_sea_water_velocity"              - zonal velocity, u component of ocean surface current
    !   field_id(8)  represents "northward_sea_water_velocity"             - meridional velocity, v component of ocean surface current
    !   field_id(9)  represents "ocean_sea_ice_bundle"                     - ice thickness, snow thickness, ice concentration
    !   field_id(12) represents "co2_flux"                                 - ocean co2 flux
    !
    !-------------------------------------------------------------------------


!! WRITE(*,*) 'ocean_interface 1 '

!  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
!  Send fields from atmosphere to ocean
!  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****

    write_coupler_restart = .FALSE.

    ! Calculate fractional ocean mask 
    ! evaporation over ice-free and ice-covered water fraction, of whole ocean part, without land part
    !  - lake part is included in land part, must be subtracted as well
    IF ( mask_checksum > 0 ) THEN

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jc, i_startidx, i_endidx) ICON_OMP_GUIDED_SCHEDULE
      DO jb = i_startblk, i_endblk
  
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         & i_startidx, i_endidx, rl_start, rl_end)
  
        DO jc = i_startidx, i_endidx
          IF ( lseaice ) THEN
            frac_oce(jc,jb) = ext_data%atm%frac_t(jc,jb,isub_water)+ ext_data%atm%frac_t(jc,jb,isub_seaice) ! sea ice
          ELSE
            frac_oce(jc,jb) = ext_data%atm%frac_t(jc,jb,isub_water)                                         ! open sea
          END IF
        ENDDO

      ENDDO
!ICON_OMP_END_DO
!ICON_OMP_END_PARALLEL

    ELSE

!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(jb, jc, i_startidx, i_endidx) ICON_OMP_RUNTIME_SCHEDULE
      DO jb = i_startblk, i_endblk
  
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)
  
        DO jc = i_startidx, i_endidx
          frac_oce(jc,jb) = 1.0
        ENDDO

      ENDDO
!ICON_OMP_END_DO
!ICON_OMP_END_PARALLEL

    ENDIF


!------------------------------------------------
!  Send zonal wind stress bundle
!    field_id(1) represents "surface_downward_eastward_stress" bundle 
!    - zonal wind stress component over ice and water
!------------------------------------------------

    ncount = 0
!ICON_OMP_PARALLEL_DO PRIVATE(jb, jc, i_startidx, i_endidx, ncount) ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                       & i_startidx, i_endidx, rl_start, rl_end)
      DO jc = i_startidx, i_endidx
        ncount = ncount + 1
        buffer(ncount,1) = prm_diag%umfl_s_t(jc,jb,isub_water)      ! open sea
        buffer(ncount,2) = prm_diag%umfl_s_t(jc,jb,isub_seaice)     ! sea ice
      ENDDO
    ENDDO
!ICON_OMP_END_PARALLEL_DO

!  The nwp_ocean_interface and YAC has not yet been adapted for nested ICON setups.  
!  It only works for global configuations.
    IF (ncount /= nbr_hor_cells) THEN
      WRITE(*,*) 'ncount, nbr_hor_cells', ncount, nbr_hor_cells
      CALL finish (TRIM(str_module), 'unequal ncount and nbr_hor_cells')
    ENDIF

    IF (ltimer) CALL timer_start(timer_coupling_put)

    no_arr = 2
    CALL yac_fput ( field_id(1), nbr_hor_cells, no_arr, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    IF ( info > COUPLING .AND. info < OUT_OF_BOUND ) write_coupler_restart = .TRUE.
    IF ( info == OUT_OF_BOUND ) &
         & CALL warning('nwp_couple_ocean', 'YAC says fput called after end of run - id=1, u-stress')

    IF (ltimer) CALL timer_stop(timer_coupling_put)


!------------------------------------------------
!  Send meridional wind stress bundle
!    field_id(2) represents "surface_downward_northward_stress" bundle 
!    - meridional wind stress component over ice and water
!------------------------------------------------

    ncount = 0
!ICON_OMP_PARALLEL_DO PRIVATE(jb, jc, i_startidx, i_endidx, ncount) ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                       & i_startidx, i_endidx, rl_start, rl_end)
      DO jc = i_startidx, i_endidx
        ncount = ncount + 1
        buffer(ncount,1) = prm_diag%vmfl_s_t(jc,jb,isub_water)      ! open sea
        buffer(ncount,2) = prm_diag%vmfl_s_t(jc,jb,isub_seaice)     ! sea ice 
      ENDDO
    ENDDO
!ICON_OMP_END_PARALLEL_DO

    IF (ltimer) CALL timer_start(timer_coupling_put)

    no_arr = 2
    CALL yac_fput ( field_id(2), nbr_hor_cells, no_arr, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    IF ( info > COUPLING .AND. info < OUT_OF_BOUND ) write_coupler_restart = .TRUE.
    IF ( info == OUT_OF_BOUND ) &
         & CALL warning('nwp_couple_ocean', 'YAC says fput called after end of run - id=2, v-stress')

    IF (ltimer) CALL timer_stop(timer_coupling_put)


!------------------------------------------------
!  Send surface fresh water flux bundle
!    field_id(3) represents "surface_fresh_water_flux" bundle 
!    - liquid rain, snowfall, evaporation
!
!    Note: the evap_tile should be properly updated and added;
!          as long as evaporation over sea-ice is not used in ocean thermodynamics, the evaporation over the
!          whole ocean part of grid-cell is passed to the ocean
!          for pre04 a preliminary solution for evaporation in ocean model is to exclude the land fraction
!          evap.oce = (evap.wtr*frac.wtr + evap.ice*frac.ice)/(1-frac.lnd)
!------------------------------------------------

    IF ( idbg_mxmn >= 1 .OR. idbg_val >=1 )  &
      scr(:,:) = 0.0_wp

    ncount = 0
!ICON_OMP_PARALLEL_DO PRIVATE(jb, jc, i_startidx, i_endidx, ncount) ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                       & i_startidx, i_endidx, rl_start, rl_end)
      DO jc = i_startidx, i_endidx
        ncount = ncount + 1
    
        ! total rates of rain and snow over whole cell
        buffer(ncount,1) = prm_diag%rain_con_rate(jc,jb) + prm_diag%rain_gsp_rate   (jc,jb)
        buffer(ncount,2) = prm_diag%snow_con_rate(jc,jb) + prm_diag%snow_gsp_rate   (jc,jb)!+ &
!             &            prm_diag%hail_gsp_rate(jc,jb) + prm_diag%graupel_gsp_rate(jc,jb) + &  ! when available
!             &            prm_diag%ice_gsp_rate (jc,jb)                                         ! when available

        ! evaporation over ice-free and ice-covered water fraction, of whole ocean part, without land part
        IF (frac_oce(jc,jb) <= 0.0_wp) THEN
          ! ocean part is zero
          buffer(ncount,3) = 0.0_wp
        ELSE  ! attention: this will underestimate fluxes for fractional coast
          buffer(ncount,3) = prm_diag%qhfl_s_t(jc,jb,isub_water)  * ext_data%atm%frac_t(jc,jb,isub_water) + &
            &                prm_diag%qhfl_s_t(jc,jb,isub_seaice) * ext_data%atm%frac_t(jc,jb,isub_seaice)
        ENDIF
        IF ( idbg_mxmn >= 1 .OR. idbg_val >=1 ) scr(jc,jb) = buffer(ncount,3)
      ENDDO
    ENDDO
!ICON_OMP_END_PARALLEL_DO

    IF ( idbg_mxmn >= 1 .OR. idbg_val >=1 )  &
      &  CALL dbg_print('NWPOce: evapo-cpl',scr,str_module,3,in_subset=p_patch%cells%owned)

    IF (ltimer) CALL timer_start(timer_coupling_put)

    no_arr = 3
    CALL yac_fput ( field_id(3), nbr_hor_cells, no_arr, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    IF ( info > COUPLING .AND. info < OUT_OF_BOUND ) write_coupler_restart = .TRUE.
    IF ( info == OUT_OF_BOUND )                  &
      & CALL warning('nwp_couple_ocean', 'YAC says fput called after end of run - id=3, fresh water flux')

    IF (ltimer) CALL timer_stop(timer_coupling_put)


!------------------------------------------------
!  Send total heat flux bundle
!    field_id(4) represents "total heat flux" bundle 
!    - short wave, long wave, sensible, latent heat flux
!------------------------------------------------

    ncount = 0
!ICON_OMP_PARALLEL_DO PRIVATE(jb, jc, i_startidx, i_endidx, ncount) ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                       & i_startidx, i_endidx, rl_start, rl_end)
      DO jc = i_startidx, i_endidx
        ncount = ncount + 1
        buffer(ncount,1) = prm_diag%swflxsfc_t (jc,jb,isub_water)
        buffer(ncount,2) = prm_diag%lwflxsfc_t (jc,jb,isub_water)
        buffer(ncount,3) = prm_diag%shfl_s_t   (jc,jb,isub_water)
        buffer(ncount,4) = prm_diag%lhfl_s_t   (jc,jb,isub_water)
      ENDDO
    ENDDO
!ICON_OMP_END_PARALLEL_DO

    IF (ltimer) CALL timer_start(timer_coupling_put)

    no_arr = 4
    CALL yac_fput ( field_id(4), nbr_hor_cells, no_arr, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    IF ( info > COUPLING .AND. info < OUT_OF_BOUND ) write_coupler_restart = .TRUE.
    IF ( info == OUT_OF_BOUND ) &
         & CALL warning('nwp_couple_ocean', 'YAC says fput called after end of run - id=4, heat flux')

    IF (ltimer) CALL timer_stop(timer_coupling_put)


!------------------------------------------------
!  Send sea ice flux bundle
!    field_id(5) represents "atmosphere_sea_ice_bundle" 
!    - sea ice surface and bottom melt potentials Qtop, Qbot (conductive heat flux)
!------------------------------------------------

    ncount = 0
!ICON_OMP_PARALLEL_DO PRIVATE(jb, jc, i_startidx, i_endidx, ncount) ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                       & i_startidx, i_endidx, rl_start, rl_end)
      DO jc = i_startidx, i_endidx
        ncount = ncount + 1
        buffer(ncount,1) = prm_diag%shfl_s_t(jc,jb,isub_seaice) + prm_diag%swflxsfc_t(jc,jb,isub_seaice) &
                       & + prm_diag%lhfl_s_t(jc,jb,isub_seaice) + prm_diag%lwflxsfc_t(jc,jb,isub_seaice) 
        buffer(ncount,2) = lnd_diag%condhf_ice(jc,jb)  !  melt potential via conductive heat flux at bottom of sea-ice
      ENDDO
    ENDDO
!ICON_OMP_END_PARALLEL_DO

    IF (ltimer) CALL timer_start(timer_coupling_put)

    no_arr = 2
    CALL yac_fput ( field_id(5), nbr_hor_cells, no_arr, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    IF ( info > COUPLING .AND. info < OUT_OF_BOUND ) write_coupler_restart = .TRUE.
    IF ( info == OUT_OF_BOUND )  &
      & CALL warning('nwp_couple_ocean', 'YAC says fput called after end of run - id=5, atmos sea ice')

    IF (ltimer) CALL timer_stop(timer_coupling_put)

    IF ( write_coupler_restart ) THEN
      CALL message('nwp_couple_ocean', 'YAC says it is put for restart - ids 1 to 5, atmosphere fields')
    ENDIF


!------------------------------------------------
!  Send 10m wind speed
!    field_id(10) represents "10m_wind_speed" 
!    - atmospheric wind speed
!------------------------------------------------

    ncount = 0
!ICON_OMP_PARALLEL_DO PRIVATE(jb, jc, i_startidx, i_endidx, ncount) ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                       & i_startidx, i_endidx, rl_start, rl_end)
      DO jc = i_startidx, i_endidx
        ncount = ncount + 1
        ! attention: using the grid-point mean of 10m wind instead of ocean wind
        buffer(ncount,1) = prm_diag%sp_10m(jc,jb)
       ENDDO
    ENDDO
!ICON_OMP_END_PARALLEL_DO

    IF (ltimer) CALL timer_start(timer_coupling_put)

    no_arr = 1
    CALL yac_fput ( field_id(10), nbr_hor_cells, no_arr, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    IF ( info > COUPLING .AND. info < OUT_OF_BOUND ) THEN
      write_coupler_restart = .TRUE.
    ELSE
      write_coupler_restart = .FALSE.
    ENDIF

    IF ( info == OUT_OF_BOUND )  &
       & CALL warning('nwp_couple_ocean', 'YAC says fput called after end of run - id=10, wind speed')

    IF (ltimer) CALL timer_stop(timer_coupling_put)

    IF ( write_coupler_restart ) THEN
      CALL message('nwp_couple_ocean', 'YAC says it is put for restart - ids 10, wind speed')
    ENDIF


!------------------------------------------------
!  Send sea level pressure
!    field_id(13) represents "pres_msl" 
!    - atmospheric sea level pressure
!    - pres_sfc is used insted of pres_msl because
!      * it is available at each fast physics timestep
!      * it calculated the hydrostatic surface pressure with less noise
!------------------------------------------------

    ncount = 0
!ICON_OMP_PARALLEL_DO PRIVATE(jb, jc, i_startidx, i_endidx, ncount) ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                       & i_startidx, i_endidx, rl_start, rl_end)
      DO jc = i_startidx, i_endidx
        ncount = ncount + 1
        buffer(ncount,1) = pt_diag%pres_sfc(jc,jb)
       ENDDO
    ENDDO
!ICON_OMP_END_PARALLEL_DO

    IF (ltimer) CALL timer_start(timer_coupling_put)

    no_arr = 1
    CALL yac_fput ( field_id(13), nbr_hor_cells, no_arr, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    IF ( info > COUPLING .AND. info < OUT_OF_BOUND ) THEN
       write_coupler_restart = .TRUE.
    ELSE
       write_coupler_restart = .FALSE.
    ENDIF

    IF ( info == OUT_OF_BOUND )   &
       & CALL warning('nwp_couple_ocean', &
       &              'YAC says fput called after end of run - id=13, sea level pressure')

    IF (ltimer) CALL timer_stop(timer_coupling_put)

    IF ( write_coupler_restart ) THEN
       CALL message('nwp_couple_ocean', 'YAC says it is put for restart - ids 13, sea level pressure')
    ENDIF


!------------------------------------------------
!  Send co2 mixing ratio
!    field_id(11) represents "co2_mixing_ratio" 
!    - CO2 mixing ratio in ppmv
!------------------------------------------------

#ifndef __NO_ICON_OCEAN__
    IF (ccycle_config(jg)%iccycle /= 0) THEN

      ncount = 0
!ICON_OMP_PARALLEL_DO PRIVATE(jb, jc, i_startidx, i_endidx, ncount) ICON_OMP_RUNTIME_SCHEDULE
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         & i_startidx, i_endidx, rl_start, rl_end)

        SELECT CASE (ccycle_config(jg)%iccycle)
        CASE (1) ! c-cycle with interactive atm. co2 concentration, qtrc in kg/kg
          DO jc = i_startidx, i_endidx
            ncount = ncount + 1
!ECHAM      buffer(ncount,1)     =  amd/amco2 * 1.0e6_wp * prm_field(jg)%qtrc(n,nlev,i_blk,ico2)
!NWP:  prognostic CO2 not yet available in NWP physics
            buffer(ncount,1)    =  0.0_wp
          END DO
        CASE (2) ! c-cycle with prescribed  atm. co2 concentration
          SELECT CASE (ccycle_config(jg)%ico2conc)
          CASE (2) ! constant  co2 concentration, vmr_co2 in m3/m3
            DO jc = i_startidx, i_endidx
              ncount = ncount + 1
              buffer(ncount,1) = 1.0e6_wp * ccycle_config(jg)%vmr_co2
            END DO
          CASE (4) ! transient co2 concentration, ghg_co2mmr in kg/kg
            DO jc = i_startidx, i_endidx
              ncount = ncount + 1
              buffer(ncount,1) = amd/amco2 * 1.0e6_wp * ghg_co2mmr
            END DO
          END SELECT
        END SELECT

      ENDDO
!ICON_OMP_END_PARALLEL_DO

      IF (ltimer) CALL timer_start(timer_coupling_put)

      no_arr = 1
      CALL yac_fput ( field_id(11), nbr_hor_cells, no_arr, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
      IF ( info > COUPLING .AND. info < OUT_OF_BOUND ) THEN
        write_coupler_restart = .TRUE.
      ELSE
        write_coupler_restart = .FALSE.
      ENDIF

      IF ( info == OUT_OF_BOUND )  &
         & CALL warning('nwp_couple_ocean', 'YAC says fput called after end of run - id=11, co2 mr')

      IF (ltimer) CALL timer_stop(timer_coupling_put)

      IF ( write_coupler_restart ) THEN
        CALL message('nwp_couple_ocean', 'YAC says it is put for restart - id=11, co2 mr')
      ENDIF

    ENDIF
#endif



!  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
!  Receive fields from ocean to atmosphere
!  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
!
!  Receive fields, only assign values if something was received ( info > 0 )
!   - ocean fields have undefined values on land, which are not sent to the atmosphere,
!     therefore buffer is set to zero to avoid unintended usage of ocean values over land

    buffer(:,:) = 0.0_wp

    ! exclude nest boundary and halo points
    rl_start   = grf_bdywidth_c+1          
    rl_end     = min_rlcell_int
    i_nchdom   = MAX(1,p_patch%n_childdom)

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)


!------------------------------------------------
!  Receive SST
!    field_id(6) represents "sea_surface_temperature" 
!    - SST
!------------------------------------------------

    IF ( .NOT. lyac_very_1st_get ) THEN
      IF (ltimer) CALL timer_start(timer_coupling_1stget)
    ENDIF

    CALL yac_fget ( field_id(6), nbr_hor_cells, 1, buffer(1:nbr_hor_cells,1:1), info, ierror )
    IF ( info > COUPLING .AND. info < OUT_OF_BOUND ) &
         & CALL message('nwp_couple_ocean', 'YAC says it is get for restart - id=6, SST')
    IF ( info == OUT_OF_BOUND ) &
         & CALL warning('nwp_couple_ocean', 'YAC says fget called after end of run - id=6, SST')

    IF ( .NOT. lyac_very_1st_get ) THEN
       IF (ltimer) CALL timer_stop(timer_coupling_1stget)
    ENDIF

    lyac_very_1st_get = .FALSE.

    IF ( info > 0 .AND. info < 7 ) THEN

      ! new lists are calculated in process_sst_and_seaice in mo_nwp_sfc_utils.f90

      scr(:,:) = 285.0_wp  !  value over land - for dbg_print

!ICON_OMP_PARALLEL_DO PRIVATE(jb, ic, jc, nn, scr) ICON_OMP_RUNTIME_SCHEDULE
      DO jb = i_startblk, i_endblk
        nn = (jb-1)*nproma                                            ! translation to 1-d buffer fields

        DO ic = 1, ext_data%atm%list_sea%ncount(jb)                   ! number of ocean points (open water & sea ice)
          jc = ext_data%atm%list_sea%idx(ic,jb)                       ! index list of ocean points

          lnd_diag%t_seasfc(jc,jb) = buffer(nn+jc,1)

          !  for dbg_print only
          IF ( idbg_mxmn >= 1 .OR. idbg_val >=1 ) THEN
            scr(jc,jb) = buffer(nn+jc,1)
          ENDIF

        END DO
      ENDDO
!ICON_OMP_END_PARALLEL_DO

      IF ( idbg_mxmn >= 1 .OR. idbg_val >=1 )  &
        &  CALL dbg_print('NWPOce: SSToce-cpl',scr,str_module,3,in_subset=p_patch%cells%owned)

      CALL sync_patch_array(sync_c, p_patch, lnd_diag%t_seasfc(:,:) )

    END IF


!------------------------------------------------
!  Receive zonal velocity
!    field_id(7) represents "eastward_sea_water_velocity" 
!    - zonal velocity, u component of ocean surface current
!    RR: not used in NWP so far, not activated for exchange in coupling.xml
!------------------------------------------------

    IF (ltimer) CALL timer_start(timer_coupling_get)
    CALL yac_fget ( field_id(7), nbr_hor_cells, 1, buffer(1:nbr_hor_cells,1:1), info, ierror )
    IF ( info > COUPLING .AND. info < OUT_OF_BOUND ) &
         & CALL message('nwp_couple_ocean', 'YAC says it is get for restart - id=7, u velocity')
    IF ( info == OUT_OF_BOUND ) &
         & CALL warning('nwp_couple_ocean', 'YAC says fget called after end of run - id=7, u velocity')
    IF (ltimer) CALL timer_stop(timer_coupling_get)

!--- no ocean current in turbtran yet
!     IF ( info > 0 .AND. info < 7 ) THEN
! 
! !ICON_OMP_PARALLEL_DO PRIVATE(jb, ic, jc, nn) ICON_OMP_RUNTIME_SCHEDULE
!       DO jb = i_startblk, i_endblk
!         nn = (jb-1)*nproma                                            ! translation to 1-d buffer fields
! 
!         DO ic = 1, ext_data%atm%list_sea%ncount(jb)                   ! number of ocean points (open water & sea ice)
!           jc = ext_data%atm%list_sea%idx(ic,jb)                       ! index list of ocean points
!          !ECHAM  prm_field(jg)%ocu(n,i_blk) = buffer(nn+jc,1)
!          !NWP: no ocean current in turbtran yet
!         END DO
! 
!       ENDDO
! !ICON_OMP_END_PARALLEL_DO
! 
!       !ECHAM CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%ocu(:,:))
!       !NWP: no ocean current in turbtran yet
! 
!     END IF


!------------------------------------------------
!  Receive meridional velocity
!    field_id(8) represents "northward_sea_water_velocity" 
!    - meridional velocity, v component of ocean surface current
!    RR: not used in NWP so far, not activated for exchange in coupling.xml
!------------------------------------------------

    IF (ltimer) CALL timer_start(timer_coupling_get)
    CALL yac_fget ( field_id(8), nbr_hor_cells, 1, buffer(1:nbr_hor_cells,1:1), info, ierror )
    IF ( info > COUPLING .AND. info < OUT_OF_BOUND ) &
         & CALL message('nwp_couple_ocean', 'YAC says it is get for restart - id=8, v velocity')
    IF ( info == OUT_OF_BOUND ) &
         & CALL warning('nwp_couple_ocean', 'YAC says fget called after end of run - id=8, v velocity')
    IF (ltimer) CALL timer_stop(timer_coupling_get)

!--- no ocean current in turbtran yet
!     IF ( info > 0 .AND. info < 7 ) THEN
! 
! !ICON_OMP_PARALLEL_DO PRIVATE(jb, ic, jc, nn) ICON_OMP_RUNTIME_SCHEDULE
!       DO jb = i_startblk, i_endblk
!         nn = (jb-1)*nproma                                            ! translation to 1-d buffer fields
! 
!         DO ic = 1, ext_data%atm%list_sea%ncount(jb)                   ! number of ocean points (open water & sea ice)
!           jc = ext_data%atm%list_sea%idx(ic,jb)                       ! index list of ocean points
!          !ECHAM  prm_field(jg)%ocv(n,i_blk) = buffer(nn+jc,1)
!          !NWP: no ocean current in turbtran yet
!         END DO
! 
!       ENDDO
! !ICON_OMP_END_PARALLEL_DO
! 
!       !ECHAM CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%ocv(:,:))
!       !NWP: no ocean current in turbtran yet
!     END IF


!------------------------------------------------
!  Receive sea ice bundle
!    field_id(9) represents "ocean_sea_ice_bundle" 
!    - ice thickness, snow thickness, ice concentration
!------------------------------------------------

    IF (ltimer) CALL timer_start(timer_coupling_get)

    no_arr = 3
    CALL yac_fget ( field_id(9), nbr_hor_cells, no_arr, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    IF ( info > COUPLING .AND. info < OUT_OF_BOUND ) &
         & CALL message('nwp_couple_ocean', 'YAC says it is get for restart - id=9, sea ice')
    IF ( info == OUT_OF_BOUND ) &
         & CALL warning('nwp_couple_ocean', 'YAC says fget called after end of run - id=9, sea ice')

    IF (ltimer) CALL timer_stop(timer_coupling_get)

    IF ( info > 0 .AND. info < 7 ) THEN

! --- Here we loop only over ocean points, because fr_seaice and h_ice are also used for oceans and lakes.

!ICON_OMP_PARALLEL_DO PRIVATE(jb, ic, jc, nn) ICON_OMP_RUNTIME_SCHEDULE
      DO jb = i_startblk, i_endblk
        nn = (jb-1)*nproma                                            ! translation to 1-d buffer fields

        DO ic = 1, ext_data%atm%list_sea%ncount(jb)                   ! number of ocean points (open water & sea ice)
          jc = ext_data%atm%list_sea%idx(ic,jb)                       ! index list of ocean points
          lnd_diag%fr_seaice    (jc,jb) = buffer(nn+jc,3) 
          wtr_prog_new%h_ice    (jc,jb) = buffer(nn+jc,1)             ! overwrite new ice thickness only
! --- Limiting h_ice from ocean to hice_max used in atmospheric sea-ice.  Note that the ocean uses 
!     seaice_limit*dzlev_m(1), a fractional thickness of the ocean top layer.  The user is responsible
! --- for consistency of these three namelist parameters.  There is no ocean/atmo consistency check in ICON. 
          wtr_prog_new%h_ice    (jc,jb) = MIN (wtr_prog_new%h_ice(jc,jb), hice_max-csmall)
!         wtr_prog_now%h_ice    (jc,jb) = buffer(nn+jc,1)
!         wtr_prog_now%h_snow_si(jc,jb) = buffer(nn+jc,2)             ! Dmitrii's seaice doesn't do snow.
!         wtr_prog_new%h_snow_si(jc,jb) = buffer(nn+jc,2)             ! ...
        END DO

      ENDDO
!ICON_OMP_END_PARALLEL_DO

      CALL sync_patch_array(sync_c, p_patch, lnd_diag%fr_seaice    (:,:) ) 
      CALL sync_patch_array(sync_c, p_patch, wtr_prog_new%h_ice    (:,:) )
!     CALL sync_patch_array(sync_c, p_patch, wtr_prog_now%h_ice    (:,:) )
!     CALL sync_patch_array(sync_c, p_patch, wtr_prog_now%h_snow_si(:,:) )
!     CALL sync_patch_array(sync_c, p_patch, wtr_prog_new%h_snow_si(:,:) )

    END IF


!------------------------------------------------
!  Receive co2 flux
!    field_id(12) represents "co2_flux" 
!    - ocean co2 flux
!------------------------------------------------

    IF (ccycle_config(jg)%iccycle /= 0) THEN

      IF (ltimer) CALL timer_start(timer_coupling_get)

      buffer(:,:) = 0.0_wp ! needs to be checked if this is necessary

      CALL yac_fget ( field_id(12), nbr_hor_cells, 1, buffer(1:nbr_hor_cells,1:1), info, ierror )
      IF ( info > COUPLING .AND. info < OUT_OF_BOUND ) &
           & CALL message('nwp_couple_ocean', 'YAC says it is get for restart - id=12, CO2 flux')
      IF ( info == OUT_OF_BOUND )                      &
           & CALL warning('nwp_couple_ocean', 'YAC says fget called after end of run - id=12, CO2 flux')

      IF (ltimer) CALL timer_stop(timer_coupling_get)

!--- no prognostic CO2 in NWP physics
!       IF ( info > 0 .AND. info < 7 ) THEN
! 
! !ICON_OMP_PARALLEL_DO PRIVATE(jb, jc, nn, nlen) ICON_OMP_RUNTIME_SCHEDULE
!         DO jb = 1, p_patch%nblks_c
!           nn = (jb-1)*nproma
!           IF (jb /= p_patch%nblks_c) THEN
!             nlen = nproma
!           ELSE
!             nlen = p_patch%npromz_c
!           END IF
!           DO jc = 1, nlen
!             IF ( nn+jc > nbr_inner_cells ) THEN
!               !ECHAM  prm_field(jg)%co2_flux_tile(n,i_blk,iwtr) = 0.0_wp
!               !NWP: no prognostic CO2 in NWP physics
!             ELSE
!               !ECHAM  prm_field(jg)%co2_flux_tile(n,i_blk,iwtr) = buffer(nn+jc,1)
!               !NWP: no prognostic CO2 in NWP physics
!             ENDIF
!           ENDDO
!         ENDDO
! !ICON_OMP_END_PARALLEL_DO
! 
!         !ECHAM  CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%co2_flux_tile(:,:,iwtr))
!         !NWP: no prognostic CO2 in NWP physics
!       ENDIF

    END IF


!---------DEBUG DIAGNOSTICS-------------------------------------------

    ! calculations for debug print output for namelist debug-values >0 only
    IF ( idbg_mxmn >= 1 .OR. idbg_val >=1 ) THEN

      ! u/v-stress on ice and water
      CALL dbg_print('NWPOce: u_stress wtr',prm_diag%umfl_s_t(:,:,isub_water) ,str_module,3,in_subset=p_patch%cells%owned)
      CALL dbg_print('NWPOce: u_stress ice',prm_diag%umfl_s_t(:,:,isub_seaice),str_module,3,in_subset=p_patch%cells%owned)
      CALL dbg_print('NWPOce: v_stress wtr',prm_diag%vmfl_s_t(:,:,isub_water) ,str_module,4,in_subset=p_patch%cells%owned)
      CALL dbg_print('NWPOce: v_stress ice',prm_diag%vmfl_s_t(:,:,isub_seaice),str_module,4,in_subset=p_patch%cells%owned)

      ! rain, snow, evaporation
!      scr(:,:) = prm_diag%rain_gsp_rate(:,:) + prm_diag%rain_con_rate(:,:)
!      CALL dbg_print('NWPOce: total rain  ',scr,str_module,3,in_subset=p_patch%cells%owned)
!      scr(:,:) = prm_diag%snow_gsp_rate(:,:) + prm_diag%snow_con_rate(:,:)
!      CALL dbg_print('NWPOce: total snow  ',scr,str_module,4,in_subset=p_patch%cells%owned)
!      scr(:,:) = prm_diag%hail_gsp_rate   (:,:) &
!        &      + prm_diag%graupel_gsp_rate(:,:) &
!        &      + prm_diag%ice_gsp_rate    (:,:)
!      CALL dbg_print('NWPOce: ice, hail,..',scr,str_module,4,in_subset=p_patch%cells%owned)
      scr(:,:) = prm_diag%qhfl_s_t(:,:,isub_water)  * ext_data%atm%frac_t(:,:,isub_water) + &
        &        prm_diag%qhfl_s_t(:,:,isub_seaice) * ext_data%atm%frac_t(:,:,isub_seaice)
      CALL dbg_print('NWPOce: evaporation ',scr,str_module,4,in_subset=p_patch%cells%owned)

      ! total: short wave, long wave, sensible, latent heat flux sent
      scr(:,:) = prm_diag%swflxsfc_t(:,:,isub_water) &
        &      + prm_diag%lwflxsfc_t(:,:,isub_water) &
        &      + prm_diag%shfl_s_t  (:,:,isub_water) &
        &      + prm_diag%lhfl_s_t  (:,:,isub_water)
      CALL dbg_print('NWPOce: totalhfx.wtr',scr,str_module,2,in_subset=p_patch%cells%owned)
      scr(:,:) = prm_diag%swflxsfc_t(:,:,isub_water)
      CALL dbg_print('NWPOce: swflxsfc.wtr',scr,str_module,3,in_subset=p_patch%cells%owned)
      scr(:,:) = prm_diag%lwflxsfc_t(:,:,isub_water)
      CALL dbg_print('NWPOce: lwflxsfc.wtr',scr,str_module,4,in_subset=p_patch%cells%owned)
      scr(:,:) =  prm_diag%shfl_s_t (:,:,isub_water)
      CALL dbg_print('NWPOce: shflx.wtr   ',scr,str_module,4,in_subset=p_patch%cells%owned)
      scr(:,:) =  prm_diag%lhfl_s_t (:,:,isub_water)
      CALL dbg_print('NWPOce: lhflx.wtr   ',scr,str_module,4,in_subset=p_patch%cells%owned)

      ! Qtop and Qbot, windspeed sent
      scr(:,:) = prm_diag%shfl_s_t(:,:,isub_seaice) + prm_diag%swflxsfc_t(:,:,isub_seaice) &
             & + prm_diag%lhfl_s_t(:,:,isub_seaice) + prm_diag%lwflxsfc_t(:,:,isub_seaice)
      CALL dbg_print('NWPOce: ice-Qtop    ',scr                ,str_module,4,in_subset=p_patch%cells%owned)
      CALL dbg_print('NWPOce: ice-Qbot    ',lnd_diag%condhf_ice,str_module,3,in_subset=p_patch%cells%owned)
      CALL dbg_print('NWPOce: sp_10m      ',prm_diag%sp_10m    ,str_module,3,in_subset=p_patch%cells%owned)

      ! SST, sea ice, ocean velocity received
      CALL dbg_print('NWPOce: t_seasfc    ',lnd_diag%t_seasfc(:,:)              ,str_module,2,in_subset=p_patch%cells%owned)
      CALL dbg_print('NWPOce: new%h_ice   ',wtr_prog_new%h_ice(:,:)             ,str_module,4,in_subset=p_patch%cells%owned)
      CALL dbg_print('NWPOce: now%h_ice   ',wtr_prog_now%h_ice(:,:)             ,str_module,4,in_subset=p_patch%cells%owned)
!     CALL dbg_print('NWPOce: h_snow_si   ',wtr_prog_now%h_snow_si(:,:)         ,str_module,4,in_subset=p_patch%cells%owned)
!     CALL dbg_print('NWPOce: siced       ',prm_field(jg)%siced       ,str_module,3,in_subset=p_patch%cells%owned)
!     CALL dbg_print('NWPOce: seaice      ',prm_field(jg)%seaice      ,str_module,4,in_subset=p_patch%cells%owned)
!     CALL dbg_print('NWPOce: ocu         ',prm_field(jg)%ocu         ,str_module,4,in_subset=p_patch%cells%owned)
!     CALL dbg_print('NWPOce: ocv         ',prm_field(jg)%ocv         ,str_module,4,in_subset=p_patch%cells%owned)

      !error?
      !CALL dbg_print('NWPOce: ts_tile.iwtr:iwtr',prm_field(jg)%ts_tile(:,:,iwtr:iwtr),str_module,2,in_subset=p_patch%cells%owned)

      ! Fraction of tiles:
      CALL dbg_print('NWPOce: lnd%fr_seaic',lnd_diag%fr_seaice(:,:)   ,str_module,4,in_subset=p_patch%cells%owned)
      scr(:,:) = ext_data%atm%frac_t(:,:,isub_seaice)
      CALL dbg_print('NWPOce: ext%frac.si  ',scr                      ,str_module,3,in_subset=p_patch%cells%owned)
      scr(:,:) = ext_data%atm%frac_t(:,:,isub_water)
      CALL dbg_print('NWPOce: ext%frac.wtr ',scr                      ,str_module,3,in_subset=p_patch%cells%owned)
      scr(:,:) = 1.0_wp - ext_data%atm%frac_t(:,:,isub_water) - ext_data%atm%frac_t(:,:,isub_seaice)
      CALL dbg_print('NWPOce: ext%frac.land',scr                      ,str_module,4,in_subset=p_patch%cells%owned)
      CALL dbg_print('NWPOce: frac_oce     ',frac_oce                 ,str_module,3,in_subset=p_patch%cells%owned)
      CALL dbg_print('NWPOce: ext:frac.lake',ext_data%atm%frac_t(:,:,isub_lake),str_module,4,in_subset=p_patch%cells%owned)
    ENDIF

!---------------------------------------------------------------------

    DEALLOCATE(buffer)

!YAC_coupling
#endif

  END SUBROUTINE nwp_couple_ocean

  
END MODULE mo_nwp_ocean_interface
