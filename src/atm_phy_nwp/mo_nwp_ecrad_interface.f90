!>
!! This module is the interface between ICON:nwp_radiation to the radiation scheme ecRad
!!
!! - There are two interfaces within this module: nwp_ecrad_radiation and
!!   nwp_ecrad_radiation_reduced. The former provides the interface to ecrad on the full
!!   ICON dynamics grid. The latter one provides the interface on a grid with a reduced
!!   spatial resolution.
!! - The decision which of the two interface routines is used, is done via the namelist
!!   switch lredgrid_phys. Based on the value of lredgrid_phys, the correct interface
!!   is called by mo_nwp_rad_interface:nwp_radiation.
!! - The interfaces have to fill the different ecRad input types (ecrad_aerosol,
!!   ecrad_single_level, ecrad_thermodynamics, ecrad_gas, ecrad_cloud) with data from
!!   ICON variables. Then, the ecRad radiation code is called. At the end, the fluxes
!!   calculated by ecRad and stored in the ecrad_flux structure are copied to ICON variables.
!! - The difference between nwp_ecrad_radiation and nwp_ecrad_radiation_reduced is mostly
!!   an upscaling at the beginning and a downscaling at the end of the interface.
!! - The transfer of data from ICON to ecRad and vice versa is performed within
!!   routines from mo_nwp_ecrad_utilities and mo_nwp_ecrad_prep_aerosol, independent of
!!   the choice to use a reduced radiation grid or not.
!!
!! @author Daniel Rieger, Deutscher Wetterdienst, Offenbach
!!
!! @par Revision History
!! Initial release by Daniel Rieger, Deutscher Wetterdienst, Offenbach (2019-01-31)
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

MODULE mo_nwp_ecrad_interface

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: finish, message
  USE mo_math_constants,         ONLY: pi
  USE mo_model_domain,           ONLY: t_patch, p_patch_local_parent
  USE mo_impl_constants,         ONLY: min_rlcell_int
  USE mo_impl_constants_grf,     ONLY: grf_bdywidth_c, grf_ovlparea_start_c, grf_fbk_start_c
  USE mo_fortran_tools,          ONLY: init, assert_acc_device_only
  USE mo_parallel_config,        ONLY: nproma, nproma_sub, nblocks_sub
  USE mo_loopindices,            ONLY: get_indices_c
  USE mo_grid_config,            ONLY: l_limited_area, nexlevs_rrg_vnest
  USE mo_ext_data_types,         ONLY: t_external_data
  USE mo_nwp_lnd_types,          ONLY: t_lnd_prog
  USE mo_nonhydro_types,         ONLY: t_nh_prog, t_nh_diag
  USE mo_nwp_phy_types,          ONLY: t_nwp_phy_diag
  USE mo_physical_constants,     ONLY: rhoh2o
  USE mo_run_config,             ONLY: msg_level, iqv, iqi, iqc, iqr, iqs, iqg
  USE mo_atm_phy_nwp_config,     ONLY: atm_phy_nwp_config
  USE mo_radiation_config,       ONLY: irad_aero, ssi_radt,                                   &
                                   &   iRadAeroNone, iRadAeroConst, iRadAeroTegen,            &
                                   &   iRadAeroART, iRadAeroConstKinne, iRadAeroKinne,        &
                                   &   iRadAeroVolc, iRadAeroKinneVolc,  iRadAeroKinneVolcSP, &
                                   &   iRadAeroKinneSP
  USE mo_phys_nest_utilities,    ONLY: t_upscale_fields, upscale_rad_input, downscale_rad_output
  USE mtime,                     ONLY: datetime
#ifdef __ECRAD
  USE mo_ecrad,                  ONLY: ecrad, ecrad_ssi_default,                 &
                                   &   t_ecrad_conf, t_ecrad_aerosol_type,       &
                                   &   t_ecrad_single_level_type,                &
                                   &   t_ecrad_thermodynamics_type,              &
                                   &   t_ecrad_gas_type, t_ecrad_flux_type,      &
                                   &   t_ecrad_cloud_type, t_opt_ptrs
  USE mo_nwp_ecrad_prep_aerosol, ONLY: nwp_ecrad_prep_aerosol
  USE mo_nwp_ecrad_utilities,    ONLY: ecrad_set_single_level,                   &
                                   &   ecrad_set_thermodynamics,                 &
                                   &   ecrad_set_clouds,                         &
                                   &   ecrad_set_gas,                            &
                                   &   ecrad_store_fluxes, add_3D_diffuse_rad,   &
                                   &   get_indices_rad_subblock
#ifndef __ECRAD_ACC
  USE mo_nwp_ecrad_utilities,    ONLY: ecrad_acc_allocation,                   &
                                   &   ecrad_acc_deallocation,                 &
                                   &   update_host_pre_ecrad,                  &
                                   &   update_device_post_ecrad
#endif
#endif


  IMPLICIT NONE

  PRIVATE
#ifdef __ECRAD
  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_nwp_ecrad_interface'


  PUBLIC :: nwp_ecrad_radiation
  PUBLIC :: nwp_ecrad_radiation_reduced


CONTAINS


  !---------------------------------------------------------------------------------------
  !>
  !! SUBROUTINE nwp_ecrad_radiation:
  !! Interface to ecRad on full grid. This routine
  !!  ... allocates the ecRad data types
  !!  ... fills the ecRad data types with current atmospheric and external data
  !!  ... saves the output to ICON physics data structure
  !!
  !! @par Revision History
  !! Initial release by Daniel Rieger, Deutscher Wetterdienst, Offenbach (2019-01-31)
  !! Add separate (diffuse) near-IR, visible, and PAR fluxes by Roland Wirth, Deutscher
  !!     Wetterdienst, Offenbach (2021-09-15)
  !!
  SUBROUTINE nwp_ecrad_radiation ( current_datetime, pt_patch, ext_data,      &
    &  zaeq1, zaeq2, zaeq3, zaeq4, zaeq5, od_lw, od_sw, ssa_sw,               &
    &  g_sw, pt_diag, prm_diag, pt_prog, lnd_prog, zsct, ecrad_conf, lacc )

    CHARACTER(len=*), PARAMETER:: routine = modname//'::nwp_ecrad_radiation'

    TYPE(datetime), POINTER, INTENT(in)    :: current_datetime !< Current date and time

    TYPE(t_patch), TARGET,   INTENT(in)    :: pt_patch         !< Current domain info
    TYPE(t_external_data),   INTENT(in)    :: ext_data         !< External data container

    REAL(wp),                INTENT(in)    ::             &
      & zaeq1(nproma,pt_patch%nlev,pt_patch%nblks_c),     & !< Climatological aerosol (Tegen)
      & zaeq2(nproma,pt_patch%nlev,pt_patch%nblks_c),     & !< Climatological aerosol (Tegen)
      & zaeq3(nproma,pt_patch%nlev,pt_patch%nblks_c),     & !< Climatological aerosol (Tegen)
      & zaeq4(nproma,pt_patch%nlev,pt_patch%nblks_c),     & !< Climatological aerosol (Tegen)
      & zaeq5(nproma,pt_patch%nlev,pt_patch%nblks_c)        !< Climatological aerosol (Tegen)
    REAL(wp), TARGET, ALLOCATABLE, INTENT(in) ::          &
      & od_lw (:,:,:,:)                             ,     & !< LW aerosol optical thickness
      & od_sw (:,:,:,:)                             ,     & !< SW aerosol optical thickness
      & g_sw  (:,:,:,:)                             ,     & !< SW aerosol asymmetry factor
      & ssa_sw(:,:,:,:)                                     !< SW aerosol single scattering albedo

    TYPE(t_nh_diag), TARGET, INTENT(in)         :: pt_diag       !< ICON diagnostic variables
    TYPE(t_nwp_phy_diag), TARGET, INTENT(inout) :: prm_diag      !< ICON physics diagnostics
    TYPE(t_nh_prog), TARGET, INTENT(in)         :: pt_prog        !< ICON dyn prog vars
    TYPE(t_lnd_prog),        INTENT(inout)      :: lnd_prog      !< ICON prognostic land state
    
    REAL(wp),                INTENT(in)         ::   zsct        !< Time-dependent solar constant

    TYPE(t_ecrad_conf),      INTENT(in)         :: ecrad_conf    !< ecRad configuration object
    LOGICAL,                 INTENT(IN), OPTIONAL:: lacc
! Local variables
    TYPE(t_ecrad_aerosol_type)        :: &
      &  ecrad_aerosol                     !< ecRad aerosol information (input)
    TYPE(t_ecrad_single_level_type)   :: &
      &  ecrad_single_level                !< ecRad single level information (input)
    TYPE(t_ecrad_thermodynamics_type) :: &
      &  ecrad_thermodynamics              !< ecRad thermodynamics information (input)
    TYPE(t_ecrad_gas_type)            :: &
      &  ecrad_gas                         !< ecRad gas information (input)
    TYPE(t_ecrad_cloud_type)          :: &
      &  ecrad_cloud                       !< ecRad cloud information (input)
    TYPE(t_ecrad_flux_type)           :: &
      &  ecrad_flux                        !< ecRad flux information (output)
    TYPE(t_opt_ptrs),ALLOCATABLE      :: &
      &  opt_ptrs_lw(:), opt_ptrs_sw(:)    !< Contains pointers to aerosol optical properties
    REAL(wp)                 :: &
      &  fact_reffc               !< Factor in the calculation of cloud droplet effective radius
    INTEGER                  :: &
      &  jc, jb, jk, jw,        & !< Loop indices
      &  jg,                    & !< Domain index
      &  nlev, nlevp1,          & !< Number of vertical levels (full, half)
      &  rl_start, rl_end,      & !<
      &  i_startblk, i_endblk,  & !< blocks
      &  i_startidx, i_endidx,  & !< slices
      &  jb_rad,                & !< index of subblock
      &  i_startidx_sub,        & !< start index of subblock in nproma block
      &  i_endidx_sub,          & !< end index of subblock in nproma block
      &  i_startidx_rad,        & !< start index of subblock in nproma_sub
      &  i_endidx_rad,          & !< end index of subblock in nproma_sub
      &  jcs, jce                 !< raw start and end index of subblock in nproma with boundaries
    LOGICAL, ALLOCATABLE     :: &
      &  cosmu0mask(:)            !< Mask if cosmu0 > 0
    REAL(wp), ALLOCATABLE    :: &
      &  zlwflx_up(:,:),        & !< longwave upward flux
      &  zlwflx_dn(:,:),        & !< longwave downward flux
      &  zswflx_up(:,:),        & !< shortwave upward flux
      &  zswflx_dn(:,:),        & !< shortwave downward flux
      &  zlwflx_up_clr(:,:),    & !< longwave upward clear-sky flux
      &  zlwflx_dn_clr(:,:),    & !< longwave downward clear-sky flux
      &  zswflx_up_clr(:,:),    & !< shortave upward clear-sky flux
      &  zswflx_dn_clr(:,:)       !< shortave downward clear-sky flux
    REAL(wp), DIMENSION(:,:),  POINTER :: &
      &  ptr_acdnc => NULL(),                                                 &
      &  ptr_qr => NULL(),      ptr_qs => NULL(),      ptr_qg => NULL(),      &
      &  ptr_reff_qc => NULL(), ptr_reff_qi => NULL(), ptr_reff_qr => NULL(), &
      &  ptr_reff_qs => NULL(), ptr_reff_qg => NULL()
    REAL(wp), DIMENSION(:),    POINTER :: &
      &  ptr_fr_glac => NULL(), ptr_fr_land => NULL()

    CALL assert_acc_device_only(routine, lacc)

    !$ACC DATA CREATE(ecrad_aerosol, ecrad_cloud, ecrad_flux, ecrad_gas, ecrad_single_level, ecrad_thermodynamics) &
    !$ACC   PRESENT(lnd_prog, prm_diag)

    nlev      = pt_patch%nlev
    nlevp1    = nlev+1
    jg        = pt_patch%id

    fact_reffc = (3.0e-9_wp/(4.0_wp*pi*rhoh2o))**(1.0_wp/3.0_wp)

    IF (msg_level >= 7) &
      &       CALL message(routine, 'ecrad radiation on full grid')

    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = pt_patch%cells%start_block(rl_start)
    i_endblk   = pt_patch%cells%end_block(rl_end)

!$OMP PARALLEL PRIVATE(cosmu0mask, opt_ptrs_lw, opt_ptrs_sw, jw,                     &
!$OMP                  zlwflx_up,     zlwflx_dn,     zswflx_up,     zswflx_dn,       &
!$OMP                  zlwflx_up_clr, zlwflx_dn_clr, zswflx_up_clr, zswflx_dn_clr,   &
!$OMP                  ecrad_aerosol,ecrad_single_level, ecrad_thermodynamics,       &
!$OMP                  ecrad_gas, ecrad_cloud, ecrad_flux)

    ALLOCATE( cosmu0mask   (nproma_sub)     )
    ALLOCATE( zlwflx_up    (nproma_sub,nlevp1), zlwflx_dn    (nproma_sub,nlevp1) )
    ALLOCATE( zswflx_up    (nproma_sub,nlevp1), zswflx_dn    (nproma_sub,nlevp1) )
    ALLOCATE( zlwflx_up_clr(nproma_sub,nlevp1), zlwflx_dn_clr(nproma_sub,nlevp1) )
    ALLOCATE( zswflx_up_clr(nproma_sub,nlevp1), zswflx_dn_clr(nproma_sub,nlevp1) )
    ALLOCATE( opt_ptrs_lw(ecrad_conf%n_bands_lw), opt_ptrs_sw(ecrad_conf%n_bands_sw) )
    !$ACC ENTER DATA CREATE(cosmu0mask, zlwflx_up, zlwflx_dn, zswflx_up, zswflx_dn, zlwflx_up_clr, zlwflx_dn_clr) &
    !$ACC   CREATE(zswflx_up_clr, zswflx_dn_clr) ASYNC(1)

    CALL ecrad_single_level%allocate(nproma_sub, 2, 1, .true.) !< use_sw_albedo_direct, 2 bands
    ecrad_single_level%solar_irradiance = 1._wp            !< Obtain normalized fluxes which corresponds to the
                                                           !< transmissivity needed in the following
    !$ACC UPDATE DEVICE(ecrad_single_level%solar_irradiance) ASYNC(1)

    IF (ecrad_conf%use_spectral_solar_scaling) THEN
      ALLOCATE(ecrad_single_level%spectral_solar_scaling(ecrad_conf%n_bands_sw))
      ecrad_single_level%spectral_solar_scaling = ssi_radt / ecrad_ssi_default
      !$ACC ENTER DATA COPYIN(ecrad_single_level%spectral_solar_scaling) ASYNC(1)
    ENDIF

    CALL ecrad_thermodynamics%allocate(nproma_sub, nlev, use_h2o_sat=.false., rrtm_pass_temppres_fl=.true.)

    CALL ecrad_gas%allocate(nproma_sub, nlev)

    CALL ecrad_cloud%allocate(nproma_sub, nlev)
    ! Currently hardcoded values for FSD
    !$ACC WAIT
    CALL ecrad_cloud%create_fractional_std(nproma_sub, nlev, 1._wp)

    IF ( ecrad_conf%use_aerosols ) THEN
      ! Allocate aerosol container
      CALL ecrad_aerosol%allocate_direct(ecrad_conf, nproma_sub, 1, nlev)
    ENDIF

    CALL ecrad_flux%allocate(ecrad_conf, 1, nproma_sub, nlev)

#ifndef __ECRAD_ACC
    ! The current master of libecrad submodule doesnot support OpenACC.
    ! If a libecrad submodule is used that is ported with OpenACC __ECRAD_ACC
    ! has to be defined to disable this functionality.
    ! Once the libecrad submodule master includes the OpenACC port, everything
    ! that is covered in this #ifndef __ECRAD_ACC can be removed.
    ! In the subroutine ecrad_acc_allocation, the arrays are allocated on the GPU that are allocated on the CPU in:
    !   - ecrad_single_level%allocate
    !   - ecrad_thermodynamics%allocate
    !   - ecrad_gas%allocate
    !   - ecrad_cloud%allocate
    !   - ecrad_cloud%create_fractional_std
    !   - ecrad_aerosol%allocate_direct
    !   - ecrad_flux%allocate
    call ecrad_acc_allocation(ecrad_conf, ecrad_single_level, ecrad_thermodynamics, ecrad_gas, ecrad_cloud, &
      &                       ecrad_aerosol, ecrad_flux)
#endif

!$OMP DO PRIVATE(jb, jc, i_startidx, i_endidx,                   &
!$OMP            jb_rad, jcs, jce, i_startidx_sub, i_endidx_sub, &
!$OMP            i_startidx_rad, i_endidx_rad,                   &
!$OMP            ptr_acdnc, ptr_fr_land, ptr_fr_glac,            &
!$OMP            ptr_reff_qc, ptr_reff_qi, ptr_qr, ptr_reff_qr,  &
!$OMP            ptr_qs, ptr_reff_qs, ptr_qg, ptr_reff_qg),      &
!$OMP ICON_OMP_GUIDED_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, rl_start, rl_end)

      ! It may happen that an MPI patch contains only nest boundary points
      ! In this case, no action is needed
      IF (i_startidx > i_endidx) CYCLE

      DO jb_rad = 1, nblocks_sub
        CALL get_indices_rad_subblock(i_startidx, i_endidx, nproma_sub, jb_rad, jcs, jce, &
          &  i_startidx_rad, i_endidx_rad, i_startidx_sub=i_startidx_sub, i_endidx_sub=i_endidx_sub)

        IF (i_startidx_rad > i_endidx_rad) CYCLE

        NULLIFY(ptr_acdnc,ptr_fr_land,ptr_fr_glac,ptr_reff_qc,ptr_reff_qi,    &
                ptr_qr, ptr_reff_qr, ptr_qs, ptr_reff_qs,ptr_qg, ptr_reff_qg)

        IF (atm_phy_nwp_config(jg)%icpl_rad_reff == 0) THEN  ! Own calculation of reff inside ecrad_set_clouds()
          ptr_acdnc   => prm_diag%acdnc(jcs:jce,:,jb)
          ptr_fr_land => ext_data%atm%fr_land(jcs:jce,jb)
          ptr_fr_glac => ext_data%atm%fr_glac(jcs:jce,jb)
        ELSE
          ptr_reff_qc => prm_diag%reff_qc(jcs:jce,:,jb)
          ptr_reff_qi => prm_diag%reff_qi(jcs:jce,:,jb)
        ENDIF

        IF (atm_phy_nwp_config(jg)%icpl_rad_reff == 2) THEN ! Option to use all hydrometeors reff individually
          ! Set extra hydropmeteors
          ptr_qr => pt_prog%tracer(jcs:jce,:,jb,iqr)
          ptr_qs => pt_prog%tracer(jcs:jce,:,jb,iqs)
          IF (iqg > 0) ptr_qg => pt_prog%tracer(jcs:jce,:,jb,iqg)
          ptr_reff_qr => prm_diag%reff_qr(jcs:jce,:,jb)
          ptr_reff_qs => prm_diag%reff_qs(jcs:jce,:,jb)
          IF (iqg > 0) ptr_reff_qg => prm_diag%reff_qg(jcs:jce,:,jb)
        END IF

        !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
        !$ACC LOOP GANG VECTOR
        do jc = i_startidx_sub, i_endidx_sub
          prm_diag%tsfctrad(jc,jb) = lnd_prog%t_g(jc,jb)
        end do
        !$ACC END PARALLEL

        !$ACC DATA PRESENT(cosmu0mask)
        !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
        !$ACC LOOP GANG VECTOR
        do jc = 1, nproma_sub
          cosmu0mask(jc) = .FALSE.
        end do
        !$ACC END PARALLEL
        !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
        !$ACC LOOP GANG VECTOR
        DO jc = i_startidx_rad, i_endidx_rad
          IF ( prm_diag%cosmu0(jc+nproma_sub*(jb_rad-1),jb) > 0._wp ) THEN
            cosmu0mask(jc) = .TRUE.
          ENDIF
        ENDDO
        !$ACC END PARALLEL
        !$ACC END DATA

! Fill single level configuration type
        !$ACC WAIT
        CALL ecrad_set_single_level(ecrad_single_level, current_datetime, pt_patch%cells%center(jcs:jce,jb),           &
          &                         prm_diag%cosmu0(jcs:jce,jb), prm_diag%tsfctrad(jcs:jce,jb),                        &
          &                         prm_diag%albvisdif(jcs:jce,jb), prm_diag%albnirdif(jcs:jce,jb),                    &
          &                         prm_diag%albvisdir(jcs:jce,jb), prm_diag%albnirdir(jcs:jce,jb),                    &
          &                         prm_diag%lw_emiss(jcs:jce,jb), i_startidx_rad, i_endidx_rad)

! Fill thermodynamics configuration type
        CALL ecrad_set_thermodynamics(ecrad_thermodynamics, pt_diag%temp(jcs:jce,:,jb), pt_diag%pres(jcs:jce,:,jb),    &
          &                           pt_diag%pres_ifc(jcs:jce,:,jb), nlev, nlevp1, i_startidx_rad, i_endidx_rad, lacc=.TRUE.)

! Fill gas configuration type
        CALL ecrad_set_gas(ecrad_gas, ecrad_conf, ext_data%atm%o3(jcs:jce,:,jb), prm_diag%tot_cld(jcs:jce,:,jb,iqv), &
          &                pt_diag%pres(jcs:jce,:,jb), i_startidx_rad, i_endidx_rad, nlev)

! Fill clouds configuration type
        CALL ecrad_set_clouds(ecrad_cloud, ecrad_thermodynamics, prm_diag%tot_cld(jcs:jce,:,jb,iqc),               &
          &                   prm_diag%tot_cld(jcs:jce,:,jb,iqi), prm_diag%clc(jcs:jce,:,jb),                      &
          &                   pt_diag%temp(jcs:jce,:,jb), pt_diag%pres(jcs:jce,:,jb),                              &
          &                   ptr_acdnc, ptr_fr_glac, ptr_fr_land,                                                 &
          &                   ptr_qr, ptr_qs, ptr_qg, ptr_reff_qc, ptr_reff_qi,                                    &
          &                   ptr_reff_qr, ptr_reff_qs, ptr_reff_qg,                                               &
          &                   atm_phy_nwp_config(jg)%icpl_rad_reff,                                                &
          &                   fact_reffc, ecrad_conf%cloud_fraction_threshold, nlev, i_startidx_rad, i_endidx_rad)
        ! $ACC WAIT

! Fill aerosol configuration type
        SELECT CASE (irad_aero)
          CASE(iRadAeroNone)
            ! No aerosol, nothing to do
          CASE(iRadAeroConst)
            !         Arguments can be added to fill ecrad_aerosol with actual values. For the time being,
            !         we stay consistent with RRTM where irad_aero=2 does not add any aerosol
            CALL nwp_ecrad_prep_aerosol(1, nlev, i_startidx_rad, i_endidx_rad, &
              &                         ecrad_conf, ecrad_aerosol, lacc=.TRUE.)
          CASE(iRadAeroTegen)
            ! Fill aerosol configuration type with Tegen aerosol
            CALL nwp_ecrad_prep_aerosol(1, nlev, i_startidx_rad, i_endidx_rad,     &
              &                         zaeq1(jcs:jce,:,jb), zaeq2(jcs:jce,:,jb),  &
              &                         zaeq3(jcs:jce,:,jb), zaeq4(jcs:jce,:,jb),  &
              &                         zaeq5(jcs:jce,:,jb),                       &
              &                         ecrad_conf, ecrad_aerosol, lacc=.TRUE.)
          CASE(iRadAeroART)
            ! Use ART aerosol
            CALL nwp_ecrad_prep_aerosol(1, nlev, i_startidx_rad, i_endidx_rad, jb, jg, nproma,  &
              &                         zaeq1(jcs:jce,:,jb), zaeq2(jcs:jce,:,jb),               &
              &                         zaeq3(jcs:jce,:,jb), zaeq4(jcs:jce,:,jb),               &
              &                         zaeq5(jcs:jce,:,jb),                                    &
              &                         ecrad_conf, ecrad_aerosol, lacc=.TRUE.)

          CASE(iRadAeroConstKinne,iRadAeroKinne,iRadAeroVolc,iRadAeroKinneVolc,iRadAeroKinneVolcSP,iRadAeroKinneSP)
            DO jw = 1, ecrad_conf%n_bands_lw
              opt_ptrs_lw(jw)%ptr_od  => od_lw(jcs:jce,:,jb,jw)
            ENDDO
            DO jw = 1, ecrad_conf%n_bands_sw
              opt_ptrs_sw(jw)%ptr_od  => od_sw(jcs:jce,:,jb,jw)
              opt_ptrs_sw(jw)%ptr_ssa => ssa_sw(jcs:jce,:,jb,jw)
              opt_ptrs_sw(jw)%ptr_g   => g_sw(jcs:jce,:,jb,jw)
            ENDDO
            CALL nwp_ecrad_prep_aerosol(1, nlev, i_startidx_rad, i_endidx_rad,   &
              &                         opt_ptrs_lw, opt_ptrs_sw,                &
              &                         ecrad_conf, ecrad_aerosol, lacc=.TRUE.)
          CASE DEFAULT
            CALL finish(routine, 'irad_aero not valid for ecRad')
        END SELECT

        !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
        !$ACC LOOP GANG VECTOR
        DO jc = 1, nproma_sub
          ecrad_flux%cloud_cover_sw(jc) = 0._wp
          ecrad_flux%cloud_cover_lw(jc) = 0._wp
        END DO
        !$ACC END PARALLEL

!---------------------------------------------------------------------------------------
! Call the radiation scheme ecRad
!---------------------------------------------------------------------------------------
#ifndef __ECRAD_ACC
        ! The current master of libecrad submodule doesnot support OpenACC.
        ! If a libecrad submodule is used that is ported with OpenACC __ECRAD_ACC
        ! has to be defined to disable this functionality.
        ! Once the libecrad submodule master includes the OpenACC port, everything
        ! that is covered in this #ifndef __ECRAD_ACC can be removed.
        ! In the subroutine update_host_pre_ecrad, all data is updated on the
        ! CPU that is needed for the ecrad computation.
        CALL update_host_pre_ecrad(ecrad_conf, ecrad_single_level, ecrad_thermodynamics, ecrad_gas, &
          &                        ecrad_cloud, ecrad_aerosol, ecrad_flux)
#endif
        !$ACC WAIT
        CALL ecrad(nproma_sub, nlev,                        & !< Array and loop bounds (input)
          &        i_startidx_rad, i_endidx_rad,            & !< Array and loop bounds (input)
          &        ecrad_conf,                              & !< General ecRad configuration object (input)
          &        ecrad_single_level,                      & !< ecRad single level configuration object (input)
          &        ecrad_thermodynamics,                    & !< ecRad thermodynamics configuration object (input)
          &        ecrad_gas,                               & !< ecRad gas configuration object (input)
          &        ecrad_cloud,                             & !< ecRad cloud configuration object (input)
          &        ecrad_aerosol,                           & !< ecRad aerosol configuration object (input)
          &        ecrad_flux                               ) !< ecRad fluxes in the longwave BUT flux/solar constant in the shortwave (output)
#ifndef __ECRAD_ACC
        ! The current master of libecrad submodule doesnot support OpenACC.
        ! If a libecrad submodule is used that is ported with OpenACC __ECRAD_ACC
        ! has to be defined to disable this functionality.
        ! Once the libecrad submodule master includes the OpenACC port, everything
        ! that is covered in this #ifndef __ECRAD_ACC can be removed.
        ! In the subroutine update_device_post_ecrad, all data is updated on the
        ! GPU that result from the ecrad computation.
        CALL update_device_post_ecrad(ecrad_conf, ecrad_flux)
#endif

!---------------------------------------------------------------------------------------

! Update ICON variables with fluxes from ecRad
        CALL ecrad_store_fluxes(jg, ecrad_flux, prm_diag%cosmu0(jcs:jce,jb), prm_diag%trsolall       (jcs:jce,:,jb),  &
          &                     prm_diag%trsol_up_toa          (jcs:jce,jb), prm_diag%trsol_up_sfc     (jcs:jce,jb),  &
          &                     prm_diag%trsol_nir_sfc         (jcs:jce,jb), prm_diag%trsol_vis_sfc    (jcs:jce,jb),  &
          &                     prm_diag%trsol_par_sfc         (jcs:jce,jb), prm_diag%fr_nir_sfc_diff  (jcs:jce,jb),  &
          &                     prm_diag%fr_vis_sfc_diff       (jcs:jce,jb), prm_diag%fr_par_sfc_diff  (jcs:jce,jb),  &
          &                     prm_diag%trsol_dn_sfc_diff     (jcs:jce,jb),                                          &
          &                     prm_diag%trsolclr_sfc          (jcs:jce,jb), prm_diag%lwflxall       (jcs:jce,:,jb),  &
          &                     prm_diag%lwflx_up_sfc_rs       (jcs:jce,jb), prm_diag%lwflxclr_sfc     (jcs:jce,jb),  &
          &                     zlwflx_up    (:,:), zlwflx_dn       (:,:), zswflx_up    (:,:), zswflx_dn    (:,:),    &
          &                     zlwflx_up_clr(:,:), zlwflx_dn_clr   (:,:), zswflx_up_clr(:,:), zswflx_dn_clr(:,:),    &
          &                     cosmu0mask, zsct, i_startidx_rad, i_endidx_rad, nlevp1)

        IF (atm_phy_nwp_config(jg)%l_3d_rad_fluxes) THEN
          !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) PRESENT(zlwflx_up, zlwflx_dn, zswflx_up, zswflx_dn, zlwflx_up_clr) &
          !$ACC   PRESENT(zlwflx_dn_clr, zswflx_up_clr, zswflx_dn_clr)
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          DO jk = 1, nlevp1
            DO jc = i_startidx_rad, i_endidx_rad
              prm_diag%lwflx_up    (jc+nproma_sub*(jb_rad-1),jk,jb) = zlwflx_up    (jc,jk)
              prm_diag%lwflx_dn    (jc+nproma_sub*(jb_rad-1),jk,jb) = zlwflx_dn    (jc,jk)
              prm_diag%swflx_up    (jc+nproma_sub*(jb_rad-1),jk,jb) = zswflx_up    (jc,jk)
              prm_diag%swflx_dn    (jc+nproma_sub*(jb_rad-1),jk,jb) = zswflx_dn    (jc,jk)
              prm_diag%lwflx_up_clr(jc+nproma_sub*(jb_rad-1),jk,jb) = zlwflx_up_clr(jc,jk)
              prm_diag%lwflx_dn_clr(jc+nproma_sub*(jb_rad-1),jk,jb) = zlwflx_dn_clr(jc,jk)
              prm_diag%swflx_up_clr(jc+nproma_sub*(jb_rad-1),jk,jb) = zswflx_up_clr(jc,jk)
              prm_diag%swflx_dn_clr(jc+nproma_sub*(jb_rad-1),jk,jb) = zswflx_dn_clr(jc,jk)
            ENDDO
          ENDDO
          !$ACC END PARALLEL
        ENDIF

        ! Add 3D contribution to diffuse radiation
        !$ACC WAIT
        CALL add_3D_diffuse_rad(ecrad_flux, prm_diag%clc(jcs:jce,:,jb), pt_diag%pres(jcs:jce,:,jb),              &
          &                     pt_diag%temp(jcs:jce,:,jb), prm_diag%cosmu0(jcs:jce,jb),                         &
          &                     prm_diag%fr_nir_sfc_diff(jcs:jce,jb), prm_diag%fr_vis_sfc_diff(jcs:jce,jb),      &
          &                     prm_diag%fr_par_sfc_diff(jcs:jce,jb), prm_diag%trsol_dn_sfc_diff(jcs:jce,jb),    &
          &                     i_startidx_rad, i_endidx_rad, nlev)

      ENDDO ! jb_rad
    ENDDO ! jb
!$OMP END DO

! CLEANUP
    !$ACC WAIT
#ifndef __ECRAD_ACC
    ! The current master of libecrad submodule doesnot support OpenACC.
    ! If a libecrad submodule is used that is ported with OpenACC __ECRAD_ACC
    ! has to be defined to disable this functionality.
    ! Once the libecrad submodule master includes the OpenACC port, everything
    ! that is covered in this #ifndef __ECRAD_ACC can be removed.
    ! In the subroutine ecrad_acc_deallocation the derived ecrad types are dealocated on the GPU.
    CALL ecrad_acc_deallocation(ecrad_conf, ecrad_single_level, ecrad_thermodynamics, ecrad_gas, ecrad_cloud, &
      &                         ecrad_aerosol, ecrad_flux)
#endif
    CALL ecrad_single_level%deallocate()
    CALL ecrad_thermodynamics%deallocate()
    CALL ecrad_gas%deallocate()
    CALL ecrad_cloud%deallocate()
    IF ( ecrad_conf%use_aerosols ) CALL ecrad_aerosol%deallocate()
    CALL ecrad_flux%deallocate()
    !$ACC EXIT DATA DELETE(cosmu0mask, zlwflx_up, zlwflx_dn, zswflx_up, zswflx_dn, zlwflx_up_clr, zlwflx_dn_clr) &
    !$ACC   DELETE(zswflx_up_clr, zswflx_dn_clr)
    DEALLOCATE( cosmu0mask )
    DEALLOCATE( zlwflx_up,     zlwflx_dn,     zswflx_up,    zswflx_dn     )
    DEALLOCATE( zlwflx_up_clr, zlwflx_dn_clr, zswflx_up_clr,zswflx_dn_clr )
    DO jw = 1, ecrad_conf%n_bands_lw
      CALL opt_ptrs_lw(jw)%finalize()
    ENDDO
    DO jw = 1, ecrad_conf%n_bands_sw
      CALL opt_ptrs_sw(jw)%finalize()
    ENDDO
    DEALLOCATE( opt_ptrs_lw, opt_ptrs_sw )
!$OMP END PARALLEL

    !$ACC END DATA

  END SUBROUTINE nwp_ecrad_radiation
  !---------------------------------------------------------------------------------------


  !---------------------------------------------------------------------------------------
  !>
  !! SUBROUTINE nwp_ecrad_radiation_reduced:
  !! Interface to ecRad on reduced radiation grid. This routine
  !!  ... allocates the ecRad data types
  !!  ... fills the ecRad data types with current atmospheric and external data
  !!  ... saves the output to ICON physics data structure
  !!
  !! @par Revision History
  !! Initial release by Daniel Rieger, Deutscher Wetterdienst, Offenbach (2019-01-31)
  !! Add separate (diffuse) near-IR, visible, and PAR fluxes by Roland Wirth, Deutscher
  !!     Wetterdienst, Offenbach (2021-09-15)
  !! Open TODOs: dust_tunefac not considered so far
  !!
  SUBROUTINE nwp_ecrad_radiation_reduced (current_datetime, pt_patch, pt_par_patch, ext_data,  &
    &                                     zaeq1,zaeq2,zaeq3,zaeq4,zaeq5,                       &
    &                                     od_lw, od_sw, ssa_sw, g_sw,                          &
    &                                     pt_diag,prm_diag,pt_prog, lnd_prog, zsct, ecrad_conf, lacc )

    CHARACTER(len=*), PARAMETER :: &
      &  routine = modname//'::nwp_ecrad_radiation_reduced'

    TYPE(datetime), POINTER, INTENT(in)    :: current_datetime !< Current date and time

    TYPE(t_patch), TARGET,   INTENT(in)    :: pt_patch         !< Current domain info
    TYPE(t_patch), TARGET,   INTENT(in)    :: pt_par_patch     !< Parent domain info
    TYPE(t_external_data),   INTENT(in)    :: ext_data         !< External data container

    REAL(wp),                INTENT(in)    ::             &
      & zaeq1(nproma,pt_patch%nlev,pt_patch%nblks_c),     & !< Climatological aerosol (Tegen)
      & zaeq2(nproma,pt_patch%nlev,pt_patch%nblks_c),     & !< Climatological aerosol (Tegen)
      & zaeq3(nproma,pt_patch%nlev,pt_patch%nblks_c),     & !< Climatological aerosol (Tegen)
      & zaeq4(nproma,pt_patch%nlev,pt_patch%nblks_c),     & !< Climatological aerosol (Tegen)
      & zaeq5(nproma,pt_patch%nlev,pt_patch%nblks_c)        !< Climatological aerosol (Tegen)
    REAL(wp), ALLOCATABLE,   INTENT(in)    ::             &
      & od_lw (:,:,:,:)                             ,     & !< LW aerosol optical thickness
      & od_sw (:,:,:,:)                             ,     & !< SW aerosol optical thickness
      & g_sw  (:,:,:,:)                             ,     & !< SW aerosol asymmetry factor
      & ssa_sw(:,:,:,:)                                     !< SW aerosol single scattering albedo

    TYPE(t_nh_diag), TARGET, INTENT(in)    :: pt_diag       !< ICON diagnostic variables
    TYPE(t_nwp_phy_diag),    INTENT(inout) :: prm_diag      !< ICON physics diagnostics
    TYPE(t_nh_prog), TARGET, INTENT(in)    :: pt_prog        !< ICON dyn prog vars
    TYPE(t_lnd_prog),        INTENT(inout) :: lnd_prog      !< ICON prognostic land state
    REAL(wp),                INTENT(in)    :: zsct        !< Time-dependent solar constant

    TYPE(t_ecrad_conf),      INTENT(in)    :: ecrad_conf    !< ecRad configuration object
    LOGICAL, OPTIONAL,       INTENT(in)    :: lacc ! If true, use openacc
! Local variables
    TYPE(t_patch), POINTER            :: &
      &  ptr_pp                            !< Pointer to parent patch of current domain
    TYPE(t_ecrad_aerosol_type)        :: &
      &  ecrad_aerosol                     !< ecRad aerosol information (input)
    TYPE(t_ecrad_single_level_type)   :: &
      &  ecrad_single_level                !< ecRad single level information (input)
    TYPE(t_ecrad_thermodynamics_type) :: &
      &  ecrad_thermodynamics              !< ecRad thermodynamics information (input)
    TYPE(t_ecrad_gas_type)            :: &
      &  ecrad_gas                         !< ecRad gas information (input)
    TYPE(t_ecrad_cloud_type)          :: &
      &  ecrad_cloud                       !< ecRad cloud information (input)
    TYPE(t_ecrad_flux_type)           :: &
      &  ecrad_flux                        !< ecRad flux information (output)
    TYPE(t_upscale_fields)            :: &
      & input_extra_flds,  input_extra_2D, input_extra_reff !< pointer array for input in upscale routine
    REAL(wp)                 :: &
      &  fact_reffc               !< Factor in the calculation of cloud droplet effective radius
    INTEGER                  :: &
      &  nblks_par_c,           & !< nblks for reduced grid (parent domain)
      &  nblks_lp_c,            & !< nblks for reduced grid (local parent)
      &  jb, jc, jw,            & !< loop indices
      &  jg,                    & !< domain id
      &  nlev,                  & !< number of full levels
      &  nlev_rg, nlev_rgp1,    & !< number of full and half levels at reduced grid
      &  rl_start, rl_end,      & !<
      &  i_startblk, i_endblk,  & !< blocks
      &  i_startidx, i_endidx,  & !< slices
      &  np, nl,                & !< dimension variables for allocation (3d fluxes)
      &  jb_rad,                & !< index of subblock
      &  i_startidx_rad,        & !< start index of subblock in nproma_sub
      &  i_endidx_rad,          & !< end index of subblock in nproma_sub
      &  jcs, jce,              & !< raw start and end index of subblock in nproma with boundaries
      &  jnps, jnpe               !< raw start and end index of subblock in nproma with boundaries for potential empty arrays

    ! For radiation on reduced grid
    ! These fields need to be allocatable because they have different dimensions for
    ! the global grid and nested grids, and for runs with/without MPI parallelization
    ! Input fields
    REAL(wp), ALLOCATABLE, TARGET :: &
      &  zrg_cosmu0(:,:),            & !< Cosine of solar zenith angle on reduced grid
      &  zrg_tsfc(:,:),              & !< Surface temperature on reduced grid
      &  zrg_emis_rad(:,:),          & !< Longwave surface emissivity on reduced grid
      &  zrg_albvisdir(:,:),         & !< Surface albedo for visible range (direct) on reduced grid
      &  zrg_albnirdir(:,:),         & !< Surface albedo for near IR range (direct) on reduced grid
      &  zrg_albvisdif(:,:),         & !< Surface albedo for visible range (diffuse) on reduced grid
      &  zrg_albnirdif(:,:),         & !< Surface albedo for near IR range (diffuse) on reduced grid
      &  zrg_pres(:,:,:),            & !< Pressure at full levels
      &  zrg_pres_ifc(:,:,:),        & !< Pressure at half levels
      &  zrg_temp(:,:,:),            & !< Temperature at full levels
      &  zrg_o3(:,:,:),              & !< Ozone mass mixing ratio on reduced grid
      &  zrg_tot_cld(:,:,:,:),       & !< Mass mixing ratio of water vapor, cloud water and cloud ice on reduced grid
      &  zrg_clc(:,:,:),             & !< Cloud cover on reduced grid
      &  zrg_aeq1(:,:,:),            & !< Climatological aerosol on reduced grid
      &  zrg_aeq2(:,:,:),            & !< Climatological aerosol on reduced grid
      &  zrg_aeq3(:,:,:),            & !< Climatological aerosol on reduced grid
      &  zrg_aeq4(:,:,:),            & !< Climatological aerosol on reduced grid
      &  zrg_aeq5(:,:,:),            & !< Climatological aerosol on reduced grid
      &  zrg_trsolall(:,:,:),        & !< solar transmissivity, all sky, net down on reduced grid
      &  zrg_lwflxall(:,:,:),        & !< Terrestrial flux, all sky, net down on reduced grid
      &  zrg_lwflx_up_sfc(:,:),      & !< Longwave upward flux at surface on reduced grid
      &  zrg_trsol_up_toa(:,:),      & !< Upward solar transmissivity at TOA on reduced grid
      &  zrg_trsol_up_sfc(:,:),      & !< Upward solar transmissivity at surface on reduced grid
      &  zrg_trsol_dn_sfc_diff(:,:), & !< Downward diffuse solar transmissivity at surface on reduced grid
      &  zrg_trsol_clr_sfc(:,:),     & !< Clear-sky net transmissvity at surface on reduced grid
      &  zrg_aclcov(:,:),            & !< Cloud cover on reduced grid
      &  zrg_trsol_nir_sfc(:,:),     & !< Near-infrared radiation
      &  zrg_trsol_vis_sfc(:,:),     & !< Visible radiation
      &  zrg_trsol_par_sfc(:,:),     & !< Photosynthetically active radiation
      &  zrg_fr_nir_sfc_diff(:,:),   & !< Diffuse fraction of near-infrared radiation
      &  zrg_fr_vis_sfc_diff(:,:),   & !< Diffuse fraction of visible radiation
      &  zrg_fr_par_sfc_diff(:,:),   & !< Diffuse fraction of photosynthetically active radiation
      &  zrg_lwflx_clr_sfc(:,:),     & !< clear-sky net LW flux at surface
      &  zrg_lwflx_up    (:,:,:),    & !< longwave  3D upward   flux
      &  zrg_lwflx_dn    (:,:,:),    & !< longwave  3D downward flux
      &  zrg_swflx_up    (:,:,:),    & !< shortwave 3D upward   flux
      &  zrg_swflx_dn    (:,:,:),    & !< shortwave 3D downward flux
      &  zrg_lwflx_up_clr(:,:,:),    & !< longwave  3D upward   flux clear-sky
      &  zrg_lwflx_dn_clr(:,:,:),    & !< longwave  3D downward flux clear-sky
      &  zrg_swflx_up_clr(:,:,:),    & !< shortwave 3D upward   flux clear-sky
      &  zrg_swflx_dn_clr(:,:,:),    & !< shortwave 3D downward flux clear-sky
      &  zrg_reff_liq(:,:,:),        & !< Effective radius (m) of liquid phase on reduced grid
      &  zrg_reff_frz(:,:,:),        & !< Effective radius (m) of frozen phase on reduced grid
      &  zrg_extra_flds(:,:,:,:),    & !< Extra fields for the upscaling routine (indices by irg_)
      &  zrg_extra_2D(:,:,:),        & !< Extra 2D fields for the upscaling routine (indices by irg_)
      &  zrg_extra_reff(:,:,:,:)       !< Extra effective radius (indices by irg_)

    ! Indices and pointers of extra (optional) fields that are needed by radiation
    ! and therefore have be aggregated to the radiation grid
    INTEGER :: irg_acdnc, irg_fr_glac, irg_fr_land,  irg_qr, irg_qs, irg_qg,  &
      &        irg_reff_qr, irg_reff_qs, irg_reff_qg
    INTEGER, DIMENSION (ecrad_conf%n_bands_lw) :: irg_od_lw
    INTEGER, DIMENSION (ecrad_conf%n_bands_sw) :: irg_od_sw, irg_ssa_sw, irg_g_sw
    REAL(wp), DIMENSION(:,:),  POINTER :: &
      &  ptr_acdnc => NULL(),                                                 &
      &  ptr_qr => NULL(),      ptr_qs => NULL(),      ptr_qg => NULL(),      &
      &  ptr_reff_qc => NULL(), ptr_reff_qi => NULL(), ptr_reff_qr => NULL(), &
      &  ptr_reff_qs => NULL(), ptr_reff_qg => NULL()

    TYPE(t_opt_ptrs),ALLOCATABLE      :: &
      &  opt_ptrs_lw(:), opt_ptrs_sw(:)    !< Contains pointers to aerosol optical properties

    REAL(wp), DIMENSION(:),    POINTER :: &
      &  ptr_fr_glac => NULL(), ptr_fr_land => NULL()

    ! Some unused variables to be up- and downscaled (to not change the interface to up- and downscale)
    REAL(wp), ALLOCATABLE, TARGET :: &
      &  aclcov(:,:),                & !< Cloud cover
      &  zrg_albdif(:,:),            &
      &  zrg_rtype(:,:),             &
      &  zlp_pres_ifc(:,:,:),        &
      &  zlp_tot_cld(:,:,:,:)
    LOGICAL, ALLOCATABLE          :: &
      &  cosmu0mask(:)                 !< Mask if cosmu0 > 0

    jg         = pt_patch%id
    nlev       = pt_patch%nlev

    CALL assert_acc_device_only(routine, lacc)

    fact_reffc = (3.0e-9_wp/(4.0_wp*pi*rhoh2o))**(1.0_wp/3.0_wp)

    IF (msg_level >= 7) &
      &       CALL message(routine, 'ecrad radiation on reduced grid')

    !$ACC DATA CREATE(ecrad_aerosol, ecrad_single_level, ecrad_thermodynamics) &
    !$ACC   CREATE(ecrad_gas, ecrad_cloud, ecrad_flux) PRESENT(prm_diag, lnd_prog)

    IF (jg == 1 .AND. .NOT. l_limited_area) THEN
      ptr_pp      => pt_par_patch
      nblks_par_c =  pt_par_patch%nblks_c
      nblks_lp_c  =  p_patch_local_parent(jg)%nblks_c
    ELSE ! Nested domain with MPI parallelization
      ptr_pp      => p_patch_local_parent(jg)
      nblks_par_c =  ptr_pp%nblks_c
      nblks_lp_c  =  ptr_pp%nblks_c
    ENDIF

    ! Add extra layer for atmosphere above model top if requested
    IF (atm_phy_nwp_config(jg)%latm_above_top) THEN
      IF (jg == 1 .OR. pt_patch%nshift == 0) THEN
        nlev_rg = nlev + 1
      ELSE ! add a specified number levels up to the top of the parent domain in case of vertical nesting
        nlev_rg = MIN(nlev+nexlevs_rrg_vnest, pt_par_patch%nlev)
      ENDIF
    ELSE
      nlev_rg = nlev
    ENDIF
    nlev_rgp1 = nlev_rg+1

    ! Allocate for reduced radiation grid
    ALLOCATE(zrg_cosmu0           (nproma,nblks_par_c),     &
      &      zrg_tsfc             (nproma,nblks_par_c),     &
      &      zrg_emis_rad         (nproma,nblks_par_c),     &
      &      zrg_albvisdir        (nproma,nblks_par_c),     &
      &      zrg_albnirdir        (nproma,nblks_par_c),     &
      &      zrg_albvisdif        (nproma,nblks_par_c),     &
      &      zrg_albnirdif        (nproma,nblks_par_c),     &
      &      zrg_aclcov           (nproma,nblks_par_c),     &
      &      zrg_lwflx_up_sfc     (nproma,nblks_par_c),     &
      &      zrg_trsol_up_toa     (nproma,nblks_par_c),     &
      &      zrg_trsol_up_sfc     (nproma,nblks_par_c),     &
      &      zrg_trsol_nir_sfc     (nproma,nblks_par_c),    &
      &      zrg_trsol_vis_sfc     (nproma,nblks_par_c),    &
      &      zrg_trsol_par_sfc    (nproma,nblks_par_c),     &
      &      zrg_fr_nir_sfc_diff   (nproma,nblks_par_c),    &
      &      zrg_fr_vis_sfc_diff   (nproma,nblks_par_c),    &
      &      zrg_fr_par_sfc_diff   (nproma,nblks_par_c),    &
      &      zrg_trsol_dn_sfc_diff (nproma,nblks_par_c),    &
      &      zrg_trsol_clr_sfc    (nproma,nblks_par_c),     &
      &      zrg_lwflx_clr_sfc    (nproma,nblks_par_c),     &
      &      aclcov               (nproma,pt_patch%nblks_c))
      !$ACC ENTER DATA CREATE(zrg_cosmu0, zrg_tsfc, zrg_emis_rad, zrg_albvisdir) &
      !$ACC   CREATE(zrg_albnirdir, zrg_albvisdif, zrg_albnirdif, zrg_aclcov) &
      !$ACC   CREATE(zrg_lwflx_up_sfc, zrg_trsol_up_toa, zrg_trsol_up_sfc) &
      !$ACC   CREATE(zrg_trsol_nir_sfc, zrg_trsol_vis_sfc, zrg_trsol_par_sfc) &
      !$ACC   CREATE(zrg_fr_nir_sfc_diff, zrg_fr_vis_sfc_diff, zrg_fr_par_sfc_diff) &
      !$ACC   CREATE(zrg_trsol_dn_sfc_diff, zrg_trsol_clr_sfc) &
      !$ACC   CREATE(zrg_lwflx_clr_sfc, aclcov) ASYNC(1)

    ! Set dimensions for 3D radiative flux variables
    IF (atm_phy_nwp_config(jg)%l_3d_rad_fluxes) THEN
       np = nproma
       nl = nlev_rgp1
    ELSE
       np = 1
       nl = 1
    END IF

    ALLOCATE(zrg_pres_ifc    (nproma,nlev_rgp1,nblks_par_c),&
      &      zrg_lwflxall    (nproma,nlev_rgp1,nblks_par_c),&
      &      zrg_trsolall    (nproma,nlev_rgp1,nblks_par_c),&
      &      zrg_lwflx_up    (np, nl, nblks_par_c),&
      &      zrg_lwflx_dn    (np, nl, nblks_par_c),&
      &      zrg_swflx_up    (np, nl, nblks_par_c),&
      &      zrg_swflx_dn    (np, nl, nblks_par_c),&
      &      zrg_lwflx_up_clr(np, nl, nblks_par_c),&
      &      zrg_lwflx_dn_clr(np, nl, nblks_par_c),&
      &      zrg_swflx_up_clr(np, nl, nblks_par_c),&
      &      zrg_swflx_dn_clr(np, nl, nblks_par_c) )
    !$ACC ENTER DATA CREATE(zrg_pres_ifc, zrg_lwflxall, zrg_trsolall) &
    !$ACC   CREATE(zrg_lwflx_up, zrg_lwflx_dn, zrg_swflx_up, zrg_swflx_dn) &
    !$ACC   CREATE(zrg_lwflx_up_clr, zrg_lwflx_dn_clr, zrg_swflx_up_clr) &
    !$ACC   CREATE(zrg_swflx_dn_clr) ASYNC(1)

    ALLOCATE(zrg_pres     (nproma,nlev_rg  ,nblks_par_c),   &
      &      zrg_temp     (nproma,nlev_rg  ,nblks_par_c),   &
      &      zrg_o3       (nproma,nlev_rg  ,nblks_par_c),   &
      &      zrg_aeq1     (nproma,nlev_rg  ,nblks_par_c),   &
      &      zrg_aeq2     (nproma,nlev_rg  ,nblks_par_c),   &
      &      zrg_aeq3     (nproma,nlev_rg  ,nblks_par_c),   &
      &      zrg_aeq4     (nproma,nlev_rg  ,nblks_par_c),   &
      &      zrg_aeq5     (nproma,nlev_rg  ,nblks_par_c),   &
      &      zrg_clc      (nproma,nlev_rg  ,nblks_par_c))
    !$ACC ENTER DATA CREATE(zrg_pres, zrg_temp, zrg_o3, zrg_aeq1, zrg_aeq2) &
    !$ACC   CREATE(zrg_aeq3, zrg_aeq4, zrg_aeq5, zrg_clc) ASYNC(1)

    IF (atm_phy_nwp_config(jg)%icpl_rad_reff > 0) THEN
      ALLOCATE(zrg_reff_liq (nproma,nlev_rg,nblks_par_c),   &
        &      zrg_reff_frz (nproma,nlev_rg,nblks_par_c))
      !$ACC ENTER DATA CREATE(zrg_reff_liq, zrg_reff_frz) ASYNC(1)
    ENDIF

    ALLOCATE(zrg_tot_cld  (nproma,nlev_rg  ,nblks_par_c,3))
    !$ACC ENTER DATA CREATE(zrg_tot_cld) ASYNC(1)

    ! Unused variables, allocated to not change the upscale_rad_input interface
    ALLOCATE(zrg_albdif   (nproma,nblks_par_c),             &
      &      zrg_rtype    (nproma,nblks_par_c),             &
      &      zlp_pres_ifc (nproma,nlev_rgp1,nblks_lp_c ),   &
      &      zlp_tot_cld  (nproma,nlev_rg  ,nblks_lp_c,3))
    !$ACC ENTER DATA CREATE(zrg_albdif, zrg_rtype, zlp_pres_ifc, zlp_tot_cld) &
    !$ACC   ASYNC(1)


    ! Set indices for extra fields in the upscaling routine
    irg_acdnc     = 0
    irg_fr_land   = 0
    irg_fr_glac   = 0
    irg_qr        = 0
    irg_qs        = 0
    irg_qg        = 0
    irg_reff_qr   = 0
    irg_reff_qs   = 0
    irg_reff_qg   = 0
    irg_od_lw     = 0
    irg_od_sw     = 0
    irg_ssa_sw    = 0
    irg_g_sw      = 0

    CALL input_extra_flds%construct(nlev_rg)  ! Extra fields in upscaling routine: 3D fields with nlev_rg
    CALL input_extra_2D%construct(1)          ! Extra fields in upscaling routine: 2D fields
    CALL input_extra_reff%construct(nlev_rg)  ! Extra fields in upscaling routine: extra Reff

    SELECT CASE (atm_phy_nwp_config(jg)%icpl_rad_reff)
    CASE (0)  ! Own calculation of reff inside ecrad_set_clouds()
      CALL input_extra_flds%assign(prm_diag%acdnc, irg_acdnc)
      CALL input_extra_2D%assign(ext_data%atm%fr_land, irg_fr_land)
      CALL input_extra_2D%assign(ext_data%atm%fr_glac, irg_fr_glac)
    CASE (2) ! Option to use all hydrometeors reff individually
      ! Set extra hydrometeors
      CALL input_extra_flds%assign(pt_prog%tracer(:,:,:,iqr), irg_qr)
      CALL input_extra_flds%assign(pt_prog%tracer(:,:,:,iqs), irg_qs)
      IF (iqg >0) CALL input_extra_flds%assign(pt_prog%tracer(:,:,:,iqg), irg_qg)

      ! Set extra effective radius (in different array due to different interpolation)
      CALL input_extra_reff%assign(prm_diag%reff_qr(:,:,:), irg_reff_qr, assoc_hyd = irg_qr )
      CALL input_extra_reff%assign(prm_diag%reff_qs(:,:,:), irg_reff_qs, assoc_hyd = irg_qs )
      IF (iqg >0) CALL input_extra_reff%assign(prm_diag%reff_qg(:,:,:), irg_reff_qg, assoc_hyd = irg_qg )
    END SELECT

    IF (ANY( irad_aero == (/iRadAeroConstKinne,iRadAeroKinne,iRadAeroVolc,  &
      &                     iRadAeroKinneVolc,iRadAeroKinneVolcSP,iRadAeroKinneSP/) )) THEN
      ! Aerosol extra fields
      DO jw = 1, ecrad_conf%n_bands_lw
        CALL input_extra_flds%assign(od_lw(:,:,:,jw), irg_od_lw(jw))
      ENDDO
      DO jw = 1, ecrad_conf%n_bands_sw
        CALL input_extra_flds%assign(od_sw(:,:,:,jw), irg_od_sw(jw))
        CALL input_extra_flds%assign(ssa_sw(:,:,:,jw), irg_ssa_sw(jw))
        CALL input_extra_flds%assign(g_sw(:,:,:,jw), irg_g_sw(jw))
      ENDDO
    END IF

    !$ACC DATA COPYIN(input_extra_flds, input_extra_2D, input_extra_reff)
    CALL input_extra_flds%acc_attach()
    CALL input_extra_2D%acc_attach()
    CALL input_extra_reff%acc_attach()

    ! Allocate output arrays
    IF ( input_extra_flds%ntot > 0 )  THEN
      ALLOCATE( zrg_extra_flds(nproma,input_extra_flds%nlev_rg,nblks_par_c,input_extra_flds%ntot) )
      !$ACC ENTER DATA CREATE(zrg_extra_flds) ASYNC(1)
    END IF
    IF ( input_extra_2D%ntot > 0 )  THEN
      ALLOCATE( zrg_extra_2D(nproma,nblks_par_c,input_extra_2D%ntot) )
      !$ACC ENTER DATA CREATE(zrg_extra_2D) ASYNC(1)
    END IF
    IF ( input_extra_reff%ntot  > 0 ) THEN
      ALLOCATE( zrg_extra_reff(nproma,input_extra_reff%nlev_rg,nblks_par_c,input_extra_reff%ntot) )
      !$ACC ENTER DATA CREATE(zrg_extra_reff) ASYNC(1)
    END IF


    rl_start = 1 ! SR radiation is not set up to handle boundaries of nested domains
    rl_end   = min_rlcell_int
    i_startblk = pt_patch%cells%start_block(rl_start)
    i_endblk   = pt_patch%cells%end_block(rl_end)

!$OMP PARALLEL

    ! Initialize output fields
    CALL init(zrg_trsolall(:,:,:), 0._wp )
    CALL init(zrg_lwflxall(:,:,:), 0._wp )
    CALL init(zrg_trsol_up_toa(:,:)     )
    CALL init(zrg_trsol_up_sfc(:,:)     )
    CALL init(zrg_trsol_nir_sfc(:,:)     )
    CALL init(zrg_trsol_vis_sfc(:,:)     )
    CALL init(zrg_trsol_par_sfc(:,:)    )
    CALL init(zrg_fr_nir_sfc_diff(:,:)  )
    CALL init(zrg_fr_vis_sfc_diff(:,:)  )
    CALL init(zrg_fr_par_sfc_diff(:,:)  )
    CALL init(zrg_trsol_dn_sfc_diff(:,:) )
    CALL init(zrg_trsol_clr_sfc(:,:)    )
    CALL init(zrg_lwflx_up_sfc(:,:)     )
    CALL init(zrg_lwflx_clr_sfc(:,:)    )
    CALL init(zrg_aclcov(:,:)           )

    IF (atm_phy_nwp_config(jg)%l_3d_rad_fluxes) THEN
      CALL init(zrg_lwflx_up    (:,:,:), 0._wp)
      CALL init(zrg_lwflx_dn    (:,:,:), 0._wp)
      CALL init(zrg_swflx_up    (:,:,:), 0._wp)
      CALL init(zrg_swflx_dn    (:,:,:), 0._wp)
      CALL init(zrg_lwflx_up_clr(:,:,:), 0._wp)
      CALL init(zrg_lwflx_dn_clr(:,:,:), 0._wp)
      CALL init(zrg_swflx_up_clr(:,:,:), 0._wp)
      CALL init(zrg_swflx_dn_clr(:,:,:), 0._wp)
    END IF

!$OMP DO PRIVATE(jb, i_startidx, i_endidx), ICON_OMP_GUIDED_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
        &                       i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(NONE) ASYNC(1)
      DO jc = i_startidx, i_endidx
        prm_diag%tsfctrad(jc,jb) = lnd_prog%t_g(jc,jb)
      ENDDO ! jc
      !$ACC END PARALLEL
    ENDDO ! jb
!$OMP END DO NOWAIT

!$OMP END PARALLEL

! Upscale ICON input fields from full grid to reduced radiation grid

    !$ACC WAIT
    CALL upscale_rad_input(pt_patch%id, pt_par_patch%id,                                 &
      &                    nlev_rg,                                                      &
      &                    prm_diag%lw_emiss, prm_diag%cosmu0,                           &
      &                    prm_diag%albvisdir, prm_diag%albnirdir, prm_diag%albvisdif,   &
      &                    prm_diag%albnirdif, prm_diag%albdif, prm_diag%tsfctrad,       &
      &                    prm_diag%ktype, pt_diag%pres_ifc, pt_diag%pres,               &
      &                    pt_diag%temp, prm_diag%tot_cld, prm_diag%clc,                 &
      &                    ext_data%atm%o3, zaeq1, zaeq2, zaeq3, zaeq4, zaeq5,           &
      &                    zrg_emis_rad,                                                 &
      &                    zrg_cosmu0, zrg_albvisdir, zrg_albnirdir, zrg_albvisdif,      &
      &                    zrg_albnirdif, zrg_albdif, zrg_tsfc, zrg_rtype, zrg_pres_ifc, &
      &                    zrg_pres, zrg_temp,                                           &
      &                    zrg_tot_cld, zrg_clc, zrg_o3,                                 &
      &                    zrg_aeq1, zrg_aeq2, zrg_aeq3, zrg_aeq4, zrg_aeq5,             &
      &                    zlp_pres_ifc, zlp_tot_cld, prm_diag%buffer_rrg,               &
      &                    atm_phy_nwp_config(jg)%icpl_rad_reff,                         &
      &                    prm_diag%reff_qc, prm_diag%reff_qi,                           &
      &                    zrg_reff_liq, zrg_reff_frz,                                   &
      &                    input_extra_flds, zrg_extra_flds,                             &
      &                    input_extra_2D, zrg_extra_2D,                                 &
      &                    input_extra_reff, zrg_extra_reff, lacc=.TRUE.)

! Set indices for reduced grid loop
    IF (jg == 1 .AND. l_limited_area) THEN
      rl_start = grf_fbk_start_c
    ELSE
      rl_start = grf_ovlparea_start_c
    ENDIF
    rl_end     = min_rlcell_int
    i_startblk = ptr_pp%cells%start_block(rl_start)
    i_endblk   = ptr_pp%cells%end_block(rl_end)

!$OMP PARALLEL PRIVATE(cosmu0mask, opt_ptrs_lw, opt_ptrs_sw, jw,                &
!$OMP                  ecrad_aerosol, ecrad_single_level, ecrad_thermodynamics, &
!$OMP                  ecrad_gas, ecrad_cloud, ecrad_flux)

    CALL ecrad_single_level%allocate(nproma_sub, 2, 1, .true.) !< use_sw_albedo_direct, 2 bands
    ecrad_single_level%solar_irradiance = 1._wp            !< Obtain normalized fluxes which corresponds to the
                                                           !< transmissivity needed in the following
    !$ACC UPDATE DEVICE(ecrad_single_level%solar_irradiance) ASYNC(1)

    IF (ecrad_conf%use_spectral_solar_scaling) THEN
      ALLOCATE(ecrad_single_level%spectral_solar_scaling(ecrad_conf%n_bands_sw))
      ecrad_single_level%spectral_solar_scaling = ssi_radt / ecrad_ssi_default
      !$ACC ENTER DATA COPYIN(ecrad_single_level%spectral_solar_scaling) ASYNC(1)
    ENDIF

    CALL ecrad_thermodynamics%allocate(nproma_sub, nlev_rg, use_h2o_sat=.false., rrtm_pass_temppres_fl=.true.)

    CALL ecrad_gas%allocate(nproma_sub, nlev_rg)

    CALL ecrad_cloud%allocate(nproma_sub, nlev_rg)
    ! Currently hardcoded values for FSD
    !$ACC WAIT
    CALL ecrad_cloud%create_fractional_std(nproma_sub, nlev_rg, 1._wp)

    IF ( ecrad_conf%use_aerosols ) THEN
      ! Allocate aerosol container
      CALL ecrad_aerosol%allocate_direct(ecrad_conf, nproma_sub, 1, nlev_rg)
    ENDIF

    CALL ecrad_flux%allocate(ecrad_conf, 1, nproma_sub, nlev_rg)

#ifndef __ECRAD_ACC
    ! The current master of libecrad submodule doesnot support OpenACC.
    ! If a libecrad submodule is used that is ported with OpenACC __ECRAD_ACC
    ! has to be defined to disable this functionality.
    ! Once the libecrad submodule master includes the OpenACC port, everything
    ! that is covered in this #ifndef __ECRAD_ACC can be removed.
    ! In the subroutine ecrad_acc_allocation, the arrays are allocated on the GPU that are allocated on the CPU in:
    !   - ecrad_single_level%allocate
    !   - ecrad_thermodynamics%allocate
    !   - ecrad_gas%allocate
    !   - ecrad_cloud%allocate
    !   - ecrad_cloud%create_fractional_std
    !   - ecrad_aerosol%allocate_direct
    !   - ecrad_flux%allocate
    call ecrad_acc_allocation(ecrad_conf, ecrad_single_level, ecrad_thermodynamics, ecrad_gas, ecrad_cloud, &
      &                       ecrad_aerosol, ecrad_flux)
#endif

    ALLOCATE(cosmu0mask(nproma_sub))
    !$ACC ENTER DATA CREATE(cosmu0mask) ASYNC(1)
    ALLOCATE(opt_ptrs_lw(ecrad_conf%n_bands_lw))
    ALLOCATE(opt_ptrs_sw(ecrad_conf%n_bands_sw))

    !$ACC DATA PRESENT(cosmu0mask, zrg_cosmu0, zrg_tsfc) &
    !$ACC   PRESENT(zrg_albvisdif, zrg_albnirdif, zrg_albvisdir, zrg_albnirdir) &
    !$ACC   PRESENT(zrg_emis_rad, zrg_pres_ifc, zrg_temp, zrg_pres, zrg_o3) &
    !$ACC   PRESENT(zrg_tot_cld, zrg_clc, zrg_aeq1, zrg_aeq2, zrg_aeq3, zrg_aeq4) &
    !$ACC   PRESENT(zrg_aeq5)


!$OMP DO PRIVATE(jb, jc, i_startidx, i_endidx,                  &
!$OMP            jb_rad, jcs, jce, jnps, jnpe,                  &
!$OMP            i_startidx_rad,i_endidx_rad,                   &
!$OMP            ptr_acdnc, ptr_fr_land, ptr_fr_glac,           &
!$OMP            ptr_reff_qc, ptr_reff_qi, ptr_qr, ptr_reff_qr, &
!$OMP            ptr_qs, ptr_reff_qs, ptr_qg, ptr_reff_qg),     &
!$OMP ICON_OMP_GUIDED_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(ptr_pp, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, rl_start, rl_end)

      ! It may happen that an MPI patch contains only nest boundary points
      ! In this case, no action is needed
      IF (i_startidx > i_endidx) CYCLE

      DO jb_rad = 1, nblocks_sub
        CALL get_indices_rad_subblock(i_startidx, i_endidx, nproma_sub, jb_rad, jcs, jce, &
          &  i_startidx_rad, i_endidx_rad, l_3d_rad_fluxes=atm_phy_nwp_config(jg)%l_3d_rad_fluxes, &
          &  jnps=jnps, jnpe=jnpe)

        IF (i_startidx_rad > i_endidx_rad) CYCLE

        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(NONE) ASYNC(1)
        DO jc = 1,nproma_sub
          cosmu0mask(jc) = .FALSE.
        ENDDO
        !$ACC END PARALLEL
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(NONE) ASYNC(1)
        DO jc = i_startidx_rad, i_endidx_rad
          IF ( zrg_cosmu0(jc+nproma_sub*(jb_rad-1),jb) > 0._wp ) THEN
            cosmu0mask(jc) = .TRUE.
          ENDIF
        ENDDO
        !$ACC END PARALLEL

        IF (atm_phy_nwp_config(jg)%icpl_rad_reff > 0) THEN
          ptr_reff_qc => zrg_reff_liq(jcs:jce,:,jb)
          ptr_reff_qi => zrg_reff_frz(jcs:jce,:,jb)
        ENDIF
        IF ( irg_acdnc   > 0 ) ptr_acdnc   => zrg_extra_flds(jcs:jce,:,jb,irg_acdnc)
        IF ( irg_qr      > 0 ) ptr_qr      => zrg_extra_flds(jcs:jce,:,jb,irg_qr)
        IF ( irg_qs      > 0 ) ptr_qs      => zrg_extra_flds(jcs:jce,:,jb,irg_qs)
        IF ( irg_qg      > 0 ) ptr_qg      => zrg_extra_flds(jcs:jce,:,jb,irg_qg)
        IF ( irg_fr_land > 0 ) ptr_fr_land => zrg_extra_2D(jcs:jce,jb,irg_fr_land)
        IF ( irg_fr_glac > 0 ) ptr_fr_glac => zrg_extra_2D(jcs:jce,jb,irg_fr_glac)
        IF ( irg_reff_qr > 0 ) ptr_reff_qr => zrg_extra_reff(jcs:jce,:,jb,irg_reff_qr)
        IF ( irg_reff_qs > 0 ) ptr_reff_qs => zrg_extra_reff(jcs:jce,:,jb,irg_reff_qs)
        IF ( irg_reff_qg > 0 ) ptr_reff_qg => zrg_extra_reff(jcs:jce,:,jb,irg_reff_qg)
        IF ( ALL(irg_od_lw(:)  > 0) ) THEN
          DO jw = 1, ecrad_conf%n_bands_lw
            opt_ptrs_lw(jw)%ptr_od  => zrg_extra_flds(jcs:jce,:,jb,irg_od_lw(jw))
          ENDDO
        ENDIF
        IF ( ALL(irg_od_sw(:)  > 0) ) THEN
          DO jw = 1, ecrad_conf%n_bands_sw
            opt_ptrs_sw(jw)%ptr_od  => zrg_extra_flds(jcs:jce,:,jb,irg_od_sw(jw))
          ENDDO
        ENDIF
        IF ( ALL(irg_ssa_sw(:) > 0) ) THEN
          DO jw = 1, ecrad_conf%n_bands_sw
            opt_ptrs_sw(jw)%ptr_ssa => zrg_extra_flds(jcs:jce,:,jb,irg_ssa_sw(jw))
          ENDDO
        ENDIF
        IF ( ALL(irg_g_sw(:)   > 0) ) THEN
          DO jw = 1, ecrad_conf%n_bands_sw
            opt_ptrs_sw(jw)%ptr_g   => zrg_extra_flds(jcs:jce,:,jb,irg_g_sw(jw))
          ENDDO
        ENDIF

! Fill single level configuration type
        !$ACC WAIT
        CALL ecrad_set_single_level(ecrad_single_level, current_datetime, ptr_pp%cells%center(jcs:jce,jb),            &
          &                         zrg_cosmu0(jcs:jce,jb), zrg_tsfc(jcs:jce,jb), zrg_albvisdif(jcs:jce,jb),          &
          &                         zrg_albnirdif(jcs:jce,jb), zrg_albvisdir(jcs:jce,jb), zrg_albnirdir(jcs:jce,jb),  &
          &                         zrg_emis_rad(jcs:jce,jb), i_startidx_rad, i_endidx_rad)

! Fill thermodynamics configuration type
        CALL ecrad_set_thermodynamics(ecrad_thermodynamics, zrg_temp(jcs:jce,:,jb), zrg_pres(jcs:jce,:,jb),     &
          &                           zrg_pres_ifc(jcs:jce,:,jb), nlev_rg, nlev_rgp1, i_startidx_rad, i_endidx_rad, lacc=.TRUE.)

! Fill gas configuration type
        CALL ecrad_set_gas(ecrad_gas, ecrad_conf, zrg_o3(jcs:jce,:,jb), zrg_tot_cld(jcs:jce,:,jb,iqv), &
          &                zrg_pres(jcs:jce,:,jb), i_startidx_rad, i_endidx_rad, nlev_rg)

! Fill clouds configuration type
        CALL ecrad_set_clouds(ecrad_cloud, ecrad_thermodynamics, zrg_tot_cld(jcs:jce,:,jb,iqc), &
          &                   zrg_tot_cld(jcs:jce,:,jb,iqi), zrg_clc(jcs:jce,:,jb),             &
          &                   zrg_temp(jcs:jce,:,jb), zrg_pres(jcs:jce,:,jb), ptr_acdnc,        &
          &                   ptr_fr_glac, ptr_fr_land,                                         &
          &                   ptr_qr, ptr_qs, ptr_qg, ptr_reff_qc, ptr_reff_qi,                 &
          &                   ptr_reff_qr, ptr_reff_qs, ptr_reff_qg,                            &
          &                   atm_phy_nwp_config(jg)%icpl_rad_reff,                             &
          &                   fact_reffc, ecrad_conf%cloud_fraction_threshold, nlev_rg,         &
          &                   i_startidx_rad, i_endidx_rad)

! Fill aerosol configuration type
        SELECT CASE (irad_aero)
          CASE(iRadAeroNone)
            ! No aerosol, nothing to do
          CASE(iRadAeroConst)
            !         Arguments can be added to fill ecrad_aerosol with actual values. For the time being,
            !         we stay consistent with RRTM where irad_aero=2 does not add any aerosol
            CALL nwp_ecrad_prep_aerosol(1, nlev_rg, i_startidx_rad, i_endidx_rad, &
              &                         ecrad_conf, ecrad_aerosol, lacc=.TRUE.)
          CASE(iRadAeroTegen)
            ! Fill aerosol configuration type with Tegen aerosol
            CALL nwp_ecrad_prep_aerosol(1, nlev_rg, i_startidx_rad, i_endidx_rad,         &
              &                         zrg_aeq1(jcs:jce,:,jb), zrg_aeq2(jcs:jce,:,jb),   &
              &                         zrg_aeq3(jcs:jce,:,jb), zrg_aeq4(jcs:jce,:,jb),   &
              &                         zrg_aeq5(jcs:jce,:,jb),                           &
              &                         ecrad_conf, ecrad_aerosol, lacc=.TRUE.)
          CASE(iRadAeroConstKinne,iRadAeroKinne,iRadAeroVolc,iRadAeroKinneVolc,iRadAeroKinneVolcSP,iRadAeroKinneSP)
#ifdef _OPENACC
            CALL finish(routine, 'irad_aero not valid for OpenACC ecrad')
#endif
            CALL nwp_ecrad_prep_aerosol(1, nlev_rg, i_startidx_rad, i_endidx_rad,         &
              &                         opt_ptrs_lw, opt_ptrs_sw,                         &
              &                         ecrad_conf, ecrad_aerosol, lacc=.TRUE.)
          CASE DEFAULT
            CALL finish(routine, 'irad_aero not valid for ecRad')
        END SELECT

        ! Reset output values
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(NONE) ASYNC(1)
        DO jc = 1, nproma_sub
          ecrad_flux%cloud_cover_sw(jc) = 0._wp
          ecrad_flux%cloud_cover_lw(jc) = 0._wp
        END DO
        !$ACC END PARALLEL

!---------------------------------------------------------------------------------------
! Call the radiation scheme ecRad
!---------------------------------------------------------------------------------------
#ifndef __ECRAD_ACC
        ! The current master of libecrad submodule doesnot support OpenACC.
        ! If a libecrad submodule is used that is ported with OpenACC __ECRAD_ACC
        ! has to be defined to disable this functionality.
        ! Once the libecrad submodule master includes the OpenACC port, everything
        ! that is covered in this #ifndef __ECRAD_ACC can be removed.
        ! In the subroutine update_host_pre_ecrad, all data is updated on the
        ! CPU that is needed for the ecrad computation.
        CALL update_host_pre_ecrad(ecrad_conf, ecrad_single_level, ecrad_thermodynamics, ecrad_gas, &
          &                        ecrad_cloud, ecrad_aerosol, ecrad_flux)
#endif
        !$ACC WAIT
        CALL ecrad(nproma_sub, nlev_rg,                     & !< Array and loop bounds (input)
          &        i_startidx_rad, i_endidx_rad,            & !< Array and loop bounds (input)
          &        ecrad_conf,                              & !< General ecRad configuration object (input)
          &        ecrad_single_level,                      & !< ecRad single level configuration object (input)
          &        ecrad_thermodynamics,                    & !< ecRad thermodynamics configuration object (input)
          &        ecrad_gas,                               & !< ecRad gas configuration object (input)
          &        ecrad_cloud,                             & !< ecRad cloud configuration object (input)
          &        ecrad_aerosol,                           & !< ecRad aerosol configuration object (input)
          &        ecrad_flux                               ) !< ecRad fluxes (output)
#ifndef __ECRAD_ACC
        ! The current master of libecrad submodule doesnot support OpenACC.
        ! If a libecrad submodule is used that is ported with OpenACC __ECRAD_ACC
        ! has to be defined to disable this functionality.
        ! Once the libecrad submodule master includes the OpenACC port, everything
        ! that is covered in this #ifndef __ECRAD_ACC can be removed.
        ! In the subroutine update_device_post_ecrad, all data is updated on the
        ! GPU that result from the ecrad computation.
        CALL update_device_post_ecrad(ecrad_conf, ecrad_flux)
#endif
!---------------------------------------------------------------------------------------

! Update ICON variables with fluxes from ecRad
        CALL ecrad_store_fluxes(jg, ecrad_flux, zrg_cosmu0(jcs:jce,jb), zrg_trsolall       (jcs:jce,:,jb),    &
          &                     zrg_trsol_up_toa          (jcs:jce,jb), zrg_trsol_up_sfc     (jcs:jce,jb),    &
          &                     zrg_trsol_nir_sfc         (jcs:jce,jb), zrg_trsol_vis_sfc    (jcs:jce,jb),    &
          &                     zrg_trsol_par_sfc         (jcs:jce,jb), zrg_fr_nir_sfc_diff  (jcs:jce,jb),    &
          &                     zrg_fr_vis_sfc_diff       (jcs:jce,jb), zrg_fr_par_sfc_diff  (jcs:jce,jb),    &
          &                     zrg_trsol_dn_sfc_diff     (jcs:jce,jb),                                       &
          &                     zrg_trsol_clr_sfc         (jcs:jce,jb), zrg_lwflxall       (jcs:jce,:,jb),    &
          &                     zrg_lwflx_up_sfc          (jcs:jce,jb), zrg_lwflx_clr_sfc    (jcs:jce,jb),    &
          &                     zrg_lwflx_up          (jnps:jnpe,:,jb), zrg_lwflx_dn     (jnps:jnpe,:,jb),    &
          &                     zrg_swflx_up          (jnps:jnpe,:,jb), zrg_swflx_dn     (jnps:jnpe,:,jb),    &
          &                     zrg_lwflx_up_clr      (jnps:jnpe,:,jb), zrg_lwflx_dn_clr (jnps:jnpe,:,jb),    &
          &                     zrg_swflx_up_clr      (jnps:jnpe,:,jb), zrg_swflx_dn_clr (jnps:jnpe,:,jb),    &
          &                     cosmu0mask, zsct, i_startidx_rad, i_endidx_rad, nlev_rgp1)

        ! Add 3D contribution to diffuse radiation
        !$ACC WAIT
        CALL add_3D_diffuse_rad(ecrad_flux, zrg_clc(jcs:jce,:,jb), zrg_pres(jcs:jce,:,jb), zrg_temp(jcs:jce,:,jb),       &
          &                     zrg_cosmu0            (jcs:jce,jb), zrg_fr_nir_sfc_diff  (jcs:jce,jb),         &
          &                     zrg_fr_vis_sfc_diff   (jcs:jce,jb), zrg_fr_par_sfc_diff  (jcs:jce,jb),         &
          &                     zrg_trsol_dn_sfc_diff (jcs:jce,jb), i_startidx_rad, i_endidx_rad, nlev_rg)

      ENDDO !jb_rad
    ENDDO !jb
!$OMP END DO

! CLEANUP
    !$ACC WAIT
    !$ACC END DATA
#ifndef __ECRAD_ACC
    ! The current master of libecrad submodule doesnot support OpenACC.
    ! If a libecrad submodule is used that is ported with OpenACC __ECRAD_ACC
    ! has to be defined to disable this functionality.
    ! Once the libecrad submodule master includes the OpenACC port, everything
    ! that is covered in this #ifndef __ECRAD_ACC can be removed.
    ! In the subroutine ecrad_acc_deallocation the derived ecrad types are dealocated on the GPU.
    CALL ecrad_acc_deallocation(ecrad_conf, ecrad_single_level, ecrad_thermodynamics, ecrad_gas, ecrad_cloud, &
      &                         ecrad_aerosol, ecrad_flux)
#endif
    CALL ecrad_single_level%deallocate()
    CALL ecrad_thermodynamics%deallocate()
    CALL ecrad_gas%deallocate()
    CALL ecrad_cloud%deallocate()
    IF ( ecrad_conf%use_aerosols ) CALL ecrad_aerosol%deallocate()
    CALL ecrad_flux%deallocate()
    !$ACC EXIT DATA DELETE(cosmu0mask)
    DEALLOCATE( cosmu0mask )
    DO jw = 1, ecrad_conf%n_bands_lw
      CALL opt_ptrs_lw(jw)%finalize()
    ENDDO
    DO jw = 1, ecrad_conf%n_bands_sw
      CALL opt_ptrs_sw(jw)%finalize()
    ENDDO
    DEALLOCATE( opt_ptrs_lw, opt_ptrs_sw )
!$OMP END PARALLEL

! Downscale radiative fluxes from reduced radiation grid to full grid
    !$ACC WAIT
    CALL downscale_rad_output(pt_patch%id, pt_par_patch%id,                                         &
      &  nlev_rg, zrg_aclcov, zrg_lwflxall, zrg_trsolall, zrg_trsol_clr_sfc, zrg_lwflx_clr_sfc,     &
      &  zrg_lwflx_up_sfc, zrg_trsol_up_toa, zrg_trsol_up_sfc, zrg_trsol_nir_sfc, zrg_trsol_vis_sfc,&
      &  zrg_trsol_par_sfc, zrg_fr_nir_sfc_diff, zrg_fr_vis_sfc_diff, zrg_fr_par_sfc_diff,          &
      &  zrg_trsol_dn_sfc_diff, zrg_tsfc, zrg_albdif, zrg_emis_rad, zrg_cosmu0, zrg_tot_cld,        &
      &  zlp_tot_cld, zrg_pres_ifc, zlp_pres_ifc, prm_diag%tsfctrad, prm_diag%albdif, aclcov,       &
      &  prm_diag%lwflxall, prm_diag%trsolall, prm_diag%lwflx_up_sfc_rs, prm_diag%trsol_up_toa,     &
      &  prm_diag%trsol_up_sfc, prm_diag%trsol_nir_sfc, prm_diag%trsol_vis_sfc,                     &
      &  prm_diag%trsol_par_sfc, prm_diag%fr_nir_sfc_diff, prm_diag%fr_vis_sfc_diff,                &
      &  prm_diag%fr_par_sfc_diff, prm_diag%trsol_dn_sfc_diff, prm_diag%trsolclr_sfc,               &
      &  prm_diag%lwflxclr_sfc,                                                                     &
      &  zrg_lwflx_up         , zrg_lwflx_dn         , zrg_swflx_up         , zrg_swflx_dn,         &
      &  zrg_lwflx_up_clr     , zrg_lwflx_dn_clr     , zrg_swflx_up_clr     , zrg_swflx_dn_clr,     &
      &  prm_diag%lwflx_up    , prm_diag%lwflx_dn    , prm_diag%swflx_up    , prm_diag%swflx_dn,    &
      &  prm_diag%lwflx_up_clr, prm_diag%lwflx_dn_clr, prm_diag%swflx_up_clr, prm_diag%swflx_dn_clr,&
      &  lacc=.TRUE. )

    !$ACC WAIT
    !$ACC END DATA
    !$ACC EXIT DATA DELETE(zrg_cosmu0, zrg_tsfc, zrg_emis_rad, zrg_albvisdir) &
    !$ACC   DELETE(zrg_albnirdir, zrg_albvisdif, zrg_albnirdif, zrg_pres_ifc, zrg_o3) &
    !$ACC   DELETE(zrg_aeq1, zrg_aeq2, zrg_aeq3, zrg_clc, zrg_aeq4, zrg_aeq5) &
    !$ACC   DELETE(zrg_tot_cld, zrg_pres, zrg_temp, zrg_trsolall, zrg_lwflxall) &
    !$ACC   DELETE(zrg_lwflx_up_sfc, zrg_trsol_up_toa, zrg_trsol_up_sfc) &
    !$ACC   DELETE(zrg_trsol_dn_sfc_diff, zrg_trsol_clr_sfc, zrg_aclcov) &
    !$ACC   DELETE(zrg_trsol_nir_sfc, zrg_trsol_vis_sfc, zrg_trsol_par_sfc) &
    !$ACC   DELETE(zrg_fr_nir_sfc_diff, zrg_fr_vis_sfc_diff, zrg_fr_par_sfc_diff) &
    !$ACC   DELETE(aclcov, zrg_albdif, zrg_rtype, zlp_pres_ifc) &
    !$ACC   DELETE(zlp_tot_cld, zrg_lwflx_clr_sfc, zrg_lwflx_up, zrg_lwflx_dn) &
    !$ACC   DELETE(zrg_swflx_up, zrg_swflx_dn, zrg_lwflx_up_clr, zrg_lwflx_dn_clr) &
    !$ACC   DELETE(zrg_swflx_up_clr, zrg_swflx_dn_clr)
    DEALLOCATE (zrg_cosmu0, zrg_tsfc, zrg_emis_rad, zrg_albvisdir, zrg_albnirdir, zrg_albvisdif,   &
      &         zrg_albnirdif, zrg_pres_ifc, zrg_o3, zrg_aeq1, zrg_aeq2, zrg_aeq3, zrg_clc,        &
      &         zrg_aeq4, zrg_aeq5, zrg_tot_cld, zrg_pres, zrg_temp, zrg_trsolall, zrg_lwflxall,   &
      &         zrg_lwflx_up_sfc, zrg_trsol_up_toa, zrg_trsol_up_sfc, zrg_trsol_dn_sfc_diff,       &
      &         zrg_trsol_clr_sfc, zrg_aclcov, zrg_trsol_nir_sfc, zrg_trsol_vis_sfc,               &
      &         zrg_trsol_par_sfc, zrg_fr_nir_sfc_diff, zrg_fr_vis_sfc_diff,                       &
      &         zrg_fr_par_sfc_diff, aclcov,                                                       &
      &         zrg_albdif, zrg_rtype, zlp_pres_ifc, zlp_tot_cld, zrg_lwflx_clr_sfc,               &
      &         zrg_lwflx_up    , zrg_lwflx_dn    , zrg_swflx_up    , zrg_swflx_dn,                &
      &         zrg_lwflx_up_clr, zrg_lwflx_dn_clr, zrg_swflx_up_clr, zrg_swflx_dn_clr             )

    IF (atm_phy_nwp_config(jg)%icpl_rad_reff > 0) THEN
      !$ACC EXIT DATA DELETE(zrg_reff_liq, zrg_reff_frz)
      DEALLOCATE(zrg_reff_liq, zrg_reff_frz)
    ENDIF
    !$ACC EXIT DATA DELETE(zrg_extra_flds) IF(input_extra_flds%ntot>0)
    IF (input_extra_flds%ntot > 0 ) DEALLOCATE(zrg_extra_flds)
    !$ACC EXIT DATA DELETE(zrg_extra_2D) IF(input_extra_2D%ntot>0)
    IF (input_extra_2D%ntot   > 0 ) DEALLOCATE(zrg_extra_2D)
    !$ACC EXIT DATA DELETE(zrg_extra_reff) IF(input_extra_reff%ntot>0)
    IF (input_extra_reff%ntot > 0 ) DEALLOCATE(zrg_extra_reff)

    CALL input_extra_flds%destruct()
    CALL input_extra_2D%destruct()
    CALL input_extra_reff%destruct()

    !$ACC END DATA

  END SUBROUTINE nwp_ecrad_radiation_reduced
  !---------------------------------------------------------------------------------------

#endif
END MODULE mo_nwp_ecrad_interface
