!
! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_ocean_atmo_coupling

  USE mo_kind,                ONLY: wp
  USE mo_parallel_config,     ONLY: nproma
  USE mo_exception,           ONLY: warning, message
  USE mo_impl_constants,      ONLY: max_char_length
  USE mo_mpi,                 ONLY: p_comm_work, p_sum
  USE mo_physical_constants,  ONLY: tmelt, rhoh2o
  USE mo_run_config,          ONLY: ltimer
  USE mo_dynamics_config,     ONLY: nnew
  USE mo_timer,               ONLY: timer_start, timer_stop, timer_coupling, &
       &                            timer_coupling_put, timer_coupling_get,  &
       &                            timer_coupling_1stget
  USE mo_sync,                ONLY: sync_c, sync_patch_array
  USE mo_util_dbg_prnt,       ONLY: dbg_print
  USE mo_model_domain,        ONLY: t_patch, t_patch_3d

  USE mo_ocean_types
  USE mo_sea_ice_types,       ONLY: t_sea_ice, t_atmos_fluxes
  USE mo_ocean_surface_types, ONLY: t_atmos_for_ocean

  !-------------------------------------------------------------
  ! For the coupling
  !
  USE mo_coupling,            ONLY: lyac_very_1st_get
  USE mo_coupling_utils,      ONLY: def_field
  USE mo_parallel_config,     ONLY: nproma
#ifdef YAC_coupling
  USE mo_yac_finterface,      ONLY: yac_fput, yac_fget, yac_fdef_mask,      &
    &                               YAC_LOCATION_CELL,                      &
    &                               YAC_ACTION_COUPLING,                    &
    &                               YAC_ACTION_OUT_OF_BOUND
#endif
  USE mo_coupling_config,     ONLY: is_coupled_run
  USE mo_hamocc_nml,          ONLY: l_cpl_co2

  !-------------------------------------------------------------

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: str_module = 'mo_ocean_atmo_coupling'  ! Output of module for 1 line debug

  PUBLIC :: construct_ocean_atmo_coupling, couple_ocean_toatmo_fluxes

  INTEGER, TARGET :: field_id_umfl
  INTEGER, TARGET :: field_id_vmfl
  INTEGER, TARGET :: field_id_freshflx
  INTEGER, TARGET :: field_id_heatflx
  INTEGER, TARGET :: field_id_seaice_atm
  INTEGER, TARGET :: field_id_sst
  INTEGER, TARGET :: field_id_oce_u
  INTEGER, TARGET :: field_id_oce_v
  INTEGER, TARGET :: field_id_seaice_oce
  INTEGER, TARGET :: field_id_sp10m
  INTEGER, TARGET :: field_id_co2_vmr
  INTEGER, TARGET :: field_id_co2_flx
  INTEGER, TARGET :: field_id_pres_msl
  INTEGER :: field_id_freshflx_runoff

  INTEGER, SAVE :: nbr_inner_cells

CONTAINS

  !--------------------------------------------------------------------------

  SUBROUTINE construct_ocean_atmo_coupling( &
    patch_3d, comp_id, grid_id, cell_point_id, timestepstring, &
    nbr_inner_cells_)

    TYPE(t_patch_3d ), TARGET, INTENT(in) :: patch_3d
    INTEGER, INTENT(IN) :: comp_id
    INTEGER, INTENT(IN) :: grid_id
    INTEGER, INTENT(IN) :: cell_point_id
    CHARACTER(LEN=*), INTENT(IN) :: timestepstring
    INTEGER, INTENT(IN) :: nbr_inner_cells_

    INTEGER                :: patch_no
    TYPE(t_patch), POINTER :: patch_horz

    INTEGER :: cell_mask_id

    INTEGER :: mask_checksum
    INTEGER :: blockNo, cell_index, i

    LOGICAL, ALLOCATABLE  :: is_valid(:)

    TYPE field_id_ptr
      INTEGER, POINTER :: p
    END TYPE field_id_ptr

    INTEGER, PARAMETER :: no_of_fields = 13
    CHARACTER(LEN=max_char_length) :: field_name(no_of_fields)
    INTEGER                        :: collection_size(no_of_fields)
    TYPE(field_id_ptr)             :: field_ids(no_of_fields)

    CHARACTER(LEN=*), PARAMETER   :: &
      routine = str_module // ':construct_ocean_atmo_coupling'

#ifndef YAC_coupling
    CALL finish(routine, 'built without coupling support.')
#else

    patch_no = 1
    patch_horz => patch_3d%p_patch_2d(patch_no)

    nbr_inner_cells = nbr_inner_cells_

    field_name(1) = "surface_downward_eastward_stress"   ! bundled field containing two components
    collection_size(1) = 2
    field_ids(1)%p => field_id_umfl
    field_name(2) = "surface_downward_northward_stress"  ! bundled field containing two components
    collection_size(2) = 2
    field_ids(2)%p => field_id_vmfl
    field_name(3) = "surface_fresh_water_flux"           ! bundled field containing three components
    collection_size(3) = 3
    field_ids(3)%p => field_id_freshflx
    field_name(4) = "total_heat_flux"                    ! bundled field containing four components
    collection_size(4) = 4
    field_ids(4)%p => field_id_heatflx
    field_name(5) = "atmosphere_sea_ice_bundle"          ! bundled field containing two components
    collection_size(5) = 2
    field_ids(5)%p => field_id_seaice_atm
    field_name(6) = "sea_surface_temperature"
    collection_size(6) = 1
    field_ids(6)%p => field_id_sst
    field_name(7) = "eastward_sea_water_velocity"
    collection_size(7) = 1
    field_ids(7)%p => field_id_oce_u
    field_name(8) = "northward_sea_water_velocity"
    collection_size(8) = 1
    field_ids(8)%p => field_id_oce_v
    field_name(9) = "ocean_sea_ice_bundle"               ! bundled field containing three components
    collection_size(9) = 3
    field_ids(9)%p => field_id_seaice_oce
    field_name(10) = "10m_wind_speed"
    collection_size(10) = 1
    field_ids(10)%p => field_id_sp10m
    field_name(11) = "co2_mixing_ratio"
    collection_size(11) = 1
    field_ids(11)%p => field_id_co2_vmr
    field_name(12) = "co2_flux"
    collection_size(12) = 1
    field_ids(12)%p => field_id_co2_flx
    field_name(13) = "sea_level_pressure"
    collection_size(13) = 1
    field_ids(13)%p => field_id_pres_msl

    !
    ! mask generation : ... not yet defined ...
    !
    ! We could use the patch_horz%cells%decomp_info%owner_local information
    ! e.g. to mask out halo points. We do we get the info about what is local and what
    ! is remote.
    !
    ! The integer land-sea mask:
    !          -2: inner ocean
    !          -1: boundary ocean
    !           1: boundary land
    !           2: inner land
    !
    ! This integer mask for the ocean is available in patch_3D%surface_cell_sea_land_mask(:,:)
    ! The logical mask for the coupler is set to .FALSE. for land points to exclude them from mapping by yac.
    ! These points are not touched by yac.

    mask_checksum = 0
!ICON_OMP_PARALLEL_DO PRIVATE(blockNo,cell_index) REDUCTION(+:mask_checksum) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = 1, patch_horz%nblks_c
      DO cell_index = 1, nproma
        mask_checksum = mask_checksum + ABS(patch_3d%surface_cell_sea_land_mask(cell_index, blockNo))
      ENDDO
    ENDDO
!ICON_OMP_END_PARALLEL_DO
    mask_checksum = p_sum(mask_checksum, comm=p_comm_work)

    IF ( mask_checksum > 0 ) THEN

     ALLOCATE(is_valid(nproma*patch_horz%nblks_c))

!ICON_OMP_PARALLEL_DO PRIVATE(blockNo, cell_index) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = 1, patch_horz%nblks_c
        DO cell_index = 1, nproma
          IF ( patch_3d%surface_cell_sea_land_mask(cell_index, blockNo) < 0 ) THEN
            ! ocean and ocean-coast is valid (-2, -1)
            is_valid((blockNo-1)*nproma+cell_index) = .TRUE.
          ELSE
            ! land is undef (1, 2)
            is_valid((blockNo-1)*nproma+cell_index) = .FALSE.
          ENDIF
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO

      CALL yac_fdef_mask (          &
        & grid_id,                  &
        & patch_horz%n_patch_cells, &
        & YAC_LOCATION_CELL,        &
        & is_valid,                 &
        & cell_mask_id )

      DO i = 1, no_of_fields

        CALL def_field( &
          comp_id, cell_point_id, cell_mask_id, timestepstring, &
          field_name(i), collection_size(i), field_ids(i)%p)

      END DO

      ! Define mask for runoff: all ocean points are valid.
      ! Define field for river runoff using same mask as for other fields
      !   (all ocean points are valid)
      !!slo! Define mask for runoff: ocean coastal points only are valid.
      CALL def_field( &
        comp_id, cell_point_id, cell_mask_id, timestepstring, &
        "river_runoff", 1, field_id_freshflx_runoff)

    ELSE ! no mask

      DO i = 1, no_of_fields

        CALL def_field( &
          comp_id, cell_point_id, timestepstring, &
          TRIM(field_name(i)), collection_size(i), field_ids(i)%p)

      END DO

      CALL def_field( &
        comp_id, cell_point_id, timestepstring, &
        "river_runoff", 1, field_id_freshflx_runoff)

    ENDIF

! YAC_coupling
#endif

  END SUBROUTINE construct_ocean_atmo_coupling

  SUBROUTINE couple_ocean_toatmo_fluxes(patch_3d, ocean_state, ice, atmos_fluxes, atmos_forcing)

    TYPE(t_patch_3d ),TARGET, INTENT(in)        :: patch_3d
    TYPE(t_hydro_ocean_state)                   :: ocean_state
    TYPE(t_sea_ice)                             :: ice
    TYPE(t_atmos_fluxes)                        :: atmos_fluxes
    TYPE(t_atmos_for_ocean)                     :: atmos_forcing

    ! Local declarations for coupling:
    LOGICAL :: write_coupler_restart
    INTEGER :: nbr_hor_cells  ! = inner and halo points
    INTEGER :: cell_index     ! nproma loop count
    INTEGER :: nn             ! block offset
    INTEGER :: blockNo        ! block loop count
    INTEGER :: nlen           ! nproma/npromz
    INTEGER :: no_arr         ! no of arrays in bundle for put/get calls
    TYPE(t_patch), POINTER:: patch_horz

    INTEGER                             :: info, ierror   !< return values from cpl_put/get calls
    REAL(wp), PARAMETER                 :: dummy = 0.0_wp
    REAL(wp)                            :: total_rain

    REAL(wp), ALLOCATABLE :: buffer(:,:)

    CHARACTER(LEN=*), PARAMETER   :: routine = str_module // ':couple_ocean_toatmo_fluxes'

#ifndef YAC_coupling
    CALL finish(routine, 'built without coupling support.')
#else

    IF (.NOT. is_coupled_run() ) RETURN

    IF (ltimer) CALL timer_start(timer_coupling)

    patch_horz   => patch_3D%p_patch_2D(1)

    nbr_hor_cells = patch_horz%n_patch_cells

    ! adjust size if larger bundles are used (no_arr > 4 below)

    ALLOCATE(buffer(nproma*patch_horz%nblks_c,4))

    ! As YAC does not touch masked data an explicit initialisation
    ! is required as some compilers are asked to initialise with NaN
    ! iand as we sometimes loop over the full array

    buffer(:,:) = 0.0_wp
    !
    !  Receive fields from atmosphere
    !   "surface_downward_eastward_stress" bundle  - zonal wind stress component over ice and water
    !   "surface_downward_northward_stress" bundle - meridional wind stress component over ice and water
    !   "surface_fresh_water_flux" bundle          - liquid rain, snowfall, evaporation
    !   "total heat flux" bundle                   - short wave, long wave, sensible, latent heat flux
    !   "atmosphere_sea_ice_bundle"                - sea ice surface and bottom melt potentials
    !   "10m_wind_speed"                           - atmospheric wind speed
    !   "co2_mixing_ratio"                         - co2 mixing ratio
    !
    !  Receive field from HD-model:
    !   "river_runoff"                             - river discharge into the ocean
    !
    !  Send fields to atmosphere:
    !   "sea_surface_temperature"                  - SST
    !   "eastward_sea_water_velocity"              - zonal velocity, u component of ocean surface current
    !   "northward_sea_water_velocity"             - meridional velocity, v component of ocean surface current
    !   "ocean_sea_ice_bundle"                     - ice thickness, snow thickness, ice concentration
    !
    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !  Send fields from ocean to atmosphere
    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !
    write_coupler_restart = .FALSE.
    !
    ! ------------------------------
    !  Send SST
    !   "sea_surface_temperature" - SST
    !
!ICON_OMP_PARALLEL_DO PRIVATE(blockNo, cell_index, nn, nlen) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = 1, patch_horz%nblks_c
      nn = (blockNo-1)*nproma
      IF (blockNo /= patch_horz%nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = patch_horz%npromz_c
      END IF
      DO cell_index = 1, nlen
        ! use SST updated by ocean dynamics - coupling is done at end of timestep
        buffer(nn+cell_index,1) = ocean_state%p_prog(nnew(1))%tracer(cell_index,1,blockNo,1) + tmelt
      ENDDO
    ENDDO
!ICON_OMP_END_PARALLEL_DO
    !    
    IF (ltimer) CALL timer_start(timer_coupling_put)

    CALL yac_fput ( field_id_sst, nbr_hor_cells, 1, buffer(1:nbr_hor_cells,1:1), info, ierror )
    IF ( info > YAC_ACTION_COUPLING .AND. info < YAC_ACTION_OUT_OF_BOUND ) write_coupler_restart = .TRUE.
    IF ( info == YAC_ACTION_OUT_OF_BOUND ) &
         CALL warning('couple_ocean_toatmo_fluxes', &
                      'YAC says fput called after end of run - id=6, SST')

    IF (ltimer) CALL timer_stop(timer_coupling_put)

    !
    ! ------------------------------
    !  Send zonal velocity
    !   "eastward_sea_water_velocity" - zonal velocity, u component of ocean surface current
    !
!ICON_OMP_PARALLEL_DO PRIVATE(blockNo, cell_index, nn, nlen) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = 1, patch_horz%nblks_c
      nn = (blockNo-1)*nproma
      IF (blockNo /= patch_horz%nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = patch_horz%npromz_c
      END IF
      DO cell_index = 1, nlen
        buffer(nn+cell_index,1) = (1.0_wp - ice%conc(cell_index,1,blockNo))*ocean_state%p_diag%u(cell_index,1,blockNo) &
            & + (ice%conc(cell_index,1,blockNo)*ice%u(cell_index,blockNo))
      ENDDO
    ENDDO
!ICON_OMP_END_PARALLEL_DO
    !
    IF (ltimer) CALL timer_start(timer_coupling_put)

    CALL yac_fput ( field_id_oce_u, nbr_hor_cells, 1, buffer(1:nbr_hor_cells,1:1), info, ierror )
    IF ( info > YAC_ACTION_COUPLING .AND. info < YAC_ACTION_OUT_OF_BOUND ) write_coupler_restart = .TRUE.
    IF ( info == YAC_ACTION_OUT_OF_BOUND ) &
         CALL warning('couple_ocean_toatmo_fluxes', &
                      'YAC says fput called after end of run - id=7, u velocity')

    IF (ltimer) CALL timer_stop(timer_coupling_put)

    !
    ! ------------------------------
    !  Send meridional velocity
    !   "northward_sea_water_velocity" - meridional velocity, v component of ocean surface current
    !
!ICON_OMP_PARALLEL_DO PRIVATE(blockNo, cell_index, nn, nlen) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = 1, patch_horz%nblks_c
      nn = (blockNo-1)*nproma
      IF (blockNo /= patch_horz%nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = patch_horz%npromz_c
      END IF
      DO cell_index = 1, nlen
        buffer(nn+cell_index,1) = (1.0_wp - ice%conc(cell_index,1,blockNo))*ocean_state%p_diag%v(cell_index,1,blockNo) &
             & + (ice%conc(cell_index,1,blockNo)*ice%v(cell_index,blockNo))
      ENDDO
    ENDDO
!ICON_OMP_END_PARALLEL_DO
    !
    IF (ltimer) CALL timer_start(timer_coupling_put)

    CALL yac_fput ( field_id_oce_v, nbr_hor_cells, 1, buffer(1:nbr_hor_cells,1:1), info, ierror )
    IF ( info > YAC_ACTION_COUPLING .AND. info < YAC_ACTION_OUT_OF_BOUND ) write_coupler_restart = .TRUE.
    IF ( info == YAC_ACTION_OUT_OF_BOUND ) &
         CALL warning('couple_ocean_toatmo_fluxes', &
                      'YAC says fput called after end of run - id=8, v velocity')

    IF (ltimer) CALL timer_stop(timer_coupling_put)

    !
    ! ------------------------------
    !  Send sea ice bundle
    !   "ocean_sea_ice_bundle" - ice thickness, snow thickness, ice concentration
    !
!ICON_OMP_PARALLEL_DO PRIVATE(blockNo, cell_index, nn, nlen) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = 1, patch_horz%nblks_c
      nn = (blockNo-1)*nproma
      IF (blockNo /= patch_horz%nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = patch_horz%npromz_c
      END IF
      DO cell_index = 1, nlen
        buffer(nn+cell_index,1) = ice%hi  (cell_index,1,blockNo)
        buffer(nn+cell_index,2) = ice%hs  (cell_index,1,blockNo)
        buffer(nn+cell_index,3) = ice%conc(cell_index,1,blockNo)
      ENDDO
    ENDDO
!ICON_OMP_END_PARALLEL_DO
    !
    IF (ltimer) CALL timer_start(timer_coupling_put)

    no_arr = 3
    CALL yac_fput ( field_id_seaice_oce, nbr_hor_cells, no_arr, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    IF ( info > YAC_ACTION_COUPLING .AND. info < YAC_ACTION_OUT_OF_BOUND ) write_coupler_restart = .TRUE.
    IF ( info == YAC_ACTION_OUT_OF_BOUND ) &
         CALL warning('couple_ocean_toatmo_fluxes', &
                      'YAC says fput called after end of run - id=8, sea ice bundle')

    IF (ltimer) CALL timer_stop(timer_coupling_put)

    IF ( write_coupler_restart ) THEN
       CALL message('couple_ocean_toatmo_fluxes', 'YAC says it is put for restart - ids = 6 to 9, ocean fields')
    ENDIF


    IF (l_cpl_co2) THEN
      !
      ! ------------------------------
      !  Send CO2 flux
      !   "co2_flux" - co2flux
      !
      !ICON_OMP_PARALLEL_DO PRIVATE(blockNo, cell_index, nn, nlen) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = 1, patch_horz%nblks_c
        nn = (blockNo-1)*nproma
        IF (blockNo /= patch_horz%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = patch_horz%npromz_c
        END IF
        DO cell_index = 1, nlen
          buffer(nn+cell_index,1) = atmos_forcing%co2flx(cell_index,blockNo) 
        ENDDO
      ENDDO
      !ICON_OMP_END_PARALLEL_DO
      !    
      IF (ltimer) CALL timer_start(timer_coupling_put)

      CALL yac_fput ( field_id_co2_flx, nbr_hor_cells, 1, buffer(1:nbr_hor_cells,1:1), info, ierror )
      IF ( info > YAC_ACTION_COUPLING .AND. info < YAC_ACTION_OUT_OF_BOUND ) THEN
        CALL message('couple_ocean_toatmo_fluxes', &
                     'YAC says it is put for restart - id=13, CO2 flux')
      ENDIF
      IF ( info == YAC_ACTION_OUT_OF_BOUND ) &
         CALL warning('couple_ocean_toatmo_fluxes', &
                      'YAC says fput called after end of run - id=13, CO2 flux')

      IF (ltimer) CALL timer_stop(timer_coupling_put)

    ENDIF

    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !  Receive fields from atmosphere
    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !
    !

    !
    ! ------------------------------
    !  Receive zonal wind stress bundle
    !   "surface_downward_eastward_stress" bundle - zonal wind stress component over ice and water
    !
    IF ( .NOT. lyac_very_1st_get ) THEN
      IF (ltimer) CALL timer_start(timer_coupling_1stget)
    ENDIF

    no_arr = 2
    CALL yac_fget ( field_id_umfl, nbr_hor_cells, no_arr, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    IF ( info > YAC_ACTION_COUPLING .AND. info < YAC_ACTION_OUT_OF_BOUND ) &
         CALL message('couple_ocean_toatmo_fluxes', &
                      'YAC says it is get for restart - id=1, u-stress')
    IF ( info == YAC_ACTION_OUT_OF_BOUND ) &
         CALL warning('couple_ocean_toatmo_fluxes', &
                     'YAC says fget called after end of run - id=1, u-stress')

    IF ( .NOT. lyac_very_1st_get ) THEN
       IF (ltimer) CALL timer_stop(timer_coupling_1stget)
    ENDIF

    lyac_very_1st_get = .FALSE.

    !
    IF (info > 0 .AND. info < 7 ) THEN
      !
!ICON_OMP_PARALLEL_DO PRIVATE(blockNo, cell_index, nn, nlen) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = 1, patch_horz%nblks_c
        nn = (blockNo-1)*nproma
        IF (blockNo /= patch_horz%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = patch_horz%npromz_c
        END IF
        DO cell_index = 1, nlen
          IF ( nn+cell_index > nbr_inner_cells ) THEN
            atmos_fluxes%stress_xw(cell_index,blockNo) = dummy
            atmos_fluxes%stress_x (cell_index,blockNo) = dummy
          ELSE
            atmos_fluxes%stress_xw(cell_index,blockNo) = buffer(nn+cell_index,1)
            atmos_fluxes%stress_x (cell_index,blockNo) = buffer(nn+cell_index,2)
          ENDIF
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO
      !
      CALL sync_patch_array(sync_c, patch_horz, atmos_fluxes%stress_xw(:,:))
      CALL sync_patch_array(sync_c, patch_horz, atmos_fluxes%stress_x (:,:))
    ENDIF

    !
    ! ------------------------------
    !  Receive meridional wind stress bundle
    !   "surface_downward_northward_stress" bundle - meridional wind stress component over ice and water
    !
    IF (ltimer) CALL timer_start(timer_coupling_get)

    no_arr = 2
    CALL yac_fget ( field_id_vmfl, nbr_hor_cells, no_arr, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    IF ( info > YAC_ACTION_COUPLING .AND. info < YAC_ACTION_OUT_OF_BOUND ) &
         CALL message('couple_ocean_toatmo_fluxes', &
                      'YAC says it is get for restart - id=2, v-stress')
    IF ( info == YAC_ACTION_OUT_OF_BOUND ) &
         CALL warning('couple_ocean_toatmo_fluxes', &
                      'YAC says fget called after end of run - id=2, v-stress')

    IF (ltimer) CALL timer_stop(timer_coupling_get)
    !
    IF (info > 0 .AND. info < 7 ) THEN
      !
!ICON_OMP_PARALLEL_DO PRIVATE(blockNo, cell_index, nn, nlen) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = 1, patch_horz%nblks_c
        nn = (blockNo-1)*nproma
        IF (blockNo /= patch_horz%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = patch_horz%npromz_c
        END IF
        DO cell_index = 1, nlen
          IF ( nn+cell_index > nbr_inner_cells ) THEN
            atmos_fluxes%stress_yw(cell_index,blockNo) = dummy
            atmos_fluxes%stress_y (cell_index,blockNo) = dummy
          ELSE
            atmos_fluxes%stress_yw(cell_index,blockNo) = buffer(nn+cell_index,1)
            atmos_fluxes%stress_y (cell_index,blockNo) = buffer(nn+cell_index,2)
          ENDIF
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO
      !
      CALL sync_patch_array(sync_c, patch_horz, atmos_fluxes%stress_yw(:,:))
      CALL sync_patch_array(sync_c, patch_horz, atmos_fluxes%stress_y (:,:))
    ENDIF

    !
    ! ------------------------------
    !  Receive surface fresh water flux bundle
    !   "surface_fresh_water_flux" bundle - liquid rain, snowfall, evaporation

    ! Note: freshwater fluxes are received in kg/m^2/s and are converted to m/s by division by rhoh2o below.
    ! Note: precipitation is the sum of rain and snowfall
    !
    IF (ltimer) CALL timer_start(timer_coupling_get)

    no_arr = 3
    CALL yac_fget ( field_id_freshflx, nbr_hor_cells, no_arr, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    IF ( info > YAC_ACTION_COUPLING .AND. info < YAC_ACTION_OUT_OF_BOUND ) &
         CALL message('couple_ocean_toatmo_fluxes', &
                      'YAC says it is get for restart - id=3, surface fresh water flux')
    IF ( info == YAC_ACTION_OUT_OF_BOUND ) &
         CALL warning('couple_ocean_toatmo_fluxes', &
                      'YAC says fget called after end of run - id=3, surface fresh water flux')

    IF (ltimer) CALL timer_stop(timer_coupling_get)
    !
    IF (info > 0 .AND. info < 7 ) THEN
      !
!ICON_OMP_PARALLEL_DO PRIVATE(blockNo, cell_index, nn, nlen, total_rain) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = 1, patch_horz%nblks_c
        nn = (blockNo-1)*nproma
        IF (blockNo /= patch_horz%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = patch_horz%npromz_c
        END IF
        DO cell_index = 1, nlen
          IF ( nn+cell_index > nbr_inner_cells ) THEN
            total_rain                                              = dummy
            atmos_fluxes%FrshFlux_SnowFall     (cell_index,blockNo) = dummy
            atmos_fluxes%FrshFlux_Evaporation  (cell_index,blockNo) = dummy
            atmos_fluxes%FrshFlux_Precipitation(cell_index,blockNo) = dummy
          ELSE
            total_rain                                              = buffer(nn+cell_index,1) / rhoh2o
            atmos_fluxes%FrshFlux_SnowFall     (cell_index,blockNo) = buffer(nn+cell_index,2) / rhoh2o
            atmos_fluxes%FrshFlux_Evaporation  (cell_index,blockNo) = buffer(nn+cell_index,3) / rhoh2o
            atmos_fluxes%FrshFlux_Precipitation(cell_index,blockNo) = total_rain + &
              &  atmos_fluxes%FrshFlux_SnowFall(cell_index,blockNo)
          ENDIF
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO
      !
      CALL sync_patch_array(sync_c, patch_horz, atmos_fluxes%FrshFlux_Precipitation(:,:))
      CALL sync_patch_array(sync_c, patch_horz, atmos_fluxes%FrshFlux_SnowFall     (:,:))
      CALL sync_patch_array(sync_c, patch_horz, atmos_fluxes%FrshFlux_Evaporation  (:,:))
    END IF

    !
    ! ------------------------------
    !  Receive total heat flux bundle
    !   "total heat flux" bundle - short wave, long wave, sensible, latent heat flux
    !
    ! atmos_fluxes%swflx(:,:)  ocean short wave heat flux                              [W/m2]
    ! atmos_fluxes%lwflx(:,:)  ocean long  wave heat fluxe                             [W/m2]
    ! atmos_fluxes%ssflx(:,:)  ocean sensible heat fluxes                              [W/m2]
    ! atmos_fluxes%slflx(:,:)  ocean latent heat fluxes                                [W/m2]
    !
    IF (ltimer) CALL timer_start(timer_coupling_get)

    no_arr = 4
    CALL yac_fget ( field_id_heatflx, nbr_hor_cells, no_arr, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    IF ( info > YAC_ACTION_COUPLING .AND. info < YAC_ACTION_OUT_OF_BOUND ) &
         CALL message('couple_ocean_toatmo_fluxes', &
                      'YAC says it is get for restart - id=4, heat flux')
    IF ( info == YAC_ACTION_OUT_OF_BOUND ) &
         CALL warning('couple_ocean_toatmo_fluxes', &
                      'YAC says fget called after end of run - id=4, heat flux')

    IF (ltimer) CALL timer_stop(timer_coupling_get)
    !
    IF (info > 0 .AND. info < 7 ) THEN
      !
!ICON_OMP_PARALLEL_DO PRIVATE(blockNo, cell_index, nn, nlen) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = 1, patch_horz%nblks_c
        nn = (blockNo-1)*nproma
        IF (blockNo /= patch_horz%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = patch_horz%npromz_c
        END IF
        DO cell_index = 1, nlen
          IF ( nn+cell_index > nbr_inner_cells ) THEN
            atmos_fluxes%HeatFlux_ShortWave(cell_index,blockNo) = dummy
            atmos_fluxes%HeatFlux_LongWave (cell_index,blockNo) = dummy
            atmos_fluxes%HeatFlux_Sensible (cell_index,blockNo) = dummy
            atmos_fluxes%HeatFlux_Latent   (cell_index,blockNo) = dummy
          ELSE
            atmos_fluxes%HeatFlux_ShortWave(cell_index,blockNo) = buffer(nn+cell_index,1)
            atmos_fluxes%HeatFlux_LongWave (cell_index,blockNo) = buffer(nn+cell_index,2)
            atmos_fluxes%HeatFlux_Sensible (cell_index,blockNo) = buffer(nn+cell_index,3)
            atmos_fluxes%HeatFlux_Latent   (cell_index,blockNo) = buffer(nn+cell_index,4)
          ENDIF
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO
      !
      CALL sync_patch_array(sync_c, patch_horz, atmos_fluxes%HeatFlux_ShortWave(:,:))
      CALL sync_patch_array(sync_c, patch_horz, atmos_fluxes%HeatFlux_LongWave (:,:))
      CALL sync_patch_array(sync_c, patch_horz, atmos_fluxes%HeatFlux_Sensible (:,:))
      CALL sync_patch_array(sync_c, patch_horz, atmos_fluxes%HeatFlux_Latent   (:,:))

      ! sum of fluxes for ocean boundary condition
!ICON_OMP_PARALLEL_DO PRIVATE(blockNo, cell_index, nlen) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = 1, patch_horz%nblks_c
        IF (blockNo /= patch_horz%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = patch_horz%npromz_c
        END IF
        DO cell_index = 1, nlen
          atmos_fluxes%HeatFlux_Total(cell_index,blockNo) = atmos_fluxes%HeatFlux_ShortWave(cell_index,blockNo) &
        &                                                 + atmos_fluxes%HeatFlux_LongWave (cell_index,blockNo) &
        &                                                 + atmos_fluxes%HeatFlux_Sensible (cell_index,blockNo) &
        &                                                 + atmos_fluxes%HeatFlux_Latent   (cell_index,blockNo)
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO
    ENDIF

    !
    ! ------------------------------
    !  Receive sea ice flux bundle
    !   "atmosphere_sea_ice_bundle" - sea ice surface and bottom melt potentials
    !
    ! ice%Qtop(:,:)         Surface melt potential of ice                           [W/m2]
    ! ice%Qbot(:,:)         Bottom melt potential of ice                            [W/m2]
    !
    IF (ltimer) CALL timer_start(timer_coupling_get)

    no_arr = 2
    CALL yac_fget ( field_id_seaice_atm, nbr_hor_cells, no_arr, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    IF ( info > YAC_ACTION_COUPLING .AND. info < YAC_ACTION_OUT_OF_BOUND ) &
         CALL message('couple_ocean_toatmo_fluxes', &
                      'YAC says it is get for restart - id=5, sea ice')
    IF ( info == YAC_ACTION_OUT_OF_BOUND ) &
         CALL warning('couple_ocean_toatmo_fluxes', &
                      'YAC says fget called after end of run - id=5, sea ice')

    IF (ltimer) CALL timer_stop(timer_coupling_get)
    !
    IF (info > 0 .AND. info < 7 ) THEN
      !
!ICON_OMP_PARALLEL_DO PRIVATE(blockNo, cell_index, nn, nlen) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = 1, patch_horz%nblks_c
        nn = (blockNo-1)*nproma
        IF (blockNo /= patch_horz%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = patch_horz%npromz_c
        END IF
        DO cell_index = 1, nlen
          IF ( nn+cell_index > nbr_inner_cells ) THEN
            ice%qtop(cell_index,1,blockNo) = dummy
            ice%qbot(cell_index,1,blockNo) = dummy
          ELSE
            ice%qtop(cell_index,1,blockNo) = buffer(nn+cell_index,1)
            ice%qbot(cell_index,1,blockNo) = buffer(nn+cell_index,2)
          ENDIF
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO
      !
      CALL sync_patch_array(sync_c, patch_horz, ice%qtop(:,1,:))
      CALL sync_patch_array(sync_c, patch_horz, ice%qbot(:,1,:))
    END IF

    !
    ! ------------------------------
    !  Receive 10m wind speed
    !   "10m_wind_speed" - atmospheric wind speed
    !
    IF (ltimer) CALL timer_start(timer_coupling_get)

    no_arr = 1
    CALL yac_fget ( field_id_sp10m, nbr_hor_cells, no_arr, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    IF ( info > YAC_ACTION_COUPLING .AND. info < YAC_ACTION_OUT_OF_BOUND ) &
         CALL message('couple_ocean_toatmo_fluxes', &
                      'YAC says it is get for restart - id=10, wind speed')
    IF ( info == YAC_ACTION_OUT_OF_BOUND ) &
         CALL warning('couple_ocean_toatmo_fluxes', &
                      'YAC says fget called after end of run - id=10, wind speed')

    IF (ltimer) CALL timer_stop(timer_coupling_get)
    !
    IF (info > 0 .AND. info < 7 ) THEN
      !
!!ICON_OMP_PARALLEL_DO PRIVATE(blockNo, cell_index, nn, nlen) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = 1, patch_horz%nblks_c
        nn = (blockNo-1)*nproma
        IF (blockNo /= patch_horz%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = patch_horz%npromz_c
        END IF
        DO cell_index = 1, nlen
          IF ( nn+cell_index > nbr_inner_cells ) THEN
            atmos_forcing%fu10(cell_index,blockNo) = dummy
          ELSE
            atmos_forcing%fu10(cell_index,blockNo) = buffer(nn+cell_index,1)
            IF ( atmos_forcing%fu10(cell_index,blockNo) < 0.0_wp ) atmos_forcing%fu10(cell_index,blockNo) = 0.0_wp
          ENDIF
        ENDDO
      ENDDO
!!ICON_OMP_END_PARALLEL_DO
      !
      CALL sync_patch_array(sync_c, patch_horz, atmos_forcing%fu10(:,:))
    END IF

    !  Receive slp
    !   atmospheric sea level pressure
    !
    IF (ltimer) CALL timer_start(timer_coupling_get)

    no_arr = 1
    CALL yac_fget ( field_id_pres_msl, nbr_hor_cells, no_arr, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    IF ( info > YAC_ACTION_COUPLING .AND. info < YAC_ACTION_OUT_OF_BOUND ) &
         CALL message('couple_ocean_toatmo_fluxes', &
                      'YAC says it is get for restart - id=14, sea level pressure')
    IF ( info == YAC_ACTION_OUT_OF_BOUND ) &
         CALL warning('couple_ocean_toatmo_fluxes', &
                      'YAC says fget called after end of run - id=14, sea level pressure')

    IF (ltimer) CALL timer_stop(timer_coupling_get)
    !
    IF (info > 0 .AND. info < 7 ) THEN
      !
!!ICON_OMP_PARALLEL_DO PRIVATE(blockNo, cell_index, nn, nlen) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = 1, patch_horz%nblks_c
        nn = (blockNo-1)*nproma
        IF (blockNo /= patch_horz%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = patch_horz%npromz_c
        END IF
        DO cell_index = 1, nlen
          IF ( nn+cell_index > nbr_inner_cells ) THEN
            atmos_forcing%pao(cell_index,blockNo) = dummy
          ELSE
            atmos_forcing%pao(cell_index,blockNo) = buffer(nn+cell_index,1)
          ENDIF
        ENDDO
      ENDDO
!!ICON_OMP_END_PARALLEL_DO
      !
      CALL sync_patch_array(sync_c, patch_horz, atmos_forcing%pao(:,:))
    END IF
    IF(l_cpl_co2)then
    !
    ! ------------------------------
    !  Receive co2 mixing ratio
    !   "co2 mixing ratio" - atmospheric co2 mixing ratio
    !
    IF (ltimer) CALL timer_start(timer_coupling_get)

    no_arr = 1
    CALL yac_fget ( field_id_co2_vmr, nbr_hor_cells, no_arr, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    IF ( info > YAC_ACTION_COUPLING .AND. info < YAC_ACTION_OUT_OF_BOUND ) &
         CALL message('couple_ocean_toatmo_fluxes', &
                      'YAC says it is get for restart - id=12, co2 mr')
    IF ( info == YAC_ACTION_OUT_OF_BOUND ) &
         CALL warning('couple_ocean_toatmo_fluxes', &
                      'YAC says fget called after end of run - id=12, co2 mixing ratio')

    IF (ltimer) CALL timer_stop(timer_coupling_get)
    !
    IF (info > 0 .AND. info < 7 ) THEN
      !
!!ICON_OMP_PARALLEL_DO PRIVATE(blockNo, cell_index, nn, nlen) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = 1, patch_horz%nblks_c
        nn = (blockNo-1)*nproma
        IF (blockNo /= patch_horz%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = patch_horz%npromz_c
        END IF
        DO cell_index = 1, nlen
          IF ( nn+cell_index > nbr_inner_cells ) THEN
            atmos_forcing%co2(cell_index,blockNo) = dummy
          ELSE
            atmos_forcing%co2(cell_index,blockNo) = buffer(nn+cell_index,1)
            IF ( atmos_forcing%co2(cell_index,blockNo) < 0.0_wp ) atmos_forcing%co2(cell_index,blockNo) = 0.0_wp
          ENDIF
        ENDDO
      ENDDO
!!ICON_OMP_END_PARALLEL_DO
      !
      CALL sync_patch_array(sync_c, patch_horz, atmos_forcing%co2(:,:))
    END IF
    END IF !l_cpl_co2

    !
    ! ------------------------------
    !  Receive river runoff
    !   "river_runoff" - river discharge into the ocean
    !
    ! Note: river runoff fluxes are received in m^3/s and are converted to m/s by division by whole grid area
    !
    IF (ltimer) CALL timer_start(timer_coupling_get)

    buffer(:,:) = 0.0_wp  ! mandatory as river ruoff comes with a different mask!

    CALL yac_fget ( field_id_freshflx_runoff, nbr_hor_cells, 1, buffer(1:nbr_hor_cells,1:1), info, ierror )
    IF ( info > YAC_ACTION_COUPLING .AND. info < YAC_ACTION_OUT_OF_BOUND )&
         CALL message('couple_ocean_toatmo_fluxes', &
                      'YAC says it is get for restart - id=11, runoff')
    IF ( info == YAC_ACTION_OUT_OF_BOUND ) &
         CALL warning('couple_ocean_toatmo_fluxes', &
                      'YAC says fget called after end of run - id=11, runoff')

    IF (ltimer) CALL timer_stop(timer_coupling_get)
    !
    IF (info > 0 .AND. info < 7 ) THEN
      !
!ICON_OMP_PARALLEL_DO PRIVATE(blockNo, cell_index, nn, nlen) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = 1, patch_horz%nblks_c
        nn = (blockNo-1)*nproma
        IF (blockNo /= patch_horz%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = patch_horz%npromz_c
        END IF
        DO cell_index = 1, nlen
          IF ( nn+cell_index > nbr_inner_cells ) THEN
            atmos_fluxes%FrshFlux_Runoff(cell_index,blockNo) = dummy
          ELSE
    ! !!! Note: freshwater fluxes are received in kg/m^2/s and are converted to m/s by division by rhoh2o below.
    ! !!!   atmos_fluxes%FrshFlux_Runoff(cell_index,blockNo) = buffer(nn+cell_index,1) / rhoh2o
    ! discharge_ocean is in m3/s
            atmos_fluxes%FrshFlux_Runoff(cell_index,blockNo) = buffer(nn+cell_index,1) / patch_horz%cells%area(cell_index,blockNo)
          ENDIF
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO
      !
      CALL sync_patch_array(sync_c, patch_horz, atmos_fluxes%FrshFlux_Runoff(:,:))
    END IF

    !---------DEBUG DIAGNOSTICS-------------------------------------------

    CALL dbg_print('toatmo: AtmFluxStress_x  ',atmos_fluxes%stress_x              ,str_module,3,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: AtmFluxStress_xw ',atmos_fluxes%stress_xw             ,str_module,3,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: AtmFluxStress_y  ',atmos_fluxes%stress_y              ,str_module,4,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: AtmFluxStress_yw ',atmos_fluxes%stress_yw             ,str_module,4,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: FrshFluxPrecip   ',atmos_fluxes%FrshFlux_Precipitation,str_module,3,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: FrshFluxEvapo    ',atmos_fluxes%FrshFlux_Evaporation  ,str_module,3,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: FrshFluxSnowFall ',atmos_fluxes%FrshFlux_SnowFall     ,str_module,3,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: HeatFluxTotal    ',atmos_fluxes%HeatFlux_Total        ,str_module,2,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: HeatFluxShortwave',atmos_fluxes%HeatFlux_ShortWave    ,str_module,3,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: HeatFluxLongwave ',atmos_fluxes%HeatFlux_Longwave     ,str_module,4,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: HeatFluxSensible ',atmos_fluxes%HeatFlux_Sensible     ,str_module,4,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: HeatFluxLatent   ',atmos_fluxes%HeatFlux_Latent       ,str_module,4,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: ice%Qtop         ',ice%qtop                           ,str_module,4,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: ice%Qbot         ',ice%qbot                           ,str_module,3,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: FrshFluxRunoff   ',atmos_fluxes%FrshFlux_Runoff       ,str_module,3,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: 10m_wind_speed   ',atmos_forcing%fu10                 ,str_module,3,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: sea_level_pressure',atmos_forcing%pao                 ,str_module,3,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: tracer(1): SST   ', ocean_state%p_prog(nnew(1))%tracer(:,1,:,1) + tmelt &
      &                                                                           ,str_module,2,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: p_diag%u         ',ocean_state%p_diag%u(:,1,:)        ,str_module,4,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: p_diag%v         ',ocean_state%p_diag%v(:,1,:)        ,str_module,4,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: ice%hi           ',ice%hi(:,1,:)                      ,str_module,3,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: ice%hs           ',ice%hs(:,1,:)                      ,str_module,4,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: ice%conc         ',ice%conc(:,1,:)                    ,str_module,4,in_subset=patch_horz%cells%owned)

    !---------------------------------------------------------------------

    DEALLOCATE(buffer)

    IF (ltimer) CALL timer_stop(timer_coupling)

! YAC_coupling
#endif

  END SUBROUTINE couple_ocean_toatmo_fluxes
  !--------------------------------------------------------------------------

END MODULE mo_ocean_atmo_coupling

