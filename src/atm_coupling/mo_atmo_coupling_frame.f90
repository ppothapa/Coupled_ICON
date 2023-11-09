! @brief Initialisation of atmosphere coupling
!
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

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_atmo_coupling_frame

  USE mo_kind                ,ONLY: wp
  USE mo_model_domain        ,ONLY: t_patch
  USE mo_ext_data_state      ,ONLY: ext_data
#ifndef __NO_AES__
  USE mo_aes_phy_memory      ,ONLY: prm_field
#endif
  USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config

  USE mo_parallel_config     ,ONLY: nproma

  USE mo_run_config          ,ONLY: iforcing, ltimer, modelTimeStep
  USE mo_timer,               ONLY: timer_start, timer_stop, timer_coupling_init

  USE mo_impl_constants      ,ONLY: MAX_CHAR_LENGTH, inwp, iaes, LSS_JSBACH

#if !defined(__NO_JSBACH__) && !defined(__NO_JSBACH_HD__)
  USE mo_interface_hd_ocean  ,ONLY: jsb_fdef_hd_fields
#endif

  USE mo_master_control      ,ONLY: get_my_process_name

  USE mo_mpi                 ,ONLY: p_pe_work, p_comm_work, p_sum
  USE mo_parallel_config     ,ONLY: nproma

  USE mo_coupling_config     ,ONLY: is_coupled_run
  USE mo_aes_rad_config      ,ONLY: aes_rad_config
  USE mo_aes_phy_config      ,ONLY: aes_phy_config
  USE mo_time_config         ,ONLY: time_config
  USE mo_util_dbg_prnt       ,ONLY: dbg_print

  USE mo_exception           ,ONLY: finish, message

  USE mo_yac_finterface      ,ONLY: yac_fget_version, yac_fdef_comp,        &
    &                               yac_fdef_datetime, yac_fdef_grid,       &
    &                               yac_fdef_points, yac_fset_global_index, &
    &                               yac_fset_core_mask, yac_fdef_mask,      &
    &                               yac_fdef_field_mask, yac_fenddef,       &
    &                               YAC_LOCATION_CELL, YAC_TIME_UNIT_ISO_FORMAT, &
    &                               yac_fget_field_collection_size

  USE mtime                  ,ONLY: datetimeToString, MAX_DATETIME_STR_LEN

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: str_module = 'mo_atmo_coupling_frame' ! Output of module for debug

  INTEGER :: min_no_of_fields = 13
  INTEGER :: max_no_of_fields = 24

  ! These constants are only valid AFTER calling construct_atmo_coupling !
  INTEGER, PARAMETER :: CPF_UMFL = 1 !< Surface zonal stress for water and ice.
  INTEGER, PARAMETER :: CPF_VMFL = 2 !< Surface meridional stress for water and ice.
  INTEGER, PARAMETER :: CPF_FRESHFLX = 3 !< Fresh water flux (rain, snow, evap).
  INTEGER, PARAMETER :: CPF_HEATFLX = 4 !< Heat flux (SW net, LW net, sensible, latent).
  INTEGER, PARAMETER :: CPF_SEAICE_ATM = 5 !< Sea-ice fluxes (SW+LW+H+L, conductive at bottom).
  INTEGER, PARAMETER :: CPF_SST = 6 !< Sea-surface temperature.
  INTEGER, PARAMETER :: CPF_OCE_U = 7 !< Zonal ocean surface velocity.
  INTEGER, PARAMETER :: CPF_OCE_V = 8 !< Meridional ocean surface velocity.
  INTEGER, PARAMETER :: CPF_SEAICE_OCE = 9 !< Sea-ice state (h_ice, h_snow, fr_seaice).
  INTEGER, PARAMETER :: CPF_SP10M = 10 !< 10m wind speed.
  INTEGER, PARAMETER :: CPF_CO2_VMR = 11 !< CO2 volume mixing ratio at surface in PPM.
  INTEGER, PARAMETER :: CPF_CO2_FLX = 12 !< CO2 flux.
  INTEGER, PARAMETER :: CPF_PRES_MSL = 13 !< Sea-level pressure.

  PUBLIC :: construct_atmo_coupling
  PUBLIC :: lyac_very_1st_get, nbr_inner_cells, mask_checksum, field_id

  PUBLIC :: CPF_UMFL, CPF_VMFL, CPF_FRESHFLX, CPF_HEATFLX, CPF_SEAICE_ATM, CPF_SST, &
      & CPF_OCE_U, CPF_OCE_V, CPF_SEAICE_OCE, CPF_SP10M, CPF_CO2_VMR, CPF_CO2_FLX, &
      & CPF_PRES_MSL

  INTEGER, ALLOCATABLE  :: field_id(:), tmp_field_id(:)

  INTEGER, SAVE         :: nbr_inner_cells
  INTEGER, SAVE         :: mask_checksum = -1
  LOGICAL, SAVE         :: lyac_very_1st_get

CONTAINS

  !>
  !! SUBROUTINE construct_atmo_coupling -- the initialisation for the coupling
  !! of atmosphere and the ocean, through a coupler

  SUBROUTINE construct_atmo_coupling (p_patch)

    TYPE(t_patch), TARGET, INTENT(IN) :: p_patch(:)

    INTEGER :: error_status

    TYPE(t_patch), POINTER :: patch_horz

    !---------------------------------------------------------------------
    ! 11. Do the setup for the coupled run
    !
    ! For the time being this could all go into a subroutine which is
    ! common to atmo and ocean. Does this make sense if the setup deviates
    ! too much in future.
    !---------------------------------------------------------------------

    CHARACTER(LEN=max_char_length) :: grid_name
    CHARACTER(LEN=max_char_length) :: comp_name

    INTEGER :: comp_id
    INTEGER :: comp_ids(1)
    INTEGER :: cell_point_ids(1)
    INTEGER :: cell_mask_ids(1)
    INTEGER :: grid_id

    INTEGER :: jg
    INTEGER :: nblks
    INTEGER :: jb, jc, nn
    INTEGER :: nbr_vertices_per_cell

    REAL(wp), ALLOCATABLE :: buffer_lon(:)
    REAL(wp), ALLOCATABLE :: buffer_lat(:)
    INTEGER,  ALLOCATABLE :: buffer_c(:,:)

    LOGICAL,  ALLOCATABLE :: is_valid(:)

    REAL(wp), ALLOCATABLE :: lsmnolake(:,:)

    REAL(wp), PARAMETER :: eps = 1.E-10_wp

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: startdatestring
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: stopdatestring

    CHARACTER(len=40), ALLOCATABLE :: field_name(:)
    INTEGER :: collection_size(max_no_of_fields)
    INTEGER :: intermediate_no_of_fields

    LOGICAL :: needs_o3_coupling = .FALSE.
    LOGICAL :: needs_areo_coupling = .FALSE.

    ! Skip time measurement of the very first yac_fget
    ! as this will measure mainly the wait time caused
    ! by the initialisation of the model components
    ! and does not tell us much about the load balancing
    ! in subsequent calls.

    lyac_very_1st_get = .TRUE.

    IF ( .NOT. is_coupled_run() ) RETURN

    IF (ltimer) CALL timer_start (timer_coupling_init)

    ! setup default coupling fields {{{
    ALLOCATE(field_name(max_no_of_fields))
    ALLOCATE(field_id(min_no_of_fields))
    field_name(CPF_UMFL)             = "surface_downward_eastward_stress"
    collection_size(CPF_UMFL)        = 2

    field_name(CPF_VMFL)             = "surface_downward_northward_stress"
    collection_size(CPF_VMFL)        = 2

    field_name(CPF_FRESHFLX)         = "surface_fresh_water_flux"
    collection_size(CPF_FRESHFLX)    = 3

    field_name(CPF_HEATFLX)          = "total_heat_flux"
    collection_size(CPF_HEATFLX)     = 4

    field_name(CPF_SEAICE_ATM)       = "atmosphere_sea_ice_bundle"
    collection_size(CPF_SEAICE_ATM)  = 2

    field_name(CPF_SST)              = "sea_surface_temperature"
    collection_size(CPF_SST)         = 1

    field_name(CPF_OCE_U)            = "eastward_sea_water_velocity"
    collection_size(CPF_OCE_U)       = 1

    field_name(CPF_OCE_V)            = "northward_sea_water_velocity"
    collection_size(CPF_OCE_V)       = 1

    field_name(CPF_SEAICE_OCE)       = "ocean_sea_ice_bundle"
    collection_size(CPF_SEAICE_OCE)  = 3

    field_name(CPF_SP10M)            = "10m_wind_speed"
    collection_size(CPF_SP10M)       = 1

    field_name(CPF_CO2_VMR)          = "co2_mixing_ratio"
    collection_size(CPF_CO2_VMR)     = 1

    field_name(CPF_CO2_FLX)          = "co2_flux"
    collection_size(CPF_CO2_FLX)     = 1

    field_name(CPF_PRES_MSL)         = "sea_level_pressure"
    collection_size(CPF_PRES_MSL)    = 1
    !}}}

    jg = 1
    patch_horz => p_patch(jg)

    ! Print the YAC version
    CALL message('Running ICON atmosphere in coupled mode with YAC version ', TRIM(yac_fget_version()) )

    ! Overwrite job start and end date with component data
    CALL datetimeToString(time_config%tc_startdate, startdatestring)
    CALL datetimeToString(time_config%tc_stopdate, stopdatestring)

    CALL yac_fdef_datetime ( start_datetime = TRIM(startdatestring), &
         &                   end_datetime   = TRIM(stopdatestring)   )

    ! Inform the coupler about what we are
    comp_name = TRIM(get_my_process_name())
    CALL yac_fdef_comp ( TRIM(comp_name), comp_id )
    comp_ids(1) = comp_id

    ! Announce one grid (patch) to the coupler
    grid_name = "icon_atmos_grid"

    ! Extract cell information {{{
    !
    ! cartesian coordinates of cell vertices are stored in
    ! patch_horz%verts%cartesian(:,:)%x(1:3)
    ! Here we use the longitudes and latitudes in rad.

    nblks = max(patch_horz%nblks_c,patch_horz%nblks_v)

    ALLOCATE(buffer_lon(nproma*nblks))
    ALLOCATE(buffer_lat(nproma*nblks))
    ALLOCATE(buffer_c(3,nproma*nblks))

    ALLOCATE(lsmnolake(nproma,nblks))

    nbr_vertices_per_cell = 3

    !ICON_OMP_PARALLEL
    !ICON_OMP_DO PRIVATE(jb, jc, nn) ICON_OMP_RUNTIME_SCHEDULE
    DO jb = 1, patch_horz%nblks_v
      DO jc = 1, nproma
        nn = (jb-1)*nproma+jc
        buffer_lon(nn) = patch_horz%verts%vertex(jc,jb)%lon
        buffer_lat(nn) = patch_horz%verts%vertex(jc,jb)%lat
      ENDDO
    ENDDO
    !ICON_OMP_END_DO NOWAIT

    !ICON_OMP_DO PRIVATE(jb, jc, nn) ICON_OMP_RUNTIME_SCHEDULE
    DO jb = 1, patch_horz%nblks_c
      DO jc = 1, nproma
        nn = (jb-1)*nproma+jc
        buffer_c(1,nn) = (patch_horz%cells%vertex_blk(jc,jb,1)-1)*nproma + &
          &               patch_horz%cells%vertex_idx(jc,jb,1)
        buffer_c(2,nn) = (patch_horz%cells%vertex_blk(jc,jb,2)-1)*nproma + &
          &               patch_horz%cells%vertex_idx(jc,jb,2)
        buffer_c(3,nn) = (patch_horz%cells%vertex_blk(jc,jb,3)-1)*nproma + &
                          patch_horz%cells%vertex_idx(jc,jb,3)
      ENDDO
    ENDDO
    !ICON_OMP_END_DO
    !ICON_OMP_END_PARALLEL

    ! Description of elements, here as unstructured grid
    CALL yac_fdef_grid(           &
      & TRIM(grid_name),          &
      & patch_horz%n_patch_verts, &
      & patch_horz%n_patch_cells, &
      & nbr_vertices_per_cell,    &
      & buffer_lon,               &
      & buffer_lat,               &
      & buffer_c,                 &
      & grid_id)

    !
    ! Define cell center points (location = 0)
    !
    ! cartesian coordinates of cell centers are stored in
    ! patch_horz%cells%cartesian_center(:,:)%x(1:3)
    ! Here we use the longitudes and latitudes.

    !ICON_OMP_PARALLEL_DO PRIVATE(jb, jc, nn) ICON_OMP_RUNTIME_SCHEDULE
    DO jb = 1, patch_horz%nblks_c
      DO jc = 1, nproma
        nn = (jb-1)*nproma+jc
        buffer_lon(nn) = patch_horz%cells%center(jc,jb)%lon
        buffer_lat(nn) = patch_horz%cells%center(jc,jb)%lat
      ENDDO
    ENDDO
    !ICON_OMP_END_PARALLEL_DO

    ! center points in cells (needed e.g. for patch recovery and nearest neighbour interpolation)
    CALL yac_fdef_points (        &
      & grid_id,                  &
      & patch_horz%n_patch_cells, &
      & YAC_LOCATION_CELL,        &
      & buffer_lon,               &
      & buffer_lat,               &
      & cell_point_ids(1) )

    DEALLOCATE (buffer_lon, buffer_lat, buffer_c)

    CALL yac_fset_global_index (                &
      & patch_horz%cells%decomp_info%glb_index, &
      & YAC_LOCATION_CELL,                      &
      & grid_id )

    ALLOCATE(is_valid(nproma*patch_horz%nblks_c))

    nbr_inner_cells = 0
!ICON_OMP_PARALLEL_DO PRIVATE(jc) REDUCTION(+:nbr_inner_cells) ICON_OMP_RUNTIME_SCHEDULE
    DO jc = 1, patch_horz%n_patch_cells
       IF ( p_pe_work == patch_horz%cells%decomp_info%owner_local(jc) ) THEN
         is_valid(jc) = .TRUE.
         nbr_inner_cells = nbr_inner_cells + 1
       ELSE
         is_valid(jc) = .FALSE.
       ENDIF
    ENDDO
!ICON_OMP_END_PARALLEL_DO

    CALL yac_fset_core_mask ( &
      & is_valid,             &
      & YAC_LOCATION_CELL,    &
      & grid_id )

    ! }}}
    ! The integer land-sea mask:
    !          -2: inner ocean
    !          -1: boundary ocean
    !           1: boundary land
    !           2: inner land
    !
    ! The (fractional) mask which is used in the AES physics is prm_field(1)%lsmask(:,:).
    !
    ! The logical mask for the coupler must be generated from the fractional mask by setting
    !   only those gridpoints to land that have no ocean part at all (lsf<1 is ocean).
    ! The logical mask is then set to .FALSE. for land points to exclude them from mapping by yac.
    ! These points are not touched by yac.
    !

    SELECT CASE( iforcing ) !{{{

      CASE ( inwp )

        !ICON_OMP_PARALLEL PRIVATE(jb,jc)
          !ICON_OMP_WORKSHARE
          is_valid(:) = .FALSE.
          !ICON_OMP_END_WORKSHARE

          !ICON_OMP_DO ICON_OMP_DEFAULT_SCHEDULE
          DO jb = 1, patch_horz%nblks_c
            DO jc = 1, ext_data(jg)%atm%list_sea%ncount(jb)
              is_valid((jb-1)*nproma + ext_data(jg)%atm%list_sea%idx(jc,jb)) = .TRUE.
            END DO
          END DO
          !ICON_OMP_END_DO
        !ICON_OMP_END_PARALLEL

        CALL dbg_print('AtmFrame: fr_land',ext_data(jg)%atm%fr_land,str_module,3,in_subset=patch_horz%cells%owned)
        CALL dbg_print('AtmFrame: fr_lake',ext_data(jg)%atm%fr_lake,str_module,3,in_subset=patch_horz%cells%owned)

      CASE ( iaes )
#ifdef __NO_AES__
        CALL finish ('mo_atmo_coupling_frame:construct_atmo_coupling', &
            & 'coupled model needs aes; remove --disable-aes and reconfigure')
#else
        !ICON_OMP_PARALLEL_DO PRIVATE(jb,jc) ICON_OMP_RUNTIME_SCHEDULE
        DO jb = 1, patch_horz%nblks_c
            DO jc = 1, nproma
              !  slo: caution - lsmask includes alake, must be added to refetch pure lsm:
              lsmnolake(jc, jb) = prm_field(1)%lsmask(jc,jb) + prm_field(1)%alake(jc,jb)
            ENDDO
        ENDDO
        !ICON_OMP_END_PARALLEL_DO

        mask_checksum = 0
        !ICON_OMP_PARALLEL_DO PRIVATE(jb,jc) REDUCTION(+:mask_checksum) ICON_OMP_RUNTIME_SCHEDULE
        DO jb = 1, patch_horz%nblks_c
          DO jc = 1, nproma
            mask_checksum = mask_checksum + ABS( lsmnolake(jc,jb))
          ENDDO
        ENDDO
        !ICON_OMP_END_PARALLEL_DO

        mask_checksum = p_sum(mask_checksum, comm=p_comm_work)

        !
        ! Define cell_mask_ids(1): all ocean and coastal points are valid
        !   This is the standard for the coupling of atmospheric fields listed below
        !
        IF ( mask_checksum > 0 ) THEN
          !ICON_OMP_PARALLEL_DO PRIVATE(jb, jc, nn) ICON_OMP_RUNTIME_SCHEDULE
          DO jb = 1, patch_horz%nblks_c
            DO jc = 1, nproma

              IF ( lsmnolake(jc, jb) .LT. (1.0_wp - eps) ) THEN
                ! ocean point (fraction of ocean is >0., lsmnolake .lt. 1.) is valid
                is_valid((jb-1)*nproma+jc) = .TRUE.
              ELSE
                ! land point (fraction of land is one, no sea water, lsmnolake=1.) is undef
                is_valid((jb-1)*nproma+jc) = .FALSE.
              ENDIF

            ENDDO
          ENDDO
          !ICON_OMP_END_PARALLEL_DO
        ELSE
          !ICON_OMP_PARALLEL_DO PRIVATE(jb, jc, nn) ICON_OMP_RUNTIME_SCHEDULE
          DO jc = 1, patch_horz%nblks_c * nproma
            is_valid(jc) = .TRUE.
          ENDDO
          !ICON_OMP_END_PARALLEL_DO
        ENDIF
#endif
       CASE DEFAULT

          CALL finish ('Please mask handling for new forcing in ' &
            & //'src/coupling/mo_atmo_coupling_frame: construct_atmo_coupling. Thank you!')

    END SELECT !}}}

    CALL yac_fdef_mask (          &
      & grid_id,                  &
      & patch_horz%n_patch_cells, &
      & YAC_LOCATION_CELL,        &
      & is_valid,                 &
      & cell_mask_ids(1) )

    CALL yac_fsync_def()

    DO jc = 1, min_no_of_fields
      CALL yac_fdef_field_mask ( &
        & TRIM(field_name(jc)),  &
        & comp_id,               &
        & cell_point_ids,        &
        & cell_mask_ids(1),      &
        & 1,                     &
        & collection_size(jc),   &
        & modelTimeStep,         &
        & YAC_TIME_UNIT_ISO_FORMAT, &
        & field_id(jc) )
    ENDDO

    ! add Ozone data field if needed
    needs_o3_coupling = (aes_rad_config(jg)%lrad_yac     .AND.  &
                        & (aes_rad_config(jg)%irad_o3 == 5 .OR. &
                        &  aes_rad_config(jg)%irad_o3 == 6))
    IF  (needs_o3_coupling) THEN
      ! enlarge the field_id array because it's public and ppl might want to use
      ! its shape info for looping
       ALLOCATE(tmp_field_id(min_no_of_fields+1))
       tmp_field_id(1:min_no_of_fields) = field_id(1:min_no_of_fields)
       CALL move_alloc(tmp_field_id, field_id)

       field_name(min_no_of_fields+1) = "o3"
       collection_size(min_no_of_fields+1) = &
         & yac_fget_field_collection_size("o3_provider", "o3_grid", field_name(min_no_of_fields+1) )
       CALL yac_fdef_field( &
         & TRIM(field_name(min_no_of_fields+1)),  &
         & comp_id,               &
         & cell_point_ids,        &
         & 1,                     &
         & collection_size(min_no_of_fields+1),   &
         & TRIM(aes_phy_config(jg)%dt_rad), &
         & YAC_TIME_UNIT_ISO_FORMAT, &
         & field_id(min_no_of_fields+1) )

    END IF

    needs_areo_coupling =  (aes_rad_config(jg)%lrad_yac        .AND.  &
                           & (aes_rad_config(jg)%irad_aero == 12 .OR. &
                           &  aes_rad_config(jg)%irad_aero == 13 .OR. &
                           &  aes_rad_config(jg)%irad_aero == 15 .OR. &
                           &  aes_rad_config(jg)%irad_aero == 18 .OR. &
                           &  aes_rad_config(jg)%irad_aero == 19))
    IF (needs_areo_coupling) THEN
      intermediate_no_of_fields = MERGE(min_no_of_fields+1, min_no_of_fields, needs_o3_coupling)
      field_name(intermediate_no_of_fields+1) = "aod_lw_b16_coa" !TODO 15, 16, ..., 24
      field_name(intermediate_no_of_fields+2) = "ssa_lw_b16_coa"
      field_name(intermediate_no_of_fields+3) = "aer_lw_b16_coa"
      field_name(intermediate_no_of_fields+4) = "aod_sw_b14_coa"
      field_name(intermediate_no_of_fields+5) = "ssa_sw_b14_coa"
      field_name(intermediate_no_of_fields+6) = "asy_sw_b14_coa"
      field_name(intermediate_no_of_fields+7) = "aod_sw_b14_fin"
      field_name(intermediate_no_of_fields+8) = "ssa_sw_b14_fin"
      field_name(intermediate_no_of_fields+9) = "asy_sw_b14_fin"
      field_name(intermediate_no_of_fields+10)= "aer_sw_b14_fin"

      ALLOCATE(tmp_field_id(max_no_of_fields))
      tmp_field_id(1:intermediate_no_of_fields) = field_id(1:intermediate_no_of_fields)
      CALL move_alloc(tmp_field_id, field_id)

      DO jc = intermediate_no_of_fields+1, max_no_of_fields
        collection_size(jc) = yac_fget_field_collection_size("aero_provider", "aero_grid", TRIM(field_name(jc)))
        CALL yac_fdef_field( &
          & TRIM(field_name(jc)),  &
          & comp_id,               &
          & cell_point_ids,        &
          & 1,                     &
          & collection_size(jc),   &
          & TRIM(aes_phy_config(jg)%dt_rad), &
          & YAC_TIME_UNIT_ISO_FORMAT, &
          & field_id(jc) )
      ENDDO

    END IF

#if !defined(__NO_JSBACH__) && !defined(__NO_JSBACH_HD__)

    ! Define coupling of runoff if HD model is present and interface is coded
    !  - discrimination between Proto2 (no HD) and Proto3 (with HD) is needed
    ! preliminary: coupling to jsbach/hd is active
    IF (iforcing /= INWP .OR. atm_phy_nwp_config(jg)%inwp_surface == LSS_JSBACH) THEN

!     !
!     ! Attention: needs to be checked with Roland Wirth and JSBACH users in case a second runoff mask is used
!     !            if a separate runoff mask is used, the dimension of cell_mask_ids is (2)
!     !
!     ! ! Define cell_mask_ids(2) for runoff:
!     ! !slo old!   Ocean coastal points with respect to HDmodel mask only are valid.
!     ! !slo old!   The integer mask for the HDmodel is ext_data(1)%atm%lsm_hd_c(:,:).
!     ! !slo old!   Caution: jg=1 is only valid for coupling to ocean
!     ! !
!     IF ( mask_checksum > 0 ) THEN
!
! !ICON_OMP_PARALLEL_DO PRIVATE(jb, jc, nn) ICON_OMP_RUNTIME_SCHEDULE
!         DO jb = 1, patch_horz%nblks_c
!           DO jc = 1, nproma
!
!              IF ( lsmnolake(jc, jb) .LT. (1.0_wp - eps) ) THEN
!                ! ocean point (fraction of ocean is >0., lsmnolake .lt. 1.) is valid
!                ! eps necessary because lsmnolake=lsm_ctr_c has input values of 0.999999 instead of 1.0)
!                is_valid((jb-1)*nproma+jc) = .TRUE.
!              ELSE
!                ! land point (fraction of land is one, lsmnolake=1.) is undef
!                is_valid((jb-1)*nproma+jc) = .FALSE.
!              ENDIF
!
!           ENDDO
!         ENDDO
! !ICON_OMP_END_PARALLEL_DO
!     ELSE
! !ICON_OMP_PARALLEL_DO PRIVATE(jb, jc, nn) ICON_OMP_RUNTIME_SCHEDULE
!        DO jc = 1,patch_horz%nblks_c * nproma
!           is_valid(jc) = .TRUE.
!        ENDDO
! !ICON_OMP_END_PARALLEL_DO
!
!     ENDIF
!
!     CALL yac_fdef_mask (          &
!       & grid_id,                  &
!       & patch_horz%n_patch_cells, &
!       & YAC_LOCATION_CELL,        &
!       & is_valid,                 &
!       & cell_mask_ids(2) )
!
!     ! Define additional coupling field(s) for JSBACH/HD
!     ! Utilize mask field for runoff
!     ! cell_mask_ids(2) shall contain ocean coast points only for source point mapping (source_to_target_map)
!     ! Currently it is the same mask as for the rest.

      ! change of call to be checked - no second mask used for this call
      !CALL jsb_fdef_hd_fields(comp_id, cell_point_ids, cell_mask_ids(2:2) )

      CALL jsb_fdef_hd_fields(comp_id, cell_point_ids, grid_id, patch_horz%n_patch_cells)

    ENDIF

#endif

    DEALLOCATE (is_valid)

    DEALLOCATE (lsmnolake)

    ! End definition of coupling fields and search

    CALL yac_fenddef ( )

    IF (ltimer) CALL timer_stop(timer_coupling_init)

  END SUBROUTINE construct_atmo_coupling

END MODULE mo_atmo_coupling_frame

!vim:fdm=marker
