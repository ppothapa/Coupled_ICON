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
  USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config

  USE mo_parallel_config     ,ONLY: nproma

  USE mo_run_config          ,ONLY: iforcing, ltimer
  USE mo_timer,               ONLY: timer_start, timer_stop, timer_coupling_init
  USE mo_impl_constants      ,ONLY: MAX_CHAR_LENGTH, inwp, LSS_JSBACH, SUCCESS

#if !defined(__NO_JSBACH__) && !defined(__NO_JSBACH_HD__)
  USE mo_interface_hd_ocean  ,ONLY: jsb_fdef_hd_fields
#endif

  USE mo_master_control      ,ONLY: get_my_process_name

  USE mo_mpi                 ,ONLY: p_pe_work

  USE mo_coupling_config     ,ONLY: is_coupled_run, is_coupled_to_ocean, &
    &                               is_coupled_to_hydrodisc, &
    &                               is_coupled_to_waves, is_coupled_to_output
  USE mo_aes_rad_config      ,ONLY: aes_rad_config
  USE mo_aes_phy_config      ,ONLY: aes_phy_config
  USE mo_time_config         ,ONLY: time_config

  USE mo_atmo_wave_coupling  ,ONLY: construct_atmo_wave_coupling
  USE mo_atmo_ocean_coupling ,ONLY: construct_atmo_ocean_coupling
  USE mo_atmo_o3_provider_coupling,ONLY: &
    construct_atmo_o3_provider_coupling_post_sync
  USE mo_atmo_aero_provider_coupling,ONLY: &
    construct_atmo_aero_provider_coupling_post_sync
  USE mo_nwp_hydrodisc_coupling,ONLY: construct_nwp_hydrodisc_coupling

  USE mo_exception           ,ONLY: finish, message

#ifdef YAC_coupling
  USE mo_yac_finterface      ,ONLY: yac_fget_version, yac_fdef_comp,        &
    &                               yac_fdef_comps,                         &
    &                               yac_fdef_datetime, yac_fdef_grid,       &
    &                               yac_fdef_points, yac_fset_global_index, &
    &                               yac_fset_core_mask, yac_fsync_def,      &
    &                               yac_fenddef, YAC_LOCATION_CELL,         &
    &                               YAC_LOCATION_CORNER
#endif

  USE mtime                  ,ONLY: datetimeToString, MAX_DATETIME_STR_LEN, &
    &                               timedeltaToString, MAX_TIMEDELTA_STR_LEN

#ifndef __NO_ICON_COMIN__
  USE comin_host_interface, ONLY: EP_ATM_YAC_DEFCOMP_BEFORE,       &
       &                          EP_ATM_YAC_DEFCOMP_AFTER,        &
       &                          EP_ATM_YAC_SYNCDEF_BEFORE,       &
       &                          EP_ATM_YAC_SYNCDEF_AFTER,        &
       &                          EP_ATM_YAC_ENDDEF_BEFORE,        &
       &                          EP_ATM_YAC_ENDDEF_AFTER,         &
       &                          COMIN_DOMAIN_OUTSIDE_LOOP
  USE mo_comin_adapter,     ONLY: icon_call_callback
#endif

  USE mo_output_coupling     ,ONLY: construct_output_coupling, winnow_field_list

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: str_module = 'mo_atmo_coupling_frame' ! Output of module for debug

  PUBLIC :: construct_atmo_coupling
  PUBLIC :: nbr_inner_cells

  INTEGER, SAVE         :: nbr_inner_cells

CONTAINS

  !>
  !! SUBROUTINE construct_atmo_coupling -- the initialisation for the coupling
  !! of atmosphere and the ocean, through a coupler

  SUBROUTINE construct_atmo_coupling (p_patch)

    TYPE(t_patch), TARGET, INTENT(IN) :: p_patch(:)

    TYPE(t_patch), POINTER :: patch_horz

    CHARACTER(LEN=max_char_length) :: comp_name

    !---------------------------------------------------------------------
    ! 11. Do the setup for the coupled run
    !
    ! For the time being this could all go into a subroutine which is
    ! common to atmo and ocean. Does this make sense if the setup deviates
    ! too much in future.
    !---------------------------------------------------------------------

    CHARACTER(LEN=max_char_length) :: grid_name

    INTEGER :: comp_id, output_comp_id, comp_ids(2)
    INTEGER :: cell_point_id, vertex_point_id
    INTEGER :: grid_id

    INTEGER :: jg
    INTEGER :: nblks
    INTEGER :: jb, jc, nn
    INTEGER :: nbr_vertices_per_cell
    INTEGER :: error

    REAL(wp), ALLOCATABLE :: buffer_lon(:)
    REAL(wp), ALLOCATABLE :: buffer_lat(:)
    INTEGER,  ALLOCATABLE :: buffer_c(:,:)

    LOGICAL,  ALLOCATABLE :: is_valid(:)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: startdatestring
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: stopdatestring
    CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN) :: timestepstring

    CHARACTER(LEN=*), PARAMETER   :: routine = str_module // ':construct_atmo_coupling'

#ifndef YAC_coupling
    CALL finish(routine, 'built without coupling support.')
#else

    IF ( .NOT. is_coupled_run() ) RETURN

    IF (ltimer) CALL timer_start (timer_coupling_init)

    comp_name = TRIM(get_my_process_name())

    jg = 1
    patch_horz => p_patch(jg)

    ! Print the YAC version
    CALL message('Running ICON atmosphere in coupled mode with YAC version ', TRIM(yac_fget_version()) )

    ! Overwrite job start and end date with component data
    CALL datetimeToString(time_config%tc_startdate, startdatestring)
    CALL datetimeToString(time_config%tc_stopdate, stopdatestring)
    CALL timedeltaToString(time_config%tc_dt_model, timestepstring)

    CALL yac_fdef_datetime ( start_datetime = TRIM(startdatestring), &
         &                   end_datetime   = TRIM(stopdatestring)   )

#ifndef __NO_ICON_COMIN__
    CALL icon_call_callback(EP_ATM_YAC_DEFCOMP_BEFORE, COMIN_DOMAIN_OUTSIDE_LOOP)
#endif

    ! Inform the coupler about what we are
    IF( is_coupled_to_output() ) THEN
      CALL yac_fdef_comps (          &
        [TRIM(comp_name)//"       ", &
         TRIM(comp_name)//"_output"], 2, comp_ids )
      comp_id = comp_ids(1)
      output_comp_id = comp_ids(2)
    ELSE
      CALL yac_fdef_comp ( TRIM(comp_name), comp_id )
    ENDIF

#ifndef __NO_ICON_COMIN__
    CALL icon_call_callback(EP_ATM_YAC_DEFCOMP_AFTER, COMIN_DOMAIN_OUTSIDE_LOOP)
#endif

    ! Announce one grid (patch) to the coupler
    grid_name = "icon_atmos_grid"

    ! Extract cell information {{{
    !
    ! cartesian coordinates of cell vertices are stored in
    ! patch_horz%verts%cartesian(:,:)%x(1:3)
    ! Here we use the longitudes and latitudes in rad.

    nblks = max(patch_horz%nblks_c,patch_horz%nblks_v)

    ALLOCATE(buffer_lon(nproma*nblks), STAT = error)
    IF(error /= SUCCESS) CALL finish(str_module, "memory allocation failure for buffer_lon")
    ALLOCATE(buffer_lat(nproma*nblks), STAT = error)
    IF(error /= SUCCESS) CALL finish(str_module, "memory allocation failure for buffer_lat")
    ALLOCATE(buffer_c(3,nproma*nblks), STAT = error)
    IF(error /= SUCCESS) CALL finish(str_module, "memory allocation failure for buffer_c")

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

    ! vertex points (needed for output_coupling)
    CALL yac_fdef_points (        &
      & grid_id,                  &
      & patch_horz%n_patch_verts, &
      & YAC_LOCATION_CORNER,      &
      & buffer_lon(1:patch_horz%n_patch_verts), &
      & buffer_lat(1:patch_horz%n_patch_verts), &
      & vertex_point_id )

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
      & cell_point_id )

    DEALLOCATE (buffer_lon, buffer_lat, buffer_c, STAT = error)
    IF(error /= SUCCESS) CALL finish(str_module, "Deallocation failed for buffer_lon, buffer_lat, buffer_c")

    CALL yac_fset_global_index (                &
      & patch_horz%cells%decomp_info%glb_index, &
      & YAC_LOCATION_CELL,                      &
      & grid_id )

    ALLOCATE(is_valid(nproma*patch_horz%nblks_c), STAT = error)
    IF(error /= SUCCESS) CALL finish(str_module, "memory allocation failure for is_valid")

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

    DEALLOCATE (is_valid, STAT = error)
    IF(error /= SUCCESS) CALL finish(str_module, "Deallocation failed for is_valid")

    IF( is_coupled_to_output() ) THEN

      CALL message(str_module, 'Constructing the coupling frame atmosphere-output.')

      CALL construct_output_coupling ( &
        p_patch, output_comp_id, cell_point_id, vertex_point_id, timestepstring)

    END IF

    IF ( is_coupled_to_ocean() ) THEN

      CALL message(str_module, 'Constructing the coupling frame atmosphere-ocean.')

      CALL construct_atmo_ocean_coupling( &
        p_patch, comp_id, grid_id, cell_point_id, timestepstring)

#if !defined(__NO_JSBACH__) && !defined(__NO_JSBACH_HD__)

      ! Define coupling of runoff if HD model is present and interface is coded
      !  - discrimination between Proto2 (no HD) and Proto3 (with HD) is needed
      ! preliminary: coupling to jsbach/hd is active
      IF ( iforcing /= INWP .OR. ( atm_phy_nwp_config(jg)%inwp_surface == LSS_JSBACH .AND. .NOT. is_coupled_to_hydrodisc() ) ) THEN

        ! Construct coupling frame for atmosphere/JSBACH-hydrological discharge
        CALL message(str_module, 'Constructing the coupling frame atmosphere/JSBACH-hydrological discharge.')

        CALL jsb_fdef_hd_fields(comp_id, (/cell_point_id/), grid_id, patch_horz%n_patch_cells)

      ENDIF

#endif
    ENDIF   ! Construct coupling frame for atmosphere-ocean

    IF ( is_coupled_to_hydrodisc() ) THEN

      ! Construct coupling frame for atmosphere-hydrological discharge
      CALL message(str_module, 'Constructing the coupling frame atmosphere-hydrological discharge.')


      CALL construct_nwp_hydrodisc_coupling( &
        p_patch, comp_id, grid_id, cell_point_id, timestepstring)

    ENDIF ! Construct coupling frame for atmosphere-hydrological discharge

    IF ( is_coupled_to_waves() ) THEN

      CALL message(str_module, 'Constructing the coupling frame atmosphere-wave.')

      CALL construct_atmo_wave_coupling( &
        comp_id, cell_point_id, timestepstring)

    END IF

#ifndef __NO_ICON_COMIN__
    CALL icon_call_callback(EP_ATM_YAC_SYNCDEF_BEFORE, COMIN_DOMAIN_OUTSIDE_LOOP)
#endif

    ! Synchronize all definitions until this point with other components
    CALL yac_fsync_def()

#ifndef __NO_ICON_COMIN__
    CALL icon_call_callback(EP_ATM_YAC_SYNCDEF_AFTER, COMIN_DOMAIN_OUTSIDE_LOOP)
#endif

    ! add Ozone data field if needed
    IF  ((aes_rad_config(jg)%lrad_yac     .AND.  &
          (aes_rad_config(jg)%irad_o3 == 5 .OR. &
          aes_rad_config(jg)%irad_o3 == 6))) THEN

      ! Construct coupling frame for atmosphere-o3 provider
      CALL message(str_module, 'Constructing the coupling frame atmosphere-o3 provider.')

      CALL construct_atmo_o3_provider_coupling_post_sync( &
        comp_id, cell_point_id, TRIM(aes_phy_config(jg)%dt_rad))

    END IF

    IF ((aes_rad_config(jg)%lrad_yac        .AND.  &
        (aes_rad_config(jg)%irad_aero == 12 .OR. &
          aes_rad_config(jg)%irad_aero == 13 .OR. &
          aes_rad_config(jg)%irad_aero == 15 .OR. &
          aes_rad_config(jg)%irad_aero == 18 .OR. &
          aes_rad_config(jg)%irad_aero == 19))) THEN

      ! Construct coupling frame for atmosphere-aero provider
      CALL message(str_module, 'Constructing the coupling frame atmosphere-aero provider.')

      CALL construct_atmo_aero_provider_coupling_post_sync( &
        comp_id, cell_point_id, TRIM(aes_phy_config(jg)%dt_rad))

    END IF

    ! End definition of coupling fields and search

#ifndef __NO_ICON_COMIN__
    CALL icon_call_callback(EP_ATM_YAC_ENDDEF_BEFORE, COMIN_DOMAIN_OUTSIDE_LOOP)
#endif

    CALL yac_fenddef ( )

#ifndef __NO_ICON_COMIN__
    CALL icon_call_callback(EP_ATM_YAC_ENDDEF_AFTER, COMIN_DOMAIN_OUTSIDE_LOOP)
#endif

    IF( is_coupled_to_output() ) CALL winnow_field_list()

    IF (ltimer) CALL timer_stop(timer_coupling_init)

! YAC_coupling
#endif

  END SUBROUTINE construct_atmo_coupling

END MODULE mo_atmo_coupling_frame

!vim:fdm=marker
