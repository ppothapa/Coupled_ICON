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

MODULE mo_ocean_coupling_frame

  USE mo_master_control,      ONLY: get_my_process_name
  USE mo_kind,                ONLY: wp
  USE mo_parallel_config,     ONLY: nproma
  USE mo_exception,           ONLY: warning, message, finish
  USE mo_impl_constants,      ONLY: max_char_length
  USE mo_mpi,                 ONLY: p_pe_work
  USE mo_run_config,          ONLY: ltimer
  USE mo_timer,               ONLY: timer_start, timer_stop, &
       &                            timer_coupling_init
  USE mo_model_domain,        ONLY: t_patch, t_patch_3d
  USE mtime,                  ONLY: datetimeToString, MAX_DATETIME_STR_LEN, &
       &                            timedeltaToString, MAX_TIMEDELTA_STR_LEN

  !-------------------------------------------------------------
  ! For the coupling
  !
  USE mo_math_constants,      ONLY: pi
  USE mo_parallel_config,     ONLY: nproma
#ifdef YAC_coupling
  USE mo_yac_finterface,      ONLY: yac_fdef_comp, yac_fget_version,        &
    &                               yac_fdef_datetime, yac_fdef_grid,       &
    &                               yac_fdef_points, yac_fset_global_index, &
    &                               yac_fset_core_mask, yac_fenddef,        &
    &                               YAC_LOCATION_CELL
#endif
  USE mo_ocean_atmo_coupling, ONLY: construct_ocean_atmo_coupling
  USE mo_coupling_config,     ONLY: is_coupled_run, is_coupled_to_atmo
  USE mo_time_config,         ONLY: time_config 

  !-------------------------------------------------------------

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: str_module = 'mo_ocean_coupling_frame' ! Output of module for debug

  PUBLIC :: construct_ocean_coupling, destruct_ocean_coupling
  PUBLIC :: nbr_inner_cells

  INTEGER, SAVE :: nbr_inner_cells

CONTAINS

  !--------------------------------------------------------------------------
  ! Prepare the coupling
  !
  ! For the time being this could all go into a subroutine which is
  ! common to atmo and ocean. Does this make sense if the setup deviates
  ! too much in future.
  !------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE construct_ocean_coupling(patch_3d)
    TYPE(t_patch_3d ), TARGET, INTENT(in)    :: patch_3d

    INTEGER                :: patch_no
    TYPE(t_patch), POINTER :: patch_horz

    CHARACTER(LEN=max_char_length) :: grid_name
    CHARACTER(LEN=max_char_length) :: comp_name

    INTEGER :: comp_id
    INTEGER :: cell_point_id
    INTEGER :: grid_id
    INTEGER :: nbr_vertices_per_cell

    INTEGER :: nblks
    INTEGER :: blockNo, cell_index, nn

    REAL(wp), ALLOCATABLE :: buffer_lon(:)
    REAL(wp), ALLOCATABLE :: buffer_lat(:)
    INTEGER, ALLOCATABLE  :: buffer_c(:,:)
    LOGICAL, ALLOCATABLE  :: is_valid(:)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: startdatestring
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: stopdatestring
    CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN) :: timestepstring

    CHARACTER(LEN=*), PARAMETER   :: routine = str_module // ':construct_ocean_coupling'

#ifndef YAC_coupling
    CALL finish(routine, 'built without coupling support.')
#else

    IF (.NOT. is_coupled_run()) RETURN

    IF (ltimer) CALL timer_start(timer_coupling_init)

    comp_name = TRIM(get_my_process_name())

    patch_no = 1
    patch_horz => patch_3d%p_patch_2d(patch_no)

    ! Inform the coupler about what we are
    CALL yac_fdef_comp ( TRIM(comp_name), comp_id )

    ! Print the YAC version
    CALL message('Running ICON ocean in coupled mode with YAC version ', TRIM(yac_fget_version()) )

    ! Overwrite job start and end date with component data
    CALL datetimeToString(time_config%tc_startdate, startdatestring)
    CALL datetimeToString(time_config%tc_stopdate, stopdatestring)
    CALL timedeltaToString(time_config%tc_dt_model, timestepstring)

    CALL yac_fdef_datetime ( start_datetime = TRIM(startdatestring), &
         &                   end_datetime   = TRIM(stopdatestring)   )

    ! Announce one subdomain (patch) to the coupler
    grid_name = "icon_ocean_grid"

    ! Extract cell information
    !
    ! cartesian coordinates of cell vertices are stored in
    ! patch_horz%verts%cartesian(:,:)%x(1:3)
    ! Here we use the longitudes and latitudes.

    nblks = max(patch_horz%nblks_c,patch_horz%nblks_v)

    ALLOCATE(buffer_lon(nproma*nblks))
    ALLOCATE(buffer_lat(nproma*nblks))
    ALLOCATE(buffer_c(3,nproma*nblks))

    nbr_vertices_per_cell = 3

!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(blockNo, cell_index, nn) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = 1, patch_horz%nblks_v
      DO cell_index = 1, nproma
        nn = (blockNo-1)*nproma+cell_index
        buffer_lon(nn) = patch_horz%verts%vertex(cell_index,blockNo)%lon
        buffer_lat(nn) = patch_horz%verts%vertex(cell_index,blockNo)%lat
      ENDDO
    ENDDO
!ICON_OMP_END_DO NOWAIT

!ICON_OMP_DO PRIVATE(blockNo, cell_index, nn) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = 1, patch_horz%nblks_c
      DO cell_index = 1, nproma
        nn = (blockNo-1)*nproma+cell_index
        buffer_c(1,nn) = (patch_horz%cells%vertex_blk(cell_index,blockNo,1)-1)*nproma + &
          &               patch_horz%cells%vertex_idx(cell_index,blockNo,1)
        buffer_c(2,nn) = (patch_horz%cells%vertex_blk(cell_index,blockNo,2)-1)*nproma + &
          &               patch_horz%cells%vertex_idx(cell_index,blockNo,2)
        buffer_c(3,nn) = (patch_horz%cells%vertex_blk(cell_index,blockNo,3)-1)*nproma + &
                          patch_horz%cells%vertex_idx(cell_index,blockNo,3)
      ENDDO
    ENDDO
!ICON_OMP_END_DO
!ICON_OMP_END_PARALLEL

    ! Description of elements, here as unstructured grid
    CALL yac_fdef_grid (          &
      & TRIM(grid_name),          &
      & patch_horz%n_patch_verts, &
      & patch_horz%n_patch_cells, &
      & nbr_vertices_per_cell,    &
      & buffer_lon,               &
      & buffer_lat,               &
      & buffer_c,                 &
      & grid_id )

    ! Can we have two fdef_point calls for the same subdomain, i.e.
    ! one single set of cells?
    !
    ! Define cell center points (location = 0)
    !
    ! cartesian coordinates of cell centers are stored in
    ! patch_horz%cells%cartesian_center(:,:)%x(1:3)
    ! Here we use the longitudes and latitudes.

!ICON_OMP_PARALLEL_DO PRIVATE(blockNo, cell_index, nn) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = 1, patch_horz%nblks_c
      DO cell_index = 1, nproma
        nn = (blockNo-1)*nproma+cell_index
        buffer_lon(nn) = patch_horz%cells%center(cell_index,blockNo)%lon
        buffer_lat(nn) = patch_horz%cells%center(cell_index,blockNo)%lat
      ENDDO
    ENDDO
!ICON_OMP_END_PARALLEL_DO

    ! center points in cells (needed e.g. for patch recovery and nearest neighbour)
    CALL yac_fdef_points (        &
      & grid_id,                  &
      & patch_horz%n_patch_cells, &
      & YAC_LOCATION_CELL,        &
      & buffer_lon,               &
      & buffer_lat,               &
      & cell_point_id )

    DEALLOCATE (buffer_lon, buffer_lat, buffer_c)

    CALL yac_fset_global_index (                &
      & patch_horz%cells%decomp_info%glb_index, &
      & YAC_LOCATION_CELL,                      &
      & grid_id )

    ALLOCATE(is_valid(nproma*patch_horz%nblks_c))

    nbr_inner_cells = 0
!ICON_OMP_PARALLEL DO PRIVATE(cell_index) REDUCTION(+:nbr_inner_cells) ICON_OMP_DEFAULT_SCHEDULE
    DO cell_index = 1, patch_horz%n_patch_cells
       IF ( p_pe_work == patch_horz%cells%decomp_info%owner_local(cell_index) ) THEN
         is_valid(cell_index) =.TRUE.
         nbr_inner_cells = nbr_inner_cells + 1
       ELSE
         is_valid(cell_index) = .FALSE.
       ENDIF
    ENDDO
!ICON_OMP_END_PARALLEL_DO

    CALL yac_fset_core_mask ( &
      & is_valid,             &
      & YAC_LOCATION_CELL,    &
      & grid_id )

    DEALLOCATE(is_valid)

    IF ( is_coupled_to_atmo() ) THEN

      ! Construct coupling frame for ocean-atmosphere
      CALL message(str_module, 'Constructing the coupling frame ocean-atmosphere.')

      CALL construct_ocean_atmo_coupling( &
        patch_3d, comp_id, grid_id, cell_point_id, timestepstring, &
        nbr_inner_cells)

    END IF

    CALL yac_fenddef ( )

    IF (ltimer) CALL timer_stop(timer_coupling_init)

! YAC_coupling
#endif

  END SUBROUTINE construct_ocean_coupling

  !--------------------------------------------------------------------------

!<Optimize:inUse>
  SUBROUTINE destruct_ocean_coupling()

    IF (.NOT. is_coupled_run()) RETURN

  END SUBROUTINE destruct_ocean_coupling

  !--------------------------------------------------------------------------

END MODULE mo_ocean_coupling_frame

