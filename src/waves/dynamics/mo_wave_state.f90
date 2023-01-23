!>
!!        Contains the variables to set up the wave  model.
!=============================================================================================
!!
!! @author Mikhail Dobrynin, DWD, 11.06.2019
!!
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!
!=============================================================================================

MODULE mo_wave_state

  USE mo_exception,            ONLY: message, finish
  USE mo_parallel_config,      ONLY: nproma
  USE mo_wave_types,           ONLY: t_wave_diag, t_wave_state, t_wave_state_lists
  USE mo_model_domain,         ONLY: t_patch
  USE mo_grid_config,          ONLY: n_dom
  USE mo_impl_constants,       ONLY: max_char_length, SUCCESS
  USE mo_var_list,             ONLY: add_var, t_var_list_ptr
  USE mo_var_list_register,    ONLY: vlr_add, vlr_del
  USE mo_cdi_constants,        ONLY: GRID_UNSTRUCTURED_CELL, GRID_CELL
  USE mo_cdi,                  ONLY: DATATYPE_FLT32, DATATYPE_FLT64, GRID_UNSTRUCTURED, &
       &                             DATATYPE_PACK16
  USE mo_zaxis_type,           ONLY: ZA_SURFACE
  USE mo_cf_convention,        ONLY: t_cf_var
  USE mo_grib2,                ONLY: t_grib2_var, grib2_var
  USE mo_io_config,            ONLY: lnetcdf_flt64_output
  USE mo_mpi,                  ONLY: get_my_mpi_work_id


  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_wave_state'

  PUBLIC :: construct_wave_state    ! Constructor for the wave state
  PUBLIC :: destruct_wave_state     ! Destructor

  PUBLIC :: p_wave_state            ! state vector of wave variables
  PUBLIC :: p_wave_state_lists      ! lists for state vector of wave variables

  TYPE(t_wave_state),       TARGET, ALLOCATABLE :: p_wave_state(:)
  TYPE(t_wave_state_lists), TARGET, ALLOCATABLE :: p_wave_state_lists(:)

CONTAINS

  SUBROUTINE construct_wave_state(p_patch)

    TYPE(t_patch),      INTENT(IN) :: p_patch(:)

    CHARACTER(len=max_char_length) :: listname
    CHARACTER(len=*), PARAMETER :: routine = modname//'::construct_wave_state'
    INTEGER :: jg
    INTEGER :: ist

    ALLOCATE (p_wave_state(n_dom), stat=ist)
    IF (ist /= success) THEN
      CALL finish(TRIM(routine),'allocation for p_wave_state failed')
    END IF

    ALLOCATE (p_wave_state_lists(n_dom), stat=ist)
    IF (ist /= success) THEN
      CALL finish(TRIM(routine),'allocation for p_wave_state_lists failed')
    END IF

    DO jg = 1, n_dom
      WRITE(listname,'(a,i2.2)') 'wave_state_diag_of_domain_',jg
      CALL new_wave_state_diag_list(&
            p_patch(jg), &
            p_wave_state(jg)%diag, &
            p_wave_state_lists(jg)%diag_list, &
            listname)

      p_wave_state(jg)%diag%process_id(:,:) = get_my_mpi_work_id()
    END DO

    CALL message(routine, 'done.')

  END SUBROUTINE construct_wave_state




  SUBROUTINE new_wave_state_diag_list(p_patch, p_diag, p_diag_list, listname)

    TYPE(t_patch),         INTENT(IN)    :: p_patch
    TYPE(t_wave_diag),     INTENT(INOUT) :: p_diag
    TYPE(t_var_list_ptr),  INTENT(INOUT) :: p_diag_list
    CHARACTER(len=*),      INTENT(IN)    :: listname

    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: ibits         !< "entropy" of horizontal slice
    INTEGER :: datatype_flt  !< floating point accuracy in NetCDF output
    INTEGER :: nblks_c
    INTEGER :: shape2d_c(2)

    nblks_c = p_patch%nblks_c

    shape2d_c = (/nproma, nblks_c/)

    ibits = DATATYPE_PACK16   ! "entropy" of horizontal slice

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    CALL vlr_add(p_diag_list, TRIM(listname), patch_id=p_patch%id, lrestart=.TRUE. )

    cf_desc    = t_cf_var('process_id', 'i','process ID', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(&
         p_diag_list, 'process_id', p_diag%process_id, &
         GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         ldims=shape2d_c)

  END SUBROUTINE new_wave_state_diag_list


  !>
  !! Destruction of wave-specific variable lists and memory deallocation
  !!
  !! @par Revision History
  !! Initial revision by Mikhail Dobrynin, DWD (2022-12-20)
  !!
  SUBROUTINE destruct_wave_state ()

    INTEGER :: ist
    INTEGER :: jg
    CHARACTER(len=*), PARAMETER :: routine = modname//'::destruct_wave_state'

    DO jg = 1, n_dom

      ! delete diagnostics state list elements
      CALL vlr_del(p_wave_state_lists(jg)%diag_list)
    ENDDO

    ! destruct state list and state
    DEALLOCATE(p_wave_state, p_wave_state_lists, stat=ist)
    IF (ist /= success) THEN
      CALL finish(TRIM(routine),'deallocation for p_wave_state failed')
    END IF

  END SUBROUTINE destruct_wave_state


END MODULE mo_wave_state
