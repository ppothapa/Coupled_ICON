!>
!! This module contains often-used Fortran language constructs.
!!
!! The small functions and subroutines in this module should depend
!! only on most basic types and should not call other model-specific
!! subroutines.
!!
!! @par Revision History
!!    Initial implementation : F. Prill, DWD (2012-07-04)
!!    moved routines from "mo_var_list.f90"
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
MODULE mo_fortran_tools

  USE mo_kind,                    ONLY: wp, sp, vp, dp, ik4 => i4
  USE mo_exception,               ONLY: finish
  USE mo_impl_constants,          ONLY: SUCCESS, vname_len
#ifdef _OPENACC
  USE mo_mpi,                     ONLY: i_am_accel_node
  USE openacc
#endif
  USE iso_c_binding,              ONLY: c_ptr, c_f_pointer, c_loc, c_null_ptr
  USE mo_util_stride,             ONLY: util_stride_1d, util_stride_2d

  IMPLICIT NONE

  PUBLIC :: assign_if_present
  PUBLIC :: t_ptr_2d3d, t_ptr_2d3d_vp
  PUBLIC :: assign_if_present_allocatable
  PUBLIC :: t_ptr_1d, t_ptr_1d_sp, t_ptr_1d_int
  PUBLIC :: t_ptr_1d_ptr_1d
  PUBLIC :: t_ptr_2d, t_ptr_2d_sp, t_ptr_2d_int
  PUBLIC :: t_ptr_3d, t_ptr_3d_sp, t_ptr_3d_int
  PUBLIC :: t_ptr_i2d3d
  PUBLIC :: t_ptr_4d, t_ptr_4d_sp, t_ptr_4d_int
  PUBLIC :: t_ptr_tracer
  PUBLIC :: copy, init, swap, negative2zero
  PUBLIC :: var_scale, var_add
  PUBLIC :: init_zero_contiguous_dp, init_zero_contiguous_sp
  PUBLIC :: init_contiguous_dp, init_contiguous_sp
  PUBLIC :: init_contiguous_i4, init_contiguous_l
  PUBLIC :: minval_1d
  PUBLIC :: resize_arr_c1d
  PUBLIC :: DO_DEALLOCATE
  PUBLIC :: DO_PTR_DEALLOCATE
  PUBLIC :: insert_dimension
  PUBLIC :: assert_acc_host_only
  PUBLIC :: assert_acc_device_only
  PUBLIC :: assert_lacc_equals_i_am_accel_node
  PUBLIC :: set_acc_host_or_device

  PRIVATE

  TYPE t_ptr_1d
    REAL(wp), POINTER :: p(:)  ! pointer to 1D (spatial) array
  END TYPE t_ptr_1d

  TYPE t_ptr_1d_sp
    REAL(sp), POINTER :: p(:)  ! pointer to 1D (spatial) array
  END TYPE t_ptr_1d_sp

  TYPE t_ptr_1d_int
    INTEGER, POINTER :: p(:)  ! pointer to 1D (spatial) array
  END TYPE t_ptr_1d_int

  TYPE t_ptr_1d_ptr_1d
    TYPE(t_ptr_1d), POINTER :: p(:)  ! pointer to a 1D array of pointers to 1D (spatial) arrays
  END TYPE t_ptr_1d_ptr_1d

  TYPE t_ptr_2d
    REAL(dp), POINTER :: p(:,:)  ! pointer to 2D (spatial) array
  END TYPE t_ptr_2d

  TYPE t_ptr_2d_sp
    REAL(sp), POINTER :: p(:,:)  ! pointer to 2D (spatial) array
  END TYPE t_ptr_2d_sp

  TYPE t_ptr_2d_int
    INTEGER, POINTER :: p(:,:)  ! pointer to 2D (spatial) array
  END TYPE t_ptr_2d_int

  TYPE t_ptr_3d
    REAL(dp), POINTER :: p(:,:,:)  ! pointer to 3D (spatial) array
  END TYPE t_ptr_3d

  TYPE t_ptr_3d_sp
    REAL(sp), POINTER :: p(:,:,:)  ! pointer to 3D (spatial) array
  END TYPE t_ptr_3d_sp

  TYPE t_ptr_3d_int
    INTEGER, POINTER :: p(:,:,:)  ! pointer to 3D (spatial) array
  END TYPE t_ptr_3d_int

  TYPE t_ptr_4d
    REAL(dp), POINTER :: p(:,:,:,:)  ! pointer to 3D (spatial) array
  END TYPE t_ptr_4d

  TYPE t_ptr_4d_sp
    REAL(sp), POINTER :: p(:,:,:,:)  ! pointer to 3D (spatial) array
  END TYPE t_ptr_4d_sp

  TYPE t_ptr_4d_int
    INTEGER, POINTER :: p(:,:,:,:)  ! pointer to 3D (spatial) array
  END TYPE t_ptr_4d_int

  TYPE t_ptr_2d3d
    REAL(wp),POINTER :: p_3d(:,:,:)  ! REAL pointer to 3D (spatial) array
    REAL(wp),POINTER :: p_2d(:,:)    ! REAL pointer to 2D (spatial) array
  END TYPE t_ptr_2d3d

  TYPE t_ptr_2d3d_vp
    REAL(vp),POINTER :: p_3d(:,:,:)  ! REAL pointer to 3D (spatial) array
    REAL(vp),POINTER :: p_2d(:,:)    ! REAL pointer to 2D (spatial) array
  END TYPE t_ptr_2d3d_vp

  TYPE t_ptr_i2d3d
    INTEGER,POINTER :: p_3d(:,:,:)  ! INTEGER pointer to 3D (spatial) array
    INTEGER,POINTER :: p_2d(:,:)    ! INTEGER pointer to 2D (spatial) array
  END TYPE t_ptr_i2d3d

  ! Type to pass pointer arrays to convection and turbulent diffusion subroutines
  TYPE t_ptr_tracer
    REAL(wp), POINTER :: ptr(:,:)
    INTEGER           :: idx_tracer
  END TYPE t_ptr_tracer

  INTERFACE assign_if_present
    MODULE PROCEDURE assign_if_present_character
    MODULE PROCEDURE assign_if_present_logical
    MODULE PROCEDURE assign_if_present_logicals
    MODULE PROCEDURE assign_if_present_integer
    MODULE PROCEDURE assign_if_present_integers
    MODULE PROCEDURE assign_if_present_real
    MODULE PROCEDURE assign_if_present_real_sp
  END INTERFACE assign_if_present

  INTERFACE assign_if_present_allocatable
    MODULE PROCEDURE assign_if_present_logical_allocatable_1d
    MODULE PROCEDURE assign_if_present_integer_allocatable
    MODULE PROCEDURE assign_if_present_integer_allocatable_1d
    MODULE PROCEDURE assign_if_present_real_allocatable
    MODULE PROCEDURE assign_if_present_real_allocatable_1d
    MODULE PROCEDURE assign_if_present_character_allocatable
  END INTERFACE assign_if_present_allocatable

  !> `copy(b, a)` is meant to make it easier for compilers to circumvent
  !! temporaries as are too often created in a(:, :, :) = b(:, :, :)
  !!
  !! `copy` uses openMP orphaning, i.e. it must be called inside an
  !! OMP PARALLEL region. However, it must not be called inside another
  !! OMP DO region. 
  INTERFACE copy
    MODULE PROCEDURE copy_2d_dp
    MODULE PROCEDURE copy_3d_dp
    MODULE PROCEDURE copy_4d_dp
    MODULE PROCEDURE copy_5d_dp
    MODULE PROCEDURE copy_5d_sp
    MODULE PROCEDURE copy_2d_spdp
    MODULE PROCEDURE copy_3d_spdp
    MODULE PROCEDURE copy_4d_spdp
    MODULE PROCEDURE copy_5d_spdp
    MODULE PROCEDURE copy_2d_i4
    MODULE PROCEDURE copy_3d_i4
    MODULE PROCEDURE copy_5d_i4
    MODULE PROCEDURE copy_5d_l
  END INTERFACE copy

  !> `init` uses openMP orphaning (explanation see `copy`)
  INTERFACE init
    MODULE PROCEDURE init_zero_1d_dp
    MODULE PROCEDURE init_zero_1d_sp
    MODULE PROCEDURE init_zero_2d_dp
    MODULE PROCEDURE init_zero_2d_i4
    MODULE PROCEDURE init_zero_3d_dp
    MODULE PROCEDURE init_zero_3d_sp
    MODULE PROCEDURE init_zero_3d_i4
    MODULE PROCEDURE init_zero_4d_i4
    MODULE PROCEDURE init_zero_4d_dp
    MODULE PROCEDURE init_zero_4d_sp
    MODULE PROCEDURE init_2d_dp
    MODULE PROCEDURE init_3d_dp
    MODULE PROCEDURE init_3d_spdp
    MODULE PROCEDURE init_5d_dp
    MODULE PROCEDURE init_5d_sp
    MODULE PROCEDURE init_5d_i4
    MODULE PROCEDURE init_5d_l
  END INTERFACE init

  INTERFACE negative2zero
    MODULE PROCEDURE negative2zero_4d_dp
  END INTERFACE negative2zero

  INTERFACE var_scale
    MODULE PROCEDURE var_scale_3d_dp
  END INTERFACE var_scale

  INTERFACE var_add
    MODULE PROCEDURE var_addc_3d_dp
  END INTERFACE var_add

  INTERFACE swap
    MODULE PROCEDURE swap_int
  END INTERFACE swap

  ! auxiliary routines
  INTERFACE DO_DEALLOCATE
    MODULE PROCEDURE DO_DEALLOCATE_r4D
    MODULE PROCEDURE DO_DEALLOCATE_r3D
    MODULE PROCEDURE DO_DEALLOCATE_r2D
    MODULE PROCEDURE DO_DEALLOCATE_r1D
    MODULE PROCEDURE DO_DEALLOCATE_i3D
    MODULE PROCEDURE DO_DEALLOCATE_i2D
    MODULE PROCEDURE DO_DEALLOCATE_i1D
  END INTERFACE DO_DEALLOCATE

  INTERFACE DO_PTR_DEALLOCATE
    MODULE PROCEDURE DO_PTR_DEALLOCATE_r3D
    MODULE PROCEDURE DO_PTR_DEALLOCATE_r2D
    MODULE PROCEDURE DO_PTR_DEALLOCATE_dp1D
    MODULE PROCEDURE DO_PTR_DEALLOCATE_sp1D
    MODULE PROCEDURE DO_PTR_DEALLOCATE_int1D
  END INTERFACE DO_PTR_DEALLOCATE


  CHARACTER(LEN = *), PARAMETER :: modname = "mo_fortran_tools"

  INTERFACE insert_dimension
    MODULE PROCEDURE insert_dimension_r_dp_3_2, insert_dimension_r_dp_3_2_s
    MODULE PROCEDURE insert_dimension_r_sp_3_2, insert_dimension_r_sp_3_2_s
    MODULE PROCEDURE insert_dimension_i4_3_2, insert_dimension_i4_3_2_s
    MODULE PROCEDURE insert_dimension_l_3_2, insert_dimension_l_3_2_s
    MODULE PROCEDURE insert_dimension_r_dp_6_5, insert_dimension_r_dp_6_5_s
    MODULE PROCEDURE insert_dimension_r_sp_6_5, insert_dimension_r_sp_6_5_s
    MODULE PROCEDURE insert_dimension_i4_6_5, insert_dimension_i4_6_5_s
  END INTERFACE insert_dimension

CONTAINS

  !
  ! private helper OpenACC function
  !
  SUBROUTINE acc_wait_if_requested(acc_async_queue, opt_acc_async)
    INTEGER, INTENT(IN) :: acc_async_queue
    LOGICAL, INTENT(IN), OPTIONAL :: opt_acc_async

#ifdef _OPENACC
    IF ( PRESENT(opt_acc_async) ) THEN
      IF ( .NOT. opt_acc_async ) THEN
        !$ACC WAIT(acc_async_queue)
      END IF
    ELSE
      !$ACC WAIT(acc_async_queue)
    END IF
#endif
  END SUBROUTINE acc_wait_if_requested

  ! routines to assign values if actual parameters are present
  !
  SUBROUTINE assign_if_present_character (y,x)
    CHARACTER(len=*), INTENT(inout)        :: y
    CHARACTER(len=*), INTENT(in) ,OPTIONAL :: x
    IF (.NOT. PRESENT(x)) RETURN
    IF ( x == ' ' )       RETURN
    y = x
  END SUBROUTINE assign_if_present_character

  SUBROUTINE assign_if_present_logical (y,x)
    LOGICAL, INTENT(inout)        :: y
    LOGICAL, INTENT(in) ,OPTIONAL :: x
    IF (PRESENT(x)) y = x
  END SUBROUTINE assign_if_present_logical

  SUBROUTINE assign_if_present_logicals (y,x)
    LOGICAL, INTENT(inout)        :: y(:)
    LOGICAL, INTENT(in) ,OPTIONAL :: x(:)
    INTEGER :: n
    IF (PRESENT(x)) THEN
      n = MIN(SIZE(x), SIZE(y))
      y(1:n) = x(1:n)
    ENDIF
  END SUBROUTINE assign_if_present_logicals

  SUBROUTINE assign_if_present_integer (y,x)
    INTEGER, INTENT(inout)        :: y
    INTEGER, INTENT(in) ,OPTIONAL :: x
    IF (.NOT. PRESENT(x)) RETURN
    IF ( x == -HUGE(x)  ) RETURN
    y = x
  END SUBROUTINE assign_if_present_integer

  SUBROUTINE assign_if_present_integers (y,x)
    INTEGER, INTENT(inout)        :: y (:)
    INTEGER, INTENT(in) ,OPTIONAL :: x (:)
    INTEGER :: n
    IF (PRESENT(x)) THEN
      n = MIN(SIZE(x), SIZE(y))
      y(1:n) = x(1:n)
    ENDIF
  END SUBROUTINE assign_if_present_integers

  SUBROUTINE assign_if_present_real (y,x)
    REAL(wp), INTENT(inout)        :: y
    REAL(wp), INTENT(in) ,OPTIONAL :: x
    IF (.NOT.PRESENT(x)) RETURN
    IF ( x == -HUGE(x) ) RETURN
    y = x
  END SUBROUTINE assign_if_present_real

  SUBROUTINE assign_if_present_real_sp (y,x)
    REAL(sp), INTENT(inout)        :: y
    REAL(sp), INTENT(in) ,OPTIONAL :: x
    IF (.NOT.PRESENT(x)) RETURN
    IF ( x == -HUGE(x) ) RETURN
    y = x
  END SUBROUTINE assign_if_present_real_sp

  SUBROUTINE assign_if_present_logical_allocatable_1d(y, x)
    LOGICAL, ALLOCATABLE, INTENT(INOUT) :: y(:)
    LOGICAL, OPTIONAL, INTENT(IN) :: x(:)

    INTEGER :: error
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":assign_if_present_logical_allocatable_1d"

    IF(.NOT.PRESENT(x)) RETURN
    IF(ALLOCATED(y)) THEN
        IF(SIZE(y) /= SIZE(x)) DEALLOCATE(y)
    END IF
    IF(.NOT.ALLOCATED(y)) THEN
        ALLOCATE(y(SIZE(x)), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
    END IF
    y(:) = x(:)
  END SUBROUTINE assign_if_present_logical_allocatable_1d

  SUBROUTINE assign_if_present_integer_allocatable(y, x)
    INTEGER, ALLOCATABLE, INTENT(INOUT) :: y
    INTEGER, OPTIONAL, INTENT(IN) :: x

    INTEGER :: error
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":assign_if_present_integer_allocatable"

    IF(.NOT.PRESENT(x)) RETURN
    IF(.NOT.ALLOCATED(y)) THEN
        ALLOCATE(y, STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
    END IF
    y = x
  END SUBROUTINE assign_if_present_integer_allocatable

  SUBROUTINE assign_if_present_integer_allocatable_1d(y, x)
    INTEGER, ALLOCATABLE, INTENT(INOUT) :: y(:)
    INTEGER, OPTIONAL, INTENT(IN) :: x(:)

    INTEGER :: error
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":assign_if_present_integer_allocatable_1d"

    IF(.NOT.PRESENT(x)) RETURN
    IF(ALLOCATED(y)) THEN
        IF(SIZE(y) /= SIZE(x)) DEALLOCATE(y)
    END IF
    IF(.NOT.ALLOCATED(y)) THEN
        ALLOCATE(y(SIZE(x)), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
    END IF
    y(:) = x(:)
  END SUBROUTINE assign_if_present_integer_allocatable_1d

  SUBROUTINE assign_if_present_real_allocatable(y, x)
    REAL(wp), ALLOCATABLE, INTENT(INOUT) :: y
    REAL(wp), OPTIONAL, INTENT(IN) :: x

    INTEGER :: error
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":assign_if_present_real_allocatable"

    IF(.NOT.PRESENT(x)) RETURN
    IF(.NOT.ALLOCATED(y)) THEN
        ALLOCATE(y, STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
    END IF
    y = x
  END SUBROUTINE assign_if_present_real_allocatable

  SUBROUTINE assign_if_present_real_allocatable_1d(y, x)
    REAL(wp), ALLOCATABLE, INTENT(INOUT) :: y(:)
    REAL(wp), OPTIONAL, INTENT(IN) :: x(:)

    INTEGER :: error
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":assign_if_present_real_allocatable_1d"

    IF(.NOT.PRESENT(x)) RETURN
    IF(ALLOCATED(y)) THEN
        IF(SIZE(y) /= SIZE(x)) DEALLOCATE(y)
    END IF
    IF(.NOT.ALLOCATED(y)) THEN
        ALLOCATE(y(SIZE(x)), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
    END IF
    y(:) = x(:)
  END SUBROUTINE assign_if_present_real_allocatable_1d

  SUBROUTINE assign_if_present_character_allocatable(y,x)
    CHARACTER(LEN=:), ALLOCATABLE, INTENT(INOUT) :: y
    CHARACTER(LEN=*), OPTIONAL,    INTENT(IN)    :: x

    IF (.NOT. PRESENT(x)) RETURN
    IF ( TRIM(x) == '' )  RETURN
    y = x
  END SUBROUTINE assign_if_present_character_allocatable

  SUBROUTINE alloc_int_1d(array, allocSize)
    INTEGER, ALLOCATABLE, INTENT(INOUT) :: array(:)
    INTEGER, VALUE :: allocSize

    INTEGER :: error
    CHARACTER(*), PARAMETER :: routine = modname//":alloc_int_1d"

    IF(allocSize < 1) allocSize = 1
    IF(ALLOCATED(array)) THEN
        IF(SIZE(array) == allocSize) RETURN
        DEALLOCATE(array)
    END IF
    ALLOCATE(array(allocSize), STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
  END SUBROUTINE alloc_int_1d

  SUBROUTINE alloc_double_1d(array, allocSize)
    REAL(dp), ALLOCATABLE, INTENT(INOUT) :: array(:)
    INTEGER, VALUE :: allocSize

    INTEGER :: error
    CHARACTER(*), PARAMETER :: routine = modname//":alloc_double_1d"

    IF(allocSize < 1) allocSize = 1
    IF(ALLOCATED(array)) THEN
        IF(SIZE(array) == allocSize) RETURN
        DEALLOCATE(array)
    END IF
    ALLOCATE(array(allocSize), STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
  END SUBROUTINE alloc_double_1d

  SUBROUTINE alloc_single_1d(array, allocSize)
    REAL(sp), ALLOCATABLE, INTENT(INOUT) :: array(:)
    INTEGER, VALUE :: allocSize

    INTEGER :: error
    CHARACTER(*), PARAMETER :: routine = modname//":alloc_single_1d"

    IF(allocSize < 1) allocSize = 1
    IF(ALLOCATED(array)) THEN
        IF(SIZE(array) == allocSize) RETURN
        DEALLOCATE(array)
    END IF
    ALLOCATE(array(allocSize), STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
  END SUBROUTINE alloc_single_1d

  SUBROUTINE ensureSize_dp_1d(buffer, requiredSize, do_copy_in)
    REAL(dp), POINTER, INTENT(INOUT) :: buffer(:)
    INTEGER, VALUE ::requiredSize
    LOGICAL, OPTIONAL, INTENT(IN) :: do_copy_in
    REAL(dp), POINTER :: newBuffer(:)
    INTEGER :: oldSize, error
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":ensureSize_dp_1d"
    LOGICAL :: do_copy

    IF (PRESENT(do_copy_in)) THEN
      do_copy = do_copy_in
    ELSE
      do_copy = .true.
    END IF
    IF(ASSOCIATED(buffer)) THEN
      oldSize = SIZE(buffer, 1)
      IF(oldSize >= requiredSize) RETURN  ! nothing to DO IF it's already big enough
      requiredSize = MAX(requiredSize, INT(REAL(oldSize, dp) * 1.1_dp)) ! avoid quadratic complexity
      ALLOCATE(newBuffer(requiredSize), STAT = error)
      IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
      IF (do_copy) THEN
        newBuffer(1:oldSize) = buffer(1:oldSize)
        newBuffer(oldSize + 1:requiredSize) = 0._dp
      END IF
      DEALLOCATE(buffer)
      buffer => newBuffer
      newBuffer => NULL()
    ELSE
      ALLOCATE(buffer(requiredSize), STAT = error)
      IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
      IF (do_copy) THEN
        buffer(1:requiredSize) = 0._dp
      END IF
    END IF
  END SUBROUTINE ensureSize_dp_1d

  SUBROUTINE ensureSize_sp_1d(buffer, requiredSize, do_copy_in)
    REAL(sp), POINTER, INTENT(INOUT) :: buffer(:)
    INTEGER, VALUE ::requiredSize
    LOGICAL, OPTIONAL, INTENT(IN) :: do_copy_in
    REAL(sp), POINTER :: newBuffer(:)
    INTEGER :: oldSize, error
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":ensureSize_dp_1d"
    LOGICAL :: do_copy

    IF (PRESENT(do_copy_in)) THEN
      do_copy = do_copy_in
    ELSE
      do_copy = .true.
    END IF
    IF(ASSOCIATED(buffer)) THEN
      oldSize = SIZE(buffer, 1)
      IF(oldSize >= requiredSize) RETURN  ! nothing to DO IF it's already big enough
      requiredSize = MAX(requiredSize, INT(REAL(oldSize, dp) * 1.1_dp))
      ALLOCATE(newBuffer(requiredSize), STAT = error)
      IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
      IF (do_copy) THEN
        newBuffer(1:oldSize) = buffer(1:oldSize)
        newBuffer(oldSize + 1:requiredSize) = 0._sp
      END IF
      DEALLOCATE(buffer)
      buffer => newBuffer
      newBuffer => NULL()
    ELSE
      ALLOCATE(buffer(requiredSize), STAT = error)
      IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
      IF (do_copy) THEN
          buffer(1:requiredSize) = 0._sp
      END IF
    END IF
  END SUBROUTINE ensureSize_sp_1d

  SUBROUTINE ensureSize_int_1d(buffer, requiredSize, do_copy_in)
    INTEGER, POINTER, INTENT(INOUT) :: buffer(:)
    INTEGER, VALUE ::requiredSize
    LOGICAL, OPTIONAL, INTENT(IN) :: do_copy_in
    INTEGER, POINTER :: newBuffer(:)
    INTEGER :: oldSize, error
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":ensureSize_dp_1d"
    LOGICAL :: do_copy

    IF (PRESENT(do_copy_in)) THEN
      do_copy = do_copy_in
    ELSE
      do_copy = .true.
    END IF
    IF(ASSOCIATED(buffer)) THEN
      oldSize = SIZE(buffer, 1)
      IF(oldSize >= requiredSize) RETURN  ! nothing to DO IF it's already big enough
      requiredSize = MAX(requiredSize, INT(REAL(oldSize, dp) * 1.1_dp))
      ALLOCATE(newBuffer(requiredSize), STAT = error)
      IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
      IF (do_copy) THEN
        newBuffer(1:oldSize) = buffer(1:oldSize)
        newBuffer(oldSize + 1:requiredSize) = 0
      END IF
      DEALLOCATE(buffer)
      buffer => newBuffer
      newBuffer => NULL()
    ELSE
      ALLOCATE(buffer(requiredSize), STAT = error)
      IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
      IF (do_copy) THEN
        buffer(1:requiredSize) = 0
      END IF
    END IF
  END SUBROUTINE ensureSize_int_1d

  !>
  !! Swap content of two Integers
  !!
  !! Swap content of two Integers
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2014-10-28)
  !!
  SUBROUTINE swap_int(a,b)
    INTEGER, INTENT(INOUT) :: a
    INTEGER, INTENT(INOUT) :: b

    ! local variables
    INTEGER :: temp
  !-----------------------------
    temp = a
    a    = b
    b    = temp
  END SUBROUTINE swap_int

  !>
  !! Expand array by given size
  !!
  !! Expand a 1D character array by given size.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2014-10-28)
  !!
  SUBROUTINE resize_arr_c1d(arr,nelem)
    ! GCC 4.9.0 complained about CHARACTER(:); Cray did not!
    CHARACTER(len=vname_len), ALLOCATABLE, INTENT(INOUT) :: arr(:)   ! array to be resized
    INTEGER                  , INTENT(IN)    :: nelem    ! number of elements to expand
    !
    ! local variables
    CHARACTER(len=vname_len), ALLOCATABLE :: tmp_arr(:)
    INTEGER :: istat                   ! status
    !-----------------------------

    ! If arr has not yet been allocated, do it.
    IF (.NOT. ALLOCATED(arr)) THEN
      ALLOCATE(arr(1), STAT=istat)
      IF (istat /= 0) THEN
        CALL finish('mo_fortran_tools:resize_arr_c1d', &
             &      ' initial allocation of array arr failed')
      ENDIF
    ELSE
      ! check for appropriate nelem
      IF ( nelem < 0) CALL finish('mo_fortran_tools:resize_arr_c1d', &
        &                         ' nelem must be > 0')

      ! allocate temporary array of size SIZE(arr)+nelem
      ALLOCATE(tmp_arr(1:SIZE(arr)+nelem), STAT=istat)
      IF (istat /= 0) THEN
        CALL finish('mo_fortran_tools:resize_char_arr_1d', &
             &      'allocation of array tmp_arr failed')
      ENDIF

      ! copy
      tmp_arr(1:SIZE(arr)) = arr(1:SIZE(arr))

      CALL move_alloc(tmp_arr, arr)
      ! now arr has been resized to the size of tmp_arr,
      ! and tmp_arr is deallocated.
    ENDIF

  END SUBROUTINE resize_arr_c1d

  !> copy state, omp parallel, does not wait for other threads to complete
  SUBROUTINE copy_2d_dp(src, dest, opt_acc_async)
    REAL(dp), INTENT(in) :: src(:, :)
    REAL(dp), INTENT(out) :: dest(:, :)
    LOGICAL, INTENT(IN), OPTIONAL :: opt_acc_async
    INTEGER :: i1, i2, m1, m2
    m1 = SIZE(dest, 1)
    m2 = SIZE(dest, 2)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) COLLAPSE(2) IF(i_am_accel_node)
#ifdef __INTEL_COMPILER
!$omp do private(i1,i2)
#else
!$omp do collapse(2)
#endif
    DO i2 = 1, m2
      DO i1 = 1, m1
        dest(i1, i2) = src(i1, i2)
      END DO
    END DO
!$omp end do nowait
    CALL acc_wait_if_requested(1, opt_acc_async)

  END SUBROUTINE copy_2d_dp

  !> copy state, omp parallel, does not wait for other threads to complete
  SUBROUTINE copy_3d_dp(src, dest, opt_acc_async)
    REAL(dp), INTENT(in) :: src(:, :, :)
    REAL(dp), INTENT(out) :: dest(:, :, :)
    LOGICAL, INTENT(IN), OPTIONAL :: opt_acc_async
    INTEGER :: i1, i2, i3, m1, m2, m3
    m1 = SIZE(dest, 1)
    m2 = SIZE(dest, 2)
    m3 = SIZE(dest, 3)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) COLLAPSE(3) IF(i_am_accel_node)
#if (defined(_CRAYFTN) || defined(__INTEL_COMPILER))
!$omp do private(i1,i2,i3)
#else
!$omp do collapse(3)
#endif
    DO i3 = 1, m3
      DO i2 = 1, m2
        DO i1 = 1, m1
          dest(i1, i2, i3) = src(i1, i2, i3)
        END DO
      END DO
    END DO
!$omp end do nowait

    CALL acc_wait_if_requested(1, opt_acc_async)
  END SUBROUTINE copy_3d_dp

  !> copy state, omp parallel, does not wait for other threads to complete
  SUBROUTINE copy_4d_dp(src, dest, opt_acc_async)
    REAL(dp), INTENT(in) :: src(:, :, :, :)
    REAL(dp), INTENT(out) :: dest(:, :, :, :)
    LOGICAL, INTENT(IN), OPTIONAL :: opt_acc_async
    INTEGER :: i1, i2, i3, i4, m1, m2, m3, m4
    m1 = SIZE(dest, 1)
    m2 = SIZE(dest, 2)
    m3 = SIZE(dest, 3)
    m4 = SIZE(dest, 4)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) COLLAPSE(4) IF(i_am_accel_node)
#if (defined(_CRAYFTN) || defined(__INTEL_COMPILER))
!$omp do private(i1,i2,i3,i4)
#else
!$omp do collapse(4)
#endif
    DO i4 = 1, m4
      DO i3 = 1, m3
        DO i2 = 1, m2
          DO i1 = 1, m1
            dest(i1, i2, i3, i4) = src(i1, i2, i3, i4)
          END DO
        END DO
      END DO
    END DO
!$omp end do nowait

    CALL acc_wait_if_requested(1, opt_acc_async)
  END SUBROUTINE copy_4d_dp

  !> copy state, omp parallel, does not wait for other threads to complete
  SUBROUTINE copy_5d_dp(src, dest, opt_acc_async)
    REAL(dp), INTENT(in) :: src(:, :, :, :, :)
    REAL(dp), INTENT(out) :: dest(:, :, :, :, :)
    LOGICAL, INTENT(IN), OPTIONAL :: opt_acc_async
    INTEGER :: i1, i2, i3, i4, i5, m1, m2, m3, m4, m5
    m1 = SIZE(dest, 1)
    m2 = SIZE(dest, 2)
    m3 = SIZE(dest, 3)
    m4 = SIZE(dest, 4)
    m5 = SIZE(dest, 5)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) COLLAPSE(5) IF(i_am_accel_node)
#if (defined(__INTEL_COMPILER))
!$omp do private(i1,i2,i3,i4,i5)
#else
!$omp do collapse(5)
#endif
    DO i5 = 1, m5
      DO i4 = 1, m4
        DO i3 = 1, m3
          DO i2 = 1, m2
            DO i1 = 1, m1
              dest(i1, i2, i3, i4, i5) = src(i1, i2, i3, i4, i5)
            END DO
          END DO
        END DO
      END DO
    END DO
!$omp end do nowait

    CALL acc_wait_if_requested(1, opt_acc_async)
  END SUBROUTINE copy_5d_dp

  !> copy state, omp parallel, does not wait for other threads to complete
  SUBROUTINE copy_5d_sp(src, dest, opt_acc_async)
    REAL(sp), INTENT(in) :: src(:, :, :, :, :)
    REAL(sp), INTENT(out) :: dest(:, :, :, :, :)
    LOGICAL, INTENT(IN), OPTIONAL :: opt_acc_async
    INTEGER :: i1, i2, i3, i4, i5, m1, m2, m3, m4, m5
    m1 = SIZE(dest, 1)
    m2 = SIZE(dest, 2)
    m3 = SIZE(dest, 3)
    m4 = SIZE(dest, 4)
    m5 = SIZE(dest, 5)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) COLLAPSE(5) IF(i_am_accel_node)
#if (defined(__INTEL_COMPILER))
!$omp do private(i1,i2,i3,i4,i5)
#else
!$omp do collapse(5)
#endif
    DO i5 = 1, m5
      DO i4 = 1, m4
        DO i3 = 1, m3
          DO i2 = 1, m2
            DO i1 = 1, m1
              dest(i1, i2, i3, i4, i5) = src(i1, i2, i3, i4, i5)
            END DO
          END DO
        END DO
      END DO
    END DO
!$omp end do nowait

    CALL acc_wait_if_requested(1, opt_acc_async)
  END SUBROUTINE copy_5d_sp

  !> copy state, omp parallel, does not wait for other threads to complete
  SUBROUTINE copy_2d_spdp(src, dest, opt_acc_async)
    REAL(sp), INTENT(in) :: src(:, :)
    REAL(dp), INTENT(out) :: dest(:, :)
    LOGICAL, INTENT(IN), OPTIONAL :: opt_acc_async
    INTEGER :: i1, i2, m1, m2
    m1 = SIZE(dest, 1)
    m2 = SIZE(dest, 2)
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) COLLAPSE(2) IF(i_am_accel_node)
#if (defined(__INTEL_COMPILER))
!$omp do private(i1,i2)
#else
!$omp do collapse(2)
#endif
    DO i2 = 1, m2
      DO i1 = 1, m1
        dest(i1, i2) = REAL(src(i1, i2),KIND=dp)
      END DO
    END DO
!$omp end do nowait
    CALL acc_wait_if_requested(1, opt_acc_async)

  END SUBROUTINE copy_2d_spdp

  !> copy state, omp parallel, does not wait for other threads to complete
  SUBROUTINE copy_3d_spdp(src, dest, opt_acc_async)
    REAL(sp), INTENT(in) :: src(:, :, :)
    REAL(dp), INTENT(out) :: dest(:, :, :)
    LOGICAL, INTENT(IN), OPTIONAL :: opt_acc_async
    INTEGER :: i1, i2, i3, m1, m2, m3
    m1 = SIZE(dest, 1)
    m2 = SIZE(dest, 2)
    m3 = SIZE(dest, 3)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) COLLAPSE(3) IF(i_am_accel_node)
#if (defined(__INTEL_COMPILER))
!$omp do private(i1,i2,i3)
#else
!$omp do collapse(3)
#endif
    DO i3 = 1, m3
      DO i2 = 1, m2
        DO i1 = 1, m1
          dest(i1, i2, i3) = REAL(src(i1, i2, i3),KIND=dp)
        END DO
      END DO
    END DO
!$omp end do nowait
    CALL acc_wait_if_requested(1, opt_acc_async)

  END SUBROUTINE copy_3d_spdp

  !> copy state, omp parallel, does not wait for other threads to complete
  SUBROUTINE copy_4d_spdp(src, dest, opt_acc_async)
    REAL(sp), INTENT(in) :: src(:, :, :, :)
    REAL(dp), INTENT(out) :: dest(:, :, :, :)
    LOGICAL, INTENT(IN), OPTIONAL :: opt_acc_async
    INTEGER :: i1, i2, i3, i4, m1, m2, m3, m4
    m1 = SIZE(dest, 1)
    m2 = SIZE(dest, 2)
    m3 = SIZE(dest, 3)
    m4 = SIZE(dest, 4)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) COLLAPSE(4) IF(i_am_accel_node)
#if (defined(__INTEL_COMPILER))
!$omp do private(i1,i2,i3,i4)
#else
!$omp do collapse(4)
#endif
    DO i4 = 1, m4
      DO i3 = 1, m3
        DO i2 = 1, m2
          DO i1 = 1, m1
            dest(i1, i2, i3, i4) = REAL(src(i1, i2, i3, i4),KIND=dp)
          END DO
        END DO
      END DO
    END DO
!$omp end do nowait
    CALL acc_wait_if_requested(1, opt_acc_async)

  END SUBROUTINE copy_4d_spdp

  !> copy state, omp parallel, does not wait for other threads to complete
  SUBROUTINE copy_5d_spdp(src, dest, opt_acc_async)
    REAL(sp), INTENT(in) :: src(:, :, :, :, :)
    REAL(dp), INTENT(out) :: dest(:, :, :, :, :)
    LOGICAL, INTENT(IN), OPTIONAL :: opt_acc_async
    INTEGER :: i1, i2, i3, i4, i5, m1, m2, m3, m4, m5
    m1 = SIZE(dest, 1)
    m2 = SIZE(dest, 2)
    m3 = SIZE(dest, 3)
    m4 = SIZE(dest, 4)
    m5 = SIZE(dest, 5)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) COLLAPSE(5) IF(i_am_accel_node)
#if (defined(__INTEL_COMPILER))
!$omp do private(i1,i2,i3,i4,i5)
#else
!$omp do collapse(5)
#endif
    DO i5 = 1, m5
      DO i4 = 1, m4
        DO i3 = 1, m3
          DO i2 = 1, m2
            DO i1 = 1, m1
              dest(i1, i2, i3, i4, i5) = REAL(src(i1, i2, i3, i4, i5),KIND=dp)
            END DO
          END DO
        END DO
      END DO
    END DO
!$omp end do nowait

    CALL acc_wait_if_requested(1, opt_acc_async)
  END SUBROUTINE copy_5d_spdp

  !> copy state, omp parallel, does not wait for other threads to complete
  SUBROUTINE copy_2d_i4(src, dest, opt_acc_async)
    INTEGER(ik4), INTENT(in) :: src(:, :)
    INTEGER(ik4), INTENT(out) :: dest(:, :)
    LOGICAL, INTENT(IN), OPTIONAL :: opt_acc_async
    INTEGER :: i1, i2, m1, m2
    m1 = SIZE(dest, 1)
    m2 = SIZE(dest, 2)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) COLLAPSE(2) IF(i_am_accel_node)
#if (defined(__INTEL_COMPILER))
!$omp do private(i1,i2)
#else
!$omp do collapse(2)
#endif
    DO i2 = 1, m2
      DO i1 = 1, m1
        dest(i1, i2) = src(i1, i2)
      END DO
    END DO
!$omp end do nowait

    CALL acc_wait_if_requested(1, opt_acc_async)
  END SUBROUTINE copy_2d_i4

  !> copy state, omp parallel, does not wait for other threads to complete
  SUBROUTINE copy_5d_i4(src, dest, opt_acc_async)
    INTEGER(ik4), INTENT(in) :: src(:, :, :, :, :)
    INTEGER(ik4), INTENT(out) :: dest(:, :, :, :, :)
    LOGICAL, INTENT(IN), OPTIONAL :: opt_acc_async
    INTEGER :: i1, i2, i3, i4, i5, m1, m2, m3, m4, m5
    m1 = SIZE(dest, 1)
    m2 = SIZE(dest, 2)
    m3 = SIZE(dest, 3)
    m4 = SIZE(dest, 4)
    m5 = SIZE(dest, 5)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) COLLAPSE(5) IF(i_am_accel_node)
#if (defined(__INTEL_COMPILER))
!$omp do private(i1,i2,i3,i4,i5)
#else
!$omp do collapse(5)
#endif
    DO i5 = 1, m5
      DO i4 = 1, m4
        DO i3 = 1, m3
          DO i2 = 1, m2
            DO i1 = 1, m1
              dest(i1, i2, i3, i4, i5) = src(i1, i2, i3, i4, i5)
            END DO
          END DO
        END DO
      END DO
    END DO
!$omp end do nowait

    CALL acc_wait_if_requested(1, opt_acc_async)
  END SUBROUTINE copy_5d_i4

  !> copy state, omp parallel, does not wait for other threads to complete
  SUBROUTINE copy_5d_l(src, dest, opt_acc_async)
    LOGICAL, INTENT(in) :: src(:, :, :, :, :)
    LOGICAL, INTENT(out) :: dest(:, :, :, :, :)
    LOGICAL, INTENT(IN), OPTIONAL :: opt_acc_async
    INTEGER :: i1, i2, i3, i4, i5, m1, m2, m3, m4, m5
    m1 = SIZE(dest, 1)
    m2 = SIZE(dest, 2)
    m3 = SIZE(dest, 3)
    m4 = SIZE(dest, 4)
    m5 = SIZE(dest, 5)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) COLLAPSE(5) IF(i_am_accel_node)
#if (defined(__INTEL_COMPILER))
    !$omp do private(i1,i2,i3,i4,i5)
#else
    !$omp do collapse(5)
#endif
    DO i5 = 1, m5
      DO i4 = 1, m4
        DO i3 = 1, m3
          DO i2 = 1, m2
            DO i1 = 1, m1
              dest(i1, i2, i3, i4, i5) = src(i1, i2, i3, i4, i5)
            END DO
          END DO
        END DO
      END DO
    END DO
    !$omp end do nowait

    CALL acc_wait_if_requested(1, opt_acc_async)
  END SUBROUTINE copy_5d_l

  !> copy state, omp parallel, does not wait for other threads to complete
  SUBROUTINE copy_3d_i4(src, dest, opt_acc_async)
    INTEGER(ik4), INTENT(in) :: src(:, :, :)
    INTEGER(ik4), INTENT(out) :: dest(:, :, :)
    LOGICAL, INTENT(IN), OPTIONAL :: opt_acc_async
    INTEGER :: i1, i2, i3, m1, m2, m3
    m1 = SIZE(dest, 1)
    m2 = SIZE(dest, 2)
    m3 = SIZE(dest, 3)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) COLLAPSE(3) IF(i_am_accel_node)
#if (defined(__INTEL_COMPILER))
!$omp do private(i1,i2,i3)
#else
!$omp do collapse(3)
#endif
    DO i3 = 1, m3
      DO i2 = 1, m2
        DO i1 = 1, m1
          dest(i1, i2, i3) = src(i1, i2, i3)
        END DO
      END DO
    END DO
!$omp end do nowait

    CALL acc_wait_if_requested(1, opt_acc_async)
  END SUBROUTINE copy_3d_i4

  SUBROUTINE init_zero_3d_dp(init_var, opt_acc_async)
    REAL(dp), INTENT(out) :: init_var(:, :, :)
    LOGICAL, INTENT(IN), OPTIONAL :: opt_acc_async
    INTEGER :: i1, i2, i3, m1, m2, m3

    m1 = SIZE(init_var, 1)
    m2 = SIZE(init_var, 2)
    m3 = SIZE(init_var, 3)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) COLLAPSE(3) IF(i_am_accel_node)
#if (defined(__INTEL_COMPILER) || defined(_CRAYFTN))
!$omp do private(i1,i2,i3)
#else
!$omp do collapse(3)
#endif
    DO i3 = 1, m3
      DO i2 = 1, m2
        DO i1 = 1, m1
          init_var(i1, i2, i3) = 0.0_dp
        END DO
      END DO
    END DO
!$omp end do nowait

    CALL acc_wait_if_requested(1, opt_acc_async)
  END SUBROUTINE init_zero_3d_dp

  SUBROUTINE init_zero_4d_dp(init_var, opt_acc_async)
    REAL(dp), INTENT(out) :: init_var(:, :, :, :)
    LOGICAL, INTENT(IN), OPTIONAL :: opt_acc_async
    INTEGER :: i1, i2, i3, i4, m1, m2, m3, m4

    m1 = SIZE(init_var, 1)
    m2 = SIZE(init_var, 2)
    m3 = SIZE(init_var, 3)
    m4 = SIZE(init_var, 4)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) COLLAPSE(4) IF(i_am_accel_node)
#if (defined(__INTEL_COMPILER) || defined(_CRAYFTN))
!$omp do private(i1,i2,i3,i4)
#else
!$omp do collapse(4)
#endif
    DO i4 = 1, m4
      DO i3 = 1, m3
        DO i2 = 1, m2
          DO i1 = 1, m1
            init_var(i1, i2, i3, i4) = 0.0_dp
          END DO
        END DO
      END DO
    END DO
!$omp end do nowait

    CALL acc_wait_if_requested(1, opt_acc_async)
  END SUBROUTINE init_zero_4d_dp

  SUBROUTINE init_zero_4d_sp(init_var, opt_acc_async)
    REAL(sp), INTENT(out) :: init_var(:, :, :, :)
    LOGICAL, INTENT(IN), OPTIONAL :: opt_acc_async
    INTEGER :: i1, i2, i3, i4, m1, m2, m3, m4

    m1 = SIZE(init_var, 1)
    m2 = SIZE(init_var, 2)
    m3 = SIZE(init_var, 3)
    m4 = SIZE(init_var, 4)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) COLLAPSE(4) IF(i_am_accel_node)
#if (defined(__INTEL_COMPILER) || defined(_CRAYFTN))
!$omp do private(i1,i2,i3,i4)
#else
!$omp do collapse(4)
#endif
    DO i4 = 1, m4
      DO i3 = 1, m3
        DO i2 = 1, m2
          DO i1 = 1, m1
            init_var(i1, i2, i3, i4) = 0.0_sp
          END DO
        END DO
      END DO
    END DO
!$omp end do nowait

    CALL acc_wait_if_requested(1, opt_acc_async)
  END SUBROUTINE init_zero_4d_sp

  SUBROUTINE init_zero_3d_sp(init_var, opt_acc_async)
    REAL(sp), INTENT(out) :: init_var(:, :, :)
    LOGICAL, INTENT(IN), OPTIONAL :: opt_acc_async
    INTEGER :: i1, i2, i3, m1, m2, m3

    m1 = SIZE(init_var, 1)
    m2 = SIZE(init_var, 2)
    m3 = SIZE(init_var, 3)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) COLLAPSE(3) IF(i_am_accel_node)
#if (defined(__INTEL_COMPILER))
!$omp do private(i1,i2,i3)
#else
!$omp do collapse(3)
#endif
    DO i3 = 1, m3
      DO i2 = 1, m2
        DO i1 = 1, m1
          init_var(i1, i2, i3) = 0.0_sp
        END DO
      END DO
    END DO
!$omp end do nowait

    CALL acc_wait_if_requested(1, opt_acc_async)

  END SUBROUTINE init_zero_3d_sp

  SUBROUTINE init_zero_2d_i4(init_var, opt_acc_async)
    INTEGER(ik4), INTENT(out) :: init_var(:, :)
    LOGICAL, INTENT(IN), OPTIONAL :: opt_acc_async
    INTEGER :: i1, i2, m1, m2

    m1 = SIZE(init_var, 1)
    m2 = SIZE(init_var, 2)


    !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) COLLAPSE(2) IF(i_am_accel_node)
#if (defined(__INTEL_COMPILER))
!$omp do private(i1,i2)
#else
!$omp do collapse(2)
#endif
    DO i2 = 1, m2
      DO i1 = 1, m1
        init_var(i1, i2) = 0_ik4
      END DO
    END DO
!$omp end do nowait

    CALL acc_wait_if_requested(1, opt_acc_async)
  END SUBROUTINE init_zero_2d_i4

  SUBROUTINE init_zero_3d_i4(init_var, opt_acc_async)
    INTEGER(ik4), INTENT(out) :: init_var(:, :, :)
    LOGICAL, INTENT(IN), OPTIONAL :: opt_acc_async
    INTEGER :: i1, i2, i3, m1, m2, m3

    m1 = SIZE(init_var, 1)
    m2 = SIZE(init_var, 2)
    m3 = SIZE(init_var, 3)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) COLLAPSE(3) IF(i_am_accel_node)
#if (defined(__INTEL_COMPILER))
!$omp do private(i1,i2,i3)
#else
!$omp do collapse(3)
#endif
    DO i3 = 1, m3
      DO i2 = 1, m2
        DO i1 = 1, m1
          init_var(i1, i2, i3) = 0_ik4
        END DO
      END DO
    END DO
!$omp end do nowait

    CALL acc_wait_if_requested(1, opt_acc_async)
  END SUBROUTINE init_zero_3d_i4

  SUBROUTINE init_zero_4d_i4(init_var, opt_acc_async)
    INTEGER(ik4), INTENT(out) :: init_var(:, :, :, :)
    LOGICAL, INTENT(IN), OPTIONAL :: opt_acc_async
    INTEGER :: i1, i2, i3, i4, m1, m2, m3, m4

    m1 = SIZE(init_var, 1)
    m2 = SIZE(init_var, 2)
    m3 = SIZE(init_var, 3)
    m4 = SIZE(init_var, 4)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) COLLAPSE(4) IF(i_am_accel_node)
#if (defined(__INTEL_COMPILER))
!$omp do private(i1,i2,i3,i4)
#else
!$omp do collapse(4)
#endif
    DO i4 = 1, m4
      DO i3 = 1, m3
        DO i2 = 1, m2
          DO i1 = 1, m1
            init_var(i1, i2, i3, i4) = 0_ik4
          END DO
        END DO
      END DO
    END DO
!$omp end do nowait

    CALL acc_wait_if_requested(1, opt_acc_async)
  END SUBROUTINE init_zero_4d_i4

  SUBROUTINE init_zero_1d_dp(init_var, opt_acc_async)
    REAL(dp), INTENT(out) :: init_var(:)
    LOGICAL, INTENT(IN), OPTIONAL :: opt_acc_async
    INTEGER :: i1, m1

    m1 = SIZE(init_var, 1)
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
!$omp do
    DO i1 = 1, m1
      init_var(i1) = 0.0_dp
    END DO
!$omp end do nowait

    CALL acc_wait_if_requested(1, opt_acc_async)
  END SUBROUTINE init_zero_1d_dp

  SUBROUTINE init_zero_1d_sp(init_var, opt_acc_async)
    REAL(sp), INTENT(out) :: init_var(:)
    LOGICAL, INTENT(IN), OPTIONAL :: opt_acc_async
    INTEGER :: i1, m1

    m1 = SIZE(init_var, 1)
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
    !$omp do
    DO i1 = 1, m1
      init_var(i1) = 0.0_dp
    END DO
!$omp end do nowait

    CALL acc_wait_if_requested(1, opt_acc_async)
  END SUBROUTINE init_zero_1d_sp

  SUBROUTINE init_zero_2d_dp(init_var, opt_acc_async)
    REAL(dp), INTENT(out) :: init_var(:, :)
    LOGICAL, INTENT(IN), OPTIONAL :: opt_acc_async
    INTEGER :: i1, i2, m1, m2

    m1 = SIZE(init_var, 1)
    m2 = SIZE(init_var, 2)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) COLLAPSE(2) IF(i_am_accel_node)
#if (defined(__INTEL_COMPILER))
!$omp do private(i1,i2)
#else
!$omp do collapse(2)
#endif
    DO i2 = 1, m2
      DO i1 = 1, m1
        init_var(i1, i2) = 0.0_dp
      END DO
    END DO
!$omp end do nowait

    CALL acc_wait_if_requested(1, opt_acc_async)
  END SUBROUTINE init_zero_2d_dp

  SUBROUTINE init_2d_dp(init_var, init_val, opt_acc_async)
    REAL(dp), INTENT(out) :: init_var(:, :)
    REAL(dp), INTENT(in) :: init_val
    LOGICAL, INTENT(IN), OPTIONAL :: opt_acc_async

    INTEGER :: i1, i2, m1, m2

    m1 = SIZE(init_var, 1)
    m2 = SIZE(init_var, 2)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) COLLAPSE(2) IF(i_am_accel_node)
#if (defined(__INTEL_COMPILER))
!$omp do private(i1,i2)
#else
!$omp do collapse(2)
#endif
      DO i2 = 1, m2
        DO i1 = 1, m1
          init_var(i1, i2) = init_val
        END DO
      END DO
!$omp end do nowait

    CALL acc_wait_if_requested(1, opt_acc_async)
  END SUBROUTINE init_2d_dp

  SUBROUTINE init_3d_dp(init_var, init_val, opt_acc_async)
    REAL(dp), INTENT(out) :: init_var(:, :, :)
    REAL(dp), INTENT(in) :: init_val
    LOGICAL, INTENT(IN), OPTIONAL :: opt_acc_async

    INTEGER :: i1, i2, i3, m1, m2, m3

    m1 = SIZE(init_var, 1)
    m2 = SIZE(init_var, 2)
    m3 = SIZE(init_var, 3)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) COLLAPSE(3) IF(i_am_accel_node)
#if (defined(__INTEL_COMPILER))
!$omp do private(i1,i2,i3)
#else
!$omp do collapse(3)
#endif
    DO i3 = 1, m3
      DO i2 = 1, m2
        DO i1 = 1, m1
          init_var(i1, i2, i3) = init_val
        END DO
      END DO
    END DO
!$omp end do nowait

    CALL acc_wait_if_requested(1, opt_acc_async)
  END SUBROUTINE init_3d_dp

  SUBROUTINE init_3d_spdp(init_var, init_val, opt_acc_async)
    REAL(sp), INTENT(out) :: init_var(:, :, :)
    REAL(dp), INTENT(in) :: init_val
    LOGICAL, INTENT(IN), OPTIONAL :: opt_acc_async

    INTEGER :: i1, i2, i3, m1, m2, m3

    m1 = SIZE(init_var, 1)
    m2 = SIZE(init_var, 2)
    m3 = SIZE(init_var, 3)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) COLLAPSE(3) IF(i_am_accel_node)
#if (defined(__INTEL_COMPILER))
!$omp do private(i1,i2,i3)
#else
!$omp do collapse(3)
#endif
    DO i3 = 1, m3
      DO i2 = 1, m2
        DO i1 = 1, m1
          init_var(i1, i2, i3) = REAL(init_val,KIND=sp)
        END DO
      END DO
    END DO
!$omp end do nowait

    CALL acc_wait_if_requested(1, opt_acc_async)
  END SUBROUTINE init_3d_spdp

  SUBROUTINE init_5d_dp(init_var, init_val, opt_acc_async)
    REAL(dp), INTENT(out) :: init_var(:, :, :, :, :)
    REAL(dp), INTENT(in) :: init_val
    LOGICAL, INTENT(IN), OPTIONAL :: opt_acc_async

    INTEGER :: i1, i2, i3, i4, i5, m1, m2, m3, m4, m5

    m1 = SIZE(init_var, 1)
    m2 = SIZE(init_var, 2)
    m3 = SIZE(init_var, 3)
    m4 = SIZE(init_var, 4)
    m5 = SIZE(init_var, 5)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) COLLAPSE(5) IF(i_am_accel_node)
#if (defined(__INTEL_COMPILER))
!$omp do private(i1,i2,i3,i4,i5)
#else
!$omp do collapse(5)
#endif
    DO i5 = 1, m5
      DO i4 = 1, m4
        DO i3 = 1, m3
          DO i2 = 1, m2
            DO i1 = 1, m1
              init_var(i1, i2, i3, i4, i5) = init_val
            END DO
          END DO
        END DO
      END DO
    END DO
!$omp end do nowait

    CALL acc_wait_if_requested(1, opt_acc_async)
  END SUBROUTINE init_5d_dp


  SUBROUTINE init_5d_sp(init_var, init_val, opt_acc_async)
    REAL(sp), INTENT(out) :: init_var(:, :, :, :, :)
    REAL(sp), INTENT(in) :: init_val
    LOGICAL, INTENT(IN), OPTIONAL :: opt_acc_async

    INTEGER :: i1, i2, i3, i4, i5, m1, m2, m3, m4, m5

    m1 = SIZE(init_var, 1)
    m2 = SIZE(init_var, 2)
    m3 = SIZE(init_var, 3)
    m4 = SIZE(init_var, 4)
    m5 = SIZE(init_var, 5)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) COLLAPSE(5) IF(i_am_accel_node)
#if (defined(__INTEL_COMPILER))
!$omp do private(i1,i2,i3,i4,i5)
#else
!$omp do collapse(5)
#endif
    DO i5 = 1, m5
      DO i4 = 1, m4
        DO i3 = 1, m3
          DO i2 = 1, m2
            DO i1 = 1, m1
              init_var(i1, i2, i3, i4, i5) = init_val
            END DO
          END DO
        END DO
      END DO
    END DO
!$omp end do nowait

    CALL acc_wait_if_requested(1, opt_acc_async)
  END SUBROUTINE init_5d_sp



  SUBROUTINE init_5d_i4(init_var, init_val, opt_acc_async)
    INTEGER(ik4), INTENT(out) :: init_var(:, :, :, :, :)
    INTEGER(ik4), INTENT(in) :: init_val
    LOGICAL, INTENT(IN), OPTIONAL :: opt_acc_async

    INTEGER :: i1, i2, i3, i4, i5, m1, m2, m3, m4, m5

    m1 = SIZE(init_var, 1)
    m2 = SIZE(init_var, 2)
    m3 = SIZE(init_var, 3)
    m4 = SIZE(init_var, 4)
    m5 = SIZE(init_var, 5)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) COLLAPSE(5) IF(i_am_accel_node)
#if (defined(__INTEL_COMPILER))
!$omp do private(i1,i2,i3,i4,i5)
#else
!$omp do collapse(5)
#endif
    DO i5 = 1, m5
      DO i4 = 1, m4
        DO i3 = 1, m3
          DO i2 = 1, m2
            DO i1 = 1, m1
              init_var(i1, i2, i3, i4, i5) = init_val
            END DO
          END DO
        END DO
      END DO
    END DO
!$omp end do nowait

    CALL acc_wait_if_requested(1, opt_acc_async)
  END SUBROUTINE init_5d_i4


  SUBROUTINE init_5d_l(init_var, init_val, opt_acc_async)
    LOGICAL, INTENT(out) :: init_var(:, :, :, :, :)
    LOGICAL, INTENT(in)  :: init_val
    LOGICAL, INTENT(IN), OPTIONAL :: opt_acc_async

    INTEGER :: i1, i2, i3, i4, i5, m1, m2, m3, m4, m5

    m1 = SIZE(init_var, 1)
    m2 = SIZE(init_var, 2)
    m3 = SIZE(init_var, 3)
    m4 = SIZE(init_var, 4)
    m5 = SIZE(init_var, 5)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) COLLAPSE(5) IF(i_am_accel_node)
#if (defined(__INTEL_COMPILER))
!$omp do private(i1,i2,i3,i4,i5)
#else
!$omp do collapse(5)
#endif
    DO i5 = 1, m5
      DO i4 = 1, m4
        DO i3 = 1, m3
          DO i2 = 1, m2
            DO i1 = 1, m1
              init_var(i1, i2, i3, i4, i5) = init_val
            END DO
          END DO
        END DO
      END DO
    END DO
!$omp end do nowait

    CALL acc_wait_if_requested(1, opt_acc_async)
  END SUBROUTINE init_5d_l



  SUBROUTINE var_scale_3d_dp(var, scale_val, opt_acc_async)
    REAL(dp), INTENT(inout) :: var(:, :, :)
    REAL(dp), INTENT(in) :: scale_val
    LOGICAL, INTENT(IN), OPTIONAL :: opt_acc_async

    INTEGER :: i1, i2, i3, m1, m2, m3

    m1 = SIZE(var, 1)
    m2 = SIZE(var, 2)
    m3 = SIZE(var, 3)


    !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) COLLAPSE(3) IF(i_am_accel_node)
#if (defined(__INTEL_COMPILER))
!$omp do private(i1,i2,i3)
#else
!$omp do collapse(3)
#endif
    DO i3 = 1, m3
      DO i2 = 1, m2
        DO i1 = 1, m1
          var(i1, i2, i3) = var(i1, i2, i3) * scale_val
        END DO
      END DO
    END DO
!$omp end do nowait

    CALL acc_wait_if_requested(1, opt_acc_async)
  END SUBROUTINE var_scale_3d_dp


  ! add a constant value to a 3D field
  SUBROUTINE var_addc_3d_dp(var, add_val, opt_acc_async)
    REAL(dp), INTENT(inout) :: var(:, :, :)
    REAL(dp), INTENT(in) :: add_val
    LOGICAL, INTENT(IN), OPTIONAL :: opt_acc_async

    INTEGER :: i1, i2, i3, m1, m2, m3

    m1 = SIZE(var, 1)
    m2 = SIZE(var, 2)
    m3 = SIZE(var, 3)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) COLLAPSE(3) IF(i_am_accel_node)
#if (defined(__INTEL_COMPILER))
!$omp do private(i1,i2,i3)
#else
!$omp do collapse(3)
#endif
    DO i3 = 1, m3
      DO i2 = 1, m2
        DO i1 = 1, m1
          var(i1, i2, i3) = var(i1, i2, i3) + add_val
        END DO
      END DO
    END DO
!$omp end do nowait

    CALL acc_wait_if_requested(1, opt_acc_async)
  END SUBROUTINE var_addc_3d_dp


  SUBROUTINE negative2zero_4d_dp(var, opt_acc_async)
    REAL(dp), INTENT(inout) :: var(:, :, :, :)
    LOGICAL, INTENT(IN), OPTIONAL :: opt_acc_async

    INTEGER :: i1, i2, i3, i4, m1, m2, m3, m4
    REAL(dp) :: v

    m1 = SIZE(var, 1)
    m2 = SIZE(var, 2)
    m3 = SIZE(var, 3)
    m4 = SIZE(var, 4)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) PRIVATE(v) ASYNC(1) COLLAPSE(4) IF(i_am_accel_node)
#if (defined(__INTEL_COMPILER))
!$omp do private(i1,i2,i3,i4)
#else
!$omp do collapse(4)
#endif
    DO i4 = 1, m4
      DO i3 = 1, m3
        DO i2 = 1, m2
          DO i1 = 1, m1
            v = var(i1, i2, i3, i4)
            var(i1, i2, i3, i4) = (ABS(v) + v) * 0.5_dp
          END DO
        END DO
      END DO
    END DO
!$omp end do nowait

    CALL acc_wait_if_requested(1, opt_acc_async)
  END SUBROUTINE negative2zero_4d_dp

  SUBROUTINE init_contiguous_dp(var, n, v, opt_acc_async, lacc)
    INTEGER, INTENT(in) :: n
    REAL(dp), INTENT(out) :: var(n)
    REAL(dp), INTENT(in) :: v
    LOGICAL, INTENT(in), OPTIONAL :: opt_acc_async
    LOGICAL, INTENT(in), OPTIONAL :: lacc

    INTEGER :: i
    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
!$omp do
    DO i = 1, n
      var(i) = v
    END DO
!$omp end do nowait

    CALL acc_wait_if_requested(1, opt_acc_async)
  END SUBROUTINE init_contiguous_dp

  SUBROUTINE init_zero_contiguous_dp(var, n, opt_acc_async, lacc)
    INTEGER, INTENT(in) :: n
    REAL(dp), INTENT(out) :: var(n)
    LOGICAL, INTENT(IN), OPTIONAL :: opt_acc_async
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    CALL init_contiguous_dp(var, n, 0.0_dp, opt_acc_async, lacc)
  END SUBROUTINE init_zero_contiguous_dp

  SUBROUTINE init_contiguous_sp(var, n, v, opt_acc_async, lacc)
    INTEGER, INTENT(in) :: n
    REAL(sp), INTENT(out) :: var(n)
    REAL(sp), INTENT(in) :: v
    LOGICAL, INTENT(in), OPTIONAL :: opt_acc_async
    LOGICAL, INTENT(in), OPTIONAL :: lacc

    INTEGER :: i
    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
!$omp do
    DO i = 1, n
      var(i) = v
    END DO
!$omp end do nowait
    CALL acc_wait_if_requested(1, opt_acc_async)

  END SUBROUTINE init_contiguous_sp

  SUBROUTINE init_zero_contiguous_sp(var, n, opt_acc_async, lacc)
    INTEGER, INTENT(in) :: n
    REAL(sp), INTENT(out) :: var(n)
    LOGICAL, INTENT(IN), OPTIONAL :: opt_acc_async
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    CALL init_contiguous_sp(var, n, 0.0_sp, opt_acc_async, lacc=lacc)
  END SUBROUTINE init_zero_contiguous_sp

  SUBROUTINE init_contiguous_i4(var, n, v, opt_acc_async, lacc)
    INTEGER, INTENT(in) :: n
    INTEGER(ik4), INTENT(out) :: var(n)
    INTEGER(ik4), INTENT(in) :: v
    LOGICAL, INTENT(in), OPTIONAL :: opt_acc_async
    LOGICAL, INTENT(in), OPTIONAL :: lacc

    INTEGER :: i
    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
!$omp do
    DO i = 1, n
      var(i) = v
    END DO
!$omp end do nowait

    CALL acc_wait_if_requested(1, opt_acc_async)
  END SUBROUTINE init_contiguous_i4

  SUBROUTINE init_contiguous_l(var, n, v, opt_acc_async, lacc)
    INTEGER, INTENT(in) :: n
    LOGICAL, INTENT(out) :: var(n)
    LOGICAL, INTENT(in) :: v
    LOGICAL, INTENT(IN), OPTIONAL :: opt_acc_async
    LOGICAL, INTENT(in), OPTIONAL :: lacc

    INTEGER :: i
    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
!$omp do
    DO i = 1, n
      var(i) = v
    END DO
!$omp end do nowait

    CALL acc_wait_if_requested(1, opt_acc_async)
  END SUBROUTINE init_contiguous_l

  FUNCTION minval_1d(var, lacc)
  !! Computes the MINVAL(var)
  !! This wrapper enables the use of OpenACC without using ACC-KERNELS
    INTEGER, INTENT(IN) :: var(:) ! input array
    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! if true, use OpenACC
    LOGICAL :: lzacc ! non-optional version of lacc
    INTEGER :: minval_1d, i, s1

#ifdef _OPENACC
    CALL set_acc_host_or_device(lzacc, lacc)

    s1 = SIZE(var, 1)

    minval_1d = HUGE(minval_1d)

    !$ACC PARALLEL DEFAULT(PRESENT) COPY(minval_1d) ASYNC(1) REDUCTION(MIN: minval_1d) IF(lacc)
    !$ACC LOOP GANG VECTOR
    DO i = 1, s1
      minval_1d = MIN(minval_1d, var(i)) ! The loop is equivalent to MINVAL(var(:))
    ENDDO
    !$ACC END PARALLEL
    !$ACC WAIT ! required to sync result back to CPU
#else
    minval_1d = MINVAL(var(:))
#endif

  END FUNCTION minval_1d

  SUBROUTINE insert_dimension_r_dp_3_2_s(ptr_out, ptr_in, in_shape, &
       new_dim_rank)
    INTEGER, PARAMETER :: out_rank = 3
    INTEGER, INTENT(in) :: in_shape(out_rank-1), new_dim_rank
    REAL(dp), POINTER, INTENT(out) :: ptr_out(:,:,:)
    REAL(dp), TARGET, INTENT(in) :: ptr_in
    INTEGER :: out_shape(out_rank), i
    TYPE(c_ptr) :: cptr
    out_shape(1:out_rank-1) = in_shape
    cptr = C_LOC(ptr_in)
    DO i = out_rank, new_dim_rank+1, -1
      out_shape(i) = out_shape(i-1)
    END DO
    out_shape(new_dim_rank) = 1
    CALL C_F_POINTER(cptr, ptr_out, out_shape)
  END SUBROUTINE insert_dimension_r_dp_3_2_s

  SUBROUTINE insert_dimension_r_dp_3_2(ptr_out, ptr_in, new_dim_rank)
    INTEGER, PARAMETER :: out_rank = 3
    REAL(dp), POINTER, INTENT(out) :: ptr_out(:,:,:)
    ! note: must have target attribute in caller!
    REAL(dp), TARGET, INTENT(in) :: ptr_in(:,:)
    INTEGER, INTENT(in) :: new_dim_rank
    INTEGER :: base_shape(out_rank-1), &
         in_shape(out_rank-1), in_stride(out_rank-1), &
         out_shape(out_rank), out_stride(out_rank), i
    INTEGER, PARAMETER :: elem_byte_size=8
    IF (SIZE(ptr_in) > 0) THEN
      ! reconstruct underlying array shape and corresponding stride
      in_shape = SHAPE(ptr_in)
      in_stride(1) = 1
      in_stride(2) = in_shape(1)
      IF (in_shape(1) > 1 .AND. in_shape(2) > 1) THEN
        CALL util_stride_2d(in_stride, elem_byte_size, &
             C_LOC(ptr_in(1, 1)), C_LOC(ptr_in(2, 1)), &
             C_LOC(ptr_in(1, 2)))
        base_shape(1) = in_stride(2)
      ELSE IF (in_shape(1) > 1) THEN
        CALL util_stride_1d(in_stride(1), elem_byte_size, &
             C_LOC(ptr_in(1, 1)), C_LOC(ptr_in(2, 1)))
        base_shape(1) = in_stride(1) * in_shape(1)
      ELSE IF (in_shape(2) > 1) THEN
        CALL util_stride_1d(in_stride(2), elem_byte_size, &
             C_LOC(ptr_in(1, 1)), C_LOC(ptr_in(1, 2)))
        base_shape(1) = in_stride(2)
      END IF
      base_shape(2) = in_shape(2)
      CALL insert_dimension_r_dp_3_2_s(ptr_out, ptr_in(1,1), &
           base_shape, new_dim_rank)
      IF (in_stride(1) > 1 .OR. in_stride(2) > in_shape(1) &
           .OR. base_shape(1) /= in_shape(1)) THEN
        out_stride(1) = in_stride(1)
        out_stride(2) = 1
        out_shape(1:out_rank-1) = in_shape
        DO i = out_rank, new_dim_rank+1, -1
          out_shape(i) = out_shape(i - 1)
          out_stride(i) = out_stride(i - 1)
        END DO
        out_stride(new_dim_rank) = 1
        out_shape(new_dim_rank) = 1
        out_shape = (out_shape - 1) * out_stride + 1
        ptr_out => ptr_out(:out_shape(1):out_stride(1), &
             &             :out_shape(2):out_stride(2), &
             &             :out_shape(3):out_stride(3))
      END IF
    ELSE
      out_shape(1:out_rank-1) = SHAPE(ptr_in)
      DO i = out_rank, new_dim_rank+1, -1
        out_shape(i) = out_shape(i - 1)
      END DO
      out_shape(new_dim_rank) = 1
      CALL C_F_POINTER(c_null_ptr, ptr_out, out_shape)
    END IF
  END SUBROUTINE insert_dimension_r_dp_3_2

  SUBROUTINE insert_dimension_r_sp_3_2_s(ptr_out, ptr_in, in_shape, &
       new_dim_rank)
    INTEGER, PARAMETER :: out_rank = 3
    INTEGER, INTENT(in) :: in_shape(out_rank-1), new_dim_rank
    REAL(sp), POINTER, INTENT(out) :: ptr_out(:,:,:)
    REAL(sp), TARGET, INTENT(in) :: ptr_in
    INTEGER :: out_shape(out_rank), i
    TYPE(c_ptr) :: cptr
    out_shape(1:out_rank-1) = in_shape
    cptr = C_LOC(ptr_in)
    DO i = out_rank, new_dim_rank+1, -1
      out_shape(i) = out_shape(i-1)
    END DO
    out_shape(new_dim_rank) = 1
    CALL C_F_POINTER(cptr, ptr_out, out_shape)
  END SUBROUTINE insert_dimension_r_sp_3_2_s

  SUBROUTINE insert_dimension_r_sp_3_2(ptr_out, ptr_in, new_dim_rank)
    INTEGER, PARAMETER :: out_rank = 3
    REAL(sp), POINTER, INTENT(out) :: ptr_out(:,:,:)
    ! note: must have target attribute in caller!
    REAL(sp), TARGET, INTENT(in) :: ptr_in(:,:)
    INTEGER, INTENT(in) :: new_dim_rank
    INTEGER :: base_shape(out_rank-1), &
         in_shape(out_rank-1), in_stride(out_rank-1), &
         out_shape(out_rank), out_stride(out_rank), i
    INTEGER, PARAMETER :: elem_byte_size=8
    IF (SIZE(ptr_in) > 0) THEN
      ! reconstruct underlying array shape and corresponding stride
      in_shape = SHAPE(ptr_in)
      in_stride(1) = 1
      in_stride(2) = in_shape(1)
      IF (in_shape(1) > 1 .AND. in_shape(2) > 1) THEN
        CALL util_stride_2d(in_stride, elem_byte_size, &
             C_LOC(ptr_in(1, 1)), C_LOC(ptr_in(2, 1)), &
             C_LOC(ptr_in(1, 2)))
        base_shape(1) = in_stride(2)
      ELSE IF (in_shape(1) > 1) THEN
        CALL util_stride_1d(in_stride(1), elem_byte_size, &
             C_LOC(ptr_in(1, 1)), C_LOC(ptr_in(2, 1)))
        base_shape(1) = in_stride(1) * in_shape(1)
      ELSE IF (in_shape(2) > 1) THEN
        CALL util_stride_1d(in_stride(2), elem_byte_size, &
             C_LOC(ptr_in(1, 1)), C_LOC(ptr_in(1, 2)))
        base_shape(1) = in_stride(2)
      END IF
      base_shape(2) = in_shape(2)
      CALL insert_dimension_r_sp_3_2_s(ptr_out, ptr_in(1,1), &
           base_shape, new_dim_rank)
      IF (in_stride(1) > 1 .OR. in_stride(2) > in_shape(1) &
           .OR. base_shape(1) /= in_shape(1)) THEN
        out_stride(1) = in_stride(1)
        out_stride(2) = 1
        out_shape(1:out_rank-1) = in_shape
        DO i = out_rank, new_dim_rank+1, -1
          out_shape(i) = out_shape(i - 1)
          out_stride(i) = out_stride(i - 1)
        END DO
        out_stride(new_dim_rank) = 1
        out_shape(new_dim_rank) = 1
        out_shape = (out_shape - 1) * out_stride + 1
        ptr_out => ptr_out(:out_shape(1):out_stride(1), &
             &             :out_shape(2):out_stride(2), &
             &             :out_shape(3):out_stride(3))
      END IF
    ELSE
      out_shape(1:out_rank-1) = SHAPE(ptr_in)
      DO i = out_rank, new_dim_rank+1, -1
        out_shape(i) = out_shape(i - 1)
      END DO
      out_shape(new_dim_rank) = 1
      CALL C_F_POINTER(c_null_ptr, ptr_out, out_shape)
    END IF
  END SUBROUTINE insert_dimension_r_sp_3_2

  SUBROUTINE insert_dimension_i4_3_2_s(ptr_out, ptr_in, in_shape, new_dim_rank)
    INTEGER, PARAMETER :: out_rank = 3
    INTEGER, INTENT(in) :: in_shape(out_rank-1), new_dim_rank
    INTEGER(ik4), POINTER, INTENT(out) :: ptr_out(:,:,:)
    INTEGER(ik4), TARGET, INTENT(in) :: ptr_in
    INTEGER :: out_shape(out_rank), i
    TYPE(c_ptr) :: cptr
    out_shape(1:out_rank-1) = in_shape
    cptr = C_LOC(ptr_in)
    DO i = out_rank, new_dim_rank+1, -1
      out_shape(i) = out_shape(i-1)
    END DO
    out_shape(new_dim_rank) = 1
    CALL C_F_POINTER(cptr, ptr_out, out_shape)
  END SUBROUTINE insert_dimension_i4_3_2_s

  ! insert dimension of size 1 (so that total array size remains the
  ! same but an extra dimension is inserted into the shape)
  SUBROUTINE insert_dimension_i4_3_2(ptr_out, ptr_in, new_dim_rank)
    INTEGER, PARAMETER :: out_rank = 3
    INTEGER(ik4), POINTER, INTENT(out) :: ptr_out(:,:,:)
    INTEGER(ik4), TARGET, INTENT(in) :: ptr_in(:,:)
    INTEGER, INTENT(in) :: new_dim_rank
    INTEGER :: base_shape(out_rank-1), &
         in_shape(out_rank-1), in_stride(out_rank-1), &
         out_shape(out_rank), out_stride(out_rank), i
    INTEGER, PARAMETER :: elem_byte_size=4
    IF (SIZE(ptr_in) > 0) THEN
      ! reconstruct underlying array shape and corresponding stride
      in_shape = SHAPE(ptr_in)
      in_stride(1) = 1
      in_stride(2) = in_shape(1)
      IF (in_shape(1) > 1 .AND. in_shape(2) > 1) THEN
        CALL util_stride_2d(in_stride, elem_byte_size, &
             C_LOC(ptr_in(1, 1)), C_LOC(ptr_in(2, 1)), &
             C_LOC(ptr_in(1, 2)))
        base_shape(1) = in_stride(2)
      ELSE IF (in_shape(1) > 1) THEN
        CALL util_stride_1d(in_stride(1), elem_byte_size, &
             C_LOC(ptr_in(1, 1)), C_LOC(ptr_in(2, 1)))
        base_shape(1) = in_stride(1) * in_shape(1)
      ELSE IF (in_shape(2) > 1) THEN
        CALL util_stride_1d(in_stride(2), elem_byte_size, &
             C_LOC(ptr_in(1, 1)), C_LOC(ptr_in(1, 2)))
        base_shape(1) = in_stride(2)
      END IF
      base_shape(2) = in_shape(2)
      CALL insert_dimension_i4_3_2_s(ptr_out, ptr_in(1,1), &
           base_shape, new_dim_rank)
      IF (in_stride(1) > 1 .OR. in_stride(2) > in_shape(1) &
           .OR. base_shape(1) /= in_shape(1)) THEN
        out_stride(1) = in_stride(1)
        out_stride(2) = 1
        out_shape(1:out_rank-1) = in_shape
        DO i = out_rank, new_dim_rank+1, -1
          out_shape(i) = out_shape(i - 1)
          out_stride(i) = out_stride(i - 1)
        END DO
        out_stride(new_dim_rank) = 1
        out_shape(new_dim_rank) = 1
        out_shape = (out_shape - 1) * out_stride + 1
        ptr_out => ptr_out(:out_shape(1):out_stride(1), &
             &             :out_shape(2):out_stride(2), &
             &             :out_shape(3):out_stride(3))
      END IF
    ELSE
      out_shape(1:out_rank-1) = SHAPE(ptr_in)
      DO i = out_rank, new_dim_rank+1, -1
        out_shape(i) = out_shape(i - 1)
      END DO
      out_shape(new_dim_rank) = 1
      CALL C_F_POINTER(c_null_ptr, ptr_out, out_shape)
    END IF
  END SUBROUTINE insert_dimension_i4_3_2

  SUBROUTINE insert_dimension_l_3_2_s(ptr_out, ptr_in, in_shape, new_dim_rank)
    INTEGER, PARAMETER :: out_rank = 3
    INTEGER, INTENT(in) :: in_shape(out_rank-1), new_dim_rank
    LOGICAL, POINTER, INTENT(out) :: ptr_out(:,:,:)
    LOGICAL, TARGET, INTENT(in) :: ptr_in
    INTEGER :: out_shape(out_rank), i
    out_shape(1:out_rank-1) = in_shape
    DO i = out_rank, new_dim_rank+1, -1
      out_shape(i) = out_shape(i-1)
    END DO
    out_shape(new_dim_rank) = 1
    CALL C_F_POINTER(C_LOC(ptr_in), ptr_out, out_shape)
  END SUBROUTINE insert_dimension_l_3_2_s

  ! insert dimension of size 1 (so that total array size remains the
  ! same but an extra dimension is inserted into the shape)
  SUBROUTINE insert_dimension_l_3_2(ptr_out, ptr_in, new_dim_rank)
    INTEGER, PARAMETER :: out_rank = 3
    LOGICAL, POINTER, INTENT(out) :: ptr_out(:,:,:)
    LOGICAL, TARGET, INTENT(in) :: ptr_in(:,:)
    INTEGER, INTENT(in) :: new_dim_rank
    INTEGER :: base_shape(out_rank-1), &
         in_shape(out_rank-1), in_stride(out_rank-1), &
         out_shape(out_rank), out_stride(out_rank), i
    INTEGER, PARAMETER :: elem_byte_size=4
    IF (SIZE(ptr_in) > 0) THEN
      ! reconstruct underlying array shape and corresponding stride
      in_shape = SHAPE(ptr_in)
      in_stride(1) = 1
      in_stride(2) = in_shape(1)
      IF (in_shape(1) > 1 .AND. in_shape(2) > 1) THEN
        CALL util_stride_2d(in_stride, elem_byte_size, &
             C_LOC(ptr_in(1, 1)), C_LOC(ptr_in(2, 1)), &
             C_LOC(ptr_in(1, 2)))
        base_shape(1) = in_stride(2)
      ELSE IF (in_shape(1) > 1) THEN
        CALL util_stride_1d(in_stride(1), elem_byte_size, &
             C_LOC(ptr_in(1, 1)), C_LOC(ptr_in(2, 1)))
        base_shape(1) = in_stride(1) * in_shape(1)
      ELSE IF (in_shape(2) > 1) THEN
        CALL util_stride_1d(in_stride(2), elem_byte_size, &
             C_LOC(ptr_in(1, 1)), C_LOC(ptr_in(1, 2)))
        base_shape(1) = in_stride(2)
      END IF
      base_shape(2) = in_shape(2)
      CALL insert_dimension_l_3_2_s(ptr_out, ptr_in(1,1), &
           base_shape, new_dim_rank)
      IF (in_stride(1) > 1 .OR. in_stride(2) > in_shape(1) &
           .OR. base_shape(1) /= in_shape(1)) THEN
        out_stride(1) = in_stride(1)
        out_stride(2) = 1
        out_shape(1:out_rank-1) = in_shape
        DO i = out_rank, new_dim_rank+1, -1
          out_shape(i) = out_shape(i - 1)
          out_stride(i) = out_stride(i - 1)
        END DO
        out_stride(new_dim_rank) = 1
        out_shape(new_dim_rank) = 1
        out_shape = (out_shape - 1) * out_stride + 1
        ptr_out => ptr_out(:out_shape(1):out_stride(1), &
             &             :out_shape(2):out_stride(2), &
             &             :out_shape(3):out_stride(3))
      END IF
    ELSE
      out_shape(1:out_rank-1) = SHAPE(ptr_in)
      DO i = out_rank, new_dim_rank+1, -1
        out_shape(i) = out_shape(i - 1)
      END DO
      out_shape(new_dim_rank) = 1
      CALL C_F_POINTER(c_null_ptr, ptr_out, out_shape)
    END IF
  END SUBROUTINE insert_dimension_l_3_2

  SUBROUTINE insert_dimension_r_dp_6_5_s(ptr_out, ptr_in, in_shape, &
       new_dim_rank)
    INTEGER, INTENT(in) :: in_shape(5), new_dim_rank
    REAL(dp), POINTER, INTENT(out) :: ptr_out(:,:,:,:,:,:)
    REAL(dp), TARGET, INTENT(in) :: ptr_in(in_shape(1),in_shape(2),&
         in_shape(3),in_shape(4),in_shape(5))
    INTEGER :: out_shape(6), i
    TYPE(c_ptr) :: cptr
    out_shape(1:5) = SHAPE(ptr_in)
    cptr = C_LOC(ptr_in)
    DO i = 6, new_dim_rank+1, -1
      out_shape(i) = out_shape(i-1)
    END DO
    out_shape(new_dim_rank) = 1
    CALL C_F_POINTER(cptr, ptr_out, out_shape)
  END SUBROUTINE insert_dimension_r_dp_6_5_s

  ! insert dimension of size 1 (so that total array size remains the
  ! same but an extra dimension is inserted into the shape)
  SUBROUTINE insert_dimension_r_dp_6_5(ptr_out, ptr_in, new_dim_rank)
    REAL(dp), POINTER, INTENT(out) :: ptr_out(:,:,:,:,:,:)
    REAL(dp), TARGET, INTENT(in) :: ptr_in(:,:,:,:,:)
    INTEGER, INTENT(in) :: new_dim_rank
    INTEGER :: in_shape(5)
    in_shape = SHAPE(ptr_in)
    CALL insert_dimension(ptr_out, ptr_in, in_shape, new_dim_rank)
  END SUBROUTINE insert_dimension_r_dp_6_5

  SUBROUTINE insert_dimension_r_sp_6_5_s(ptr_out, ptr_in, in_shape, new_dim_rank)
    INTEGER, INTENT(in) :: in_shape(5), new_dim_rank
    REAL(sp), POINTER, INTENT(out) :: ptr_out(:,:,:,:,:,:)
    REAL(sp), TARGET, INTENT(in) :: ptr_in(in_shape(1),in_shape(2),&
         in_shape(3),in_shape(4),in_shape(5))
    INTEGER :: out_shape(6), i
    TYPE(c_ptr) :: cptr
    out_shape(1:5) = SHAPE(ptr_in)
    cptr = C_LOC(ptr_in)
    DO i = 6, new_dim_rank+1, -1
      out_shape(i) = out_shape(i-1)
    END DO
    out_shape(new_dim_rank) = 1
    CALL C_F_POINTER(cptr, ptr_out, out_shape)
  END SUBROUTINE insert_dimension_r_sp_6_5_s

  ! insert dimension of size 1 (so that total array size remains the
  ! same but an extra dimension is inserted into the shape)
  SUBROUTINE insert_dimension_r_sp_6_5(ptr_out, ptr_in, new_dim_rank)
    REAL(sp), POINTER, INTENT(out) :: ptr_out(:,:,:,:,:,:)
    REAL(sp), TARGET, INTENT(in) :: ptr_in(:,:,:,:,:)
    INTEGER, INTENT(in) :: new_dim_rank
    INTEGER :: in_shape(5)
    in_shape = SHAPE(ptr_in)
    CALL insert_dimension(ptr_out, ptr_in, in_shape, new_dim_rank)
  END SUBROUTINE insert_dimension_r_sp_6_5

  SUBROUTINE insert_dimension_i4_6_5_s(ptr_out, ptr_in, in_shape, new_dim_rank)
    INTEGER, INTENT(in) :: in_shape(5), new_dim_rank
    INTEGER(ik4), POINTER, INTENT(out) :: ptr_out(:,:,:,:,:,:)
    INTEGER(ik4), TARGET, INTENT(in) :: ptr_in(in_shape(1),in_shape(2),&
         in_shape(3),in_shape(4),in_shape(5))
    INTEGER :: out_shape(6), i
    TYPE(c_ptr) :: cptr
    out_shape(1:5) = SHAPE(ptr_in)
    cptr = C_LOC(ptr_in)
    DO i = 6, new_dim_rank+1, -1
      out_shape(i) = out_shape(i-1)
    END DO
    out_shape(new_dim_rank) = 1
    CALL C_F_POINTER(cptr, ptr_out, out_shape)
  END SUBROUTINE insert_dimension_i4_6_5_s

  ! insert dimension of size 1 (so that total array size remains the
  ! same but an extra dimension is inserted into the shape)
  SUBROUTINE insert_dimension_i4_6_5(ptr_out, ptr_in, new_dim_rank)
    INTEGER(ik4), POINTER, INTENT(out) :: ptr_out(:,:,:,:,:,:)
    INTEGER(ik4), TARGET, INTENT(in) :: ptr_in(:,:,:,:,:)
    INTEGER, INTENT(in) :: new_dim_rank
    INTEGER :: in_shape(5)
    in_shape = SHAPE(ptr_in)
    CALL insert_dimension(ptr_out, ptr_in, in_shape, new_dim_rank)
  END SUBROUTINE insert_dimension_i4_6_5


  ! AUXILIARY ROUTINES FOR DEALLOCATION
  SUBROUTINE DO_DEALLOCATE_r4D(object)
    REAL(wp), ALLOCATABLE, INTENT(INOUT) :: object(:,:,:,:)
    INTEGER :: ierrstat
    IF (ALLOCATED(object)) THEN
      DEALLOCATE(object, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish("DO_DEALLOCATE_r4D", "DEALLOCATE failed!")
    END IF
  END SUBROUTINE DO_DEALLOCATE_R4D

  SUBROUTINE DO_DEALLOCATE_r3D(object)
    REAL(wp), ALLOCATABLE, INTENT(INOUT) :: object(:,:,:)
    INTEGER :: ierrstat
    IF (ALLOCATED(object)) THEN
      DEALLOCATE(object, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish("DO_DEALLOCATE_r3D", "DEALLOCATE failed!")
    END IF
  END SUBROUTINE DO_DEALLOCATE_R3D

  SUBROUTINE DO_DEALLOCATE_r2D(object)
    REAL(wp), ALLOCATABLE, INTENT(INOUT) :: object(:,:)
    INTEGER :: ierrstat
    IF (ALLOCATED(object)) THEN
      DEALLOCATE(object, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish("DO_DEALLOCATE_r2D", "DEALLOCATE failed!")
    END IF
  END SUBROUTINE DO_DEALLOCATE_R2D

  SUBROUTINE DO_DEALLOCATE_r1D(object)
    REAL(wp), ALLOCATABLE, INTENT(INOUT) :: object(:)
    INTEGER :: ierrstat
    IF (ALLOCATED(object)) THEN
      DEALLOCATE(object, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish("DO_DEALLOCATE_r1D", "DEALLOCATE failed!")
    END IF
  END SUBROUTINE DO_DEALLOCATE_R1D

  SUBROUTINE DO_DEALLOCATE_i3D(object)
    INTEGER, ALLOCATABLE, INTENT(INOUT) :: object(:,:,:)
    INTEGER :: ierrstat
    IF (ALLOCATED(object)) THEN
      DEALLOCATE(object, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish("DO_DEALLOCATE_i3D", "DEALLOCATE failed!")
    END IF
  END SUBROUTINE DO_DEALLOCATE_i3D

  SUBROUTINE DO_DEALLOCATE_i2D(object)
    INTEGER, ALLOCATABLE, INTENT(INOUT) :: object(:,:)
    INTEGER :: ierrstat
    IF (ALLOCATED(object)) THEN
      DEALLOCATE(object, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish("DO_DEALLOCATE_i2D", "DEALLOCATE failed!")
    END IF
  END SUBROUTINE DO_DEALLOCATE_i2D

  SUBROUTINE DO_DEALLOCATE_i1D(object)
    INTEGER, ALLOCATABLE, INTENT(INOUT) :: object(:)
    INTEGER :: ierrstat
    IF (ALLOCATED(object)) THEN
      DEALLOCATE(object, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish("DO_DEALLOCATE_i1D", "DEALLOCATE failed!")
    END IF
  END SUBROUTINE DO_DEALLOCATE_i1D

  SUBROUTINE DO_PTR_DEALLOCATE_r3D(object)
    REAL(wp), POINTER, INTENT(INOUT) :: object(:,:,:)
    INTEGER :: ierrstat
    IF (ASSOCIATED(object)) THEN
      DEALLOCATE(object, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish("DO_PTR_DEALLOCATE_r3D", "DEALLOCATE failed!")
    END IF
    NULLIFY(object)
  END SUBROUTINE DO_PTR_DEALLOCATE_R3D

  SUBROUTINE DO_PTR_DEALLOCATE_r2D(object)
    REAL(wp), POINTER, INTENT(INOUT) :: object(:,:)
    INTEGER :: ierrstat
    IF (ASSOCIATED(object)) THEN
      DEALLOCATE(object, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish("DO_PTR_DEALLOCATE_r2D", "DEALLOCATE failed!")
    END IF
    NULLIFY(object)
  END SUBROUTINE DO_PTR_DEALLOCATE_R2D

  SUBROUTINE DO_PTR_DEALLOCATE_dp1D(object)
    REAL(KIND=dp), POINTER, INTENT(INOUT) :: object(:)
    INTEGER :: ierr

    IF (ASSOCIATED(object)) THEN
      DEALLOCATE(object, STAT=ierr)
      IF (ierr /= SUCCESS) CALL finish("DO_PTR_DEALLOCATE_dp1D", "DEALLOCATE failed!")
    END IF
    NULLIFY(object)
  END SUBROUTINE DO_PTR_DEALLOCATE_dp1D

  SUBROUTINE DO_PTR_DEALLOCATE_sp1D(object)
    REAL(KIND=sp), POINTER, INTENT(INOUT) :: object(:)
    INTEGER :: ierr

    IF (ASSOCIATED(object)) THEN
      DEALLOCATE(object, STAT=ierr)
      IF (ierr /= SUCCESS) CALL finish("DO_PTR_DEALLOCATE_sp1D", "DEALLOCATE failed!")
    END IF
    NULLIFY(object)
  END SUBROUTINE DO_PTR_DEALLOCATE_sp1D

  SUBROUTINE DO_PTR_DEALLOCATE_int1D(object)
    INTEGER, POINTER, INTENT(INOUT) :: object(:)
    INTEGER :: ierr

    IF (ASSOCIATED(object)) THEN
      DEALLOCATE(object, STAT=ierr)
      IF (ierr /= SUCCESS) CALL finish("DO_PTR_DEALLOCATE_dp1D", "DEALLOCATE failed!")
    END IF
    NULLIFY(object)
  END SUBROUTINE DO_PTR_DEALLOCATE_int1D

  SUBROUTINE assert_acc_host_only(routine_name, lacc)
    CHARACTER(len=*), INTENT(in) :: routine_name
    LOGICAL, INTENT(in), OPTIONAL :: lacc

#ifdef _OPENACC
    IF ( PRESENT(lacc) .AND. lacc ) THEN
      CALL finish (routine_name, ' not supported on ACC device.')
    ENDIF
#endif
  END SUBROUTINE assert_acc_host_only

  SUBROUTINE assert_acc_device_only(routine_name, lacc)
    CHARACTER(len=*), INTENT(in) :: routine_name
    LOGICAL, INTENT(in), OPTIONAL :: lacc

#ifdef _OPENACC
    IF ( (.NOT. PRESENT(lacc)) .OR. (.NOT. lacc) ) THEN
      CALL finish (routine_name, ' not supported on ACC host.')
    ENDIF
#endif
  END SUBROUTINE assert_acc_device_only

  SUBROUTINE assert_lacc_equals_i_am_accel_node(routine_name, lacc)
    CHARACTER(len=*), INTENT(in) :: routine_name
    LOGICAL, INTENT(in), OPTIONAL :: lacc

#ifdef _OPENACC
    IF ( PRESENT(lacc) .AND. (lacc .neqv. i_am_accel_node)) THEN
      CALL finish(routine_name, 'lacc /= i_am_accel_node')
    ENDIF
#endif

  END SUBROUTINE assert_lacc_equals_i_am_accel_node

  SUBROUTINE set_acc_host_or_device(lzacc, lacc)
    LOGICAL, INTENT(out) :: lzacc
    LOGICAL, INTENT(in), OPTIONAL :: lacc

    lzacc = .FALSE.
#ifdef _OPENACC
    IF ( PRESENT(lacc) ) THEN
      lzacc = lacc
    ENDIF
#endif
  END SUBROUTINE set_acc_host_or_device

END MODULE mo_fortran_tools
