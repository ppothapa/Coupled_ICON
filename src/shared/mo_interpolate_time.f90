!----------------------------
#include "omp_definitions.inc"
!----------------------------
MODULE mo_interpolate_time

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message, message_text, finish
  USE mo_parallel_config,   ONLY: nproma
  USE mo_impl_constants,    ONLY: MAX_CHAR_LENGTH
  USE mtime,                ONLY: datetime, max_datetime_str_len,             &
       &                          julianday, getJulianDayFromDatetime,        &
       &                          datetimetostring, getDatetimeFromJulianDay, &
       &                          juliandelta, OPERATOR(<), OPERATOR(>),      &
       &                          OPERATOR(>=), OPERATOR(-), ASSIGNMENT(=),   &
       &                          no_of_ms_in_a_day

  USE mo_time_config,    ONLY: time_config
  USE mo_mpi,            ONLY: my_process_is_mpi_workroot, &
       &                       process_mpi_root_id,        &
       &                       p_comm_work, p_bcast, p_pe_work
  USE mo_reader_abstract,   ONLY: t_abstract_reader
#ifdef _OPENACC
  USE mo_mpi,            ONLY: i_am_accel_node
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_time_intp

  TYPE t_time_intp
    TYPE(julianday), ALLOCATABLE :: times(:)
    ! tidx and tidx+1 are held by dataold and datanew
    INTEGER                      :: tidx
    ! these two give access to the data
    REAL(wp),            POINTER :: dataold(:,:,:,:)
    REAL(wp),            POINTER :: datanew(:,:,:,:)
    ! these two hold the data and should not be accessed directly
    REAL(wp),        ALLOCATABLE :: dataa(:,:,:,:)
    REAL(wp),        ALLOCATABLE :: datab(:,:,:,:)

    ! Who am I? Passed to netcdf while reading.
    CHARACTER(len=MAX_CHAR_LENGTH) :: var_name

    ! Options are constant, linear, and weird linear. Default is linear.
    INTEGER :: interpolation_mode

    CLASS(t_abstract_reader), POINTER :: reader
  CONTAINS
    procedure :: init => time_intp_init
    procedure :: intp => time_intp_intp
  END TYPE t_time_intp

  INTEGER, PARAMETER :: intModeConstant    = 0
  INTEGER, PARAMETER :: intModeLinear      = 1
  INTEGER, PARAMETER :: intModeLinearWeird = 11

  CHARACTER(len=*), PARAMETER :: modname = 'mo_interpolate_time'

CONTAINS

  SUBROUTINE time_intp_init(this, reader, local_time, var_name, int_mode)
    CLASS(t_time_intp),       TARGET, INTENT(  out) :: this
    CLASS(t_abstract_reader), TARGET, INTENT(inout) :: reader
    TYPE(datetime),          POINTER, INTENT(in   ) :: local_time
    CHARACTER(*),                     INTENT(in   ) :: var_name
    INTEGER,                OPTIONAL, INTENT(in   ) :: int_mode

    INTEGER :: ntimes
    INTEGER :: i

    CHARACTER(len=max_datetime_str_len)      :: date_str1, date_str2

    TYPE(julianday) :: current_jd
    TYPE(datetime)  :: current_dt

    CHARACTER(len=*), PARAMETER :: routine = modname//"::time_intp_init"

    this%reader   => reader
    this%var_name =  var_name

    IF (PRESENT(int_mode)) THEN
      this%interpolation_mode = int_mode
    ELSE
      this%interpolation_mode = intModeLinear
    ENDIF

    CALL this%reader%get_times(this%times)

    CALL getJulianDayFromDatetime(local_time, current_jd)
    CALL datetimetostring(local_time, date_str2)

    ntimes = SIZE(this%times)
    this%tidx = -1
    DO i = 1,ntimes
      CALL getDatetimeFromJulianDay(this%times(i), current_dt)
      CALL datetimetostring(current_dt, date_str1)
      IF (current_jd >= this%times(i)) THEN
        this%tidx = i
      ENDIF
    ENDDO

    IF (this%tidx ==  -1) THEN
      CALL datetimetostring(local_time, date_str1)
      CALL finish(routine,"Time to interpolate for "//trim(date_str1)//" not covered by boundary condition data file")
    ENDIF

    CALL getDatetimeFromJulianDay(this%times(this%tidx), current_dt)
    CALL datetimetostring(current_dt, date_str1)

    check_for_time_coverage: BLOCK
      TYPE(datetime) :: first_time_in_file, last_time_in_file

      CALL getDatetimeFromJulianDay(this%times(1), first_time_in_file)
      CALL getDatetimeFromJulianDay(this%times(ntimes), last_time_in_file)

      IF (time_config%tc_startdate < first_time_in_file) THEN
        CALL datetimetostring(time_config%tc_startdate, date_str1)
        CALL datetimetostring(first_time_in_file, date_str2)
        CALL finish(routine, "Start of run ("//TRIM(date_str1)//") before start of data ("//TRIM(date_str2)//") to be read")
      ENDIF
      CALL getDatetimeFromJulianDay(this%times(ntimes), last_time_in_file)
      IF (time_config%tc_stopdate > last_time_in_file) THEN
        CALL datetimetostring(time_config%tc_stopdate, date_str1)
        CALL datetimetostring(last_time_in_file, date_str2)
        CALL finish(routine, "End of run ("//TRIM(date_str1)//") after end of data ("//TRIM(date_str2)//") to be read")
      ENDIF

    END BLOCK check_for_time_coverage

    log_output: BLOCK
      TYPE(datetime) :: load_time_in_file, load_next_time_in_file

      CALL getDatetimeFromJulianDay(this%times(this%tidx), load_time_in_file)
      CALL datetimetostring(load_time_in_file, date_str1)
      CALL getDatetimeFromJulianDay(this%times(this%tidx+1), load_next_time_in_file)
      CALL datetimetostring(load_next_time_in_file, date_str2)
      WRITE(message_text,'(a,i0,a,i0,a)') &
           &         " loading data for "//TRIM(date_str1)//" (", this%tidx, &
           &                     ") and "//TRIM(date_str2)//" (", this%tidx+1, ")"
      CALL message(TRIM(routine),message_text)
    END BLOCK log_output

    !$ACC ENTER DATA COPYIN(this)
    CALL reader%get_one_timelev(this%tidx,   this%var_name, this%dataa)
    this%dataold => this%dataa
    CALL reader%get_one_timelev(this%tidx+1, this%var_name, this%datab)
    this%datanew => this%datab

  END SUBROUTINE time_intp_init

  SUBROUTINE time_intp_intp(this, local_time, interpolated)
    CLASS(t_time_intp), TARGET, INTENT(inout) :: this
    TYPE(datetime),    POINTER, INTENT(in   ) :: local_time
    REAL(wp),      ALLOCATABLE, INTENT(inout) :: interpolated(:,:,:,:)

    TYPE(julianday)                     :: current_jd
    TYPE(juliandelta)                   :: delta_1, delta_2
    REAL(wp)                            :: ds1, ds2, weight
    INTEGER                             :: jc,jk,jb,jw
    INTEGER                             :: nlen, nblks, npromz, nlev
    CHARACTER(len=max_datetime_str_len) :: date_str

    CHARACTER(*), PARAMETER :: routine = modname//"::time_intp_intp"

    CALL getJulianDayFromDatetime(local_time, current_jd)

    IF (current_jd > this%times(this%tidx+1)) THEN
      this%tidx    = this%tidx + 1

      log_output: BLOCK
        TYPE(datetime) :: time_in_file

        CALL getDatetimeFromJulianDay(this%times(this%tidx+1), time_in_file)
        CALL datetimetostring(time_in_file, date_str)
        WRITE(message_text,'(a,i0,a)') " loading new data, with tidx :", this%tidx+1, " ("//TRIM(date_str)//")"
        CALL message(TRIM(routine),message_text)
      END BLOCK log_output

      ! FORTRAN!?!1! this should look like:
      ! DEALLOCATE(this%dataold)
      ! this%dataold => this%datanew
      ! CALL this%reader%get_one_timelev(this%tidx+1, this%var_name, this%datanew)
      ! but fails, since POINTERs are a strange hybrid of pointers and
      ! allocatables, so no chance to implicitly access the allocatable
      ! behind the pointer.
      IF(ASSOCIATED(this%dataold, this%dataa)) THEN
        this%dataold => this%datab
        CALL this%reader%get_one_timelev(this%tidx+1, this%var_name, this%dataa)
        this%datanew => this%dataa
      ELSE
        this%dataold => this%dataa
        CALL this%reader%get_one_timelev(this%tidx+1, this%var_name, this%datab)
        this%datanew => this%datab
      ENDIF
    ENDIF

    IF (this%interpolation_mode == intModeConstant) THEN
      ! A weight of 0 makes the interpolation return dataold.
      weight = 0.0_wp
    ELSE IF (this%interpolation_mode == intModeLinear) THEN
      delta_1 = current_jd - this%times(this%tidx)
      delta_2 = this%times(this%tidx+1) - this%times(this%tidx)
      ds1 = 1.0e-3_wp * (no_of_ms_in_a_day * delta_1%day + delta_1%ms)
      ds2 = 1.0e-3_wp * (no_of_ms_in_a_day * delta_2%day + delta_2%ms)
      weight = ds1/ds2
    ELSE IF (this%interpolation_mode == intModeLinearWeird) THEN
      ! This should
      ! a) be renamed. (But what is this?)
      ! b) done similar to calculate_time_interpolation_weights in
      !    shared/mo_bcs_time_interpolation.f90 . Since I do not know how this
      !    works and whether it can be generalized to data with non-monthly
      !    intervals, I leave as this for now.
      CALL finish(routine, "You are weird")
    ENDIF

    ! ATTENTION: This is a trivial Fortran 2008 feature. Does explode with Intel!
    !ALLOCATE(interpolated, MOLD=this%dataa)

    IF (ALLOCATED(interpolated)) THEN
      IF (.NOT. ALL( SHAPE(interpolated) .EQ. SHAPE(this%dataa) )) THEN
        !$ACC EXIT DATA DELETE(interpolated)
        DEALLOCATE(interpolated)
      END IF
    END IF
    IF (.NOT. ALLOCATED(interpolated)) THEN
      ALLOCATE(interpolated(size(this%dataa,1), size(this%dataa,2), size(this%dataa,3), size(this%dataa,4)))
      !$ACC ENTER DATA CREATE(interpolated)
    END IF

    ! DA: Need to list this%dataxxx in the PRESENT section for attach
    !$ACC DATA PRESENT(interpolated,this)
    !$ACC DATA PRESENT(this%dataold, this%datanew)
    !$ACC KERNELS DEFAULT(NONE) IF(i_am_accel_node)
    interpolated(:,:,:,:) = 0.d0
    !$ACC END KERNELS

    nblks  = this%reader%get_nblks()
    npromz = this%reader%get_npromz()
    nlev   = size(interpolated,2)

    ! we need this mess, since npromz == nproma is not garantueed
    DO jw = 1,size(interpolated,4)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,nlen,jc) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1,nblks
        nlen = MERGE(nproma, npromz, jb /= nblks)
        ! DA: Need to list this%dataxxx in the PRESENT section for attach
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) IF(i_am_accel_node) &
        !$ACC   DEFAULT(PRESENT)
        DO jk = 1,nlev
          DO jc = 1,nlen
            interpolated(jc,jk,jb,jw) = (1.0_wp-weight) * this%dataold(jc,jk,jb,jw) &
              &                                +weight  * this%datanew(jc,jk,jb,jw)
          ENDDO
        ENDDO
        !$ACC END PARALLEL LOOP
      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    ENDDO

    !$ACC END DATA
    !$ACC END DATA

!    if (my_process_is_mpi_workroot()) THEN
!      print *, "blubba weight", weight, interpolated(1,1,1,1), &
!      this%datanew(1,1,1,1), this%dataold(1,1,1,1)
!      print *, (1-weight) * this%dataold(1,1,1,1) + weight*this%datanew(1,1,1,1)
!    ENDIF
  END SUBROUTINE time_intp_intp

END MODULE mo_interpolate_time
