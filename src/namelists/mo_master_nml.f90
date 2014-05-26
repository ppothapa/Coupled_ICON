!>
!! @brief Master namelist.
!!        
!! @par Revision History
!! Created by Rene Redler (2011-03-22)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_master_nml

  USE mo_impl_constants, ONLY: MAX_CHAR_LENGTH, SUCCESS
  USE mo_exception,      ONLY: warning, finish
  USE mo_io_units,       ONLY: filename_max, nnml
  USE mo_namelist,       ONLY: open_nml, position_nml, POSITIONED
  USE mo_util_string,    ONLY: t_keyword_list, associate_keyword, with_keywords
  USE mo_nml_annotate,   ONLY: temp_defaults, temp_settings
  USE mo_mpi,            ONLY: my_process_is_stdio

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_master_namelist, lrestart, no_of_models
  PUBLIC :: t_master_nml, master_nml_array
  PUBLIC :: model_base_dir


  ! Component models
  !--------------------------------------------------------------
  ! TYPE definitions
  !> Holds a list of integers
  TYPE t_master_nml
    CHARACTER(len=132) :: model_name

    CHARACTER(len=filename_max) :: model_namelist_filename
    CHARACTER(len=filename_max) :: model_restart_info_filename

    INTEGER :: model_type
    INTEGER :: model_min_rank
    INTEGER :: model_max_rank
    INTEGER :: model_inc_rank
  END TYPE t_master_nml

  !> base directory (for file names)
  CHARACTER(len=99)  :: model_base_dir

  INTEGER, PARAMETER :: max_no_of_models=10
  INTEGER :: no_of_models
  TYPE(t_master_nml) :: master_nml_array(max_no_of_models)
  
    !-------------------------------------------------------------------------
    ! Namelist variables
    !-------------------------------------------------------------------------
    LOGICAL :: lrestart, nml_debug_coupler

CONTAINS
  !>
  !! Initialization of variables that contain general information
  !! about the coupled model run. The configuration is read from
  !! namelist 'master_nml'.
  !!
  !! @par Revision History
  !!
  INTEGER FUNCTION read_master_namelist(namelist_filename)
    
    CHARACTER(LEN=*), INTENT(in) :: namelist_filename
    !
    ! Local variables
    !
    !-------------------------------------------------------------------------
    ! Namelist variables
    !-------------------------------------------------------------------------
    CHARACTER(len=132) :: model_name
    CHARACTER(len=filename_max) :: model_namelist_filename
    CHARACTER(len=filename_max) :: model_restart_info_filename
    INTEGER :: model_type
    INTEGER :: model_min_rank
    INTEGER :: model_max_rank
    INTEGER :: model_inc_rank


    NAMELIST /master_model_nml/    &
      model_name,                  &
      model_namelist_filename,     &
      model_restart_info_filename, &
      model_type,                  &
      model_min_rank,              &
      model_max_rank,              &
      model_inc_rank               
    NAMELIST /master_nml/          &
      lrestart,                    &
      model_base_dir

    INTEGER :: istat
    INTEGER :: iunit
    LOGICAL :: rewnd
    CHARACTER(len=*), PARAMETER :: routine = 'mo_master_nml:read_master_namelist'
    TYPE (t_keyword_list), POINTER :: keywords         => NULL()
    TYPE (t_keyword_list), POINTER :: keywords_restart => NULL()

    !------------------------------------------------------------------
    ! Read  master_nml (done so far by all MPI processes)
    !------------------------------------------------------------------
    lrestart      = .FALSE.
    model_base_dir=''
    OPEN( nnml, FILE=TRIM(namelist_filename), IOSTAT=istat, &
        & STATUS='old', ACTION='read', DELIM='apostrophe')
    IF (istat/=0) THEN
      CALL warning(namelist_filename,"not found")
      read_master_namelist=-1
      RETURN
    ENDIF
     CALL position_nml('master_nml',STATUS=istat)
     IF (istat==POSITIONED) THEN
       READ (nnml, master_nml)
     ENDIF        
    !------------------------------------------------------------------
    ! Read  master_model_nml (done so far by all MPI processes)
    !------------------------------------------------------------------
    no_of_models = 0
    rewnd = .true.
    DO
      CALL position_nml('master_model_nml', lrewind=rewnd, status=istat)
      IF ( istat /= POSITIONED ) EXIT
      IF (no_of_models >= max_no_of_models) THEN
        CALL finish(routine, 'no_of_models >= max_no_of_models')
      ENDIF
      rewnd=.false.
      
      ! default values
      model_name=''
      model_namelist_filename=''
      model_restart_info_filename=''
      model_type=-1
      model_min_rank=0
      model_max_rank=-1 
      model_inc_rank=1
      
      IF (my_process_is_stdio()) THEN
        iunit = temp_defaults()
        WRITE(iunit, master_model_nml)  ! write defaults to temporary text file
      END IF
      READ (nnml, master_model_nml)                                        ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, master_model_nml)  ! write settings to temporary text file
      END IF

      no_of_models=no_of_models+1
      master_nml_array(no_of_models)%model_name              = model_name

      CALL associate_keyword("<path>", TRIM(model_base_dir), keywords)
      master_nml_array(no_of_models)%model_namelist_filename = &
        &  TRIM(with_keywords(keywords, model_namelist_filename))

      CALL associate_keyword("<path>", TRIM(model_base_dir), keywords_restart)
      master_nml_array(no_of_models)%model_restart_info_filename=&
        & TRIM(with_keywords(keywords_restart, model_restart_info_filename))

      master_nml_array(no_of_models)%model_type              = model_type
      master_nml_array(no_of_models)%model_min_rank          = model_min_rank
      master_nml_array(no_of_models)%model_max_rank          = model_max_rank
      master_nml_array(no_of_models)%model_inc_rank          = model_inc_rank

    ENDDO
      
    CLOSE (nnml, IOSTAT=istat)
   !------------------------------------------------------------------
   read_master_namelist=SUCCESS

  END FUNCTION read_master_namelist

END MODULE mo_master_nml
