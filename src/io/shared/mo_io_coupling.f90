!>
!! @brief Routines for the initialisation of YAC on IO processes
!!
!! @author
!!  Rene Redler (MPI-M)
!!
!! The purpose of routines construct_io_coupler and destruct_io_coupler is
!! to initialise YAC on asynchronous IO processes. The YAC component definition
!! and search are collective operations in the MPI sense. Thus all MPI processes
!! must call yac_fdef_comp(s) and yac_fsearch.
!! Currently, the IO processes do not receive any fields via YAC for output.
!! Therefore, yac_fsearch is called with an empty list of fields.
!! 
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_io_coupling

  USE mo_coupling_config,           ONLY: is_coupled_run
  USE mo_impl_constants,            ONLY: MAX_CHAR_LENGTH
  USE mo_exception,                 ONLY: message, message_text
  USE mo_yac_finterface,            ONLY: yac_fdef_comp, yac_fsearch

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_io_coupler
  PUBLIC :: destruct_io_coupler
 
CONTAINS

  SUBROUTINE construct_io_coupler ( comp_name )

    CHARACTER(LEN=*), INTENT(IN) :: comp_name
    CHARACTER(*), PARAMETER :: routine = "mo_io_coupling:construct_io_coupler"

    INTEGER :: comp_ids(1)
    INTEGER :: error_status

    IF ( is_coupled_run() ) THEN

      WRITE(message_text,*) "YAC initialisation for asynchronous I/O processes of type ", &
        &                   TRIM(comp_name)
      CALL message(routine, message_text)

      ! Inform YAC about what we are
      CALL yac_fdef_comp ( TRIM(comp_name), comp_ids(1) )

      ! IO processes need to participate in the YAC search
      CALL yac_fsearch ( error_status )
    ENDIF

  END SUBROUTINE construct_io_coupler

  SUBROUTINE destruct_io_coupler ( comp_name )

    CHARACTER(LEN=*), INTENT(IN) :: comp_name
 
    CHARACTER(*), PARAMETER :: routine = "mo_io_coupling:construct_io_coupler"

    IF ( is_coupled_run() ) THEN

      WRITE(message_text,*) "YAC termination of I/O process of type ", &
        &                   TRIM(comp_name)
      CALL message(routine, message_text)

    ENDIF

  END SUBROUTINE destruct_io_coupler

END MODULE mo_io_coupling
