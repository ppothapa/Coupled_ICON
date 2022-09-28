!----------------------------------
!>
!! Routines to access basic coupler functionality
!!
!! @par Revision History
!! First version by Moritz Hanke,  DKRZ, May 2022.
!!
!! @par
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
#include "icon_definitions.inc"
!----------------------------
MODULE mo_coupling

  USE mo_io_units, ONLY: nerr
#if !defined NOMPI && defined YAC_coupling
  USE mo_yac_finterface, ONLY: yac_finit_comm, yac_ffinalize, &
                               yac_finit_comm_group, yac_fget_groupcomm
  USE mpi
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_coupler
  PUBLIC :: finalize_coupler
  PUBLIC :: coupler_config_files_exist

  CHARACTER(LEN=*), PARAMETER :: xml_filename = "coupling.xml"
  CHARACTER(LEN=*), PARAMETER :: xsd_filename = "coupling.xsd"

  LOGICAL :: config_files_have_been_checked = .FALSE.
  LOGICAL :: config_files_exist = .FALSE.
  LOGICAL :: yac_is_initialised = .FALSE.

  CHARACTER(*), PARAMETER :: modname = "mo_coupling"

CONTAINS

  LOGICAL FUNCTION coupler_config_files_exist()

    LOGICAL :: xml_exists, xsd_exists

    IF (config_files_have_been_checked) THEN

      coupler_config_files_exist = config_files_exist

    ELSE

      INQUIRE(FILE=TRIM(ADJUSTL(xml_filename)), EXIST=xml_exists)
      INQUIRE(FILE=TRIM(ADJUSTL(xsd_filename)), EXIST=xsd_exists)

      config_files_have_been_checked = .TRUE.
      config_files_exist = xml_exists .AND. xsd_exists
      coupler_config_files_exist = config_files_exist

    END IF

  END FUNCTION

  SUBROUTINE init_coupler(world_communicator, global_name)

    INTEGER, INTENT(INOUT) :: world_communicator
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: global_name

#if !defined NOMPI && defined YAC_coupling

    CHARACTER(*), PARAMETER :: routine = modname//":init_coupler"

    INTEGER :: group_comm, group_rank
    INTEGER :: result, ierror

    IF (coupler_config_files_exist()) THEN

      yac_is_initialised = .TRUE.

      IF (PRESENT(global_name)) THEN

        CALL yac_finit_comm_group( &
          TRIM(xml_filename), TRIM(xsd_filename), world_communicator, global_name)

        CALL yac_fget_groupcomm(group_comm)

        CALL mpi_comm_compare(world_communicator, group_comm, result, ierror)

        IF ((result /= MPI_IDENT) .AND. (result /= MPI_CONGRUENT)) THEN

          world_communicator = group_comm

          CALL mpi_comm_rank(group_comm, group_rank, ierror)

          IF (group_rank == 0) THEN

            CALL print_info_stderr( &
              routine, &
              'splitting the set of available MPI ranks (MPI_COMM_WORLD communicator).')
            CALL print_info_stderr( &
              routine, &
              '"' // global_name // '" operates on a proper subset of MPI_COMM_WORLD')
          END IF

        ELSE

          CALL mpi_comm_free(group_comm, ierror)

        END IF

      ELSE

        CALL yac_finit_comm( &
          TRIM(xml_filename), TRIM(xsd_filename), world_communicator)
      END IF
    END IF
#endif

  END SUBROUTINE init_coupler

  SUBROUTINE finalize_coupler

#if !defined NOMPI && defined YAC_coupling
    IF (yac_is_initialised) CALL yac_ffinalize
#endif

  END SUBROUTINE finalize_coupler

  SUBROUTINE print_info_stderr (name, text)
    CHARACTER (len=*), INTENT(in) :: name
    CHARACTER (len=*), INTENT(in) :: text

    WRITE (nerr,'(4a)') " ", TRIM(name), ": ", TRIM(text)

  END SUBROUTINE print_info_stderr

END MODULE mo_coupling
