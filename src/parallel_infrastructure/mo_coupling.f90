! Routines to access basic coupler functionality
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
#include "icon_definitions.inc"
!----------------------------
MODULE mo_coupling

  USE mo_io_units, ONLY: nerr
#if !defined NOMPI && defined YAC_coupling
  USE mo_yac_finterface, ONLY: yac_finit_comm, yac_ffinalize, &
                               yac_fmpi_handshake,            &
                               YAC_MAX_CHARLEN
  USE mpi
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_coupler
  PUBLIC :: finalize_coupler
  PUBLIC :: coupler_config_files_exist

  CHARACTER(LEN=*), PARAMETER :: yaml_filename = "coupling.yaml"

  LOGICAL :: config_files_have_been_checked = .FALSE.
  LOGICAL :: config_files_exist = .FALSE.
  LOGICAL :: yac_is_initialised = .FALSE.

  CHARACTER(*), PARAMETER :: modname = "mo_coupling"

CONTAINS

  LOGICAL FUNCTION coupler_config_files_exist()

    LOGICAL :: yaml_exists

    IF (config_files_have_been_checked) THEN

      coupler_config_files_exist = config_files_exist

    ELSE

      INQUIRE(FILE=TRIM(ADJUSTL(yaml_filename)), EXIST=yaml_exists)

      config_files_have_been_checked = .TRUE.
      config_files_exist = yaml_exists
      coupler_config_files_exist = config_files_exist

    END IF

  END FUNCTION

  SUBROUTINE init_coupler(world_communicator, global_name)

    INTEGER, INTENT(INOUT) :: world_communicator
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: global_name

#if !defined NOMPI && defined YAC_coupling

    CHARACTER(*), PARAMETER :: routine = modname//":init_coupler"

    INTEGER :: ierror

    INTEGER :: global_rank
    INTEGER :: yac_comm
    INTEGER :: group_comms(2)
    CHARACTER(len=YAC_MAX_CHARLEN) :: group_names(2)

    IF (coupler_config_files_exist()) THEN

      yac_is_initialised = .TRUE.

      IF (PRESENT(global_name)) THEN

         group_names(1) = "yac"
         group_names(2) = TRIM(global_name)

         CALL yac_fmpi_handshake( MPI_COMM_WORLD, group_names, group_comms)
         yac_comm = group_comms(1)
         world_communicator = group_comms(2)
         CALL yac_finit_comm(yac_comm)

         CALL MPI_COMM_RANK ( world_communicator, global_rank, ierror )
         IF ( global_rank == 0 ) &
            CALL yac_fread_config_yaml( TRIM(yaml_filename) )

         CALL mpi_comm_free(yac_comm, ierror)

      ELSE

        CALL yac_finit_comm( world_communicator )
        CALL MPI_COMM_RANK ( world_communicator, global_rank, ierror )
        IF ( global_rank == 0 ) &
           CALL yac_fread_config_yaml( TRIM(yaml_filename) )
        
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
