!! Initializes and controls the time stepping in the wave model.
!!
!!
!! @author Mikhail Dobrynin, DWD, 24.06.19
!!
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.


!----------------------------
MODULE mo_wave_stepping

  USE mo_exception,                ONLY: message
  USE mo_run_config,               ONLY: output_mode
  USE mo_name_list_output,         ONLY: write_name_list_output
  USE mtime,                       ONLY: datetime

  IMPLICIT NONE

  PRIVATE

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_wave_stepping'

  PUBLIC :: perform_wave_stepping

CONTAINS

  !-------------------------------------------------------------------------
  !>
  !! Organizes wave time stepping
  !!
  !! @par Revision History
  !! Mikhail Dobrynin, DWD, 24.06.19
  !!
  SUBROUTINE perform_wave_stepping (mtime_current)

    TYPE(datetime), INTENT(IN) :: mtime_current

    CHARACTER(len=*), PARAMETER :: routine = modname//':perform_wave_stepping'


    IF (output_mode%l_nml) THEN
      CALL write_name_list_output(jstep=0)
    END IF

    CALL message(TRIM(routine),'finished')

  END SUBROUTINE perform_wave_stepping

END MODULE mo_wave_stepping
