!>
!! Main program for the ICON ocean waves model
!!
!! @author Mikhail Dobrynin, DWD, 26.04.2019
!!
!! @par Revision History
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_wave_model

  USE mo_exception,               ONLY: message

  PUBLIC :: wave_model

CONTAINS

  SUBROUTINE wave_model(wave_namelist_filename,shr_namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: wave_namelist_filename
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename

    CHARACTER(*), PARAMETER :: routine = "mo_wave_model:wave_model"

    CALL message(routine, 'ICON-waves work in progress... Stop')

  END SUBROUTINE wave_model

END MODULE mo_wave_model
