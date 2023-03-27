!!  Module determines some constants to be used by wave model.
!!
!! @author Mikhail Dobrynin, DWD, 11.10.2019
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
MODULE mo_wave_constants

  USE mo_kind,                ONLY: wp

  IMPLICIT NONE

  PUBLIC

  REAL(wp), PARAMETER :: ALPHA = 0.0075_wp  !! MINIMUM CHARNOCK CONSTANT (ECMWF CY45R1).
                                            !! 0.0060, if LE 30 frequencies changed !@waves TODO
                                            !! to 0.0075 in subroutine INITMDL !@waves TODO

  REAL(wp), PARAMETER :: EPS1  = 0.00001_wp !! SMALL NUMBER TO MAKE SURE THAT A
                                            !! SOLUTION IS OBTAINED IN ITERATION
                                            !! WITH TAU>TAUW.

  REAL(wp), PARAMETER :: EMIN = 1.0E-12_wp  !! REPLACES THE INTRINSIC TINY
  INTEGER,  PARAMETER :: EX_TAIL = -5       !! TAIL FREQUENCY EXPONENT.



END MODULE mo_wave_constants
