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
  USE mo_math_constants,      ONLY: pi2
  USE mo_physical_constants,  ONLY: grav


  IMPLICIT NONE

  PUBLIC

  REAL(wp), PARAMETER :: EPS1  = 0.00001_wp !! small number to make sure that a
                                            !! solution is obtained in iteration
                                            !! with tau>tauw.

  REAL(wp), PARAMETER :: EMIN = 1.0E-12_wp  !! replaces the intrinsic tiny
  INTEGER,  PARAMETER :: EX_TAIL = -5       !! tail frequency exponent.

  REAL(wp), PARAMETER :: CDIS = 1.33_wp !! dissipation constant
  REAL(wp), PARAMETER :: DELTA = 0.5_wp !! weight linear, quadratic part
  REAL(wp), PARAMETER :: CONSD = -CDIS * pi2**9_wp / grav**4_wp !! dissipation constant for deep water
  REAL(wp), PARAMETER :: CONSS = -CDIS * pi2 !! dissipation constant for shallow water


END MODULE mo_wave_constants
