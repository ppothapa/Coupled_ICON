!>
!! Waves config.
!! The content is mostly adopted from the WAM 4.5.4.
!!
!! @author Mikhail Dobrynin, DWD
!!
!! @par Revision History
!! Initial revision by Mikhail Dobrynin, DWD (21.12.2022)

!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_wave_config

  USE mo_kind,                ONLY: wp

  IMPLICIT NONE

  PUBLIC
  !
  ! TO DO: create t_wave_config type which includes the following switches
  !

  INTEGER  :: ndirs    ! NUMBER OF DIRECTIONS.
  INTEGER  :: nfreqs   ! NUMBER OF FREQUENCIES.

  REAL(wp) :: fr1      ! FIRST FREQUENCY [HZ].
  REAL(wp) :: CO       ! FREQUENCY RATIO

  REAL(wp) :: ALPHA      ! PHILLIPS' PARAMETER  (NOT USED IF IOPTI = 1)
  REAL(wp) :: FM         ! PEAK FREQUENCY (HZ) AND/OR MAXIMUM FREQUENCY
  REAL(wp) :: GAMMA_wave ! OVERSHOOT FACTOR
  REAL(wp) :: SIGMA_A    ! LEFT PEAK WIDTH
  REAL(wp) :: SIGMA_B    ! RIGHT PEAK WIDTH
  REAL(wp) :: THETAQ     ! WAVE DIRECTION (DEG) (NOT USED IF IOPTI = 1)
  REAL(wp) :: FETCH      ! FETCH IN METRES (IF ZERO THEN 0.5 OF THE LATITUDE INCREMENT IS USED.).

  REAL(wp) :: ROAIR   ! AIR DENSITY
  REAL(wp) :: RNUAIR  ! KINEMATIC AIR VISCOSITY
  REAL(wp) :: RNUAIRM ! KINEMATIC AIR VISCOSITY FOR MOMENTUM TRANSFER
  REAL(wp) :: ROWATER ! WATER DENSITY
  REAL(wp) :: XEPS
  REAL(wp) :: XINVEPS

  REAL(wp) :: XKAPPA  ! VON KARMAN CONSTANT.
  REAL(wp) :: XNLEV   ! WINDSPEED REF. LEVEL.
  REAL(wp) :: BETAMAX ! PARAMETER FOR WIND INPUT (ECMWF CY45R1).
  REAL(wp) :: ZALP    ! SHIFTS GROWTH CURVE (ECMWF CY45R1).

  LOGICAL :: coldstart           ! if .TRUE. start from initialisation without restart file
  INTEGER :: iforc_waves         ! 1 - test case
                                         ! 2 - forcing from coupled atmosphete
                                         ! 3 - forcing from data reader

  LOGICAL :: linput_sf1      ! if .TRUE., calculate wind input source function term, first call
  LOGICAL :: linput_sf2      ! if .TRUE., calculate wind input source function term, second call
  LOGICAL :: ldissip_sf      ! if .TRUE., calculate dissipation source function term
  LOGICAL :: lnon_linear_sf  ! if .TRUE., calculate non linear source function term
  LOGICAL :: lbottom_fric_sf ! if .TRUE., calculate bottom_friction source function term
  LOGICAL :: lwave_stress1   ! if .TRUE., calculate wave stress, first call
  LOGICAL :: lwave_stress2   ! if .TRUE., calculate wave stress, second call
  LOGICAL :: lgrid_refr      ! if .TRUE., calculate grid refraction

END MODULE mo_wave_config
