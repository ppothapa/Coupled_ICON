!>
!! Namelist setup for waves.
!! The content of namelists is mostly adopted from the WAM 4.5.4.
!!
!! @author Mikhail Dobrynin, DWD
!!
!! @par Revision History
!! Initial revision by Mikhail Dobrynin, DWD (10.04.2019)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_wave_nml

  USE mo_kind,                ONLY: wp
  USE mo_impl_constants,      ONLY: max_dom
  USE mo_namelist,            ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_master_control,      ONLY: use_restart_namelists
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_restart_nml_and_att, ONLY: open_tmpfile, store_and_close_namelist,     &
    &                               open_and_restore_namelist, close_tmpfile
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings
  USE mo_wave_config,         ONLY: wave_config


  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: read_wave_namelist

CONTAINS
  !>
  !!
  SUBROUTINE read_wave_namelist(filename)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, iunit, funit
    INTEGER :: jg

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_wave_nml:read_wave_namelist'

    INTEGER  :: ndirs    ! NUMBER OF DIRECTIONS.
    INTEGER  :: nfreqs   ! NUMBER OF FREQUENCIES.
    INTEGER  :: IREF     ! FREQUENCY BIN NUMBER OF REFERENCE FREQUENCY

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
    !  REAL(wp)  :: ALPHA      ! MINIMUM CHARNOCK CONSTANT (ECMWF CY45R1).
    ! if LE 30 frequencies changed
    ! to 0.0075 in subroutine INITMDL

    INTEGER :: dt_wave             ! PROPAGATION TIMESTEP !@waves: add units, s?
    INTEGER :: dt_fastphy          ! time step for fast physics processes !@waves: replace dt_wave by dt_fstphy?
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


    NAMELIST /wave_nml/ &
         coldstart, iforc_waves,                            &
         ndirs, nfreqs, fr1, CO, IREF,                      &
         ALPHA, FM, GAMMA_wave, SIGMA_A, SIGMA_B, THETAQ, FETCH,    &
         dt_wave,  dt_fastphy, ROAIR, RNUAIR, RNUAIRM, ROWATER, XEPS, XINVEPS, &
         XKAPPA, XNLEV, BETAMAX, ZALP, &
         linput_sf1, linput_sf2, ldissip_sf, lnon_linear_sf, lbottom_fric_sf, &
         lwave_stress1, lwave_stress2, lgrid_refr

    !-----------------------------------------------------------
    ! 1. default settings
    !-----------------------------------------------------------
    coldstart  = .TRUE.         !! TRUE IF COLDSTART
    ndirs      = 24             !! NUMBER OF DIRECTIONS.
    nfreqs     = 25             !! NUMBER OF FREQUENCIES.
    fr1        = 0.04177248_wp  !! FIRST FREQUENCY [HZ].
    CO         = 1.1_wp         !! FREQUENCY RATIO
    IREF       = 1              !! FREQUENCY BIN NUMBER OF REFERENCE FREQUENCY

    ALPHA      = 0.018_wp       !! PHILLIPS PARAMETER.
    FM         = 0.2_wp         !! PEAK FREQUENCY (HZ) AND/OR MAXIMUM FREQUENCY.
    GAMMA_wave = 3.0_wp         !! OVERSHOOT FACTOR.
    SIGMA_A    = 0.07_wp        !! LEFT PEAK WIDTH.
    SIGMA_B    = 0.09_wp        !! RIGHT PEAK WIDTH.
    THETAQ     = 0.0_wp         !! WAVE DIRECTION (DEG) (NOT USED IF IOPTI = 1).

    FETCH      = 300000._wp     !! FETCH IN METRES.

    ROAIR      = 1.225_wp       !! AIR DENSITY
    RNUAIR     = 1.5E-5_wp      !! KINEMATIC AIR VISCOSITY
    RNUAIRM    = 0.11_wp*RNUAIR !! KINEMATIC AIR VISCOSITY FOR MOMENTUM TRANSFER

    ROWATER    = 1000._wp       !! WATER DENSITY
    XEPS       = ROAIR/ROWATER
    XINVEPS    = 1./XEPS

    BETAMAX    = 1.20_wp        !! PARAMETER FOR WIND INPUT (ECMWF CY45R1).
    ZALP       = 0.0080_wp      !! SHIFTS GROWTH CURVE (ECMWF CY45R1).

    XKAPPA     = 0.40_wp        !! VON KARMAN CONSTANT.
    XNLEV      = 10.0_wp        !! WINDSPEED REF. LEVEL.

    dt_wave    = 600            !! PROPAGATION TIMESTEP, s

    iforc_waves = 1         !! 1 - test case, 2 - forcing from coupled  atmosphere

    linput_sf1 =       .TRUE. !< if .TRUE., calculate wind input source function term, first call
    linput_sf2 =       .TRUE. !< if .TRUE., calculate wind input source function term, second call
    ldissip_sf =       .TRUE. !< if .TRUE., calculate dissipation source function term
    lnon_linear_sf =   .TRUE. !< if .TRUE., calculate non linear source function term
    lbottom_fric_sf =  .TRUE. !< if .TRUE., calculate bottom_friction source function term
    lwave_stress1  =   .TRUE. !< if .TRUE., calculate wave stress, first call
    lwave_stress2  =   .TRUE. !< if .TRUE., calculate wave stress, second call
    lgrid_refr =       .TRUE. !< if .TRUE., calculate grid refraction


    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('wave_nml')
      READ(funit,NML=wave_nml)
      CALL close_tmpfile(funit)
    END IF


    !--------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, wave_nml)   ! write defaults to temporary text file
    END IF

    CALL open_nml(TRIM(filename))
    CALL position_nml ('wave_nml', STATUS=istat)

    SELECT CASE (istat)
    CASE (positioned)
       READ (nnml, wave_nml)                                        ! overwrite default settings
       IF (my_process_is_stdio()) THEN
          iunit = temp_settings()
          WRITE(iunit, wave_nml)    ! write settings to temporary text file
       END IF
    END SELECT
    CALL close_nml

    !----------------------------------------------------
    ! 5. Fill the configuration state
    !----------------------------------------------------

    DO jg=1,max_dom
      wave_config(jg)%ndirs           = ndirs
      wave_config(jg)%nfreqs          = nfreqs
      wave_config(jg)%fr1             = fr1
      wave_config(jg)%CO              = CO
      wave_config(jg)%IREF            = IREF
      wave_config(jg)%ALPHA           = ALPHA
      wave_config(jg)%FM              = FM
      wave_config(jg)%GAMMA_wave      = GAMMA_wave
      wave_config(jg)%SIGMA_A         = SIGMA_A
      wave_config(jg)%SIGMA_B         = SIGMA_B
      wave_config(jg)%THETAQ          = THETAQ
      wave_config(jg)%FETCH           = FETCH
      wave_config(jg)%ROAIR           = ROAIR
      wave_config(jg)%RNUAIR          = RNUAIR
      wave_config(jg)%RNUAIRM         = RNUAIRM
      wave_config(jg)%ROWATER         = ROWATER
      wave_config(jg)%XEPS            = XEPS
      wave_config(jg)%XINVEPS         = XINVEPS
      wave_config(jg)%XKAPPA          = XKAPPA
      wave_config(jg)%XNLEV           = XNLEV
      wave_config(jg)%BETAMAX         = BETAMAX
      wave_config(jg)%ZALP            = ZALP
      wave_config(jg)%coldstart       = coldstart
      wave_config(jg)%iforc_waves     = iforc_waves
      wave_config(jg)%linput_sf1      = linput_sf1
      wave_config(jg)%linput_sf2      = linput_sf2
      wave_config(jg)%ldissip_sf      = ldissip_sf
      wave_config(jg)%lnon_linear_sf  = lnon_linear_sf
      wave_config(jg)%lbottom_fric_sf = lbottom_fric_sf
      wave_config(jg)%lwave_stress1   = lwave_stress1
      wave_config(jg)%lwave_stress2   = lwave_stress2
      wave_config(jg)%lgrid_refr      = lgrid_refr
    ENDDO


    !-----------------------------------------------------
    ! 6. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=wave_nml)
      CALL store_and_close_namelist(funit, 'wave_nml')
    ENDIF

    !------------------------------------------------------------------
    ! 7. Write the namelist to an ASCII file
    !------------------------------------------------------------------
    IF ( my_process_is_stdio() ) WRITE(nnml_output,nml=wave_nml)

  END SUBROUTINE read_wave_namelist

END MODULE mo_wave_nml
