!>
!! Namelist setup for waves. 
!! The content of namelists is mostly adopted from the WAM 4.5.4.
!!
!! @author <Mikhail Dobrynin, DWD, 10.04.2019>
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
  USE mo_namelist,            ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_io_units,            ONLY: nnml, nnml_output, find_next_free_unit
  USE mo_master_control,      ONLY: use_restart_namelists
  !  USE mo_restart_namelist,    ONLY: open_tmpfile, store_and_close_namelist, &
  !                                  & open_and_restore_namelist, close_tmpfile
  USE mo_exception,           ONLY: finish
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC  :: read_wave_namelist
  
  CHARACTER(LEN=80)  :: HEADER            !! HEADER OF MODEL RUN.
  
  INTEGER, PUBLIC  :: ndirs    ! NUMBER OF DIRECTIONS.
  INTEGER, PUBLIC  :: nfreqs   ! NUMBER OF FREQUENCIES.
  
  REAL(wp), PUBLIC :: fr1      ! FIRST FREQUENCY [HZ].
  REAL(wp), PUBLIC :: CO       ! FREQUENCY RATIO
  
  REAL(wp), PUBLIC :: ALPHA      ! PHILLIPS' PARAMETER  (NOT USED IF IOPTI = 1)
  REAL(wp), PUBLIC :: FM         ! PEAK FREQUENCY (HZ) AND/OR MAXIMUM FREQUENCY
  REAL(wp), PUBLIC :: GAMMA_wave ! OVERSHOOT FACTOR
  REAL(wp), PUBLIC :: SIGMA_A    ! LEFT PEAK WIDTH
  REAL(wp), PUBLIC :: SIGMA_B    ! RIGHT PEAK WIDTH
  REAL(wp), PUBLIC :: THETAQ     ! WAVE DIRECTION (DEG) (NOT USED IF IOPTI = 1)
  REAL(wp), PUBLIC :: FETCH      ! FETCH IN METRES (IF ZERO THEN 0.5 OF THE LATITUDE INCREMENT IS USED.).

  REAL(wp), PUBLIC :: ROAIR   ! AIR DENSITY
  REAL(wp), PUBLIC :: RNUAIR  ! KINEMATIC AIR VISCOSITY
  REAL(wp), PUBLIC :: RNUAIRM ! KINEMATIC AIR VISCOSITY FOR MOMENTUM TRANSFER
  REAL(wp), PUBLIC :: ROWATER ! WATER DENSITY
  REAL(wp), PUBLIC :: XEPS
  REAL(wp), PUBLIC :: XINVEPS

  REAL(wp), PUBLIC :: XKAPPA  ! VON KARMAN CONSTANT.
  REAL(wp), PUBLIC :: XNLEV   ! WINDSPEED REF. LEVEL.
  REAL(wp), PUBLIC :: BETAMAX ! PARAMETER FOR WIND INPUT (ECMWF CY45R1).
  REAL(wp), PUBLIC :: ZALP    ! SHIFTS GROWTH CURVE (ECMWF CY45R1).
!  REAL(wp)  :: ALPHA      ! MINIMUM CHARNOCK CONSTANT (ECMWF CY45R1).
                           ! if LE 30 frequencies changed
                           ! to 0.0075 in subroutine INITMDL

  INTEGER, PUBLIC :: dt_wave             ! PROPAGATION TIMESTEP !@waves: add units, s?
  INTEGER, PUBLIC :: dt_fastphy          ! time step for fast physics processes !@waves: replace dt_wave by dt_fstphy?
  LOGICAL, PUBLIC :: coldstart           ! if .TRUE. start from initialisation without restart file
  INTEGER, PUBLIC :: iforc_waves         ! 1 - test case
                                 ! 2 - forcing from coupled atmosphete
                                 ! 3 - forcing from data reader

  LOGICAL, PUBLIC :: linput_sf1      ! if .TRUE., calculate wind input source function term, first call
  LOGICAL, PUBLIC :: linput_sf2      ! if .TRUE., calculate wind input source function term, second call
  LOGICAL, PUBLIC :: ldissip_sf      ! if .TRUE., calculate dissipation source function term
  LOGICAL, PUBLIC :: lnon_linear_sf  ! if .TRUE., calculate non linear source function term
  LOGICAL, PUBLIC :: lbottom_fric_sf ! if .TRUE., calculate bottom_friction source function term
  LOGICAL, PUBLIC :: lwave_stress1   ! if .TRUE., calculate wave stress, first call
  LOGICAL, PUBLIC :: lwave_stress2   ! if .TRUE., calculate wave stress, second call
  LOGICAL, PUBLIC :: lgrid_refr      ! if .TRUE., calculate grid refraction


  NAMELIST /wave_nml/ &
       HEADER, coldstart, iforc_waves,                            &
       ndirs, nfreqs, fr1, CO,                                    &
       ALPHA, FM, GAMMA_wave, SIGMA_A, SIGMA_B, THETAQ, FETCH,    &
       dt_wave,  dt_fastphy, ROAIR, RNUAIR, RNUAIRM, ROWATER, XEPS, XINVEPS, &
       XKAPPA, XNLEV, BETAMAX, ZALP, &
       linput_sf1, linput_sf2, ldissip_sf, lnon_linear_sf, lbottom_fric_sf, &
       lwave_stress1, lwave_stress2, lgrid_refr
  
CONTAINS
  !>
  !!
  SUBROUTINE read_wave_namelist(filename)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit, iunit
    
    CHARACTER(len=*), PARAMETER ::  &
         &  routine = 'mo_wave_nml:read_wave_namelist'
    
    !------------------------------------------------------------------
    ! Set default values
    !------------------------------------------------------------------
    HEADER = ' '                !! HEADER OF MODEL RUN.
    coldstart  = .TRUE.         !! TRUE IF COLDSTART
    ndirs      = 24             !! NUMBER OF DIRECTIONS.
    nfreqs     = 25             !! NUMBER OF FREQUENCIES.
    fr1        = 0.04177248_wp  !! FIRST FREQUENCY [HZ].
    CO         = 1.1_wp         !! FREQUENCY RATIO
    ALPHA      = 0.018_wp       !! PHILLIPS PARAMETER.
    FM         = 0.2_wp         !! PEAK FREQUENCY (HZ) AND/OR MAXIMUM FREQUENCY.
    GAMMA_wave = 3.0_wp         !! OVERSHOOT FACTOR.
    SIGMA_A    = 0.07_wp        !! LEFT PEAK WIDTH.
    SIGMA_B    = 0.09_wp        !! RIGHT PEAK WIDTH.
    THETAQ     = 0.0_wp         !! WAVE DIRECTION (DEG) (NOT USED IF IOPTI = 1).
    FETCH      = 300000._wp      !! FETCH IN METRES.
    
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

    CALL open_nml(TRIM(filename))
    CALL position_nml ('wave_nml', STATUS=istat)
    IF (my_process_is_stdio()) THEN
       iunit = temp_defaults()
       WRITE(iunit, wave_nml)    ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (positioned)
       READ (nnml, wave_nml)                                        ! overwrite default settings
       IF (my_process_is_stdio()) THEN
          iunit = temp_settings()
          WRITE(iunit, wave_nml)    ! write settings to temporary text file
       END IF
    END SELECT
    CALL close_nml
    
    !------------------------------------------------------------------
    ! Write the namelist to an ASCII file
    !------------------------------------------------------------------
    IF ( my_process_is_stdio() ) WRITE(nnml_output,nml=wave_nml)
    
  END SUBROUTINE read_wave_namelist

END MODULE mo_wave_nml
