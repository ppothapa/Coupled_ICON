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

  USE mo_kind,                 ONLY: wp
  USE mo_exception,            ONLY: finish, message, message_text
  USE mo_impl_constants,       ONLY: max_dom, SUCCESS
  USE mo_math_constants,       ONLY: pi2, rad2deg
  USE mo_physical_constants,   ONLY: grav, rhoh2o
  USE mo_wave_constants,       ONLY: EX_TAIL
  USE mo_fortran_tools,        ONLY: DO_DEALLOCATE
  USE mo_idx_list,             ONLY: t_idx_list1D

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_wave_config'

  ! types
  PUBLIC :: t_wave_config

  ! objects
  PUBLIC :: wave_config

  ! subroutines
  PUBLIC :: configure_wave

  TYPE t_wave_config
    INTEGER  :: ndirs    ! NUMBER OF DIRECTIONS.
    INTEGER  :: nfreqs   ! NUMBER OF FREQUENCIES.

    REAL(wp) :: fr1      ! FIRST FREQUENCY [HZ].
    REAL(wp) :: CO       ! FREQUENCY RATIO
    INTEGER  :: IREF     ! FREQUENCY BIN NUMBER OF REFERENCE FREQUENCY

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
                                   ! 2 - forcing from coupled atmosphere
                                   ! 3 - forcing from data reader

    LOGICAL :: linput_sf1      ! if .TRUE., calculate wind input source function term, first call
    LOGICAL :: linput_sf2      ! if .TRUE., calculate wind input source function term, second call
    LOGICAL :: ldissip_sf      ! if .TRUE., calculate dissipation source function term
    LOGICAL :: lnon_linear_sf  ! if .TRUE., calculate non linear source function term
    LOGICAL :: lbottom_fric_sf ! if .TRUE., calculate bottom_friction source function term
    LOGICAL :: lwave_stress1   ! if .TRUE., calculate wave stress, first call
    LOGICAL :: lwave_stress2   ! if .TRUE., calculate wave stress, second call
    LOGICAL :: lgrid_refr      ! if .TRUE., calculate grid refraction


    ! derived variables and fields
    !
    REAL(wp) ::            &
      &  DELTH,            & ! ANGULAR INCREMENT OF SPECTRUM [RAD].
      &  MO_TAIL,          & ! MO  TAIL FACTOR.
      &  MM1_TAIL,         & ! M-1 TAIL FACTOR.
      &  MP1_TAIL,         & ! M+1 TAIL FACTOR.
      &  MP2_TAIL            ! M+2 TAIL FACTOR.

    REAL(wp), ALLOCATABLE :: &
      &  freqs(:),         & ! frequencies (1:nfreqs) of wave spectrum [Hz]
      &  dfreqs(:),        & ! FREQUENCY INTERVAL (1:nfreqs)
      &  dfreqs_freqs(:),  & ! dfreqs * freqs
      &  dfreqs_freqs2(:), & ! dfreqs * freqs * freqs
      &  dirs(:),          & ! directions (1:ndirs) of wave spectrum [rad]
      &  DFIM(:),          & ! MO  INTEGRATION WEIGHTS.
      &  DFIMOFR(:),       & ! M-1 INTEGRATION WEIGHTS.
      &  DFIM_FR(:),       & ! M+1 INTEGRATION WEIGHTS.
      &  DFIM_FR2(:),      & ! M+2 INTEGRATION WEIGHTS.
      &  RHOWG_DFIM(:)       ! MOMENTUM AND ENERGY FLUX WEIGHTS.

    INTEGER, ALLOCATABLE :: &
      &  freq_ind(:),      & ! index of frequency (1:ntracer=ndirs*nfreq)
      &  dir_ind(:)          ! index of direction (1:ntracer=ndirs*nfreq)

  CONTAINS
    !
    ! get tracer index as a function of direction and frequency index
    PROCEDURE :: get_tracer_id       => wave_config_get_tracer_id
    !
    ! get list of tracer ids for a given frequency index
    PROCEDURE :: get_tracer_ids_freq => wave_config_get_tracer_ids_freq
    !
    ! destruct wave_config object
    PROCEDURE :: destruct            => wave_config_destruct
  END type t_wave_config

  TYPE(t_wave_config), TARGET:: wave_config(max_dom)

CONTAINS

  !>
  !! return tracer index as a function of direction and frequency index
  !!
  !! @par Revision History
  !! Initial revision by Mikhail Dybrynin, DWD (2023-01-30)
  !!
  FUNCTION wave_config_get_tracer_id(me, jdir, jfreq) RESULT(tracer_id)
    CLASS(t_wave_config) :: me
    INTEGER, INTENT(IN)  :: jdir
    INTEGER, INTENT(IN)  :: jfreq

    INTEGER :: tracer_id

    tracer_id = (jfreq-1) * me%ndirs + jdir

  END FUNCTION wave_config_get_tracer_id


  !>
  !! return list of tracer ids for a given frequency index
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2023-01-30)
  !!
  SUBROUTINE wave_config_get_tracer_ids_freq(me, jfreq, list)

    CLASS(t_wave_config)              :: me
    INTEGER,            INTENT(IN)    :: jfreq   !< frequency index
    TYPE(t_idx_list1D), INTENT(INOUT) :: list

    INTEGER :: jt      ! loop index
    INTEGER :: cnt

    IF (.NOT.ALLOCATED(list%idx)) THEN
      CALL list%construct(SIZE(me%freq_ind))
    ENDIF

    cnt = 0
    DO jt=1,SIZE(me%freq_ind)
      IF (me%freq_ind(jt) == jfreq) THEN
        cnt = cnt + 1
        list%idx(cnt) = jt
      ENDIF
    ENDDO

    list%ncount = cnt
  END SUBROUTINE wave_config_get_tracer_ids_freq


  !>
  !! deallocate memory used by object of type t_wave_config
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2023-02-02)
  !!
  SUBROUTINE wave_config_destruct(me)
    CLASS(t_wave_config) :: me

    CALL DO_DEALLOCATE(me%freqs)
    CALL DO_DEALLOCATE(me%dfreqs)
    CALL DO_DEALLOCATE(me%dfreqs_freqs)
    CALL DO_DEALLOCATE(me%dfreqs_freqs2)
    CALL DO_DEALLOCATE(me%dirs)
    CALL DO_DEALLOCATE(me%DFIM)
    CALL DO_DEALLOCATE(me%DFIMOFR)
    CALL DO_DEALLOCATE(me%DFIM_FR)
    CALL DO_DEALLOCATE(me%DFIM_FR2)
    CALL DO_DEALLOCATE(me%RHOWG_DFIM)
    CALL DO_DEALLOCATE(me%freq_ind)
    CALL DO_DEALLOCATE(me%dir_ind)

  END SUBROUTINE wave_config_destruct


  !>
  !! setup the waves model
  !!
  !! Setup of additional waves control variables and constant fields
  !! which depend on the waves-NAMELIST and potentially other namelists.
  !! This routine is called, after all namelists have been read and a
  !! synoptic consistency check has been done.
  !!
  !! @par Revision History
  !! Initial revision by Mikhail Dobrynin, DWD (2023-02-01)
  !!
  SUBROUTINE configure_wave(n_dom, ntracer)

    INTEGER, INTENT(IN) :: n_dom    !< number of domains
    INTEGER, INTENT(IN) :: ntracer  !< number of tracers

    INTEGER :: jd, jf, jfjd   ! loop index
    INTEGER :: jg             ! patch ID
    INTEGER :: ist            ! error status
    TYPE(t_wave_config), POINTER :: wc =>NULL()     ! convenience pointer

    REAL(wp):: CO1

    CHARACTER(*), PARAMETER :: routine = modname//'::configure_waves'


    DO jg=1,n_dom

      ! convenience pointer
      wc => wave_config(jg)

      ALLOCATE(wc%dirs         (wc%ndirs),  &
        &      wc%freqs        (wc%nfreqs), &
        &      wc%dfreqs       (wc%nfreqs), &
        &      wc%dfreqs_freqs (wc%nfreqs), &
        &      wc%dfreqs_freqs2(wc%nfreqs), &
        &      wc%dfim         (wc%nfreqs), &
        &      wc%dfimofr      (wc%nfreqs), &
        &      wc%dfim_fr      (wc%nfreqs), &
        &      wc%dfim_fr2     (wc%nfreqs), &
        &      wc%rhowg_dfim   (wc%nfreqs), &
        &      stat=ist)
      IF (ist/=SUCCESS) CALL finish(routine, "allocation for fields of type REAL failed")

      ALLOCATE(wc%freq_ind     (ntracer),  &
        &      wc%dir_ind      (ntracer),  &
        &      stat=ist)
      IF (ist/=SUCCESS) CALL finish(routine, "allocation for fields of type INTEGER failed")


      !
      ! configuration of spectral setup
      !
      CALL message ('  ','')
      CALL message (':----------------------------------------------------------','')
      WRITE(message_text,'(a,i4)') 'grid ', jg
      CALL message ('Frequencies and directions of wave spectrum',message_text)
      WRITE(message_text,'(a,i4)') 'Number of directions  = ', wc%ndirs
      CALL message ('  ',message_text)
      WRITE(message_text,'(a,i4)') 'Number of frequencies = ', wc%nfreqs
      CALL message ('  ',message_text)
      CALL message (':----------------------------------------------------------','')

      wave_config(jg)%DELTH    = pi2 / wave_config(jg)%ndirs !! ANGULAR INCREMENT OF SPECTRUM [RAD].

      ! calculate directions for wave spectrum
      !
      CALL message ('  ','Directions [Degree]: ')
      DO jd = 1,wc%ndirs
        wc%dirs(jd) = REAL(jd-1) *  wc%DELTH + 0.5_wp * wc%DELTH
        WRITE(message_text,'(i3,f10.5)') jd, wc%dirs(jd)*rad2deg
        CALL message ('  ',message_text)
      END DO
      CALL message ('  ','')

      ! calculate frequencies for wave spectrum
      !
      CALL message ('  ','Frequencies, [Hz]: ')
      WRITE(message_text,'(i3,f10.5)') 1, wc%FR1
      CALL message ('  ',message_text)
      wc%freqs(1) = wc%CO**(-wc%iref + 1._wp) * wc%FR1
      DO jf = 2,wc%nfreqs
        wc%freqs(jf) = wc%CO * wc%freqs(jf-1)
        WRITE(message_text,'(i3,f10.5)') jf, wc%freqs(jf)
        CALL message ('  ',message_text)
      END DO
      CALL message (':----------------------------------------------------------','')


      ! calculate frequency intervals
      !
      CO1 = 0.5_wp * (wc%CO - 1.0_wp)
      wc%dfreqs(1) = CO1 * wc%freqs(1)
      wc%dfreqs(2:wc%nfreqs-1) = CO1 &
        &                      * (wc%freqs(2:wc%nfreqs-1) &
        &                      + (wc%freqs(1:wc%nfreqs-2)))
      wc%dfreqs(wc%nfreqs) = CO1 * wc%freqs(wc%nfreqs-1)
      !
      wc%dfreqs_freqs  = wc%dfreqs * wc%freqs
      wc%dfreqs_freqs2 = wc%dfreqs_freqs * wc%freqs

      ! index of frequency and directions in ntracer loop
      jfjd = 1
      DO jf = 1,wc%nfreqs
        DO jd = 1,wc%ndirs
          wc%freq_ind(jfjd) = jf
          wc%dir_ind(jfjd)  = jd

          jfjd = jfjd + 1
        END DO
      END DO

      ! MO  TAIL FACTOR.
      wc%MO_TAIL  = - wc%DELTH / (EX_TAIL + 1.0_wp) * wc%freqs(wc%nfreqs)

      ! M-1 TAIL FACTOR.
      wc%MM1_TAIL = - wc%DELTH / EX_TAIL

      ! M+1 TAIL FACTOR
      wc%MP1_TAIL = - wc%DELTH / (EX_TAIL + 2.0_wp) * wc%freqs(wc%nfreqs)**2.0_wp

      ! M+2 TAIL FACTOR.
      wc%MP2_TAIL = - wc%DELTH / (EX_TAIL + 3.0_wp) * wc%freqs(wc%nfreqs)**3.0_wp


      ! calculate freqs_dirs parameters
      ! MO INTEGRATION WEIGHTS
      wc%DFIM = wc%dfreqs * wc%DELTH

      ! M-1 INTEGRATION WEIGHTS.
      wc%DFIMOFR = wc%DFIM / wc%freqs

      ! M+1 INTEGRATION WEIGHTS.
      wc%DFIM_FR = wc%dfreqs_freqs * wc%DELTH

      ! M+2 INTEGRATION WEIGHTS.
      wc%DFIM_FR2 = wc%dfreqs_freqs2 * wc%DELTH

      ! MOMENTUM AND ENERGY FLUX WEIGHTS.
      wc%RHOWG_DFIM(:) = rhoh2o * grav * wc%DELTH * LOG(wc%CO) * wc%freqs(:)
      !
      wc%RHOWG_DFIM(1)         = 0.5_wp * wc%RHOWG_DFIM(1)
      wc%RHOWG_DFIM(wc%nfreqs) = 0.5_wp * wc%RHOWG_DFIM(wc%nfreqs)

    ENDDO

  END SUBROUTINE configure_wave

END MODULE mo_wave_config
