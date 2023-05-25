!! Description:  Contains the data structures
!!  for initialisation of the physical model state and other auxiliary variables
!!  in order to run wave physics.
!!
!! @author Mikhail Dobrynin, DWD, 05.09.2019
!!
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
!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_init_wave_physics

  USE mo_kind,                 ONLY: wp
  USE mo_exception,            ONLY: message, message_text
  USE mo_model_domain,         ONLY: t_patch
  USE mo_impl_constants,       ONLY: MAX_CHAR_LENGTH, min_rlcell
  USE mo_physical_constants,   ONLY: grav
  USE mo_math_constants,       ONLY: pi2, rpi_2, deg2rad
  USE mo_loopindices,          ONLY: get_indices_c

  USE mo_wave_types,           ONLY: t_wave_prog, t_wave_diag
  USE mo_wave_forcing_types,   ONLY: t_wave_forcing
  USE mo_wave_config,          ONLY: t_wave_config
  USE mo_wave_ext_data_types,  ONLY: t_external_wave
  USE mo_wave_constants,       ONLY: EMIN
  USE mo_wave_physics,         ONLY: wave_group_velocity_c,   &
    &                                wave_group_velocity_e,   &
    &                                wave_group_velocity_nt,  &
    &                                wave_group_velocity_bnd, &
    &                                wave_number_c, wave_number_e

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_wave_phy

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_init_wave_physics'

CONTAINS

  !>
  !! Initialisation of the wave physics
  !!
  !! 1. Calculation of initial spectrum from the fetch law
  !! and 1D JONSWAP spectum. The minimum of wave energy
  !! is limited to FLMIN.
  !! 2. tba
  !!
  !! @par Revision History
  !! Initial revision by Mikhail Dobrynin, DWD (2019-09-05)
  !!
  SUBROUTINE init_wave_phy(p_patch, wave_config, p_prog, p_diag, wave_ext_data, p_forcing)

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    TYPE(t_wave_prog),           INTENT(INOUT) :: p_prog
    TYPE(t_wave_diag),           INTENT(INOUT) :: p_diag
    TYPE(t_external_wave),       INTENT(IN)    :: wave_ext_data
    TYPE(t_wave_forcing),        INTENT(IN)    :: p_forcing

    TYPE(t_wave_config), POINTER :: wc => NULL()

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = modname//'::init_wave_phy'

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc,jb,jf
    
    ! save some paperwork
    wc => wave_config

    CALL FETCH_LAW(p_patch, wave_config, p_diag, p_forcing)

    ! Set minimum values of energy allowed in the spectrum
    CALL JONSWAP(p_patch,           &
         wc%freqs,                  &
         p_diag%ALPHAJ*0.01_wp,     &
         wc%GAMMA_wave, wc%SIGMA_A, wc%SIGMA_B, &
         p_diag%FP,                 &
         p_diag%FLMINFR)           !OUT

    ! Set JONSWAP spectrum
    CALL JONSWAP(p_patch,           &
         wc%freqs,                  &
         p_diag%ALPHAJ,             &
         wc%GAMMA_wave, wc%SIGMA_A, wc%SIGMA_B, &
         p_diag%FP,                 &
         p_diag%ET)           !OUT

    CALL init_wave_spectrum(p_patch, wc, p_diag, p_forcing, p_prog%tracer)

    ! get wave number as a function of circular frequency and water depth
    ! at cell center
    CALL wave_number_c(p_patch     = p_patch,                    & !IN
      &                wave_config = wave_config,                & !IN
      &                depth       = wave_ext_data%bathymetry_c, & !IN
      &                wave_num_c  = p_diag%wave_num_c)            !OUT

    ! get wave number as a function of circular frequency and water depth
    ! at edge midpoint
    CALL wave_number_e(p_patch     = p_patch,                    & !IN
      &                wave_config = wave_config,                & !IN
      &                depth       = wave_ext_data%bathymetry_e, & !IN
      &                wave_num_e  = p_diag%wave_num_e)            !OUT


    ! compute absolute value of group velocity at cell centers
    !
    CALL wave_group_velocity_c(p_patch, wc, &
         p_diag%wave_num_c,          &  ! IN
         wave_ext_data%bathymetry_c, &  ! IN
         p_diag%gv_c)                   !INOUT

    ! compute absolute value of group velocity at edge midpoints
    !
    CALL wave_group_velocity_e(p_patch, wc, &
         p_diag%wave_num_e,          &  ! IN
         wave_ext_data%bathymetry_e, &  ! IN
         p_diag%gv_e)                   !INOUT

    ! compute normal and tangential components of group velociy vector
    ! at edge midpoints
    CALL wave_group_velocity_nt(p_patch, wc, &
         p_diag%gv_e, &
         p_diag%gvn_e, & !INOUT
         p_diag%gvt_e)   !INOUT

    CALL wave_group_velocity_bnd(p_patch, wc, &
         p_diag%gvn_e)  !INOUT

  END SUBROUTINE init_wave_phy


  !>
  !! Initialisation of the wave spectrum
  !!
  !! Calculation of wind dependent initial spectrum from
  !! the fetch law and from the 1D JONSWAP spectum. The minimum
  !! of wave energy is limited to FLMIN.
  !!
  !! @par Revision History
  !! Initial revision by Mikhail Dobrynin, DWD (2019-09-05)
  !! Vectorization by Daniel Reinert, DWD (2023-01-31)
  !!
  SUBROUTINE init_wave_spectrum(p_patch, wave_config, p_diag, p_forcing, tracer)
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
         &  routine = 'init_wave_spectrum'

    TYPE(t_patch),        INTENT(IN)    :: p_patch
    TYPE(t_wave_config),  INTENT(IN)    :: wave_config
    TYPE(t_wave_diag),    INTENT(IN)    :: p_diag
    TYPE(t_wave_forcing), INTENT(IN)    :: p_forcing
    REAL(wp),             INTENT(INOUT) :: tracer(:,:,:,:)

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc,jb,jd,jf,jk,jt
    REAL(wp):: st

    ! halo points must be included
    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)
    jk         = p_patch%nlev


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jf,jd,jt,jc,i_startidx,i_endidx,st)
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
           &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jf = 1,wave_config%nfreqs
        DO jd = 1,wave_config%ndirs
          !
          jt = wave_config%get_tracer_id(jd,jf)
          !
        DO jc = i_startidx, i_endidx
            st = rpi_2*MAX(0._wp, COS(wave_config%dirs(jd)-p_forcing%dir10m(jc,jb)*deg2rad) )**2
            IF (st < 0.1E-08_wp) st = 0._wp

            ! WAM initialisation
            tracer(jc,jk,jb,jt) = p_diag%ET(jc,jb,jf) * st
            tracer(jc,jk,jb,jt) = MAX(tracer(jc,jk,jb,jt),EMIN)
          END DO  !jc
        END DO  !jd
      END DO  !jf
    END DO  !jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE init_wave_spectrum

  !>
  !! Calculation of the JONSWAP spectrum according to
  !! Hasselmann et al. 1973. Adaptation of WAM 4.5
  !! subroutine JONSWAP.
  !!
  !! @par Revision History
  !! Initial revision by Mikhail Dobrynin, DWD (2019-09-05)
  !! Vectorization by Daniel Reinert, DWD (2023-01-31)
  !!
  SUBROUTINE JONSWAP (p_patch, freqs, ALPHAJ, GAMMA, SA, SB, FP, ET)
    TYPE(t_patch), INTENT(IN)  :: p_patch
    REAL(wp),      INTENT(IN)  :: freqs(:)      !! FREQUENCiIES.
    REAL(wp),      INTENT(IN)  :: ALPHAJ(:,:)   !! OVERALL ENERGY LEVEL OF JONSWAP SPECTRA.
    REAL(wp),      INTENT(IN)  :: GAMMA         !! OVERSHOOT FACTOR.
    REAL(wp),      INTENT(IN)  :: SA            !! LEFT PEAK WIDTH.
    REAL(wp),      INTENT(IN)  :: SB            !! RIGHT PEAK WIDTH.
    REAL(wp),      INTENT(IN)  :: FP(:,:)       !! PEAK FREQUENCIES.
    REAL(wp),      INTENT(OUT) :: ET(:,:,:)     !! JONSWAP SPECTRA.

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
         &  routine = 'JONSWAP'

    REAL(wp), PARAMETER :: FLMIN = 0.000001_wp !! absolute minimum energy in spectral bins
    
    REAL(wp) :: ARG, sigma, G2ZPI4FRH5M

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc,jb,jf

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jf,jc,i_startidx,i_endidx,G2ZPI4FRH5M,sigma,ARG) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
           &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jf = 1,SIZE(freqs)

        G2ZPI4FRH5M = grav**2._wp / pi2**4._wp * freqs(jf)**(-5._wp)

        DO jc = i_startidx, i_endidx

          sigma = MERGE(sb,sa, freqs(jf)>fp(jc,jb))
          ET(jc,jb,jf) = 0._wp

          ARG = 1.25_wp*(FP(jc,jb)/freqs(jf))**4._wp
          IF (ARG.LT.50.0_wp) THEN
            ET(jc,jb,jf) = ALPHAJ(jc,jb) * G2ZPI4FRH5M * EXP(-ARG)
          END IF

          ARG = 0.5_wp*((freqs(jf)-FP(jc,jb)) / (sigma*FP(jc,jb)))**2._wp
          IF (ARG.LT.99._wp) THEN
            ET(jc,jb,jf) = ET(jc,jb,jf)*exp(log(GAMMA)*EXP(-ARG))
          END IF

          ! Avoid too small numbers of p_diag%FLMINFR          
          ET(jc,jb,jf) = MAX(ET(jc,jb,jf),FLMIN)
          
        END DO  !jc
      END DO  !jf
    END DO  !jb
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE JONSWAP

  !>
  !! Calculation of JONSWAP parameters.
  !!
  !! Calculate the peak frequency from a fetch law
  !! and the JONSWAP alpha.
  !!
  !! Developted by S. Hasselmann (July 1990) and H. Guenther (December 1990).
  !! K.HASSELMAN,D.B.ROOS,P.MUELLER AND W.SWELL. A parametric wave prediction
  !! model. Journal of physical oceanography, Vol. 6, No. 2, March 1976.
  !!
  !! Adopted from WAM 4.5.
  !!
  !! @par Revision History
  !! Initial revision by Mikhail Dobrynin, DWD (2019-09-05)
  !! Vectorization by Daniel Reinert, DWD (2023-01-31)
  !!
  SUBROUTINE FETCH_LAW (p_patch, wave_config, p_diag, p_forcing)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
         &  routine = 'FETCH_LAW'
    !
    TYPE(t_patch),        INTENT(IN)    :: p_patch
    TYPE(t_wave_config),  INTENT(IN)    :: wave_config
    TYPE(t_wave_diag),    INTENT(INOUT) :: p_diag
    TYPE(t_wave_forcing), INTENT(IN)    :: p_forcing

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc, jb

    !@waves move to nml?
    REAL(wp), PARAMETER :: A = 2.84_wp,  D = -(3._wp/10._wp) !! PEAK FREQUENCY FETCH LAW CONSTANTS
    REAL(wp), PARAMETER :: B = 0.033_wp, E = 2._wp/3._wp     !! ALPHA-PEAK FREQUENCY LAW CONSTANTS

    REAL(wp) :: UG
    REAL(wp) :: fetch
    REAL(wp) :: fm

    ! halo points must be included !
    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

    fetch = wave_config%fetch
    fm    = wave_config%fm

    ! ---------------------------------------------------------------------------- !
    !                                                                              !
    !     1. COMPUTE VALUES FROM FETCH LAWS.                                       !
    !        -------------------------------
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,UG)
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
           &                 i_startidx, i_endidx, i_rlstart, i_rlend)
      DO jc = i_startidx, i_endidx
        IF (p_forcing%sp10m(jc,jb).GT.0.1E-08_wp) THEN
          UG = grav / p_forcing%sp10m(jc,jb)
          p_diag%FP(jc,jb) = MAX(0.13_wp, A*((grav*fetch)/(p_forcing%sp10m(jc,jb)**2))**D)

          p_diag%FP(jc,jb) = MIN(p_diag%FP(jc,jb), fm/UG)
!!! set min of ALPHAJ to 0 (was 0.0081), otherwise it produces Hs > 1m with wind=0 !!!
          p_diag%ALPHAJ(jc,jb) = MAX(0.0081_wp, B * p_diag%FP(jc,jb)**E)
          p_diag%FP(jc,jb) = p_diag%FP(jc,jb) * UG
        ELSE
          p_diag%ALPHAJ(jc,jb) = 0.0081_wp
          p_diag%FP(jc,jb) = fm
        END IF
      END DO

    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE FETCH_LAW

END MODULE mo_init_wave_physics
