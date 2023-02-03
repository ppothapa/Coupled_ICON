!! Contains the subroutines with wave physics parametrisation
!!
!! @author Mikhail Dobrynin
!!
!! @par Revision History
!! Initial revision by Mikhail Dobrynin, DWD  (2019.09.05)
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

MODULE mo_wave_physics

  USE mo_model_domain,        ONLY: t_patch
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, min_rlcell
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_kind,                ONLY: wp
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: dtime
  USE mo_physical_constants,  ONLY: grav
  USE mo_math_constants,      ONLY: deg2rad, dbl_eps
  USE mo_idx_list,            ONLY: t_idx_list1D

  USE mo_wave_types,          ONLY: t_wave_prog, t_wave_diag
  USE mo_wave_config,         ONLY: t_wave_config
  USE mo_wave_constants,      ONLY: ALPHA, EPS1, EMIN

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: air_sea
  PUBLIC :: new_spectrum
  PUBLIC :: total_energy
  PUBLIC :: mean_frequency_energy


  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_wave_physics'

CONTAINS

  !>
  !! Calculation of total stress and sea surface roughness
  !!
  !! Adaptation of WAM 4.5 algorithm and code for calculation of total stress
  !! and sea surface roughness for ICON-waves (P.A.E.M. Janssen, 1990).
  !!
  !! @par Revision History
  !! Initial revision by Mikhail Dobrynin, DWD (2019-09-05)
  !! Vectorization by Daniel Reinert, DWD (2023-01-30)
  !!
  SUBROUTINE air_sea(p_patch, wave_config, p_diag ) !(UTOP, TAUW, USTAR, Z0)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'air_sea'

    TYPE(t_patch),       INTENT(IN)    :: p_patch
    TYPE(t_wave_diag),   INTENT(INOUT) :: p_diag
    TYPE(t_wave_config), INTENT(IN)    :: wave_config

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc,jb,iter

    REAL(wp), PARAMETER :: TWOXMP1 = 3.0_wp
    REAL(wp), PARAMETER :: EPSUS   = 1.0E-6_wp

    !     *ACD*       COEFFICIENTS FOR SIMPLE CD(U10) RELATION
    !     *BCD*       CD = ACD + BCD*U10
    REAL(wp), PARAMETER :: ACD = 8.0E-4_wp
    REAL(wp), PARAMETER :: BCD = 8.0E-5_wp

    INTEGER, PARAMETER :: NITER = 15

    REAL(wp):: xkutop(nproma), tauold(nproma), ustm1(nproma), z0ch(nproma)
    REAL(wp):: xlogxl, alphaog, xologz0
    REAL(wp):: ustold, taunew, x, f, delf
    REAL(wp):: z0tot, z0vis, zz
    LOGICAL :: l_converged(nproma)

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

    xlogxl  = LOG(wave_config%XNLEV)
    alphaog = ALPHA / grav

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,iter,i_startidx,i_endidx,xkutop,ustold,tauold,ustm1, &
!$OMP            l_converged,x,z0ch,z0vis,z0tot,xologz0,f,zz,delf,taunew) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      ! initialization
      DO jc = i_startidx, i_endidx
        xkutop(jc) = wave_config%XKAPPA * p_diag%sp_10m(jc,jb) !xkutop = XKAPPA*UTOP(IJ)
        ustold     = p_diag%sp_10m(jc,jb) * SQRT(ACD + BCD * p_diag%sp_10m(jc,jb))
        tauold(jc) = MAX(ustold**2.0_wp,p_diag%tauw(jc,jb)+EPS1) !MAX(ustold**2,TAUW(IJ)+EPS1)
        p_diag%ustar(jc,jb) = SQRT(tauold(jc))
        ustm1(jc) = 1.0_wp/MAX(p_diag%ustar(jc,jb),EPSUS)
        !
        l_converged(jc) = .FALSE.
      ENDDO

      DO iter = 1,NITER
        DO jc = i_startidx, i_endidx

          x        = p_diag%tauw(jc,jb) / tauold(jc)
          z0ch(jc) = alphaog * tauold(jc) / SQRT(MAX(1.0_wp-x,EPS1))
          z0vis    = wave_config%RNUAIRM * ustm1(jc)
          z0tot    = z0ch(jc) + z0vis

          xologz0 = 1.0 / (xlogxl - LOG(z0tot))
          f = p_diag%ustar(jc,jb) - xkutop(jc) * xologz0
          zz = ustm1(jc) &
               * (z0ch(jc) * (2.0_wp-TWOXMP1*X) / (1.0_wp -x) - z0vis) &
               / z0tot
          delf = 1.0_wp - xkutop(jc) * xologz0**2.0_wp * zz

          p_diag%ustar(jc,jb) = p_diag%ustar(jc,jb) - f / delf

          taunew = MAX(p_diag%ustar(jc,jb)**2.0_wp,p_diag%tauw(jc,jb) + EPS1)

          p_diag%ustar(jc,jb) = SQRT(taunew)

          IF (ABS(taunew-tauold(jc))<= dbl_eps) l_converged(jc) = .TRUE.
          ustm1(jc)  = 1.0_wp/MAX(p_diag%ustar(jc,jb),EPSUS)
          tauold(jc) = taunew
        END DO  !jc
        !
        IF ( ALL(l_converged(i_startidx:i_endidx)) ) EXIT
        !
      END DO  !iter
      !
      DO jc = i_startidx, i_endidx
        p_diag%z0(jc,jb) = z0ch(jc)
      ENDDO
    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE air_sea




  !>
  !! Calculation of mean frequency energy
  !!
  !! Integration over frequencies for calculation of
  !! of mean frequency energy. Adaptation of WAM 4.5 code
  !! of the subroutine FEMEAN developed by S.D. HASSELMANN,
  !! optimized by L. Zambresky and H. Guenther, GKSS, 2001                              !
  !!
  !! @par Revision History
  !! Initial revision by Mikhail Dobrynin, DWD (2019-09-05)
  !! Optimization and vectorization by Daniel Reinert, DWD (2023-01-30)
  !!
  SUBROUTINE mean_frequency_energy(p_patch, wave_config, p_prog, p_diag)
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: routine = 'mean_frequency'

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    TYPE(t_wave_prog),           INTENT(IN)    :: p_prog
    TYPE(t_wave_diag),           INTENT(INOUT) :: p_diag

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc,jb,jf,jk
    INTEGER :: jt                       !< tracer index
    INTEGER :: n                        !< loop index
    TYPE(t_idx_list1D) :: list_tr       !< list of tracer ids

    REAL(wp) :: temp(nproma,wave_config%nfreqs), temp_1(nproma,wave_config%nfreqs)
    TYPE(t_wave_config), POINTER :: wc => NULL()

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)
    jk         = p_patch%nlev

    ! save some paperwork
    wc => wave_config


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jf,n,jt,i_startidx,i_endidx,list_tr,temp,temp_1) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jf = 1,wc%nfreqs

        CALL list_tr%construct(SIZE(wc%freq_ind))
!!$        CALL get_tracer_idx_list(wc%freq_ind, jf, list_tr)
        CALL wc%get_tracer_ids_freq(jf,list_tr)

        DO jc = i_startidx, i_endidx
          temp(jc,jf)   = 0._wp
          temp_1(jc,jf) = 0._wp
        ENDDO

        DO n=1,list_tr%ncount
          jt = list_tr%idx(n)
          DO jc = i_startidx, i_endidx
            temp(jc,jf) = temp(jc,jf) + p_prog%tracer(jc,jk,jb,jt)
            IF (p_diag%LLWS(jc,jb,jt)==1) THEN
              temp_1(jc,jf) = temp_1(jc,jf) + p_prog%tracer(jc,jk,jb,jt)
            ENDIF
          ENDDO
        ENDDO  ! n

        CALL list_tr%finalize()
      END DO  ! jf

      DO jc = i_startidx, i_endidx
        p_diag%femean(jc,jb)   = wc%MM1_TAIL * temp(jc,wc%nfreqs)
        p_diag%femeanws(jc,jb) = wc%MM1_TAIL * temp_1(jc,wc%nfreqs)
      ENDDO

      DO jf = 1,wc%nfreqs
        DO jc = i_startidx, i_endidx
          p_diag%femean(jc,jb)   = p_diag%femean(jc,jb) + temp(jc,jf) * wc%DFIMOFR(jf)
          p_diag%femeanws(jc,jb) = p_diag%femeanws(jc,jb) + temp_1(jc,jf) * wc%DFIMOFR(jf)
        END DO
      END DO

      DO jc = i_startidx, i_endidx
        p_diag%femean(jc,jb)   = p_diag%emean(jc,jb) / MAX(p_diag%femean(jc,jb),EMIN)
        p_diag%femeanws(jc,jb) = p_diag%emeanws(jc,jb) / MAX(p_diag%femeanws(jc,jb),EMIN)
      END DO
    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL
  END SUBROUTINE mean_frequency_energy

  !>
  !! Calculation of total energy
  !!
  !! Calculation of total energy by integtation over directions and frequencies.
  !! A tail correction is added.
  !! Adaptation of WAM 4.5 code of the subroutine TOTAL_ENERGY
  !! developed by S.D. HASSELMANN, optimized by L. Zambresky
  !! and H. Guenther, GKSS, 2001
  !!
  !! @par Revision History
  !! Initial revision by Mikhail Dobrynin, DWD (2019-09-05)
  !! Optimization and vectorization by Daniel Reinert, DWD (2023-01-30)
  !!
  SUBROUTINE total_energy(p_patch, wave_config, p_prog, p_diag)
    ! in  p_patch, p_prog%tracer, wave_config%freq_ind, wave_config%MO_TAIL
    !out  p_diag%emean(jc,jb), p_diag%emeanws(jc,jb)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = modname//'total_energy'

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    TYPE(t_wave_prog),           INTENT(IN)    :: p_prog
    TYPE(t_wave_diag),           INTENT(INOUT) :: p_diag

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc,jb,jf,jk
    INTEGER :: jt                       !< tracer index
    INTEGER :: n                        !< loop index
    TYPE(t_idx_list1D) :: list_tr       !< list of tracer ids

    REAL(wp):: sum1(nproma,wave_config%nfreqs), sum2(nproma,wave_config%nfreqs)
    TYPE(t_wave_config), POINTER :: wc => NULL()

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)
    jk         = p_patch%nlev

    ! save some paperwork
    wc => wave_config

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jf,jt,n,i_startidx,i_endidx,list_tr,sum1,sum2) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)


      ! compute sum of all tracers that match a specific frequency
      DO jf = 1,wc%nfreqs

        CALL list_tr%construct(SIZE(wc%freq_ind))
!!$        CALL get_tracer_idx_list(wc%freq_ind, jf, list_tr)
        CALL wc%get_tracer_ids_freq(jf, list_tr)

        ! initialization
        DO jc = i_startidx, i_endidx
          sum1(jc,jf) = 0._wp
          sum1(jc,jf) = 0._wp
        ENDDO

        DO n=1,list_tr%ncount
          jt = list_tr%idx(n)
          DO jc = i_startidx, i_endidx
            sum1(jc,jf) = sum1(jc,jf) + p_prog%tracer(jc,jk,jb,jt)
            IF (p_diag%LLWS(jc,jb,jt) == 1) THEN
              sum2(jc,jf) = sum2(jc,jf) + p_prog%tracer(jc,jk,jb,jt)
            ENDIF
          ENDDO
        ENDDO  ! n

        CALL list_tr%finalize()
      ENDDO  ! jf

      ! initialization
      DO jc = i_startidx, i_endidx
        p_diag%emean(jc,jb)   = wc%MO_TAIL * sum1(jc,wc%nfreqs)
        p_diag%emeanws(jc,jb) = wc%MO_TAIL * sum2(jc,wc%nfreqs)
      ENDDO

      ! sum over all frequencies
      DO jf = 1,wc%nfreqs
        DO jc = i_startidx, i_endidx
          p_diag%emean(jc,jb)   = p_diag%emean(jc,jb)   + sum1(jc,jf) * wc%DFIM(jf)
          p_diag%emeanws(jc,jb) = p_diag%emeanws(jc,jb) + sum2(jc,jf) * wc%DFIM(jf)
        END DO
      ENDDO  ! jf

      ! clipping
      DO jc = i_startidx, i_endidx
        p_diag%emean(jc,jb)   = MAX(p_diag%emean(jc,jb),EMIN)
        p_diag%emeanws(jc,jb) = MAX(p_diag%emeanws(jc,jb),EMIN)
      END DO
    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE total_energy


  !>
  !! Calculation of new spectrum.
  !!
  !! Increase of spectrum in a time step and limit to a finite
  !! fraction of a typical F**(-4) equilibrium spectrum.
  !! Adaptation of WAM 4.5 code.
  !!
  !! @par Revision History
  !! Initial revision by Mikhail Dobrynin, DWD (2019-09-05)
  !! Optimization and vectorization by Daniel Reinert, DWD (2023-01-30)
  !!
  SUBROUTINE new_spectrum(p_patch, wave_config, p_prog, p_diag)
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'new_spectrum'

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    TYPE(t_wave_prog),           INTENT(INOUT) :: p_prog
    TYPE(t_wave_diag),           INTENT(IN)    :: p_diag

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb,jc,jf,jd,jt,jk

    REAL(wp) :: temp_1, temp_2, temp_3
    REAL(wp) :: sprd
    REAL(wp) :: delfl
    TYPE(t_wave_config), POINTER :: wc => NULL()

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)
    jk         = p_patch%nlev

    ! save some paperwork
    wc => wave_config

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jb,jf,jd,jt,i_startidx,i_endidx,delfl,sprd,temp_1,temp_2,temp_3) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jf = 1,wc%nfreqs
        DO jd = 1,wc%ndirs
          !
          jt = wc%get_tracer_id(jd,jf)
          !
          DO jc = i_startidx, i_endidx

            delfl = 5.0E-07_wp * grav / wc%freqs(jf)**4 * dtime

            temp_2 = p_diag%ustar(jc,jb) * delfl &
                 &  * MAX(p_diag%femeanws(jc,jb),p_diag%femean(jc,jb))

            sprd = MAX(0._wp, COS(wc%dirs(jd)-p_diag%dir_10m(jc,jb)*deg2rad) )**2

            temp_1 = dtime * p_diag%sl(jc,jb,jt) &
                 / MAX(1._wp, 1._wp -  dtime * p_diag%fl(jc,jb,jt))

            temp_3 = MIN(ABS(temp_1),temp_2)

            p_prog%tracer(jc,jk,jb,jt) = p_prog%tracer(jc,jk,jb,jt)  + SIGN(temp_3,temp_1)
            p_prog%tracer(jc,jk,jb,jt) = MAX(p_prog%tracer(jc,jk,jb,jt), p_diag%FLMINFR(jc,jb,jf)*sprd)

          END DO
        END DO
      END DO
    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE new_spectrum

END MODULE mo_wave_physics
