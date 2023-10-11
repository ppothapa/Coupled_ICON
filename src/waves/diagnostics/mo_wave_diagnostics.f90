!>
!! Computes diagnostic parameters and some diagnostics in the wave model
!!
!!
!! @par Revision History
!! Initial release by Mikhail Dorbynin, DWD (2021-02-18)
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
!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_wave_diagnostics
  USE mo_kind,                ONLY: wp
  USE mo_model_domain,        ONLY: t_patch
  USE mo_wave_config,         ONLY: t_wave_config
  USE mo_wave_types,          ONLY: t_wave_diag
  USE mo_impl_constants,      ONLY: min_rlcell
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_physical_constants,  ONLY: grav
  USE mo_math_constants,      ONLY: pi2
  USE mo_parallel_config,     ONLY: nproma
  USE mo_kind,                ONLY: wp
  USE mo_fortran_tools,       ONLY: init

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname   = 'mo_wave_diagnostics'

  PUBLIC :: calculate_output_diagnostics


CONTAINS
  !>
  !! Calculation of diagnostic parameters
  !!
  !! @par Revision History
  !! Initial revision by Mikhail Dobrynin, DWD (2019-09-22)
  !!
  SUBROUTINE calculate_output_diagnostics(p_patch, wave_config, tracer, depth, p_diag)

    TYPE(t_patch),         INTENT(IN)    :: p_patch
    TYPE(t_wave_config),   INTENT(IN)    :: wave_config
    REAL(wp),              INTENT(IN)    :: tracer(:,:,:,:) !energy spectral bins
    REAL(wp),              INTENT(IN)    :: depth(:,:)
    TYPE(t_wave_diag),     INTENT(INOUT) :: p_diag

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = modname//':calculate_diagnostics'

      CALL significant_wave_height(p_patch = p_patch, &
        &                          emean   = p_diag%emean(:,:), &
        &                          hs      = p_diag%hs(:,:))

      CALL stokes_drift(p_patch = p_patch, &
        &           wave_config = wave_config, &
        &            wave_num_c = p_diag%wave_num_c, &
        &                 depth = depth, &
        &                tracer = tracer, &
        &              u_stokes = p_diag%u_stokes, & ! OUT
        &              v_stokes = p_diag%v_stokes)   ! OUT

  END SUBROUTINE calculate_output_diagnostics


  !>
  !! TO DO: Description
  !!
  !! @par Revision History
  !! Initial revision by Mikhail Dobrynin, DWD (2019-09-22)
  !!
  SUBROUTINE significant_wave_height(p_patch, emean, hs)

    TYPE(t_patch),     INTENT(IN)    :: p_patch
    REAL(wp),          INTENT(IN)    :: emean(:,:)  !< total energy [m^2]
    REAL(wp),          INTENT(INOUT) :: hs(:,:)     !< significant wave height [m]

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = modname//':significant_wave_height'

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc,jb


    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jb,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)
      DO jc = i_startidx, i_endidx
        hs(jc,jb) = 4.0_wp * SQRT(emean(jc,jb))
      END DO
    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL
  END SUBROUTINE significant_wave_height

  !>
  !! Calculation of Stokes drift components
  !!
  !! Adaptation of WAM 4.5 code of the subroutine STOKES_DRIFT
  !! developed by M.REISTAD, O.SAETRA, and H.GUNTHER
  !!
  !! Reference
  !! Kern E. Kenton, JGR, Vol 74 NO 28, 1969
  !!
  !! @par Revision History
  !! Initial revision by Mikhail Dobrynin, DWD (2019-09-18)
  !!
  SUBROUTINE stokes_drift(p_patch, wave_config, wave_num_c, depth, tracer, u_stokes, v_stokes)

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = modname//':stokes_drift'

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    REAL(wp),                    INTENT(IN)    :: wave_num_c(:,:,:)  !< wave number (1/m)
    REAL(wp),                    INTENT(IN)    :: depth(:,:)
    REAL(wp),                    INTENT(IN)    :: tracer(:,:,:,:) !energy spectral bins
    REAL(wp),                    INTENT(INOUT) :: u_stokes(:,:)
    REAL(wp),                    INTENT(INOUT) :: v_stokes(:,:)

    TYPE(t_wave_config), POINTER :: wc => NULL()

    REAL(wp) :: ak, akd, fact, tailfac
    REAL(wp) :: si(nproma), ci(nproma)

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc,jb,jf,jd,jk,jt

    wc => wave_config

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)
    jk         = p_patch%nlev

!$OMP PARALLEL
    CALL init(u_stokes)
    CALL init(v_stokes)
!$OMP BARRIER
!$OMP DO PRIVATE(jb,jc,jf,jd,jt,i_startidx,i_endidx,ak,akd,si,ci,fact,tailfac) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
           &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      ! initialisation of si, ci
      DO jc = i_startidx, i_endidx
        si(jc) = 0._wp
        ci(jc) = 0._wp
      END DO

      freqs:DO jf = 1,wc%nfreqs
        DO jd = 1, wc%ndirs
          jt =  wc%tracer_ind(jd,jf)
          DO jc = i_startidx, i_endidx
            si(jc) = si(jc) + tracer(jc,jk,jb,jt) * SIN(wc%dirs(jd))
            ci(jc) = ci(jc) + tracer(jc,jk,jb,jt) * COS(wc%dirs(jd))
          END DO
        END DO

        DO jc = i_startidx, i_endidx
          ak = wave_num_c(jc,jb,jf)
          akd = ak * depth(jc,jb)
          fact = 2._wp*grav*ak**2/(pi2*wc%freqs(jf)*TANH(2._wp*akd)) * wc%DFIM(jf)
          si(jc) = fact * si(jc)
          ci(jc) = fact * ci(jc)
          u_stokes(jc,jb) = u_stokes(jc,jb) + si(jc)
          v_stokes(jc,jb) = v_stokes(jc,jb) + ci(jc)
        END DO

      END DO freqs

      DO jc = i_startidx, i_endidx
        tailfac = wc%freqs(wc%nfreqs)**2 / &
          &      (wc%dfreqs(wc%nfreqs) * (wc%freqs(wc%nfreqs)+0.5_wp*wc%dfreqs(wc%nfreqs)))
        u_stokes(jc,jb) = u_stokes(jc,jb) + tailfac * si(jc)
        v_stokes(jc,jb) = v_stokes(jc,jb) + tailfac * ci(jc)
      END DO

    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE stokes_drift

END MODULE mo_wave_diagnostics
