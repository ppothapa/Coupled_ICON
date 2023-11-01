! Contains the subroutines with wave physics parametrisation
!
!
! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_wave_physics

  USE mo_model_domain,        ONLY: t_patch
  USE mo_grid_config,         ONLY: grid_sphere_radius
  USE mo_exception,           ONLY: message
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, min_rlcell, min_rledge
  USE mo_loopindices,         ONLY: get_indices_c, get_indices_e
  USE mo_kind,                ONLY: wp, vp
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: dtime
  USE mo_physical_constants,  ONLY: grav
  USE mo_math_constants,      ONLY: deg2rad, dbl_eps, pi, pi2
  USE mo_fortran_tools,       ONLY: init

  USE mo_wave_types,          ONLY: t_wave_diag
  USE mo_wave_config,         ONLY: t_wave_config
  USE mo_wave_constants,      ONLY: EPS1, EMIN, DELTA, CONSS

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: air_sea, bottom_friction, nonlinear_transfer
  PUBLIC :: input_source_function, dissipation_source_function
  PUBLIC :: last_prog_freq_ind
  PUBLIC :: impose_high_freq_tail
  PUBLIC :: tm1_tm2_periods
  PUBLIC :: wm1_wm2_wavenumber
  PUBLIC :: new_spectrum
  PUBLIC :: total_energy
  PUBLIC :: wave_stress
  PUBLIC :: mean_frequency_energy, set_energy2emin
  PUBLIC :: wave_group_velocity_c
  PUBLIC :: wave_group_velocity_e
  PUBLIC :: wave_group_velocity_nt
  PUBLIC :: wave_group_velocity_bnd
  PUBLIC :: wave_number_c
  PUBLIC :: wave_number_e
  PUBLIC :: wave_refraction


  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_wave_physics'

CONTAINS

  !>
  !! Calculation of wave group velocity
  !!
  !! Calculation of shallow water cell centered
  !! wave group velocity
  !!
  SUBROUTINE wave_group_velocity_c(p_patch, p_config, wave_num_c, bathymetry_c, gv_c)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
         &  routine = modname//':wave_group_velocity_c'

    TYPE(t_patch),       INTENT(IN)   :: p_patch
    TYPE(t_wave_config), INTENT(IN)   :: p_config
    REAL(wp),            INTENT(IN)   :: wave_num_c(:,:,:) !< wave number (1/m)
    REAL(wp),            INTENT(IN)   :: bathymetry_c(:,:) !< bathymetric height at cell centers (nproma,nblks_c) ( m )
    REAL(wp),            INTENT(INOUT):: gv_c(:,:,:)       !< group velocity (nproma,nlev,nblks_c,ndirs*nfreqs)  ( m/s )

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb,jf,jc
    INTEGER :: nfreqs

    REAL(wp) :: gh, ak, akd, gv

    gh = grav / (4.0_wp * pi)

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

    nfreqs = p_config%nfreqs

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jf,jc,i_startidx,i_endidx,ak,akd,gv) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jf = 1, nfreqs

        DO jc = i_startidx, i_endidx
          ! shallow water group velocity
          ak = wave_num_c(jc,jb,jf)
          akd = ak * bathymetry_c(jc,jb)

          IF (akd <= 10.0_wp) THEN
            gv = 0.5_wp * SQRT(grav * TANH(akd)/ak) * (1.0_wp + 2.0_wp*akd/SINH(2.0_wp*akd))
          ELSE
            gv = gh / p_config%freqs(jf)
          END IF
          gv_c(jc,jb,jf) = gv
        END DO
      END DO
    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL
  END SUBROUTINE wave_group_velocity_c



  !>
  !! Calculation of wave group velocity
  !!
  !! Calculation of shallow water edges centered
  !! wave group velocities
  !!
  SUBROUTINE wave_group_velocity_e(p_patch, p_config, wave_num_e, bathymetry_e, gv_e)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
         &  routine = modname//':wave_group_velocity_e'

    TYPE(t_patch),       INTENT(IN)   :: p_patch
    TYPE(t_wave_config), INTENT(IN)   :: p_config
    REAL(wp),            INTENT(IN)   :: wave_num_e(:,:,:) !< wave number (1/m)
    REAL(wp),            INTENT(IN)   :: bathymetry_e(:,:) !< bathymetric height at cell centers (nproma,nblks_c) ( m )
    REAL(wp),            INTENT(INOUT):: gv_e(:,:,:)     !< group velocity (nproma,nlev,nblks_c,ndirs*nfreqs)  ( m/s )

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb,jf,je
    INTEGER :: nfreqs

    REAL(wp) :: gh, ak, akd, gv

    gh = grav / (4.0_wp * pi)

    i_rlstart  = 1
    i_rlend    = min_rledge
    i_startblk = p_patch%edges%start_block(i_rlstart)
    i_endblk   = p_patch%edges%end_block(i_rlend)

    nfreqs = p_config%nfreqs

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jf,je,i_startidx,i_endidx,ak,akd,gv) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_e( p_patch, jb, i_startblk, i_endblk,           &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jf = 1, nfreqs

        DO je = i_startidx, i_endidx
          ! shallow water group velocity
          ak = wave_num_e(je,jb,jf)
          akd = ak * bathymetry_e(je,jb)

          IF (akd <= 10.0_wp) THEN
            gv = 0.5_wp * SQRT(grav * TANH(akd)/ak) * (1.0_wp + 2.0_wp*akd/SINH(2.0_wp*akd))
          ELSE
            gv = gh / p_config%freqs(jf)
          END IF
          gv_e(je,jb,jf) = gv
        END DO
      END DO
    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE wave_group_velocity_e

  !>
  !! Calculation of wave group velocity
  !!
  !! Calculation of shallow water
  !! edge-normal and -tangential projections of
  !! wave group velocities using spectral directions
  !!
  SUBROUTINE wave_group_velocity_nt(p_patch, p_config, gv_e, gvn_e, gvt_e)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
         &  routine = modname//':wave_group_velocity_nt'

    TYPE(t_patch),       INTENT(IN)   :: p_patch
    TYPE(t_wave_config), INTENT(IN)   :: p_config
    REAL(wp),            INTENT(IN)   :: gv_e(:,:,:) !< group velocity (nproma,nblks_e,nfreqs)  ( m/s )
    REAL(wp),            INTENT(INOUT):: gvn_e(:,:,:)!< normal group velocity (nproma,nblks_e,dirs*nfreqs)  ( m/s )
    REAL(wp),            INTENT(INOUT):: gvt_e(:,:,:)!< tangential group velocity (nproma,nblks_e,ndirs*nfreqs)  ( m/s )

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb,jf,jd,jt,je
    INTEGER :: nfreqs, ndirs

    REAL(wp) :: gvu, gvv

    i_rlstart  = 1
    i_rlend    = min_rledge
    i_startblk = p_patch%edges%start_block(i_rlstart)
    i_endblk   = p_patch%edges%end_block(i_rlend)

    nfreqs = p_config%nfreqs
    ndirs  = p_config%ndirs

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jf,jd,jt,je,i_startidx,i_endidx,gvu,gvv) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_e( p_patch, jb, i_startblk, i_endblk,           &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)
      DO jf = 1, nfreqs
        DO jd = 1, ndirs
          jt = p_config%tracer_ind(jd,jf)
          DO je = i_startidx, i_endidx
            gvu = gv_e(je,jb,jf) * SIN(p_config%dirs(jd))
            gvv = gv_e(je,jb,jf) * COS(p_config%dirs(jd))

            gvn_e(je,jb,jt) = &
                 gvu * p_patch%edges%primal_normal(je,jb)%v1 + &
                 gvv * p_patch%edges%primal_normal(je,jb)%v2

            gvt_e(je,jb,jt) = &
                 gvu * p_patch%edges%dual_normal(je,jb)%v1 + &
                 gvv * p_patch%edges%dual_normal(je,jb)%v2
          END DO
        END DO
      END DO
    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL
  END SUBROUTINE wave_group_velocity_nt



  !>
  !! Correction of normal to edge group velocity,
  !! avoiding of wave energy propagation from land
  !! and insuring full "outflow" of wave energy towards land.
  !!
  !! Set the wave group velocity to zero at the boundary edge
  !! in case of wave energy propagation towards the ocean,
  !! and set gn = deep water group velocity otherwise
  !!
  !! we make use of the fact that the edge-normal velocity vector points
  !! * towards the coast, if
  !!   cells%edge_orientation > 0 .AND. vn > 0
  !!   OR
  !!   cells%edge_orientation < 0 .AND. vn < 0
  !! * towards the sea, if
  !!   cells%edge_orientation > 0 .AND. vn < 0
  !!   OR
  !!   cells%edge_orientation < 0 .AND. vn > 0
  !!
  !! hence:
  !! * towards the coast, if (cells%edge_orientation * vn) > 0
  !! * towards the sea,   if (cells%edge_orientation * vn) < 0
  !!
  SUBROUTINE wave_group_velocity_bnd(p_patch, p_config, gvn_e)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = modname//':wave_group_velocity_bnd'

    TYPE(t_patch),       INTENT(IN)   :: p_patch
    TYPE(t_wave_config), INTENT(IN)   :: p_config
    REAL(wp),            INTENT(INOUT):: gvn_e(:,:,:)!< normal group velocity (nproma,nblks_e,dirs*nfreqs)  ( m/s )

    ! local variables
    REAL(wp):: gv
    INTEGER :: jb, jc, jf, jd, jt, ic
    INTEGER :: jce                   !< loop index for cell edges
    INTEGER :: eidx, eblk            !< edge index and block
    INTEGER :: jje, jjb              !< line and block index of boundary edge
    INTEGER :: i_rlstart_c, i_rlend_c
    INTEGER :: i_startblk_c, i_endblk_c
    INTEGER :: i_startidx_c, i_endidx_c
    INTEGER :: nfreqs, ndirs
    INTEGER :: ile(nproma), ibe(nproma)
    REAL(wp):: e_orient(nproma)
    INTEGEr :: cnt
    LOGICAL :: is_towards_coastline  !< TRUE if normal component of group velocity vector
                                     !  at land-sea boundary points towards coastline

    nfreqs = p_config%nfreqs
    ndirs  = p_config%ndirs

    ! set up loop over boundary cells (refine_c_ctrl==1)
    i_rlstart_c  = 1
    i_rlend_c    = 1
    i_startblk_c = p_patch%cells%start_block(i_rlstart_c)
    i_endblk_c   = p_patch%cells%end_block(i_rlend_c)

    DO jb = i_startblk_c, i_endblk_c
      CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c, &
        &                i_startidx_c, i_endidx_c, i_rlstart_c, i_rlend_c)

      cnt = 0
      !
      DO jc = i_startidx_c, i_endidx_c

        ! build list of coastline edges (refin_e_ctrl==1)
        ! and store edge orientation.
        !
        DO jce =1,3
          eidx = p_patch%cells%edge_idx(jc,jb,jce)
          eblk = p_patch%cells%edge_blk(jc,jb,jce)

          IF (p_patch%edges%refin_ctrl(eidx,eblk) == 1) THEN
            ! coastline edge found
            cnt = cnt + 1
            ile(cnt) = eidx
            ibe(cnt) = eblk
            e_orient(cnt) = p_patch%cells%edge_orientation(jc,jb,jce)
          ENDIF
        ENDDO

      ENDDO  !jc

      ! Correction of normal to edge group velocity, avoiding of wave energy propagation from land
      ! and insuring full "outflow" of wave energy towards land.
      !
      DO jf = 1, nfreqs
        ! deep water group velocity
        gv = grav / (2.0_wp * pi2 * p_config%freqs(jf))
        DO jd = 1, ndirs
          jt = p_config%tracer_ind(jd,jf)
!$NEC ivdep
          DO ic = 1, cnt
            jje = ile(ic)
            jjb = ibe(ic)
            is_towards_coastline = (e_orient(ic) * gvn_e(jje,jjb,jt)) > 0._wp
            gvn_e(jje,jjb,jt) = MERGE(gv, 0.0_wp, is_towards_coastline)
          ENDDO  !jc
        ENDDO  !jd
      ENDDO  !jf

    ENDDO  !jb

  END SUBROUTINE wave_group_velocity_bnd

!!$  !>
!!$  !! Correction of normal to edge group velocity,
!!$  !! avoiding of wave energy propagation from land
!!$  !! and insuring full "outflow" of wave energy towards land.
!!$  !!
!!$  !! Set the wave group velocity to zero at the boundary edge
!!$  !! in case of wave energy propagation towards the ocean,
!!$  !! and set gn = deep water group velocity otherwise
!!$  !!
!!$  SUBROUTINE wave_group_velocity_bnd(p_patch, p_config, gvn_e)
!!$
!!$    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
!!$         &  routine = modname//':wave_group_velocity_bnd'
!!$
!!$    TYPE(t_patch),       INTENT(IN)   :: p_patch
!!$    TYPE(t_wave_config), INTENT(IN)   :: p_config
!!$    REAL(wp),            INTENT(INOUT):: gvn_e(:,:,:)!< normal group velocity (nproma,nblks_e,dirs*nfreqs)  ( m/s )
!!$
!!$    INTEGER :: i_rlstart_c, i_rlend_c
!!$    INTEGER :: i_rlstart_e, i_rlend_e
!!$    INTEGER :: i_startblk_c, i_endblk_c
!!$    INTEGER :: i_startblk_e, i_endblk_e
!!$    INTEGER :: i_startidx_c, i_endidx_c
!!$    INTEGER :: i_startidx_e, i_endidx_e
!!$    INTEGER :: jb,jf,jd,jt,je,jc,ie_local,iel
!!$    INTEGER :: nfreqs, ndirs
!!$
!!$    REAL(wp) :: gv
!!$
!!$    nfreqs = p_config%nfreqs
!!$    ndirs  = p_config%ndirs
!!$
!!$    !TO DO: make index for edges with corrected gn in base/mo_wave_toolbox.f90
!!$    !TO DO: make function for deep_water_gv(freq) (gv = group velocity)
!!$    !                gv = grav / (2.0_wp * pi2 * freqs))
!!$
!!$    i_rlstart_e  = 1
!!$    i_rlend_e    = min_rledge
!!$    i_startblk_e   = p_patch%edges%start_block(i_rlstart_e)
!!$    i_endblk_e     = p_patch%edges%end_block(i_rlend_e)
!!$
!!$    i_rlstart_c = 1
!!$    i_rlend_c   = min_rlcell
!!$    i_startblk_c = p_patch%cells%start_block(i_rlstart_c)
!!$    i_endblk_c   = p_patch%cells%end_block(i_rlend_c)
!!$
!!$    ie_local = 1
!!$    DO jb = i_startblk_e, i_endblk_e
!!$      CALL get_indices_e(p_patch, jb, i_startblk_e, i_endblk_e, &
!!$           &                i_startidx_e, i_endidx_e, i_rlstart_e, i_rlend_e)
!!$      CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c, &
!!$           &                i_startidx_c, i_endidx_c, i_rlstart_c, i_rlend_c)
!!$
!!$      DO je = i_startidx_e, i_endidx_e
!!$        IF (p_patch%edges%refin_ctrl(je,jb).eq.1) THEN !or .le. 3
!!$          !IF (ANY(p_patch%edges%cell_idx(je,jb,1:2) .le. 0)) THEN
!!$          DO jc = i_startidx_c, i_endidx_c
!!$            DO iel = 1, 3
!!$              IF ((p_patch%cells%edge_idx(jc,jb,iel) == je) &
!!$                   .AND. (p_patch%cells%edge_blk(jc,jb,iel) == jb)) THEN
!!$                ie_local = iel
!!$                DO jf = 1, nfreqs
!!$                  !deep water group velocity
!!$                  gv = grav / (2.0_wp * pi2 * p_config%freqs(jf))
!!$                  DO jd = 1, ndirs
!!$                    jt = p_config%get_tracer_id(jd,jf)
!!$                    If (p_patch%cells%edge_orientation(jc,jb,ie_local).gt.0) THEN ! > 0 - outwards
!!$                      IF (gvn_e(je,jb,jt).ge.0.0_wp) THEN ! > 0 - outwards
!!$                        ! towards coastline
!!$                        gvn_e(je,jb,jt) = gv
!!$                      ELSE
!!$                        ! towards ocean
!!$                        gvn_e(je,jb,jt) = 0.0_wp
!!$                      END IF
!!$                    ELSE IF (p_patch%cells%edge_orientation(jc,jb,ie_local).lt.0) THEN ! < 0 - inwards
!!$                      IF (gvn_e(je,jb,jt).ge.0.0_wp) THEN !  0 - inwards
!!$                        ! towards ocean
!!$                        gvn_e(je,jb,jt) = 0.0_wp
!!$                      ELSE
!!$                        ! towards coastline
!!$                        gvn_e(je,jb,jt) = gv
!!$                      END IF
!!$                    END IF
!!$                  END DO
!!$                END DO
!!$              END IF
!!$            END DO
!!$          END DO
!!$        END IF
!!$      END DO
!!$    END DO
!!$  END SUBROUTINE wave_group_velocity_bnd

  !>
  !! Calculation of total stress and sea surface roughness
  !!
  !! Adaptation of WAM 4.5 algorithm and code for calculation of total stress
  !! and sea surface roughness for ICON-waves (P.A.E.M. Janssen, 1990).
  !!
  SUBROUTINE air_sea(p_patch, wave_config, wsp10m, tauw, ustar, z0)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
         & routine =  modname//'air_sea'

    TYPE(t_patch),      INTENT(IN)  :: p_patch
    TYPE(t_wave_config),INTENT(IN)  :: wave_config
    REAL(wp),           INTENT(IN)  :: wsp10m(:,:)!10m wind speed (nproma,nblks_c) ( m/s )
    REAL(wp),           INTENT(IN)  :: tauw(:,:)  !wave stress (nproma,nblks_c) ( (m/s)^2 )
    REAL(wp),           INTENT(OUT) :: ustar(:,:) !friction velocity (nproma,nblks_c) ( m/s )
    REAL(wp),           INTENT(OUT) :: z0(:,:)    !roughness length (nproma,nblks_c) ( m )

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
    alphaog = wave_config%ALPHA_CH / grav

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,iter,i_startidx,i_endidx,xkutop,ustold,tauold,ustm1, &
!$OMP            l_converged,x,z0ch,z0vis,z0tot,xologz0,f,zz,delf,taunew) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      ! initialization
      DO jc = i_startidx, i_endidx
        xkutop(jc) = wave_config%XKAPPA * wsp10m(jc,jb)
        ustold     = wsp10m(jc,jb) * SQRT(ACD + BCD * wsp10m(jc,jb))
        tauold(jc) = MAX(ustold**2,tauw(jc,jb)+EPS1)
        ustar(jc,jb) = SQRT(tauold(jc))
        ustm1(jc) = 1.0_wp/MAX(ustar(jc,jb),EPSUS)

        l_converged(jc) = .FALSE.
      END DO

      DO iter = 1,NITER
        DO jc = i_startidx, i_endidx
          x        = tauw(jc,jb) / tauold(jc)
          z0ch(jc) = alphaog * tauold(jc) / SQRT(MAX(1.0_wp-x,EPS1))
          z0vis    = wave_config%RNUAIRM * ustm1(jc)
          z0tot    = z0ch(jc) + z0vis

          xologz0 = 1.0 / (xlogxl - LOG(z0tot))
          f = ustar(jc,jb) - xkutop(jc) * xologz0
          zz = ustm1(jc) &
               * (z0ch(jc) * (2.0_wp-TWOXMP1*X) / (1.0_wp -x) - z0vis) &
               / z0tot
          delf = 1.0_wp - xkutop(jc) * xologz0**2 * zz

          ustar(jc,jb) = ustar(jc,jb) - f / delf

          taunew = MAX(ustar(jc,jb)**2,tauw(jc,jb) + EPS1)

          ustar(jc,jb) = SQRT(taunew)

          IF (ABS(taunew-tauold(jc))<= dbl_eps) l_converged(jc) = .TRUE.

          ustm1(jc)  = 1.0_wp/MAX(ustar(jc,jb),EPSUS)
          tauold(jc) = taunew
        END DO  !jc
        !
        IF ( ALL(l_converged(i_startidx:i_endidx)) ) EXIT

      END DO  !iter
      !
      DO jc = i_startidx, i_endidx
        z0(jc,jb) = z0ch(jc)
      END DO
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
  SUBROUTINE mean_frequency_energy(p_patch, wave_config, tracer, llws, emean, emeanws, femean, femeanws)
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
         & routine =  modname//'mean_frequency_energy'

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    REAL(wp), INTENT(IN)  :: tracer(:,:,:,:) !energy spectral bins (nproma,nlev,nblks_c,ntracer)
    INTEGER,  INTENT(IN)  :: llws(:,:,:)   !=1 where wind_input is positive (nproma,nblks_c,ntracer)
    REAL(wp), INTENT(IN)  :: emean(:,:)    !total energy (nproma,nblks_c)
    REAL(wp), INTENT(IN)  :: emeanws(:,:)  !total windsea energy (nproma,nblks_c)
    REAL(wp), INTENT(INOUT) :: femean(:,:)   !mean frequency energy (nproma,nblks_c)
    REAL(wp), INTENT(INOUT) :: femeanws(:,:) !mean windsea frequency energy (nproma,nblks_c)

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc,jb,jf,jk
    INTEGER :: jt                       !< tracer index
    INTEGER :: n                        !< loop index

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
!$OMP DO PRIVATE(jb,jc,jf,n,jt,i_startidx,i_endidx,temp,temp_1) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jf = 1,wc%nfreqs

        DO jc = i_startidx, i_endidx
          temp(jc,jf)   = 0._wp
          temp_1(jc,jf) = 0._wp
        ENDDO

        DO n=1,SIZE(wc%list_tr(jf)%p)
          jt = wc%list_tr(jf)%p(n)
          DO jc = i_startidx, i_endidx
            temp(jc,jf) = temp(jc,jf) + tracer(jc,jk,jb,jt)
            IF (llws(jc,jb,jt)==1) THEN
              temp_1(jc,jf) = temp_1(jc,jf) + tracer(jc,jk,jb,jt)
            ENDIF
          ENDDO
        ENDDO  ! n

      END DO  ! jf

      DO jc = i_startidx, i_endidx
        femean(jc,jb)   = wc%MM1_TAIL * temp(jc,wc%nfreqs)
        femeanws(jc,jb) = wc%MM1_TAIL * temp_1(jc,wc%nfreqs)
      END DO

      DO jf = 1,wc%nfreqs
        DO jc = i_startidx, i_endidx
          femean(jc,jb)   = femean(jc,jb) + temp(jc,jf) * wc%DFIMOFR(jf)
          femeanws(jc,jb) = femeanws(jc,jb) + temp_1(jc,jf) * wc%DFIMOFR(jf)
        END DO
      END DO

      DO jc = i_startidx, i_endidx
        femean(jc,jb)   = emean(jc,jb) / MAX(femean(jc,jb),EMIN)
        femeanws(jc,jb) = emeanws(jc,jb) / MAX(femeanws(jc,jb),EMIN)
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
  SUBROUTINE total_energy(p_patch, wave_config, tracer, llws, emean, emeanws)
    ! in  p_patch, p_prog%tracer, wave_config%freq_ind, wave_config%MO_TAIL
    !out  p_diag%emean(jc,jb), p_diag%emeanws(jc,jb)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
         &  routine = modname//'total_energy'

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    REAL(wp), INTENT(IN)    :: tracer(:,:,:,:) !energy spectral bins (nproma,nlev,nblks_c,ntracer)
    INTEGER,  INTENT(IN)    :: llws(:,:,:)   !=1 where wind_input is positive (nproma,nblks_c,ntracer)
    REAL(wp), INTENT(INOUT) :: emean(:,:)    !total energy (nproma,nblks_c)
    REAL(wp), INTENT(INOUT) :: emeanws(:,:)  !total windsea energy (nproma,nblks_c)


    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc,jb,jf,jk
    INTEGER :: jt                       !< tracer index
    INTEGER :: n                        !< loop index

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
!$OMP DO PRIVATE(jb,jc,jf,jt,n,i_startidx,i_endidx,sum1,sum2) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)


      ! compute sum of all tracers that match a specific frequency
      DO jf = 1,wc%nfreqs

        ! initialization
        DO jc = i_startidx, i_endidx
          sum1(jc,jf) = 0._wp
          sum2(jc,jf) = 0._wp
        ENDDO

        DO n=1,SIZE(wc%list_tr(jf)%p)
          jt = wc%list_tr(jf)%p(n)
          DO jc = i_startidx, i_endidx
            sum1(jc,jf) = sum1(jc,jf) + tracer(jc,jk,jb,jt)
            IF (llws(jc,jb,jt) == 1) THEN
              sum2(jc,jf) = sum2(jc,jf) + tracer(jc,jk,jb,jt)
            ENDIF
          ENDDO
        ENDDO  ! n

      ENDDO  ! jf

      ! initialization
      DO jc = i_startidx, i_endidx
        emean(jc,jb)   = wc%MO_TAIL * sum1(jc,wc%nfreqs)
        emeanws(jc,jb) = wc%MO_TAIL * sum2(jc,wc%nfreqs)
      ENDDO

      ! sum over all frequencies
      DO jf = 1,wc%nfreqs
        DO jc = i_startidx, i_endidx
          emean(jc,jb)   = emean(jc,jb)   + sum1(jc,jf) * wc%DFIM(jf)
          emeanws(jc,jb) = emeanws(jc,jb) + sum2(jc,jf) * wc%DFIM(jf)
        END DO
      ENDDO  ! jf

      ! clipping
      DO jc = i_startidx, i_endidx
        emean(jc,jb)   = MAX(emean(jc,jb),EMIN)
        emeanws(jc,jb) = MAX(emeanws(jc,jb),EMIN)
      END DO
    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE total_energy


  !>
  !! Calculation of wave stress.
  !!
  !! Compute normalized wave stress from input source function
  !!
  !! Adaptation of WAM 4.5 code.
  !! STRESSO
  !!     H. GUNTHER      GKSS/ECMWF  NOVEMBER  1989 CODE MOVED FROM SINPUT
  !!     P.A.E.M. JANSSEN      KNMI  AUGUST    1990
  !!     J. BIDLOT             ECMWF FEBRUARY  1996-97
  !!     H. GUENTHER   GKSS  FEBRUARY 2002       FT 90
  !!     J. BIDLOT             ECMWF           2007  ADD MIJ
  !!     P.A.E.M. JANSSEN     ECMWF            2011  ADD FLUX CALULATIONS
  !!
  !! Reference
  !!       R SNYDER ET AL,1981.
  !!       G. KOMEN, S. HASSELMANN AND K. HASSELMANN, JPO, 1984.
  !!       P. JANSSEN, JPO, 1985
  !!
  SUBROUTINE wave_stress(p_patch, wave_config, p_diag, dir10m, tracer)
     CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
          &  routine = modname//'wave_stress'

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    TYPE(t_wave_diag),           INTENT(INOUT) :: p_diag
    REAL(wp),                    INTENT(IN)    :: dir10m(:,:)
    REAL(wp),                    INTENT(IN)    :: tracer(:,:,:,:)

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc,jb,jf,jd,jtd,jk

    REAL(wp) :: gm1, const, sinplus, cosw
    REAL(wp) :: cmrhowgdfth
    REAL(wp) :: const1(nproma), const2(nproma)
    REAL(wp) :: rhowgdfth(nproma,wave_config%nfreqs)
    REAL(wp) :: cm(nproma,wave_config%nfreqs)
    REAL(wp) :: xstress(nproma), xstress_tot
    REAL(wp) :: ystress(nproma), ystress_tot
    REAL(wp) :: temp1(nproma), temp2(nproma)
    REAL(wp) :: sumt(nproma), sumx(nproma), sumy(nproma)
    REAL(wp) :: roair    ! air density

    TYPE(t_wave_config), POINTER :: wc => NULL()

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)
    jk         = p_patch%nlev

    ! save some paperwork
    wc => wave_config

    gm1   = 1.0_wp/grav
    const = wc%delth*(pi2)**4*gm1
    roair = MAX(wc%roair,1._wp)


!$OMP PARALLEL
    CALL init(p_diag%phiaw)
    CALL init(p_diag%tauw)
!$OMP BARRIER
!$OMP DO PRIVATE(jb,jc,jf,jd,jtd,i_startidx,i_endidx,cm,const1,const2,         &
!$OMP            rhowgdfth,sinplus,sumt,sumx,sumy,cmrhowgdfth,xstress,ystress, &
!$OMP            xstress_tot,ystress_tot,cosw,temp1,temp2) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
           &                 i_startidx, i_endidx, i_rlstart, i_rlend)
      DO jf = 1,wc%nfreqs
        DO jc = i_startidx, i_endidx
          cm(jc,jf) = p_diag%wave_num_c(jc,jb,jf) * 1.0_wp/(pi2*wc%freqs(jf))
        ENDDO
      END DO


      DO jc = i_startidx, i_endidx

        rhowgdfth(jc,1:p_diag%last_prog_freq_ind(jc,jb)) = &
             wc%rhowg_dfim(1:p_diag%last_prog_freq_ind(jc,jb))

        IF (p_diag%last_prog_freq_ind(jc,jb).NE.wc%nfreqs)  &
             rhowgdfth(jc,p_diag%last_prog_freq_ind(jc,jb)) = &
             0.5_wp * rhowgdfth(jc,p_diag%last_prog_freq_ind(jc,jb))

        rhowgdfth(jc,p_diag%last_prog_freq_ind(jc,jb)+1:wc%nfreqs) = 0.0_wp

        !initialisation
        xstress(jc) = 0._wp
        ystress(jc) = 0._wp

      END DO


      !sum
      DO jf = 1, MAXVAL(p_diag%last_prog_freq_ind(i_startidx:i_endidx,jb))
        DO jc = i_startidx, i_endidx
          sumt(jc) = 0._wp
          sumx(jc) = 0._wp
          sumy(jc) = 0._wp
        END DO

        DO jd = 1, wc%ndirs
          jtd = wc%tracer_ind(jd,jf)
          DO jc = i_startidx, i_endidx
            sinplus = MAX(p_diag%sl(jc,jb,jtd),0._wp)
            sumt(jc) = sumt(jc) + sinplus
            sumx(jc) = sumx(jc) + sinplus * SIN(wc%dirs(jd))
            sumy(jc) = sumy(jc) + sinplus * COS(wc%dirs(jd))
          END DO
        END DO

        DO jc = i_startidx, i_endidx
          p_diag%phiaw(jc,jb) =  p_diag%phiaw(jc,jb) + sumt(jc)*rhowgdfth(jc,jf)
          cmrhowgdfth = cm(jc,jf) * rhowgdfth(jc,jf)
          xstress(jc) = xstress(jc) + sumx(jc)*cmrhowgdfth
          ystress(jc) = ystress(jc) + sumy(jc)*cmrhowgdfth
        END DO
      END DO  ! jf



      ! calculate high-frequency contribution to stress
      !
      DO jc = i_startidx, i_endidx
        temp1(jc)  = 0._wp
        temp2(jc)  = 0._wp
        const1(jc) = const * wc%freqs(p_diag%last_prog_freq_ind(jc,jb))**5 * gm1
        const2(jc) = roair * const * wc%freqs(p_diag%last_prog_freq_ind(jc,jb))**5
      ENDDO

      DO jd = 1, wc%ndirs
        DO jc = i_startidx, i_endidx
          jtd = wc%tracer_ind(jd,p_diag%last_prog_freq_ind(jc,jb))

          cosw = MAX(COS(wc%dirs(jd)-dir10m(jc,jb)*deg2rad),0.0_wp)
          temp1(jc) = temp1(jc) + tracer(jc,jk,jb,jtd) * cosw**3
          temp2(jc) = temp2(jc) + tracer(jc,jk,jb,jtd) * cosw**2
        END DO
      END DO

      CALL high_frequency_stress(wave_config        = wave_config,                     & !IN
        &                        i_startidx         = i_startidx,                      & !IN
        &                        i_endidx           = i_endidx,                        & !IN
        &                        last_prog_freq_ind = p_diag%last_prog_freq_ind(:,jb), & !IN
        &                        ustar              = p_diag%ustar(:,jb),              & !IN
        &                        z0                 = p_diag%z0(:,jb),                 & !IN
        &                        xlevtail           = p_diag%xlevtail(:,jb),           & !IN
        &                        tauhf1             = p_diag%tauhf1(:,jb),             & !INOUT
        &                        phihf1             = p_diag%phihf1(:,jb) )              !INOUT

      DO jc = i_startidx, i_endidx
        p_diag%tauhf(jc,jb) = const1(jc)*temp1(jc)*p_diag%tauhf1(jc,jb)
        p_diag%phihf(jc,jb) = const2(jc)*temp2(jc)*p_diag%phihf1(jc,jb)

        p_diag%phiaw(jc,jb) = p_diag%phiaw(jc,jb) + p_diag%phihf(jc,jb)

        xstress_tot = xstress(jc)/roair + p_diag%tauhf(jc,jb)*SIN(dir10m(jc,jb)*deg2rad)
        ystress_tot = ystress(jc)/roair + p_diag%tauhf(jc,jb)*COS(dir10m(jc,jb)*deg2rad)

        p_diag%tauw(jc,jb) = SQRT(xstress_tot**2+ystress_tot**2)
        p_diag%tauw(jc,jb) = MIN(p_diag%tauw(jc,jb),p_diag%ustar(jc,jb)**2 - EPS1)
        p_diag%tauw(jc,jb) = MAX(p_diag%tauw(jc,jb),0.0_wp)

      END DO

    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE wave_stress


  !>
  !! Calculation of high frequency stress.
  !!
  !! Adaptation of WAM 4.5 code.
  !!
  SUBROUTINE high_frequency_stress(wave_config, i_startidx, i_endidx, last_prog_freq_ind, &
    &                              ustar, z0, xlevtail, tauhf1, phihf1)
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
         &  routine = modname//'high_frequency_stress'

    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    INTEGER,                     INTENT(IN)    :: i_startidx, i_endidx
    INTEGER,                     INTENT(IN)    :: last_prog_freq_ind(:) !<last frequency index of the prognostic range
    REAL(wp),                    INTENT(IN)    :: ustar(:)              !< friction velocity
    REAL(wp),                    INTENT(IN)    :: z0(:)                 !< roughness length
    REAL(wp),                    INTENT(IN)    :: xlevtail(:)           !< tail level
    REAL(wp),                    INTENT(INOUT) :: tauhf1(:)             !< high-frequency stress
    REAL(wp),                    INTENT(INOUT) :: phihf1(:)             !< high-frequency energy flux into ocean

    INTEGER :: jc, j
    REAL(wp):: gm1, x0g
    REAL(wp):: OMEGA, OMEGAC, OMEGACC
    REAL(wp):: UST(nproma), UST0, TAUW0, tauw(nproma)
    REAL(wp):: YC, Y, CM1, ZX, ZARG, ZLOG, ZBETA
    REAL(wp):: DELZ(nproma), ZINF(nproma)
    REAL(wp):: FNC2, SQRTZ0OG(nproma), GZ0, SQRTGZ0(nproma), XLOGGZ0(nproma)

    TYPE(t_wave_config), POINTER :: wc => NULL()
    REAL(wp), PARAMETER :: ZSUP = 0.0_wp  !  LOG(1.)

    wc => wave_config

    gm1 = 1.0_wp/grav
    x0g = wc%X0TAUHF * grav

    DO jc = i_startidx, i_endidx
      OMEGAC      = pi2 * wc%freqs(last_prog_freq_ind(jc))
      UST0        = ustar(jc)
      TAUW0       = UST0**2.0_wp
      GZ0         = grav * z0(jc)
      OMEGACC     = MAX(OMEGAC,X0G/UST0)

      XLOGGZ0(jc) = LOG(GZ0)
      SQRTZ0OG(jc)= SQRT(z0(jc)*GM1)
      SQRTGZ0(jc) = 1.0_wp / SQRTZ0OG(jc)
      YC          = OMEGACC * SQRTZ0OG(jc)
      ZINF(jc)    = LOG(YC)
      DELZ(jc)    = MAX((ZSUP-ZINF(jc))/REAL(wc%JTOT_TAUHF-1),0.0_wp)

      tauw(jc)    = TAUW0
      UST(jc)     = UST0

      tauhf1(jc)  = 0.0_wp
      phihf1(jc)  = 0.0_wp
    ENDDO


    ! Integrals are integrated following a change of variable : Z=LOG(Y)
    DO J = 1, wc%JTOT_TAUHF
      DO jc = i_startidx, i_endidx
        Y         = EXP(ZINF(jc)+REAL(J-1)*DELZ(jc))
        OMEGA     = Y * SQRTGZ0(jc)
        CM1       = OMEGA * GM1
        ZX        = UST(jc) * CM1 + wc%zalp
        ZARG      = wc%XKAPPA/ZX
        ZLOG      = XLOGGZ0(jc) + 2.0_wp * LOG(CM1)+ZARG
        ZLOG      = MIN(ZLOG,0.0_wp)
        ZBETA     = EXP(ZLOG) * ZLOG**4

        FNC2      = ZBETA * tauw(jc) * wc%WTAUHF(J) * DELZ(jc)
        tauw(jc)  = MAX(tauw(jc)-xlevtail(jc) * FNC2,0.0_wp)
        UST(jc)   = SQRT(tauw(jc))

        tauhf1(jc) = tauhf1(jc) + FNC2
        phihf1(jc) = phihf1(jc) + FNC2/Y
      END DO
    END DO

    DO jc = i_startidx, i_endidx
      phihf1(jc) = SQRTZ0OG(jc) * phihf1(jc)
    END DO

  END SUBROUTINE high_frequency_stress


  !>
  !! Calculation of TM1 and TM2 periods.
  !!
  !! Adaptation of WAM 4.5 code.
  !! TM1_TM2_PERIODS_B
  !! Integration of spectra and adding of tail factors.
  !!
  SUBROUTINE tm1_tm2_periods(p_patch, wave_config, tracer, emean, tm1, tm2, f1mean)
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
         &  routine = modname//'tm1_tm2_period'

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    REAL(wp),                    INTENT(IN)    :: tracer(:,:,:,:) !< energy spectral bins (nproma,nlev,nblks_c,ntracer)
    REAL(wp),                    INTENT(IN)    :: emean(:,:)      !< total energy (nproma,nblks_c)
    REAL(wp),                    INTENT(INOUT) :: tm1(:,:)        !< tm1 period (nproma,nblks_c)
    REAL(wp),                    INTENT(INOUT) :: tm2(:,:)        !< tm2 period (nproma,nblks_c)
    REAL(wp),                    INTENT(INOUT) :: f1mean(:,:)     !< tm1 frequency (nproma,nblks_c)

    ! local
    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc,jb,jf,jk
    INTEGER :: jt                       !< tracer index
    INTEGER :: n                        !< loop index
    REAL(wp):: temp(nproma,wave_config%nfreqs)

    TYPE(t_wave_config), POINTER :: wc => NULL()

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)
    jk         = p_patch%nlev

    ! save some paperwork
    wc => wave_config

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jf,jt,n,i_startidx,i_endidx,temp) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
           &                 i_startidx, i_endidx, i_rlstart, i_rlend)
      ! compute sum of all tracers that match a specific frequency
      DO jf = 1,wc%nfreqs

        ! initialization
        DO jc = i_startidx, i_endidx
          temp(jc,jf) = 0._wp
        END DO

        DO n=1,SIZE(wc%list_tr(jf)%p)
          jt = wc%list_tr(jf)%p(n)
          DO jc = i_startidx, i_endidx
            temp(jc,jf) = temp(jc,jf) + tracer(jc,jk,jb,jt)
          END DO
        END DO  ! n

      END DO  ! jf

      ! tail part
      DO jc = i_startidx, i_endidx
        tm1(jc,jb) = wc%MP1_TAIL * temp(jc,wc%nfreqs)
        tm2(jc,jb) = wc%MP2_TAIL * temp(jc,wc%nfreqs)
      END DO

      ! add all other frequencies
      DO jf = 1,wc%nfreqs
        DO jc = i_startidx, i_endidx
          tm1(jc,jb) = tm1(jc,jb) + temp(jc,jf) * wc%dfim_fr(jf)
          tm2(jc,jb) = tm2(jc,jb) + temp(jc,jf) * wc%dfim_fr2(jf)
        END DO
      END DO

      ! clipping
      DO jc = i_startidx, i_endidx
        IF (emean(jc,jb).gt.EMIN) THEN
          tm1(jc,jb) = emean(jc,jb) / tm1(jc,jb)
          tm2(jc,jb) = SQRT(emean(jc,jb) / tm2(jc,jb))
        ELSE
          tm1(jc,jb) =  1.0_wp
          tm2(jc,jb) =  1.0_wp
        END IF
        f1mean(jc,jb) = 1.0_wp / tm1(jc,jb)
      END DO
    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE tm1_tm2_periods

  !>
  !! Calculation of new spectrum.
  !!
  !! Increase of spectrum in a time step and limit to a finite
  !! fraction of a typical F**(-4) equilibrium spectrum.
  !! Adaptation of WAM 4.5 code.
  !!
  SUBROUTINE new_spectrum(p_patch, wave_config, p_diag, dir10m, tracer)
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
         &  routine = 'new_spectrum'

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    TYPE(t_wave_diag),           INTENT(INOUT) :: p_diag
    REAL(wp),                    INTENT(IN)    :: dir10m(:,:)
    REAL(wp),                    INTENT(INOUT) :: tracer(:,:,:,:)

    TYPE(t_wave_config), POINTER :: wc => NULL()

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb,jc,jf,jd,jt,jk

    REAL(wp) :: temp_1, temp_2, temp_3
    REAL(wp) :: sprd
    REAL(wp) :: delfl

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
          jt = wc%tracer_ind(jd,jf)
          !
          DO jc = i_startidx, i_endidx
            delfl = 5.0E-07_wp * grav / wc%freqs(jf)**4 * dtime
            temp_2 = p_diag%ustar(jc,jb) * delfl &
                 &  * MAX(p_diag%femeanws(jc,jb),p_diag%femean(jc,jb))
            sprd = MAX(0._wp, COS(wc%dirs(jd)-dir10m(jc,jb)*deg2rad) )**2

            temp_1 = dtime * p_diag%sl(jc,jb,jt) &
                 / MAX(1._wp, 1._wp -  dtime * p_diag%fl(jc,jb,jt))

            temp_3 = MIN(ABS(temp_1),temp_2)

            tracer(jc,jk,jb,jt) = tracer(jc,jk,jb,jt)  + SIGN(temp_3,temp_1)
            tracer(jc,jk,jb,jt) = MAX(tracer(jc,jk,jb,jt), p_diag%FLMINFR(jc,jb,jf)*sprd)
          END DO
        END DO
      END DO
    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE new_spectrum

  !>
  !! Returns the last frequency index of prognostic part of spectrum
  !!
  !! Compute last frequency index of prognostic part of spectrum.
  !! Frequencies le MAX(tailfactor*max(fmnws,fm),tailfactor_pm*fpm),
  !! where fpm is the Pierson-Moskowitz frequency based on friction
  !! velocity. (fpm=g/(fric*zpi*ustar))
  !!
  !! Adaptation of WAM 4.5 code.
  !! FRCUTINDEX
  !!
  !! Initial revision by Mikhail Dobrynin, DWD (2019-10-10)
  !!
  SUBROUTINE last_prog_freq_ind(p_patch, wave_config, femeanws, femean, ustar, lpfi)
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = modname//'last_prog_freq_ind'

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    REAL(wp),                    INTENT(IN)    :: femeanws(:,:) !< mean frequency wind sea wave energy
    REAL(wp),                    INTENT(IN)    :: femean(:,:)   !< mean frequency wave energy
    REAL(wp),                    INTENT(IN)    :: ustar(:,:)    !< friction velocity
    INTEGER,                     INTENT(INOUT) :: lpfi(:,:)     !< last frequency index
                                                                !  of prognostic part of spectrum
    ! local
    TYPE(t_wave_config), POINTER :: wc => NULL()
    !
    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb, jc
    !
    REAL(wp) :: fpmh, fppm, fm2, fpm, fpm4, inv_log_co
    !
    REAL(wp), PARAMETER :: epsus = 1.0e-6_wp
    REAL(wp), PARAMETER :: fric = 28.0_wp
    REAL(wp), PARAMETER :: tailfactor = 2.5_wp
    REAL(wp), PARAMETER :: tailfactor_pm = 3.0_wp

    wc => wave_config

    fpmh = tailfactor / wc%freqs(1)
    fppm = tailfactor_pm * grav / (fric * pi2 * wc%freqs(1))
    inv_log_co = 1.0_wp / LOG10(wc%co)

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jb,i_startidx,i_endidx,fm2,fpm,fpm4) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jc = i_startidx, i_endidx
        fm2 = MAX(femeanws(jc,jb),femean(jc,jb)) * fpmh
        fpm = fppm / MAX(ustar(jc,jb),epsus)
        fpm4 = MAX(fm2,fpm)
        lpfi(jc,jb) = NINT(LOG10(fpm4)*inv_log_co)+1
        lpfi(jc,jb) = MIN(MAX(1,lpfi(jc,jb)),wc%nfreqs)
      END DO

    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE last_prog_freq_ind


  !>
  !! Impose high frequency tail to the spectrum
  !!
  !! Adaptation of WAM 4.5 code.
  !! IMPHFTAIL
  !!
  SUBROUTINE impose_high_freq_tail(p_patch, wave_config, wave_num_c, depth, last_prog_freq_ind, tracer)
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
         & routine =  modname//'impose_high_freq_tail'

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    REAL(wp),                    INTENT(IN)    :: wave_num_c(:,:,:)  !< wave number (1/m)
    REAL(wp),                    INTENT(IN)    :: depth(:,:)
    INTEGER,                     INTENT(IN)    :: last_prog_freq_ind(:,:)
    REAL(wp),                    INTENT(INOUT) :: tracer(:,:,:,:)


    TYPE(t_wave_config), POINTER :: wc => NULL()

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb,jf,jc,jd,jt,jtl,jk

    REAL(wp) :: gh, ak, akd, tcgond, akm1
    REAL(wp) :: temp(nproma, wave_config%nfreqs)
    REAL(wp) :: tfac(nproma)

    wc => wave_config
    gh = grav / (4.0_wp * pi)

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)
    jk         = p_patch%nlev

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jf,jc,jd,jt,jtl,i_startidx,i_endidx,ak,akd,tcgond,akm1,temp,tfac) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
           &                 i_startidx, i_endidx, i_rlstart, i_rlend)
      DO jf = 1, wc%nfreqs
        DO jc = i_startidx, i_endidx
          ak = wave_num_c(jc,jb,jf)
          IF (jf >=last_prog_freq_ind(jc,jb)) THEN
            akd = ak * depth(jc,jb)
            IF (akd.le.10.0_wp) THEN
              tcgond = 0.5_wp * SQRT(grav * TANH(akd)/ak) * (1.0_wp + 2.0_wp*akd/SINH(2.0_wp*akd))
            ELSE
              tcgond = gh / wc%freqs(jf)
            END IF
            akm1 = 1.0_wp/ak
            temp(jc,jf) = akm1**3/tcgond
          END IF
        END DO
      END DO

      DO jc = i_startidx, i_endidx
        DO jf = last_prog_freq_ind(jc,jb)+1, wc%nfreqs
          temp(jc,jf) = temp(jc,jf) / temp(jc,last_prog_freq_ind(jc,jb))
        END DO
      END DO

      DO jd = 1, wc%ndirs
        DO jc = i_startidx, i_endidx
          jtl = wc%tracer_ind(jd,last_prog_freq_ind(jc,jb))
          tfac(jc) = tracer(jc,jk,jb,jtl)
        ENDDO
        !
        DO jf = 1, wc%nfreqs
          jt = wc%tracer_ind(jd,jf)
          DO jc = i_startidx, i_endidx
            IF (jf >=last_prog_freq_ind(jc,jb)+1) THEN
              tracer(jc,jk,jb,jt) = temp(jc,jf) * tfac(jc)
            END IF
          END DO  !jc
        END DO  !jf
      END DO  !jd

    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE impose_high_freq_tail

  !>
  !! Calculation of input source function
  !!
  !! Compute input source function and the functional derivative of input
  !! source function.
  !!
  !! Adaptation of WAM 4.5 code.
  !! SINPUT
  !!
  !! Reference
  !! P. Janssen, JPO, 1989.
  !! P. Janssen, JPO., 1991.
  !!
  SUBROUTINE input_source_function(p_patch, wave_config, dir10m, tracer, p_diag)
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
         & routine =  modname//'input_source_function'

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    REAL(wp),                    INTENT(IN)    :: dir10m(:,:)
    REAL(wp),                    INTENT(INOUT) :: tracer(:,:,:,:)
    TYPE(t_wave_diag),           INTENT(INOUT) :: p_diag

    TYPE(t_wave_config), POINTER :: wc => NULL()

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb,jc,jf,jd,jt,jk

    REAL(wp) :: fac, const, const3, xk, cm, ucn, zcn, sh, cnsn
    REAL(wp) :: xv1d, temp, zbeta1, x, zlog, zlog2x, ufac

    wc => wave_config

    const3   = 2.0_wp * wc%xkappa / wc%betamax * wc%xkappa**2

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)
    jk         = p_patch%nlev

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jb,jf,jd,jt,i_startidx,i_endidx,fac,const,xk,cm,  &
!$OMP            ucn,zcn,sh,cnsn,xv1d,temp,zbeta1,x,zlog,zlog2x,ufac) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
           &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      FRE:DO jf = 1,wc%nfreqs
        fac = pi2 * wc%freqs(jf)
        const = fac * wc%xeps * wc%betamax / (wc%xkappa*wc%xkappa)

        DIR:DO jd = 1,wc%ndirs
          jt = wc%tracer_ind(jd,jf)

          DO jc = i_startidx, i_endidx
            xk = p_diag%wave_num_c(jc,jb,jf)
            cm = xk / fac
            ucn = p_diag%ustar(jc,jb) * cm + wc%zalp
            zcn = LOG(xk * p_diag%z0(jc,jb))
            sh = fac*fac / (grav * xk)
            cnsn = const * sh
            xv1d = -1.0_wp / (p_diag%ustar(jc,jb) / wc%xkappa * zcn * cm)
            temp = COS(wc%dirs(jd) - dir10m(jc,jb)*deg2rad)
            zbeta1 = const3 * (temp - xv1d) * ucn*ucn

            IF (temp > 0.01_wp) THEN
              x = temp * ucn
              zlog = zcn + wc%xkappa / x
              IF (zlog < 0._wp) THEN
                zlog2x = zlog*zlog * x
                ufac = EXP(zlog) * zlog2x*zlog2x + zbeta1
                p_diag%llws(jc,jb,jt) = 1
              ELSE
                ufac = zbeta1
                p_diag%llws(jc,jb,jt) = 0
              END IF
            ELSE
              ufac = zbeta1
              p_diag%llws(jc,jb,jt) = 0
            END IF

            p_diag%fl(jc,jb,jt) = cnsn * ufac
            p_diag%sl(jc,jb,jt) = tracer(jc,jk,jb,jt) * p_diag%fl(jc,jb,jt) !SL
          END DO
        END DO DIR
      END DO FRE
    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL
  END SUBROUTINE input_source_function

  !>
  !! Calculation of WM1 and WM2 wavenumbers
  !!
  !! Adaptation of WAM 4.5 code.
  !! WM1_WM2_WAVENUMBER_B
  !!
  !! C.Schneggenburger 08/97.
  !!
  !! Integration of spectra and adding of tail factors.
  !!   WM1 IS SQRT(1/K)*F
  !!   WM2 IS SQRT(K)*F
  !!
  SUBROUTINE wm1_wm2_wavenumber(p_patch, wave_config, wave_num_c, tracer, emean, akmean, xkmean)
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
         & routine =  modname//'wm1_wm2_wavenumber'

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    REAL(wp),                    INTENT(IN)    :: wave_num_c(:,:,:) !< wave number (1/m)
    REAL(wp),                    INTENT(IN)    :: tracer(:,:,:,:)
    REAL(wp),                    INTENT(IN)    :: emean(:,:)    !< total wave energy
    REAL(wp),                    INTENT(INOUT) :: akmean(:,:)   !< mean wavenumber based on SQRT(1/K)-moment, wm1
    REAL(wp),                    INTENT(INOUT) :: xkmean(:,:)   !< mean wavenumber based on SQRT(K)-moment, wm2

    TYPE(t_wave_config), POINTER :: wc => NULL()

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb,jc,jf,jt,jk
    INTEGER :: n

    REAL(wp) :: temp(nproma, wave_config%nfreqs)
    REAL(wp) :: temp2(nproma, wave_config%nfreqs)

    wc => wave_config

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)
    jk         = p_patch%nlev

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jb,jf,jt,i_startidx,i_endidx,temp,temp2) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
           &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jf = 1,wc%nfreqs

        DO jc = i_startidx, i_endidx
          temp(jc,jf) = 0._wp
          temp2(jc,jf)= SQRT(wave_num_c(jc,jb,jf))
        ENDDO
        !
        ! sum
        DO n=1,SIZE(wc%list_tr(jf)%p)
          jt = wc%list_tr(jf)%p(n)
          DO jc = i_startidx, i_endidx
            temp(jc,jf) = temp(jc,jf) + tracer(jc,jk,jb,jt)
          END DO
        ENDDO
      END DO  !jf

      !initialisation
      DO jc = i_startidx, i_endidx
        akmean(jc,jb) = wc%MM1_TAIL * SQRT(grav)/pi2 * temp(jc,wc%nfreqs)
        xkmean(jc,jb) = wc%MM1_TAIL * SQRT(grav)/pi2 * temp(jc,wc%nfreqs)
      END DO

      DO jf = 1,wc%nfreqs
        DO jc = i_startidx, i_endidx
          akmean(jc,jb) = akmean(jc,jb) + temp(jc,jf) / temp2(jc,jf) * wc%DFIM(jf)
          xkmean(jc,jb) = xkmean(jc,jb) + temp(jc,jf) * temp2(jc,jf) * wc%DFIM(jf)
        END DO
      END DO

      DO jc = i_startidx, i_endidx
        IF (emean(jc,jb) > EMIN) THEN
          akmean(jc,jb) = ( emean(jc,jb) / akmean(jc,jb))**2
          xkmean(jc,jb) = ( xkmean(jc,jb) / emean(jc,jb))**2
        ELSE
          akmean(jc,jb) = 1._wp
          xkmean(jc,jb) = 1._wp
        END IF
      END DO

    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE wm1_wm2_wavenumber


!!$  !>
!!$  !! Calculation of wave number.
!!$  !!
!!$  !! Wave number as a function of circular frequency and water depth.
!!$  !! Newtons method to solve the dispersion relation in shallow water.
!!$  !! G. KOMEN, P. JANSSEN   KNMI              01/06/1986
!!$  !! Adaptation of WAM 4.5 code, function AKI
!!$  !!
!!$  FUNCTION wave_number(OM, depth) RESULT(wave_num)
!!$    REAL(wp), INTENT(IN) :: OM                  !< CIRCULAR FREQUENCY 2*pi*freq (nfreqs)
!!$    REAL(wp), INTENT(IN) :: depth(:)            !< bathymetric height at cell centers (nproma)
!!$
!!$    REAL(wp), PARAMETER  :: EBS = 0.0001_wp     !< RELATIVE ERROR LIMIT OF NEWTON'S METHOD.
!!$    REAL(wp), PARAMETER  :: DKMAX = 40.0_wp     !< MAXIMUM VALUE OF DEPTH*WAVENUMBER.
!!$    REAL(wp)             :: BO, TH, STH
!!$    INTEGER              :: jc                  !< loop index for cells
!!$    INTEGER              :: ncell               !< loop bounds
!!$
!!$    LOGICAL              :: l_converged(SIZE(depth)), all_converged
!!$    REAL(wp)             :: wave_num(SIZE(depth))
!!$    REAL(wp)             :: AKP
!!$    ! ---------------------------------------------------------------------------- !
!!$    !     1. START WITH MAXIMUM FROM DEEP AND EXTREM SHALLOW WATER WAVE NUMBER.    !
!!$
!!$    ncell = SIZE(depth)
!!$    AKP   = 10000.0_wp
!!$
!!$    DO jc=1,ncell
!!$      !
!!$      ! initialization
!!$      wave_num(jc) = MAX( OM**2/(4.0_wp* grav), OM/(2.0_wp*SQRT(grav*depth(jc))) )
!!$
!!$      l_converged(jc) = (ABS(AKP-wave_num(jc)) .le. EBS*wave_num(jc))
!!$    ENDDO
!!$
!!$    all_converged = ALL(l_converged(1:ncell).EQV..TRUE.)
!!$
!!$    ! ---------------------------------------------------------------------------- !
!!$    !     2. ITERATION LOOP.                                                       !
!!$    DO WHILE (all_converged.EQV..FALSE.)
!!$      DO jc=1,ncell
!!$        IF (.NOT.l_converged(jc)) THEN
!!$          BO = depth(jc)*wave_num(jc)
!!$          IF (BO .gt. DKMAX) THEN
!!$            wave_num(jc) = OM**2/grav
!!$            l_converged(jc) = .TRUE.
!!$          ELSE
!!$            AKP = wave_num(jc)
!!$            TH = grav*wave_num(jc)*TANH(BO)
!!$            STH = SQRT(TH)
!!$            wave_num(jc) = wave_num(jc) &
!!$                 &                + (OM-STH)*STH*2.0_wp / (TH/wave_num(jc) + grav*BO/COSH(BO)**2)
!!$            ! check for converged solution
!!$            l_converged(jc) = (ABS(AKP-wave_num(jc)) .le. EBS*wave_num(jc))
!!$          END IF
!!$        END IF  ! l_converged
!!$      ENDDO  !jc
!!$      all_converged = ALL(l_converged(1:ncell).EQV..TRUE.)
!!$    ENDDO !while
!!$  END FUNCTION wave_number




  !>
  !! Calculation of wave number at cell centers.
  !!
  !! Wave number as a function of circular frequency and water depth.
  !! Newtons method to solve the dispersion relation in shallow water.
  !! G. KOMEN, P. JANSSEN   KNMI              01/06/1986
  !! Adaptation of WAM 4.5 code, function AKI
  !!
  SUBROUTINE wave_number_c(p_patch, wave_config, depth, wave_num_c)

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    REAL(wp),                    INTENT(IN)    :: depth(:,:)        !< bathymetric height at cell centers (nproma)
    REAL(wp),                    INTENT(INOUT) :: wave_num_c(:,:,:) !< Wave number as a function of circular frequency
    ! local
    REAL(wp), PARAMETER  :: EBS = 0.0001_wp     !< RELATIVE ERROR LIMIT OF NEWTON'S METHOD.
    REAL(wp), PARAMETER  :: DKMAX = 40.0_wp     !< MAXIMUM VALUE OF DEPTH*WAVENUMBER.
    REAL(wp)             :: BO, TH, STH
    INTEGER              :: jb, jc, jf          !< loop indices
    REAL(wp)             :: OM                  !< CIRCULAR FREQUENCY 2*pi*freq (nfreqs)
    LOGICAL              :: l_converged(nproma), all_converged
    REAL(wp)             :: wave_num(nproma)
    REAL(wp)             :: AKP
    INTEGER              :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER              :: i_startidx, i_endidx
    TYPE(t_wave_config), POINTER :: wc => NULL()
    !  ---------------------------------------------------------------------------- !
    !     1. START WITH MAXIMUM FROM DEEP AND EXTREM SHALLOW WATER WAVE NUMBER.    !

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

    wc  => wave_config
    AKP = 10000.0_wp

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jb,jf,i_startidx,i_endidx,OM,wave_num,l_converged,all_converged,BO,AKP,TH,STH) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,      &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jf = 1,wc%nfreqs
        OM = pi2 * wc%freqs(jf)
        DO jc=i_startidx, i_endidx
          !
          ! initialization
          wave_num(jc) = MAX( OM**2/(4.0_wp* grav), OM/(2.0_wp*SQRT(grav*depth(jc,jb))) )

          l_converged(jc) = (ABS(AKP-wave_num(jc)) .le. EBS*wave_num(jc))
        ENDDO  ! jc
        all_converged = ALL(l_converged(i_startidx:i_endidx).EQV..TRUE.)

        ! ---------------------------------------------------------------------------- !
        !     2. ITERATION LOOP.                                                       !

        DO WHILE (all_converged.EQV..FALSE.)
          DO jc = i_startidx, i_endidx
            IF (.NOT.l_converged(jc)) THEN
              BO = depth(jc,jb)*wave_num(jc)
              IF (BO > DKMAX) THEN
                wave_num(jc) = OM**2/grav
                l_converged(jc) = .TRUE.
              ELSE
                AKP = wave_num(jc)
                TH  = grav*wave_num(jc)*TANH(BO)
                STH = SQRT(TH)
                wave_num(jc) = wave_num(jc) &
                  &                + (OM-STH)*STH*2.0_wp / (TH/wave_num(jc) + grav*BO/COSH(BO)**2)
                ! check for converged solution
                l_converged(jc) = (ABS(AKP-wave_num(jc)) .le. EBS*wave_num(jc))
              END IF
            END IF  ! l_converged
          ENDDO  !jc
          all_converged = ALL(l_converged(i_startidx:i_endidx).EQV..TRUE.)
        ENDDO !while
        !
        DO jc = i_startidx, i_endidx
          wave_num_c(jc,jb,jf) = wave_num(jc)
        ENDDO
        !
      ENDDO ! jf
    ENDDO  ! jb
!$OMP END DO
!$OMP END PARALLEL
  END SUBROUTINE wave_number_c

  !>
  !! Calculation of wave number at edge midpoints.
  !!
  !! Wave number as a function of circular frequency and water depth.
  !! Newtons method to solve the dispersion relation in shallow water.
  !! G. KOMEN, P. JANSSEN   KNMI              01/06/1986
  !! Adaptation of WAM 4.5 code, function AKI
  !!
  SUBROUTINE wave_number_e(p_patch, wave_config, depth, wave_num_e)

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    REAL(wp),                    INTENT(IN)    :: depth(:,:)        !< bathymetric height at cell centers (nproma)
    REAL(wp),                    INTENT(INOUT) :: wave_num_e(:,:,:) !< Wave number as a function of circular frequency
    ! local
    REAL(wp), PARAMETER  :: EBS = 0.0001_wp     !< RELATIVE ERROR LIMIT OF NEWTON'S METHOD.
    REAL(wp), PARAMETER  :: DKMAX = 40.0_wp     !< MAXIMUM VALUE OF DEPTH*WAVENUMBER.
    REAL(wp)             :: BO, TH, STH
    INTEGER              :: jb, je, jf          !< loop indices
    REAL(wp)             :: OM                  !< CIRCULAR FREQUENCY 2*pi*freq (nfreqs)
    LOGICAL              :: l_converged(nproma), all_converged
    REAL(wp)             :: wave_num(nproma)
    REAL(wp)             :: AKP
    INTEGER              :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER              :: i_startidx, i_endidx
    TYPE(t_wave_config), POINTER :: wc => NULL()
    !  ---------------------------------------------------------------------------- !
    !     1. START WITH MAXIMUM FROM DEEP AND EXTREM SHALLOW WATER WAVE NUMBER.    !

    i_rlstart  = 1
    i_rlend    = min_rledge
    i_startblk = p_patch%edges%start_block(i_rlstart)
    i_endblk   = p_patch%edges%end_block(i_rlend)

    wc  => wave_config
    AKP = 10000.0_wp

!$OMP PARALLEL
!$OMP DO PRIVATE(je,jb,jf,i_startidx,i_endidx,OM,wave_num,l_converged,all_converged,BO,AKP,TH,STH) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_e( p_patch, jb, i_startblk, i_endblk,      &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jf = 1,wc%nfreqs
        OM = pi2 * wc%freqs(jf)
        DO je = i_startidx, i_endidx
          !
          ! initialization
          wave_num(je) = MAX( OM**2/(4.0_wp* grav), OM/(2.0_wp*SQRT(grav*depth(je,jb))) )

          l_converged(je) = (ABS(AKP-wave_num(je)) .le. EBS*wave_num(je))
        ENDDO  ! je
        all_converged = ALL(l_converged(i_startidx:i_endidx).EQV..TRUE.)

        ! ---------------------------------------------------------------------------- !
        !     2. ITERATION LOOP.                                                       !

        DO WHILE (all_converged.EQV..FALSE.)
          DO je = i_startidx, i_endidx
            IF (.NOT.l_converged(je)) THEN
              BO = depth(je,jb)*wave_num(je)
              IF (BO > DKMAX) THEN
                wave_num(je) = OM**2/grav
                l_converged(je) = .TRUE.
              ELSE
                AKP = wave_num(je)
                TH  = grav*wave_num(je)*TANH(BO)
                STH = SQRT(TH)
                wave_num(je) = wave_num(je) &
                  &                + (OM-STH)*STH*2.0_wp / (TH/wave_num(je) + grav*BO/COSH(BO)**2)
                ! check for converged solution
                l_converged(je) = (ABS(AKP-wave_num(je)) .le. EBS*wave_num(je))
              END IF
            END IF  ! l_converged
          ENDDO  !je
          all_converged = ALL(l_converged(i_startidx:i_endidx).EQV..TRUE.)
        ENDDO !while
        !
        DO je = i_startidx, i_endidx
          wave_num_e(je,jb,jf) = wave_num(je)
        ENDDO
        !
      ENDDO ! jf
    ENDDO  ! jb
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE wave_number_e

  !>
  !! Calculation of dissipation source function
  !!
  !! Compute dissipation source function and store additively into
  !! net source function array. Also compute functional derivative
  !! of dissipation source function.
  !!
  !! Adaptation of WAM 4.5 code SDISSIP
  !!     S.D.HASSELMANN.
  !!     MODIFIED TO SHALLOW WATER : G. KOMEN , P. JANSSEN
  !!     OPTIMIZATION : L. ZAMBRESKY
  !!     J. BIDLOT   ECMWF  FEBRUARY 1997   ADD SL IN SUBROUTINE CALL
  !!     H. GUENTHER GKSS  FEBRUARY 2002       FT 90
  !!     J. BIDLOT   ECMWF  NOVEMBER 2004  REFORMULATION BASED ON AKMEAN
  !!                                       AND FMEAN.
  !! Reference
  !! G.Komen, S. Hasselmann And K. Hasselmann, On The Existence
  !!          Of A Fully Developed Windsea Spectrum, JGR, 1984.
  !!
  SUBROUTINE dissipation_source_function(p_patch, wave_config, wave_num_c, tracer, p_diag)
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
         & routine =  modname//'dissipation_source_function'

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    REAL(wp),                    INTENT(IN)    :: wave_num_c(:,:,:) !< wave number (1/m)
    REAL(wp),                    INTENT(INOUT) :: tracer(:,:,:,:)
    TYPE(t_wave_diag),           INTENT(INOUT) :: p_diag

    TYPE(t_wave_config), POINTER :: wc => NULL()

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb,jk,jc,jf,jd,jt

    REAL(wp) :: sds, temp, sdiss

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)
    jk         = p_patch%nlev

    ! save some paperwork
    wc => wave_config

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jf,jd,jt,jc,sds,temp,sdiss,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
           &                 i_startidx, i_endidx, i_rlstart, i_rlend)
      DO jf = 1,wc%nfreqs
        DO jd = 1, wc%ndirs
          jt = wc%tracer_ind(jd,jf)
          DO jc = i_startidx, i_endidx

            sds = CONSS * p_diag%f1mean(jc,jb) * p_diag%emean(jc,jb)**2 * p_diag%xkmean(jc,jb)**4
            temp = wave_num_c(jc,jb,jf) / p_diag%xkmean(jc,jb)
            temp = sds * ((1.0_wp - DELTA) * temp +  DELTA * temp**2)
            sdiss = temp * tracer(jc,jk,jb,jt)

            p_diag%sl(jc,jb,jt) = p_diag%sl(jc,jb,jt) + sdiss
            p_diag%fl(jc,jb,jt) = p_diag%fl(jc,jb,jt) + temp
          END DO
        END DO
      END DO
    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL
  END SUBROUTINE dissipation_source_function

  !>
  !! Set wave spectrum to absolute allowed minimum
  !! Adaptation of WAM 4.5 code.
  !!
  SUBROUTINE set_energy2emin(p_patch, wave_config, tracer)
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
         &  routine = 'set_energy2emin'

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    REAL(wp),                    INTENT(INOUT) :: tracer(:,:,:,:)

    TYPE(t_wave_config), POINTER :: wc => NULL()

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb,jc,jf,jd,jt,jk
    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)
    jk         = p_patch%nlev

    ! save some paperwork
    wc => wave_config
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jf,jd,jt,jc,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
           &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jf = 1,wc%nfreqs
        DO jd = 1,wc%ndirs
          !
          jt = wc%tracer_ind(jd,jf)
          !
          DO jc = i_startidx, i_endidx
            tracer(jc,jk,jb,jt) = MAX(tracer(jc,jk,jb,jt),EMIN)
          END DO
        END DO
      END DO
    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL
  END SUBROUTINE set_energy2emin


  !>
  !! Calculation of bottom friction
  !!
  !! Compute dissipation of wave energy due to bottom friction
  !!
  !! Adaptation of WAM 4.5 code
  !!     G.J.KOMEN AND Q.D.GAO
  !!     OPTIMIZED BY L.F. ZAMBRESKY
  !!     H. GUENTHER   GKSS  FEBRUARY 2002       FT 90
  !!     E. MYKLEBUST        FEBRUARY 2005       OPTIMIZATION
  !!
  !! Reference
  !!  HASSELMANN ET AL, D. HYDR. Z SUPPL A12(1973) (JONSWAP)
  !!  BOUWS AND KOMEN, JPO 13(1983)1653-1658
  !!
  SUBROUTINE bottom_friction(p_patch, wave_config, wave_num_c, depth, tracer, p_diag)
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
         & routine =  modname//'bottom_friction'

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    REAL(wp),                    INTENT(IN)    :: wave_num_c(:,:,:) !< wave number (1/m)
    REAL(wp),                    INTENT(IN)    :: depth(:,:)
    REAL(wp),                    INTENT(INOUT) :: tracer(:,:,:,:)
    TYPE(t_wave_diag),           INTENT(INOUT) :: p_diag

    TYPE(t_wave_config), POINTER :: wc => NULL()

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb,jk,jc,jf,jd,jt

    REAL(wp) :: const, sbo

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)
    jk         = p_patch%nlev

    ! save some paperwork
    wc => wave_config

    const = -2.0_wp*0.038_wp/grav

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jf,jd,jt,jc,sbo,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
           &                 i_startidx, i_endidx, i_rlstart, i_rlend)
      DO jf = 1,wc%nfreqs
        DO jd = 1, wc%ndirs
          jt = wc%tracer_ind(jd,jf)
          DO jc = i_startidx, i_endidx
            sbo = MIN(2.0_wp * depth(jc,jb) * wave_num_c(jc,jb,jf),50.0_wp)
            sbo = const * wave_num_c(jc,jb,jf) / SINH(sbo)
            p_diag%sl(jc,jb,jt) = p_diag%sl(jc,jb,jt) + sbo*tracer(jc,jk,jb,jt)
            p_diag%fl(jc,jb,jt) = p_diag%fl(jc,jb,jt) + sbo
          END DO
        END DO
      END DO
    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL
  END SUBROUTINE bottom_friction




  !>
  !!
  !! Computation of nonlinear transfer rate and its
  !! functional derivative (diagonal terms only) and
  !! addition to corresponding net expressions.
  !!
  !! Adaptation of WAM 4.5 SNONLIN
  !! S.D. Hasselmann.  MPI
  !! G. Komen, P. Janssen   KNMI        modified to shallow water
  !! H. Guenther, L. Zambresky          optimized
  !! H. Guenther  GKSS/ECMWF  June 1991 interactions between diag-
  !!                                    and prognostic part.
  !! H. Guenther  GKSS  February 2002   FT 90
  !! E. Myklebust       February 2005   optimization
  !!
  SUBROUTINE nonlinear_transfer(p_patch, wave_config, depth, tracer, p_diag)
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
         & routine =  modname//'nonlinear_transfer'

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    REAL(wp),                    INTENT(IN)    :: depth(:,:)
    REAL(wp),                    INTENT(INOUT) :: tracer(:,:,:,:)
    TYPE(t_wave_diag),           INTENT(INOUT) :: p_diag

    ! local
    TYPE(t_wave_config), POINTER :: wc             => NULL()
    INTEGER,             POINTER:: tr_idx(:,:,:,:) => NULL()

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc,jk,jb,jf

    INTEGER :: nfreqs, ndirs
    INTEGER :: MP, MP1, MM, MM1, IC, IP, IP1, IM, IM1, KH, K

    REAL(wp) :: FFACP, FFACP1, FFACM1, FTAIL, FKLAMP, FKLAMP1, GW1, GW2, GW3, GW4
    REAL(wp) :: FKLAMPA, FKLAMPB, FKLAMP2, FKLAPA2, FKLAPB2, FKLAP12, FKLAP22
    REAL(wp) :: FKLAMM, FKLAMM1, GW5, GW6, GW7, GW8, FKLAMMA, FKLAMMB, FKLAMM2
    REAL(wp) :: FKLAMA2, FKLAMB2, FKLAM12, FKLAM22
    REAL(wp) :: SAP, SAM, FIJ, FAD1, FAD2, FCEN

    REAL(wp) :: AD
    REAL(wp) :: DELAD, DELAP, DELAM
    REAL(wp) :: FTEMP, ENHFR


    ! convenience pointers
    wc     => wave_config
    tr_idx => p_diag%non_lin_tr_ind(:,:,:,:)
    nfreqs = wc%nfreqs
    ndirs  = wc%ndirs

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)
    jk         = p_patch%nlev

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jf,i_startidx,i_endidx,ENHFR,FTEMP,                              &
!$OMP           MP,MP1,MM,MM1,IC,IP,IP1,IM,IM1,KH,K,FFACP,                               &
!$OMP           FFACP1,FFACM1,FTAIL,FKLAMP,FKLAMP1, GW1, GW2, GW3, GW4,FKLAMPA, FKLAMPB, &
!$OMP           FKLAMP2, FKLAPA2, FKLAPB2, FKLAP12, FKLAP22, FKLAMM, FKLAMM1, GW5, GW6, &
!$OMP           GW7, GW8, FKLAMMA, FKLAMMB, FKLAMM2,FKLAMA2, FKLAMB2, FKLAM12, FKLAM22, &
!$OMP           SAP, SAM, FIJ, FAD1, FAD2, FCEN, AD, DELAD, DELAP, DELAM ) ICON_OMP_DEFAULT_SCHEDULE
    ljb: DO jb = i_startblk, i_endblk
       CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
            &                 i_startidx, i_endidx, i_rlstart, i_rlend)
       DO jc = i_startidx, i_endidx
          ENHFR = MAX(0.75_wp*depth(jc,jb)*p_diag%AKMEAN(jc,jb), 0.5_wp)
          ENHFR = 1.0_wp + (5.5_wp/ENHFR) * (1.0_wp-0.833_wp*ENHFR) &
               &                * EXP(-1.25_wp*ENHFR)
          p_diag%ENH(jc,jb) = ENHFR
       END DO

      FRE4: DO jf = 1,nfreqs+4
        MP  = p_diag%IKP (jf)
        MP1 = p_diag%IKP1(jf)
        MM  = p_diag%IKM (jf)
        MM1 = p_diag%IKM1(jf)
        FFACP  = 1._wp
        FFACP1 = 1._wp
        FFACM1 = 1._wp
        FTAIL  = 1._wp
        IC  = jf
        IP  = MP
        IP1 = MP1
        IM  = MM
        IM1 = MM1
        IF (IP1.GT.nfreqs) THEN
          FFACP1 = wc%FRH(IP1-nfreqs+1)
          IP1 = nfreqs
          IF (IP.GT.nfreqs) THEN
            FFACP  = wc%FRH(IP-nfreqs+1)
            IP  = nfreqs
            IF (IC.GT.nfreqs) THEN
              FTAIL  = wc%FRH(IC-nfreqs+1)
              IC  = nfreqs
              IF (IM1.GT.nfreqs) THEN
                FFACM1 = wc%FRH(IM1-nfreqs+1)
                IM1 = nfreqs
              END IF
            END IF
          END IF
        END IF

        FKLAMP  = p_diag%FKLAP(jf)
        FKLAMP1 = p_diag%FKLAP1(jf)
        GW2 = FKLAMP1*FFACP*wc%DAL1
        GW1 = GW2*wc%CL11
        GW2 = GW2*wc%ACL1
        GW4 = FKLAMP*FFACP1*wc%DAL1
        GW3 = GW4*wc%CL11
        GW4 = GW4*wc%ACL1
        FKLAMPA = FKLAMP*wc%CL11
        FKLAMPB = FKLAMP*wc%ACL1
        FKLAMP2 = FKLAMP1*wc%ACL1
        FKLAMP1 = FKLAMP1*wc%CL11
        FKLAPA2 = FKLAMPA**2
        FKLAPB2 = FKLAMPB**2
        FKLAP12 = FKLAMP1**2
        FKLAP22 = FKLAMP2**2
        FKLAMM  = p_diag%FKLAM(jf)
        FKLAMM1 = p_diag%FKLAM1(jf)
        GW6 = FKLAMM1*wc%DAL2
        GW5 = GW6*wc%CL21
        GW6 = GW6*wc%ACL2
        GW8 = FKLAMM*FFACM1*wc%DAL2
        GW7 = GW8*wc%CL21
        GW8 = GW8*wc%ACL2
        FKLAMMA = FKLAMM*wc%CL21
        FKLAMMB = FKLAMM*wc%ACL2
        FKLAMM2 = FKLAMM1*wc%ACL2
        FKLAMM1 = FKLAMM1*wc%CL21
        FKLAMA2 = FKLAMMA**2
        FKLAMB2 = FKLAMMB**2
        FKLAM12 = FKLAMM1**2
        FKLAM22 = FKLAMM2**2

        ! shallow water case
        IF (jf.GT.4) THEN! UNTIL 7.
          IF (MM1.LE.nfreqs) THEN! UNTIL 6.
            IF (jf .LE.nfreqs) THEN! UNTIL 5.
              IF (MP .LE.nfreqs) THEN! UNTIL 4.
                IF (MP1.LE.nfreqs) THEN! UNTIL 3.
                  !     2.1.1   ANGULAR LOOP.                                     !
                  DIR2: DO K = 1,ndirs !DIR2

                    DO jc = i_startidx, i_endidx
                      !     2.1.1.1 LOOP OVER GRIDPOINTS.. NONLINEAR TRANSFER AND     !
                      !             DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.         !
                      !             ----------------------------------------------    !

                      !     2.1 LOOP FOR ANLULAR SYMMETRY.                            !
                      MIR2: DO KH = 1,2

                         SAP = &
                              GW1*tracer(jc,jk,jb,tr_idx(1,jf,KH,K)) + &
                              GW2*tracer(jc,jk,jb,tr_idx(2,jf,KH,K)) + &
                              GW3*tracer(jc,jk,jb,tr_idx(3,jf,KH,K)) + &
                              GW4*tracer(jc,jk,jb,tr_idx(4,jf,KH,K))
                         SAM = &
                              GW5*tracer(jc,jk,jb,tr_idx(5,jf,KH,K)) + &
                              GW6*tracer(jc,jk,jb,tr_idx(6,jf,KH,K)) + &
                              GW7*tracer(jc,jk,jb,tr_idx(7,jf,KH,K)) + &
                              GW8*tracer(jc,jk,jb,tr_idx(8,jf,KH,K))

                         FTEMP = p_diag%AF11(jf) * p_diag%ENH(jc,jb)
                         FIJ = tracer(jc,jk,jb,tr_idx(9,jf,KH,K))*FTAIL
                         FAD1 = FIJ*(SAP+SAM)
                         FAD2 = FAD1-2._wp*SAP*SAM
                         FAD1 = FAD1+FAD2
                         FCEN = FTEMP*FIJ
                         AD = FAD2*FCEN
                         DELAD = FAD1*FTEMP
                         DELAP = (FIJ-2._wp*SAM)*wc%DAL1*FCEN
                         DELAM = (FIJ-2._wp*SAP)*wc%DAL2*FCEN

                         p_diag%sl(jc,jb,tr_idx(10,jf,KH,K)) = p_diag%sl(jc,jb,tr_idx(10,jf,KH,K)) + AD*FKLAMM1 !SL
                         p_diag%sl(jc,jb,tr_idx(11,jf,KH,K)) = p_diag%sl(jc,jb,tr_idx(11,jf,KH,K)) + AD*FKLAMM2 !SL
                         p_diag%fl(jc,jb,tr_idx(10,jf,KH,K)) = p_diag%fl(jc,jb,tr_idx(10,jf,KH,K)) + DELAM*FKLAM12 !FL
                         p_diag%fl(jc,jb,tr_idx(11,jf,KH,K)) = p_diag%fl(jc,jb,tr_idx(11,jf,KH,K)) + DELAM*FKLAM22 !FL

                         p_diag%sl(jc,jb,tr_idx(12,jf,KH,K)) = p_diag%sl(jc,jb,tr_idx(12,jf,KH,K)) + AD*FKLAMMA
                         p_diag%sl(jc,jb,tr_idx(13,jf,KH,K)) = p_diag%sl(jc,jb,tr_idx(13,jf,KH,K)) + AD*FKLAMMB
                         p_diag%fl(jc,jb,tr_idx(12,jf,KH,K)) = p_diag%fl(jc,jb,tr_idx(12,jf,KH,K)) + DELAM*FKLAMA2
                         p_diag%fl(jc,jb,tr_idx(13,jf,KH,K)) = p_diag%fl(jc,jb,tr_idx(13,jf,KH,K)) + DELAM*FKLAMB2

                         p_diag%sl(jc,jb,tr_idx(14,jf,KH,K)) = p_diag%sl(jc,jb,tr_idx(14,jf,KH,K)) - 2._wp*AD
                         p_diag%fl(jc,jb,tr_idx(14,jf,KH,K)) = p_diag%fl(jc,jb,tr_idx(14,jf,KH,K)) - 2._wp*DELAD

                         p_diag%sl(jc,jb,tr_idx(15,jf,KH,K)) = p_diag%sl(jc,jb,tr_idx(15,jf,KH,K)) + AD*FKLAMP1
                         p_diag%sl(jc,jb,tr_idx(16,jf,KH,K)) = p_diag%sl(jc,jb,tr_idx(16,jf,KH,K)) + AD*FKLAMP2
                         p_diag%fl(jc,jb,tr_idx(15,jf,KH,K)) = p_diag%fl(jc,jb,tr_idx(15,jf,KH,K)) + DELAP*FKLAP12
                         p_diag%fl(jc,jb,tr_idx(16,jf,KH,K)) = p_diag%fl(jc,jb,tr_idx(16,jf,KH,K)) + DELAP*FKLAP22

                         p_diag%sl(jc,jb,tr_idx(17,jf,KH,K)) = p_diag%sl(jc,jb,tr_idx(17,jf,KH,K)) + AD*FKLAMPA
                         p_diag%sl(jc,jb,tr_idx(18,jf,KH,K)) = p_diag%sl(jc,jb,tr_idx(18,jf,KH,K)) + AD*FKLAMPB
                         p_diag%fl(jc,jb,tr_idx(17,jf,KH,K)) = p_diag%fl(jc,jb,tr_idx(17,jf,KH,K)) + DELAP*FKLAPA2
                         p_diag%fl(jc,jb,tr_idx(18,jf,KH,K)) = p_diag%fl(jc,jb,tr_idx(18,jf,KH,K)) + DELAP*FKLAPB2
                      END DO MIR2
                    END DO  ! jc
                  END DO DIR2
                ELSE!IF (MP1.LE.ML) THEN
                  !     3.1.1   ANGULAR LOOP.                                   !
                  DIR3: DO K = 1, ndirs

                    DO jc = i_startidx, i_endidx
                      !     3.1.1.1 LOOP OVER GRIDPOINTS.. NONLINEAR TRANSFER AND   !
                      !             DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.       !
                      !             ----------------------------------------------  !

                      !     3.1 LOOP FOR ANGULAR SYMMETRY.                          !
                      MIR3: DO KH = 1,2

                        SAP = &
                             GW1*tracer(jc,jk,jb,tr_idx(1,jf,KH,K)) + &
                             GW2*tracer(jc,jk,jb,tr_idx(2,jf,KH,K)) + &
                             GW3*tracer(jc,jk,jb,tr_idx(3,jf,KH,K)) + &
                             GW4*tracer(jc,jk,jb,tr_idx(4,jf,KH,K))
                        SAM = &
                             GW5*tracer(jc,jk,jb,tr_idx(5,jf,KH,K)) + &
                             GW6*tracer(jc,jk,jb,tr_idx(6,jf,KH,K)) + &
                             GW7*tracer(jc,jk,jb,tr_idx(7,jf,KH,K)) + &
                             GW8*tracer(jc,jk,jb,tr_idx(8,jf,KH,K))

                        FTEMP = p_diag%AF11(jf) * p_diag%ENH(jc,jb)
                        FIJ = tracer(jc,jk,jb,tr_idx(9,jf,KH,K))*FTAIL
                        FAD1 = FIJ*(SAP+SAM)
                        FAD2 = FAD1-2._wp*SAP*SAM
                        FAD1 = FAD1+FAD2
                        FCEN = FTEMP*FIJ
                        AD = FAD2*FCEN
                        DELAD = FAD1*FTEMP
                        DELAP = (FIJ-2._wp*SAM)*wc%DAL1*FCEN
                        DELAM = (FIJ-2._wp*SAP)*wc%DAL2*FCEN

                        p_diag%sl(jc,jb,tr_idx(10,jf,KH,K)) = p_diag%sl(jc,jb,tr_idx(10,jf,KH,K)) + AD*FKLAMM1
                        p_diag%sl(jc,jb,tr_idx(11,jf,KH,K)) = p_diag%sl(jc,jb,tr_idx(11,jf,KH,K)) + AD*FKLAMM2
                        p_diag%fl(jc,jb,tr_idx(10,jf,KH,K)) = p_diag%fl(jc,jb,tr_idx(10,jf,KH,K)) + DELAM*FKLAM12
                        p_diag%fl(jc,jb,tr_idx(11,jf,KH,K)) = p_diag%fl(jc,jb,tr_idx(11,jf,KH,K)) + DELAM*FKLAM22

                        p_diag%sl(jc,jb,tr_idx(12,jf,KH,K)) = p_diag%sl(jc,jb,tr_idx(12,jf,KH,K)) + AD*FKLAMMA
                        p_diag%sl(jc,jb,tr_idx(13,jf,KH,K)) = p_diag%sl(jc,jb,tr_idx(13,jf,KH,K)) + AD*FKLAMMB
                        p_diag%fl(jc,jb,tr_idx(12,jf,KH,K)) = p_diag%fl(jc,jb,tr_idx(12,jf,KH,K)) + DELAM*FKLAMA2
                        p_diag%fl(jc,jb,tr_idx(13,jf,KH,K)) = p_diag%fl(jc,jb,tr_idx(13,jf,KH,K)) + DELAM*FKLAMB2

                        p_diag%sl(jc,jb,tr_idx(14,jf,KH,K)) = p_diag%sl(jc,jb,tr_idx(14,jf,KH,K)) - 2._wp*AD
                        p_diag%fl(jc,jb,tr_idx(14,jf,KH,K)) = p_diag%fl(jc,jb,tr_idx(14,jf,KH,K)) - 2._wp*DELAD

                        p_diag%sl(jc,jb,tr_idx(15,jf,KH,K)) = p_diag%sl(jc,jb,tr_idx(15,jf,KH,K)) + AD*FKLAMP1
                        p_diag%sl(jc,jb,tr_idx(16,jf,KH,K)) = p_diag%sl(jc,jb,tr_idx(16,jf,KH,K)) + AD*FKLAMP2
                        p_diag%fl(jc,jb,tr_idx(15,jf,KH,K)) = p_diag%fl(jc,jb,tr_idx(15,jf,KH,K)) + DELAP*FKLAP12
                        p_diag%fl(jc,jb,tr_idx(16,jf,KH,K)) = p_diag%fl(jc,jb,tr_idx(16,jf,KH,K)) + DELAP*FKLAP22
                      END DO MIR3
                    END DO !jc
                  END DO DIR3!  BRANCH BACK TO 3.1.1 FOR NEXT DIRECTION.
                END IF
              ELSE!IF (MP .LE.ML) THEN
                !     4.1.1   ANGULAR LOOP.                                   !
                DIR4: DO K = 1, ndirs!

                  DO jc = i_startidx, i_endidx
                    !     4.1.1.1 LOOP OVER GRIDPOINTS.. NONLINEAR TRANSFER AND   !
                    !             DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.       !
                    !             ----------------------------------------------  !

                    !     4.1 LOOP FOR ANGULAR SYMMETRY.
                    MIR4: DO KH = 1,2

                      SAP = &
                           GW1*tracer(jc,jk,jb,tr_idx(1,jf,KH,K)) + &
                           GW2*tracer(jc,jk,jb,tr_idx(2,jf,KH,K)) + &
                           GW3*tracer(jc,jk,jb,tr_idx(3,jf,KH,K)) + &
                           GW4*tracer(jc,jk,jb,tr_idx(4,jf,KH,K))
                      SAM = &
                           GW5*tracer(jc,jk,jb,tr_idx(5,jf,KH,K)) + &
                           GW6*tracer(jc,jk,jb,tr_idx(6,jf,KH,K)) + &
                           GW7*tracer(jc,jk,jb,tr_idx(7,jf,KH,K)) + &
                           GW8*tracer(jc,jk,jb,tr_idx(8,jf,KH,K))

                      FTEMP = p_diag%AF11(jf) * p_diag%ENH(jc,jb)
                      FIJ = tracer(jc,jk,jb,tr_idx(9,jf,KH,K))*FTAIL
                      FAD1 = FIJ*(SAP+SAM)
                      FAD2 = FAD1-2._wp*SAP*SAM
                      FAD1 = FAD1+FAD2
                      FCEN = FTEMP*FIJ
                      AD = FAD2*FCEN
                      DELAD = FAD1*FTEMP
                      DELAP = (FIJ-2._wp*SAM)*wc%DAL1*FCEN
                      DELAM = (FIJ-2._wp*SAP)*wc%DAL2*FCEN

                      p_diag%sl(jc,jb,tr_idx(10,jf,KH,K)) = p_diag%sl(jc,jb,tr_idx(10,jf,KH,K)) + AD*FKLAMM1
                      p_diag%sl(jc,jb,tr_idx(11,jf,KH,K)) = p_diag%sl(jc,jb,tr_idx(11,jf,KH,K)) + AD*FKLAMM2
                      p_diag%fl(jc,jb,tr_idx(10,jf,KH,K)) = p_diag%fl(jc,jb,tr_idx(10,jf,KH,K)) + DELAM*FKLAM12
                      p_diag%fl(jc,jb,tr_idx(11,jf,KH,K)) = p_diag%fl(jc,jb,tr_idx(11,jf,KH,K)) + DELAM*FKLAM22

                      p_diag%sl(jc,jb,tr_idx(12,jf,KH,K)) = p_diag%sl(jc,jb,tr_idx(12,jf,KH,K)) + AD*FKLAMMA
                      p_diag%sl(jc,jb,tr_idx(13,jf,KH,K)) = p_diag%sl(jc,jb,tr_idx(13,jf,KH,K)) + AD*FKLAMMB
                      p_diag%fl(jc,jb,tr_idx(12,jf,KH,K)) = p_diag%fl(jc,jb,tr_idx(12,jf,KH,K)) + DELAM*FKLAMA2
                      p_diag%fl(jc,jb,tr_idx(13,jf,KH,K)) = p_diag%fl(jc,jb,tr_idx(13,jf,KH,K)) + DELAM*FKLAMB2

                      p_diag%sl(jc,jb,tr_idx(14,jf,KH,K)) = p_diag%sl(jc,jb,tr_idx(14,jf,KH,K)) - 2._wp*AD
                      p_diag%fl(jc,jb,tr_idx(14,jf,KH,K)) = p_diag%fl(jc,jb,tr_idx(14,jf,KH,K)) - 2._wp*DELAD
                    END DO MIR4
                  END DO !jc
                END DO DIR4
              END IF
            ELSE
              !     5.1.1   ANGULAR LOOP.                                         !
              DIR5: DO K = 1, ndirs

                DO jc = i_startidx, i_endidx

                  !     5.1 LOOP FOR ANLULAR SYMMETRY.                                !
                  MIR5: DO KH = 1,2

                    SAP = &
                         GW1*tracer(jc,jk,jb,tr_idx(1,jf,KH,K)) + &
                         GW2*tracer(jc,jk,jb,tr_idx(2,jf,KH,K)) + &
                         GW3*tracer(jc,jk,jb,tr_idx(3,jf,KH,K)) + &
                         GW4*tracer(jc,jk,jb,tr_idx(4,jf,KH,K))
                    SAM = &
                         GW5*tracer(jc,jk,jb,tr_idx(5,jf,KH,K)) + &
                         GW6*tracer(jc,jk,jb,tr_idx(6,jf,KH,K)) + &
                         GW7*tracer(jc,jk,jb,tr_idx(7,jf,KH,K)) + &
                         GW8*tracer(jc,jk,jb,tr_idx(8,jf,KH,K))

                    FTEMP = p_diag%AF11(jf) * p_diag%ENH(jc,jb)
                    FIJ = tracer(jc,jk,jb,tr_idx(9,jf,KH,K))*FTAIL
                    FAD1 = FIJ*(SAP+SAM)
                    FAD2 = FAD1-2._wp*SAP*SAM
                    FAD1 = FAD1+FAD2
                    FCEN = FTEMP*FIJ
                    AD = FAD2*FCEN
                    DELAD = FAD1*FTEMP
                    DELAP = (FIJ-2._wp*SAM)*wc%DAL1*FCEN
                    DELAM = (FIJ-2._wp*SAP)*wc%DAL2*FCEN

                    p_diag%sl(jc,jb,tr_idx(10,jf,KH,K)) = p_diag%sl(jc,jb,tr_idx(10,jf,KH,K)) + AD*FKLAMM1
                    p_diag%sl(jc,jb,tr_idx(11,jf,KH,K)) = p_diag%sl(jc,jb,tr_idx(11,jf,KH,K)) + AD*FKLAMM2
                    p_diag%fl(jc,jb,tr_idx(10,jf,KH,K)) = p_diag%fl(jc,jb,tr_idx(10,jf,KH,K)) + DELAM*FKLAM12
                    p_diag%fl(jc,jb,tr_idx(11,jf,KH,K)) = p_diag%fl(jc,jb,tr_idx(11,jf,KH,K)) + DELAM*FKLAM22

                    p_diag%sl(jc,jb,tr_idx(12,jf,KH,K)) = p_diag%sl(jc,jb,tr_idx(12,jf,KH,K)) + AD*FKLAMMA
                    p_diag%sl(jc,jb,tr_idx(13,jf,KH,K)) = p_diag%sl(jc,jb,tr_idx(13,jf,KH,K)) + AD*FKLAMMB
                    p_diag%fl(jc,jb,tr_idx(12,jf,KH,K)) = p_diag%fl(jc,jb,tr_idx(12,jf,KH,K)) + DELAM*FKLAMA2
                    p_diag%fl(jc,jb,tr_idx(13,jf,KH,K)) = p_diag%fl(jc,jb,tr_idx(13,jf,KH,K)) + DELAM*FKLAMB2
                  END DO MIR5
                END DO !jc
              END DO DIR5
            END IF
          ELSE
            !     6.1.1   ANGULAR LOOP.                                        !
            DIR6: DO K = 1,ndirs

              DO jc = i_startidx, i_endidx

                !     6.1 LOOP FOR ANGULAR SYMMETRY.                                !
                MIR6: DO KH = 1,2

                  SAP = &
                       GW1*tracer(jc,jk,jb,tr_idx(1,jf,KH,K)) + &
                       GW2*tracer(jc,jk,jb,tr_idx(2,jf,KH,K)) + &
                       GW3*tracer(jc,jk,jb,tr_idx(3,jf,KH,K)) + &
                       GW4*tracer(jc,jk,jb,tr_idx(4,jf,KH,K))
                  SAM = &
                       GW5*tracer(jc,jk,jb,tr_idx(5,jf,KH,K)) + &
                       GW6*tracer(jc,jk,jb,tr_idx(6,jf,KH,K)) + &
                       GW7*tracer(jc,jk,jb,tr_idx(7,jf,KH,K)) + &
                       GW8*tracer(jc,jk,jb,tr_idx(8,jf,KH,K))

                  FTEMP = p_diag%AF11(jf) * p_diag%ENH(jc,jb)
                  FIJ = tracer(jc,jk,jb,tr_idx(9,jf,KH,K))*FTAIL
                  FAD1 = FIJ*(SAP+SAM)
                  FAD2 = FAD1-2._wp*SAP*SAM
                  FAD1 = FAD1+FAD2
                  FCEN = FTEMP*FIJ
                  AD = FAD2*FCEN
                  DELAD = FAD1*FTEMP
                  DELAP = (FIJ-2._wp*SAM)*wc%DAL1*FCEN
                  DELAM = (FIJ-2._wp*SAP)*wc%DAL2*FCEN

                  p_diag%sl(jc,jb,tr_idx(10,jf,KH,K)) = p_diag%sl(jc,jb,tr_idx(10,jf,KH,K)) + AD*FKLAMM1
                  p_diag%sl(jc,jb,tr_idx(11,jf,KH,K)) = p_diag%sl(jc,jb,tr_idx(11,jf,KH,K)) + AD*FKLAMM2
                  p_diag%fl(jc,jb,tr_idx(10,jf,KH,K)) = p_diag%fl(jc,jb,tr_idx(10,jf,KH,K)) + DELAM*FKLAM12
                  p_diag%fl(jc,jb,tr_idx(11,jf,KH,K)) = p_diag%fl(jc,jb,tr_idx(11,jf,KH,K)) + DELAM*FKLAM22
                END DO MIR6
              END DO !jc
            END DO DIR6
          END IF
        ELSE
          !     7.1.1   ANGULAR LOOP.                                           !
          DIR7: DO K = 1,ndirs

            DO jc = i_startidx, i_endidx

              !     7.1 LOOP FOR ANGULAR SYMMETRY.                                     !
              MIR7: DO KH = 1,2

                SAP = &
                     GW1*tracer(jc,jk,jb,tr_idx(1,jf,KH,K)) + &
                     GW2*tracer(jc,jk,jb,tr_idx(2,jf,KH,K)) + &
                     GW3*tracer(jc,jk,jb,tr_idx(3,jf,KH,K)) + &
                     GW4*tracer(jc,jk,jb,tr_idx(4,jf,KH,K))

                FTEMP = p_diag%AF11(jf) * p_diag%ENH(jc,jb)
                FIJ = tracer(jc,jk,jb,tr_idx(9,jf,KH,K))
                FAD2 = FIJ*SAP
                FAD1 = 2._wp*FAD2
                FCEN = FTEMP*FIJ
                AD = FAD2*FCEN
                DELAD = FAD1*FTEMP
                DELAP = FIJ*wc%DAL1*FCEN

                p_diag%sl(jc,jb,tr_idx(14,jf,KH,K)) = p_diag%sl(jc,jb,tr_idx(14,jf,KH,K)) - 2._wp*AD
                p_diag%sl(jc,jb,tr_idx(15,jf,KH,K)) = p_diag%sl(jc,jb,tr_idx(15,jf,KH,K)) + AD*FKLAMP1
                p_diag%sl(jc,jb,tr_idx(16,jf,KH,K)) = p_diag%sl(jc,jb,tr_idx(16,jf,KH,K)) + AD*FKLAMP2
                p_diag%sl(jc,jb,tr_idx(17,jf,KH,K)) = p_diag%sl(jc,jb,tr_idx(17,jf,KH,K)) + AD*FKLAMPA
                p_diag%sl(jc,jb,tr_idx(18,jf,KH,K)) = p_diag%sl(jc,jb,tr_idx(18,jf,KH,K)) + AD*FKLAMPB

                p_diag%fl(jc,jb,tr_idx(14,jf,KH,K)) = p_diag%fl(jc,jb,tr_idx(14,jf,KH,K)) - 2._wp*DELAD
                p_diag%fl(jc,jb,tr_idx(15,jf,KH,K)) = p_diag%fl(jc,jb,tr_idx(15,jf,KH,K)) + DELAP*FKLAP12
                p_diag%fl(jc,jb,tr_idx(16,jf,KH,K)) = p_diag%fl(jc,jb,tr_idx(16,jf,KH,K)) + DELAP*FKLAP22
                p_diag%fl(jc,jb,tr_idx(17,jf,KH,K)) = p_diag%fl(jc,jb,tr_idx(17,jf,KH,K)) + DELAP*FKLAPA2
                p_diag%fl(jc,jb,tr_idx(18,jf,KH,K)) = p_diag%fl(jc,jb,tr_idx(18,jf,KH,K)) + DELAP*FKLAPB2
              END DO MIR7
            END DO !jc
          END DO DIR7
        END IF
      END DO FRE4
    END DO ljb
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE nonlinear_transfer


  !>
  !! Calculate grid, depth and current refraction
  !!
  !! Based on WAM shallow water with depth and current refraction
  !! P_SPHER_SHALLOW_CURR
  !!
  SUBROUTINE wave_refraction(p_patch, wave_config, wave_num_c, gv_c, depth, depth_grad, tracer)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
         routine = modname//':wave_refraction'

    TYPE(t_patch),               INTENT(IN)  :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)  :: wave_config
    REAL(wp),                    INTENT(IN)  :: wave_num_c(:,:,:)
    REAL(wp),                    INTENT(IN)  :: gv_c(:,:,:)         ! group velocity at cell centers
    REAL(wp),                    INTENT(IN)  :: depth(:,:)
    REAL(vp),                    INTENT(IN)  :: depth_grad(:,:,:,:) ! bathymetry gradient (2,jc,jk,jb)
    REAL(wp),                    INTENT(INOUT) :: tracer(:,:,:,:)


    TYPE(t_wave_config), POINTER :: wc => NULL()

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc,jb,jf,jd,jk
    INTEGER :: jt,jtm1,jtp1                       !< tracer index

    REAL(wp) :: DELTHR, DELTH, DELTR, DELTH0, sm, sp, ak, akd, DTP, DTM, dDTC, temp, tsihkd
    REAL(wp) :: thdd(nproma,wave_config%ndirs)
    REAL(wp) :: delta_ref(nproma,wave_config%nfreqs*wave_config%ndirs)

    wc => wave_config

    DELTH = 2.0_wp*pi/REAL(wc%ndirs,wp)
    DELTR = DELTH * grid_sphere_radius
    DELTH0 = 0.5_wp * dtime / DELTR
    DELTHR = 0.5_wp * dtime / DELTH

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)
    jk         = p_patch%nlev


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jf,jd,jc,jt,i_startidx,i_endidx,temp,ak,akd,tsihkd,thdd, &
!$OMP            sm,sp,jtm1,jtp1,dtp,dtm,dDTC,delta_ref) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
           &                 i_startidx, i_endidx, i_rlstart, i_rlend)
      DO jf = 1,wc%nfreqs
        DO jd = 1,wc%ndirs
          DO jc = i_startidx, i_endidx

            temp = (SIN(wc%dirs(jd)) + SIN(wc%dirs(wc%dir_neig_ind(2,jd)))) * depth_grad(2,jc,jk,jb) &
                 - (COS(wc%dirs(jd)) + COS(wc%dirs(wc%dir_neig_ind(2,jd)))) * depth_grad(1,jc,jk,jb)

            ak = wave_num_c(jc,jb,jf)
            akd = ak * depth(jc,jb)

            IF (akd <= 10.0_wp) THEN
              tsihkd = (pi2 * wc%freqs(jf))/SINH(2.0_wp*akd)
            ELSE
              tsihkd = 0.0_wp
            END IF

            thdd(jc,jd) = temp * tsihkd

          END DO !jc
        END DO !jd


        DO jd = 1,wc%ndirs

          sm = DELTH0 * (SIN(wc%dirs(jd)) + SIN(wc%dirs(wc%dir_neig_ind(1,jd)))) !index of direction - 1
          sp = DELTH0 * (SIN(wc%dirs(jd)) + SIN(wc%dirs(wc%dir_neig_ind(2,jd)))) !index of direction + 1

          jt   = wc%tracer_ind(jd,jf)
          jtm1 = wc%tracer_ind(wc%dir_neig_ind(1,jd),jf)
          jtp1 = wc%tracer_ind(wc%dir_neig_ind(2,jd),jf)

          DO jc = i_startidx, i_endidx

            DTP = SIN(p_patch%cells%center(jc,jb)%lat) / COS(p_patch%cells%center(jc,jb)%lat) * gv_c(jc,jb,jf)

            DTM = DTP * SM + thdd(jc,wc%dir_neig_ind(1,jd)) * DELTHR
            DTP = DTP * SP + thdd(jc,jd) * DELTHR

            dDTC = -MAX(0._wp , DTP) + MIN(0._wp , DTM)
            DTP  = -MIN(0._wp , DTP)
            DTM  =  MAX(0._wp , DTM)

            delta_ref(jc,jt) = dDTC * tracer(jc,jk,jb,jt) &
                 + DTM * tracer(jc,jk,jb,jtm1)  &
                 + DTP * tracer(jc,jk,jb,jtp1)
          END DO !jc
        END DO !jd
      END DO !jf

      DO jf = 1,wc%nfreqs
        DO jd = 1,wc%ndirs
          jt = wc%tracer_ind(jd,jf)
          DO jc = i_startidx, i_endidx
            tracer(jc,jk,jb,jt) = tracer(jc,jk,jb,jt) + delta_ref(jc,jt)
          END DO !jc
        END DO !jd
      END DO !jf

    END DO !jb
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL
  END SUBROUTINE wave_refraction

END MODULE mo_wave_physics
