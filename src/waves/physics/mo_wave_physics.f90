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
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, min_rlcell, min_rledge
  USE mo_loopindices,         ONLY: get_indices_c, get_indices_e
  USE mo_kind,                ONLY: wp
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: dtime
  USE mo_physical_constants,  ONLY: grav
  USE mo_math_constants,      ONLY: deg2rad, dbl_eps, pi, pi2
  USE mo_idx_list,            ONLY: t_idx_list1D

  USE mo_wave_types,          ONLY: t_wave_diag
  USE mo_wave_config,         ONLY: t_wave_config
  USE mo_wave_forcing_types,  ONLY: t_wave_forcing
  USE mo_wave_constants,      ONLY: ALPHA, EPS1, EMIN

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: air_sea
  PUBLIC :: new_spectrum
  PUBLIC :: total_energy
  PUBLIC :: mean_frequency_energy
  PUBLIC :: wave_group_velocity_c, wave_group_velocity_e, wave_group_velocity_nt
  PUBLIC :: wave_group_velocity_bnd

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_wave_physics'

CONTAINS

  !>
  !! Calculation of wave group velocity
  !!
  !! Calculation of shallow water cell centered
  !! wave group velocity
  !!
  !! @par Revision History
  !! Initial revision by Mikhail Dobrynin, DWD (2023-02-01)
  !! Vectorization by Daniel Reinert, DWD (2023-02-13)
  !!
  SUBROUTINE wave_group_velocity_c(p_patch, p_config, bathymetry_c, gv_c)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
         &  routine = modname//':wave_group_velocity_c'

    TYPE(t_patch),       INTENT(IN)   :: p_patch
    TYPE(t_wave_config), INTENT(IN)   :: p_config
    REAL(wp),            INTENT(IN)   :: bathymetry_c(:,:) !< bathymetric height at cell centers (nproma,nblks_c) ( m )
    REAL(wp),            INTENT(INOUT):: gv_c(:,:,:)       !< group velocity (nproma,nlev,nblks_c,ndirs*nfreqs)  ( m/s )

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb,jf,jc
    INTEGER :: nfreqs

    REAL(wp) :: gh, ak, akd, gv
    REAL(wp) :: wave_num(SIZE(bathymetry_c,1))

    gh = grav / (4.0_wp * pi)

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

    nfreqs = p_config%nfreqs

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jf,jc,i_startidx,i_endidx,ak,akd,gv,wave_num) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jf = 1, nfreqs

        ! get wave number as a function of circular frequency and water depth
        wave_num(i_startidx:i_endidx) = wave_number(pi2 * p_config%freqs(jf), bathymetry_c(i_startidx:i_endidx,jb))

        DO jc = i_startidx, i_endidx
          ! shallow water group velocity
          ak = wave_num(jc)
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
  !! @par Revision History
  !! Initial revision by Mikhail Dobrynin, DWD (2023-02-01)
  !! Vectorization by Daniel Reinert, DWD (2023-02-13)
  !!
 SUBROUTINE wave_group_velocity_e(p_patch, p_config, bathymetry_e, gv_e)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
         &  routine = modname//':wave_group_velocity_e'

    TYPE(t_patch),       INTENT(IN)   :: p_patch
    TYPE(t_wave_config), INTENT(IN)   :: p_config
    REAL(wp),            INTENT(IN)   :: bathymetry_e(:,:) !< bathymetric height at cell centers (nproma,nblks_c) ( m )
    REAL(wp),            INTENT(INOUT):: gv_e(:,:,:)     !< group velocity (nproma,nlev,nblks_c,ndirs*nfreqs)  ( m/s )

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb,jf,je
    INTEGER :: nfreqs

    REAL(wp) :: gh, ak, akd, gv
    REAL(wp) :: wave_num(SIZE(bathymetry_e,1))

    gh = grav / (4.0_wp * pi)

    i_rlstart  = 1
    i_rlend    = min_rledge
    i_startblk = p_patch%edges%start_block(i_rlstart)
    i_endblk   = p_patch%edges%end_block(i_rlend)

    nfreqs = p_config%nfreqs

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jf,je,i_startidx,i_endidx,ak,akd,gv,wave_num) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_e( p_patch, jb, i_startblk, i_endblk,           &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jf = 1, nfreqs

        ! get wave number as a function of circular frequency and water depth
        wave_num(i_startidx:i_endidx) = wave_number(pi2 * p_config%freqs(jf), bathymetry_e(i_startidx:i_endidx,jb))

        DO je = i_startidx, i_endidx
          ! shallow water group velocity
          ak = wave_num(je)
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
  !! @par Revision History
  !! Initial revision by Mikhail Dobrynin, DWD (2023-02-01)
  !! Vectorization by Daniel Reinert, DWD (2023-02-13)
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
          jt = p_config%get_tracer_id(jd,jf)
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
  !! @par Revision History
  !! Initial revision by Mikhail Dobrynin, DWD (2023-02-15)
  !! Vectorization by Daniel Reinert, DWD (2023-02-16)
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
    INTEGER :: ile(nproma), ibe(nproma), e_orient(nproma)
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
          jt = p_config%get_tracer_id(jd,jf)
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
!!$  !! @par Revision History
!!$  !! Initial revision by Mikhail Dobrynin, DWD (2023-02-15)
!!$  !! Vectorization by Daniel Reinert, DWD (2023-02-xx)
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
  !! @par Revision History
  !! Initial revision by Mikhail Dobrynin, DWD (2019-09-05)
  !! Vectorization by Daniel Reinert, DWD (2023-01-30)
  !!
  SUBROUTINE air_sea(p_patch, wave_config, wsp10m, tauw, ustar, z0)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'air_sea'

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
    alphaog = ALPHA / grav

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,iter,i_startidx,i_endidx,xkutop,ustold,tauold,ustm1, &
!$OMP            l_converged,x,z0ch,z0vis,z0tot,xologz0,f,zz,delf,taunew) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      ! initialization
      DO jc = i_startidx, i_endidx
        xkutop(jc) = wave_config%XKAPPA * wsp10m(jc,jb) !xkutop = XKAPPA*UTOP(IJ)
        ustold     = wsp10m(jc,jb) * SQRT(ACD + BCD * wsp10m(jc,jb))
        tauold(jc) = MAX(ustold**2.0_wp,tauw(jc,jb)+EPS1) !MAX(ustold**2,TAUW(IJ)+EPS1)
        ustar(jc,jb) = SQRT(tauold(jc))
        ustm1(jc) = 1.0_wp/MAX(ustar(jc,jb),EPSUS)
        !
        l_converged(jc) = .FALSE.
      ENDDO

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
          delf = 1.0_wp - xkutop(jc) * xologz0**2.0_wp * zz

          ustar(jc,jb) = ustar(jc,jb) - f / delf

          taunew = MAX(ustar(jc,jb)**2.0_wp,tauw(jc,jb) + EPS1)

          ustar(jc,jb) = SQRT(taunew)

          IF (ABS(taunew-tauold(jc))<= dbl_eps) l_converged(jc) = .TRUE.
          ustm1(jc)  = 1.0_wp/MAX(ustar(jc,jb),EPSUS)
          tauold(jc) = taunew
        END DO  !jc
        !
        IF ( ALL(l_converged(i_startidx:i_endidx)) ) EXIT
        !
      END DO  !iter
      !
      DO jc = i_startidx, i_endidx
        z0(jc,jb) = z0ch(jc)
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
  SUBROUTINE mean_frequency_energy(p_patch, wave_config, tracer, llws, emean, emeanws, femean, femeanws)
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: routine = 'mean_frequency'

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    REAL(wp), INTENT(IN)  :: tracer(:,:,:,:) !energy spectral bins (nproma,nlev,nblks_c,ntracer)
    INTEGER,  INTENT(IN)  :: llws(:,:,:)   !=1 where wind_input is positive (nproma,nblks_c,ntracer)
    REAL(wp), INTENT(IN)  :: emean(:,:)    !total energy (nproma,nblks_c)
    REAL(wp), INTENT(IN)  :: emeanws(:,:)  !total windsea energy (nproma,nblks_c)
    REAL(wp), INTENT(OUT) :: femean(:,:)   !mean frequency energy (nproma,nblks_c)
    REAL(wp), INTENT(OUT) :: femeanws(:,:) !mean windsea frequency energy (nproma,nblks_c)

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
        CALL wc%get_tracer_ids_freq(jf,list_tr)

        DO jc = i_startidx, i_endidx
          temp(jc,jf)   = 0._wp
          temp_1(jc,jf) = 0._wp
        ENDDO

        DO n=1,list_tr%ncount
          jt = list_tr%idx(n)
          DO jc = i_startidx, i_endidx
            temp(jc,jf) = temp(jc,jf) + tracer(jc,jk,jb,jt)
            IF (llws(jc,jb,jt)==1) THEN
              temp_1(jc,jf) = temp_1(jc,jf) + tracer(jc,jk,jb,jt)
            ENDIF
          ENDDO
        ENDDO  ! n

        CALL list_tr%finalize()
      END DO  ! jf

      DO jc = i_startidx, i_endidx
        femean(jc,jb)   = wc%MM1_TAIL * temp(jc,wc%nfreqs)
        femeanws(jc,jb) = wc%MM1_TAIL * temp_1(jc,wc%nfreqs)
      ENDDO

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
  !! @par Revision History
  !! Initial revision by Mikhail Dobrynin, DWD (2019-09-05)
  !! Optimization and vectorization by Daniel Reinert, DWD (2023-01-30)
  !!
  SUBROUTINE total_energy(p_patch, wave_config, tracer, llws, emean, emeanws)
    ! in  p_patch, p_prog%tracer, wave_config%freq_ind, wave_config%MO_TAIL
    !out  p_diag%emean(jc,jb), p_diag%emeanws(jc,jb)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = modname//'total_energy'

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    REAL(wp), INTENT(IN)  :: tracer(:,:,:,:) !energy spectral bins (nproma,nlev,nblks_c,ntracer)
    INTEGER,  INTENT(IN)  :: llws(:,:,:)   !=1 where wind_input is positive (nproma,nblks_c,ntracer)
    REAL(wp), INTENT(OUT) :: emean(:,:)    !total energy (nproma,nblks_c)
    REAL(wp), INTENT(OUT) :: emeanws(:,:)  !total windsea energy (nproma,nblks_c)


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
        CALL wc%get_tracer_ids_freq(jf, list_tr)

        ! initialization
        DO jc = i_startidx, i_endidx
          sum1(jc,jf) = 0._wp
          sum2(jc,jf) = 0._wp
        ENDDO

        DO n=1,list_tr%ncount
          jt = list_tr%idx(n)
          DO jc = i_startidx, i_endidx
            sum1(jc,jf) = sum1(jc,jf) + tracer(jc,jk,jb,jt)
            IF (llws(jc,jb,jt) == 1) THEN
              sum2(jc,jf) = sum2(jc,jf) + tracer(jc,jk,jb,jt)
            ENDIF
          ENDDO
        ENDDO  ! n

        CALL list_tr%finalize()
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
  SUBROUTINE new_spectrum(p_patch, wave_config, p_diag, p_forcing, tracer)
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'new_spectrum'

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    TYPE(t_wave_diag),           INTENT(IN)    :: p_diag
    TYPE(t_wave_forcing),        INTENT(IN)    :: p_forcing
    REAL(wp),                    INTENT(INOUT) :: tracer(:,:,:,:)

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

            sprd = MAX(0._wp, COS(wc%dirs(jd)-p_forcing%dir10m(jc,jb)*deg2rad) )**2

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
  !! Calculation of wave number.
  !!
  !! Wave number as a function of circular frequency and water depth
  !! Adaptation of WAM 4.5 code.
  !!
  !! @par Revision History
  !! Initial revision by Mikhail Dobrynin, DWD (2019-09-05)
  !! Optimization and vectorization by Daniel Reinert, DWD (2023-02-13)
  !!
  FUNCTION wave_number(OM, depth) RESULT(wave_num)
    ! ---------------------------------------------------------------------------- !
    !   AKI - FUNCTION TO COMPUTE WAVE NUMBER.                                     !
    !                                                                              !
    !     G. KOMEN, P. JANSSEN   KNMI              01/06/1986                       !
    !     Mikhail Dobrynin, DWD, adjusted for ICON 19/09/2021                       !
    !                                                                              !
    !     PURPOSE.                                                                 !
    !     -------                                                                  !
    !       WAVE NUMBER AS FUNCTION OF CIRCULAR FREQUENCY AND WATER DEPTH.         !
    !                                                                              !
    !     METHOD.                                                                  !
    !     -------                                                                  !
    !       NEWTONS METHOD TO SOLVE THE DISPERSION RELATION IN SHALLOW WATER.      !

    REAL(wp)          , INTENT(IN)   :: OM       !! CIRCULAR FREQUENCY 2*pi*freq (nfreqs)
    REAL(wp)          , INTENT(IN)   :: depth(:) !! bathymetric height at cell centers (nproma)

    REAL(wp), PARAMETER :: EBS = 0.0001_wp     !! RELATIVE ERROR LIMIT OF NEWTON'S METHOD.
    REAL(wp), PARAMETER :: DKMAX = 40.0_wp     !! MAXIMUM VALUE OF DEPTH*WAVENUMBER.
    REAL(wp)            :: BO, TH, STH
    INTEGER             :: jc                  ! loop index for cells
    INTEGER             :: ncell               ! loop bounds

    LOGICAL             :: l_converged(SIZE(depth))
    REAL(wp)            :: wave_num(SIZE(depth))
    REAL(wp)            :: AKP(SIZE(depth))
    ! ---------------------------------------------------------------------------- !
    !     1. START WITH MAXIMUM FROM DEEP AND EXTREM SHALLOW WATER WAVE NUMBER.    !

    ncell = SIZE(depth)

    DO jc=1,ncell
      wave_num(jc) = MAX( OM**2/(4.0_wp* grav), OM/(2.0_wp*SQRT(grav*depth(jc))) )
      !
      ! initialization
      AKP(jc)         = 10000.0_wp
      l_converged(jc) = .FALSE.
    ENDDO

    ! ---------------------------------------------------------------------------- !
    !     2. ITERATION LOOP.                                                       !

    DO WHILE (ANY(l_converged(1:ncell)).EQV..FALSE.)

      DO jc=1,ncell
        IF (.NOT.l_converged(jc)) THEN
          BO = depth(jc)*wave_num(jc)
          IF (BO > DKMAX) THEN
            wave_num(jc) = OM**2/grav
            l_converged(jc) = .TRUE.
          ELSE
            AKP(jc) = wave_num(jc)
            TH = grav*wave_num(jc)*TANH(BO)
            STH = SQRT(TH)
            wave_num(jc) = wave_num(jc) &
              &                + (OM-STH)*STH*2.0_wp / (TH/wave_num(jc) + grav*BO/COSH(BO)**2)
            ! check for converged solution
            l_converged(jc) = (ABS(AKP(jc)-wave_num(jc)) <= EBS*wave_num(jc))
          END IF
        END IF  ! l_converged
      ENDDO  !jc

    ENDDO !while

  END FUNCTION wave_number



!!$  !>
!!$  !! Calculation of wave number.
!!$  !!
!!$  !! Wave number as a function of circular frequency and water depth
!!$  !! Adaptation of WAM 4.5 code.
!!$  !!
!!$  !! @par Revision History
!!$  !! Initial revision by Mikhail Dobrynin, DWD (2019-09-05)
!!$  !! Optimization and vectorization by Daniel Reinert, DWD (2023-XX-XX)
!!$  !!
!!$  REAL(wp) FUNCTION wave_number(OM, depth)
!!$    ! ---------------------------------------------------------------------------- !
!!$    !   AKI - FUNCTION TO COMPUTE WAVE NUMBER.                                     !
!!$    !                                                                              !
!!$    !     G. KOMEN, P. JANSSEN   KNMI              01/06/1986                       !
!!$    !     Mikhail Dobrynin, DWD, adjusted for ICON 19/09/2021                       !
!!$    !                                                                              !
!!$    !     PURPOSE.                                                                 !
!!$    !     -------                                                                  !
!!$    !       WAVE NUMBER AS FUNCTION OF CIRCULAR FREQUENCY AND WATER DEPTH.         !
!!$    !                                                                              !
!!$    !     METHOD.                                                                  !
!!$    !     -------                                                                  !
!!$    !       NEWTONS METHOD TO SOLVE THE DISPERSION RELATION IN SHALLOW WATER.      !
!!$
!!$    REAL(wp)          , INTENT(IN) :: OM    !! CIRCULAR FREQUENCY 2*pi*freq
!!$    REAL(wp)          , INTENT(IN) :: depth
!!$
!!$    REAL(wp), PARAMETER :: EBS = 0.0001_wp     !! RELATIVE ERROR LIMIT OF NEWTON'S METHOD.
!!$    REAL(wp), PARAMETER :: DKMAX = 40.0_wp     !! MAXIMUM VALUE OF DEPTH*WAVENUMBER.
!!$    REAL(wp)            :: AKP, BO, TH, STH
!!$
!!$    ! ---------------------------------------------------------------------------- !
!!$    !     1. START WITH MAXIMUM FROM DEEP AND EXTREM SHALLOW WATER WAVE NUMBER.    !
!!$
!!$    wave_number = MAX( OM**2/(4.0_wp* grav), OM/(2.0_wp*SQRT(grav*depth)) )
!!$
!!$    ! ---------------------------------------------------------------------------- !
!!$    !     2. ITERATION LOOP.                                                       !
!!$
!!$    AKP = 10000.0_wp
!!$    DO WHILE (ABS(AKP-wave_number).gt.EBS*wave_number)
!!$       BO = depth*wave_number
!!$       IF (BO.GT.DKMAX) THEN
!!$          wave_number = OM**2/grav
!!$          EXIT
!!$       ELSE
!!$          AKP = wave_number
!!$          TH = grav*wave_number*TANH(BO)
!!$          STH = SQRT(TH)
!!$          wave_number = wave_number + (OM-STH)*STH*2.0_wp / (TH/wave_number + grav*BO/COSH(BO)**2)
!!$       END IF
!!$    END DO
!!$
!!$  END FUNCTION wave_number


END MODULE mo_wave_physics
