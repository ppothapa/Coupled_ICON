!>
!!
!! @author Hui Wan, MPI-M
!!
!! @par Revision History
!! First version by Hui Wan (2010-09)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_vdiff_solver

  USE mo_kind,              ONLY: wp
  USE mo_impl_constants,    ONLY: SUCCESS
  USE mo_exception,         ONLY: message, message_text, finish
  USE mo_physical_constants,ONLY: grav, rgrav, cpd, cpv
  USE mo_echam_vdiff_params,ONLY: tke_min, &
    &                             tpfac1, tpfac2, tpfac3, cchar, z0m_min
  USE mo_echam_phy_config,  ONLY: phy_config => echam_phy_config, get_lebudget
  USE mo_echam_phy_memory,  ONLY: cdimissval

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: init_vdiff_solver         !< subroutine
  PUBLIC :: cleanup_vdiff_solver      !< subroutine
  PUBLIC :: matrix_setup_elim         !< subroutine
  PUBLIC :: matrix_to_richtmyer_coeff !< subroutine
  PUBLIC :: rhs_setup, rhs_elim       !< subroutines
  PUBLIC :: rhs_bksub                 !< subroutines
  PUBLIC :: vdiff_tendencies          !< subroutine
  PUBLIC :: nvar_vdiff, nmatrix       !< parameters
  PUBLIC :: ih,iqv,iu,iv              !< parameters
  PUBLIC :: imh,imqv, imuv            !< parameters

  ! Module variables

  INTEGER :: nvar_vdiff    !< total number of variables affected by turbulent mixing
  INTEGER :: iu, iv, ih, iqv
  INTEGER :: ixl, ixi, ixv
  INTEGER :: itke, ithv
  INTEGER :: itrc_start
  INTEGER :: nmatrix
  INTEGER :: imh, imqv, imuv

  INTEGER, ALLOCATABLE :: matrix_idx(:)    !< shape: (nvar_vdiff)
  INTEGER, ALLOCATABLE :: ibtm_var  (:)    !< shape: (nvar_vdiff)
  INTEGER, ALLOCATABLE :: ibtm_mtrx (:)    !< shape: (nmatrix)

  CHARACTER(len=*), PARAMETER :: thismodule = 'mo_vdiff_solver'

CONTAINS
  !>
  !!
  !! In this prototype it is assumed that the following variables are subject to
  !! turbulent mixing:
  !!
  !!    variables                              |  # of variables
  !! -----------------------------------------------------------
  !!   u, v, T, qv                             |  4
  !!   all hydrometeors                        |  khydromet
  !!   variance of cloud droplet concentration |  1
  !!   TKE                                     |  1
  !!   variance of theta_v                     |  1
  !!   additional tracers                      |  ktrac
  !! -----------------------------------------------------------
  !!
  SUBROUTINE init_vdiff_solver( khydromet, ktrac, klev )

    INTEGER,INTENT(IN) :: khydromet, ktrac
    INTEGER,INTENT(IN) :: klev
    INTEGER :: ist

    !------------------------------------------
    ! Set up index for prognostic variables
    !------------------------------------------

    nvar_vdiff = 7 + khydromet + ktrac

    iu   = 1;   iv   = 2
    ixl  = 3;   ixi  = 4;  ixv = 5
    itke = 6;   ithv = 7
    ih   = 8;   iqv  = 9

    !>KF suggestion
    IF((7 + khydromet) > iqv ) &
    CALL finish( TRIM(thismodule),'matrix for vdiff is not properly defined')

    IF(ktrac > 0)  THEN
      itrc_start = 7 + khydromet +1
    ELSE
      itrc_start = 7 + khydromet
    ENDIF
    !<KF

    !-------------------------------------------------------------------
    ! # of vertical levels on which the prognostic equations are solved
    !-------------------------------------------------------------------

    ALLOCATE( ibtm_var(nvar_vdiff),STAT=ist )
    IF (ist/=SUCCESS) CALL finish(TRIM(thismodule),&
      & 'Allocation of ibtm_var failed')

    ! momentum, heat, water substances and tracers are solved on        
    ! klev full levels 

    ibtm_var(:)    = klev

    ! TKE and the variance of $\theta_v$ are solved at klev-1 half levels.
    ! The upper and lower boundaries of the atmosphere are excluded.  

    ibtm_var(itke) = klev -1
    ibtm_var(ithv) = klev -1

    !------------------------------------------
    ! Set up matrix indices
    !------------------------------------------

    ALLOCATE( matrix_idx(nvar_vdiff),STAT=ist )
    IF (ist/=SUCCESS) CALL finish(TRIM(thismodule),&
      & 'Allocation of matrix_idx failed')

    matrix_idx(iu)   = 1  ; imuv = 1
    matrix_idx(iv)   = 1  ! u and v share the same exchange coeff.
    matrix_idx(ixl)  = 2
    matrix_idx(ixi)  = 2  ! cloud water and ice share the same exchange coeff.
    matrix_idx(ixv)  = 3
    matrix_idx(itke) = 4
    matrix_idx(ithv) = 5
    matrix_idx(ih)   = 6 ; imh  = 6
    matrix_idx(iqv)  = 7 ; imqv = 7

    IF (ktrac>0) matrix_idx(nvar_vdiff-ktrac+1:nvar_vdiff) = 2

    nmatrix = 7    ! total number of matrices

    !---------------------------------------------------------------------------
    ! # of vertical levels on which elimination of the coefficient will be done
    !---------------------------------------------------------------------------

    ALLOCATE( ibtm_mtrx(nmatrix),STAT=ist )
    IF (ist/=SUCCESS) CALL finish(TRIM(thismodule),&
    & 'Allocation of ibtm_mtrx failed')

    ibtm_mtrx(:)                = klev
    ibtm_mtrx(matrix_idx(itke)) = klev -1
    ibtm_mtrx(matrix_idx(ithv)) = klev -1

  END SUBROUTINE init_vdiff_solver
  !-------------
  !>
  !!
  SUBROUTINE cleanup_vdiff_solver

    INTEGER :: ist

    DEALLOCATE( matrix_idx,ibtm_mtrx,ibtm_var, STAT=ist)
    IF (ist/=SUCCESS) CALL finish('cleanup_vdiff_solver','Deallocation failed')

  END SUBROUTINE cleanup_vdiff_solver
  !-------------
  !>
  !!
  !! Set up coeffient matrix of the linear algebraic system and 
  !! perform Gauss elimination. For moisture, the last row of the 
  !! matrix (aa_btm) can not be finished yet because the evapotranspiration 
  !! coefficients "cair" and "csat" are not yet available. Thus for this variable 
  !! elimination is performed only till level klev-1. 
  !! 
  SUBROUTINE matrix_setup_elim( kproma, kbdim, klev, klevm1,  &! in
                              & ksfc_type, itop,              &! in
                              & pcfm, pcfh, pcfh_tile, pcfv,  &! in
                              & pcftke, pcfthv,               &! in
                              & pprfac, prdpm, prdph,         &! in
                              & prmairm, prmairh, prmdrym,    &! in
                              & aa, aa_btm                    )! out
    ! Arguments

    INTEGER, INTENT(IN) :: kproma, kbdim, klev, klevm1, ksfc_type
    INTEGER, INTENT(IN) :: itop

    REAL(wp),INTENT(IN) :: pcfm     (kbdim,klev)      !< exchange coeff. for u, v
    REAL(wp),INTENT(IN) :: pcfh     (kbdim,klevm1)    !< exchange coeff. for heat and tracers
    REAL(wp),INTENT(IN) :: pcfh_tile(kbdim,ksfc_type) !< exchange coeff. for heat and qv, at surface
    REAL(wp),INTENT(IN) :: pcfv     (kbdim,klev)      !< exchange coeff. for total water variance
    REAL(wp),INTENT(IN) :: pcftke   (kbdim,klev)      !< exchange coeff. for TKE
    REAL(wp),INTENT(IN) :: pcfthv   (kbdim,klev)      !< exchange coeff. for variance of theta_v
    REAL(wp),INTENT(IN) :: pprfac   (kbdim,klev)      !< prefactor for the exchange coefficients
    REAL(wp),INTENT(IN) :: prdpm    (kbdim,klev)      !< reciprocal of layer thickness, full levels
    REAL(wp),INTENT(IN) :: prdph    (kbdim,klevm1)    !< reciprocal of layer thickness, half levels
    REAL(wp),INTENT(IN) :: prmairm  (kbdim,klev)      !< reciprocal of layer air mass, full levels
    REAL(wp),INTENT(IN) :: prmairh  (kbdim,klevm1)    !< reciprocal of layer dry aor mass, half levels
    REAL(wp),INTENT(IN) :: prmdrym  (kbdim,klev)      !< reciprocal of layer dry air mass, full levels

    REAL(wp),INTENT(INOUT) :: aa    (kbdim,klev,3,nmatrix) !< exchange coeff. matrices    out
    REAL(wp),INTENT(INOUT) :: aa_btm(kbdim,3,ksfc_type,imh:imqv)   !  out
                                     !< last (the klev-th) row of the coeff. matrices 
                                     !< of dry static energy and moisture

    ! Local variables

    REAL(wp) :: zkstar (kbdim,itop-1:klev)     !< scaled exchange coeff on half-levels
    REAL(wp) :: zkh    (kbdim,itop-1:klevm1)   !< scaled exchange doeff on full-levels, 
                                               !< for TKE and variance of theta_v
    INTEGER  :: im             !< index of coefficient matrix
    INTEGER  :: jc, jk, jsfc   !< loop indices
    INTEGER  :: jkm1, jmax

    !-----------------------------------------------------------------------
    ! For all prognostic variables: no turbulent flux at the upper boundary
    !-----------------------------------------------------------------------

    zkstar(1:kproma,itop-1) = 0._wp

    !-----------------------------------------------------------------------
    ! For momentum: surface flux is considered
    !-----------------------------------------------------------------------
    im = matrix_idx(iu)    ! also = matrix_idx(iv)
    zkstar(1:kproma,itop:klev) = pprfac(1:kproma,itop:klev)  &
                               &  *pcfm(1:kproma,itop:klev)

    DO jk = itop,klev
      DO jc = 1,kproma
        aa(jc,jk,1,im) = -zkstar(jc,jk-1)*prmairm(jc,jk)  ! -K*_{k-1/2}/dm_k
        aa(jc,jk,3,im) = -zkstar(jc,jk  )*prmairm(jc,jk)  ! -K*_{k+1/2}/dm_k
        aa(jc,jk,2,im) = 1._wp - aa(jc,jk,1,im) - aa(jc,jk,3,im)
      ENDDO
    ENDDO

    !---------------------------------------------------------------------
    ! Dry static energy: surface fluxes on different surface types 
    ! are handled separately. 
    !---------------------------------------------------------------------
    im = imh
    zkstar(1:kproma,itop:klevm1) =  pprfac(1:kproma,itop:klevm1) &
                                 &   *pcfh(1:kproma,itop:klevm1)
    DO jk = itop,klevm1
      DO jc = 1,kproma
        aa(jc,jk,1,im) = -zkstar(jc,jk-1)*prmairm(jc,jk)  ! -K*_{k-1/2}/dm_k
        aa(jc,jk,3,im) = -zkstar(jc,jk  )*prmairm(jc,jk)  ! -K*_{k+1/2}/dm_k
        aa(jc,jk,2,im) = 1._wp - aa(jc,jk,1,im) - aa(jc,jk,3,im)
      ENDDO
    ENDDO

    ! Set the bottom row of the coeff matrix. The same formula applies 
    ! for all surface types (land, water, ice).

    jk = klev
    DO jsfc = 1,ksfc_type
      DO jc = 1,kproma
        aa_btm(jc,1,jsfc,im) = -zkstar(jc,jk-1)*prmairm(jc,jk)    ! -K*_{k-1/2}/dm_k
        aa_btm(jc,3,jsfc,im) = -pcfh_tile(jc,jsfc)*pprfac(jc,jk)*prmairm(jc,jk)
        aa_btm(jc,2,jsfc,im) = 1._wp - aa_btm(jc,1,jsfc,im) - aa_btm(jc,3,jsfc,im)
      ENDDO
    ENDDO

    !---------------------------------------------------------------------
    ! Moisture: different surface types are handled separately.
    !---------------------------------------------------------------------
    im = imqv
    zkstar(1:kproma,itop:klevm1) =  pprfac(1:kproma,itop:klevm1) &
                                 &   *pcfh(1:kproma,itop:klevm1)
    DO jk = itop,klevm1
      DO jc = 1,kproma
        aa(jc,jk,1,im) = -zkstar(jc,jk-1)*prmdrym(jc,jk)  ! -K*_{k-1/2}/dm_k
        aa(jc,jk,3,im) = -zkstar(jc,jk  )*prmdrym(jc,jk)  ! -K*_{k+1/2}/dm_k
        aa(jc,jk,2,im) = 1._wp - aa(jc,jk,1,im) - aa(jc,jk,3,im)
      ENDDO
    ENDDO

    ! Bottom row of the matrix: finish the setup over water and ice;
    ! do part of the computation for land. Later in subroutine
    ! matrix_to_richtmyer_coeff, aa_btm(:,3,idx_land,imqv) will be 
    ! modified, and aa_btm(:,2,idx_land,imqv) re-computed.

    jk = klev
    DO jsfc = 1,ksfc_type
      DO jc = 1,kproma
        aa_btm(jc,1,jsfc,im) = -zkstar(jc,jk-1)*prmdrym(jc,jk)    ! -K*_{k-1/2}/dm_k
        aa_btm(jc,3,jsfc,im) = -pcfh_tile(jc,jsfc)*pprfac(jc,jk)*prmdrym(jc,jk)
        aa_btm(jc,2,jsfc,im) = 1._wp - aa_btm(jc,1,jsfc,im) - aa_btm(jc,3,jsfc,im)
      ENDDO
    ENDDO

    !----------------------------------------------------------------------
    ! For all advected tracers except water vapour: no turbulent flux at 
    ! the surface.
    !----------------------------------------------------------------------
    im = matrix_idx(ixl)
    zkstar(1:kproma,klev) = 0._wp  ! lower boundary, no turbulent flux

    DO jk = itop,klev
      DO jc = 1,kproma
        aa(jc,jk,1,im) = -zkstar(jc,jk-1)*prmdrym(jc,jk)  ! -K*_{k-1/2}/dm_k
        aa(jc,jk,3,im) = -zkstar(jc,jk  )*prmdrym(jc,jk)  ! -K*_{k+1/2}/dm_k
        aa(jc,jk,2,im) = 1._wp - aa(jc,jk,1,im) - aa(jc,jk,3,im)
      ENDDO
    ENDDO

    !----------------------------------------------------------------------
    ! For total water variance: no surface flux. The exchange coefficient
    ! pcfv has been set to to zero in subroutine sfc_exchange_coeff, which
    ! automatically leads to zkstar(:,klev) = 0._wp, thus no additional 
    ! attention is needed here.
    !----------------------------------------------------------------------
    im = matrix_idx(ixv)
    zkstar(1:kproma,itop:klev) = pprfac(1:kproma,itop:klev) &
                               &  *pcfv(1:kproma,itop:klev)
    DO jk = itop,klev
      DO jc = 1,kproma
        aa(jc,jk,1,im) = -zkstar(jc,jk-1)*prmdrym(jc,jk)
        aa(jc,jk,3,im) = -zkstar(jc,jk  )*prmdrym(jc,jk)
        aa(jc,jk,2,im) = 1._wp - aa(jc,jk,1,im) - aa(jc,jk,3,im)
      ENDDO
    ENDDO

    !------------------------------------------------------------------------
    ! For TKE: Note that
    ! - Vertical averaging is needed to convert exchange coefficient from
    !   half to full levels, because TKE equation is solved on half levels.
    ! - TKE equation is solved only till array subscript klevm1, which
    !   corresponds to half level (klev - 1/2), i.e., the lowest
    !   interface above surface. Surface value of TKE is (already)
    !   computed in subroutine "sfc_exchange_coeff".
    !------------------------------------------------------------------------
    im = matrix_idx(itke)

    zkstar(1:kproma,itop:klev) =  pprfac(1:kproma,itop:klev) &
                               & *pcftke(1:kproma,itop:klev)
    DO jk = itop,klevm1
      zkh(1:kproma,jk) = 0.5_wp*(zkstar(1:kproma,jk)+zkstar(1:kproma,jk+1))
    ENDDO
    zkh(1:kproma,itop-1) = 0._wp  ! upper boundary, no flux

    DO jk = itop,klevm1
      DO jc = 1,kproma
        aa(jc,jk,1,im) = -zkh(jc,jk-1)*prmairh(jc,jk)
        aa(jc,jk,3,im) = -zkh(jc,jk  )*prmairh(jc,jk)
        aa(jc,jk,2,im) = 1._wp - aa(jc,jk,1,im) - aa(jc,jk,3,im)
      ENDDO
    ENDDO

    !------------------------------------------------
    ! For the variance of theta_v (similar to TKE)
    !------------------------------------------------
    im = matrix_idx(ithv)
    zkstar(1:kproma,itop:klev) =  pprfac(1:kproma,itop:klev) &
                               & *pcfthv(1:kproma,itop:klev)
    DO jk = itop,klevm1
      zkh(1:kproma,jk) = 0.5_wp*(zkstar(1:kproma,jk)+zkstar(1:kproma,jk+1))
    ENDDO

    DO jk = itop,klevm1
      DO jc = 1,kproma
        aa(jc,jk,1,im) = -zkh(jc,jk-1)*prmairh(jc,jk)
        aa(jc,jk,3,im) = -zkh(jc,jk  )*prmairh(jc,jk)
        aa(jc,jk,2,im) = 1._wp - aa(jc,jk,1,im) - aa(jc,jk,3,im)
      ENDDO
    ENDDO

    !-----------------------------------------------------------------------------
    ! Gauss elimination for the coefficient matrices at
    ! - vertical levels [itop,klev-2], for TKE and variance of theta_v;
    ! - vertical levels [itop,klev-1], for all the other variables.
    !-----------------------------------------------------------------------------

    DO im = 1,nmatrix
      aa(1:kproma,itop,3,im) = aa(1:kproma,itop,3,im)/aa(1:kproma,itop,2,im)

      jmax = ibtm_mtrx(im) - 1
      DO jk = itop+1,jmax
        jkm1 = jk - 1
        aa(1:kproma,jk,2,im) =  aa(1:kproma,jk,2,im)                       &
                             & -aa(1:kproma,jk,1,im)*aa(1:kproma,jkm1,3,im)
        aa(1:kproma,jk,3,im) =  aa(1:kproma,jk,3,im)/aa(1:kproma,jk,2,im)
      ENDDO
    END DO

    ! Translation for developers who prefer to think in terms of 
    ! the Richtmyer-Morthon formula and are familiar with the paper by
    ! Polcher et al (1998): after this elimination, 
    !  aa(:,1:ibtm_mtrx(im)-1,2,:) becomes C  (Eqn. 17),
    !  aa(:,1:ibtm_mtrx(im)-1,3,:) becomes -A (Eqn. 19).
    ! See subroutine matrix_to_richtmyer_coeff.

  END SUBROUTINE matrix_setup_elim
  !--------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------
  !>
  SUBROUTINE rhs_setup( kproma, kbdim, itop, klev, klevm1,   &! in
                      & ksfc_type, ktrac, ptpfac2, pstep_len,&! in
                      & pum1, pvm1, pcptgz, pqm1,            &! in
                      & pxlm1, pxim1, pxvar, pxtm1, pxt_emis,&! in
                      & prdpm, prmdrym, ptkevn, pzthvvar, aa,&! in
                      & bb, bb_btm                           )! out

    ! Arguments

    INTEGER, INTENT(IN) :: kproma, kbdim, itop, klev, klevm1
    INTEGER, INTENT(IN) :: ksfc_type, ktrac
    REAL(wp),INTENT(IN) :: ptpfac2, pstep_len

    REAL(wp),INTENT(IN) :: pum1     (kbdim,klev)
    REAL(wp),INTENT(IN) :: pvm1     (kbdim,klev)
    REAL(wp),INTENT(IN) :: pcptgz   (kbdim,klev)
    REAL(wp),INTENT(IN) :: pqm1     (kbdim,klev)
    REAL(wp),INTENT(IN) :: pxlm1    (kbdim,klev)
    REAL(wp),INTENT(IN) :: pxim1    (kbdim,klev)
    REAL(wp),INTENT(IN) :: pxvar    (kbdim,klev)
    REAL(wp),INTENT(IN) :: pxtm1    (kbdim,klev,ktrac)
    REAL(wp),INTENT(IN) :: pxt_emis (kbdim,ktrac)
   !REAL(wp),INTENT(IN) :: pxt_emis (kbdim,klev,ktrac) ! backup for later use
    REAL(wp),INTENT(IN) :: ptkevn   (kbdim,klev)
    REAL(wp),INTENT(IN) :: pzthvvar (kbdim,klev)
    REAL(wp),INTENT(IN) :: prdpm    (kbdim,klev)
    REAL(wp),INTENT(IN) :: prmdrym  (kbdim,klev)
    REAL(wp),INTENT(IN) :: aa       (kbdim,klev,3,nmatrix)

    REAL(wp),INTENT(INOUT) :: bb    (kbdim,klev,nvar_vdiff)  ! OUT
    REAL(wp),INTENT(INOUT) :: bb_btm(kbdim,ksfc_type,ih:iqv)  ! OUT

    ! Local variables

    REAL(wp) :: ztmp(kbdim,klev)
    INTEGER  :: jsfc, jt, irhs, im

    !-------------------------------------------------------------------
    ! First handle variables that are defined on full levels
    !-------------------------------------------------------------------
    ! u and v

    bb(1:kproma,itop:klev,iu) = pum1(1:kproma,itop:klev)
    bb(1:kproma,itop:klev,iv) = pvm1(1:kproma,itop:klev)

    ! Hydrometeors and the variance of cloud droplets

    bb(1:kproma,itop:klev,ixl) = pxlm1(1:kproma,itop:klev)
    bb(1:kproma,itop:klev,ixi) = pxim1(1:kproma,itop:klev)
    bb(1:kproma,itop:klev,ixv) = pxvar(1:kproma,itop:klev)

    ! Other tracers

    DO jt = 1,ktrac
       irhs = jt - 1 + itrc_start
       bb(1:kproma,itop:klev,irhs) =  pxtm1(1:kproma,itop:klev,jt)
    ENDDO

    ! Heat and moisture

    bb(1:kproma,itop:klevm1,ih ) = pcptgz(1:kproma,itop:klevm1)
    bb(1:kproma,itop:klevm1,iqv) = pqm1  (1:kproma,itop:klevm1)

    DO jsfc = 1,ksfc_type
       bb_btm(1:kproma,jsfc,ih)  = pcptgz(1:kproma,klev)
       bb_btm(1:kproma,jsfc,iqv) =   pqm1(1:kproma,klev)
    ENDDO

    !-------------------------------------------------------------------
    ! TKE and the variance of theta_v:
    ! These variables are defined at half levels. Array index jk
    ! correspond to half level k+1/2. Thus klev correspond to the
    ! lower boundary. The linear solver only solves till index klevm1.
    !-------------------------------------------------------------------
    im = matrix_idx(itke)
    bb(1:kproma,itop:klevm1,itke) =  ptkevn(1:kproma,itop:klevm1)
    bb(1:kproma,     klevm1,itke) =  bb(1:kproma,klevm1,itke)   &
                                  & -aa(1:kproma,klevm1,3,im)   &
                                  & *ptkevn(1:kproma,klev)

    im = matrix_idx(ithv)
    bb(1:kproma,itop:klevm1,ithv) =  pzthvvar(1:kproma,itop:klevm1)
    bb(1:kproma,     klevm1,ithv) =  bb(1:kproma,klevm1,ithv)     &
                                  & -aa(1:kproma,klevm1,3,im)     &
                                  & *pzthvvar(1:kproma,klev)

    !--------------------------------------------------------------------
    ! Apply the implicitness factor
    !--------------------------------------------------------------------
    !bb     = ptpfac2*bb
    !bb_btm = ptpfac2*bb_btm

     bb(1:kproma,1:klev,  1:itke-1) = ptpfac2*bb(1:kproma,1:klev,  1:itke-1)
     bb(1:kproma,1:klevm1,itke:iqv) = ptpfac2*bb(1:kproma,1:klevm1,itke:iqv)

     IF (ktrac>0) THEN
       bb(1:kproma,1:klev,itrc_start:) = ptpfac2*bb(1:kproma,1:klev,itrc_start:)
     ENDIF

     bb_btm(1:kproma,:,:) = ptpfac2*bb_btm(1:kproma,:,:)

    !--------------------------------------------------------------------
    ! Add tracer emissions
    !--------------------------------------------------------------------
    ! Currently we follow ECHAM in which only the surface emission
    ! is treated in "vdiff".

    ztmp(1:kproma,klev) = prmdrym(1:kproma,klev)*pstep_len

    DO jt = 1,ktrac
       irhs = jt - 1 + itrc_start
       bb(1:kproma,klev,irhs) =         bb(1:kproma,klev,irhs) &
                              & + pxt_emis(1:kproma,jt)        &
                              &      *ztmp(1:kproma,klev)
    ENDDO

    ! Later we may consider treating emission on all vertical levels
    ! in the same way.
    !
    !ztmp(1:kproma,itop:klev) = prmdrym(1:kproma,itop:klev)*pstep_len
    !
    !DO jt = 1,ktrac
    !   irhs = jt - 1 + itrc_start
    !   bb(1:kproma,itop:klev,irhs) =         bb(1:kproma,itop:klev,irhs) &
    !                               & + pxt_emis(1:kproma,itop:klev,jt)   &
    !                               &      *ztmp(1:kproma,itop:klev)
    !ENDDO

  END SUBROUTINE rhs_setup
  !--------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------
  !>
  !!
  !! Gauss elimination of the right-hand-side vector 
  !! using coefficients obtained in subroutine "matrix_setup_elim".
  !!
  !! Translation for developers who prefer to think in terms of 
  !! the Richtmyer-Morthon formula and are familiar with the paper by
  !! Polcher et al (1998): after the elimination at the end of 
  !! subroutine matrix_setup_elim, aa(:,1:ibtm_mtrx(im)-1,2,:) 
  !! became the coeff C defined by Eqn. 17 of Polcher et al (1998). 
  !! It is used in this subroutine to convert the variable bb 
  !! into the Richtmyer coeff B (cf Eqn. 19 of Polcher et al 1998).
  !!
  SUBROUTINE rhs_elim( kproma, kbdim, itop, klev, klevm1, &! in
                     & aa, bb                             )! in, inout

    INTEGER, INTENT(IN)    :: kproma, kbdim, itop, klev, klevm1
    REAL(wp),INTENT(IN)    :: aa(kbdim,klev,3,nmatrix)
    REAL(wp),INTENT(INOUT) :: bb(kbdim,klev,nvar_vdiff)

    REAL(wp) :: znum(kbdim), zden(kbdim)
    INTEGER  :: jvar, im, jk, jkm1, jmax

    ! 1. Vertical levels [itop+1,klev-2] for TKE and variance of theta_v;
    !    [itop+1,klev-1] for all the other variables.

    DO jvar = 1,nvar_vdiff
      im = matrix_idx(jvar)  ! Index of coefficient matrix
      bb(1:kproma,itop,jvar) =  bb(1:kproma,itop,jvar)/aa(1:kproma,itop,2,im)

      jmax = ibtm_var(jvar) - 1
      DO jk = itop+1,jmax
         jkm1 = jk - 1
         znum(1:kproma) =  bb(1:kproma,jk  ,jvar)                     &
                        & -bb(1:kproma,jkm1,jvar)*aa(1:kproma,jk,1,im)
         bb(1:kproma,jk,jvar) = znum(1:kproma)/aa(1:kproma,jk,2,im)
      ENDDO
    ENDDO !jvar: variable loop

    ! 2. Bottom level for all variables except u, v, dry static energy
    !    and moisture. After this step the array bb contains the
    !    solution of the linear system.

    DO jvar = 1,nvar_vdiff

      IF (jvar==iu.OR.jvar==iv.OR.jvar==ih.OR.jvar==iqv ) THEN
         CYCLE
      ELSE

        im   = matrix_idx(jvar)  ! Index of coefficient matrix
        jk   = ibtm_var(jvar)    ! Bottom level index
        jkm1 = jk - 1
        zden(1:kproma) =  aa(1:kproma,jk,2,im)                      &
                       & -aa(1:kproma,jk,1,im)*aa(1:kproma,jkm1,3,im)
        znum(1:kproma) =  bb(1:kproma,jk,jvar)                      &
                       & -aa(1:kproma,jk,1,im)*bb(1:kproma,jkm1,jvar)
        bb(1:kproma,jk,jvar) = znum(1:kproma)/zden(1:kproma)

      END IF
    ENDDO !jvar: variable loop

    ! Note that for TKE and the variance of theta_v, klev-1 is the lowest
    ! level above surface. Now set boundary condition for the variance 
    ! of theta_v.

    bb(1:kproma,klev,ithv) = bb(1:kproma,klevm1,ithv)

  END SUBROUTINE rhs_elim
  !--------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------
  !>
  !!
  !! Prepare the Richtmyer-Morton coeffcients for dry static energy and 
  !! moisture, to be used by the surface models (ocean, sea-ice, land).
  !!
  SUBROUTINE matrix_to_richtmyer_coeff( kproma, kbdim, klev, ksfc_type, idx_lnd, &! in
                                      & aa, bb,                                  &! in
                                      & aa_btm, bb_btm,                          &! inout
                                      & pen_h, pfn_h, pen_qv, pfn_qv,            &! out
                                      & pcair,                                   &! in
                                      & pcsat)                                    ! in

    INTEGER,INTENT(IN)     :: kproma, kbdim, klev, ksfc_type, idx_lnd
    REAL(wp),INTENT(IN)    :: aa    (kbdim,klev,3,imh:imqv)
    REAL(wp),INTENT(IN)    :: bb    (kbdim,klev,ih:iqv)
    REAL(wp),INTENT(INOUT) :: aa_btm(kbdim,3,ksfc_type,imh:imqv)
    REAL(wp),INTENT(INOUT) :: bb_btm(kbdim,ksfc_type,ih:iqv)

    REAL(wp),INTENT(INOUT) :: pen_h (kbdim,ksfc_type)  ! OUT
    REAL(wp),INTENT(INOUT) :: pfn_h (kbdim,ksfc_type)  ! OUT
    REAL(wp),INTENT(INOUT) :: pen_qv(kbdim,ksfc_type)  ! OUT
    REAL(wp),INTENT(INOUT) :: pfn_qv(kbdim,ksfc_type)  ! OUT

    REAL(wp),OPTIONAL,INTENT(IN)    :: pcair(kbdim)
    REAL(wp),OPTIONAL,INTENT(IN)    :: pcsat(kbdim)

    INTEGER  :: jsfc, klevm1

    klevm1 = klev - 1

    !---------------------------------------------------------
    ! Matrix setup and bottom level elimination for moisture
    !---------------------------------------------------------
    ! Evapotranspiration has to be considered over land 

    IF (phy_config%ljsbach .AND. idx_lnd<=ksfc_type) THEN

      jsfc = idx_lnd

      aa_btm(1:kproma,2,jsfc,imqv) =           1._wp - aa_btm(1:kproma,1,jsfc,imqv) &
                                   & - pcair(1:kproma)*aa_btm(1:kproma,3,jsfc,imqv)
      aa_btm(1:kproma,3,jsfc,imqv) =   pcsat(1:kproma)*aa_btm(1:kproma,3,jsfc,imqv)

    END IF ! ljsbach

    ! Bottom level elimination for all surface types

    DO jsfc = 1,ksfc_type

      aa_btm(1:kproma,2,jsfc,imqv) =  aa_btm(1:kproma,2,jsfc,imqv)  &
                                   & -aa_btm(1:kproma,1,jsfc,imqv)  &
                                   & *aa    (1:kproma,klevm1,3,imqv)

      aa_btm(1:kproma,3,jsfc,imqv) =  aa_btm(1:kproma,3,jsfc,imqv)  &
                                   & /aa_btm(1:kproma,2,jsfc,imqv)

      bb_btm(1:kproma,jsfc,iqv)    = (bb_btm(1:kproma,jsfc,iqv)    &          
                                   & -aa_btm(1:kproma,1,jsfc,imqv) &
                                   & *bb    (1:kproma,klevm1,iqv) )&
                                   & /aa_btm(1:kproma,2,jsfc,imqv)
    END DO

    !---------------------------------------------------------
    ! Bottom level elimination for dry static energy
    !---------------------------------------------------------
    DO jsfc = 1,ksfc_type

      aa_btm(1:kproma,2,jsfc,imh) =  aa_btm(1:kproma,2,jsfc,imh) &
                                  & -aa_btm(1:kproma,1,jsfc,imh) &
                                  & *aa    (1:kproma,klevm1,3,imh)

      aa_btm(1:kproma,3,jsfc,imh) =  aa_btm(1:kproma,3,jsfc,imh) &
                                  & /aa_btm(1:kproma,2,jsfc,imh)

      bb_btm(1:kproma,jsfc,ih)    = (bb_btm(1:kproma,jsfc,ih)    &          
                                  & -aa_btm(1:kproma,1,jsfc,imh) &
                                  & *bb    (1:kproma,klevm1,ih) )&
                                  & /aa_btm(1:kproma,2,jsfc,imh)
    END DO

    !---------------------------------------------------------
    ! Convert matrix entries to Richtmyer-Morton coefficients
    !---------------------------------------------------------

    pen_h (1:kproma,1:ksfc_type) = -aa_btm(1:kproma,3,1:ksfc_type,imh)
    pen_qv(1:kproma,1:ksfc_type) = -aa_btm(1:kproma,3,1:ksfc_type,imqv)

    pfn_h (1:kproma,1:ksfc_type) =  bb_btm(1:kproma,1:ksfc_type,ih )*tpfac1
    pfn_qv(1:kproma,1:ksfc_type) =  bb_btm(1:kproma,1:ksfc_type,iqv)*tpfac1

  END SUBROUTINE matrix_to_richtmyer_coeff
  !--------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------
  !>
  !!
  !! Do Back-substitution to get the solution of the linear system.
  !!
  !! Translation for developers who prefer to think in terms of 
  !! the Richtmyer-Morthon formula and are familiar with the paper by
  !! Polcher et al (1998): on entry bb contains the solution 
  !! at the bottom level and the coeff B on upper levels;
  !! On exit it becomes the solution of the linear system.
  !! aa(:,:,3,:) used here corresponds to -A in the Appendix of 
  !! Polcher et al (1998).
  !! Note that VDIFF uses the implicit time stepping as in IFS
  !! in contrast to Polcher et al (1998). Thus the solution is 
  !! not yet the new value at time step t+dt. 
  !!
  SUBROUTINE rhs_bksub( kproma, kbdim, itop, klev, aa, bb )

    INTEGER, INTENT(IN)   :: kproma, kbdim, itop, klev
    REAL(wp),INTENT(IN)   :: aa(kbdim,klev,3,nmatrix)
    REAL(wp),INTENT(INOUT):: bb(kbdim,klev,nvar_vdiff)

    INTEGER  :: jvar, im, jk, jkp1

    DO jvar = 1,nvar_vdiff
      im = matrix_idx(jvar)
      DO jk = ibtm_var(jvar)-1,itop,-1
         jkp1 = jk + 1
         bb(1:kproma,jk,jvar) =  bb(1:kproma,jk  ,jvar) &
                              & -bb(1:kproma,jkp1,jvar) &
                              & *aa(1:kproma,jk  ,3,im)
      ENDDO
    ENDDO

  END SUBROUTINE rhs_bksub
  !-------------
  !>
  !!
  SUBROUTINE vdiff_tendencies( kproma, kbdim, itop, klev, klevm1,          &! in
                             & ktrac, ksfc_type, idx_wtr,                  &! in
                             & pdtime,                                     &! in
                             & pum1, pvm1, ptm1,                           &! in
                             & pmair, pmdry,                               &! in
                             & pqm1, pxlm1, pxim1, pxtm1,                  &! in
                             & pgeom1, pcptgz,                             &! in
                             & pztkevn, pzthvvar,                          &! in
                             & pcfm_tile, pfrc, bb,                        &! in
                             & pkedisp,                                    &! out
                             & pxvar, pz0m_tile,                           &! inout
                             & pute_vdf, pvte_vdf, pq_vdf,                 &! out
                             & pqte_vdf, pxlte_vdf, pxite_vdf, pxtte_vdf,  &! out
                             & pz0m, ptke, pthvvar,                        &! out
                             & psh_vdiff,pqv_vdiff                         )! out

    INTEGER, INTENT(IN) :: kproma, kbdim, itop, klev, klevm1, ktrac !!$, klevp1
    INTEGER, INTENT(IN) :: ksfc_type, idx_wtr
    REAL(wp),INTENT(IN) :: pdtime

    REAL(wp),INTENT(IN)  :: pum1   (kbdim,klev)
    REAL(wp),INTENT(IN)  :: pvm1   (kbdim,klev)
    REAL(wp),INTENT(IN)  :: ptm1   (kbdim,klev)
    REAL(wp),INTENT(IN)  :: pmair  (kbdim,klev)  !< moist air mass [kg/m2]
    REAL(wp),INTENT(IN)  :: pmdry  (kbdim,klev)  !< dry   air mass [kg/m2]
    REAL(wp),INTENT(IN)  :: pqm1   (kbdim,klev)
    REAL(wp),INTENT(IN)  :: pxlm1  (kbdim,klev)
    REAL(wp),INTENT(IN)  :: pxim1  (kbdim,klev)
    REAL(wp),INTENT(IN)  :: pxtm1  (kbdim,klev,ktrac)
    REAL(wp),INTENT(IN)  :: pgeom1 (kbdim,klev)
    REAL(wp),INTENT(IN)  :: pcptgz (kbdim,klev)
    REAL(wp),INTENT(IN)  :: pztkevn (kbdim,klev)
    REAL(wp),INTENT(IN)  :: pzthvvar(kbdim,klev)
    REAL(wp),INTENT(IN)  :: pcfm_tile     (kbdim,ksfc_type)
    REAL(wp),INTENT(IN)  :: pfrc          (kbdim,ksfc_type)
    REAL(wp),INTENT(IN)  :: bb            (kbdim,klev,nvar_vdiff)

    REAL(wp),INTENT(OUT) :: pkedisp(kbdim) !< vertically integrated dissipation
                                           !  of kinetic energy [W/m2]

    REAL(wp),INTENT(INOUT) :: pxvar    (kbdim,klev)
    REAL(wp),INTENT(INOUT) :: pz0m_tile(kbdim,ksfc_type)

    REAL(wp),INTENT(OUT) :: pute_vdf (kbdim,klev)
    REAL(wp),INTENT(OUT) :: pvte_vdf (kbdim,klev)
    REAL(wp),INTENT(OUT) :: pq_vdf   (kbdim,klev)
    REAL(wp),INTENT(OUT) :: pqte_vdf (kbdim,klev)
    REAL(wp),INTENT(OUT) :: pxlte_vdf(kbdim,klev)
    REAL(wp),INTENT(OUT) :: pxite_vdf(kbdim,klev)
    REAL(wp),INTENT(OUT) :: pxtte_vdf(kbdim,klev,ktrac)

    REAL(wp),INTENT(OUT) :: pz0m     (kbdim)
    REAL(wp),INTENT(OUT) :: ptke     (kbdim,klev)
    REAL(wp),INTENT(OUT) :: pthvvar  (kbdim,klev)
    REAL(wp),INTENT(OUT) :: psh_vdiff(kbdim)
    REAL(wp),INTENT(OUT) :: pqv_vdiff(kbdim)

    REAL(wp) :: ztest, zrdt
    REAL(wp) :: zunew, zvnew, zqnew, zsnew, zhnew
    REAL(wp) :: zcp
    REAL(wp) :: zdis  (kbdim,klev)

    INTEGER  :: jk, jl, jt, irhs, jsfc


    zrdt   = 1._wp/pdtime

    pute_vdf (:,:)   = 0._wp
    pvte_vdf (:,:)   = 0._wp
    pq_vdf   (:,:)   = 0._wp
    pqte_vdf (:,:)   = 0._wp
    pxlte_vdf(:,:)   = 0._wp
    pxite_vdf(:,:)   = 0._wp
    pxtte_vdf(:,:,:) = 0._wp

    ptke     (:,:)   = 0._wp
    pthvvar  (:,:)   = 0._wp
    pz0m     (:)     = 0._wp

    !-------------------------------------------------------------------
    ! Compute TKE at the new time step.
    !-------------------------------------------------------------------

    DO jk = itop,klevm1
      ztest = 0._wp
      DO jl = 1,kproma
        ptke(jl,jk) = bb(jl,jk,itke) + tpfac3*pztkevn(jl,jk)
        ztest = ztest+MERGE(1._wp,0._wp,ptke(jl,jk)<0._wp)
      END DO
      IF(ztest.NE.0._wp) THEN
        WRITE(message_text,'(a,I4,2E15.5)') 'level, MIN TKE components = ',&
             & jk, MINVAL(bb(:,jk,itke)),MINVAL(pztkevn(:,jk))
        CALL message('', TRIM(message_text))
        CALL finish('vdiff_tendencies','TKE IS NEGATIVE')
      ENDIF
    END DO
    ptke(1:kproma,klev) = pztkevn(1:kproma,klev)


    !-------------------------------------------------------------
    ! Variance of virtual potential temperature
    !-------------------------------------------------------------
    DO jk = itop,klev
      DO jl = 1,kproma
        pthvvar(jl,jk) = bb(jl,jk,ithv) + tpfac3*pzthvvar(jl,jk)
        pthvvar(jl,jk) = MAX(tke_min,pthvvar(jl,jk))
      END DO
    END DO

    !-------------------------------------------------------------
    ! Tendency of velocity; kinetic energy dissipation
    !-------------------------------------------------------------
    pkedisp(:) = 0._wp   ! initilize the vertical integral

    DO jk = itop,klev
      DO jl = 1,kproma
        pute_vdf(jl,jk) = (bb(jl,jk,iu)-tpfac2*pum1(jl,jk))*zrdt
        pvte_vdf(jl,jk) = (bb(jl,jk,iv)-tpfac2*pvm1(jl,jk))*zrdt

        zunew = bb(jl,jk,iu) + tpfac3*pum1(jl,jk)
        zvnew = bb(jl,jk,iv) + tpfac3*pvm1(jl,jk)

        zdis(jl,jk) = 0.5_wp*( pum1(jl,jk)**2 - zunew**2 &
                    &         +pvm1(jl,jk)**2 - zvnew**2 )
        pkedisp(jl)  = pkedisp(jl) + zdis(jl,jk)*pmair(jl,jk)*zrdt
      END DO
    END DO

    !-------------------------------------------------------------
    ! Tendency of T and qv, ql, qi; xvar at the new time step
    !-------------------------------------------------------------
    DO jk=itop,klev
      DO jl=1,kproma

        zqnew = bb(jl,jk,iqv) + tpfac3*pqm1(jl,jk)
        pqte_vdf(jl,jk) = (zqnew-pqm1(jl,jk))*zrdt


        ! The computation of the new temperature must be consistent with the computation of
        ! the static energy pcptgz in the subroutine mo_turbulence_diag:atm_exchange_coeff.
        ! The same specific heat must be used.
        !
        zsnew = bb(jl,jk,ih) + tpfac3*pcptgz(jl,jk)
        zhnew = (zsnew + zdis(jl,jk) - pgeom1(jl,jk))
        zcp   = cpd+(cpv-cpd)*pqm1(jl,jk) ! cp of moist air
        !
        ! Now derive the heating for constant pressure conditions
        ! as needed for provisional updating in the physics.
        ! 
        pq_vdf(jl,jk)   = (zhnew - ptm1(jl,jk)*zcp)*zrdt*pmair(jl,jk)

        pxlte_vdf(jl,jk) = (bb(jl,jk,ixl) - tpfac2*pxlm1(jl,jk))*zrdt
        pxite_vdf(jl,jk) = (bb(jl,jk,ixi) - tpfac2*pxim1(jl,jk))*zrdt

        pxvar(jl,jk) = bb(jl,jk,ixv) + tpfac3*pxvar(jl,jk)
      END DO
    END DO

    IF ( get_lebudget() ) THEN
      psh_vdiff(:) = 0._wp
      pqv_vdiff(:) = 0._wp
      DO jk=itop,klev
        ! compute heat budget diagnostic
        psh_vdiff(1:kproma) = psh_vdiff(1:kproma) + pmdry(1:kproma,jk) * &
        & (bb(1:kproma,jk,ih)  + (tpfac3 - 1._wp)*pcptgz(1:kproma,jk)) * zrdt
        ! compute moisture budget diagnostic
        ! ? zdis appears to be dissipation, probably we don't need this for qv??
        pqv_vdiff(1:kproma) = pqv_vdiff(1:kproma) + pmdry(1:kproma,jk)* &
        & (bb(1:kproma,jk,iqv) + (tpfac3 - 1._wp)*pqm1(1:kproma,jk)) * zrdt
      END DO
    END IF
    !-------------------------------------------------------------
    ! Tendency of tracers
    !-------------------------------------------------------------
!   IF (trlist% anyvdiff /= 0) THEN   ! ECHAM
!     DO 577 jt=1,trlist% ntrac       ! ECHAM
        DO jt = 1,ktrac
          irhs = itrc_start + jt - 1
!         IF (trlist% ti(jt)% nvdiff /= 1) CYCLE  ! ECHAM
          DO jk = itop,klev
            DO jl = 1,kproma
              pxtte_vdf(jl,jk,jt) = (bb(jl,jk,irhs)-tpfac2*pxtm1(jl,jk,jt))*zrdt
            ENDDO
          ENDDO
        ENDDO
!577  ENDDO
!     END IF

    !----------------------------------------------------------------------------
    ! Update roughness height over open water, then update the grid-box mean
    !----------------------------------------------------------------------------
    IF (idx_wtr<=ksfc_type) THEN  ! water surface exists in the simulation
      DO jl = 1,kproma
        IF(pfrc(jl,idx_wtr).GT.0._wp) THEN
          pz0m_tile(jl,idx_wtr) = tpfac1*SQRT( bb(jl,klev,iu)**2+bb(jl,klev,iv)**2 ) &
                                & *pcfm_tile(jl,idx_wtr)*cchar*rgrav
          pz0m_tile(jl,idx_wtr) = MAX(z0m_min,pz0m_tile(jl,idx_wtr))
        ELSE
          pz0m_tile(jl,idx_wtr) = cdimissval
        ENDIF
      ENDDO
    ENDIF

    ! Compute grid-box mean 

    pz0m(:) = 0._wp
    DO jsfc = 1,ksfc_type
       pz0m(1:kproma) = pz0m(1:kproma) + pfrc(1:kproma,jsfc)*pz0m_tile(1:kproma,jsfc)
    ENDDO

  END SUBROUTINE vdiff_tendencies
  !-------------



END MODULE mo_vdiff_solver
