!OPTION! -cont -msg o
! this command should fix the problem of copying arrays in a subroutine call
!
! This module is the interface between nwp_nh_interface and the
! turbulence parameterisations:
! inwp_turb == icosmo == turbulence scheme by from COSMO (Raschendorfer)
! inwp_turb == iedmf  == eddy-diffusivity mass-flux DUALM scheme from IFS (Koehler, Neggers)
! inwp_turb == igme   == turbulence scheme imported from the GME
! This module handles the atmospheric turbulent above the surface layer diffusion.
!
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

MODULE mo_nwp_turbdiff_interface

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: message, message_text, finish, warning
  USE mo_model_domain,           ONLY: t_patch
  USE mo_impl_constants,         ONLY: min_rlcell_int, igme, icosmo, iedmf, &
    &                                  max_ntracer
  USE mo_impl_constants_grf,     ONLY: grf_bdywidth_c
  USE mo_loopindices,            ONLY: get_indices_c
  USE mo_physical_constants,     ONLY: alv, grav, vtmpc1, rd, cpd, cvd
  USE mo_ext_data_types,         ONLY: t_external_data
  USE mo_nonhydro_types,         ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_nwp_phy_types,          ONLY: t_nwp_phy_diag, t_nwp_phy_tend
  USE mo_nwp_phy_state,          ONLY: phy_params 
  USE mo_nwp_lnd_types,          ONLY: t_lnd_prog, t_wtr_prog, t_lnd_diag
  USE mo_parallel_config,        ONLY: nproma
  USE mo_run_config,             ONLY: msg_level, iqv, iqc, iqi, iqnc, iqni, iqtke, &
    &                                  iqs, iqns, iqtvar, lart, ltestcase
  USE mo_atm_phy_nwp_config,     ONLY: atm_phy_nwp_config
  USE mo_nonhydrostatic_config,  ONLY: kstart_moist, kstart_tracer
  USE turb_data,                 ONLY: get_turbdiff_param, lsflcnd, modvar, ndim
  USE turb_diffusion,            ONLY: turbdiff
  USE turb_vertdiff,             ONLY: vertdiff
  USE mo_gme_turbdiff,           ONLY: partura, progimp_turb

  USE mo_art_config,             ONLY: art_config
  USE mo_advection_config,       ONLY: advection_config
  USE mo_turbdiff_config,        ONLY: turbdiff_config
#ifdef __ICON_ART
  USE mo_art_turbdiff_interface, ONLY: art_turbdiff_interface
#endif
#ifndef __NO_ICON_EDMF__
  USE mo_edmf_param,             ONLY: ntiles_edmf
  USE mo_vdfouter,               ONLY: vdfouter
#endif
  USE mo_lnd_nwp_config,         ONLY: nlev_soil, nlev_snow, ntiles_total, ntiles_water
  USE mo_grid_config,            ONLY: l_scm_mode
  USE mo_scm_nml,                ONLY: scm_sfc_mom, scm_sfc_temp ,scm_sfc_qv
  USE mo_nh_torus_exp,           ONLY: set_scm_bnd
  USE mo_timer
  USE mo_run_config,             ONLY: timers_level
  USE mo_fortran_tools,          ONLY: assert_acc_device_only
  IMPLICIT NONE

  PRIVATE



  PUBLIC  ::  nwp_turbdiff

CONTAINS
  !!
  !!-------------------------------------------------------------------------
  !!
SUBROUTINE nwp_turbdiff  ( tcall_turb_jg,                     & !>in
                          & p_patch,                          & !>in
                          & p_metrics,                        & !>in
                          & ext_data,                         & !>in
                          & p_prog,                           & !>in
                          & p_prog_now_rcf,                   & !>in
                          & p_prog_rcf,                       & !>inout
                          & p_diag ,                          & !>inout
                          & prm_diag, prm_nwp_tend,           & !>inout
                          & wtr_prog_now,                     & !>in 
                          & lnd_prog_now,                     & !>in 
                          & lnd_diag,                         & !>in
                          & lacc                              ) !>in


  TYPE(t_patch),        TARGET,INTENT(in)   :: p_patch        !!<grid/patch info.
  TYPE(t_external_data),       INTENT(in)   :: ext_data        !< external data
  TYPE(t_nh_metrics)          ,INTENT(in)   :: p_metrics
  TYPE(t_nh_prog),      TARGET,INTENT(inout):: p_prog          !<the prog vars
  TYPE(t_nh_prog),      TARGET,INTENT(in)   :: p_prog_now_rcf  !<progs with red.
  TYPE(t_nh_prog),      TARGET,INTENT(inout):: p_prog_rcf      !<call freq
  TYPE(t_nh_diag),      TARGET,INTENT(inout):: p_diag          !<the diag vars
  TYPE(t_nwp_phy_diag),        INTENT(inout):: prm_diag        !< atm phys vars
  TYPE(t_nwp_phy_tend), TARGET,INTENT(inout):: prm_nwp_tend    !< atm tend vars
  TYPE(t_wtr_prog),            INTENT(in)   :: wtr_prog_now    !< prog vars for wtr
  TYPE(t_lnd_prog),            INTENT(inout):: lnd_prog_now    !< prog vars for sfc
  TYPE(t_lnd_diag),            INTENT(inout):: lnd_diag        !< diag vars for sfc
  REAL(wp),                    INTENT(in)   :: tcall_turb_jg   !< time interval for 
                                                               !< turbulence
  LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc


  ! Local array bounds

  INTEGER :: rl_start, rl_end
  INTEGER :: i_startblk, i_endblk    !> blocks
  INTEGER :: i_startidx, i_endidx    !< slices

  ! Local scalars:

  INTEGER :: jc,jk,jb,jg      !loop indices

  ! local variables for turbdiff

  REAL(wp), PARAMETER :: cpd_o_cvd = cpd/cvd

  INTEGER :: nzprv=1

  INTEGER  :: nlev, nlevp1, nlevcm                  !< number of full, half and canopy levels

  REAL(wp) :: tke_inc_ic(nproma)                    !< TKE increment at half levels

  REAL(wp) :: l_hori(nproma)                        !< horizontal length scale

  REAL(wp) :: zvari(nproma,p_patch%nlevp1,ndim), &  !< to pass values from turbdiff to vertdiff
              zrhon(nproma,p_patch%nlevp1)

  REAL(wp) :: z_tvs(nproma,p_patch%nlevp1,1)        !< aux turbulence velocity scale [m/s]
  REAL(wp) :: tempv_sfc(nproma)                     !< surface virtual temperature [K]
  REAL(wp) :: rho_sfc(nproma)                       !< surface density [Pa]

  ! type structure to hand over additional tracers to turbdiff
  TYPE(modvar) :: ptr(max_ntracer)

  INTEGER :: nturb                               !< Loop counter
  INTEGER :: ncloud_offset                       !< offset for ptr-indexing in ART
                                                 !< interface due to additionally
                                                 !< diffused cloud fields
  INTEGER, ALLOCATABLE :: idx_nturb_tracer(:)    !< indices of the turbulent tracers in the prognostic list

  LOGICAL  :: ltwomoment                    !< using 2mom microphysics?
  REAL(wp), TARGET      :: & 
    &  ddt_turb_qnc(nproma,p_patch%nlev), & !< tendency field for qnc
    &  ddt_turb_qni(nproma,p_patch%nlev), & !< tendendy field for qni
    &  ddt_turb_qs (nproma,p_patch%nlev), & !< tendency field for qs 
    &  ddt_turb_qns(nproma,p_patch%nlev)    !< tendendy field for qns

  ! EDMF variables:
  ! (attention: if no block index - p_patch%nblks_c - variables need to be declared private)

  INTEGER  :: jt
  INTEGER, PARAMETER :: itrac_vdf = 0
  INTEGER  :: khpbln(nproma)      , kvartop(nproma)     , kpbltype(nproma)
  REAL(wp) :: pdifts(nproma,p_patch%nlev+1) , pdiftq(nproma,p_patch%nlev+1)   , &
    &         pdiftl(nproma,p_patch%nlev+1) , pdifti(nproma,p_patch%nlev+1)   , &
    &         pstrtu(nproma,p_patch%nlev+1) , pstrtv(nproma,p_patch%nlev+1)   , &
    &         pkh   (nproma,p_patch%nlev)   , pkm   (nproma,p_patch%nlev)
  REAL(wp) :: z_omega_p(nproma,p_patch%nlev), zchar(nproma)                   , &
    &         zucurr(nproma)                , zvcurr(nproma)                  , &
    &         zsoteu(nproma,p_patch%nlev)   , zsotev(nproma,p_patch%nlev)     , &
    &         zsobeta(nproma,p_patch%nlev)  , zz0m(nproma),  zz0h(nproma)     , &
    &         zae(nproma,p_patch%nlev)      , ztskrad(nproma)                 , &
    &         zsigflt(nproma)               ,                                   &
    &         shfl_s_t(nproma,ntiles_total+ntiles_water)                      , &
    &         evap_s_t(nproma,ntiles_total+ntiles_water)                      , &
    &         tskin_t (nproma,ntiles_total+ntiles_water)                      , &
    &         ustr_s_t(nproma,ntiles_total+ntiles_water)                      , &
    &         vstr_s_t(nproma,ntiles_total+ntiles_water)

  REAL(wp) :: ut_sso(nproma, p_patch%nlev), vt_sso(nproma, p_patch%nlev)
  REAL(wp),DIMENSION(nproma,p_patch%nlev+1) :: tet_flux, vap_flux, liq_flux
 
  INTEGER, SAVE :: nstep_turb = 0

!--------------------------------------------------------------

  IF (msg_level >= 15) CALL message('mo_nwp_turbdiff:', 'turbulence')

  CALL assert_acc_device_only("mo_nwp_turbdiff", lacc)

  jg = p_patch%id

  ! number of vertical levels
  nlev   = p_patch%nlev
  nlevp1 = p_patch%nlevp1

  nlevcm = nlevp1
   
  ! local variables related to the blocking
  !
  ! exclude boundary interpolation zone of nested domains
  rl_start = grf_bdywidth_c+1
  rl_end   = min_rlcell_int

  i_startblk = p_patch%cells%start_block(rl_start)
  i_endblk   = p_patch%cells%end_block(rl_end)

  
  IF ( atm_phy_nwp_config(jg)%inwp_turb == icosmo ) THEN
     CALL get_turbdiff_param(jg)
  ENDIF

  IF ( lart .AND. art_config(jg)%nturb_tracer > 0 ) THEN
     ALLOCATE(idx_nturb_tracer(art_config(jg)%nturb_tracer))
  END IF

  ! logical for SB two-moment scheme
  ltwomoment = atm_phy_nwp_config(jg)%l2moment

  !$ACC DATA CREATE(khpbln, kvartop, kpbltype, pdifts, pdiftq, pdiftl, pdifti, pstrtu, pstrtv, pkh, pkm, z_omega_p) &
  !$ACC   CREATE(zchar, zucurr, zvcurr, zsoteu, zsotev, zsobeta, zz0m, zz0h, zae, ztskrad, zsigflt, shfl_s_t) &
  !$ACC   CREATE(evap_s_t, tskin_t, ustr_s_t, vstr_s_t, ddt_turb_qnc, ddt_turb_qni, ddt_turb_qs, ddt_turb_qns) &
  !$ACC   CREATE(tke_inc_ic, l_hori, zvari, zrhon, z_tvs, tempv_sfc, rho_sfc, ut_sso, vt_sso)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,tke_inc_ic,z_tvs, &
!$OMP            ncloud_offset,ptr,nzprv,l_hori,zvari,zrhon,                               &
!$OMP            jt       , khpbln  , kvartop , kpbltype,                                  &
!$OMP            pdifts   , pdiftq  , pdiftl  , pdifti  , pstrtu  , pstrtv , pkh , pkm ,   &
!$OMP            z_omega_p, zchar   , zucurr  , zvcurr  , zsoteu  , zsotev , zsobeta   ,   &
!$OMP            zz0m     , zz0h    , zae     , ztskrad , zsigflt ,                        &
!$OMP            shfl_s_t , evap_s_t, tskin_t , ustr_s_t, vstr_s_t, rho_sfc, tempv_sfc,    &
!$OMP            ddt_turb_qnc, ddt_turb_qni, ddt_turb_qs, ddt_turb_qns, ut_sso, vt_sso)  ICON_OMP_GUIDED_SCHEDULE

  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
      & i_startidx, i_endidx, rl_start, rl_end)

   !-------------------------------------------------------------------------
   !<  turbulent diffusion
   !-------------------------------------------------------------------------

   !<  NOTE: since  turbulence is a fast process it is
   !!        allowed to do a sequential updating except for wind speed
   !!        because back-and-forth interpolation would cause too large errors
   !!  (GZ, 2011-08-29): Nevertheless, tendency fields are now passed to turbdiff
   !!        to have them available for extended diagnostic output


    SELECT CASE (atm_phy_nwp_config(jg)%inwp_turb)

    CASE(icosmo)

!-------------------------------------------------------------------------
!< COSMO turbulence scheme by M. Raschendorfer  
!-------------------------------------------------------------------------


      !
      ! convert TKE to the turbulence velocity scale SQRT(2*TKE) as required by turbdiff
      ! INPUT to turbdiff is timestep now
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) ASYNC(1) DEFAULT(PRESENT)
      DO jk=1, nlevp1
        DO jc=i_startidx, i_endidx
          z_tvs(jc,jk,1) = SQRT(2._wp* (p_prog_now_rcf%tke(jc,jk,jb))) 
         ENDDO
      ENDDO
      !$ACC END PARALLEL LOOP



      IF (advection_config(jg)%iadv_tke > 0) THEN
        ! Interpolate advective tvs tendency from full levels to half levels
        ! Note that both the advective TKE tendency and ddt_tke actually carry time tendencies
        ! of tvs; attempts to horizontally advect TKE failed because of numerical instability
        !
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
        !$ACC LOOP GANG
        DO jk=2, nlev
          !$ACC LOOP VECTOR
          DO jc=i_startidx, i_endidx

            tke_inc_ic(jc) = p_metrics%wgtfac_c(jc,jk,jb) * p_diag%ddt_tracer_adv(jc,jk,jb,iqtke) &
              &             + (1._wp - p_metrics%wgtfac_c(jc,jk,jb)) &
              &             * p_diag%ddt_tracer_adv(jc,jk-1,jb,iqtke)

            ! add advective TKE (actually tvs) tendency to ddt_tke, which is provided to turbdiff as input
            prm_nwp_tend%ddt_tke(jc,jk,jb) = prm_nwp_tend%ddt_tke(jc,jk,jb) + tke_inc_ic(jc)

          ENDDO  ! jc
        ENDDO  ! jk
        !$ACC END PARALLEL
        !
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
        !$ACC LOOP GANG VECTOR
        DO jc=i_startidx, i_endidx

          ! zero gradient assumption for TKE increment at model bottom (top level not needed)
          prm_nwp_tend%ddt_tke(jc,nlevp1,jb) = prm_nwp_tend%ddt_tke(jc,nlevp1,jb) + p_diag%ddt_tracer_adv(jc,nlev,jb,iqtke)

        ENDDO  ! jc
        !$ACC END PARALLEL
      ENDIF



      !KF tendencies  have to be set to zero
      !GZ: this should be replaced by an appropriate switch in turbdiff
      !$ACC KERNELS ASYNC(1) DEFAULT(PRESENT)
      prm_nwp_tend%ddt_u_turb(:,:,jb) = 0._wp
      prm_nwp_tend%ddt_v_turb(:,:,jb) = 0._wp
      prm_nwp_tend%ddt_temp_turb(:,:,jb) = 0._wp
      prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqv) = 0._wp
      prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqc) = 0._wp
      !$ACC END KERNELS

      ! offset for ptr-indexing in ART-Interface
      ncloud_offset = 0

      IF (ltwomoment) THEN
        ! register cloud droplet number for turbulent diffusion
        ncloud_offset = ncloud_offset+1
        !$ACC KERNELS ASYNC(1) DEFAULT(PRESENT)
        ddt_turb_qnc(:,:) = 0.0_wp
        !$ACC END KERNELS
        ptr(ncloud_offset)%av     => p_prog_rcf%tracer(:,:,jb,iqnc)
        ptr(ncloud_offset)%at     => ddt_turb_qnc(:,:)
        ptr(ncloud_offset)%sv     => NULL()
        ptr(ncloud_offset)%kstart =  kstart_moist(jg)
      ENDIF ! ltwomoment

      IF (turbdiff_config(jg)%ldiff_qi) THEN
        ! register cloud ice for turbulent diffusion
        ncloud_offset = ncloud_offset + 1
        !$ACC KERNELS ASYNC(1) DEFAULT(PRESENT)
        prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqi) = 0.0_wp
        !$ACC END KERNELS
        ptr(ncloud_offset)%av     => p_prog_rcf%tracer(:,:,jb,iqi)
        ptr(ncloud_offset)%at     => prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqi)
        ptr(ncloud_offset)%sv     => NULL()
        ptr(ncloud_offset)%kstart =  kstart_moist(jg)
        IF (ltwomoment) THEN
          ! register cloud ice number for turbulent diffusion
          ncloud_offset = ncloud_offset + 1
          !$ACC KERNELS ASYNC(1) DEFAULT(PRESENT)
          ddt_turb_qni(:,:) = 0.0_wp
          !$ACC END KERNELS
          ptr(ncloud_offset)%av     => p_prog_rcf%tracer(:,:,jb,iqni)
          ptr(ncloud_offset)%at     => ddt_turb_qni(:,:)
          ptr(ncloud_offset)%sv     => NULL()
          ptr(ncloud_offset)%kstart =  kstart_moist(jg)
        ENDIF ! ltwomoment
      ENDIF ! turbdiff_config(jg)%ldiff_qi

      IF (turbdiff_config(jg)%ldiff_qs) THEN
        ! register snow mass for turbulent diffusion
        ncloud_offset = ncloud_offset + 1
        !$ACC KERNELS ASYNC(1) DEFAULT(PRESENT)
        ddt_turb_qs (:,:) = 0.0_wp
        !$ACC END KERNELS
        ptr(ncloud_offset)%av     => p_prog_rcf%tracer(:,:,jb,iqs )
        ptr(ncloud_offset)%at     => ddt_turb_qs (:,:)
        ptr(ncloud_offset)%sv     => NULL()
        ptr(ncloud_offset)%kstart =  kstart_moist(jg)
        IF (ltwomoment) THEN
          ! register snow number for turbulent diffusion
          ncloud_offset = ncloud_offset + 1
          !$ACC KERNELS ASYNC(1) DEFAULT(PRESENT)
          ddt_turb_qns(:,:) = 0.0_wp
          !$ACC END KERNELS
          ptr(ncloud_offset)%av     => p_prog_rcf%tracer(:,:,jb,iqns)
          ptr(ncloud_offset)%at     => ddt_turb_qns(:,:)
          ptr(ncloud_offset)%sv     => NULL()
          ptr(ncloud_offset)%kstart =  kstart_moist(jg)
        ENDIF ! ltwomoment
      ENDIF ! turbdiff_config(jg)%ldiff_qs

#ifdef __ICON_ART
      IF ( lart .AND. art_config(jg)%nturb_tracer > 0 ) THEN
         CALL art_turbdiff_interface( 'setup_ptr', p_patch, p_prog_rcf, prm_nwp_tend,  &
           &                          ncloud_offset=ncloud_offset,                     &
           &                          ptr=ptr(:), dt=tcall_turb_jg,                    &
           &                          p_rho=p_prog%rho(:,:,:),                         &
           &                          p_metrics=p_metrics,                             &
           &                          p_diag=p_diag, prm_diag=prm_diag,                &
           &                          jb=jb, idx_nturb_tracer=idx_nturb_tracer,        &
           &                          lacc=lacc )
         DO nturb = 1, art_config(jg)%nturb_tracer
           ptr(ncloud_offset+nturb)%kstart = kstart_tracer(jg,idx_nturb_tracer(nturb))
         ENDDO
      ENDIF
#endif
  
      IF ( ltestcase .AND. l_scm_mode .AND. &
        &  ((scm_sfc_mom .GE. 2) .OR. (scm_sfc_temp .GE. 2) .OR. (scm_sfc_qv .GE. 2)) ) THEN
        CALL set_scm_bnd( nvec=nproma, ivstart=i_startidx, ivend=i_endidx, &
          & u_s          = p_diag%u(:,nlev,jb),                            & !in
          & v_s          = p_diag%v(:,nlev,jb),                            & !in
          & th_b         = p_diag%temp(:,nlev,jb)/p_prog%exner(:,nlev,jb), & !in
          & qv_b         = p_prog_rcf%tracer(:,nlev,jb,iqv),               & !in
          & pres_sfc     = p_diag%pres_sfc(:,jb),                          & !in
          & dz_bs=p_metrics%z_mc(:,nlev,jb)-p_metrics%z_ifc(:,nlevp1,jb),  & !in
          & z0m=prm_diag%gz0(:,jb)/grav,                                   & !in
          !for noq z0m is assumed to be equal to z0h - GABLS1
          & z0h=prm_diag%gz0(:,jb)/grav,                                   & !in
          & prm_nwp_tend = prm_nwp_tend,                                   & !in 
          & tvm          = prm_diag%tvm(:,jb),                             & !inout
          & tvh          = prm_diag%tvh(:,jb),                             & !inout
          & shfl_s       = prm_diag%shfl_s(:,jb),                          & !out
          & qhfl_s       = prm_diag%qhfl_s(:,jb),                          & !out
          & lhfl_s       = prm_diag%lhfl_s(:,jb),                          & !out
          & umfl_s       = prm_diag%umfl_s(:,jb),                          & !out
          & vmfl_s       = prm_diag%vmfl_s(:,jb),                          & !out
          & qv_s         = lnd_diag%qv_s(:,jb),                            & !out
          & t_g          = lnd_prog_now%t_g(:,jb) )                          !out
      ENDIF

      !should be dependent on location in future!
      !$ACC KERNELS ASYNC(1) DEFAULT(PRESENT)
      l_hori(i_startidx:i_endidx)=phy_params(jg)%mean_charlen
      !$ACC END KERNELS

      nzprv = 1

      !$ACC KERNELS ASYNC(1) DEFAULT(PRESENT)
      ut_sso(:,:)=REAL(prm_nwp_tend%ddt_u_sso(:,:,jb), wp)
      vt_sso(:,:)=REAL(prm_nwp_tend%ddt_v_sso(:,:,jb), wp)
      !$ACC END KERNELS

      IF (timers_level > 9) CALL timer_start(timer_nwp_turbdiff)

      ! turbdiff
      CALL turbdiff( &
        &  iini=0,                         & !atmosph. turbulence and vertical diffusion
        &  ltkeinp=.FALSE.,        & !
        &  lstfnct=.TRUE. ,        & !
        &  l3dturb=.FALSE.,        & !
        &  lrunsso=(atm_phy_nwp_config(jg)%inwp_sso > 0),        & ! running COSMO SSO scheme
        &  lruncnv=(atm_phy_nwp_config(jg)%inwp_convection > 0), & ! running convection
        &  lrunscm=.FALSE.,                                      & ! no single column model
        &  lsfluse=lsflcnd,        & !
!
        &  dt_var=tcall_turb_jg,                                                      & !in
        &  dt_tke=tcall_turb_jg,                                                      & !in
        &  nprv=nzprv, ntur=1, ntim=1,                                                & !in
        &  nvec=nproma, ke=nlev, ke1=nlevp1, kcm=nlevcm, iblock=jb,                   & !in
        &  ivstart=i_startidx, ivend=i_endidx,                                        & !in
!
        &  l_hori=l_hori,                                                             & !in
        &  hhl=p_metrics%z_ifc(:,:,jb),                                               & !in
        &  dp0=p_diag%dpres_mc(:,:,jb),                                               & !in
        &  trop_mask=prm_diag%tropics_mask(:,jb),                                     & !in
        &  innertrop_mask=prm_diag%innertropics_mask(:,jb),                           & !in
        &  gz0=prm_diag%gz0(:,jb),                                                    & !in
        &  l_pat = ext_data%atm%l_pat(:,jb),                                          & !in
        &  t_g=lnd_prog_now%t_g(:,jb),                                                & !in
        &  qv_s=lnd_diag%qv_s(:,jb),                                                  & !in
        &  ps=p_diag%pres_sfc(:,jb),                                                  & !in
        &  u=p_diag%u(:,:,jb),                                                        & !in
        &  v=p_diag%v(:,:,jb),                                                        & !in
        &  w=p_prog%w(:,:,jb),                                                        & !in
        &  t=p_diag%temp(:,:,jb),                                                     & !in
        &  qv=p_prog_rcf%tracer(:,:,jb,iqv),                                          & !in
        &  qc=p_prog_rcf%tracer(:,:,jb,iqc),                                          & !in
        &  prs=p_diag%pres(:,:,jb),                                                   & !in
        &  rhoh=p_prog%rho(:,:,jb),                                                   & !in
        &  rhon=zrhon(:,:),                                                           & !out
        &  epr=p_prog%exner(:,:,jb),                                                  & !in
        &  impl_weight=turbdiff_config(jg)%impl_weight,                               & !in
        &  tvm=prm_diag%tvm(:,jb),                                                    & !inout
        &  tvh=prm_diag%tvh(:,jb),                                                    & !inout
        &  tfm=prm_diag%tfm(:,jb),                                                    & !inout
        &  tfh=prm_diag%tfh(:,jb),                                                    & !inout
        &  tkred_sfc=prm_diag%tkred_sfc(:,jb),                                        & !inout
        &  tkred_sfc_h=prm_diag%tkred_sfc_h(:,jb),                                    & !inout
        &  tke=z_tvs(:,:,:),                                                          & !inout
        &  tkvm=prm_diag%tkvm(:,:,jb),                                                & !inout
        &  tkvh=prm_diag%tkvh(:,:,jb),                                                & !inout
        &  rcld=prm_diag%rcld(:,:,jb),                                                & !inout
        &  hdef2=p_diag%hdef_ic(:,:,jb),                                              & !in
        &  hdiv=p_diag%div_ic(:,:,jb),                                                & !in
        &  dwdx=p_diag%dwdx(:,:,jb),                                                  & !in
        &  dwdy=p_diag%dwdy(:,:,jb),                                                  & !in
        &  tket_conv=prm_nwp_tend%ddt_tke_pconv(:,:,jb),                              & !in
        &  tket_hshr=prm_nwp_tend%ddt_tke_hsh(:,:,jb),                                & !out
        &  u_tens=prm_nwp_tend%ddt_u_turb(:,:,jb),                                    & !inout
        &  v_tens=prm_nwp_tend%ddt_v_turb(:,:,jb),                                    & !inout
        &  t_tens=prm_nwp_tend%ddt_temp_turb(:,:,jb),                                 & !inout
        &  qv_tens=prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqv),                          & !inout
        &  qc_tens=prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqc),                          & !inout
        &  tketens=prm_nwp_tend%ddt_tke(:,:,jb),                                      & !inout
        &  ut_sso=ut_sso(:,:),                                                        & !in
        &  vt_sso=vt_sso(:,:),                                                        & !in
        &  shfl_s=prm_diag%shfl_s(:,jb),                                              & !in
        &  qvfl_s=prm_diag%qhfl_s(:,jb),                                              & !out
        &  zvari=zvari(:,:,:),                                                        & !out
        &  tet_flux=tet_flux,                                      & !out 
        &  vap_flux=vap_flux,                                      & !out
        &  liq_flux=liq_flux,                                      & !out
        &  l_3d_turb_fluxes=atm_phy_nwp_config(jg)%l_3d_turb_fluxes)                    !in

      IF (timers_level > 9) CALL timer_stop(timer_nwp_turbdiff)

      ! vertdiff
      CALL vertdiff( &

        &  iini=0, lturatm=.TRUE. ,        & !atmosph. turbulence and vertical diffusion


        &  itnd=0, lum_dif=.TRUE. ,        & ! with vertical diffusion
        &          lvm_dif=.TRUE. ,        & !
        &          lscadif=.TRUE.,         & !
        &          lsfluse=lsflcnd,        & !and thus (implicitly) neither surface-layer diagnosis
        &          lqvcrst=.FALSE.,        & !nor surface-flux calculation (both in turbtran)
        &          lrunscm=.FALSE.,        & !single column model
!
        &  dt_var=tcall_turb_jg,                                  & !in
        &  nvec=nproma,                                           & !in
        &  ke=nlev, ke1=nlevp1,                                   & !in
        &  kcm=nlevcm, kstart_cloud=kstart_moist(jg),             & !in
        &  iblock=jb,                                             & !in
        &  ivstart=i_startidx, ivend=i_endidx,                    & !in
        &  hhl       = p_metrics%z_ifc(:,:,jb),                   & !in
        &  zvari     = zvari(:,:,:),                              & !out
        &  t_g       = lnd_prog_now%t_g(:,jb),                    & !in
        &  qv_s      = lnd_diag%qv_s(:,jb),                       & !in
        &  ps        = p_diag%pres_sfc(:,jb),                     & !in
        &  u         = p_diag%u(:,:,jb),                          & !inout
        &  v         = p_diag%v(:,:,jb),                          & !inout
        &  t         = p_diag%temp(:,:,jb),                       & !inout
        &  qv        = p_prog_rcf%tracer(:,:,jb,iqv),             & !inout
        &  qc        = p_prog_rcf%tracer(:,:,jb,iqc),             & !inout
        &  prs       = p_diag%pres(:,:,jb),                       & !in
        &  rhoh      = p_prog%rho(:,:,jb),                        & !in
        &  rhon      = zrhon(:,:),                                & !in
        &  epr       = p_prog%exner(:,:,jb),                      & !in
        &  impl_weight=turbdiff_config(jg)%impl_weight,           & !in
        &  ptr       = ptr(:),                                    & !inout
        &  ndtr      = art_config(jg)%nturb_tracer+ncloud_offset, & !in: diffusion of additional tracer variables!
        &  tvm       = prm_diag%tvm(:,jb),                        & !inout
        &  tvh       = prm_diag%tvh(:,jb),                        & !inout
        &  tkvm      = prm_diag%tkvm(:,:,jb),                     & !inout
        &  tkvh      = prm_diag%tkvh(:,:,jb),                     & !inout
        &  u_tens    = prm_nwp_tend%ddt_u_turb(:,:,jb),           & !inout
        &  v_tens    = prm_nwp_tend%ddt_v_turb(:,:,jb),           & !inout
        &  t_tens    = prm_nwp_tend%ddt_temp_turb(:,:,jb),        & !inout
        &  qv_tens   = prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqv),  & !inout
        &  qc_tens   = prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqc),  & !inout
!          qv_conv: missing
        &  shfl_s    = prm_diag%shfl_s(:,jb),                     & !inout
        &  qvfl_s    = prm_diag%qhfl_s(:,jb)                      & !inout
!          umfl_s: missing                                          !inout
!          vmfl_s: missing                                          !inout
        )


      ! re-diagnose turbulent deposition fluxes for qc and qi (positive downward)
      ! So far these fluxes only serve diagnostic purposes. I.e. they 
      ! must be taken into account when checking the atmospheric water mass balance.
      !
      ! ToDo: In the midterm, these fluxes should rather be computed by turbtran and 
      !       being treated analogous to qhfl_s. I.e. they should also be passed to 
      !       the soil/surface scheme TERRA.
      !
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO jc = i_startidx, i_endidx
        tempv_sfc(jc) = lnd_prog_now%t_g(jc,jb) * (1._wp + vtmpc1*lnd_diag%qv_s(jc,jb))
        rho_sfc(jc)   = p_diag%pres_sfc(jc,jb)/(rd*tempv_sfc(jc))
        prm_diag%qcfl_s(jc,jb) = rho_sfc(jc) * prm_diag%tvh(jc,jb) * p_prog_rcf%tracer(jc,nlev,jb,iqc)
      ENDDO
      !$ACC END PARALLEL

      IF (turbdiff_config(jg)%ldiff_qi) THEN
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR
        DO jc = i_startidx, i_endidx
          prm_diag%qifl_s(jc,jb) = rho_sfc(jc) * prm_diag%tvh(jc,jb) * p_prog_rcf%tracer(jc,nlev,jb,iqi)
        ENDDO
        !$ACC END PARALLEL
      ENDIF

      IF (atm_phy_nwp_config(jg)%l_3d_turb_fluxes) THEN
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = 1, nlevp1
          DO jc = i_startidx, i_endidx
                    prm_diag%tetfl_turb(jc,jk,jb) = tet_flux(jc,jk)
                    prm_diag%vapfl_turb(jc,jk,jb) = vap_flux(jc,jk)
                    prm_diag%liqfl_turb(jc,jk,jb) = liq_flux(jc,jk)
          END DO
        END DO
        !$ACC END PARALLEL
      ENDIF

!DR If accumulated deposition fluxes are required ...
!!$      DO jc = i_startidx, i_endidx
!!$        p_diag%extra_2d(jc,jb,1) = p_diag%extra_2d(jc,jb,1) + tcall_turb_jg*prm_diag%qcfl_s(jc,jb)
!!$        p_diag%extra_2d(jc,jb,2) = p_diag%extra_2d(jc,jb,2) + tcall_turb_jg*prm_diag%qifl_s(jc,jb)
!!$      ENDDO
!DR End Test



      ! preparation for concentration boundary condition. Usually inactive for standard ICON runs.
      IF ( .NOT. lsflcnd ) THEN
        !$ACC KERNELS ASYNC(1) DEFAULT(PRESENT)
        prm_diag%lhfl_s(i_startidx:i_endidx,jb) = &
          &  prm_diag%qhfl_s(i_startidx:i_endidx,jb) * alv
        !$ACC END KERNELS
      END IF

#ifdef __ICON_ART
      IF ( lart .AND. art_config(jg)%nturb_tracer > 0 ) THEN
         CALL art_turbdiff_interface( 'update_ptr', p_patch, p_prog_rcf, prm_nwp_tend,  &
           &                          ncloud_offset=ncloud_offset,                      &
           &                          ptr=ptr(:), dt=tcall_turb_jg,                     &
           &                          i_st=i_startidx, i_en=i_endidx,                   &
           &                          lacc=lacc )
      ENDIF
#endif

      ! transform updated turbulent velocity scale back to TKE
      ! Note: ddt_tke is purely diagnostic and has already been added to z_tvs
      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
      !$ACC LOOP GANG
      DO jk=1, nlevp1
        !$ACC LOOP VECTOR
        DO jc=i_startidx, i_endidx
          p_prog_rcf%tke(jc,jk,jb) = 0.5_wp*(z_tvs(jc,jk,1))**2
        ENDDO
      ENDDO
      !$ACC END PARALLEL


      ! Interpolate updated TKE (actually tvs) back to main levels
      ! Note that TKE at lowest main level is re-computed in nwp_turbtrans, after surface TKE
      ! has been updated.
      IF (advection_config(jg)%iadv_tke > 0) THEN
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
        !$ACC LOOP GANG
        DO jk=1, nlev
          !$ACC LOOP VECTOR
          DO jc=i_startidx, i_endidx
            p_prog_rcf%tracer(jc,jk,jb,iqtke) = 0.5_wp* ( z_tvs(jc,jk,1) + z_tvs(jc,jk+1,1) )
          ENDDO
        ENDDO
        !$ACC END PARALLEL
      ENDIF

#ifndef __NO_ICON_EDMF__
    CASE(iedmf)

!-------------------------------------------------------------------------
!> EDMF DUALM turbulence scheme (eddy-diffusivity/mass-flux dual mass-flux)
!-------------------------------------------------------------------------

#ifdef _OPENACC
      CALL finish('GPU:mo_nwp_turbdiff_interface:nwp_turbdiff', 'iedmf unsupported. Only cosmo turbulence is supported on GPU!')
#endif

!     Calculate vertical velocity in p-system

      DO jk = 1,p_patch%nlev
        DO jc = i_startidx,i_endidx
          z_omega_p(jc,jk) = - p_prog%w(jc,jk,jb) * p_prog%rho(jc,jk,jb) * grav
        ENDDO
      ENDDO

!     Various variables for VDFOUTER

      DO jc = i_startidx, i_endidx
        zchar  (jc) = 0.018_wp                 ! default value from IFS if no wave model
        zucurr (jc) = 0.0_wp
        zvcurr (jc) = 0.0_wp
        zsigflt(jc) = 0.0_wp                   ! just for testing (standard dev. of filtered orogrphy)
        zz0m   (jc) = prm_diag%gz0(jc,jb) /grav! diagnostic z0,m from turbtrans_interface
        zz0h   (jc) = prm_diag%gz0(jc,jb) /grav! guess of z0,h = z0,m (not used in EDMF)
        ztskrad(jc) = lnd_prog_now%t_g (jc,jb) ! skin temperature at last radiation step - used in SURFSEB
      ENDDO

      DO jk = 1,nlev
        DO jc = i_startidx, i_endidx
          zsoteu (jc,jk) = 0.0_wp
          zsotev (jc,jk) = 0.0_wp
          zsobeta(jc,jk) = 0.0_wp
          zae    (jc,jk) = 0.0_wp              ! cloud tendency ???
        ENDDO
      ENDDO

!     TERRA tile fluxes: input to EDMF

      DO jt = 1,ntiles_total+ntiles_water
        DO jc = i_startidx, i_endidx
          shfl_s_t(jc,jt) = prm_diag%shfl_s_t (jc,jb,jt) 
          evap_s_t(jc,jt) = prm_diag%lhfl_s_t (jc,jb,jt) / alv ! evaporation [kg/(m2 s)]
         !tskin_t (jc,jt) = lnd_prog_now%t_g_t(jc,jb,jt) 
          tskin_t (jc,jt) = 0.0_wp                       ! not needed as TSK is transfered in t_g_ex
          ustr_s_t(jc,jt) = prm_diag%umfl_s_t (jc,jb,jt) ! prognostic surface stress U (sfc momentum flux)  
          vstr_s_t(jc,jt) = prm_diag%vmfl_s_t (jc,jb,jt) ! prognostic surface stress V (sfc momentum flux)
        ENDDO
      ENDDO

!     Tendencies are set to include dynamics and radiation
!     Question: should SSO tendendies be included in "ddt_u_turb"?
!     ATTENTION: currently for simplicity all input tendencies = 0
!                when updated after vdfouter difference needed (see convection)

      prm_nwp_tend%ddt_u_turb     (:,:,jb)     = 0._wp
      prm_nwp_tend%ddt_v_turb     (:,:,jb)     = 0._wp
      prm_nwp_tend%ddt_temp_turb  (:,:,jb)     = 0._wp
      prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqv) = 0._wp
      prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqc) = 0._wp
      prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqi) = 0._wp

      CALL vdfouter ( &
        & KIDIA   = i_startidx                                 ,&! (IN)
        & KFDIA   = i_endidx                                   ,&! (IN)
        & KLON    = nproma                                     ,&! (IN)
        & KLEV    = p_patch%nlev                               ,&! (IN)
        & KLEVS   = nlev_soil-1                                ,&! (IN)  skip lowermost (climat.) layer
        & KSTEP   = nstep_turb                                 ,&! (IN)  ??? used in surfexdriver!!
        & KTILES  = ntiles_edmf                                ,&! (IN)
        & KTRAC   = itrac_vdf                                  ,&! (IN)  default 0 (itrac?)
        & KLEVSN  = nlev_snow                                  ,&! (IN)  # snow layers (1!)
        & KLEVI   = 1                                          ,&! (IN)  # sea ice layers
        & PTSPHY  = tcall_turb_jg                              ,&! (IN)
!xmk ?  & PSIGFLT = ext_data%atm%sso_stdh(:,jb)                ,&! (IN)  input for TOFD (needs to be passed down!!!)
        & PSIGFLT = zsigflt                                    ,&! (IN)  input for TOFD (needs to be passed down!!!)
        & PUM1    = p_diag%u(:,:,jb)                           ,&! (IN)
        & PVM1    = p_diag%v(:,:,jb)                           ,&! (IN)
        & PTM1    = p_diag%temp(:,:,jb)                        ,&! (IN)
        & PQM1    = p_prog_rcf%tracer(:,:,jb,iqv)              ,&! (IN)
        & PLM1    = p_prog_rcf%tracer(:,:,jb,iqc)              ,&! (IN)
        & PIM1    = p_prog_rcf%tracer(:,:,jb,iqi)              ,&! (IN)
        & PAM1    = prm_diag%clc (:,:,jb)                      ,&! (IN)
        & PAPHM1  = p_diag%pres_ifc         (:,:,jb)           ,&! (IN)
        & PAPM1   = p_diag%pres             (:,:,jb)           ,&! (IN)
        & PGEOM1  = p_metrics%geopot_agl    (:,:,jb)           ,&! (IN)
        & PGEOH   = p_metrics%geopot_agl_ifc(:,:,jb)           ,&! (IN)
        & PTSKM1M = lnd_prog_now%t_g(:,jb)                     ,&! (IN)  T,skin
        & PSSRFL  = prm_diag%swflxsfc(:,jb)                    ,&! (IN)
        & PSLRFL  = prm_diag%lwflxsfc(:,jb)                    ,&! (IN)
        & PEMIS   = ext_data%atm%emis_rad(:,jb)                ,&! (IN)
        & PHRLW   = prm_nwp_tend%ddt_temp_radlw(:,:,jb)        ,&! (IN)
        & PHRSW   = prm_nwp_tend%ddt_temp_radsw(:,:,jb)        ,&! (IN)
        & PTSNOW  = lnd_prog_now%t_snow_t(:,jb,1)              ,&! (IN)  T,snow - unused (attention: tile 1????)
        & PCHAR   = zchar                                      ,&! (IN)  Charnock parameter (for z0 over ocean)
        & PUCURR  = zucurr                                     ,&! (IN)  Ocean current x
        & PVCURR  = zvcurr                                     ,&! (IN)  Ocean current y
        & PTSKRAD = ztskrad                                    ,&! (IN)  unused: T,skin at last radiation step
        & PSOTEU  = zsoteu                                     ,&! (IN)  unused: Explicit part of U-tendency from SSO
        & PSOTEV  = zsotev                                     ,&! (IN)  unused: Explicit part of V-tendency from SSO
        & PSOBETA = zsobeta                                    ,&! (IN)  unused: Implicit part of subgrid orography
        & PVERVEL = z_omega_p                                  ,&! (IN)
        & PZ0M    = zz0m                                       ,&! (INOUT) z0,m (calculated in vupdz0)
        & PZ0H    = zz0h                                       ,&! (INOUT) z0,h (should be diagnostic output ???)
        & PVAR    = p_prog_rcf%tracer(:,:,jb,iqtvar)           ,&! (INOUT) qt,variance - prognostic advected tracer
        & KHPBLN  = khpbln                                     ,&! (OUT) optional out: PBL top level
        & KVARTOP = kvartop                                    ,&! (OUT) optional out: top level of predictied qt,var
        & KPBLTYPE= kpbltype                                   ,&! (OUT) optional out: PBL type
        & PFPLVL  = prm_diag%rain_edmf_rate_3d(:,:,jb)         ,&! (OUT) optional out: PBL rain flux
        & PFPLVN  = prm_diag%snow_edmf_rate_3d(:,:,jb)         ,&! (OUT) optional out: PBL snow flux
        & KLEVX   = p_patch%nlev                               ,&! (IN) 
        & KFLDX   = 0                                          ,&! (IN) 
        & KFLDX2  = 0                                          ,&! (IN)
        & LLDIAG  = .FALSE.                                    ,&! (IN)
        & PTE     = prm_nwp_tend%ddt_temp_turb  (:,:,jb)       ,&! (INOUT)
        & PQE     = prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqv)   ,&! (INOUT)
        & PLE     = prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqc)   ,&! (INOUT)
        & PIE     = prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqi)   ,&! (INOUT)
        & PAE     = zae                                        ,&! (INOUT)
        & PVOM    = prm_nwp_tend%ddt_u_turb(:,:,jb)            ,&! (INOUT)
        & PVOL    = prm_nwp_tend%ddt_v_turb(:,:,jb)            ,&! (INOUT)
        & PUSTRTI = ustr_s_t                                   ,&! (INOUT) tile u stress
        & PVSTRTI = vstr_s_t                                   ,&! (INOUT) tile v stress
        & PAHFSTI = shfl_s_t                                   ,&! (INOUT) tile sensible heat flux
        & PEVAPTI = evap_s_t                                   ,&! (INOUT) tile latent heat flux
        & PTSKTI  = tskin_t                                    ,&! (INOUT) currently unused!!!
        & PDIFTS  = pdifts                                     ,&! (OUT)  optional out: turbulent heat flux
        & PDIFTQ  = pdiftq                                     ,&! (OUT)  optional out: turbulent moisture flux
        & PDIFTL  = pdiftl                                     ,&! (OUT)  optional out: turbulent liquid water flux
        & PDIFTI  = pdifti                                     ,&! (OUT)  optional out: turbulent ice water flux
        & PSTRTU  = pstrtu                                     ,&! (OUT)  optional out: turbulent U flux
        & PSTRTV  = pstrtv                                     ,&! (OUT)  optional out: turbulent V flux
        & PKH     = pkh                                        ,&! (OUT)  used by TURBTRAN !!
        & PKM     = pkm                                        ,&! (OUT)  used by TURBTRAN !!
        & LDLAND  = ext_data%atm%llsm_atm_c(:,jb)              ,&! (IN)   logical for land
! surface fluxes from TERRA for EDMF atmospheric transport:
        & SHFL_S  = prm_diag%shfl_s(:,jb)                      ,&! (IN)    sensible heat flux from TERRA
        & LHFL_S  = prm_diag%lhfl_s(:,jb)                      ,&! (IN)    latent   heat from from TERRA
        & UMFL_S  = prm_diag%umfl_s(:,jb)                      ,&! (IN)    u-flux from TURBTRAN & SFCinterface
        & VMFL_S  = prm_diag%vmfl_s(:,jb)                      ,&! (IN)    v-flux from TURBTRAN & SFCinterface
! TERRA data
        & frac_t  = ext_data%atm%frac_t (:,jb,:)               ,&! (IN)    tile fraction
        & t_g_t   = lnd_prog_now%t_g_t  (:,jb,:)               ,&! (IN)    tile surface temperature (TERRA)
        & qv_s    = lnd_diag%qv_s       (:,jb)                 ,&! (IN) ?? for testcase only
        & t_ice   = wtr_prog_now%t_ice  (:,jb)                  )! (IN) ?? used in vdfdifh for edmf_conf=1


      prm_diag%tkvh(:,2:nlev,jb) = pkh(:,1:nlev-1)     ! K,surface calculated by TURBTRAN
      prm_diag%tkvm(:,2:nlev,jb) = pkm(:,1:nlev-1)             
 
! EDMF precipitation: accumulate in turbulence time steps, in convection add to con_precip and zero out!!!
!     ????????????????????????
!     prm_diag%rain_con_rate(:,jb) = prm_diag%rain_con_rate(:,jb) + prm_diag%rain_edmf_rate_3d(:,nlev+1,jb)
!     prm_diag%rain_edmf_rate_3d(:,:,jb)  PBL rain flux
!     prm_diag%snow_edmf_rate_3d(:,:,jb)  PBL snow flux

! Turbulence updating strategy:
! * Update T, prognostic QV, QC, QI and diagnostic CC with turbulence tendencies
! * Set diagnostic QV, QC, QI equal to prognostic values
! * Give U, V tendencies to dynamics (but see below)

! turn off shallow convection

      DO jc = i_startidx, i_endidx
        !allow Tiedtke shallow convection when not DUALM shallow convection or strcu
        IF (KPBLTYPE(jc) .EQ. 2 .OR. KPBLTYPE(jc) .EQ. 3) THEN
          prm_diag%ldshcv(jc,jb) = .FALSE.
        ELSE
          prm_diag%ldshcv(jc,jb) = .TRUE.
        ENDIF
      ENDDO
#endif

   CASE(igme)

!-------------------------------------------------------------------------
!> GME turbulence scheme 
!-------------------------------------------------------------------------

#ifdef _OPENACC
      CALL finish('GPU:mo_nwp_turbdiff_interface:nwp_turbdiff', 'igme unsupported. Only cosmo turbulence is supported on GPU!')
#endif

      ! turbulent diffusion coefficients in atmosphere
      CALL partura( zh=p_metrics%z_ifc(:,:,jb), zf=p_metrics%z_mc(:,:,jb),                 & !in
        &           u=p_diag%u(:,:,jb),         v=p_diag%v(:,:,jb), t=p_diag%temp(:,:,jb), & !in
        &           qv=p_prog_rcf%tracer(:,:,jb,iqv), qc=p_prog_rcf%tracer(:,:,jb,iqc),    & !in
        &           ph=p_diag%pres_ifc(:,:,jb), pf=p_diag%pres(:,:,jb),                    & !in
        &           ie=nproma, ke=nlev, ke1=nlevp1,                                        & !in
        &           i_startidx=i_startidx, i_endidx=i_endidx,                              & !in
        &           tkvm=prm_diag%tkvm(:,2:nlev,jb), tkvh=prm_diag%tkvh(:,2:nlev,jb)       ) !inout


      ! tendencies from turbulent diffusion
      CALL progimp_turb( t=p_diag%temp(:,:,jb), qv=p_prog_rcf%tracer(:,:,jb,iqv),      & !in
        &                qc=p_prog_rcf%tracer(:,:,jb,iqc),                             & !in
        &                u=p_diag%u(:,:,jb),    v=p_diag%v(:,:,jb),                    & !in
        &                zh=p_metrics%z_ifc(:,:,jb), zf=p_metrics%z_mc(:,:,jb),        & !in
        &                rho=p_prog%rho(:,:,jb), ps=p_diag%pres_ifc(:,nlevp1,jb),      & !in
        &                tkvm=prm_diag%tkvm(:,2:nlev,jb),                              & !in
        &                tkvh=prm_diag%tkvh(:,2:nlev,jb),                              & !in
        &                t_g=lnd_prog_now%t_g(:,jb), qv_s=lnd_diag%qv_s(:,jb),         & !in
        &                h_ice=wtr_prog_now%h_ice(:,jb),                               & !in
        &                tcm=prm_diag%tcm(:,jb), tch=prm_diag%tch(:,jb),               & !in
        &                ie=nproma, ke=nlev, ke1=nlevp1,                               & !in
        &                i_startidx=i_startidx, i_endidx=i_endidx, dt=tcall_turb_jg,   & !in
        &                du_turb=prm_nwp_tend%ddt_u_turb(:,:,jb),                      & !out
        &                dv_turb=prm_nwp_tend%ddt_v_turb(:,:,jb),                      & !out
        &                dt_turb=prm_nwp_tend%ddt_temp_turb(:,:,jb),                   & !out
        &                dqv_turb=prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqv),            & !out
        &                dqc_turb=prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqc))              !out
!       &                shfl_s=prm_diag%shfl_s_t(:,jb,1),                             & !out
!       &                lhfl_s=prm_diag%lhfl_s_t(:,jb,1),                             & !out
!       &                qhfl_s=prm_diag%qhfl_s_t(:,jb,1),                             & !out
!       &                umfl_s=prm_diag%umfl_s(:,jb), vmfl_s=prm_diag%vmfl_s(:,jb)    ) !out

!DR NOTE: computation of sensible and latent heat fluxes (over non-land points) should be 
!DR moved either to the turbtran interface, or the surface interface!!

!     DO jc = i_startidx, i_endidx
!       prm_diag%shfl_s(jc,jb) = prm_diag%shfl_s_t(jc,jb,1)
!       prm_diag%lhfl_s(jc,jb) = prm_diag%lhfl_s_t(jc,jb,1)
!       prm_diag%qhfl_s(jc,jb) = prm_diag%qhfl_s_t(jc,jb,1)
!     ENDDO


    END SELECT !inwp_turb



    ! Update wind speed, QV and temperature with turbulence tendencies
    ! Note: wind speed is updated here, in order to pass u and v at the correct time level
    ! to turbtran and the convection scheme. However, the prognostic variable vn is updated
    ! at the end of the NWP interface by first interpolating the u/v tendencies to the 
    ! velocity points (in order to minimize interpolation errors) and then adding the tendencies
    ! to vn (for efficiency reasons)
    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = 1, nlev
!DIR$ IVDEP
      DO jc = i_startidx, i_endidx
      
        prm_nwp_tend%ddt_u_turb(jc,jk,jb) = &
          SIGN(MIN(0.1_wp,ABS(prm_nwp_tend%ddt_u_turb(jc,jk,jb))),prm_nwp_tend%ddt_u_turb(jc,jk,jb))
        prm_nwp_tend%ddt_v_turb(jc,jk,jb) = &
          SIGN(MIN(0.1_wp,ABS(prm_nwp_tend%ddt_v_turb(jc,jk,jb))),prm_nwp_tend%ddt_v_turb(jc,jk,jb))

        p_prog_rcf%tracer(jc,jk,jb,iqv) =MAX(0._wp, p_prog_rcf%tracer(jc,jk,jb,iqv) &
             &           + tcall_turb_jg*prm_nwp_tend%ddt_tracer_turb(jc,jk,jb,iqv))
        p_diag%temp(jc,jk,jb) = p_diag%temp(jc,jk,jb)  &
         &  + cpd_o_cvd * tcall_turb_jg*prm_nwp_tend%ddt_temp_turb(jc,jk,jb)
        p_diag%u(jc,jk,jb) = p_diag%u(jc,jk,jb) + tcall_turb_jg*prm_nwp_tend%ddt_u_turb(jc,jk,jb)
        p_diag%v(jc,jk,jb) = p_diag%v(jc,jk,jb) + tcall_turb_jg*prm_nwp_tend%ddt_v_turb(jc,jk,jb)
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    ! QC is updated only in that part of the model domain where moisture physics is active
    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = kstart_moist(jg), nlev
!DIR$ IVDEP
      DO jc = i_startidx, i_endidx
        p_prog_rcf%tracer(jc,jk,jb,iqc) =MAX(0._wp, p_prog_rcf%tracer(jc,jk,jb,iqc) &
             &           + tcall_turb_jg*prm_nwp_tend%ddt_tracer_turb(jc,jk,jb,iqc))
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    IF (ltwomoment) THEN
      ! QNC update
      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = kstart_moist(jg), nlev
!DIR$ IVDEP
        DO jc = i_startidx, i_endidx
          p_prog_rcf%tracer(jc,jk,jb,iqnc) = MAX(0.0_wp, p_prog_rcf%tracer(jc,jk,jb,iqnc) &
            &                              + tcall_turb_jg                                & 
            &                              * ddt_turb_qnc(jc,jk))
        ENDDO
      ENDDO
      !$ACC END PARALLEL
    ENDIF ! ltwomoment

    IF (turbdiff_config(jg)%ldiff_qi) THEN
      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = kstart_moist(jg), nlev
!DIR$ IVDEP
        DO jc = i_startidx, i_endidx
          p_prog_rcf%tracer(jc,jk,jb,iqi) = MAX(0.0_wp, p_prog_rcf%tracer(jc,jk,jb,iqi) &
            &                             + tcall_turb_jg                               &
            &                             * prm_nwp_tend%ddt_tracer_turb(jc,jk,jb,iqi))
        ENDDO
      ENDDO
      !$ACC END PARALLEL
      IF (ltwomoment) THEN
        ! QNI update
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = kstart_moist(jg), nlev
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx
            p_prog_rcf%tracer(jc,jk,jb,iqni) = MAX(0.0_wp, p_prog_rcf%tracer(jc,jk,jb,iqni) &
              &                              + tcall_turb_jg                                & 
              &                              * ddt_turb_qni(jc,jk))
          ENDDO
        ENDDO
        !$ACC END PARALLEL
      ENDIF ! ltwomoment
    ENDIF ! ldiff_qi

    IF (turbdiff_config(jg)%ldiff_qs) THEN
      ! QS update
      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = kstart_moist(jg), nlev
!DIR$ IVDEP
        DO jc = i_startidx, i_endidx
          p_prog_rcf%tracer(jc,jk,jb,iqs) = MAX(0.0_wp, p_prog_rcf%tracer(jc,jk,jb,iqs) &
            &                             + tcall_turb_jg                               & 
            &                             * ddt_turb_qs(jc,jk))
        ENDDO
      ENDDO
      !$ACC END PARALLEL
      IF (ltwomoment) THEN
        ! QNS update
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = kstart_moist(jg), nlev
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx
            p_prog_rcf%tracer(jc,jk,jb,iqns) = MAX(0.0_wp, p_prog_rcf%tracer(jc,jk,jb,iqns) &
              &                              + tcall_turb_jg                                &
              &                              * ddt_turb_qns(jc,jk))
          ENDDO
        ENDDO
        !$ACC END PARALLEL
      ENDIF ! ltwomoment
    ENDIF ! ldiff_qs

    ! EDMF: diagnostic clouds are same as prognostic clouds
    IF ( atm_phy_nwp_config(jg)%inwp_turb == iedmf ) THEN
      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = 1, nlev
!DIR$ IVDEP
        DO jc = i_startidx, i_endidx
          prm_diag%tot_cld (jc,jk,jb,iqv) = p_prog_rcf%tracer(jc,jk,jb,iqv)
          prm_diag%tot_cld (jc,jk,jb,iqc) = p_prog_rcf%tracer(jc,jk,jb,iqc)
          prm_diag%tot_cld (jc,jk,jb,iqi) = p_prog_rcf%tracer(jc,jk,jb,iqi)
          prm_diag%clc     (jc,jk,jb)     = MIN(MAX(prm_diag%clc(jc,jk,jb) &
                                          & + tcall_turb_jg * zae(jc,jk), 0._wp), 1._wp)
        ENDDO
      ENDDO
      !$ACC END PARALLEL
    ENDIF

  ENDDO ! jb
  !$ACC WAIT

!$OMP END DO NOWAIT
!$OMP END PARALLEL


  IF ( atm_phy_nwp_config(jg)%inwp_turb == iedmf ) nstep_turb = nstep_turb + 1

  !$ACC END DATA

END SUBROUTINE nwp_turbdiff

END MODULE mo_nwp_turbdiff_interface
