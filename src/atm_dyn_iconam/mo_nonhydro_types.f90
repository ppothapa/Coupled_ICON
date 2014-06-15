#if (defined (__GNUC__) || defined(__SUNPRO_F95) || defined(__SX__))
#define HAVE_F95
#endif
!>

!! Type definition for the dynamical core of ICONAM.
!!
!! @author Almut Gassmann (MPI-M)
!! @author Daniel Reinert (DWD)
!! @author Guenther Zaengl (DWD)
!!
!! @par Revision History
!! Initial release by Daniel Reinert, DWD (2012-02-07)
!! - Moved here from mo_nonhydro_state to avoid circular dependencies
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
MODULE mo_nonhydro_types

  USE mo_kind,                 ONLY: wp, vp
  USE mo_fortran_tools,        ONLY: t_ptr_2d3d,t_ptr_tracer
  USE mo_linked_list,          ONLY: t_var_list


  IMPLICIT NONE

  PRIVATE


  PUBLIC :: t_nh_prog             ! state vector of prognostic variables (type)
  PUBLIC :: t_nh_diag             ! state vector of diagnostic variables (type)
                                  ! on p- and/or z-levels
  PUBLIC :: t_nh_ref              ! state vector of reference state (type)
  PUBLIC :: t_nh_metrics          ! state vector of metrics variables (type)
  PUBLIC :: t_nh_state            ! state vector of nonhydrostatic variables (type)



  ! prognostic variables state vector
  TYPE t_nh_prog

    REAL(wp), POINTER    &
#ifdef _CRAYFTN
    , CONTIGUOUS         &
#endif
      ::                 &
      w(:,:,:),          & !> orthogonal vertical wind (nproma,nlevp1,nblks_c)     [m/s]
      vn(:,:,:),         & !! orthogonal normal wind (nproma,nlev,nblks_e)         [m/s]
      rho(:,:,:),        & !! density (nproma,nlev,nblks_c)                     [kg/m^3]
      exner(:,:,:),      & !! Exner pressure (nproma,nlev,nblks_c)                   [-]
      theta_v(:,:,:),    & !! virtual potential temperature (nproma,nlev,nblks_c)    [K]
      tracer(:,:,:,:),   & !! tracer concentration (nproma,nlev,nblks_c,ntracer) [kg/kg]
      tke   (:,:,:)        !! turbulent kinetic energy                         [m^2/s^2]
                           !! (defined on half levels) with 2 time levels
    TYPE(t_ptr_2d3d),ALLOCATABLE :: tracer_ptr(:)  !< pointer array: one pointer for each tracer
    TYPE(t_ptr_tracer),ALLOCATABLE :: conv_tracer(:,:)  
    TYPE(t_ptr_tracer),ALLOCATABLE :: turb_tracer(:,:)  
  END TYPE t_nh_prog


  ! diagnostic variables state vector
  TYPE t_nh_diag

    REAL(wp), POINTER       &
#ifdef _CRAYFTN
    , CONTIGUOUS            &
#endif
    &  ::                   &
    ! a) variables needed for intermediate storage and physics-dynamics coupling
    &  u(:,:,:),            & ! zonal wind (nproma,nlev,nblks_c)               [m/s]
    &  v(:,:,:),            & ! meridional wind (nproma,nlev,nblks_c)          [m/s]
    &  vt(:,:,:),           & ! tangential wind (nproma,nlev,nblks_e)          [m/s]
    &  omega_z(:,:,:),      & ! relative vertical vorticity at dual grid
                              ! (nproma,nlev,nblks_v)                          [1/s]
    &  vor(:,:,:),          & ! relative vertical vorticity interpolated to cells
                              ! (nproma,nlev,nblks_c)                          [1/s]
    &  ddt_vn_phy(:,:,:),   & ! normal wind tendency from forcing
                              ! (nproma,nlev,nblks_e)                          [m/s^2]
    &  ddt_exner_phy(:,:,:),& ! exner pressure tendency from physical forcing 
                              ! (nproma,nlev,nblks_c)                     [1/s]
    &  ddt_temp_dyn(:,:,:), & ! rediagnosed temperature tendency from dynamics [K/s]
    &  ddt_tracer_adv(:,:,:,:), &! advective tendency of tracers          [kg/kg/s]
    &  tracer_vi(:,:,:),    & ! vertically integrated tracers(for q1,q2,q3) [kg/m**2]
    &  tracer_vi_avg(:,:,:),& ! average since last output of tracer_vi [kg/m**2]
    &  exner_old(:,:,:),    & ! exner pres from previous step (nproma,nlev,nblks_c)
    &  exner_dyn_incr(:,:,:), & ! exner pres dynamics increment (nproma,nlev,nblks_c)
    &  temp(:,:,:),         & ! temperature (nproma,nlev,nblks_c)                 [K]
    &  tempv(:,:,:),        & ! virtual temperature (nproma,nlev,nblks_c)         [K]
    &  temp_ifc(:,:,:),     & ! temperature at half levels (nproma,nlevp1,nblks_c)[K]
    &  pres(:,:,:),         & ! pressure (nproma,nlev,nblks_c)                  [Pa]
    &  pres_ifc(:,:,:),     & ! pressure at interfaces (nproma,nlevp1,nblks_c)  [Pa]
    &  pres_sfc(:,:),       & ! diagnosed surface pressure (nproma,nblks_c)     [Pa]
    &  pres_sfc_old(:,:),   & ! diagnosed surface pressure at previous timestep (nproma,nblks_c) [Pa]
    &  pres_msl(:,:),       & ! diagnosed mean sea level pressure (nproma,nblks_c)  [Pa]
    &  dpres_mc(:,:,:),     & ! pressure thickness at masspoints(nproma,nlevp,nblks_c)  [Pa]
    &  omega(:,:,:),        & ! vertical velocity ( omega=dp/dt )           [Pa/s]
    &  hfl_tracer(:,:,:,:), & ! horizontal tracer flux at edges             [kg/m/s]
                              ! (nproma,nlev,nblks_e,ntracer)
    &  vfl_tracer(:,:,:,:), & ! vertical tracer flux at cells               [kg/m/s]
                              ! (nproma,nlevp1,nblks_c,ntracer)
    &  div(:,:,:),          & ! divergence(nproma,nlev,nblks_c)     [1/s]
    &  div_ic(:,:,:),       & ! divergence at half levels(nproma,nlevp1,nblks_c)     [1/s]
    &  hdef_ic(:,:,:),      & ! horizontal wind field deformation (nproma,nlevp1,nblks_c)     [1/s^2]
    &  dwdx(:,:,:),         & ! zonal gradient of vertical wind speed (nproma,nlevp1,nblks_c)     [1/s]
    &  dwdy(:,:,:),         & ! meridional gradient of vertical wind speed (nproma,nlevp1,nblks_c)     [1/s]
    &  mass_fl_e(:,:,:),    & ! horizontal mass flux at edges (nproma,nlev,nblks_e) [kg/m/s]
    &  mass_fl_e_sv(:,:,:), & ! storage field for horizontal mass flux at edges (nproma,nlev,nblks_e) [kg/m/s]
    &  rho_ic(:,:,:),       & ! density at half levels (nproma,nlevp1,nblks_c)     [kg/m^3]
    &  theta_v_ic(:,:,:),   & ! theta_v at half levels (nproma,nlevp1,nblks_c)         [K]
    &  w_concorr_c(:,:,:),  & ! contravariant vert correction (nproma,nlevp1,nblks_c)[m/s]
    &  vn_ie(:,:,:),        & ! normal wind at half levels (nproma,nlevp1,nblks_e)   [m/s]
    &  ddt_vn_adv(:,:,:,:), & ! normal wind tendency from advection
                              ! (nproma,nlev,nblks_e,1:3)                    [m/s^2]
    &  ddt_w_adv(:,:,:,:),  & ! vert. wind tendency from advection
                              ! (nproma,nlevp1,nblks_c,1:3)                  [m/s^2]
    &  airmass_now(:,:,:),  & ! mass of air in layer at physics time step now [kg/m^2]
    &  airmass_new(:,:,:),  & ! mass of air in layer at physics time step new [kg/m^2]

    !
    ! b) variables needed for grid nesting
    &  grf_tend_vn(:,:,:),  & ! vn tendency field for use in grid refinement
                              ! (nproma,nlev,nblks_e)                        [m/s^2]
    &  grf_tend_w(:,:,:),   & ! w tendency field for use in grid refinement
                              ! (nproma,nlevp1,nblks_c)                      [m/s^2]
    &  grf_tend_rho(:,:,:), & ! rho tendency field for use in grid refinement
                              ! (nproma,nlev,nblks_c)                     [kg/m^3/s]
    &  grf_tend_mflx(:,:,:),& ! rho*vn tendency field for use in grid refinement
                              ! (nproma,nlev,nblks_e)                     [kg/m^2/s^2]
    &  grf_bdy_mflx(:,:,:),&  ! rho*vn boundary field for use in grid refinement
                              ! (nlev,npoints,2)                            [kg/m^2/s^2]
    &  grf_tend_thv(:,:,:), & ! theta_v tendency field for use in grid refinement
                              ! (nproma,nlev,nblks_c)                          [K/s]
    &  grf_tend_tracer(:,:,:,:), & ! tracer tendency field for use in grid refinement
                                   ! (nproma,nlev,nblks_c,ntracer)          [kg/kg/s]
    &  dvn_ie_int(:,:),    & ! Storage field for vertical nesting: vn at parent interface level
    &  dvn_ie_ubc(:,:),    & ! Storage field for vertical nesting: vn at child upper boundary
    &  mflx_ic_int(:,:,:), & ! Storage field for vertical nesting: mass flux at parent interface level
    &  mflx_ic_ubc(:,:,:), & ! Storage field for vertical nesting: mass flux at child upper boundary
    &  dtheta_v_ic_int(:,:,:),& ! Storage field for vertical nesting: theta at parent interface level
    &  dtheta_v_ic_ubc(:,:),& ! Storage field for vertical nesting: theta at child upper boundary
    &  dw_int(:,:,:),      & ! Storage field for vertical nesting: w at parent interface level
    &  dw_ubc(:,:),        & ! Storage field for vertical nesting: w at child upper boundary
    &  q_int(:,:,:),       & ! Storage field for vertical nesting: q at parent interface level
    &  q_ubc(:,:,:),       & ! Storage field for vertical nesting: q at child upper boundary

    !
    ! c) analysis increments (necessary for incremental analysis update)
    &  vn_incr   (:,:,:),   & ! normal velocity increment        [m/s]
    &  exner_incr(:,:,:),   & ! exner inrement                   [-]
    &  rho_incr  (:,:,:),   & ! moist density increment          [kg/m^3]
    &  qv_incr   (:,:,:)      ! specific humidity increment      [kg/kg]



    REAL(wp), POINTER ::    & !
     &  extra_2d(:,:,:)  ,  & !> extra debug output in 2d and
     &  extra_3d(:,:,:,:)     !!                       3d


    TYPE(t_ptr_2d3d),ALLOCATABLE ::   &
      &  ddt_grf_trc_ptr(:),   &  !< pointer array: one pointer for each tracer
      &  hfl_trc_ptr    (:),   &  !< pointer array: one pointer for each tracer
      &  vfl_trc_ptr    (:),   &  !< pointer array: one pointer for each tracer
      &  ddt_trc_adv_ptr(:),   &  !< pointer array: one pointer for each tracer
      &  ddt_vn_adv_ptr (:),   &  !< pointer array: one pointer for each tracer
      &  ddt_w_adv_ptr  (:),   &  !< pointer array: one pointer for each tracer
      &  q_int_ptr      (:),   &  
      &  q_ubc_ptr      (:),   &
      &  tracer_vi_ptr  (:),   &  !< pointer array: one pointer for each tracer
      &  tracer_vi_avg_ptr(:), &  !< pointer array: one pointer for each tracer
      &  extra_2d_ptr   (:),   &
      &  extra_3d_ptr   (:)

  END TYPE t_nh_diag


  TYPE t_nh_ref
    REAL(wp), POINTER ::    &
      vn_ref  (:,:,:),      & !! orthogonal normal wind (nproma,nlev,nblks_e)      [m/s]
      w_ref   (:,:,:)         !> orthogonal vertical wind (nproma,nlevp1,nblks_c)  [m/s]
  END TYPE t_nh_ref


  TYPE t_nh_metrics

    ! Variables that are always in double precision
    REAL(wp), POINTER      &
#ifdef _CRAYFTN
    , CONTIGUOUS           &
#endif
     ::                    &
     ! a) General geometric quantities
     !
     z_ifc(:,:,:)        , & ! geometric height at the vertical interface of cells (nproma,nlevp1,nblks_c)
     z_mc(:,:,:)         , & ! geometric height at full levels (nproma,nlev,nblks_c)
     ddqz_z_full(:,:,:)  , & ! functional determinant of the metrics [sqrt(gamma)] (nproma,nlev,nblks_c)
     geopot(:,:,:)       , & ! geopotential at cell center (nproma,nlev,nblks_c)
     geopot_agl(:,:,:)   , & ! geopotential above ground level at cell center (nproma,nlev,nblks_c)
     geopot_agl_ifc(:,:,:),& ! geopotential above ground level at interfaces and cell center (nproma,nlevp1,nblks_c)
     dgeopot_mc(:,:,:)   , & ! geopotential at cell center (nproma,nlev,nblks_c)
     !
     ! b) Specific fields for the dynamical core
     !
     rayleigh_w(:)       , & ! Rayleigh damping on the vertical velocity
     rayleigh_vn(:)      , & ! Rayleigh damping on the normal velocity
     enhfac_diffu(:)     , & ! Enhancement factor for nabla4 background diffusion
     vwind_expl_wgt(:,:)  , & ! explicit weight in vertical wind solver (nproma,nblks_c)
     vwind_impl_wgt(:,:)  , & ! implicit weight in vertical wind solver (nproma,nblks_c)
     vwind_impl_wgt_sv(:,:),& ! save array for the above coefficient field (nproma,nblks_c)
     !
     ! c) Fields for the reference atmosphere
     !
     theta_ref_mc(:,:,:) , & 
     theta_ref_me(:,:,:) , & 
     theta_ref_ic(:,:,:) , & 
     tsfc_ref(:,:)       , & 
     exner_ref_mc(:,:,:) , & 
     rho_ref_mc  (:,:,:) , &  
     rho_ref_me  (:,:,:) , & 
     !
     ! d) Fields for truly horizontal temperature diffusion
     !
     zd_intcoef(:,:) , & 
     zd_geofac(:,:)  , & 
     zd_e2cell(:,:)  , & 
     zd_diffcoef(:)  , & 
     !
     ! e) Fields for LES Model : Anurag Dipankar, MPIM (2013-04)
     !
     ! Vertical grid related
     inv_ddqz_z_half_e(:,:,:)  , & 
     inv_ddqz_z_full_e(:,:,:)  , & 
     inv_ddqz_z_half(:,:,:)    , & 
     inv_ddqz_z_half_v(:,:,:)  , & 
     wgtfac_v(:,:,:)           , & 
     ! Mixing length for Smagorinsky model
     mixing_length_sq(:,:,:)   , & 
     !
     ! f) Other stuff
     !
     ! Correction term needed to use perturbation density for lateral boundary nudging
     ! (note: this field is defined on the local parent grid in case of MPI parallelization)
     rho_ref_corr(:,:,:) , & 
     ! Area of subdomain for which feedback is performed; dim: (nlev)
     fbk_dom_volume(:)

    ! Variables that are in single precision when "__MIXED_PRECISION" is defined
    REAL(vp), POINTER      &
#ifdef _CRAYFTN
    , CONTIGUOUS           &
#endif
     ::                    &
     ! a) Layer thicknesses
     !
     ddxn_z_full(:,:,:)  , & ! slope of the terrain in normal direction (nproma,nlev,nblks_e)
     ddxt_z_full(:,:,:)  , & ! slope of the terrain in tangential direction (nproma,nlev,nblks_e)
     ddqz_z_full_e(:,:,:), & ! functional determinant of the metrics [sqrt(gamma)] (nproma,nlev,nblks_e)
     ddqz_z_half(:,:,:)  , & ! functional determinant of the metrics [sqrt(gamma)] (nproma,nlevp1,nblks_c)
     inv_ddqz_z_full(:,:,:),& ! Inverse layer thickness of full levels (nproma,nlev,nblks_c)
     !
     ! b) Interpolation coefficients
     !
     wgtfac_c(:,:,:)      , & ! weighting factor for interpolation from full to half levels (nproma,nlevp1,nblks_c)
     wgtfac_e(:,:,:)      , & ! weighting factor for interpolation from full to half levels (nproma,nlevp1,nblks_e)
     wgtfacq_c(:,:,:)     , & ! weighting factor for quadratic interpolation to surface (nproma,3,nblks_c)
     wgtfacq_e(:,:,:)     , & ! weighting factor for quadratic interpolation to surface (nproma,3,nblks_e)
     wgtfacq1_c(:,:,:)    , & ! weighting factor for quadratic interpolation to model top (nproma,3,nblks_c)
     wgtfacq1_e(:,:,:)    , & ! weighting factor for quadratic interpolation to model top (nproma,3,nblks_e)
     coeff_gradekin(:,:,:), & ! Coefficients for improved discretization of horizontal kinetic energy gradient (nproma,2,nblks_e)
     coeff1_dwdz(:,:,:)   , & 
     coeff2_dwdz(:,:,:)   , & ! Coefficients for second-order accurate dw/dz term (nproma,nlev,nblks_c)
     zdiff_gradp(:,:,:,:) , & ! Height differences between local edge point and neighbor cell points used for
                              ! pressure gradient computation (2,nproma,nlev,nblks_e)
     coeff_gradp(:,:,:,:) , & ! Interpolation coefficients for cubic interpolation of Exner pressure (8,nproma,nlev,nblks_e)
     exner_exfac(:,:,:)   , & ! extrapolation factor for Exner pressure (slope-dependent for stability optimization) 
     !
     ! c) Fields for reference atmosphere
     !
     d_exner_dz_ref_ic(:,:,:), & 
     d2dexdz2_fac1_mc(:,:,:) , & 
     d2dexdz2_fac2_mc(:,:,:) , &
     !
     ! d) other stuff
     !
     pg_exdist (:)           ! extrapolation distance needed for igradp_method = 3


    INTEGER, POINTER          &
#ifdef _CRAYFTN
    , CONTIGUOUS              &
#endif
     ::                       &
     vertidx_gradp(:,:,:,:) , &  ! Vertical index of neighbor points needed for Taylor-expansion-based 
                                 ! pressure gradient (2,nproma,nlev,nblks_e)
     !
     ! Fields for truly horizontal temperature diffusion
     zd_indlist(:,:) , & 
     zd_blklist(:,:) , & 
     zd_edgeidx(:,:) , & 
     zd_edgeblk(:,:) , & 
     zd_vertidx(:,:) , & 
     !
     ! Fields for igradp_method = 3
     pg_edgeidx(:) , & 
     pg_edgeblk(:) , & 
     pg_vertidx(:) , & 
     !
     ! Index lists for grid points on which lateral boundary nudging is applied
     nudge_c_idx(:) , & 
     nudge_e_idx(:) , & 
     nudge_c_blk(:) , & 
     nudge_e_blk(:) , & 
     !
     ! Index lists and mask fields needed to minimize the number of halo communications
     ! a) index lists for halo points belonging to the nest boundary region
     bdy_halo_c_idx(:) , & 
     bdy_halo_c_blk(:) , & 
     ! b) index lists for halo points belonging to a nest overlap region
     !    the additional dimension is n_childdom
     ovlp_halo_c_dim(:)   , &
     ovlp_halo_c_idx(:,:) , & 
     ovlp_halo_c_blk(:,:) , & 
     ! c) index lists for mass fluxes at lateral nest boundary (including the required halo points)
     bdy_mflx_e_idx(:) , & 
     bdy_mflx_e_blk(:)



   ! Corresponding scalar list dimensions
   INTEGER  :: zd_listdim  ! for truly horizontal temperature diffusion
   INTEGER  :: pg_listdim  ! for igradp_method = 3
   INTEGER  :: nudge_c_dim, nudge_e_dim ! for grid points on which lateral boundary nudging is applied
   INTEGER  :: bdy_halo_c_dim ! for halo points belonging to the nest boundary region
   INTEGER  :: bdy_mflx_e_dim ! for mass fluxes at lateral nest boundary


   ! Finally, a mask field that excludes boundary halo points
   LOGICAL,  POINTER :: mask_prog_halo_c(:,:)

  END TYPE t_nh_metrics


!-------------------------------------------------------------------------
!                      STATE VECTORS AND LISTS
!-------------------------------------------------------------------------
  TYPE t_nh_state

    !array of prognostic states at different timelevels
    TYPE(t_nh_prog),  ALLOCATABLE :: prog(:)       !< shape: (timelevels)
    TYPE(t_var_list), ALLOCATABLE :: prog_list(:)  !< shape: (timelevels)

    TYPE(t_nh_diag)    :: diag
    TYPE(t_var_list)   :: diag_list

    TYPE(t_nh_ref)     :: ref
    TYPE(t_var_list)   :: ref_list

    TYPE(t_nh_metrics) :: metrics
    TYPE(t_var_list)   :: metrics_list

    TYPE(t_var_list), ALLOCATABLE :: tracer_list(:) !< shape: (timelevels)

  END TYPE t_nh_state


END MODULE mo_nonhydro_types




