!+ Data module for variables of the turbulence parameterization
!------------------------------------------------------------------------------

MODULE turb_data

!------------------------------------------------------------------------------
!
! Description:
!  This module contains parameters that are used in the turbulence
!  parameterizations. With some of these parameters a tuning of the schemes
!  is possible.
!
! Current Code Owner: DWD, Matthias Raschendorfer
!  phone:  +49  69  8062 2708
!  fax:    +49  69  8062 3721
!  email:  Matthias.Raschendorfer@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! V5_4a        2016-05-10 Matthias Raschendorfer, Ulrich Schaettler
!  Initial release, based on data_turbulence with modifications for ICON
!  Note that this module is now also used for the old (not blocked) scheme
! V5_4b        2016-07-12 Ulrich Schaettler
!  Activate old defaults for some namelist variables for COSMO 
!  as long as non-blocked turbulence scheme is running per default
! V5_4c        2016-10-06 Ulrich Schaettler
!  Again use local memory if not running on GPUs 
!     (because on vectorization problem on CRAY)
! V5_4d        2016-12-12 Matthias Raschendorfer
!  Adapting some comments and reshuffling variables.
! V5_4e        2017-03-23 Ulrich Schaettler
!  Modified several Namelist defaults to the ICON defaults
! V5_4f        2017-09-01 Matthias Raschendorfer, Ulrich Schaettler
!  Unified settings for ICON and COSMO:
!   - imode_rat_sea:    1  (ICON value,  changes in COSMO)
!   - loutshs        TRUE  (COSMO value, changes in ICON)
!  To get an old-COSMO-like behaviour of the blocked turbulence scheme some 
!   hardcoded variables need different values. These are set by compiling 
!   with the pragma -DCOSMO_OLD
!  Removed wichfakt, securi from the old version (US)
! V5_4g        2017-11-13 Ulrich Schaettler
!  Implemented namelist variable loldtur (Default .FALSE.) and removed 
!  pragma -DCOSMO_OLD to be able to switch the behaviour during run time.
! V5_4h        2017-12-15 Xavier Lapillonne
!  Modifications to port turbulence scheme to GPU
! V5_5         2018-02-23 Ulrich Schaettler
!  Updated with ICON Version 7bcba73: 
!   - new (still internal) switch imode_tkesso (should replace ltkesso)
!   - modified some default values of namelist variables with ifdef ICON/COSMO
! V5_6         2019-02-27 Ulrich Schaettler
!  Updated with ICON Version d7e0252
!    (Set alpha1 = 0.75; was 1.0 before)
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Declarations:
!
! Modules used:

USE mo_kind,                ONLY: wp           ! KIND-type parameter for real variables

USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config

USE mo_turbdiff_config,     ONLY: turbdiff_config

!==============================================================================

IMPLICIT NONE

PUBLIC

!==============================================================================
! Configuration parameters:
! ----------------------------------------------------------------------------

INTEGER, PARAMETER :: &
!
    mom=1,        & ! index for a momentum variable
    sca=2           ! index for a scalar   variable

INTEGER, PARAMETER :: &
!
!  Indexgrenzen:
!
   nscal=3,     & !aktive skalare Groessen 1-ter Ordnung ('tem', 'vap', 'liq')
   ninv=2,      & !daraus abgeleitete gegenueber vertikalen (feuchtadiabatischen)
                  !Verrueckungen invarianten Groessen ('tet_l', 'h2o_g')
   nvel=2,      & !aktive Geschwindigkeitskomponenten ('u_m', 'v_m')
   nmvar=nscal+nvel, &
   nred=ninv+nvel,   &
   naux=5,      & !number of auxilary variables
   ntyp=2,      & !Anzahl von Variablentypen (mom) und (sca)
   ntmax=3,     & !max. Anzahl der Zeitebenen fuer die TKE
!
!    Zeiger fuer die Variablen :
!
   u_m=1,       & !zonale Geschw.komp. im Massenzentrum
   v_m=2,       & !meridionale  ,,      ,,     ,,
   tet_l=3,     & !feucht-potentielle Temperatur
   tem_l=tet_l, & !Fluessigwasser-Temperatur
   h2o_g=4,     & !Gesamtwasseergehalt
   liq=5,       & !Fluessigwasser  ,,
   w_m=6,       & !vertikale Geschw.komp. im Massenzentrum
!
   tet=tet_l,   & !pot.Temperatur
   tem=tet,     & !Temperatur
   vap=h2o_g,   & !Wasserdampfmischungsverh.

   ndim=MAX(nmvar,naux)

!  Beachte:     u_m,v_m muessen in [1,nvel] liegen;
!         aber: tet_l,h2o_g in [nvel+1,nred]
!         und   tem (tet),vap,liq in [nvel+1,nmvar]


!==============================================================================
! Parameters that may be used for tuning and special configurations:
! ----------------------------------------------------------------------------

! Attention:
! The given initializations are default settings of the boundary layer
! parameters. Some of these initial parameter values may be changed afterwards
! by model input NAMELISTs!

! 1. Numerical parameters:
!-----------------------------

REAL (KIND=wp)     ::        &
  impl_s       =  1.20_wp,   & ! implicit weight near the surface (maximal value)
  impl_t       =  0.75_wp,   & ! implicit weight near top of the atmosphere (maximal value)

  ! Minimal diffusion coefficients in [m^2/s] for vertical
  tkhmin       =  0.75_wp,   & ! scalar (heat) transport
  tkmmin       =  0.75_wp,   & ! momentum transport
  tkhmin_strat =  0.75_wp,   & ! scalar (heat) transport, enhanced value for stratosphere
  tkmmin_strat =  4.00_wp,   & ! momentum transport,      enhanced value for stratosphere

  ditsmot      =  0.00_wp,   & ! smoothing factor for direct time-step iteration

  tndsmot      =  0.00_wp,   & ! vertical smoothing factor for diffusion tendencies
  frcsmot      =  0.00_wp,   & ! vertical smoothing factor for TKE forcing (in ICON only in the tropics)
  tkesmot      =  0.15_wp,   & ! time smoothing factor for TKE and diffusion coefficients
  stbsmot      =  0.00_wp,   & ! time smoothing factor for stability function
  frcsecu      =  1.00_wp,   & ! security factor for TKE-forcing       (<=1)
  tkesecu      =  1.00_wp,   & ! security factor in  TKE equation      (out of [0; 1])
  stbsecu      =  0.00_wp,   & ! security factor in stability function (out of [0; 1])

  epsi         =  1.0E-6_wp    ! relative limit of accuracy for comparison of numbers

INTEGER            ::        &

  it_end       =  1            ! number of initialization iterations (>=0)

! 2. Parameters describing physical properties of the lower boundary 
!    of the atmosphere:
!------------------------------------------

REAL (KIND=wp)     ::        &
  rlam_mom     =  0.0_wp,    & ! scaling factor of the laminar boundary layer for momentum
  rlam_heat    = 10.0_wp,    & ! scaling factor of the laminar boundary layer for heat

  rat_lam      =  1.0_wp,    & ! ratio of laminar scaling factors for vapour and heat
  rat_can      =  1.0_wp,    & ! ratio of canopy height over z0m
  rat_sea      =  0.8_wp,    & ! ratio of laminar scaling factors for heat over sea and land
  rat_glac     =  3.0_wp,    & ! ratio of laminar scaling factors for heat for glaciers

  z0m_dia      =  0.2_wp,    & ! roughness length of a typical synoptic station [m]

  alpha0       =  0.0123_wp, & ! Charnock-parameter
  alpha0_max   =  0.0335_wp, & ! upper limit of velocity-dependent Charnock-parameter
  alpha0_pert  =  0.0_wp,    & ! additive ensemble perturbation of Charnock-parameter

  alpha1       =  0.7500_wp    ! parameter scaling the molecular roughness of water waves


  !$ACC DECLARE COPYIN(alpha0, alpha0_max, alpha0_pert)

! 3. Parameters that should be external parameter fields being not yet 
!    available:
!------------------------------------------

REAL (KIND=wp)     ::        &
  c_lnd        = 2.0_wp,     & ! surface area density of the roughness elements over land
  c_sea        = 1.5_wp,     & ! surface area density of the waves over sea
  c_soil       = 1.0_wp,     & ! surface area density of the (evaporative) soil surface
  e_surf       = 1.0_wp        ! exponent to get the effective surface area


! 4. Parameters that should be dynamical fields being not yet available:
!------------------------------------------

REAL (KIND=wp)     ::        &
  z0_ice       =  0.001_wp     !roughness length of sea ice


! 5. Parameters for modelling turbulent diffusion:
!------------------------------------------

REAL (KIND=wp)     ::        &
  tur_len      = 500.0_wp,   & ! asymptotic maximal turbulent distance [m]
  pat_len      = 100.0_wp,   & ! effective length scale of subscale surface patterns over land [m]
                               ! (should be dependent on location)
  len_min      =  1.0E-6_wp, & ! minimal turbulent length scale [m]

  vel_min      =  0.01_wp,   & ! minimal velocity scale [m/s]

  akt          =  0.4_wp,    & ! von Karman-constant

  ! Length scale factors for pressure destruction of turbulent
  a_heat       =  0.74_wp,   & ! scalar (heat) transport
  a_mom        =  0.92_wp,   & ! momentum transport

  ! Length scale factors for dissipation of
  d_heat       =  10.1_wp,   & ! scalar (temperature) variance
  d_mom        =  16.6_wp,   & ! momentum variance

  ! Length scale factors for turbulent transport (vertical diffusion)
  c_diff       =  0.20_wp,   & ! of TKE

  ! Length scale factor for separate horizontal shear production
  a_hshr       =  1.00_wp,   & ! of TKE

  ! Length scale factor for the stability correction
  a_stab       =  0.00_wp,   & ! no stability correction so far

  ! Dimensionless parameters used in the sub grid scale condensation scheme
  ! (statistical cloud scheme):
  clc_diag     =  0.5_wp,    & !cloud cover at saturation
  q_crit       =  1.6_wp,    & !critical value for normalized over-saturation
  c_scld       =  1.0_wp       !factor for liquid water flux density in sub grid scale clouds

!==============================================================================
! Switches controlling the turbulence model, turbulent transfer and diffusion:
! ----------------------------------------------------------------------------

LOGICAL :: &

  loldtur       =.FALSE., & ! use settings to simulate old ijk turbulence version
                            ! if .TRUE.: new ICON-like settings are used
  ltkesso       =.TRUE.,  & ! calculation SSO-wake turbulence production for TKE
  ltkecon       =.FALSE., & ! consider convective buoyancy production for TKE
  ltkeshs       =.TRUE. , & ! consider separ. horiz. shear production for TKE
  loutshs       =.TRUE. , & ! consider separ. horiz. shear production of TKE for output
  lnonloc       =.FALSE., & ! nonlocal calculation of vertical gradients used for turbul. diff.
  lprfcor       =.FALSE., & ! using the profile values of the lowest main level instead of
                            ! the mean value of the lowest layer for surface flux calulations

  ltmpcor       =.FALSE., & ! consideration of thermal TKE-sources in the enthalpy budget
  lcpfluc       =.FALSE., & ! consideration of fluctuations of the heat capacity of air

  lexpcor       =.FALSE., & ! explicit corrections of the implicit calculated turbul. diff.

! for semi-implicit vertical diffusion: 
  lsflcnd       =.TRUE. , & ! lower flux condition for vertical diffusion calculation
  ldynimp       =.FALSE., & ! dynamical calculation of implicit weights
  lprecnd       =.FALSE., & ! preconditioning of tridiagonal matrix
  lfreeslip     =.FALSE.    ! free-slip lower boundary condition (use for idealized runs only!)

! Notice that the following switches are provided by the parameter-list of 
! SUB 'turb_diffusion' or 'turb_transfer':

! lstfnct                   :calculation of stability function required
! lnsfdia                   :calculation of (synoptical) near-surface variables required
! lmomdif                   :calculation of complete gradient diffusion of horizontal momenum
! lscadif                   :calculation of complete gradient diffusion of scalar properties
! lturatm                   :running turbulence model between atmosph. layers (updating diffusion coefficients)
! ltursrf                   :running turbulence model at the surface layer (updating transfer coefficients
! lsfluse                   :use explicit heat flux densities at the surface
! ltkeinp                   :TKE present as input (at level k=ke1 for current time level 'ntur')
! lgz0inp                   :gz0 present as input

!==============================================================================
! Selectors controlling the turbulence model, turbulent transfer and diffusion:
! ----------------------------------------------------------------------------

INTEGER :: &

  imode_tran    =0, & ! mode of TKE-equation in transfer scheme             (compare 'imode_turb')
  imode_turb    =1, & ! mode of TKE-equation in turbulence scheme
                      !  0: diagnostic equation
                      !  1: prognostic equation (default)
                      !  2: prognostic equation (implicitly positive definit)
  icldm_tran    =2, & ! mode of cloud representation in transfer parametr.  (compare 'icldm_turb')
  icldm_turb    =2, & ! mode of cloud representation in turbulence parametr.
                      ! -1: ignoring cloud water completely (pure dry scheme)
                      !  0: no clouds considered (all cloud water is evaporated)
                      !  1: only grid scale condensation possible
                      !  2: also sub grid (turbulent) condensation considered
  itype_wcld    =2, & ! type of water cloud diagnosis within the turbulence scheme:
                      ! 1: employing a scheme based on relative humitidy
                      ! 2: employing a statistical saturation adjustment
  itype_sher    =0    ! type of shear production for TKE
                      ! 0: only vertical shear of horizontal wind
                      ! 1: previous plus horizontal shear correction
                      ! 2: previous plus shear from vertical velocity
INTEGER :: &
  itype_diag_t2m=1    ! type of diagnostics of 2m-temperature and -dewpoint
                      ! 1: Considering a fictive surface roughness of a SYNOP lawn
                      ! 2: Considering the mean surface roughness of a grid box
                      !    and using an exponential roughness layer profile

! To reproduce the old ijk turbulence settings as good as possible, all these switches 
! have to be set to 1. If loldtur=.TRUE., the re-setting is done in organize_physics.

! These are the settings for the ICON-like setup of the physics
INTEGER :: &

  imode_stbcorr =1, & ! mode of correcting the stability function (related to 'stbsecu')
                      ! 1: always for strict.-non-stb. strat. using a restr. 'gama' in terms of prev. forc.
                      ! 2: only to avoid non-physic. solution or if current 'gama' is too large
  imode_tkemini =1, & ! mode of fixing a lower limit of q=2TKE**2
                      ! 1: by using 'vel_min' only
                      ! 2: by adapting to minimal diffusion coefficients
  imode_lamdiff =1, & ! mode of considering laminar diffusion at surface layer
                      ! 1: only when calculating the profile functions
                      ! 2: surface-layer diffusion coeff. always at least at laminar value
  ilow_def_cond =2, & ! type of the default condition at the lower boundary
                      ! 1: zero surface flux density
                      ! 2: zero surface value
  imode_calcirc =2, & ! mode of treating the circulation term (related to 'pat_len', imode_pat_len')
                      ! 1: explicit calculation of the flux convergence
                      ! 2: quasi implicit treatment by calculation of effective TKE-gradients
  imode_pat_len =2, & ! mode of determining the length scale of surface patterns (related to 'pat_len')
                      ! 1: by the constant value 'pat_len' only
                      ! 2: and the std. deviat. of SGS orography as a lower limit (only if 'd_pat' is pres.)
  imode_frcsmot =2, & ! if "frcsmot>0", apply smoothing of TKE source terms 
                      ! 1: globally or 
                      ! 2: in the tropics only (if 'trop_mask' is present) 
  imode_shshear =2, & ! mode of calculat. the separated horizontal shear mode (related to 'ltkeshs', 'a_hshr')
                      ! 0: with a constant lenght scale and based on 3D-shear and incompressibility
                      ! 1: with a constant lenght scale and considering the trace constraint for the 2D-strain tensor
                      ! 2: with a Ri-depend. length sclale correct. and the trace constraint for the 2D-strain tensor
  imode_tkesso =1, &  ! mode of calculat. the SSO source term for TKE production
                      ! 1: original implementation
                      ! 2: with a Ri-dependent reduction factor for Ri>1
  imode_tkvmini =2, & ! mode of calculating the minimal turbulent diff. coeffecients
                      ! 1: with a constant value
                      ! 2: with a stability dependent correction
  imode_rat_sea =1, & ! mode of scaling the laminar resistance for heat over sea (related to 'rat_sea')
                      ! 1: constant ratio 'rat_sea' to land 
                      ! 2: with a correction for a strongly overheated SST
  imode_vel_min =2, & ! mode of calculating the minimal turbulent velocity scale (in the surface layer only)
                      ! 1: with a constant value
                      ! 2: with a stability dependent correction
  imode_charpar =2    ! mode of estimating the Charnock-Parameter
                      ! 1: use a constant value 
                      ! 2: use a wind-dependent value with a constant lower bound
                      ! 3: as 2, but with reduction at wind speeds above 25 m/s for more realistic TC wind speeds
  !$ACC DECLARE COPYIN(imode_charpar)

INTEGER :: &

  imode_syndiag =2, & ! mode of diagnostics at the synoptic near surface levels (related to 'itype_diag_t2m')
                      ! 1: direct interpolation of temperature and specific humidity
                      ! 2: interpol. of conserved quantities and subsequent statistical saturation adjustm.,
                      !    allowing particularly for the diagnostic of cloud water at the 2m-level (fog)
  imode_qvsatur =2, & ! mode of calculating the saturat. humidity
                      ! 1: old version using total pressure
                      ! 2: new version using partial pressure of dry air
  imode_stadlim =2, & ! mode of limitting statist. saturation adjustment (SUB 'turb_cloud')
                      ! 1: only absolut upper limit of stand. dev. of local oversatur. (sdsd)
                      ! 2: relative limit of sdsd and upper limit of cloud-water
  imode_trancnf =2, & ! mode of configuring the transfer-scheme (SUB 'turbtran')
                      ! 1: old version: start. with lamin. diffus.; with a lamin. correct. for profile-funct.;
                      !    interpol. T_s rather then Tet_l onto zero-level; calcul. only approx. Tet_l-grads.;
                      !    using an upper bound for TKE-forcing; without transmit. skin-layer depth to turbul.
                      ! 2: 1-st ConSAT: start. with estim. Ustar, without a laminar correct. for prof.-funct.;
                      !    interpol. Tet_l onto zero-level; calcul. Tet_l-gradients directly; 
                      !    without an upper bound for TKE-forcing; with transmit. skin-layer depth to turbul.
                      ! 3: 2-nd ConSAT: as "2", but with a hyperbolic interpol. of profile function
                      !    for stable stratification
                      ! 4: 3-rd ConSAT: as "3", but without using an upper interpolation node
  imode_tkediff =2, & ! mode of implicit TKE-Diffusion (related to 'c_diff')
                      ! 1: in terms of q=SQRT(2*TKE)) 
                      ! 2; in terms of TKE=0.5*TKE**2
  imode_adshear =2    ! mode of considering addit. shear by scale interaction (realt. to 'ltkesso', 'ltkeshs',
                      ! 'ltkecon')
                      ! 1: not consid. for stability functions
                      ! 2:  considered for stability functions

! Notice that the following selectors are provided by the parameter-list of 
! SUB 'turb_diffusion' or 'turb_transfer':

! iini                :type of initialization (0: no, 1: separate before the time loop
!                                                   , 2: within the first time step)
! itnd                :type of tendency cons. (0: no, 1: in implicit vertical diffusion equation
!                                                     2: by adding to current profile before vertical diffusion
!                                                     3: by using corrected virtual vertical profiles

!==============================================================================
! Declarations of utility variables:

! Turbulence parameters which are computed during model run
!-------------------------------------------------------------------------------

REAL (KIND=wp)     ::        &
  ! do we need it as TARGET?
  ! these variables are set in SR turb_param
  c_tke,tet_g,rim, &
  c_m,c_h, b_m,b_h,  sm_0, sh_0, &
  d_0,d_1,d_2,d_3,d_4,d_5,d_6, &
  a_3,a_5,a_6,                 &

  ! these parameters are physical constants which are either taken as
  ! they are or set to 0.0 in the turbulence for special applications
  tur_rcpv,         & ! cp_v/cp_d - 1
  tur_rcpl            ! cp_l/cp_d - 1 (where cp_l=cv_l)

! Definition of used data types
!-------------------------------

TYPE modvar !model variable
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
     REAL (KIND=wp), POINTER, CONTIGUOUS     ::         &
#else
     REAL (KIND=wp), POINTER                 ::         &
#endif
             av(:,:) => NULL(), & !atmospheric values
             sv(:)   => NULL(), & !surface     values (concentration of flux density)
             at(:,:) => NULL()    !atmospheric time tendencies
     LOGICAL                                 ::         &
             fc                !surface values are flux densities
END TYPE modvar

TYPE turvar !turbulence variables
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
     REAL (KIND=wp), POINTER, CONTIGUOUS     ::         &
#else
     REAL (KIND=wp), POINTER                 ::         &
#endif
             tkv(:,:) => NULL(), & !turbulent coefficient for vert. diff.
             tsv(:)   => NULL()    !turbulent velocity at the surface
END TYPE turvar

TYPE varprf !variable profile
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
     REAL (KIND=wp), POINTER, CONTIGUOUS     ::         &
#else
     REAL (KIND=wp), POINTER                 ::         &
#endif
             bl(:,:), & !variable at boundary model levels
             ml(:,:)    !variable at main     model levels
END TYPE varprf


! 7. Switches from the COSMO-Model: must only be defined for ICON
! ---------------------------------------------------------------

REAL(KIND=wp), POINTER :: &
  impl_weight(:)  ! implicit weights for tridiagonal solver

! Switches controlling turbulent diffusion:
! ------------------------------------------

INTEGER  :: &

    itype_tran   =2,       & ! type of surface-atmosphere transfer
    imode_circ   =2,       & ! mode of treating the circulation term
    ilow_dcond   =1          ! type of the default condition at the lower boundary

!==============================================================================

CONTAINS

!==============================================================================

! MR: Uli, the following SUB 'get_turbdiff' should be ported to the ICON-interface 
!     for TURBDIFF!!

! this subroutines sets the switches defined above during an ICON run

SUBROUTINE get_turbdiff_param (jg)

   INTEGER, INTENT(IN) :: jg !patch index

   impl_weight => turbdiff_config(jg)%impl_weight

   imode_tran   = turbdiff_config(jg)%imode_tran
   icldm_tran   = turbdiff_config(jg)%icldm_tran
   imode_turb   = turbdiff_config(jg)%imode_turb
   icldm_turb   = turbdiff_config(jg)%icldm_turb
   itype_sher   = turbdiff_config(jg)%itype_sher
   imode_frcsmot= turbdiff_config(jg)%imode_frcsmot
   imode_tkesso = turbdiff_config(jg)%imode_tkesso

   ltkesso      = turbdiff_config(jg)%ltkesso
   ltkecon      = turbdiff_config(jg)%ltkecon
   ltkeshs      = turbdiff_config(jg)%ltkeshs
   lexpcor      = turbdiff_config(jg)%lexpcor
   ltmpcor      = turbdiff_config(jg)%ltmpcor
   lprfcor      = turbdiff_config(jg)%lprfcor
   lnonloc      = turbdiff_config(jg)%lnonloc
   lfreeslip    = turbdiff_config(jg)%lfreeslip
   lcpfluc      = turbdiff_config(jg)%lcpfluc
   lsflcnd      = turbdiff_config(jg)%lsflcnd

   itype_wcld   = turbdiff_config(jg)%itype_wcld

   tur_len      = turbdiff_config(jg)%tur_len
   pat_len      = turbdiff_config(jg)%pat_len
   a_stab       = turbdiff_config(jg)%a_stab
   tkhmin       = turbdiff_config(jg)%tkhmin
   tkmmin       = turbdiff_config(jg)%tkmmin
   tkhmin_strat = turbdiff_config(jg)%tkhmin_strat
   tkmmin_strat = turbdiff_config(jg)%tkmmin_strat

   imode_charpar = turbdiff_config(jg)%imode_charpar
   alpha0       = turbdiff_config(jg)%alpha0
   alpha0_max   = turbdiff_config(jg)%alpha0_max
   alpha0_pert  = turbdiff_config(jg)%alpha0_pert

   c_diff       = turbdiff_config(jg)%c_diff
   rlam_heat    = turbdiff_config(jg)%rlam_heat
   rlam_mom     = turbdiff_config(jg)%rlam_mom
   rat_sea      = turbdiff_config(jg)%rat_sea
   tkesmot      = turbdiff_config(jg)%tkesmot
   frcsmot      = turbdiff_config(jg)%frcsmot
   impl_s       = turbdiff_config(jg)%impl_s
   impl_t       = turbdiff_config(jg)%impl_t
   q_crit       = turbdiff_config(jg)%q_crit

   loutshs      = ltkeshs .OR. itype_sher > 0

END SUBROUTINE get_turbdiff_param


!==============================================================================

END MODULE turb_data
