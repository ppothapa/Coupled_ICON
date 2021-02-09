!>
!! Computation of cloud cover and grid mean cloud liquid water and cloud ice
!!
!! This routine takes information from turbulence, convection and grid-scale
!! to produce cloud properties used in radiation (and microphysics).
!!
!! Possible future options
!! - simple diagnostic (from turbulence, convection and grid scale)
!! - prognostic total water variance AND prognostic ice
!!
!!
!! @author Martin Koehler, DWD
!!
!!
!! @par Revision History
!! Developed by Martin Koehler  (starting 2010-08-23)
!! Modification by Martin Koehler, DWD (2010-11-20)
!! - options 0,1,2,3 ready
!!
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
#include "consistent_fma.inc"
MODULE mo_cover_koe

  USE mo_kind,               ONLY: wp, i4

  USE mo_physical_constants, ONLY: rdv    , & !! r_d / r_v
                                   rv     , & !! Rv
                                   tmelt  , & !! melting temperature of ice/snow
                                   grav   , & !! gravitational acceleration
                                   alvdcp     !! lh_v / cp_d

  USE mo_convect_tables,     ONLY: c1es   , & !! constants for computing the sat. vapour
                                   c3les  , & !! pressure over water (l) and ice (i)
                                   c3ies  , & !!               -- " --
                                   c4les  , & !!               -- " --
                                   c4ies, c5les

  USE mo_cufunctions,        ONLY: foealfcu    !! liquid fraction as in Tiedtke/Bechtold convection

  USE mo_cloud_diag,         ONLY: cloud_diag

  USE mo_cover_cosmo,        ONLY: cover_cosmo

  USE mo_impl_constants,     ONLY: iedmf

  USE mo_nwp_tuning_config,  ONLY: tune_box_liq, tune_box_liq_asy, tune_thicklayfac, tune_sgsclifac, icpl_turb_clc

  USE mo_ensemble_pert_config, ONLY: box_liq_sv, thicklayfac_sv, box_liq_asy_sv

  USE mo_impl_constants,      ONLY: max_dom  

  IMPLICIT NONE

  PRIVATE


  PUBLIC :: cover_koe, t_cover_koe_config, cover_koe_config


!  Cloud cover derived type with physics configuration options
   
  TYPE t_cover_koe_config
    INTEGER(KIND=i4)        ::     icldscheme    ! cloud cover option 
    INTEGER(KIND=i4)        ::     inwp_turb     ! turbulence scheme number
    INTEGER(KIND=i4)        ::     inwp_cpl_re   ! coupling reff (for qs altering qi)
    INTEGER(KIND=i4)        ::     inwp_reff     ! reff option (for qs altering qi)
  END TYPE t_cover_koe_config

!-------------------------------------------------------------------------

  TYPE(t_cover_koe_config)  ::     cover_koe_config(max_dom)  ! This should be in interface module

CONTAINS

!-------------------------------------------------------------------------
!
!  Cloud cover and cloud water/ice calculation.
!
!
!  @par Revision History
!  Initial version by Martin Koehler, DWD (2010-08-23)
!
!  Options:
!  (0) no clouds
!  (1) diagnostic cloud cover
!  (2) prognostic total water variance (not yet started)
!  (3) clouds as in COSMO
!  (4) clouds as in turbulence
!  (5) grid-scale cloud cover [1 or 0]
!
!-------------------------------------------------------------------------

SUBROUTINE cover_koe( &
  & kidia, kfdia, klon, kstart, klev, & ! in:    dimensions (turn off physics above kstart)
  & cover_koe_config                , & ! in:    configure state
  & tt                              , & ! in:    temperature (main levels)
  & pp                              , & ! in:    pressure (")
  & ps                              , & ! in:    surface pressure
  & t_g                             , & ! in;    surface temperature
  & pgeo, deltaz                    , & ! in:    geopotential above ground, layer thickness
  & rho                             , & ! in:    density
  & rcld                            , & ! inout: standard deviation of saturation deficit
  & ldland                          , & ! in:    land/sea mask
  & ldcum, kcbot, kctop, ktype      , & ! in:    convection: on/off, bottom, top, type
  & pmfude_rate                     , & ! in:    convection: updraft detrainment rate
  & plu                             , & ! in:    convection: updraft condensate
  & rhoc_tend                       , & ! in:    convective rhoc tendency
  & qv, qc, qi, qs, qtvar           , & ! inout: prognostic cloud variables
  & cc_tot, qv_tot, qc_tot, qi_tot )    ! out:   cloud output diagnostic


!! Subroutine arguments:
!! --------------------

INTEGER(KIND=i4), INTENT(IN) ::  &
  & kidia            , & ! horizontal start index
  & kfdia            , & ! horizontal end   index
  & klon             , & ! horizontal dimension
  & kstart           , & ! vertical start index (turn off physics above)
  & klev                 ! vertical dimension

TYPE(t_cover_koe_config), INTENT(IN) :: cover_koe_config ! configure state

REAL(KIND=wp), DIMENSION(klon,klev), INTENT(IN) ::  &
  & pp               , & ! full pressure                                 (  Pa )
  & pgeo             , & ! geopotential (above ground)                   (m2/s2)
  & deltaz           , & ! layer thickness                               (m)
  & rho                  ! density                                       (kg/m3)

REAL(KIND=wp), DIMENSION(klon,klev), INTENT(IN) ::  &
  & tt               , & ! temperature                                   (  K  )
  & qv               , & ! specific water vapor content                  (kg/kg)
  & qc               , & ! specific cloud water content                  (kg/kg)
  & qi               , & ! specific cloud ice   content                  (kg/kg)
  & qs               , & ! specific snow        content                  (kg/kg)
  & qtvar                ! total water variance (qt'2)                   (kg2/kg2)

REAL(KIND=wp), DIMENSION(klon), INTENT(IN) ::  &
  & ps               , & ! surface pressure
  & t_g                  ! surface temperature

REAL(KIND=wp), DIMENSION(klon,klev+1), INTENT(INOUT) ::  &
  & rcld                 ! standard deviation of saturation deficit

LOGICAL, DIMENSION(klon), INTENT(IN) ::  &
  & ldland           , & ! true for land points
  & ldcum                ! true for convection points

INTEGER(KIND=i4), DIMENSION(klon), INTENT(IN) ::  &
  & ktype            , & ! convection type
  & kcbot            , & ! convective cloud base level (klev: bottom level, -1: no conv)
  & kctop                ! convective cloud top level

REAL(KIND=wp), DIMENSION(klon,klev), INTENT(IN) ::  &
  & pmfude_rate      , & ! convective updraft detrainment rate           (kg/(m3*s))
  & plu              , & ! updraft condensate                            (kg/kg)
  & rhoc_tend            ! convective rho_c tendency                     (kg/(m3*s))

REAL(KIND=wp), DIMENSION(klon,klev), INTENT(INOUT) ::   &
  & cc_tot           , & ! cloud cover diagnostic
  & qv_tot           , & ! specific water vapor content diagnostic       (kg/kg)
  & qc_tot           , & ! specific cloud water content diagnostic       (kg/kg)
  & qi_tot               ! specific cloud ice   content diagnostic       (kg/kg)


!! Local variables:
!! ----------------

LOGICAL ::  &
  & lprog_qi, l_addsnow

INTEGER (KIND=i4) :: &
  & jl, jk, jkp1,         &
  & itype_wcld, icldm_rad

REAL(KIND=wp), DIMENSION(klon,klev)  :: &
  & cc_turb, qc_turb, qi_turb, &
  & cc_conv, qc_conv, qi_conv, &
  & cc_turb_liq, cc_turb_ice , &
  & p0


REAL(KIND=wp) :: &
  & fgew   , fgee   , fgqs   , dqsdt,   & !fgqv   , & ! name of statement functions
  & ztt    , zzpv   , zzpa   , zzps   , zqs, &
  & zf_ice , deltaq , qisat_grid, zdeltaq, zrcld, thicklay_fac, tfac, satdef_fac, rhcrit_sgsice, &
  & vap_pres, zaux, zqisat_m50, zqisat_m25, qi_mod, par1, qcc, box_liq_asy, fac_aux, fac_sfc, &
  & rcld_asyfac, dq1, dq2, dq3, tfmax

REAL(KIND=wp), DIMENSION(klon,klev)  :: &
  zqlsat , zqisat, zagl_lim, zdqlsat_dT

!! Local parameters:
!! -----------------

REAL(KIND=wp), PARAMETER  :: &
  & zcldlim  = 1.0e-8_wp, & ! threshold of cloud water/ice for cloud cover  (kg/kg)
  & taudecay = 1500.0_wp, & ! decay time scale of convective anvils
  & box_ice  = 0.05_wp  , & ! box width scale ice clouds
  & tm10     = tmelt - 10.0_wp, &
  & tm40     = tmelt - 40.0_wp

  REAL(kind=wp), PARAMETER :: grav_i = 1._wp/grav

!-----------------------------------------------------------------------

! statement function to calculate saturation vapour pressure over water
  fgew(ztt)            = c1es * EXP( c3les*(ztt - tmelt)/(ztt - c4les) )  ! ztt: temperature

! statement function to calculate saturation vapour pressure over ice
  fgee(ztt)            = c1es * EXP( c3ies*(ztt - tmelt)/(ztt - c4ies) )  ! ztt: temperature

! statement function to calculate specific humitdity
! fgqv(zzpv,zzpa)      = rdv * zzpv / (zzpa - (1._wp-rdv)*zzpv)           ! zzpv: vapour pressure

! statement function to calculate saturation specific humidities from RH=esat/e (proper and safe)
  fgqs(zzps,zzpv,zzpa) = rdv * zzps / (zzpa - (1._wp-rdv)*zzpv)           ! zzps: saturation vapour pressure

! statement function for dq_sat_dT
  dqsdt(ztt,zqs) = c5les * (1._wp-zqs) * zqs / (ztt-c4les)**2
!-----------------------------------------------------------------------

! saturation mixing ratio at -50 C and 200 hPa
zqisat_m50 = fgqs ( fgee(223.15_wp), 0._wp, 20000._wp )

! saturation mixing ratio at -25 C and 700 hPa
zqisat_m25 = fgqs ( fgee(248.15_wp), 0._wp, 70000._wp )

IF (icpl_turb_clc == 1) THEN
  rcld_asyfac = 0._wp
ELSE
  rcld_asyfac = 2._wp
ENDIF
! auxiliary factors depending on ensemble perturbations in order to increase spread
tfmax = 0.6_wp  + 100._wp*(tune_thicklayfac-thicklayfac_sv)
dq1   = 0.8_wp  + 100._wp*(tune_box_liq-box_liq_sv)*(tune_box_liq_asy-box_liq_asy_sv)
dq2   = 1._wp   + 400._wp*(tune_thicklayfac-thicklayfac_sv)*(tune_box_liq_asy-box_liq_asy_sv)
dq3   = 0.25_wp + 2500._wp*(tune_box_liq-box_liq_sv)*(tune_thicklayfac-thicklayfac_sv)


! Snow is added to qi_dia in three cases: 
! 1) No coupling of reff with radiation
! 2) No param for reff 
! 3) Using the original RRTM parameterization for reff 
l_addsnow = (cover_koe_config%inwp_cpl_re == 0) .OR. (cover_koe_config%inwp_reff == 0) .OR. &
            (cover_koe_config%inwp_reff == 101)

! Set cloud fields for stratospheric levels to zero
DO jk = 1,kstart-1
  DO jl = kidia,kfdia
    qv_tot(jl,jk) = qv(jl,jk)
    qc_tot(jl,jk) = 0.0_wp
    qi_tot(jl,jk) = 0.0_wp
    cc_tot(jl,jk) = 0.0_wp
  ENDDO
ENDDO

!-----------------------------------------------------------------------
! Calculate water vapour saturation mixing ratios of over water
! and over ice (from mo_cover_cosmo.f90)
!-----------------------------------------------------------------------

DO jk = kstart,klev
  DO jl = kidia,kfdia
    vap_pres = qv(jl,jk) * rho(jl,jk) * rv * tt(jl,jk)
    ! specific humidity at saturation over water (zqlsat) and ice (zqisat)
    zqlsat (jl,jk) = fgqs ( fgew(tt(jl,jk)), vap_pres, pp(jl,jk) )
    zqisat (jl,jk) = fgqs ( fgee(tt(jl,jk)), vap_pres, pp(jl,jk) )
    ! derivative of qsat_w w.r.t. temperature
    zdqlsat_dT(jl,jk) = dqsdt(tt(jl,jk), zqlsat(jl,jk))
    ! limit on box width near the surface, reaches 0.05 at 500 m AGL
    zagl_lim(jl,jk) = 0.025_wp + 5.e-5_wp * pgeo(jl,jk)*grav_i
  ENDDO
ENDDO

!-----------------------------------------------------------------------
! Select desired cloud cover framework
!-----------------------------------------------------------------------

SELECT CASE( cover_koe_config%icldscheme )

!-----------------------------------------------------------------------

! no clouds
CASE( 0 )

  DO jk = kstart,klev
    DO jl = kidia,kfdia
      qc_tot(jl,jk) = 0.0_wp
      qi_tot(jl,jk) = 0.0_wp
      cc_tot(jl,jk) = 0.0_wp
    ENDDO
  ENDDO

!-----------------------------------------------------------------------

! diagnostic cloud cover
CASE( 1 )

  DO jk = kstart,klev
    jkp1 = MIN(jk+1,klev)
    DO jl = kidia,kfdia

! stratiform cloud
!  liquid cloud
     !
     ! quadratic increase of cloud cover from 0 to 1 between RH = (1 - box_liq_asy*tune_box_liq) and (1 + tune_box_liq);
     ! the box width is defined as 4*rcld (standard dev. of saturation deficit) with upper and lower limits;
     ! in addition, sub-grid scale moisture variations in the vertical are parameterized depending on vertical resolution
     ! Diagnosed cloud water is proportional to clcov**2
     !
      thicklay_fac = MIN(tfmax,MAX(0._wp,tune_thicklayfac*(deltaz(jl,jk)-150._wp))) ! correction for thick model layers
      zdeltaq = MIN(tune_box_liq*(1._wp+0.5_wp*thicklay_fac), zagl_lim(jl,jk)) * zqlsat(jl,jk)
      zrcld = 0.5_wp*(rcld(jl,jk)+rcld(jl,jk+1))
      IF (icpl_turb_clc == 1) THEN
        deltaq = MAX(dq1*zdeltaq,(4._wp*dq2+thicklay_fac)*zrcld)
      ELSE
        deltaq = dq1*zdeltaq+(dq2+dq3*thicklay_fac)*zrcld
      ENDIF
      deltaq = MIN(deltaq,2._wp*zdeltaq)
      IF ( ( qv(jl,jk) + qc(jl,jk) - deltaq ) > zqlsat(jl,jk) ) THEN
        cc_turb_liq(jl,jk) = 1.0_wp
        qc_turb  (jl,jk)   = qv(jl,jk) + qc(jl,jk) - zqlsat(jl,jk)
      ELSE
        ! asymmetry factor for water clouds and derived parameters;
        ! the asymmetry factor is reduced to 2 close to the surface because 'long tails' are unrealistic (and detrimental) in this case
        fac_sfc = MAX(0._wp,zqlsat(jl,jk)-(qv(jl,jk)+deltaq))/deltaq*MAX(0._wp,(tune_box_liq-zagl_lim(jl,jk))/zagl_lim(jl,jk))
        fac_sfc = MIN(1._wp,fac_sfc)
        box_liq_asy = tune_box_liq_asy*(1._wp+0.5_wp*thicklay_fac+rcld_asyfac*zrcld/zqlsat(jl,jk))*(1._wp-fac_sfc) + 2._wp*fac_sfc
        par1 = box_liq_asy+1._wp
        !
        zaux = qv(jl,jk) + qc(jl,jk) + box_liq_asy*deltaq - zqlsat(jl,jk)
        cc_turb_liq(jl,jk) = SIGN((zaux/(par1*deltaq))**2,zaux)
        ! compensating reduction of cloud water content if the thick-layer correction is active
        fac_aux = 1._wp + (alvdcp*zdqlsat_dT(jl,jk)+thicklay_fac)*MIN(1._wp,2._wp*(1._wp-cc_turb_liq(jl,jk)))
        IF ( cc_turb_liq(jl,jk) > 0.0_wp ) THEN
          qc_turb  (jl,jk) = deltaq*cc_turb_liq(jl,jk)**2/fac_aux
        ELSE
          qc_turb  (jl,jk) = 0.0_wp
        ENDIF
      ENDIF

!  ice cloud
      rhcrit_sgsice = 1._wp - 0.25_wp*tune_sgsclifac*MAX(0._wp,0.75_wp-zqisat(jl,jk)/zqlsat(jl,jk))
      fac_aux = 1._wp - MIN(1._wp,MAX(0._wp,tt(jl,jk)-tm40)/15._wp)

      qi_mod = MERGE( MAX(qi(jl,jk), 0.1_wp*(qi(jl,jk)+qs(jl,jk))), qi(jl,jk), l_addsnow) + &
               fac_aux*MIN(1._wp,tune_sgsclifac*zrcld/(box_ice*zqisat(jl,jk)))* &
               MAX(0._wp,qv(jl,jk)-rhcrit_sgsice*zqisat(jl,jk))

     !ice cloud: assumed box distribution, width 0.1 qisat, saturation above qv 
     !           (qv is microphysical threshold for ice as seen by grid scale microphysics)
      IF ( qi_mod > zcldlim ) THEN
        deltaq     = box_ice * MIN(zqisat_m25, zqisat(jl,jk))  ! box width = 2*deltaq
        qisat_grid = MAX( qv(jl,jk), zqisat(jl,jk) )           ! qsat grid-scale
        IF ( ( qv(jl,jk) + qi_mod - deltaq) > qisat_grid ) THEN
          cc_turb_ice(jl,jk) = 1.0_wp
          qi_turb    (jl,jk) = qi_mod
        ELSE
          zaux = MIN(2._wp*(qi_mod + MAX(0._wp, qv(jl,jk) - zqisat(jl,jk))), &
                 qv(jl,jk) + qi_mod + deltaq - qisat_grid )
          cc_turb_ice(jl,jk) = zaux / (2._wp*deltaq)
          IF ( cc_turb_ice(jl,jk) > 0.0_wp ) THEN
            qi_turb  (jl,jk) = MAX(MIN(qi_mod,zqisat(jl,jk)*cc_turb_ice(jl,jk)), zaux**2/(4._wp*deltaq))
          ELSE
            qi_turb  (jl,jk) = 0.0_wp
          ENDIF
        ENDIF
      ELSE
        cc_turb_ice(jl,jk) = 0.0_wp
        qi_turb    (jl,jk) = 0.0_wp
      ENDIF

      ! reduce cloud cover fraction of very thin ice clouds, defined as clouds with a mixing ratio
      ! of less than 5% of the saturation mixing ratio w.r.t. ice at -50 deg C
      cc_turb_ice(jl,jk) = MIN(cc_turb_ice(jl,jk),qi_turb(jl,jk)/(box_ice*zqisat_m50))

      cc_turb(jl,jk) = max( cc_turb_liq(jl,jk), cc_turb_ice(jl,jk) )          ! max overlap liq/ice
      cc_turb(jl,jk) = min(max(0.0_wp,cc_turb(jl,jk)),1.0_wp)
      qc_turb(jl,jk) =     max(0.0_wp,qc_turb(jl,jk))
      qi_turb(jl,jk) =     max(0.0_wp,qi_turb(jl,jk))


! convective cloud
      tfac = foealfcu(tt(jl,jk)) ! foealfa = liquid/(liquid+ice); controls partitioning between cloud water and cloud ice
      ! reduction of decay time scale depending on saturation deficit
      satdef_fac = 1._wp - MIN(0.9_wp,125._wp*( tfac*(zqlsat(jl,jk)-qv(jl,jk)) + (1._wp-tfac)*(zqisat(jl,jk)-qv(jl,jk)) ))
      cc_conv(jl,jk) = ( pmfude_rate(jl,jk) / rho(jl,jk) ) &  ! cc = detrainment/rho / (Du/rho + 1/tau,decay)
                   & / ( pmfude_rate(jl,jk) / rho(jl,jk) + 1.0_wp / (taudecay*satdef_fac) )
      ! detrainment water is defined as detrainment rate * updraft liquid water in layer below
      qc_conv(jl,jk) = cc_conv(jl,jk) * plu(jl,jkp1)*       tfac ! ql up
      qi_conv(jl,jk) = cc_conv(jl,jk) * plu(jl,jkp1)*(1._wp-tfac)! qi up
      ! alternative formulation of source term for liquid convective clouds depending on detrained cloud water and RH;
      ! as most important difference, it uses the same clcov-qc relationship as turbulent clouds but is restricted to low mixing ratios
      qcc = MAX(0._wp, MIN(0.075_wp*tune_box_liq*zqlsat(jl,jk), (rhoc_tend(jl,jk)/rho(jl,jk))*taudecay* &
        (1._wp - 4._wp*(1._wp-qv(jl,jk)/zqlsat(jl,jk))) ))
      cc_conv(jl,jk) = MAX(cc_conv(jl,jk),SQRT(qcc/(tune_box_liq*zqlsat(jl,jk))) )
      qc_conv(jl,jk) = MAX(qcc,qc_conv(jl,jk))
      !
      cc_conv(jl,jk) = min(max(0.0_wp,cc_conv(jl,jk)),1.0_wp)
      qc_conv(jl,jk) = min(max(0.0_wp,qc_conv(jl,jk)),0.1_wp*qv(jl,jk))       ! qc limit to 10%qv
      qi_conv(jl,jk) = min(max(0.0_wp,qi_conv(jl,jk)),0.1_wp*qv(jl,jk))       ! qi limit to 10%qv

    ENDDO
  ENDDO


  IF (cover_koe_config%inwp_turb == iedmf) THEN
    DO jk = kstart,klev
      DO jl = kidia,kfdia
!       IF (SQRT(qtvar(jl,jk)) / MAX(qv(jl,jk)+qc(jl,jk)+qi(jl,jk),0.000001_wp) > 0.01_wp) THEN
        IF (SQRT(qtvar(jl,jk)) > 0.001_wp * (qv(jl,jk)+qc(jl,jk)+qi(jl,jk)) ) THEN
! for EDMF DUALM: take values written within EDMF - only done within high qtvar grid points
          cc_tot(jl,jk) = cc_tot(jl,jk)
          qc_tot(jl,jk) = qc_tot(jl,jk)
          qi_tot(jl,jk) = qi_tot(jl,jk)
        ELSE
! combination strat/conv cloud
          cc_tot(jl,jk)  = max( cc_turb(jl,jk), cc_conv(jl,jk) )
          qc_tot(jl,jk)  = max( qc_turb(jl,jk), qc_conv(jl,jk) )
          qi_tot(jl,jk)  = max( qi_turb(jl,jk), qi_conv(jl,jk) )
        ENDIF
      ENDDO
    ENDDO
  ELSE ! use always combination of strat/conv cloud
    DO jk = kstart,klev
      DO jl = kidia,kfdia
        cc_tot(jl,jk)  = max( cc_turb(jl,jk), cc_conv(jl,jk) )
        qc_tot(jl,jk)  = max( qc_turb(jl,jk), qc_conv(jl,jk) )
        qi_tot(jl,jk)  = max( qi_turb(jl,jk), qi_conv(jl,jk) )
      ENDDO
    ENDDO
  ENDIF


!-----------------------------------------------------------------------

! prognostic total water variance
CASE( 2 )

  cc_tot = 0.0_wp
  qc_tot = 0.0_wp
  qi_tot = 0.0_wp

!-----------------------------------------------------------------------

! clouds as in COSMO
CASE( 3 )

  lprog_qi   = .true.       ! .true.: running with cloud ice
  icldm_rad  = 3            ! 3:      standard COSMO/GME
  itype_wcld = 2            ! cloud calculation in cloud_diag (only active if used)
  p0         = 0.0_wp       ! base state presssure set to zero, pp is full pressure
  !DO jl = kidia,kfdia
  !  t_g(jl)  = tt(jl,klev)  ! should be surface temperature!
  !ENDDO

  CALL cover_cosmo ( &
    klon       ,    & ! number of grid points in zonal direction
    klev       ,    & ! number of grid points in vertical direction
    kidia      ,    & ! start- and end-indices for computing the radiation (i-start)
    kfdia      ,    & !    -"- (i-end)
    kstart     ,    & ! vertical start index
    tt         ,    & ! temperature                                   (  K  )
    ps         ,    & ! surface pressure                              ( Pa  )
    p0         ,    & ! reference pressure at full levels             ( Pa  )
    pp         ,    & ! deviation from the reference pressure         ( Pa  )
    pgeo       ,    & ! geopotential                                  (m2/s2)
    qv         ,    & ! specific water vapor content                  (kg/kg)
    qc         ,    & ! specific cloud water content                  (kg/kg)
    qi         ,    & ! specific cloud ice content                    (kg/kg)
    rcld       ,    & ! standard deviation of saturation deficit
    t_g        ,    & ! weighted surface temperature                  (  K  )
    cc_turb    ,    & ! OUT: subgrid-scale stratiform cloud cover     (  1  )
    cc_conv    ,    & ! OUT: cloud cover due to convection            (  1  )
    icldm_rad  ,    & ! mode of cloud representation in radiation parametr.
    itype_wcld ,    & ! type of water cloud diagnosis
    lprog_qi   ,    & ! if .TRUE., running with cloud ice
    ldland     ,    & ! landpoint mask
    ldcum      ,    & ! in:  convection on/off
    kcbot      ,    & ! in:  convective cloud base
    kctop      ,    & ! in:  convective cloud top
    qv_tot     ,    & ! OUT: water vapour mixing ratio
    cc_tot     ,    & ! OUT: cloud cover in each layer
    qc_tot     ,    & ! OUT: liquid water mixing ratio
    qi_tot     )      ! OUT: ice mixing ratio

 !cc_tot = cc_turb   ! CC=turb for testing
 !cc_tot = cc_conv   ! CC=conv for testing


!-----------------------------------------------------------------------

! cloud cover as in turbulence
CASE( 4 )

  itype_wcld = 2     ! 2: Gaussian calculation; 1: Sundqvist type
  p0         = 0.0_wp! base state presssure set to zero, pp is full pressure

  CALL cloud_diag ( cc_tot, qc_tot,                   &
                    kidia, kfdia, 1, 1, kstart, klev, &
                    klon , 1 , klev,                  &
                    tt, qv, qc, pp, p0, rcld, ps,     &
                    itype_wcld )

  qi_tot     = 0.0_wp

!-----------------------------------------------------------------------

! grid-scale cloud cover [1 or 0]
CASE( 5 )

  DO jk = kstart,klev
    DO jl = kidia,kfdia
      IF ( qc(jl,jk) + qi(jl,jk) > zcldlim ) THEN
        cc_tot(jl,jk) = 1.0_wp
      ELSE
        cc_tot(jl,jk) = 0.0_wp
      ENDIF
      qc_tot(jl,jk) = qc(jl,jk)
      qi_tot(jl,jk) = qi(jl,jk)
    ENDDO
  ENDDO

!-----------------------------------------------------------------------

END SELECT


! total water vapor by conservation of grid-scale total water

!PREVENT_INCONSISTENT_IFORT_FMA
DO jk = kstart,klev
  DO jl = kidia,kfdia
    qv_tot(jl,jk) = qv(jl,jk) + qc(jl,jk) + qi(jl,jk) - qc_tot(jl,jk) - qi_tot(jl,jk)

!    IF ( (qc_tot(jl,jk)+qi_tot(jl,jk)) > 0.0_wp ) THEN
!      zf_ice = MIN( 1.0_wp, MAX( 0.0_wp, qi_tot(jl,jk) / (qc_tot(jl,jk)+qi_tot(jl,jk)) ) )
!    ELSE
!      zf_ice = 0.0_wp
!    ENDIF

! GZ: except for the no-clouds option (0), it is impossible to violate this sanity check criterion;
!     in this case, however, we'd like to have really no clouds
! sanity check 1: qv_tot < 1.1 q,sat,liq (to allow for small qc_tot diagnostic)
!    IF (qv_tot(jl,jk) .GT. 1.1_wp * zqlsat(jl,jk) ) THEN
!      qc_tot(jl,jk) = qc_tot(jl,jk) + ( qv_tot(jl,jk) - 1.1_wp * zqlsat (jl,jk) ) *(1.0_wp-zf_ice)
!      qi_tot(jl,jk) = qi_tot(jl,jk) + ( qv_tot(jl,jk) - 1.1_wp * zqlsat (jl,jk) ) * zf_ice
!      qv_tot(jl,jk) = 1.1_wp * zqlsat (jl,jk)
!    ENDIF

! sanity check 2: qv_tot > 0.1 qv (take qc, qi if available)
    IF ( qv_tot(jl,jk) < 0.1_wp * qv(jl,jk) ) THEN
      zf_ice = qi_tot(jl,jk) / MAX(zcldlim,qc_tot(jl,jk)+qi_tot(jl,jk))
      qc_tot(jl,jk) = MAX(0._wp, qc_tot(jl,jk) + ( qv_tot(jl,jk) - 0.1_wp * qv(jl,jk) ) * (1.0_wp - zf_ice) )
      qi_tot(jl,jk) = MAX(0._wp, qi_tot(jl,jk) + ( qv_tot(jl,jk) - 0.1_wp * qv(jl,jk) ) * zf_ice )
      qv_tot(jl,jk) = 0.1_wp * qv(jl,jk)
    ENDIF

! sanity check 3: remove 'epsilon' clouds
    IF ( qc_tot(jl,jk) + qi_tot(jl,jk) < zcldlim ) THEN
      qv_tot(jl,jk) = qv_tot(jl,jk) + qc_tot(jl,jk) + qi_tot(jl,jk)
      qc_tot(jl,jk) = 0.0_wp
      qi_tot(jl,jk) = 0.0_wp
      cc_tot(jl,jk) = 0.0_wp
    ENDIF

  ENDDO
ENDDO

END SUBROUTINE cover_koe


END MODULE mo_cover_koe
