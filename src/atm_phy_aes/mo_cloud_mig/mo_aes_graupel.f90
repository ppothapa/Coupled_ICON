!>
!! @brief Module containing thermodynamic functions used by the AES department in MPI-M
!!
!! @contact B. Stevens, MPI-M, 2022-06, clode slightly refactored from gspc_graupel which
!! is owned and maintained by A. Seifert, DWD, and has older origins as indicated in the
!! original documentation.
!!
!! @par Copyright and License
!!
!! Copyright 2022 Max Planck Institute for Meteorology
!! This program is free software: you can restribute it and/or modify it under the terms of the
!! BSD-3-clause license
!!

MODULE mo_aes_graupel

USE mo_kind,               ONLY: wp         , &
                                 i4
USE mo_physical_constants, ONLY: alv   , & !! latent heat of vapourization
                                 als   , & !! latent heat of sublimation
                                 rv    , & !! vapor gas constant
                                 cpv   , & !! isobaric specific heat of vapor
                                 cvd   , & !! isometric specific heat of dry air
                                 cvv   , & !! isometric specific heat of vapor
                                 clw   , & !! specific heat of liquid water
                                 tmelt     !! melting temperature of ice/snow

USE mo_aes_thermo,         ONLY: sat_pres_water,      &  !! saturation vapor pressure over water
                                 vaporization_energy, &
                                 sublimation_energy,  &
                                 qsat_rho,            &  !! sat_pres_water (rho constant)
                                 qsat_ice_rho            !! sat_pres_ice (rho constant)
USE mo_exception,          ONLY: message, message_text

USE gscp_data, ONLY:                                                     &
    ccsrim,    ccsagg,    ccsdep,    ccsvel,    ccsvxp,    ccslam,       &
    ccslxp,    ccsaxp,    ccsdxp,    ccshi1,    ccdvtp,    ccidep,       &
    ccswxp,    zconst,    zcev,      zbev,      zcevxp,    zbevxp,       &
    zvzxp,     zvz0r,                                                    &
    v0snow,                                                              &
    x13o8,     x1o2,      x27o16,    x3o4,      x7o4,      x7o8,         &
    zbvi,      zcac,      zccau,     zciau,     zcicri,                  &
    zcrcri,    zcrfrz,    zcrfrz1,   zcrfrz2,   zeps,      zkcac,        &
    zkphi1,    zkphi2,    zkphi3,    zmi0,      zmimax,    zmsmin,       &
    zn0s0,     zn0s1,     zn0s2,     znimax_thom,          zqmin,        &
    zrho0,     zthet,     zthn,      ztmix,     ztrfrz,                  &
    zvz0i,     icesedi_exp,          zams => zams_gr,                    &
    iautocon,  isnow_n0temp, dist_cldtop_ref,   reduce_dep_ref,          &
    tmin_iceautoconv,     zceff_fac, zceff_min,                          &
    mma, mmb, v_sedi_rain_min, v_sedi_snow_min, v_sedi_graupel_min

IMPLICIT NONE
PRIVATE
PUBLIC :: graupel

LOGICAL, PARAMETER :: &
  lrain        = .TRUE.  , & ! switch for disabling rain
  lautocnv     = .TRUE.  , & ! switch for disabling auto-conversion of cloud water
  lcold        = .TRUE.      ! switch for disabling freezing processes

REAL(wp), PARAMETER :: &
   ci = 2108._wp ! specific heat of ice


CONTAINS

#ifdef _OPENACC
! GPU code can't flush to zero double precision denormals
! So to avoid CPU-GPU differences we'll do it manually
FUNCTION make_normalized(v)
  !$ACC ROUTINE SEQ
  REAL(wp) :: v, make_normalized

  IF (ABS(v) <= 2.225073858507201e-308_wp) THEN
    make_normalized = 0.0_wp
  ELSE
    make_normalized = v
  END IF
END FUNCTION
#else
FUNCTION make_normalized(v)
  REAL(wp) :: v, make_normalized
    make_normalized = v
END FUNCTION
#endif

SUBROUTINE graupel                 ( &
  nvec,ke,                           & !> array dimensions
  ivstart,ivend, kstart,             & !! start/end indicies
  zdt, dz,                           & !! numerics parameters
  t,p,rho,qv,qc,qi,qr,qs,qg,qnc,     & !! prognostic variables
  prr_gsp,pri_gsp,prs_gsp,prg_gsp,   & !! surface precipitation rates
  qrsflux                          )  !  total precipitation flux


  INTEGER, INTENT(IN) ::  &
    nvec      ,    & !> number of horizontal points
    ke        ,    & !! number of grid points in vertical direction
    ivstart   ,    & !! start index for horizontal direction
    ivend     ,    & !! end index   for horizontal direction
    kstart           !! start index for the vertical 

  REAL(KIND=wp), INTENT(IN) :: &
    zdt                    !> time step for integration of microphysics     (  s  )

  REAL(KIND=wp), DIMENSION(:,:), INTENT(IN) ::      &   ! (ie,ke)
    dz              ,    & !> layer thickness of full levels                (  m  )
    rho             ,    & !! density of moist air                          (kg/m3)
    p                      !! pressure                                      ( Pa  )

  REAL(KIND=wp), DIMENSION(:,:), INTENT(INOUT) ::   &   ! dim (ie,ke)
    t               ,    & !> temperature                                   (  K  )
    qv              ,    & !! specific water vapor content                  (kg/kg)
    qc              ,    & !! specific cloud water content                  (kg/kg)
    qi              ,    & !! specific cloud ice   content                  (kg/kg)
    qr              ,    & !! specific rain content                         (kg/kg)
    qs              ,    & !! specific snow content                         (kg/kg)
    qg              ,    & !! specific graupel content                      (kg/kg)
    qrsflux                ! total precipitation flux (nudg)

  REAL(KIND=wp), DIMENSION(:), INTENT(INOUT) ::   &   ! dim (ie)
    prr_gsp,             & !> precipitation rate of rain, grid-scale        (kg/(m2*s))
    pri_gsp,             & !> precipitation rate of ice, grid-scale        (kg/(m2*s))
    prs_gsp,             & !! precipitation rate of snow, grid-scale        (kg/(m2*s))
    prg_gsp,             & !! precipitation rate of graupel, grid-scale     (kg/(m2*s))
    qnc                    !! cloud number concentration

  REAL    (KIND=wp   ), PARAMETER ::  &
    qi0=0.0_wp,           & !threshold for cloud ice 
    qc0=0.0_wp,           & !threshold for cloud liquid 
    zcsg=0.5_wp,          & !coefficient for snow-graupel conversion by riming
    zcrim_g=4.43_wp,      & !
    zrimexp_g=0.94878_wp, & !
    zcagg_g = 2.46_wp ,   & !
    zasmel= 2.95E3_wp ,   & ! DIFF*lh_v*RHO/LHEAT
    zexpsedg=0.217_wp,    & ! exponent for graupel sedimentation
    zvz0g = 12.24_wp  ,   & ! coefficient of sedimentation velocity for graupel
    ztcrit=3339.5_wp        ! factor in calculation of critical temperature
  
  INTEGER (KIND=i4)    ::  &
    iv, k               !> loop indices

  REAL    (KIND=wp   ) ::  &
    nnr,               & ! 
    fpvsw,             & ! name of statement function
    fxna_cooper ,      & ! statement function for ice crystal number, Cooper(1986) 
    ztx  ,             & ! dummy argument for statement functions
    znimax,            & ! maximum number of cloud ice crystals
    znimix,            & ! number of ice crystals at ztmix -> threshold temp for mixed-phase clouds 
    zpvsw0,            & ! sat.vap. pressure at melting temperature
    zqvsw0,            & ! sat.specific humidity at melting temperature
    zqvsw0diff,        & ! qv-zqvsw0  
    zdtr ,             & ! reciprocal of timestep for integration
    zscsum, zscmax, zcorr,  & ! terms for limiting  total cloud water depletion
    zsrsum,            & ! terms for limiting  total rain water depletion
    zsssum,            & ! terms for limiting snow depletion
    znin,              & ! number of cloud ice crystals at nucleation
    fnuc,              & !FR: coefficient needed for Forbes (2012) SLW layer parameterization 
    znid,              & ! number of cloud ice crystals for deposition
    zmi ,              & ! mass of a cloud ice crystal
    zsvidep, zsvisub,  & ! deposition, sublimation of cloud ice
    zsimax , zsisum , zsvmax,& ! terms for limiting total cloud ice depletion
    zqvsw,             & ! sat. specitic humidity at ice and water saturation
    zqvsidiff,         & ! qv-zqvsi
    ztfrzdiff,         & ! ztrfrz-t  
    zztau, zxfac, zx1,  ztt, qliq, qice, qtot, cv,  &   ! some help variables
    ztau, zphi, zhi, zdvtp, ztc, zeff, zlog_10

  REAL    (KIND=wp   ) ::  &
    zqct   ,& ! layer tendency of cloud water
    zqvt   ,& ! layer tendency of water vapour
    zqit   ,& ! layer tendency of cloud ice
    zqrt   ,& ! layer tendency of rain
    zqst   ,& ! layer tendency of snow
    zqgt      ! layer tendency of graupel

  REAL    (KIND=wp   ) ::  &
    zlnqrk,zlnqsk,zlnqik,     & !
    zlnlogmi,zlnqgk,ccswxp_ln1o2,zvzxp_ln1o2,zbvi_ln1o2,zexpsedg_ln1o2, &
    qcg,tg,qvg,qrg,qsg,qgg,qig,rhog,ppg,alf,bet,m2s,m3s,hlp,            &
    qcgk_1,maxevap,temp_c

  LOGICAL :: &
    llqs,llqc,llqi,llqg,llqr  !   switch for existence of qr, qs, qc, qi

  REAL(KIND=wp), DIMENSION(nvec,ke) ::   &
    t_in               ,    & !> temperature                                   (  K  )
    qv_in              ,    & !! specific water vapor content                  (kg/kg)
    qc_in              ,    & !! specific cloud water content                  (kg/kg)
    qi_in              ,    & !! specific cloud ice   content                  (kg/kg)
    qr_in              ,    & !! specific rain content                         (kg/kg)
    qs_in              ,    & !! specific snow content                         (kg/kg)
    qg_in                     !! specific graupel content                      (kg/kg)

  REAL    (KIND=wp   ) ::  &
    zqvsi             ,     & !> sat. specitic humidity at ice and water saturation
    zvzr        (nvec),     & !
    zvzs        (nvec),     & !
    zvzg        (nvec),     & ! 
    zvzi        (nvec),     & ! terminal fall velocity of ice
    zpkr        (nvec),     & !
    zpks        (nvec),     & !
    zpkg        (nvec),     & ! 
    zpki        (nvec),     & ! precipitation flux of ice
    zprvr       (nvec),     & !
    zprvs       (nvec),     & !
    zprvg       (nvec),     & !
    zprvi       (nvec),     & !
    zqvsw_up    (nvec),     & ! sat. specitic humidity at ice and water saturation
    zcsdep            ,     & !
    zcidep            

 REAL (KIND = wp )  ::  zlhv(nvec), zlhs(nvec)
    
 REAL    (KIND=wp   ) ::  &    
    zsrmax            ,     & !
    zssmax            ,     & !
    zsgmax            ,     & !
    zvz0s             ,     & !
    zcrim             ,     & !
    zcagg             ,     & !
    zbsdep            ,     & !
    zcslam            ,     & !
    zn0s              ,     & !
    zimr              ,     & !
    zims              ,     & !
    zimg              ,     & !
    zimi              ,     & !
    zzar              ,     & !
    zzas              ,     & !
    zzag              ,     & !
    zzai              ,     & !
    zqrk              ,     & !
    zqsk              ,     & !
    zqgk              ,     & !
    zqik              ,     & !
    zdtdh             ,     & !
    z1orhog           ,     & ! 1/rhog
    zrho1o2           ,     & ! (rho0/rhog)**1/2
    zrhofac_qi        ,     & ! (rho0/rhog)**icesedi_exp
    zeln7o8qrk        ,     & !
    zeln7o4qrk        ,     & ! FR new  
    zeln27o16qrk      ,     & !
    zeln13o8qrk       ,     & !
    zeln3o4qsk        ,     & ! 
    zeln6qgk          ,     & !
    zeln8qsk          ,     & !
    zelnrimexp_g      


  REAL    (KIND=wp   ) ::  &
    scau   , & ! transfer rate due to autoconversion of cloud water
    scac   , & ! transfer rate due to accretion of cloud water
    snuc   , & ! transfer rate due nucleation of cloud ice
    scfrz  , & ! transfer rate due homogeneous freezing of cloud water
    simelt , & ! transfer rate due melting of cloud ice
    sidep  , & ! transfer rate due depositional growth of cloud ice
    ssdep  , & ! transfer rate due depositional growth of snow
    sgdep  , & ! transfer rate due depositional growth of graupel
    sdau   , & ! transfer rate due depositional cloud ice autoconversion
    srim   , & ! transfer rate due riming of snow
    srim2  , & ! transfer rate due riming of graupel
    sconsg , & ! transfer rate due to conversion from snow to graupel by riming  
    sshed  , & ! transfer rate due shedding
    sicri  , & ! transfer rate due cloud ice collection by rain (sink qi)
    srcri  , & ! transfer rate due cloud ice collection by rain (sink qr)
    sagg   , & ! transfer rate due aggregation of snow and cloud ice
    sagg2  , & ! transfer rate due aggregation of snow and cloud ice
    siau   , & ! transfer rate due autoconversion of cloud ice
    ssmelt , & ! transfer rate due melting of snow
    sgmelt , & ! transfer rate due melting of graupel
    sev    , & ! transfer rate due evaporation of rain
    sconr  , & ! transfer rate due to condensation on melting snow/graupel
    srfrz  , & ! transfer rate due to rainwater freezing
    reduce_dep,&!FR: coefficient: reduce deposition at cloud top (Forbes 2012)
    dist_cldtop(nvec) !FR: distance from cloud top layer 


! saturation vapour pressure over water (fpvsw), over ice (fpvsi)
! and specific humidity at vapour saturation (fqvs)
  fpvsw(ztx)     = sat_pres_water(ztx)

! Number of activate ice crystals;  ztx is temperature
  fxna_cooper(ztx) = 5.0E+0_wp * EXP(0.304_wp * (tmelt - ztx))   ! FR: Cooper (1986) used by Greg Thompson(2008)

  ! Input data
  !$ACC DATA                                                     &
  !$ACC PRESENT( dz, t, p, rho, qv, qc, qi, qr, qs, qg, qnc )    &
  !$ACC PRESENT( prr_gsp, pri_gsp, prs_gsp, prg_gsp, qrsflux )   &
  ! automatic arrays
  !$ACC CREATE( zvzr, zvzs, zvzg, zvzi )                         &
  !$ACC CREATE( zpkr, zpks, zpkg, zpki )                         &
  !$ACC CREATE( zprvr, zprvs, zprvi, zqvsw_up, zprvg )           &
  !$ACC CREATE( dist_cldtop, zlhv, zlhs )

  znimax = znimax_Thom
  znimix = fxna_cooper(ztmix) ! number of ice crystals at temp threshold for mixed-phase clouds

  zpvsw0 = fpvsw(tmelt)  ! sat. vap. pressure for t = tmelt
  zlog_10 = LOG(10._wp) ! logarithm of 10
  
  ! Precomputations for optimization
  ccswxp_ln1o2   = EXP (ccswxp * LOG (0.5_wp))
  zvzxp_ln1o2    = EXP (zvzxp * LOG (0.5_wp))
  zbvi_ln1o2     = EXP (zbvi * LOG (0.5_wp))
  zexpsedg_ln1o2 = EXP (zexpsedg * LOG (0.5_wp))

! timestep for calculations
  zdtr  = 1.0_wp / zdt

  ! Delete precipitation fluxes from previous timestep
  !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
  !$ACC LOOP GANG VECTOR
  DO iv = ivstart, ivend
    prr_gsp (iv) = 0.0_wp
    pri_gsp (iv) = 0.0_wp
    prs_gsp (iv) = 0.0_wp
    prg_gsp (iv) = 0.0_wp
    zpkr(iv)     = 0.0_wp
    zpks(iv)     = 0.0_wp
    zpkg(iv)     = 0.0_wp
    zpki(iv)     = 0.0_wp
    zprvr(iv)    = 0.0_wp
    zprvs(iv)    = 0.0_wp
    zprvg(iv)    = 0.0_wp
    zprvi(iv)    = 0.0_wp
    zvzr(iv)     = 0.0_wp
    zvzs(iv)     = 0.0_wp
    zvzg(iv)     = 0.0_wp
    zvzi(iv)     = 0.0_wp
    dist_cldtop(iv) = 0.0_wp
    zqvsw_up(iv) = 0.0_wp
  END DO 
  !$ACC END PARALLEL

! *********************************************************************
! Loop from the top of the model domain to the surface to calculate the
! transfer rates  and sedimentation terms
! *********************************************************************

  !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
  !$ACC LOOP SEQ
  DO  k = kstart, ke  ! loop over levels

    !$ACC LOOP GANG(STATIC:1) VECTOR PRIVATE( alf, bet, fnuc, hlp, llqc, llqg, llqi, llqr, &
    !$ACC                           llqs, m2s, m3s, maxevap, nnr, ppg, qcg,      &
    !$ACC                           qcgk_1, qgg, qig, qrg, qsg, qvg, reduce_dep, &
    !$ACC                           rhog, sagg, sagg2, scac, scau, scfrz, sconr, &
    !$ACC                           sconsg, sdau, sev, sgdep, sgmelt, siau,      &
    !$ACC                           sicri, sidep, simelt, snuc, srcri, srfrz,    &
    !$ACC                           srim, srim2, ssdep, sshed, ssmelt, temp_c,   &
    !$ACC                           tg, z1orhog, zbsdep, zcagg, zcidep, zcorr,   &
    !$ACC                           zcrim, zcsdep, zcslam, zdtdh, zdvtp, zeff,   &
    !$ACC                           zeln13o8qrk, zeln27o16qrk, zeln3o4qsk,       &
    !$ACC                           zeln6qgk, zeln7o4qrk, zeln7o8qrk, zeln8qsk,  &
    !$ACC                           zelnrimexp_g, zhi, zimg, zimi, zimr, zims,   &
    !$ACC                           zlnlogmi, zlnqgk, zlnqik, zlnqrk, zlnqsk,    &
    !$ACC                           zmi, zn0s, znid, znin, zphi, zqct, zqgk,     &
    !$ACC                           zqgt, zqik, zqit, zqrk, zqrt, zqsk, zqst,    &
    !$ACC                           zqvsi, zqvsidiff, zqvsw, zqvsw0, zqvsw0diff, &
    !$ACC                           zqvt, zrho1o2, zrhofac_qi, zscmax, zscsum,   &
    !$ACC                           zsgmax, zsimax, zsisum, zsrmax, zsrsum,      &
    !$ACC                           zssmax, zsssum, zsvidep, zsvisub, zsvmax,    &
    !$ACC                           ztau, ztc, ztfrzdiff, ztt, zvz0s, zx1,       &
    !$ACC                           qliq, qice, qtot, cv,                        &
    !$ACC                           zxfac, zzag, zzai, zzar, zzas, zztau )
    DO iv = ivstart, ivend  !loop over horizontal domain

      !----------------------------------------------------------------------------
      ! Section 2: Check for existence of rain and snow
      !            Initialize microphysics and sedimentation scheme
      !----------------------------------------------------------------------------

      zcrim  = 0.0_wp
      zcagg  = 0.0_wp
      zbsdep = 0.0_wp
      zvz0s  = 0.0_wp
      zn0s   = zn0s0
      reduce_dep = 1.0_wp  !FR: Reduction coeff. for dep. growth of rain and ice  

      !----------------------------------------------------------------------------
      ! 2.1: Preparations for computations and to check the different conditions
      !----------------------------------------------------------------------------

      qrsflux(iv,k) = 0.0_wp

      qrg  = make_normalized(qr(iv,k))
      qsg  = make_normalized(qs(iv,k))
      qgg  = make_normalized(qg(iv,k))
      qvg  = make_normalized(qv(iv,k))
      qcg  = make_normalized(qc(iv,k))
      qig  = make_normalized(qi(iv,k))

      tg   = t(iv,k) 
      ppg  = p(iv,k)
      rhog = rho(iv,k)

      !..for density correction of fall speeds
      z1orhog = 1.0_wp/rhog
      hlp     = LOG(zrho0*z1orhog)
      zrho1o2 = EXP(hlp*x1o2)
      zrhofac_qi = EXP(hlp*icesedi_exp)

      zqrk = qrg * rhog
      zqsk = qsg * rhog
      zqgk = qgg * rhog
      zqik = qig * rhog

      llqr = zqrk > zqmin
      llqs = zqsk > zqmin
      llqg = zqgk > zqmin
      llqi = zqik > zqmin

      zdtdh = 0.5_wp * zdt / dz(iv,k)

      zzar = zqrk/zdtdh + zprvr(iv) + zpkr(iv)
      zzas = zqsk/zdtdh + zprvs(iv) + zpks(iv)
      zzag = zqgk/zdtdh + zprvg(iv) + zpkg(iv)
      zzai = zqik/zdtdh + zprvi(iv) + zpki(iv)

      zpkr(iv) = 0.0_wp
      zpks(iv) = 0.0_wp
      zpkg(iv) = 0.0_wp
      zpki(iv) = 0.0_wp

      !-------------------------------------------------------------------------
      ! qs_prepare:
      !-------------------------------------------------------------------------
      IF (llqs) THEN
        ! Calculate n0s using the temperature-dependent moment
        ! relations of Field et al. (2005)
        ztc = tg - tmelt
        ztc = MAX(MIN(ztc,0.0_wp),-40.0_wp)

        nnr  = 3._wp
        hlp = mma(1) + mma(2)*ztc + mma(3)*nnr + mma(4)*ztc*nnr &
            + mma(5)*ztc**2 + mma(6)*nnr**2 + mma(7)*ztc**2*nnr &
            + mma(8)*ztc*nnr**2 + mma(9)*ztc**3 + mma(10)*nnr**3
        alf = EXP(hlp*zlog_10) ! 10.0_wp**hlp
        bet = mmb(1) + mmb(2)*ztc + mmb(3)*nnr + mmb(4)*ztc*nnr &
            + mmb(5)*ztc**2 + mmb(6)*nnr**2 + mmb(7)*ztc**2*nnr &
            + mmb(8)*ztc*nnr**2 + mmb(9)*ztc**3 + mmb(10)*nnr**3

        ! Here is the exponent bms=2.0 hardwired! not ideal! (Uli Blahak)
        m2s = qsg * rhog / zams   ! UB rho added as bugfix
        m3s = alf*EXP(bet*LOG(m2s))

        hlp  = zn0s1*EXP(zn0s2*ztc)
        zn0s = 13.50_wp * m2s * (m2s / m3s)**3
        zn0s = MAX(zn0s,0.5_wp*hlp)
        zn0s = MIN(zn0s,1.0E2_wp*hlp)
        zn0s = MIN(zn0s,1.0E9_wp)
        zn0s = MAX(zn0s,1.0E6_wp)

        zcrim  = ccsrim*zn0s
        zcagg  = ccsagg*zn0s
        zbsdep = ccsdep*SQRT(v0snow)
        zvz0s  = ccsvel*EXP(ccsvxp * LOG(zn0s))
        IF (lrain) THEN
          zlnqsk = zvz0s * EXP (ccswxp * LOG (zqsk)) * zrho1o2
          ! Prevent terminal fall speed of snow from being zero at the surface level
          IF ( k == ke ) zlnqsk = MAX( zlnqsk, v_sedi_snow_min )
          zpks(iv) = zqsk * zlnqsk
          IF (zvzs(iv) == 0.0_wp) THEN
            zvzs(iv) = zlnqsk * ccswxp_ln1o2
          ENDIF
        ENDIF 
      ENDIF ! qs_prepare
    
      ! sedimentation fluxes

      !-------------------------------------------------------------------------
      ! qr_sedi:
      !-------------------------------------------------------------------------

      IF (llqr .AND. lrain) THEN
        zlnqrk = zvz0r * EXP (zvzxp * LOG (zqrk)) * zrho1o2
        ! Prevent terminal fall speed of rain from being zero at the surface level
        IF ( k == ke ) zlnqrk = MAX( zlnqrk, v_sedi_rain_min )
        zpkr(iv) = zqrk * zlnqrk
        IF (zvzr(iv) == 0.0_wp) THEN
          zvzr(iv) = zlnqrk * zvzxp_ln1o2
        ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      ! qg_sedi:
      !-------------------------------------------------------------------------

      IF (llqg .AND. lrain) THEN
        zlnqgk = zvz0g * EXP (zexpsedg * LOG (zqgk)) * zrho1o2
        ! Prevent terminal fall speed of graupel from being zero at the surface level
        IF ( k == ke ) zlnqgk = MAX( zlnqgk, v_sedi_graupel_min )
        zpkg(iv) = zqgk * zlnqgk
        IF (zvzg(iv) == 0.0_wp) THEN
          zvzg(iv) = zlnqgk * zexpsedg_ln1o2
        ENDIF
      ENDIF ! qg_sedi

      !-------------------------------------------------------------------------
      ! qi_sedi:
      !-------------------------------------------------------------------------

      IF (llqi .AND. lrain) THEN
        zlnqik = zvz0i * EXP (zbvi * LOG (zqik)) * zrhofac_qi
        zpki(iv) = zqik * zlnqik
        IF (zvzi(iv) == 0.0_wp) THEN
          zvzi(iv) = zlnqik * zbvi_ln1o2
        ENDIF
      ENDIF  ! qi_sedi

      ! Prevent terminal fall speeds of precip hydrometeors from being zero at the surface level
      IF ( k == ke .AND. lrain ) THEN
        zvzr(iv) = MAX( zvzr(iv), v_sedi_rain_min )
        zvzs(iv) = MAX( zvzs(iv), v_sedi_snow_min )
        zvzg(iv) = MAX( zvzg(iv), v_sedi_graupel_min )
      ENDIF

      !--------------------------------------------------------------------------
      ! 2.3: Second part of preparations
      !--------------------------------------------------------------------------

      zeln7o8qrk    = 0.0_wp
      zeln7o4qrk    = 0.0_wp
      zeln27o16qrk  = 0.0_wp
      zeln13o8qrk   = 0.0_wp
      zeln3o4qsk    = 0.0_wp
      zeln8qsk      = 0.0_wp
      zeln6qgk      = 0.0_wp
      zelnrimexp_g  = 0.0_wp
      zsrmax        = 0.0_wp
      zssmax        = 0.0_wp
      zsgmax        = 0.0_wp

      zcsdep        = 3.367E-2_wp
      zcidep        = 1.3E-5_wp
      zcslam        = 1e10_wp

      scau          = 0.0_wp
      scac          = 0.0_wp
      snuc          = 0.0_wp
      scfrz         = 0.0_wp
      simelt        = 0.0_wp
      sidep         = 0.0_wp
      ssdep         = 0.0_wp
      sgdep         = 0.0_wp
      sdau          = 0.0_wp
      srim          = 0.0_wp
      srim2         = 0.0_wp
      sshed         = 0.0_wp
      sicri         = 0.0_wp
      srcri         = 0.0_wp
      sagg          = 0.0_wp
      sagg2         = 0.0_wp
      siau          = 0.0_wp
      ssmelt        = 0.0_wp
      sgmelt        = 0.0_wp
      sev           = 0.0_wp
      sconr         = 0.0_wp
      sconsg        = 0.0_wp
      srfrz         = 0.0_wp

      zpkr(iv)   = MIN( zpkr(iv) , zzar )
      zpks(iv)   = MIN( zpks(iv) , zzas )
      zpkg(iv)   = MIN( zpkg(iv) , MAX(0._wp,zzag) )
      zpki(iv)   = MIN( zpki(iv) , zzai )

      zzar   = zdtdh * (zzar-zpkr(iv))
      zzas   = zdtdh * (zzas-zpks(iv))
      zzag   = zdtdh * (zzag-zpkg(iv))
      zzai   = zdtdh * (zzai-zpki(iv))

      zimr   = 1.0_wp / (1.0_wp + zvzr(iv) * zdtdh)
      zims   = 1.0_wp / (1.0_wp + zvzs(iv) * zdtdh)
      zimg   = 1.0_wp / (1.0_wp + zvzg(iv) * zdtdh)
      zimi   = 1.0_wp / (1.0_wp + zvzi(iv) * zdtdh)

      zqrk   = zzar*zimr
      zqsk   = zzas*zims
      zqgk   = zzag*zimg
      zqik   = zzai*zimi

      zqvsi = qsat_ice_rho(tg,rhog)

      llqr = zqrk > zqmin
      llqs = zqsk > zqmin
      llqg = zqgk > zqmin
      llqc =  qcg > zqmin
      llqi =  qig > zqmin

      !!----------------------------------------------------------------------------
      !! 2.4: IF (llqr): ic1
      !!----------------------------------------------------------------------------

      IF (llqr) THEN
        zlnqrk   = LOG (zqrk)
        zsrmax   = zzar/rhog*zdtr  ! GZ: shifting this computation ahead of the IF condition changes results!
        IF ( qig+qcg > zqmin ) THEN
          zeln7o8qrk   = EXP (x7o8   * zlnqrk)
        ENDIF
        IF ( tg < ztrfrz ) THEN
          zeln7o4qrk   = EXP (x7o4   * zlnqrk) !FR new
          zeln27o16qrk = EXP (x27o16 * zlnqrk)
        ENDIF
        IF (llqi) THEN
          zeln13o8qrk  = EXP (x13o8  * zlnqrk)
        ENDIF
      ENDIF

      !!----------------------------------------------------------------------------
      !! 2.5: IF (llqs): ic2
      !!----------------------------------------------------------------------------

! ** GZ: the following computation differs substantially from the corresponding code in cloudice **
      IF (llqs) THEN
        zlnqsk   = LOG (zqsk)
        zssmax   = zzas / rhog*zdtr  ! GZ: shifting this computation ahead of the IF condition changes results!
        IF (qig+qcg > zqmin) THEN
          zeln3o4qsk = EXP (x3o4 *zlnqsk)
        ENDIF
        zeln8qsk = EXP (0.8_wp *zlnqsk)
      ENDIF

      !!----------------------------------------------------------------------------
      !! 2.6: IF (llqg): ic3
      !!----------------------------------------------------------------------------

      IF (zqgk > zqmin) THEN
        zlnqgk   = LOG (zqgk)
        zsgmax   = zzag / rhog*zdtr
        IF (qig+qcg > zqmin) THEN
          zelnrimexp_g = EXP (zrimexp_g * zlnqgk)
        ENDIF
        zeln6qgk = EXP (0.6_wp *zlnqgk)
      ENDIF

      !!----------------------------------------------------------------------------
      !! 2.7:  slope of snow PSD and coefficients for depositional growth (llqi,llqs)
      !!----------------------------------------------------------------------------    

      IF ((qig > zqmin) .OR. (zqsk > zqmin)) THEN
        zdvtp  = ccdvtp * EXP(1.94_wp * LOG(tg)) / ppg
        zhi    = ccshi1*zdvtp*rhog*zqvsi/(tg*tg)
        hlp    = zdvtp / (1.0_wp + zhi)
        zcidep = ccidep * hlp

        IF (llqs) THEN
          zcslam = EXP(ccslxp * LOG(ccslam * zn0s / zqsk ))
          zcslam = MIN(zcslam,1.0E15_wp)
          zcsdep = 4.0_wp * zn0s * hlp
        ENDIF
      ENDIF

      !!----------------------------------------------------------------------------
      !! 2.8: Deposition nucleation for low temperatures below a threshold (llqv)
      !!----------------------------------------------------------------------------    

      IF (( tg < zthet .AND. qvg >  8.E-6_wp .AND. qig <= 0.0_wp )) THEN
        IF( qvg > zqvsi ) THEN
          znin  = MIN( fxna_cooper(tg), znimax )
          snuc = zmi0 * z1orhog * znin * zdtr
        ENDIF
      ENDIF

      !!--------------------------------------------------------------------------
      !! Section 3: Search for cloudy grid points with cloud water and
      !!            calculation of the conversion rates involving qc (ic6)
      !!--------------------------------------------------------------------------

      IF (qcg > zqmin .AND. lautocnv) THEN
        llqs = zqsk > zqmin

        zscmax = qcg*zdtr
        IF( tg > zthn ) THEN
          IF (iautocon == 0) THEN
            ! Kessler (1969) autoconversion rate
            scau = zccau * MAX( qcg - qc0, 0.0_wp )
            scac = zcac  * qcg * zeln7o8qrk
          ELSEIF (iautocon == 1) THEN
            ! Seifert and Beheng (2001) autoconversion rate
            ! with constant cloud droplet number concentration qnc
            IF (qcg > 1.0E-6_wp) THEN
              ztau = MIN(1.0_wp-qcg/(qcg+qrg),0.9_wp)
              ztau = MAX(ztau,1.E-30_wp)
              hlp  = EXP(zkphi2*LOG(ztau))
              zphi = zkphi1 * hlp * (1.0_wp - hlp)**3
              scau = zconst * qcg*qcg*qcg*qcg/(qnc(iv)*qnc(iv)) &
                   * (1.0_wp + zphi/(1.0_wp - ztau)**2)
              zphi = (ztau/(ztau+zkphi3))**4
              scac = zkcac * qcg * qrg * zphi
            ELSE
              scau = 0.0_wp
              scac = 0.0_wp
            ENDIF
          ENDIF
          IF (llqr) THEN
            ! Calculation of in-cloud rainwater freezing
            IF ( tg < ztrfrz .AND. qrg > 0.1_wp*qcg ) THEN
              srfrz = zcrfrz1*(EXP(zcrfrz2*(ztrfrz-tg))-1.0_wp ) * zeln7o4qrk
            ENDIF
          ENDIF
          IF (llqs) THEN
            srim = zcrim * qcg *  EXP(ccsaxp * LOG(zcslam))
          ENDIF
          srim2 = zcrim_g * qcg * zelnrimexp_g
          IF( tg >= tmelt ) THEN
            sshed = srim+srim2
            srim  = 0.0_wp
            srim2 = 0.0_wp
          ELSE
            IF (qcg >= qc0) THEN
              sconsg = zcsg * qcg * zeln3o4qsk
            ENDIF
          ENDIF
          ! Check for maximum depletion of cloud water and adjust the
          ! transfer rates accordingly
          zscsum = scau + scac + srim + srim2 + sshed
          zcorr  = zscmax / MAX( zscmax, zscsum )
          scau   = zcorr*scau
          scac   = zcorr*scac
          srim   = zcorr*srim
          srim2  = zcorr*srim2
          sshed  = zcorr*sshed
          sconsg = MIN (sconsg, srim+zssmax)
        ELSE !tg >= tg: ! hom. freezing of cloud and rain water
          scfrz = zscmax
          srfrz = zsrmax
        ENDIF
        ! Calculation of heterogeneous nucleation of cloud ice.
        ! This is done in this section, because we require water saturation
        ! for this process (i.e. the existence of cloud water) to exist.
        ! Heterogeneous nucleation is assumed to occur only when no
        ! cloud ice is present and the temperature is below a nucleation
        ! threshold.

        IF( tg <= 267.15_wp .AND. .NOT.llqi ) THEN   
          znin  = MIN( fxna_cooper(tg), znimax )
          snuc = zmi0 * z1orhog * znin * zdtr
        ENDIF
        ! Calculation of reduction of depositional growth at cloud top (Forbes 2012)
        IF( k>kstart .AND. k<ke ) THEN
          znin = MIN(fxna_cooper(tg), znimax )
          fnuc = MIN(znin/znimix, 1.0_wp)

          qcgk_1 = qi(iv,k-1) + qs(iv,k-1) + qg(iv,k-1)

          !! distance from cloud top
          IF( qv(iv,k-1) + qc(iv,k-1) < zqvsw_up(iv) .AND. qcgk_1 .LT. zqmin ) THEN      ! upper cloud layer
            dist_cldtop(iv) = 0.0_wp    ! reset distance to upper cloud layer
          ELSE
            dist_cldtop(iv) = dist_cldtop(iv) + dz(iv,k)
          END IF

          ! with asymptotic behaviour dz -> 0 (xxx)
          !        reduce_dep = MIN(fnuc + (1.0_wp-fnuc)*(reduce_dep_ref + &
          !                             dist_cldtop(iv)/dist_cldtop_ref + &
          !                             (1.0_wp-reduce_dep_ref)*(zdh/dist_cldtop_ref)**4), 1.0_wp)

          ! without asymptotic behaviour dz -> 0
          reduce_dep = MIN(fnuc + (1.0_wp-fnuc)*(reduce_dep_ref + &
                        dist_cldtop(iv)/dist_cldtop_ref), 1.0_wp)

        END IF ! Reduction of dep. growth of snow/ice 

      ENDIF

      !------------------------------------------------------------------------
      ! Section 4: Search for cold grid points with cloud ice and/or snow and
      !            calculation of the conversion rates involving qi, qs and qg
      !------------------------------------------------------------------------

      IF ( (qig > zqmin) .OR. (zqsk > zqmin) .OR. zqgk > zqmin ) THEN
        llqs =  zqsk > zqmin
        llqi =   qig > zqmin

        IF (tg<=tmelt) THEN           ! cold case 

          zqvsidiff = qvg-zqvsi
          zsvmax    = zqvsidiff * zdtr
          IF (llqi) THEN

            znin   = MIN( fxna_cooper(tg), znimax )
            ! Change in sticking efficiency needed in case of cloud ice sedimentation  (GZ with addition of ice sed)
            zeff     = MAX(MIN(EXP(0.09_wp*(tg-tmelt)),1.0_wp), zceff_min, zceff_fac*(tg-tmin_iceautoconv)) 
            sagg      = zeff * qig * zcagg * EXP(ccsaxp*LOG(zcslam))
            sagg2     = zeff * qig * zcagg_g * zelnrimexp_g
            siau      = zeff * zciau * MAX( qig - qi0, 0.0_wp )
            zmi       = MIN( rhog*qig/znin, zmimax )
            zmi       = MAX( zmi0, zmi )
            znid      = rhog * qig/zmi
            zlnlogmi  = LOG (zmi)
            sidep     = zcidep * znid * EXP(0.33_wp * zlnlogmi) * zqvsidiff
            zsvidep   = 0.0_wp
            zsvisub   = 0.0_wp
            ! for sedimenting quantities the maximum 
            ! allowed depletion is determined by the predictor value. 
            zsimax  = zzai*z1orhog*zdtr
            IF( sidep > 0.0_wp ) THEN
              sidep = sidep * reduce_dep  !FR new: depositional growth reduction
              zsvidep = MIN( sidep, zsvmax )
            ELSEIF ( sidep < 0.0_wp ) THEN
              zsvisub  =   MAX (   sidep,  zsvmax)
              zsvisub  = - MAX ( zsvisub, -zsimax)
            ENDIF
            zlnlogmi   = LOG  (zmsmin/zmi)
            zztau      = 1.5_wp*( EXP(0.66_wp*zlnlogmi) - 1.0_wp)
            sdau       = zsvidep/zztau
            sicri      = zcicri * qig * zeln7o8qrk
            IF (qsg > 1.e-7_wp) srcri = zcrcri * (qig/zmi) * zeln13o8qrk
          ELSE
            zsimax    =  0.0_wp
            zsvidep   =  0.0_wp
            zsvisub   =  0.0_wp
          ENDIF

          zxfac = 1.0_wp + zbsdep * EXP(ccsdxp*LOG(zcslam))
          ssdep = zcsdep * zxfac * zqvsidiff / (zcslam+zeps)**2
          ! FR new: depositional growth reduction
          ! GZ: This limitation, which was missing in the original graupel scheme,
          ! is crucial for numerical stability in the tropics!
          IF (ssdep > 0.0_wp) ssdep = MIN(ssdep*reduce_dep, zsvmax-zsvidep)
          ! Suppress depositional growth of snow if the existing amount is too small for a
          ! a meaningful distiction between cloud ice and snow
          IF (qsg <= 1.e-7_wp) ssdep = MIN(ssdep, 0.0_wp)
          ! ** GZ: this numerical fit should be replaced with a physically more meaningful formulation **
          sgdep = (0.398561_wp-0.00152398_wp*tg                 &
                   + 2554.99_wp/ppg+ 2.6531E-7_wp*ppg) *        &
                 zqvsidiff * zeln6qgk
          ! Check for maximal depletion of cloud ice
          ! No check is done for depositional autoconversion (sdau) because
          ! this is a always a fraction of the gain rate due to
          ! deposition (i.e the sum of this rates is always positive)
          zsisum = siau + sagg + sagg2 + sicri + zsvisub
          zcorr  = 0.0_wp
          IF( zsimax > 0.0_wp ) zcorr  = zsimax / MAX( zsimax, zsisum )
          sidep  = zsvidep - zcorr*zsvisub
          siau   = zcorr*siau
          sagg   = zcorr*sagg
          sagg2  = zcorr*sagg2
          sicri  = zcorr*sicri
          IF ( zqvsidiff < 0.0_wp ) THEN
            ssdep = MAX(ssdep, - zssmax)
            sgdep = MAX(sgdep, - zsgmax)
          ENDIF

        ELSE ! tg > 0 - warm case

          !------------------------------------------------------------------------
          ! Section 5: Search for warm grid points with cloud ice and/or snow and
          !            calculation of the melting rates of qi and ps
          !------------------------------------------------------------------------

          ! cloud ice melts instantaneously
          simelt = zzai*z1orhog*zdtr

          zqvsw0     = zpvsw0/(rhog * rv *tmelt)
          zqvsw0diff = qvg-zqvsw0
          ! ** GZ: several numerical fits in this section should be replaced with physically more meaningful formulations **
          IF ( tg > (tmelt-ztcrit*zqvsw0diff) ) THEN
            !calculate melting rate
            zx1         = (tg - tmelt) + zasmel*zqvsw0diff
            ssmelt = (79.6863_wp/ppg+0.612654E-3_wp)* zx1 * zeln8qsk
            ssmelt = MIN (ssmelt,zssmax)
            sgmelt = (12.31698_wp/ppg+7.39441e-05_wp)* zx1 * zeln6qgk
            sgmelt = MIN (sgmelt, zsgmax)
            !deposition + melting, ice particle temperature: tmelt
            !calculation without howell-factor!
            ssdep  = (31282.3_wp/ppg+0.241897_wp)       &
                    * zqvsw0diff * zeln8qsk
            sgdep  = (0.153907_wp-ppg*7.86703e-07_wp)  &
                    * zqvsw0diff * zeln6qgk
            IF (zqvsw0diff < 0.0_wp) THEN
              !melting + evaporation of snow/graupel
              ssdep = MAX (-zssmax,ssdep)
              sgdep = MAX (-zsgmax,sgdep)
              !melt water evaporates
              ssmelt = ssmelt+ssdep
              sgmelt = sgmelt+sgdep
              ssmelt = MAX( ssmelt, 0.0_wp )
              sgmelt = MAX( sgmelt, 0.0_wp )
            ELSE
              !deposition on snow/graupel is interpreted as increase
              !in rain water ( qv --> qr, sconr)
              !therefore,  sconr=(zssdep+zsgdep)
              sconr=ssdep+sgdep
              ssdep=0.0_wp
              sgdep=0.0_wp
            ENDIF
          ELSE
            !if t<t_crit
            !no melting, only evaporation of snow/graupel
            zqvsw      = qsat_rho(tg,rhog)
            zqvsw_up(iv) = zqvsw
            zqvsidiff  = qvg-zqvsw
            ssdep = (0.28003_wp-ppg*0.146293E-6_wp) &
                     * zqvsidiff * zeln8qsk
            sgdep = (0.0418521_wp-ppg*4.7524E-8_wp) &
                     * zqvsidiff *zeln6qgk
            ssdep = MAX(-zssmax ,ssdep )
            sgdep = MAX(-zsgmax ,sgdep )
          ENDIF !t_crit
        ENDIF !tg
      ENDIF

      !--------------------------------------------------------------------------
      ! Section 6: Search for grid points with rain in subsaturated areas
      !            and calculation of the evaporation rate of rain
      !--------------------------------------------------------------------------

      zqvsw        = qsat_rho(tg,rhog)
      zqvsw_up(iv) = zqvsw

      IF( (llqr) .AND. (qvg+qcg <= zqvsw)) THEN

        zlnqrk   = LOG (zqrk)
        zx1      = 1.0_wp + zbev * EXP (zbevxp  * zlnqrk)
        !sev  = zcev*zx1*(zqvsw - qvg) * EXP (zcevxp  * zlnqrk)
        ! Limit evaporation rate in order to avoid overshoots towards supersaturation
        ! the pre-factor approximates (esat(T_wb)-e)/(esat(T)-e) at temperatures between 0 degC and 30 degC
        temp_c = tg - tmelt
        maxevap     = (0.61_wp-0.0163_wp*temp_c+1.111e-4_wp*temp_c**2)*(zqvsw-qvg)/zdt
        sev    = MIN(zcev*zx1*(zqvsw - qvg) * EXP (zcevxp  * zlnqrk), maxevap)

        IF( tg > zthn ) THEN
          ! Calculation of below-cloud rainwater freezing
          IF ( tg < ztrfrz ) THEN
            !FR new: reduced rain freezing rate
            srfrz = zcrfrz1*(EXP(zcrfrz2*(ztrfrz-tg))-1.0_wp ) * zeln7o4qrk
          ENDIF
        ELSE ! Hom. freezing of rain water
          srfrz = zsrmax
        ENDIF
      ENDIF

      !--------------------------------------------------------------------------
      ! Section 7: Calculate the total tendencies of the prognostic variables.
      !            Update the prognostic variables in the interior domain.
      !--------------------------------------------------------------------------

      zsrsum = sev + srfrz + srcri
      zcorr  = 1.0_wp
      IF(zsrsum > 0._wp) THEN
        zcorr  = zsrmax / MAX( zsrmax, zsrsum )
      ENDIF
      sev   = zcorr*sev
      srfrz = zcorr*srfrz
      srcri = zcorr*srcri

      ! limit snow depletion in order to avoid negative values of qs
      zcorr  = 1.0_wp
      IF (ssdep <= 0._wp) THEN
        zsssum = ssmelt + sconsg - ssdep
        IF(zsssum > 0._wp) THEN
          zcorr  = zssmax / MAX( zssmax, zsssum )
        ENDIF
        ssmelt = zcorr * ssmelt
        sconsg = zcorr * sconsg
        ssdep  = zcorr * ssdep
      ELSE
        zsssum = ssmelt + sconsg
        IF(zsssum > 0._wp) THEN
          zcorr  = zssmax / MAX( zssmax, zsssum )
        ENDIF
        ssmelt = zcorr * ssmelt
        sconsg = zcorr * sconsg
      ENDIF

      IF (lcold) THEN
        zqct =  -scau   - scac   - sshed  + simelt - scfrz  - srim   - srim2 
        zqrt =   scau   + scac   + sshed  - sev    + sconr  + ssmelt + sgmelt - srcri  - srfrz 
        zqit =   snuc   + scfrz  - simelt - sicri  + sidep  - sdau   - sagg   - sagg2  - siau
        zqst =   siau   + sdau   - ssmelt + srim   + ssdep  + sagg   - sconsg
        zqgt =   sagg2  - sgmelt + sicri  + srcri  + sgdep  + srfrz  + srim2  + sconsg
      ELSE 
        zqct = - scau   - scac   - sshed
        zqrt =   scau   + sshed  + scac   - sev    + sconr
        zqit =   0.0_wp
        zqst =   0.0_wp
        zqgt =   0.0_wp
      END IF
      zqvt =  -zqct   - zqrt   - zqit   - zqst   - zqgt

      qice  =  qig + qsg + qgg + (zqit+zqst+zqgt)*zdt
      qliq  =  qcg + qrg + (zqct+zqrt)*zdt
      qtot  =  qvg + qcg + qrg + qig + qsg + qgg
      cv    =  cvd + (cvv-cvd)*qtot + (clw-cvv)*qliq + (ci-cvv)*qice
      ztt   =  (-zqvt*tg*cvv                                   &
                +(zqct+zqrt     )*(alv-(cpv-clw)*tmelt-clw*tg) &
                +(zqit+zqst+zqgt)*(als-(cpv-ci )*tmelt-ci *tg) &
               ) / cv

      ! Update variables and add qi to qrs for water loading
      qig = MAX ( 0.0_wp, (zzai*z1orhog + zqit*zdt)*zimi)
      qrg = MAX ( 0.0_wp, (zzar*z1orhog + zqrt*zdt)*zimr)
      qsg = MAX ( 0.0_wp, (zzas*z1orhog + zqst*zdt)*zims)
      qgg = MAX ( 0.0_wp, (zzag*z1orhog + zqgt*zdt)*zimg)
      
      !----------------------------------------------------------------------
      ! Section 10: Complete time step
      !----------------------------------------------------------------------

      IF ( k /= ke) THEN
        ! Store precipitation fluxes and sedimentation velocities 
        ! for the next level
        zprvr(iv) = qrg*rhog*zvzr(iv)
        zprvs(iv) = qsg*rhog*zvzs(iv)
        zprvg(iv) = qgg*rhog*zvzg(iv)
        zprvi(iv) = qig*rhog*zvzi(iv)
        IF (zprvr(iv) .LE. zqmin) zprvr(iv)=0.0_wp
        IF (zprvs(iv) .LE. zqmin) zprvs(iv)=0.0_wp
        IF (zprvg(iv) .LE. zqmin) zprvg(iv)=0.0_wp
        IF (zprvi(iv) .LE. zqmin) zprvi(iv)=0.0_wp

        qrsflux(iv,k) = zprvr(iv)+zprvs(iv)+zprvg(iv)+zprvi(iv)
        qrsflux(iv,k) = 0.5_wp*(qrsflux(iv,k)+zpkr(iv)+zpks(iv)+zpkg(iv)+zpki(iv))

        IF (qrg+qr(iv,k+1) <= zqmin .OR. .NOT. lrain) THEN
          zvzr(iv)= 0.0_wp
        ELSE
          zvzr(iv)= zvz0r * EXP(zvzxp*LOG((qrg+qr(iv,k+1))*0.5_wp*rhog)) * zrho1o2
        ENDIF
        IF (qsg+qs(iv,k+1) <= zqmin .OR. .NOT. lrain) THEN
          zvzs(iv)= 0.0_wp
        ELSE
          zvzs(iv)= zvz0s * EXP(ccswxp*LOG((qsg+qs(iv,k+1))*0.5_wp*rhog)) * zrho1o2
        ENDIF
        IF (qgg+qg(iv,k+1) <= zqmin .OR. .NOT. lrain) THEN
          zvzg(iv)= 0.0_wp
        ELSE
          zvzg(iv)=zvz0g * EXP(zexpsedg*LOG((qgg+qg(iv,k+1))*0.5_wp*rhog)) * zrho1o2
        ENDIF
        IF (qig+qi(iv,k+1) <= zqmin .OR. .NOT. lrain) THEN
          zvzi(iv)= 0.0_wp
        ELSE
          zvzi(iv)= zvz0i * EXP(zbvi*LOG((qig+qi(iv,k+1))*0.5_wp*rhog)) * zrhofac_qi
        ENDIF
          
      ELSE
        ! Precipitation fluxes at the ground
        IF (lrain) THEN
           prr_gsp(iv) = 0.5_wp * (qrg*rhog*zvzr(iv) + zpkr(iv))
           pri_gsp(iv) = 0.5_wp * (qig*rhog*zvzi(iv) + zpki(iv))
           prs_gsp(iv) = 0.5_wp * (qsg*rhog*zvzs(iv) + zpks(iv))
           prg_gsp(iv) = 0.5_wp * (qgg*rhog*zvzg(iv) + zpkg(iv))
        END IF
        qrsflux(iv,k) = prr_gsp(iv)+pri_gsp(iv)+prs_gsp(iv)+prg_gsp(iv)
      ENDIF

      ! Update of prognostic variables or tendencies
      qr (iv,k) = MAX ( 0.0_wp, qrg )
      qs (iv,k) = MAX ( 0.0_wp, qsg )
      qi (iv,k) = MAX ( 0.0_wp, qig )
      qg (iv,k) = MAX ( 0.0_wp, qgg )
      t  (iv,k) = t (iv,k) + ztt*zdt 
      qv (iv,k) = MAX ( 0.0_wp, qv(iv,k) + zqvt*zdt )
      qc (iv,k) = MAX ( 0.0_wp, qc(iv,k) + zqct*zdt )

    ENDDO  !loop over iv

  END DO ! loop over levels
  !$ACC END PARALLEL
  !$ACC WAIT
  !$ACC END DATA

END SUBROUTINE graupel

END MODULE mo_aes_graupel
