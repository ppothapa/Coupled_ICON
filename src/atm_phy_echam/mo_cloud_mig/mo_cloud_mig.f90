!>
!! @brief Subroutine cloud_mig calls the saturation adjustment and cloud scheme
!!
!! @author Marco Giorgetta, MPI-M, 2020-01
!!
!! @par Revision History
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

MODULE mo_cloud_mig

  USE mo_kind                ,ONLY: wp

  USE mo_timer               ,ONLY: ltimer, timer_start, timer_stop, &
       &                            timer_sat, timer_grp

  USE mo_echam_phy_config    ,ONLY: echam_phy_config
  USE mo_cloud_mig_config    ,ONLY: cloud_mig_config
  USE mo_physical_constants  ,ONLY: cvd
  USE mo_satad               ,ONLY: satad_v_3d
  USE mo_satad_v_el          ,ONLY: satad_v_el
  USE mo_satad_v_1col        ,ONLY: satad_v_1col
  USE mo_satad_v_1cell       ,ONLY: satad_v_1cell
  USE gscp_data              ,ONLY: cloud_num
  USE gscp_graupel           ,ONLY: graupel
  USE gscp_graupel_1col      ,ONLY: graupel_1col

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: cloud_mig

  INTERFACE cloud_mig
     MODULE PROCEDURE cloud_mig_1col
     MODULE PROCEDURE cloud_mig_ncol
  END INTERFACE cloud_mig

CONTAINS

  SUBROUTINE cloud_mig_ncol( jg         ,&
       &                     jcs, jce   ,&
       &                     msg_level  ,&
       &                     pdtime     ,&
       &                     dz         ,&
       &                     rho        ,&
       &                     pf         ,&
       &                     cpair      ,&
       &                     ta         ,&
       &                     qv         ,&
       &                     qc         ,&
       &                     qi         ,&
       &                     qr         ,&
       &                     qs         ,&
       &                     qg         ,&
       &                     tend_ta    ,&
       &                     tend_qv    ,&
       &                     tend_qc    ,&
       &                     tend_qi    ,&
       &                     tend_qr    ,&
       &                     tend_qs    ,&
       &                     tend_qg    ,&
       &                     pr_rain    ,&
       &                     pr_snow    ,&
       &                     pr_grpl    )

    ! Arguments
    !
    INTEGER , INTENT(in)  :: jg            !< grid index
    INTEGER , INTENT(in)  :: jcs, jce      !< column index range
    INTEGER , INTENT(in)  :: msg_level     !< message level
    REAL(wp), INTENT(in)  :: pdtime        !< timestep
    !
    REAL(wp), INTENT(in)  :: dz      (:,:) !< vertical layer thickness
    REAL(wp), INTENT(in)  :: rho     (:,:) !< density
    REAL(wp), INTENT(in)  :: pf      (:,:) !< pressure
    REAL(wp), INTENT(in)  :: cpair   (:,:) !< specific heat of air
    !
    REAL(wp), INTENT(in)  :: ta      (:,:) !< temperature
    REAL(wp), INTENT(in)  :: qv      (:,:) !< sp humidity
    REAL(wp), INTENT(in)  :: qc      (:,:) !< cloud water
    REAL(wp), INTENT(in)  :: qi      (:,:) !< ice
    REAL(wp), INTENT(in)  :: qr      (:,:) !< rain
    REAL(wp), INTENT(in)  :: qs      (:,:) !< snow
    REAL(wp), INTENT(in)  :: qg      (:,:) !< graupel
    !
    REAL(wp), INTENT(out) :: tend_ta (:,:) !< tendency of temperature
    REAL(wp), INTENT(out) :: tend_qv (:,:) !< tendency of water vapor
    REAL(wp), INTENT(out) :: tend_qc (:,:) !< tendency of cloud water
    REAL(wp), INTENT(out) :: tend_qi (:,:) !< tendency of cloud ice
    REAL(wp), INTENT(out) :: tend_qr (:,:) !< tendency of rain
    REAL(wp), INTENT(out) :: tend_qs (:,:) !< tendency of snow
    REAL(wp), INTENT(out) :: tend_qg (:,:) !< tendency of graupel

    REAL(wp), INTENT(out) :: pr_rain (:)   !< precip rate rain
    REAL(wp), INTENT(out) :: pr_snow (:)   !< precip rate snow
    REAL(wp), INTENT(out) :: pr_grpl (:)   !< precip rate graupel

    ! Local variables
    !
    INTEGER  :: nproma
    INTEGER  :: jc
    INTEGER  :: jk, jks, jke
    !
    REAL(wp) :: qi0, qc0
    REAL(wp) :: zqnc(SIZE(dz,1))
    !
    REAL(wp) :: zta(SIZE(dz,1),SIZE(dz,2))
    REAL(wp) :: zqv(SIZE(dz,1),SIZE(dz,2))
    REAL(wp) :: zqc(SIZE(dz,1),SIZE(dz,2))
    REAL(wp) :: zqi(SIZE(dz,1),SIZE(dz,2))
    REAL(wp) :: zqr(SIZE(dz,1),SIZE(dz,2))
    REAL(wp) :: zqs(SIZE(dz,1),SIZE(dz,2))
    REAL(wp) :: zqg(SIZE(dz,1),SIZE(dz,2))
    !
    REAL(wp) :: zqrsflux(SIZE(dz,1),SIZE(dz,2))
    !
    REAL(wp) :: zdtr ! reciprocal of timestep

    !$ACC DATA PRESENT( dz, rho, pf, cpair, ta, qv, qc, qi, qr, qs, qg,                &
    !$ACC               tend_ta, tend_qv, tend_qc, tend_qi, tend_qr, tend_qs, tend_qg, &
    !$ACC               pr_rain, pr_snow, pr_grpl )                                    &
    !$ACC       CREATE( zqnc, zta, zqv, zqc, zqi, zqr, zqs, zqg, zqrsflux )

    nproma = SIZE(dz,1)

    jks = 1
    jke = (SIZE(dz,2))

    qi0 = cloud_mig_config(jg)% qi0
    qc0 = cloud_mig_config(jg)% qc0

    !$ACC PARALLEL DEFAULT(NONE)
    !$ACC LOOP GANG VECTOR
    DO jc = jcs,jce
       zqnc(jc) = cloud_num
    END DO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(NONE)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = jks,jke
       DO jc = jcs,jce
          zta(jc,jk) = ta(jc,jk)
          zqv(jc,jk) = qv(jc,jk)
          zqc(jc,jk) = qc(jc,jk)
          zqi(jc,jk) = qi(jc,jk)
          zqr(jc,jk) = qr(jc,jk)
          zqs(jc,jk) = qs(jc,jk)
          zqg(jc,jk) = qg(jc,jk)
       END DO
    END DO
    !$ACC END PARALLEL

    zdtr = 1._wp/pdtime

    ! Initial saturation adjustment
    !
    IF (ltimer) call timer_start(timer_sat)
    !
    SELECT CASE(echam_phy_config(jg)%if_mig)
    CASE(1)
       CALL satad_v_3d( maxiter  = 10              ,& !> in
            &           idim     = nproma          ,& !> in
            &           kdim     = jke             ,& !> in
            &           ilo      = jcs             ,& !> in
            &           iup      = jce             ,& !> in
            &           klo      = jks             ,& !> in
            &           kup      = jke             ,& !> in
            &           tol      = 1.e-3_wp        ,& !> in
            &           te       = zta       (:,:) ,& !> inout
            &           qve      = zqv       (:,:) ,& !> inout
            &           qce      = zqc       (:,:) ,& !> inout
            &           rhotot   = rho       (:,:) )
    CASE(11)
       CALL satad_v_el( rho(jcs:jce,jks:jke) ,& !> in
            &           zta(jcs:jce,jks:jke) ,& !> inout
            &           zqv(jcs:jce,jks:jke) ,& !> inout
            &           zqc(jcs:jce,jks:jke) )  !> inout
    END SELECT
          !
    IF (ltimer) call timer_stop(timer_sat)

    ! Single moment cloud microphyiscs for water vapor,
    ! cloud water, cloud ice, rain, snow and graupel
    !
    IF (ltimer) call timer_start(timer_grp)
    !
    CALL graupel( nvec    = nproma        ,& !< in
         &        ke      = jke           ,& !< in
         &        ivstart = jcs           ,& !< in
         &        ivend   = jce           ,& !< in
         &        kstart  = jks           ,& !< in
         &        idbg    = msg_level     ,& !< in   : message level 
         &        l_cv    = .TRUE.        ,& !< in   : if temp. is changed for const. volumes
!!!      & ithermo_water=atm_phy_nwp_config(jg)%ithermo_water) ,& !< in: latent heat choice
!!!                Not yet in use in sapphire physics part
         &        zdt     = pdtime        ,& !< in   : timestep
         &        qi0     = qi0           ,& !< in   : cloud ice threshold for autoconversion
         &        qc0     = qc0           ,& !< in   : cloud water threshold for autoconversion
         &        qnc     = zqnc    (:)   ,& !< in
         !
         &        dz      = dz      (:,:) ,& !< in   : vertical layer thickness
         &        rho     = rho     (:,:) ,& !< in   : density
         &        p       = pf      (:,:) ,& !< in   : pressure
         !
         &        t       = zta     (:,:) ,& !< inout: temp
         &        qv      = zqv     (:,:) ,& !< inout: sp humidity
         &        qc      = zqc     (:,:) ,& !< inout: cloud water
         &        qi      = zqi     (:,:) ,& !< inout: ice
         &        qr      = zqr     (:,:) ,& !< inout: rain
         &        qs      = zqs     (:,:) ,& !< inout: snow
         &        qg      = zqg     (:,:) ,& !< inout: graupel
         &        qrsflux = zqrsflux(:,:) ,& !<   out: precip flux in atmosphere
         !
         &        prr_gsp = pr_rain (:)   ,& !<   out: precip rate rain
         &        prs_gsp = pr_snow (:)   ,& !<   out: precip rate snow
         &        prg_gsp = pr_grpl (:)   )  !<   out: precip rate graupel
    !
    IF (ltimer) call timer_stop(timer_grp)

    ! Final saturation adjustment
    !
    IF (ltimer) call timer_start(timer_sat)
    !
    SELECT CASE(echam_phy_config(jg)%if_mig)
    CASE(1)
       CALL satad_v_3d( maxiter  = 10              ,& !> in
            &           idim     = nproma          ,& !> in
            &           kdim     = jke             ,& !> in
            &           ilo      = jcs             ,& !> in
            &           iup      = jce             ,& !> in
            &           klo      = jks             ,& !> in
            &           kup      = jke             ,& !> in
            &           tol      = 1.e-3_wp        ,& !> in
            &           te       = zta       (:,:) ,& !> inout
            &           qve      = zqv       (:,:) ,& !> inout
            &           qce      = zqc       (:,:) ,& !> inout
            &           rhotot   = rho       (:,:) )
    CASE(11)
       CALL satad_v_el( rho(jcs:jce,jks:jke) ,& !> in
            &           zta(jcs:jce,jks:jke) ,& !> inout
            &           zqv(jcs:jce,jks:jke) ,& !> inout
            &           zqc(jcs:jce,jks:jke) )  !> inout
    END SELECT
    !
    IF (ltimer) call timer_stop(timer_sat)

    ! Calculate tendencies and convert temperature tendency, as computed
    ! in satad/graupel for constant volume to constant pressure
    !
    !$ACC PARALLEL DEFAULT(NONE)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = jks,jke
       DO jc = jcs,jce
          tend_ta(jc,jk) =     (zta(jc,jk)-ta(jc,jk))*zdtr*cvd/cpair(jc,jk)
          tend_qv(jc,jk) = MAX((zqv(jc,jk)-qv(jc,jk))*zdtr,-qv(jc,jk)*zdtr)
          tend_qc(jc,jk) = MAX((zqc(jc,jk)-qc(jc,jk))*zdtr,-qc(jc,jk)*zdtr)
          tend_qi(jc,jk) = MAX((zqi(jc,jk)-qi(jc,jk))*zdtr,-qi(jc,jk)*zdtr)
          tend_qr(jc,jk) = MAX((zqr(jc,jk)-qr(jc,jk))*zdtr,-qr(jc,jk)*zdtr)
          tend_qs(jc,jk) = MAX((zqs(jc,jk)-qs(jc,jk))*zdtr,-qs(jc,jk)*zdtr)
          tend_qg(jc,jk) = MAX((zqg(jc,jk)-qg(jc,jk))*zdtr,-qg(jc,jk)*zdtr)
       END DO
    END DO
    !$ACC END PARALLEL

    !$ACC END DATA

  END SUBROUTINE cloud_mig_ncol

  !-------------------------

  SUBROUTINE cloud_mig_1col( jg         ,&
       &                     msg_level  ,&
       &                     pdtime     ,&
       &                     dz         ,&
       &                     rho        ,&
       &                     pf         ,&
       &                     cpair      ,&
       &                     ta         ,&
       &                     qv         ,&
       &                     qc         ,&
       &                     qi         ,&
       &                     qr         ,&
       &                     qs         ,&
       &                     qg         ,&
       &                     tend_ta    ,&
       &                     tend_qv    ,&
       &                     tend_qc    ,&
       &                     tend_qi    ,&
       &                     tend_qr    ,&
       &                     tend_qs    ,&
       &                     tend_qg    ,&
       &                     pr_rain    ,&
       &                     pr_snow    ,&
       &                     pr_grpl    )

    INTEGER , INTENT(in)  :: jg          !< grid index
    INTEGER , INTENT(in)  :: msg_level   !< message level
    REAL(wp), INTENT(in)  :: pdtime      !< timestep
    !
    REAL(wp), INTENT(in)  :: dz      (:) !< vertical layer thickness
    REAL(wp), INTENT(in)  :: rho     (:) !< density
    REAL(wp), INTENT(in)  :: pf      (:) !< pressure
    REAL(wp), INTENT(in)  :: cpair   (:) !< specific heat of air
    !
    REAL(wp), INTENT(in)  :: ta      (:) !< temperature
    REAL(wp), INTENT(in)  :: qv      (:) !< sp humidity
    REAL(wp), INTENT(in)  :: qc      (:) !< cloud water
    REAL(wp), INTENT(in)  :: qi      (:) !< ice
    REAL(wp), INTENT(in)  :: qr      (:) !< rain
    REAL(wp), INTENT(in)  :: qs      (:) !< snow
    REAL(wp), INTENT(in)  :: qg      (:) !< graupel
    !
    REAL(wp), INTENT(out) :: tend_ta (:) !< tendency of temperature
    REAL(wp), INTENT(out) :: tend_qv (:) !< tendency of water vapor
    REAL(wp), INTENT(out) :: tend_qc (:) !< tendency of cloud water
    REAL(wp), INTENT(out) :: tend_qi (:) !< tendency of cloud ice
    REAL(wp), INTENT(out) :: tend_qr (:) !< tendency of rain
    REAL(wp), INTENT(out) :: tend_qs (:) !< tendency of snow
    REAL(wp), INTENT(out) :: tend_qg (:) !< tendency of graupel
    !
    REAL(wp), INTENT(out) :: pr_rain     !< precip rate rain
    REAL(wp), INTENT(out) :: pr_snow     !< precip rate snow
    REAL(wp), INTENT(out) :: pr_grpl     !< precip rate graupel

    ! Local variables
    !
    INTEGER  :: jk, jks, jke
    !
    REAL(wp) :: qi0, qc0
    !
    REAL(wp) :: zta(SIZE(dz))
    REAL(wp) :: zqv(SIZE(dz))
    REAL(wp) :: zqc(SIZE(dz))
    REAL(wp) :: zqi(SIZE(dz))
    REAL(wp) :: zqr(SIZE(dz))
    REAL(wp) :: zqs(SIZE(dz))
    REAL(wp) :: zqg(SIZE(dz))
    REAL(wp) :: zqnc(SIZE(dz))
    !
    REAL(wp) :: zqrsflux(SIZE(dz))
    !
    REAL(wp) :: zdtr ! reciprocal of timestep

    jks = 1
    jke = (SIZE(dz))

    qi0 = cloud_mig_config(jg)% qi0
    qc0 = cloud_mig_config(jg)% qc0

    DO jk = jks,jke
       zta(jk) = ta(jk)
       zqv(jk) = qv(jk)
       zqc(jk) = qc(jk)
       zqi(jk) = qi(jk)
       zqr(jk) = qr(jk)
       zqs(jk) = qs(jk)
       zqg(jk) = qg(jk)
       zqnc(jk) = cloud_num
    END DO

    zdtr = 1._wp/pdtime

    ! Initial saturation adjustment
    !
    IF (ltimer) call timer_start(timer_sat)
    !
    SELECT CASE(echam_phy_config(jg)%if_mig)
    CASE(2,3)
       CALL satad_v_1col( rho(:) ,& !> in
            &             zta(:) ,& !> inout
            &             zqv(:) ,& !> inout
            &             zqc(:) )  !> inout
    CASE(4)
       DO jk = jks,jke
          CALL satad_v_1cell( rho(jk) ,& !> in
               &              zta(jk) ,& !> inout
               &              zqv(jk) ,& !> inout
               &              zqc(jk) )  !> inout
       END DO
    CASE(12,13)
       CALL satad_v_el( rho(:) ,& !> in
            &           zta(:) ,& !> inout
            &           zqv(:) ,& !> inout
            &           zqc(:) )  !> inout
    CASE(14)
       DO jk = jks,jke
          CALL satad_v_el( rho(jk) ,& !> in
               &           zta(jk) ,& !> inout
               &           zqv(jk) ,& !> inout
               &           zqc(jk) )  !> inout
       END DO
    END SELECT
    !
    IF (ltimer) call timer_stop(timer_sat)

    ! Single moment cloud microphyiscs for water vapor,
    ! cloud water, cloud ice, rain, snow and graupel
    !
    IF (ltimer) call timer_start(timer_grp)
    !
    CALL graupel_1col( msg_level   ,& !< in   : message level 
         &             pdtime      ,& !< in   : timestep
         &             qi0         ,& !< in   : cloud ice threshold for autoconversion
         &             qc0         ,& !< in   : cloud water threshold for autoconversion
         &             .TRUE.      ,& !< in   : if temp. is changed for const. volumes
!!!      & ithermo_water=atm_phy_nwp_config(jg)%ithermo_water) ,& !< in: latent heat choice
!!!                Not yet in use in sapphire physics part
         !
         &             dz      (:) ,& !< in   : vertical layer thickness
         &             rho     (:) ,& !< in   : density
         &             pf      (:) ,& !< in   : pressure
         !
         &             zta     (:) ,& !< inout: temp
         &             zqv     (:) ,& !< inout: sp humidity
         &             zqc     (:) ,& !< inout: cloud water
         &             zqi     (:) ,& !< inout: ice
         &             zqr     (:) ,& !< inout: rain
         &             zqs     (:) ,& !< inout: snow
         &             zqg     (:) ,& !< inout: graupel
         &             zqnc    (:) ,& !< in
         &             zqrsflux(:) ,& !<   out: precip flux in atmosphere
         !
         &             pr_rain     ,& !<   out: precip rate rain
         &             pr_snow     ,& !<   out: precip rate snow
         &             pr_grpl     )  !<   out: precip rate graupel
    !
    IF (ltimer) call timer_stop(timer_grp)

    ! Final saturation adjustment
    !
    IF (ltimer) call timer_start(timer_sat)
    !
    SELECT CASE(echam_phy_config(jg)%if_mig)
    CASE(2,3)
       CALL satad_v_1col( rho(:) ,& !> in
            &             zta(:) ,& !> inout
            &             zqv(:) ,& !> inout
            &             zqc(:) )  !> inout
    CASE(4)
       DO jk = jks,jke
          CALL satad_v_1cell( rho(jk) ,& !> in
               &              zta(jk) ,& !> inout
               &              zqv(jk) ,& !> inout
               &              zqc(jk) )  !> inout
       END DO
    CASE(12,13)
       CALL satad_v_el( rho(:) ,& !> in
            &           zta(:) ,& !> inout
            &           zqv(:) ,& !> inout
            &           zqc(:) )  !> inout
    CASE(14)
       DO jk = jks,jke
          CALL satad_v_el( rho(jk) ,& !> in
               &           zta(jk) ,& !> inout
               &           zqv(jk) ,& !> inout
               &           zqc(jk) )  !> inout
       END DO
    END SELECT
    !
    IF (ltimer) call timer_stop(timer_sat)

    ! Calculate tendencies and convert temperature tendency, as computed
    ! in satad/graupel for constant volume to constant pressure
    !
    DO jk = jks,jke
       tend_ta(jk) =     (zta(jk)-ta(jk))*zdtr*cvd/cpair(jk)
       tend_qv(jk) = MAX((zqv(jk)-qv(jk))*zdtr,-qv(jk)*zdtr)
       tend_qc(jk) = MAX((zqc(jk)-qc(jk))*zdtr,-qc(jk)*zdtr)
       tend_qi(jk) = MAX((zqi(jk)-qi(jk))*zdtr,-qi(jk)*zdtr)
       tend_qr(jk) = MAX((zqr(jk)-qr(jk))*zdtr,-qr(jk)*zdtr)
       tend_qs(jk) = MAX((zqs(jk)-qs(jk))*zdtr,-qs(jk)*zdtr)
       tend_qg(jk) = MAX((zqg(jk)-qg(jk))*zdtr,-qg(jk)*zdtr)
    END DO

  END SUBROUTINE cloud_mig_1col

END MODULE mo_cloud_mig
