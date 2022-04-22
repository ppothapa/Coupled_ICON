!>
!! @brief Subroutine interface_cloud_mig calls NWP graupel scheme
!!
!! @author Monika Esch, MPI-M
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

MODULE mo_interface_cloud_mig_3

  USE mo_kind                ,ONLY: wp
  USE mtime                  ,ONLY: datetime
  USE mo_copy                ,ONLY: copy

  USE mo_run_config          ,ONLY: iqv, iqc, iqi, iqr, iqs, iqg, msg_level

  USE mo_aes_phy_config      ,ONLY: aes_phy_config
  USE mo_aes_phy_memory      ,ONLY: t_aes_phy_field, t_aes_phy_tend, &
       &                                    prm_field,         prm_tend

  USE mo_cloud_mig_types     ,ONLY: t_cloud_mig_input, t_cloud_mig_output
  USE mo_cloud_mig_memory    ,ONLY:   cloud_mig_input,   cloud_mig_output
  USE mo_cloud_mig           ,ONLY:   cloud_mig

  USE mo_timer               ,ONLY: ltimer, timer_start, timer_stop, &
       &                            timer_mig, timer_cld_mig

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: interface_cloud_mig_3

CONTAINS

  SUBROUTINE interface_cloud_mig_3(jg,jb,jc             ,&
       &                           nproma,nlev,ntracer  ,& 
       &                           is_in_sd_ed_interval ,&
       &                           is_active            ,&
       &                           datetime_old         ,&
       &                           pdtime               )

    ! Arguments
    !
    INTEGER                 ,INTENT(in) :: jg,jb,jc
    INTEGER                 ,INTENT(in) :: nproma,nlev,ntracer
    LOGICAL                 ,INTENT(in) :: is_in_sd_ed_interval
    LOGICAL                 ,INTENT(in) :: is_active
    TYPE(datetime)          ,POINTER    :: datetime_old
    REAL(wp)                ,INTENT(in) :: pdtime

    ! Pointers
    !
    ! to aes_phy_memory
    TYPE(t_aes_phy_field) ,POINTER    :: field
    TYPE(t_aes_phy_tend)  ,POINTER    :: tend
    !
    ! to cloud_mig_memory
    TYPE(t_cloud_mig_input ),POINTER    :: input
    TYPE(t_cloud_mig_output),POINTER    :: output

    ! Local variables
    !
    LOGICAL  :: lparamcpl
    INTEGER  :: fc_mig
    !
    REAL(wp) :: tend_ta_mig(nlev)
    REAL(wp) :: tend_qv_mig(nlev)
    REAL(wp) :: tend_qc_mig(nlev)
    REAL(wp) :: tend_qi_mig(nlev)
    REAL(wp) :: tend_qr_mig(nlev)
    REAL(wp) :: tend_qs_mig(nlev)
    REAL(wp) :: tend_qg_mig(nlev)
    !
    REAL(wp) :: pr_rain
    REAL(wp) :: pr_snow
    REAL(wp) :: pr_grpl

    INTEGER  :: jk, jks, jke

    IF (ltimer) call timer_start(timer_mig)

    ! associate pointers
    !
    ! to aes_phy_memory
    field  => prm_field(jg)
    tend   => prm_tend (jg)
    !
    ! to cloud_mig memory
    input  => cloud_mig_input (jg)
    output => cloud_mig_output(jg)

    lparamcpl = aes_phy_config(jg)% lparamcpl
    fc_mig    = aes_phy_config(jg)% fc_mig

    jks       = aes_phy_config(jg)% jks_cloudy
    jke       = nlev

    ! store input in memory for output
    !
    ! input parameters
    !
    IF (ASSOCIATED(input% jg        )) CALL copy(jg       , input% jg        (jc,jb))
    IF (ASSOCIATED(input% msg_level )) CALL copy(msg_level, input% msg_level (jc,jb))
    IF (ASSOCIATED(input% pdtime    )) CALL copy(pdtime   , input% pdtime    (jc,jb))
    !
    ! input fields
    !
    IF (ASSOCIATED(input% dz    )) CALL copy(jks,jke, field% dz        (jc,:,jb)    , input% dz    (jc,:,jb))
    IF (ASSOCIATED(input% rho   )) CALL copy(jks,jke, field% rho       (jc,:,jb)    , input% rho   (jc,:,jb))
    IF (ASSOCIATED(input% pf    )) CALL copy(jks,jke, field% pfull     (jc,:,jb)    , input% pf    (jc,:,jb))
    IF (ASSOCIATED(input% cpair )) CALL copy(jks,jke, field% cpair     (jc,:,jb)    , input% cpair (jc,:,jb))
    IF (ASSOCIATED(input% ta    )) CALL copy(jks,jke, field% ta        (jc,:,jb)    , input% ta    (jc,:,jb))
    IF (ASSOCIATED(input% qv    )) CALL copy(jks,jke, field% qtrc      (jc,:,jb,iqv), input% qv    (jc,:,jb))
    IF (ASSOCIATED(input% qc    )) CALL copy(jks,jke, field% qtrc      (jc,:,jb,iqc), input% qc    (jc,:,jb))
    IF (ASSOCIATED(input% qi    )) CALL copy(jks,jke, field% qtrc      (jc,:,jb,iqi), input% qi    (jc,:,jb))
    IF (ASSOCIATED(input% qr    )) CALL copy(jks,jke, field% qtrc      (jc,:,jb,iqr), input% qr    (jc,:,jb))
    IF (ASSOCIATED(input% qs    )) CALL copy(jks,jke, field% qtrc      (jc,:,jb,iqs), input% qs    (jc,:,jb))
    IF (ASSOCIATED(input% qg    )) CALL copy(jks,jke, field% qtrc      (jc,:,jb,iqg), input% qg    (jc,:,jb))

    IF ( is_in_sd_ed_interval ) THEN
       !
       IF ( is_active ) THEN

          IF (ltimer) CALL timer_start(timer_cld_mig)
          !
          CALL cloud_mig( jg                                   ,& !< in : grid index
               &          msg_level                            ,& !< in : message level 
               &          pdtime                               ,& !< in : timestep
               &          field% dz        (jc,jks:jke,jb)     ,& !< in : vertical layer thickness
               &          field% rho       (jc,jks:jke,jb)     ,& !< in : density
               &          field% pfull     (jc,jks:jke,jb)     ,& !< in : pressure
               &          field% cpair     (jc,jks:jke,jb)     ,& !< in : specific heat of air
               &          field% ta        (jc,jks:jke,jb)     ,& !< in : temperature
               &          field% qtrc      (jc,jks:jke,jb,iqv) ,& !< in : sp humidity
               &          field% qtrc      (jc,jks:jke,jb,iqc) ,& !< in : cloud water
               &          field% qtrc      (jc,jks:jke,jb,iqi) ,& !< in : ice
               &          field% qtrc      (jc,jks:jke,jb,iqr) ,& !< in : rain
               &          field% qtrc      (jc,jks:jke,jb,iqs) ,& !< in : snow
               &          field% qtrc      (jc,jks:jke,jb,iqg) ,& !< in : graupel
               &          tend_ta_mig         (jks:jke)        ,& !< out: tendency of temperature
               &          tend_qv_mig         (jks:jke)        ,& !< out: tendency of water vapor
               &          tend_qc_mig         (jks:jke)        ,& !< out: tendency of cloud water
               &          tend_qi_mig         (jks:jke)        ,& !< out: tendency of cloud ice
               &          tend_qr_mig         (jks:jke)        ,& !< out: tendency of rain
               &          tend_qs_mig         (jks:jke)        ,& !< out: tendency of snow
               &          tend_qg_mig         (jks:jke)        ,& !< out: tendency of graupel 
               &          pr_rain                              ,& !& out: precip rate rain
               &          pr_snow                              ,& !& out: precip rate snow
               &          pr_grpl                              )  !& out: precip rate graupel
          !
          IF (ltimer) CALL timer_stop(timer_cld_mig)

          field% rsfl(jc,jb) = pr_rain
          field% ssfl(jc,jb) = pr_snow + pr_grpl

          ! store output in memory for output or recycling
          !
          IF (ASSOCIATED(output% tend_ta_mig )) CALL copy(jks,jke, tend_ta_mig (:), output% tend_ta_mig (jc,:,jb))
          IF (ASSOCIATED(output% tend_qv_mig )) CALL copy(jks,jke, tend_qv_mig (:), output% tend_qv_mig (jc,:,jb))
          IF (ASSOCIATED(output% tend_qc_mig )) CALL copy(jks,jke, tend_qc_mig (:), output% tend_qc_mig (jc,:,jb))
          IF (ASSOCIATED(output% tend_qi_mig )) CALL copy(jks,jke, tend_qi_mig (:), output% tend_qi_mig (jc,:,jb))
          IF (ASSOCIATED(output% tend_qr_mig )) CALL copy(jks,jke, tend_qr_mig (:), output% tend_qr_mig (jc,:,jb))
          IF (ASSOCIATED(output% tend_qs_mig )) CALL copy(jks,jke, tend_qs_mig (:), output% tend_qs_mig (jc,:,jb))
          IF (ASSOCIATED(output% tend_qg_mig )) CALL copy(jks,jke, tend_qg_mig (:), output% tend_qg_mig (jc,:,jb))
          !
          IF (ASSOCIATED(output% pr_rain     )) CALL copy(pr_rain, output% pr_rain (jc,jb))
          IF (ASSOCIATED(output% pr_snow     )) CALL copy(pr_snow, output% pr_snow (jc,jb))
          IF (ASSOCIATED(output% pr_grpl     )) CALL copy(pr_grpl, output% pr_grpl (jc,jb))
          !
       ELSE    ! is_active
          !
          ! retrieve output from memory for recycling
          !
          IF (ASSOCIATED(output% tend_ta_mig )) CALL copy(jks,jke, output% tend_ta_mig (jc,:,jb), tend_ta_mig (:))
          IF (ASSOCIATED(output% tend_qv_mig )) CALL copy(jks,jke, output% tend_qv_mig (jc,:,jb), tend_qv_mig (:))
          IF (ASSOCIATED(output% tend_qc_mig )) CALL copy(jks,jke, output% tend_qc_mig (jc,:,jb), tend_qc_mig (:))
          IF (ASSOCIATED(output% tend_qi_mig )) CALL copy(jks,jke, output% tend_qi_mig (jc,:,jb), tend_qi_mig (:))
          IF (ASSOCIATED(output% tend_qr_mig )) CALL copy(jks,jke, output% tend_qr_mig (jc,:,jb), tend_qr_mig (:))
          IF (ASSOCIATED(output% tend_qs_mig )) CALL copy(jks,jke, output% tend_qs_mig (jc,:,jb), tend_qs_mig (:))
          IF (ASSOCIATED(output% tend_qg_mig )) CALL copy(jks,jke, output% tend_qg_mig (jc,:,jb), tend_qg_mig (:))
          !
       END IF  ! is_active
       !
       ! accumulate tendencies for later updating the model state
       SELECT CASE(fc_mig)
       CASE(0)
          ! diagnostic, do not use tendency
       CASE(1)
          DO jk = jks,jke
             ! use tendency to update the model state
             tend%   ta_phy(jc,jk,jb)      = tend%   ta_phy(jc,jk,jb)     + tend_ta_mig(jk)
             tend% qtrc_phy(jc,jk,jb,iqv)  = tend% qtrc_phy(jc,jk,jb,iqv) + tend_qv_mig(jk)
             tend% qtrc_phy(jc,jk,jb,iqc)  = tend% qtrc_phy(jc,jk,jb,iqc) + tend_qc_mig(jk)
             tend% qtrc_phy(jc,jk,jb,iqi)  = tend% qtrc_phy(jc,jk,jb,iqi) + tend_qi_mig(jk)
             tend% qtrc_phy(jc,jk,jb,iqr)  = tend% qtrc_phy(jc,jk,jb,iqr) + tend_qr_mig(jk)
             tend% qtrc_phy(jc,jk,jb,iqs)  = tend% qtrc_phy(jc,jk,jb,iqs) + tend_qs_mig(jk)
             tend% qtrc_phy(jc,jk,jb,iqg)  = tend% qtrc_phy(jc,jk,jb,iqg) + tend_qg_mig(jk)
          END DO
!!$       CASE(2)
!!$          ! use tendency as forcing in the dynamics
!!$          ...
       END SELECT
       !
       ! update physics state for input to the next physics process
       IF (lparamcpl) THEN
         SELECT CASE(fc_mig)
         CASE(0)
            ! diagnostic, do not use tendency
         CASE(1)
             DO jk = jks,jke
                field%   ta(jc,jk,jb)      = field%   ta(jc,jk,jb)      + tend_ta_mig(jk)*pdtime
                field% qtrc(jc,jk,jb,iqv)  = field% qtrc(jc,jk,jb,iqv)  + tend_qv_mig(jk)*pdtime
                field% qtrc(jc,jk,jb,iqc)  = field% qtrc(jc,jk,jb,iqc)  + tend_qc_mig(jk)*pdtime
                field% qtrc(jc,jk,jb,iqi)  = field% qtrc(jc,jk,jb,iqi)  + tend_qi_mig(jk)*pdtime
                field% qtrc(jc,jk,jb,iqr)  = field% qtrc(jc,jk,jb,iqr)  + tend_qr_mig(jk)*pdtime
                field% qtrc(jc,jk,jb,iqs)  = field% qtrc(jc,jk,jb,iqs)  + tend_qs_mig(jk)*pdtime
                field% qtrc(jc,jk,jb,iqg)  = field% qtrc(jc,jk,jb,iqg)  + tend_qg_mig(jk)*pdtime
             END DO
             !
!!$         CASE(2)
!!$            ! use tendency as forcing in the dynamics
!!$            ...
         END SELECT
       END IF
       !
    ELSE       ! is_in_sd_ed_interval
       !
       ! initialize
       !
       IF (ASSOCIATED(output% tend_ta_mig )) CALL copy(jks,jke, 0._wp, output% tend_ta_mig (jc,:,jb))
       IF (ASSOCIATED(output% tend_qv_mig )) CALL copy(jks,jke, 0._wp, output% tend_qv_mig (jc,:,jb))
       IF (ASSOCIATED(output% tend_qc_mig )) CALL copy(jks,jke, 0._wp, output% tend_qc_mig (jc,:,jb))
       IF (ASSOCIATED(output% tend_qi_mig )) CALL copy(jks,jke, 0._wp, output% tend_qi_mig (jc,:,jb))
       IF (ASSOCIATED(output% tend_qr_mig )) CALL copy(jks,jke, 0._wp, output% tend_qr_mig (jc,:,jb))
       IF (ASSOCIATED(output% tend_qs_mig )) CALL copy(jks,jke, 0._wp, output% tend_qs_mig (jc,:,jb))
       IF (ASSOCIATED(output% tend_qg_mig )) CALL copy(jks,jke, 0._wp, output% tend_qg_mig (jc,:,jb))
       !
       IF (ASSOCIATED(output% pr_rain )) CALL copy(0._wp, output% pr_rain (jc,jb))
       IF (ASSOCIATED(output% pr_snow )) CALL copy(0._wp, output% pr_snow (jc,jb))
       IF (ASSOCIATED(output% pr_grpl )) CALL copy(0._wp, output% pr_grpl (jc,jb))
       !
    END IF     ! is_in_sd_ed_interval

    ! disassociate pointers
    NULLIFY(field, tend)
    NULLIFY(input, output)

    IF (ltimer) call timer_stop(timer_mig)

  END SUBROUTINE interface_cloud_mig_3

END MODULE mo_interface_cloud_mig_3
