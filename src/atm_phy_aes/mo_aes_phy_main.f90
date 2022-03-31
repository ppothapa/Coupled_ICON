!>
!! @brief Subroutine aes_phy_main calls all the parameterization schemes
!!
!! @author Hui Wan, MPI-M
!! @author Marco Giorgetta, MPI-M
!!
!! @par Revision History
!!  Original version from ECHAM6 (revision 2028)
!!  Modified for ICOHAM by Hui Wan and Marco Giorgetta (2010)
!!  Modified for ICONAM by Marco Giorgetta (2014)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

#if defined __xlC__ && !defined NOXLFPROCESS
@PROCESS HOT
@PROCESS SPILLSIZE(5000)
#endif
!OCL NOALIAS

MODULE mo_aes_phy_main

  USE mo_kind                ,ONLY: wp
  USE mo_exception           ,ONLY: message, warning
  USE mtime                  ,ONLY: datetime, isCurrentEventActive, &
       &                            OPERATOR(<=), OPERATOR(>)

  USE mo_model_domain        ,ONLY: t_patch

  USE mo_omp_loop            ,ONLY: omp_loop_cell, omp_loop_cell_3

  USE mo_aes_phy_config      ,ONLY: aes_phy_config, aes_phy_tc, dt_zero
  USE mo_aes_phy_diag        ,ONLY: surface_fractions, &
    &                               droplet_number,    &
    &                               cpair_cvair_qconv, &
    &                               initialize,        &
    &                               finalize

  USE mo_aes_diagnostics     ,ONLY: aes_global_diagnostics, aes_diag_output_minmax_micro
#if defined( _OPENACC )
  USE mo_var_list_gpu        ,ONLY: gpu_update_var_list
#endif
  USE mo_run_config          ,ONLY: msg_level


  USE mo_interface_aes_cov   ,ONLY: interface_aes_cov
  USE mo_interface_aes_wmo   ,ONLY: interface_aes_wmo
  USE mo_interface_aes_rad   ,ONLY: interface_aes_rad
  USE mo_interface_aes_rht   ,ONLY: interface_aes_rht
  USE mo_interface_aes_vdf   ,ONLY: interface_aes_vdf
  USE mo_interface_aes_car   ,ONLY: interface_aes_car
  !
  ! experimental interfaces selected by aes_phy_config(:)%if_mig = 1,2!!$,3
  USE mo_interface_cloud_mig_1 ,ONLY: interface_cloud_mig_1 ! use 2d interface, jcs:jce loop in cloud_mig
  USE mo_interface_cloud_mig_2 ,ONLY: interface_cloud_mig_2 ! use 2d interface, jcs:jce loop in interface
  USE mo_interface_cloud_mig_3 ,ONLY: interface_cloud_mig_3 ! use 1d interface, jcs:jce loop in omp_loop_cell_1d
  !
  USE mo_interface_cloud_two ,ONLY: interface_cloud_two

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: aes_phy_main

CONTAINS

  !>
  !!
  SUBROUTINE aes_phy_main(patch         ,&
    &                     datetime_old  ,&
    &                     pdtime        )


    ! Arguments
    !
    TYPE(t_patch)  ,TARGET ,INTENT(INOUT) :: patch
    TYPE(datetime)         ,POINTER       :: datetime_old
    REAL(wp)               ,INTENT(IN)    :: pdtime

    ! Local variables
    !
    INTEGER  :: jg                                         !< grid level/domain index

    LOGICAL  :: is_in_sd_ed_interval                       !< time is in process interval [sd,ed[
    LOGICAL  :: is_active                                  !< process is active

    jg = patch%id

    !-------------------------------------------------------------------
    ! Initialize (diagnostic)
    !-------------------------------------------------------------------
    !
    CALL omp_loop_cell(patch,initialize)

    !-------------------------------------------------------------------
    ! Specific heat of moist air (diagnostic)
    !-------------------------------------------------------------------
    !
    CALL omp_loop_cell(patch,cpair_cvair_qconv)

    !-------------------------------------------------------------------
    ! Calculate surface fraction (diagnostic)
    !-------------------------------------------------------------------
    !
    CALL omp_loop_cell(patch,surface_fractions)
 
    !-------------------------------------------------------------------
    ! Cloud cover (diagnostic)
    !-------------------------------------------------------------------
    !
    CALL omp_loop_cell(patch,interface_aes_cov)

    !---------------------------------------------------------------------
    ! 3.9 Determine tropopause height (diagnostic)
    !---------------------------------------------------------------------
    !
    CALL omp_loop_cell(patch,interface_aes_wmo)

    !---------------------------------------------------------------------
    ! Cloud droplet number concentration (diagnostic)
    ! used in radiation and cloud
    !---------------------------------------------------------------------
    !
    CALL omp_loop_cell(patch,droplet_number)

    !-------------------------------------------------------------------
    ! Radiation (one interface for LW+SW)
    !-------------------------------------------------------------------
    !
    IF ( aes_phy_tc(jg)%dt_rad > dt_zero ) THEN
       !
       is_in_sd_ed_interval =          (aes_phy_tc(jg)%sd_rad <= datetime_old) .AND. &
            &                          (aes_phy_tc(jg)%ed_rad >  datetime_old)
       is_active = isCurrentEventActive(aes_phy_tc(jg)%ev_rad,   datetime_old)
       !
       ! RTE-RRTMGP
       !
       CALL message_forcing_action('LW and SW radiation (rad:fluxes )' ,&
            &                      is_in_sd_ed_interval, is_active )
       !
       ! radiative fluxes
       CALL omp_loop_cell     (patch, interface_aes_rad        ,&
            &                  is_in_sd_ed_interval, is_active ,&
            &                  datetime_old, pdtime            )

       !
       ! always compute radiative heating
       is_active = .TRUE.
       !
       CALL message_forcing_action('LW and SW radiation (rht:heating)' ,&
            &                      is_in_sd_ed_interval, is_active )
       !
       ! radiative heating
       CALL omp_loop_cell(patch, interface_aes_rht        ,&
            &             is_in_sd_ed_interval, is_active ,&
            &             datetime_old, pdtime            )
       !
    END IF


    !-------------------------------------------------------------------
    ! Vertical diffusion, boundary layer and surface
    !-------------------------------------------------------------------
    !
    IF ( aes_phy_tc(jg)%dt_vdf > dt_zero ) THEN
       !
       is_in_sd_ed_interval =          (aes_phy_tc(jg)%sd_vdf <= datetime_old) .AND. &
            &                          (aes_phy_tc(jg)%ed_vdf >  datetime_old)
       is_active = isCurrentEventActive(aes_phy_tc(jg)%ev_vdf,   datetime_old)
       !
       CALL message_forcing_action('vertical diffusion (vdf)'      ,&
            &                      is_in_sd_ed_interval, is_active )
       !
       CALL interface_aes_vdf(patch ,&
            &                   is_in_sd_ed_interval, is_active ,&
            &                   datetime_old, pdtime            )
       !
    END IF

    !-------------------------------------------------------------------
    ! Linearized ozone chemistry of Cariolle
    !-------------------------------------------------------------------
    !
    IF ( aes_phy_tc(jg)%dt_car > dt_zero ) THEN
#if defined( _OPENACC )
       CALL warning('GPU:aes_car_main','GPU host synchronization should be removed when port is done!')
       CALL gpu_update_var_list('prm_field_D', .false., jg)
       CALL gpu_update_var_list('prm_tend_D', .false., jg)
#endif
       !
       is_in_sd_ed_interval =          (aes_phy_tc(jg)%sd_car <= datetime_old) .AND. &
            &                          (aes_phy_tc(jg)%ed_car >  datetime_old)
       is_active = isCurrentEventActive(aes_phy_tc(jg)%ev_car,   datetime_old)
       !
       CALL message_forcing_action('lin. Cariolle ozone chem. (car)' ,&
            &                      is_in_sd_ed_interval, is_active   )
       !
       CALL omp_loop_cell(patch, interface_aes_car        ,&
            &             is_in_sd_ed_interval, is_active ,&
            &             datetime_old, pdtime            )
       !
#if defined( _OPENACC )
       CALL warning('GPU:aes_car_main','GPU device synchronization should be removed when port is done!')
       CALL gpu_update_var_list('prm_field_D', .true., jg)
       CALL gpu_update_var_list('prm_tend_D', .true., jg)
#endif
    END IF

    !-------------------------------------------------------------------
    ! Graupel (microphysics) processes
    !-------------------------------------------------------------------
    !
    IF ( aes_phy_tc(jg)%dt_mig > dt_zero ) THEN
       !
       is_in_sd_ed_interval =          (aes_phy_tc(jg)%sd_mig <= datetime_old) .AND. &
            &                          (aes_phy_tc(jg)%ed_mig >  datetime_old)
       is_active = isCurrentEventActive(aes_phy_tc(jg)%ev_mig,   datetime_old)
       !
       CALL message_forcing_action('graupel microphysics (mig)',    &
            &                      is_in_sd_ed_interval, is_active)
       !
       SELECT CASE(aes_phy_config(jg)%if_mig)
       CASE(1,11)
          CALL omp_loop_cell  (patch, interface_cloud_mig_1    ,&
               &               is_in_sd_ed_interval, is_active ,&
               &               datetime_old, pdtime            )
       CASE(2,12)
          CALL omp_loop_cell  (patch, interface_cloud_mig_2    ,&
               &               is_in_sd_ed_interval, is_active ,&
               &               datetime_old, pdtime            )
       CASE(3,4,13,14)
          CALL omp_loop_cell_3(patch, interface_cloud_mig_3    ,&
               &               is_in_sd_ed_interval, is_active ,&
               &               datetime_old, pdtime            )
       END SELECT
       !
    END IF

    !--------------------------------------------------------------------
    ! two-moment bulk microphysics by Seifert and Beheng (2006) processes
    !--------------------------------------------------------------------
    !
    IF ( aes_phy_tc(jg)%dt_two > dt_zero ) THEN
       !
       is_in_sd_ed_interval =          (aes_phy_tc(jg)%sd_two <= datetime_old) .AND. &
            &                          (aes_phy_tc(jg)%ed_two >  datetime_old)
       is_active = isCurrentEventActive(aes_phy_tc(jg)%ev_two,   datetime_old)
       !
       CALL message_forcing_action('two-moment bulk microphysics (two)',    &
            &                      is_in_sd_ed_interval, is_active)
       !
       ! Preliminary: Some run time diagnostics (can also be used for other schemes)
       IF (msg_level>14) THEN
          CALL aes_diag_output_minmax_micro(patch,.TRUE.)
       END IF

       CALL omp_loop_cell(patch, interface_cloud_two      ,&
            &             is_in_sd_ed_interval, is_active ,&
            &             datetime_old, pdtime            )
       !
       ! Preliminary: Some run time diagnostics (can also be used for other schemes)
       IF (msg_level>14) THEN
          CALL aes_diag_output_minmax_micro(patch,.FALSE.)
       END IF

    END IF

    !-------------------------------------------------------------------
    ! Global output diagnostics
    !-------------------------------------------------------------------
    !
    CALL aes_global_diagnostics(patch)


    !-------------------------------------------------------------------
    ! Finalize (diagnostic)
    !-------------------------------------------------------------------
    !
    CALL omp_loop_cell(patch,finalize)

  END SUBROUTINE aes_phy_main
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  SUBROUTINE message_forcing_action(process, is_in_sd_ed_interval, is_active)
    CHARACTER(LEN=*) ,INTENT(in) :: process
    LOGICAL          ,INTENT(in) :: is_in_sd_ed_interval
    LOGICAL          ,INTENT(in) :: is_active

    IF (is_in_sd_ed_interval) THEN
       IF (is_active) THEN
          CALL message('aes_phy_main','compute forcing by '//process)
       ELSE
          CALL message('aes_phy_main','recycle forcing by '//process)
       END IF
    ELSE
       CALL    message('aes_phy_main','no      forcing by '//process)
    END IF

  END SUBROUTINE message_forcing_action
  !---------------------------------------------------------------------


END MODULE mo_aes_phy_main
