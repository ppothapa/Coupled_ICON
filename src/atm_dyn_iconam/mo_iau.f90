!>
!! Routines and functions for the Incremental Analysis Update (IAU)
!!
!! This module contains routines and functions for the Incremental Analysis 
!! Update (IAU), including the iterative IAU. The characteristics of the 
!! iterative IAU are described in the following using the example of a global 
!! NWP forecast run:
!!
!! IAU iteration
!!
!!                     input
!!                       /
!!                      /
!!                     /
!!          ........../
!!         /
!!        /
!!       /
!!      /
!!     /
!!  -90min               0min              90min
!! ---|------------------|------------------|------------->
!!    |//////////////////| - - - - - - - - - - - - - - - ->
!!                               free forecast (iteration = false)
!!    \________IAU_______/
!!                       |
!!                       /
!!                      /
!!                     /
!!          ........../
!!         /   reset
!!        /
!!       /
!!      /
!!     /
!!  -90min               0min              90min
!! ---|------------------|------------------|------------->
!!    |//////////////////|//////////////////| free forecast
!!
!!    \_________________IAU________________/
!!
!!
!!
!! @author Daniel Reinert, DWD
!!
!! @Literature:
!! Bloom, S. C., Takacs, L. L., da Silva, A. M., & Ledvina, D. (1996).
!! Data Assimilation Using Incremental Analysis Updates, Monthly Weather Review,
!! 124(6), 1256-1271
!! Polavarapu, S., Ren, S., Clayton, A. M., Sankey, D., & Rochon, Y. (2004).
!! On the Relationship between Incremental Analysis Updating and
!! Incremental Digital Filtering, Monthly Weather Review, 132(10), 2495-2502
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2022-12-28)
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
MODULE mo_iau

  USE mo_kind,                  ONLY: wp
  USE mo_impl_constants,        ONLY: REAL_T, SINGLE_T, INT_T, BOOL_T, TLEV_NNOW, TLEV_NNOW_RCF, &
    &                                 LSS_JSBACH, nclass_aero
  USE mo_exception,             ONLY: finish
  USE mo_math_constants,        ONLY: pi
  USE mo_initicon_types,        ONLY: t_saveinit_state, t_saved_field
  USE mo_model_domain,          ONLY: t_patch
  USE mo_nonhydro_types,        ONLY: t_nh_state, t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_nwp_phy_types,         ONLY: t_nwp_phy_diag, t_nwp_phy_tend, t_nwp_phy_stochconv
  USE mo_nwp_lnd_types,         ONLY: t_lnd_state, t_lnd_prog, t_lnd_diag, t_wtr_prog
  USE mo_ext_data_types,        ONLY: t_external_data
  USE mo_radar_data_types,      ONLY: t_lhn_diag
  USE mo_parallel_config,       ONLY: nproma
  USE mo_grid_config,           ONLY: n_dom
  USE mo_dynamics_config,       ONLY: nnow, nnow_rcf
  USE mo_atm_phy_nwp_config,    ONLY: atm_phy_nwp_config, iprog_aero
  USE mo_lnd_nwp_config,        ONLY: ntiles_total, ntiles_water, lmulti_snow, nlev_soil, nlev_snow, &
    &                                 l2lay_rho_snow, lsnowtile, lprog_albsi, itype_trvg, itype_snowevap
  USE mo_io_config,             ONLY: var_in_output
  USE mo_initicon_config,       ONLY: type_iau_wgt, is_iau_active, iau_wgt_dyn, iau_wgt_adv, &
    &                                 ltile_coldstart, qcana_mode, qiana_mode, qrsgana_mode
  USE mo_run_config,            ONLY: ntracer, iqv, iqc, iqi, iqr, iqs, iqg, iqh, &
    &                                 iqm_max, iqni, iqnc, iqnr, iqns, iqng, iqnh
  USE mo_advection_config,      ONLY: advection_config
  USE mo_upatmo_config,         ONLY: upatmo_config
  USE mo_nonhydrostatic_config, ONLY: kstart_moist
  USE mo_fortran_tools,         ONLY: init, copy, assert_acc_device_only, assert_acc_host_only
  USE mo_var,                   ONLY: t_var
  USE mo_var_groups,            ONLY: var_groups_dyn
  USE mo_var_list_register,     ONLY: t_vl_register_iter
  USE mo_hash_table,            ONLY: t_HashTable
  USE mo_var_metadata,          ONLY: get_var_timelevel, get_var_name
  USE mo_satad,                 ONLY: qsat_rho
  USE mo_nh_diagnose_pres_temp, ONLY: diag_pres, diag_temp

#if defined(__PGI) || defined(__FLANG)
  USE mo_util_texthash,         ONLY: t_char_workaround
#endif

#ifdef _OPENACC
  USE mo_mpi,                   ONLY: i_am_accel_node
#endif

  IMPLICIT NONE

  PRIVATE

  ! subroutines
  !
  PUBLIC :: save_initial_state, restore_initial_state
  PUBLIC :: compute_iau_wgt
  PUBLIC :: iau_update_tracer

  TYPE(t_saveinit_state), ALLOCATABLE  :: saveinit(:)

CONTAINS



  !----------------------------------------------------------------------------
  !>
  !! Saves the initial state of NWP applications for the IAU iteration mode.
  !!
  !! @par Revision History
  !! Initial release by Guenther Zaengl, DWD, (2016-06-17)
  !!
  SUBROUTINE save_initial_state(p_patch, p_nh, prm_diag, prm_stochconv, p_lnd, ext_data)

    TYPE(t_patch),             INTENT(IN) :: p_patch(:)
    TYPE(t_nh_state),          INTENT(IN) :: p_nh(:)
    TYPE(t_nwp_phy_diag),      INTENT(IN) :: prm_diag(:)
    TYPE(t_nwp_phy_stochconv), INTENT(IN) :: prm_stochconv(:)
    TYPE(t_lnd_state), TARGET, INTENT(IN) :: p_lnd(:)
    TYPE(t_external_data),     INTENT(IN) :: ext_data(:)

    INTEGER :: jg, ntl, ntw, nlev, nlevp1, nblks_c, nblks_e

    TYPE(t_lnd_prog), POINTER :: lnd_prog
    TYPE(t_lnd_diag), POINTER :: lnd_diag
    TYPE(t_wtr_prog), POINTER :: wtr_prog

#ifdef _OPENACC
    if (i_am_accel_node) CALL finish('save_initial_state', 'This should be called in CPU mode only.')
#endif

    ntl = ntiles_total
    ntw = ntiles_total+ntiles_water

    ALLOCATE(saveinit(n_dom))

    DO jg = 1, n_dom

      IF(.NOT. p_patch(jg)%ldom_active) CYCLE

      CALL saveinit(jg)%init

      nlev    = p_patch(jg)%nlev
      nlevp1  = p_patch(jg)%nlevp1
      nblks_c = p_patch(jg)%nblks_c
      nblks_e = p_patch(jg)%nblks_e

      lnd_prog => p_lnd(jg)%prog_lnd(nnow_rcf(jg))
      lnd_diag => p_lnd(jg)%diag_lnd
      wtr_prog => p_lnd(jg)%prog_wtr(nnow_rcf(jg))

      ALLOCATE (saveinit(jg)%fr_seaice(nproma,nblks_c), saveinit(jg)%t_ice(nproma,nblks_c),    &
                saveinit(jg)%h_ice(nproma,nblks_c),     saveinit(jg)%gz0(nproma,nblks_c),      &
                saveinit(jg)%t_mnw_lk(nproma,nblks_c),  saveinit(jg)%t_wml_lk(nproma,nblks_c), &
                saveinit(jg)%h_ml_lk(nproma,nblks_c),   saveinit(jg)%t_bot_lk(nproma,nblks_c), &
                saveinit(jg)%c_t_lk(nproma,nblks_c),    saveinit(jg)%t_b1_lk(nproma,nblks_c),  &
                saveinit(jg)%h_b1_lk(nproma,nblks_c) )

      ALLOCATE (saveinit(jg)%theta_v(nproma,nlev,nblks_c), &
                saveinit(jg)%rho(nproma,nlev,nblks_c),     &
                saveinit(jg)%exner(nproma,nlev,nblks_c),   &
                saveinit(jg)%w(nproma,nlevp1,nblks_c),     &
                saveinit(jg)%tke(nproma,nlevp1,nblks_c),   &
                saveinit(jg)%vn(nproma,nlev,nblks_e),      &
                saveinit(jg)%gz0_t(nproma,nblks_c,ntw),    &
                saveinit(jg)%t_g_t(nproma,nblks_c,ntw),    &
                saveinit(jg)%t_sk_t(nproma,nblks_c,ntw),   &
                saveinit(jg)%qv_s_t(nproma,nblks_c,ntw),   &
                saveinit(jg)%freshsnow_t(nproma,nblks_c,ntl), &
                saveinit(jg)%snowfrac_t(nproma,nblks_c,ntl), &
                saveinit(jg)%snowfrac_lc_t(nproma,nblks_c,ntl), &
                saveinit(jg)%w_snow_t(nproma,nblks_c,ntl), &
                saveinit(jg)%w_i_t(nproma,nblks_c,ntl),    &
                saveinit(jg)%h_snow_t(nproma,nblks_c,ntl), &
                saveinit(jg)%t_snow_t(nproma,nblks_c,ntl), &
                saveinit(jg)%rho_snow_t(nproma,nblks_c,ntl), &
                saveinit(jg)%snowtile_flag_t(nproma,nblks_c,ntl), &
                saveinit(jg)%idx_lst_t(nproma,nblks_c,ntl), &
                saveinit(jg)%frac_t(nproma,nblks_c,ntw),    &
                saveinit(jg)%gp_count_t(nblks_c,ntl)        )

      ALLOCATE (saveinit(jg)%tracer(nproma,nlev,nblks_c,ntracer),      &
                saveinit(jg)%w_so_t(nproma,nlev_soil,nblks_c,ntl),     &
                saveinit(jg)%w_so_ice_t(nproma,nlev_soil,nblks_c,ntl), &
                saveinit(jg)%t_so_t(nproma,nlev_soil+1,nblks_c,ntl)    )

      IF (lmulti_snow) THEN
        ALLOCATE (saveinit(jg)%t_snow_mult_t(nproma,nlev_snow+1,nblks_c,ntl), &
                  saveinit(jg)%rho_snow_mult_t(nproma,nlev_snow,nblks_c,ntl), &
                  saveinit(jg)%wtot_snow_t(nproma,nlev_snow,nblks_c,ntl),     &
                  saveinit(jg)%wliq_snow_t(nproma,nlev_snow,nblks_c,ntl),     &
                  saveinit(jg)%dzh_snow_t(nproma,nlev_snow,nblks_c,ntl)       )
      ELSE IF (l2lay_rho_snow) THEN
        ALLOCATE (saveinit(jg)%rho_snow_mult_t(nproma,nlev_snow,nblks_c,ntl))
      ENDIF

      IF (atm_phy_nwp_config(jg)%lstoch_deep) THEN
        ALLOCATE (saveinit(jg)%clnum_d(nproma,nblks_c), &
                  saveinit(jg)%clmf_d(nproma,nblks_c)   )
      ENDIF
      IF (atm_phy_nwp_config(jg)%lstoch_sde) THEN
        ALLOCATE (saveinit(jg)%clnum_a(nproma,nblks_c), &
                  saveinit(jg)%clmf_a(nproma,nblks_c),  &
                  saveinit(jg)%clnum_p(nproma,nblks_c), &
                  saveinit(jg)%clmf_p(nproma,nblks_c)   )
      ENDIF

      IF (iprog_aero >= 1)     ALLOCATE (saveinit(jg)%aerosol(nproma,nclass_aero,nblks_c))
      IF (lprog_albsi)         ALLOCATE (saveinit(jg)%alb_si(nproma,nblks_c))
      IF (itype_trvg == 3)     ALLOCATE (saveinit(jg)%plantevap_t(nproma,nblks_c,ntl))
      IF (itype_snowevap == 3) ALLOCATE (saveinit(jg)%hsnow_max(nproma,nblks_c),saveinit(jg)%h_snow(nproma,nblks_c),&
                                         saveinit(jg)%snow_age(nproma,nblks_c))

!$OMP PARALLEL
      CALL copy(lnd_diag%fr_seaice, saveinit(jg)%fr_seaice)
      CALL copy(wtr_prog%t_ice, saveinit(jg)%t_ice)
      CALL copy(wtr_prog%h_ice, saveinit(jg)%h_ice)
      CALL copy(prm_diag(jg)%gz0, saveinit(jg)%gz0)
      CALL copy(wtr_prog%t_mnw_lk, saveinit(jg)%t_mnw_lk)
      CALL copy(wtr_prog%t_wml_lk, saveinit(jg)%t_wml_lk)
      CALL copy(wtr_prog%h_ml_lk, saveinit(jg)%h_ml_lk)
      CALL copy(wtr_prog%t_bot_lk, saveinit(jg)%t_bot_lk)
      CALL copy(wtr_prog%c_t_lk, saveinit(jg)%c_t_lk)
      CALL copy(wtr_prog%t_b1_lk, saveinit(jg)%t_b1_lk)
      CALL copy(wtr_prog%h_b1_lk, saveinit(jg)%h_b1_lk)

      CALL copy(p_nh(jg)%prog(nnow(jg))%theta_v, saveinit(jg)%theta_v)
      CALL copy(p_nh(jg)%prog(nnow(jg))%rho, saveinit(jg)%rho)
      CALL copy(p_nh(jg)%prog(nnow(jg))%exner, saveinit(jg)%exner)
      CALL copy(p_nh(jg)%prog(nnow(jg))%w, saveinit(jg)%w)
      CALL copy(p_nh(jg)%prog(nnow_rcf(jg))%tke, saveinit(jg)%tke)
      CALL copy(p_nh(jg)%prog(nnow(jg))%vn, saveinit(jg)%vn)
      CALL copy(p_nh(jg)%prog(nnow_rcf(jg))%tracer, saveinit(jg)%tracer)

      CALL copy(prm_diag(jg)%gz0_t, saveinit(jg)%gz0_t)
      CALL copy(lnd_prog%t_g_t, saveinit(jg)%t_g_t)
      CALL copy(lnd_prog%t_sk_t, saveinit(jg)%t_sk_t)
      CALL copy(lnd_diag%qv_s_t, saveinit(jg)%qv_s_t)
      CALL copy(lnd_diag%freshsnow_t, saveinit(jg)%freshsnow_t)
      CALL copy(lnd_diag%snowfrac_t, saveinit(jg)%snowfrac_t)
      CALL copy(lnd_diag%snowfrac_lc_t, saveinit(jg)%snowfrac_lc_t)
      CALL copy(lnd_prog%w_snow_t, saveinit(jg)%w_snow_t)
      CALL copy(lnd_prog%w_i_t, saveinit(jg)%w_i_t)
      CALL copy(lnd_diag%h_snow_t, saveinit(jg)%h_snow_t)
      CALL copy(lnd_prog%t_snow_t, saveinit(jg)%t_snow_t)
      CALL copy(lnd_prog%rho_snow_t, saveinit(jg)%rho_snow_t)
      CALL copy(lnd_prog%w_so_t, saveinit(jg)%w_so_t)
      CALL copy(lnd_prog%w_so_ice_t, saveinit(jg)%w_so_ice_t)
      CALL copy(lnd_prog%t_so_t, saveinit(jg)%t_so_t)

      IF (ntiles_total > 1 .AND. lsnowtile .AND. .NOT. ltile_coldstart) THEN
        CALL copy(ext_data(jg)%atm%snowtile_flag_t, saveinit(jg)%snowtile_flag_t)
        CALL copy(ext_data(jg)%atm%idx_lst_t, saveinit(jg)%idx_lst_t)
        CALL copy(ext_data(jg)%atm%frac_t, saveinit(jg)%frac_t)
        CALL copy(ext_data(jg)%atm%gp_count_t, saveinit(jg)%gp_count_t)
      ENDIF

      IF (lmulti_snow) THEN
        CALL copy(lnd_prog%t_snow_mult_t, saveinit(jg)%t_snow_mult_t)
        CALL copy(lnd_prog%rho_snow_mult_t, saveinit(jg)%rho_snow_mult_t)
        CALL copy(lnd_prog%wtot_snow_t, saveinit(jg)%wtot_snow_t)
        CALL copy(lnd_prog%wliq_snow_t, saveinit(jg)%wliq_snow_t)
        CALL copy(lnd_prog%dzh_snow_t, saveinit(jg)%dzh_snow_t)
      ELSE IF (l2lay_rho_snow) THEN
        CALL copy(lnd_prog%rho_snow_mult_t, saveinit(jg)%rho_snow_mult_t)
      ENDIF

      IF (atm_phy_nwp_config(jg)%lstoch_deep) THEN
        CALL copy(prm_stochconv(jg)%clnum_d, saveinit(jg)%clnum_d)
        CALL copy(prm_stochconv(jg)%clmf_d,  saveinit(jg)%clmf_d)
      ENDIF
      IF (atm_phy_nwp_config(jg)%lstoch_sde) THEN
        CALL copy(prm_stochconv(jg)%clnum_a, saveinit(jg)%clnum_a)
        CALL copy(prm_stochconv(jg)%clmf_a,  saveinit(jg)%clmf_a)
        CALL copy(prm_stochconv(jg)%clnum_p, saveinit(jg)%clnum_p)
        CALL copy(prm_stochconv(jg)%clmf_p,  saveinit(jg)%clmf_p)
      ENDIF
      
      IF (iprog_aero >= 1)  CALL copy(prm_diag(jg)%aerosol, saveinit(jg)%aerosol)
      IF (lprog_albsi)      CALL copy(wtr_prog%alb_si, saveinit(jg)%alb_si)
      IF (itype_trvg == 3)  CALL copy(lnd_diag%plantevap_t, saveinit(jg)%plantevap_t)
      IF (itype_snowevap == 3) THEN
        CALL copy(lnd_diag%hsnow_max, saveinit(jg)%hsnow_max)
        CALL copy(lnd_diag%h_snow, saveinit(jg)%h_snow)
        CALL copy(lnd_diag%snow_age, saveinit(jg)%snow_age)
      ENDIF

!$OMP END PARALLEL

    ENDDO

    IF (ANY(atm_phy_nwp_config(:)%inwp_surface == LSS_JSBACH)) THEN
      CALL save_var_group_state('jsb_init_vars', p_patch(:))
    END IF

  END SUBROUTINE save_initial_state

  !----------------------------------------------------------------------------
  !>
  !! Restores the initial state of NWP applications for the IAU iteration mode.
  !!
  !! @par Revision History
  !! Initial release by Guenther Zaengl, DWD, (2016-06-17)
  !!
  SUBROUTINE restore_initial_state(p_patch, p_nh, prm_diag, prm_tend, prm_stochconv, p_lnd, ext_data, lhn_fields)

    TYPE(t_patch),             INTENT(IN)    :: p_patch(:)
    TYPE(t_nh_state),          INTENT(INOUT) :: p_nh(:)
    TYPE(t_nwp_phy_diag),      INTENT(INOUT) :: prm_diag(:)
    TYPE(t_nwp_phy_tend),      INTENT(INOUT) :: prm_tend(:)
    TYPE(t_nwp_phy_stochconv), INTENT(INOUT) :: prm_stochconv(:)
    TYPE(t_lnd_state), TARGET, INTENT(INOUT) :: p_lnd(:)
    TYPE(t_external_data),     INTENT(INOUT) :: ext_data(:)
    TYPE(t_lhn_diag),          INTENT(INOUT) :: lhn_fields(:)

    INTEGER :: jg, ic, je, jb

    TYPE(t_lnd_prog), POINTER :: lnd_prog
    TYPE(t_lnd_diag), POINTER :: lnd_diag
    TYPE(t_wtr_prog), POINTER :: wtr_prog

#ifdef _OPENACC
    ! make sure that accelerator is not used (this check is necessary as long as copy/init use i_am_accel_node)
    CALL assert_acc_host_only("restore_initial_state", i_am_accel_node)
#endif

    IF (ANY(atm_phy_nwp_config(:)%inwp_surface == LSS_JSBACH)) THEN
      CALL restore_var_group_state('jsb_init_vars', p_patch(:))
    END IF

    DO jg = 1, n_dom

      IF(.NOT. p_patch(jg)%ldom_active) CYCLE

      lnd_prog => p_lnd(jg)%prog_lnd(nnow_rcf(jg))
      lnd_diag => p_lnd(jg)%diag_lnd
      wtr_prog => p_lnd(jg)%prog_wtr(nnow_rcf(jg))

!$OMP PARALLEL
      CALL copy(saveinit(jg)%fr_seaice, lnd_diag%fr_seaice)
      CALL copy(saveinit(jg)%t_ice, wtr_prog%t_ice)
      CALL copy(saveinit(jg)%h_ice, wtr_prog%h_ice)
      CALL copy(saveinit(jg)%gz0, prm_diag(jg)%gz0)
      CALL copy(saveinit(jg)%t_mnw_lk, wtr_prog%t_mnw_lk)
      CALL copy(saveinit(jg)%t_wml_lk, wtr_prog%t_wml_lk)
      CALL copy(saveinit(jg)%h_ml_lk, wtr_prog%h_ml_lk)
      CALL copy(saveinit(jg)%t_bot_lk, wtr_prog%t_bot_lk)
      CALL copy(saveinit(jg)%c_t_lk, wtr_prog%c_t_lk)
      CALL copy(saveinit(jg)%t_b1_lk, wtr_prog%t_b1_lk)
      CALL copy(saveinit(jg)%h_b1_lk, wtr_prog%h_b1_lk)

      CALL copy(saveinit(jg)%theta_v, p_nh(jg)%prog(nnow(jg))%theta_v)
      CALL copy(saveinit(jg)%rho, p_nh(jg)%prog(nnow(jg))%rho)
      CALL copy(saveinit(jg)%exner, p_nh(jg)%prog(nnow(jg))%exner)
      CALL copy(saveinit(jg)%w, p_nh(jg)%prog(nnow(jg))%w)
      CALL copy(saveinit(jg)%tke, p_nh(jg)%prog(nnow_rcf(jg))%tke)
      CALL copy(saveinit(jg)%vn, p_nh(jg)%prog(nnow(jg))%vn)
      CALL copy(saveinit(jg)%tracer, p_nh(jg)%prog(nnow_rcf(jg))%tracer)

      CALL copy(saveinit(jg)%gz0_t, prm_diag(jg)%gz0_t)
      CALL copy(saveinit(jg)%t_g_t, lnd_prog%t_g_t)
      CALL copy(saveinit(jg)%t_sk_t, lnd_prog%t_sk_t)
      CALL copy(saveinit(jg)%qv_s_t, lnd_diag%qv_s_t)
      CALL copy(saveinit(jg)%freshsnow_t, lnd_diag%freshsnow_t)
      CALL copy(saveinit(jg)%snowfrac_t, lnd_diag%snowfrac_t)
      CALL copy(saveinit(jg)%snowfrac_lc_t, lnd_diag%snowfrac_lc_t)
      CALL copy(saveinit(jg)%w_snow_t, lnd_prog%w_snow_t)
      CALL copy(saveinit(jg)%w_i_t, lnd_prog%w_i_t)
      CALL copy(saveinit(jg)%h_snow_t, lnd_diag%h_snow_t)
      CALL copy(saveinit(jg)%t_snow_t, lnd_prog%t_snow_t)
      CALL copy(saveinit(jg)%rho_snow_t, lnd_prog%rho_snow_t)
      CALL copy(saveinit(jg)%w_so_t, lnd_prog%w_so_t)
      CALL copy(saveinit(jg)%w_so_ice_t, lnd_prog%w_so_ice_t)
      CALL copy(saveinit(jg)%t_so_t, lnd_prog%t_so_t)

      IF (ntiles_total > 1 .AND. lsnowtile .AND. .NOT. ltile_coldstart) THEN
        CALL copy(saveinit(jg)%snowtile_flag_t, ext_data(jg)%atm%snowtile_flag_t)
        CALL copy(saveinit(jg)%idx_lst_t, ext_data(jg)%atm%idx_lst_t)
        CALL copy(saveinit(jg)%frac_t, ext_data(jg)%atm%frac_t)
        CALL copy(saveinit(jg)%gp_count_t, ext_data(jg)%atm%gp_count_t)
      ENDIF

      IF (lmulti_snow) THEN
        CALL copy(saveinit(jg)%t_snow_mult_t, lnd_prog%t_snow_mult_t)
        CALL copy(saveinit(jg)%rho_snow_mult_t, lnd_prog%rho_snow_mult_t)
        CALL copy(saveinit(jg)%wtot_snow_t, lnd_prog%wtot_snow_t)
        CALL copy(saveinit(jg)%wliq_snow_t, lnd_prog%wliq_snow_t)
        CALL copy(saveinit(jg)%dzh_snow_t, lnd_prog%dzh_snow_t)
      ELSE IF (l2lay_rho_snow) THEN
        CALL copy(saveinit(jg)%rho_snow_mult_t, lnd_prog%rho_snow_mult_t)
      ENDIF

      IF (atm_phy_nwp_config(jg)%lstoch_deep) THEN
        CALL copy(saveinit(jg)%clnum_d, prm_stochconv(jg)%clnum_d)
        CALL copy(saveinit(jg)%clmf_d,  prm_stochconv(jg)%clmf_d)
      ENDIF
      IF (atm_phy_nwp_config(jg)%lstoch_sde) THEN
        CALL copy(saveinit(jg)%clnum_a, prm_stochconv(jg)%clnum_a)
        CALL copy(saveinit(jg)%clmf_a,  prm_stochconv(jg)%clmf_a)
        CALL copy(saveinit(jg)%clnum_p, prm_stochconv(jg)%clnum_p)
        CALL copy(saveinit(jg)%clmf_p,  prm_stochconv(jg)%clmf_p)
      ENDIF

      IF (iprog_aero >= 1) THEN
        CALL copy(saveinit(jg)%aerosol, prm_diag(jg)%aerosol)
      ENDIF
      IF (lprog_albsi) THEN
        CALL copy(saveinit(jg)%alb_si, wtr_prog%alb_si)
      ENDIF
      IF (itype_trvg == 3) THEN
        CALL copy(saveinit(jg)%plantevap_t, lnd_diag%plantevap_t)
      ENDIF
      IF (itype_snowevap == 3) THEN
        CALL copy(saveinit(jg)%hsnow_max, lnd_diag%hsnow_max)
        CALL copy(saveinit(jg)%h_snow, lnd_diag%h_snow)
        CALL copy(saveinit(jg)%snow_age, lnd_diag%snow_age)
      ENDIF

      ! Fields that need to be reset to zero in order to obtain identical results
      CALL init (p_nh(jg)%diag%ddt_vn_phy)
      CALL init (p_nh(jg)%diag%ddt_tracer_adv)
      CALL init (prm_tend(jg)%ddt_tracer_turb)
      CALL init (prm_tend(jg)%ddt_temp_radsw)
      CALL init (prm_tend(jg)%ddt_temp_radlw)
      CALL init (prm_tend(jg)%ddt_temp_turb)
      CALL init (p_nh(jg)%diag%ddt_temp_dyn)
      CALL init (p_nh(jg)%diag%exner_dyn_incr)
      CALL init (prm_diag(jg)%rain_gsp_rate)
      CALL init (prm_diag(jg)%snow_gsp_rate)
      IF (ANY((/1,2,4,5,6,7/) == atm_phy_nwp_config(jg)%inwp_gscp)) THEN
        CALL init (prm_diag(jg)%ice_gsp_rate)
      ENDIF
      CALL init (prm_diag(jg)%shfl_s_t)
      CALL init (prm_diag(jg)%qhfl_s_t)
      CALL init (lnd_diag%runoff_s_t)
      CALL init (lnd_diag%runoff_g_t)
      CALL init (lnd_diag%runoff_s_inst_t)
      CALL init (lnd_diag%runoff_g_inst_t)
      IF (var_in_output(jg)%res_soilwatb) THEN
        CALL init (lnd_diag%resid_wso_t)
        CALL init (lnd_diag%resid_wso_inst_t)
      ENDIF

      ! If explicit stochastic scheme is selected, set all cloud ensemble
      ! variables to zero at init stage. Required to get identical results
      ! for IAU/iterative IAU.
      IF (atm_phy_nwp_config(jg)%lstoch_expl) THEN
        CALL init (prm_stochconv(jg)%mf_i)
        CALL init (prm_stochconv(jg)%time_i)
        CALL init (prm_stochconv(jg)%life_i)
        CALL init (prm_stochconv(jg)%area_i)
        CALL init (prm_stochconv(jg)%type_i)
        CALL init (prm_stochconv(jg)%ktype_i)
        CALL init (prm_stochconv(jg)%depth_i)
        CALL init (prm_stochconv(jg)%base_i)
        CALL init (prm_stochconv(jg)%used_cell)
      ENDIF

!$OMP END PARALLEL


      ! For the limited-area mode and one-way nesting, we also need to reset grf_tend_vn on the nudging points
      !
      DO ic = 1, p_nh(jg)%metrics%nudge_e_dim
        je = p_nh(jg)%metrics%nudge_e_idx(ic)
        jb = p_nh(jg)%metrics%nudge_e_blk(ic)
        p_nh(jg)%diag%grf_tend_vn(je,:,jb) = 0._wp
      ENDDO

      ! deallocate
      CALL saveinit(jg)%finalize

    ENDDO

    DEALLOCATE(saveinit)

  END SUBROUTINE restore_initial_state

  !>
  !! Helper function to extract a value from the hash table.
  FUNCTION get_val (tbl, key) RESULT(val)
    TYPE(t_HashTable), INTENT(IN) :: tbl !< Hash table.
    CHARACTER(*), TARGET, INTENT(IN) :: key !< Key string.
    CLASS(*), POINTER :: val !< Value associated with `key` or NULL

    CLASS(*), POINTER :: ptr

    ptr => key
    val => tbl%getEntry(ptr)
  END FUNCTION

  !>
  !! Save state of a variable group for iterative IAU.
  SUBROUTINE save_var_group_state (group_name, p_patch)
    CHARACTER(len=*), INTENT(IN) :: group_name !< Name of the group to save.
    TYPE(t_patch), INTENT(IN) :: p_patch(:) !< Domain patches.

    CHARACTER(*), PARAMETER :: routine = 'save_var_group_state'

    TYPE(t_vl_register_iter) :: iter
    CLASS(*), POINTER :: key
    CLASS(*), POINTER :: val
    TYPE(t_saved_field), POINTER :: field
    INTEGER :: group_id

#if defined(__PGI) || defined(__FLANG)
    TYPE(t_char_workaround), POINTER :: key_p
#endif

    INTEGER :: i, jg, tl

#ifdef _OPENACC
    CALL finish('save_var_group_state', 'not tested with OpenACC')
#endif

    group_id = var_groups_dyn%group_id(TRIM(group_name))

    DO WHILE (iter%next())
      jg = iter%cur%p%patch_id
      ASSOCIATE (vars => iter%cur%p%vl, nvars => iter%cur%p%nvars, fields => saveinit(jg)%fields)

        IF (.NOT. p_patch(jg)%ldom_active) CYCLE

        DO i = 1, nvars
          IF (.NOT. vars(i)%p%info%in_group(group_id)) CYCLE

          tl = get_var_timelevel(vars(i)%p%info%name)

          IF (tl > 0) THEN
            IF (vars(i)%p%info%tlev_source == TLEV_NNOW .AND. tl /= nnow(jg)) CYCLE
            IF (vars(i)%p%info%tlev_source == TLEV_NNOW_RCF .AND. tl /= nnow_rcf(jg)) CYCLE
          END IF

          IF (ASSOCIATED(get_val(fields, TRIM(get_var_name(vars(i)%p%info))))) THEN
            CALL finish(routine, 'Variable ' // TRIM(vars(i)%p%info%name) // ' already saved!')
          END IF

#if defined(__PGI) || defined(__FLANG)
          ALLOCATE(key_p)
          key_p%c = TRIM(get_var_name(vars(i)%p%info))
          key => key_p
#else
          ALLOCATE(key, SOURCE=TRIM(get_var_name(vars(i)%p%info)))
#endif
          ALLOCATE(field)

          SELECT CASE (vars(i)%p%info%data_type)
          CASE (REAL_T)
            CALL field%put(vars(i)%p%r_ptr)
          CASE (SINGLE_T)
            CALL field%put(vars(i)%p%s_ptr)
          CASE (INT_T)
            CALL field%put(vars(i)%p%i_ptr)
          CASE (BOOL_T)
            CALL field%put(vars(i)%p%l_ptr)
          END SELECT

          val => field
          CALL fields%setEntry(key, val)
        END DO
      END ASSOCIATE
    END DO

  END SUBROUTINE save_var_group_state

  !>
  !! Restore state of a variable group for iterative IAU.
  SUBROUTINE restore_var_group_state (group_name, p_patch)
    CHARACTER(len=*), INTENT(IN) :: group_name !< Name of the group to save.
    TYPE(t_patch), INTENT(IN) :: p_patch(:) !< Domain patches.

    CHARACTER(*), PARAMETER :: routine = 'restore_var_group_state'

    TYPE(t_vl_register_iter) :: iter
    INTEGER :: group_id

    INTEGER :: i, jg

#ifdef _OPENACC
    CALL finish('restore_var_group_state', 'not tested with OpenACC')
#endif

    group_id = var_groups_dyn%group_id(TRIM(group_name))

    DO WHILE (iter%next())
      jg = iter%cur%p%patch_id
      IF (.NOT. p_patch(jg)%ldom_active) CYCLE

      ! The selection of variables is carried out by all threads redundantly. The actual work-
      ! sharing takes place inside `field%get`, which calls `copy`. The call to the nested
      ! subroutine is a workaround to some compilers (NEC) refusing CLASS(*) pointers in private
      ! clauses.
      !$OMP PARALLEL PRIVATE(i)
      DO i = 1, iter%cur%p%nvars
        CALL restore(saveinit(jg)%fields, iter%cur%p%vl(i)%p)
      END DO
      !$OMP END PARALLEL
    END DO

  CONTAINS

    SUBROUTINE restore(fields, var)
      TYPE(t_HashTable), INTENT(IN) :: fields
      TYPE(t_var), INTENT(INOUT) :: var

      CLASS(*), POINTER :: val
      TYPE(t_saved_field), POINTER :: field
      INTEGER :: tl

      IF (.NOT. var%info%in_group(group_id)) RETURN

      tl = get_var_timelevel(var%info%name)

      IF (tl > 0) THEN
        IF (var%info%tlev_source == TLEV_NNOW .AND. tl /= nnow(jg)) RETURN
        IF (var%info%tlev_source == TLEV_NNOW_RCF .AND. tl /= nnow_rcf(jg)) RETURN
      END IF

      field => NULL()
      val => get_val(fields, TRIM(get_var_name(var%info)))

      IF (.NOT. ASSOCIATED(val)) THEN
        CALL finish(routine, 'Variable ' // TRIM(var%info%name) // ' has not been saved!')
      END IF

      SELECT TYPE (val)
      TYPE IS (t_saved_field)
        field => val
      END SELECT

      SELECT CASE (var%info%data_type)
      CASE (REAL_T)
        CALL field%get(var%r_ptr)
      CASE (SINGLE_T)
        CALL field%get(var%s_ptr)
      CASE (INT_T)
        CALL field%get(var%i_ptr)
      CASE (BOOL_T)
        CALL field%get(var%l_ptr)
      END SELECT

    END SUBROUTINE restore

  END SUBROUTINE restore_var_group_state


  !>
  !! Compute weights for incremental analysis update
  !!
  !! Compute weights for incremental analysis update.
  !! 2 weights are provided:
  !! - iau_wgt_dyn can be used for all fields that need to be updated
  !!   every (fast) dynamics time step
  !! - iau_wgt_adv can be used for all fields that need to be updated
  !!   every (slow) advection time step.
  !!
  !! @par Revision History
  !! Initial revision by daniel Reinert, DWD (2014-01-29)
  !!
  SUBROUTINE compute_iau_wgt(sim_time, dt, dt_iau, lreset_wgt_adv)

    REAL(wp)        , INTENT(IN)  :: sim_time          !< Simulation time since model
                                                       !< start
    REAL(wp)        , INTENT(IN)  :: dt                !< time step
    REAL(wp)        , INTENT(IN)  :: dt_iau            !< width of IAU window
    LOGICAL         , INTENT(IN)  :: lreset_wgt_adv    !< If true, reset the accumulated weight for the advective time step

    ! local variables
    REAL(wp)  :: time_iau_elapsed                      !< elapsed time since IAU start [s]
    REAL(wp)  :: fct_eval                              !< result of top-hat or sin2 function

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_nh_init_utils:compute_iau_wgt'
    !-------------------------------------------------------------------------

    ! initialize
    IF (lreset_wgt_adv) iau_wgt_adv   = 0._wp

    ! compute elapsed time (in s) since IAU start
    !
    ! trivial so far, however will be changed to mtime when the functionality of
    ! computing the timedelta between two dates becomes available.
    time_iau_elapsed = sim_time


    IF (time_iau_elapsed <= dt_iau) THEN
      is_iau_active = .TRUE.

      SELECT CASE (type_iau_wgt)
        CASE(1)  ! top-hat function
          fct_eval = iau_top_hat(dt_iau,time_iau_elapsed)

        CASE(2)  ! sin2 function
          fct_eval = iau_sin2   (dt_iau,time_iau_elapsed)

        CASE(3)  ! sin function
          fct_eval = iau_sin    (dt_iau,time_iau_elapsed)

        CASE default
          CALL finish(routine,&
                      'Invalid IAU weighting function. Must be 1, 2 or 3.')
      END SELECT

      ! compute weights by multiplying with the time step
      iau_wgt_dyn = fct_eval * dt
      iau_wgt_adv = iau_wgt_adv + iau_wgt_dyn

    ELSE
      is_iau_active = .FALSE.
      iau_wgt_dyn   = 0._wp
      iau_wgt_adv   = 0._wp
    ENDIF

!!$write(0,*) "sim_time, is_iau_active, iau_wgt_dyn, iau_wgt_adv: ", &
!!$  & sim_time, is_iau_active, iau_wgt_dyn, iau_wgt_adv

  END SUBROUTINE compute_iau_wgt


  !>
  !! Evaluates top-hat function at a particular point in time
  !!
  !! Evaluates top-hat function at a particular point in time
  !! Top-hat function is non-zero for 0<=t<=dt and is normalized such that
  !! \int_{t=0}^{t=dt} f(t)\,dt=1
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2014-01-29)
  !!
  FUNCTION iau_top_hat (dt,cur_time)  RESULT (fct_eval)

    REAL(wp), INTENT(IN) :: dt                 ! time interval [s]
    REAL(wp), INTENT(in) :: cur_time           ! current time  [s]

    REAL(wp) :: fct_eval
    !-------------------------------------------------------------------------

    IF (cur_time <= dt) THEN
      fct_eval = 1._wp/dt
    ELSE
      fct_eval = 0._wp
    ENDIF

  END FUNCTION iau_top_hat


  !>
  !! Evaluates SIN2 function at a particular point in time
  !!
  !! Evaluates SIN2 function at a particular point in time
  !! SIN2 function is non-zero for 0<=t<=dt and is normalized such that
  !! \int_{t=0}^{t=dt} f(t)\,dt=1
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2014-01-29)
  !!
  FUNCTION iau_sin2 (dt,cur_time)  RESULT (fct_eval)

    REAL(wp), INTENT(IN) :: dt                 ! time interval [s]
    REAL(wp), INTENT(in) :: cur_time           ! current time  [s]

    REAL(wp) :: fct_eval
    !-------------------------------------------------------------------------

    IF (cur_time <= dt) THEN
      fct_eval = (2._wp/dt) * SIN(pi*cur_time/dt)**2
    ELSE
      fct_eval = 0._wp
    ENDIF

  END FUNCTION iau_sin2

  !>
  !! Evaluates SIN function at a particular point in time
  !!
  !! Evaluates SIN function at a particular point in time
  !! SIN function is non-zero for 0<=t<=dt and is normalized such that
  !! \int_{t=0}^{t=dt} f(t)\,dt=1
  !!
  !! @par Revision History
  !! Initial revision by Harald Anlauf, DWD (2014-04-03)
  !!
  FUNCTION iau_sin (dt, cur_time)  RESULT (fct_eval)

    REAL(wp), INTENT(IN) :: dt                 ! time interval [s]
    REAL(wp), INTENT(in) :: cur_time           ! current time  [s]

    REAL(wp) :: fct_eval
    !-------------------------------------------------------------------------

    IF (cur_time <= dt) THEN
      fct_eval = ((PI/2._wp)/dt) * SIN(PI*cur_time/dt)
    ELSE
      fct_eval = 0._wp
    ENDIF

  END FUNCTION iau_sin


  !
  ! Add IAU increment to qv during IAU phase
  !
  ! Add analysis increments from data assimilation to qv
  !
  ! Initial revision by Daniel Reinert, DWD (2018-05-18)
  ! Previously, this code snippet was part of nh_update_tracer_phy
  ! 
  SUBROUTINE iau_update_tracer( pt_prog, p_metrics, pt_diag, pt_prog_rcf, &
    &                     jg, jb, i_startidx, i_endidx, kend, lacc )

    TYPE(t_nh_prog)    ,INTENT(IN)   :: pt_prog      !< NH prog state at dynamic time step
    TYPE(t_nh_metrics) ,INTENT(IN)   :: p_metrics    !< NH metrics variables
    TYPE(t_nh_diag)    ,INTENT(INOUT):: pt_diag      !< the diagnostic variables
    TYPE(t_nh_prog)    ,INTENT(INOUT):: pt_prog_rcf  !< the tracer field at
                                                      !< reduced calling frequency
    INTEGER            ,INTENT(IN)   :: jg           !< domain ID
    INTEGER            ,INTENT(IN)   :: jb           !< block index
    INTEGER            ,INTENT(IN)   :: i_startidx   !< hor. start idx
    INTEGER            ,INTENT(IN)   :: i_endidx     !< hor. end idx
    INTEGER            ,INTENT(IN)   :: kend         !< vert. end idx
    LOGICAL, OPTIONAL  ,INTENT(IN)   :: lacc         ! If true, use openacc

    ! Local variables
    INTEGER  :: jk,jc
    REAL(wp) :: zqin
    REAL(wp) :: zrhw(nproma, kend) ! relative humidity w.r.t. water


    CALL assert_acc_device_only("iau_update_tracer", lacc)

    ! add analysis increments from data assimilation to qv
    !
    ! Diagnose pressure and temperature for subsequent calculations
    CALL diag_temp (pt_prog, pt_prog_rcf, advection_config(jg)%trHydroMass%list, pt_diag, &
                    jb, i_startidx, i_endidx, 1, kstart_moist(jg), kend)
    CALL diag_pres (pt_prog, pt_diag, p_metrics, jb, i_startidx, i_endidx, 1, kend, &
      &             opt_lconstgrav=upatmo_config(jg)%nwp_phy%l_constgrav)

    !$ACC DATA CREATE(zrhw) PRESENT(pt_prog, p_metrics, pt_diag, pt_prog_rcf, atm_phy_nwp_config)

    ! Compute relative humidity w.r.t. water
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP SEQ
    DO jk = 1, kend
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jc = i_startidx, i_endidx
        zrhw(jc,jk) = pt_prog_rcf%tracer(jc,jk,jb,iqv)/qsat_rho(pt_diag%temp(jc,jk,jb),pt_prog%rho(jc,jk,jb))
      ENDDO
    ENDDO

    ! GZ: This loop needs to be split for correct vectorization because rhoc_incr is allocated for qcana_mode >= 1 only;
    !     otherwise, the NEC runs into a segfault. Likewise, the remaining case selections need to be done outside the
    !     vectorized loops in order to avoid invalid memory accesses.
    !$ACC LOOP SEQ
    DO jk = 1, kend
      IF (qcana_mode >= 1) THEN
          !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zqin)
          DO jc = i_startidx, i_endidx
          IF (qcana_mode == 2 .AND. pt_prog_rcf%tracer(jc,jk,jb,iqc) > 0._wp) THEN
            pt_prog_rcf%tracer(jc,jk,jb,iqv) = pt_prog_rcf%tracer(jc,jk,jb,iqv) + &
              iau_wgt_adv*pt_diag%rhov_incr(jc,jk,jb)/pt_prog%rho(jc,jk,jb)
            pt_prog_rcf%tracer(jc,jk,jb,iqc) = MAX(0._wp,pt_prog_rcf%tracer(jc,jk,jb,iqc) + &
              iau_wgt_adv*pt_diag%rhoc_incr(jc,jk,jb)/pt_prog%rho(jc,jk,jb))
          ELSE 
            zqin = (pt_diag%rhov_incr(jc,jk,jb)+pt_diag%rhoc_incr(jc,jk,jb))/pt_prog%rho(jc,jk,jb)
            ! DA increments of humidity are limited to positive values if p > 150 hPa and RH < 2% or QV < 5.e-7
            IF (pt_diag%pres(jc,jk,jb) > 15000._wp .AND. zrhw(jc,jk) < 0.02_wp .OR. &
              pt_prog_rcf%tracer(jc,jk,jb,iqv) < 5.e-7_wp) zqin = MAX(0._wp, zqin)
            pt_prog_rcf%tracer(jc,jk,jb,iqv) = MAX(0._wp,pt_prog_rcf%tracer(jc,jk,jb,iqv) + iau_wgt_adv*zqin)
          ENDIF
        ENDDO
      ELSE
          !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zqin)
          DO jc = i_startidx, i_endidx
          zqin = pt_diag%rhov_incr(jc,jk,jb)/pt_prog%rho(jc,jk,jb)
          ! DA increments of humidity are limited to positive values if p > 150 hPa and RH < 2% or QV < 5.e-7
          IF (pt_diag%pres(jc,jk,jb) > 15000._wp .AND. zrhw(jc,jk) < 0.02_wp .OR. &
            pt_prog_rcf%tracer(jc,jk,jb,iqv) < 5.e-7_wp) zqin = MAX(0._wp, zqin)
          pt_prog_rcf%tracer(jc,jk,jb,iqv) = MAX(0._wp,pt_prog_rcf%tracer(jc,jk,jb,iqv) + iau_wgt_adv*zqin)
        ENDDO
      ENDIF

      IF (qiana_mode > 0) THEN
          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO jc = i_startidx, i_endidx
          pt_prog_rcf%tracer(jc,jk,jb,iqi) = MAX(0._wp,pt_prog_rcf%tracer(jc,jk,jb,iqi) + &
            iau_wgt_adv*pt_diag%rhoi_incr(jc,jk,jb)/pt_prog%rho(jc,jk,jb))
        ENDDO
      ENDIF

      IF (qrsgana_mode > 0) THEN
          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO jc = i_startidx, i_endidx
          pt_prog_rcf%tracer(jc,jk,jb,iqr) = MAX(0._wp,pt_prog_rcf%tracer(jc,jk,jb,iqr) + &
            iau_wgt_adv * pt_diag%rhor_incr(jc,jk,jb)/pt_prog%rho(jc,jk,jb))
          pt_prog_rcf%tracer(jc,jk,jb,iqs) = MAX(0._wp,pt_prog_rcf%tracer(jc,jk,jb,iqs) + &
            iau_wgt_adv * pt_diag%rhos_incr(jc,jk,jb)/pt_prog%rho(jc,jk,jb))
        ENDDO
      ENDIF

      IF (qrsgana_mode > 0 .AND. iqg <= iqm_max) THEN
          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO jc = i_startidx, i_endidx
          pt_prog_rcf%tracer(jc,jk,jb,iqg) = MAX(0._wp,pt_prog_rcf%tracer(jc,jk,jb,iqg) + &
            iau_wgt_adv * pt_diag%rhog_incr(jc,jk,jb)/pt_prog%rho(jc,jk,jb))
        ENDDO
      ENDIF

      IF (atm_phy_nwp_config(jg)%l2moment) THEN
          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO jc = i_startidx, i_endidx
          IF (qcana_mode > 0) THEN
            pt_prog_rcf%tracer(jc,jk,jb,iqnc) = MAX(0._wp,pt_prog_rcf%tracer(jc,jk,jb,iqnc) + &
                 iau_wgt_adv * pt_diag%rhonc_incr(jc,jk,jb)/pt_prog%rho(jc,jk,jb))
          END IF
          IF (qiana_mode > 0) THEN
            pt_prog_rcf%tracer(jc,jk,jb,iqni) = MAX(0._wp,pt_prog_rcf%tracer(jc,jk,jb,iqni) + &
                 iau_wgt_adv * pt_diag%rhoni_incr(jc,jk,jb)/pt_prog%rho(jc,jk,jb))
          END IF
          IF (qrsgana_mode > 0) THEN
            pt_prog_rcf%tracer(jc,jk,jb,iqh) = MAX(0._wp,pt_prog_rcf%tracer(jc,jk,jb,iqh) + &
                 iau_wgt_adv * pt_diag%rhoh_incr(jc,jk,jb)/pt_prog%rho(jc,jk,jb))
            pt_prog_rcf%tracer(jc,jk,jb,iqnr) = MAX(0._wp,pt_prog_rcf%tracer(jc,jk,jb,iqnr) + &
                 iau_wgt_adv * pt_diag%rhonr_incr(jc,jk,jb)/pt_prog%rho(jc,jk,jb))
            pt_prog_rcf%tracer(jc,jk,jb,iqns) = MAX(0._wp,pt_prog_rcf%tracer(jc,jk,jb,iqns) + &
                 iau_wgt_adv * pt_diag%rhons_incr(jc,jk,jb)/pt_prog%rho(jc,jk,jb))
            pt_prog_rcf%tracer(jc,jk,jb,iqng) = MAX(0._wp,pt_prog_rcf%tracer(jc,jk,jb,iqng) + &
                 iau_wgt_adv * pt_diag%rhong_incr(jc,jk,jb)/pt_prog%rho(jc,jk,jb))
            pt_prog_rcf%tracer(jc,jk,jb,iqnh) = MAX(0._wp,pt_prog_rcf%tracer(jc,jk,jb,iqnh) + &
                 iau_wgt_adv * pt_diag%rhonh_incr(jc,jk,jb)/pt_prog%rho(jc,jk,jb))
          END IF
        ENDDO
      ENDIF

    ENDDO
    !$ACC END PARALLEL

    !$ACC WAIT
    !$ACC END DATA

  END SUBROUTINE iau_update_tracer


END MODULE mo_iau
