MODULE mo_nwp_gpu_util

  USE mo_ext_data_types,          ONLY: t_external_data
  USE mo_nwp_phy_types,           ONLY: t_nwp_phy_diag
  USE mo_model_domain,            ONLY: t_patch
  USE mo_dynamics_config,         ONLY: nnow, nnew, nnow_rcf, nnew_rcf
  USE mo_intp_data_strc,          ONLY: t_int_state
  USE mo_nwp_parameters,          ONLY: t_phy_params
  USE mo_nonhydrostatic_config,   ONLY: kstart_moist, kstart_tracer
  USE mo_run_config,              ONLY: ldass_lhn
  USE mo_atm_phy_nwp_config,      ONLY: t_atm_phy_nwp_config
  USE mo_fortran_tools,           ONLY: assert_acc_device_only
#ifdef _OPENACC
  USE mo_var_list_gpu,            ONLY: gpu_update_var_list
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: gpu_d2h_nh_nwp, gpu_h2d_nh_nwp, devcpy_nwp, hostcpy_nwp

  CONTAINS

  SUBROUTINE gpu_d2h_nh_nwp(pt_patch, prm_diag, ext_data, p_int, phy_params, atm_phy_nwp_config, lacc)

    TYPE(t_patch), TARGET, INTENT(in) :: pt_patch
    TYPE(t_nwp_phy_diag), INTENT(inout) :: prm_diag
    TYPE(t_external_data), OPTIONAL, INTENT(inout):: ext_data
    TYPE(t_int_state), OPTIONAL, INTENT(inout) :: p_int
    TYPE(t_phy_params), OPTIONAL, INTENT(inout) :: phy_params
    TYPE(t_atm_phy_nwp_config), OPTIONAL, TARGET, INTENT(inout) :: atm_phy_nwp_config
    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc

    TYPE(t_atm_phy_nwp_config), POINTER :: a

    INTEGER :: jg

    CALL assert_acc_device_only("gpu_d2h_nh_nwp", lacc)

    !$ACC UPDATE HOST(prm_diag%qrs_flux)

    !$ACC UPDATE HOST(ext_data%atm%gp_count_t, ext_data%atm%lp_count_t, ext_data%atm%list_seaice%ncount) &
    !$ACC   HOST(ext_data%atm%list_seaice%idx, ext_data%atm%list_lake%ncount, ext_data%atm%list_lake%idx) &
    !$ACC   HOST(ext_data%atm%list_land%ncount, ext_data%atm%list_land%idx) &
    !$ACC   HOST(ext_data%atm%list_seawtr%ncount, ext_data%atm%list_seawtr%idx, ext_data%atm%emis_rad) &
    !$ACC   HOST(ext_data%atm%z0_lcc, ext_data%atm%z0_lcc_min, ext_data%atm%plcovmax_lcc) &
    !$ACC   HOST(ext_data%atm%laimax_lcc, ext_data%atm%rootdmax_lcc, ext_data%atm%stomresmin_lcc) &
    !$ACC   HOST(ext_data%atm%snowalb_lcc, ext_data%atm%snowtile_lcc, ext_data%atm%t_cl) &
    !$ACC   IF(PRESENT(ext_data))

    !$ACC UPDATE HOST(p_int%lsq_high, p_int%lsq_lin) &
    !$ACC   HOST(p_int%c_bln_avg, p_int%c_lin_e, p_int%cells_aw_verts) &
    !$ACC   HOST(p_int%e_bln_c_s, p_int%e_flx_avg, p_int%geofac_div) &
    !$ACC   HOST(p_int%geofac_grdiv, p_int%geofac_grg, p_int%geofac_n2s) &
    !$ACC   HOST(p_int%geofac_rot, p_int%lsq_high%lsq_blk_c) &
    !$ACC   HOST(p_int%lsq_high%lsq_dim_stencil, p_int%lsq_high%lsq_idx_c) &
    !$ACC   HOST(p_int%lsq_high%lsq_moments, p_int%lsq_high%lsq_moments_hat) &
    !$ACC   HOST(p_int%lsq_high%lsq_pseudoinv, p_int%lsq_high%lsq_qtmat_c) &
    !$ACC   HOST(p_int%lsq_high%lsq_rmat_utri_c, p_int%lsq_high%lsq_weights_c) &
    !$ACC   HOST(p_int%lsq_lin%lsq_blk_c) &
    !$ACC   HOST(p_int%lsq_lin%lsq_dim_stencil, p_int%lsq_lin%lsq_idx_c) &
    !$ACC   HOST(p_int%lsq_lin%lsq_moments, p_int%lsq_lin%lsq_moments_hat) &
    !$ACC   HOST(p_int%lsq_lin%lsq_pseudoinv, p_int%lsq_lin%lsq_qtmat_c) &
    !$ACC   HOST(p_int%lsq_lin%lsq_rmat_utri_c, p_int%lsq_lin%lsq_weights_c) &
    !$ACC   HOST(p_int%nudgecoeff_c, p_int%nudgecoeff_e, p_int%pos_on_tplane_e) &
    !$ACC   HOST(p_int%rbf_c2grad_blk, p_int%rbf_c2grad_idx, p_int%rbf_c2grad_coeff) &
    !$ACC   HOST(p_int%rbf_vec_blk_c, p_int%rbf_vec_idx_c, p_int%rbf_vec_coeff_c) &
    !$ACC   HOST(p_int%rbf_vec_blk_e, p_int%rbf_vec_idx_e, p_int%rbf_vec_coeff_e) &
    !$ACC   HOST(p_int%rbf_vec_blk_v, p_int%rbf_vec_idx_v, p_int%rbf_vec_coeff_v) &
    !$ACC   HOST(p_int%verts_aw_cells) &
    !$ACC   IF(PRESENT(p_int))

#ifdef _OPENACC
    jg = pt_patch%id

    ! Update NWP fields
    CALL gpu_update_var_list('prm_diag_of_domain_', .false., domain=jg, lacc=.TRUE.)
    CALL gpu_update_var_list('prm_tend_of_domain_', .false., domain=jg, lacc=.TRUE.)
    CALL gpu_update_var_list('lnd_prog_of_domain_', .false., domain=jg, substr='_and_timelev_', timelev=nnow(jg), lacc=.TRUE.)
    CALL gpu_update_var_list('lnd_prog_of_domain_', .false., domain=jg, substr='_and_timelev_', timelev=nnew(jg), lacc=.TRUE.)
    CALL gpu_update_var_list('lnd_diag_of_domain_', .false., domain=jg, lacc=.TRUE.)
    CALL gpu_update_var_list('wtr_prog_of_domain_', .false., domain=jg, substr='_and_timelev_', timelev=nnow(jg), lacc=.TRUE.)
    CALL gpu_update_var_list('wtr_prog_of_domain_', .false., domain=jg, substr='_and_timelev_', timelev=nnew(jg), lacc=.TRUE.)
    CALL gpu_update_var_list('ext_data_atm_td_D',.false.,domain=jg, lacc=.TRUE.)
    CALL gpu_update_var_list('ext_data_atm_D', .false., domain=jg, lacc=.TRUE.)

    IF(ldass_lhn) THEN
        ! Update radar data fields
        CALL gpu_update_var_list('radar_data_ct_dom_', .false., domain=jg, lacc=.TRUE.)
        CALL gpu_update_var_list('radar_data_td_dom_', .false., domain=jg, lacc=.TRUE.)
    ENDIF

    ! Update dynamics fields
    CALL gpu_update_var_list('nh_state_metrics_of_domain_', .false., domain=jg, lacc=.TRUE.)
    CALL gpu_update_var_list('nh_state_diag_of_domain_', .false., domain=jg, lacc=.TRUE.)
    CALL gpu_update_var_list('nh_state_prog_of_domain_', .false., domain=jg, substr='_and_timelev_', timelev=nnow(jg), lacc=.TRUE.)
    CALL gpu_update_var_list('nh_state_prog_of_domain_', .false., domain=jg, substr='_and_timelev_', timelev=nnew(jg), lacc=.TRUE.) !p_prog
    CALL gpu_update_var_list('nh_state_prog_of_domain_', .false., domain=jg, substr='_and_timelev_', timelev=nnow_rcf(jg), lacc=.TRUE.) !p_prog_now_rcf
    CALL gpu_update_var_list('nh_state_prog_of_domain_', .false., domain=jg, substr='_and_timelev_', timelev=nnew_rcf(jg), lacc=.TRUE.) !p_prog_rcf
    CALL gpu_update_var_list('prepadv_of_domain_', .false., domain=jg, lacc=.TRUE.)
#endif

    IF (PRESENT(phy_params)) THEN
      ! This is save as long as all t_phy_params components are scalars.
      !$ACC UPDATE HOST(phy_params)
    ENDIF

    IF (PRESENT(atm_phy_nwp_config)) THEN
      a => atm_phy_nwp_config
      ! a%phyProc* are not allocated on GPU yet.
      !$ACC UPDATE HOST(a%inwp_gscp, a%inwp_satad, a%inwp_convection, a%lshallowconv_only, a%lgrayzone_deepconv) &
      !$ACC   HOST(a%ldetrain_conv_prec, a%inwp_radiation, a%inwp_sso, a%inwp_gwd, a%inwp_cldcover, a%inwp_turb) &
      !$ACC   HOST(a%inwp_surface, a%itype_z0, a%dt_conv, a%dt_ccov, a%dt_rad, a%dt_sso, a%dt_gwd, a%dt_fastphy) &
      !$ACC   HOST(a%mu_rain, a%mu_snow, a%rain_n0_factor, a%qi0, a%qc0, a%icpl_aero_gscp, a%ustart_raylfric) &
      !$ACC   HOST(a%efdt_min_raylfric, a%latm_above_top, a%icalc_reff, a%icpl_rad_reff, a%ithermo_water) &
      !$ACC   HOST(a%lupatmo_phy, a%lenabled, a%lcall_phy, a%lcalc_acc_avg) &
      !$ACC   HOST(a%lcalc_extra_avg, a%lhave_graupel, a%l2moment, a%lhydrom_read_from_fg, a%lhydrom_read_from_ana) &
#ifndef __NO_ICON_LES__
      !$ACC   HOST(a%is_les_phy) &
#endif
      !$ACC   HOST(a%nclass_gscp, a%l_3d_rad_fluxes, a%l_3d_turb_fluxes, a%fac_ozone, a%shapefunc_ozone) &
      !$ACC   HOST(a%ozone_maxinc)
    ENDIF

  END SUBROUTINE gpu_d2h_nh_nwp

  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------

  SUBROUTINE gpu_h2d_nh_nwp(pt_patch, prm_diag, ext_data, p_int, phy_params, atm_phy_nwp_config, lacc)

    TYPE(t_patch), TARGET, INTENT(in) :: pt_patch
    TYPE(t_nwp_phy_diag), INTENT(inout) :: prm_diag
    TYPE(t_external_data), OPTIONAL, INTENT(inout):: ext_data
    TYPE(t_int_state), OPTIONAL, INTENT(inout) :: p_int
    TYPE(t_phy_params), OPTIONAL, INTENT(inout) :: phy_params
    TYPE(t_atm_phy_nwp_config), OPTIONAL, TARGET, INTENT(inout) :: atm_phy_nwp_config
    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc

    TYPE(t_atm_phy_nwp_config), POINTER :: a

    INTEGER :: jg

    CALL assert_acc_device_only("gpu_d2h_nh_nwp", lacc)

    !$ACC UPDATE DEVICE(prm_diag%qrs_flux)

    !$ACC UPDATE DEVICE(ext_data%atm%gp_count_t, ext_data%atm%lp_count_t, ext_data%atm%list_seaice%ncount) &
    !$ACC   DEVICE(ext_data%atm%list_seaice%idx, ext_data%atm%list_lake%ncount, ext_data%atm%list_lake%idx) &
    !$ACC   DEVICE(ext_data%atm%list_land%ncount, ext_data%atm%list_land%idx) &
    !$ACC   DEVICE(ext_data%atm%list_seawtr%ncount, ext_data%atm%list_seawtr%idx, ext_data%atm%emis_rad) &
    !$ACC   DEVICE(ext_data%atm%z0_lcc, ext_data%atm%z0_lcc_min, ext_data%atm%plcovmax_lcc) &
    !$ACC   DEVICE(ext_data%atm%laimax_lcc, ext_data%atm%rootdmax_lcc, ext_data%atm%stomresmin_lcc) &
    !$ACC   DEVICE(ext_data%atm%snowalb_lcc, ext_data%atm%snowtile_lcc, ext_data%atm%t_cl) &
    !$ACC   IF(PRESENT(ext_data))

    !$ACC UPDATE DEVICE(p_int%lsq_high, p_int%lsq_lin) &
    !$ACC   DEVICE(p_int%c_bln_avg, p_int%c_lin_e, p_int%cells_aw_verts) &
    !$ACC   DEVICE(p_int%e_bln_c_s, p_int%e_flx_avg, p_int%geofac_div) &
    !$ACC   DEVICE(p_int%geofac_grdiv, p_int%geofac_grg, p_int%geofac_n2s) &
    !$ACC   DEVICE(p_int%geofac_rot, p_int%lsq_high%lsq_blk_c) &
    !$ACC   DEVICE(p_int%lsq_high%lsq_dim_stencil, p_int%lsq_high%lsq_idx_c) &
    !$ACC   DEVICE(p_int%lsq_high%lsq_moments, p_int%lsq_high%lsq_moments_hat) &
    !$ACC   DEVICE(p_int%lsq_high%lsq_pseudoinv, p_int%lsq_high%lsq_qtmat_c) &
    !$ACC   DEVICE(p_int%lsq_high%lsq_rmat_utri_c, p_int%lsq_high%lsq_weights_c) &
    !$ACC   DEVICE(p_int%lsq_lin%lsq_blk_c) &
    !$ACC   DEVICE(p_int%lsq_lin%lsq_dim_stencil, p_int%lsq_lin%lsq_idx_c) &
    !$ACC   DEVICE(p_int%lsq_lin%lsq_moments, p_int%lsq_lin%lsq_moments_hat) &
    !$ACC   DEVICE(p_int%lsq_lin%lsq_pseudoinv, p_int%lsq_lin%lsq_qtmat_c) &
    !$ACC   DEVICE(p_int%lsq_lin%lsq_rmat_utri_c, p_int%lsq_lin%lsq_weights_c) &
    !$ACC   DEVICE(p_int%nudgecoeff_c, p_int%nudgecoeff_e, p_int%pos_on_tplane_e) &
    !$ACC   DEVICE(p_int%rbf_c2grad_blk, p_int%rbf_c2grad_idx, p_int%rbf_c2grad_coeff) &
    !$ACC   DEVICE(p_int%rbf_vec_blk_c, p_int%rbf_vec_idx_c, p_int%rbf_vec_coeff_c) &
    !$ACC   DEVICE(p_int%rbf_vec_blk_e, p_int%rbf_vec_idx_e, p_int%rbf_vec_coeff_e) &
    !$ACC   DEVICE(p_int%rbf_vec_blk_v, p_int%rbf_vec_idx_v, p_int%rbf_vec_coeff_v) &
    !$ACC   DEVICE(p_int%verts_aw_cells) &
    !$ACC   IF(PRESENT(p_int))
    
#ifdef _OPENACC
    jg = pt_patch%id

    ! Update NWP fields
    CALL gpu_update_var_list('prm_diag_of_domain_', .true., domain=jg, lacc=.TRUE.)
    CALL gpu_update_var_list('prm_tend_of_domain_', .true., domain=jg, lacc=.TRUE.)
    CALL gpu_update_var_list('lnd_prog_of_domain_', .true., domain=jg, substr='_and_timelev_', timelev=nnow(jg), lacc=.TRUE.)
    CALL gpu_update_var_list('lnd_prog_of_domain_', .true., domain=jg, substr='_and_timelev_', timelev=nnew(jg), lacc=.TRUE.)
    CALL gpu_update_var_list('lnd_diag_of_domain_', .true., domain=jg, lacc=.TRUE.)
    CALL gpu_update_var_list('wtr_prog_of_domain_', .true., domain=jg, substr='_and_timelev_', timelev=nnow(jg), lacc=.TRUE.)
    CALL gpu_update_var_list('wtr_prog_of_domain_', .true., domain=jg, substr='_and_timelev_', timelev=nnew(jg), lacc=.TRUE.)
    CALL gpu_update_var_list('lnd_diag_of_domain_', .true., domain=jg, lacc=.TRUE.)
    CALL gpu_update_var_list('ext_data_atm_td_D',.true.,domain=jg, lacc=.TRUE.)
    CALL gpu_update_var_list('ext_data_atm_D', .true., domain=jg, lacc=.TRUE.)

    IF(ldass_lhn) THEN
        ! Update radar data fields
        CALL gpu_update_var_list('radar_data_ct_dom_', .true., domain=jg, lacc=.TRUE.)
        CALL gpu_update_var_list('radar_data_td_dom_', .true., domain=jg, lacc=.TRUE.)
    ENDIF

    ! Update dynamics fields
    CALL gpu_update_var_list('nh_state_metrics_of_domain_', .true., domain=jg, lacc=.TRUE.)
    CALL gpu_update_var_list('nh_state_diag_of_domain_', .true., domain=jg, lacc=.TRUE.)
    CALL gpu_update_var_list('nh_state_prog_of_domain_', .true., domain=jg, substr='_and_timelev_', timelev=nnow(jg), lacc=.TRUE.)
    CALL gpu_update_var_list('nh_state_prog_of_domain_', .true., domain=jg, substr='_and_timelev_', timelev=nnew(jg), lacc=.TRUE.) !p_prog
    CALL gpu_update_var_list('nh_state_prog_of_domain_', .true., domain=jg, substr='_and_timelev_', timelev=nnow_rcf(jg), lacc=.TRUE.) !p_prog_now_rcf
    CALL gpu_update_var_list('nh_state_prog_of_domain_', .true., domain=jg, substr='_and_timelev_', timelev=nnew_rcf(jg), lacc=.TRUE.) !p_prog_new_rcf
    CALL gpu_update_var_list('prepadv_of_domain_', .true., domain=jg, lacc=.TRUE.)
#endif

    IF (PRESENT(phy_params)) THEN
      ! This is save as long as all t_phy_params components are scalars.
      !$ACC UPDATE DEVICE(phy_params)
    END IF

    IF (PRESENT(atm_phy_nwp_config)) THEN
      a => atm_phy_nwp_config
      ! a%phyProc* are not allocated on GPU yet.
      !$ACC UPDATE DEVICE(a%inwp_gscp, a%inwp_satad, a%inwp_convection, a%lshallowconv_only, a%lgrayzone_deepconv) &
      !$ACC   DEVICE(a%ldetrain_conv_prec, a%inwp_radiation, a%inwp_sso, a%inwp_gwd, a%inwp_cldcover, a%inwp_turb) &
      !$ACC   DEVICE(a%inwp_surface, a%itype_z0, a%dt_conv, a%dt_ccov, a%dt_rad, a%dt_sso, a%dt_gwd, a%dt_fastphy) &
      !$ACC   DEVICE(a%mu_rain, a%mu_snow, a%rain_n0_factor, a%qi0, a%qc0, a%icpl_aero_gscp, a%ustart_raylfric) &
      !$ACC   DEVICE(a%efdt_min_raylfric, a%latm_above_top, a%icalc_reff, a%icpl_rad_reff, a%ithermo_water) &
      !$ACC   DEVICE(a%lupatmo_phy, a%lenabled, a%lcall_phy, a%lcalc_acc_avg) &
      !$ACC   DEVICE(a%lcalc_extra_avg, a%lhave_graupel, a%l2moment, a%lhydrom_read_from_fg, a%lhydrom_read_from_ana) &
#ifndef __NO_ICON_LES__
      !$ACC   DEVICE(a%is_les_phy) &
#endif
      !$ACC   DEVICE(a%nclass_gscp, a%l_3d_rad_fluxes, a%l_3d_turb_fluxes, a%fac_ozone, a%shapefunc_ozone) &
      !$ACC   DEVICE(a%ozone_maxinc)
    ENDIF

  END SUBROUTINE gpu_h2d_nh_nwp

  SUBROUTINE devcpy_nwp(lacc)
    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc

    CALL assert_acc_device_only("devcpy_nwp", lacc)

    !$ACC ENTER DATA COPYIN(kstart_moist, kstart_tracer)

  END SUBROUTINE devcpy_nwp

  SUBROUTINE hostcpy_nwp(lacc)
    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc

    CALL assert_acc_device_only("hostcpy_nwp", lacc)

    !$ACC EXIT DATA DELETE(kstart_moist, kstart_tracer)

  END SUBROUTINE hostcpy_nwp

END MODULE mo_nwp_gpu_util
