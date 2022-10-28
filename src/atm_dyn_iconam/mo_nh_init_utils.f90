!>
!! This module contains utility routines needed for the initialization of the
!! NH model
!!
!! @author Guenther Zaengl, DWD
!!
!!
!! @par Revision History
!! First version by Guenther Zaengl, DWD (2011-06-29)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_nh_init_utils

  USE mo_kind,                  ONLY: wp
  USE mo_model_domain,          ONLY: t_patch
  USE mo_nonhydro_types,        ONLY: t_nh_state
  USE mo_nwp_phy_types,         ONLY: t_nwp_phy_diag, t_nwp_phy_tend, t_nwp_phy_stochconv
  USE mo_nwp_lnd_types,         ONLY: t_lnd_state, t_lnd_prog, t_lnd_diag, t_wtr_prog
  USE mo_ext_data_types,        ONLY: t_external_data
  USE mo_parallel_config,       ONLY: nproma
  USE mo_run_config,            ONLY: ntracer, ldass_lhn
  USE mo_grid_config,           ONLY: l_limited_area, n_dom
  USE mo_dynamics_config,       ONLY: nnow, nnow_rcf
  USE mo_physical_constants,    ONLY: grav, cpd, rd, cvd_o_rd, p0ref
  USE mo_vertical_coord_table,  ONLY: vct_b
  USE mo_impl_constants,        ONLY: nclass_aero, min_rlcell, REAL_T, SINGLE_T, INT_T, BOOL_T, &
    &                                 TLEV_NNOW, TLEV_NNOW_RCF, vname_len, LSS_JSBACH
  USE mo_math_constants,        ONLY: pi
  USE mo_exception,             ONLY: finish
  USE mo_sync,                  ONLY: sync_patch_array, SYNC_C
  USE mo_intp_data_strc,        ONLY: t_int_state
  USE mo_intp,                  ONLY: edges2cells_scalar
  USE mo_math_laplace,          ONLY: nabla2_scalar
  USE mo_math_gradients,        ONLY: grad_fd_norm
  USE mo_loopindices,           ONLY: get_indices_c, get_indices_e
  USE mo_io_config,             ONLY: var_in_output
  USE mo_initicon_types,        ONLY: t_saveinit_state, t_init_state, t_saved_field
  USE mo_initicon_config,       ONLY: type_iau_wgt, is_iau_active, &
    &                                 iau_wgt_dyn, iau_wgt_adv, ltile_coldstart
  USE mo_util_phys,             ONLY: virtual_temp
  USE mo_atm_phy_nwp_config,    ONLY: iprog_aero,atm_phy_nwp_config
  USE mo_lnd_nwp_config,        ONLY: ntiles_total, l2lay_rho_snow, ntiles_water, lmulti_snow, &
                                      nlev_soil, nlev_snow, lsnowtile, lprog_albsi, itype_trvg,&
                                      itype_snowevap
  USE mo_fortran_tools,         ONLY: init, copy
  USE mo_ifs_coord,             ONLY: geopot
  USE mo_radar_data_types,      ONLY : t_lhn_diag
  USE mo_var,                   ONLY: t_var
  USE mo_var_groups,            ONLY: var_groups_dyn
  USE mo_var_list_register,     ONLY: t_vl_register_iter
  USE mo_var_metadata,          ONLY: get_var_timelevel, get_var_name
  USE mo_hash_table,            ONLY: t_HashTable

#ifdef __PGI
  USE mo_util_texthash,         ONLY: t_char_workaround
#endif

#ifdef _OPENACC
  USE mo_mpi,                   ONLY: i_am_accel_node
#endif

  IMPLICIT NONE

  PRIVATE


  TYPE(t_saveinit_state), ALLOCATABLE  :: saveinit(:)

  ! subroutines
  !
  PUBLIC :: compute_smooth_topo 
  PUBLIC :: interp_uv_2_vn
  PUBLIC :: init_w, adjust_w
  PUBLIC :: convert_thdvars
  PUBLIC :: convert_omega2w
  PUBLIC :: save_initial_state, restore_initial_state
  PUBLIC :: compute_iau_wgt
  PUBLIC :: compute_input_pressure_and_height

CONTAINS

  !-------------
  !> Compute pressure and height of input data.
  !
  !  OUT: initicon%const%z_mc_in
  !       initicon%atm_in%pres
  !
  SUBROUTINE compute_input_pressure_and_height(p_patch, psfc, phi_sfc, initicon, opt_lmask)
    TYPE(t_patch),          INTENT(IN)       :: p_patch
    REAL(wp),               INTENT(INOUT)    :: psfc(:,:)
    REAL(wp),               INTENT(INOUT)    :: phi_sfc(:,:)
    CLASS(t_init_state),    INTENT(INOUT)    :: initicon
    LOGICAL, OPTIONAL,      INTENT(IN)       :: opt_lmask(:,:)
    ! LOCAL VARIABLES
    INTEGER :: jb, nlen, nlev_in
    INTEGER :: jc, jc1, jb1
    REAL(wp), DIMENSION(nproma,initicon%atm_in%nlev  ) :: delp, rdelp, rdlnpr, rdalpha, geop_mc
    REAL(wp), DIMENSION(nproma,initicon%atm_in%nlev,p_patch%nblks_c) :: temp_v_in
    REAL(wp), DIMENSION(nproma,initicon%atm_in%nlev+1) :: pres_ic, lnp_ic, geop_ic

    nlev_in = initicon%atm_in%nlev

    ! Compute virtual temperature of input data
    CALL virtual_temp(p_patch, initicon%atm_in%temp, initicon%atm_in%qv, initicon%atm_in%qc, &
                      initicon%atm_in%qi, initicon%atm_in%qr, initicon%atm_in%qs,            &
                      temp_v=temp_v_in)


    ! 1. Compute pressure and height of input data, using the IFS routines
    ! (note: no deep-atmosphere modifications are applied to these subroutines, 
    ! because they are assumed to be still on an "IFS-model-internal" level)

    ! If mask field is provided, fill data-void points (mask=.FALSE.) 
    ! with dummy value.
    IF (PRESENT(opt_lmask)) THEN
      !
      ! Detect first grid point for which the mask field is .true.
      outer: DO jb = 1, p_patch%nblks_c
        IF (jb /= p_patch%nblks_c) THEN
          nlen = nproma
        ELSE
         nlen = p_patch%npromz_c
        ENDIF
        inner: DO jc = 1, nlen
          IF (opt_lmask(jc,jb)) THEN
            jc1 = jc
            jb1 = jb
            EXIT outer
          ENDIF
        ENDDO inner
      ENDDO outer

      ! Do filling for psfc, phi_sfc, temp_v_in
!$OMP PARALLEL DO PRIVATE(jb,jc,nlen)
      DO jb = 1,p_patch%nblks_c

        IF (jb /= p_patch%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = p_patch%npromz_c
        ENDIF

        DO jc = 1, nlen
          IF (.NOT. opt_lmask(jc,jb)) THEN
            psfc(jc,jb)                = psfc(jc1,jb1)
            phi_sfc(jc,jb)             = phi_sfc(jc1,jb1)
            temp_v_in(jc,1:nlev_in,jb) = temp_v_in(jc1,1:nlev_in,jb1)
          ENDIF
        ENDDO
      ENDDO  ! jb
!$OMP END PARALLEL DO

    ENDIF


!$OMP PARALLEL
!$OMP DO PRIVATE(jb, nlen, pres_ic, lnp_ic, geop_ic, delp, rdelp, rdlnpr, &
!$OMP            rdalpha, geop_mc) ICON_OMP_DEFAULT_SCHEDULE

    DO jb = 1,p_patch%nblks_c

      IF (jb /= p_patch%nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = p_patch%npromz_c
      ENDIF
      
      ! Check if psfc is really psfc or LOG(psfc)
      IF (MAXVAL(psfc(1:nlen,jb)) <= 100._wp) THEN
        psfc(1:nlen,jb) = EXP(psfc(1:nlen,jb))
      ENDIF
      
      CALL initicon%const%vct%half_level_pressure(psfc(:,jb), nproma, nlen, nlev_in, pres_ic)
      
      CALL initicon%const%vct%full_level_pressure(pres_ic,nproma, nlen, nlev_in, initicon%atm_in%pres(:,:,jb))
      
      CALL initicon%const%vct%auxhyb(pres_ic, nproma, nlen, nlev_in,     & ! in
        delp, rdelp, lnp_ic, rdlnpr, rdalpha) ! out
      
      CALL geopot(temp_v_in(:,:,jb), rdlnpr, rdalpha, phi_sfc(:,jb), & ! in
        nproma, 1, nlen, nlev_in, geop_mc, geop_ic ) ! inout
      
      ! Compute 3D height coordinate field
      initicon%const%z_mc_in(1:nlen,1:nlev_in,jb) = geop_mc(1:nlen,1:nlev_in)/grav
      
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE compute_input_pressure_and_height




  !-------------
  !>
  !! SUBROUTINE convert_thdvars
  !! Converts the hydrostatic set of thermodynamic variables into the nonhydrostatic one
  !!
  !! Required input fields: pressure, virtual temperature
  !! Output: density, Exner pressure, virtual potential temperature
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-06-29)
  !!
  !!
  !!
  SUBROUTINE convert_thdvars(p_patch, pres, temp_v, &
                             rho, exner, theta_v    )


    TYPE(t_patch), INTENT(IN) :: p_patch

    ! Input fields - all defined at full model levels
    REAL(wp), INTENT(IN) :: pres  (:,:,:) ! pressure (Pa)
    REAL(wp), INTENT(IN) :: temp_v(:,:,:) ! virtual temperature (K)

    ! Output fields (prognostic model variables) - all defined at full model levels
    REAL(wp), INTENT(OUT) :: rho(:,:,:)        ! density (kg/m**3)
    REAL(wp), INTENT(OUT) :: exner(:,:,:)      ! Exner pressure
    REAL(wp), INTENT(OUT) :: theta_v(:,:,:)    ! virtual potential temperature (K)

    ! LOCAL VARIABLES
    INTEGER :: jb, jk, jc
    INTEGER :: nlen, nlev

    nlev = p_patch%nlev

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,jc) ICON_OMP_DEFAULT_SCHEDULE

    DO jb = 1, p_patch%nblks_c
      IF (jb /= p_patch%nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = p_patch%npromz_c
      ENDIF

      DO jk = 1, nlev
        DO jc = 1, nlen
          exner(jc,jk,jb)   = (pres(jc,jk,jb)/p0ref)**(rd/cpd)
          theta_v(jc,jk,jb) = temp_v(jc,jk,jb)/exner(jc,jk,jb)
          rho(jc,jk,jb)     = exner(jc,jk,jb)**cvd_o_rd*p0ref/rd/theta_v(jc,jk,jb)
        ENDDO
      ENDDO

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE convert_thdvars

  !-------------
  !>
  !! SUBROUTINE convert_omega2w
  !! Converts the hydrostatic vertical velocity (omega, Pa/s)
  !! into physical vertical velocity (m/s)
  !! Note: this routine has to be called on the input grid,
  !! where omega, pressure and temperature are not vertically staggered
  !!
  !! Required input fields: omega, pressure, temperature
  !! Output: vertical wind speed
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-20)
  !!
  !!
  !!
  SUBROUTINE convert_omega2w(omega, w, pres, temp, nblks, npromz, nlev, opt_lmask)


    ! Input fields
    REAL(wp), INTENT(IN) :: omega (:,:,:) ! omega (Pa/s)
    REAL(wp), INTENT(IN) :: pres  (:,:,:) ! pressure (Pa)
    REAL(wp), INTENT(IN) :: temp  (:,:,:) ! virtual temperature (K)

    ! Output
    REAL(wp), INTENT(OUT) :: w(:,:,:)  ! vertical velocity (m/s)

    ! Input dimension parameters
    INTEGER , INTENT(IN) :: nblks      ! Number of blocks
    INTEGER , INTENT(IN) :: npromz     ! Length of last block
    INTEGER , INTENT(IN) :: nlev       ! Number of model levels
    
    LOGICAL , INTENT(IN), OPTIONAL :: opt_lmask(:,:) ! logical mask of points to process

    ! LOCAL VARIABLES
    INTEGER :: jb, jk, jc
    INTEGER :: nlen

    IF(PRESENT(opt_lmask)) THEN

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,jc) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1, nblks
        IF (jb /= nblks) THEN
          nlen = nproma
        ELSE
          nlen = npromz
        ENDIF

        DO jk = 1, nlev
          DO jc = 1, nlen
            IF (opt_lmask(jc,jb)) THEN
              w(jc,jk,jb) = -rd*omega(jc,jk,jb)*temp(jc,jk,jb)/(grav*pres(jc,jk,jb))
            ELSE ! fill with dummy value
              w(jc,jk,jb) = 0._wp
            ENDIF
          ENDDO
        ENDDO

      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ELSE ! not present opt_lmask

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,jc) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1, nblks
        IF (jb /= nblks) THEN
          nlen = nproma
        ELSE
          nlen = npromz
        ENDIF

        DO jk = 1, nlev
          DO jc = 1, nlen
            w(jc,jk,jb) = -rd*omega(jc,jk,jb)*temp(jc,jk,jb)/(grav*pres(jc,jk,jb))
          ENDDO
        ENDDO

      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ENDIF ! PRESENT(opt_lmask)

  END SUBROUTINE convert_omega2w


  !-------------
  !>
  !! SUBROUTINE interp_uv_2_vn
  !! Interpolates u and v on cell points to vn on edge points
  !!
  !! Required input fields: u and v on cell points
  !! Output: vn on edge points
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-11)
  !!
  !!
  !!
  SUBROUTINE interp_uv_2_vn(p_patch, p_int, u, v, vn )


    TYPE(t_patch), TARGET, INTENT(IN)   :: p_patch
    TYPE(t_int_state),     INTENT(IN)   :: p_int

    ! Input fields - all defined at full model levels
    REAL(wp), INTENT(IN) :: u(:,:,:) ! zonal wind component on cell points (m/s)
    REAL(wp), INTENT(IN) :: v(:,:,:) ! meridional wind component on cell points (m/s)

    ! Output field (prognostic model variable) - defined at full model levels
    ! Intent (INOUT) because lateral nest boundaries cannot be filled here
    REAL(wp), INTENT(INOUT) :: vn(:,:,:)  ! edge-normal wind component (m/s)

    ! LOCAL VARIABLES
    INTEGER :: jb, jk, je
    INTEGER :: nlev, nblks_e, i_startblk,i_endblk, i_startidx,i_endidx
    REAL(wp) :: z_u, z_v

    INTEGER, DIMENSION(:,:,:), POINTER :: iidx, iblk

    nlev = p_patch%nlev
    nblks_e = p_patch%nblks_e

    iidx => p_patch%edges%cell_idx
    iblk => p_patch%edges%cell_blk

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

    IF (l_limited_area .OR. p_patch%id > 1) THEN ! Fill outermost nest boundary

      i_startblk = p_patch%edges%start_blk(1,1)
      i_endblk   = p_patch%edges%end_blk(1,1)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk,z_u,z_v) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, 1, 1)

        DO je = i_startidx, i_endidx
          IF (iidx(je,jb,1) >= 1 .AND. iblk(je,jb,1) >= 1) THEN
            DO jk = 1, nlev
              z_u = u(iidx(je,jb,1),jk,iblk(je,jb,1))
              z_v = v(iidx(je,jb,1),jk,iblk(je,jb,1))
              vn(je,jk,jb) =  z_u*p_patch%edges%primal_normal(je,jb)%v1 + &
                              z_v*p_patch%edges%primal_normal(je,jb)%v2
            END DO
          ELSE IF (iidx(je,jb,2) >= 1 .AND. iblk(je,jb,2) >= 1) THEN
            DO jk = 1, nlev
              z_u = u(iidx(je,jb,2),jk,iblk(je,jb,2))
              z_v = v(iidx(je,jb,2),jk,iblk(je,jb,2))
              vn(je,jk,jb) =  z_u*p_patch%edges%primal_normal(je,jb)%v1 + &
                              z_v*p_patch%edges%primal_normal(je,jb)%v2
            END DO
          ENDIF
        END DO

      END DO
!$OMP END DO
    ENDIF

    i_startblk = p_patch%edges%start_blk(2,1)

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE

      DO jb = i_startblk, nblks_e

        CALL get_indices_e(p_patch, jb, i_startblk, nblks_e, &
                           i_startidx, i_endidx, 2)

#ifdef __LOOP_EXCHANGE
        DO je = i_startidx, i_endidx
          DO jk = 1, nlev
#else
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
#endif

            vn(je,jk,jb) = p_int%c_lin_e(je,1,jb)                                              &
              *(u(iidx(je,jb,1),jk,iblk(je,jb,1))*p_patch%edges%primal_normal_cell(je,jb,1)%v1 &
              + v(iidx(je,jb,1),jk,iblk(je,jb,1))*p_patch%edges%primal_normal_cell(je,jb,1)%v2)&
              +            p_int%c_lin_e(je,2,jb)                                              &
              *(u(iidx(je,jb,2),jk,iblk(je,jb,2))*p_patch%edges%primal_normal_cell(je,jb,2)%v1 &
              + v(iidx(je,jb,2),jk,iblk(je,jb,2))*p_patch%edges%primal_normal_cell(je,jb,2)%v2 )

          ENDDO
        ENDDO
      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE interp_uv_2_vn

  !-------------
  !>
  !! SUBROUTINE init_w
  !! Initializes the vertical wind field based on the lower boundary condition
  !! w = v grad h and an empirical vertical decay function
  !! The discretization used here is simpler than that used in the dynamical core
  !! but is sufficient to avoid excessive generation of sound waves during the start phase
  !!
  !! Required input fields: vn, z_ifc
  !! Output: w
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-11)
  !!
  !!
  !!
  SUBROUTINE init_w(p_patch, p_int, vn, z_ifc, w)


    TYPE(t_patch), TARGET, INTENT(IN)   :: p_patch
    TYPE(t_int_state),     INTENT(IN)   :: p_int

    ! Input fields
    REAL(wp), INTENT(IN) :: vn(:,:,:)    ! edge-normal wind component (m/s)
    REAL(wp), INTENT(IN) :: z_ifc(:,:,:) ! height of half levels (m)

    ! Output field - defined at half model levels
    ! Intent (INOUT) because lateral nest boundaries cannot be filled here
    REAL(wp), INTENT(INOUT) :: w(:,:,:)  ! vertical wind component (m/s)

    ! LOCAL VARIABLES
    INTEGER :: jb, jk, je, jc, ktop
    INTEGER :: nlev, nlevp1, nblks_e, nblks_c, nshift, i_startblk, i_startidx, i_endidx

    REAL(wp) :: z_wsfc_e(nproma,1,p_patch%nblks_e) ! w at surface (edge points)
    REAL(wp) :: z_wsfc_c(nproma,1,p_patch%nblks_c) ! w at surface (cell points)
    REAL(wp) :: z_slope_e(nproma,p_patch%nlevp1,p_patch%nblks_e) ! slope at edges

    nlev    = p_patch%nlev
    nlevp1  = p_patch%nlevp1
    nblks_e = p_patch%nblks_e
    nblks_c = p_patch%nblks_c
    nshift  = p_patch%nshift_total

    ! In order to initialize w(1) = 0 except for vertical nesting
    IF (nshift == 0) THEN
      ktop = 2
    ELSE
      ktop = 1
    ENDIF

    ! Compute slope at edges
    CALL grad_fd_norm (z_ifc, p_patch, z_slope_e, 1, nlevp1)

    ! slope cannot be computed at outer boundary edges
    i_startblk = p_patch%edges%start_blk(2,1)

    DO jb = i_startblk, nblks_e

      CALL get_indices_e(p_patch, jb, i_startblk, nblks_e, &
                         i_startidx, i_endidx, 2)

      ! Extrapolation of vn to half levels is neglected here
      DO je = i_startidx, i_endidx
        z_wsfc_e(je,1,jb) = vn(je,nlev,jb)*z_slope_e(je,nlevp1,jb)
      ENDDO
    ENDDO

    CALL edges2cells_scalar(z_wsfc_e,p_patch,p_int%e_inn_c,z_wsfc_c,&
                            1,1,opt_rlstart=2)

    i_startblk = p_patch%cells%start_blk(2,1)

!$OMP PARALLEL
    ! First, initialize w with zero in order to avoid undefined nest boundary points
    CALL init(w(:,:,:))
!$OMP BARRIER

    ! specify a reasonable initial vertical wind speed
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, nblks_c

      CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, &
                         i_startidx, i_endidx, 2)

      DO jc = i_startidx, i_endidx
        w(jc,nlevp1,jb) = z_wsfc_c(jc,1,jb)
      ENDDO
      DO jk = nlev, ktop, -1
        DO jc = i_startidx, i_endidx
          w(jc,jk,jb) = z_wsfc_c(jc,1,jb)*vct_b(jk+nshift)
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE init_w


  !-------------
  !>
  !! SUBROUTINE adjust_w
  !! Computes the lower boundary condition for w in a similar way as init_w,
  !! but the result is then merged with an already available vertical wind
  !! field provided by an external data source.
  !!
  !! Required input fields: vn, w, z_ifc
  !! Output: w
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-18)
  !!
  !!
  !!
  SUBROUTINE adjust_w(p_patch, p_int, vn, z_ifc, w)


    TYPE(t_patch), TARGET, INTENT(IN)   :: p_patch
    TYPE(t_int_state),     INTENT(IN)   :: p_int

    ! Input fields
    REAL(wp), INTENT(IN) :: vn(:,:,:)    ! edge-normal wind component (m/s)
    REAL(wp), INTENT(IN) :: z_ifc(:,:,:) ! height of half levels (m)

    ! INOUT field - defined at half model levels
    REAL(wp), INTENT(INOUT) :: w(:,:,:)  ! vertical wind component (m/s)

    ! LOCAL VARIABLES
    INTEGER :: jb, jk, je, jc, ktop
    INTEGER :: nlev, nlevp1, nblks_e, nblks_c, nshift, i_startblk, i_startidx, i_endidx
    REAL(wp):: wfac

    REAL(wp) :: z_wsfc_e(nproma,1,p_patch%nblks_e) ! w at surface (edge points)
    REAL(wp) :: z_wsfc_c(nproma,1,p_patch%nblks_c) ! w at surface (cell points)
    REAL(wp) :: z_slope_e(nproma,p_patch%nlevp1,p_patch%nblks_e) ! slope at edges

    nlev    = p_patch%nlev
    nlevp1  = p_patch%nlevp1
    nblks_e = p_patch%nblks_e
    nblks_c = p_patch%nblks_c
    nshift  = p_patch%nshift_total

    ! In order to initialize w(1) = 0 except for vertical nesting
    IF (nshift == 0) THEN
      ktop = 2
    ELSE
      ktop = 1
    ENDIF

    ! Compute slope at edges
    CALL grad_fd_norm (z_ifc, p_patch, z_slope_e, 1, nlevp1)

    ! slope cannot be computed at outer boundary edges
    i_startblk = p_patch%edges%start_blk(2,1)

    DO jb = i_startblk, nblks_e

      CALL get_indices_e(p_patch, jb, i_startblk, nblks_e, &
                         i_startidx, i_endidx, 2)

      ! Extrapolation of vn to half levels is neglected here
      DO je = i_startidx, i_endidx
        z_wsfc_e(je,1,jb) = vn(je,nlev,jb)*z_slope_e(je,nlevp1,jb)
      ENDDO
    ENDDO

    CALL edges2cells_scalar(z_wsfc_e,p_patch,p_int%e_inn_c,z_wsfc_c,&
                            1,1,opt_rlstart=2)

    i_startblk = p_patch%cells%start_blk(2,1)

    ! specify lower boundary condition and merge with w field provided on input
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,wfac) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, nblks_c

      CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, &
                         i_startidx, i_endidx, 2)

      ! Lower boundary condition
      DO jc = i_startidx, i_endidx
        w(jc,nlevp1,jb) = z_wsfc_c(jc,1,jb)
      ENDDO

      ! Merging of lower boundary condition with interpolated data
      DO jk = nlev, ktop, -1
        wfac = vct_b(jk+nshift)**2
        DO jc = i_startidx, i_endidx
          w(jc,jk,jb) = (1._wp-wfac)*w(jc,jk,jb) + wfac*z_wsfc_c(jc,1,jb)
        ENDDO
      ENDDO

      IF (nshift == 0) THEN
        ! Upper boundary condition and smooth transition below
        ! if domain is not vertically nested
        DO jc = i_startidx, i_endidx
          w(jc,1,jb) = 0._wp
          w(jc,2,jb) = 0.33_wp*w(jc,2,jb)
          w(jc,3,jb) = 0.66_wp*w(jc,3,jb)
        ENDDO
      ENDIF

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE adjust_w


  !---------------------------------------------------------------------------
  !>
  !! Computes the smoothed topography needed for the SLEVE coordinate.
  !! May be bypassed once an option for reading the smooth topography from data
  !! is available
  !!
  !! @par Revision History
  !! Initial release by Guenther Zaengl, DWD, (2010-07-21)
  !!
  SUBROUTINE compute_smooth_topo(p_patch, p_int, topo_c, topo_smt_c)

    TYPE(t_patch),TARGET,INTENT(INOUT) :: p_patch
    TYPE(t_int_state), INTENT(IN) :: p_int

    ! Input fields: topography on cells
    REAL(wp), INTENT(IN) :: topo_c(:,:)

    ! Output fields: smooth topography on cells
    REAL(wp), INTENT(OUT) :: topo_smt_c(:,:)

    INTEGER  :: jb, jc, iter, niter
    INTEGER  :: i_startblk, nblks_c, i_startidx, i_endidx
    REAL(wp) :: z_topo(nproma,1,p_patch%nblks_c),nabla2_topo(nproma,1,p_patch%nblks_c)

    !-------------------------------------------------------------------------

    niter = 25 ! 25 smoothing iterations (do we need this to be a namelist variable?)

    ! Initialize auxiliary fields for topography with data and nullify nabla2 field
    z_topo(:,1,:)      = topo_c(:,:)
    nabla2_topo(:,1,:) = 0._wp

    i_startblk = p_patch%cells%start_blk(2,1)
    nblks_c    = p_patch%nblks_c

    CALL sync_patch_array(SYNC_C,p_patch,z_topo)

    ! Apply nabla2-diffusion niter times to create smooth topography
    DO iter = 1, niter

      CALL nabla2_scalar(z_topo, p_patch, p_int, nabla2_topo, &
        &                 slev=1, elev=1, rl_start=2, rl_end=min_rlcell )

      DO jb = i_startblk,nblks_c

        CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, &
                           i_startidx, i_endidx, 2)

        DO jc = i_startidx, i_endidx
          z_topo(jc,1,jb) = z_topo(jc,1,jb) + 0.125_wp*nabla2_topo(jc,1,jb) &
            &                               * p_patch%cells%area(jc,jb)
        ENDDO
      ENDDO

      CALL sync_patch_array(SYNC_C,p_patch,z_topo)

    ENDDO

    ! Store smooth topography on output fields
    topo_smt_c(:,:) = z_topo(:,1,:)

  END SUBROUTINE compute_smooth_topo



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
    if (.not. i_am_accel_node) CALL finish('restore_initial_state', 'This should be called in GPU mode only.')
#endif

    IF (ANY(atm_phy_nwp_config(:)%inwp_surface == LSS_JSBACH)) THEN
      CALL restore_var_group_state('jsb_init_vars', p_patch(:))
    END IF

    DO jg = 1, n_dom

      IF(.NOT. p_patch(jg)%ldom_active) CYCLE

      lnd_prog => p_lnd(jg)%prog_lnd(nnow_rcf(jg))
      lnd_diag => p_lnd(jg)%diag_lnd
      wtr_prog => p_lnd(jg)%prog_wtr(nnow_rcf(jg))

      !$ACC DATA COPYIN(saveinit(jg)%fr_seaice, saveinit(jg)%t_ice, saveinit(jg)%h_ice, saveinit(jg)%gz0) &
      !$ACC   COPYIN(saveinit(jg)%t_mnw_lk, saveinit(jg)%t_wml_lk, saveinit(jg)%h_ml_lk, saveinit(jg)%t_bot_lk) &
      !$ACC   COPYIN(saveinit(jg)%c_t_lk, saveinit(jg)%t_b1_lk, saveinit(jg)%h_b1_lk, saveinit(jg)%theta_v) &
      !$ACC   COPYIN(saveinit(jg)%rho, saveinit(jg)%exner, saveinit(jg)%w, saveinit(jg)%tke, saveinit(jg)%vn) &
      !$ACC   COPYIN(saveinit(jg)%tracer, saveinit(jg)%gz0_t, saveinit(jg)%t_g_t, saveinit(jg)%t_sk_t) &
      !$ACC   COPYIN(saveinit(jg)%qv_s_t, saveinit(jg)%freshsnow_t, saveinit(jg)%snowfrac_t) &
      !$ACC   COPYIN(saveinit(jg)%snowfrac_lc_t, saveinit(jg)%w_snow_t, saveinit(jg)%w_i_t, saveinit(jg)%h_snow_t) &
      !$ACC   COPYIN(saveinit(jg)%t_snow_t, saveinit(jg)%rho_snow_t, saveinit(jg)%w_so_t, saveinit(jg)%w_so_ice_t) &
      !$ACC   COPYIN(saveinit(jg)%t_so_t)

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

      !$ACC WAIT
      !$ACC END DATA

      IF (ntiles_total > 1 .AND. lsnowtile .AND. .NOT. ltile_coldstart) THEN
        !$ACC DATA COPYIN(saveinit(jg)%snowtile_flag_t, saveinit(jg)%idx_lst_t, saveinit(jg)%frac_t) &
        !$ACC   COPYIN(saveinit(jg)%gp_count_t)
        CALL copy(saveinit(jg)%snowtile_flag_t, ext_data(jg)%atm%snowtile_flag_t)
        CALL copy(saveinit(jg)%idx_lst_t, ext_data(jg)%atm%idx_lst_t)
        CALL copy(saveinit(jg)%frac_t, ext_data(jg)%atm%frac_t)
        CALL copy(saveinit(jg)%gp_count_t, ext_data(jg)%atm%gp_count_t)
        !$ACC WAIT
        !$ACC END DATA
      ENDIF

      IF (lmulti_snow) THEN
#ifdef _OPENACC
        CALL finish('restore_initial_state', 'lmulti_snow is not supported/tested with OpenACC')
#endif
        CALL copy(saveinit(jg)%t_snow_mult_t, lnd_prog%t_snow_mult_t)
        CALL copy(saveinit(jg)%rho_snow_mult_t, lnd_prog%rho_snow_mult_t)
        CALL copy(saveinit(jg)%wtot_snow_t, lnd_prog%wtot_snow_t)
        CALL copy(saveinit(jg)%wliq_snow_t, lnd_prog%wliq_snow_t)
        CALL copy(saveinit(jg)%dzh_snow_t, lnd_prog%dzh_snow_t)
      ELSE IF (l2lay_rho_snow) THEN
        CALL copy(saveinit(jg)%rho_snow_mult_t, lnd_prog%rho_snow_mult_t)
      ENDIF

      IF (atm_phy_nwp_config(jg)%lstoch_deep) THEN
#ifdef _OPENACC
        CALL finish('restore_initial_state', 'lstoch_deep is not supported/tested with OpenACC')
#endif
        CALL copy(saveinit(jg)%clnum_d, prm_stochconv(jg)%clnum_d)
        CALL copy(saveinit(jg)%clmf_d,  prm_stochconv(jg)%clmf_d)
      ENDIF
      IF (atm_phy_nwp_config(jg)%lstoch_sde) THEN
#ifdef _OPENACC
        CALL finish('restore_initial_state', 'lstoch_sde is not supported/tested with OpenACC')
#endif
        CALL copy(saveinit(jg)%clnum_a, prm_stochconv(jg)%clnum_a)
        CALL copy(saveinit(jg)%clmf_a,  prm_stochconv(jg)%clmf_a)
        CALL copy(saveinit(jg)%clnum_p, prm_stochconv(jg)%clnum_p)
        CALL copy(saveinit(jg)%clmf_p,  prm_stochconv(jg)%clmf_p)
      ENDIF

      IF (iprog_aero >= 1) THEN
#ifdef _OPENACC
        CALL finish('restore_initial_state', 'iprog_aero >= 1 is not supported/tested with OpenACC')
#endif        
        CALL copy(saveinit(jg)%aerosol, prm_diag(jg)%aerosol)
      ENDIF
      IF (lprog_albsi) THEN
        !$ACC DATA COPYIN(saveinit(jg)%alb_si)
        CALL copy(saveinit(jg)%alb_si, wtr_prog%alb_si)
        !$ACC WAIT
        !$ACC END DATA
      ENDIF
      IF (itype_trvg == 3) THEN
        !$ACC DATA COPYIN(saveinit(jg)%plantevap_t)
        CALL copy(saveinit(jg)%plantevap_t, lnd_diag%plantevap_t)
        !$ACC WAIT
        !$ACC END DATA
      ENDIF
      IF (itype_snowevap == 3) THEN
        !$ACC DATA COPYIN(saveinit(jg)%hsnow_max, saveinit(jg)%h_snow, saveinit(jg)%snow_age)
        CALL copy(saveinit(jg)%hsnow_max, lnd_diag%hsnow_max)
        CALL copy(saveinit(jg)%h_snow, lnd_diag%h_snow)
        CALL copy(saveinit(jg)%snow_age, lnd_diag%snow_age)
        !$ACC WAIT
        !$ACC END DATA
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
      CALL init (prm_diag(jg)%ice_gsp_rate)
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
      IF (ldass_lhn) CALL init (lhn_fields(jg)%brightband(:,:), -1._wp)

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
      !$ACC PARALLEL PRESENT(p_nh) ASYNC(1)
      !$ACC LOOP GANG VECTOR PRIVATE(je, jb)
      DO ic = 1, p_nh(jg)%metrics%nudge_e_dim
        je = p_nh(jg)%metrics%nudge_e_idx(ic)
        jb = p_nh(jg)%metrics%nudge_e_blk(ic)
        p_nh(jg)%diag%grf_tend_vn(je,:,jb) = 0._wp
      ENDDO
      !$ACC END PARALLEL

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

#ifdef __PGI
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

#ifdef __PGI
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

END MODULE mo_nh_init_utils
