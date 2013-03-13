!>
!! mo_sgs_turbulence
!!
!! Calculates subgrid-scale viscosity and diffusivity in the nonhydrostatic model
!! 
!! Additional preprocessing calculations need to be done in prep_les
!! 
!! @author Anurag Dipankar, MPI-M
!!
!!
!! @par Revision History
!! Initial release by Anurag Dipankar, MPI-M (2013-02-20)
!!
!! @par Copyright
!! 2002-2013 by DWD and MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!!    violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!    copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!    an according license agreement with DWD and MPI-M.
!! </ol>
!!
!! @par Warranty
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!!

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_sgs_turbulence

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, finish,message_text
  USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_model_domain,        ONLY: t_patch
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_intp_rbf,            ONLY: rbf_vec_interpol_vertex, rbf_vec_interpol_edge
  USE mo_intp,                ONLY: verts2edges_scalar, edges2verts_scalar, &
                                    cells2verts_scalar, cells2edges_scalar, &
                                    edges2cells_scalar, verts2cells_scalar, &
                                    edges2cells_vector
  USE mo_intp_rbf,            ONLY: rbf_vec_interpol_cell
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: ltimer
  USE mo_loopindices,         ONLY: get_indices_e, get_indices_c, get_indices_v
  USE mo_impl_constants    ,  ONLY: min_rledge, min_rlcell, min_rlvert, &
                                    min_rledge_int, min_rlcell_int, min_rlvert_int
  USE mo_math_constants,      ONLY: dbl_eps, pi
  USE mo_sync,                ONLY: SYNC_E, SYNC_C, SYNC_V, sync_patch_array
  USE mo_physical_constants,  ONLY: cpd, rcvd, p0ref, grav, rcpd, alv
 ! USE mo_timer,               ONLY: timer_nh_sgs_turbulence, timer_start, timer_stop
  USE mo_nwp_lnd_types,       ONLY: t_lnd_prog, t_wtr_prog, t_lnd_diag 
  USE mo_surface_les,         ONLY: rkarman, min_wind, surface_conditions 
  USE mo_nwp_phy_types,       ONLY: t_nwp_phy_diag, t_nwp_phy_tend

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  PUBLIC :: drive_subgrid_diffusion

  !Variables for the module
  REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: D_11c, D_12v, visc_smag_v, visc_smag_ie,  &
                                             diff_smag_e,  D_13ie, visc_smag_c, rho_e

  REAL(wp), PARAMETER :: inv_turb_prandtl = 3._wp
  REAL(wp), PARAMETER :: smag_const = 0.23_wp

  CONTAINS

  !>
  !! drive_subgrid_diffusion
  !!------------------------------------------------------------------------
  !! Driver for computing the sgs viscosity and diffusivity using Smagorisnky model
  !! and evaluating diffusion term on triangles
  !!------------------------------------------------------------------------
  !! @par Revision History
  !! Initial release by Anurag Dipankar, MPI-M (2013-03-05)
 
  SUBROUTINE drive_subgrid_diffusion(p_nh_prog, p_nh_prog_rcf, p_nh_diag, p_nh_metrics, p_patch, &
                                     p_int, p_prog_lnd_now, p_diag_lnd, prm_diag, prm_nwp_tend)

    TYPE(t_nh_prog),   INTENT(in)        :: p_nh_prog     !< single nh prognostic state
    TYPE(t_nh_prog),   INTENT(in)        :: p_nh_prog_rcf !< rcf nh prognostic state 
    TYPE(t_nh_diag),   INTENT(inout)     :: p_nh_diag     !< single nh diagnostic state
    TYPE(t_nh_metrics),INTENT(in),TARGET :: p_nh_metrics  !< single nh metric state
    TYPE(t_patch),     INTENT(in),TARGET :: p_patch       !< single patch
    TYPE(t_int_state), INTENT(in),TARGET :: p_int         !< single interpolation state
    TYPE(t_lnd_prog),  INTENT(inout)     :: p_prog_lnd_now!<land prog state 
    TYPE(t_lnd_diag),  INTENT(inout)     :: p_diag_lnd    !<land diag state 
    TYPE(t_nwp_phy_diag),   INTENT(inout):: prm_diag      !< atm phys vars
    TYPE(t_nwp_phy_tend), TARGET,INTENT(inout):: prm_nwp_tend    !< atm tend vars

    REAL(wp), ALLOCATABLE :: theta(:,:,:), theta_v(:,:,:)

    INTEGER :: nlev
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
    INTEGER :: rl_start, rl_end
    INTEGER :: jk, jb, jc

    !IF (ltimer) CALL timer_start(timer_nh_sgs_turbulence)

    nlev = p_patch%nlev
    i_nchdom   = MAX(1,p_patch%n_childdom)
    
    ALLOCATE( D_11c(nproma,nlev,p_patch%nblks_c),          &
              D_12v(nproma,nlev,p_patch%nblks_v),          &
              D_13ie(nproma,nlev+1,p_patch%nblks_e),       &                              
              visc_smag_v(nproma,nlev,p_patch%nblks_v),    &
              visc_smag_c(nproma,nlev,p_patch%nblks_c),    &
              visc_smag_ie(nproma,nlev+1,p_patch%nblks_e), &                              
              diff_smag_e(nproma,nlev,p_patch%nblks_e),    &
              theta(nproma,nlev,p_patch%nblks_c),          &
              theta_v(nproma,nlev,p_patch%nblks_c),        &
              rho_e(nproma,nlev,p_patch%nblks_e)           &
             )

    !Initialize
    D_11c       = 0._wp; D_12v       = 0._wp; visc_smag_ie = 0._wp
    visc_smag_v = 0._wp; diff_smag_e = 0._wp; visc_smag_c  = 0._wp
    D_13ie      = 0._wp; rho_e       = 0._wp; theta        = 0._wp
    theta_v     = 0._wp

    !Convert temperature to potential temperature: all routines within use theta
    !Assuming all prog variables are synced
    rl_start = 2
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk,i_endblk
       CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
       DO jc = i_startidx, i_endidx
         theta(jc,1:nlev,jb) = p_nh_diag%temp(jc,1:nlev,jb) / &
                               p_nh_prog%exner(jc,1:nlev,jb)

         theta_v(jc,1:nlev,jb) = p_nh_diag%tempv(jc,1:nlev,jb) / &
                                 p_nh_prog%exner(jc,1:nlev,jb)
       END DO
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    CALL sync_patch_array(SYNC_C, p_patch, theta)
    CALL sync_patch_array(SYNC_C, p_patch, theta_v)
         
    !CALL message('mo_sgs_turbulence:', 'CALL surface_conditions')
    CALL surface_conditions(p_nh_metrics, p_patch, p_nh_diag, p_int, &
                            p_prog_lnd_now, p_diag_lnd, prm_diag, theta)

    !CALL message('mo_sgs_turbulence:', 'CALL smagorinsky_model')
    CALL smagorinsky_model(p_nh_prog, p_nh_diag, p_nh_metrics, p_patch, p_int, prm_diag, &
                           theta_v, p_prog_lnd_now%t_g, p_diag_lnd%qv_s)

    !CALL message('mo_sgs_turbulence:', 'CALL diffuse_velocity')
    CALL diffuse_velocity(p_nh_prog, p_nh_diag, p_nh_metrics, p_patch, p_int, prm_diag, &
                          prm_nwp_tend%ddt_u_turb, prm_nwp_tend%ddt_v_turb)

    !CALL message('mo_sgs_turbulence:', 'CALL diffuse_scalar')
    CALL diffuse_scalar(theta, p_nh_metrics, p_patch, p_int, p_nh_diag,  &
                        prm_nwp_tend%ddt_temp_turb, p_nh_prog%exner,     &
                        prm_diag, p_nh_prog%rho, 'theta')

    !For qv and qc
!$OMP PARALLEL WORKSHARE
    prm_nwp_tend%ddt_tracer_turb = 0._wp
!$OMP END PARALLEL WORKSHARE

   DEALLOCATE(D_11c, D_12v, visc_smag_v, visc_smag_ie, diff_smag_e, &
              D_13ie, theta, visc_smag_c, rho_e, theta_v)

  !IF (ltimer) CALL timer_stop(timer_nh_sgs_turbulence)

  END SUBROUTINE drive_subgrid_diffusion

 
  !>
  !! smagorisnky_model
  !!------------------------------------------------------------------------
  !! Computes the sgs viscosity and diffusivity using Smagorisnky model
  !! \tau_ij = KD_ij where D_ij = du_i/dx_j+du_j/dx_i
  !!
  !! and  K = cs*\Delta**2*D/sqrt(2) where D = sqrt(D_ijD_ij)
  !!
  !! and D**2 = D_11**2 + D_22**2 + D_33**2 + 2D_12**2 + 2D_13**2 + 2D_23**2
  !!
  !! where, D_11 = 2du_1/dx_1
  !!        D_22 = 2d_u2/dx_2
  !!        D_33 = 2d_u3/dx_3
  !!        D_12 = du_1/dx_2 + du_2/dx_1
  !!        D_13 = du_1/dx_3 + du_3/dx_1
  !!        D_23 = du_2/dx_3 + du_3/dx_2
  !! For triangles: 1=normal, 2=tangential, and 3 = z directions
  !!------------------------------------------------------------------------
  !! @par Revision History
  !! Initial release by Anurag Dipankar, MPI-M (2013-02-20)
  SUBROUTINE  smagorinsky_model(p_nh_prog, p_nh_diag, p_nh_metrics, p_patch, p_int, &
                                prm_diag, theta_v, temp_sfc, qv_sfc)

    TYPE(t_patch),     INTENT(in),TARGET :: p_patch    !< single patch
    TYPE(t_int_state), INTENT(in),TARGET :: p_int      !< single interpolation state
    TYPE(t_nh_prog),   INTENT(in)        :: p_nh_prog  !< single nh prognostic state
    TYPE(t_nh_diag),   INTENT(in)        :: p_nh_diag  !< single nh diagnostic state
    TYPE(t_nh_metrics),INTENT(in),TARGET :: p_nh_metrics  !< single nh metric state
    REAL(wp),          INTENT(in)        :: theta_v(:,:,:)!potential temperature
    REAL(wp),          INTENT(in)        :: temp_sfc(:,:) !surface temperature
    REAL(wp),          INTENT(in)        :: qv_sfc(:,:)   !surface qv
    TYPE(t_nwp_phy_diag),   INTENT(inout):: prm_diag      !< atm phys vars

    ! local variables
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: u_vert, v_vert, w_vert, w_ie, vt_ie,     & 
                                               theta_v_ie, DD, D_11, D_12, D_13, z_me,  &
                                               theta_v_e, visc_smag_e, temp_sfc_e  
    REAL(wp), POINTER :: diff_smag_ic(:,:,:)

    REAL(wp) :: tmp(nproma,p_patch%nlev,p_patch%nblks_c)
    REAL(wp) :: vn_vert1, vn_vert2, vn_vert3, vn_vert4, dvt_norm, dvt_tang, w_full_c1
    REAL(wp) :: w_full_c2, w_full_v1, w_full_v2, brunt_vaisala_frq, les_filter
    REAL(wp) :: mixing_length_sq, kh_k, kh_km1
    REAL(wp) :: D_22, D_23, D_33

    INTEGER  :: nlev              !< number of full levels
    INTEGER,  DIMENSION(:,:,:), POINTER :: ividx, ivblk, iecidx, iecblk, ieidx, ieblk
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
    INTEGER :: rl_start, rl_end
    INTEGER :: jk, jb, jc, je, ic

    !--------------------------------------------------------------------------

    ! number of vertical levels
    nlev = p_patch%nlev
    i_nchdom   = MAX(1,p_patch%n_childdom)

    !Some constants
    diff_smag_ic => prm_diag%tkvh

    !Allocation
    ALLOCATE( u_vert(nproma,nlev,p_patch%nblks_v),   v_vert(nproma,nlev,p_patch%nblks_v),      &
              w_vert(nproma,nlev+1,p_patch%nblks_v), w_ie(nproma,nlev+1,p_patch%nblks_e),      &
              vt_ie(nproma,nlev+1,p_patch%nblks_e),  theta_v_ie(nproma,nlev+1,p_patch%nblks_e),&
              D_11(nproma,nlev,p_patch%nblks_e),     D_12(nproma,nlev,p_patch%nblks_e),        &
              DD(nproma,nlev,p_patch%nblks_e),       D_13(nproma,nlev,p_patch%nblks_e),        &
              theta_v_e(nproma,nlev,p_patch%nblks_e),  visc_smag_e(nproma,nlev,p_patch%nblks_e), &
              z_me(nproma,nlev,p_patch%nblks_e),     temp_sfc_e(nproma,1,p_patch%nblks_e) )
 
    !Initialize
    u_vert     = 0._wp; v_vert = 0._wp; w_vert = 0._wp; w_ie = 0._wp; vt_ie = 0._wp
    theta_v_ie = 0._wp; D_13   = 0._wp; D_11   = 0._wp; D_12 = 0._wp; DD    = 0._wp    
    theta_v_e  = 0._wp; z_me   = 0._wp; visc_smag_e = 0._wp; temp_sfc_e = 0._wp

    !--------------------------------------------------------------------------
    !1) Interpolate velocities at desired locations- mostly around the quadrilateral
    !
    !<It assumes that prog values are all synced at this stage while diag values might not>
    !--------------------------------------------------------------------------

    ! Include halo vertices because they might be used later on in the loop over edge
    CALL cells2verts_scalar(p_nh_prog%w, p_patch, p_int%cells_aw_verts, w_vert, &
                            opt_rlend=min_rlvert_int-1) 

    ! No need to include halo edge points for loop over interior edges
    CALL cells2edges_scalar(p_nh_prog%w, p_patch, p_int%c_lin_e, w_ie, opt_rlend=min_rledge_int)

    ! RBF reconstruction of velocity at vertices: include halos
    CALL rbf_vec_interpol_vertex( p_nh_prog%vn, p_patch, p_int, &
                                  u_vert, v_vert, opt_rlend=min_rlvert_int-1 )

    ! RBF reconstruction of tangential wind at level interface
    ! p_nh_diag%vn_ie is used that is computed in mo_solve_nonhydro
    ! therefore, sgs_turbulence better be called after that!
    CALL rbf_vec_interpol_edge(p_nh_diag%vn_ie, p_patch, p_int, vt_ie, opt_rlend=min_rledge_int)
    CALL sync_patch_array(SYNC_E, p_patch, vt_ie)

    !--------------------------------------------------------------------------
    !2) Compute horizontal strain rate tensor at full levels
    !--------------------------------------------------------------------------
    ividx => p_patch%edges%vertex_idx
    ivblk => p_patch%edges%vertex_blk

    iecidx => p_patch%edges%cell_idx
    iecblk => p_patch%edges%cell_blk

    ieidx => p_patch%cells%edge_idx
    ieblk => p_patch%cells%edge_blk

!$OMP PARALLEL PRIVATE(rl_start, rl_end, i_startblk,i_endblk)
    rl_start = 2
    rl_end   = min_rledge_int

    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,vn_vert1,vn_vert2,vn_vert3,vn_vert4,&
!$OMP            dvt_norm,dvt_tang,w_full_c1,w_full_c2,w_full_v1,w_full_v2,D_22,  &
!$OMP            D_23,D_33), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk,i_endblk

       CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
       DO jk = 1, nlev
         DO je = i_startidx, i_endidx

            vn_vert1 = u_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) * &
                       p_patch%edges%primal_normal_vert(je,jb,1)%v1 + &
                       v_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) * &
                       p_patch%edges%primal_normal_vert(je,jb,1)%v2

            vn_vert2 = u_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) * &
                       p_patch%edges%primal_normal_vert(je,jb,2)%v1 + &
                       v_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) * &
                       p_patch%edges%primal_normal_vert(je,jb,2)%v2

            dvt_tang = p_patch%edges%system_orientation(je,jb)* (   &
                       u_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) * &
                       p_patch%edges%dual_normal_vert(je,jb,2)%v1 + &
                       v_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) * &
                       p_patch%edges%dual_normal_vert(je,jb,2)%v2 - &
                      (u_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) * &
                       p_patch%edges%dual_normal_vert(je,jb,1)%v1 + &
                       v_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) * &
                       p_patch%edges%dual_normal_vert(je,jb,1)%v2) )

            vn_vert3 = u_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) * &
                       p_patch%edges%primal_normal_vert(je,jb,3)%v1 + &
                       v_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) * &
                       p_patch%edges%primal_normal_vert(je,jb,3)%v2

            vn_vert4 = u_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) * &
                       p_patch%edges%primal_normal_vert(je,jb,4)%v1 + &
                       v_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) * &
                       p_patch%edges%primal_normal_vert(je,jb,4)%v2

            dvt_norm = u_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) * &
                       p_patch%edges%dual_normal_vert(je,jb,4)%v1 + &
                       v_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) * &
                       p_patch%edges%dual_normal_vert(je,jb,4)%v2 - &
                      (u_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) * &
                       p_patch%edges%dual_normal_vert(je,jb,3)%v1 + &
                       v_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) * &
                       p_patch%edges%dual_normal_vert(je,jb,3)%v2)

            !W at full levels
            w_full_c1  = 0.5_wp * (                                        &
                         p_nh_prog%w(iecidx(je,jb,1),jk,iecblk(je,jb,1)) + &
                         p_nh_prog%w(iecidx(je,jb,1),jk+1,iecblk(je,jb,1)) ) 
 
            w_full_c2  = 0.5_wp * (                                        &
                         p_nh_prog%w(iecidx(je,jb,2),jk,iecblk(je,jb,2)) + &
                         p_nh_prog%w(iecidx(je,jb,2),jk+1,iecblk(je,jb,2)) ) 

            !W at full levels vertices from w at vertices at interface levels
            w_full_v1  = 0.5_wp * (                                    &
                         w_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) +    &
                         w_vert(ividx(je,jb,1),jk+1,ivblk(je,jb,1)) ) 

            w_full_v2  = 0.5_wp * (                                    &
                         w_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) +    &
                         w_vert(ividx(je,jb,2),jk+1,ivblk(je,jb,2)) ) 

            !Strain rates at edge center
            D_11(je,jk,jb) = 2._wp * (vn_vert4-vn_vert3) * &       
                             p_patch%edges%inv_vert_vert_length(je,jb)

            D_12(je,jk,jb) = p_patch%edges%system_orientation(je,jb) *  &
                            (vn_vert2-vn_vert1) *                       &
                             p_patch%edges%inv_primal_edge_length(je,jb)&
                             + dvt_norm *                               &
                             p_patch%edges%inv_vert_vert_length(je,jb)   
                    
            D_13(je,jk,jb) = (p_nh_diag%vn_ie(je,jk,jb)-p_nh_diag%vn_ie(je,jk+1,jb))/ &
                              p_nh_metrics%ddqz_z_full_e(je,jk,jb)  +                 &
                             (w_full_c2 - w_full_c1) *                                &
                              p_patch%edges%inv_dual_edge_length(je,jb)   

            D_22           = 2._wp * dvt_tang *  &
                             p_patch%edges%inv_primal_edge_length(je,jb)

            D_23           =(vt_ie(je,jk,jb) - vt_ie(je,jk+1,jb)) /          &
                             p_nh_metrics%ddqz_z_full_e(je,jk,jb)  +         &
                             p_patch%edges%system_orientation(je,jb) *       &
                            (w_full_v2 - w_full_v1) *                        &
                             p_patch%edges%inv_primal_edge_length(je,jb)   

            D_33           = 2._wp * (w_ie(je,jk,jb) - w_ie(je,jk+1,jb))/    &
                             p_nh_metrics%ddqz_z_full_e(je,jk,jb)               
           
            DD(je,jk,jb)   = D_11(je,jk,jb)**2 + D_22**2 + D_33**2  +        &
                   2._wp * ( D_12(je,jk,jb)**2 + D_13(je,jk,jb)**2 + D_23**2 )                  
         ENDDO
       ENDDO
    ENDDO
!$OMP END DO 


    !--------------------------------------------------------------------------
    !3) Classical Smagorinsky model with stability correction due to Lilly 1962
    !--------------------------------------------------------------------------

    ! 3(a)Calculate theta_v at half levels for gradient Richardson number and 
    ! some additional interpolations. 
    
    CALL cells2edges_scalar(theta_v, p_patch, p_int%c_lin_e, theta_v_e, opt_rlend=min_rledge_int)
    
    tmp(:,1,:) = temp_sfc(:,:) * ( 1._wp + 0.61_wp * qv_sfc(:,:) )
    CALL cells2edges_scalar(tmp, p_patch, p_int%c_lin_e, temp_sfc_e, opt_slev=1, opt_elev=1, &
                            opt_rlend=min_rledge_int)    
   
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx) ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk,i_endblk
       CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
       DO jk = 2 , nlev
         DO je = i_startidx, i_endidx
           theta_v_ie(je,jk,jb) = p_nh_metrics%wgtfac_e(je,jk,jb) * theta_v_e(je,jk,jb) + &
                        (1._wp - p_nh_metrics%wgtfac_e(je,jk,jb)) * theta_v_e(je,jk-1,jb)
         END DO
       END DO     
    END DO      
!$OMP END DO 
   
    !At top boundary: Extrapolating to save computation
    theta_v_ie(:,1,:) = 2._wp * theta_v_e(:,1,:) - theta_v_ie(:,2,:)
    
    !At surface
    theta_v_ie(:,nlev+1,:) = temp_sfc_e(:,1,:) 
 
    !additional interpolations    
    CALL cells2edges_scalar(p_nh_prog%rho, p_patch, p_int%c_lin_e, rho_e, opt_rlend=min_rledge_int)
    CALL cells2edges_scalar(p_nh_metrics%z_mc, p_patch, p_int%c_lin_e, z_me, opt_rlend=min_rledge_int)

    !3(b) Calculate stability corrected turbulent viscosity
    ! visc = mixing_length_sq * SQRT(DD/2) * SQRT(1-Ri/Pr) where Ri = (g/theta)*d_theta_dz/(DD/2), 
    ! where Brunt_vaisala_freq = (g/theta)*d_theta_dz. After simplification: 
    ! visc = mixing_length_sq * SQRT[ DD/2 - (Brunt_vaisala_frq/Pr) ]

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,les_filter,mixing_length_sq, &
!$OMP            brunt_vaisala_frq) ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk,i_endblk

       CALL get_indices_e(p_patch, jb, i_startblk, i_endblk,       &
                          i_startidx, i_endidx, rl_start, rl_end)
       DO jk = 1 , nlev 
         DO je = i_startidx, i_endidx

           !move calculation of mixing length to mo_vertical_grid
           les_filter = (0.5_wp*p_patch%edges%quad_area(je,jb)* &
                         p_nh_metrics%ddqz_z_full_e(je,jk,jb))**(1._wp/3._wp) 

           mixing_length_sq = (les_filter*smag_const*z_me(je,jk,jb)*rkarman)**2    &
                      / ((les_filter*smag_const)**2+(z_me(je,jk,jb)*rkarman)**2)

           brunt_vaisala_frq = grav * (theta_v_ie(je,jk,jb)-theta_v_ie(je,jk+1,jb)) / &
                              (theta_v_e(je,jk,jb)*p_nh_metrics%ddqz_z_full_e(je,jk,jb)) 

           visc_smag_e(je,jk,jb) = rho_e(je,jk,jb) * MAX( 0.1_wp, mixing_length_sq * &
                SQRT(MAX(0._wp, DD(je,jk,jb)*0.5_wp-inv_turb_prandtl*brunt_vaisala_frq)) ) 

         END DO
       END DO
    END DO
!$OMP END DO 

    CALL sync_patch_array(SYNC_E, p_patch, visc_smag_e)
    CALL sync_patch_array(SYNC_E, p_patch, D_11)
    CALL sync_patch_array(SYNC_E, p_patch, D_12)

    !--------------------------------------------------------------------------
    !4) Interpolate viscosity and strain rates to different locations: calculate them for 
    !   halos also because they will be used later in diffusion
    !--------------------------------------------------------------------------

    !visc and D_11 from edge to cell
    rl_start = 2
    rl_end   = min_rlcell_int-2

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx) ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk,i_endblk
       CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
       DO jk = 1 , nlev
         DO jc = i_startidx, i_endidx
      
           visc_smag_c(jc,jk,jb) =                                                       &
              (visc_smag_e(ieidx(jc,jb,1),jk,ieblk(jc,jb,1))*p_int%e_bln_c_s(jc,1,jb)  + &
               visc_smag_e(ieidx(jc,jb,2),jk,ieblk(jc,jb,2))*p_int%e_bln_c_s(jc,2,jb)  + &
               visc_smag_e(ieidx(jc,jb,3),jk,ieblk(jc,jb,3))*p_int%e_bln_c_s(jc,3,jb))           
      
           D_11c(jc,jk,jb) =                                                      &
              (D_11(ieidx(jc,jb,1),jk,ieblk(jc,jb,1))*p_int%e_bln_c_s(jc,1,jb)  + &
               D_11(ieidx(jc,jb,2),jk,ieblk(jc,jb,2))*p_int%e_bln_c_s(jc,2,jb)  + &
               D_11(ieidx(jc,jb,3),jk,ieblk(jc,jb,3))*p_int%e_bln_c_s(jc,3,jb))           

           !temporary fix
           tmp(jc,jk,jb) =                                                        &
              (D_12(ieidx(jc,jb,1),jk,ieblk(jc,jb,1))*p_int%e_bln_c_s(jc,1,jb)  + &
               D_12(ieidx(jc,jb,2),jk,ieblk(jc,jb,2))*p_int%e_bln_c_s(jc,2,jb)  + &
               D_12(ieidx(jc,jb,3),jk,ieblk(jc,jb,3))*p_int%e_bln_c_s(jc,3,jb))           

         END DO
       END DO
    END DO
!$OMP END DO

    !now interpolate visc and D_12 to vertices: TEMPORARY FIX
    CALL cells2verts_scalar(visc_smag_c, p_patch, p_int%cells_aw_verts, visc_smag_v, &
                           opt_rlend=min_rlvert_int-1) 
    CALL cells2verts_scalar(tmp, p_patch, p_int%cells_aw_verts, D_12v, opt_rlend=min_rlvert_int-1) 

     
    !Now calculate visc_smag at half levels at edge

    !vertical derivative of visc_smag_ie is then calculated- therefore no need to get
    !its values on Halos
    rl_start = 2
    rl_end   = min_rledge_int 

    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)

    !viscosity at interface levels
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx) ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk,i_endblk
       CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
       DO jk = 2 , nlev
         DO je = i_startidx, i_endidx
           !Arithmetic mean
           visc_smag_ie(je,jk,jb) = p_nh_metrics%wgtfac_e(je,jk,jb) * visc_smag_e(je,jk,jb) + &
                          (1._wp - p_nh_metrics%wgtfac_e(je,jk,jb)) * visc_smag_e(je,jk-1,jb)
            
           !Harmonic mean
           !visc_smag_ie(je,jk,jb) = visc_smag_e(je,jk,jb) * visc_smag_e(je,jk-1,jb) / &
           !             ( p_nh_metrics%wgtfac_e(je,jk,jb) * visc_smag_e(je,jk,jb) +   &
           !              (1._wp - p_nh_metrics%wgtfac_e(je,jk,jb)) * visc_smag_e(je,jk-1,jb) )

           !D_13ie(je,jk,jb) = p_nh_metrics%wgtfac_e(je,jk,jb) * D_13(je,jk,jb) + &
           !               (1._wp - p_nh_metrics%wgtfac_e(je,jk,jb)) * D_13(je,jk-1,jb)
         END DO
       END DO     
    END DO      
!$OMP END DO 

    !At the TOP boundary
    visc_smag_ie(:,1,:) = 0._wp
    !D_13ie(:,1,:)       = D_13ie(:,2,:) 
    
 
    !At the bottom: temporary fix because we will use surface flux directly while solving
    !the diffusion equation. Need to get proper surface scheme and derive visc_smag_ie
    !from there
    visc_smag_ie(:,nlev+1,:) = visc_smag_ie(:,nlev,:) 
    !D_13ie(:,nlev+1,:)       = D_13ie(:,nlev,:) 

    !--------------------------------------------------------------------------
    !5) Calculate turbulent diffusivity
    !--------------------------------------------------------------------------

    !Turbulent diffusivity at edge at full levels
!$OMP WORKSHARE
    diff_smag_e(:,:,:) = visc_smag_e(:,:,:) * inv_turb_prandtl
!$OMP END WORKSHARE

    !Turbulent diffusivity at cell center at half levels and like visc_smag_ie 
    !it is also calculated for interior points only
    rl_start = 2
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,kh_k,kh_km1) ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk,i_endblk
       CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
       DO jk = 2 , nlev
         DO jc = i_startidx, i_endidx

           kh_k    = visc_smag_c(jc,jk,jb)   * inv_turb_prandtl
           kh_km1  = visc_smag_c(jc,jk-1,jb) * inv_turb_prandtl        
       
           !Harmonic mean
           !diff_smag_ic(jc,jk,jb) = kh_k * kh_km1 / (                               & 
           !                         p_nh_metrics%wgtfac_c(jc,jk,jb) * kh_k +        &
           !                        (1._wp - p_nh_metrics%wgtfac_c(jc,jk,jb)) * kh_km1 )

           !Arithmetic mean
           diff_smag_ic(jc,jk,jb) = p_nh_metrics%wgtfac_c(jc,jk,jb) * kh_k +        &
                                   (1._wp - p_nh_metrics%wgtfac_c(jc,jk,jb)) * kh_km1

         END DO
       END DO     
    END DO      
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !At the TOP boundary
    diff_smag_ic(:,1,:) = 0._wp

    !At the bottom: temporary fix because we will use surface flux directly while solving
    !the diffusion equation. Need to get proper surface scheme and derive diff_smag_ie
    !from there
    diff_smag_ic(:,nlev+1,:) = diff_smag_ic(:,nlev,:) 

    !DEALLOCATE variables
    DEALLOCATE( u_vert, v_vert, w_vert, w_ie, vt_ie, theta_v_ie, D_11, D_12, DD, D_13, &
                visc_smag_e, theta_v_e, temp_sfc_e, z_me )
  
  END SUBROUTINE smagorinsky_model

  !>
  !! diffuse_velocity
  !!------------------------------------------------------------------------
  !! Calculate the SGS diffusion term for normal velocity component
  !! - Uses the forward Euler time scheme in time split (sequential) manner adopted
  !!   in the NH version.
  !! - The Euler scheme is generally known to be unconditionaly unstable for diffusion term
  !!   but it is not certain when implemented in time split manner
  !! - This needs to be investigated and if found unstable we will switch to MacKormack 
  !!   scheme which is known to be stable
  !!  
  !! d_vn/d_t =  d_tau_11/d_x1 + d_tau_12/d_x2 + d_tau_13/d_x3
  !!------------------------------------------------------------------------
  !! @par Revision History
  !! Initial release by Anurag Dipankar, MPI-M (2013-02-05)
  SUBROUTINE diffuse_velocity(p_nh_prog, p_nh_diag, p_nh_metrics, p_patch, p_int, &
                              prm_diag, ddt_u, ddt_v)

    TYPE(t_nh_prog),   INTENT(in)        :: p_nh_prog    !< single nh prognostic state
    TYPE(t_nh_diag),   INTENT(in)        :: p_nh_diag    !< single nh diagnostic state
    TYPE(t_nh_metrics),INTENT(in),TARGET :: p_nh_metrics !< single nh metric state
    TYPE(t_patch), TARGET, INTENT(in)    :: p_patch      !< single patch
    TYPE(t_int_state), INTENT(in),TARGET :: p_int        !< single interpolation state
    TYPE(t_nwp_phy_diag),  INTENT(in)    :: prm_diag     !< atm phys vars
    REAL(wp),   TARGET, INTENT(out)      :: ddt_u(:,:,:) !< u tendency
    REAL(wp),   TARGET, INTENT(out)      :: ddt_v(:,:,:) !< v tendency

    REAL(wp) :: flux_up, flux_dn, d_tau_11_d_x1, d_tau_12_d_x2, d_tau_13_d_x3
    REAL(wp) :: stress_uc, stress_vc, stress_c1n, stress_c2n, inv_mwind
    REAL(wp) :: inv_rhoe(nproma,p_patch%nlev,p_patch%nblks_e)
    REAL(wp) :: tot_tend(nproma,p_patch%nlev,p_patch%nblks_e)
    REAL(wp) :: dz_half_e(nproma,p_patch%nlev+1,p_patch%nblks_e)

    INTEGER,  DIMENSION(:,:,:), POINTER :: ividx, ivblk, iecidx, iecblk
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
    INTEGER :: rl_start, rl_end
    INTEGER :: jk, jb, je, jcn, jbn
    INTEGER  :: nlev              

    ! number of vertical levels
    nlev     = p_patch%nlev
    i_nchdom = MAX(1,p_patch%n_childdom)
   
    ividx  => p_patch%edges%vertex_idx
    ivblk  => p_patch%edges%vertex_blk

    iecidx => p_patch%edges%cell_idx
    iecblk => p_patch%edges%cell_blk
   
    CALL cells2edges_scalar(p_nh_metrics%ddqz_z_half, p_patch, p_int%c_lin_e, dz_half_e, &
                            opt_rlend=min_rledge_int)

    !Some initializations

    !Inverse of density (global rho_e in this module)
    inv_rhoe = 1._wp / rho_e

    !total tendency
    tot_tend = 0._wp

    ! 1) First get the horizontal tendencies

    rl_start = 2
    rl_end   = min_rledge_int

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,flux_up,flux_dn,d_tau_11_d_x1, &
!$OMP            d_tau_12_d_x2), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk,i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)
      DO jk = 1 , nlev
       DO je = i_startidx, i_endidx
        
         !tendency in normal direction
         flux_up = visc_smag_c(iecidx(je,jb,2),jk,iecblk(je,jb,2)) * &
                   D_11c(iecidx(je,jb,2),jk,iecblk(je,jb,2))        

         flux_dn = visc_smag_c(iecidx(je,jb,1),jk,iecblk(je,jb,1)) * &
                   D_11c(iecidx(je,jb,1),jk,iecblk(je,jb,1))                

         d_tau_11_d_x1 = (flux_up - flux_dn) * p_patch%edges%inv_dual_edge_length(je,jb)   

         !tendency in tangential direction
         flux_up = visc_smag_v(ividx(je,jb,2),jk,ivblk(je,jb,2)) * &
                   D_12v(ividx(je,jb,2),jk,ivblk(je,jb,2))        

         flux_dn = visc_smag_v(ividx(je,jb,1),jk,ivblk(je,jb,1)) * &
                   D_12v(ividx(je,jb,1),jk,ivblk(je,jb,1))                

         d_tau_12_d_x2 = p_patch%edges%system_orientation(je,jb) * &
                        (flux_up - flux_dn) * p_patch%edges%inv_primal_edge_length(je,jb)   

         tot_tend(je,jk,jb) = tot_tend(je,jk,jb) + &
                              (d_tau_11_d_x1 + d_tau_12_d_x2) * inv_rhoe(je,jk,jb)

       END DO
      END DO
    END DO                         
!$OMP END DO 

   ! 2) Vertical tendency

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,flux_up,flux_dn, &
!$OMP            d_tau_13_d_x3), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk,i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)
      DO jk = 2 , nlev-1
       DO je = i_startidx, i_endidx

         !tendency in vertical direction
         flux_up = visc_smag_ie(je,jk,jb) *                             &
                ( (p_nh_prog%vn(je,jk-1,jb) - p_nh_prog%vn(je,jk,jb)) / &
                   dz_half_e(je,jk,jb) +                                &
                  (p_nh_prog%w(iecidx(je,jb,2),jk,iecblk(je,jb,2)) -    &
                   p_nh_prog%w(iecidx(je,jb,1),jk,iecblk(je,jb,1))) *   &
                   p_patch%edges%inv_dual_edge_length(je,jb) )
                   
         !flux_up = visc_smag_ie(je,jk,jb) * D_13ie(je,jk,jb) 
                         
         flux_dn = visc_smag_ie(je,jk+1,jb) *                           &
                ( (p_nh_prog%vn(je,jk,jb) - p_nh_prog%vn(je,jk+1,jb)) / &
                   dz_half_e(je,jk+1,jb) +                              &
                  (p_nh_prog%w(iecidx(je,jb,2),jk+1,iecblk(je,jb,2)) -  &
                   p_nh_prog%w(iecidx(je,jb,1),jk+1,iecblk(je,jb,1))) * &
                   p_patch%edges%inv_dual_edge_length(je,jb) )
                  
         !flux_dn = visc_smag_ie(je,jk+1,jb) * D_13ie(je,jk+1,jb) 

         d_tau_13_d_x3 = (flux_up - flux_dn) / p_nh_metrics%ddqz_z_full_e(je,jk,jb)   

         tot_tend(je,jk,jb) = tot_tend(je,jk,jb) + d_tau_13_d_x3 * inv_rhoe(je,jk,jb)

         !Advance velocity
         !p_nh_prog%vn(je,jk,jb) = p_nh_prog%vn(je,jk,jb) + dtime * tot_tend(je,jk,jb)

       END DO
      END DO
    END DO  
!$OMP END DO 

   !-----------------------------------------------------------------
   !3) Boundary treatment in vertical: use surface fluxes from 
   !   surface_conditions
   !-----------------------------------------------------------------
   
   !-----------------------------------------------------------------
   jk = 1
   !-----------------------------------------------------------------

!$OMP DO PRIVATE(jb,je,i_startidx,i_endidx,flux_up,flux_dn, &
!$OMP            d_tau_13_d_x3), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk,i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)
       DO je = i_startidx, i_endidx

         flux_up = 0._wp
                   
         flux_dn = visc_smag_ie(je,jk+1,jb) *                           &
                ( (p_nh_prog%vn(je,jk,jb) - p_nh_prog%vn(je,jk+1,jb)) / &
                   dz_half_e(je,jk+1,jb) +                              &
                  (p_nh_prog%w(iecidx(je,jb,2),jk+1,iecblk(je,jb,2)) -  &
                   p_nh_prog%w(iecidx(je,jb,1),jk+1,iecblk(je,jb,1))) * &
                   p_patch%edges%inv_dual_edge_length(je,jb) )
                  
         !flux_dn = visc_smag_ie(je,jk+1,jb) * D_13ie(je,jk+1,jb) 

         d_tau_13_d_x3 = (flux_up - flux_dn) / p_nh_metrics%ddqz_z_full_e(je,jk,jb)   

         tot_tend(je,jk,jb) = tot_tend(je,jk,jb) + d_tau_13_d_x3 * inv_rhoe(je,jk,jb)

         !Advance velocity
         !p_nh_prog%vn(je,jk,jb) = p_nh_prog%vn(je,jk,jb) + dtime * tot_tend(je,jk,jb)

       END DO
    END DO                         
!$OMP END DO 

   !-----------------------------------------------------------------
    jk = nlev
   !-----------------------------------------------------------------

!$OMP DO PRIVATE(jb,je,i_startidx,i_endidx,flux_up,flux_dn,jcn,jbn,inv_mwind,&
!$OMP            stress_uc,stress_vc,stress_c1n,stress_c2n,d_tau_13_d_x3),   &
!$OMP            ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk,i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)
       DO je = i_startidx, i_endidx

         flux_up = visc_smag_ie(je,jk,jb) *                             &
                ( (p_nh_prog%vn(je,jk-1,jb) - p_nh_prog%vn(je,jk,jb)) / &
                   dz_half_e(je,jk,jb) +                                &
                  (p_nh_prog%w(iecidx(je,jb,2),jk,iecblk(je,jb,2)) -    &
                   p_nh_prog%w(iecidx(je,jb,1),jk,iecblk(je,jb,1))) *   &
                   p_patch%edges%inv_dual_edge_length(je,jb) )
                  
         !flux_up = visc_smag_ie(je,jk,jb) * D_13ie(je,jk,jb) 

         !Get net shear stress in the direction of vn at surface

         !shear stress in normal direction from cell 1 
         jcn = iecidx(je,jb,1)
         jbn = iecblk(je,jb,1)        

         inv_mwind = 1._wp / MAX( min_wind, &
                     SQRT(p_nh_diag%u(jcn,jk,jbn)**2+p_nh_diag%v(jcn,jk,jbn)**2))        

         stress_uc = p_nh_diag%u(jcn,jk,jbn) * prm_diag%umfl_s(jcn,jbn) * inv_mwind 
         stress_vc = p_nh_diag%v(jcn,jk,jbn) * prm_diag%umfl_s(jcn,jbn) * inv_mwind 

         stress_c1n = stress_uc*p_patch%edges%primal_normal_cell(je,jb,1)%v1 + &
                      stress_vc*p_patch%edges%primal_normal_cell(je,jb,1)%v2

         !shear stress in normal direction from cell 2 
         jcn = iecidx(je,jb,2)
         jbn = iecblk(je,jb,2)        

         inv_mwind = 1._wp / MAX( min_wind, &
                     SQRT(p_nh_diag%u(jcn,jk,jbn)**2+p_nh_diag%v(jcn,jk,jbn)**2))        

         stress_uc = p_nh_diag%u(jcn,jk,jbn) * prm_diag%umfl_s(jcn,jbn) * inv_mwind 
         stress_vc = p_nh_diag%v(jcn,jk,jbn) * prm_diag%umfl_s(jcn,jbn) * inv_mwind 

         stress_c2n = stress_uc*p_patch%edges%primal_normal_cell(je,jb,2)%v1 + &
                      stress_vc*p_patch%edges%primal_normal_cell(je,jb,2)%v2

         !Net stress at the edge
         flux_dn    = stress_c1n * p_int%c_lin_e(je,1,jb) + &
                      stress_c2n * p_int%c_lin_e(je,2,jb) 

         d_tau_13_d_x3 = (flux_up - flux_dn) / p_nh_metrics%ddqz_z_full_e(je,jk,jb)   

         tot_tend(je,jk,jb) = tot_tend(je,jk,jb) + d_tau_13_d_x3 * inv_rhoe(je,jk,jb)

         !Advance velocity
         !p_nh_prog%vn(je,jk,jb) = p_nh_prog%vn(je,jk,jb) + dtime * tot_tend(je,jk,jb)

       END DO
    END DO                         
!$OMP END DO NOWAIT        
!$OMP END PARALLEL

    !Get turbulent tendency at cell center
    CALL sync_patch_array(SYNC_E, p_patch, tot_tend)
    CALL rbf_vec_interpol_cell(tot_tend, p_patch, p_int, ddt_u, ddt_v, &
                               opt_rlend=min_rlcell_int-1)

  END SUBROUTINE diffuse_velocity

  !>
  !! diffuse_scalar
  !!------------------------------------------------------------------------
  !! Calculate the SGS diffusion term for cell based scalars
  !! - Uses the forward Euler time scheme in time split (sequential) manner adopted
  !!   in the NH version.
  !! - The Euler scheme is generally known to be unconditionaly unstable for diffusion term
  !!   but it is not certain when implemented in time split manner
  !! - This needs to be investigated and if found unstable we will switch to MacKormack 
  !!   scheme which is known to be stable
  !!  
  !!------------------------------------------------------------------------
  !! @par Revision History
  !! Initial release by Anurag Dipankar, MPI-M (2013-02-05)
  SUBROUTINE diffuse_scalar(var, p_nh_metrics, p_patch, p_int, p_nh_diag, tot_tend,  &
                            exner, prm_diag, rho, scalar_name)

    REAL(wp),          INTENT(in)        :: var(:,:,:)   !input scalar
    TYPE(t_nh_metrics),INTENT(in),TARGET :: p_nh_metrics !< single nh metric state
    TYPE(t_patch),     INTENT(in),TARGET :: p_patch      !< single patch
    TYPE(t_int_state), INTENT(in),TARGET :: p_int        !< single interpolation state
    TYPE(t_nh_diag),   INTENT(in)        :: p_nh_diag    !< single nh diagnostic state
    REAL(wp),          INTENT(out),TARGET:: tot_tend(:,:,:)!<total tendency
    REAL(wp),          INTENT(in)        :: exner(:,:,:)   !
    REAL(wp),          INTENT(in)        :: rho(:,:,:)     !density at cell center
    TYPE(t_nwp_phy_diag),  INTENT(in)    :: prm_diag     !< atm phys vars
    CHARACTER(*), INTENT(in)             :: scalar_name    

    REAL(wp), POINTER :: diff_smag_ic(:,:,:)
    REAL(wp) :: flux_up, flux_dn, tend_hor, tend_ver
    REAL(wp) :: nabla2_e(nproma,p_patch%nlev,p_patch%nblks_e)
    REAL(wp) :: fac(nproma,p_patch%nlev,p_patch%nblks_c)
    REAL(wp) :: sflux(nproma,p_patch%nblks_c)
     
    INTEGER,  DIMENSION(:,:,:), POINTER :: iecidx, iecblk, ieidx, ieblk
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
    INTEGER :: rl_start, rl_end
    INTEGER :: jk, jb, je, jc
    INTEGER  :: nlev              

    ! number of vertical levels
    nlev = p_patch%nlev
    i_nchdom   = MAX(1,p_patch%n_childdom)

    iecidx => p_patch%edges%cell_idx
    iecblk => p_patch%edges%cell_blk
   
    ieidx => p_patch%cells%edge_idx
    ieblk => p_patch%cells%edge_blk

    diff_smag_ic => prm_diag%tkvh

    !Special treatment for different scalars: includes
    !boundary treatment also
    IF(TRIM(scalar_name)=='theta')THEN
      !multiply by exner to convert from theta tend to temp tend

!$OMP PARALLEL WORKSHARE
      fac(:,:,:) = cpd * rcvd * exner(:,:,:) / rho(:,:,:)
      sflux(:,:) = prm_diag%shfl_s(:,:) * rcpd
!$OMP END PARALLEL WORKSHARE

    ELSE

!$OMP PARALLEL WORKSHARE
      fac(:,:,:) = 1._wp / rho(:,:,:)
      sflux(:,:) = prm_diag%lhfl_s(:,:) / alv
!$OMP END PARALLEL WORKSHARE

    END IF

    ! use conservative discretization div(k*grad(var))-horizontal part  
    ! following mo_nh_diffusion

!$OMP PARALLEL PRIVATE(rl_start, rl_end, i_startblk,i_endblk)

    !include halo points because these values will be used in next loop
    rl_start = 2
    rl_end   = min_rledge_int-1 

    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jk,je,jb,i_startidx,i_endidx), ICON_OMP_RUNTIME_SCHEDULE
        DO jb = i_startblk,i_endblk

          CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)

          ! compute kh_smag_e * grad(var) 
          DO jk = 1, nlev
            DO je = i_startidx, i_endidx
                nabla2_e(je,jk,jb) = diff_smag_e(je,jk,jb) *          &
                         p_patch%edges%inv_dual_edge_length(je,jb)*   &
                        (var(iecidx(je,jb,2),jk,iecblk(je,jb,2)) -    &
                         var(iecidx(je,jb,1),jk,iecblk(je,jb,1)))
            ENDDO
          ENDDO
        ENDDO
!$OMP END DO

        ! now compute the divergence of the quantity above
        rl_start = 2
        rl_end   = min_rlcell_int

        i_startblk = p_patch%cells%start_blk(rl_start,1)
        i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jc,jb,i_startidx,i_endidx,tend_hor,flux_up,flux_dn,tend_ver), ICON_OMP_RUNTIME_SCHEDULE
        DO jb = i_startblk,i_endblk

          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)
          DO jk = 2, nlev-1
            DO jc = i_startidx, i_endidx
              !horizontal tendency
              tend_hor           =                                                     &
                nabla2_e(ieidx(jc,jb,1),jk,ieblk(jc,jb,1))*p_int%geofac_div(jc,1,jb) + &
                nabla2_e(ieidx(jc,jb,2),jk,ieblk(jc,jb,2))*p_int%geofac_div(jc,2,jb) + &
                nabla2_e(ieidx(jc,jb,3),jk,ieblk(jc,jb,3))*p_int%geofac_div(jc,3,jb)

              !vertical tendency
              flux_up = diff_smag_ic(jc,jk,jb) * (var(jc,jk-1,jb) - var(jc,jk,jb)) / &
                        p_nh_metrics%ddqz_z_half(jc,jk,jb) 

              flux_dn = diff_smag_ic(jc,jk+1,jb) * (var(jc,jk,jb) - var(jc,jk+1,jb)) / &
                        p_nh_metrics%ddqz_z_half(jc,jk+1,jb) 

              tend_ver = (flux_up - flux_dn) * p_nh_metrics%inv_ddqz_z_full(jc,jk,jb) 
             
              !Advance scalar
              tot_tend(jc,jk,jb) = (tend_hor + tend_ver) * fac(jc,jk,jb)
              !var(jc,jk,jb) = var(jc,jk,jb) + dtime * tot_tend(jc,jk,jb)
            ENDDO
          ENDDO
        ENDDO
!$OMP END DO

        !Boundary treatment
       
        !--------------------------------------------------------    
        jk = 1        
        !--------------------------------------------------------    
!$OMP DO PRIVATE(jc,jb,i_startidx,i_endidx,tend_hor,flux_up,flux_dn,tend_ver), ICON_OMP_RUNTIME_SCHEDULE
        DO jb = i_startblk,i_endblk

          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)
            DO jc = i_startidx, i_endidx
              !horizontal tendency
              tend_hor           =                                                     &
                nabla2_e(ieidx(jc,jb,1),jk,ieblk(jc,jb,1))*p_int%geofac_div(jc,1,jb) + &
                nabla2_e(ieidx(jc,jb,2),jk,ieblk(jc,jb,2))*p_int%geofac_div(jc,2,jb) + &
                nabla2_e(ieidx(jc,jb,3),jk,ieblk(jc,jb,3))*p_int%geofac_div(jc,3,jb)

              !vertical tendency
              flux_up = 0._wp

              flux_dn = diff_smag_ic(jc,jk+1,jb) * (var(jc,jk,jb) - var(jc,jk+1,jb)) / &
                        p_nh_metrics%ddqz_z_half(jc,jk+1,jb) 

              tend_ver = (flux_up - flux_dn) * p_nh_metrics%inv_ddqz_z_full(jc,jk,jb) 
             
              !Advance scalar
              tot_tend(jc,jk,jb) = (tend_hor + tend_ver) * fac(jc,jk,jb)
              !var(jc,jk,jb) = var(jc,jk,jb) + dtime * tot_tend(jc,jk,jb)
            ENDDO
        ENDDO
!$OMP END DO

        !--------------------------------------------------------    
        jk = nlev
        !--------------------------------------------------------    

!$OMP DO PRIVATE(jc,jb,i_startidx,i_endidx,tend_hor,flux_up,flux_dn,tend_ver),ICON_OMP_RUNTIME_SCHEDULE
        DO jb = i_startblk,i_endblk

          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)
            DO jc = i_startidx, i_endidx
              !horizontal tendency
              tend_hor           =                                                     &
                nabla2_e(ieidx(jc,jb,1),jk,ieblk(jc,jb,1))*p_int%geofac_div(jc,1,jb) + &
                nabla2_e(ieidx(jc,jb,2),jk,ieblk(jc,jb,2))*p_int%geofac_div(jc,2,jb) + &
                nabla2_e(ieidx(jc,jb,3),jk,ieblk(jc,jb,3))*p_int%geofac_div(jc,3,jb)

              !vertical tendency
              flux_up = diff_smag_ic(jc,jk,jb) * (var(jc,jk-1,jb) - var(jc,jk,jb)) / &
                        p_nh_metrics%ddqz_z_half(jc,jk,jb) 

              flux_dn  = -sflux(jc,jb)

              tend_ver = (flux_up - flux_dn) * p_nh_metrics%inv_ddqz_z_full(jc,jk,jb) 
 
              !Advance scalar
              tot_tend(jc,jk,jb) = (tend_hor + tend_ver) * fac(jc,jk,jb)
              !var(jc,jk,jb) = var(jc,jk,jb) + dtime * tot_tend(jc,jk,jb)
            
            ENDDO
        ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL


  END SUBROUTINE diffuse_scalar

!---------------------------------------------------------------------------------

END MODULE mo_sgs_turbulence
