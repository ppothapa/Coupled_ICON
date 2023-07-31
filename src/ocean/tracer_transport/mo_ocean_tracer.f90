!>
!! Contains the implementation of the tracer transport routines for the ICON ocean model.
!! This comprises advection and diffusion in horizontal and vertical direction.
!!
!!
!! @par Revision History
!!  Developed  by Peter Korn,       MPI-M (2011/01)
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
#include "icon_definitions.inc"
#include "iconfor_dsl_definitions.inc"
!----------------------------
MODULE mo_ocean_tracer
  !-------------------------------------------------------------------------
  USE mo_kind,                      ONLY: wp
  USE mo_ocean_nml,                 ONLY: n_zlev,                         &
    & threshold_min_t, threshold_max_t, threshold_min_s, threshold_max_s, &
    & iswm_oce,                 use_none,                                 &
    & l_with_vert_tracer_diffusion, l_with_vert_tracer_advection,         &
    & GMRedi_configuration,                                               &
    & Cartesian_Mixing, tracer_threshold_min, tracer_threshold_max,       &
    & tracer_update_mode, l_with_horz_tracer_diffusion,                   &
    & vert_mix_type,vmix_kpp
  USE mo_util_dbg_prnt,             ONLY: dbg_print
  USE mo_parallel_config,           ONLY: nproma
  USE mo_run_config,                ONLY: dtime, ltimer, debug_check_level
  USE mo_model_domain,              ONLY: t_patch, t_patch_3d
  USE mo_exception,                 ONLY: finish !, message_text, message
  USE mo_ocean_tracer_diffusion,    ONLY: tracer_diffusion_vertical_implicit
  USE mo_ocean_tracer_transport_horz, ONLY: advect_horz, diffuse_horz
  USE mo_ocean_tracer_transport_vert, ONLY: advect_flux_vertical
  USE mo_operator_ocean_coeff_3d,   ONLY: t_operator_coeff
  USE mo_grid_subset,               ONLY: t_subset_range, get_index_range
  USE mo_sync,                      ONLY: sync_c, sync_patch_array
  USE mo_timer,                     ONLY: timer_start, timer_stop, timers_level, timer_extra30
  USE mo_statistics,                ONLY: global_minmaxmean, print_value_location
  USE mo_ocean_tracer_transport_types,  ONLY: t_ocean_tracer, t_tracer_collection, t_ocean_transport_state

  IMPLICIT NONE

  PRIVATE

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_ocean_tracer'

  CHARACTER(LEN=12) :: str_module = 'oceTracer   '  ! Output of module for 1 line debug
  INTEGER :: idt_src    = 1               ! Level of detail for 1 line debug

  PUBLIC :: advect_ocean_tracers
!   PUBLIC :: prepare_tracer_transport

CONTAINS

  !-------------------------------------------------------------------------
  !>
  !! !  SUBROUTINE advects the tracers present in the ocean model.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2010).
  !!
!<Optimize:inUse>
  SUBROUTINE advect_ocean_tracers(old_tracers, new_tracers, transport_state, operators_coeff, use_acc)
    TYPE(t_tracer_collection), INTENT(inout)      :: old_tracers
    TYPE(t_tracer_collection), INTENT(inout)      :: new_tracers
    TYPE(t_ocean_transport_state), TARGET         :: transport_state
    TYPE(t_operator_coeff), INTENT(in) :: operators_coeff
    LOGICAL, INTENT(in), OPTIONAL :: use_acc

    !Local variables
    TYPE(t_patch_3d ), POINTER     :: patch_3d
    INTEGER :: tracer_index
    LOGICAL :: lacc
    !-------------------------------------------------------------------------------
    patch_3d => transport_state%patch_3d

    IF (PRESENT(use_acc)) THEN
      lacc = use_acc
    ELSE
      lacc = .FALSE.
    END IF

    DO tracer_index = 1, old_tracers%no_of_tracers
      IF ( old_tracers%tracer(tracer_index)%is_advected) THEN
        CALL advect_diffuse_individual_tracer( patch_3d,    &
          & old_tracers%tracer(tracer_index),               &
          & transport_state, operators_coeff,               &
          & new_tracers%tracer(tracer_index),               &
          & old_tracers%typeOfTracers,                      &
          & use_acc=lacc)
      ENDIF
    END DO
    
!     IF ( old_tracers%no_of_tracers > 2) THEN
!     
!       CALL dbg_print('h_old all'       , transport_state%h_old, "transport", 1,  patch_3d%p_patch_2d(1)%cells%all)
!       CALL dbg_print('h_new all'       , transport_state%h_new, "transport", 1,  patch_3d%p_patch_2d(1)%cells%all)
!       CALL dbg_print('w all'           , transport_state%w, "transport", 1,  patch_3d%p_patch_2d(1)%cells%all)
!       CALL dbg_print('vn all'          , transport_state%vn, "transport", 1,  patch_3d%p_patch_2d(1)%edges%all)
!       CALL dbg_print('flux all'        , transport_state%mass_flux_e, "transport", 1,  patch_3d%p_patch_2d(1)%edges%all)
! 
!       CALL dbg_print('hor diff all'    , old_tracers%tracer(1)%hor_diffusion_coeff, "transport", 1,  &
!         & patch_3d%p_patch_2d(1)%edges%all)
!       CALL dbg_print('ver diff all'    , old_tracers%tracer(1)%ver_diffusion_coeff, "transport", 1,  &
!         & patch_3d%p_patch_2d(1)%cells%all)
!       
!       CALL dbg_print('trac 1 old all'    , old_tracers%tracer(1)%concentration, "transport", 1, &
!         patch_3d%p_patch_2d(1)%cells%all)
!       CALL dbg_print('trac 1 new all'    , new_tracers%tracer(1)%concentration, "transport", 1, &
!         patch_3d%p_patch_2d(1)%cells%all)
! 
!     ENDIF


  END SUBROUTINE advect_ocean_tracers
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE copy_individual_tracer_ab(patch_3d, old_tracer, new_tracer, use_acc)

    TYPE(t_patch_3d ),TARGET, INTENT(inout)   :: patch_3d
    TYPE(t_ocean_tracer), TARGET :: old_tracer
    TYPE(t_ocean_tracer), TARGET :: new_tracer
    LOGICAL, INTENT(in), OPTIONAL :: use_acc

    INTEGER :: jc,level,jb
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp), POINTER :: trac_old(:,:,:), trac_new(:,:,:) ! temporary pointers to the concentration arrays
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER :: patch_2D
    INTEGER :: start_block, end_block
    LOGICAL :: lacc
    CHARACTER(len=*), PARAMETER :: routine = modname//':copy_individual_tracer_ab'
    ! CHARACTER(len=max_char_length), PARAMETER :: &
    !        & routine = ('mo_tracer_advection:advect_individual_tracer')
    !-------------------------------------------------------------------------------

#ifdef _OPENACC
    CALL finish(routine, 'OpenACC version currently not tested/validated')
#endif

    trac_old => old_tracer%concentration
    trac_new => new_tracer%concentration
    patch_2D => patch_3d%p_patch_2d(1)
    all_cells => patch_2D%cells%all
    start_block = all_cells%start_block
    end_block = all_cells%end_block

    IF (PRESENT(use_acc)) THEN
      lacc = use_acc
    ELSE
      lacc = .FALSE.
    END IF

    !$ACC DATA COPYIN(trac_old) &
    !$ACC   COPY(trac_new) IF(lacc)

    DO jb = start_block, end_block
      CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
      !$ACC LOOP GANG VECTOR
      DO jc = start_cell_index, end_cell_index
        DO level = 1, n_zlev
          trac_new(jc,level,jb) = trac_old(jc,level,jb)
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
    !$ACC WAIT(1)
    !$ACC END DATA
  END SUBROUTINE copy_individual_tracer_ab
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! !  SUBROUTINE advects the tracers present in the ocean model.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2010).
  !!
!<Optimize:inUse>
  SUBROUTINE advect_diffuse_individual_tracer(patch_3d, old_tracer,       &
    & transport_state, operators_coeff, new_tracer, typeOfTracers, use_acc)

    TYPE(t_patch_3d ),TARGET, INTENT(inout)   :: patch_3d
    TYPE(t_ocean_tracer), TARGET :: old_tracer
    TYPE(t_ocean_tracer), TARGET :: new_tracer
    TYPE(t_ocean_transport_state), TARGET :: transport_state
    TYPE(t_operator_coeff),INTENT(in) :: operators_coeff
    CHARACTER(LEN=*), INTENT(in)         :: typeOfTracers
    LOGICAL, INTENT(in), OPTIONAL :: use_acc
!     REAL(wp), INTENT(inout), OPTIONAL :: horizontally_diffused_tracer(:,:,:)

    !Local variables
    TYPE(t_subset_range), POINTER :: cells_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    LOGICAL :: lacc
    ! CHARACTER(len=max_char_length), PARAMETER :: &
    !        & routine = ('mo_tracer_advection:advect_diffuse_tracer')
    !-------------------------------------------------------------------------------_

    patch_2D => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2D%cells%in_domain

    IF (PRESENT(use_acc)) THEN
      lacc = use_acc
    ELSE
      lacc = .FALSE.
    END IF

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=2  ! output print level (1-5, fix)
    CALL dbg_print('on entry: IndTrac: trac_old',old_tracer%concentration(:,:,:), &
      & str_module,idt_src, in_subset=patch_2D%cells%owned)
    !---------------------------------------------------------------------
    IF(tracer_update_mode == use_none ) THEN
      CALL copy_individual_tracer_ab( patch_3d,            &
        & old_tracer,new_tracer,use_acc=lacc)
      RETURN
    ENDIF
   
    !Shallow water is done with horizontal advection
    IF(iswm_oce == 1) THEN
      CALL advect_diffuse_SW_tracer(patch_3d, old_tracer,       &
        & transport_state, operators_coeff,                      &
        & old_tracer%hor_diffusion_coeff,        &
        & new_tracer, use_acc=lacc)
        
    !The 3D-case
    ELSE ! IF( iswm_oce /= 1) THEN

      CALL advect_diffuse_tracer( patch_3d, &
        & old_tracer,                       &
        & transport_state,                  &
        & operators_coeff,                  &
        & old_tracer%hor_diffusion_coeff,   &
        & old_tracer%ver_diffusion_coeff,   &
        & new_tracer,                       &
        & typeOfTracers,                    &
        & use_acc=lacc)

    ENDIF

  END SUBROUTINE advect_diffuse_individual_tracer
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! !  SUBROUTINE advects the tracers for shallow water case
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2010).
  !!
!<Optimize:inUse>
  SUBROUTINE advect_diffuse_SW_tracer(patch_3d, old_tracer,       &
    & transport_state, operators_coeff,                      &
    & k_h,                                   &
    & new_tracer, use_acc)

    TYPE(t_patch_3d ),TARGET, INTENT(inout)   :: patch_3d
    TYPE(t_ocean_tracer), TARGET :: old_tracer
    TYPE(t_ocean_transport_state), TARGET :: transport_state
    TYPE(t_operator_coeff),INTENT(in) :: operators_coeff
    REAL(wp), INTENT(in)                 :: k_h(:,:,:)       !horizontal mixing coeff
    TYPE(t_ocean_tracer), TARGET :: new_tracer
    LOGICAL, INTENT(in), OPTIONAL :: use_acc
 
    !Local variables
    REAL(wp) :: delta_t, delta_z,delta_z_new
    REAL(wp) :: div_adv_flux_horz(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: div_adv_flux_vert(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)    
    REAL(wp) :: div_diff_flux_horz(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: flux_vert(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: div_diff_flx(nproma, n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)

    REAL(wp), POINTER :: trac_old(:,:,:), trac_new(:,:,:) ! temporary pointers to the concentration arrays

    INTEGER :: jc,level,jb,ic,ib,alloc_cell_blocks,start_block,end_block
    INTEGER :: z_dolic
    INTEGER :: start_cell_index, end_cell_index
    LOGICAL :: lacc
    TYPE(t_subset_range), POINTER :: cells_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    REAL(wp), DIMENSION(:,:,:), POINTER :: new_tracer_concentration, old_tracer_concentration
    REAL(wp), DIMENSION(:,:), POINTER :: transport_state_h_old, transport_state_h_new
    CHARACTER(len=*), PARAMETER :: routine = modname//':advect_diffuse_SW_tracer'
    ! CHARACTER(len=max_char_length), PARAMETER :: &
    !        & routine = ('mo_tracer_advection:advect_diffuse_tracer')
    !-------------------------------------------------------------------------------

#ifdef _OPENACC
    CALL finish(routine, 'OpenACC version currently not tested/validated')
#endif

    trac_old => old_tracer%concentration
    trac_new => new_tracer%concentration

    patch_2D => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2D%cells%in_domain
    old_tracer_concentration => old_tracer%concentration
    new_tracer_concentration => new_tracer%concentration
    transport_state_h_old => transport_state%h_old
    transport_state_h_new => transport_state%h_new

    delta_t = dtime
    alloc_cell_blocks = patch_2D%alloc_cell_blocks
    start_block = cells_in_domain%start_block
    end_block = cells_in_domain%end_block

    IF (PRESENT(use_acc)) THEN
      lacc = use_acc
    ELSE
      lacc = .FALSE.
    END IF

    !$ACC DATA PRESENT(nproma, n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks) &
    !$ACC   COPYIN(k_h, transport_state_h_old, transport_state_h_new) &
    !$ACC   COPY(old_tracer_concentration) &
    !$ACC   CREATE(div_adv_flux_horz, div_adv_flux_vert, div_diff_flux_horz) IF(lacc)

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=2  ! output print level (1-5, fix)
    CALL dbg_print('on entry: IndTrac: trac_old',trac_old(:,:,:) ,str_module,idt_src, in_subset=patch_2D%cells%owned)
    !---------------------------------------------------------------------
    
    !Shallow water is done with horizontal advection
    !$ACC KERNELS DEFAULT(PRESENT) IF(lacc)
    div_adv_flux_horz   (1:nproma,1:n_zlev,1:alloc_cell_blocks) = 0.0_wp
    div_adv_flux_vert   (1:nproma,1:n_zlev,1:alloc_cell_blocks) = 0.0_wp
    ! div_diff_flux_horz  (1:nproma,1:n_zlev,1:alloc_cell_blocks) = 0.0_wp
    !$ACC END KERNELS

    !---------------------------------------------------------------------
    CALL advect_horz( patch_3d,        &
      & old_tracer_concentration,      &
      & transport_state,               &
      & operators_coeff,               &
      & k_h,                           &
      & transport_state_h_old,         &
      & transport_state_h_new,         &
      & div_adv_flux_horz,             &
      & div_adv_flux_vert)

    IF ( l_with_horz_tracer_diffusion) THEN

      CALL diffuse_horz( patch_3d,       &
        & old_tracer_concentration,      &
        & transport_state,               &
        & operators_coeff,               &
        & k_h,                           &
        & transport_state_h_old,         &
        & transport_state_h_new,         &
        & div_diff_flux_horz)
    ELSE
      div_diff_flux_horz  (1:nproma,1:n_zlev,1:alloc_cell_blocks) = 0.0_wp
    ENDIF

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('DivAdvFlx',div_adv_flux_horz,str_module,idt_src, in_subset=cells_in_domain)
    CALL dbg_print('DivDiffFlx',div_diff_flux_horz,str_module,idt_src, in_subset=cells_in_domain)
    !---------------------------------------------------------------------

    !$ACC DATA COPYIN(div_adv_flux_horz, div_diff_flux_horz) &
    !$ACC   COPYIN(old_tracer_concentration, transport_state_h_new) &
    !$ACC   COPY(new_tracer_concentration) IF(lacc)

    !level=1
    DO jb = start_block, end_block
      CALL get_index_range(cells_in_domain, jb, start_cell_index, end_cell_index)
      
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
      !$ACC LOOP GANG VECTOR
      DO jc = start_cell_index, end_cell_index
        DO level = 1, 1!MIN(dolic_c(jc,jb),1)  ! this at most should be 1

          !delta_z =  del_zlev_m(1)
          delta_z = transport_state_h_new(jc,jb) !- bathymetry_c(jc,jb)
          new_tracer_concentration(jc,level,jb)= old_tracer_concentration(jc,level,jb) - &
            & (delta_t/delta_z) * (div_adv_flux_horz(jc,level,jb)-div_diff_flux_horz(jc,level,jb))

        END DO
      END DO
      !$ACC END PARALLEL
    END DO

    CALL sync_patch_array(sync_c, patch_2D, new_tracer_concentration)
    !$ACC WAIT(1)
    !$ACC END DATA

    !$ACC END DATA
  END SUBROUTINE advect_diffuse_SW_tracer
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! !  SUBROUTINE advects the tracers present in the ocean model.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2010).
  !!
!<Optimize:inUse>
  SUBROUTINE advect_diffuse_tracer(        &
    & patch_3d, old_tracer,                &
    & transport_state, operators_coeff,    &
    & k_h, a_v,                            &
    & new_tracer, typeOfTracers, use_acc)!,        &
    ! & horizontally_diffused_tracer        )

    TYPE(t_patch_3d ),TARGET, INTENT(inout)   :: patch_3d
    TYPE(t_ocean_tracer), TARGET :: old_tracer
    TYPE(t_ocean_transport_state), TARGET :: transport_state
    TYPE(t_operator_coeff),INTENT(in) :: operators_coeff
    REAL(wp), INTENT(in)                 :: k_h(:,:,:)       !horizontal mixing coeff
    REAL(wp), INTENT(inout)              :: a_v(:,:,:)       !vertical mixing coeff, in
    TYPE(t_ocean_tracer), TARGET         :: new_tracer
    CHARACTER(LEN=*), INTENT(in)         :: typeOfTracers
    LOGICAL, INTENT(in), OPTIONAL        :: use_acc
!     REAL(wp), INTENT(inout), OPTIONAL :: horizontally_diffused_tracer(:,:,:)

    !Local variables
    REAL(wp) :: delta_t, delta_z,delta_z_new
    REAL(wp) :: div_adv_flux_horz(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: div_diff_flux_horz(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: flux_horz(nproma,n_zlev, patch_3d%p_patch_2D(1)%nblks_e)
    REAL(wp) :: div_adv_flux_vert(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    INTEGER :: jc,level,jb, je
    INTEGER :: z_dolic, max_dolic_c
    INTEGER :: start_cell_index, end_cell_index, start_block, end_block
    TYPE(t_subset_range), POINTER :: cells_in_domain, edges_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    REAL(wp) :: top_bc(nproma)
    LOGICAL :: lacc

    INTEGER, DIMENSION(:,:), POINTER :: dolic_c
    REAL(wp), DIMENSION(:,:,:), POINTER :: prism_thick_flat_sfc_c, new_tracer_concentration, &
      old_tracer_concentration, prism_thick_c, old_transport_tendencies
    REAL(wp), DIMENSION(:,:), POINTER :: transport_state_h_old, transport_state_h_new, old_tracer_top_bc
    CHARACTER(len=*), PARAMETER :: method_name = 'mo_ocean_tracer:advect_diffuse_tracer'
    !-------------------------------------------------------------------------------
    patch_2D        => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2D%cells%in_domain
    edges_in_domain => patch_2D%edges%in_domain
    delta_t = dtime
  
    dolic_c => patch_3d%p_patch_1D(1)%dolic_c
    transport_state_h_old => transport_state%h_old
    transport_state_h_new => transport_state%h_new
    old_tracer_top_bc => old_tracer%top_bc
    prism_thick_flat_sfc_c => patch_3d%p_patch_1D(1)%prism_thick_flat_sfc_c
    new_tracer_concentration => new_tracer%concentration
    old_tracer_concentration => old_tracer%concentration
    prism_thick_c => patch_3d%p_patch_1D(1)%prism_thick_c
    old_transport_tendencies => old_tracer%vertical_trasnport_tendencies
    start_block = cells_in_domain%start_block
    end_block = cells_in_domain%end_block

    IF (PRESENT(use_acc)) THEN
      lacc = use_acc
    ELSE
      lacc = .FALSE.
    END IF

    !$ACC DATA CREATE(div_adv_flux_horz, div_adv_flux_vert, div_diff_flux_horz, top_bc) IF(lacc)

    !---------------------------------------------------------------------
#ifndef _OPENACC
    !$ACC KERNELS DEFAULT(PRESENT) IF(lacc)
    ! these are probably not necessary
    div_adv_flux_vert = 0.0_wp
    div_adv_flux_horz = 0.0_wp
    div_diff_flux_horz = 0.0_wp
    !$ACC END KERNELS
#endif
    !---------------------------------------------------------------------
    IF ( l_with_vert_tracer_advection ) THEN

      CALL advect_flux_vertical( patch_3d,&
        & old_tracer_concentration, &
        & transport_state,                           &
        & operators_coeff,                     &
        & div_adv_flux_vert,        &
        & use_acc=lacc)

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src=3  ! output print level (1-5, fix)
      CALL dbg_print('aft. AdvFluxVert:divfluxvert',div_adv_flux_vert          ,str_module,idt_src, in_subset=cells_in_domain)
      !---------------------------------------------------------------------

    ENDIF  ! l_with_vert_tracer_advection

    !---------------------------------------------------------------------
    CALL advect_horz( patch_3d,       &
      & old_tracer_concentration, &
      & transport_state,                           &
      & operators_coeff,                     &
      & k_h,                            &
      & transport_state_h_old,         &
      & transport_state_h_new,         &
      & div_adv_flux_horz,              &
      & div_adv_flux_vert,             &
      & use_acc=lacc)
    !---------------------------------------------------------------------

    IF(GMRedi_configuration==Cartesian_Mixing)THEN
      !horizontal diffusion, vertical is handled implicitely below
      CALL diffuse_horz( patch_3d,      &
      & old_tracer_concentration, &
      & transport_state,                           &
      & operators_coeff,                     &
      & k_h,                            &
      & transport_state_h_old,         &
      & transport_state_h_new,         &
      & div_diff_flux_horz,            &
      & use_acc=lacc)

    ELSE
      CALL finish(method_name, "wrong GMredi call")
    ENDIF  
      
    !Case: Implicit Vertical diffusion

    !Calculate preliminary tracer value out of horizontal advective and
    !diffusive fluxes and vertical advective fluxes, plus surface forcing.
    !Surface forcing applied as volume forcing at rhs, i.e.part of explicit term
    !in tracer (and also momentum) eqs. In this case, top boundary condition of
    !vertical Laplacians are homogeneous

!ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index, end_cell_index, jc, level, &
!ICON_OMP delta_z, delta_z_new, top_bc) ICON_OMP_DEFAULT_SCHEDULE

    DO jb = start_block, end_block
      CALL get_index_range(cells_in_domain, jb, start_cell_index, end_cell_index)
      IF (ASSOCIATED(old_tracer_top_bc)) THEN
!         CALL dbg_print('top_bc all'       , old_tracer_top_bc, "transport", 1,  patch_3d%p_patch_2d(1)%cells%all)
!         CALL finish("ASSOCIATED(old_tracer_top_bc)","")
        !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lacc)
        top_bc(:) = old_tracer_top_bc(:,jb)
        !$ACC END KERNELS
      ELSE
        !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lacc)
        top_bc(:) = 0.0_wp
        !$ACC END KERNELS
      ENDIF
      
#ifdef __LVECTOR__
      level = 1
      max_dolic_c = -1
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) REDUCTION(MAX: max_dolic_c) IF(lacc)
      DO jc = start_cell_index, end_cell_index
        max_dolic_c = MAX(max_dolic_c, dolic_c(jc,jb))
      END DO
      !$ACC END PARALLEL LOOP

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
      !$ACC LOOP GANG VECTOR
      DO jc = start_cell_index, end_cell_index
        IF (dolic_c(jc,jb) < level) CYCLE
#else
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
      !$ACC LOOP GANG VECTOR
      DO jc = start_cell_index, end_cell_index
        !TODO check algorithm: inv_prism_thick_c vs. del_zlev_m | * vs. /
        DO level = 1, MIN(dolic_c(jc,jb),1)  ! this at most should be 1
#endif

          delta_z     = prism_thick_flat_sfc_c(jc,level,jb)+transport_state_h_old(jc,jb)
          delta_z_new = prism_thick_flat_sfc_c(jc,level,jb)+transport_state_h_new(jc,jb)

          new_tracer_concentration(jc,level,jb)= &
            & (old_tracer_concentration(jc,level,jb) * delta_z &
            & - delta_t * (&
            &  div_adv_flux_horz(jc,level,jb) +div_adv_flux_vert(jc,level,jb)&
            & -div_diff_flux_horz(jc,level,jb) - top_bc(jc))) / delta_z_new

#ifndef __LVECTOR__
        ENDDO
#endif
      ENDDO
      !$ACC END PARALLEL

      IF (vert_mix_type .EQ. vmix_kpp .and. typeOfTracers == "ocean") THEN
        !by_Oliver: account for nonlocal transport term for heat and scalar
        !(salinity) if KPP scheme is used
#ifdef __LVECTOR__
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO level = 2, max_dolic_c
          DO jc = start_cell_index, end_cell_index
            IF (dolic_c(jc,jb) < level) CYCLE
#else
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
        !$ACC LOOP GANG VECTOR
        DO jc = start_cell_index, end_cell_index
          DO level = 2, dolic_c(jc,jb)
#endif 
            new_tracer_concentration(jc,level,jb) =                          &
              &  old_tracer_concentration(jc,level,jb)                       &
              &  - (delta_t / prism_thick_c(jc,level,jb))    &
              &  * (  div_adv_flux_horz(jc,level,jb)  &
              &     + div_adv_flux_vert(jc,level,jb)  &
              &     - div_diff_flux_horz(jc,level,jb) &
              &     - old_transport_tendencies(jc,level,jb))
          END DO
        END DO
        !$ACC END PARALLEL

      ELSE

#ifdef __LVECTOR__
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO level = 2, max_dolic_c
          DO jc = start_cell_index, end_cell_index
            IF (dolic_c(jc,jb) < level) CYCLE
#else
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
        !$ACC LOOP GANG VECTOR
        DO jc = start_cell_index, end_cell_index
          DO level = 2, dolic_c(jc,jb)
#endif
            new_tracer_concentration(jc,level,jb) =                          &
              &  old_tracer_concentration(jc,level,jb)                       &
              &  - (delta_t / prism_thick_c(jc,level,jb))    &
              &  * (  div_adv_flux_horz(jc,level,jb)  &
              &     + div_adv_flux_vert(jc,level,jb)  &
              &     - div_diff_flux_horz(jc,level,jb) )
          END DO
        END DO
        !$ACC END PARALLEL

      END IF

    END DO
!ICON_OMP_END_PARALLEL_DO


    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('BefImplDiff: div_adv_flux_vert',div_adv_flux_vert, str_module,idt_src, in_subset=cells_in_domain)
    CALL dbg_print('BefImplDiff: trac_inter', new_tracer_concentration,  str_module,idt_src, in_subset=cells_in_domain)
    !---------------------------------------------------------------------

    !calculate vert diffusion impicit: result is stored in trac_out
    ! no sync because of columnwise computation
    IF ( l_with_vert_tracer_diffusion ) THEN
          
      !Vertical mixing: implicit and with coefficient a_v
      !that is the sum of PP-coeff and implicit part of Redi-scheme      
      CALL tracer_diffusion_vertical_implicit( &
          & patch_3d,                   &
          & new_tracer,                 &
          & a_v,                        &
          & transport_state_h_new,      &
          & use_acc=lacc)
          
    ENDIF!IF ( l_with_vert_tracer_diffusion )

    CALL sync_patch_array(sync_c, patch_2D, new_tracer_concentration)

    !$ACC WAIT(1)
    !$ACC END DATA

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('aft. AdvIndivTrac: trac_old', old_tracer_concentration, str_module, 3, in_subset=cells_in_domain)
    CALL dbg_print('aft. AdvIndivTrac: trac_new', new_tracer_concentration, str_module, 3, in_subset=cells_in_domain)
    !---------------------------------------------------------------------

  END SUBROUTINE advect_diffuse_tracer
  !-------------------------------------------------------------------------
 
  !-------------------------------------------------------------------------
  SUBROUTINE check_min_max_tracer(info_text, tracer, min_tracer, max_tracer, tracer_name, in_subset)
    CHARACTER(*) :: info_text
    REAL(wp), POINTER :: tracer(:,:,:)
    REAL(wp), INTENT(in) :: min_tracer, max_tracer
    CHARACTER(*) :: tracer_name
    TYPE(t_subset_range), POINTER :: in_subset

!     INTEGER  :: level
    REAL(wp) :: minmaxmean(3)
!     REAL(wp) :: lon, lat

      minmaxmean(:) = global_minmaxmean(values = tracer(:,:,:), in_subset=in_subset)
        IF (minmaxmean(1) < min_tracer) THEN
          WRITE(0,*) TRIM(tracer_name), ' too low:', minmaxmean(1)
          CALL print_value_location(tracer(:,:,:), minmaxmean(1), in_subset)
          CALL finish(TRIM(info_text), 'tracer below threshold')
        ENDIF

        IF (minmaxmean(2) > max_tracer) THEN
          WRITE(0,*) TRIM(tracer_name), ' too high:', minmaxmean(2)
          CALL print_value_location(tracer(:,:,:), minmaxmean(2), in_subset)
          CALL finish(TRIM(info_text), 'tracer above threshold')
        ENDIF

  END SUBROUTINE check_min_max_tracer
  !-------------------------------------------------------------------------

END MODULE mo_ocean_tracer


