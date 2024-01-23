! Interface between NWP physics and the hydrological discharge model, through a coupler.
! Based on the ocean coupling interface
!
! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------
!
!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_nwp_hydrodisc_coupling

  USE mo_kind                ,ONLY: wp
  USE mo_model_domain        ,ONLY: t_patch
  USE mo_ext_data_state      ,ONLY: ext_data
  USE mo_nwp_lnd_types       ,ONLY: t_lnd_diag
  USE mo_nwp_phy_types       ,ONLY: t_nwp_phy_diag
  USE mo_lnd_nwp_config      ,ONLY: ntiles_total, isub_lake
  USE mo_ext_data_types      ,ONLY: t_external_data
  USE mo_fortran_tools       ,ONLY: init
  USE mo_parallel_config     ,ONLY: nproma
  USE mo_atm_phy_nwp_config  ,ONLY: atm_phy_nwp_config
  USE mo_impl_constants      ,ONLY: min_rlcell, LSS_TERRA, SUCCESS
  USE mo_loopindices         ,ONLY: get_indices_c

  USE mo_run_config          ,ONLY: ltimer, dtime
  USE mo_timer               ,ONLY: timer_start, timer_stop, timer_coupling_put

  USE mo_coupling_config     ,ONLY: is_coupled_to_hydrodisc
  USE mo_coupling_utils      ,ONLY: def_field
#ifdef YAC_coupling
  USE mo_yac_finterface      ,ONLY: yac_fdef_mask, yac_fput, yac_dble_ptr, &
                                    YAC_ACTION_COUPLING, YAC_ACTION_OUT_OF_BOUND, &
                                    YAC_LOCATION_CELL
#endif

  USE mo_exception           ,ONLY: warning, message, finish

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_nwp_hydrodisc_coupling, nwp_couple_hydrodisc

  CHARACTER(len=*), PARAMETER :: str_module = 'mo_nwp_hydrodisc_coupling' ! Output of module for debug

  INTEGER :: field_id_runoffs, field_id_runoffg

CONTAINS

  SUBROUTINE construct_nwp_hydrodisc_coupling( &
    p_patch, comp_id, grid_id, cell_point_id, timestepstring)

    TYPE(t_patch), TARGET, INTENT(IN) :: p_patch(:)
    INTEGER, INTENT(IN) :: comp_id
    INTEGER, INTENT(IN) :: grid_id
    INTEGER, INTENT(IN) :: cell_point_id
    CHARACTER(LEN=*), INTENT(IN) :: timestepstring

    TYPE(t_patch), POINTER :: patch_horz

    INTEGER :: cell_mask_id

    INTEGER :: jg, jb, jc, error

    LOGICAL,  ALLOCATABLE :: is_valid(:)

    CHARACTER(LEN=*), PARAMETER   :: routine = str_module // ':construct_nwp_hydrodisc_coupling'

#ifndef YAC_coupling
    CALL finish(routine, 'built without coupling support.')
#else

    jg = 1
    patch_horz => p_patch(jg)

    ALLOCATE(is_valid(nproma*patch_horz%nblks_c), STAT = error)
    IF(error /= SUCCESS) CALL finish(str_module, "memory allocation failure for is_valid")

    !ICON_OMP_PARALLEL PRIVATE(jb,jc)
      !ICON_OMP_WORKSHARE
      is_valid(:) = .FALSE.
      !ICON_OMP_END_WORKSHARE

      !ICON_OMP_DO ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1, patch_horz%nblks_c
        DO jc = 1, nproma
          IF ( ext_data(jg)%atm%fr_land(jc,jb)+ext_data(jg)%atm%fr_lake(jc,jb) .GE. 0.05_wp ) THEN
            is_valid((jb-1)*nproma+jc) = .TRUE.
          END IF
        END DO
      END DO
      !ICON_OMP_END_DO
    !ICON_OMP_END_PARALLEL

    CALL yac_fdef_mask (          &
      & grid_id,                  &
      & patch_horz%n_patch_cells, &
      & YAC_LOCATION_CELL,        &
      & is_valid,                 &
      & cell_mask_id )

    DEALLOCATE (is_valid, STAT = error)
    IF(error /= SUCCESS) CALL finish(str_module, "Deallocation failed for is_valid")

    CALL def_field( &
      comp_id, cell_point_id, cell_mask_id, timestepstring, &
      "surface_water_runoff", 1, field_id_runoffs)

    CALL def_field( &
      comp_id, cell_point_id, cell_mask_id, timestepstring, &
      "soil_water_runoff", 1, field_id_runoffg)
! YAC_coupling
#endif

  END SUBROUTINE construct_nwp_hydrodisc_coupling

  !>
  !! SUBROUTINE nwp_couple_hydrodisc -- the interface between
  !! NWP physics and the hydrodisc, through a coupler
  !!
  !! This subroutine is called from nwp_nh_interface.

  SUBROUTINE nwp_couple_hydrodisc( p_patch, lnd_diag, prm_diag, ext_data )

    ! Arguments

    TYPE(t_patch),   TARGET, INTENT(INOUT)  :: p_patch
    TYPE(t_lnd_diag),        INTENT(INOUT)  :: lnd_diag
    TYPE(t_nwp_phy_diag),    INTENT(INOUT)  :: prm_diag
    TYPE(t_external_data),   INTENT(INOUT)  :: ext_data

    ! Local variables

    LOGICAL               :: write_coupler_restart
    INTEGER               :: nbr_hor_cells         ! inner points
    INTEGER               :: jg                    ! patch ID
    INTEGER               :: jb                    ! block loop count
    INTEGER               :: jc                    ! nproma loop count
    INTEGER               :: error
    INTEGER               :: info, ierror          ! return values from cpl_put/get calls
    INTEGER               :: rl_start, rl_end
    INTEGER               :: i_startblk, i_endblk  ! blocks
    INTEGER               :: i_startidx, i_endidx  ! slices
    INTEGER               :: isubs                 ! tile index
    REAL(wp), TARGET, ALLOCATABLE :: buffer(:,:)   ! buffer transferred to YAC coupler
    CHARACTER(LEN=*), PARAMETER   :: routine = str_module // ':nwp_couple_hydrodisc'

#ifndef YAC_coupling
    CALL finish(routine, 'built without coupling support.')
#else

    TYPE(yac_dble_ptr)    :: ptr(1,1)

    ALLOCATE(buffer(nproma, p_patch%nblks_c), STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")

    jg          = p_patch%id

    ! include boundary interpolation zone of nested domains and halo points
    rl_start = 1
    rl_end   = min_rlcell

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

    nbr_hor_cells = p_patch%n_patch_cells

    !-------------------------------------------------------------------------
    ! If running in atm-hydrological discharge coupled mode, exchange information 
    !-------------------------------------------------------------------------
    !
    ! Possible fields that contain information to be sent to the hydrological
    ! discharge model include
    !
    ! 1. lnd_diag%runoff_s_inst_t(:,:,:)    averaged surface water runoff   [kg/m2/s] 
    ! 2. lnd_diag%runoff_g_inst_t(:,:,:)    averaged soil water runoff      [kg/m2/s] 
    !
    !-------------------------------------------------------------------------

    !------------------------------------------------
    !  Send surface water runoff
    !    "surface water runoff"
    !------------------------------------------------

    !$OMP PARALLEL
    CALL init(buffer(:,:))
    !$OMP END PARALLEL

!ICON_OMP_PARALLEL_DO PRIVATE(jb, jc, i_startidx, i_endidx, isubs) ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                       & i_startidx, i_endidx, rl_start, rl_end)
      IF ( atm_phy_nwp_config(jg)%inwp_surface == LSS_TERRA ) THEN

        ! aggregate over tiles

        DO isubs = 1, ntiles_total
          IF ( isubs == isub_lake ) THEN
            DO jc = i_startidx, i_endidx
              ! Take P-E over the lake as runoff
              buffer(jc,jb) = buffer(jc,jb) + ( prm_diag%tot_prec_rate(jc,jb) + prm_diag%qhfl_s_t(jc,jb,isubs) )  &
                &           * ext_data%atm%frac_t(jc,jb,isubs)
            ENDDO
          ELSE
            DO jc = i_startidx, i_endidx
              buffer(jc,jb) = buffer(jc,jb) + lnd_diag%runoff_s_inst_t(jc,jb,isubs) / dtime                       &
                &           * ext_data%atm%frac_t(jc,jb,isubs)
            ENDDO
          ENDIF
        ENDDO  ! isubs
        
      ELSE ! JSBACH

        DO jc = i_startidx, i_endidx
          buffer(jc,jb) = lnd_diag%runoff_s_inst_t(jc,jb,1) / dtime
        ENDDO

      ENDIF
    ENDDO ! jb

!ICON_OMP_END_PARALLEL_DO

    ptr(1,1)%p(1:nbr_hor_cells) => buffer

    IF (ltimer) CALL timer_start(timer_coupling_put)
    CALL yac_fput ( field_id_runoffs, SIZE(ptr(:,1:1),1), 1, ptr(:,1:1), info, ierror )
    IF ( info > YAC_ACTION_COUPLING .AND. info < YAC_ACTION_OUT_OF_BOUND ) THEN
      write_coupler_restart = .TRUE.
    ELSE
      write_coupler_restart = .FALSE.
    ENDIF

    IF ( info == YAC_ACTION_OUT_OF_BOUND ) &
         & CALL warning('nwp_couple_hydrodisc', 'YAC says fput called after end of run - id=1, surface water runoff')

    IF (ltimer) CALL timer_stop(timer_coupling_put)

    IF ( write_coupler_restart ) THEN
      CALL message('nwp_couple_hydrodisc', 'YAC says it is put for restart - ids 1, surface water runoff')
    ENDIF

    !------------------------------------------------
    !  Send soil water runoff
    !    "soil water runoff"
    !------------------------------------------------

    !$OMP PARALLEL
    CALL init(buffer(:,:))
    !$OMP END PARALLEL

!ICON_OMP_PARALLEL_DO PRIVATE(jb, jc, i_startidx, i_endidx, isubs) ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                       & i_startidx, i_endidx, rl_start, rl_end)

      IF ( atm_phy_nwp_config(jg)%inwp_surface == LSS_TERRA ) THEN

        ! aggregate over tiles
        DO isubs = 1, ntiles_total
          DO jc = i_startidx, i_endidx
            buffer (jc,jb) = buffer(jc,jb) + lnd_diag%runoff_g_inst_t(jc,jb,isubs) / dtime  &
                &            * ext_data%atm%frac_t(jc,jb,isubs)
          ENDDO
        ENDDO  ! isubs
        
      ELSE ! JSBACH

        DO jc = i_startidx, i_endidx
          buffer(jc,jb) = lnd_diag%runoff_g_inst_t(jc,jb,1) / dtime
        ENDDO

      ENDIF
    ENDDO ! jb
!ICON_OMP_END_PARALLEL_DO

    ptr(1,1)%p(1:nbr_hor_cells) => buffer

    IF (ltimer) CALL timer_start(timer_coupling_put)
    CALL yac_fput ( field_id_runoffg, SIZE(ptr(:,1:1),1), 1, ptr(:,1:1), info, ierror )
    IF ( info > YAC_ACTION_COUPLING .AND. info < YAC_ACTION_OUT_OF_BOUND ) THEN
      write_coupler_restart = .TRUE.
    ELSE
      write_coupler_restart = .FALSE.
    ENDIF

    IF ( info == YAC_ACTION_OUT_OF_BOUND ) &
         & CALL warning('nwp_couple_hydrodisc', 'YAC says fput called after end of run - id=2, ground water runoff')

    IF (ltimer) CALL timer_stop(timer_coupling_put)

    IF ( write_coupler_restart ) THEN
      CALL message('nwp_couple_hydrodisc', 'YAC says it is put for restart - ids 2, ground water runoff')
    ENDIF

    DEALLOCATE(buffer, STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")

!YAC_coupling
#endif

  END SUBROUTINE nwp_couple_hydrodisc

END MODULE mo_nwp_hydrodisc_coupling
