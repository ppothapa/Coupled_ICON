!>
!! Core subroutine used for SPPT
!! (Stochastic Perturbation of Physics Tendencies)
!!
!! @author Sascha Bellaire, MCH
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
!<

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_sppt_core

  USE mtime,                      ONLY: datetime, timedelta, OPERATOR(-), getTotalSecondsTimedelta
  USE mo_kind,                    ONLY: wp
  USE mo_model_domain,            ONLY: t_patch
  USE mo_impl_constants,          ONLY: min_rlcell_int
  USE mo_impl_constants_grf,      ONLY: grf_bdywidth_c
  USE mo_loopindices,             ONLY: get_indices_c
  USE mo_sppt_config,             ONLY: t_sppt_config
  USE mo_nonhydro_types,          ONLY: t_nh_prog
  USE mo_nonhydrostatic_config,   ONLY: kstart_moist
  USE mo_atm_phy_nwp_config,      ONLY: atm_phy_nwp_config
  USE mo_sppt_types,              ONLY: t_sppt
  USE mo_nwp_phy_types,           ONLY: t_nwp_phy_tend
  USE mo_run_config,              ONLY: iqv, iqc, iqi, iqr, iqs, iqg
  USE mo_fortran_tools,           ONLY: negative2zero

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: time_interpol_rn ! time interpolation/tapering of random numbers
  PUBLIC :: calc_tend        ! calculate tendencies
  PUBLIC :: pert_tend        ! perturb tendencies
  PUBLIC :: apply_tend       ! apply tendencies
  PUBLIC :: save_state       ! save prognostic state for key variables

  CONTAINS


  !>--------------------------------------------------------------------
  !! Time interpolation and tapering of random numbers
  !!
  !! @par Revision History
  !!
  !!<-------------------------------------------------------------------

  SUBROUTINE time_interpol_rn (p_patch, sppt_config, rn_2d_now, rn_2d_new, &
    &                          mtime_datetime, rn_3d)

    TYPE(t_patch),                INTENT(IN)        :: p_patch          !< patch variables
    TYPE(t_sppt_config),          INTENT(IN)        :: sppt_config
    REAL(wp),                     INTENT(IN)        :: rn_2d_now(:,:)   !< 2d utility field
    REAL(wp),                     INTENT(IN)        :: rn_2d_new(:,:)   !< 2d utility field

    TYPE(datetime),  POINTER,     INTENT(IN)        :: mtime_datetime   !< Date/time information
    REAL(wp),                     INTENT(INOUT)     :: rn_3d(:,:,:)     !< 3d field of random numbers

    ! Local variables
    INTEGER  :: jk, jc, jb

    REAL(wp) :: int_fact_1                    ! time interpolation weight
    Real(wp) :: int_fact_2                    ! time interpolation weight

    TYPE(timedelta) :: td                     !< timedelta elapsed since last rapa event
    REAL(wp) :: elapsed_time                  !< timedelta converted to seconds

    LOGICAL, PARAMETER  :: ltaper = .TRUE.    !< if TRUE tapering is performed
    LOGICAL, PARAMETER  :: linter = .TRUE.    !< if TRUE interpolation is performed

    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_endidx, i_startidx

    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

    !<---------------------------------------------------------------------
    ! Interpolation
    !>---------------------------------------------------------------------

    IF(linter) THEN

      ! Calculate interpolation factors

      ! get timedelta elapsed since last rapa event
      td = sppt_config%validity_date_rn_2d_new - mtime_datetime
      ! get elapsed time in seconds
      elapsed_time = getTotalSecondsTimeDelta(td, mtime_datetime)

      int_fact_1 = elapsed_time / (sppt_config%hinc_rn*3600.0_wp)  ! assumes that hinc_rn is given in hours
      int_fact_2 = 1.0_wp - int_fact_1

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,  &
                           i_startidx, i_endidx, rl_start, rl_end)

        ! Interpolate random numbers in time
        !
        DO jk = 1, p_patch%nlev
          DO jc = i_startidx, i_endidx

            rn_3d(jc,jk,jb)  = int_fact_1*rn_2d_now(jc,jb) + int_fact_2*rn_2d_new(jc,jb)

          ENDDO
        ENDDO

        !<---------------------------------------------------------------------
        ! Tapering
        !>---------------------------------------------------------------------

        IF(ltaper) THEN

          ! Apply tapering factor
          DO jk = 1, p_patch%nlev
            DO jc = i_startidx, i_endidx

              rn_3d(jc,jk,jb)  = rn_3d(jc,jk,jb) * sppt_config%taper(jk)

            ENDDO ! end of jc
          ENDDO ! end of jk

        ENDIF ! end of ltaper

      ENDDO  ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ENDIF ! end of linter


  END SUBROUTINE time_interpol_rn



  !>--------------------------------------------------------------------
  !! Calculate tendencies
  !!
  !! @par Revision History
  !!
  !!<-------------------------------------------------------------------

  SUBROUTINE calc_tend (i_startidx, i_endidx, nlev, ddt, var_new, var_now, dt)

    INTEGER,                  INTENT(IN)        :: i_startidx
    INTEGER,                  INTENT(IN)        :: i_endidx
    INTEGER,                  INTENT(IN)        :: nlev

    REAL(wp),                 INTENT(INOUT)     :: ddt(:,:)       !< tendencies to be calculated
    REAL(wp),                 INTENT(IN)        :: var_new(:,:)   !< current state of variable of interest
    REAL(wp),                 INTENT(IN)        :: var_now(:,:)   !< current state of variable of interest

    REAL(wp),                 INTENT(IN)        :: dt             !< (advective) time step applicable to local grid level

    ! Local variables
    INTEGER  :: jk,jc

    !<---------------------------------------------------------------------
    ! Calculate tendencies
    !>---------------------------------------------------------------------

    DO jk = 1, nlev
      DO jc = i_startidx, i_endidx

        ddt(jc,jk) = (var_new(jc,jk)  - var_now(jc,jk)) / dt

      ENDDO
    ENDDO

  END SUBROUTINE calc_tend


  !>--------------------------------------------------------------------
  !! Perturb tendencies
  !!
  !! @par Revision History
  !!
  !!<-------------------------------------------------------------------

  SUBROUTINE pert_tend(jb, jg, i_startidx, i_endidx, nlev, sppt,  &
                        prm_nwp_tend, rho_atm,                &
                        ddt_temp, ddt_u_tot, ddt_v_tot) 

    INTEGER,                     INTENT(IN)        :: jb
    INTEGER,                     INTENT(IN)        :: jg             !< patch ID
    INTEGER,                     INTENT(IN)        :: i_startidx
    INTEGER,                     INTENT(IN)        :: i_endidx
    INTEGER,                     INTENT(IN)        :: nlev 

    TYPE(t_sppt),                INTENT(INOUT)     :: sppt 
    TYPE(t_nwp_phy_tend),        INTENT(IN)        :: prm_nwp_tend

    REAL(wp),                    INTENT(IN)        :: rho_atm(:,:)

    REAL(wp),                    INTENT(INOUT)     :: ddt_temp(:,:)
    REAL(wp),                    INTENT(INOUT)     :: ddt_u_tot(:,:)
    REAL(wp),                    INTENT(INOUT)     :: ddt_v_tot(:,:)

    ! Local variables
    INTEGER  :: jk, jc

    !<---------------------------------------------------------------------
    ! Perturb tendencies
    !>---------------------------------------------------------------------

    DO jk = 1, nlev
      DO jc = i_startidx, i_endidx

        ! air temperature
        ddt_temp(jc,jk)  = (1._wp + sppt%rn_3d(jc,jk,jb)) *  ddt_temp(jc,jk)    &
          &              + sppt%ddt_temp_fast(jc,jk,jb)   *  sppt%rn_3d(jc,jk,jb)

        ! wind u- and v-components
        ddt_u_tot(jc,jk)  = (1._wp + sppt%rn_3d(jc,jk,jb)) * ddt_u_tot(jc,jk)   &
          &               + sppt%ddt_u_fast(jc,jk,jb) *  sppt%rn_3d(jc,jk,jb)

        ddt_v_tot(jc,jk)  = (1._wp + sppt%rn_3d(jc,jk,jb)) * ddt_v_tot(jc,jk)   &
          &               + sppt%ddt_v_fast(jc,jk,jb) *  sppt%rn_3d(jc,jk,jb)


        ! tracer -  Special treatment for water tracers here only the perturbed
        ! tendencies are used to avoid double-counting.

        ! water vapor
        sppt%ddt_qv(jc,jk,jb) = ( ( prm_nwp_tend%ddt_tracer_pconv(jc,jk,jb,iqv)/rho_atm(jc,jk) ) + sppt%ddt_qv_fast(jc,jk,jb) ) &
          &                    * sppt%rn_3d(jc,jk,jb)

        ! cloud water
        sppt%ddt_qc(jc,jk,jb) = ( ( prm_nwp_tend%ddt_tracer_pconv(jc,jk,jb,iqc)/rho_atm(jc,jk) ) + sppt%ddt_qc_fast(jc,jk,jb) ) &
          &                    * sppt%rn_3d(jc,jk,jb)

        ! cloud ice
        sppt%ddt_qi(jc,jk,jb) = ( ( prm_nwp_tend%ddt_tracer_pconv(jc,jk,jb,iqi)/rho_atm(jc,jk) ) + sppt%ddt_qi_fast(jc,jk,jb) ) &
          &                    * sppt%rn_3d(jc,jk,jb)

        IF ( iqg /= 0 ) THEN
          ! graupel -  no contribution from the slow physics
          sppt%ddt_qg(jc,jk,jb)  =  sppt%ddt_qg_fast(jc,jk,jb)  * sppt%rn_3d(jc,jk,jb)
        ENDIF

      ENDDO  ! end of jc
    ENDDO ! end of jk


    ! note that convective tendencies for qr and qs exist only if ldetrain_conv_prec=.TRUE.
    IF (atm_phy_nwp_config(jg)%ldetrain_conv_prec) THEN

      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
          ! rain
          sppt%ddt_qr(jc,jk,jb) = ( ( prm_nwp_tend%ddt_tracer_pconv(jc,jk,jb,iqr)/rho_atm(jc,jk) ) + sppt%ddt_qr_fast(jc,jk,jb) ) &
            &                    * sppt%rn_3d(jc,jk,jb)

          ! snow
          sppt%ddt_qs(jc,jk,jb) = ( ( prm_nwp_tend%ddt_tracer_pconv(jc,jk,jb,iqs)/rho_atm(jc,jk) ) + sppt%ddt_qs_fast(jc,jk,jb) ) &
            &                    * sppt%rn_3d(jc,jk,jb)
        ENDDO  ! end of jc
      ENDDO ! end of jk

    ELSE

      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
          ! rain
          sppt%ddt_qr(jc,jk,jb) = sppt%ddt_qr_fast(jc,jk,jb) * sppt%rn_3d(jc,jk,jb)

          ! snow
          sppt%ddt_qs(jc,jk,jb) = sppt%ddt_qs_fast(jc,jk,jb) * sppt%rn_3d(jc,jk,jb)
        ENDDO  ! end of jc
      ENDDO ! end of jk

    ENDIF


  END SUBROUTINE pert_tend



  !>--------------------------------------------------------------------
  !! Apply tendencies
  !!
  !! @par Revision History
  !!
  !!<-------------------------------------------------------------------

  SUBROUTINE apply_tend(pt_patch, sppt, pt_prog_rcf)

    ! Subroutine Arguments ((in/out/inout)
    TYPE(t_patch),   INTENT(IN)        :: pt_patch       !< grid/patch info

    TYPE(t_sppt),    INTENT(IN)        :: sppt

    TYPE(t_nh_prog), INTENT(INOUT)     :: pt_prog_rcf


    ! Local variables
    INTEGER  :: jb, jk, jc

    INTEGER  :: nlev
    INTEGER  :: rl_start, rl_end
    INTEGER  :: i_startblk, i_endblk
    INTEGER  :: i_startidx, i_endidx

    ! Number of vertical levels
    nlev = pt_patch%nlev

    ! Loop boundaries for prognostic domain.
    rl_start   = grf_bdywidth_c + 1
    rl_end     = min_rlcell_int

    i_startblk = pt_patch%cells%start_block(rl_start)
    i_endblk   = pt_patch%cells%end_block(rl_end)

    !<---------------------------------------------------------------------
    ! Calculate tendencies
    !>---------------------------------------------------------------------

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk,  &
        &                i_startidx, i_endidx, rl_start, rl_end)


      DO jk = kstart_moist(pt_patch%id), nlev
        DO jc = i_startidx, i_endidx

          pt_prog_rcf%tracer(jc,jk,jb,iqv) = sppt%ddt_qv(jc,jk,jb) + pt_prog_rcf%tracer (jc,jk,jb,iqv) ! water vapor
          pt_prog_rcf%tracer(jc,jk,jb,iqi) = sppt%ddt_qi(jc,jk,jb) + pt_prog_rcf%tracer (jc,jk,jb,iqi) ! ice
          pt_prog_rcf%tracer(jc,jk,jb,iqr) = sppt%ddt_qr(jc,jk,jb) + pt_prog_rcf%tracer (jc,jk,jb,iqr) ! rain
          pt_prog_rcf%tracer(jc,jk,jb,iqs) = sppt%ddt_qs(jc,jk,jb) + pt_prog_rcf%tracer (jc,jk,jb,iqs) ! snow
          pt_prog_rcf%tracer(jc,jk,jb,iqc) = sppt%ddt_qc(jc,jk,jb) + pt_prog_rcf%tracer (jc,jk,jb,iqc) ! cloud ice

          IF ( iqg /= 0 ) THEN
            pt_prog_rcf%tracer(jc,jk,jb,iqg) = sppt%ddt_qg(jc,jk,jb) + pt_prog_rcf%tracer (jc,jk,jb,iqg) ! graupel
          ENDIF

        ENDDO ! jc
      ENDDO ! jk
    ENDDO ! jb
!$OMP END DO

    ! perform clipping of negative tracer concentrations
    CALL negative2zero(pt_prog_rcf%tracer(:,:,:,:))
!$OMP END PARALLEL

  END SUBROUTINE apply_tend

  !>--------------------------------------------------------------------
  !! Save prognostic state for key variables
  !!
  !! @par Revision History
  !!
  !!<-------------------------------------------------------------------

  SUBROUTINE save_state(jb, i_startidx, i_endidx, nlev,        &
                         temp, tracer, sppt)

    INTEGER,             INTENT(IN)        :: jb
    INTEGER,             INTENT(IN)        :: i_startidx
    INTEGER,             INTENT(IN)        :: i_endidx
    INTEGER,             INTENT(IN)        :: nlev

    REAL(wp),            INTENT(IN)        :: temp(:,:,:)

    REAL(wp),            INTENT(IN)        :: tracer(:,:,:,:)

    TYPE(t_sppt),        INTENT(INOUT)     :: sppt

    ! Local variables
    INTEGER  :: jk, jc

    !<---------------------------------------------------------------------
    ! Save state
    !>---------------------------------------------------------------------

    DO jk = 1, nlev
      DO jc = i_startidx, i_endidx

        sppt%temp_now(jc,jk,jb) = temp(jc,jk,jb)

        sppt%qv_now(jc,jk,jb) = tracer(jc,jk,jb,iqv)
        sppt%qi_now(jc,jk,jb) = tracer(jc,jk,jb,iqi)
        sppt%qr_now(jc,jk,jb) = tracer(jc,jk,jb,iqr)
        sppt%qs_now(jc,jk,jb) = tracer(jc,jk,jb,iqs)
        sppt%qc_now(jc,jk,jb) = tracer(jc,jk,jb,iqc)

        IF ( iqg /= 0 ) THEN
          sppt%qg_now(jc,jk,jb) = tracer(jc,jk,jb,iqg)
        ENDIF

      ENDDO
    ENDDO

  END SUBROUTINE save_state


END MODULE mo_sppt_core	








