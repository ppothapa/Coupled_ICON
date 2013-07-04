!>
!! @brief branch for the non-hydrostatic ICON workflow
!!
!! @author Kristina Froehlich, MPI-M (2011-07-19)
!!
!! @par Copyright
!! 2002-2011 by DWD and MPI-M
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
MODULE mo_atmo_nonhydrostatic

USE mo_kind,                 ONLY: wp
USE mo_exception,            ONLY: message, finish
USE mo_impl_constants,       ONLY: SUCCESS, max_dom, inwp
USE mo_mpi,                  ONLY: my_process_is_stdio
USE mo_timer,                ONLY: print_timer, timers_level, timer_start, &
  &                                timer_stop, timer_model_init
USE mo_master_control,       ONLY: is_restart_run
USE mo_output,               ONLY: init_output_files, close_output_files
USE mo_var_list,             ONLY: print_var_list
USE mo_time_config,          ONLY: time_config      ! variable
USE mo_io_restart,           ONLY: read_restart_files
USE mo_io_restart_attributes,ONLY: get_restart_attribute
USE mo_io_config,            ONLY: dt_file,dt_diag,dt_checkpoint
USE mo_parallel_config,      ONLY: nproma
USE mo_nh_pzlev_config,      ONLY: configure_nh_pzlev
USE mo_advection_config,     ONLY: configure_advection
USE mo_art_config,           ONLY: configure_art
USE mo_run_config,           ONLY: dtime, dtime_adv,     & !    namelist parameter
  &                                ltestcase,            &
  &                                nsteps,               & !    :
  &                                ltimer,               & !    :
  &                                iforcing,             & !    namelist parameter
  &                                output_mode,          &
  &                                msg_level,            & !    namelist parameter
  &                                lvert_nest, ntracer,  &
  &                                iqc, iqi, iqr, iqs
USE mo_dynamics_config,      ONLY: nnow, nnow_rcf, iequations
! Horizontal grid
USE mo_model_domain,         ONLY: p_patch
USE mo_grid_config,          ONLY: n_dom, start_time, is_plane_torus
! to break circular dependency KF???
USE mo_intp_data_strc,       ONLY: p_int_state
USE mo_grf_intp_data_strc,   ONLY: p_grf_state
! NH-namelist state
USE mo_nonhydrostatic_config,ONLY: iadv_rcf, kstart_moist, kend_qvsubstep, l_open_ubc

USE mo_atm_phy_nwp_config,   ONLY: configure_atm_phy_nwp, atm_phy_nwp_config
! NH-Model states
USE mo_nonhydro_state,       ONLY: p_nh_state, construct_nh_state, destruct_nh_state
USE mo_opt_diagnostics,      ONLY: construct_opt_diag, destruct_opt_diag
USE mo_nwp_phy_state,        ONLY: prm_diag, prm_nwp_diag_list, prm_nwp_tend_list, &
  &                                construct_nwp_phy_state,                        &
  &                                destruct_nwp_phy_state
USE mo_nwp_lnd_state,        ONLY: p_lnd_state, construct_nwp_lnd_state,       &
  &                                destruct_nwp_lnd_state
USE mo_nh_diagnose_pres_temp,ONLY: diagnose_pres_temp
! Time integration
USE mo_nh_stepping,          ONLY: prepare_nh_integration, perform_nh_stepping
! Initialization with real data
USE mo_nh_initicon,         ONLY: init_icon
USE mo_ext_data_state,      ONLY: ext_data, init_index_lists
! meteogram output
USE mo_meteogram_output,    ONLY: meteogram_init, meteogram_finalize
USE mo_meteogram_config,    ONLY: meteogram_output_config
USE mo_name_list_output_config,   ONLY: first_output_name_list, &
  &                               is_any_output_nml_active
USE mo_name_list_output_init, ONLY: init_name_list_output, &
  &                               parse_variable_groups
USE mo_name_list_output,    ONLY: close_name_list_output
USE mo_name_list_output_init, ONLY: output_file
USE mo_pp_scheduler,        ONLY: pp_scheduler_init, pp_scheduler_finalize
USE mo_intp_lonlat,         ONLY: compute_lonlat_area_weights

!-------------------------------------------------------------------------

IMPLICIT NONE
PRIVATE

PUBLIC :: atmo_nonhydrostatic
PUBLIC :: construct_atmo_nonhydrostatic, destruct_atmo_nonhydrostatic

! module data
INTEGER :: jfile, n_file, n_diag, n_chkpt
LOGICAL :: l_have_output

CONTAINS

  !---------------------------------------------------------------------
  SUBROUTINE atmo_nonhydrostatic
    
    CHARACTER(*), PARAMETER :: routine = "mo_atmo_nonhydrostatic"

    CALL construct_atmo_nonhydrostatic()
    
    !------------------------------------------------------------------
    ! Now start the time stepping:
    ! The special initial time step for the three time level schemes
    ! is executed within process_grid_level
    !------------------------------------------------------------------

    CALL perform_nh_stepping(  &
      &                       time_config%cur_datetime,       &
      &                       n_file, jfile, n_chkpt, n_diag, l_have_output  )
 
    IF (ltimer) CALL print_timer

    !---------------------------------------------------------------------
    ! 6. Integration finished. Clean up.
    !---------------------------------------------------------------------
    CALL destruct_atmo_nonhydrostatic()

    
  END SUBROUTINE atmo_nonhydrostatic
  !---------------------------------------------------------------------
  
  !---------------------------------------------------------------------
  SUBROUTINE construct_atmo_nonhydrostatic
    
    CHARACTER(*), PARAMETER :: routine = "construct_atmo_nonhydrostatic"
    

    INTEGER :: jg, ist
    LOGICAL :: l_realcase
    INTEGER :: pat_level(n_dom)
    LOGICAL :: l_pres_msl(n_dom) !< Flag. TRUE if computation of mean sea level pressure desired
    LOGICAL :: l_rh(n_dom)       !< Flag. TRUE if computation of relative humidity desired

    IF (timers_level > 3) CALL timer_start(timer_model_init)

    DO jg=1,n_dom
       pat_level(jg)= p_patch(jg)%level
    ENDDO

    IF(iforcing == inwp) THEN
      !
      ! - generate index lists for tiles (land, ocean, lake)
      ! index lists for ice-covered and non-ice covered ocean points 
      ! are initialized in init_nwp_phy
      CALL init_index_lists (p_patch(1:), ext_data)

      CALL configure_atm_phy_nwp(n_dom, pat_level(:), ltestcase, dtime_adv )

     ! initialize number of chemical tracers for convection 
     DO jg = 1, n_dom
       CALL configure_art(jg)
     ENDDO

    ENDIF

    IF (.NOT. ltestcase .AND. iforcing == inwp) THEN
      l_realcase = .TRUE.
    ELSE
      l_realcase = .FALSE.
    ENDIF
 
    ! initialize ldom_active flag
    DO jg=1, n_dom
      IF (jg > 1 .AND. start_time(jg) > 0._wp) THEN
        p_patch(jg)%ldom_active = .FALSE. ! domain not active from the beginning
      ELSE
        p_patch(jg)%ldom_active = .TRUE.
      ENDIF
    ENDDO

    !---------------------------------------------------------------------
    ! 4.c Non-Hydrostatic / NWP
    !---------------------------------------------------------------------

    ALLOCATE (p_nh_state(n_dom), stat=ist)
    IF (ist /= success) THEN
      CALL finish(TRIM(routine),'allocation for p_nh_state failed')
    ENDIF

    ! Note(GZ): Land state now needs to be allocated even if physics is turned
    ! off because ground temperature is included in feedback since r8133
    ! However, setting inwp_surface = 0 effects that only a few 2D fields are allocated
    ALLOCATE (p_lnd_state(n_dom), stat=ist)
    IF (ist /= success) THEN
      CALL finish(TRIM(routine),'allocation for p_lnd_state failed')
    ENDIF

    IF(iforcing /= inwp) atm_phy_nwp_config(:)%inwp_surface = 0



    ! Now allocate memory for the states
    DO jg=1,n_dom
      l_pres_msl(jg) = is_any_output_nml_active(first_output_name_list, &
        &                                 dtime=dtime, iadv_rcf=iadv_rcf, var_name="pres_msl")
    END DO
    CALL construct_nh_state(p_patch(1:), p_nh_state, n_timelevels=2, &
      &                     l_pres_msl=l_pres_msl)

    ! Add optional diagnostic variable lists (might remain empty)
    CALL construct_opt_diag(p_patch(1:), .TRUE.)

    IF(iforcing == inwp) THEN
      DO jg=1,n_dom
        l_rh(jg) = is_any_output_nml_active(first_output_name_list, &
          &                                 dtime=dtime, iadv_rcf=iadv_rcf, var_name="rh")
      END DO
      CALL construct_nwp_phy_state( p_patch(1:), l_rh )
    ENDIF

    CALL construct_nwp_lnd_state( p_patch(1:),p_lnd_state,n_timelevels=2 )


    ! Due to the required ability to overwrite advection-Namelist settings 
    ! via add_ref/add_tracer_ref for ICON-ART, configure_advection is called 
    ! AFTER the nh_state is created. Otherwise, potential modifications of the 
    ! advection-Namelist can not be taken into account properly.
    ! Unfortunatley this conflicts with our trying to call the config-routines 
    ! as early as possible. 
    DO jg =1,n_dom
     CALL configure_advection( jg, p_patch(jg)%nlev, p_patch(1)%nlev,      &
       &                      iequations, iforcing, iqc, iqi, iqr, iqs,    &
       &                      kstart_moist(jg), kend_qvsubstep(jg),        &
       &                      lvert_nest, l_open_ubc, ntracer ) 
    ENDDO


    !---------------------------------------------------------------------
    ! 5. Perform time stepping
    !---------------------------------------------------------------------
      !------------------------------------------------------------------
      ! Prepare for time integration
      !------------------------------------------------------------------


    ! CALL prepare_nh_integration(p_patch(1:), p_nh_state, p_int_state(1:), p_grf_state(1:))

    CALL prepare_nh_integration( )

    !
    ! Read restart files (if necessary)
    !
    IF (is_restart_run()) THEN
      ! This is a resumed integration. Read model state from restart file(s).

#ifdef NOMPI
      CALL read_restart_files
#else
      jg = 1
      !DO jg = n_dom_start,n_dom
      CALL read_restart_files( p_patch(jg) )
      !END DO
#endif      
      CALL message(TRIM(routine),'normal exit from read_restart_files')

    ENDIF


    !---------------------------------------------------------------------
    !     Setup of meteogram output
    !---------------------------------------------------------------------

    IF ((output_mode%l_nml) .OR. (output_mode%l_vlist)) THEN
      DO jg =1,n_dom
        IF (meteogram_output_config(jg)%lenabled) THEN
          ! For dry test cases: do not sample moist variables
          ! (but allow for TORUS moist runs; see also mo_mtgrm_output.F90)
          IF (ltestcase .and. .not. is_plane_torus) THEN
             CALL meteogram_init(meteogram_output_config(jg), jg, p_patch(jg), &
              &                ext_data(jg), p_nh_state(jg),                  &
              &                p_lnd_state=p_lnd_state(jg), iforcing=iforcing)
          ELSE
            CALL meteogram_init(meteogram_output_config(jg), jg, p_patch(jg), &
              &                ext_data(jg), p_nh_state(jg), prm_diag(jg),    &
              &                p_lnd_state(jg), iforcing)
          END IF
        END IF
      END DO
    END IF



    !------------------------------------------------------------------
    ! Prepare initial conditions for time integration.
    !------------------------------------------------------------------
    !
    ! Initialize model with real atmospheric data if appropriate switches are set
    !
    IF (l_realcase .AND. .NOT. is_restart_run()) THEN

      CALL init_icon (p_patch(1:), p_nh_state(1:), prm_diag(1:), p_lnd_state(1:), &
        &             p_int_state(1:), p_grf_state(1:), ext_data(1:))

    ENDIF


    !---------------------------------------------------------
    ! The most primitive event handling algorithm: 
    ! compute time step interval for taking a certain action
    !--------------------------------------------------------- 
 
    ! writing output is now controlled via 'istime4output'
    n_file  = NINT(dt_file/dtime)        ! trigger new output file
    n_chkpt = NINT(dt_checkpoint/dtime)  ! write restart files
    n_diag  = MAX(1,NINT(dt_diag/dtime)) ! diagnose of total integrals


    !------------------------------------------------------------------
    ! Prepare output file
    !------------------------------------------------------------------
    !
    ! if output on z and/or p-levels is required do some config
    DO jg = 1, n_dom
      CALL configure_nh_pzlev(jg, nproma, p_patch(jg)%npromz_c,  &
        &                     p_patch(jg)%nblks_c)
    ENDDO

    ! Map the variable groups given in the output namelist onto the
    ! corresponding variable subsets:
    IF (output_mode%l_nml) THEN
      CALL parse_variable_groups()
    END IF

    IF (.NOT.is_restart_run()) THEN
      ! Initialize the first output file which will contain also the 
      ! initial conditions.

      jfile = 1
      CALL init_output_files(jfile, lclose=.FALSE.)

    ELSE
    ! No need to write out the initial condition, thus no output
    ! during the first integration step. This run will produce
    ! output if n_io <= integration_length. 

      CALL get_restart_attribute('next_output_file',jfile)

!DR may need to be fine-tuned later on
!DR      IF (n_io.le.(nsteps-1)) THEN
         CALL init_output_files(jfile, lclose=.FALSE.)
         l_have_output = .TRUE.  
!DR      ELSE
!DR         l_have_output = .FALSE.
!DR      END IF

    END IF
    
    ! Add a special metrics variable containing the area weights of
    ! the regular lon-lat grid.
    CALL compute_lonlat_area_weights()

    ! setup of post-processing job queue, e.g. setup of optional
    ! diagnostic quantities like pz-level interpolation
    CALL pp_scheduler_init( (.NOT. ltestcase .OR. iforcing == inwp) )

    ! If async IO is in effect, init_name_list_output is a collective call
    ! with the IO procs and effectively starts async IO
    IF (output_mode%l_nml) THEN
      CALL init_name_list_output(isample=iadv_rcf)
    END IF


    ! for debug purpose: print var lists
    IF ( msg_level >=20 .AND. my_process_is_stdio() .AND. .NOT. ltestcase) THEN
      CALL print_var_list (p_nh_state(1)%prog_list(1))
      CALL print_var_list (p_nh_state(1)%diag_list)
      CALL print_var_list (p_nh_state(1)%metrics_list)
      CALL print_var_list (prm_nwp_diag_list(1))
      CALL print_var_list (prm_nwp_tend_list(1))
      CALL print_var_list (p_lnd_state(1)%lnd_prog_nwp_list(1))
      CALL print_var_list (p_lnd_state(1)%lnd_diag_nwp_list)
      CALL print_var_list (ext_data(1)%atm_list)
      CALL print_var_list (ext_data(1)%atm_td_list)
    ENDIF


    IF (timers_level > 3) CALL timer_stop(timer_model_init)

  END SUBROUTINE construct_atmo_nonhydrostatic
  
  !---------------------------------------------------------------------
  SUBROUTINE destruct_atmo_nonhydrostatic
    
    CHARACTER(*), PARAMETER :: routine = "destruct_atmo_nonhydrostatic"


    INTEGER :: jg, ist

    !---------------------------------------------------------------------
    ! 6. Integration finished. Clean up.
    !---------------------------------------------------------------------

    CALL message(TRIM(routine),'start to clean up')

    ! Destruction of post-processing job queue
    CALL pp_scheduler_finalize()

    ! Delete optional diagnostics
    CALL destruct_opt_diag()
   
    ! Delete state variables

    CALL destruct_nh_state( p_nh_state )
    DEALLOCATE (p_nh_state, STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for p_nh_state failed')
    ENDIF

    IF (iforcing == inwp) THEN
      CALL destruct_nwp_phy_state
    ENDIF
    CALL destruct_nwp_lnd_state( p_lnd_state )

    ! Delete output variable lists
    IF (output_mode%l_vlist) THEN
      IF (l_have_output) CALL close_output_files
    END IF
    IF (output_mode%l_nml) THEN
      CALL close_name_list_output
    END IF

    ! finalize meteogram output
    IF ((output_mode%l_nml) .OR. (output_mode%l_vlist)) THEN
      DO jg = 1, n_dom
        IF (meteogram_output_config(jg)%lenabled) THEN
          CALL meteogram_finalize(jg)
        END IF
      END DO
      DO jg = 1, max_dom
        DEALLOCATE(meteogram_output_config(jg)%station_list)
      END DO
    END IF

    CALL message(TRIM(routine),'clean-up finished')
    
  END SUBROUTINE destruct_atmo_nonhydrostatic
  !---------------------------------------------------------------------

END MODULE mo_atmo_nonhydrostatic

