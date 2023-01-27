!>
!! Main program for the ICON ocean waves model
!!
!! @author Mikhail Dobrynin, DWD, 26.04.2019
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
MODULE mo_wave_model

  USE mo_exception,               ONLY: message, finish
  USE mo_mpi,                     ONLY: set_mpi_work_communicators, process_mpi_pref_size, &
       &                                my_process_is_io, my_process_is_pref, my_process_is_mpi_test, &
       &                                stop_mpi, my_process_is_work, process_mpi_io_size, &
       &                                my_process_is_stdio
  USE mo_timer,                   ONLY: init_timer, timer_start, timer_stop, &
       &                                timers_level,timer_model_init, &
       &                                timer_domain_decomp, print_timer
  USE mo_master_config,           ONLY: isRestart
  USE mo_master_control,          ONLY: wave_process
  USE mo_intp_lonlat_types,       ONLY: lonlat_grids
  USE mo_impl_constants,          ONLY: success, pio_type_async, pio_type_cdipio
  USE mo_run_config,              ONLY: configure_run, &
       &                                ltimer, dtime, &
       &                                nshift, num_lev, output_mode, msg_level
  USE mo_time_config,             ONLY: time_config
  USE mo_io_config,               ONLY: restartWritingParameters, configure_io
  USE mo_load_restart,            ONLY: read_restart_header
  USE mo_restart,                 ONLY: detachRestartProcs
  USE mo_name_list_output,        ONLY: name_list_io_main_proc, &
       &                                close_name_list_output
  USE mo_icon_output_tools,       ONLY: init_io_processes
  USE mo_wave_read_namelists,     ONLY: read_wave_namelists
  USE mo_wave_crosscheck,         ONLY: wave_crosscheck
  USE mo_parallel_config,         ONLY: p_test_run, num_test_pe, l_test_openmp, num_io_procs, &
       &                                proc0_shift, num_prefetch_proc, pio_type, num_io_procs_radar, &
       &                                ignore_nproma_use_nblocks_c, nproma, update_nproma_for_io_procs
  USE mo_grid_config,             ONLY: n_dom, n_dom_start
  USE mo_build_decomposition,     ONLY: build_decomposition
  USE mo_zaxis_type,              ONLY: zaxisTypeList, t_zaxisTypeList
  USE mo_wave_state,              ONLY: p_wave_state, p_wave_state_lists, construct_wave_state
  USE mo_model_domain,            ONLY: p_patch
  USE mo_name_list_output_config, ONLY: use_async_name_list_io
  USE mo_name_list_output,        ONLY: write_name_list_output
  USE mo_name_list_output_init,   ONLY: init_name_list_output, parse_variable_groups, &
       &                                output_file, create_vertical_axes
  USE mo_wave,                    ONLY: wave
  USE mo_wave_ext_data_state,     ONLY: wave_ext_data, wave_ext_data_list, destruct_wave_ext_data_state
  USE mo_wave_ext_data_init,      ONLY: init_wave_ext_data
  USE mo_alloc_patches,           ONLY: destruct_patches

  USE mo_intp_data_strc,          ONLY: p_int_state
  USE mo_intp_state,              ONLY: construct_2d_interpol_state, destruct_2d_interpol_state
  USE mo_icon_comm_interface,     ONLY: construct_icon_communication, destruct_icon_communication
  USE mo_interpol_config,         ONLY: configure_interpolation
  USE mo_complete_subdivision,    ONLY: setup_phys_patches


  PUBLIC :: wave_model

CONTAINS
  !-------------------------------------------------------------------
  !>
  SUBROUTINE wave_model(wave_namelist_filename,shr_namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: wave_namelist_filename
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename

    CHARACTER(*), PARAMETER :: routine = "mo_wave_model:wave_model"

    !---------------------------------------------------------------------
    ! construct the wave model
    CALL construct_wave_model(wave_namelist_filename,shr_namelist_filename)

    CALL wave()

    CALL destruct_wave_model()

    CALL message(routine, 'finished')

    ! print performance timers:
    IF (ltimer) CALL print_timer

  END SUBROUTINE wave_model
  !-------------------------------------------------------------------

  !-------------------------------------------------------------------
  !>
  SUBROUTINE construct_wave_model(wave_namelist_filename,shr_namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: wave_namelist_filename
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename

    CHARACTER(*), PARAMETER :: routine = "mo_wave_model:construct_wave_model"
    INTEGER                 :: dedicatedRestartProcs
    LOGICAL                 :: radar_flag_doms_model(n_dom)
    INTEGER :: error_status
    ! initialize global registry of lon-lat grids
    CALL lonlat_grids%init()

    !---------------------------------------------------------------------
    ! 0. If this is a resumed or warm-start run...
    !---------------------------------------------------------------------
    IF (isRestart()) THEN
      CALL message('','Read restart file meta data ...')
      CALL read_restart_header("wave")
    ENDIF

    !---------------------------------------------------------------------
    ! 1.1 Read namelists (newly) specified by the user; fill the
    !     corresponding sections of the configuration states.
    !---------------------------------------------------------------------

    CALL read_wave_namelists(wave_namelist_filename,shr_namelist_filename)

    !---------------------------------------------------------------------
    ! 1.2 Cross-check namelist setups
    !---------------------------------------------------------------------

    CALL wave_crosscheck

    !---------------------------------------------------------------------
    ! 2. Call configure_run to finish filling the run_config state.
    !    This needs to be done very early (but anyway after
    !    wave_crosscheck)
    !    because some components of the state, e.g., num_lev, may be
    !    modified in this subroutine which affects the following
    !    CALLs.
    !---------------------------------------------------------------------

    CALL configure_run()

    !-------------------------------------------------------------------
    ! 3.1 Initialize the mpi work groups
    !-------------------------------------------------------------------
    CALL restartWritingParameters(opt_dedicatedProcCount = dedicatedRestartProcs)

    radar_flag_doms_model(1:n_dom) = .FALSE.
    CALL set_mpi_work_communicators(p_test_run, l_test_openmp,                    &
         &                          num_io_procs, dedicatedRestartProcs,          &
         &                          wave_process, num_prefetch_proc, num_test_pe, &
         &                          pio_type, &
         &                          num_io_procs_radar=num_io_procs_radar,        &
         &                          radar_flag_doms_model=radar_flag_doms_model,  &
         &                          num_dio_procs=proc0_shift)

    !-------------------------------------------------------------------
    ! 3.2 Initialize various timers
    !-------------------------------------------------------------------
    IF (ltimer) CALL init_timer
    IF (timers_level > 1) CALL timer_start(timer_model_init)

    !-------------------------------------------------------------------
    ! initialize dynamic list of vertical axes
    !-------------------------------------------------------------------

    zaxisTypeList = t_zaxisTypeList()


    IF (timers_level > 4) CALL timer_start(timer_domain_decomp)

    CALL build_decomposition(num_lev, nshift, is_ocean_decomposition = .true.)

    IF (timers_level > 4) CALL timer_stop(timer_domain_decomp)

    CALL init_io_processes()

    !--------------------------------------------------------------------------------
    ! 6. Construct interpolation state, compute interpolation coefficients.
    !--------------------------------------------------------------------------------

    CALL configure_interpolation( n_dom, p_patch(1:)%level, &
                                  p_patch(1:)%geometry_info )

    ALLOCATE( p_int_state(n_dom_start:n_dom), &
            & STAT=error_status)
    IF (error_status /= SUCCESS) THEN
      CALL finish(routine, 'allocation for ptr_int_state failed')
    ENDIF

    ! Construct interpolation state
    ! Please note that for parallel runs the divided state is constructed here
    CALL construct_2d_interpol_state(p_patch, p_int_state)

    CALL construct_icon_communication(p_patch, n_dom)

    !------------------------------------------------------------------
    ! Setup the information for the physical patches
    !------------------------------------------------------------------
    CALL setup_phys_patches

    !------------------------------------------------------------------
    ! Create and optionally read external data fields
    !------------------------------------------------------------------
    CALL init_wave_ext_data (p_patch(1:), wave_ext_data, wave_ext_data_list)

    CALL message(routine, 'finished.')

  END SUBROUTINE construct_wave_model
  !-------------------------------------------------------------------

  SUBROUTINE destruct_wave_model()

    CHARACTER(*), PARAMETER :: routine = "mo_wave_model:destruct_wave_model"

    INTEGER                 :: error_status


    ! Deallocate wave external data state and variable lists
    CALL destruct_wave_ext_data_state()

    ! Deallocate interpolation fields
    CALL destruct_2d_interpol_state( p_int_state )
    IF (msg_level>5) CALL message(routine,'destruct_2d_interpol_state is done')

    DEALLOCATE (p_int_state, STAT=error_status)
    IF (error_status /= SUCCESS) THEN
       CALL finish(routine, 'deallocation for ptr_int_state failed')
    ENDIF

    ! Deallocate global registry for lon-lat grids
    CALL lonlat_grids%finalize()

    ! Deallocate grid patches
    CALL destruct_patches( p_patch )
    IF (msg_level>5) CALL message(routine, 'destruct_patches is done')

    DEALLOCATE( p_patch, STAT=error_status )
    IF (error_status/=SUCCESS) THEN
       CALL finish(routine, 'deallocate for patch array failed')
    ENDIF

    ! Delete output variable lists
    IF (output_mode%l_nml) THEN
       CALL message(routine, 'delete output variable lists')
       CALL close_name_list_output
       CALL message(routine, 'finish statistics streams')
    END IF

    CALL destruct_icon_communication()

    CALL message(routine, 'clean-up finished')

  END SUBROUTINE destruct_wave_model
  !-------------------------------------------------------------------

END MODULE mo_wave_model
