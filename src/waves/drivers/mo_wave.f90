!>
!! @brief branch for the non-hydrostatic ICON workflow
!!
!! @author Mikhail Dobrynin, DWD, 20.06.19
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_wave
  USE mo_kind,                  ONLY: wp
  USE mo_exception,             ONLY: message
  USE mo_model_domain,          ONLY: p_patch
  USE mo_grid_config,           ONLY: n_dom, start_time, end_time
  USE mo_wave_state,            ONLY: construct_wave_state, destruct_wave_state
  USE mo_time_config,           ONLY: time_config
  USE mo_util_mtime,            ONLY: getElapsedSimTimeInSeconds
  USE mo_output_event_types,    ONLY: t_sim_step_info
  USE mo_name_list_output_init, ONLY: init_name_list_output, &
       &                              output_file, create_vertical_axes
  USE mo_run_config,            ONLY: dtime, output_mode
  USE mo_io_config,             ONLY: configure_io
  USE mo_wave_stepping,         ONLY: perform_wave_stepping

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: wave


CONTAINS

  SUBROUTINE wave

    CHARACTER(*), PARAMETER :: routine = "mo_wave:wave"

    CALL construct_wave()

    CALL perform_wave_stepping(time_config%tc_current_date)

    CALL destruct_wave()

    CALL message(TRIM(routine),'finished')
  END SUBROUTINE wave


  SUBROUTINE construct_wave()

    CHARACTER(*), PARAMETER :: routine = "construct_wave"

    TYPE(t_sim_step_info) :: sim_step_info

    INTEGER :: jg
    REAL(wp):: sim_time


    ! calculate elapsed simulation time in seconds
    sim_time = getElapsedSimTimeInSeconds(time_config%tc_current_date)

    DO jg=1, n_dom
      p_patch(jg)%ldom_active &
           =        (jg <= 1 .OR. start_time(jg) <= sim_time) &
           &  .AND. end_time(jg) > sim_time
    END DO

    CALL construct_wave_state(p_patch(1:))

    CALL configure_io()   ! set n_chkpt and n_diag, which control
                          ! writing of restart files and tot_int diagnostics.

    ! If async IO is in effect, init_name_list_output is a collective call
    ! with the IO procs and effectively starts async IO
    IF (output_mode%l_nml) THEN
      ! compute sim_start, sim_end
      sim_step_info%sim_start = time_config%tc_exp_startdate
      sim_step_info%sim_end = time_config%tc_exp_stopdate
      sim_step_info%run_start = time_config%tc_startdate
      sim_step_info%restart_time = time_config%tc_stopdate

      sim_step_info%dtime      = dtime
      sim_step_info%jstep0 = 0

      CALL init_name_list_output(sim_step_info)

      CALL create_vertical_axes(output_file)
    END IF

    CALL message(TRIM(routine),'finished')

  END SUBROUTINE construct_wave


  SUBROUTINE destruct_wave()

    CHARACTER(*), PARAMETER :: routine = "destruct_wave"


    CALL destruct_wave_state( )

    CALL message(TRIM(routine),'finished')

  END SUBROUTINE destruct_wave


END MODULE mo_wave
