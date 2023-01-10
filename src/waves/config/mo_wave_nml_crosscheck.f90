!>
!! This module checks the read-in namelist parameters and, in case of
!! inconsistencies, it tries to correct these.
!!
!!
!! @author Mikhail Dobrynin, DWD, 04.06.2019
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
MODULE mo_wave_crosscheck

  !USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message!, finish, warning
!  USE mo_grid_config,       ONLY: init_grid_configuration
  USE mo_parallel_config,   ONLY: check_parallel_configuration!, p_test_run, l_fast_sum, &
!       &                          use_dp_mpi2io
  USE mo_run_config,        ONLY: nsteps!, dtime, nlev
  USE mo_time_config,       ONLY: time_config, dt_restart
!  USE mo_io_config,         ONLY: dt_checkpoint, write_initial_state, lnetcdf_flt64_output
!  USE mo_grid_config,       ONLY: grid_rescale_factor, use_duplicated_connectivity
!  USE mo_wave_nml
!  USE mo_master_config,     ONLY: isRestart
  USE mo_time_management,   ONLY: compute_timestep_settings,                        &
       &                          compute_restart_settings,                         &
       &                          compute_date_settings
!  USE mo_event_manager,     ONLY: initEventManager

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: wave_crosscheck

  CHARACTER(LEN = *), PARAMETER :: modname = "mo_wave_nml_crosscheck"

CONTAINS

  SUBROUTINE wave_crosscheck

    CHARACTER(len=*), PARAMETER :: routine =   modname//'::wave_crosscheck'


    !--------------------------------------------------------------------
    ! Compute date/time/time step settings
    ! and initialize the event manager
    !--------------------------------------------------------------------
    !
    ! Note that the ordering of the following three calls must not be
    ! changed, since they rely on previous results:
    !
    CALL compute_timestep_settings()
    CALL compute_restart_settings()
    CALL compute_date_settings("wave", dt_restart, nsteps)
    !
    ! Create an event manager, ie. a collection of different events
    !
!    CALL initEventManager(time_config%tc_exp_refdate)

    !--------------------------------------------------------------------
    ! Parallelization
    !--------------------------------------------------------------------
    CALL check_parallel_configuration()

    !--------------------------------------------------------------------
    ! @WAVES TO DO: check the wave model settings
    !--------------------------------------------------------------------
    !
    !ADD HERE ...

    CALL message(routine,'done.')

  END SUBROUTINE wave_crosscheck

END MODULE mo_wave_crosscheck
