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

  USE mo_exception,         ONLY: message, finish
  USE mo_parallel_config,   ONLY: check_parallel_configuration
  USE mo_coupling_config,   ONLY: is_coupled_run
  USE mo_run_config,        ONLY: nsteps, ldynamics, ltransport, ntracer, num_lev, ltestcase
  USE mo_grid_config,       ONLY: n_dom
  USE mo_time_config,       ONLY: time_config, dt_restart
  USE mo_time_management,   ONLY: compute_timestep_settings,                        &
       &                          compute_restart_settings,                         &
       &                          compute_date_settings
  USE mo_event_manager,     ONLY: initEventManager
  USE mo_wave_config,       ONLY: wave_config

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: wave_crosscheck

  CHARACTER(LEN = *), PARAMETER :: modname = "mo_wave_nml_crosscheck"

CONTAINS

  SUBROUTINE wave_crosscheck

    CHARACTER(len=*), PARAMETER :: routine =   modname//'::wave_crosscheck'
    INTEGER :: jg

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
    CALL initEventManager(time_config%tc_exp_refdate)

    !--------------------------------------------------------------------
    ! Parallelization
    !--------------------------------------------------------------------
    CALL check_parallel_configuration()

    IF (.not.ldynamics) THEN
      CALL finish(TRIM(routine),'Error: ldynamics must be TRUE')
    END IF

    IF (.not.ltransport) THEN
      CALL finish(TRIM(routine),'Error: ltransport must be TRUE')
    END IF

    DO jg=1,n_dom
      IF (ntracer /= wave_config(jg)%nfreqs*wave_config(jg)%ndirs) THEN
        CALL finish(TRIM(routine),'Error: ntracer must be equal to nfreqs*ndirs')
      END IF
    ENDDO

    ! For standalone runs reading of forcing data from file is mandatory.
    ! The only exception are testcase runs (ltestcase=.TRUE).
    IF (.NOT. is_coupled_run()) THEN
      DO jg=1,n_dom
        IF (TRIM(wave_config(jg)%forc_file_prefix) == '' .AND. .NOT.ltestcase) THEN
          CALL finish(TRIM(routine),'Error: For standalone runs reading of forcing data from file is mandatory')
        ENDIF
      ENDDO
    ENDIF

    IF (ANY(num_lev(1:n_dom).ne.1)) THEN
      CALL finish(TRIM(routine),'Error: num_lev must be 1')
    END IF

    CALL message(routine,'finished.')

  END SUBROUTINE wave_crosscheck

END MODULE mo_wave_crosscheck
