!! Initializes and controls the time stepping in the wave model.
!!
!!
!! @author Mikhail Dobrynin, DWD, 24.06.19
!!
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.


!----------------------------
MODULE mo_wave_stepping
  USE mo_exception,                ONLY: message, message_text, finish
  USE mo_run_config,               ONLY: output_mode, ltestcase
  USE mo_name_list_output,         ONLY: write_name_list_output
  USE mo_parallel_config,          ONLY: proc0_offloading
  USE mo_time_config,              ONLY: t_time_config
  USE mtime,                       ONLY: datetime, timedelta,       &
    &                                    OPERATOR(+), OPERATOR(>=)
  USE mo_util_mtime,               ONLY: mtime_utils, FMT_DDHHMMSS_DAYSEP, is_event_active
  USE mo_model_domain,             ONLY: p_patch
  USE mo_grid_config,              ONLY: n_dom, nroot
  USE mo_io_units,                 ONLY: filename_max
  USE mo_master_config,            ONLY: getModelBaseDir
  USE mo_initicon_config,          ONLY: generate_filename

  USE mo_wave_adv_exp,             ONLY: init_wind_adv_test
  USE mo_init_wave_physics,        ONLY: init_wave_phy
  USE mo_wave_state,               ONLY: p_wave_state
  USE mo_wave_ext_data_state,      ONLY: wave_ext_data
  USE mo_wave_forcing_state,       ONLY: wave_forcing_state
  USE mo_wave_diagnostics,         ONLY: significant_wave_height
  USE mo_wave_physics,             ONLY: new_spectrum, total_energy, &
                                         mean_frequency_energy, air_sea
  USE mo_wave_config,              ONLY: wave_config
  USE mo_wave_forcing_state,       ONLY: wave_forcing_state
  USE mo_wave_events,              ONLY: create_wave_events, dummyWaveEvent

  IMPLICIT NONE

  PRIVATE

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_wave_stepping'

  PUBLIC :: perform_wave_stepping

CONTAINS

  !-------------------------------------------------------------------------
  !>
  !! Organizes wave time stepping
  !!
  !! @par Revision History
  !! Initial revision by Mikhail Dobrynin, DWD, (2019-06-24)
  !!
  SUBROUTINE perform_wave_stepping (time_config)
    TYPE(t_time_config), INTENT(IN) :: time_config  !< information for time control

    TYPE(datetime),  POINTER :: mtime_current     => NULL() !< current datetime
    TYPE(timedelta), POINTER :: model_time_step   => NULL()
    INTEGER                  :: jstep                       !< time step number
    LOGICAL                  :: lprint_timestep             !< print current datetime information
    INTEGER                  :: jg, jt, jlev

    CHARACTER(len=*), PARAMETER :: routine = modname//':perform_wave_stepping'


    CHARACTER(LEN=filename_max) :: wave_forc_wind_fn(n_dom) ! forc_file_prefix+'_wind' for U and V 10 meter wind (m/s)
    CHARACTER(LEN=filename_max) :: wave_forc_ice_fn(n_dom) ! forc_file_prefix+'_ice'  for sea ice concentration (fraction of 1)
    CHARACTER(LEN=filename_max) :: wave_forc_slh_fn(n_dom) ! forc_file_prefix+'_slh'  for sea level height (m
    CHARACTER(LEN=filename_max) :: wave_forc_osc_fn(n_dom) ! forc_file_prefix+'_osc'  for U and V ocean surface currents (m/s)

    DO jg = 1, n_dom
      IF (TRIM(wave_config(jg)%forc_file_prefix) /= '') THEN

        ! TO DO generate name of input forcing files 
        jlev = p_patch(jg)%level

        wave_forc_wind_fn(jg) = " "
        wave_forc_ice_fn(jg) = " "
        wave_forc_slh_fn(jg) = " "
        wave_forc_osc_fn(jg) = " "

        wave_forc_wind_fn(jg) = generate_filename(TRIM(wave_config(jg)%forc_file_prefix)//"_wind", getModelBaseDir(), &
             &                  nroot, jlev, jg)
        wave_forc_ice_fn(jg) = generate_filename(TRIM(wave_config(jg)%forc_file_prefix)//"_ice", getModelBaseDir(), &
             &                  nroot, jlev, jg)
        wave_forc_slh_fn(jg) = generate_filename(TRIM(wave_config(jg)%forc_file_prefix)//"_slh", getModelBaseDir(), &
             &                  nroot, jlev, jg)
        wave_forc_osc_fn(jg) = generate_filename(TRIM(wave_config(jg)%forc_file_prefix)//"_osc", getModelBaseDir(), &
             &                  nroot, jlev, jg)

        CALL message(routine,'10m wind from: '//wave_forc_wind_fn(jg))
        CALL message(routine,'ice concentration from: '//wave_forc_ice_fn(jg))
        CALL message(routine,'sea level height from: '//wave_forc_slh_fn(jg))
        CALL message(routine,'ocean surface currents from: '//wave_forc_osc_fn(jg))

        CALL finish(routine,'read of forcing data from '//TRIM(wave_config(jg)%forc_file_prefix)//'* files is not implemented yet')


      END IF
    END DO

    IF (ltestcase) THEN
      !-----------------------------------------------------------------------
      ! advection experiment
      CALL message(routine,'test case run: advection experiment')

      DO jg = 1, n_dom
        !initialisation of 10 meter wind
        CALL init_wind_adv_test(p_patch(jg), wave_forcing_state(jg))

        DO jt = 1, 2 !time levels
          CALL init_wave_phy(p_patch(jg), &
               wave_config(jg), &
               p_wave_state(jg)%prog(jt), &
               p_wave_state(jg)%diag, &
               wave_ext_data(jg), &
               wave_forcing_state(jg))

          CALL new_spectrum(p_patch(jg), &
               wave_config(jg), &
               p_wave_state(jg)%diag, &
               wave_forcing_state(jg), &
               p_wave_state(jg)%prog(jt)%tracer)

          CALL total_energy(p_patch(jg), &
               wave_config(jg), &
               p_wave_state(jg)%prog(jt)%tracer, &
               p_wave_state(jg)%diag%llws,&
               p_wave_state(jg)%diag%emean, &
               p_wave_state(jg)%diag%emeanws)

          CALL mean_frequency_energy(p_patch(jg), &
               wave_config(jg), &
               p_wave_state(jg)%prog(jt)%tracer, &
               p_wave_state(jg)%diag%llws,&
               p_wave_state(jg)%diag%emean, &
               p_wave_state(jg)%diag%emeanws, &
               p_wave_state(jg)%diag%femean, &
               p_wave_state(jg)%diag%femeanws)
        END DO

        CALL air_sea(p_patch(jg), &
             wave_config(jg), &
             wave_forcing_state(jg)%sp10m, &
             p_wave_state(jg)%diag%tauw, &
             p_wave_state(jg)%diag%ustar, &
             p_wave_state(jg)%diag%z0)

        CALL significant_wave_height(p_patch = p_patch(jg), &
             &                       emean   = p_wave_state(jg)%diag%emean(:,:), &
             &                       hs      = p_wave_state(jg)%diag%hs(:,:))

      END DO


      ! advection experiment
      !-----------------------------------------------------------------------
    ELSE
      CALL message(routine,'normal run: work in progress...')
    END IF

    ! create wave events
    CALL create_wave_events(time_config)

    ! convenience pointer
    mtime_current   => time_config%tc_current_date
    model_time_step => time_config%tc_dt_model

    ! TODO: write logical function such that the timestep information is
    !       prnted only under certain conditions (see atmospheric code)
    lprint_timestep = .TRUE.

    ! initialize time step counter
    jstep = 0

    ! output at initial time
    IF (output_mode%l_nml) THEN
      CALL write_name_list_output(jstep=jstep)
    END IF


    TIME_LOOP: DO

      ! update model date and time
      mtime_current = mtime_current + model_time_step
      jstep = jstep + 1

      IF (lprint_timestep) THEN
        CALL message('','')

        WRITE(message_text,'(a,i8,a,i0,a,5(i2.2,a),i3.3,a,a)') &
          &             'Time step: ', jstep, ', model time: ',                              &
          &             mtime_current%date%year,   '-', mtime_current%date%month,    '-',    &
          &             mtime_current%date%day,    ' ', mtime_current%time%hour,     ':',    &
          &             mtime_current%time%minute, ':', mtime_current%time%second,   '.',    &
          &             mtime_current%time%ms, ' forecast time ',                            &
          &             TRIM(mtime_utils%ddhhmmss(time_config%tc_exp_startdate, &
          &                                       mtime_current, FMT_DDHHMMSS_DAYSEP))

        CALL message('',message_text)
      ENDIF

      IF (is_event_active(dummyWaveEvent, mtime_current, proc0_offloading)) THEN
        WRITE(message_text,'(a)') "dummyWaveEvent is active"
        CALL message('',message_text)

        !----------
        !TEST change sea level height
        wave_forcing_state(1)%sea_level_c = wave_forcing_state(1)%sea_level_c + 5.0
        wave_ext_data(1)%bathymetry_c = wave_ext_data(1)%bathymetry_c + wave_forcing_state(1)%sea_level_c
        !TEST
        !-----------
      ENDIF

      IF (output_mode%l_nml) THEN
        CALL write_name_list_output(jstep=jstep)
      END IF

      IF (mtime_current >= time_config%tc_stopdate) THEN
        ! leave time loop
        EXIT TIME_LOOP
      END IF

    ENDDO TIME_LOOP


    CALL message(routine,'finished')
  END SUBROUTINE perform_wave_stepping

END MODULE mo_wave_stepping
