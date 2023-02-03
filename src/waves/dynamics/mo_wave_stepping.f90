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
  USE mo_exception,                ONLY: message
  USE mo_run_config,               ONLY: output_mode, ltestcase
  USE mo_name_list_output,         ONLY: write_name_list_output
  USE mtime,                       ONLY: datetime
  USE mo_model_domain,             ONLY: p_patch
  USE mo_grid_config,              ONLY: n_dom

  USE mo_wave_adv_exp,             ONLY: init_wind_adv_test
  USE mo_init_wave_physics,        ONLY: init_wave_phy
  USE mo_wave_state,               ONLY: p_wave_state
  USE mo_wave_diagnostics,         ONLY: significant_wave_height
  USE mo_wave_physics,             ONLY: new_spectrum, total_energy, &
                                         mean_frequency_energy, air_sea
  USE mo_wave_config,              ONLY: wave_config

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
  SUBROUTINE perform_wave_stepping (mtime_current)
    TYPE(datetime), INTENT(IN) :: mtime_current

    CHARACTER(len=*), PARAMETER :: routine = modname//':perform_wave_stepping'

    INTEGER :: jg, jt

    IF (ltestcase) THEN
      !-----------------------------------------------------------------------
      ! advection experiment
      CALL message(TRIM(routine),'test case run: advection experiment')

      DO jg = 1, n_dom
        !initialisation of 10 meter wind
        CALL init_wind_adv_test(p_patch(jg), p_wave_state(jg)%diag)
        DO jt = 1, 2 ! @wave initialisation for the time level 1 only?
          CALL init_wave_phy(p_patch(jg), &
               wave_config(jg), &
               p_wave_state(jg)%prog(jt), &
               p_wave_state(jg)%diag)

          CALL new_spectrum(p_patch(jg), &
               wave_config(jg), &
               p_wave_state(jg)%prog(jt), &
               p_wave_state(jg)%diag)

          CALL total_energy(p_patch(jg), &
               wave_config(jg), &
               p_wave_state(jg)%prog(jt), &
               p_wave_state(jg)%diag)

          CALL mean_frequency_energy(p_patch(jg), &
               wave_config(jg), &
               p_wave_state(jg)%prog(jt), &
               p_wave_state(jg)%diag)

          CALL air_sea(p_patch(jg), wave_config(jg), p_wave_state(jg)%diag )

          CALL significant_wave_height(p_patch = p_patch(jg), &
               &                          emean   = p_wave_state(jg)%diag%emean(:,:), &
               &                          hs      = p_wave_state(jg)%diag%hs(:,:))

        END DO
      END DO
      ! advection experiment
      !-----------------------------------------------------------------------
    ELSE
      CALL message(TRIM(routine),'normal run: work in progress...')
    END IF

    IF (output_mode%l_nml) THEN
      CALL write_name_list_output(jstep=0)
    END IF

    CALL message(TRIM(routine),'finished')
  END SUBROUTINE perform_wave_stepping

END MODULE mo_wave_stepping
