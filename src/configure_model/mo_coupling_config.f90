!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
MODULE mo_coupling_config

  IMPLICIT NONE
  PRIVATE

  !>
  !! Namelist input to steer the coupling modes
  !! note that default is potentially overwritten in corresponding Namelist routine(s)
  !!
  LOGICAL :: config_coupled_to_ocean = .FALSE.
  LOGICAL :: config_coupled_to_waves = .FALSE.
  LOGICAL :: config_coupled_to_atmo  = .FALSE.


  ! variables
  PUBLIC :: config_coupled_to_ocean
  PUBLIC :: config_coupled_to_waves
  PUBLIC :: config_coupled_to_atmo

  ! functions
  PUBLIC :: is_coupled_run
  PUBLIC :: is_coupled_to_ocean
  PUBLIC :: is_coupled_to_waves
  PUBLIC :: is_coupled_to_atmo

CONTAINS

  !------------------------------------------------------------------------
  LOGICAL FUNCTION is_coupled_run()

    is_coupled_run = config_coupled_to_ocean .OR.  &
      &              config_coupled_to_waves .OR.  &
      &              config_coupled_to_atmo

  END FUNCTION is_coupled_run

  !------------------------------------------------------------------------
  LOGICAL FUNCTION is_coupled_to_ocean()

    is_coupled_to_ocean = config_coupled_to_ocean

  END FUNCTION is_coupled_to_ocean

  !------------------------------------------------------------------------
  LOGICAL FUNCTION is_coupled_to_waves()

    is_coupled_to_waves = config_coupled_to_waves

  END FUNCTION is_coupled_to_waves

  !------------------------------------------------------------------------
  LOGICAL FUNCTION is_coupled_to_atmo()

    is_coupled_to_atmo = config_coupled_to_atmo

  END FUNCTION is_coupled_to_atmo

  !------------------------------------------------------------------------

END MODULE mo_coupling_config
