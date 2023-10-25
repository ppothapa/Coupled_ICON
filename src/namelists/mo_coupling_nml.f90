! Contains the variables to set up the coupling.
!
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

MODULE mo_coupling_nml

  !-------------------------------------------------------------------------
  !
  !    ProTeX FORTRAN source: Style 2
  !    modified for ICON project, DWD/MPI-M 2006
  !
  !-------------------------------------------------------------------------

  USE mo_impl_constants,  ONLY: max_char_length
  USE mo_io_units,        ONLY: nnml
  USE mo_namelist,        ONLY: open_nml, close_nml, position_nml, POSITIONED
  USE mo_exception,       ONLY: finish
  USE mo_coupling_config, ONLY: config_coupled_to_ocean, config_coupled_to_waves,       &
    &                           config_coupled_to_atmo, config_use_sens_heat_flux_hack, &
    &                           config_suppress_sens_heat_flux_hack_over_ice
  USE mo_coupling,        ONLY: coupler_config_files_exist

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_coupling_namelist

CONTAINS

  !!  Initialization of variables that contain general information.
  !!
  !!               Initialization of variables that contain general information
  !!               about the coupled model run. The configuration is read from
  !!               namelist 'icon_cpl'.
  !!

  SUBROUTINE read_coupling_namelist (namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: namelist_filename

    !
    ! Local variables
    !
    LOGICAL :: coupled_to_ocean, coupled_to_waves, coupled_to_atmo
    LOGICAL :: coupled_mode
    LOGICAL :: use_sens_heat_flux_hack
    LOGICAL :: suppress_sens_heat_flux_hack_over_ice
    INTEGER :: istat

    CHARACTER(len=max_char_length), PARAMETER :: &
         &   routine = 'mo_coupling_nml:read_coupling_namelist'

    NAMELIST /coupling_mode_nml/ coupled_to_ocean, coupled_to_waves, coupled_to_atmo, &
         use_sens_heat_flux_hack, suppress_sens_heat_flux_hack_over_ice

    !--------------------------------------------------------------------
    ! 1. Set default values
    !--------------------------------------------------------------------

    coupled_to_ocean = .FALSE.
    coupled_to_waves = .FALSE.
    coupled_to_atmo  = .FALSE.
    use_sens_heat_flux_hack = .FALSE.
    suppress_sens_heat_flux_hack_over_ice = .FALSE.

    !--------------------------------------------------------------------
    ! 2. Read user's (new) specifications (done so far by all MPI processes)
    !--------------------------------------------------------------------

#ifdef YAC_coupling

    CALL open_nml (TRIM(namelist_filename))

    CALL position_nml('coupling_mode_nml',STATUS=istat)
    IF (istat==POSITIONED) THEN
      READ (nnml, coupling_mode_nml)
    ENDIF

    CALL close_nml

#endif

    config_coupled_to_ocean = coupled_to_ocean
    config_coupled_to_waves = coupled_to_waves
    config_coupled_to_atmo  = coupled_to_atmo
    config_use_sens_heat_flux_hack = use_sens_heat_flux_hack
    config_suppress_sens_heat_flux_hack_over_ice = suppress_sens_heat_flux_hack_over_ice


    !----------------------------------------------------
    ! 3. Sanity checks
    !----------------------------------------------------

    coupled_mode = ANY((/coupled_to_ocean,coupled_to_waves,coupled_to_atmo/))

    IF (coupled_mode .AND. .NOT. coupler_config_files_exist()) THEN
      CALL finish( &
        routine, &
        'run is configured to be coupled, but coupler configuration files are not available')
    END IF

  END SUBROUTINE read_coupling_namelist

END MODULE mo_coupling_nml
