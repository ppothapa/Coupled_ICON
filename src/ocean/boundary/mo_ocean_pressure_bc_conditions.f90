!! ===========================================================================================================================
!! Implementation of tides by computation of the Sun's and Moon's full tidal potential
!! This will be used in the pressure gradient calculation
!!
!! Authors: Kai Logemann, Helmholtz-Zentrum Geesthacht, and Leonidas Linardakis, Max-Planck-Institute for Meteorology, Hamburg
!!
!! ===========================================================================================================================

!----------------------------
#include "icon_definitions.inc"
#include "omp_definitions.inc"
#include "iconfor_dsl_definitions.inc"
!----------------------------
MODULE mo_ocean_pressure_bc_conditions
  !-------------------------------------------------------------------------
  USE mo_kind,                   ONLY: wp, dp
   USE mtime,                    ONLY: datetime
  USE mo_exception,              ONLY: finish
  USE mo_model_domain,           ONLY: t_patch_3d
  USE mo_ocean_nml,              ONLY: use_tides,   &
    & use_tides_SAL, atm_pressure_included_in_ocedyn,       &
    & OceanReferenceDensity_inv, vert_cor_type
  USE mo_physical_constants,     ONLY: grav
 ! USE mo_grid_subset,            ONLY: t_subset_range, get_index_range
 ! USE mo_parallel_config,        ONLY: nproma
  USE mo_impl_constants,         ONLY: sea_boundary
  USE mo_util_dbg_prnt,          ONLY: dbg_print, debug_printValue
  USE mo_ocean_tides,            ONLY: calculate_tides_potential
  USE mo_ocean_surface_types,    ONLY: t_atmos_for_ocean
  USE mo_ocean_types,            ONLY: t_hydro_ocean_state
  USE mo_sea_ice_types,          ONLY: t_sea_ice
  USE mo_dynamics_config,        ONLY: nold


  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: create_pressure_bc_conditions

  CHARACTER(LEN=12)  :: str_module = 'mo_ocean_pressure_bc_conditions'  ! Output of module for 1 line debug
  !-------------------------------------------------------------------------

CONTAINS

  !-------------------------------------------------------------------------
  SUBROUTINE create_pressure_bc_conditions(patch_3d,ocean_state, p_as,sea_ice, current_time, use_acc)
    TYPE(t_patch_3d ),TARGET, INTENT(in)             :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state
    TYPE(t_atmos_for_ocean), TARGET, INTENT(in)      :: p_as
    TYPE (t_sea_ice), TARGET, INTENT(in)             :: sea_ice
    TYPE(datetime), POINTER, INTENT(in)              :: current_time
    LOGICAL, INTENT(IN), OPTIONAL                    :: use_acc

    REAL(wp) :: switch_atm_pressure, switch_vert_cor_type, switch_tides
    LOGICAL  :: lacc
    CHARACTER(len=*), PARAMETER :: routine = 'create_pressure_bc_conditions'

    IF (PRESENT(use_acc)) THEN
      lacc = use_acc
    ELSE
      lacc = .FALSE.
    END IF

    !------------------------------------------------------------------------
    IF (use_tides .OR. use_tides_SAL) THEN
#ifdef _OPENACC
      CALL finish(routine, 'use_tides not ported')
#endif
      ! compute tidal potential
      CALL calculate_tides_potential(patch_3d,current_time,ocean_state%p_diag%rho, ocean_state%p_prog(nold(1))%h, &
           ocean_state%p_aux%bc_tides_potential, ocean_state%p_aux%bc_SAL_potential)
    ENDIF

    !------------------------------------------------------------------------
    ! total top potential (individiual terms)

    IF (atm_pressure_included_in_ocedyn) THEN
      switch_atm_pressure=1.0_wp
    ELSE
      switch_atm_pressure=0.0_wp
    ENDIF

    IF (vert_cor_type .EQ. 1 ) THEN
      switch_vert_cor_type=1.0_wp
    ELSE
      switch_vert_cor_type=0.0_wp
    ENDIF

    IF (use_tides .OR. use_tides_SAL) THEN
      switch_tides=1.0_wp
    ELSE
      switch_tides=0.0_wp
    ENDIF

    !------------------------------------------------------------------------
    ! total top potential

    IF ( (switch_atm_pressure + switch_vert_cor_type + switch_tides) > 0.0_wp ) THEN

      !$ACC KERNELS DEFAULT(PRESENT) IF(lacc)
      ocean_state%p_aux%bc_total_top_potential =  &
        & ocean_state%p_aux%bc_tides_potential &
        & + ocean_state%p_aux%bc_SAL_potential &
        & + p_as%pao * OceanReferenceDensity_inv * switch_atm_pressure & ! add acceleration by air pressure
        & + grav * sea_ice%draftave * switch_vert_cor_type ! only zstar: pressure of sea ice on top of the first layer (divided by rhoref to create an acceleration)
      !$ACC END KERNELS

    ENDIF

     IF (use_tides .OR. use_tides_SAL) THEN
      CALL dbg_print('tides_potential',  ocean_state%p_aux%bc_tides_potential, &
           str_module, 3, in_subset=patch_3d%p_patch_2d(1)%cells%owned)
      CALL dbg_print('tides_SAL',      ocean_state%p_aux%bc_SAL_potential, &
           str_module, 3, in_subset=patch_3d%p_patch_2d(1)%cells%owned)
    ENDIF

  END SUBROUTINE create_pressure_bc_conditions
  !-------------------------------------------------------------------------


END MODULE mo_ocean_pressure_bc_conditions
!=============================================================================
