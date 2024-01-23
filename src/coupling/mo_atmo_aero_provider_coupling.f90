! @brief Interface between atmosphere and the aero provider, through a coupler
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

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_atmo_aero_provider_coupling

  USE mo_exception,      ONLY: finish
#ifdef YAC_coupling
  USE mo_yac_finterface, ONLY: yac_fget_field_collection_size
#endif
  USE mo_coupling_utils, ONLY: def_field

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: str_module = 'mo_atmo_aero_provider_coupling' ! Output of module for debug

  PUBLIC :: construct_atmo_aero_provider_coupling_post_sync

  PUBLIC :: field_id_aod_c_f, field_id_ssa_c_f, field_id_z_km_aer_c_mo, &
            field_id_aod_c_s, field_id_ssa_c_s, field_id_asy_c_s, &
            field_id_aod_f_s, field_id_ssa_f_s, field_id_asy_f_s, &
            field_id_z_km_aer_f_mo
  PUBLIC:: nblw_aero_provider, nbsw_aero_provider, nlev_aero_provider

  INTEGER :: field_id_aod_c_f
  INTEGER :: field_id_ssa_c_f
  INTEGER :: field_id_z_km_aer_c_mo
  INTEGER :: field_id_aod_c_s
  INTEGER :: field_id_ssa_c_s
  INTEGER :: field_id_asy_c_s
  INTEGER :: field_id_aod_f_s
  INTEGER :: field_id_ssa_f_s
  INTEGER :: field_id_asy_f_s
  INTEGER :: field_id_z_km_aer_f_mo

  INTEGER :: nblw_aero_provider
  INTEGER :: nbsw_aero_provider
  INTEGER :: nlev_aero_provider

CONTAINS

  SUBROUTINE construct_atmo_aero_provider_coupling_post_sync( &
    comp_id, cell_point_id, timestepstring)

    INTEGER, INTENT(IN) :: comp_id
    INTEGER, INTENT(IN) :: cell_point_id
    CHARACTER(LEN=*), INTENT(IN) :: timestepstring

    CHARACTER(LEN=*), PARAMETER   :: &
      routine = str_module // ':construct_atmo_aero_provider_coupling_post_sync'

#ifdef YAC_coupling
    nblw_aero_provider = &
      yac_fget_field_collection_size( &
        "aero_provider", "aero_grid", "aod_lw_b16_coa")
    nbsw_aero_provider = &
      yac_fget_field_collection_size( &
        "aero_provider", "aero_grid", "aod_sw_b14_coa")
    nlev_aero_provider = &
      yac_fget_field_collection_size( &
        "aero_provider", "aero_grid", "aer_lw_b16_coa")
#else
    CALL finish(routine, 'built without coupling support.')
#endif

    CALL def_field( &
      comp_id, cell_point_id, timestepstring, &
      "aod_lw_b16_coa", nblw_aero_provider, field_id_aod_c_f)
    CALL def_field( &
      comp_id, cell_point_id, timestepstring, &
      "ssa_lw_b16_coa", nblw_aero_provider, field_id_ssa_c_f)
    CALL def_field( &
      comp_id, cell_point_id, timestepstring, &
      "aer_lw_b16_coa", nlev_aero_provider, field_id_z_km_aer_c_mo)
    CALL def_field( &
      comp_id, cell_point_id, timestepstring, &
      "aod_sw_b14_coa", nbsw_aero_provider, field_id_aod_c_s)
    CALL def_field( &
      comp_id, cell_point_id, timestepstring, &
      "ssa_sw_b14_coa", nbsw_aero_provider, field_id_ssa_c_s)
    CALL def_field( &
      comp_id, cell_point_id, timestepstring, &
      "asy_sw_b14_coa", nbsw_aero_provider, field_id_asy_c_s)
    CALL def_field( &
      comp_id, cell_point_id, timestepstring, &
      "aod_sw_b14_fin", nbsw_aero_provider, field_id_aod_f_s)
    CALL def_field( &
      comp_id, cell_point_id, timestepstring, &
      "ssa_sw_b14_fin", nbsw_aero_provider, field_id_ssa_f_s)
    CALL def_field( &
      comp_id, cell_point_id, timestepstring, &
      "asy_sw_b14_fin", nbsw_aero_provider, field_id_asy_f_s)
    CALL def_field( &
      comp_id, cell_point_id, timestepstring, &
      "aer_sw_b14_fin", nlev_aero_provider, field_id_z_km_aer_f_mo)

  END SUBROUTINE construct_atmo_aero_provider_coupling_post_sync
END MODULE mo_atmo_aero_provider_coupling
