! @brief Interface between atmosphere and the O3 provider, through a coupler
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

MODULE mo_atmo_o3_provider_coupling

  USE mo_exception,       ONLY: finish
#ifdef YAC_coupling
  USE mo_yac_finterface, ONLY: yac_fget_field_collection_size
#endif
  USE mo_coupling_utils, ONLY: def_field

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: str_module = 'mo_atmo_o3_provider_coupling' ! Output of module for debug

  PUBLIC :: construct_atmo_o3_provider_coupling_post_sync

  PUBLIC :: field_id_o3, nplev_o3_provider

  INTEGER :: field_id_o3
  INTEGER :: nplev_o3_provider

CONTAINS

  SUBROUTINE construct_atmo_o3_provider_coupling_post_sync( &
    comp_id, cell_point_id, timestepstring)

    INTEGER, INTENT(IN) :: comp_id
    INTEGER, INTENT(IN) :: cell_point_id
    CHARACTER(LEN=*), INTENT(IN) :: timestepstring

    CHARACTER(LEN=*), PARAMETER   :: &
      routine = str_module // ':construct_atmo_o3_provider_coupling_post_sync'

#ifndef YAC_coupling
    CALL finish(routine, 'built without coupling support.')
#else

    nplev_o3_provider = &
      yac_fget_field_collection_size( &
        "o3_provider", "o3_grid", "o3")

    CALL def_field( &
      comp_id, cell_point_id, timestepstring, &
      "o3", nplev_o3_provider, field_id_o3)

! YAC_coupling
#endif

  END SUBROUTINE construct_atmo_o3_provider_coupling_post_sync
END MODULE mo_atmo_o3_provider_coupling
