! Set of routines shared by various coupling related modules
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

MODULE mo_coupling_utils

#if !defined NOMPI && defined YAC_coupling
  USE mo_yac_finterface,   ONLY: yac_fdef_field, yac_fdef_field_mask, &
    &                            YAC_TIME_UNIT_ISO_FORMAT
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: def_field

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_coupling_utils'

  INTERFACE def_field
    MODULE PROCEDURE def_field_no_mask
    MODULE PROCEDURE def_field_mask
  END INTERFACE def_field

CONTAINS

  SUBROUTINE def_field_no_mask( &
    comp_id, cell_point_id, timestepstring, &
    field_name, collection_size, field_id)

    INTEGER, INTENT(IN) :: comp_id
    INTEGER, INTENT(IN) :: cell_point_id
    CHARACTER(LEN=*), INTENT(IN) :: timestepstring
    CHARACTER(LEN=*), INTENT(IN) :: field_name
    INTEGER, INTENT(IN) :: collection_size
    INTEGER, INTENT(OUT) :: field_id

#if !defined NOMPI && defined YAC_coupling
    CALL yac_fdef_field (                           &
      & field_name      = TRIM(field_name),         & !in
      & component_id    = comp_id,                  & !in
      & point_ids       = (/cell_point_id/),        & !in
      & num_pointsets   = 1,                        & !in
      & collection_size = collection_size,          & !in
      & timestep        = timestepstring,           & !in
      & time_unit       = YAC_TIME_UNIT_ISO_FORMAT, & !in
      & field_id        = field_id )                  !out
#endif

  END SUBROUTINE def_field_no_mask

  SUBROUTINE def_field_mask( &
    comp_id, cell_point_id, cell_mask_id, timestepstring, &
    field_name, collection_size, field_id)

    INTEGER, INTENT(IN) :: comp_id
    INTEGER, INTENT(IN) :: cell_point_id
    INTEGER, INTENT(IN) :: cell_mask_id
    CHARACTER(LEN=*), INTENT(IN) :: timestepstring
    CHARACTER(LEN=*), INTENT(IN) :: field_name
    INTEGER, INTENT(IN) :: collection_size
    INTEGER, INTENT(OUT) :: field_id

#if !defined NOMPI && defined YAC_coupling
    CALL yac_fdef_field_mask (                      &
      & field_name      = TRIM(field_name),         & !in
      & component_id    = comp_id,                  & !in
      & point_ids       = (/cell_point_id/),        & !in
      & mask_ids        = (/cell_mask_id/),         & !in
      & num_pointsets   = 1,                        & !in
      & collection_size = collection_size,          & !in
      & timestep        = timestepstring,           & !in
      & time_unit       = YAC_TIME_UNIT_ISO_FORMAT, & !in
      & field_id        = field_id )                  !out
#endif

  END SUBROUTINE def_field_mask

END MODULE mo_coupling_utils
