! Specification of vertical axes for the wave model
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

MODULE mo_waves_vertical_axes

  USE mo_kind,                              ONLY: dp
  USE mo_zaxis_type,                        ONLY: ZA_SURFACE, ZA_HEIGHT_10M
  USE mo_name_list_output_zaxes_types,      ONLY: t_verticalAxisList
  USE mo_name_list_output_zaxes,            ONLY: single_level_axis

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: setup_zaxes_waves

CONTAINS


  SUBROUTINE setup_zaxes_waves(verticalAxisList)
    TYPE(t_verticalAxisList), INTENT(INOUT) :: verticalAxisList

    ! --------------------------------------------------------------------------------------
    ! Definitions for single levels --------------------------------------------------------
    ! --------------------------------------------------------------------------------------

    ! surface level
    CALL verticalAxisList%append(single_level_axis(ZA_surface))

    ! Specified height level above ground: 10m
    CALL verticalAxisList%append(single_level_axis(ZA_height_10m, opt_level_value=10._dp))

  END SUBROUTINE setup_zaxes_waves

END MODULE mo_waves_vertical_axes

