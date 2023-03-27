!>
!! Specification of vertical axes for the wave model
!!
!!
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2023-01-24)
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

