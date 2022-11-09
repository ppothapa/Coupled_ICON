!>
!!
!! This module provides data structures for reading a NetCDF file in a distributed way.
!!
!! @par Revision History
!! Initial version by Roland Wirth (DWD), October 2022
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!
MODULE mo_read_netcdf_types

  USE mo_communication_types, ONLY: t_comm_pattern

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_distrib_read_data

  TYPE t_distrib_read_data
    INTEGER :: basic_data_index = -1
    CLASS(t_comm_pattern), POINTER :: pat => NULL()
  END TYPE t_distrib_read_data

END MODULE mo_read_netcdf_types