!>
!! Provides the interface for ART-routines to ICON parts
!!
!! This module provides an interface to ICON parts.
!! The interface is written in such a way, that ICON will compile and run
!! properly, even if the ART-routines are not available at compile time.
!!
!!
!! @author Sven Werchner, KIT
!!
!! @par Revision History
!! Initial revision by Sven Werchner, KIT (2022-09-27)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_art_general_interface
  
  USE mo_cdi,                           ONLY: DATATYPE_FLT32, DATATYPE_FLT64
  USE mo_io_config,                     ONLY: lnetcdf_flt64_output

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: routine = 'mo_art_general_interface'

  PUBLIC  :: getNetcdfPrecision 

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
FUNCTION getNetcdfPrecision() RESULT(flt_prec)

  INTEGER :: flt_prec

  flt_prec = MERGE(DATATYPE_FLT64,DATATYPE_FLT32,lnetcdf_flt64_output)

END FUNCTION getNetcdfPrecision 
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_general_interface
