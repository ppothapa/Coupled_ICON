!>
!! Read configuration parameters as Fortran namelist from an external file. 
!!
!! @author Marco, Giorgetta, MPI-M, 2017-12
!!
!! @par Revision History
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_aes_vdf_nml

  USE mo_aes_vdf_config   ,ONLY: aes_vdf_config, init_aes_vdf_config
  USE mo_process_nml      ,ONLY: process_nml
  
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: process_aes_vdf_nml

  NAMELIST /aes_vdf_nml/ aes_vdf_config

CONTAINS

  SUBROUTINE nml_read(funit)
    INTEGER, INTENT(in) :: funit
    READ(funit, NML=aes_vdf_nml)
  END SUBROUTINE nml_read
  !
  SUBROUTINE nml_write(funit)
    INTEGER, INTENT(in) :: funit
    WRITE(funit, NML=aes_vdf_nml)
  END SUBROUTINE nml_write
  !
  SUBROUTINE process_aes_vdf_nml(filename)
    !
    CHARACTER(LEN=*), INTENT(in) :: filename
    !
    CALL init_aes_vdf_config
    !
    CALL process_nml(filename, 'aes_vdf_nml', nml_read, nml_write)
    !$ACC ENTER DATA COPYIN(aes_vdf_config)
    !
  END SUBROUTINE process_aes_vdf_nml

END MODULE mo_aes_vdf_nml
