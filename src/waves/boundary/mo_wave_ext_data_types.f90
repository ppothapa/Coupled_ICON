!>
!! Data type defintion for wave external data state
!!
!! Defines the data type for storing wave-specific external parameter 
!! fields such as bathymetry.
!!
!! @author Mikhail Dobrynin, DWD
!!
!!
!! @par Revision History
!! Initial Revision by Mikhail Dobrynin, DWD (2023-01-16)
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
MODULE mo_wave_ext_data_types

  USE mo_kind,               ONLY: wp

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: t_external_wave


  ! wave-specific external data type
  !
  TYPE :: t_external_wave

    ! ocean topography <=> bathymetric height used in the ocean 
    ! cell centers and edges only
    !
    REAL(wp), POINTER ::   &  !<  bathymetric height at cell centers [m]
      &  bathymetry_c(:,:)    !  index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::   &  !< topographic height at cell edges    [m]
      &  bathymetry_e(:,:)    !  index1=1,nproma, index2=1,nblks_e

    REAL(wp), POINTER ::   &  !< topographic height at cell vertices [m]
      &  bathymetry_v(:,:)    !  index1=1,nproma, index2=1,nblks_v

    ! *** Land-Sea-Mask ***
    INTEGER, POINTER  ::   &  !< land-sea-mask for cell centers          [ ]
      &  lsm_wave_ctr_c(:,:)       !  index1=1,nproma, index2=1,nblks_c
    INTEGER, POINTER ::    &  !< land-sea-mask for cell edges
      &  lsm_wave_ctr_e(:,:)       !  index1=1,nproma, index2=1,nblks_e

  END TYPE t_external_wave

END MODULE mo_wave_ext_data_types
