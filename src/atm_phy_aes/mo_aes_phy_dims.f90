!>
!! Dimensions for the AES physics package.
!!
!! @author Marco Giorgetta, MPI-M
!!
!! @par Revision History
!!     First version by Marco Giorgetta, MPI-M (2023-01)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_aes_phy_dims

  USE mo_impl_constants,  ONLY: max_dom
  USE mo_parallel_config, ONLY: nproma
  USE mo_run_config,      ONLY: num_lev, ntracer

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: aes_phy_dims, init_aes_phy_dims

  !>
  !! Dimensions for physics, for multiple domains/grids.
  !!
  TYPE t_aes_phy_dims
     !
     INTEGER :: nproma  !< size of cells dimension
     INTEGER :: nlev    !< size of levels dimension
     INTEGER :: ntracer !< size of tracers dimension
     !
  END type t_aes_phy_dims

  TYPE(t_aes_phy_dims)  , TARGET :: aes_phy_dims   (max_dom)

CONTAINS

  !----

  !>
  !! Initialize the dimensions
  !!
  SUBROUTINE init_aes_phy_dims
    !
    aes_phy_dims(:)%nproma  = nproma
    aes_phy_dims(:)%ntracer = ntracer
    aes_phy_dims(:)%nlev    = num_lev(:)
    !
  END SUBROUTINE init_aes_phy_dims

  !----

END MODULE mo_aes_phy_dims
