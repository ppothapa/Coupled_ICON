!>
!! @author Hui Wan, MPI-M
!!
!! @par Revision History
!! First version by Hui Wan, 2010-08-16
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_aes_phy_cleanup

  USE mtime,                       ONLY: OPERATOR(>)
  USE mo_grid_config,              ONLY: n_dom
  USE mo_aes_phy_memory,           ONLY: destruct_aes_phy_memory
  USE mo_cloud_mig_memory,         ONLY: destruct_cloud_mig_memory
  USE mo_radiation_forcing_memory, ONLY: destruct_radiation_forcing_list
  USE mo_aes_phy_config,           ONLY: aes_phy_tc, dt_zero
  USE mo_turb_vdiff,               ONLY: vdiff_cleanup

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: cleanup_aes_phy

CONTAINS
  !>
  !! Top-level routine for the cleanup for ECHAM6 physics.
  !! It calls a series of subroutines to deallocate parameter arrays with
  !! "allocatable" attribute.
  !!
  !! @par Revision History
  !! Initial version by Hui Wan, MPI-M (2010-08)
  !!
  SUBROUTINE cleanup_aes_phy

    LOGICAL :: lany
    INTEGER :: jg

    lany=.FALSE.
    DO jg = 1,n_dom
       lany = lany .OR. (aes_phy_tc(jg)%dt_vdf > dt_zero)
    END DO
    IF (lany) CALL vdiff_cleanup             ! deallocate array "matrix_idx"

    lany=.FALSE.
    DO jg = 1,n_dom
       lany = lany .OR. (aes_phy_tc(jg)%dt_rad > dt_zero)
    END DO
    IF (lany) CALL destruct_radiation_forcing_list
   
    lany=.FALSE.
    DO jg = 1,n_dom
       lany = lany .OR. (aes_phy_tc(jg)%dt_mig > dt_zero)
    END DO
    IF (lany) CALL destruct_cloud_mig_memory

    CALL destruct_aes_phy_memory

  END SUBROUTINE cleanup_aes_phy
  !-------------

END MODULE mo_aes_phy_cleanup

