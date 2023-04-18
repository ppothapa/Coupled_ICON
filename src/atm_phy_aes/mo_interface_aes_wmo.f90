!>
!! @brief Subroutine aes_phy_main calls all the parameterization schemes
!!
!! @author Hui Wan, MPI-M
!! @author Marco Giorgetta, MPI-M
!!
!! @par Revision History
!!  Original version from ECHAM6 (revision 2028)
!!  Modified for ICOHAM by Hui Wan and Marco Giorgetta (2010)
!!  Modified for ICONAM by Marco Giorgetta (2014)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

MODULE mo_interface_aes_wmo

  USE mo_aes_phy_dims        ,ONLY: aes_phy_dims
  USE mo_aes_phy_memory      ,ONLY: t_aes_phy_field, prm_field

  USE mo_timer               ,ONLY: ltimer, timer_start, timer_stop, timer_wmo

  USE mo_tropopause          ,ONLY: WMO_tropopause

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: interface_aes_wmo

CONTAINS

  SUBROUTINE interface_aes_wmo(jg, jb, jcs, jce)

    INTEGER, INTENT(in)     :: jg, jb, jcs, jce

    ! Local variables
    !
    INTEGER  :: nlev
    INTEGER  :: nproma
    !
    TYPE(t_aes_phy_field), POINTER :: field

    IF (ltimer) call timer_start(timer_wmo)

    nlev   = aes_phy_dims(jg)%nlev
    nproma = aes_phy_dims(jg)%nproma

    ! associate pointers
    field  => prm_field(jg)

    CALL WMO_tropopause( jg,                       &! in
                       & jcs, jce, nproma, nlev,   &! in
                       & field% ta(:,:,jb),        &! in
                       & field% pfull(:,:,jb),     &! in
                       & field% ptp(:,jb)          )! inout for diagnostics
    !

    IF (ltimer) call timer_stop(timer_wmo)

  END SUBROUTINE interface_aes_wmo
  !-------------------------------------------------------------------

END MODULE mo_interface_aes_wmo
