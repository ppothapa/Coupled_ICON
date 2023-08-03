!>
!! @brief Subroutine interface_aes_cov calls the cloud cover scheme.
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

MODULE mo_interface_aes_cov

  USE mo_kind                ,ONLY: wp

  USE mo_aes_phy_dims        ,ONLY: aes_phy_dims
  USE mo_aes_phy_config      ,ONLY: aes_phy_config
  USE mo_aes_phy_memory      ,ONLY: t_aes_phy_field, prm_field
  USE mo_aes_sfc_indices     ,ONLY: nsfc_type, iwtr, iice
  
  USE mo_timer               ,ONLY: ltimer, timer_start, timer_stop, timer_cov

  USE mo_run_config          ,ONLY: iqv, iqc, iqi
  USE mo_cover               ,ONLY: cover
  
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: interface_aes_cov

CONTAINS

  SUBROUTINE interface_aes_cov(jg, jb, jcs, jce) 

    ! Arguments
    !
    INTEGER, INTENT(in)     :: jg, jb, jcs, jce

    ! Pointers
    !
    TYPE(t_aes_phy_field), POINTER :: field

    ! Local variables
    !
    INTEGER  :: nlev
    INTEGER  :: nproma
    !
    INTEGER  :: nlevp1, jc, jks
    REAL(wp) :: zfrw (aes_phy_dims(jg)%nproma) !< cell area fraction of open water
    REAL(wp) :: zfri (aes_phy_dims(jg)%nproma) !< cell area fraction of ice covered water

    IF (ltimer) call timer_start(timer_cov)

    nlev   = aes_phy_dims(jg)%nlev
    nproma = aes_phy_dims(jg)%nproma

    field  => prm_field(jg)

    jks    = aes_phy_config(jg)%jks_cloudy
    nlevp1 = nlev+1

    !$ACC DATA PRESENT(field%frac_tile) &
    !$ACC   CREATE(zfrw, zfri)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR
    DO jc = 1, nproma
      IF (iwtr.LE.nsfc_type) THEN
         zfrw(jc) = field%frac_tile(jc,jb,iwtr)
      ELSE
         zfrw(jc) = 0.0_wp
      END IF

      IF (iice.LE.nsfc_type) THEN
         zfri(jc) = field%frac_tile(jc,jb,iice)
      ELSE
         zfri(jc) = 0.0_wp
      END IF
    END DO
    !$ACC END PARALLEL

    CALL cover( jg,                            &! in
         &      jb,                            &! in
         &      jcs, jce, nproma,              &! in
         &      jks, nlev, nlevp1,             &! in
         &      zfrw(:),                       &! in
         &      zfri(:),                       &! in
         &      field% zf(:,:,jb),             &! in
         &      field% phalf(:,:,jb),          &! in
         &      field% pfull(:,:,jb),          &! in
         &      field%  ta(:,:,jb),            &! in    tm1
         &      field%  qtrc_phy(:,:,jb,iqv),  &! in    qm1
         &      field%  qtrc_phy(:,:,jb,iqc),  &! in    xlm1
         &      field%  qtrc_phy(:,:,jb,iqi),  &! in    xim1
         &      field%  aclc(:,:,jb)          ) ! out   (for "radiation" and "vdiff_down")

    !$ACC WAIT(1)
    !$ACC END DATA

    NULLIFY(field)

    IF (ltimer) call timer_stop(timer_cov)

  END SUBROUTINE interface_aes_cov

END MODULE mo_interface_aes_cov
