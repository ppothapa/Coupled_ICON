!>
!! @brief Subroutine interface_aes_art calls the art reaction interface.
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

MODULE mo_interface_aes_art

  USE mo_model_domain           ,ONLY: t_patch

  USE mo_kind                   ,ONLY: wp
  USE mtime                     ,ONLY: t_datetime => datetime

  USE mo_aes_phy_config         ,ONLY: aes_phy_tc
  USE mo_aes_phy_memory         ,ONLY: t_aes_phy_field, prm_field

  USE mo_nonhydro_state         ,ONLY: p_nh_state_lists
  USE mo_dynamics_config        ,ONLY: nnew_rcf

#ifdef __ICON_ART
  USE mo_art_reaction_interface ,ONLY: art_reaction_interface
#endif

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: interface_aes_art

CONTAINS

  SUBROUTINE interface_aes_art(patch)

    ! Arguments
    !
    TYPE(t_patch), TARGET, INTENT(in) :: patch

    ! Pointers
    !
    TYPE(t_aes_phy_field), POINTER    :: field

    ! Local variables
    !
    TYPE(t_datetime), POINTER :: datetime
    REAL(wp) :: pdtime
    LOGICAL  :: is_in_sd_ed_interval
    LOGICAL  :: is_active
    !
    INTEGER  :: jg

    jg  = patch%id

    datetime             => aes_phy_tc(jg)%datetime
    pdtime               =  aes_phy_tc(jg)%dt_phy_sec
    is_in_sd_ed_interval =  aes_phy_tc(jg)%is_in_sd_ed_interval_art
    is_active            =  aes_phy_tc(jg)%is_active_art

    ! associate pointers
    field     => prm_field(jg)

    IF ( is_in_sd_ed_interval ) THEN
       !
       IF ( is_active ) THEN
          !
#ifdef __ICON_ART
          CALL art_reaction_interface(jg,                                           & !> in
               &                      datetime,                                     & !> in
               &                      pdtime,                                       & !> in
               &                      p_nh_state_lists(jg)%prog_list(nnew_rcf(jg)), & !> in
               &                      field%qtrc_phy)
#endif
          !
       END IF
       !
    END IF

    ! disassociate pointers
    NULLIFY(datetime, field)

  END SUBROUTINE interface_aes_art

END MODULE mo_interface_aes_art
