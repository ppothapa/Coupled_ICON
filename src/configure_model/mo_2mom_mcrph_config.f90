!>
!! @brief Setup for 2-moment cloud microphysics scheme
!!
!! configuration setup for synthetic radar data on the model grid
!!
!! @author Ulrich Blahak, DWD
!!
!!
!! @par Revision History
!! Initial revision by Ulrich Blahak, DWD (2022-07-06)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_2mom_mcrph_config
  
  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PUBLIC

  !--------------------------------------------------------------------------
  ! Namelist parameters
  !--------------------------------------------------------------------------

  ! Container for some configuration parameters of the Seifert-Beheng 2-moment cloud microphysical scheme,
  ! which may be changed by namelist:
  TYPE t_cfg_2mom
    INTEGER  :: i2mom_solver ! 0) explicit (1) semi-implicit solver
    INTEGER  :: ccn_type     ! if not set by namelist, the ccn_type_gscp4 or ccn_type_gscp5 will win
    REAL(wp) :: alpha_spacefilling  !..factor involved in the conversion of ice/snow to graupel by riming
    REAL(wp) :: D_conv_ii    ! D-threshold for conversion to snow ice_selfcollection 
    REAL(wp) :: D_rainfrz_ig ! rain --> ice oder graupel
    REAL(wp) :: D_rainfrz_gh ! rain --> graupel oder hail
    REAL(wp) :: rain_cmu0    ! asymptotic mue-value for small D_m in the mu-Dm-Relation of Seifert (2008)
    REAL(wp) :: rain_cmu1    ! asymptotic mue-value for large D_m in the mu-Dm-Relation of Seifert (2008)
    REAL(wp) :: rain_cmu3    ! D_br: equilibrium diameter for breakup and selfcollection
  END TYPE t_cfg_2mom


END MODULE mo_2mom_mcrph_config
