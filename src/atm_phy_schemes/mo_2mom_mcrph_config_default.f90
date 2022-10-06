!>
!! @brief Setup for 2-moment cloud microphysics scheme
!!
!! default configuration setup for synthetic radar data on the model grid
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
MODULE mo_2mom_mcrph_config_default
  
  USE mo_kind, ONLY: wp
  USE mo_2mom_mcrph_config, ONLY: t_cfg_2mom

  IMPLICIT NONE
  PUBLIC

  !--------------------------------------------------------------------------
  ! Namelist parameters
  !--------------------------------------------------------------------------

  !.. Type instance to hold the defaults for the config params:
  TYPE(t_cfg_2mom), PARAMETER :: cfg_2mom_default = t_cfg_2mom ( &
       &            1, &          ! i2mom_solver: 0) explicit (1) semi-implicit solver
       &           -1, &          ! ccn_type: 6,7,8,9; if not set by namelist, the ccn_type_gscp4 or ccn_type_gscp5 will win
       &            0.01_wp, &    ! alpha_spacefilling
       &            75.0e-6_wp, & ! D-threshold for conversion to snow ice_selfcollection
       &            0.50e-3_wp, & ! D_rainfrz_ig
       &            1.25e-3_wp, & ! D_rainfrz_gh
       &            6.0_wp, &     ! rain_cmu0
       &            30.0_wp, &    ! rain_cmu1
       &            1.1e-3_wp &   ! rain_cmu3 = D_br
       &            )

END MODULE mo_2mom_mcrph_config_default
