!>
!! Configuration structure for the VDIFF turbulence scheme.
!!
!! @author Marco Giorgetta, MPI-M
!! @author Roland Wirth, DWD
!!
!! @par Revision History
!! 2017-04, Marco Giorgetta (MPI-M), First revision
!! 2021-07, Roland Wirth (DWD), Adapted from ICON-ECHAM
!!
!! References:
!!     Angevine, W. M., Jiang, H., & Mauritsen T. (2010).
!!           Performance of an eddy diffusivity mass flux scheme for shallow cumulus boundary layers.
!!           Monthly Weather Review, 138(7), 2895-2912. https://doi.org/10.1175/2010MWR3142.1
!!     Mauritsen, T., & Svensson, G. (2007).
!!           Observations of stably stratified shear-driven atmospheric turbulence at low and high Richardson numbers.
!!           Journal of the Atmospheric Sciences, 64(2), 645-655. https://doi.org/10.1175/JAS3856.1
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_turb_vdiff_config
  USE mo_kind, ONLY: wp
  USE mo_exception, ONLY: message, finish
  USE mo_impl_constants, ONLY: max_dom
  USE mo_turb_vdiff_params, ONLY: VDIFF_TURB_TTE

  IMPLICIT NONE
  PRIVATE

  !> Configuration for the VDIFF TTE and 3D Smagorinsky schemes.
  !! The TTE scheme is described in Mauritsen et al. J Atm Sci 64, 4113 (2007).
  TYPE t_vdiff_config
    LOGICAL  :: lsfc_mom_flux   !< switch on/off surface momentum flux
    LOGICAL  :: lsfc_heat_flux  !< switch on/off surface heat flux (sensible AND latent)

    REAL(wp) :: pr0             !< neutral limit Prandtl number, can be varied from about 0.6 to 1.0
    REAL(wp) :: f_tau0          !< neutral non-dimensional stress factor
    REAL(wp) :: f_theta0        !< neutral non-dimensional heat flux factor

    REAL(wp) :: f_tau_limit_fraction   !< Limit of f_tau/f_tau0 in the stable limit (Mauritsen: 0.25).
    REAL(wp) :: f_theta_limit_fraction !< Limit of f_theta/f_theta0 in the stable limit (Mauritsen: 0.).

    REAL(wp) :: f_tau_decay     !< Decay rate of f_tau/f_tau0 towards its limit (Mauritsen: 4.).
    REAL(wp) :: f_theta_decay   !< Decay rate of f_theta/f_theta0 towards its limit (Mauritsen: 4.).

    REAL(wp) :: ek_ep_ratio_stable   !< Ek/Ep ratio in the stable limit (Mauritsen: 1/(0.3+-0.1) - 1).
    REAL(wp) :: ek_ep_ratio_unstable !< Ek/Ep ratio in the unstable limit (Mauritsen: 1).

    REAL(wp) :: c_f             !< mixing length: coriolis term tuning parameter
    REAL(wp) :: c_n             !< mixing length: stability term tuning parameter
    REAL(wp) :: c_e             !< dissipation coefficient (=f_tau0^(3/2))
    REAL(wp) :: wmc             !< ratio of typical horizontal velocity to wstar at free convection
    REAL(wp) :: fsl             !< fraction of first-level height at which surface fluxes
    !                              are nominally evaluated, tuning param for sfc stress
    REAL(wp) :: fbl             !< 1/fbl: fraction of BL height at which lmix hat its max
    REAL(wp) :: lmix_max        !< Maximum mixing length [m].
    REAL(wp) :: z0m_min         !< Minimum roughness length for momentum [m].
    REAL(wp) :: z0m_ice         !< Roughness length for momentum over ice [m].
    REAL(wp) :: z0m_oce         !< Roughness length for momentum over ocean [m].
                                !! See mo_surface_ocean.f90 of ECHAM6.

    INTEGER :: turb             !< turbulence scheme: VDIFF_TURB_TTE or VDIFF_TURB_3DSMAGORINSKY.
    REAL(wp) :: smag_constant
    REAL(wp) :: turb_prandtl    !< Turbulent Prandtl number
    REAL(wp) :: rturb_prandtl   !< inverse turbulent prandtl number
    REAL(wp) :: km_min          !< min mass weighted turbulent viscosity
    REAL(wp) :: max_turb_scale  !< max turbulence length scale
    REAL(wp) :: min_sfc_wind    !< min sfc wind in free convection limit
  END TYPE t_vdiff_config

  !---------------------
  ! Types and constants
  !---------------------

  PUBLIC :: t_vdiff_config

  !------------
  ! Procedures
  !------------

  INTERFACE vdiff_config_check
    MODULE PROCEDURE vdiff_config_check_0d
    MODULE PROCEDURE vdiff_config_check_1d
  END INTERFACE vdiff_config_check

  PUBLIC :: vdiff_config_init
  PUBLIC :: vdiff_config_update
  PUBLIC :: vdiff_config_check

  !-----------
  ! Variables
  !-----------

  !> Config variables for each domain.
  TYPE(t_vdiff_config) :: vdiff_config(max_dom)
  PUBLIC :: vdiff_config

CONTAINS

  !>
  !! Initialize the configuration state structure.
  !!
  ELEMENTAL SUBROUTINE vdiff_config_init(config)
    TYPE(t_vdiff_config), INTENT(OUT) :: config

    config%lsfc_mom_flux  = .TRUE.
    config%lsfc_heat_flux = .TRUE.

    config%pr0      =  1.0_wp
    config%f_tau0   =  0.17_wp

    config%f_tau_decay = 4._wp
    config%f_theta_decay = 4._wp
    config%f_tau_limit_fraction = 0.25_wp
    config%f_theta_limit_fraction = 0._wp

    config%ek_ep_ratio_stable = 3._wp
    config%ek_ep_ratio_unstable = 2._wp

    config%c_f      =  0.185_wp
    config%c_n      =  2.0_wp
    config%wmc      =  0.5_wp
    config%fsl      =  0.4_wp
    config%fbl      =  3._wp
    config%z0m_min  =  1.5e-5_wp
    config%z0m_ice  =  1e-3_wp
    config%z0m_oce  =  1e-3_wp
    config%lmix_max =  150._wp
    config%turb     =  VDIFF_TURB_TTE
    config%smag_constant = 0.23_wp
    config%max_turb_scale= 300._wp
    config%turb_prandtl  = 0.33333333333_wp
    config%km_min        =  0.001_wp
    config%min_sfc_wind  =  1._wp

    CALL vdiff_config_update(config)

  END SUBROUTINE vdiff_config_init

  !>
  !! Update derived parameters.
  !!
  ELEMENTAL SUBROUTINE vdiff_config_update(config)
    TYPE(t_vdiff_config), INTENT(INOUT) :: config

    config%f_theta0 = -SQRT(config%f_tau0**2/2.0_wp/config%pr0)
    config%c_e      =  SQRT(config%f_tau0**3)
    config%rturb_prandtl = 1/config%turb_prandtl
  END SUBROUTINE vdiff_config_update

  !<
  !! Check validity of configuration (scalar version).
  !!
  !! Prints a warning or aborts if chosen settings are not recommended or unsupported.
  !! Use via the generic interface mo_turb_vdiff_config::vdiff_config_check.
  SUBROUTINE vdiff_config_check_0d(config)
    TYPE(t_vdiff_config), INTENT(IN) :: config

    CHARACTER(LEN=*), PARAMETER :: routine = 'vdiff_config_check'

    ! check range of pr0
    IF ( config%pr0 < 0.6_wp ) CALL message(routine//' WARNING:',' config%pr0 < 0.6 ')
    IF ( config%pr0 > 1.0_wp ) CALL message(routine//' WARNING:',' config%pr0 > 1.0 ')

  END SUBROUTINE vdiff_config_check_0d

  !<
  !! Check validity of configuration (array version).
  !!
  !! Prints a warning or aborts if chosen settings are not recommended or unsupported.
  !! Use via the generic interface mo_turb_vdiff_config::vdiff_config_check.
  SUBROUTINE vdiff_config_check_1d(config)
    TYPE(t_vdiff_config), INTENT(IN) :: config(:)

    CHARACTER(LEN=*), PARAMETER :: routine = 'vdiff_config_check'

    ! check range of pr0
    IF ( ANY(config(:)%pr0 < 0.6_wp ) ) CALL message(routine//' WARNING:',' config(:)%pr0 < 0.6 ')
    IF ( ANY(config(:)%pr0 > 1.0_wp ) ) CALL message(routine//' WARNING:',' config(:)%pr0 > 1.0 ')

  END SUBROUTINE vdiff_config_check_1d

END MODULE mo_turb_vdiff_config
