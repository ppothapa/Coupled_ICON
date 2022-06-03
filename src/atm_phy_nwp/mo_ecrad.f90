!>
!! This module is the interface between ICON and ECMWFs radiation code
!! ecRad which is provided as a library. 
!! - Modules, variables and data types from ecRad are used and provided
!!   to ICON. Possibly ambiguous names get an ecrad_ prefix.
!! - By this approach, no other ICON module has to make a direct use
!!   statement on an ecRad module.
!! - This module also holds a few ecRad configuration objects:
!!   ecrad_conf, nweight_par_ecrad, iband_par_ecrad, weight_par_ecrad
!!
!! Literature references:
!! Coddington et al. (2016) - Coddington, O., Lean, J. L., Pilewskie, P., Snow, M., & Lindholm, D. (2016).
!!                            A solar irradiance climate data record.
!!                            Bulletin of the American Meteorological Society, 97(7), 1265-1282.
!!
!! @author Daniel Rieger, DWD, Offenbach
!!
!! @par Revision History
!! Initial release by Daniel Rieger, DWD, Offenbach (2018-06-01)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_ecrad

  USE mo_kind,                    ONLY: wp
#ifdef __ECRAD
  USE radiation_config,           ONLY: t_ecrad_conf=>config_type,                        &
                                    &   ISolverHomogeneous, ISolverMcICA,                 &
                                    &   ISolverSpartacus, ISolverTripleclouds,            &
                                    &   ISolverCloudless,                                 &
                                    &   IGasModelMonochromatic, IGasModelIFSRRTMG,        &
                                    &   ILiquidModelMonochromatic, ILiquidModelSlingo,    &
                                    &   ILiquidModelSOCRATES,                             &
                                    &   IOverlapMaximumRandom, IOverlapExponentialRandom, &
                                    &   IOverlapExponential,                              &
                                    &   IIceModelMonochromatic,                           &
                                    &   IIceModelFu, IIceModelBaran, IIceModelBaran2016,  &
                                    &   IIceModelBaran2017, IIceModelYi
  USE radiation_single_level,     ONLY: t_ecrad_single_level_type=>single_level_type
  USE radiation_thermodynamics,   ONLY: t_ecrad_thermodynamics_type=>thermodynamics_type
  USE radiation_gas,              ONLY: t_ecrad_gas_type=>gas_type,                       &
                                    &   IMassMixingRatio, IVolumeMixingRatio,             &
                                    &   ecRad_IH2O=>IH2O,     ecRad_ICO2=>ICO2,           &
                                    &   ecRad_IO3=>IO3,       ecRad_IN2O=>IN2O,           &
                                    &   ecRad_ICO=>ICO,       ecRad_ICH4=>ICH4,           &
                                    &   ecRad_IO2=>IO2,       ecRad_ICFC11=>ICFC11,       &
                                    &   ecRad_ICFC12=>ICFC12, ecRad_IHCFC22=>IHCFC22,     &
                                    &   ecRad_ICCl4=>ICCl4
  USE radiation_flux,             ONLY: t_ecrad_flux_type=>flux_type
  USE radiation_cloud,            ONLY: t_ecrad_cloud_type=>cloud_type
  USE radiation_aerosol,          ONLY: t_ecrad_aerosol_type=>aerosol_type
  
  USE radiation_interface,        ONLY: ecrad_setup=>setup_radiation,                     &
                                    &   ecrad_set_gas_units=>set_gas_units,               &
                                    &   ecrad=>radiation
#endif

  IMPLICIT NONE

  PRIVATE

#ifdef __ECRAD
! ecRad subroutines
  PUBLIC :: ecrad_setup, ecrad_set_gas_units, ecrad

! ecRad configuration types
  PUBLIC :: t_ecrad_conf
  PUBLIC :: t_ecrad_single_level_type
  PUBLIC :: t_ecrad_thermodynamics_type
  PUBLIC :: t_ecrad_gas_type
  PUBLIC :: t_ecrad_flux_type
  PUBLIC :: t_ecrad_cloud_type
  PUBLIC :: t_ecrad_aerosol_type
! ecRad configuration state
  PUBLIC :: ecrad_conf

! Aerosol optical properties
  PUBLIC :: t_opt_ptrs

! ecRad enumerators
  ! Solver
  PUBLIC :: ISolverHomogeneous, ISolverMcICA, ISolverSpartacus, ISolverTripleclouds
  ! Gas model
  PUBLIC :: IGasModelMonochromatic, IGasModelIFSRRTMG
  ! Liquid hydrometeor scattering
  PUBLIC :: ILiquidModelMonochromatic, ILiquidModelSlingo, ILiquidModelSOCRATES
  ! Ice scattering
  PUBLIC :: IIceModelMonochromatic, IIceModelFu, IIceModelBaran, IIceModelBaran2016, IIceModelBaran2017, IIceModelYi
  ! Cloud overlap
  PUBLIC :: IOverlapMaximumRandom, IOverlapExponentialRandom, IOverlapExponential
  ! Gas units
  PUBLIC :: IMassMixingRatio, IVolumeMixingRatio
  ! Gas indices
  PUBLIC :: ecRad_IH2O, ecRad_ICO2, ecRad_IO3, ecRad_IN2O, ecRad_ICO, ecRad_ICH4
  PUBLIC :: ecRad_IO2, ecRad_ICFC11, ecRad_ICFC12, ecRad_IHCFC22, ecRad_ICCl4
  ! Photosynthetically active radiation weightings
  PUBLIC :: nweight_par_ecrad, iband_par_ecrad, weight_par_ecrad
  ! Spectral Solar Insolation
  PUBLIC :: ecrad_ssi_default, ecrad_ssi_coddington

! ----------------------------------------------------
! Configuration state

  TYPE(t_ecrad_conf) :: ecrad_conf

! Photosynthetically active radiation weightings
  INTEGER            :: nweight_par_ecrad
  INTEGER            :: iband_par_ecrad(100)
  REAL(KIND=wp)      :: weight_par_ecrad(100)

  !$ACC DECLARE COPYIN(iband_par_ecrad, weight_par_ecrad)

! Pointers to aerosol optical properties
  TYPE t_opt_ptrs
    REAL(wp), POINTER, DIMENSION(:,:) :: &
      &  ptr_od   => NULL(), &
      &  ptr_ssa  => NULL(), &
      &  ptr_g    => NULL()
    CONTAINS
! WARNING: Call finalize only if ptr_od, ptr_ssa, ptr_g are associated to a
!          target. If they were allocated, this might cause a memory leak
      PROCEDURE :: finalize => del_opt_ptrs
  END TYPE t_opt_ptrs

  REAL(wp) :: ecrad_ssi_default(14) = (/ 12.045647_wp, 20.257584_wp, 23.604472_wp , 22.309308_wp , 55.332985_wp , &
    &                                   102.388219_wp, 24.165380_wp, 343.917494_wp, 217.035256_wp, 345.359642_wp, &
    &                                   128.811472_wp, 49.887519_wp, 3.063681_wp  , 12.821340_wp /)
  ! New measurements of SSI from Coddington et al. 2016
  REAL(wp) :: ecrad_ssi_coddington(14) = (/ 12.045647_wp  , 20.257584_wp  , 23.604472_wp  , 23.37569292_wp, 57.56843759_wp, &
                                            105.6339255_wp, 24.72360028_wp, 345.7746485_wp, 213.5909065_wp, 344.8864993_wp, &
                                            128.6916773_wp, 45.19260459_wp, 2.825112161_wp, 12.821340_wp/)

CONTAINS

  SUBROUTINE del_opt_ptrs(self)
    CLASS(t_opt_ptrs),INTENT(inout) :: self

    IF (ASSOCIATED(self%ptr_od) ) THEN
      NULLIFY(self%ptr_od)
    ENDIF
    IF (ASSOCIATED(self%ptr_ssa) ) THEN
      NULLIFY(self%ptr_ssa)
    ENDIF
    IF (ASSOCIATED(self%ptr_g) ) THEN
      NULLIFY(self%ptr_g)
    ENDIF
  END SUBROUTINE del_opt_ptrs

#endif


END MODULE mo_ecrad
