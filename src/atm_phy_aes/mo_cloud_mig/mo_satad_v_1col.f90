!>
!! @brief Single column subroutine for saturation adjustment at constant volume
!!
!! @author Marco Giorgetta, MPI-M, 2020-01
!!
!! @par Revision History
!!
!! This code is derived from the satad_v_3d code of module mo_satad,
!! which processes 2-dimensional arrays (columns x levels) and originates from
!! the COSMO model, see comments in the original code.
!!
!! This code is designed for ICON only, to explore ways to
!! formulate phyiscs at the natural dimensionality.
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

MODULE mo_satad_v_1col

  USE mo_kind,               ONLY: wp

  USE mo_physical_constants, ONLY: rv    ,&  !> gas constant for water vapour
!!$       &                           cpv   ,&  !> specific heat of water vapor at constant pressure
       &                           cvd   ,&  !> specific heat of dry air at constant volume
       &                           clw   ,&  !> specific heat of water
       &                           alv   ,&  !> latent heat for vapourization
       &                           tmelt     !> melting temperature of ice/snow

  USE mo_satad,              ONLY: qsat_rho   ,&!> function for qsat at constant total density 
       &                           dqsatdt_rho  !> function for d(qsat_rho)/dt at const. tot. density

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: satad_v_1col

  ! Control parameters
  INTEGER , PARAMETER :: maxiter = 10        ! Max. number of iterations
  REAL(wp), PARAMETER :: tol     = 1.e-3_wp  ! Target temperature accuracy for iterations
  REAL(wp), PARAMETER :: zqwmin  = 1.e-20_wp ! Minimum cloud water content for adjustment

  ! Physical constant
  ! ?? why not cpv from mo_physical_constants ??
  REAL(wp), PARAMETER :: cpv     = 1850._wp  ! Specific heat of water vapor J at constant pressure
  ! (Landolt-Bornstein)

CONTAINS

  SUBROUTINE satad_v_1col (rho, ta, qv, qc)

    IMPLICIT NONE

    !-------------------------------------------------------------------------------
    !
    ! Description:
    !   This routine corrects the temperature (ta), the specific humidity (qv),
    !   and the cloud water content (qc) for condensation/evaporation.
    !   Pressure adapts itself in ICON.
    !
    ! Method:
    !   Saturation adjustment at constant volume  assuming chemical equilibrium
    !   of water and vapor.
    !   For the heat capacity  of the total system (dry air, vapor, and hydrometeors)
    !   the value of dry air is taken, which is a common approximation and introduces
    !   only a small error.
    !
    ! Subroutine arguments
    REAL(wp), INTENT(in)    :: rho (:) ! total density          [kg/m3]
    REAL(wp), INTENT(inout) :: ta  (:) ! temperature            [K]
    REAL(wp), INTENT(inout) :: qv  (:) ! specific vapor content [kg/kg]
    REAL(wp), INTENT(inout) :: qc  (:) ! specific cloud content [kg/kg]
    !
    !-------------------------------------------------------------------------------

    ! Local variables
    REAL(wp) :: alvocvd         ! (temperature-dependent) latent heat of vaporization over cv
    REAL(wp) :: twork, tworkold ! temperature
    REAL(wp) :: qw              ! specific vapor+water content
    REAL(wp) :: qsat            ! specific vapor content at saturation
    REAL(wp) :: ft, dft         ! Newton method
    INTEGER  :: count           ! iteration counter

    INTEGER  :: jk, jks, jke

    !-------------------------------------------------------------------------------

    jks = 1
    jke = SIZE(rho)

    DO jk = jks,jke

       ! total content of the species which are changed by the adjustment:
       qw = qv(jk) + qc(jk)

       ! check, which points will still be subsaturated even
       ! if all the cloud water would have been evaporated.
       ! At such points, the Newton iteration is not necessary and the
       ! adjusted values of T, p, qv and qc can be obtained directly.
       alvocvd = ( alv + (cpv-clw)*(ta(jk)-tmelt) - rv*ta(jk) )/cvd
       twork   = ta(jk) - alvocvd*qc(jk)
       qsat    = qsat_rho(twork,rho(jk))

       IF (qw <= qsat ) THEN

          ! In this case, all the cloud water evaporates and there is still (sub)saturation.
          ! The resulting state depends only on the available cloud water and is
          ! not saturated, which enables direct computation of the adjusted variables:
          ta(jk) = twork
          qv(jk) = qw
          qc(jk) = 0.0_wp

       ELSE

          ! In this case, the Newton interation method is needed.
          ! Field for the iterated temperature, here set the starting value for the
          ! iteration below to the "old" temperature.
          ! (As an alternative, the arithmetic mean between the "old" temperature and the
          ! dew point temperature has been tested, but did not significantly increase the
          ! convergence speed of the iteration):
          twork = ta(jk)
          !
          ! And this is the storage variable for the "old" values in the iteration below:
          ! Add some nonsense increment to the starting, which is sufficient to trigger the
          ! iteration below:
          tworkold = ta(jk) + 2.0_wp*tol
          !
          ! Newton iteration
          count = 0
          DO WHILE (ABS(twork-tworkold)>tol .AND. count<maxiter)
             ! Here we still have to iterate ...
             tworkold = twork
             qsat     = qsat_rho(twork,rho(jk))
             ! Newton:
             ft    = twork - ta(jk) + alvocvd*(qsat-qv(jk))
             dft   = 1.0_wp         + alvocvd*dqsatdt_rho(qsat,twork,rho(jk))
             twork = twork - ft/dft;
             count = count + 1
          END DO
          !
          ! We disregard here the extrapolation of qsat from the second-last iteration
          ! step, which is done in the original routine to exactly preserve the internal energy.
          ! This introduces a small error (see the COSMO-Documentation, Part III).
          qsat = qsat_rho(twork,rho(jk))
          !
          ta(jk) = twork
          qv(jk) = qsat
          qc(jk) = MAX(qw-qsat,zqwmin)

       END IF

    END DO

  END SUBROUTINE satad_v_1col

END MODULE mo_satad_v_1col

