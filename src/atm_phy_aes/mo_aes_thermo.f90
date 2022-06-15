!>
!>
!! @brief Module aes_thermo thermodyanmic routines for moist atmosphere
!!
!! @author Bjorn Stevens, MPI-M, 2020-01
!!
!! @par Copyright and License
!!
!!   Code distributed under BSD-3Clause license
!!
MODULE mo_aes_thermo

USE mo_kind,               ONLY: wp     , &
                                 i4
USE mo_physical_constants, ONLY: rv    , & !> gas constant for water vapour
                                 rd    , & !! rd
                                 vtmpc1, & !! rv/rd-1._wp
                                 cvd   , & !!
                                 cpv   , & !!
                                 cvv   , & !!
                                 clw   , & !! specific heat of water
                                 alv   , & !! latent heat of vaporization
                                 als   , & !! latent heat of sublimation
                                 tmelt

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: saturation_adjustment ! partitions water mass to maintain saturation
  PUBLIC  :: qsat_rho              ! sat. vapor pres. (over liquid) at constant density
  PUBLIC  :: qsat_ice_rho          ! sat. vapor pres. (over ice) at constant density
  PUBLIC  :: dqsatdT_rho           ! d(qsat_rho)/dT
  PUBLIC  :: vaporization_energy   ! internal energy of vaporization
  PUBLIC  :: sublimation_energy    ! internal energy of sublimation
  PUBLIC  :: sat_pres_water        ! saturation pressure over water
  
  REAL (KIND=wp), PARAMETER ::     &
       ci  = 2108.0_wp,            & !! specific heat of ice
       lvc = alv-(cpv-clw)*tmelt,  & !! invariant part of vaporization enthalpy    
       c1es  = 610.78_wp,          & !! constants for saturation vapor pressure
       c2es  = c1es*rd/rv,         & !!
       c3les = 17.269_wp,          & !!
       c3ies = 21.875_wp,          & !!
       c4les = 35.86_wp,           & !!
       c4ies = 7.66_wp,            & !!
       c5les = c3les*(tmelt-c4les),& !!
       c5ies = c3ies*(tmelt-c4ies)

CONTAINS

SUBROUTINE saturation_adjustment ( idim, kdim,  ilo,  iup,  klo,  kup, & 
                                     te,  qve,  qce,  qre,  qti,  rho  ) 
  !-------------------------------------------------------------------------------
  !
  ! Description:
  !   This routine performs the saturation adjustment to find the combination
  !   of cloud water and temperature that is in equilibirum at the same internal
  !   energy as the initial fields.
  !
  ! Method:
  !   Saturation adjustment in the presence of non-zero cloud water requires
  !   solving a non-linear equation, which is done using a Newton-Raphson
  !   method.  The procedure first checks for the special case of sub-saturation
  !   in which case the solution can be calculated directly.   If not then the
  !   the solver looks for the zero of the function f(T) denoted fT, whereby T
  !   is temperature and f is the difference between the internal energy at T
  !   and the internal energy at the initial T and qc, denoted ue.
  !
  !-------------------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER (KIND=i4), INTENT (IN) ::  &
       idim, kdim,              & !  Dimension of I/O-fields
       ilo, iup, klo, kup         !  start- and end-indices for the computations

  REAL    (KIND=wp),    INTENT (INOUT), DIMENSION(:,:) ::  &  !  dim (idim,kdim)
       te      , & ! temperature on input/ouput
       qve     , & ! specific humidity on input/output
       qce         ! specific cloud water content on input/output

  REAL    (KIND=wp),    INTENT (IN   ), DIMENSION(:,:) ::  &  !  dim (idim,kdim)
       qre     , & ! specific rain water
       qti     , & ! specific mass of all ice species (total-ice)
       rho         ! density containing dry air and water constituents

  INTEGER (KIND=i4) ::  &
       i, k, iter     ! loop indices

  REAL    (KIND=wp   ) ::  &
       Tx,  & ! Test temperature variable
       qx,  & ! Test saturation vapor mixing ratio (at Tx)
       qcx, & ! Test cloud water, can be negative
       qt,  & ! Total water specific humidity
       cvc, & ! contribution to cv that is constant (not varying with fast condensation)
       cv,  & ! isometric specific heat of moist system with condensate
       ue,  & ! partial (that which varies with condensation) internal energy at qce and Te
       ux,  & ! partial (that which varies with condensation) internal energy at qcx and Tx
       dqx, & ! change in saturation vapor pressure at Tx 
       dux    ! derivative of ux wrt Tx

  !------------ End of header ----------------------------------------------------

!!!=============================================================================================

  !$ACC PARALLEL DEFAULT(PRESENT)
  !$ACC LOOP GANG VECTOR TILE(128,1)
  DO k = klo, kup
    DO i = ilo , iup
      qt      = qve(i,k) + qce(i,k) + qre(i,k) + qti(i,k)
      cvc     = cvd*(1.0_wp-qt) + clw*qre(i,k) + ci*qti(i,k)

      cv      = cvc + cvv*qve(i,k) + clw*qce(i,k)
      ue      = cv*te(i,k) - qce(i,k)*lvc
      Tx      = ue / (cv + qce(i,k)*(cvv-clw))
      qx      = qsat_rho(Tx, rho(i,k))
      !
      ! If subsaturated upon evaporating all cloud water, T can be diagnosed explicitly,
      ! so test for this.  If not then T needs to be solved for iteratively.
      !
      IF (qve(i,k)+qce(i,k) <= qx ) THEN 
        qve(i,k)  = qve(i,k)+qce(i,k)
        qce(i,k)  = 0.0_wp
      ELSE 
      !$ACC LOOP SEQ
        DO iter = 1, 6 
           qx   = qsat_rho(Tx, rho(i,k))
           dqx  = dqsatdT_rho(qx, Tx, rho(i,k) )
           qcx  = qve(i,k)+qce(i,k) - qx
           cv   = cvc + cvv*qx + clw*qcx 
           ux   = cv*Tx -qcx*lvc
           dux  = cv + dqx*(lvc + (cvv-clw)*Tx)
           Tx   = Tx - (ux-ue) / dux 
        END DO 
        qx       = qsat_rho(Tx, rho(i,k))
        qce(i,k) = MAX( qve(i,k) + qce(i,k) - qx, 0.0_wp)
        qve(i,k) = qx
      ENDIF
      te (i,k) = Tx
    ENDDO !i
  ENDDO !k
  !$ACC END PARALLEL

!!!=============================================================================================

END SUBROUTINE saturation_adjustment

  !-------------------------------------------------------------------------------
  !
  ! Description:
  !   Below are thermodynamic functions, the take SI units (Kelvin) for
  !   temperature or unitless for specific masses. They include
  !      * sat_pres_water :saturation pressure over planar liquid
  !      * sat_pres_icei  :saturation pressure over planar ic
  !      * qsat_rho            ! sat. vapor pressure (over liquid) at constant density
  !      * qsat_ice_rho        ! sat. vapor pressure (over ice) at constant density
  !      * dqsatdT_rho         ! d(qsat_rho)/dT
  !      * vaporization_energy ! internal energy of vaporization
  !      * sublimation_energy  ! internal energy of sublimation
  
  !
  ! Method (also GPU directives):
  !   Most functions are elemental.  However functions that use the directive 
  !   ACC ROUTINE SEQ are conditionally elemental as this directive is not
  !   compatible with an elemental function
  !
  !-------------------------------------------------------------------------------

#ifndef _OPENACC
ELEMENTAL &
#endif
FUNCTION sat_pres_water(TK)

  IMPLICIT NONE
  REAL (KIND=wp)              :: sat_pres_water
  REAL (KIND=wp), INTENT(IN)  :: TK

  !$ACC ROUTINE SEQ
  sat_pres_water = c1es*EXP( c3les*(TK-tmelt)/(TK-c4les) )

END FUNCTION sat_pres_water

!!!=============================================================================================

#ifndef _OPENACC
ELEMENTAL &
#endif
FUNCTION sat_pres_ice(TK)

  IMPLICIT NONE
  REAL (KIND=wp)              :: sat_pres_ice
  REAL (KIND=wp), INTENT(IN)  :: TK

  !$ACC ROUTINE SEQ
  sat_pres_ice  = c1es*EXP( c3ies*(TK-tmelt)/(TK-c4ies) )

END FUNCTION sat_pres_ice

!!!=============================================================================================

#ifndef _OPENACC
ELEMENTAL &
#endif
FUNCTION qsat_rho(TK, rho)

  IMPLICIT NONE
  REAL (KIND=wp)             :: qsat_rho
  REAL (KIND=wp), INTENT(IN) :: TK, rho

  !$ACC ROUTINE SEQ
  qsat_rho   = sat_pres_water(TK) / (rho * rv * TK)

END FUNCTION qsat_rho

!!!=============================================================================================

#ifndef _OPENACC
ELEMENTAL &
#endif
FUNCTION qsat_ice_rho(TK, rho)

  IMPLICIT NONE
  REAL (KIND=wp)             :: qsat_ice_rho
  REAL (KIND=wp), INTENT(IN) :: TK, rho

  !$ACC ROUTINE SEQ
  qsat_ice_rho   = sat_pres_ice(TK) / (rho * rv * TK)

END FUNCTION qsat_ice_rho

!!!=============================================================================================

#ifndef _OPENACC
ELEMENTAL &
#endif
FUNCTION dqsatdT_rho(qs, TK, rho)

  IMPLICIT NONE
  REAL (KIND=wp)            :: dqsatdT_rho
  REAL (KIND=wp), INTENT(IN):: qs, TK, rho

  !$ACC ROUTINE SEQ
  dqsatdT_rho = qs * (c5les/(TK-c4les)**2_i4 - 1.0_wp / TK)
  
END FUNCTION dqsatdT_rho

!!!=============================================================================================
#ifndef _OPENACC
ELEMENTAL &
#endif
FUNCTION dqsatdT (qs, TK)

  IMPLICIT NONE
  REAL (KIND=wp)            :: dqsatdT
  REAL (KIND=wp), INTENT(IN):: qs, TK

  !$ACC ROUTINE SEQ
  dqsatdT =     c5les * ( 1.0_wp + vtmpc1*qs ) * qs / (TK-c4les)**2

END FUNCTION dqsatdT

!!!=============================================================================================

#ifndef _OPENACC
ELEMENTAL &
#endif
FUNCTION dqsatdT_ice (qs, TK)

  IMPLICIT NONE
  REAL (KIND=wp)            :: dqsatdT_ice
  REAL (KIND=wp), INTENT(IN):: qs, TK

  !$ACC ROUTINE SEQ
  dqsatdT_ice = c5ies * ( 1.0_wp + vtmpc1*qs ) * qs / (TK-c4ies)**2

END FUNCTION dqsatdT_ice

!!!=============================================================================================

#ifndef _OPENACC
ELEMENTAL &
#endif
FUNCTION vaporization_energy(TK)

  IMPLICIT NONE
  REAL(KIND=wp)             :: vaporization_energy
  REAL(KIND=wp), INTENT(IN) :: TK

  !$ACC ROUTINE SEQ
  vaporization_energy = lvc + (cvv - clw)*TK
  
END FUNCTION vaporization_energy

!!!=============================================================================================

#ifndef _OPENACC
ELEMENTAL &
#endif
FUNCTION sublimation_energy(TK)

  IMPLICIT NONE
  REAL(KIND=wp)             :: sublimation_energy
  REAL(KIND=wp), INTENT(IN) :: TK

  !$ACC ROUTINE SEQ
  sublimation_energy = als + (cpv - ci)*(TK-tmelt) -rv*TK
  
END FUNCTION sublimation_energy

END MODULE mo_aes_thermo

