!>
!> Data module for variables of the turbulence parameterization
!------------------------------------------------------------------------------
!!
!
!> Description:
!!  This module contains variables that are used in the turbulence
!!  parameterizations. With these variables a tuning of the scheme is
!!  pssible.
!!
!! @author DWD, Ulrich.Schaettler
!!  phone:  +49  69  8062 2739
!!  fax:    +49  69  8062 3721
!!  email:  ulrich.schaettler@dwd.de
!!
!! 
!! @par Revision History  
!> ---------- ---------- ----
!! 3.21       2006/12/04 Ulrich Schaettler
!!  Initial release
!! V3_23        2007/03/30 Matthias Raschendorfer
!!  Importing 'rat_lam' from data_soil
!!  and 'clc_diag', 'q_crit', 'akt' from data_constants.
!!  Introduction of some parameters from turb_param.incf.
!!  Initialization of all parameters with default values.
!! V4_10        2009/09/11 Matthias Raschendorfer
!!  Introduction of 'a_hshr' and 'a_stab'.
!! V4_13        2010/05/11 Michael Gertz
!! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
!!
!!  Increased solar albedo of snow by T. Reinhardt / J. Helmert (2011-09-21)
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!!    violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!    copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!    an according license agreement with DWD and MPI-M.
!! </ol>
!!
!! @par Warranty
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.

MODULE mo_phyparam_soil
 
! Modules used:

#ifdef __COSMO__

USE data_parameters, ONLY :   &
    ireals,    & ! KIND-type parameter for real variables
    iintegers    ! KIND-type parameter for standard integer variables

#endif

#ifdef __ICON__

USE mo_kind, ONLY:     &
    ireals   =>wp    , &
    iintegers=>i4

#endif

!==============================================================================

IMPLICIT NONE
PUBLIC           ! All constants and variables in this module are public

INTEGER, PARAMETER :: nclass_lu = 23

!==============================================================================
! 1. Data arrays for properties of different soil types (array index)     
! -------------------------------------------------------------------

    REAL (ireals) ::          &
         c_lnd      = 2.0_ireals,     & ! surface area density of the roughness elements over land
         c_sea      = 1.5_ireals,     & ! surface area density of the waves over sea
         c_soil     = 1.0_ireals,     & ! surface area density of the (evaporative) soil surface
         e_surf     = 1.0_ireals,     & ! exponent to get the effective surface area
         red_fac                        ! reduction factor for effective area indeces multiplication

    REAL (ireals) :: zplcmxc_lu(nclass_lu)  = (/ &     !< lookup table landuse class to maximal plant cover
         &            0.80, &       ! 1 evergreen broadleaf forest   
         &            0.90, &       ! 2 deciduous broadleaf closed forest
         &            0.80, &       ! 3 deciduous broadleaf open   forest
         &            0.80, &       ! 4 evergreen needleleaf forest   
         &            0.90, &       ! 5 deciduous needleleaf forest
         &            0.90, &       ! 6 mixed leaf trees            
         &            0.80, &       ! 7 fresh water flooded trees
         &            0.80, &       ! 8 saline water flooded trees
         &            0.80, &       ! 9 mosaic tree / natural vegetation
         &            0.50, &       ! 10 burnt tree cover
         &            0.80, &       ! 11 evergreen shrubs closed-open
         &            0.80, &       ! 12 decidous shrubs closed-open
         &            0.90, &       ! 13 herbaceous vegetation closed-open
         &            0.50, &       ! 14 sparse herbaceous or grass 
         &            0.80, &       ! 15 flooded shrubs or herbaceous
         &            0.90, &       ! 16 cultivated & managed areas
         &            0.80, &       ! 17 mosaic crop / tree / natural vegetation
         &            0.90, &       ! 18 mosaic crop / shrub / grass
         &            0.05, &       ! 19 bare areas                       
         &            0.00, &       ! 20 water
         &            0.00, &       ! 21 snow & ice 
         &            0.20, &       ! 22 artificial surface   
         &            0.       /)   ! 23 undefined

    REAL (ireals) :: zlaimxc_lu(nclass_lu)  = (/ &      !< lookup table landuse class to maximal leaf area index
         &            5.00, &       ! evergreen broadleaf forest   
         &            6.00, &       ! deciduous broadleaf closed forest
         &            4.00, &       ! deciduous broadleaf open   forest
         &            5.00, &       ! evergreen needleleaf forest   
         &            5.00, &       ! deciduous needleleaf forest
         &            5.00, &       ! mixed leaf trees            
         &            5.00, &       ! fresh water flooded trees
         &            5.00, &       ! saline water flooded trees
         &            2.50, &       ! mosaic tree / natural vegetation
         &            0.60, &       ! burnt tree cover
         &            3.00, &       ! evergreen shrubs closed-open
         &            1.50, &       ! decidous shrubs closed-open
         &            3.10, &       ! herbaceous vegetation closed-open
         &            0.60, &       ! sparse herbaceous or grass 
         &            2.00, &       ! flooded shrubs or herbaceous
         &            3.30, &       ! cultivated & managed areas
         &            3.00, &       ! mosaic crop / tree / natural vegetation
         &            3.50, &       ! mosaic crop / shrub / grass
         &            0.60, &       ! bare areas                       
         &            0.00, &       ! water
         &            0.00, &       ! snow & ice 
         &            1.00, &       ! artificial surface   
         &            0.      /)    ! undefined


    REAL (ireals) :: zrd_lu(nclass_lu)   = (/ &         !< lookup table landuse class to root depth [m]
         &          1.00,       &       ! evergreen broadleaf forest   
         &          1.00,       &       ! deciduous broadleaf closed forest
         &          2.00,       &       ! deciduous broadleaf open   forest
         &          0.60,       &       ! evergreen needleleaf forest   
         &          0.60,       &       ! deciduous needleleaf forest
         &          0.80,       &       ! mixed leaf trees            
         &          1.00,       &       ! fresh water flooded trees
         &          1.00,       &       ! saline water flooded trees
         &          1.00,       &       ! mosaic tree / natural vegetation
         &          0.30,       &       ! burnt tree cover
         &          1.00,       &       ! evergreen shrubs closed-open
         &          2.00,       &       ! decidous shrubs closed-open
         &          0.60,       &       ! herbaceous vegetation closed-open
         &          0.30,       &       ! sparse herbaceous or grass 
         &          0.40,       &       ! flooded shrubs or herbaceous
         &          1.00,       &       ! cultivated & managed areas
         &          1.00,       &       ! mosaic crop / tree / natural vegetation
         &          1.00,       &       ! mosaic crop / shrub / grass
         &          0.30,       &       ! bare areas                       
         &          0.00,       &       ! water
         &          0.00,       &       ! snow & ice 
         &          0.60,       &       ! artificial surface   
         &          0.         /)       ! undefined

    REAL (ireals) :: zrs_min_lu(nclass_lu) =(/ &
         &          250.0,      &       ! evergreen broadleaf forest   
         &          150.0,      &       ! deciduous broadleaf closed forest
         &          150.0,      &       ! deciduous broadleaf open   forest
         &          150.0,      &       ! evergreen needleleaf forest   
         &          150.0,      &       ! deciduous needleleaf forest
         &          150.0,      &       ! mixed leaf trees            
         &          150.0,      &       ! fresh water flooded trees
         &          150.0,      &       ! saline water flooded trees
         &          150.0,      &       ! mosaic tree / natural vegetation
         &          150.0,      &       ! burnt tree cover
         &          120.0,      &       ! evergreen shrubs closed-open
         &          120.0,      &       ! decidous shrubs closed-open
         &          40.0,      &       ! herbaceous vegetation closed-open
         &          40.0,      &       ! sparse herbaceous or grass 
         &          40.0,      &       ! flooded shrubs or herbaceous
         &          120.0,      &       ! cultivated & managed areas
         &          120.0,      &       ! mosaic crop / tree / natural vegetation
         &          100.0,      &       ! mosaic crop / shrub / grass
         &          120.0,      &       ! bare areas                       
         &          120.0,      &       ! water
         &          120.0,      &       ! snow & ice 
         &          120.0,     &       ! artificial surface   
         &          0.         /)       ! undefined
 
  REAL  (KIND=ireals) ::  &
!   a) parameters describing the soil water budget
    cporv (10), &  !  pore volume (fraction of volume)
    cfcap (10), &  !  field capacity (fraction of volume)
    cpwp  (10), &  !  plant wilting point (fraction of volume)
    cadp  (10), &  !  air dryness point (fraction of volume)
    cik2  (10), &  !  minimum infiltration rate (kg/s*m**2)
    ckw0  (10), &  !  parameter for determination of hydr. conductivity (m/s)
    ckw1  (10), &  !  parameter for determination of hydr. conductivity (1)
    cdw0  (10), &  !  parameter for determination of hydr. diffusivity (m**2/s)
    cdw1  (10), &  !  parameter for determination of hydr. diffusivity (1)
    crock (10), &  !  rock/ice/water indicator (hydrological calculations 
                   !  only for crock=1)
!   b) parameters describing the soil heat budget
    cdz1  (10), &  !  top layer thickness (EFR-method)
    crhoc (10), &  !  soil heat capacity  (J/K*m**3)
    cala0 (10), &  !  parameters for the determination of
    cala1 (10), &  !      the soil heat conductivity (W/(K*m))
    csalb (10), &  !  solar albedo for dry soil                            
    csalbw(10), &  !  slope of solar albedo with respect to soil water content     
!   c) additional parameters for the BATS scheme (Dickinson)
    ck0di (10), &  !  (m/s)
    cbedi (10), &  !  (1)
    clgk0 (10), &  !  auxiliary variable
!   d) additional parameters for soil water content dependent freezing/melting
    csandf(10), &  !  mean fraction of sand (weight percent)
    cclayf(10)     !  mean fraction of clay (weight percent)
 

  ! Initialization of soil type parameters except cdz1 
  ! (being calculated during execution)

  ! soil type  / ice           , rock          , sand          , sandyloam     , loam          , &
  !            & clayloam      , clay          , peat          , sea water     , sea ice         /

  DATA  cporv  / 1.E-10_ireals , 1.E-10_ireals , 0.364_ireals  , 0.445_ireals  , 0.455_ireals  , &
               & 0.475_ireals  , 0.507_ireals  , 0.863_ireals  , 1.E-10_ireals , 1.E-10_ireals   /

  DATA  cfcap  / 1.E-10_ireals , 1.E-10_ireals , 0.196_ireals  , 0.260_ireals  , 0.340_ireals  , &
               & 0.370_ireals  , 0.463_ireals  , 0.763_ireals  , 1.E-10_ireals , 1.E-10_ireals   /

  DATA  cpwp   / 0.000_ireals  , 0.000_ireals  , 0.042_ireals  , 0.100_ireals  , 0.110_ireals  , &
               & 0.185_ireals  , 0.257_ireals  , 0.265_ireals  , 0.000_ireals  , 0.000_ireals    /

  DATA  cadp   / 0.000_ireals  , 0.000_ireals  , 0.012_ireals  , 0.030_ireals  , 0.035_ireals  , &
               & 0.060_ireals  , 0.065_ireals  , 0.098_ireals  , 0.000_ireals  , 0.000_ireals    /

  DATA  crhoc  / 1.92E6_ireals , 2.10E6_ireals , 1.28E6_ireals , 1.35E6_ireals , 1.42E6_ireals , &
               & 1.50E6_ireals , 1.63E6_ireals , 0.58E6_ireals , 4.18E6_ireals , 1.92E6_ireals   /

  DATA  cik2   / 0.0000_ireals , 0.0000_ireals , 0.0035_ireals , 0.0023_ireals , 0.0010_ireals , &
               & 0.0006_ireals , 0.0001_ireals , 0.0002_ireals , 0.0000_ireals , 0.0000_ireals   /

  DATA  ckw0   / 0.00_ireals   , 0.00_ireals   , 479.E-7_ireals, 943.E-8_ireals, 531.E-8_ireals, &
               & 764.E-9_ireals,  17.E-9_ireals, 58.E-9_ireals , 0.00_ireals   , 0.00_ireals     /

  DATA  ckw1   / 0.00_ireals   , 0.00_ireals   , -19.27_ireals , -20.86_ireals , -19.66_ireals , &
               & -18.52_ireals , -16.32_ireals , -16.48_ireals , 0.00_ireals   , 0.00_ireals     /

  DATA  cdw0   / 0.00_ireals   , 0.00_ireals   , 184.E-7_ireals, 346.E-8_ireals, 357.E-8_ireals, &
               & 118.E-8_ireals, 442.E-9_ireals, 106.E-9_ireals, 0.00_ireals   , 0.00_ireals     /

  DATA  cdw1   / 0.00_ireals   , 0.00_ireals   , -8.45_ireals  , -9.47_ireals  , -7.44_ireals  , &
               & -7.76_ireals  , -6.74_ireals  , -5.97_ireals  , 0.00_ireals   , 0.00_ireals     /

  DATA  crock  / 0.00_ireals   , 0.00_ireals   , 1.00_ireals   , 1.00_ireals   , 1.00_ireals   , &
               & 1.00_ireals   , 1.00_ireals   , 1.00_ireals   , 0.00_ireals   , 0.00_ireals     /

  DATA  cala0  / 2.26_ireals   , 2.41_ireals   , 0.30_ireals   , 0.28_ireals   , 0.25_ireals   , &
               & 0.21_ireals   , 0.18_ireals   , 0.06_ireals   , 1.0_ireals    , 2.26_ireals     /

  DATA  cala1  / 2.26_ireals   , 2.41_ireals   , 2.40_ireals   , 2.40_ireals   , 1.58_ireals   , &
               & 1.55_ireals   , 1.50_ireals   , 0.50_ireals   , 1.0_ireals    , 2.26_ireals     /

  DATA  csalb  / 0.70_ireals   , 0.30_ireals   , 0.30_ireals   , 0.25_ireals   , 0.25_ireals   , &
               & 0.25_ireals   , 0.25_ireals   , 0.20_ireals   , 0.07_ireals   , 0.70_ireals     /

  DATA  csalbw / 0.00_ireals   , 0.00_ireals   , 0.44_ireals   , 0.27_ireals   , 0.24_ireals   , &
               & 0.23_ireals   , 0.22_ireals   , 0.10_ireals   , 0.00_ireals   , 0.00_ireals     /

  DATA  ck0di  / 1.E-4_ireals  , 1.E-4_ireals  , 2.E-4_ireals  , 2.E-5_ireals  , 6.E-6_ireals  , &
               & 2.E-6_ireals  , 1.E-6_ireals  , 1.5E-6_ireals , 0.00_ireals   , 0.00_ireals     /

  DATA  cbedi  / 1.0_ireals    , 1.00_ireals   , 3.5_ireals    , 4.8_ireals    , 6.1_ireals    , &
               & 8.6_ireals    , 10.00_ireals  , 9.0_ireals    , 0.00_ireals   , 0.00_ireals     /

  DATA  csandf / 0.00_ireals   , 0.00_ireals   , 90.00_ireals  , 65.0_ireals   , 40.0_ireals   , &
               & 35.0_ireals   , 15.0_ireals   , 90.00_ireals  , 0.00_ireals   , 0.00_ireals     /

  DATA  cclayf / 0.00_ireals   , 0.00_ireals   , 5.0_ireals    , 10.0_ireals   , 20.0_ireals   , &
               & 35.0_ireals   , 70.0_ireals   , 5.0_ireals    , 0.00_ireals   , 0.00_ireals     /
 

!==============================================================================

! 2. Additional parameters for the soil model                             
! -------------------------------------------------------------------

  REAL  (KIND=ireals) ::  &
!==============================================================================

    csalb_p    = 0.15_ireals  , & !  solar albedo of ground covered by plants
    csalb_snow = 0.70_ireals  , & !  solar albedo of ground covered by snow
    !csalb_snow_min = 0.400_ireals, &
    !                      ! min. solar albedo of snow for forest free surfaces
    !csalb_snow_max = 0.700_ireals, &
    !                       ! max. solar albedo of snow for forest free surfaces
    ! T.R. 2011-09-21 csalb_snow_min/max set to values used in GME
    csalb_snow_min = 0.500_ireals, &
                           ! min. solar albedo of snow for forest free surfaces
    csalb_snow_max = 0.850_ireals, &
                           ! max. solar albedo of snow for forest free surfaces
  ! for possible later use:
    !csalb_snow_fe  = 0.200_ireals , &  ! solar albedo of snow for surfaces with evergreen forest
    !csalb_snow_fd  = 0.200_ireals , &  ! solar albedo of snow for surfaces with deciduous forest
    ! T.R. 2011-09-21 snow albedos for forests set to values used in GME
    csalb_snow_fe  = 0.270_ireals , &  ! solar albedo of snow for surfaces with evergreen forest
    csalb_snow_fd  = 0.320_ireals , &  ! solar albedo of snow for surfaces with deciduous forest
    ctalb      = 0.004_ireals , & !  thermal albedo ( of all soil types )   
    cf_snow    = 0.0150_ireals, & !  parameter for the calculation of the 
                                  !  fractional snow coverage
  ! for the multi-layer soil model
    cwhc       = 0.04_ireals,   & !  water holding capacity of snow ()
    chcond     = 0.01_ireals,   & !  saturation hydraulic conductivity of snow ()
    ca2        = 6.6E-07_ireals,& !  activation energy (for snow metamorphosis) (J)
    csigma     = 75._ireals,    & !  snow metamorphosis, Pa

  ! cf_w changed from 0.0004 to 0.0010 (in agreement with GME)
    cf_w       = 0.0010_ireals, & !  parameter for the calculation of the
                                  !  fractional water coverage

    csvoro     = 1.0000_ireals, & !  parameter to estimate the subgrid-scale 
                                  !  variation of orography
    cik1       = 0.0020_ireals, & !  parameter for the determination of the 
                                  !  maximum infiltaration
    cwimax     = 0.0005_ireals, & !  parameter for the determination of the 
    cwimax_ml  = 1.E-6_ireals,  & !  maximum interception water content
    ctau_i     = 1000.0_ireals, & !  time constatant for the drainage from the 
                                  !  interception storeage 
    cakw       = 0.8000_ireals, & !  parameter for averaging the water contents
                                  !  of the top and middle soil water layers to 
                                  !  calculate the hydraulic diffusivity and 
                                  !  conductiviy

    ctau1      = 1.0000_ireals, & !  first adjustment time period in EFR-method
    ctau2      = 5.0000_ireals, & !  second adjustment time period in EFR-method
    chc_i      = 2100.0_ireals, & !  heat capacity of ice     
    chc_w      = 4180.0_ireals, & !  heat capacity of water     

    cdzw12     = 0.1000_ireals, & !  thickness of upper soil water layer in 
                                  !  two-layer model         
    cdzw22     = 0.9000_ireals, & !  thickness of lower soil water layer in 
                                  !  two-layer model      
    cdzw13     = 0.0200_ireals, & !  thickness of upper soil water layer in 
                                  !  three-layer model
    cdzw23     = 0.0800_ireals, & !  thickness of middle soil water layer in 
                                  !  three-layer model 
    cdzw33     = 0.9000_ireals    !  thickness of lower soil water layer in 
                                  !  three-layer model

  REAL  (KIND=ireals) ::  &
    cdsmin     = 0.0100_ireals, & !  minimum snow depth
    crhosmin   = 500.00_ireals, & !  minimum density of snow
    crhosmax   = 800.00_ireals, & !  maximum density of snow
    crhosmin_ml=  50.00_ireals, & !  minimum density of snow
    crhosmax_ml= 400.00_ireals, & !  maximum density of snow
    crhosminf  =  50.00_ireals, & !  minimum density of fresh snow
    crhosmaxf  = 150.00_ireals, & !  maximum density of fresh snow
    crhosmint  =   0.20_ireals, & !  minimum value of time constant for ageing 
                                  !  of snow
    crhosmaxt  =   0.40_ireals, & !  maximum value of time constant for ageing 
                                  !  of snow
    csnow_tmin = 258.15_ireals, & !  lower threshold temperature of snow for 
                                  !  ageing and fresh snow density computation 
                                  !  ( = 273.15-15.0)
    crhos_dw   = 300.00_ireals, & !  change of snow density with water content
    calasmin   = 0.2000_ireals, & !  minimum heat conductivity of snow (W/m K)
    calasmax   = 1.5000_ireals, & !  maximum heat conductivity of snow (W/m K)
    calas_dw   = 1.3000_ireals, & !  change of snow heat conductivity with
                                  !  water content                (W/(m**2) K)
   
    crhowm     =    0.8_ireals    , & !  BATS (1)
    cdmin      =    0.25E-9_ireals, & !  BATS (m**2/s)
    cfinull    =    0.2_ireals    , & !  BATS (m)
    ckrdi      =    1.0E-5_ireals , & !  BATS (m/s)
    cdash      =    0.05_ireals   , & !  BATS ((m/s)**1/2)
    clai       =    3.0_ireals    , & !  BATS
    cparcrit   =  100.0_ireals    , & !  BATS (W/m**2)
    ctend      =  313.15_ireals   , & !  BATS (K)
    csatdef    = 4000.0_ireals    , & !  BATS (Pa)

    !Minimum and maximum value of tomatal resistance (s/m)
    !used by the Pen.-Mont. method for vegetation transpiration
    !(itype_trvg=2):
    crsmin     = 150.0_ireals     , & !  BATS (s/m)
    crsmax     = 4000.0_ireals    , & !  BATS (s/m)

! crsmax increased from 1000 to 4000 s/m (to reduce latent heat flux).

!   Determine constants clgk0 for BATS-scheme ! from TERRA_ML

!>JH to transfered into TERRA
!!$    DO jb       = 1, 10
!!$      clgk0(jb) = LOG10(MAX(zepsi,ck0di(jb)/ckrdi))
!!$    END DO

    b3=273.16_ireals
!==============================================================================

END MODULE mo_phyparam_soil
