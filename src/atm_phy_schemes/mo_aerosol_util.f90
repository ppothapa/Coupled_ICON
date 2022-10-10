!>
!! <Short description of module for listings and indices>
!!
!! <Describe the concepts of the procedures and algorithms used in the module.>
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author Thorsten Reinhardt, AGeoBw, Offenbach
!!
!!
!! @par Revision History
!! Initial release by Thorsten Reinhardt, AGeoBw, Offenbach (2011-09-06)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_aerosol_util

  USE mo_impl_constants,         ONLY: min_rlcell_int, iss, iorg, ibc, iso4, idu, nclass_aero
  USE mo_impl_constants_grf,     ONLY: grf_bdywidth_c
  USE mo_math_constants,         ONLY: rad2deg
  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: finish
  USE mo_loopindices,            ONLY: get_indices_c
  USE mo_lrtm_par,               ONLY: jpband => nbndlw
  USE mo_model_domain,           ONLY: t_patch
  USE mo_intp_data_strc,         ONLY: t_int_state
  USE mo_srtm_config,            ONLY: jpsw
  USE mo_lnd_nwp_config,         ONLY: ntiles_lnd, dzsoil, isub_water
  USE mo_nwp_tuning_config,      ONLY: tune_dust_abs
  USE mo_aerosol_sources_types,  ONLY: p_dust_source_const
  USE mo_aerosol_sources,        ONLY: aerosol_dust_aod_source, aerosol_ssa_aod_source
  USE mo_math_laplace,           ONLY: nabla2_scalar

  IMPLICIT NONE

  PRIVATE

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_aerosol_util'

  !RRTM
  REAL  (wp)              ::           &
  zaea_rrtm(jpsw+jpband,5), &  ! ratio of optical thickness for the absorption in spectral
                               ! interval jpspec  and total optical thickness at 0.55m*1.E-06 
                               ! for an aerosoltyp specified by second array index
  zaes_rrtm(jpsw+jpband,5), &  ! analog for the optical thickness of scattering 
  zaeg_rrtm(jpsw+jpband,5)!, zaef_rrtm(jpsw+jpband,5)

  !ecRad
  TYPE t_tegen_scal_factors
    ! Total number of wavelength bands
    INTEGER :: n_bands
    ! Scaling factors from 550nm to wavelengths bands
    REAL(wp), ALLOCATABLE :: &
      &  absorption(:,:),            & !< Dim [n_bands, nspecies=5]
      &  scattering(:,:),            & !< Dim [n_bands, nspecies=5] 
      &  asymmetry(:,:)                !< Dim [n_bands, nspecies=5]
    CONTAINS
      PROCEDURE :: init     => init_tegen_scal_factors
      PROCEDURE :: finalize => finalize_tegen_scal_factors
  END TYPE t_tegen_scal_factors

  TYPE(t_tegen_scal_factors), TARGET :: tegen_scal_factors

  PUBLIC :: zaea_rrtm, zaes_rrtm, zaeg_rrtm
  PUBLIC :: aerdis
  PUBLIC :: init_aerosol_props_tegen_rrtm, tune_dust
  PUBLIC :: prog_aerosol_2D, aerosol_2D_diffusion
  PUBLIC :: tegen_scal_factors, init_aerosol_props_tegen_ecrad

  !$ACC DECLARE CREATE(zaea_rrtm, zaes_rrtm, zaeg_rrtm)

CONTAINS

  !!  Subroutine aerdis is simplified version from COSMO model (version 4.16).
  !!
  !! @par Revision History
  !! Initial Release by Thorsten Reinhardt, AGeoBw, Offenbach (2011-02-28)
  !! Transferred to mo_aerosol_util Sophia Schaefer, DWD (2021-06-21)
   
  SUBROUTINE aerdis ( klevp1, kbdim, jcs, jce, petah,  pvdaes, pvdael, pvdaeu, pvdaed, lacc )
    
    !------------------------------------------------------------------------------
    !
    ! Description:
    !
    ! The module procedure aerdis provides parameters for the vertical distribution
    ! of aerosols (based on the original code of J.F. Geleyn (ECMWF, 4.11.82).
    !
    ! The routine computes the values PVDAE* (* = s, l, u or d for sea, land
    ! urban or desert) of a surfach-normalised vertical distribution of aerosols'
    ! optical depth from the argument petah (vertical coordinate) at klevp1 levels.
    ! It also sets values for non-geograpically weighted total optical depths (at
    ! 55 micrometer wavelength) paeopn for the same four types and similar optical
    ! depths diveded by pressure for bachground well-mixed aerosols of three types
    ! p**bga (** = tr, vo or st for tropospheric, volcanic (stratosperic ashes) or
    ! stratosperic (sulfuric type)). It finally sets values for the power to be
    ! applied to a temperature ratio smaller than two in order to obtain an index
    ! one in the stratosphere and zero in the troposphere with a relatively smooth
    ! transistion (ptrpt), as well as for adsorption coefficients fo water to the
    ! three type of troposperic aerosols (paeadk) with a minimum value ( in the 
    ! whole atmosphere) for the sum of the products paeadk by the optical depths
    ! divided by pressure thickness: paeadm. 
    !
    ! Method:
    !
    ! Straightforward, equivalent heights are given in meters (8434 for the
    ! atmosphere) and tropospheric and stratospheric pressure boundary values
    ! are set at 101325 and 19330 Pascal. 
    !
    !------------------------------------------------------------------------------
    
    ! Subroutine arguments:
    ! --------------------
    
    ! Input data
    ! ----------
    INTEGER, INTENT (IN) ::  &
      & klevp1,         &           ! number of model layer interfaces
      & kbdim,          &
      & jcs,            &
      & jce

    REAL    (wp), INTENT (IN) ::  &
      petah(kbdim,klevp1)    ! normalized vertical coordinate at half levels

    LOGICAL, OPTIONAL, INTENT(IN) :: lacc

    ! Output data
    ! -----------
    REAL    (wp), INTENT (OUT) ::  &
      pvdaes(kbdim,klevp1), & ! normalized vertical distribution (sea)
      pvdael(kbdim,klevp1), & ! normalized vertical distribution (land)
      pvdaeu(kbdim,klevp1), & ! normalized vertical distribution (urban)
      pvdaed(kbdim,klevp1)    ! normalized vertical distrubution (desert)

    ! Local parameters:
    ! -------------
    REAL (wp), PARAMETER  ::  &
      zhss = 8434.0_wp/1000.0_wp ,  & !
      zhsl = 8434.0_wp/1000.0_wp ,  & !
      zhsu = 8434.0_wp/1000.0_wp ,  & !
      zhsd = 8434.0_wp/3000.0_wp      !

    INTEGER :: jc,jk
    REAL(wp) :: log_eta
    LOGICAL :: lzacc

    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    ! Begin Subroutine aerdis              
    !------------------------------------------------------------------------------

    IF(PRESENT(lacc)) THEN
        lzacc = lacc
    ELSE
        lzacc = .FALSE.
    ENDIF
    !$acc data present(petah,pvdaes,pvdael,pvdaeu,pvdaed) if (lzacc)

    ! default data present
    !$acc parallel default(none) ASYNC(1) if (lzacc)
    !$acc loop gang vector
    DO jc=jcs,jce
      pvdaes(jc,1) = 0.0_wp
      pvdael(jc,1) = 0.0_wp
      pvdaeu(jc,1) = 0.0_wp
      pvdaed(jc,1) = 0.0_wp
    ENDDO
    !$acc end parallel

!!$  IF(petah(1).NE.0._wp) THEN
!!$     pvdaes(1) = petah(1)**zhss
!!$     pvdael(1) = petah(1)**zhsl
!!$     pvdaeu(1) = petah(1)**zhsu
!!$     pvdaed(1) = petah(1)**zhsd
!!$  END IF

    !$acc parallel default(none) ASYNC(1) if (lzacc)
    !$acc loop gang vector collapse(2) PRIVATE(log_eta)
    DO jk=2,klevp1
      DO jc=jcs,jce
        log_eta       = LOG(petah(jc,jk))
        pvdaes(jc,jk) = EXP(zhss*log_eta) ! petah(jc,jk)**zhss
        pvdael(jc,jk) = pvdaes(jc,jk)     ! petah(jc,jk)**zhsl; zhsl is the same as zhss
        pvdaeu(jc,jk) = pvdaes(jc,jk)     ! petah(jc,jk)**zhsu; zhsu is the same as zhss
        pvdaed(jc,jk) = EXP(zhsd*log_eta) ! petah(jc,jk)**zhsd
      ENDDO
    ENDDO
    !$acc end parallel

    !$acc end data

  END SUBROUTINE aerdis

  
  SUBROUTINE init_aerosol_props_tegen_rrtm

  ! the following aerosol types (second array index) are considered:
  ! 1. continental, 2. maritime, 3. desert, 4. urban, 5. stratospheric background (SB)

   !absorption
   zaea_rrtm=RESHAPE( (/ &
     &0.0304_wp,0.0367_wp,0.0462_wp,0.0566_wp,0.0496_wp,0.0336_wp,0.0355_wp,0.0456_wp,&
     &0.0272_wp,0.0264_wp,0.0290_wp,0.0156_wp,0.0165_wp,0.0157_wp,0.0138_wp,0.0401_wp,&
     &0.0401_wp,0.0760_wp,0.0214_wp,0.0227_wp,0.0295_wp,0.0394_wp,0.0431_wp,0.0519_wp,&
     &0.0611_wp,0.0774_wp,0.1012_wp,0.1412_wp,0.2632_wp,0.0324_wp,                    &
     &0.1096_wp,0.1614_wp,0.2294_wp,0.2506_wp,0.2242_wp,0.1190_wp,0.0680_wp,0.0664_wp,&
     &0.0656_wp,0.0749_wp,0.1250_wp,0.0425_wp,0.0498_wp,0.0425_wp,0.0259_wp,0.1619_wp,&
     &0.1619_wp,0.2152_wp,0.0139_wp,0.0119_wp,0.0046_wp,0.0036_wp,0.0020_wp,0.0016_wp,&
     &0.0012_wp,0.0013_wp,0.0016_wp,0.0035_wp,0.0147_wp,0.0882_wp,                    &
     &0.0974_wp,0.1529_wp,0.1643_wp,0.1373_wp,0.1753_wp,0.1923_wp,0.2804_wp,0.2426_wp,&
     &0.1263_wp,0.1321_wp,0.0979_wp,0.0664_wp,0.0360_wp,0.0311_wp,0.0325_wp,0.0833_wp,&
     &0.0833_wp,0.1170_wp,0.0739_wp,0.0631_wp,0.0604_wp,0.0628_wp,0.0645_wp,0.0677_wp,&
     &0.0843_wp,0.1328_wp,0.2224_wp,0.3022_wp,0.3579_wp,0.1820_wp,                    &
     &0.0267_wp,0.0329_wp,0.0420_wp,0.0515_wp,0.0461_wp,0.0332_wp,0.0354_wp,0.0447_wp,&
     &0.0303_wp,0.0306_wp,0.0342_wp,0.0248_wp,0.0274_wp,0.0276_wp,0.0271_wp,0.0526_wp,&
     &0.0526_wp,0.0903_wp,0.0450_wp,0.0492_wp,0.0596_wp,0.0754_wp,0.0842_wp,0.1082_wp,&
     &0.1429_wp,0.1926_wp,0.2595_wp,0.3379_wp,0.4761_wp,0.0340_wp,                    &
     &0.0060_wp,0.0117_wp,0.0269_wp,0.0222_wp,0.0195_wp,0.0398_wp,0.0733_wp,0.1091_wp,&     ! SB
     &0.1124_wp,0.0415_wp,0.0424_wp,0.0495_wp,0.0451_wp,0.0484_wp,0.0540_wp,0.0735_wp,&     ! SB
     &0.0735_wp,0.0188_wp,0.0021_wp,0.0014_wp,0.0007_wp,0.0002_wp,0.0000_wp,0.0000_wp,&     ! SB
     &0.0000_wp,0.0000_wp,0.0000_wp,0.0000_wp,0.0000_wp,0.0628_wp/),(/jpsw+jpband,5/))      ! SB

   !scattering
   zaes_rrtm=RESHAPE( (/ &
     &0.0060_wp,0.0107_wp,0.0134_wp,0.0150_wp,0.0152_wp,0.0200_wp,0.0232_wp,0.0211_wp,&
     &0.0112_wp,0.0186_wp,0.0128_wp,0.0260_wp,0.0339_wp,0.0368_wp,0.0409_wp,0.0527_wp,&
     &0.0527_wp,0.0621_wp,0.0715_wp,0.0929_wp,0.1276_wp,0.1895_wp,0.2350_wp,0.3930_wp,&
     &0.6641_wp,0.9834_wp,1.3737_wp,1.7160_wp,1.9115_wp,0.0198_wp,                    &
     &0.0188_wp,0.0421_wp,0.0576_wp,0.0547_wp,0.0430_wp,0.0367_wp,0.0806_wp,0.1209_wp,&
     &0.1681_wp,0.2257_wp,0.2440_wp,0.3622_wp,0.4540_wp,0.5026_wp,0.5765_wp,0.5986_wp,&
     &0.5986_wp,0.5225_wp,0.7420_wp,0.8311_wp,0.8970_wp,0.9444_wp,0.9637_wp,0.9763_wp,&
     &0.9855_wp,1.0034_wp,1.0337_wp,1.0640_wp,1.0795_wp,0.1312_wp,                    &
     &0.0458_wp,0.0823_wp,0.0667_wp,0.0642_wp,0.1080_wp,0.1471_wp,0.2422_wp,0.1216_wp,&
     &0.0717_wp,0.1616_wp,0.2027_wp,0.3042_wp,0.4045_wp,0.4369_wp,0.4685_wp,0.5043_wp,&
     &0.5043_wp,0.5782_wp,0.6898_wp,0.7477_wp,0.7926_wp,0.8320_wp,0.8503_wp,0.8736_wp,&
     &0.8874_wp,0.8737_wp,0.8278_wp,0.7857_wp,0.7571_wp,0.1714_wp,                    &
     &0.0048_wp,0.0085_wp,0.0107_wp,0.0119_wp,0.0121_wp,0.0160_wp,0.0185_wp,0.0170_wp,&
     &0.0090_wp,0.0150_wp,0.0103_wp,0.0210_wp,0.0274_wp,0.0298_wp,0.0332_wp,0.0430_wp,&
     &0.0430_wp,0.0485_wp,0.0593_wp,0.0776_wp,0.1073_wp,0.1610_wp,0.2008_wp,0.3398_wp,&
     &0.5809_wp,0.8701_wp,1.2309_wp,1.5535_wp,1.7368_wp,0.0159_wp,                    &
     &0.0000_wp,0.0000_wp,0.0000_wp,0.0000_wp,0.0001_wp,0.0003_wp,0.0006_wp,0.0008_wp,&     ! SB
     &0.0005_wp,0.0003_wp,0.0008_wp,0.0013_wp,0.0024_wp,0.0030_wp,0.0040_wp,0.0059_wp,&     ! SB
     &0.0059_wp,0.0123_wp,0.0236_wp,0.0384_wp,0.0651_wp,0.1246_wp,0.1801_wp,0.3807_wp,&     ! SB
     &0.7105_wp,1.0514_wp,1.3754_wp,1.5334_wp,1.5495_wp,0.0009_wp/),(/jpsw+jpband,5/))      ! SB

   !asymmetry factor
   zaeg_rrtm=RESHAPE( (/ &
     &0.4388_wp,0.5396_wp,0.6191_wp,0.6535_wp,0.6876_wp,0.6718_wp,0.6493_wp,0.6782_wp,&
     &0.7958_wp,0.7537_wp,0.7757_wp,0.7821_wp,0.7583_wp,0.7487_wp,0.7351_wp,0.6917_wp,&
     &0.6917_wp,0.6989_wp,0.6982_wp,0.6726_wp,0.6426_wp,0.6294_wp,0.6337_wp,0.6582_wp,&
     &0.6850_wp,0.7061_wp,0.7212_wp,0.7306_wp,0.7417_wp,0.6978_wp,                    &
     &0.4062_wp,0.4507_wp,0.4878_wp,0.5302_wp,0.5850_wp,0.6962_wp,0.7242_wp,0.7293_wp,&
     &0.7414_wp,0.7484_wp,0.7607_wp,0.7785_wp,0.7805_wp,0.7785_wp,0.7724_wp,0.7690_wp,&
     &0.7690_wp,0.8348_wp,0.8316_wp,0.8170_wp,0.8074_wp,0.7990_wp,0.7954_wp,0.7897_wp,&
     &0.7884_wp,0.7927_wp,0.8001_wp,0.8057_wp,0.8076_wp,0.7462_wp,                    &
     &0.4219_wp,0.3928_wp,0.5306_wp,0.6229_wp,0.5544_wp,0.5454_wp,0.4353_wp,0.5736_wp,&
     &0.7502_wp,0.6957_wp,0.7038_wp,0.6881_wp,0.6740_wp,0.6739_wp,0.6784_wp,0.6969_wp,&
     &0.6969_wp,0.7068_wp,0.6965_wp,0.6918_wp,0.6904_wp,0.6911_wp,0.6915_wp,0.6952_wp,&
     &0.7080_wp,0.7326_wp,0.7689_wp,0.8000_wp,0.8206_wp,0.5788_wp,                    &
     &0.4387_wp,0.5394_wp,0.6187_wp,0.6531_wp,0.6871_wp,0.6712_wp,0.6482_wp,0.6756_wp,&
     &0.7930_wp,0.7498_wp,0.7685_wp,0.7766_wp,0.7520_wp,0.7419_wp,0.7277_wp,0.6828_wp,&
     &0.6828_wp,0.6875_wp,0.6872_wp,0.6622_wp,0.6333_wp,0.6209_wp,0.6250_wp,0.6479_wp,&
     &0.6725_wp,0.6912_wp,0.7043_wp,0.7129_wp,0.7254_wp,0.6956_wp,                    &
     &0.0021_wp,0.0039_wp,0.0061_wp,0.0078_wp,0.0109_wp,0.0161_wp,0.0201_wp,0.0206_wp,&     ! SB
     &0.0217_wp,0.0320_wp,0.0428_wp,0.0583_wp,0.0773_wp,0.0856_wp,0.0985_wp,0.1310_wp,&     ! SB
     &0.1310_wp,0.1906_wp,0.2625_wp,0.3154_wp,0.3869_wp,0.4787_wp,0.5279_wp,0.6272_wp,&     ! SB
     &0.6941_wp,0.7286_wp,0.7358_wp,0.7177_wp,0.6955_wp,0.0616_wp/),(/jpsw+jpband,5/))      ! SB

    !$ACC UPDATE DEVICE(zaea_rrtm, zaes_rrtm, zaeg_rrtm)

  END SUBROUTINE init_aerosol_props_tegen_rrtm

  !---------------------------------------------------------------------------------------
  !>
  !! SUBROUTINE init_aerosol_props_tegen_ecrad
  !! Initializes scaling factors used to scale 550nm AOD to different wavelength bands
  !! of ecRad.
  !! In case the RRTM gas model is chosen, the previously used lookup tables are copied
  !! into the new data structure.
  !!
  !! @par Revision History
  !! Initial release by Daniel Rieger, Deutscher Wetterdienst, Offenbach (2021-10-20)
  !!
  !---------------------------------------------------------------------------------------
  SUBROUTINE init_aerosol_props_tegen_ecrad(n_bands_sw,n_bands_lw, l_rrtm_gas_model)
    INTEGER, INTENT(in) ::      &
      &  n_bands_sw, n_bands_lw   !< Number of SW/LW bands
    LOGICAL, INTENT(in) ::      &
      &  l_rrtm_gas_model         !< Use RRTM gas model (mimic legacy behavior)
    ! Local variables
    CHARACTER(len=*), PARAMETER :: routine = modname//'::init_aerosol_props_tegen_ecrad'

    CALL tegen_scal_factors%init(n_bands_sw+n_bands_lw)

    IF (l_rrtm_gas_model) THEN ! Use look up tables from legacy routine
      ! Sanity check
      IF ((tegen_scal_factors%n_bands) /= 30) &
        &  CALL finish(routine,'ecRad wavelength bands / gas model mismatch.')
      CALL init_aerosol_props_tegen_rrtm()
      tegen_scal_factors%absorption(:,:) = zaea_rrtm(:,:)
      tegen_scal_factors%scattering(:,:) = zaes_rrtm(:,:)
      tegen_scal_factors%asymmetry (:,:) = zaeg_rrtm(:,:)
      !$ACC UPDATE DEVICE(tegen_scal_factors%absorption, &
      !$ACC   tegen_scal_factors%scattering, tegen_scal_factors%asymmetry)
    ELSE
      ! This part will be used for ecckd in the future.
      ! Here, the number of bands is flexible and the
      ! scaling coefficients must be read from a file
      ! (or new lookup tables for heavily-used configurations)
      CALL finish(routine,'Currently only RRTM Gas Optics implemented.')
    ENDIF
    
  END SUBROUTINE init_aerosol_props_tegen_ecrad
  !---------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------
  !>
  !! SUBROUTINE init_tegen_scal_factors
  !! Constructor for t_tegen_scal_factors
  !!
  !! @par Revision History
  !! Initial release by Daniel Rieger, Deutscher Wetterdienst, Offenbach (2021-10-20)
  !!
  !---------------------------------------------------------------------------------------
  SUBROUTINE init_tegen_scal_factors(this, n_bands)
    CLASS(t_tegen_scal_factors), INTENT(inout) :: &
      &  this      !< Scaling factor information
    INTEGER, INTENT(in) :: &
      &  n_bands   !< Total number of wavelength bands

    this%n_bands = n_bands

    ALLOCATE(this%absorption(this%n_bands,5))
    ALLOCATE(this%scattering(this%n_bands,5))
    ALLOCATE(this%asymmetry (this%n_bands,5))
    !$ACC ENTER DATA CREATE(this%absorption, this%scattering, this%asymmetry)
  END SUBROUTINE
  !---------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------
  !>
  !! SUBROUTINE finalize_tegen_scal_factors
  !! Destructor for t_tegen_scal_factors
  !!
  !! @par Revision History
  !! Initial release by Daniel Rieger, Deutscher Wetterdienst, Offenbach (2021-10-20)
  !!
  !---------------------------------------------------------------------------------------
  SUBROUTINE finalize_tegen_scal_factors(this)
    CLASS(t_tegen_scal_factors), INTENT(inout) :: &
      &  this      !< Scaling factor information

    this%n_bands = 0

    !$ACC EXIT DATA DELETE(this%absorption) IF(ALLOCATED(this%absorption))
    IF (ALLOCATED(this%absorption)) &
      DEALLOCATE(this%absorption)
    !$ACC EXIT DATA DELETE(this%scattering) IF(ALLOCATED(this%scattering))
    IF (ALLOCATED(this%scattering)) &
      DEALLOCATE(this%scattering)
    !$ACC EXIT DATA DELETE(this%asymmetry) IF(ALLOCATED(this%asymmetry))
    IF (ALLOCATED(this%asymmetry)) &
      DEALLOCATE(this%asymmetry)
  END SUBROUTINE
  !---------------------------------------------------------------------------------------

  ! Very simple parameterization of source and sink terms for prognostic 2D aerosol fields
  !
  SUBROUTINE prog_aerosol_2D (jcs, jce, jg, dtime, iprog_aero, aerosol,               &
    &                         aercl_ss,aercl_or,aercl_bc,aercl_su,aercl_du,           &
    &                         rr_gsp,sr_gsp,rr_con,sr_con,                            &
    &                         soiltype,plcov_t,frac_t,w_so_t, w_so_ice_t, h_snow_t,   &
    &                         t_seasfc, lc_class_t, rho, tcm_t, u, v, sp_10m,         &
    &                         idx_lst_t, gp_count_t , i_count_sea, idx_sea)
    REAL(wp), INTENT(in)            :: &
      &  dtime,                        & !< Time step (s)
      &  aercl_ss(:), aercl_or(:),     & !< AOD climatology (sea salt, organic)
      &  aercl_bc(:), aercl_su(:),     & !< AOD climatology (black carbon, sulfate)
      &  aercl_du(:),                  & !< AOD climatology (dust)
      &  rr_gsp(:),sr_gsp(:),          & !< Grid-scale rain & snow rate
      &  rr_con(:),sr_con(:),          & !< Convective rain & snow rate
      &  plcov_t(:,:),                 & !< Plant cover (tiled)
      &  frac_t(:,:),                  & !< Tile fraction
      &  w_so_t(:,:), w_so_ice_t(:,:), & !< Soil water & ice (tiled)
      &  h_snow_t(:,:),                & !< Snow height (tiled)
      &  t_seasfc(:),                  & !< Sea surface temperature
      &  rho(:),                       & !< Air density
      &  tcm_t(:,:),                   & !< Transfer coefficient for momentum
      &  u(:), v(:),                   & !< Wind vector components
      &  sp_10m(:)                       !< Wind speed in 10m
    INTEGER,  INTENT(in) :: &
      &  jcs, jce,          & !< Start and end index of nproma loop
      &  jg,                & !< Domain index
      &  iprog_aero,        & !< Prognostic aerosol mode: 1 only dust, 2 all
      &  soiltype(:),       & !< Soil type index (dim: nproma)
      &  lc_class_t(:,:),   & !< Land use class index (dim: nproma, ntiles)
      &  idx_lst_t(:,:),    & !< Tiled index list to loop over land points (dim: nproma,ntiles)
      &  gp_count_t(:),     & !< Returns number of local grid points per tile (dim: ntiles)
      &  i_count_sea,       & !< Number of open water points in current block
      &  idx_sea(:)           !< Indices of open water points in current block
    REAL(wp), INTENT(inout) :: &
      &  aerosol(:,:)         !< Aerosol Optical Depth (AOD)
    ! Local variables
    REAL(wp) ::                 &
      &  relax_scale(2),        & !< Target values for relaxation
      &  relax_fac(2),          & !< Relaxation time scales
      &  od_clim(nclass_aero),  & !< AOD offsets for climatology-based source terms
      &  minfrac,               & !< minimum allowed fraction of climatological AOD
      &  ts_orgsrc,             & !< Time scales for sources
      &  ts_bcsrc, ts_susrc,    & !< Time scales for sources
      &  washout, washout_scale,& !< Washout and washout scale for dust
      &  dust_flux,             & !< Flux of dust
      &  aod_flux                 !< Source function for aerosol optical depth
    INTEGER ::              &
      &  jc, jt, jcl,       & !< Loop indices
      &  i_count_lnd          !< Number of land grid points in current block

    relax_scale(1) = 0.7_wp ! tuned to approximately balance 
                            ! the source terms for non-dust aerosol classes
    relax_scale(2) = 1._wp
    relax_fac(1)   = 1._wp/(3._wp*86400._wp)  ! 3 days for aerosols with shallow vertical extent
    relax_fac(2)   = 1._wp/(8._wp*86400._wp)  ! 8 days for aerosols with deep vertical extent (dust)
    ts_orgsrc      = 1._wp/(2.5_wp*86400._wp) ! 2.5  days for organic aerosol
    ts_bcsrc       = 1._wp/(2.5_wp*86400._wp) ! 2.5  days for black carbon
    ts_susrc       = 1._wp/(2.5_wp*86400._wp) ! 2.5  days for sulfate aerosol
    washout_scale  = 1._wp/7.5_wp             ! e-folding scale 7.5 mm WE precipitation
    minfrac        = 0.025_wp
    od_clim(iorg)  = 0.015_wp
    od_clim(ibc)   = 0.002_wp
    od_clim(iso4)  = 0.015_wp

    ! Prediction of mineral dust; other aerosol classes are treated prognostically only if iprog_aero=2

    ! Relaxation to scaled climatology
    DO jc = jcs, jce
      aerosol(jc,idu)  = aerosol(jc,idu)  + dtime*relax_fac(2)*(relax_scale(2)*aercl_du(jc)-aerosol(jc,idu))
    ENDDO

    DO jt = 1, ntiles_lnd
      i_count_lnd = gp_count_t(jt)
      IF (i_count_lnd == 0) CYCLE ! skip loop if the index list for the given tile is empty
!$NEC ivdep
      DO jcl = 1, i_count_lnd
        jc = idx_lst_t(jcl,jt)
        ! dust_flux is not used here, but could be used for more sophisticated aerosol modules
        CALL aerosol_dust_aod_source (p_dust_source_const(jg), dzsoil(1), w_so_t(jc,jt), h_snow_t(jc,jt), &
          &                           w_so_ice_t(jc,jt), soiltype(jc), plcov_t(jc,jt), lc_class_t(jc,jt), &
          &                           rho(jc), tcm_t(jc,jt), u(jc), v(jc), aod_flux, dust_flux)
        ! Update AOD field with tendency from aod_flux
        aerosol(jc,idu) = aerosol(jc,idu) + aod_flux * frac_t(jc,jt) * dtime
      ENDDO ! jcl
    ENDDO !jt

    DO jc = jcs, jce
      ! Calculate washout everywhere (not only above land)
      washout = dtime*washout_scale*(rr_gsp(jc)+sr_gsp(jc)+0.5_wp*(rr_con(jc)+sr_con(jc)))*aerosol(jc,idu)
      aerosol(jc,idu)  = aerosol(jc,idu) - washout
      ! Ensure that the aerosol optical depth does not fall below 2.5% of the climatological value
      aerosol(jc,idu)  = MAX(aerosol(jc,idu),  minfrac*aercl_du(jc))
    ENDDO

    IF (iprog_aero == 2) THEN

      ! Relaxation to scaled climatology
      DO jc = jcs, jce
        aerosol(jc,iss)  = aerosol(jc,iss)  + dtime*relax_fac(1)*(relax_scale(2)*aercl_ss(jc)-aerosol(jc,iss))
        aerosol(jc,iorg) = aerosol(jc,iorg) + dtime*relax_fac(1)*(relax_scale(1)*aercl_or(jc)-aerosol(jc,iorg))
        aerosol(jc,ibc)  = aerosol(jc,ibc)  + dtime*relax_fac(1)*(relax_scale(1)*aercl_bc(jc)-aerosol(jc,ibc))
        aerosol(jc,iso4) = aerosol(jc,iso4) + dtime*relax_fac(1)*(relax_scale(1)*aercl_su(jc)-aerosol(jc,iso4))
      ENDDO

      ! Sources are specified where either the climatological value exceeds the global average
      ! sources of sea salt and black carbon are restricted to water and land surfaces, respectively
      DO jc = jcs, jce
        aerosol(jc,iso4) = aerosol(jc,iso4) + dtime*ts_susrc*MAX(0._wp,aercl_su(jc)-MAX(od_clim(iso4),0.3_wp*aerosol(jc,iso4)))
        aerosol(jc,iorg) = aerosol(jc,iorg) + dtime*ts_orgsrc*MAX(0._wp,aercl_or(jc)-MAX(od_clim(iorg),0.3_wp*aerosol(jc,iorg)))
        IF (soiltype(jc) /= 9) THEN ! land
          aerosol(jc,ibc)  = aerosol(jc,ibc)  + dtime*ts_bcsrc* MAX(0._wp,aercl_bc(jc)-od_clim(ibc))
        ENDIF
      ENDDO
      
      ! New SSA source
!$NEC ivdep
      DO jcl = 1, i_count_sea
        jc = idx_sea(jcl)
        CALL aerosol_ssa_aod_source (t_seasfc(jc), sp_10m(jc), aod_flux)
        aerosol(jc,iss) = aerosol(jc,iss) + aod_flux * frac_t(jc,isub_water) * dtime
      ENDDO

      DO jc = jcs, jce
        ! Washout (for seasalt only)
        washout = dtime*washout_scale*(rr_gsp(jc)+sr_gsp(jc)+0.5_wp*(rr_con(jc)+sr_con(jc)))*aerosol(jc,iss)
        aerosol(jc,iss)  = aerosol(jc,iss) - washout
        ! Ensure that the aerosol optical depth does not fall below 2.5% of the climatological value
        aerosol(jc,iss)  = MAX(aerosol(jc,iss),  minfrac*aercl_ss(jc))
        aerosol(jc,iorg) = MAX(aerosol(jc,iorg), minfrac*aercl_or(jc))
        aerosol(jc,ibc)  = MAX(aerosol(jc,ibc),  minfrac*aercl_bc(jc))
        aerosol(jc,iso4) = MAX(aerosol(jc,iso4), minfrac*aercl_su(jc))
      ENDDO

    ENDIF

  END SUBROUTINE prog_aerosol_2D


  ! Tuning of longwave absorption coefficient of mineral dust in order to reduce cold bias in the Saharan region
  !
  SUBROUTINE tune_dust (lat,lon,iend,tunefac)

    REAL(wp), INTENT(in) :: lat(:), lon(:)
    INTEGER,  INTENT(in) :: iend

    REAL(wp), INTENT(out) :: tunefac(:,:)

    INTEGER :: jc, jb
    REAL(wp) :: maxfac

    DO jb = 1, jpband
      maxfac = tune_dust_abs*5._wp*(jpband-MAX(8,jb))/REAL(jpband-8,wp)
      DO jc = 1, iend
        tunefac(jc,jb) = 1._wp + maxfac*(1._wp - MIN(1._wp,((rad2deg*lat(jc)-15._wp)/20._wp)**4)) * &
         (1._wp - MIN(1._wp,((rad2deg*lon(jc)-20._wp)/50._wp)**4))
      ENDDO
    ENDDO


  END SUBROUTINE tune_dust


  SUBROUTINE aerosol_2D_diffusion( p_patch, p_int_state, nproma, aerosol )
    TYPE(t_patch), INTENT(in)     :: &
      &  p_patch                       !< Current patch
    TYPE(t_int_state), INTENT(in) :: &
      &  p_int_state                   !< interpolation state
    INTEGER,  INTENT(in)          :: &
      &  nproma
    REAL(wp), INTENT(inout)       :: &
      &  aerosol(:,:,:)                !< Aerosol container
    ! Local variables
    REAL(wp)                      :: &
      &  diff_coeff,                 & !< Diffusion coefficient
      &  nabla2_aero(nproma,nclass_aero,p_patch%nblks_c) !< Laplacian of aerosol(:,:,:)
    INTEGER                       :: &
      &  jb, jc,                     &
      &  i_rlstart, i_rlend,         &
      &  i_startblk, i_endblk,       & 
      &  i_startidx, i_endidx

    diff_coeff = 0.125_wp

    CALL nabla2_scalar(aerosol(:,:,:),          &
      &                p_patch, p_int_state,    &
      &                nabla2_aero(:,:,:),      &
      &                iss, idu, grf_bdywidth_c+1, min_rlcell_int)

    i_rlstart  = grf_bdywidth_c+1
    i_rlend    = min_rlcell_int
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, i_rlstart, i_rlend)
      DO jc = i_startidx, i_endidx
        aerosol(jc,idu,jb) = MAX(0.0_wp, aerosol(jc,idu,jb) + diff_coeff *           &
                                 p_patch%cells%area(jc,jb) * nabla2_aero(jc,idu,jb))
        aerosol(jc,iss,jb) = MAX(0.0_wp, aerosol(jc,iss,jb) + diff_coeff *           &
                                 p_patch%cells%area(jc,jb) * nabla2_aero(jc,iss,jb))
      ENDDO !jc
    ENDDO !jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE aerosol_2D_diffusion

END MODULE mo_aerosol_util

