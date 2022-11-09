!>
!! <Short description of module for listings and indices>
!!
!! <Describe the concepts of the procedures and algorithms used in the module.>
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author <name, affiliation>
!! @author <name, affiliation>
!!
!!
!! @par Revision History
!! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_oem_config

  USE mo_kind,           ONLY: wp
  USE mo_io_units,       ONLY: filename_max
  USE mo_impl_constants, ONLY: MAX_CHAR_LENGTH

  IMPLICIT NONE
  PUBLIC

  !--------------------------------------------------------------------------
  ! Basic configuration setup for the online emission module
  !--------------------------------------------------------------------------

    !------------------------------------------------------------------------
    ! oemctrl_nml:
    !------------------------------------------------------------------------
    CHARACTER(LEN=filename_max) :: vertical_profile_nc,   & !< name of the oae vertical profile
                              &    hour_of_day_nc,        & !< name of the oae hour of day file
                              &    day_of_week_nc,        & !< name of the oae day of week file
                              &    month_of_year_nc,      & !< name of the oae month of year file
                              &    hour_of_year_nc,       & !< name of the oae hour of year file
                              &    gridded_emissions_nc,  & !< name of the oae gridded emission file
                              &    ens_reg_nc,            & !< name of file with ensemble-regions
                              &    ens_lambda_nc,         & !< name of file with ensemble-lambdas
                              &    vegetation_indices_nc    !< name of file with MODIS reflectances
    REAL(wp), DIMENSION(8) ::      vprm_par,              & !< VPRM parameter values for PAR_0
                              &    vprm_lambda,           & !< VPRM parameter values for lambda
                              &    vprm_alpha,            & !< VPRM parameter values for alpha
                              &    vprm_beta,             & !< VPRM parameter values for beta
                              &    vprm_tmin,             & !< VPRM parameter values for T_min
                              &    vprm_tmax,             & !< VPRM parameter values for T_max
                              &    vprm_topt,             & !< VPRM parameter values for T_opt
                              &    vprm_tlow                !< VPRM parameter values for T_low
    LOGICAL ::                     lcut_area                !< Switch to turn on/off to select an
                                                            !< area where no fluxes are applied
    REAL(wp) ::                    lon_cut_start,         & !< longitude start coordinate 
                              &    lon_cut_end,           & !< longitude end coordinate
                              &    lat_cut_start,         & !< latitude start coordinate
                              &    lat_cut_end              !< latitude end coordinate
 



END MODULE mo_oem_config

