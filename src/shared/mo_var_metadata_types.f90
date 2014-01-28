MODULE mo_var_metadata_types

  USE mo_kind,           ONLY: dp, wp
  USE mo_grib1,          ONLY: t_grib1_var
  USE mo_grib2,          ONLY: t_grib2_var
  USE mo_action_types,   ONLY: t_var_action
  USE mo_cf_convention,  ONLY: t_cf_var

  IMPLICIT NONE

  PRIVATE

  ! ---------------------------------------------------------------
  ! CONSTANTS
  ! ---------------------------------------------------------------


  ! maximum string length for variable names
  INTEGER, PARAMETER :: VARNAME_LEN = 32

  ! list of variable groups
  ! 
  ! A variable can have any combination of this which means that it is
  ! part of each of these different variable sets.
  CHARACTER(len=VARNAME_LEN), PARAMETER :: var_groups(33) = &
    (/ "ALL                   ",  &
    &  "ATMO_ML_VARS          ",  &
    &  "ATMO_PL_VARS          ",  &
    &  "ATMO_ZL_VARS          ",  &
    &  "NH_PROG_VARS          ",  &
    &  "ATMO_DERIVED_VARS     ",  &
    &  "RAD_VARS              ",  &
    &  "PRECIP_VARS           ",  &
    &  "CLOUD_DIAG            ",  &
    &  "PBL_VARS              ",  &
    &  "PHYS_TENDENCIES       ",  &
    &  "LAND_VARS             ",  &
    &  "LAND_TILE_VARS        ",  &
    &  "MULTISNOW_VARS        ",  &
    &  "ADDITIONAL_PRECIP_VARS",  &
    &  "SNOW_VARS             ",  &
    &  "DWD_FG_ATM_VARS       ",  &  ! DWD First Guess (atmosphere) 
    &  "DWD_FG_SFC_VARS       ",  &  ! DWD First Guess (surface/soil)
    &  "MODE_DWD_FG_IN        ",  &  ! Input first guess fields for MODE_DWD
    &  "MODE_DWD_ANA_IN       ",  &  ! Input analysis fields for MODE_DWD
    &  "MODE_COMBINED_IN      ",  &  ! Input fields for MODE_COMBINED
    &  "MODE_COSMODE_IN       ",  &  ! Input fields for MODE_COSMODE
    &  "OCE_PROG              ",  &
    &  "OCE_DIAG              ",  &
    &  "OCE_DEFAULT           ",  &
    &  "oce_essentials        ",  &
    &  "oce_force_essentials  ",  &
    &  "OCE_AUX               ",  &
    &  "OCE_GEOMETRY          ",  &
    &  "OCE_PHYSICS           ",  &
    &  "OCE_COEFFS            ",  &
    &  "ICE_DEFAULT           ",  &
    &  "ICE_DIAG              "/)


  ! list of vertical interpolation types
  ! 
  ! A variable can have any combination of this which means that it
  ! can be interpolated vertically in these different ways.
  CHARACTER(len=VARNAME_LEN), PARAMETER :: VINTP_TYPE_LIST(3) = &
    (/ "Z                     ",  &
    &  "P                     ",  &
    &  "I                     " /)


  ! list of available post-op's (small arithmetic operations on
  ! fields). The implementation is placed in "mo_post_op.f90".o
  INTEGER, PARAMETER, PUBLIC   :: POST_OP_NONE      = -1  !< trivial post-op ("do nothing")
  INTEGER, PARAMETER, PUBLIC   :: POST_OP_SCALE     =  1  !< multiply by scalar factor "arg1"
  INTEGER, PARAMETER, PUBLIC   :: POST_OP_RHO       =  2  !< multiply by rho to get densities instead

  ! ---------------------------------------------------------------
  ! TYPE DEFINITIONS
  ! ---------------------------------------------------------------


  TYPE t_union_vals
    REAL(dp) :: rval
    INTEGER  :: ival
    LOGICAL  :: lval
  END type t_union_vals


  TYPE t_tracer_meta
    !
    LOGICAL :: lis_tracer         ! this is a tracer field (TRUE/FALSE)
    CHARACTER(len=VARNAME_LEN) :: tracer_class ! type of tracer
    !  
    INTEGER :: ihadv_tracer       ! method for horizontal transport
    INTEGER :: ivadv_tracer       ! method for vertical transport
    !
    LOGICAL :: lturb_tracer       ! turbulent transport (TRUE/FALSE)
    LOGICAL :: lsed_tracer        ! sedimentation (TRUE/FALSE)
    LOGICAL :: ldep_tracer        ! dry deposition (TRUE/FALSE)  
    LOGICAL :: lconv_tracer       ! convection  (TRUE/FALSE)
    LOGICAL :: lwash_tracer       ! washout (TRUE/FALSE)
    !
    REAL(wp) :: rdiameter_tracer   ! particle diameter in m
    REAL(wp) :: rrho_tracer        ! particle density in kg m^-3
    !
    REAL(wp) :: halflife_tracer   ! radioactive half-life in s^-1
    INTEGER  :: imis_tracer     ! IMIS number
    !
  END TYPE t_tracer_meta


  !> data specific for pz-level interpolation.
  TYPE t_vert_interp_meta
    ! meta data containing the groups to which a variable belongs
    LOGICAL  :: vert_intp_type(SIZE(VINTP_TYPE_LIST))
    INTEGER  :: vert_intp_method
    LOGICAL  :: l_hires_intp, l_restore_fricred, l_loglin, &
         &      l_extrapol, l_satlimit, l_restore_pbldev,  &
         &      l_pd_limit, l_restore_sfcinv, l_hires_corr
    REAL(wp) :: lower_limit, extrapol_dist
  END TYPE t_vert_interp_meta


  !> data specific for horizontal interpolation.
  TYPE t_hor_interp_meta
    INTEGER :: hor_intp_type ! NONE/LONLAT
    INTEGER :: lonlat_id     ! lon-lat grid (ID in global list)
  END TYPE t_hor_interp_meta


  !> This type defines small arithmetic operations ("post-ops") as
  !  post-processing tasks.
  !
  !  These post-processing tasks are restricted to point-wise
  !  operations (no halo synchronization) of a single field, like
  !  value scaling.
  !
  !  @note The "post-ops" are performed at output time and DO NOT
  !        MODIFY THE FIELD ITSELF.
  !
  TYPE t_post_op_meta
    INTEGER                    :: ipost_op_type         !< type of post-processing operation
    !
    LOGICAL                    :: lnew_cf
    TYPE(t_cf_var)             :: new_cf                !< CF information of modified field
    LOGICAL                    :: lnew_grib2
    TYPE(t_grib2_var)          :: new_grib2             !< GRIB2 information of modified field
    !
    REAL(wp)                   :: arg1                  !< post-op argument (e.g. scaling factor)
  END TYPE t_post_op_meta


  TYPE t_var_metadata
    !
    INTEGER                    :: key                   ! hash value of name
    CHARACTER(len=VARNAME_LEN) :: name                  ! variable name  
    !
    TYPE(t_cf_var)             :: cf                    ! CF convention information 
    TYPE(t_grib2_var)          :: grib2                 ! GRIB2 related information
    !
    LOGICAL                    :: allocated             ! allocation status
    INTEGER                    :: ndims                 ! number of dimensions used
    INTEGER                    :: used_dimensions(5)    ! final dimensions of variable
    ! 
    LOGICAL                    :: lrestart              ! write field to restart
    LOGICAL                    :: loutput               ! write field to output
    INTEGER                    :: isteptype             ! Type of statistical processing
    !                                         
    TYPE(t_union_vals)         :: resetval              ! reset value for accumulated fields
    LOGICAL                    :: lmiss                 ! missing value flag
    TYPE(t_union_vals)         :: missval               ! missing value
    LOGICAL                    :: lrestart_cont         ! continue if not in restart file     
    LOGICAL                    :: lrestart_read         ! field has been set from restart file
    TYPE(t_union_vals)         :: initval               ! value if not in restart file
    !     
    LOGICAL                    :: lcontainer            ! true, if this is a container
    LOGICAL                    :: lcontained            ! true, if this is in a container
    INTEGER                    :: ncontained            ! index in container   
    !
    INTEGER                    :: hgrid                 ! CDI horizontal grid type
    INTEGER                    :: vgrid                 ! CDI vertical grid type
    !
    INTEGER                    :: tlev_source           ! Information where to find the actual
    !                                                     timelevel for timelevel dependent variables:        
    !                                                      = 0 : nnow
    !                                                      = 1 : nnow_rcf
    !                                                      ... more may follow
    !
    INTEGER                    :: cdiVarID
    INTEGER                    :: cdiVarID_2            ! for 2nd vector component in LatLon interpolation
    INTEGER                    :: cdiGridID
    INTEGER                    :: cdiZaxisID
    INTEGER                    :: cdiDataType
    !
    TYPE(t_tracer_meta)        :: tracer                ! metadata for tracer fields
    !
    ! Metadata for "post-ops" (small arithmetic operations)
    !
    TYPE(t_post_op_meta)       :: post_op               !<  "post-op" (small arithmetic operations) for this variable
    !
    ! Metadata for "actions" (regularly triggered events)
    !
    TYPE(t_var_action)         :: action_list
    !
    ! Metadata for vertical/horizontal interpolation
    !
    ! Note that setting these parameters to non-default values does
    ! not mean that interpolation is actually performed for this
    ! variables (this is controlled by namelist settings) but only
    ! that this is possible!
    !
    TYPE(t_vert_interp_meta)   :: vert_interp 
    TYPE(t_hor_interp_meta)    :: hor_interp 
    !
    ! meta data containing the groups to which a variable belongs
    LOGICAL :: in_group(SIZE(var_groups))

    ! Flag: defines, if this field is updated by the internal
    ! post-processing scheduler
    INTEGER :: l_pp_scheduler_task

  END TYPE t_var_metadata

  PUBLIC :: VINTP_TYPE_LIST
  PUBLIC :: VARNAME_LEN

  PUBLIC :: t_union_vals
  PUBLIC :: t_var_metadata
  PUBLIC :: t_tracer_meta
  PUBLIC :: t_vert_interp_meta
  PUBLIC :: t_hor_interp_meta
  PUBLIC :: t_post_op_meta

  PUBLIC :: var_groups

END MODULE mo_var_metadata_types