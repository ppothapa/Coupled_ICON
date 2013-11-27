!>
!! Namelist for the configuration of interpolation and reconstruction.
!!
!! @par Revision History
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
!!
MODULE mo_interpol_nml

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish
  USE mo_impl_constants,      ONLY: max_dom
  USE mo_mpi,                 ONLY: my_process_is_stdio 
  USE mo_master_control,      ONLY: is_restart_run
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_namelist,            ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_io_restart_namelist, ONLY: open_tmpfile, store_and_close_namelist,  &
                                  & open_and_restore_namelist, close_tmpfile

  USE mo_interpol_config,     ONLY: config_llsq_lin_consv    => llsq_lin_consv    , &
                                  & config_llsq_high_consv   => llsq_high_consv   , &
                                  & config_lsq_high_ord      => lsq_high_ord      , &
                                  & config_rbf_vec_kern_c    => rbf_vec_kern_c    , &
                                  & config_rbf_vec_scale_c   => rbf_vec_scale_c   , &
                                  & config_rbf_vec_kern_v    => rbf_vec_kern_v    , &
                                  & config_rbf_vec_scale_v   => rbf_vec_scale_v   , &
                                  & config_rbf_vec_kern_e    => rbf_vec_kern_e    , &
                                  & config_rbf_vec_scale_e   => rbf_vec_scale_e   , &
                                  & config_rbf_vec_kern_ll   => rbf_vec_kern_ll   , &
                                  & config_rbf_vec_scale_ll  => rbf_vec_scale_ll  , &
                                  & config_i_cori_method     => i_cori_method     , &
                                  & config_nudge_max_coeff   => nudge_max_coeff   , &
                                  & config_nudge_efold_width => nudge_efold_width , &
                                  & config_nudge_zone_width  => nudge_zone_width  , &
                                  & config_l_corner_vort     => l_corner_vort     , &
                                  & config_l_intp_c2l        => l_intp_c2l        , &
                                  & config_rbf_dim_c2l       => rbf_dim_c2l       , &
                                  & config_l_mono_c2l        => l_mono_c2l
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_interpol_namelist

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  !--------------------
  ! namelist variables
  !--------------------

  LOGICAL  :: llsq_lin_consv      ! conservative (TRUE) or non-conservative (FALSE)
                                  ! linear least squares reconstruction
  LOGICAL  :: llsq_high_consv     ! conservative (TRUE) or non-conservative (FALSE)
                                  ! high order least squares reconstruction 

  INTEGER  :: lsq_high_ord        ! specific order for higher order lsq

  INTEGER  :: rbf_vec_kern_c,   & ! parameter determining the type
              rbf_vec_kern_v,   & ! of vector rbf kernel
              rbf_vec_kern_e,   &
              rbf_vec_kern_ll

  ! Parameter fields determining the scale factor used by the vector rbf
  ! interpolator.
  ! Note: these fields are defined on each grid level; to allow the namelist input
  ! going from 1 to depth (rather than from start_lev to end_lev), the namelist input
  ! fields defined here differ from those used in the model

  REAL(wp) :: rbf_vec_scale_c(max_dom),  &
              rbf_vec_scale_v(max_dom),  &
              rbf_vec_scale_e(max_dom),  &
              rbf_vec_scale_ll(max_dom)

  INTEGER  :: i_cori_method       ! Identifier for the method with wich the tangential
                                  ! wind reconstruction in Coriolis force is computed,
                                  ! if the Thuburn method is used. (To be
                                  ! implemented for triangles, currently only for
                                  ! hexagons)
                                  ! i_cori_method = 1 : Almut's method for reconstruction
                                  !                     but TRSK method for PV
                                  ! i_cori_method = 2 : Thuburn/Ringler/Skamarock/Klemp
                                  ! i_cori_method = 3 : Almut's method for reconstruction
                                  !                     Almut's method also for PV
                                  ! i_cori_method = 4 : Almut's method for reconstruction, 
                                  !                     but PV on averaged on vertices

  ! Namelist variables setting up the lateral boundary nudging (applicable to limited-area
  ! runs and one-way nesting). The nudging coefficients start with nudge_max_coeff in
  ! the cell row bordering to the boundary interpolation zone, and decay exponentially
  ! with nudge_efold_width (in units of cell rows)

  REAL(wp) :: nudge_max_coeff, nudge_efold_width
  INTEGER  :: nudge_zone_width    ! total width of nudging zone in units of cell rows

  LOGICAL :: l_corner_vort        ! yields for i_cori_method>=3
                                  ! Decision wheter the hexagon vector reconstruction is
                                  ! combined with either of the two vorticities :
                                  ! .TRUE. : Three rhombi are combined to the corner
                                  !          and afterwards averaged to the hexagon center
                                  ! .FALSE.: 6 rhombi are directly averaged to the
                                  !          hexagon center (original method).  
                                  ! After the writing of the paper to be published in JCP 
                                  ! it seems that l_corner_vort=.TRUE. should be the right way.

  LOGICAL :: l_intp_c2l, l_mono_c2l
  INTEGER :: rbf_dim_c2l

  NAMELIST/interpol_nml/ llsq_lin_consv,    llsq_high_consv,     &
                       & lsq_high_ord,      rbf_vec_kern_c,      &
                       & rbf_vec_scale_c,   rbf_vec_kern_v,      &
                       & rbf_vec_scale_v,   rbf_vec_kern_e,      &
                       & rbf_vec_scale_e,   i_cori_method,       &
                       & nudge_max_coeff,   nudge_efold_width,   &
                       & nudge_zone_width,  l_corner_vort,       &
                       & l_intp_c2l, rbf_dim_c2l, l_mono_c2l,    &
                       & rbf_vec_scale_ll,  rbf_vec_kern_ll

CONTAINS
  !-------------------------------------------------------------------------
  !>
  !! Read Namelist for interpolation. 
  !!
  !! This subroutine 
  !! - reads the Namelist for interpolation
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a 
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state (partly)    
  !!
  !! @par Revision History
  !!  by Daniel Reinert, DWD (2011-07-06)
  !!  Kristina Froehlich and Hui Wan, MPI-M (2011-07)
  !!
  SUBROUTINE read_interpol_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_interpol_nml: read_interpol_namelist'

    !-----------------------
    ! 1. default settings   
    !-----------------------
    ! LSQ reconstruction at cell center
    llsq_lin_consv   = .FALSE.  ! non-conservative linear reconstruction
    llsq_high_consv  = .TRUE.   ! conservative high order reconstruction
    lsq_high_ord     = 3        ! cubic polynomial

    ! Kernals for RBF vector reconstruction
    rbf_vec_kern_c  = 1         ! Gaussian kernel for cell centers
    rbf_vec_kern_v  = 1         ! Gaussian kernel for vertices
    rbf_vec_kern_e  = 3         ! Inverse multiquadric kernel for edge midpoints
    rbf_vec_kern_ll = 1         ! Gaussian kernel for lon-lat interpolation

    ! Initialize namelist fields for scaling factors (dimension 1:depth) with dummy values
    ! A meaningful initialization follows after reading the namelist
    rbf_vec_scale_c(:)  = -1.0_wp
    rbf_vec_scale_v(:)  = -1.0_wp
    rbf_vec_scale_e(:)  = -1.0_wp
    rbf_vec_scale_ll(:) = -1.0_wp

    ! Initialize the namelist for the method for the vorticity flux term
    i_cori_method = 3
    l_corner_vort=.TRUE.

    ! Coefficients for lateral boundary nudging
    nudge_max_coeff   = 0.02_wp  ! Maximum nudging coefficient
    nudge_efold_width = 2._wp    ! e-folding width in units of cell rows
    nudge_zone_width  = 8        ! Width of nudging zone in units of cell rows

    ! direct interpolation from cell centers to lon-lat points:
    l_intp_c2l   = .TRUE.
    ! stencil size: 4  = nearest neighbor, 
    !               13 = vertex stencil,
    !               rbf_c2grad_dim = edge stencil
    rbf_dim_c2l  = 10
    ! no monotonicity cutoff by default:
    l_mono_c2l   = .TRUE.

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (is_restart_run()) THEN
      funit = open_and_restore_namelist('interpol_nml')
      READ(funit,NML=interpol_nml)
      CALL close_tmpfile(funit)
    END IF

    !-------------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !-------------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('interpol_nml', status=istat)
    IF (my_process_is_stdio()) WRITE(temp_defaults(), interpol_nml) ! write defaults to temporary text file
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, interpol_nml)                                     ! overwrite default settings
      IF (my_process_is_stdio()) WRITE(temp_settings(), interpol_nml) ! write settings to temporary text file
    END SELECT
    CALL close_nml

    !-----------------------------------------------------------------------
    ! Sanity check 
    !-----------------------------------------------------------------------
    IF ((rbf_vec_kern_c/=1 ).AND.(rbf_vec_kern_c/=3)) THEN
      CALL finish( TRIM(routine),'wrong value of rbf_vec_kern_c, must be 1 or 3')
    ENDIF

    IF (.NOT.((rbf_vec_kern_v==1 ).OR.(rbf_vec_kern_v==3))) THEN
      CALL finish( TRIM(routine),'wrong value of rbf_vec_kern_v, must be 1 or 3')
    ENDIF

    IF (.NOT.((rbf_vec_kern_e==1 ).OR.(rbf_vec_kern_e==3))) THEN
      CALL finish( TRIM(routine),'wrong value of rbf_vec_kern_e, must be 1 or 3')
    ENDIF

    IF ((rbf_vec_kern_ll/=1 ).AND.(rbf_vec_kern_ll/=3)) THEN
      CALL finish( TRIM(routine),'wrong value of rbf_vec_kern_ll, must be 1 or 3')
    ENDIF

    SELECT CASE(lsq_high_ord)
    CASE(1,2,30,3)  ! OK
    CASE DEFAULT
      CALL finish( TRIM(routine),'wrong value of lsq_high_ord, must be 1,2,30 or 3')
    END SELECT

    !----------------------------------------------------
    ! 4. Fill the configuration state
    !----------------------------------------------------
    config_llsq_lin_consv      = llsq_lin_consv
    config_llsq_high_consv     = llsq_high_consv
    config_lsq_high_ord        = lsq_high_ord 
    config_rbf_vec_kern_c      = rbf_vec_kern_c 
    config_rbf_vec_kern_v      = rbf_vec_kern_v
    config_rbf_vec_kern_e      = rbf_vec_kern_e
    config_rbf_vec_kern_ll     = rbf_vec_kern_ll
    config_rbf_vec_scale_c(:)  = rbf_vec_scale_c(:)
    config_rbf_vec_scale_v(:)  = rbf_vec_scale_v(:)
    config_rbf_vec_scale_e(:)  = rbf_vec_scale_e(:)
    config_rbf_vec_scale_ll(:) = rbf_vec_scale_ll(:)
    config_i_cori_method       = i_cori_method
    config_nudge_max_coeff     = nudge_max_coeff
    config_nudge_efold_width   = nudge_efold_width
    config_nudge_zone_width    = nudge_zone_width
    config_l_corner_vort       = l_corner_vort
    config_l_intp_c2l          = l_intp_c2l
    config_rbf_dim_c2l         = rbf_dim_c2l
    config_l_mono_c2l          = l_mono_c2l

    !-----------------------------------------------------
    ! 5. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=interpol_nml)
      CALL store_and_close_namelist(funit, 'interpol_nml')
    ENDIF

    ! 6. write the contents of the namelist to an ASCII file
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=interpol_nml)

  END SUBROUTINE read_interpol_namelist

END MODULE mo_interpol_nml
