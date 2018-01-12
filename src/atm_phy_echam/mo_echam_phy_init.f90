!>
!! @brief Contains subroutines for initializing the ECHAM physics package.
!!
!! @author Hui Wan, MPI-M
!!
!! @par Revision History
!! First version by Hui Wan, 2010-07-20
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

MODULE mo_echam_phy_init

  ! infrastructure
  USE mo_kind,                 ONLY: wp
  USE mo_exception,            ONLY: finish
  USE mtime,                   ONLY: datetime, OPERATOR(>)
  USE mo_io_config,            ONLY: default_read_method
  USE mo_read_interface,       ONLY: openInputFile, closeFile, read_2D, &
    &                                t_stream_id, on_cells
  USE mo_timer,                ONLY: timers_level, timer_start, timer_stop, &
    &                                timer_prep_echam_phy

  ! run configuration
  USE mo_run_config,           ONLY: nlev, iqv, iqt, ico2, io3, &
    &                                ntracer, ltestcase, lart

  ! horizontal grid and indices
  USE mo_model_domain,         ONLY: t_patch
  USE mo_loopindices,          ONLY: get_indices_c
  USE mo_grid_config,          ONLY: n_dom

  ! vertical grid
  USE mo_vertical_coord_table, ONLY: vct
  
  ! test cases
  USE mo_nh_testcases_nml,     ONLY: ape_sst_case, th_cbl, tpe_temp
  USE mo_ape_params,           ONLY: ape_sst
  USE mo_physical_constants,   ONLY: tmelt, Tf, albi, albedoW

  ! echam phyiscs
  USE mo_echam_phy_config,     ONLY: eval_echam_phy_config, eval_echam_phy_tc, print_echam_phy_config, &
    &                                echam_phy_config, echam_phy_tc, dt_zero
  USE mo_echam_phy_memory,     ONLY: construct_echam_phy_state,    &
    &                                prm_field, t_echam_phy_field, &
    &                                prm_tend,  t_echam_phy_tend

  ! radiation
  USE mo_echam_rad_config,     ONLY: eval_echam_rad_config, print_echam_rad_config, echam_rad_config
  USE mo_psrad_interface,      ONLY: setup_psrad
  USE mo_psrad_cloud_optics,   ONLY: setup_cloud_optics  

  ! subgrid scale orographic effects
  USE mo_echam_sso_config,     ONLY: eval_echam_sso_config, print_echam_sso_config

  ! atmospheric gravity wave drag
  USE mo_echam_gwd_config,     ONLY: eval_echam_gwd_config, print_echam_gwd_config

  ! vertical diffusion
  USE mo_echam_vdf_config,     ONLY: eval_echam_vdf_config, print_echam_vdf_config
  USE mo_echam_vdiff_params,   ONLY: init_vdiff_params
  USE mo_vdiff_solver,         ONLY: init_vdiff_solver

#ifndef __NO_JSBACH__
  ! land surface
  USE mo_master_control,       ONLY: master_namelist_filename
  USE mo_jsb_base,             ONLY: jsbach_init_base => init_base
  USE mo_jsb_model_init,       ONLY: jsbach_init_model => init_model
#endif

  ! cumulus convection
  USE mo_echam_cnv_config,     ONLY: alloc_echam_cnv_config, eval_echam_cnv_config, print_echam_cnv_config
  USE mo_convect_tables,       ONLY: init_convect_tables
  USE mo_echam_convect_tables, ONLY: init_echam_convect_tables => init_convect_tables

  ! cloud microphysics
  USE mo_echam_cld_config,     ONLY: print_echam_cld_config

  ! Cariolle interactive ozone scheme
  USE mo_lcariolle_externals,  ONLY: read_bcast_real_3d_wrap, &
    &                                read_bcast_real_1d_wrap, &
    &                                closeFile_wrap, openInputFile_wrap, &
    &                                get_constants

  ! water vapour production by methane oxidation
  ! and destruction by photolysis
  USE mo_methox,               ONLY: init_methox

  ! air-sea-land interface
  USE mo_echam_sfc_indices,    ONLY: nsfc_type, iwtr, iice, ilnd, init_sfc_indices

  ! for coupling
  USE mo_coupling_config,      ONLY: is_coupled_run

  ! for AMIP boundary conditions
  USE mo_bcs_time_interpolation, ONLY: t_time_interpolation_weights, calculate_time_interpolation_weights
  USE mo_bc_sst_sic,           ONLY: read_bc_sst_sic, bc_sst_sic_time_interpolation
  USE mo_bc_greenhouse_gases,  ONLY: read_bc_greenhouse_gases, bc_greenhouse_gases_time_interpolation, &
    &                                bc_greenhouse_gases_file_read, ghg_co2mmr
  USE mo_bc_aeropt_splumes,    ONLY: setup_bc_aeropt_splumes

  ! radiative forcing diagnostics
  USE mo_psrad_memory,         ONLY: construct_psrad_forcing_list

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: init_echam_phy, initcond_echam_phy
  PUBLIC  :: additional_restart_init

CONTAINS
  !>
  !! Top-level routine for initialization of ECHAM6 physics.
  !! It calls a series of subroutines to initialize tunable parameters,
  !! lookup tables, and the physics state vectors "prm_field" and "prm_tend".
  !!
  !! @par Revision History
  !! Initial version by Hui Wan, MPI-M (2010-07)
  !! name change to init_echam_phy by Levi Silvers
  !!
  SUBROUTINE init_echam_phy( p_patch, ctest_name, nlev, mtime_current)

    TYPE(t_patch), TARGET, INTENT(in) :: p_patch(:)
    CHARACTER(LEN=*),INTENT(in) :: ctest_name
    INTEGER,         INTENT(in) :: nlev
    TYPE(datetime),  INTENT(in), POINTER    :: mtime_current !< Date and time information

    INTEGER :: khydromet, ktrac
    INTEGER :: jg
    TYPE(t_stream_id) :: stream_id

    CHARACTER(len=*), PARAMETER :: land_frac_fn = 'bc_land_frac.nc'
    CHARACTER(len=*), PARAMETER :: land_phys_fn = 'bc_land_phys.nc'
    CHARACTER(len=*), PARAMETER :: land_sso_fn  = 'bc_land_sso.nc'

    TYPE(t_time_interpolation_weights) :: current_time_interpolation_weights

    LOGICAL :: lany

    ! Shortcuts to components of echam_rad_config
    !
    INTEGER, POINTER :: ighg(:), irad_aero(:)
    !
    ighg      => echam_rad_config(1:n_dom)% ighg
    irad_aero => echam_rad_config(1:n_dom)% irad_aero
    
    IF (timers_level > 1) CALL timer_start(timer_prep_echam_phy)

    !-------------------------------------------------------------------
    ! Initialize parameters and lookup tables
    !-------------------------------------------------------------------
    ! Evaluate the ECHAM physics configuration variables echam_phy_config(:)
    ! and the derived time control variables echam_phy_tc(:) on all grids
    ! and for all controled processes.

    CALL  eval_echam_phy_config
    CALL  eval_echam_phy_tc
    CALL print_echam_phy_config

    ! For radiation:

    lany=.FALSE.
    DO jg = 1,n_dom
       lany = lany .OR. (echam_phy_tc(jg)%dt_rad > dt_zero)
    END DO
    IF (lany) THEN
      !
      CALL  eval_echam_rad_config
      CALL print_echam_rad_config
      !
      CALL setup_psrad
      CALL setup_cloud_optics
      !
    END IF

    ! For surface processes:
    ! nsfc_type, iwtr, etc. are set in this subroutine.
    ! See mo_sfc_indices.f90 for further details.

    CALL init_sfc_indices( ctest_name )

    ! For turbulent mixing:
    ! Allocate memory for the tri-diagonal solver needed by the implicit
    ! time stepping scheme; Compute time-independent parameters.
    !
    lany=.FALSE.
    DO jg = 1,n_dom
       lany = lany .OR. (echam_phy_tc(jg)%dt_vdf > dt_zero)
    END DO
    IF (lany) THEN
      !
      CALL  eval_echam_vdf_config
      CALL print_echam_vdf_config
      !
      CALL init_vdiff_params( nlev, nlev+1, nlev+1, vct )
      !
      ! Currently the tracer indices are sorted such that we count
      ! the water substances first, and then other species like
      ! aerosols and their precursors. "ntracer" is the total number
      ! of tracers (including water substances) handled in the model;
      ! "iqt" is the starting index for non-water species.
      ! Before more sophisticated meta-data structure becomes available,
      ! it is assumed here that all tracers are subject to turbulent mixing.
      !
      khydromet = iqt - 2        ! # of hydrometeors
      ktrac = ntracer - iqt + 1  ! # of non-water species
      !
      CALL init_vdiff_solver( khydromet, ktrac, nlev )
      !
    ENDIF

    ! For cumulus convection:
    !
    lany=.FALSE.
    DO jg = 1,n_dom
       lany = lany .OR. (echam_phy_tc(jg)%dt_cnv > dt_zero)
    END DO
    IF (lany) THEN
      CALL alloc_echam_cnv_config 
      CALL  eval_echam_cnv_config
      CALL print_echam_cnv_config
    END IF

    lany=.FALSE.
    DO jg = 1,n_dom
       lany = lany .OR. (echam_phy_tc(jg)%dt_cld > dt_zero)
    END DO
    IF (lany) THEN
      CALL print_echam_cld_config
    END IF

    ! Lookup tables for saturation vapour pressure
    !
    lany=.FALSE.
    DO jg = 1,n_dom
       lany = lany .OR. (echam_phy_tc(jg)%dt_cnv > dt_zero)
       lany = lany .OR. (echam_phy_tc(jg)%dt_cld > dt_zero)
       lany = lany .OR. (echam_phy_tc(jg)%dt_vdf > dt_zero)
    END DO
    IF (lany) THEN
       CALL init_convect_tables
       CALL init_echam_convect_tables 
    END IF

    ! For atmospheric gravity wave drag
    !
    lany=.FALSE.
    DO jg = 1,n_dom
       lany = lany .OR. (echam_phy_tc(jg)%dt_gwd > dt_zero)
    END DO
    IF (lany) THEN
       CALL  eval_echam_gwd_config
       CALL print_echam_gwd_config
    END IF

    ! For subgrid scale orographic effects
    !
    lany=.FALSE.
    DO jg = 1,n_dom
       lany = lany .OR. (echam_phy_tc(jg)%dt_sso > dt_zero)
    END DO
    IF (lany) THEN
       CALL  eval_echam_sso_config
       CALL print_echam_sso_config
    END IF


    !-------------------------------------------------------------------
    ! Allocate memory for the state vectors "prm_field" and "prm_tend"
    !-------------------------------------------------------------------
    CALL construct_echam_phy_state( ntracer, p_patch )

    DO jg= 1,n_dom

      IF (ilnd <= nsfc_type) THEN

        ! read time-constant boundary conditions from files

        ! land, glacier and lake masks
        stream_id = openInputFile(land_frac_fn, p_patch(jg), default_read_method)
        CALL read_2D(stream_id=stream_id, location=on_cells,&
             &          variable_name='land',               &
             &          fill_array=prm_field(jg)%lsmask(:,:))
        CALL read_2D(stream_id=stream_id, location=on_cells, &
             &          variable_name='glac',               &
             &          fill_array=prm_field(jg)% glac(:,:))
        CALL read_2D(stream_id=stream_id, location=on_cells, &
             &          variable_name='lake',               &
             &          fill_array=prm_field(jg)% alake(:,:))
        CALL closeFile(stream_id)
        !
        ! At this point, %lsmask is the fraction of land (incl. glacier, but not lakes) in the grid box.
        ! If running without lakes, add lake mask to %lsmask to remove lakes and set %alake to zero.
        IF (.NOT. echam_phy_config(jg)%llake) THEN
          prm_field(jg)%lsmask(:,:) = prm_field(jg)%lsmask(:,:) + prm_field(jg)%alake(:,:)
          prm_field(jg)%alake (:,:) = 0._wp
        END IF

        ! roughness length and background albedo
        stream_id = openInputFile(land_phys_fn, p_patch(jg), default_read_method)

        IF (echam_phy_tc(jg)%dt_vdf > dt_zero) THEN
          CALL read_2D(stream_id=stream_id, location=on_cells, &
                &       variable_name='roughness_length',      &
                &       fill_array=prm_field(jg)% z0m(:,:))
        END IF

        IF (echam_phy_tc(jg)%dt_rad > dt_zero) THEN
          CALL read_2D(stream_id=stream_id, location=on_cells, &
               &       variable_name='albedo',                &
               &       fill_array=prm_field(jg)% alb(:,:))
        END IF
     
        CALL closeFile(stream_id)

        ! orography
        IF (echam_phy_tc(jg)%dt_sso > dt_zero) THEN
          stream_id = openInputFile(land_sso_fn, p_patch(jg), default_read_method)
          CALL read_2D(stream_id=stream_id, location=on_cells, &
               &       variable_name='oromea',                &
               &       fill_array=prm_field(jg)% oromea(:,:))
          CALL read_2D(stream_id=stream_id, location=on_cells, &
             &         variable_name='orostd',                &
             &         fill_array=prm_field(jg)% orostd(:,:))
          CALL read_2D(stream_id=stream_id, location=on_cells, &
             &         variable_name='orosig',                &
             &         fill_array=prm_field(jg)% orosig(:,:))
          CALL read_2D(stream_id=stream_id, location=on_cells, &
             &         variable_name='orogam',                &
             &         fill_array=prm_field(jg)% orogam(:,:))
          CALL read_2D(stream_id=stream_id, location=on_cells, &
             &         variable_name='orothe',                &
             &         fill_array=prm_field(jg)% orothe(:,:))
          CALL read_2D(stream_id=stream_id, location=on_cells, &
             &         variable_name='oropic',                &
             &         fill_array=prm_field(jg)% oropic(:,:))
          CALL read_2D(stream_id=stream_id, location=on_cells, &
             &         variable_name='oroval',                &
             &         fill_array=prm_field(jg)% oroval(:,:))
          CALL closeFile(stream_id)
        END IF

      ELSE

        prm_field(jg)%lsmask(:,:) = 0._wp
        prm_field(jg)%glac  (:,:) = 0._wp
        prm_field(jg)%alake (:,:) = 0._wp

      END IF ! (ilnd <= nsfc_type)

    END DO ! jg

    ! read time-dependent boundary conditions from file

    ! well mixed greenhouse gases, horizontally constant
    IF (ANY(ighg(:) > 0)) THEN
      ! read annual means
      IF (.NOT. bc_greenhouse_gases_file_read) THEN
        CALL read_bc_greenhouse_gases
      END IF
      ! interpolate to the current date and time, placing the annual means at
      ! the mid points of the current and preceding or following year, if the
      ! current date is in the 1st or 2nd half of the year, respectively.
      CALL bc_greenhouse_gases_time_interpolation(mtime_current)
      !
      ! IF a CO2 tracer exists, then copy the time interpolated scalar ghg_co2mmr
      ! to the 3-dimensional tracer field.
      IF ( iqt <= ico2 .AND. ico2 <= ntracer .AND. .NOT. lart) THEN
        DO jg = 1,n_dom
          prm_field(jg)%qtrc(:,:,:,ico2) = ghg_co2mmr
        END DO
      END IF
      !
    ENDIF

    ! interpolation weights for linear interpolation
    ! of monthly means onto the actual integration time step
    current_time_interpolation_weights = calculate_time_interpolation_weights(mtime_current)

!    IF (.NOT. ctest_name(1:3) == 'TPE') THEN

    IF (iice <= nsfc_type .AND. iwtr > nsfc_type) THEN
      CALL finish('','ice tile and no wtr tile not supported yet!')
    END IF
    IF (iice > nsfc_type .AND. iwtr > nsfc_type .AND. ctest_name(1:3) /= 'TPE') THEN
      CALL finish('','only lnd tile present: must use TPE* testcase!')
    END IF

    ! construct stream for radiative forcing diagnostics
    CALL construct_psrad_forcing_list ( p_patch )

    ! read data for simple plumes of aerosols

    IF (ANY(irad_aero(:) == 18)) THEN
      CALL setup_bc_aeropt_splumes
    END IF

    DO jg= 1,n_dom

      ! Read AMIP SST and SIC data
      ! Note: For coupled runs, this is only used for initialization of surface temperatures
      IF (echam_phy_config(jg)%lamip .OR.                   &
          (is_coupled_run() .AND. .NOT. ltestcase) ) THEN
        !
        ! sea surface temperature, sea ice concentration and depth
        CALL read_bc_sst_sic(mtime_current%date%year, p_patch(1))
        !
        CALL bc_sst_sic_time_interpolation(current_time_interpolation_weights, &
             &                             prm_field(jg)%lsmask(:,:) + prm_field(jg)%alake(:,:) > 1._wp - 10._wp*EPSILON(1._wp), &
             &                             prm_field(jg)%ts_tile(:,:,iwtr)   , &
             &                             prm_field(jg)%seaice(:,:)         , &
             &                             prm_field(jg)%siced(:,:)          , &
             &                             p_patch(1)                        )
        !

      ELSE

        prm_field(jg)%seaice(:,:) = 0._wp

      END IF

    END DO ! jg

    lany=.FALSE.
    DO jg = 1,n_dom
       lany = lany .OR. (echam_phy_tc(jg)%dt_car > dt_zero)
    END DO
    IF (lany) THEN
      IF(io3 > ntracer) THEN
        CALL finish('init_echam_phy: mo_echam_phy_init.f90', &
                   &'cannot find an ozone tracer - abort')
      END IF
      IF(n_dom > 1) THEN
        CALL finish('init_echam_phy: mo_echam_phy_init.f90', &
                   &'Cariolle initialization not ready for n_dom>1')
      END IF
      CALL lcariolle_init(                                     &
         & openInputFile_wrap,       closeFile_wrap,           &
         & read_bcast_real_3d_wrap,  read_bcast_real_1d_wrap,  &
         & get_constants                                       )
    END IF

    ! water vapour production by methane oxidation
    lany=.FALSE.
    DO jg = 1,n_dom
      lany = lany .OR. (echam_phy_tc(jg)%dt_mox > dt_zero)
    END DO
    IF (lany) THEN
      CALL init_methox
    END IF
   
#ifndef __NO_JSBACH__
    IF (ilnd <= nsfc_type .AND. ANY(echam_phy_config(:)%ljsb)) THEN

      ! Do basic initialization of JSBACH
      CALL jsbach_init_base(master_namelist_filename)

      ! Now continue initialization of JSBACH for the different grids
      DO jg=1,n_dom
        IF (echam_phy_config(jg)%ljsb) THEN 
          CALL jsbach_init_model( jg, p_patch(jg)) !< in
        END IF
      END DO ! jg

    END IF ! 
#endif


    IF (timers_level > 1) CALL timer_stop(timer_prep_echam_phy)

  END SUBROUTINE init_echam_phy


  !-------------
  !>
  !! Loop over all grid levels and give proper values to some components
  !! of the state vectors "prm_field" and "prm_tend".
  !! This subroutine plays a role similar to "init_g3" in ECHAM6.
  !!
  !! @par Revision History
  !! Initial version by Hui Wan, MPI-M (2010-07)
  !!
  SUBROUTINE initcond_echam_phy( p_patch        ,&
    &                            z_ifc          ,&
    &                            z_mc           ,&
    &                            ddqz_z_full    ,&
    &                            geopot_agl_ifc ,&
    &                            geopot_agl     ,&
    &                            temp           ,&
    &                            qv             ,&
    &                            ctest_name      )

    TYPE(t_patch)    ,INTENT(in) :: p_patch
    REAL(wp)         ,INTENT(in) :: z_ifc         (:,:,:)
    REAL(wp)         ,INTENT(in) :: z_mc          (:,:,:)
    REAL(wp)         ,INTENT(in) :: ddqz_z_full   (:,:,:)
    REAL(wp)         ,INTENT(in) :: geopot_agl_ifc(:,:,:)
    REAL(wp)         ,INTENT(in) :: geopot_agl    (:,:,:)
    REAL(wp)         ,INTENT(in) :: temp          (:,:,:)
    REAL(wp)         ,INTENT(in) :: qv            (:,:,:)
    CHARACTER(LEN=*), INTENT(in) :: ctest_name

    ! local variables and pointers

    INTEGER  :: jg, nblks_c, jb, jbs, jc, jcs, jce
    REAL(wp) :: zlat

    TYPE(t_echam_phy_field),POINTER :: field => NULL()
    TYPE(t_echam_phy_tend) ,POINTER :: tend  => NULL()
    !----

      jg = p_patch%id
    
      field => prm_field(jg)
      tend  => prm_tend (jg)

      nblks_c = p_patch%nblks_c
      jbs     = p_patch%cells%start_blk(2,1)

      ! Assign initial values for some components of the "field" and
      ! "tend" state vectors.

!$OMP PARALLEL
!$OMP WORKSHARE
      !
      ! constant-in-time fields
      field%      clon(:,  :) = p_patch% cells% center(:,:)% lon
      field%      clat(:,  :) = p_patch% cells% center(:,:)% lat
      field% areacella(:,  :) = p_patch% cells%   area(:,:)
      field%    coriol(:,  :) = p_patch% cells%    f_c(:,:)
      !
      field%        zh(:,:,:) =          z_ifc(:,:,:)
      field%        zf(:,:,:) =           z_mc(:,:,:)
      field%        dz(:,:,:) =    ddqz_z_full(:,:,:)
      !
      field%      geoi(:,:,:) = geopot_agl_ifc(:,:,:)
      field%      geom(:,:,:) =     geopot_agl(:,:,:)
      !
      ! initial conditions
      field% qtrc (:,:,:,iqv) = qv(:,:,:)
      field% xvar (:,:,:)     = qv(:,:,:)*0.1_wp
 
      field% swflxsfc_tile(:,:,:) = 0._wp
      field% lwflxsfc_tile(:,:,:) = 0._wp

!$OMP END WORKSHARE
!$OMP END PARALLEL

!$OMP PARALLEL WORKSHARE
        field% ocu   (:,:)   = 0._wp
        field% ocv   (:,:)   = 0._wp
!$OMP END PARALLEL WORKSHARE

      IF (echam_phy_tc(jg)%dt_vdf > dt_zero) THEN
        IF (iwtr<=nsfc_type) field% z0m_tile(:,:,iwtr) = 1e-3_wp !see init_surf in echam (or z0m_oce?)
        IF (iice<=nsfc_type) field% z0m_tile(:,:,iice) = 1e-3_wp !see init_surf in echam (or z0m_ice?)
        IF (ilnd<=nsfc_type) THEN
          field% z0m_tile(:,:,ilnd) = field%z0m(:,:) ! or maybe a larger value?
          field% z0h_lnd(:,:)       = field%z0m(:,:) ! or maybe a larger value?
        END IF
      ENDIF

      ! Initialize some variables for water, ice and land tiles
      ! This can be overridden by the testcases below

      IF (iwtr <= nsfc_type) THEN
        prm_field(jg)% albvisdir_tile(:,:,iwtr) = albedoW ! albedo in the visible range for direct radiation
        prm_field(jg)% albnirdir_tile(:,:,iwtr) = albedoW ! albedo in the NIR range for direct radiation
        prm_field(jg)% albvisdif_tile(:,:,iwtr) = albedoW ! albedo in the visible range for diffuse radiation
        prm_field(jg)% albnirdif_tile(:,:,iwtr) = albedoW ! albedo in the NIR range for diffuse radiation
        prm_field(jg)% albedo_tile   (:,:,iwtr) = albedoW
      END IF

      IF (ilnd <= nsfc_type) THEN

        IF (echam_phy_config(jg)%lamip .OR. (is_coupled_run() .AND. .NOT. ltestcase)) THEN
          prm_field(jg)%ts_tile(:,:,ilnd) = prm_field(jg)%ts_tile(:,:,iwtr)
        END IF

        prm_field(jg)% albvisdir_tile(:,:,ilnd) = prm_field(jg)%alb(:,:)    ! albedo in the visible range for direct radiation
        prm_field(jg)% albnirdir_tile(:,:,ilnd) = prm_field(jg)%alb(:,:)    ! albedo in the NIR range for direct radiation
        prm_field(jg)% albvisdif_tile(:,:,ilnd) = prm_field(jg)%alb(:,:)    ! albedo in the visible range for diffuse radiation
        prm_field(jg)% albnirdif_tile(:,:,ilnd) = prm_field(jg)%alb(:,:)    ! albedo in the NIR range for diffuse radiation
        prm_field(jg)% albedo_tile   (:,:,ilnd) = prm_field(jg)%alb(:,:)

      END IF

      IF (iice <= nsfc_type) THEN

        prm_field(jg)%ts_tile(:,:,iice) = prm_field(jg)%ts_tile(:,:,iwtr)
        !
        prm_field(jg)% albvisdir_tile(:,:,iice) = albi    ! albedo in the visible range for direct radiation
        prm_field(jg)% albnirdir_tile(:,:,iice) = albi    ! albedo in the NIR range for direct radiation
        prm_field(jg)% albvisdif_tile(:,:,iice) = albi    ! albedo in the visible range for diffuse radiation
        prm_field(jg)% albnirdif_tile(:,:,iice) = albi    ! albedo in the NIR range for diffuse radiation
        prm_field(jg)% albedo_tile   (:,:,iice) = albi
        !
        ! The ice model should be able to handle different thickness classes,
        ! but for AMIP we ONLY USE one ice class.
        prm_field(jg)% albvisdir_ice(:,:,:) = albi ! albedo in the visible range for direct radiation
        prm_field(jg)% albnirdir_ice(:,:,:) = albi ! albedo in the NIR range for direct radiation
        prm_field(jg)% albvisdif_ice(:,:,:) = albi ! albedo in the visible range for diffuse radiation
        prm_field(jg)% albnirdif_ice(:,:,:) = albi ! albedo in the NIR range for diffuse radiation
        prm_field(jg)% Tsurf(:,:,:) = Tf
        prm_field(jg)% T1   (:,:,:) = Tf
        prm_field(jg)% T2   (:,:,:) = Tf
        WHERE (prm_field(jg)%seaice(:,:) > 0.0_wp)
           prm_field(jg)% hs   (:,1,:) = 0.1_wp       ! set initial snow depth on sea ice
        ELSEWHERE
           prm_field(jg)% hs   (:,1,:) = 0.0_wp
        ENDWHERE
        prm_field(jg)% hi   (:,1,:) = prm_field(jg)%siced(:,:)
        prm_field(jg)% conc (:,1,:) = prm_field(jg)%seaice(:,:)

      END IF

      ! For idealized test cases

      SELECT CASE (ctest_name)
      CASE('APE','APE_echam','RCEhydro','RCE_glb') !Note that there is only one surface type in this case

!$OMP PARALLEL DO PRIVATE(jb,jc,jcs,jce,zlat) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = jbs,nblks_c
          CALL get_indices_c( p_patch, jb,jbs,nblks_c, jcs,jce, 2)
          DO jc = jcs,jce
            zlat = p_patch%cells%center(jc,jb)%lat
            field% ts_tile(jc,jb,iwtr) = ape_sst(ape_sst_case,zlat)
          END DO
          field% lsmask(jcs:jce,jb) = 0._wp   ! zero land fraction
          field% alake (jcs:jce,jb) = 0._wp   ! zero lake fraction
          field% glac  (jcs:jce,jb) = 0._wp   ! zero glacier fraction
          field% seaice(jcs:jce,jb) = 0._wp   ! zeor sea ice fraction
        END DO
!$OMP END PARALLEL DO

        IF ( is_coupled_run() ) CALL finish('ERROR: Use testcase APEc or APEc_nh for a coupled run')

      CASE('RCE') !Note that there is only one surface type in this case

!$OMP PARALLEL DO PRIVATE(jb,jc,jcs,jce,zlat) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = jbs,nblks_c
          CALL get_indices_c( p_patch, jb,jbs,nblks_c, jcs,jce, 2)
          DO jc = jcs,jce
            zlat = p_patch%cells%center(jc,jb)%lat
            field% ts_tile(jc,jb,iwtr) = th_cbl(1)
          END DO
          field% lsmask(jcs:jce,jb) = 0._wp   ! zero land fraction
          field% alake (jcs:jce,jb) = 0._wp   ! zero lake fraction
          field% glac  (jcs:jce,jb) = 0._wp   ! zero glacier fraction
          field% seaice(jcs:jce,jb) = 0._wp   ! zeor sea ice fraction
        END DO
!$OMP END PARALLEL DO

      CASE('APEi')
        ! The same as APE, except that whenever SST reaches tmelt, we put
        ! 1m-thick ice with a concentration of 0.9 on top

!$OMP PARALLEL DO PRIVATE(jb,jc,jcs,jce,zlat) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = jbs,nblks_c
          CALL get_indices_c( p_patch, jb,jbs,nblks_c, jcs,jce, 2)
          DO jc = jcs,jce
            zlat = p_patch%cells%center(jc,jb)%lat
            ! SST must reach Tf where there's ice. It may be better to modify ape_sst it self.
            field% ts_tile    (jc,jb,iwtr) = ape_sst(ape_sst_case,zlat) + Tf
            ! Initialise the ice - Tsurf, T1 & T2 must be in degC
            field% ts_tile    (jc,jb,iice) = Tf + tmelt
            field% Tsurf      (jc,1, jb  ) = Tf
            field% T1         (jc,1, jb  ) = Tf
            field% T2         (jc,1, jb  ) = Tf
            field% hs         (jc,1, jb  ) = 0._wp
            IF ( field%ts_tile(jc,jb,iwtr) <= Tf + tmelt ) THEN
              field%Tsurf (jc,1,jb) = field% ts_tile(jc,jb,iice) - tmelt
              field%conc  (jc,1,jb) = 0.9_wp
              field%hi    (jc,1,jb) = 1.0_wp
              field%seaice(jc,  jb) = field%conc(jc,1,jb)
            ELSE
              field%conc  (jc,1,jb) = 0._wp
              field%hi    (jc,1,jb) = 0._wp
              field%seaice(jc,  jb) = field%conc(jc,1,jb)
            ENDIF
          END DO
          field% lsmask(jcs:jce,jb) = 0._wp   ! zero land fraction
          field% alake (jcs:jce,jb) = 0._wp   ! zero lake fraction
          field% glac  (jcs:jce,jb) = 0._wp   ! zero glacier fraction
        END DO
!$OMP END PARALLEL DO
        field% albvisdir_ice(:,:,:) = albi    ! albedo in the visible range for direct radiation
        field% albnirdir_ice(:,:,:) = albi    ! albedo in the NIR range for direct radiation
        field% albvisdif_ice(:,:,:) = albi    ! albedo in the visible range for diffuse radiation
        field% albnirdif_ice(:,:,:) = albi    ! albedo in the NIR range for diffuse radiation

      CASE('APEc','APEc_nh')
        ! The same as APEi, except we initialize with no ice and don't modify the surface
        ! temperature. This is meant for a coupled run.

!$OMP PARALLEL DO PRIVATE(jb,jc,jcs,jce,zlat) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = jbs,nblks_c
          CALL get_indices_c( p_patch, jb,jbs,nblks_c, jcs,jce, 2)
          DO jc = jcs,jce
            zlat = p_patch%cells%center(jc,jb)%lat
            field% ts_tile(jc,jb,iwtr) = ape_sst(ape_sst_case,zlat)
            ! Initialise the ice - Tsurf, T1 & T2 must be in degC
            field% ts_tile    (jc,jb,iice) = Tf + tmelt
            field% Tsurf      (jc,1, jb  ) = Tf
            field% T1         (jc,1, jb  ) = Tf
            field% T2         (jc,1, jb  ) = Tf
            field% hs         (jc,1, jb  ) = 0._wp
            field%conc  (jc,1,jb) = 0._wp
            field%hi    (jc,1,jb) = 0._wp
            field%seaice(jc,  jb) = field%conc(jc,1,jb)
          END DO
          field% lsmask(jcs:jce,jb) = 0._wp   ! zero land fraction
          field% alake (jcs:jce,jb) = 0._wp   ! zero lake fraction
          field% glac  (jcs:jce,jb) = 0._wp   ! zero glacier fraction
        END DO
!$OMP END PARALLEL DO
        field% albvisdir_ice(:,:,:) = albi    ! albedo in the visible range for direct radiation
        field% albnirdir_ice(:,:,:) = albi    ! albedo in the NIR range for direct radiation
        field% albvisdif_ice(:,:,:) = albi    ! albedo in the visible range for diffuse radiation
        field% albnirdif_ice(:,:,:) = albi    ! albedo in the NIR range for diffuse radiation

      CASE('TPEc', 'TPEo') !Note that there is only one surface type (ilnd) in this case

!$OMP PARALLEL DO PRIVATE(jb,jc,jcs,jce,zlat) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = jbs,nblks_c
          CALL get_indices_c( p_patch, jb,jbs,nblks_c, jcs,jce, 2)
          field% lsmask(jcs:jce,jb) = 1._wp   ! land fraction = 1
          field% alake (jcs:jce,jb) = 0._wp   ! zero lake fraction
          field% glac  (jcs:jce,jb) = 0._wp   ! zero glacier fraction
          field% seaice(jcs:jce,jb) = 0._wp   ! zeor sea ice fraction

          field% ts_tile(jcs:jce,jb,ilnd) = tpe_temp
        END DO
!$OMP END PARALLEL DO

      CASE('JWw-Moist','LDF-Moist','jabw_m')

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jcs,jce) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = jbs,nblks_c
          CALL get_indices_c( p_patch, jb,jbs,nblks_c, jcs,jce, 2)

          ! Set the surface temperature to the same value as the lowest model
          ! level above surface. For this test case, currently we assume
          ! there is no land or sea ice.

          field% ts_tile(jcs:jce,jb,iwtr) = temp(jcs:jce,nlev,jb)

          field% lsmask(jcs:jce,jb) = 0._wp   ! zero land fraction
          field% alake (jcs:jce,jb) = 0._wp   ! zero lake fraction
          field% glac  (jcs:jce,jb) = 0._wp   ! zero glacier fraction
          field% seaice(jcs:jce,jb) = 0._wp   ! zero sea ice fraction
        END DO
!$OMP END DO  NOWAIT
!$OMP END PARALLEL

      END SELECT

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jcs,jce) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = jbs,nblks_c
        CALL get_indices_c( p_patch, jb,jbs,nblks_c, jcs,jce, 2)

        ! Initialize the flag lfland (.TRUE. if the fraction of land+lake in
        ! a grid box is larger than zero). In ECHAM a local array
        ! is initialized in each call of the subroutine "physc"
        DO jc = jcs,jce
          field%lfland(jc,jb) = (field%lsmask(jc,jb) + field%alake(jc,jb)) > 0._wp
          field%lfglac(jc,jb) = field%glac  (jc,jb) > 0._wp
        END DO

      END DO      !jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      ! Settings for total surface
      ! (after tile masks and variables potentially have been overwritten by testcases above)

      IF (iwtr <= nsfc_type) THEN
        prm_field(jg)%ts       (:,:) = prm_field(jg)%ts_tile(:,:,iwtr)
        prm_field(jg)%albvisdir(:,:) = albedoW
        prm_field(jg)%albvisdif(:,:) = albedoW
        prm_field(jg)%albnirdir(:,:) = albedoW
        prm_field(jg)%albnirdif(:,:) = albedoW
        prm_field(jg)%albedo   (:,:) = albedoW
      ELSE
        prm_field(jg)%ts       (:,:) = prm_field(jg)%ts_tile(:,:,ilnd)
        prm_field(jg)%albvisdir(:,:) = prm_field(jg)%alb(:,:)
        prm_field(jg)%albvisdif(:,:) = prm_field(jg)%alb(:,:)
        prm_field(jg)%albnirdir(:,:) = prm_field(jg)%alb(:,:)
        prm_field(jg)%albnirdif(:,:) = prm_field(jg)%alb(:,:)
        prm_field(jg)%albedo   (:,:) = prm_field(jg)%alb(:,:)
      END IF

      prm_field(jg)%ts_rad     (:,:) = prm_field(jg)%ts(:,:)
      prm_field(jg)%ts_rad_rt  (:,:) = prm_field(jg)%ts(:,:)

      NULLIFY( field,tend )

  END SUBROUTINE initcond_echam_phy
  !-------------
  !>
  !!
  SUBROUTINE additional_restart_init( p_patch, ctest_name )

    TYPE(t_patch),   INTENT(IN) :: p_patch(:)
    CHARACTER(LEN=*),INTENT(IN) :: ctest_name

    INTEGER :: nblks_c, jg, jb, jbs, jc, jcs, jce
    REAL(wp):: zlat

    TYPE(t_echam_phy_field),POINTER :: field => NULL()

!!$    CHARACTER(LEN=*),PARAMETER :: routine = 'additional_restart_init'

    !-------------------------
    ! Loop over all domains
    !-------------------------
    DO jg = 1,n_dom

      field => prm_field(jg)

      nblks_c = p_patch(jg)%nblks_c
      jbs     = p_patch(jg)%cells%start_blk(2,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jcs,jce,zlat) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = jbs,nblks_c
        CALL get_indices_c( p_patch(jg), jb,jbs,nblks_c, jcs,jce, 2)

        !---------------------------------------------------------------------
        ! Re-initialize SST, sea ice and glacier if necessary
        !---------------------------------------------------------------------
          SELECT CASE (ctest_name)
          CASE('APE','APE_echam','RCEhydro')
          ! For an aqua-planet experiment, re-initialization is necessary if
          ! the restart file in use was generated during a differently configured
          ! experiment (e.g., an APE exp with a different SST setup, or
          ! a real-world simulation such as AMIP, etc).

            DO jc = jcs,jce
              zlat = p_patch(jg)%cells%center(jc,jb)%lat
              field% ts_tile(jc,jb,iwtr) = ape_sst(ape_sst_case,zlat)   ! SST
              field% ts     (jc,     jb) = field% ts_tile(jc,jb,iwtr)
            END DO
            field% lsmask(jcs:jce,jb) = 0._wp   ! zero land fraction
            field% alake (jcs:jce,jb) = 0._wp   ! zero lake fraction
            field% glac  (jcs:jce,jb) = 0._wp   ! zero glacier fraction
            field% seaice(jcs:jce,jb) = 0._wp   ! zeor sea ice fraction

          END SELECT

        !--------------------------------------------------------------------
        ! Initialize the flag lfland (.TRUE. if the fraction of land+lake in
        ! a grid box is larger than zero). In ECHAM a local array
        ! is initialized in each call of the subroutine "physc".
        ! Note that this initialization is needed for all resumed integrations
        ! regardless of the choice of "ctest_name", because
        ! logical variables can not yet be stored in restart files.
        !--------------------------------------------------------------------

        DO jc = jcs,jce
          field%lfland(jc,jb) = (field%lsmask(jc,jb) + field%alake(jc,jb)) > 0._wp
          field%lfglac(jc,jb) = field%glac  (jc,jb) > 0._wp
        ENDDO !jc
      ENDDO   !jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      NULLIFY( field )
    ENDDO !jg

  END SUBROUTINE additional_restart_init

END MODULE mo_echam_phy_init
