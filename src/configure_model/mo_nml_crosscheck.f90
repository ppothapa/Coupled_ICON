!>
!! This module checks the read-in namelist parameters and, in case of
!! inconsistencies, it tries to correct these.
!!
!!
!! @author Kristina Froehlich, MPI-M (2011-07-12)
!! @author Hui Wan, MPI-M (2011-07-12)
!!
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_nml_crosscheck

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message, message_text, finish
  USE mo_impl_constants,     ONLY: ildf_echam, inwp, iheldsuarez,                    &
    &                              ildf_dry, inoforcing, ihs_atm_temp,               &
    &                              ihs_atm_theta, tracer_only, inh_atmosphere,       &
    &                              ishallow_water, LEAPFROG_EXPL, LEAPFROG_SI,       &
    &                              NO_HADV, UP, MIURA, MIURA3, FFSL, FFSL_HYB,       &
    &                              MCYCL, MIURA_MCYCL, MIURA3_MCYCL,                 &
    &                              FFSL_MCYCL, FFSL_HYB_MCYCL, iecham,               &
    &                              RAYLEIGH_CLASSIC,                                 &
    &                              iedmf, icosmo, MODE_IAU, MODE_IAU_OLD
  USE mo_time_config,        ONLY: time_config, dt_restart
  USE mo_extpar_config,      ONLY: itopo                                             
  USE mo_io_config,          ONLY: dt_checkpoint, lflux_avg,inextra_2d, inextra_3d,  &
    &                              lnetcdf_flt64_output
  USE mo_parallel_config,    ONLY: check_parallel_configuration,                     &
    &                              num_io_procs, itype_comm,                         &
    &                              num_prefetch_proc, use_dp_mpi2io
  USE mo_limarea_config,     ONLY: latbc_config, LATBC_TYPE_CONST
  USE mo_master_config,      ONLY: isRestart
  USE mo_run_config,         ONLY: nsteps, dtime, iforcing,                          &
    &                              ltransport, ntracer, nlev, ltestcase,             &
    &                              nqtendphy, iqtke, iqv, iqc, iqi,                  &
    &                              iqs, iqr, iqt, iqtvar, ltimer,             &
    &                              ico2, ich4, in2o, io3,                     &
    &                              iqni, iqni_nuc, iqg, iqm_max,                     &
    &                              iqh, iqnr, iqns, iqng, iqnh, iqnc,                & 
    &                              inccn, ininact, ininpot,                          &
    &                              activate_sync_timers, timers_level, lart
  USE mo_dynamics_config,    ONLY: iequations, lshallow_water, ltwotime
  USE mo_advection_config,   ONLY: advection_config

  USE mo_nonhydrostatic_config, ONLY: itime_scheme_nh => itime_scheme,               &
    &                                 lhdiff_rcf, rayleigh_type,                     &
    &                                 ivctype, ndyn_substeps
  USE mo_ha_dyn_config,      ONLY: ha_dyn_config
  USE mo_diffusion_config,   ONLY: diffusion_config
  USE mo_atm_phy_nwp_config, ONLY: atm_phy_nwp_config, icpl_aero_conv, iprog_aero
  USE mo_lnd_nwp_config,     ONLY: ntiles_lnd, lsnowtile
  USE mo_echam_phy_config,   ONLY: echam_phy_config
  USE mo_radiation_config,   ONLY: irad_o3, irad_aero
  USE mo_turbdiff_config,    ONLY: turbdiff_config
  USE mo_initicon_config,    ONLY: init_mode, dt_iau, ltile_coldstart, timeshift
  USE mo_nh_testcases_nml,   ONLY: nh_test_name
  USE mo_ha_testcases,       ONLY: ctest_name, ape_sst_case

  USE mo_meteogram_config,   ONLY: check_meteogram_configuration
  USE mo_grid_config,        ONLY: lplane, n_dom, l_limited_area, start_time,        &
    &                              nroot, is_plane_torus

  USE mo_art_config,         ONLY: art_config
  USE mo_time_management,    ONLY: compute_timestep_settings,                        &
    &                              compute_restart_settings,                         &
    &                              compute_date_settings
  USE mo_event_manager,      ONLY: initEventManager
  USE mtime,                 ONLY: getTotalMilliSecondsTimeDelta, datetime,          &
    &                              newDatetime, deallocateDatetime
  USE mo_gridref_config,     ONLY: grf_intmethod_e
  USE mo_interpol_config
  USE mo_sleve_config
  USE mo_nudging_nml,        ONLY: check_nudging

#ifdef __ICON_ART
  USE mo_grid_config,        ONLY: lredgrid_phys
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: atm_crosscheck

  CHARACTER(LEN = *), PARAMETER :: modname = "mo_nml_crosscheck"


CONTAINS

  SUBROUTINE atm_crosscheck

    INTEGER  :: jg
    INTEGER  :: jt   ! tracer loop index
    INTEGER  :: z_go_tri(11)  ! for crosscheck
    CHARACTER(len=*), PARAMETER :: routine =  modname//'::atm_crosscheck'
    REAL(wp) :: restart_time
    TYPE(datetime), POINTER :: reference_dt
    
    !--------------------------------------------------------------------
    ! Compute date/time/time step settings
    ! and initialize the event manager
    !--------------------------------------------------------------------
    !
    ! Note that the ordering of the following three calls must not be
    ! changed, since they rely on previous results:
    !
    CALL compute_timestep_settings()
    CALL compute_restart_settings()
    CALL compute_date_settings("atm", dt_restart, nsteps)
    !
    ! Create an event manager, ie. a collection of different events
    !
    CALL initEventManager(time_config%tc_exp_refdate)

    !--------------------------------------------------------------------
    ! Parallelization
    !--------------------------------------------------------------------
    CALL check_parallel_configuration()


    !--------------------------------------------------------------------
    ! Limited Area Mode and LatBC read-in:
    !--------------------------------------------------------------------


    IF (lplane) CALL finish(routine,&
      'Currently a plane version is not available')

    ! Reset num_prefetch_proc to zero if the model does not run in limited-area mode
    ! or if there are no lateral boundary data to be read
    IF (.NOT. l_limited_area .OR. latbc_config%itype_latbc == 0) num_prefetch_proc = 0

    ! If LatBC data is unavailable: Idle-wait-and-retry
    !
    IF (latbc_config%nretries > 0) THEN
      !
      ! ... only supported for prefetching LatBC mode
      IF (.NOT. l_limited_area .OR. (num_prefetch_proc == 0)) THEN
        CALL finish(routine, "LatBC: Idle-wait-and-retry only supported for prefetching LatBC mode!")
      END IF
      !
      ! ... only supported for MPI parallel configuration
#ifdef NOMPI
      CALL finish(routine, "LatBC: Idle-wait-and-retry requires MPI!")
#endif
    END IF

    ! Limited area mode must not be enabled for torus grid:
    IF (is_plane_torus .AND. l_limited_area) THEN
      CALL finish(routine, 'Plane torus grid requires l_limited_area = .FALSE.!')
    END IF

    ! Root bisection "0" does not make sense for limited area mode; it
    ! is more likely that the user tried to use a torus grid here:
    IF (l_limited_area .AND. (nroot == 0)) THEN
      CALL finish(routine, "Root bisection 0 does not make sense for limited area mode; did you try to use a torus grid?")
    END IF


    !--------------------------------------------------------------------
    ! Grid and dynamics
    !--------------------------------------------------------------------

    IF (lplane) CALL finish( routine,&
      'Currently a plane version is not available')

    SELECT CASE (iequations)
    CASE(IHS_ATM_TEMP,IHS_ATM_THETA)         ! hydrostatic atm model

      SELECT CASE(iforcing)
      CASE(INOFORCING,IHELDSUAREZ,ILDF_DRY)  ! without moist processes
        ha_dyn_config%ldry_dycore = .TRUE.
      CASE(IECHAM,ILDF_ECHAM)                ! with ECHAM physics
        CALL finish(routine, 'Hydrostatic dynamics cannot be used with ECHAM physics')
      END SELECT

    END SELECT

    lshallow_water = (iequations==ISHALLOW_WATER)

    SELECT CASE (iequations)
    CASE (IHS_ATM_TEMP,IHS_ATM_THETA,ISHALLOW_WATER)

      ltwotime = (ha_dyn_config%itime_scheme/=LEAPFROG_EXPL).AND. &
                 (ha_dyn_config%itime_scheme/=LEAPFROG_SI)

    END SELECT


    ! Limited area mode must not be enabled for torus grid:
    IF (is_plane_torus .AND. l_limited_area) THEN
      CALL finish(routine, 'Plane torus grid requires l_limited_area = .FALSE.!')
    END IF

    ! Root bisection "0" does not make sense for limited area mode; it
    ! is more likely that the user tried to use a torus grid here:
    IF (l_limited_area .AND. (nroot == 0)) THEN
      CALL finish(routine, "Root bisection 0 does not make sense for limited area mode; did you try to use a torus grid?")
    END IF

    !--------------------------------------------------------------------
    ! If ltestcase is set to .FALSE. in run_nml set testcase name to empty
    ! (in case it is still set in the run script)
    IF (.NOT. ltestcase) THEN
      ctest_name = ''
      nh_test_name = ''
    END IF
    !--------------------------------------------------------------------

    !--------------------------------------------------------------------
    ! Testcases (hydrostatic)
    !--------------------------------------------------------------------


    IF ((ctest_name=='GW') .AND. (nlev /= 20)) THEN
      CALL finish(routine,'nlev MUST be 20 for the gravity-wave test case')
    ENDIF

    IF ((ctest_name=='SV') .AND. ntracer /= 2 ) THEN

      CALL finish(routine, &
        & 'ntracer MUST be 2 for the stationary vortex test case')
    ENDIF

    IF ((ctest_name=='DF1') .AND. ntracer == 1 ) THEN
      CALL finish(routine, &
        & 'ntracer MUST be >=2 for the deformational flow test case 1')
    ENDIF

    IF ((ctest_name=='DF2') .AND. ntracer == 1 ) THEN
      CALL finish(routine, &
        & 'ntracer MUST be >=2 for the deformational flow test case 2')
    ENDIF

    IF ((ctest_name=='DF3') .AND. ntracer == 1 ) THEN
      CALL finish(routine, &
        & 'ntracer MUST be >=2 for the deformational flow test case 3')
    ENDIF

    IF ((ctest_name=='DF4') .AND. ntracer == 1 ) THEN
      CALL finish(routine, &
        & 'ntracer MUST be >=2 for the deformational flow test case 4')
    ENDIF

    IF ((ctest_name=='APE') .AND. (ape_sst_case=='sst_ice')  ) THEN
      IF (.NOT. lflux_avg)&
      CALL finish(routine, &
        & 'lflux_avg must be set true to run this setup')
    ENDIF

    !--------------------------------------------------------------------
    ! Testcases (nonhydrostatic)
    !--------------------------------------------------------------------
    IF (.NOT. ltestcase .AND. rayleigh_type == RAYLEIGH_CLASSIC) THEN
      CALL finish(routine, &
        & 'rayleigh_type = RAYLEIGH_CLASSIC not applicable to real case runs.')
    ENDIF

    IF ( ( nh_test_name=='APE_nwp'.OR. nh_test_name=='dcmip_tc_52' ) .AND.  &
      &  ( ANY(atm_phy_nwp_config(:)%inwp_surface == 1 ) ) .AND.                       &
      &  ( ANY(atm_phy_nwp_config(:)%inwp_turb    /= iedmf ) ) ) THEN
      CALL finish(routine, &
        & 'surface scheme must be switched off, when running the APE test')
    ENDIF

    !--------------------------------------------------------------------
    ! Shallow water
    !--------------------------------------------------------------------
    IF (iequations==ISHALLOW_WATER.AND.ha_dyn_config%lsi_3d) THEN
      CALL message( routine, 'lsi_3d = .TRUE. not applicable to shallow water model')
    ENDIF

    IF ((iequations==ISHALLOW_WATER).AND.(nlev/=1)) &
    CALL finish(routine,'Multiple vertical level specified for shallow water model')

    !--------------------------------------------------------------------
    ! Hydrostatic atm
    !--------------------------------------------------------------------
    IF (iequations==IHS_ATM_THETA) ha_dyn_config%ltheta_dyn = .TRUE.

    !--------------------------------------------------------------------
    ! Nonhydrostatic atm
    !--------------------------------------------------------------------
    IF (lhdiff_rcf .AND. (itype_comm == 3)) CALL finish(routine, &
      'lhdiff_rcf is available only for idiv_method=1 and itype_comm<=2')

    IF (grf_intmethod_e >= 5 .AND. iequations /= INWP .AND. n_dom > 1) THEN
      grf_intmethod_e = 4
      CALL message( routine, 'grf_intmethod_e has been reset to 4')
    ENDIF

    !--------------------------------------------------------------------
    ! Atmospheric physics, general
    !--------------------------------------------------------------------
    IF ((iforcing==INWP).AND.(iequations/=INH_ATMOSPHERE)) &
    CALL finish( routine, 'NWP physics only implemented in the '//&
               'nonhydrostatic atm model')

    !--------------------------------------------------------------------
    ! NWP physics
    !--------------------------------------------------------------------
    IF (iforcing==inwp) THEN

      DO jg =1,n_dom

        IF (atm_phy_nwp_config(1)%inwp_turb /= iedmf) THEN
          IF( atm_phy_nwp_config(jg)%inwp_satad == 0       .AND. &
          & ((atm_phy_nwp_config(jg)%inwp_convection >0 ) .OR. &
          &  (atm_phy_nwp_config(jg)%inwp_gscp > 0      )    ) ) &
          &  CALL finish( routine,'satad has to be switched on')
        ENDIF

        IF( (atm_phy_nwp_config(jg)%inwp_gscp==0) .AND. &
          & (atm_phy_nwp_config(jg)%inwp_convection==0) .AND.&
          & (atm_phy_nwp_config(jg)%inwp_radiation==0) .AND.&
          & (atm_phy_nwp_config(jg)%inwp_sso==0)  .AND. &
          & (atm_phy_nwp_config(jg)%inwp_surface == 0) .AND.&
          & (atm_phy_nwp_config(jg)%inwp_turb> 0) )   &
        CALL message(routine,' WARNING! NWP forcing set but '//&
                    'only turbulence selected!')


        IF (( atm_phy_nwp_config(jg)%inwp_turb == icosmo ) .AND. &
          & (turbdiff_config(jg)%lconst_z0) ) THEN
          CALL message(routine,' WARNING! NWP forcing set but '//  &
                      'idealized (horizontally homogeneous) roughness '//&
                      'length z0 selected!')
        ENDIF

        IF (.NOT. ltestcase .AND. atm_phy_nwp_config(jg)%inwp_surface == 0) THEN
          CALL finish( routine,'Real-data applications require using a surface scheme!')
        ENDIF

        ! check radiation scheme in relation to chosen ozone and irad_aero=6 to itopo

        IF ( (atm_phy_nwp_config(jg)%inwp_radiation > 0) )  THEN

          SELECT CASE (irad_o3)
          CASE (0) ! ok
            CALL message(routine,'radiation is used without ozone')
          CASE (2,4,6,7,8,9,79,97) ! ok
            CALL message(routine,'radiation is used with ozone')
          CASE (10) ! ok
            CALL message(routine,'radiation is used with ozone calculated from ART')
            IF ( .NOT. lart ) THEN
              CALL finish(routine,'irad_o3 currently is 10 but lart is false.')
            ENDIF
          CASE default
            CALL finish(routine,'irad_o3 currently has to be 0, 2, 4, 6, 7, 8 or 9.')
          END SELECT

          ! Tegen aerosol and itopo (Tegen aerosol data have to be read from external data file)
          IF ( ( irad_aero == 6 ) .AND. ( itopo /=1 ) ) THEN
            CALL finish(routine,'irad_aero=6 requires itopo=1')
          ENDIF

          IF ( ( irad_aero /= 6 .AND. irad_aero /= 9 ) .AND.  &
            &  ( atm_phy_nwp_config(jg)%icpl_aero_gscp > 0 .OR. icpl_aero_conv > 0 ) ) THEN
            CALL finish(routine,'aerosol-precipitation coupling requires irad_aero=6 or =9')
          ENDIF
        ELSE

          SELECT CASE (irad_o3)
          CASE(0) ! ok
          CASE default
            irad_o3 = 0
            CALL message(routine,'running without radiation => irad_o3 reset to 0')
          END SELECT

        ENDIF !inwp_radiation

        !! check microphysics scheme
        IF (  atm_phy_nwp_config(jg)%mu_rain < 0.0   .OR. &
          &   atm_phy_nwp_config(jg)%mu_rain > 5.0)  THEN
          CALL finish(routine,'mu_rain requires: 0 < mu_rain < 5')
        END IF

        IF (  atm_phy_nwp_config(jg)%mu_snow < 0.0   .OR. &
          &   atm_phy_nwp_config(jg)%mu_snow > 5.0)  THEN
          CALL finish(routine,'mu_snow requires: 0 < mu_snow < 5')
        END IF ! microphysics

        IF (atm_phy_nwp_config(jg)%inwp_surface == 0 .AND. ntiles_lnd > 1) THEN
          ntiles_lnd = 1
          CALL message(routine,'Warning: ntiles reset to 1 because the surface scheme is turned off')
        ENDIF

      ENDDO
    END IF


    !--------------------------------------------------------------------
    ! Tracers and diabatic forcing
    !--------------------------------------------------------------------

    !ATTENTION: note that in the following we make use of jg without and jg-loop.
    !           Curently it is right for the wrong reasons.

    !
    ! Check settings of ntracer
    !
    ! Set tracer indices
    !
    SELECT CASE(iforcing)
    CASE (IECHAM,ILDF_ECHAM)  ! iforcing

      IF (ntracer < 3) CALL finish(routine,'ECHAM physics needs at least 3 tracers')

      ! 0 indicates that this tracer is not (yet) used by ECHAM  physics
      iqv    = 1     !> water vapour
      iqc    = 2     !! cloud water
      iqi    = 3     !! ice
      iqr    = 0     !! rain water
      iqs    = 0     !! snow
      iqm_max= 3     !! end index of water species mixing ratios
      iqt    = 4     !! starting index of non-water species
      io3    = 4     !! O3
      ico2   = 5     !! CO2
      ich4   = 6     !! CH4
      in2o   = 7     !! N2O
      nqtendphy = 0  !! number of water species for which convective and turbulent
                     !! tendencies are stored

    CASE (INWP) ! iforcing

      ! ** NWP physics section ** 
      !
      ! IMPORTANT: For NWP physics, five microphysics tracers (QV, QC, QI, QR and QS) must always be
      !            defined because their presence is implicitly assumed in the physics-dynamics interface
      !            when converting between temperature and virtual temperature.
      !
      !            Any additional mass-related microphysics tracers must be numbered in consecutive order
      !            after iqs = 5. The parameter "iqm_max" must signify the highest tracer index carrying a
      !            moisture mixing ratio. Additional tracers for hydrometeor number concentrations (in
      !            case of a two-moment scheme) or other purposes (aerosols, qt variance or anything else)
      !            can be numbered freely after iqm_max. The parameter "iqt", denoting the start index
      !            of tracers not related at all to moisture, is used in configure_advection to specify
      !            the index range of tracers for which advection is turned off in the stratosphere
      !            (i.e. all cloud and precipitation variables including number concentrations)
      !
      !            Note also that the namelist parameter "ntracer" is reset automatically to the correct
      !            value when NWP physics is used in order to avoid multiple namelist changes when playing
      !            around with different physics schemes.
      !
      ico2      = 0     !> co2, 0: not to be used with NWP physics
      !
      ! Default settings valid for all microphysics options
      !
      iqv       = 1     !> water vapour
      iqc       = 2     !! cloud water
      iqi       = 3     !! ice 
      iqr       = 4     !! rain water
      iqs       = 5     !! snow
      nqtendphy = 3     !! number of water species for which convective and turbulent tendencies are stored
      !
      ! The following parameters may be reset depending on the selected physics scheme
      !
      iqm_max   = 5     !! end index of water species mixing ratios
      iqt       = 6     !! start index of other tracers not related at all to moisture
      !
      ntracer   = 5     !! total number of tracers
      !
      ! dummy settings
      iqni     = ntracer+100    !! cloud ice number
      iqni_nuc = ntracer+100    !! activated ice nuclei  
      iqg      = ntracer+100    !! graupel
      iqtvar   = ntracer+100    !! qt variance (for EDMF turbulence)
      iqh      = ntracer+100
      iqnr     = ntracer+100  
      iqns     = ntracer+100
      iqng     = ntracer+100
      iqnh     = ntracer+100
      iqnc     = ntracer+100
      inccn    = ntracer+100
      ininpot  = ntracer+100
      ininact  = ntracer+100

      !
      ! Taking the 'highest' microphysics option in some cases allows using more
      ! complex microphysics schemes in nested domains than in the global domain
      ! However, a clean implementation would require a domain-dependent 'ntracer' dimension
      SELECT CASE (MAXVAL(atm_phy_nwp_config(1:n_dom)%inwp_gscp))
        

      CASE(2)  ! COSMO-DE (3-cat ice: snow, cloud ice, graupel)

       ! CALL finish('mo_atm_nml_crosscheck', 'Graupel scheme not implemented.')
        
        iqg     = 6       !! graupel
        iqm_max = iqg
        iqt     = iqt + 1

        ntracer = ntracer + 1  !! increase total number of tracers by 1

 
      CASE(3)  ! improved ice nucleation scheme C. Koehler (note: iqm_max does not change!)

        iqni     = 6     !! cloud ice number
        iqni_nuc = 7     !! activated ice nuclei  
        iqt     = iqt + 2

        ntracer = ntracer + 2  !! increase total number of tracers by 2

      CASE(4)  ! two-moment scheme 
      
        iqg  = 6
        iqh  = 7
        iqni = 8        
        iqnr = 9        
        iqns = 10        
        iqng = 11        
        iqnh = 12
        iqnc = 13
        ininact = 14

        nqtendphy = 3     !! number of water species for which convective and turbulent tendencies are stored
        iqm_max   = 7     !! end index of water species mixing ratios
        iqt       = 15    !! start index of other tracers not related at all to moisture
       
        ntracer = 14

      CASE(5)  ! two-moment scheme with CCN and IN budgets
      
        iqg  = 6
        iqh  = 7
        iqni = 8        
        iqnr = 9        
        iqns = 10        
        iqng = 11        
        iqnh = 12
        iqnc = 13
        ininact = 14
        inccn   = 15
        ininpot = 16

        nqtendphy = 3     !! number of water species for which convective and turbulent tendencies are stored
        iqm_max   = 7     !! end index of water species mixing ratios
        iqt       = 17    !! start index of other tracers not related at all to moisture
       
        ntracer = 16
        
      CASE(6)
      
        iqg  = 6
        iqh  = 7
        iqni = 8        
        iqnr = 9        
        iqns = 10        
        iqng = 11        
        iqnh = 12
        iqnc = 13
        ininact = 14
        
        nqtendphy = 3     !! number of water species for which convective and turbulent tendencies are stored
        iqm_max   = 7     !! end index of water species mixing ratios
        iqt       = 14    !! start index of other tracers not related at all to moisture
        
        ntracer = 14
        
      END SELECT ! microphysics schemes


      IF (atm_phy_nwp_config(1)%inwp_turb == iedmf) THEN ! EDMF turbulence

        iqtvar = iqt       !! qt variance
        iqt    = iqt + 1   !! start index of other tracers than hydrometeors

        ntracer = ntracer + 1  !! increase total number of tracers by 1

        DO jg =1,n_dom
          turbdiff_config(jg)%ldiff_qi = .TRUE.  !! turbulent diffusion of QI on (by EDMF)
        ENDDO

!dmk    ntiles_lnd = 5     !! EDMF currently only works with 5 land tiles - consistent with TESSEL
                           !! even if the land model is inactive ntiles_lnd should be 5

      ENDIF

      IF ( (advection_config(1)%iadv_tke) > 0 ) THEN
        IF ( atm_phy_nwp_config(1)%inwp_turb == icosmo ) THEN
          iqtke = iqt        !! TKE
 
          ! Note that iqt is not increased, since TKE does not belong to the hydrometeor group.

          ntracer = ntracer + 1  !! increase total number of tracers by 1

          WRITE(message_text,'(a,i3)') 'Attention: TKE is advected, '//&
                                       'ntracer is increased by 1 to ',ntracer
          CALL message(routine,message_text)
        ELSE
          WRITE(message_text,'(a,i2)') 'TKE advection not supported for inwp_turb= ', &
            &                          atm_phy_nwp_config(1)%inwp_turb

          CALL finish(routine, message_text )
        ENDIF
      ENDIF

      ! Note: Indices for additional tracers are assigned automatically
      ! via add_tracer_ref in mo_nonhydro_state.

      WRITE(message_text,'(a,i3)') 'Attention: NWP physics is used, '//&
                                   'ntracer is automatically reset to ',ntracer
      CALL message(routine,message_text)


      IF (lart) THEN
        
        ntracer = ntracer + art_config(1)%iart_ntracer
        
        WRITE(message_text,'(a,i3,a,i3)') 'Attention: transport of ART tracers is active, '//&
                                     'ntracer is increased by ',art_config(1)%iart_ntracer, &
                                     ' to ',ntracer
        CALL message(routine,message_text)

      ENDIF

      ! set the nclass_gscp variable for land-surface scheme to number of hydrometeor mixing ratios
      DO jg = 1, n_dom
        atm_phy_nwp_config(jg)%nclass_gscp = iqm_max
      ENDDO

    CASE default ! iforcing

        iqv    = 1     !> water vapour
        iqc    = 2     !! cloud water
        iqi    = 3     !! ice
        iqr    = 0     !! 0: no rain water
        iqs    = 0     !! 0: no snow
        ico2   = 4     !! CO2
        iqm_max= 3     !! end index of water species mixing ratios
        iqt    = 4     !! starting index of non-water species
        nqtendphy = 0  !! number of water species for which convective and turbulent
                       !! tendencies are stored

    END SELECT ! iforcing

    ! take into account additional passive tracers, if present
    ! no need to update iqt, since passive tracers do not belong to the hydrometeor group.
    ! ATTENTION: as long as ntracer is not domain specific, we set jg=1
    IF ( advection_config(1)%npassive_tracer > 0) THEN
      ntracer = ntracer + advection_config(1)%npassive_tracer
      WRITE(message_text,'(a,i3,a,i3)') 'Attention: passive tracers have been added, '//&
                                   'ntracer is increased by ',advection_config(1)%npassive_tracer, &
                                   ' to ',ntracer
      CALL message(routine,message_text)
    ENDIF


    IF (ltransport) THEN

      DO jg = 1,n_dom

        SELECT CASE ( iforcing )
        CASE ( INWP )
        !...........................................................
        ! in NWP physics
        !...........................................................


          ! Force settings for tracer iqtke, if TKE advection is performed
          !
          IF ( advection_config(jg)%iadv_tke > 0 ) THEN

            ! force monotonous slope limiter for vertical advection
            advection_config(jg)%itype_vlimit(iqtke) = 2

            ! force positive definite flux limiter for horizontal advection
            advection_config(jg)%itype_hlimit(iqtke) = 4

            SELECT CASE (advection_config(jg)%iadv_tke)
            CASE (1)
              ! switch off horizontal advection
              advection_config(jg)%ihadv_tracer(iqtke) = 0
            CASE (2)
              ! check whether horizontal substepping is switched on 
              IF (ALL( (/22,32,42,52/) /= advection_config(jg)%ihadv_tracer(iqtke)) ) THEN
                ! choose Miura with substepping
                advection_config(jg)%ihadv_tracer(iqtke) = 22
              ENDIF
            END SELECT

          ENDIF

        END SELECT ! iforcing

      END DO ! jg = 1,n_dom
    END IF ! ltransport


    !--------------------------------------------------------------------
    ! Tracer transport
    !--------------------------------------------------------------------
    ! General

    SELECT CASE (iequations)
    CASE (INH_ATMOSPHERE)

      IF ((itime_scheme_nh==tracer_only) .AND. (.NOT.ltransport)) THEN
        WRITE(message_text,'(A,i2,A)') &
          'nonhydrostatic_nml:itime_scheme set to ', tracer_only, &
          '(TRACER_ONLY), but ltransport to .FALSE.'
        CALL finish( routine,message_text)
      END IF

    CASE (IHS_ATM_TEMP,IHS_ATM_THETA,ISHALLOW_WATER)

      IF ( (ha_dyn_config%itime_scheme==tracer_only).AND. &
           (.NOT.ltransport)) THEN
        WRITE(message_text,'(A,i2,A)') &
          'ha_dyn_nml:itime_scheme set to ', tracer_only, &
          '(TRACER_ONLY), but ltransport to .FALSE.'
        CALL finish( routine,message_text)
      END IF

    END SELECT

    IF (ltransport) THEN
    DO jg = 1,n_dom


      !----------------------------------------------
      ! Flux computation methods - consistency check

        z_go_tri(1:11)=(/NO_HADV,UP,MIURA,MIURA3,FFSL,FFSL_HYB,MCYCL,       &
          &              MIURA_MCYCL,MIURA3_MCYCL,FFSL_MCYCL,FFSL_HYB_MCYCL/)
        DO jt=1,ntracer
          IF ( ALL(z_go_tri /= advection_config(jg)%ihadv_tracer(jt)) ) THEN
            CALL finish( routine,                                       &
              &  'incorrect settings for TRI-C grid ihadv_tracer. Must be '// &
              &  '0,1,2,3,4,5,6,20,22,32,42 or 52 ')
          ENDIF
        ENDDO

    END DO ! jg = 1,n_dom
    END IF ! ltransport


    !--------------------------------------------------------------------
    ! Horizontal diffusion
    !--------------------------------------------------------------------

    DO jg =1,n_dom

      SELECT CASE( diffusion_config(jg)%hdiff_order )
      CASE(-1)
        WRITE(message_text,'(a,i2.2)') 'Horizontal diffusion '//&
                                       'switched off for domain ', jg
        CALL message(routine,message_text)

      CASE(2,3,4,5)
        CONTINUE

      CASE(24,42)
        IF (.NOT.( iequations==IHS_ATM_TEMP)) CALL finish(routine, &
        ' hdiff_order = 24 or 42 only implemented for the hydrostatic atm model')

      CASE DEFAULT
        CALL finish(routine,                       &
          & 'Error: Invalid choice for  hdiff_order. '// &
          & 'Choose from -1, 2, 3, 4, 5, 24, and 42.')
      END SELECT

      IF ( diffusion_config(jg)%hdiff_efdt_ratio<=0._wp) THEN
        CALL message(routine,'No horizontal background diffusion is used')
      ENDIF

      IF (lshallow_water)  diffusion_config(jg)%lhdiff_temp=.FALSE.

      IF (itype_comm == 3 .AND. diffusion_config(jg)%hdiff_order /= 5)  &
        CALL finish(routine, 'itype_comm=3 requires hdiff_order = 5')

      IF (itype_comm == 3 .AND. (diffusion_config(jg)%itype_vn_diffu > 1 .OR. &
        diffusion_config(jg)%itype_t_diffu > 1) )                             &
        CALL finish(routine, 'itype_comm=3 requires itype_t/vn_diffu = 1')

    ENDDO

    !--------------------------------------------------------------------
    ! checking the meanings of the io settings
    !--------------------------------------------------------------------


    IF (lnetcdf_flt64_output) THEN
       CALL message(routine,'NetCDF output of floating point variables will be in 64-bit accuracy')
       IF (.NOT. use_dp_mpi2io) THEN
          use_dp_mpi2io = .TRUE.
          CALL message(routine,'--> use_dp_mpi2io is changed to .TRUE. to allow 64-bit accuracy in the NetCDF output.')
       END IF
    ELSE
       CALL message(routine,'NetCDF output of floating point variables will be in 32-bit accuracy')
    END IF

    SELECT CASE(iforcing)
    CASE ( iecham, ildf_echam )
      inextra_2d   = 0
      inextra_3d   = 0

    CASE DEFAULT
    END SELECT


    IF (activate_sync_timers .AND. .NOT. ltimer) THEN
      activate_sync_timers = .FALSE.
      WRITE (message_text,*) &
        & "warning: namelist parameter 'activate_sync_timers' has been set to .FALSE., ", &
        & "because global 'ltimer' flag is disabled."
      CALL message(routine, message_text)
    END IF
    IF (timers_level > 9 .AND. .NOT. activate_sync_timers) THEN
      activate_sync_timers = .TRUE.
      WRITE (message_text,*) &
        & "warning: namelist parameter 'activate_sync_timers' has been set to .TRUE., ", &
        & "because global 'timers_level' is > 9."
      CALL message(routine, message_text)
    END IF


    !--------------------------------------------------------------------
    ! Realcase runs
    !--------------------------------------------------------------------

    IF ( ANY((/MODE_IAU,MODE_IAU_OLD/) == init_mode) ) THEN  ! start from dwd analysis with incremental update

      ! check analysis update window
      !
      IF ( (dt_iau > 0._wp) .AND. (dt_iau < dtime)) THEN
        ! If dt_iau is chosen to be larger than 0, it must be >= dtime at least.
        dt_iau = dtime
        WRITE (message_text,'(a,a,f6.2)') "Wrong value for dt_iau. ", &
          &   "If >0 then at least equal to advective/phys tstep ",dtime
        CALL finish('initicon_nml:', message_text)
      ENDIF 

      reference_dt => newDatetime("1980-06-01T00:00:00.000")

      IF (dt_checkpoint > 0._wp) THEN
        restart_time = MIN(dt_checkpoint, &
          &                0.001_wp * getTotalMilliSecondsTimeDelta(time_config%tc_dt_restart, reference_dt))
      ELSE
        restart_time = 0.001_wp * getTotalMilliSecondsTimeDelta(time_config%tc_dt_restart, reference_dt)
      ENDIF

      CALL deallocateDatetime(reference_dt)
      IF (.NOT. isRestart() .AND. (restart_time <= dt_iau+timeshift%dt_shift)) THEN
        WRITE (message_text,'(a)') "Restarting is not allowed within the IAU phase"
        CALL finish('atm_crosscheck:', message_text)
      ENDIF

      DO jg = 2, n_dom
        IF (start_time(jg) > timeshift%dt_shift .AND. start_time(jg) < dt_iau+timeshift%dt_shift) THEN
          WRITE (message_text,'(a)') "Starting a nest is not allowed within the IAU phase"
          CALL finish('atm_crosscheck:', message_text)
        ENDIF
      ENDDO

      ! IAU modes MODE_IAU_OLD cannot be combined with snowtiles
      ! when performing snowtile warmstart.
      IF ((ntiles_lnd > 1) .AND. (.NOT. ltile_coldstart) .AND. (lsnowtile)) THEN
        IF ( init_mode == MODE_IAU_OLD ) THEN
          WRITE (message_text,'(a,i2)') "lsnowtile=.TRUE. not allowed for IAU-Mode ", init_mode   
          CALL finish(routine, message_text)
        ENDIF
      ENDIF

    ENDIF


    ! check meteogram configuration
    CALL check_meteogram_configuration(num_io_procs)

    IF (iforcing==iecham) CALL land_crosscheck()

    CALL art_crosscheck()

    CALL check_nudging( n_dom, iequations, iforcing, ivctype, top_height,                     &
      &                 l_limited_area, latbc_config%lsparse_latbc, latbc_config%itype_latbc, & 
      &                 latbc_config%nudge_hydro_pres, LATBC_TYPE_CONST, is_plane_torus,      &
      &                 lart, ndyn_substeps                                                   )

  END  SUBROUTINE atm_crosscheck
  !---------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------
  SUBROUTINE land_crosscheck

    INTEGER  :: jg
    CHARACTER(len=*), PARAMETER :: routine =  modname//'::land_crosscheck'

#ifdef __NO_JSBACH__
    IF (ANY(echam_phy_config(:)%ljsb)) THEN
      CALL finish(routine, "This version was compiled without jsbach. Compile with __JSBACH__, or set ljsb=.FALSE.")
    ENDIF
#else
    DO jg=1,n_dom
      IF (.NOT.echam_phy_config(jg)%ljsb) THEN
         IF (echam_phy_config(jg)%llake) THEN
            CALL message(routine, 'Setting llake = .FALSE. since ljsb = .FALSE.')
            echam_phy_config(jg)%llake = .FALSE.
         END IF
      END IF
    END DO
#endif

  END SUBROUTINE land_crosscheck
  !---------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------
  SUBROUTINE art_crosscheck
  
    CHARACTER(len=*), PARAMETER :: routine =  modname//'::art_crosscheck'
#ifdef __ICON_ART
    INTEGER  :: &
      &  jg
#endif
    
#ifndef __ICON_ART
    IF (lart) THEN
        CALL finish( routine,'run_nml: lart is set .TRUE. but ICON was compiled without -D__ICON_ART')
    ENDIF
#endif
    
    IF (.NOT. lart .AND. irad_aero == 9 ) THEN
      CALL finish(routine,'irad_aero=9 needs lart = .TRUE.')
    END IF

    IF ( ( irad_aero == 9 ) .AND. ( iprog_aero /= 0 ) ) THEN
      CALL finish(routine,'irad_aero=9 requires iprog_aero=0')
    ENDIF
    
#ifdef __ICON_ART
    IF ( ( irad_aero == 9 ) .AND. ( itopo /=1 ) ) THEN
      CALL finish(routine,'irad_aero=9 requires itopo=1')
    ENDIF
    
    DO jg= 1,n_dom
      IF(lredgrid_phys(jg) .AND. irad_aero == 9) THEN
        CALL finish(routine,'irad_aero=9 does not work with a reduced radiation grid')
      ENDIF
      IF(art_config(jg)%iart_ari == 0 .AND. irad_aero == 9) THEN
        CALL finish(routine,'irad_aero=9 needs iart_ari > 0')
      ENDIF
      IF(art_config(jg)%iart_ari > 0  .AND. irad_aero /= 9) THEN
        CALL finish(routine,'iart_ari > 0 requires irad_aero=9')
      ENDIF
    ENDDO
    
    ! XML specification checks
    
    DO jg= 1,n_dom
      IF(art_config(jg)%lart_aerosol .AND. art_config(jg%cart_aerosol_xml)=='') THEN
        CALL finish(routine,'lart_aerosol=.TRUE. but no cart_aerosol_xml specified')
      ENDIF
      IF(art_config(jg)%lart_chem .AND. art_config(jg%cart_chemistry_xml)=='') THEN
        CALL finish(routine,'lart_chem=.TRUE. but no cart_chemistry_xml specified')
      ENDIF
      IF(art_config(jg)%lart_passive .AND. art_config(jg%cart_passive_xml)=='') THEN
        CALL finish(routine,'lart_passive=.TRUE. but no cart_passive_xml specified')
      ENDIF
    ENDDO
    
#endif
  END SUBROUTINE art_crosscheck
  !---------------------------------------------------------------------------------------
END MODULE mo_nml_crosscheck
