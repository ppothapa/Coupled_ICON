
MODULE mo_construct_icon_hamocc


  USE mo_kind, ONLY           : wp
  USE mo_hamocc_nml, ONLY     : l_init_bgc
  USE mo_exception, ONLY      : message
  USE mo_control_bgc, ONLY    : dtb, dtbgc, inv_dtbgc, ndtdaybgc,  &
       &                        ldtrunbgc, bgc_zlevs, bgc_nproma, &
       &                        inv_dtb
!#ifdef AVFLUX
!  USE mo_avflux, ONLY         : avflux_ini
!#endif


  USE mo_memory_bgc, ONLY      : alloc_mem_biomod, n90depth,n1000depth,n2000depth,&
  &                              alloc_mem_carbch,totalarea
  USE mo_param1_bgc, ONLY      : set_tracer_indices
  USE mo_bgc_icon_comm, ONLY  : ini_bgc_regions, initial_update_icon, hamocc_state, &
      &                         print_bgc_parameters,print_wpoc, update_bgc

  USE mo_sedmnt, ONLY         : alloc_mem_sedmnt, ini_bottom, sediment_bottom
  USE mo_ini_bgc, ONLY        : ini_aquatic_tracers,            &
       &                        ini_pore_water_tracers,         &
       &                        ini_atmospheric_concentrations, &
       &                        set_parameters_bgc, &
       &                        ini_continental_carbon_input, &
       &                        ini_wpoc, bgc_param_conv_unit

  USE mo_model_domain,        ONLY: t_patch, t_patch_3D
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_run_config,          ONLY: dtime
  USE mo_impl_constants,      ONLY: max_char_length
  USE mo_ocean_types,         ONLY: t_hydro_ocean_state
  USE mo_ocean_nml,           ONLY: n_zlev, dzlev_m
  USE mo_dynamics_config,     ONLY: nold,nnew
  USE mo_parallel_config,     ONLY: nproma
  USE mo_sync,                ONLY: global_sum_array
  USE mo_math_utilities,      ONLY: set_zlev
  USE mo_ocean_surface_types, ONLY: t_atmos_for_ocean 
  USE mo_end_bgc,             ONLY: cleanup_hamocc
  USE mo_bgc_bcond,           ONLY: construct_bgc_ext_data, destruct_bgc_ext_data, ext_data_bgc
  USE mo_timer,               ONLY: timer_start, timer_stop, timer_bgc_ini, ltimer
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_hamocc_output,       ONLY: construct_hamocc_var_lists, construct_hamocc_state, &
    &                                destruct_hamocc_state         
  USE mo_master_config,       ONLY: isRestart
  USE mo_ocean_hamocc_couple_state, ONLY: t_ocean_to_hamocc_state, t_hamocc_to_ocean_state, &
    & t_hamocc_ocean_state

  IMPLICIT NONE
  PRIVATE
  
  
  PUBLIC:: construct_icon_hamocc, destruct_icon_hamocc, init_icon_hamocc 
CONTAINS

!------------------------------------------------------------
SUBROUTINE destruct_icon_hamocc()

    CALL cleanup_hamocc
    CALL destruct_bgc_ext_data
    CALL destruct_hamocc_state(hamocc_state)

END SUBROUTINE destruct_icon_hamocc
!------------------------------------------------------------

!------------------------------------------------------------
SUBROUTINE construct_icon_hamocc(patch_3d, ext_data)
   TYPE(t_patch_3D),             TARGET,INTENT(IN)    :: patch_3D
   TYPE(t_external_data), TARGET, INTENT(inout) :: ext_data
   
    CALL construct_hamocc_var_lists(patch_3d%p_patch_2d(1))
    CALL construct_hamocc_state(patch_3d, hamocc_state)
    CALL construct_bgc_ext_data(patch_3d%p_patch_2d(1), ext_data,ext_data_bgc)

END SUBROUTINE construct_icon_hamocc
!------------------------------------------------------------

!------------------------------------------------------------
SUBROUTINE init_icon_hamocc(hamocc_ocean_state)
  TYPE(t_hamocc_ocean_state), TARGET   :: hamocc_ocean_state

  if(ltimer)call timer_start(timer_bgc_ini)
  
  CALL ini_bgc_icon(hamocc_ocean_state, isRestart())
  if(ltimer)call timer_stop(timer_bgc_ini)
   
END SUBROUTINE init_icon_hamocc
!------------------------------------------------------------
 
!------------------------------------------------------------
!------------------------------------------------------------
SUBROUTINE INI_BGC_ICON(hamocc_ocean_state,l_is_restart)
  TYPE(t_hamocc_ocean_state), TARGET   :: hamocc_ocean_state
 ! USE mo_hamocc_diagnostics,  ONLY: get_inventories

  !! Arguments
  TYPE(t_patch_3D), POINTER    :: p_patch_3D
  LOGICAL, INTENT(in):: l_is_restart


  !! Local variables

  CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'ini_bgc_icon'
  TYPE(t_patch),POINTER    :: p_patch
  TYPE(t_subset_range), POINTER :: all_cells

  INTEGER :: jc,  jb
  INTEGER :: start_index, end_index
  INTEGER  :: levels(nproma)

  INTEGER, POINTER              :: regions(:,:)

  REAL(wp),ALLOCATABLE:: dlevels_m(:),dlevels_i(:)
 
  CALL message(TRIM(routine), 'start')
  p_patch_3D => hamocc_ocean_state%patch_3D
  !
  !----------------------------------------------------------------------
  !
  ! Set control constants ( mo_control_bgc )
  !
  dtbgc = dtime                         !  time step length [sec]
  inv_dtbgc = 1.0_wp / dtbgc
  ndtdaybgc = NINT(86400._wp / dtbgc) !  time steps per day [no.]
  dtb = 1._wp / REAL(ndtdaybgc, wp)   !  time step length [days]
  inv_dtb = 1._wp/ dtb


   ! determine size of arrays 
   p_patch           => p_patch_3D%p_patch_2D(1)

   bgc_zlevs = n_zlev
   bgc_nproma = nproma

   all_cells => p_patch%cells%ALL

   regions => p_patch_3D%regio_c

  !
  ! Initialize time step counter of run.
  !
  ldtrunbgc = 0

  
  CALL message(TRIM(routine), 'set_parameters_bgc' )

  CALL set_tracer_indices
  CALL set_parameters_bgc




  ! region indices
  CALL message(TRIM(routine), 'ini_bgc_regions' )

  CALL INI_BGC_REGIONS


  !
  ! Allocate memory : biology
  !
  CALL message(TRIM(routine), 'alloc_mem_biomod')
  CALL ALLOC_MEM_BIOMOD

  !
  ! Allocate memory : sediment
  !
  CALL message(TRIM(routine), 'alloc_mem_sedmnt' )
  CALL ALLOC_MEM_SEDMNT

  !
  ! Allocate memory : inorganic carbon cycle
  !
  CALL message(TRIM(routine), 'alloc_mem_carbch')
  CALL ALLOC_MEM_CARBCH

  !

  !
  ! Initialize sediment layering
  !
  CALL message(TRIM(routine), 'sediment_bottom')
  CALL sediment_bottom
  
  ! convert 1/d to 1/ts
  CALL bgc_param_conv_unit

  CALL print_bgc_parameters 

  CALL message(TRIM(routine), 'bgc_param_conv_unit')

  CALL message(TRIM(routine), 'ini weathering fluxes')
  totalarea = 0._wp
  DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, start_index, end_index)
        DO jc=start_index, end_index
           totalarea = totalarea + p_patch%cells%area(jc,jb) * &
 &             p_patch_3d%wet_halo_zero_c(jc,1,jb)

!       write(0,*) "cell area:", p_patch%cells%area(jc,jb), " wet_halo_zero_c:", p_patch_3d%wet_halo_zero_c(jc,1,jb)
        ENDDO
  ENDDO
  totalarea     = global_sum_array(totalarea)
  CALL ini_continental_carbon_input(totalarea)

! set level for 90m, 1000m, 2000m for diagnostic output
  ALLOCATE(dlevels_m(n_zlev))
  ALLOCATE(dlevels_i(n_zlev+1))

  CALL set_zlev(dlevels_i, dlevels_m, n_zlev, dzlev_m)

  n90depth   = maxloc(dlevels_i,DIM=1,MASK=(dlevels_i < 90._wp))
  n1000depth = maxloc(dlevels_i,DIM=1,MASK=(dlevels_i < 1000._wp))
  n2000depth = maxloc(dlevels_i,DIM=1,MASK=(dlevels_i < 2000._wp))

! Initialize POC sinking speed
  CALL ini_wpoc(dlevels_i)
  CALL print_wpoc

  DEALLOCATE(dlevels_m)
  DEALLOCATE(dlevels_i)
   
!DIR$ INLINE
 DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, start_index, end_index)
        !  tracer 1: potential temperature
        !  tracer 2: salinity
        
        levels(start_index:end_index) = p_patch_3d%p_patch_1d(1)%dolic_c(start_index:end_index,jb)


        CALL ini_bottom(start_index,end_index,levels,p_patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,jb))

        IF(l_is_restart.and..not.l_init_bgc)then

          CALL update_bgc(start_index,end_index,levels,&
             & p_patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,jb),&  ! cell thickness
             &jb, hamocc_state%p_prog(nnew(1))%tracer(:,:,jb,:)&
             &,hamocc_ocean_state%ocean_to_hamocc_state%co2_mixing_ratio(:,jb)&
             & ,hamocc_state%p_diag,hamocc_state%p_sed, hamocc_state%p_tend)

         ELSE

          CALL ini_aquatic_tracers(start_index,end_index,levels,regions(:,jb))

          CALL ini_pore_water_tracers(start_index,end_index)

         ENDIF 

        CALL ini_atmospheric_concentrations

        CALL initial_update_icon(start_index,end_index,levels, &
              &               p_patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,jb),&  ! cell thickness
              &              jb, &
              &              hamocc_state%p_prog(nold(1))%tracer(:,:,jb,:),&
              &              hamocc_state%p_sed, hamocc_state%p_diag,&
              &              hamocc_ocean_state%hamocc_to_ocean_state%co2_flux(:,jb))
     

  ENDDO


 CALL message(TRIM(routine), 'end ini bgc')
END SUBROUTINE INI_BGC_ICON


END MODULE mo_construct_icon_hamocc
!------------------------------------------------------------
