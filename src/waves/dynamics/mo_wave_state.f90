!>
!!        Contains the variables to set up the wave  model.
!=============================================================================================
!!
!! @author Mikhail Dobrynin, DWD, 11.06.2019
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
!!
!=============================================================================================

MODULE mo_wave_state

  USE mo_exception,            ONLY: message, finish, message_text
  USE mo_parallel_config,      ONLY: nproma
  USE mo_model_domain,         ONLY: t_patch
  USE mo_grid_config,          ONLY: n_dom, l_limited_area, ifeedback_type
  USE mo_impl_constants,       ONLY: success, max_char_length, VNAME_LEN, TLEV_NNOW_RCF
  USE mo_var_list,             ONLY: add_var, add_ref, t_var_list_ptr
  USE mo_var_list_register,    ONLY: vlr_add, vlr_del
  USE mo_var_list_register_utils, ONLY: vlr_add_vref
  USE mo_var_groups,           ONLY: groups
  USE mo_cdi_constants,        ONLY: GRID_UNSTRUCTURED_CELL, GRID_CELL
  USE mo_cdi,                  ONLY: DATATYPE_FLT32, DATATYPE_FLT64, GRID_UNSTRUCTURED, &
       &                             DATATYPE_PACK16, DATATYPE_INT
  USE mo_zaxis_type,           ONLY: ZA_SURFACE, ZA_HEIGHT_10M
  USE mo_cf_convention,        ONLY: t_cf_var
  USE mo_grib2,                ONLY: t_grib2_var, grib2_var
  USE mo_io_config,            ONLY: lnetcdf_flt64_output
  USE mo_dynamics_config,      ONLY: nsav1
  USE mo_mpi,                  ONLY: get_my_mpi_work_id
  USE mo_run_config,           ONLY: ntracer

  USE mo_var_metadata,         ONLY: get_timelevel_string
  USE mo_var_metadata_types,   ONLY: t_var_metadata,t_var_metadata_dynamic
  USE mo_advection_config,     ONLY: t_advection_config, advection_config
  USE mo_tracer_metadata,      ONLY: create_tracer_metadata
 
  USE mo_wave_types,           ONLY: t_wave_prog, t_wave_diag, &
       &                             t_wave_state, t_wave_state_lists
  USE mo_wave_config,          ONLY: t_wave_config, wave_config

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_wave_state'

  PUBLIC :: construct_wave_state    ! Constructor for the wave state
  PUBLIC :: destruct_wave_state     ! Destructor

  PUBLIC :: p_wave_state            ! state vector of wave variables
  PUBLIC :: p_wave_state_lists      ! lists for state vector of wave variables

  TYPE(t_wave_state),       TARGET, ALLOCATABLE :: p_wave_state(:)
  TYPE(t_wave_state_lists), TARGET, ALLOCATABLE :: p_wave_state_lists(:)

CONTAINS

  SUBROUTINE construct_wave_state(p_patch, n_timelevels)

    TYPE(t_patch), INTENT(IN) :: p_patch(:)
    INTEGER,       INTENT(IN) :: n_timelevels

    CHARACTER(len=max_char_length) :: listname, varname_prefix
    CHARACTER(len=*), PARAMETER :: routine = modname//'::construct_wave_state'

    INTEGER :: ntl,  &! local number of timelevels
         ist,        &! status
         jg,         &! grid level counter
         jt           ! time level counter
    LOGICAL :: l_extra_timelev

    CALL message (routine, 'Construction of wave state started')

    ALLOCATE (p_wave_state(n_dom),p_wave_state_lists(n_dom), stat=ist)
    IF (ist /= success) THEN
       CALL finish(TRIM(routine),'allocation for wave state failed')
    END IF

    DO jg = 1, n_dom

       ntl = n_timelevels

       ! As grid nesting is not called at every dynamics time step, an extra time
       ! level is needed for full-field interpolation and boundary-tendency calculation
       IF (n_dom > 1) THEN
          ntl = ntl + 1
          nsav1(jg) = ntl
       END IF

       IF (ifeedback_type == 1 .AND. jg > 1 .OR. l_limited_area .AND. jg == 1) ntl = ntl + 1

       ALLOCATE(p_wave_state(jg)%prog(1:ntl), STAT=ist)
       IF (ist/=SUCCESS) CALL finish(routine,                                   &
            'allocation of prognostic state array failed')

       ! create state list
       ALLOCATE(p_wave_state_lists(jg)%prog_list(1:ntl), STAT=ist)
       IF (ist/=SUCCESS) CALL finish(routine,                                   &
            'allocation of prognostic state list array failed')

       ! create tracer list (no extra timelevels)
       ALLOCATE(p_wave_state_lists(jg)%tracer_list(1:n_timelevels), STAT=ist)
       IF (ist/=SUCCESS) CALL finish(routine,                                   &
            'allocation of prognostic tracer list array failed')


       !
       ! Build lists for every timelevel
       !
       DO jt = 1, ntl

          ! Tracer fields do not need extra time levels because feedback is not incremental
          ! and the nest-call frequency is always synchronized with the advection time step
          l_extra_timelev = jt > n_timelevels

          WRITE(listname,'(a,i2.2,a,i2.2)') 'wave_state_prog_of_domain_',jg, &
               &                               '_and_timelev_',jt

          ! Build prog state list
          ! includes memory allocation
          varname_prefix = ''
          CALL new_wave_state_prog_list(p_patch(jg), p_wave_state(jg)%prog(jt), &
               & p_wave_state_lists(jg)%prog_list(jt), &
               & listname, TRIM(varname_prefix), &
               & l_extra_timelev, jt)

          WRITE(listname,'(a,i2.2,a,i2.2)') 'wave_state_tracer_of_domain_',jg, &
               & '_and_timelev_',jt
          varname_prefix = ''
          CALL new_wave_state_tracer_list(p_patch(jg), p_wave_state_lists(jg)%prog_list(jt), &
               p_wave_state_lists(jg)%tracer_list(jt), listname )

       END DO

       ! Build diag state list
       ! includes memory allocation
       WRITE(listname,'(a,i2.2)') 'wave_state_diag_of_domain_',jg
       CALL new_wave_state_diag_list(&
            p_patch(jg), &
            p_wave_state(jg)%diag, &
            p_wave_state_lists(jg)%diag_list, &
            listname)

       p_wave_state(jg)%diag%process_id(:,:) = get_my_mpi_work_id()

    END DO

    CALL message (routine, 'wave state construction completed')

  END SUBROUTINE construct_wave_state


  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Creates tracer var list.
  !!
  !! Creates tracer var list containing references to all prognostic tracer
  !! fields.
  !!
  !! @par Revision History
  !! Initial release by Daniel Reinert, DWD (2012-02-02)
  !!
  SUBROUTINE new_wave_state_tracer_list (p_patch, from_var_list, p_tracer_list, listname)
    TYPE(t_patch), INTENT(IN) :: p_patch ! current patch
    TYPE(t_var_list_ptr), INTENT(IN) :: from_var_list ! source list to be referenced
    TYPE(t_var_list_ptr), INTENT(INOUT) :: p_tracer_list ! new tracer list (containing all tracers)
    CHARACTER(*), INTENT(IN) :: listname
    TYPE (t_var_metadata), POINTER :: from_info
    TYPE (t_var_metadata_dynamic), POINTER :: from_info_dyn
    INTEGER :: iv

    ! Register a field list and apply default settings
    CALL vlr_add(p_tracer_list, TRIM(listname), patch_id=p_patch%id, &
         &       lrestart=.FALSE., loutput =.FALSE.)
    ! add references to all tracer fields of the source list (prognostic state)
    DO iv = 1, from_var_list%p%nvars
       ! retrieve information from actual linked list element
       from_info => from_var_list%p%vl(iv)%p%info
       from_info_dyn => from_var_list%p%vl(iv)%p%info_dyn
       ! Only add tracer fields to the tracer list
       IF (from_info_dyn%tracer%lis_tracer .AND. .NOT.from_info%lcontainer) &
            &   CALL vlr_add_vref(p_tracer_list, from_info%name, from_var_list, in_group=groups())
    END DO
  END SUBROUTINE new_wave_state_tracer_list


  SUBROUTINE new_wave_state_prog_list(p_patch, p_prog, p_prog_list, listname, &
       &                              vname_prefix, l_extra_timelev, timelev)

    TYPE(t_patch), TARGET, INTENT(IN)    :: p_patch
    TYPE(t_wave_prog),     INTENT(INOUT) :: p_prog
    TYPE(t_var_list_ptr),  INTENT(INOUT) :: p_prog_list !< current prognostic state list
    CHARACTER(len=*),      INTENT(IN)    :: listname, vname_prefix
    INTEGER,               INTENT(IN)    :: timelev
    LOGICAL,               INTENT(IN)    :: l_extra_timelev

    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: nblks_c    !< number of cell blocks to allocate

    INTEGER :: ibits         !< "entropy" of horizontal slice
    INTEGER :: datatype_flt  !< floating point accuracy in NetCDF output

    CHARACTER(len=4)         :: suffix
    CHARACTER(len=VNAME_LEN) :: freq_ind_str, dir_ind_str
    CHARACTER(len=14) :: tracer_name

    TYPE(t_advection_config), POINTER :: advconf
    TYPE(t_wave_config),      POINTER :: wc

    INTEGER :: shape3d_c(3), shape4d_c(4)

    INTEGER :: jt, nlev

    !determine size of arrays
    nblks_c = p_patch%nblks_c

    ! pointer to advection_config(jg) to save some paperwork
    advconf => advection_config(p_patch%id)
    ! pointer to wave_config(jg) to save some paperwork
    wc => wave_config(p_patch%id)

    ! number of vertical levels
    nlev   = p_patch%nlev

    shape4d_c    = (/nproma, nlev, nblks_c, ntracer/)
    shape3d_c    = (/nproma, nlev, nblks_c/)

    ibits = DATATYPE_PACK16   ! "entropy" of horizontal slice

    IF (lnetcdf_flt64_output) THEN
       datatype_flt = DATATYPE_FLT64
    ELSE
       datatype_flt = DATATYPE_FLT32
    END IF

    ! Suffix (mandatory for time level dependent variables)
    suffix = get_timelevel_string(timelev)

    !
    ! Register a field list and apply default settings
    !
    CALL vlr_add(p_prog_list, TRIM(listname), patch_id=p_patch%id, lrestart=.TRUE.)

    NULLIFY(p_prog%tracer)

    cf_desc    = t_cf_var('tracer', '', 'spectral bin of wave energy', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_prog_list, 'tracer', p_prog%tracer,                          &
         &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,  &
         &           ldims=shape4d_c ,                                           &
         &           lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE. )

    ALLOCATE(p_prog%tracer_ptr(ntracer))

    DO jt = 1, ntracer
       write(freq_ind_str,'(I0.3)') wc%freq_ind(jt)
       write(dir_ind_str,'(I0.3)') wc%dir_ind(jt)

       tracer_name = TRIM(advconf%tracer_names(jt))//'_'//TRIM(freq_ind_str)//'_'//TRIM(dir_ind_str)

       CALL add_ref( p_prog_list, 'tracer',                                          &
            & TRIM(tracer_name)//suffix, p_prog%tracer_ptr(jt)%p_3d,                 &
            & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                 &
            & t_cf_var(TRIM(tracer_name), '-','spectral bin '//TRIM(tracer_name), &
            & datatype_flt),                                                      &
            & grib2_var(0, 0, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL),            &
            & ref_idx=jt,                                                         &
            & ldims=shape3d_c,                                                    &
            & loutput=.TRUE.,                                                     &
            & tlev_source=TLEV_NNOW_RCF,                                          &
            & tracer_info=create_tracer_metadata(lis_tracer=.TRUE.,               &
            &                       name        = TRIM(tracer_name)//suffix,      &
            &                       lfeedback   = .TRUE.,                         &
            &                       ihadv_tracer=advconf%ihadv_tracer(jt),        &
            &                       ivadv_tracer=advconf%ivadv_tracer(jt)),       &
            & in_group=groups("wave_spectrum"))

    END DO

  END SUBROUTINE new_wave_state_prog_list



  SUBROUTINE new_wave_state_diag_list(p_patch, p_diag, p_diag_list, listname)

    TYPE(t_patch),         INTENT(IN)    :: p_patch
    TYPE(t_wave_diag),     INTENT(INOUT) :: p_diag
    TYPE(t_var_list_ptr),  INTENT(INOUT) :: p_diag_list
    CHARACTER(len=*),      INTENT(IN)    :: listname

    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: ibits         !< "entropy" of horizontal slice
    INTEGER :: datatype_flt  !< floating point accuracy in NetCDF output
    INTEGER :: nblks_c
    INTEGER :: jg
    INTEGER :: shape2d_c(2), shape3d_freq_c(3), shape3d_tr_c(3)


    jg      = p_patch%id
    nblks_c = p_patch%nblks_c

    shape2d_c         = (/nproma, nblks_c/)
    shape3d_freq_c    = (/nproma, nblks_c, wave_config(jg)%nfreqs/)
    shape3d_tr_c      = (/nproma, nblks_c, ntracer/)

    ibits = DATATYPE_PACK16   ! "entropy" of horizontal slice

    IF ( lnetcdf_flt64_output ) THEN
       datatype_flt = DATATYPE_FLT64
    ELSE
       datatype_flt = DATATYPE_FLT32
    ENDIF

    CALL vlr_add(p_diag_list, TRIM(listname), patch_id=p_patch%id, lrestart=.TRUE. )

    cf_desc    = t_cf_var('process_id', 'i','process ID', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'process_id', p_diag%process_id, &
         GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         ldims=shape2d_c)

    !Wind 10m group
    cf_desc    = t_cf_var('u_10m', 'm s-1 ','zonal wind in 10m', datatype_flt)
    grib2_desc = grib2_var(0, 2, 2, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'u_10m', p_diag%u_10m,                         &
         & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M, cf_desc, grib2_desc,       &
         & ldims=shape2d_c, lrestart=.FALSE., in_group=groups("wind_10m") )

    cf_desc    = t_cf_var('v_10m', 'm s-1 ','meridional wind in 10m', datatype_flt)
    grib2_desc = grib2_var(0, 2, 3, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'v_10m', p_diag%v_10m,                         &
         & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M, cf_desc, grib2_desc,       &
         & ldims=shape2d_c, lrestart=.FALSE., in_group=groups("wind_10m") )

    cf_desc    = t_cf_var('sp_10m', 'm s-1 ','wind speed in 10m', datatype_flt)
    grib2_desc = grib2_var(0, 2, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'sp_10m', p_diag%sp_10m,                       &
         & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M, cf_desc, grib2_desc,       &
         & ldims=shape2d_c, lrestart=.TRUE., in_group=groups("wind_10m") )

    cf_desc    = t_cf_var('dir_10m', 'deg ','wind direction in 10m', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'dir_10m', p_diag%dir_10m,                       &
         & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M, cf_desc, grib2_desc,       &
         & ldims=shape2d_c, lrestart=.TRUE., in_group=groups("wind_10m") )

    !Wave physics group
    cf_desc    = t_cf_var('emean', 'm^2', 'total wave energy', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'emean', p_diag%emean,                 &
         & GRID_UNSTRUCTURED_CELL,  ZA_SURFACE, cf_desc, grib2_desc,       &
         & ldims=shape2d_c, in_group=groups("wave_phy") )

    cf_desc    = t_cf_var('emeanws', 'm^2', 'wind sea wave energy', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'emeanws', p_diag%emeanws,                 &
         & GRID_UNSTRUCTURED_CELL,  ZA_SURFACE, cf_desc, grib2_desc,       &
         & ldims=shape2d_c, in_group=groups("wave_phy") )

    cf_desc    = t_cf_var('femean', 'm^2', 'mean frequency wave energy', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'femean', p_diag%femean,                 &
         & GRID_UNSTRUCTURED_CELL,  ZA_SURFACE, cf_desc, grib2_desc,       &
         & ldims=shape2d_c, in_group=groups("wave_phy") )

    cf_desc    = t_cf_var('f1mean', 'm^2', 'mean frequency wave energy based on F-moment', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'f1mean', p_diag%f1mean,                 &
         & GRID_UNSTRUCTURED_CELL,  ZA_SURFACE, cf_desc, grib2_desc,       &
         & ldims=shape2d_c, in_group=groups("wave_phy") )

    cf_desc    = t_cf_var('femeanws', 'm^2', 'mean frequency wind sea wave energy', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'femeanws', p_diag%femeanws,                 &
         & GRID_UNSTRUCTURED_CELL,  ZA_SURFACE, cf_desc, grib2_desc,       &
         & ldims=shape2d_c, in_group=groups("wave_phy") )

    cf_desc    = t_cf_var('LLWS', '-', 'wind sea, swell splitting flag ', datatype_int)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'LLWS', p_diag%LLWS,                 &
         & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M, cf_desc, grib2_desc,       &
         & ldims=shape3d_tr_c )

    cf_desc    = t_cf_var('ALPHAJ', '-', 'JONSWAP alpha', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'ALPHAJ', p_diag%ALPHAJ,                 &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,       &
         & ldims=shape2d_c , in_group=groups("wave_phy"))

    cf_desc    = t_cf_var('FP', 'Hz', 'JONSWAP peak frequency', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'FP', p_diag%FP,                 &
         & GRID_UNSTRUCTURED_CELL,  ZA_SURFACE, cf_desc, grib2_desc,       &
         & ldims=shape2d_c , in_group=groups("wave_phy"))

    cf_desc    = t_cf_var('ET', '-', 'JONSWAP spectra', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'ET', p_diag%ET,                 &
         & GRID_UNSTRUCTURED_CELL,  ZA_SURFACE, cf_desc, grib2_desc,       &
         & ldims=shape3d_freq_c)

    cf_desc    = t_cf_var('FLMINFR', '-', 'THE MINIMUM VALUE IN SPECTRAL', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'FLMINFR', p_diag%FLMINFR,                 &
         & GRID_UNSTRUCTURED_CELL,  ZA_SURFACE, cf_desc, grib2_desc,       &
         & ldims=shape3d_freq_c)

    cf_desc    = t_cf_var('friction_velocity', 'm s-1', 'friction velocity', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'ustar', p_diag%ustar,                 &
         & GRID_UNSTRUCTURED_CELL,  ZA_SURFACE, cf_desc, grib2_desc,       &
         & ldims=shape2d_c , in_group=groups("wave_phy"))

    cf_desc    = t_cf_var('roughness length', 'm', 'roughness length', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'z0', p_diag%z0,                 &
         & GRID_UNSTRUCTURED_CELL,  ZA_SURFACE, cf_desc, grib2_desc,       &
         & ldims=shape2d_c , in_group=groups("wave_phy"))

    cf_desc    = t_cf_var('wave_stress', '(m/s)**2', 'wave stress', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'tauw', p_diag%tauw,                 &
         & GRID_UNSTRUCTURED_CELL,  ZA_SURFACE, cf_desc, grib2_desc,       &
         & ldims=shape2d_c , in_group=groups("wave_phy"))

    cf_desc    = t_cf_var('fl', '-', 'DIAG. MTRX OF FUNC. DERIVATIVE', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'fl', p_diag%fl,                 &
         & GRID_UNSTRUCTURED_CELL,  ZA_SURFACE, cf_desc, grib2_desc,       &
         & ldims=shape3d_tr_c)

    cf_desc    = t_cf_var('sl', '-', 'TOTAL SOURCE FUNCTION', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'sl', p_diag%sl,                 &
         & GRID_UNSTRUCTURED_CELL,  ZA_SURFACE, cf_desc, grib2_desc,       &
         & ldims=shape3d_tr_c)

    ! wave output group
    cf_desc    = t_cf_var('Hs', 'm', 'significant wave height', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'Hs', p_diag%Hs,  &
         &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,      &
         &           ldims=shape2d_c, in_group=groups("wave_short"))!,

    cf_desc    = t_cf_var('Hs_dir', 'deg', 'mean direction of total significant wave height', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'Hs_dir', p_diag%Hs_dir,  &
         &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,      &
         &           ldims=shape2d_c, in_group=groups("wave_short"))

  END SUBROUTINE new_wave_state_diag_list

  !>
  !! Destruction of wave-specific variable lists and memory deallocation
  !!
  !! @par Revision History
  !! Initial revision by Mikhail Dobrynin, DWD (2022-12-20)
  !!
  SUBROUTINE destruct_wave_state ()

    INTEGER :: ist
    INTEGER :: jg, jt
    CHARACTER(len=*), PARAMETER :: routine = modname//'::destruct_wave_state'

    DO jg = 1, n_dom
      ! delete prognostic state list elements
      DO jt = 1, SIZE(p_wave_state_lists(jg)%prog_list(:))
        CALL vlr_del(p_wave_state_lists(jg)%prog_list(jt))
      ENDDO

      ! delete tracer list elements
      DO jt = 1, SIZE(p_wave_state_lists(jg)%tracer_list(:))
        CALL vlr_del(p_wave_state_lists(jg)%tracer_list(jt))
      ENDDO

      ! delete diagnostics state list elements
      CALL vlr_del(p_wave_state_lists(jg)%diag_list)

      ! deallocate state lists and arrays
      DEALLOCATE(p_wave_state_lists(jg)%prog_list, stat=ist)
      IF (ist /= success) THEN
        CALL finish(TRIM(routine),'deallocation for prog_list array failed')
      END IF

      DEALLOCATE(p_wave_state(jg)%prog, stat=ist)
      IF (ist /= success) THEN
        CALL finish(TRIM(routine),'deallocation of prognostic state array failed')
      END IF

      DEALLOCATE(p_wave_state_lists(jg)%tracer_list, stat=ist)
      IF (ist /= success) THEN
        CALL finish(TRIM(routine),'deallocation for tracer_list array failed')
      END IF
    ENDDO

    ! deallocate states
    DEALLOCATE(p_wave_state, p_wave_state_lists, stat=ist)
    IF (ist /= success) THEN
       CALL finish(TRIM(routine),'deallocation for p_wave_state failed')
    END IF

  END SUBROUTINE destruct_wave_state


END MODULE mo_wave_state
