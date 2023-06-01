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
  USE mo_cdi_constants,        ONLY: GRID_UNSTRUCTURED_CELL, GRID_CELL, &
       &                             GRID_UNSTRUCTURED_EDGE, GRID_EDGE
  USE mo_cdi,                  ONLY: DATATYPE_FLT32, DATATYPE_FLT64, GRID_UNSTRUCTURED, &
       &                             DATATYPE_PACK16, DATATYPE_INT
  USE mo_zaxis_type,           ONLY: ZA_SURFACE
  USE mo_cf_convention,        ONLY: t_cf_var
  USE mo_grib2,                ONLY: t_grib2_var, grib2_var
  USE mo_io_config,            ONLY: lnetcdf_flt64_output
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

    CHARACTER(len=max_char_length) :: listname
    CHARACTER(len=*), PARAMETER :: routine = modname//'::construct_wave_state'

    INTEGER :: ntl,  &! local number of timelevels
         ist,        &! status
         jg,         &! grid level counter
         jt           ! time level counter

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

          WRITE(listname,'(a,i2.2,a,i2.2)') 'wave_state_prog_of_domain_',jg, &
               &                               '_and_timelev_',jt

          ! Build prog state list
          ! includes memory allocation
          CALL new_wave_state_prog_list(p_patch(jg), p_wave_state(jg)%prog(jt), &
               & p_wave_state_lists(jg)%prog_list(jt), &
               & listname, jt)

          WRITE(listname,'(a,i2.2,a,i2.2)') 'wave_state_tracer_of_domain_',jg, &
               & '_and_timelev_',jt
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
    TYPE(t_patch),        INTENT(IN) :: p_patch ! current patch
    TYPE(t_var_list_ptr), INTENT(IN) :: from_var_list ! source list to be referenced
    TYPE(t_var_list_ptr), INTENT(INOUT) :: p_tracer_list ! new tracer list (containing all tracers)
    CHARACTER(*),         INTENT(IN) :: listname
    ! local
    TYPE (t_var_metadata),         POINTER :: from_info
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


  SUBROUTINE new_wave_state_prog_list(p_patch, p_prog, p_prog_list, listname, timelev)

    TYPE(t_patch),         INTENT(IN)    :: p_patch
    TYPE(t_wave_prog),     INTENT(INOUT) :: p_prog
    TYPE(t_var_list_ptr),  INTENT(INOUT) :: p_prog_list !< current prognostic state list
    CHARACTER(len=*),      INTENT(IN)    :: listname
    INTEGER,               INTENT(IN)    :: timelev

    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: nblks_c    !< number of cell blocks to allocate

    INTEGER :: ibits         !< "entropy" of horizontal slice
    INTEGER :: datatype_flt  !< floating point accuracy in NetCDF output

    CHARACTER(len=4)         :: suffix
    CHARACTER(len=VNAME_LEN) :: freq_ind_str, dir_ind_str
    CHARACTER(LEN=VNAME_LEN) :: tracer_container_name
    CHARACTER(len=VNAME_LEN) :: tracer_name

    TYPE(t_advection_config), POINTER :: advconf
    TYPE(t_wave_config),      POINTER :: wc


    INTEGER :: shape3d_c(3), shape4d_c(4)
    INTEGER :: jt, nlev

    !determine size of arrays
    nblks_c = p_patch%nblks_c

    ! number of vertical levels
    nlev    = p_patch%nlev

    ! pointer to advection_config(jg) to save some paperwork
    advconf => advection_config(p_patch%id)
    ! pointer to wave_config(jg) to save some paperwork
    wc => wave_config(p_patch%id)

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

    tracer_container_name = 'tracer'//suffix
    cf_desc    = t_cf_var('tracer', '', 'spectral bin of wave energy', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_prog_list, tracer_container_name, p_prog%tracer,             &
         &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
         &           ldims=shape4d_c ,                                           &
         &           lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE. )

    ALLOCATE(p_prog%tracer_ptr(ntracer))

    DO jt = 1, ntracer
       write(freq_ind_str,'(I0.3)') wc%freq_ind(jt)
       write(dir_ind_str,'(I0.3)') wc%dir_ind(jt)

       tracer_name = TRIM(advconf%tracer_names(jt))//'_'//TRIM(freq_ind_str)//'_'//TRIM(dir_ind_str)//suffix

       CALL add_ref( p_prog_list, tracer_container_name,                          &
            & TRIM(tracer_name), p_prog%tracer_ptr(jt)%p_3d,                      &
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
    INTEGER :: nblks_c, nblks_e
    INTEGER :: nfreqs, ndirs
    INTEGER :: jg, jt, jf
    INTEGER :: shape2d_c(2), shape2d_e(2)
    INTEGER :: shape3d_freq_c(3), shape3d_freq_e(3)
    INTEGER :: shape3d_freq_c_p4(3)
    INTEGER :: shape3d_tr_c(3), shape3d_tr_e(3)
    INTEGER :: shape1d_freq_p4(1), shape1d_dir_2(2)
    INTEGER :: shape4d_freq_p4_2_dir_18(4)

    CHARACTER(len=3) :: freq_ind_str, dir_ind_str
    CHARACTER(len=VNAME_LEN) :: out_name, sl_name, fl_name

    TYPE(t_wave_config),      POINTER :: wc

    ! pointer to wave_config(jg) to save some paperwork
    wc => wave_config(p_patch%id)

    jg      = p_patch%id
    nblks_c = p_patch%nblks_c
    nblks_e = p_patch%nblks_e

    nfreqs =  wave_config(jg)%nfreqs
    ndirs = wave_config(jg)%ndirs

    shape1d_freq_p4   = (/nfreqs+4/)
    shape1d_dir_2     = (/ndirs, 2/)
    shape2d_c         = (/nproma, nblks_c/)
    shape2d_e         = (/nproma, nblks_e/)
    shape3d_freq_c    = (/nproma, nblks_c, nfreqs/)
    shape3d_freq_c_p4 = (/nproma, nblks_c, nfreqs+4/)
    shape3d_freq_e    = (/nproma, nblks_e, nfreqs/)
    shape3d_tr_c      = (/nproma, nblks_c, ntracer/)
    shape3d_tr_e      = (/nproma, nblks_e, ntracer/)
    shape4d_freq_p4_2_dir_18 = (/18,nfreqs+4,2,ndirs/)

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

    !wave group velocity
    cf_desc    = t_cf_var('gv_c', 'm s-1', 'group velocity at cells', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'gv_c_freq', p_diag%gv_c,    &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,  &
         & ldims=shape3d_freq_c, &
         & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

    cf_desc    = t_cf_var('gv_e', 'm s-1', 'group velocity at edges', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
    CALL add_var( p_diag_list, 'gv_e_freq', p_diag%gv_e,    &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc,  &
         & ldims=shape3d_freq_e,&
         & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

    ALLOCATE(p_diag%freq_ptr(nfreqs))

    DO jf = 1, nfreqs
      write(freq_ind_str,'(I0.3)') jf
      out_name = 'gv_c_'//TRIM(freq_ind_str)
      CALL add_ref( p_diag_list, 'gv_c_freq',                                 &
            & TRIM(out_name), p_diag%freq_ptr(jf)%p_2d,       &
            & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                            &
            & t_cf_var(TRIM(out_name), 'm/s',TRIM(out_name), &
            & datatype_flt),                                                 &
            & grib2_var(0, 0, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL),       &
            & ref_idx=jf, ldims=shape2d_c, loutput=.TRUE.,                   &
            & in_group=groups("wave_phy_ext"))

      out_name = 'gv_e_'//TRIM(freq_ind_str)
      CALL add_ref(p_diag_list, 'gv_e_freq',                                  &
            & TRIM(out_name), p_diag%freq_ptr(jf)%p_2d,       &
            & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE,                            &
            & t_cf_var(TRIM(out_name), 'm/s',TRIM(out_name), &
            & datatype_flt),                                                 &
            & grib2_var(0, 0, 0, ibits, GRID_UNSTRUCTURED, GRID_EDGE),       &
            & ref_idx=jf, ldims=shape2d_e, loutput=.TRUE.,                   &
            & in_group=groups("wave_phy_ext"))

    END DO

    cf_desc    = t_cf_var('normal_group_velocity', 'm s-1', 'group velocity normal to edge', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
    CALL add_var(p_diag_list, 'gvn_e', p_diag%gvn_e,     &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc, &
         & ldims=shape3d_tr_e)

    cf_desc    = t_cf_var('tangential_group_velocity', 'm s-1', 'group velocity tangential to edge', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
    CALL add_var( p_diag_list, 'gvt_e', p_diag%gvt_e,     &
         &           GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc, &
         &           ldims=shape3d_tr_e)

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

    cf_desc    = t_cf_var('tm1', 's', 'wave tm1 period', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'tm1', p_diag%tm1,  &
         &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,      &
         &           ldims=shape2d_c, in_group=groups("wave_phy"))

    cf_desc    = t_cf_var('tm2', 's', 'wave tm1 period', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'tm2', p_diag%tm2,  &
         &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,      &
         &           ldims=shape2d_c, in_group=groups("wave_phy"))

    cf_desc    = t_cf_var('wave_num_c', '1/m', 'wave number at cell center', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'wave_num_c', p_diag%wave_num_c,                    &
         &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,      &
         &           ldims=shape3d_freq_c,                                         &
         &           lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE. )

    ALLOCATE(p_diag%wave_num_c_ptr(nfreqs))
    !
    DO jf = 1, nfreqs
      write(freq_ind_str,'(I0.3)') jf
      out_name = 'wave_num_c_'//TRIM(freq_ind_str)
      CALL add_ref( p_diag_list, 'wave_num_c',                               &
            & TRIM(out_name), p_diag%wave_num_c_ptr(jf)%p_2d,                &
            & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                            &
            & t_cf_var(TRIM(out_name), '1/m',TRIM(out_name), &
            & datatype_flt),                                                 &
            & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL), &
            & ref_idx=jf, ldims=shape2d_c, lrestart=.TRUE., loutput=.TRUE. )
    ENDDO


    cf_desc    = t_cf_var('wave_num_e', '1/m', 'wave number at edge midpoint', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
    CALL add_var( p_diag_list, 'wave_num_e', p_diag%wave_num_e,                    &
         &           GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc,      &
         &           ldims=shape3d_freq_e,                                         &
         &           lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE. )

    ALLOCATE(p_diag%wave_num_e_ptr(nfreqs))
    !
    DO jf = 1, nfreqs
      write(freq_ind_str,'(I0.3)') jf
      out_name = 'wave_num_e_'//TRIM(freq_ind_str)
      CALL add_ref( p_diag_list, 'wave_num_e',                               &
            & TRIM(out_name), p_diag%wave_num_e_ptr(jf)%p_2d,                &
            & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE,                            &
            & t_cf_var(TRIM(out_name), '1/m',TRIM(out_name),                 &
            & datatype_flt),                                                 &
            & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE), &
            & ref_idx=jf, ldims=shape2d_e, lrestart=.TRUE., loutput=.TRUE. )
    ENDDO


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

    cf_desc    = t_cf_var('akmean', '', 'Mean wavenumber based on SQRT(1/K)-moment, wm1', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'akmean', p_diag%akmean,                 &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,       &
         & ldims=shape2d_c, in_group=groups("wave_phy"))

    cf_desc    = t_cf_var('xkmean', '', 'Mean wavenumber based on SQRT(K)-moment, wm2', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'xkmean', p_diag%xkmean,                 &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,       &
         & ldims=shape2d_c, in_group=groups("wave_phy"))

    cf_desc    = t_cf_var('LLWS', '-', 'wind sea, swell splitting flag ', datatype_int)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'LLWS', p_diag%LLWS,                 &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,       &
         & ldims=shape3d_tr_c )

    cf_desc    = t_cf_var('last_prog_freq_ind', '-', 'last frequency index of the prognostic range', datatype_int)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'last_prog_freq_ind', p_diag%last_prog_freq_ind,                 &
         & GRID_UNSTRUCTURED_CELL,  ZA_SURFACE, cf_desc, grib2_desc,       &
         & ldims=shape2d_c, in_group=groups("wave_phy") )

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

    cf_desc    = t_cf_var('integrated_energy_flux', '-', 'integrated energy flux', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'phiaw', p_diag%phiaw,  &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,      &
         & ldims=shape2d_c , in_group=groups("wave_phy"))

    cf_desc    = t_cf_var('fl', '-', 'DIAG. MTRX OF FUNC. DERIVATIVE', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'fl', p_diag%fl,                 &
         & GRID_UNSTRUCTURED_CELL,  ZA_SURFACE, cf_desc, grib2_desc,       &
         & ldims=shape3d_tr_c, &
         & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

    cf_desc    = t_cf_var('sl', '-', 'TOTAL SOURCE FUNCTION', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'sl', p_diag%sl,                 &
         & GRID_UNSTRUCTURED_CELL,  ZA_SURFACE, cf_desc, grib2_desc,       &
         & ldims=shape3d_tr_c, &
         & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

    ALLOCATE(p_diag%tr_ptr(ntracer))

    DO jt = 1, ntracer
      write(freq_ind_str,'(I0.3)') wc%freq_ind(jt)
      write(dir_ind_str,'(I0.3)') wc%dir_ind(jt)

      sl_name = 'sl_'//TRIM(freq_ind_str)//'_'//TRIM(dir_ind_str)
      CALL add_ref( p_diag_list, 'sl',                                     &
            & sl_name, p_diag%tr_ptr(jt)%p_2d,       &
            & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                            &
            & t_cf_var(sl_name, '-',sl_name, &
            & datatype_flt),                                                 &
            & grib2_var(0, 0, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL),       &
            & ref_idx=jt, ldims=shape2d_c, loutput=.TRUE.,                   &
            & in_group=groups("wave_phy_ext"))

      fl_name = 'fl_'//TRIM(freq_ind_str)//'_'//TRIM(dir_ind_str)
      CALL add_ref( p_diag_list, 'fl',                                     &
            & fl_name, p_diag%tr_ptr(jt)%p_2d,       &
            & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                            &
            & t_cf_var(fl_name, '-',fl_name, &
            & datatype_flt),                                                 &
            & grib2_var(0, 0, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL),       &
            & ref_idx=jt, ldims=shape2d_c, loutput=.TRUE.,                   &
            & in_group=groups("wave_phy_ext"))
    END DO

    cf_desc    = t_cf_var('tauhf1', '-', 'high-fequency stress', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'tauhf1', p_diag%tauhf1,                 &
         & GRID_UNSTRUCTURED_CELL,  ZA_SURFACE, cf_desc, grib2_desc,       &
         & ldims=shape2d_c, in_group=groups("wave_phy") )

    cf_desc    = t_cf_var('phihf1', '-', 'high-frequency energy flux into ocean', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'phihf1', p_diag%phihf1,                 &
         & GRID_UNSTRUCTURED_CELL,  ZA_SURFACE, cf_desc, grib2_desc,       &
         & ldims=shape2d_c, in_group=groups("wave_phy") )

    cf_desc    = t_cf_var('tauhf', '-', 'high-fequency stress', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'tauhf', p_diag%tauhf,                 &
         & GRID_UNSTRUCTURED_CELL,  ZA_SURFACE, cf_desc, grib2_desc,       &
         & ldims=shape2d_c, in_group=groups("wave_phy") )

    cf_desc    = t_cf_var('phihf', '-', 'high-frequency energy flux into ocean', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'phihf', p_diag%phihf,                 &
         & GRID_UNSTRUCTURED_CELL,  ZA_SURFACE, cf_desc, grib2_desc,       &
         & ldims=shape2d_c, in_group=groups("wave_phy") )

    cf_desc    = t_cf_var('xlevtail', '-', 'tail level', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'xlevtail', p_diag%xlevtail,                 &
         & GRID_UNSTRUCTURED_CELL,  ZA_SURFACE, cf_desc, grib2_desc,       &
         & ldims=shape2d_c, in_group=groups("wave_phy") )

    ! nonlinear transfer function coefficients for shallow water
    cf_desc    = t_cf_var('enh', '-', 'nonlinear transfer function coefficients', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'enh', p_diag%enh,                 &
         & GRID_UNSTRUCTURED_CELL,  ZA_SURFACE, cf_desc, grib2_desc,       &
         & ldims=shape2d_c)

    ! for discrete approximation of nonlinear transfer
    cf_desc    = t_cf_var('IKP', '-', 'IKP', datatype_int)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'IKP', p_diag%IKP,                 &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc,       &
         & ldims=shape1d_freq_p4 )

    cf_desc    = t_cf_var('IKP1', '-', 'IKP1', datatype_int)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'IKP1', p_diag%IKP1,                 &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc,       &
         & ldims=shape1d_freq_p4 )

    cf_desc    = t_cf_var('IKM', '-', 'IKM', datatype_int)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'IKM', p_diag%IKM,                 &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc,       &
         & ldims=shape1d_freq_p4 )

    cf_desc    = t_cf_var('IKM1', '-', 'IKM1', datatype_int)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'IKM1', p_diag%IKM1,                 &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc,       &
         & ldims=shape1d_freq_p4 )

    cf_desc    = t_cf_var('K1W', '-', 'K1W', datatype_int)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'K1W', p_diag%K1W,                 &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc,       &
         & ldims=shape1d_dir_2 )

    cf_desc    = t_cf_var('K2W', '-', 'K2W', datatype_int)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'K2W', p_diag%K2W,                 &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc,       &
         & ldims=shape1d_dir_2 )

    cf_desc    = t_cf_var('K11W', '-', 'K11W', datatype_int)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'K11W', p_diag%K11W,                 &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc,       &
         & ldims=shape1d_dir_2 )

    cf_desc    = t_cf_var('K21W', '-', 'K21W', datatype_int)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'K21W', p_diag%K21W,                 &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc,       &
         & ldims=shape1d_dir_2 )

    cf_desc    = t_cf_var('non_lin_tr_ind', '-', 'non_lin_tr_ind', datatype_int)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'non_lin_tr_ind', p_diag%non_lin_tr_ind,   &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc,       &
         & ldims=shape4d_freq_p4_2_dir_18 )

    cf_desc    = t_cf_var('JA1', '-', 'JA1', datatype_int)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'JA1', p_diag%JA1,                 &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc,       &
         & ldims=shape1d_dir_2 )

    cf_desc    = t_cf_var('JA2', '-', 'JA2', datatype_int)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'JA2', p_diag%JA2,                 &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc,       &
         & ldims=shape1d_dir_2)

    cf_desc    = t_cf_var('AF11', '-', 'AF11', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'AF11', p_diag%AF11,                 &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc,       &
         & ldims=shape1d_freq_p4 )

    cf_desc    = t_cf_var('FKLAP', '-', 'FKLAP', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'FKLAP', p_diag%FKLAP,                 &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc,       &
         & ldims=shape1d_freq_p4 )

    cf_desc    = t_cf_var('FKLAP1', '-', 'FKLAP1', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'FKLAP1', p_diag%FKLAP1,                 &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc,       &
         & ldims=shape1d_freq_p4 )

    cf_desc    = t_cf_var('FKLAM', '-', 'FKLAM', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'FKLAM', p_diag%FKLAM,                 &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc,       &
         & ldims=shape1d_freq_p4 )

    cf_desc    = t_cf_var('FKLAM1', '-', 'FKLAM1', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'FKLAM1', p_diag%FKLAM1,                 &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc,       &
         & ldims=shape1d_freq_p4 )

    ! wave output group
    cf_desc    = t_cf_var('Hs', 'm', 'significant wave height', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'Hs', p_diag%Hs,  &
         &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,      &
         &           ldims=shape2d_c, in_group=groups("wave_short"))

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
