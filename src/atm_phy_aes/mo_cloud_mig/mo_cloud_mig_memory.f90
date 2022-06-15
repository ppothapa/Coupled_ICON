!>
!! Global variables for the cloud microphsics "graupel" parameterization
!!
!! @author Marco Giorgetta (MPI-M)
!!
!! @par Revision History
!! First version by Marco Giorgetta, 2020-01-30.
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_cloud_mig_memory

  USE mo_exception               ,ONLY: message, finish
  USE mtime                      ,ONLY: OPERATOR(>)

  USE mo_model_domain            ,ONLY: t_patch
  USE mo_parallel_config         ,ONLY: nproma
  USE mo_time_config             ,ONLY: time_config
  USE mo_aes_phy_config          ,ONLY: aes_phy_tc, dt_zero
  USE mo_io_config               ,ONLY: lnetcdf_flt64_output
  USE mo_name_list_output_config ,ONLY: is_variable_in_output

  USE mo_impl_constants          ,ONLY: success, vintp_method_lin
  USE mo_cdi_constants           ,ONLY: grid_unstructured_cell, grid_cell
  USE mo_cdi                     ,ONLY: grid_unstructured,                 &
       &                                datatype_pack16,                   &
       &                                datatype_flt32,  datatype_flt64,   &
       &                                datatype_int32,                    &
       &                                tstep_instant, tstep_constant

  USE mo_var_list                ,ONLY: add_var, t_var_list_ptr
  USE mo_var_list_register       ,ONLY: vlr_add, vlr_del
  USE mo_var_metadata            ,ONLY: create_vert_interp_metadata, vintp_types
  USE mo_cf_convention           ,ONLY: t_cf_var
  USE mo_grib2                   ,ONLY: grib2_var
  USE mo_zaxis_type              ,ONLY: za_reference, za_surface

  USE mo_cloud_mig_types         ,ONLY: t_cloud_mig_input, t_cloud_mig_output

  ! include definition for "__acc_attach(ptr)"
#include "add_var_acc_macro.inc"

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: cloud_mig_input, cloud_mig_output
  PUBLIC :: cloud_mig_list
  PUBLIC :: construct_cloud_mig_memory
  PUBLIC :: destruct_cloud_mig_memory

  CHARACTER(len=*), PARAMETER                    :: thismodule = 'mo_cloud_mig_memory'
  TYPE(t_cloud_mig_input)  , ALLOCATABLE, TARGET :: cloud_mig_input(:)  !< shape: (ng)
  TYPE(t_cloud_mig_output) , ALLOCATABLE, TARGET :: cloud_mig_output(:) !< shape: (ng)
  TYPE(t_var_list_ptr)     , ALLOCATABLE         :: cloud_mig_list(:)   !< shape: (ng)

CONTAINS

  !!--------------------------------------------------------------------------
  !!                SUBROUTINES FOR BUILDING AND DELETING VARIABLE LISTS 
  !!--------------------------------------------------------------------------
  !>
  !! Top-level procedure for building the state
  !!
  SUBROUTINE construct_cloud_mig_memory( patch_array )

    TYPE(t_patch),INTENT(IN) :: patch_array(:)

    INTEGER :: ng, jg, ist
    INTEGER :: nlev, nblks

    !---

    CALL message(thismodule,'Construction of cloud_mig_list, cloud_mig_input and cloud mig_output started.')

    ! allocate pointer arrays for the pointer arrays cloud_mig_memory and cloud_mig_list

    ng = SIZE(patch_array)

    ALLOCATE( cloud_mig_list(ng), STAT=ist)
    IF (ist/=success) CALL finish(thismodule, 'allocation of cloud_mig_list(ng) failed')

    ALLOCATE( cloud_mig_input(ng), STAT=ist)
    IF (ist/=success) CALL finish(thismodule, 'allocation of cloud_mig_input(ng) failed')

    ALLOCATE( cloud_mig_output(ng), STAT=ist)
    IF (ist/=success) CALL finish(thismodule, 'allocation of cloud_mig_output(ng) failed')

    !$ACC ENTER DATA CREATE(cloud_mig_input, cloud_mig_output)

    ! build lists and allocate memory for all grids where the cloud microphysics "graupel" is used

    DO jg = 1,ng
       IF (aes_phy_tc(jg)%dt_mig > dt_zero) THEN
          !
          nlev   = patch_array(jg)%nlev
          nblks  = patch_array(jg)%nblks_c
          !
          CALL construct_cloud_mig_list( jg,                  &
               &                         nproma, nlev, nblks, &
               &                         cloud_mig_list(jg),  &
               &                         cloud_mig_input(jg), &
               &                         cloud_mig_output(jg) )
          !
       END IF
    END DO

    CALL message(thismodule,'Construction of cloud_mig_list, cloud_mig_input and cloud mig_output finished.')

  END SUBROUTINE construct_cloud_mig_memory
  !--------------------------------------------------------------------


  !--------------------------------------------------------------------
  !>
  !! Release memory used by the state variable arrays and list arrays
  !!
  SUBROUTINE destruct_cloud_mig_memory

    INTEGER :: ng   !< total number of grids
    INTEGER :: jg   !< grid index
    INTEGER :: ist  !< system status code

    !---
    CALL message(thismodule,'Destruction of cloud_mig_list, cloud_mig_input and cloud_mig_output started.')

    ng = SIZE(cloud_mig_input)

    DO jg = 1,ng
       IF (aes_phy_tc(jg)%dt_mig > dt_zero) THEN
          !
          CALL vlr_del( cloud_mig_list(jg) )
          !
       END IF
    END DO
    
    !$ACC EXIT DATA DELETE(cloud_mig_input, cloud_mig_output)

    DEALLOCATE( cloud_mig_list, STAT=ist )
    IF (ist/=success) CALL finish(thismodule, 'deallocation of cloud_mig_list failed')

    DEALLOCATE( cloud_mig_input, STAT=ist )
    IF (ist/=success) CALL finish(thismodule, 'deallocation of cloud_mig_input failed')

    DEALLOCATE( cloud_mig_output, STAT=ist )
    IF (ist/=success) CALL finish(thismodule, 'deallocation of cloud_mig_output failed')

    CALL message(thismodule,'Destruction of cloud_mig_list, cloud_mig_input and cloud mig_output finished.')

  END SUBROUTINE destruct_cloud_mig_memory
  !--------------------------------------------------------------------


  !--------------------------------------------------------------------
  !>
  SUBROUTINE construct_cloud_mig_list( jg,                  &
       &                               nproma, nlev, nblks, &
       &                               cloud_mig_list,      &
       &                               cloud_mig_input,     &
       &                               cloud_mig_output )

    INTEGER                 , INTENT(in)    :: jg                  !< grid index
    INTEGER                 , INTENT(in)    :: nproma, nlev, nblks !< size of dimensions

    TYPE(t_var_list_ptr)    , INTENT(inout) :: cloud_mig_list      !< pointers for list of variables
    TYPE(t_cloud_mig_input) , INTENT(inout) :: cloud_mig_input     !< pointers for input variables
    TYPE(t_cloud_mig_output), INTENT(inout) :: cloud_mig_output    !< pointers for output variables

    ! Local variables

    CHARACTER(len= 2) :: cg
    CHARACTER(len=20) :: listname

    INTEGER           :: shape2d(2), shape3d(3)
    INTEGER           :: datatype_grb
    INTEGER           :: datatype_flt, datatype_int

    WRITE(cg,'(i2.2)') jg
    CALL message('construct_cloud_mig_list','create list and allocate memory for jg ='//cg)

    ! number of bits for data representation in grib2
    datatype_grb = datatype_pack16

    ! number of bits for data representation in netcdf
    datatype_flt = MERGE(DATATYPE_FLT64, DATATYPE_FLT32, lnetcdf_flt64_output)
    datatype_int = datatype_int32

    shape2d  = (/nproma,       nblks/)
    shape3d  = (/nproma, nlev, nblks/)

    ! define list name
    WRITE(listname,'(a,i2.2)') 'cloud_mig_memory_D',jg

    ! register the cloud_mig_list for grid jg
    CALL vlr_add( cloud_mig_list             ,&
         &        listname                   ,&
         &        patch_id  = jg             ,&
         &        loutput   = .TRUE.         ,&
         &        lrestart  = .FALSE.        ,&
         &        linitial  = .FALSE.        )

    ! Input parameters
    ! ----------------
    !
    ! These fields are constructed only if they are requested for output
    !
    IF ( is_variable_in_output(var_name='jg_mig') ) THEN
       CALL add_var( this_list   = cloud_mig_list                                           ,&
            &        varname     = 'jg_mig'                                                 ,&
            &        ptr         = cloud_mig_input%jg                                       ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_surface                                               ,&
            &        ldims       = shape2d                                                  ,&
            &        cf          = t_cf_var ('grid_index',                                   &
            &                                '-',                                            &
            &                                'grid index (cloud_mig input)',                 &
            &                                datatype_int)                                  ,&
            &        grib2       = grib2_var(255,255,255,                                    &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_constant                                           ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_mig_input%jg)
    END IF
    !
    IF ( is_variable_in_output(var_name='jcs_mig') ) THEN
       CALL add_var( this_list   = cloud_mig_list                                           ,&
            &        varname     = 'jcs_mig'                                                ,&
            &        ptr         = cloud_mig_input%jcs                                      ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_surface                                               ,&
            &        ldims       = shape2d                                                  ,&
            &        cf          = t_cf_var ('column_start_index',                           &
            &                                '-',                                            &
            &                                'column start index (cloud_mig input)',         &
            &                                datatype_int)                                  ,&
            &        grib2       = grib2_var(255,255,255,                                    &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_constant                                           ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_mig_input%jcs)
    END IF
    !
    IF ( is_variable_in_output(var_name='jce_mig') ) THEN
       CALL add_var( this_list   = cloud_mig_list                                           ,&
            &        varname     = 'jce_mig'                                                ,&
            &        ptr         = cloud_mig_input%jce                                      ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_surface                                               ,&
            &        ldims       = shape2d                                                  ,&
            &        cf          = t_cf_var ('column_end_index',                             &
            &                                '-',                                            &
            &                                'column end index (cloud_mig input)',           &
            &                                datatype_int)                                  ,&
            &        grib2       = grib2_var(255,255,255,                                    &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_constant                                           ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_mig_input%jce)
    END IF
    !
    IF ( is_variable_in_output(var_name='msg_level_mig') ) THEN
       CALL add_var( this_list   = cloud_mig_list                                           ,&
            &        varname     = 'msg_level_mig'                                          ,&
            &        ptr         = cloud_mig_input%msg_level                                ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_surface                                               ,&
            &        ldims       = shape2d                                                  ,&
            &        cf          = t_cf_var ('message_level',                                &
            &                                '-',                                            &
            &                                'message level (cloud_mig input)',              &
            &                                datatype_int)                                  ,&
            &        grib2       = grib2_var(255,255,255,                                    &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_constant                                           ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_mig_input%msg_level)
    END IF
    !
    IF ( is_variable_in_output(var_name='pdtime_mig') ) THEN
       CALL add_var( this_list   = cloud_mig_list                                           ,&
            &        varname     = 'pdtime_mig'                                             ,&
            &        ptr         = cloud_mig_input%pdtime                                   ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_surface                                               ,&
            &        ldims       = shape2d                                                  ,&
            &        cf          = t_cf_var ('time_step',                                    &
            &                                's',                                            &
            &                                'time step (cloud_mig input)',                  &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(255,255,255,                                    &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_constant                                           ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_mig_input%pdtime)
    END IF

    ! Input fields
    ! ------------
    !
    ! These fields are constructed only if they are requested for output
    !
    IF ( is_variable_in_output(var_name='dz_mig') ) THEN
       CALL add_var( this_list   = cloud_mig_list                                           ,&
            &        varname     = 'dz_mig'                                                 ,&
            &        ptr         = cloud_mig_input%dz                                       ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('cell_thickness',                               &
            &                                'm',                                            &
            &                                'cell thickness (cloud_mig input)',             &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,3,12,                                         &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_mig_input%dz)
    END IF
    !
    IF ( is_variable_in_output(var_name='rho_mig') ) THEN
       CALL add_var( this_list   = cloud_mig_list                                           ,&
            &        varname     = 'rho_mig'                                                ,&
            &        ptr         = cloud_mig_input%rho                                      ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('air_density',                                  &
            &                                'kg/m3',                                        &
            &                                'air density (cloud_mig input)',                &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,3,10,                                         &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_mig_input%rho)
    END IF
    !
    IF ( is_variable_in_output(var_name='pf_mig') ) THEN
       CALL add_var( this_list   = cloud_mig_list                                           ,&
            &        varname     = 'pf_mig'                                                 ,&
            &        ptr         = cloud_mig_input%pf                                       ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('air_pressure',                                 &
            &                                'kg/m3',                                        &
            &                                'air pressure (cloud_mig input)',               &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,3,0,                                          &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_mig_input%pf)
    END IF
    !
    IF ( is_variable_in_output(var_name='cpair_mig') ) THEN
       CALL add_var( this_list   = cloud_mig_list                                           ,&
            &        varname     = 'cpair_mig'                                              ,&
            &        ptr         = cloud_mig_input%cpair                                    ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('air_specific_heat',                            &
            &                                'J/K/kg',                                       &
            &                                'specific heat of air at constant pressure '//  &
            &                                '(cloud_mig input)',                            &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,0,255,                                        &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_mig_input%cpair)
    END IF
    !
    IF ( is_variable_in_output(var_name='ta_mig') ) THEN
       CALL add_var( this_list   = cloud_mig_list                                           ,&
            &        varname     = 'ta_mig'                                                 ,&
            &        ptr         = cloud_mig_input%ta                                       ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('air_temperature',                              &
            &                                'K',                                            &
            &                                'air temperature (cloud_mig input)',            &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,0,0,                                          &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_mig_input%ta)
    END IF
    !
    IF ( is_variable_in_output(var_name='qv_mig') ) THEN
       CALL add_var( this_list   = cloud_mig_list                                           ,&
            &        varname     = 'qv_mig'                                                 ,&
            &        ptr         = cloud_mig_input%qv                                       ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('specific_humidity',                            &
            &                                'kg/kg',                                        &
            &                                'specific humidity (cloud_mig input)',          &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,0,                                          &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_mig_input%qv)
    END IF
    !
    IF ( is_variable_in_output(var_name='qc_mig') ) THEN
       CALL add_var( this_list   = cloud_mig_list                                           ,&
            &        varname     = 'qc_mig'                                                 ,&
            &        ptr         = cloud_mig_input%qc                                       ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('mass_fraction_of_cloud_liquid_water_in_air',   &
            &                                'kg/kg',                                        &
            &                                'mass fraction of cloud liquid water in air '// &
            &                                '(cloud_mig input)',                            &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,83,                                         &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_mig_input%qc)
    END IF
    !
    IF ( is_variable_in_output(var_name='qi_mig') ) THEN
       CALL add_var( this_list   = cloud_mig_list                                           ,&
            &        varname     = 'qi_mig'                                                 ,&
            &        ptr         = cloud_mig_input%qi                                       ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('mass_fraction_of_cloud_ice_in_air',            &
            &                                'kg/kg',                                        &
            &                                'mass fraction of cloud ice in air '//          &
            &                                '(cloud_mig input)',                            &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,84,                                         &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_mig_input%qi)
    END IF
    !
    IF ( is_variable_in_output(var_name='qr_mig') ) THEN
       CALL add_var( this_list   = cloud_mig_list                                           ,&
            &        varname     = 'qr_mig'                                                 ,&
            &        ptr         = cloud_mig_input%qr                                       ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('mass_fraction_of_rain_in_air',                 &
            &                                'kg/kg',                                        &
            &                                'mass fraction of rain in air '//               &
            &                                '(cloud_mig input)',                            &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,85,                                         &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_mig_input%qr)
    END IF
    !
    IF ( is_variable_in_output(var_name='qs_mig') ) THEN
       CALL add_var( this_list   = cloud_mig_list                                           ,&
            &        varname     = 'qs_mig'                                                 ,&
            &        ptr         = cloud_mig_input%qs                                       ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('mass_fraction_of_snow_in_air',                 &
            &                                'kg/kg',                                        &
            &                                'mass fraction of snow in air '//               &
            &                                '(cloud_mig input)',                            &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,86,                                         &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_mig_input%qs)
    END IF
    !
    IF ( is_variable_in_output(var_name='qg_mig') ) THEN
       CALL add_var( this_list   = cloud_mig_list                                           ,&
            &        varname     = 'qg_mig'                                                 ,&
            &        ptr         = cloud_mig_input%qg                                       ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('mass_fraction_of_graupel_in_air',              &
            &                                'kg/kg',                                        &
            &                                'mass fraction of graupel in air '//            &
            &                                '(cloud_mig input)',                            &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,32,                                         &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_mig_input%qg)
    END IF


    ! Output fields
    ! -------------
    !
    ! These fields are constructed only if they are requested for output
    ! or needed for recycling in the time stepping.
    !
    ! tendencies in the atmosphere
    !
    IF ( aes_phy_tc(jg)%dt_mig > time_config%tc_dt_dyn(jg) .OR.                            &
         & is_variable_in_output(var_name='tend_ta_mig') ) THEN
       CALL add_var( this_list   = cloud_mig_list                                           ,&
            &        varname     = 'tend_ta_mig'                                            ,&
            &        ptr         = cloud_mig_output%tend_ta_mig                             ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('temperature_tendency_graupel','K s-1',         &
            &                       'temperature tendency due to graupel processes (cp) '//  &
            &                                '(cloud_mig output)',                           &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,0,203,                                        &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_mig_output%tend_ta_mig)
    END IF
    !
    IF ( aes_phy_tc(jg)%dt_mig > time_config%tc_dt_dyn(jg) .OR.                            &
         & is_variable_in_output(var_name='tend_qhus_mig') ) THEN
       CALL add_var( this_list   = cloud_mig_list                                           ,&
            &        varname     = 'tend_qhus_mig'                                          ,&
            &        ptr         = cloud_mig_output%tend_qv_mig                             ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('tend_qhus_mig',  'kg kg-1 s-1',                &
            &                                'tendency of mass mixing ratio of tracer '//    &
            &                                'qhus  due to graupel processes '//             &
            &                                '(cloud_mig output)',                           &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,203,                                        &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_mig_output%tend_qv_mig)
    END IF
    !
    IF ( aes_phy_tc(jg)%dt_mig > time_config%tc_dt_dyn(jg) .OR.                            &
         & is_variable_in_output(var_name='tend_qclw_mig') ) THEN
       CALL add_var( this_list   = cloud_mig_list                                           ,&
            &        varname     = 'tend_qclw_mig'                                          ,&
            &        ptr         = cloud_mig_output%tend_qc_mig                             ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('tend_qclw_mig',  'kg kg-1 s-1'                ,&
            &                                'tendency of mass mixing ratio of tracer '//    &
            &                                'qclw  due to graupel processes '//             &
            &                                '(cloud_mig output)',                           &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,6,203,                                        &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_mig_output%tend_qc_mig)
    END IF
    !
    IF ( aes_phy_tc(jg)%dt_mig > time_config%tc_dt_dyn(jg) .OR.                            &
         & is_variable_in_output(var_name='tend_qcli_mig') ) THEN
       CALL add_var( this_list   = cloud_mig_list                                           ,&
            &        varname     = 'tend_qcli_mig'                                          ,&
            &        ptr         = cloud_mig_output%tend_qi_mig                             ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('tend_qcli_mig',  'kg kg-1 s-1'                ,&
            &                                'tendency of mass mixing ratio of tracer '//    &
            &                                'qcli  due to graupel processes '//             &
            &                                '(cloud_mig output)',                           &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,6,213,                                        &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_mig_output%tend_qi_mig)
    END IF
    !
    IF ( aes_phy_tc(jg)%dt_mig > time_config%tc_dt_dyn(jg) .OR.                            &
         & is_variable_in_output(var_name='tend_qr_mig') ) THEN
       CALL add_var( this_list   = cloud_mig_list                                           ,&
            &        varname     = 'tend_qr_mig'                                            ,&
            &        ptr         = cloud_mig_output%tend_qr_mig                             ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('tend_qr_mig',  'kg kg-1 s-1'                  ,&
            &                                'tendency of mass mixing ratio of tracer '//    &
            &                                'qr  due to graupel processes '//               &
            &                                '(cloud_mig output)',                           &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(255,255,255,                                    &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_mig_output%tend_qr_mig)
    END IF
    !
    IF ( aes_phy_tc(jg)%dt_mig > time_config%tc_dt_dyn(jg) .OR.                            &
         & is_variable_in_output(var_name='tend_qs_mig') ) THEN
       CALL add_var( this_list   = cloud_mig_list                                           ,&
            &        varname     = 'tend_qs_mig'                                            ,&
            &        ptr         = cloud_mig_output%tend_qs_mig                             ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('tend_qs_mig',  'kg kg-1 s-1'                  ,&
            &                                'tendency of mass mixing ratio of tracer '//    &
            &                                'qs  due to graupel processes '//               &
            &                                '(cloud_mig output)',                           &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(255,255,255,                                    &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_mig_output%tend_qs_mig)
    END IF
    !
    IF ( aes_phy_tc(jg)%dt_mig > time_config%tc_dt_dyn(jg) .OR.                            &
         & is_variable_in_output(var_name='tend_qg_mig') ) THEN
       CALL add_var( this_list   = cloud_mig_list                                           ,&
            &        varname     = 'tend_qg_mig'                                            ,&
            &        ptr         = cloud_mig_output%tend_qg_mig                             ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_reference                                             ,&
            &        ldims       = shape3d                                                  ,&
            &        cf          = t_cf_var ('tend_qg_mig',  'kg kg-1 s-1'                  ,&
            &                                'tendency of mass mixing ratio of tracer '//    &
            &                                'qg  due to graupel processes '//               &
            &                                '(cloud_mig output)',                           &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(255,255,255,                                    &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        vert_interp = create_vert_interp_metadata(                              &
            &                             vert_intp_type   = vintp_types("P","Z","I"),       &
            &                             vert_intp_method = vintp_method_lin)              ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_mig_output%tend_qg_mig)
    END IF
    !
    ! fluxes at the surface
    !
    IF ( is_variable_in_output(var_name='pr_rain') ) THEN
       CALL add_var( this_list   = cloud_mig_list                                           ,&
            &        varname     = 'pr_rain'                                                ,&
            &        ptr         = cloud_mig_output%pr_rain                                 ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_surface                                               ,&
            &        ldims       = shape2d                                                  ,&
            &        cf          = t_cf_var ('stratiform_rainfall_flux',                     &
            &                                'kg/m2/s',                                      &
            &                                'stratiform rainfall flux '//                   &
            &                                '(cloud_mig output)',                           &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,77,                                         &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_mig_output%pr_rain)
    END IF
    !
    IF ( is_variable_in_output(var_name='pr_ice') ) THEN
       CALL add_var( this_list   = cloud_mig_list                                           ,&
            &        varname     = 'pr_ice'                                                 ,&
            &        ptr         = cloud_mig_output%pr_ice                                  ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_surface                                               ,&
            &        ldims       = shape2d                                                  ,&
            &        cf          = t_cf_var ('stratiform_icefall_flux',                      &
            &                                'kg/m2/s',                                      &
            &                                'stratiform icefall flux '//                    &
            &                                '(cloud_mig output)',                           &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,78,                                         &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_mig_output%pr_ice)
    END IF
    !
    IF ( is_variable_in_output(var_name='pr_snow') ) THEN
       CALL add_var( this_list   = cloud_mig_list                                           ,&
            &        varname     = 'pr_snow'                                                ,&
            &        ptr         = cloud_mig_output%pr_snow                                 ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_surface                                               ,&
            &        ldims       = shape2d                                                  ,&
            &        cf          = t_cf_var ('stratiform_snowfall_flux',                     &
            &                                'kg/m2/s',                                      &
            &                                'stratiform snowfall flux '//                   &
            &                                '(cloud_mig output)',                           &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,56,                                         &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_mig_output%pr_snow)
    END IF
    !
    IF ( is_variable_in_output(var_name='pr_grpl') ) THEN
       CALL add_var( this_list   = cloud_mig_list                                           ,&
            &        varname     = 'pr_grpl'                                                ,&
            &        ptr         = cloud_mig_output%pr_grpl                                 ,&
            &        hgrid       = grid_unstructured_cell                                   ,&
            &        vgrid       = za_surface                                               ,&
            &        ldims       = shape2d                                                  ,&
            &        cf          = t_cf_var ('stratiform_graupel_flux',                      &
            &                                'kg/m2/s',                                      &
            &                                'stratiform graupel flux '//                    &
            &                                '(cloud_mig output)',                           &
            &                                datatype_flt)                                  ,&
            &        grib2       = grib2_var(0,1,75,                                         &
            &                                datatype_grb,                                   &
            &                                grid_unstructured,                              &
            &                                grid_cell)                                     ,&
            &        isteptype   = tstep_instant                                            ,&
            &        lopenacc    =.TRUE.                                                    )
       __acc_attach(cloud_mig_output%pr_grpl)
    END IF

  END SUBROUTINE construct_cloud_mig_list
  !-------------

END MODULE mo_cloud_mig_memory
