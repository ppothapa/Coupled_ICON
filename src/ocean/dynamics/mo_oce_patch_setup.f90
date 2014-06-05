!>
!!               The module <i>mo_model_import_domain</i>.
!!
!!               The module <i>mo_model_import_domain</i>
!! provides functionality to import information about the models computational
!! domain. This information is read from several files that were generated by
!! the patch_2D generator programm. The data types describing the model domain are
!! contained in <i>mo_domain_model</i>.
!!
!! @par Revision History
!! Initial version  by: Peter Korn,  MPI-M, Hamburg, June 2005
!! Modification by Thomas Heinze (2006-02-21):
!! - renamed m_modules to mo_modules
!! Modification by Thomas Heinze (2006-09-20):
!! - added routine grid_and_patch_diagnosis
!! Modification by Pilar Ripodas, DWD, (2007-01-31)
!! - addapted to the new TYPE grid_edges (system_orientation added)
!! Modification by Peter Korn,  MPI-M, (2006-12)
!! - implementation of topography and boundary treatment, i.e.
!!   initialization of the grid & patch_2D components that carry
!!   information about topography and the lateral boundaries of
!!   the domain; this is not related to patch_2D boundaries.
!!   topography can either be computed by analytical l,eans or
!!   by reading from database files.
!! Modification by Hui Wan, MPI-M, (2007-02-23)
!! - Subroutine <i>init_import</i> was changed to <i>setup_grid</i>.
!!   Namelist hierarchy_ini was renamed to <i>grid_ctl</i>, moved
!!   from <i>mo_io_utilities</i> to this module and now read from
!!   an external file in subroutine <i>setup_grid</i>.
!! - Some changes in <i>init_ocean_patch_component</i> after
!!   discussion with Peter.
!! - Calculation of the min. primal edge length was added to
!!   <i>import_patches</i>. However, shouldn't it be an array with
!!   one element for each patch_2D, rather than a scalar?
!! Modification by P. Ripodas, DWD, (2007-03-14):
!! - Now the output of "import_patches" is the min_dual_edge_lenght
!!   instead of the min_primal_edge_lenght. It will be used to set
!!   the horizontal diffusion parameter. Now it is done as it was
!!   in the prototype.
!! Modification by Almut Gassmann, MPI-M (2007-04)
!! - removed loptimize to make compatible with new grid generator
!! - removed itoa for good programming style
!! - reorganized patch_2D input to be compatible with the new patch_2D generator
!! - cleaning up "destruct_patches"
!! Modification by Almut Gassmann, MPI-M (2007-04-13)
!! - remove grid type and perform related adaptations
!!   (grid information comes now inside a patch_2D)
!! - changed subroutine name form setup_grid to setup_files
!! Modified by Hui Wan, MPI-M, (2008-04-04)
!!  - topography_file_dir renamed topo_file_dir
!!  - for the hydro_atmos, control variable testtype renamed ctest_name.
!! Modified by Almut Gassmann, MPI-M, (2008-04-23)
!!  - itopo distinguishes now shallow water (itopo=1) orography function
!!    from hydro_atmos orography function (itopo=2)
!! Modification by Jochen Foerstner, DWD, (2008-07-16)
!!  - new fields in the derived type for the edges:
!!    grid_edges%primal_cart_normal (Cartesian normal to edge),
!!    grid_edges%quad_idx, grid_edges%quad_area and grid_edges%quad_orientation
!!    (indices of edges and area of the quadrilateral formed by two adjacent cells)
!!    up to now these new fields are initialized in the new routines
!!    calculate_primal_cart_normal and init_quad_twoadjcells
!!    rather than read from a grid/patch_2D file.
!! Modification by Almut Gassmann, MPI-M, (2008-09-21)
!!  - remove reference to mask and height files, they are never used
!!  - use cell_type to distinguish cells as triangles or hexagons
!! Modification by Almut Gassmann, MPI-M (2008-10-30)
!!  - add subroutine init_coriolis to initialize Coriolis parameter
!! Modification by Constantin Junk, MPI-M (2011-04-05)
!! - ...
!! Modification by Daniel Reinert, DWD (2012-04-11)
!!  - new routine which initializes the butterfly data structure
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

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_oce_patch_setup
  !-------------------------------------------------------------------------
  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: finish, warning, message
  USE mo_model_domain,       ONLY: t_patch
  USE mo_parallel_config,    ONLY: nproma
  USE mo_grid_config,        ONLY: corio_lat, grid_angular_velocity, use_dummy_cell_closure, &
    & grid_sphere_radius, lplane
  USE mo_sync,               ONLY: sync_c, sync_e, sync_v, sync_patch_array, sync_idx
  USE mo_grid_subset,        ONLY: fill_subset,t_subset_range, get_index_range, read_subset, write_subset
  USE mo_mpi,                ONLY: work_mpi_barrier, get_my_mpi_work_id, my_process_is_mpi_seq, global_mpi_barrier, &
    & get_my_global_mpi_id

  USE mo_loopindices
  USE mo_impl_constants
  USE mo_math_types
  USE mo_math_utilities
  USE mo_grid_geometry_info, ONLY: planar_torus_geometry
  USE mo_master_control,     ONLY: get_my_process_type, ocean_process
  USE mo_model_domimp_setup, ONLY: init_coriolis
  USE mo_dynamics_config,    ONLY: lcoriolis
  USE mo_sync,               ONLY: disable_sync_checks, enable_sync_checks

  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC :: complete_ocean_patch
  
 
  
  !-------------------------------------------------------------------------
  
CONTAINS
    
  !----------------------------------------------------------------------------
  !>
!<Optimize:inUse>
  SUBROUTINE complete_ocean_patch( patch_2D)
    TYPE(t_patch), TARGET, INTENT(inout) :: patch_2D

    CALL init_coriolis( lcoriolis, lplane, patch_2D )
    CALL complete_ocean_patch_geometry( patch_2D )

  END SUBROUTINE complete_ocean_patch
  !----------------------------------------------------------------------------


  !----------------------------------------------------------------------------
  !>
  !! Computes the local orientation of the edge primal normal and dual normal.
  !!
  !! Computes the local orientation of the edge primal normal and dual normal
  !! at the location of the cell centers and vertices.
  !! Moreover, the Cartesian orientation vectors of the edge primal normals
  !! are stored for use in the RBF initialization routines, and inverse
  !! primal and dual edge lengths are computed
  !!
  !! Note: Not clear if all the included calclulations are needed
  !!
  !! @par Revision History
  !!  developed by Guenther Zaengl, 2009-03-31
  !!
!<Optimize:inUse>
  SUBROUTINE complete_ocean_patch_geometry( patch_2D)
    TYPE(t_patch), TARGET, INTENT(inout) :: patch_2D

    INTEGER :: jb, je!, jc
    INTEGER :: start_idx, end_idx
    TYPE(t_subset_range), POINTER :: all_edges, owned_edges

    INTEGER :: ilc1, ibc1, ilv1, ibv1, ilc2, ibc2, ilv2, ibv2, &
      & ilv3, ibv3, ilv4, ibv4!, ile1, ibe1

    REAL(wp) :: z_nu, z_nv, z_lon, z_lat, z_nx1(3), z_nx2(3), z_norm

    TYPE(t_cartesian_coordinates) :: cc_ev3, cc_ev4

    !-----------------------------------------------------------------------
    all_edges    => patch_2D%edges%ALL
    owned_edges  => patch_2D%edges%owned

    ! !$OMP PARALLEL  PRIVATE(rl_start,rl_end,i_startblk,i_endblk)
    !
    ! First step: compute Cartesian coordinates and Cartesian vectors on full domain
    ! this is needed to vectorize RBF initialization; the existing field carrying
    ! the Cartesian orientation vectors (primal_cart_normal) did not work for that
    ! because it is a derived data type
    ! In addition, the fields for the inverse primal and dual edge lengths are
    ! initialized here.
    !
    ! !$OMP DO PRIVATE(jb,start_idx,end_idx,je) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, jb, start_idx, end_idx)
      DO je = start_idx, end_idx

        ! compute Cartesian coordinates (needed for RBF initialization)
        patch_2D%edges%inv_primal_edge_length(je,jb) = &
          & 1._wp/patch_2D%edges%primal_edge_length(je,jb)

        ! compute inverse dual edge length (undefined for refin_ctrl=1)
        patch_2D%edges%inv_dual_edge_length(je,jb) = &
          & 1._wp/patch_2D%edges%dual_edge_length(je,jb)

      ENDDO

    END DO !block loop
    ! !$OMP END DO

    ! Second step: computed projected orientation vectors and related information

    ! Initialization of lateral boundary points
    ! !$OMP WORKSHARE
    patch_2D%edges%vertex_idx(:,:,3)            = 0
    patch_2D%edges%vertex_blk(:,:,3)            = 0
    patch_2D%edges%vertex_idx(:,:,4)            = 0
    patch_2D%edges%vertex_blk(:,:,4)            = 0
    patch_2D%edges%primal_normal_cell(:,:,:)%v1 = 0._wp
    patch_2D%edges%dual_normal_cell  (:,:,:)%v1 = 0._wp
    patch_2D%edges%primal_normal_vert(:,:,:)%v1 = 0._wp
    patch_2D%edges%dual_normal_vert  (:,:,:)%v1 = 0._wp
    patch_2D%edges%primal_normal_cell(:,:,:)%v2 = 0._wp
    patch_2D%edges%dual_normal_cell  (:,:,:)%v2 = 0._wp
    patch_2D%edges%primal_normal_vert(:,:,:)%v2 = 0._wp
    patch_2D%edges%dual_normal_vert  (:,:,:)%v2 = 0._wp
    ! !$OMP END WORKSHARE
    !
    ! loop through all patch_2D edges
    !
    ! !$OMP DO PRIVATE(jb,start_idx,end_idx,je,ilc1,ibc1,ilv1,ibv1,ilc2,ibc2,ilv2, &
    ! !$OMP            ibv2,ilv3,ibv3,ilv4,ibv4,z_nu,z_nv,z_lon,z_lat,z_nx1,z_nx2,   &
    ! !$OMP            cc_ev3,cc_ev4,z_norm) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = owned_edges%start_block, owned_edges%end_block
      CALL get_index_range(owned_edges, jb, start_idx, end_idx)
      DO je = start_idx, end_idx

        IF(.NOT.patch_2D%edges%decomp_info%owner_mask(je,jb)) CYCLE

        ! compute edge-vertex indices (and blocks) 3 and 4, which
        ! are the outer vertices of cells 1 and 2, respectively,
        ! and the inverse length bewtween vertices 3 and 4

        ilc1 = patch_2D%edges%cell_idx(je,jb,1)
        ibc1 = patch_2D%edges%cell_blk(je,jb,1)
        ilc2 = patch_2D%edges%cell_idx(je,jb,2)
        ibc2 = patch_2D%edges%cell_blk(je,jb,2)

        ilv1 = patch_2D%edges%vertex_idx(je,jb,1)
        ibv1 = patch_2D%edges%vertex_blk(je,jb,1)
        ilv2 = patch_2D%edges%vertex_idx(je,jb,2)
        ibv2 = patch_2D%edges%vertex_blk(je,jb,2)

        IF( ilc1 < 1 .or. ilc2 < 1) CYCLE
        IF( ilv1 < 1 .or. ilv2 < 1) &
          CALL finish("complete_ocean_patch_geometry","ilv1 < 1 .or. ilv2 < 1")

        IF ((patch_2D%cells%vertex_idx(ilc1,ibc1,1) /= &
          & patch_2D%edges%vertex_idx(je,jb,1) .OR.  &
          & patch_2D%cells%vertex_blk(ilc1,ibc1,1) /= &
          & patch_2D%edges%vertex_blk(je,jb,1)) .AND.  &
          & (patch_2D%cells%vertex_idx(ilc1,ibc1,1) /= &
          & patch_2D%edges%vertex_idx(je,jb,2) .OR.  &
          & patch_2D%cells%vertex_blk(ilc1,ibc1,1) /= &
          & patch_2D%edges%vertex_blk(je,jb,2)) )        THEN

          patch_2D%edges%vertex_idx(je,jb,3) = patch_2D%cells%vertex_idx(ilc1,ibc1,1)
          patch_2D%edges%vertex_blk(je,jb,3) = patch_2D%cells%vertex_blk(ilc1,ibc1,1)

        ELSE IF ((patch_2D%cells%vertex_idx(ilc1,ibc1,2) /= &
          & patch_2D%edges%vertex_idx(je,jb,1) .OR.  &
          & patch_2D%cells%vertex_blk(ilc1,ibc1,2) /= &
          & patch_2D%edges%vertex_blk(je,jb,1)) .AND.  &
          & (patch_2D%cells%vertex_idx(ilc1,ibc1,2) /= &
          & patch_2D%edges%vertex_idx(je,jb,2) .OR.  &
          & patch_2D%cells%vertex_blk(ilc1,ibc1,2) /= &
          & patch_2D%edges%vertex_blk(je,jb,2)) )        THEN

          patch_2D%edges%vertex_idx(je,jb,3) = patch_2D%cells%vertex_idx(ilc1,ibc1,2)
          patch_2D%edges%vertex_blk(je,jb,3) = patch_2D%cells%vertex_blk(ilc1,ibc1,2)

        ELSE IF ((patch_2D%cells%vertex_idx(ilc1,ibc1,3) /= &
          & patch_2D%edges%vertex_idx(je,jb,1) .OR.  &
          & patch_2D%cells%vertex_blk(ilc1,ibc1,3) /= &
          & patch_2D%edges%vertex_blk(je,jb,1)) .AND.  &
          & (patch_2D%cells%vertex_idx(ilc1,ibc1,3) /= &
          & patch_2D%edges%vertex_idx(je,jb,2) .OR.  &
          & patch_2D%cells%vertex_blk(ilc1,ibc1,3) /= &
          & patch_2D%edges%vertex_blk(je,jb,2)) )        THEN

          patch_2D%edges%vertex_idx(je,jb,3) = patch_2D%cells%vertex_idx(ilc1,ibc1,3)
          patch_2D%edges%vertex_blk(je,jb,3) = patch_2D%cells%vertex_blk(ilc1,ibc1,3)

        ENDIF

        IF ((patch_2D%cells%vertex_idx(ilc2,ibc2,1) /= &
          & patch_2D%edges%vertex_idx(je,jb,1) .OR.  &
          & patch_2D%cells%vertex_blk(ilc2,ibc2,1) /= &
          & patch_2D%edges%vertex_blk(je,jb,1)) .AND.  &
          & (patch_2D%cells%vertex_idx(ilc2,ibc2,1) /= &
          & patch_2D%edges%vertex_idx(je,jb,2) .OR.  &
          & patch_2D%cells%vertex_blk(ilc2,ibc2,1) /= &
          & patch_2D%edges%vertex_blk(je,jb,2)) )        THEN

          patch_2D%edges%vertex_idx(je,jb,4) = patch_2D%cells%vertex_idx(ilc2,ibc2,1)
          patch_2D%edges%vertex_blk(je,jb,4) = patch_2D%cells%vertex_blk(ilc2,ibc2,1)

        ELSE IF ((patch_2D%cells%vertex_idx(ilc2,ibc2,2) /= &
          & patch_2D%edges%vertex_idx(je,jb,1) .OR.  &
          & patch_2D%cells%vertex_blk(ilc2,ibc2,2) /= &
          & patch_2D%edges%vertex_blk(je,jb,1)) .AND.  &
          & (patch_2D%cells%vertex_idx(ilc2,ibc2,2) /= &
          & patch_2D%edges%vertex_idx(je,jb,2) .OR.  &
          & patch_2D%cells%vertex_blk(ilc2,ibc2,2) /= &
          & patch_2D%edges%vertex_blk(je,jb,2)) )        THEN

          patch_2D%edges%vertex_idx(je,jb,4) = patch_2D%cells%vertex_idx(ilc2,ibc2,2)
          patch_2D%edges%vertex_blk(je,jb,4) = patch_2D%cells%vertex_blk(ilc2,ibc2,2)

        ELSE IF ((patch_2D%cells%vertex_idx(ilc2,ibc2,3) /= &
          & patch_2D%edges%vertex_idx(je,jb,1) .OR.  &
          & patch_2D%cells%vertex_blk(ilc2,ibc2,3) /= &
          & patch_2D%edges%vertex_blk(je,jb,1)) .AND.  &
          & (patch_2D%cells%vertex_idx(ilc2,ibc2,3) /= &
          & patch_2D%edges%vertex_idx(je,jb,2) .OR.  &
          & patch_2D%cells%vertex_blk(ilc2,ibc2,3) /= &
          & patch_2D%edges%vertex_blk(je,jb,2)) )        THEN

          patch_2D%edges%vertex_idx(je,jb,4) = patch_2D%cells%vertex_idx(ilc2,ibc2,3)
          patch_2D%edges%vertex_blk(je,jb,4) = patch_2D%cells%vertex_blk(ilc2,ibc2,3)

        ENDIF

        ilv3 = patch_2D%edges%vertex_idx(je,jb,3)
        ibv3 = patch_2D%edges%vertex_blk(je,jb,3)
        ilv4 = patch_2D%edges%vertex_idx(je,jb,4)
        ibv4 = patch_2D%edges%vertex_blk(je,jb,4)

        cc_ev3 = patch_2D%verts%cartesian(ilv3,ibv3)
        cc_ev4 = patch_2D%verts%cartesian(ilv4,ibv4)

        ! inverse length bewtween vertices 3 and 4
        patch_2D%edges%inv_vert_vert_length(je,jb) = 1._wp / &
          & (grid_sphere_radius * arc_length(cc_ev3,cc_ev4))
        !      ENDIF

        ! next step: compute projected orientation vectors for cells and vertices
        ! bordering to each edge (incl. vertices 3 and 4 intorduced above)

        ! transform primal normal to cartesian vector z_nx1
        z_nx1(:) = patch_2D%edges%primal_cart_normal(je,jb)%x(:)
        ! get location of cell 1
        z_lon = patch_2D%cells%center(ilc1,ibc1)%lon
        z_lat = patch_2D%cells%center(ilc1,ibc1)%lat
        ! compute local primal at cell 1
        CALL cvec2gvec(z_nx1(1),z_nx1(2),z_nx1(3),z_lon,z_lat,z_nu,z_nv)
        z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)
        patch_2D%edges%primal_normal_cell(je,jb,1)%v1 = z_nu/z_norm
        patch_2D%edges%primal_normal_cell(je,jb,1)%v2 = z_nv/z_norm

        ! transform dual normal to cartesian vector z_nx2
        z_nx2(:) = patch_2D%edges%dual_cart_normal(je,jb)%x(:)
        ! compute local dual normals at cell 1
        CALL cvec2gvec(z_nx2(1),z_nx2(2),z_nx2(3),z_lon,z_lat,z_nu,z_nv)
        z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

        patch_2D%edges%dual_normal_cell(je,jb,1)%v1 = z_nu/z_norm
        patch_2D%edges%dual_normal_cell(je,jb,1)%v2 = z_nv/z_norm

        ! get location of cell 2
        z_lon = patch_2D%cells%center(ilc2,ibc2)%lon
        z_lat = patch_2D%cells%center(ilc2,ibc2)%lat

        ! compute local primal and dual normals at cell 2
        CALL cvec2gvec(z_nx1(1),z_nx1(2),z_nx1(3),z_lon,z_lat,z_nu,z_nv)
        z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

        patch_2D%edges%primal_normal_cell(je,jb,2)%v1 = z_nu/z_norm
        patch_2D%edges%primal_normal_cell(je,jb,2)%v2 = z_nv/z_norm

        CALL cvec2gvec(z_nx2(1),z_nx2(2),z_nx2(3),z_lon,z_lat,z_nu,z_nv)
        z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

        patch_2D%edges%dual_normal_cell(je,jb,2)%v1 = z_nu/z_norm
        patch_2D%edges%dual_normal_cell(je,jb,2)%v2 = z_nv/z_norm

        ! get location of vertex 1

        z_lon = patch_2D%verts%vertex(ilv1,ibv1)%lon
        z_lat = patch_2D%verts%vertex(ilv1,ibv1)%lat

        ! compute local primal and dual normals at vertex 1

        CALL cvec2gvec(z_nx1(1),z_nx1(2),z_nx1(3),z_lon,z_lat,z_nu,z_nv)
        z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

        patch_2D%edges%primal_normal_vert(je,jb,1)%v1 = z_nu/z_norm
        patch_2D%edges%primal_normal_vert(je,jb,1)%v2 = z_nv/z_norm

        CALL cvec2gvec(z_nx2(1),z_nx2(2),z_nx2(3),z_lon,z_lat,z_nu,z_nv)
        z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

        patch_2D%edges%dual_normal_vert(je,jb,1)%v1 = z_nu/z_norm
        patch_2D%edges%dual_normal_vert(je,jb,1)%v2 = z_nv/z_norm

        ! get location of vertex 2

        z_lon = patch_2D%verts%vertex(ilv2,ibv2)%lon
        z_lat = patch_2D%verts%vertex(ilv2,ibv2)%lat

        ! compute local primal and dual normals at vertex 2
        CALL cvec2gvec(z_nx1(1),z_nx1(2),z_nx1(3),z_lon,z_lat,z_nu,z_nv)
        z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

        patch_2D%edges%primal_normal_vert(je,jb,2)%v1 = z_nu/z_norm
        patch_2D%edges%primal_normal_vert(je,jb,2)%v2 = z_nv/z_norm

        CALL cvec2gvec(z_nx2(1),z_nx2(2),z_nx2(3),z_lon,z_lat,z_nu,z_nv)
        z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

        patch_2D%edges%dual_normal_vert(je,jb,2)%v1 = z_nu/z_norm
        patch_2D%edges%dual_normal_vert(je,jb,2)%v2 = z_nv/z_norm

        ! get location of vertex 3
        z_lon = patch_2D%verts%vertex(ilv3,ibv3)%lon
        z_lat = patch_2D%verts%vertex(ilv3,ibv3)%lat

        ! compute local primal and dual normals at vertex 3
        CALL cvec2gvec(z_nx1(1),z_nx1(2),z_nx1(3),z_lon,z_lat,z_nu,z_nv)
        z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

        patch_2D%edges%primal_normal_vert(je,jb,3)%v1 = z_nu/z_norm
        patch_2D%edges%primal_normal_vert(je,jb,3)%v2 = z_nv/z_norm

        CALL cvec2gvec(z_nx2(1),z_nx2(2),z_nx2(3),z_lon,z_lat,z_nu,z_nv)
        z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

        patch_2D%edges%dual_normal_vert(je,jb,3)%v1 = z_nu/z_norm
        patch_2D%edges%dual_normal_vert(je,jb,3)%v2 = z_nv/z_norm

        ! get location of vertex 4

        z_lon = patch_2D%verts%vertex(ilv4,ibv4)%lon
        z_lat = patch_2D%verts%vertex(ilv4,ibv4)%lat

        ! compute local primal and dual normals at vertex 2
        CALL cvec2gvec(z_nx1(1),z_nx1(2),z_nx1(3),z_lon,z_lat,z_nu,z_nv)
        z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

        patch_2D%edges%primal_normal_vert(je,jb,4)%v1 = z_nu/z_norm
        patch_2D%edges%primal_normal_vert(je,jb,4)%v2 = z_nv/z_norm

        CALL cvec2gvec(z_nx2(1),z_nx2(2),z_nx2(3),z_lon,z_lat,z_nu,z_nv)
        z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

        patch_2D%edges%dual_normal_vert(je,jb,4)%v1 = z_nu/z_norm
        patch_2D%edges%dual_normal_vert(je,jb,4)%v2 = z_nv/z_norm

      ENDDO

    END DO !block loop
    ! !$OMP END DO NOWAIT

    ! !$OMP END PARALLEL

    ! primal_normal_cell must be sync'd before next loop,
    ! so do a sync for all above calculated quantities

    CALL disable_sync_checks
    CALL sync_idx(sync_e,sync_v,patch_2D,patch_2D%edges%vertex_idx(:,:,3), &
      & patch_2D%edges%vertex_blk(:,:,3))
    CALL sync_idx(sync_e,sync_v,patch_2D,patch_2D%edges%vertex_idx(:,:,4), &
      & patch_2D%edges%vertex_blk(:,:,4))

    CALL sync_patch_array(sync_e,patch_2D,patch_2D%edges%inv_vert_vert_length)

    CALL sync_patch_array(sync_e,patch_2D,patch_2D%edges%primal_normal_cell(:,:,1)%v1)
    CALL sync_patch_array(sync_e,patch_2D,patch_2D%edges%primal_normal_cell(:,:,2)%v1)
    CALL sync_patch_array(sync_e,patch_2D,patch_2D%edges%primal_normal_vert(:,:,1)%v1)
    CALL sync_patch_array(sync_e,patch_2D,patch_2D%edges%primal_normal_vert(:,:,2)%v1)
    CALL sync_patch_array(sync_e,patch_2D,patch_2D%edges%primal_normal_vert(:,:,3)%v1)
    CALL sync_patch_array(sync_e,patch_2D,patch_2D%edges%primal_normal_vert(:,:,4)%v1)

    CALL sync_patch_array(sync_e,patch_2D,patch_2D%edges%dual_normal_cell(:,:,1)%v1)
    CALL sync_patch_array(sync_e,patch_2D,patch_2D%edges%dual_normal_cell(:,:,2)%v1)
    CALL sync_patch_array(sync_e,patch_2D,patch_2D%edges%dual_normal_vert(:,:,1)%v1)
    CALL sync_patch_array(sync_e,patch_2D,patch_2D%edges%dual_normal_vert(:,:,2)%v1)
    CALL sync_patch_array(sync_e,patch_2D,patch_2D%edges%dual_normal_vert(:,:,3)%v1)
    CALL sync_patch_array(sync_e,patch_2D,patch_2D%edges%dual_normal_vert(:,:,4)%v1)

    CALL sync_patch_array(sync_e,patch_2D,patch_2D%edges%primal_normal_cell(:,:,1)%v2)
    CALL sync_patch_array(sync_e,patch_2D,patch_2D%edges%primal_normal_cell(:,:,2)%v2)
    CALL sync_patch_array(sync_e,patch_2D,patch_2D%edges%primal_normal_vert(:,:,1)%v2)
    CALL sync_patch_array(sync_e,patch_2D,patch_2D%edges%primal_normal_vert(:,:,2)%v2)
    CALL sync_patch_array(sync_e,patch_2D,patch_2D%edges%primal_normal_vert(:,:,3)%v2)
    CALL sync_patch_array(sync_e,patch_2D,patch_2D%edges%primal_normal_vert(:,:,4)%v2)

    CALL sync_patch_array(sync_e,patch_2D,patch_2D%edges%dual_normal_cell(:,:,1)%v2)
    CALL sync_patch_array(sync_e,patch_2D,patch_2D%edges%dual_normal_cell(:,:,2)%v2)
    CALL sync_patch_array(sync_e,patch_2D,patch_2D%edges%dual_normal_vert(:,:,1)%v2)
    CALL sync_patch_array(sync_e,patch_2D,patch_2D%edges%dual_normal_vert(:,:,2)%v2)
    CALL sync_patch_array(sync_e,patch_2D,patch_2D%edges%dual_normal_vert(:,:,3)%v2)
    CALL sync_patch_array(sync_e,patch_2D,patch_2D%edges%dual_normal_vert(:,:,4)%v2)

    CALL enable_sync_checks


    !!$OMP PARALLEL  PRIVATE(rl_start,rl_end,i_startblk,i_endblk)


    !!$OMP END PARALLEL
  END SUBROUTINE complete_ocean_patch_geometry
  !-------------------------------------------------------------------------
    
    
END MODULE mo_oce_patch_setup
