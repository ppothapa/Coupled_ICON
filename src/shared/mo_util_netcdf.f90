!>
!! Contains routines for reading netcdf-Files of various shape.
!!
!! Contains routines for reading data from netcdf-Files of various shape.
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! Moved here from mo_ext_data by Daniel reinert, DWD (2012-02-01)
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
MODULE mo_util_netcdf

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: finish
  USE mo_communication,      ONLY: idx_no, blk_no
  USE mo_parallel_config,    ONLY: p_test_run
  USE mo_mpi,                ONLY: my_process_is_stdio, p_io, p_bcast, &
    &                              p_comm_work_test, p_comm_work

  IMPLICIT NONE

  ! required for testing/reading topography
  INCLUDE 'netcdf.inc'


  PRIVATE

  PUBLIC  :: read_netcdf_data
  PUBLIC  :: read_netcdf_lu
  PUBLIC  :: nf


  INTERFACE read_netcdf_data
    MODULE PROCEDURE read_netcdf_2d
    MODULE PROCEDURE read_netcdf_2d_int
    MODULE PROCEDURE read_netcdf_3d
    MODULE PROCEDURE read_netcdf_aero
    MODULE PROCEDURE read_netcdf_4d
  END INTERFACE read_netcdf_data


  CHARACTER(len=*), PARAMETER :: version = '$Id$'


CONTAINS

  !-------------------------------------------------------------------------
  !>
  !! Read dataset from netcdf file
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-07-14)
  !! Adapted for parallel runs by Rainer Johanni (2010-12-07)
  !!
  SUBROUTINE read_netcdf_2d (ncid, varname, glb_arr_len, loc_arr_len, glb_index, var_out)

    CHARACTER(len=*), INTENT(IN)  ::  &  !< Var name of field to be read
      &  varname

    INTEGER, INTENT(IN) :: ncid          !< id of netcdf file
    INTEGER, INTENT(IN) :: glb_arr_len   !< length of 1D field (global)
    INTEGER, INTENT(IN) :: loc_arr_len   !< length of 1D field (local)
    INTEGER, INTENT(IN) :: glb_index(:)  !< Index mapping local to global

    REAL(wp), INTENT(INOUT) :: &         !< output field
      &  var_out(:,:)

    INTEGER :: varid, mpi_comm, j, jl, jb
    REAL(wp):: z_dummy_array(glb_arr_len)!< local dummy array
  !-------------------------------------------------------------------------

    ! Get var ID
    IF( my_process_is_stdio()) CALL nf(nf_inq_varid(ncid, TRIM(varname), varid))

    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    ENDIF

    ! I/O PE reads and broadcasts data

    IF(my_process_is_stdio()) CALL nf(nf_get_var_double(ncid, varid, z_dummy_array(:)))
    CALL p_bcast(z_dummy_array, p_io, mpi_comm)

    var_out(:,:) = 0._wp

    ! Set var_out from global data
    DO j = 1, loc_arr_len

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      var_out(jl,jb) = z_dummy_array(glb_index(j))
    ENDDO

  END SUBROUTINE read_netcdf_2d



  !-------------------------------------------------------------------------
  !>
  !! Read dataset from netcdf file
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-07-14)
  !! Adapted for parallel runs by Rainer Johanni (2010-12-07)
  !!
  SUBROUTINE read_netcdf_2d_int (ncid, varname, glb_arr_len, loc_arr_len, glb_index, var_out)

    CHARACTER(len=*), INTENT(IN)  ::  &  !< Var name of field to be read
      &  varname

    INTEGER, INTENT(IN) :: ncid          !< id of netcdf file
    INTEGER, INTENT(IN) :: glb_arr_len   !< length of 1D field (global)
    INTEGER, INTENT(IN) :: loc_arr_len   !< length of 1D field (local)
    INTEGER, INTENT(IN) :: glb_index(:)  !< Index mapping local to global

    INTEGER, INTENT(INOUT) :: &          !< output field
      &  var_out(:,:)

    INTEGER :: varid, mpi_comm, j, jl, jb
    INTEGER :: z_dummy_array(glb_arr_len)!< local dummy array
  !-------------------------------------------------------------------------

    ! Get var ID
    IF( my_process_is_stdio()) CALL nf(nf_inq_varid(ncid, TRIM(varname), varid))

    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    ENDIF

    ! I/O PE reads and broadcasts data

    IF( my_process_is_stdio()) CALL nf(nf_get_var_int(ncid, varid, z_dummy_array(:)))
    CALL p_bcast(z_dummy_array, p_io, mpi_comm)

    var_out(:,:) = 0

    ! Set var_out from global data
    DO j = 1, loc_arr_len

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      var_out(jl,jb) = z_dummy_array(glb_index(j))
    ENDDO

  END SUBROUTINE read_netcdf_2d_int



 !-------------------------------------------------------------------------
  !>
  !! Read 3D (inlcuding height) dataset from netcdf file
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-07-14)
  !! Adapted for parallel runs by Rainer Johanni (2010-12-07)
  !! Adapted for 3 D by Guenther Zaengl, DWD (2011-07-11)
  !!
  SUBROUTINE read_netcdf_3d (ncid, varname, glb_arr_len, &
       &                     loc_arr_len, glb_index, &
       &                     nlevs,      var_out)

    CHARACTER(len=*), INTENT(IN)  ::  &  !< Var name of field to be read
      &  varname

    INTEGER, INTENT(IN) :: ncid          !< id of netcdf file
    INTEGER, INTENT(IN) :: nlevs         !< vertical levels of netcdf file
    INTEGER, INTENT(IN) :: glb_arr_len   !< length of 1D field (global)
    INTEGER, INTENT(IN) :: loc_arr_len   !< length of 1D field (local)
    INTEGER, INTENT(IN) :: glb_index(:)  !< Index mapping local to global

    REAL(wp), INTENT(INOUT) :: &         !< output field
      &  var_out(:,:,:)

    INTEGER :: varid, mpi_comm, j, jl, jb, jk
    REAL(wp):: z_dummy_array(glb_arr_len,nlevs)!< local dummy array
  !-------------------------------------------------------------------------

    ! Get var ID
    IF(my_process_is_stdio()) CALL nf(nf_inq_varid(ncid, TRIM(varname), varid))

    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    ENDIF

    ! I/O PE reads and broadcasts data

    IF(my_process_is_stdio()) CALL nf(nf_get_var_double(ncid, varid, z_dummy_array(:,:)))
    CALL p_bcast(z_dummy_array, p_io, mpi_comm)

    var_out(:,:,:) = 0._wp

    ! Set var_out from global data
     DO jk = 1, nlevs
       DO j = 1, loc_arr_len

         jb = blk_no(j) ! Block index in distributed patch
         jl = idx_no(j) ! Line  index in distributed patch
               
         var_out(jl,jk,jb) = z_dummy_array(glb_index(j),jk)

       ENDDO
     ENDDO

  END SUBROUTINE read_netcdf_3d



  !-------------------------------------------------------------------------

  SUBROUTINE read_netcdf_aero (ncid, ntime, varname, glb_arr_len, &
       &                     loc_arr_len, glb_index, var_out)

    CHARACTER(len=*), INTENT(IN)  ::  &  !< Var name of field to be read
      &  varname

    INTEGER, INTENT(IN) :: ncid          !< id of netcdf file
    INTEGER, INTENT(IN) :: ntime         !< time levels of netcdf file
    INTEGER, INTENT(IN) :: glb_arr_len   !< length of 1D field (global)
    INTEGER, INTENT(IN) :: loc_arr_len   !< length of 1D field (local)
    INTEGER, INTENT(IN) :: glb_index(:)  !< Index mapping local to global

    REAL(wp), INTENT(INOUT) :: &         !< output field
      &  var_out(:,:,:)

    INTEGER :: varid, mpi_comm, j, jl, jb, jt
    REAL(wp):: z_dummy_array(glb_arr_len,ntime)!< local dummy array
  !-------------------------------------------------------------------------

    ! Get var ID
    IF(my_process_is_stdio()) CALL nf(nf_inq_varid(ncid, TRIM(varname), varid))

    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    ENDIF


    ! I/O PE reads and broadcasts data

    IF(my_process_is_stdio()) CALL nf(nf_get_var_double(ncid, varid, z_dummy_array(:,:)))
    CALL p_bcast(z_dummy_array, p_io , mpi_comm)

    var_out(:,:,:) = 0._wp

    ! Set var_out from global data
    DO jt = 1, ntime
           DO j = 1, loc_arr_len

             jb = blk_no(j) ! Block index in distributed patch
             jl = idx_no(j) ! Line  index in distributed patch
               
             var_out(jl,jb,jt) = z_dummy_array(glb_index(j),jt)

          ENDDO
    ENDDO

  END SUBROUTINE read_netcdf_aero




  !-------------------------------------------------------------------------
  ! Specific read-routine for LU_CLASS_FRACTION. Probably, read_netcdf_aero 
  ! and read_netcdf_lu can be merged into a single routine in the near 
  ! future.

  SUBROUTINE read_netcdf_lu (ncid, varname, glb_arr_len,            &
       &                     loc_arr_len, glb_index, nslice, var_out)

    CHARACTER(len=*), INTENT(IN)  ::  &  !< Var name of field to be read
      &  varname

    INTEGER, INTENT(IN) :: ncid          !< id of netcdf file
    INTEGER, INTENT(IN) :: nslice        !< slices o netcdf field
    INTEGER, INTENT(IN) :: glb_arr_len   !< length of 1D field (global)
    INTEGER, INTENT(IN) :: loc_arr_len   !< length of 1D field (local)
    INTEGER, INTENT(IN) :: glb_index(:)  !< Index mapping local to global

    REAL(wp), INTENT(INOUT) :: &         !< output field
      &  var_out(:,:,:)

    INTEGER :: varid, mpi_comm, j, jl, jb, js
    REAL(wp):: z_dummy_array(glb_arr_len, nslice)!< local dummy array

  !-------------------------------------------------------------------------


    ! Get var ID
    IF(my_process_is_stdio()) CALL nf(nf_inq_varid(ncid, TRIM(varname), varid))

    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    ENDIF



    ! I/O PE reads and broadcasts data

    IF(my_process_is_stdio()) CALL nf(nf_get_var_double(ncid, varid, z_dummy_array(:,:)))
    CALL p_bcast(z_dummy_array, p_io , mpi_comm)

    var_out(:,:,:) = 0._wp


    ! Set var_out from global data
    DO js = 1, nslice
      DO j = 1, loc_arr_len

        jb = blk_no(j) ! Block index in distributed patch
        jl = idx_no(j) ! Line  index in distributed patch
               
        var_out(jl,jb,js) = z_dummy_array(glb_index(j),js)

      ENDDO
    ENDDO

  END SUBROUTINE read_netcdf_lu


 !-------------------------------------------------------------------------
  !>
  !! Read 4D (inlcuding height and time) dataset from netcdf file
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-07-14)
  !! Adapted for parallel runs by Rainer Johanni (2010-12-07)
  !! Adapted for 4 D by Kristina Froehlich, MPI-M (2011-06-16)
  !!
  SUBROUTINE read_netcdf_4d (ncid, varname, glb_arr_len, &
       &                     loc_arr_len, glb_index, &
       &                     nlevs, ntime,      var_out)

    CHARACTER(len=*), INTENT(IN)  ::  &  !< Var name of field to be read
      &  varname

    INTEGER, INTENT(IN) :: ncid          !< id of netcdf file
    INTEGER, INTENT(IN) :: nlevs         !< vertical levels of netcdf file
    INTEGER, INTENT(IN) :: ntime         !< time levels of netcdf file
    INTEGER, INTENT(IN) :: glb_arr_len   !< length of 1D field (global)
    INTEGER, INTENT(IN) :: loc_arr_len   !< length of 1D field (local)
    INTEGER, INTENT(IN) :: glb_index(:)  !< Index mapping local to global

    REAL(wp), INTENT(INOUT) :: &         !< output field
      &  var_out(:,:,:,:)

    INTEGER :: varid, mpi_comm, j, jl, jb, jk, jt
    REAL(wp):: z_dummy_array(glb_arr_len,nlevs,ntime)!< local dummy array
  !-------------------------------------------------------------------------

    ! Get var ID
    IF(my_process_is_stdio()) CALL nf(nf_inq_varid(ncid, TRIM(varname), varid))

    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    ENDIF

    ! I/O PE reads and broadcasts data

    IF(my_process_is_stdio()) CALL nf(nf_get_var_double(ncid, varid, z_dummy_array(:,:,:)))
    CALL p_bcast(z_dummy_array, p_io , mpi_comm)

    var_out(:,:,:,:) = 0._wp

    ! Set var_out from global data
    DO jt = 1, ntime
        DO jk = 1, nlevs
           DO j = 1, loc_arr_len

             jb = blk_no(j) ! Block index in distributed patch
             jl = idx_no(j) ! Line  index in distributed patch
               
             var_out(jl,jk,jb,jt) = z_dummy_array(glb_index(j),jk,jt)

          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE read_netcdf_4d


  !-------------------------------------------------------------------------


  SUBROUTINE nf(status)

    INTEGER, INTENT(in) :: status

    IF (status /= nf_noerr) THEN
      CALL finish('mo_ext_data netCDF error', nf_strerror(status))
    ENDIF

  END SUBROUTINE nf


END MODULE mo_util_netcdf
