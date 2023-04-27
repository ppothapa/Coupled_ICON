!>
!! @brief The subroutine omp_block_loop_cell provides an OpenMP parallelized
!! loop over all data blocks on the patch. For each data block, indices are
!! computed for cells in the data block belonging to the domain, and the
!! routine passed as argument is called with these indices to work on these
!! cells.
!!
!! @author Marco Giorgetta, MPI-M
!!
!! @par Revision History
!!  Original version 2017-12
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

MODULE mo_omp_block_loop

  USE mo_model_domain       ,ONLY: t_patch
  USE mo_loopindices        ,ONLY: get_indices_c
  USE mo_impl_constants     ,ONLY: min_rlcell_int
  USE mo_impl_constants_grf ,ONLY: grf_bdywidth_c

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: omp_block_loop_cell

CONTAINS

  SUBROUTINE omp_block_loop_cell(patch, routine)

    ! Arguments
    !
    TYPE(t_patch), TARGET, INTENT(in) :: patch
    !
    INTERFACE
       !
       SUBROUTINE routine(jg, jb, jcs, jce)
         !
         INTEGER, INTENT(in) :: jg, jb, jcs, jce
         !
       END SUBROUTINE routine
       !
    END INTERFACE

    ! Local variables
    !
    INTEGER  :: jg          !< grid index
    INTEGER  :: ncd         !< number of child domains of grid jg (but at least 1)
    INTEGER  :: rls, rle    !< refinement level, start and end
    INTEGER  :: jb          !< index of block loop
    INTEGER  :: jbs, jbe    !< start and end indices of block   loop
    INTEGER  :: jcs, jce    !< start and end indices of columns loop

    jg  = patch%id

    rls = grf_bdywidth_c+1
    rle = min_rlcell_int

    ncd = MAX(1,patch%n_childdom)
    jbs = patch%cells%start_blk(rls,   1)
    jbe = patch%cells%  end_blk(rle, ncd)

!$OMP PARALLEL DO PRIVATE(jb, jcs, jce)
    DO jb = jbs, jbe
       !
       CALL get_indices_c(patch, jb, jbs, jbe, jcs, jce, rls, rle)
       IF (jcs>jce) CYCLE
       !
       CALL routine(jg, jb, jcs, jce)
       !
    END DO
!$OMP END PARALLEL DO

  END SUBROUTINE omp_block_loop_cell

END MODULE mo_omp_block_loop
