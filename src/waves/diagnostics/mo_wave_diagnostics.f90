!>
!! Computes diagnostic parameters and some diagnostics in the wave model
!!
!!
!! @par Revision History
!! Initial release by Mikhail Dorbynin, DWD (2021-02-18)
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

MODULE mo_wave_diagnostics
  USE mo_kind,                ONLY: wp
  USE mo_model_domain,        ONLY: t_patch
  USE mo_impl_constants,      ONLY: min_rlcell, MAX_CHAR_LENGTH
  USE mo_loopindices,         ONLY: get_indices_c

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=MAX_CHAR_LENGTH), PARAMETER :: modname   = 'mo_wave_diagnostics'

  PUBLIC :: significant_wave_height


CONTAINS

  SUBROUTINE significant_wave_height(p_patch, emean, hs)

    TYPE(t_patch),     INTENT(IN)    :: p_patch
    REAL(wp),          INTENT(IN)    :: emean(:,:)  !< total energy [m^2]
    REAL(wp),          INTENT(INOUT) :: hs(:,:)     !< significant wave height [m]

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
         &  routine = 'significant_wave_height'

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc,jb


    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jb,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)
      DO jc = i_startidx, i_endidx
        hs(jc,jb) = 4.0_wp * SQRT(emean(jc,jb))
      END DO
    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL
  END SUBROUTINE significant_wave_height

END MODULE mo_wave_diagnostics
