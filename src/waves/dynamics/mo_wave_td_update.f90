!! Contains the subroutines for updating time-dependent
!! wave physics parameters
!!
!! @author Mikhail Dobrynin


!! Revision History
!! Initial revision by Mikhail Dobrynin, DWD  (2023.06.05)
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

MODULE mo_wave_td_update

  USE mo_kind,                ONLY: wp, vp
  USE mo_model_domain,        ONLY: t_patch
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, min_rlcell
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_fortran_tools,       ONLY: copy
  USE mo_math_gradients,      ONLY: grad_green_gauss_cell

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: update_bathymetry_gradient

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_wave_td_update'

CONTAINS

  !>
  !! calculate bathymetry gradient
  !!
  !!
  !! @par Revision History
  !! Initial revision by Mikhail Dobrynin, DWD (2023-06-05)
  !!
  SUBROUTINE update_bathymetry_gradient(p_patch, p_int_state, bathymetry_c, geo_bath_grad_c)

    TYPE(t_patch),     INTENT(IN)    :: p_patch
    TYPE(t_int_state), INTENT(IN)    :: p_int_state
    REAL(wp),          INTENT(IN)    :: bathymetry_c(:,:)
    REAL(vp),          INTENT(INOUT) :: geo_bath_grad_c(:,:,:,:)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: routine = modname//':update_bathymetry_gradient'

    REAL(wp) :: bath_c_3d(SIZE(bathymetry_c,1),1,SIZE(bathymetry_c,2))

!$OMP PARALLEL
    CALL copy(src=bathymetry_c, dest=bath_c_3d(:,1,:))
!$OMP END PARALLEL

    CALL grad_green_gauss_cell(bath_c_3d, p_patch, p_int_state, geo_bath_grad_c, &
         &                     opt_slev=1, opt_elev=1, &
         &                     opt_rlstart=2, opt_rlend=min_rlcell)

  END SUBROUTINE update_bathymetry_gradient


END MODULE mo_wave_td_update
