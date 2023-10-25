! Contains the subroutines for updating time-dependent
! wave physics parameters
!
!
! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_wave_td_update

  USE mo_kind,                ONLY: wp, vp
  USE mo_model_domain,        ONLY: t_patch
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, min_rlcell
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_fortran_tools,       ONLY: copy
  USE mo_math_gradients,      ONLY: grad_green_gauss_cell
  USE mo_math_constants,      ONLY: rad2deg, dbl_eps

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: update_bathymetry_gradient
  PUBLIC :: update_wind_speed_and_direction
  PUBLIC :: update_ice_free_mask

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_wave_td_update'

CONTAINS

  !>
  !! calculate bathymetry gradient
  !!
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

  !>
  !! calculate wind speed and direction (deg) from U and V
  !!
  SUBROUTINE update_wind_speed_and_direction(p_patch, u10, v10, wsp, wdir)

    TYPE(t_patch),     INTENT(IN)    :: p_patch
    REAL(wp),          INTENT(IN)    :: u10(:,:), v10(:,:) ! wind components at 10m above sea
    REAL(wp),          INTENT(INOUT) :: wsp(:,:)           ! wind speed at 10m above sea
    REAL(wp),          INTENT(INOUT) :: wdir(:,:)          ! wind direction

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: routine = modname//':update_wind_speed_and_direction'

    INTEGER :: jc, jb
    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    REAL(wp):: uc, vc

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,uc,vc)
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,      &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)
      DO jc = i_startidx, i_endidx

        uc = SIGN(MAX(ABS(u10(jc,jb)),dbl_eps),u10(jc,jb))
        !IF (uc.eq.0) uc = dbl_eps
        vc = SIGN(MAX(ABS(v10(jc,jb)),dbl_eps),v10(jc,jb))
        !IF (vc.eq.0) vc = dbl_eps

        wsp(jc,jb) = SQRT( uc**2 + vc**2 )
        wdir(jc,jb) = ATAN2(vc,uc)*rad2deg

      END DO
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE update_wind_speed_and_direction


  !>
  !! Calculation of ice mask
  !!
  !! Set the ice-free mask to 1
  !! if the sea ice concentration less than the threshold value trhl_ice
  !!
  SUBROUTINE update_ice_free_mask(p_patch, sea_ice_c, ice_free_mask)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = modname//':ice_mask_c'

    TYPE(t_patch),        INTENT(IN)    :: p_patch
    REAL(wp),             INTENT(IN)    :: sea_ice_c(:,:) ! sea ice concentration at centers (fraction of 1)
    INTEGER,              INTENT(INOUT) :: ice_free_mask(:,:)

    INTEGER  :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER  :: i_startidx, i_endidx
    INTEGER  :: jb,jc
    REAL(wp) :: trhl_ice

    trhl_ice = 0.5_wp ! add to nml?

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
           &              i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jc = i_startidx, i_endidx
        IF (sea_ice_c(jc,jb) < trhl_ice) THEN
          ice_free_mask(jc,jb) = 1
        ELSE
          ice_free_mask(jc,jb) = 0
        END IF
      END DO

    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE update_ice_free_mask

END MODULE mo_wave_td_update
