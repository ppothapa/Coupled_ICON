!>
!!  Subroutine to initialize the wave test case
!!
!!
!! @par Revision History
!! - Mikhail Dobrynin, DWD, 14.08.2019
!! - first version by P. Ripodas , DWD, (2011-08)
!! - main parts extracted from the original mo_nh_testcases.f90
!!
!! @par Literature
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
MODULE mo_wave_adv_exp

  USE mo_kind,                 ONLY: wp
  USE mo_model_domain,         ONLY: t_patch
  USE mo_wave_types,           ONLY: t_wave_diag
  USE mo_math_constants,       ONLY: pi, rad2deg
  USE mo_impl_constants,       ONLY: MAX_CHAR_LENGTH, min_rlcell
  USE mo_loopindices,          ONLY: get_indices_c

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_wind_adv_test

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_wave_adv_exp'

  !--------------------------------------------------------------------

CONTAINS

  SUBROUTINE init_wind_adv_test(p_patch, p_diag)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
         &  routine = modname//'::init_wind_adv_test'

    TYPE(t_patch),      INTENT(IN)    :: p_patch
    TYPE(t_wave_diag),  INTENT(INOUT) :: p_diag

    INTEGER :: jc, jb
    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    REAL(wp):: sin_tmp, cos_tmp, zlat, zlon, d1, r
    REAL(wp), PARAMETER ::                    &
      &  RR         = 1._wp/3._wp,            & ! horizontal half width divided by 'a'
      &  lambda0    = 1.5_wp*pi,              & ! center point in longitudes
      &  phi0       = 0.0_wp*pi                 ! center point in latitudes

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,zlon,zlat,sin_tmp,cos_tmp,r,d1)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,      &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jc = i_startidx, i_endidx
        !test wind from NCAR_TESTCASE
        zlon = p_patch%cells%center(jc,jb)%lon
        zlat = p_patch%cells%center(jc,jb)%lat

        sin_tmp = SIN(zlat) * SIN(phi0)
        cos_tmp = COS(zlat) * COS(phi0)
        r  = ACOS (sin_tmp + cos_tmp*COS(zlon-lambda0))       ! great circle distance without 'a'
        d1 = MIN( 1._wp, (r/RR) )

        !45 degree towards NE
        p_diag%u_10m(jc,jb) = 0.5_wp * (1._wp + COS(pi*d1)) * 17.87_wp
        p_diag%v_10m(jc,jb) = 0.5_wp * (1._wp + COS(pi*d1)) * 17.87_wp
        p_diag%sp_10m(jc,jb) = SQRT(p_diag%u_10m(jc,jb)**2 + p_diag%v_10m(jc,jb)**2)
        p_diag%dir_10m(jc,jb) = ATAN2(p_diag%v_10m(jc,jb),p_diag%u_10m(jc,jb))*rad2deg
      END DO ! cell loop
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE init_wind_adv_test

END MODULE mo_wave_adv_exp
