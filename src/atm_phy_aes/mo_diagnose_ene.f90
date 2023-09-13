!>
!! @brief Module containing subroutine computing near surface energetics
!!
!! @contact Marco A. Giorgetta, MPI-M, 2023-07
!!
!! @par Copyright and License
!!
!! Copyright 2023 Max Planck Institute for Meteorology, Marco A. Giorgetta
!!
!! Code subject to BSD-3-C, SPDX short identifier: BSD-3-Clause, see file
!! BSD-3-C-license.pdf in the license-directory
!!
MODULE mo_diagnose_ene

  USE mo_run_config         ,ONLY: num_lev, iqv
  USE mo_dynamics_config    ,ONLY: nnew_rcf

  USE mo_nonhydro_state     ,ONLY: p_nh_state
  USE mo_aes_phy_memory     ,ONLY: prm_field

  USE mo_physical_constants ,ONLY: cpd, alv

  USE mo_timer              ,ONLY: ltimer, timer_start, timer_stop, timer_ene

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: diagnose_ene

CONTAINS

  !-------------------------------------------------------------------

  SUBROUTINE diagnose_ene(jg, jb, jcs, jce)

    INTEGER , INTENT(in)    :: jg, jb, jcs, jce

    INTEGER                 :: jc, jtl_trc

    IF (ltimer) call timer_start(timer_ene)

    jtl_trc = nnew_rcf(jg)

    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1)
    DO jc = jcs, jce
      prm_field(jg)%shfl_qsa(jc,jb) = prm_field(jg)%shflx(jc,jb)/cpd * p_nh_state(jg)%prog(jtl_trc)%tracer(jc,num_lev(jg),jb,iqv)
      prm_field(jg)%evap_tsa(jc,jb) = prm_field(jg)%lhflx(jc,jb)/alv * p_nh_state(jg)%diag%temp(jc,num_lev(jg),jb)
      prm_field(jg)%rsfl_tsa(jc,jb) = prm_field(jg)%rsfl (jc,jb)     * p_nh_state(jg)%diag%temp(jc,num_lev(jg),jb)
      prm_field(jg)%ssfl_tsa(jc,jb) = prm_field(jg)%ssfl (jc,jb)     * p_nh_state(jg)%diag%temp(jc,num_lev(jg),jb)
    END DO !jc
    !$ACC END PARALLEL LOOP

    IF (ltimer) call timer_stop(timer_ene)

  END SUBROUTINE  diagnose_ene

  !-------------------------------------------------------------------

END MODULE mo_diagnose_ene
