!>
!! @brief Subroutine interface_aes_car calls Cariolle's linearized ozone scheme.
!!
!! @author Hui Wan, MPI-M
!! @author Marco Giorgetta, MPI-M
!!
!! @par Revision History
!!  Original version from ECHAM6 (revision 2028)
!!  Modified for ICOHAM by Hui Wan and Marco Giorgetta (2010)
!!  Modified for ICONAM by Marco Giorgetta (2014)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

MODULE mo_interface_aes_car

  USE mo_kind                ,ONLY: wp
  USE mtime                  ,ONLY: t_datetime => datetime

  USE mo_aes_phy_dims        ,ONLY: aes_phy_dims
  USE mo_aes_phy_config      ,ONLY: aes_phy_config, aes_phy_tc
  USE mo_aes_phy_memory      ,ONLY: t_aes_phy_field, prm_field, &
    &                               t_aes_phy_tend,  prm_tend

  USE mo_timer               ,ONLY: ltimer, timer_start, timer_stop, timer_car

  USE mo_run_config          ,ONLY: io3
  USE mo_physical_constants  ,ONLY: amd, amo3
  USE mo_bcs_time_interpolation ,ONLY: t_time_interpolation_weights, &
       &                               calculate_time_interpolation_weights
  USE mo_lcariolle     ,ONLY: t_avi, t_time_interpolation, lcariolle_do3dt 
  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: interface_aes_car

CONTAINS

  SUBROUTINE interface_aes_car(jg, jb, jcs, jce)

    ! Arguments
    !
    INTEGER, INTENT(in)     :: jg, jb, jcs, jce

    ! Pointers
    !
    TYPE(t_aes_phy_field), POINTER :: field
    TYPE(t_aes_phy_tend),  POINTER :: tend

    ! Local variables
    !
    INTEGER  :: nlev
    INTEGER  :: nproma
    !
    TYPE(t_datetime), POINTER :: datetime
    REAL(wp) :: pdtime
    LOGICAL  :: is_in_sd_ed_interval
    LOGICAL  :: is_active
    !
    INTEGER  :: fc_car
    !
    REAL(wp) :: tend_o3_car(aes_phy_dims(jg)%nproma,aes_phy_dims(jg)%nlev)
    !
    TYPE(t_time_interpolation)          :: time_interpolation
    TYPE(t_time_interpolation_weights)  :: current_time_interpolation_weights
    TYPE(t_avi)                         :: avi

    IF (ltimer) call timer_start(timer_car)

    nlev    = aes_phy_dims(jg)%nlev
    nproma  = aes_phy_dims(jg)%nproma

    datetime             => aes_phy_tc(jg)%datetime
    pdtime               =  aes_phy_tc(jg)%dt_phy_sec
    is_in_sd_ed_interval =  aes_phy_tc(jg)%is_in_sd_ed_interval_car
    is_active            =  aes_phy_tc(jg)%is_active_car

    fc_car    =  aes_phy_config(jg)%fc_car

    ! associate pointers
    field     => prm_field(jg)
    tend      => prm_tend (jg)

    IF ( is_in_sd_ed_interval ) THEN
       !
       IF ( is_active ) THEN
          !
          current_time_interpolation_weights = calculate_time_interpolation_weights(datetime)
          time_interpolation% imonth1 = current_time_interpolation_weights% month1_index
          time_interpolation% imonth2 = current_time_interpolation_weights% month2_index
          time_interpolation% weight1 = current_time_interpolation_weights% weight1
          time_interpolation% weight2 = current_time_interpolation_weights% weight2
          !
          ALLOCATE(avi%o3_vmr(nproma,nlev), avi%vmr2molm2(nproma,nlev), avi%cell_center_lat(nproma), avi%lday(nproma))
          !
          avi%ldown=.TRUE.
          avi%o3_vmr(jcs:jce,:)        =  field% qtrc_phy(jcs:jce,:,jb,io3)*amd/amo3
          avi%tmprt                    => field% ta  (:,:,jb)
          !
          ! Note: ICON has no sources and sinks in te equation for air density. This implies
          !       that the total air mass is conserved. The parameterized turbulent mass flux
          !       at the surface and precipitation have no effect on the atmospheric mass.
          !       Therefore let us use here the total air mass as dry air mass.
          !
          avi%vmr2molm2(jcs:jce,:)     =  field% mair(jcs:jce,:,jb) / amd * 1.e3_wp
          !
          avi%pres                     => field% pfull(:,:,jb)
          avi%cell_center_lat(jcs:jce) =  field% clat(jcs:jce,jb)
          avi%lday(jcs:jce)            =  field% cosmu0(jcs:jce,jb) > 1.e-3_wp
          !
          CALL lcariolle_do3dt(jcs,                   jce,             &
               &               nproma,                nlev,            &
               &               time_interpolation,                     &
               &               avi,                   tend_o3_car(:,:) )
          !
          DEALLOCATE(avi%o3_vmr, avi%vmr2molm2, avi%cell_center_lat, avi%lday)
          !
          ! store in memory for output or recycling
          !
          IF (ASSOCIATED(tend% o3_car)) THEN
             tend% o3_car(jcs:jce,:,jb) = tend_o3_car(jcs:jce,:)
          END IF
          !
       ELSE
          !
          ! retrieve from memory for recycling
          !
          IF (ASSOCIATED(tend% o3_car)) THEN
             tend_o3_car(jcs:jce,:) = tend% o3_car(jcs:jce,:,jb)
          END IF
          !
       END IF
       !
       ! accumulate tendencies for later updating the model state
       SELECT CASE(fc_car)
       CASE(1)
          !
          ! accumulate tendencies for later updating the model state
          tend%  qtrc_phy(jcs:jce,:,jb,io3) = tend%  qtrc_phy(jcs:jce,:,jb,io3) + tend_o3_car(jcs:jce,:)*amo3/amd
          !
          ! update physics state for input to the next physics process
          field% qtrc_phy(jcs:jce,:,jb,io3) = field% qtrc_phy(jcs:jce,:,jb,io3) + tend_o3_car(jcs:jce,:)*amo3/amd*pdtime
       END SELECT
       !
    ELSE
       !
       IF (ASSOCIATED(tend% o3_car)) THEN
          tend% o3_car(jcs:jce,:,jb) = 0.0_wp
       END IF
       !
    END IF
       
    ! disassociate pointers
    NULLIFY(field)
    NULLIFY(tend )
    
    IF (ltimer) call timer_stop(timer_car)

  END SUBROUTINE interface_aes_car

END MODULE mo_interface_aes_car
