!>
!! This module prepares aerosol for the use in radiation
!!
!! @author Daniel Rieger, Deutscher Wetterdienst, Offenbach
!!
!! @par Revision History
!! Initial release by Daniel Rieger, Deutscher Wetterdienst, Offenbach (2022-11-08)
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
MODULE mo_nwp_aerosol

! ICON infrastructure
  USE mo_kind,                 ONLY: wp
  USE mo_exception,            ONLY: finish, message, message_text
  USE mo_model_domain,         ONLY: t_patch
  USE mo_parallel_config,      ONLY: nproma
  USE mo_loopindices,          ONLY: get_indices_c
  USE mo_impl_constants,       ONLY: min_rlcell_int, SUCCESS
  USE mo_impl_constants_grf,   ONLY: grf_bdywidth_c
  USE mo_fortran_tools,        ONLY: set_acc_host_or_device
! External infrastruture
  USE mtime,                   ONLY: datetime, timedelta, newDatetime, newTimedelta,       &
                                 &   operator(+), deallocateTimedelta, deallocateDatetime
! Radiation-specific
  USE mo_radiation_config,     ONLY: irad_aero, iRadAeroConstKinne, iRadAeroKinne,         &
                                 &   iRadAeroVolc, iRadAeroKinneVolc, iRadAeroART,         &
                                 &   iRadAeroKinneVolcSP, iRadAeroKinneSP, iRadAeroTegen
! Aerosol-specific
  USE mo_bc_aeropt_kinne,      ONLY: read_bc_aeropt_kinne, set_bc_aeropt_kinne
  USE mo_bc_aeropt_cmip6_volc, ONLY: read_bc_aeropt_cmip6_volc, add_bc_aeropt_cmip6_volc
  USE mo_bc_aeropt_splumes,    ONLY: add_bc_aeropt_splumes

#ifdef __ICON_ART
  USE mo_aerosol_util,            ONLY: tegen_scal_factors
  USE mo_art_radiation_interface, ONLY: art_rad_aero_interface
#endif

  IMPLICIT NONE

  PRIVATE

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_nwp_aerosol'

  PUBLIC :: nwp_aerosol_interface
  PUBLIC :: nwp_aerosol_cleanup

CONTAINS

  !---------------------------------------------------------------------------------------
  !>
  !! @par Revision History
  !! Initial release by Daniel Rieger, Deutscher Wetterdienst, Offenbach (2022-11-08)
  !!
  SUBROUTINE nwp_aerosol_interface(mtime_datetime, pt_patch, zf, zh, dz, dt_rad,                   &
    &                              inwp_radiation, nbands_lw, nbands_sw, wavenum1_sw, wavenum2_sw, &
    &                              zaeq1, zaeq2, zaeq3, zaeq4, zaeq5,                              &
    &                              od_lw, od_sw, ssa_sw, g_sw, lacc)
    CHARACTER(len=*), PARAMETER :: &
      &  routine = modname//':nwp_radiation'

    TYPE(datetime), POINTER, INTENT(in) :: &
      &  mtime_datetime          !< Current datetime
    TYPE(t_patch), TARGET, INTENT(in) :: &
      &  pt_patch                !< Grid/patch info
    REAL(wp), INTENT(in) ::    &
      &  zf(:,:,:), zh(:,:,:), & !< model full/half layer height
      &  dz(:,:,:),            & !< Layer thickness
      &  dt_rad                  !< Radiation time step
    REAL(wp), POINTER, INTENT(in) :: &
      &  wavenum1_sw(:),       & !< Shortwave wavenumber lower band bounds
      &  wavenum2_sw(:)          !< Shortwave wavenumber upper band bounds
    REAL(wp), INTENT(in) ::    &
      &  zaeq1(:,:,:),         & !< Tegen optical thicknesses       1: continental
      &  zaeq2(:,:,:),         & !< relative to 550 nm, including   2: maritime
      &  zaeq3(:,:,:),         & !< a vertical profile              3: desert
      &  zaeq4(:,:,:),         & !< for 5 different                 4: urban
      &  zaeq5(:,:,:)            !< aerosol species.                5: stratospheric background
    INTEGER, INTENT(in) ::     &
      &  inwp_radiation,       & !< Radiation scheme (1=rrtmg, 4=ecrad)
      &  nbands_lw, nbands_sw    !< Number of short and long wave bands
    REAL(wp), ALLOCATABLE, INTENT(out) :: &
      &  od_lw(:,:,:,:),       & !< Longwave optical thickness
      &  od_sw(:,:,:,:),       & !< Shortwave optical thickness
      &  ssa_sw(:,:,:,:),      & !< Shortwave asymmetry factor
      &  g_sw(:,:,:,:)           !< Shortwave single scattering albedo
    LOGICAL, OPTIONAL, INTENT(in) :: lacc ! If true, use openacc

    LOGICAL :: lzacc

    ! Local variables
#ifdef __ICON_ART
    REAL(wp), ALLOCATABLE ::   &
      &  od_lw_art_vr(:,:,:),  &      !< AOD LW (vertically reversed)
      &  od_sw_art_vr(:,:,:),  &      !< AOD SW (vertically reversed)
      &  ssa_sw_art_vr(:,:,:), &      !< SSA SW (vertically reversed)
      &  g_sw_art_vr(:,:,:)           !< Assymetry parameter SW (vertically reversed)
    INTEGER :: jc, jk, jk_vr, jband   !< Loop indices
#endif
    INTEGER :: jb, rl_start, rl_end, i_startblk, i_endblk, i_startidx, i_endidx, istat

    CALL set_acc_host_or_device(lzacc, lacc)

    rl_start   = grf_bdywidth_c+1
    rl_end     = min_rlcell_int
    i_startblk = pt_patch%cells%start_block(rl_start)
    i_endblk   = pt_patch%cells%end_block(rl_end)
  
    SELECT CASE(irad_aero)
! Tegen aerosol
      CASE(iRadAeroTegen)
        ! Open TODO: Move nwp_aerosol from rrtm_interface to here
! Prognostic ART aerosol plus Tegen
      CASE(iRadAeroART)
#if defined(__ECRAD) && defined(__ICON_ART)
#ifdef _OPENACC
        IF (lzacc) CALL finish(routine, "irad_aero==iRadAeroART is not ported to openACC.")
#endif
        IF (inwp_radiation == 4) THEN
          ! Allocations
          ALLOCATE(od_lw        (nproma,pt_patch%nlev,pt_patch%nblks_c,nbands_lw), &
            &      od_sw        (nproma,pt_patch%nlev,pt_patch%nblks_c,nbands_sw), &
            &      ssa_sw       (nproma,pt_patch%nlev,pt_patch%nblks_c,nbands_sw), &
            &      g_sw         (nproma,pt_patch%nlev,pt_patch%nblks_c,nbands_sw), &
            &      od_lw_art_vr (nproma,pt_patch%nlev,                 nbands_lw), &
            &      od_sw_art_vr (nproma,pt_patch%nlev,                 nbands_lw), &
            &      ssa_sw_art_vr(nproma,pt_patch%nlev,                 nbands_lw), &
            &      g_sw_art_vr  (nproma,pt_patch%nlev,                 nbands_lw), &
            &      STAT=istat)
          IF(istat /= SUCCESS) &
            &  CALL finish(routine, 'Allocation of od_lw, od_sw, ssa_sw, g_sw plus ART variants failed')
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,                              &
!$OMP            od_lw_art_vr,od_sw_art_vr,ssa_sw_art_vr,g_sw_art_vr, &
!$OMP            jc,jk,jk_vr,jband) ICON_OMP_DEFAULT_SCHEDULE
          DO jb = i_startblk,i_endblk
            CALL get_indices_c(pt_patch,jb,i_startblk,i_endblk,i_startidx,i_endidx,rl_start,rl_end)
            IF (i_startidx>i_endidx) CYCLE

            CALL art_rad_aero_interface(zaeq1(:,:,jb),zaeq2(:,:,jb),       & !
              &                         zaeq3(:,:,jb),zaeq4(:,:,jb),       & !< Tegen aerosol
              &                         zaeq5(:,:,jb),                     & !
              &                         tegen_scal_factors%absorption,     & !
              &                         tegen_scal_factors%scattering,     & !< Tegen coefficients
              &                         tegen_scal_factors%asymmetry,      & !
              &                         pt_patch%id, jb, 1, pt_patch%nlev, & !< Indices domain, block, level
              &                         i_startidx, i_endidx,              & !< Indices nproma loop
              &                         nbands_lw,                         & !< Number of SW bands
              &                         nbands_sw,                         & !< Number of LW bands
              &                         od_lw_art_vr(:,:,:),               & !< OUT: Optical depth LW
              &                         od_sw_art_vr(:,:,:),               & !< OUT: Optical depth SW
              &                         ssa_sw_art_vr(:,:,:),              & !< OUT: SSA SW
              &                         g_sw_art_vr(:,:,:))                  !< OUT: Assymetry parameter SW

! LONGWAVE
            DO jc = i_startidx, i_endidx
              DO jk = 1, pt_patch%nlev
                jk_vr = pt_patch%nlev+1-jk
!NEC$ nointerchange
!NEC$ nounroll
                DO jband = 1, nbands_lw
                  od_lw(jc,jk,jb,jband) = od_lw_art_vr(jc,jk_vr,jband)
                ENDDO
              ENDDO
            ENDDO

! SHORTWAVE
            DO jc = i_startidx, i_endidx
              DO jk = 1, pt_patch%nlev
                jk_vr = pt_patch%nlev+1-jk
!NEC$ nointerchange
!NEC$ nounroll
                DO jband = 1, nbands_sw
                  od_sw(jc,jk,jb,jband) = od_sw_art_vr(jc,jk_vr,jband)
                  ssa_sw(jc,jk,jb,jband) = ssa_sw_art_vr(jc,jk_vr,jband)
                  g_sw(jc,jk,jb,jband) = g_sw_art_vr(jc,jk_vr,jband)
                ENDDO
              ENDDO
            ENDDO

          END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

          ! Deallocations
          DEALLOCATE(od_lw_art_vr, od_sw_art_vr, ssa_sw_art_vr, g_sw_art_vr, &
            &        STAT=istat)
          IF(istat /= SUCCESS) &
            &  CALL finish(routine, 'Deallocation of od_lw_art_vr, od_sw_art_vr, ssa_sw_art_vr, g_sw_art_vr failed')
        ENDIF
#endif

! Kinne aerosol
      CASE(iRadAeroConstKinne, iRadAeroKinne, iRadAeroVolc, iRadAeroKinneVolc, iRadAeroKinneVolcSP, iRadAeroKinneSP)
        ! Compatibility checks
#ifdef __ECRAD
        IF (inwp_radiation /= 4) THEN
          WRITE(message_text,'(a,i2,a)') 'irad_aero = ', irad_aero,' only implemented for ecrad (inwp_radiation=4).'
          CALL finish(routine, message_text)
        ENDIF
#else
        WRITE(message_text,'(a,i2,a)') 'irad_aero = ', irad_aero,' requires to compile with --enable-ecrad.'
        CALL finish(routine, message_text)
#endif
#ifdef _OPENACC
        IF (lzacc) CALL finish(routine, "irad_aero==*Kinne* is not ported to openACC.")
#endif

        ! Update Kinne aerosol from files once per day
        CALL nwp_aerosol_daily_update_kinne(mtime_datetime, pt_patch, dt_rad, inwp_radiation, &
          &                                 nbands_lw, nbands_sw)

        ! Allocations
        ALLOCATE(od_lw (nproma,pt_patch%nlev,pt_patch%nblks_c,nbands_lw)  , &
          &      od_sw (nproma,pt_patch%nlev,pt_patch%nblks_c,nbands_sw)  , &
          &      ssa_sw(nproma,pt_patch%nlev,pt_patch%nblks_c,nbands_sw)  , &
          &      g_sw  (nproma,pt_patch%nlev,pt_patch%nblks_c,nbands_sw)  , &
          &      STAT=istat)
        IF(istat /= SUCCESS) &
          &  CALL finish(routine, 'Allocation of od_lw, od_sw, ssa_sw, g_sw failed')

        IF ( .NOT. ASSOCIATED(wavenum1_sw) .OR. .NOT. ASSOCIATED(wavenum2_sw) ) &
          &  CALL finish(routine, 'wavenum1 or wavenum2 not associated')
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk,i_endblk
          CALL get_indices_c(pt_patch,jb,i_startblk,i_endblk,i_startidx,i_endidx,rl_start,rl_end)
          IF (i_startidx>i_endidx) CYCLE

          CALL nwp_aerosol_kinne(mtime_datetime, zf(:,:,jb), zh(:,:,jb), dz(:,:,jb),   &
            &                    pt_patch%id, jb, i_startidx, i_endidx, pt_patch%nlev, &
            &                    nbands_lw, nbands_sw, wavenum1_sw(:), wavenum2_sw(:), &
            &                    od_lw(:,:,jb,:), od_sw(:,:,jb,:),                     &
            &                    ssa_sw(:,:,jb,:), g_sw(:,:,jb,:)                      )

        END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      CASE DEFAULT
        ! Currently continue as not all cases are ported to nwp_aerosol_interface yet
    END SELECT

  END SUBROUTINE nwp_aerosol_interface

  !---------------------------------------------------------------------------------------
  !>
  !! @par Revision History
  !! Initial release by Daniel Rieger, Deutscher Wetterdienst, Offenbach (2022-11-10)
  !!
  SUBROUTINE nwp_aerosol_daily_update_kinne(mtime_datetime, pt_patch, dt_rad, inwp_radiation, nbands_lw, nbands_sw)
    TYPE(datetime), POINTER, INTENT(in) :: &
      &  mtime_datetime                    !< Current datetime
    TYPE(t_patch), TARGET, INTENT(in) :: &
      &  pt_patch                          !< Grid/patch info
    REAL(wp), INTENT(in) ::              &
      &  dt_rad                            !< Radiation time step
    INTEGER, INTENT(in) ::               &
      &  inwp_radiation,                 & !< Radiation scheme (1=rrtmg, 4=ecrad)
      &  nbands_lw, nbands_sw              !< Number of short and long wave bands
    ! Local variables
    TYPE(datetime), POINTER ::           &
      &  prev_radtime                      !< Datetime of previous radiation time step
    TYPE(timedelta), POINTER ::          &
      &  td_dt_rad                         !< Radiation time step

    td_dt_rad => newTimedelta('-',0,0,0,0,0, second=NINT(dt_rad), ms=0)
    prev_radtime => newDatetime(mtime_datetime + td_dt_rad)

    IF (prev_radtime%date%day /= mtime_datetime%date%day) THEN
      IF (inwp_radiation == 4) THEN
        IF (ANY(irad_aero == [iRadAeroKinne, iRadAeroKinneVolc, iRadAeroKinneVolcSP, iRadAeroKinneSP])) &
            & CALL read_bc_aeropt_kinne(mtime_datetime, pt_patch, .TRUE., nbands_lw, nbands_sw)
        IF (ANY(irad_aero == [iRadAeroVolc, iRadAeroKinneVolc, iRadAeroKinneVolcSP])) &
            & CALL read_bc_aeropt_cmip6_volc(mtime_datetime, nbands_lw, nbands_sw)
      ENDIF
    ENDIF

    CALL deallocateTimedelta(td_dt_rad)
    CALL deallocateDatetime(prev_radtime)

  END SUBROUTINE nwp_aerosol_daily_update_kinne

  !---------------------------------------------------------------------------------------
  !>
  !! @par Revision History
  !! Initial release by Daniel Rieger, Deutscher Wetterdienst, Offenbach (2022-11-09)
  !!
  SUBROUTINE nwp_aerosol_kinne(mtime_datetime, zf, zh, dz, jg, jb, i_startidx, i_endidx, nlev, &
    &                          nbands_lw, nbands_sw, wavenum1_sw, wavenum2_sw,                 &
    &                          od_lw, od_sw, ssa_sw, g_sw)
    TYPE(datetime), POINTER, INTENT(in) :: &
      &  mtime_datetime                      !< Current datetime
    REAL(wp), INTENT(in) ::                &
      &  zf(:,:), zh(:,:), dz(:,:),        & !< model full/half layer height, layer thickness
      &  wavenum1_sw(:),                   & !< Shortwave wavenumber lower band bounds
      &  wavenum2_sw(:)                      !< Shortwave wavenumber upper band bounds
    INTEGER, INTENT(in) ::                 &
      &  jg, jb,                           & !< Domain and block index
      &  i_startidx, i_endidx,             & !< Loop indices
      &  nlev,                             & !< Number of vertical levels
      &  nbands_lw, nbands_sw                !< Number of short and long wave bands
    REAL(wp), INTENT(out) ::               &
      &  od_lw(:,:,:), od_sw(:,:,:),       & !< LW/SW optical thickness
      &  ssa_sw(:,:,:), g_sw(:,:,:)          !< SW asymmetry factor, SW single scattering albedo
    ! Local variables
    REAL(wp) ::                            &
      &  od_lw_vr (nproma,nlev,nbands_lw), & !< LW optical thickness of aerosols    (vertically reversed)
      &  od_sw_vr (nproma,nlev,nbands_sw), & !< SW aerosol optical thickness        (vertically reversed)
      &  g_sw_vr  (nproma,nlev,nbands_sw), & !< SW aerosol asymmetry factor         (vertically reversed)
      &  ssa_sw_vr(nproma,nlev,nbands_sw)    !< SW aerosol single scattering albedo (vertically reversed)
    REAL(wp) ::                            &
      &  x_cdnc(nproma)                      !< Scale factor for Cloud Droplet Number Concentration (currently not used)
    INTEGER ::                             &
      &  jk                                  !< Loop index

    od_lw_vr(:,:,:)  = 0.0_wp
    od_sw_vr(:,:,:)  = 0.0_wp
    ssa_sw_vr(:,:,:) = 1.0_wp
    g_sw_vr (:,:,:)  = 0.0_wp

    ! Tropospheric Kinne aerosol
    IF (ANY( irad_aero == (/iRadAeroConstKinne,iRadAeroKinne,iRadAeroKinneVolc, &
      &                     iRadAeroKinneVolcSP,iRadAeroKinneSP/) )) THEN
      CALL set_bc_aeropt_kinne(mtime_datetime, jg, 1, i_endidx, nproma, nlev, jb, &
        &                      nbands_sw, nbands_lw, zf(:,:), dz(:,:),            &
        &                      od_sw_vr(:,:,:), ssa_sw_vr(:,:,:),                 &
        &                      g_sw_vr (:,:,:), od_lw_vr(:,:,:)                   )
    ENDIF

    ! Volcanic stratospheric aerosols for CMIP6
    IF (ANY( irad_aero == (/iRadAeroVolc,iRadAeroKinneVolc,iRadAeroKinneVolcSP/) )) THEN 
     CALL add_bc_aeropt_cmip6_volc(mtime_datetime, jg, 1, i_endidx, nproma, nlev, jb, &
       &                           nbands_sw, nbands_lw, zf(:,:), dz(:,:),            &
       &                           od_sw_vr(:,:,:), ssa_sw_vr(:,:,:),                 &
       &                           g_sw_vr (:,:,:), od_lw_vr(:,:,:)                   )
    END IF

    ! Simple plumes
    IF (ANY( irad_aero == (/iRadAeroKinneVolcSP,iRadAeroKinneSP/) )) THEN
      CALL add_bc_aeropt_splumes(jg, 1, i_endidx, nproma, nlev, jb,  &
        &                        nbands_sw, mtime_datetime,          &
        &                        zf(:,:), dz(:,:), zh(:,nlev+1),     &
        &                        wavenum1_sw(:), wavenum2_sw(:),     &
        &                        od_sw_vr(:,:,:), ssa_sw_vr(:,:,:),  &
        &                        g_sw_vr (:,:,:), x_cdnc(:)          )
    END IF

    ! Vertically reverse the fields:
    DO jk = 1, nlev
      od_lw (:,jk,:) = od_lw_vr (:,nlev-jk+1,:)
      od_sw (:,jk,:) = od_sw_vr (:,nlev-jk+1,:)
      ssa_sw(:,jk,:) = ssa_sw_vr(:,nlev-jk+1,:)
      g_sw  (:,jk,:) = g_sw_vr  (:,nlev-jk+1,:)
    ENDDO

  END SUBROUTINE nwp_aerosol_kinne

  !---------------------------------------------------------------------------------------
  !>
  !! @par Revision History
  !! Initial release by Daniel Rieger, Deutscher Wetterdienst, Offenbach (2022-11-10)
  !!
  SUBROUTINE nwp_aerosol_cleanup(od_lw, od_sw, ssa_sw, g_sw)
    CHARACTER(len=*), PARAMETER :: &
      &  routine = modname//':nwp_aerosol_cleanup'

    REAL(wp), ALLOCATABLE, INTENT(inout) :: &
      &  od_lw(:,:,:,:),       & !< Longwave optical thickness
      &  od_sw(:,:,:,:),       & !< Shortwave optical thickness
      &  ssa_sw(:,:,:,:),      & !< Shortwave asymmetry factor
      &  g_sw(:,:,:,:)           !< Shortwave single scattering albedo
    ! Local variables
    INTEGER :: istat

    IF( ALLOCATED(od_lw) ) THEN
      DEALLOCATE(od_lw, STAT=istat)
      IF(istat /= SUCCESS) CALL finish(routine, 'Deallocation of od_lw failed.')
    ENDIF
    IF( ALLOCATED(od_sw) ) THEN
      DEALLOCATE(od_sw, STAT=istat)
      IF(istat /= SUCCESS) CALL finish(routine, 'Deallocation of od_sw failed.')
    ENDIF
    IF( ALLOCATED(ssa_sw) ) THEN
      DEALLOCATE(ssa_sw, STAT=istat)
      IF(istat /= SUCCESS) CALL finish(routine, 'Deallocation of ssa_sw failed.')
    ENDIF
    IF( ALLOCATED(g_sw) ) THEN
      DEALLOCATE(g_sw, STAT=istat)
      IF(istat /= SUCCESS) CALL finish(routine, 'Deallocation of g_sw failed.')
    ENDIF

  END SUBROUTINE nwp_aerosol_cleanup

END MODULE mo_nwp_aerosol
