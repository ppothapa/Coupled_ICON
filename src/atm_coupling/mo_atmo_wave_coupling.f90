!>
!! @brief Interface between atmosphere physics and the ocean surface waves, through a coupler
!!
!! @author Mikhail Dobrynin (DWD)
!! @author Daniel Reinert (DWD)
!!
!! @par Revision History
!! Initial submission
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_atmo_wave_coupling

  USE mo_kind,             ONLY: wp
  USE mo_exception,        ONLY: warning, message_text, message
  USE mo_model_domain,     ONLY: t_patch
  USE mo_fortran_tools,    ONLY: assert_acc_host_only
  USE mo_run_config,       ONLY: ltimer
  USE mo_timer,            ONLY: timer_start, timer_stop, &
    &                            timer_coupling_put, timer_coupling_get, &
    &                            timer_coupling_1stget
  USE mo_sync,             ONLY: SYNC_C, sync_patch_array
  USE mo_yac_finterface,   ONLY: yac_fput, yac_fget, yac_dble_ptr,         &
    &                            YAC_ACTION_COUPLING, YAC_ACTION_OUT_OF_BOUND
  USE mo_atmo_wave_coupling_frame, ONLY: lyac_very_1st_get, field_id, collection_size, &
    &                                    CPF_U10M, CPF_V10M, CPF_FR_SEAICE, CPF_Z0

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: couple_atmo_to_wave

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_atmo_wave_coupling'

CONTAINS


  !>
  !! Exchange fields between atmosphere and wave model
  !!
  !! Send fields to the wave model:
  !!   "meridional_wind_in_10m"
  !!   "zonal_wind_in_10m"
  !!   "fraction_of_ocean_covered_by_sea_ice"
  !!
  !! Receive fields from the wave model:
  !!   "roughness_length"
  !!
  !! This subroutine is called from nwp_nh_interface.
  !!
  !! @par Revision History
  !! Initial revision by Mikhail Dobrynin and Daniel Reinert, DWD (2023-08-24)
  !!
  SUBROUTINE couple_atmo_to_wave(p_patch, u10m, v10m, fr_seaice, z0_waves, lacc)

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = modname//':couple_atmo_to_wave'

    TYPE(t_patch),                INTENT(IN)   :: p_patch
    REAL(wp), CONTIGUOUS, TARGET, INTENT(IN)   :: u10m(:,:)      !< zonal wind speed in 10m [m/s]
    REAL(wp), CONTIGUOUS, TARGET, INTENT(IN)   :: v10m(:,:)      !< meridional wind speed in 10m [m/s]
    REAL(wp), CONTIGUOUS, TARGET, INTENT(IN)   :: fr_seaice(:,:) !< fraction_of_ocean_covered_by_sea_ice [1]
    REAL(wp), CONTIGUOUS, TARGET, INTENT(INOUT):: z0_waves(:,:)  !< surface roughness length [m]
    LOGICAL,  OPTIONAL,           INTENT(IN)   :: lacc           ! If true, use openacc

    LOGICAL :: write_coupler_restart
    INTEGER :: info, ierror
    TYPE(yac_dble_ptr) :: yac_ptr(1,SIZE(field_id))


    CALL assert_acc_host_only('couple_atmo_to_wave', lacc)

    write_coupler_restart = .FALSE.

    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !  Send fields from atmosphere to wave
    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !
    !------------------------------------------------
    !  Send 10m zonal wind u10m
    !  field_id(CPF_U10M) = 'zonal_wind_in_10m'
    !------------------------------------------------
    !
    IF (ltimer) CALL timer_start(timer_coupling_put)
    !
    yac_ptr(1,CPF_U10M)%p(1:p_patch%n_patch_cells) => u10m(:,:)
    !
    CALL yac_fput (                                              &
      & field_id        = field_id(CPF_U10M),                    &
      & nbr_pointsets   = SIZE(yac_ptr,1),                       &
      & collection_size = collection_size(CPF_U10M),             &
      & send_field      = yac_ptr(:,CPF_U10M:CPF_U10M),          &
      & info            = info,                                  &
      & ierror          = ierror )
    !
    IF (ltimer) CALL timer_stop(timer_coupling_put)
    !
    IF ( info > YAC_ACTION_COUPLING .AND. info < YAC_ACTION_OUT_OF_BOUND ) write_coupler_restart = .TRUE.
    IF ( info == YAC_ACTION_OUT_OF_BOUND ) THEN
      WRITE(message_text, '(a,i3,a)') 'YAC says fput called after end of run - id=',CPF_U10M,', U10'
      CALL warning(routine, message_text)
    END IF

    ! ----------------------------------------------
    !  Send 10m meridional wind v10m
    !  field_id(CPF_V10M) = 'meridional_wind_in_10m'
    ! ----------------------------------------------
    !
    IF (ltimer) CALL timer_start(timer_coupling_put)
    !
    yac_ptr(1,CPF_V10M)%p(1:p_patch%n_patch_cells) => v10m(:,:)
    !
    CALL yac_fput (                                              &
      & field_id        = field_id(CPF_V10M),                    &
      & nbr_pointsets   = SIZE(yac_ptr,1),                       &
      & collection_size = collection_size(CPF_V10M),             &
      & send_field      = yac_ptr(:,CPF_V10M:CPF_V10M),          &
      & info            = info,                                  &
      & ierror          = ierror )
    !
    IF (ltimer) CALL timer_stop(timer_coupling_put)
    !
    IF ( info > YAC_ACTION_COUPLING .AND. info < YAC_ACTION_OUT_OF_BOUND ) write_coupler_restart = .TRUE.
    IF ( info == YAC_ACTION_OUT_OF_BOUND ) THEN
      WRITE(message_text, '(a,i3,a)') 'YAC says fput called after end of run - id=',CPF_V10M,', V10'
      CALL warning(routine, message_text)
    END IF

    ! ------------------------------------------------------------------
    !  Send fraction of sea ice
    !  field_id(CPF_FR_SEAICE) = 'fraction_of_ocean_covered_by_sea_ice'
    ! ------------------------------------------------------------------
    !
    !
    IF (ltimer) CALL timer_start(timer_coupling_put)
    !
    yac_ptr(1,CPF_FR_SEAICE)%p(1:p_patch%n_patch_cells) => fr_seaice(:,:)
    !
    CALL yac_fput (                                              &
      & field_id        = field_id(CPF_FR_SEAICE),               &
      & nbr_pointsets   = SIZE(yac_ptr,1),                       &
      & collection_size = collection_size(CPF_FR_SEAICE),        &
      & send_field      = yac_ptr(:,CPF_FR_SEAICE:CPF_FR_SEAICE),&
      & info            = info,                                  &
      & ierror          = ierror )
    !
    IF (ltimer) CALL timer_stop(timer_coupling_put)
    !
    IF ( info > YAC_ACTION_COUPLING .AND. info < YAC_ACTION_OUT_OF_BOUND ) write_coupler_restart = .TRUE.
    IF ( info == YAC_ACTION_OUT_OF_BOUND ) THEN
      WRITE(message_text, '(a,i3,a)') 'YAC says fput called after end of run - id=',CPF_FR_SEAICE,', fr_seaice'
      CALL warning(routine, message_text)
    END IF

    IF (write_coupler_restart) THEN
      WRITE(message_text, '(a,3i3,a)') 'YAC says it is put for restart - id=', &
        &                             CPF_U10M, CPF_V10M, CPF_FR_SEAICE,', U10, V10, fr_seaice'
      CALL message(routine, message_text)
    ENDIF


    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !  Receive fields from wave to atmosphere
    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !
    ! --------------------------------------------
    !  Receive roughness length z0 from the wave model
    !  field_name(CPF_Z0) = 'roughness_length'
    ! --------------------------------------------
    !
    IF (lyac_very_1st_get) THEN
      IF (ltimer) CALL timer_start(timer_coupling_1stget)
    ELSE
      IF (ltimer) CALL timer_start(timer_coupling_get)
    ENDIF
    !
    yac_ptr(1,CPF_Z0)%p(1:p_patch%n_patch_cells) => z0_waves(:,:)
    !
    CALL yac_fget(                                  &
      & field_id        = field_id(CPF_Z0),         &
      & collection_size = collection_size(CPF_Z0),  &
      & recv_field      = yac_ptr(:,CPF_Z0),        &
      & info            = info,                     &
      & ierror          = ierror )
    !
    IF (lyac_very_1st_get) THEN
      IF (ltimer) CALL timer_stop(timer_coupling_1stget)
      lyac_very_1st_get = .FALSE.
    ELSE
      IF (ltimer) CALL timer_stop(timer_coupling_get)
    ENDIF
    !
    IF ( info > YAC_ACTION_COUPLING .AND. info < YAC_ACTION_OUT_OF_BOUND ) THEN
      WRITE(message_text, '(a,i3,a)') 'YAC says it is get for restart - id=',CPF_Z0,', Z0'
      CALL message(routine, message_text)
    ENDIF
    IF ( info == YAC_ACTION_OUT_OF_BOUND ) THEN
      WRITE(message_text, '(a,i3,a)') 'YAC says fget called after end of run - id=',CPF_Z0,', Z0'
      CALL warning(routine, message_text)
    ENDIF

    ! halo synchronization for fields recieved from the atmosphere
    !
    CALL sync_patch_array(SYNC_C, p_patch, z0_waves(:,:), opt_varname="z0")

  END SUBROUTINE couple_atmo_to_wave

END MODULE mo_atmo_wave_coupling
