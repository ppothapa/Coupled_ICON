!>
!! Initialization/reading reading of external datasets
!!
!! This module contains read and initialization routines for the external data state.
!!
!! @author Mikhail Dobrynin, DWD, 18.06.19
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

MODULE mo_wave_ext_data_init

  USE mo_kind,                ONLY: wp
  USE mo_io_units,            ONLY: filename_max
  USE mo_io_config,           ONLY: default_read_method
  USE mo_exception,           ONLY: message
  USE mo_model_domain,        ONLY: t_patch
  USE mo_grid_config,         ONLY: n_dom, nroot
  USE mo_wave_ext_data_types, ONLY: t_external_wave
  USE mo_wave_ext_data_state, ONLY: construct_wave_ext_data_state
  USE mo_var_list,            ONLY: t_var_list_ptr
  USE mo_extpar_config,       ONLY: extpar_filename, generate_filename
  USE mo_read_interface,      ONLY: openInputFile, closeFile, t_stream_id, on_cells, read_2D
  USE mo_master_config,       ONLY: getModelBaseDir

  IMPLICIT NONE

  ! required for reading external data
  INCLUDE 'netcdf.inc'

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_wave_data_init'

  PUBLIC :: init_wave_ext_data

CONTAINS

  SUBROUTINE init_wave_ext_data (p_patch, wave_ext_data, wave_ext_data_list)

    TYPE(t_patch),                      INTENT(IN)    :: p_patch(:)
    TYPE(t_external_wave), ALLOCATABLE, INTENT(INOUT) :: wave_ext_data(:)
    TYPE(t_var_list_ptr),  ALLOCATABLE, INTENT(INOUT) :: wave_ext_data_list(:)
    CHARACTER(len=*), PARAMETER :: routine = modname//':init_wave_ext_data'

    CALL construct_wave_ext_data_state(p_patch, wave_ext_data, wave_ext_data_list)

    CALL read_ext_data_wave(p_patch, wave_ext_data)

    CALL message(TRIM(routine),'finished.')

  END SUBROUTINE init_wave_ext_data


  SUBROUTINE read_ext_data_wave(p_patch, wave_ext_data)

    TYPE(t_patch),         INTENT(IN)    :: p_patch(:)
    TYPE(t_external_wave), INTENT(INOUT) :: wave_ext_data(:)

    INTEGER :: jg

    TYPE(t_stream_id) :: stream_id

    CHARACTER(filename_max) :: extpar_file

    DO jg = 1, n_dom
      extpar_file = generate_filename(extpar_filename, getModelBaseDir(), &
        &                             TRIM(p_patch(jg)%grid_filename),    &
        &                              nroot,                             &
        &                             p_patch(jg)%level, p_patch(jg)%id)

      CALL openInputFile(stream_id, extpar_file, p_patch(jg), default_read_method)

      CALL read_2D(stream_id, on_cells, 'z', wave_ext_data(jg)%bathymetry_c)

      CALL closeFile(stream_id)

      !set minimal depth in m
      WHERE (wave_ext_data(jg)%bathymetry_c(:,:).lt.0.5_wp)
        wave_ext_data(jg)%bathymetry_c(:,:) = 0.5_wp
      END WHERE
    END DO

  END SUBROUTINE read_ext_data_wave

END MODULE mo_wave_ext_data_init
