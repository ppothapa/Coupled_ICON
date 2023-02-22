!>
!! @brief
!!  Read namelists, make sanity checks specific to each namelist and make
!!  a cross check once all namelists of a component are available.
!!
!! @author
!!  Mikhail Dobrynin, DWD, 26.04.2019
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_wave_read_namelists

  USE mo_mpi,                ONLY: my_process_is_stdio
  USE mo_namelist,           ONLY: open_nml_output, close_nml_output
  USE mo_nml_annotate,       ONLY: log_nml_settings
  USE mo_time_nml,           ONLY: read_time_namelist
  USE mo_parallel_nml,       ONLY: read_parallel_namelist
  USE mo_run_nml,            ONLY: read_run_namelist
  USE mo_gribout_nml,        ONLY: read_gribout_namelist
  USE mo_io_nml,             ONLY: read_io_namelist
  USE mo_name_list_output_init, ONLY: read_name_list_output_namelists
  USE mo_grid_nml,           ONLY: read_grid_namelist
  USE mo_grid_config,        ONLY: init_grid_configuration
  USE mo_coupling_nml,       ONLY: read_coupling_namelist
  USE mo_extpar_nml,         ONLY: read_extpar_namelist
  USE mo_advection_nml,      ONLY: read_transport_namelist
  USE mo_wave_nml,           ONLY: read_wave_namelist
  USE mo_interpol_nml,       ONLY: read_interpol_namelist
  USE mo_sleve_nml,          ONLY: read_sleve_namelist

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_wave_namelists

CONTAINS

  !---------------------------------------------------------------------
  !>
  SUBROUTINE read_wave_namelists(wave_namelist_filename,shr_namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: wave_namelist_filename
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename

    INTEGER :: tlen
    LOGICAL :: is_stdio

    is_stdio = my_process_is_stdio()
    IF (is_stdio) CALL open_nml_output('NAMELIST_ICON_output_wave')

    CALL read_time_namelist(TRIM(shr_namelist_filename))

    tlen = LEN_TRIM(wave_namelist_filename)
    CALL read_parallel_namelist       (wave_namelist_filename(1:tlen))

    CALL read_run_namelist            (wave_namelist_filename(1:tlen))

    CALL read_io_namelist             (wave_namelist_filename(1:tlen))

    CALL read_name_list_output_namelists (wave_namelist_filename(1:tlen))

    CALL read_grid_namelist           (wave_namelist_filename(1:tlen))
    CALL read_interpol_namelist       (wave_namelist_filename(1:tlen))
    CALL read_sleve_namelist          (wave_namelist_filename(1:tlen))

    CALL init_grid_configuration()

    CALL read_transport_namelist      (wave_namelist_filename(1:tlen))

    CALL read_wave_namelist           (wave_namelist_filename(1:tlen))

    CALL read_extpar_namelist         (wave_namelist_filename(1:tlen))

    CALL read_gribout_namelist        (wave_namelist_filename(1:tlen))

    CALL read_coupling_namelist       (wave_namelist_filename(1:tlen))

    !-----------------------------------------------------------------
    ! Close the file in which all the namelist variables and their
    ! actual values were stored.
    !-----------------------------------------------------------------

    IF (is_stdio) CALL close_nml_output

    ! write an annotate table of all namelist settings to a text file
    IF (is_stdio) CALL log_nml_settings("nml.wave.log")



  END SUBROUTINE read_wave_namelists

END MODULE mo_wave_read_namelists
