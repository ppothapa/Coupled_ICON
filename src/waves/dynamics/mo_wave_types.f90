!>
!!        Contains the types to set up the wave model.
!=============================================================================================
!!
!! @author Mikhail Dobrynin, DWD, 11.06.2019
!!
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
!=============================================================================================
MODULE mo_wave_types

  USE mo_kind,                ONLY: wp
  USE mo_var_list,            ONLY: t_var_list_ptr
  USE mo_fortran_tools,       ONLY: t_ptr_2d3d

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_wave_prog
  PUBLIC :: t_wave_diag
  PUBLIC :: t_wave_state
  PUBLIC :: t_wave_state_lists

  TYPE t_wave_prog
    REAL(wp), POINTER, CONTIGUOUS :: &
    tracer(:,:,:,:) => NULL()
    !! wave energy (spectral bins) over frequencies and directions (nproma, nlev, nblks_c, ntracer) [m^2 ?]
    TYPE(t_ptr_2d3d), ALLOCATABLE :: tracer_ptr(:) !< pointer array: one pointer for each tracer
  END TYPE t_wave_prog


  ! diagnostic variables state vector
  TYPE t_wave_diag
    REAL(wp), POINTER, CONTIGUOUS :: &
      &  process_id(:,:),         &
      &  u_10m (:,:),             & ! zonal wind in 10m  (nproma,nblks_c)           ( m/s )
      &  v_10m (:,:),             & ! meridional wind in 10m (nproma,nblks_c)       ( m/s )
      &  sp_10m(:,:),             & ! wind speed in 10m  (nproma,nblks_c)           ( m/s )
      &  dir_10m(:,:),            & ! wind direction in 10m (nproma,nblks_c)        ( deg )
      &  ALPHAJ(:,:),     & ! JONSWAP ALPHA (nproma,nblks_c)
      &  FP(:,:),         & ! JONSWAP PEAK FREQUENCY (nproma,nblks_c)               ( Hz )
      &  ET(:,:,:),       & ! JONSWAP SPECTRA (nproma,nblks_c,nfreqs)
      &  FLMINFR(:,:,:),  & ! THE MINIMUM VALUE IN SPECTRAL BINS FOR A GIVEN FREQUENCY (nproma,nblks_c,nfreqs)
      &  emean(:,:),      & ! TOTAL ENERGY (nproma,nblks_c)                         ( m^2 )
      &  emeanws(:,:),    & ! TOTAL WINDSEA ENERGY (nproma,nblks_c)                 ( m^2 )
      &  femean(:,:),     & ! MEAN FREQUENCY ENERGY (nproma,nblks_c)
      &  f1mean(:,:),     & ! MEAN FREQUENCY BASED ON F-MOMENT (nproma,nblks_c)
      &  femeanws(:,:),   & ! MEAN WINDSEA FREQUENCY ENERGY (nproma,nblks_c)
      &  ustar(:,:),      & ! friction velocity  (nproma,nblks_c)                   ( m/s )
      &  z0(:,:),         & ! roughness length   (nproma,nblks_c)                   ( m )
      &  tauw(:,:),       & ! wave stress        (nproma,nblks_c)                   ( (m/s)^2 )
      &  fl(:,:,:),       & ! DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE (nproma,nblks_e,ndirs*nfreqs)
      &  sl(:,:,:),       & ! TOTAL SOURCE FUNCTION                    (nproma,nblks_e,ndirs*nfreqs)
      &  Hs(:,:),         & ! significant wave height             (nproma,nblks_c)  ( m )
      &  Hs_dir(:,:)      & ! mean direction of total wave height (nproma,nblks_c)  ( deg )
      &  => NULL()

    INTEGER, POINTER ::        &
      &  LLWS(:,:,:)           & ! TRUE WHERE SINPUT IS POSITIVE, wind sea, swell splitting flag  (nproma,nblks_c,ntracer)
      &  => NULL()
  END type t_wave_diag

  TYPE t_wave_state
    !array of prognostic states at different timelevels
    TYPE(t_wave_prog), ALLOCATABLE    :: prog(:)       !< shape: (timelevels)
    TYPE(t_wave_diag)                 :: diag
  END TYPE t_wave_state

  TYPE t_wave_state_lists
    ! array of prognostic state lists at different timelevels
    TYPE(t_var_list_ptr), ALLOCATABLE :: prog_list(:)  !< shape: (timelevels)
    TYPE(t_var_list_ptr), ALLOCATABLE :: tracer_list(:)
    TYPE(t_var_list_ptr)              :: diag_list
  END TYPE t_wave_state_lists

END MODULE mo_wave_types
