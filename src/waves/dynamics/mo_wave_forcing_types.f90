!>
!! Data type defintion for wave forcing data state
!!
!! Defines the data type for storing wave-specific forcing
!! fields.
!!
!! @author Mikhail Dobrynin, DWD
!!
!!
!! @par Revision History
!! Initial Revision by Mikhail Dobrynin, DWD (2023-02-08)
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
MODULE mo_wave_forcing_types

  USE mo_kind,               ONLY: wp


  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: t_wave_forcing

  ! wave-specific forcing data type
  !
  TYPE :: t_wave_forcing

    REAL(wp), POINTER, CONTIGUOUS ::  &
      & u10m (:,:),       & ! zonal wind in 10m (nproma,nblks_c)      ( m/s )
      & v10m (:,:),       & ! meridional wind in 10m (nproma,nblks_c) ( m/s )
      & sp10m(:,:),       & ! wind speed in 10m  (nproma,nblks_c)           ( m/s )
      & dir10m(:,:),      & ! wind direction in 10m (nproma,nblks_c)        ( deg )
      & sea_ice_c(:,:),   & ! sea ice concentration at centers (fraction of 1)
      & sea_ice_e(:,:),   & ! sea ice concentration at edges (fraction of 1)
      & sea_level_c(:,:), & ! sea level height at centers (nproma,nblks_c) ( m )
      & sea_level_e(:,:), & ! sea level height at edges (nproma,nblks_e)   ( m )
      & usoce_c(:,:),     & ! zonal ocean surface current at centers    (nproma,nblks_c) ( m/s )
      & vsoce_c(:,:),     & ! medional ocean surface current at centers (nproma,nblks_c) ( m/s )
      & usoce_e(:,:),     & ! zonal ocean surface current at edges (nproma,nblks_e)    ( m/s )
      & vsoce_e(:,:)      & ! medional ocean surface current at edges (nproma,nblks_e) ( m/s )
      & => NULL()

  END TYPE t_wave_forcing

END MODULE mo_wave_forcing_types
