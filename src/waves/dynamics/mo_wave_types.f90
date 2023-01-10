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

  IMPLICIT NONE

  PRIVATE


  PUBLIC :: t_wave_diag, t_wave_state, t_wave_state_lists

   ! diagnostic variables state vector
  TYPE t_wave_diag
     REAL(wp), POINTER       &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
   , CONTIGUOUS              &
#endif
     ::                      &
     process_id(:,:)         &
     => NULL()
  END type t_wave_diag

  TYPE t_wave_state
     TYPE(t_wave_diag)                 :: diag
  END TYPE t_wave_state

  TYPE t_wave_state_lists
     TYPE(t_var_list_ptr)              :: diag_list
  END TYPE t_wave_state_lists
  
END MODULE mo_wave_types
