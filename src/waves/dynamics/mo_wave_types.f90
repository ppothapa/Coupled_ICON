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

  USE mo_kind,                ONLY: wp, vp
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
      &  gv_c(:,:,:),         & ! group velocity                    (nproma,nblks_c,nfreqs)  (m/s)
      &  gv_e(:,:,:),         & ! group velocity                    (nproma,nblks_e,nfreqs)  (m/s)
      &  gvn_e(:,:,:),        & ! orthogonal normal group velocity  (nproma,nblks_e,nfreqs)  (m/s)
      &  gvt_e(:,:,:),        & ! tangential group velocity         (nproma,nblks_e,nfreqs)  (m/s)
      &  alphaj(:,:),         & ! jonswap alpha                     (nproma,nblks_c)
      &  fp(:,:),             & ! jonswap peak frequency            (nproma,nblks_c)         (hz)
      &  et(:,:,:),           & ! jonswap spectra                   (nproma,nblks_c,nfreqs)
      &  flminfr(:,:,:),      & ! the minimum value in spectral bins for a given frequency (nproma,nblks_c,nfreqs)
      &  emean(:,:),          & ! total energy                      (nproma,nblks_c)         (m^2)
      &  emeanws(:,:),        & ! total windsea energy              (nproma,nblks_c)         (m^2)
      &  femean(:,:),         & ! mean frequency energy             (nproma,nblks_c)         (m^2)
      &  f1mean(:,:),         & ! mean frequency based on f-moment  (nproma,nblks_c)
      &  tm1(:,:),            & ! wave tm1 period                   (nproma,nblks_c)         (s)
      &  tm2(:,:),            & ! wave tm2 period                   (nproma,nblks_c)         (s)
      &  wave_num_c(:,:,:),   & ! wave number at cell centers as a function of
                                ! circular frequency and water depth (nproma,nblks_c,nfreqs) (1/m)
      &  wave_num_e(:,:,:),   & ! wave number at cell edges as a function of
                                ! circular frequency and water depth (nproma,nblks_e,nfreqs) (1/m)
      &  akmean(:,:),         & ! mean wavenumber based on sqrt(1/k)-moment (nproma,nblks_c)
      &  xkmean(:,:),         & ! mean wavenumber based on sqrt(k)-moment   (nproma,nblks_c)
      &  femeanws(:,:),       & ! mean windsea frequency energy             (nproma,nblks_c)
      &  ustar(:,:),          & ! friction velocity                 (nproma,nblks_c)         (m/s)
      &  z0(:,:),             & ! roughness length                  (nproma,nblks_c)         (m)
      &  tauhf1(:,:),         & ! init high-frequency stress        (nproma,nblks_c)
      &  phihf1(:,:),         & ! init high-frequency energy flux into ocean (nproma,nblks_c)
      &  tauhf(:,:),          & ! high-frequency stress                 (nproma,nblks_c)
      &  phihf(:,:),          & ! high-frequency energy flux into ocean (nproma,nblks_c)
      &  xlevtail(:,:),       & ! tail level                            (nproma,nblks_c)
      &  tauw(:,:),           & ! wave stress                           (nproma,nblks_c)     (m/s)^2
      &  phiaw(:,:),          & ! energy flux from wind into waves integrated over the full frequency range  (nproma,nblks_c)
      &  fl(:,:,:),           & ! diagonal matrix of functional derivative (nproma,nblks_c,ntracer)
      &  sl(:,:,:),           & ! total source function                    (nproma,nblks_c,ntracer)
      &  hs(:,:),             & ! significant wave height                  (nproma,nblks_c)  (m)
      &  hs_dir(:,:),         & ! mean direction of total wave height      (nproma,nblks_c)  (deg)
      &  enh(:,:),            & ! nonlinear transfer function coefficients for shallow water
      &  AF11(:),             & ! for discrete approximation of nonlinear transfer
      &  FKLAP(:), FKLAP1(:), & ! --//--
      &  FKLAM(:), FKLAM1(:), & ! --//--
      &  u_stokes(:,:),       & ! U-component of surface Stokes drift
      &  v_stokes(:,:)        & ! V-component of surface Stokes drift
      &  => NULL()

    REAL(vp), POINTER, CONTIGUOUS :: &
      &  geo_bath_grad_c(:,:,:,:) & ! bathymetry geographical gradient (2,nproma,nlev,nblks_c) (m/m)
      &  => NULL()

    INTEGER, POINTER, CONTIGUOUS ::  &
      &  last_prog_freq_ind(:,:), & ! last frequency index of the prognostic range
      &  llws(:,:,:),             & ! .TRUE. where sinput is positive, wind sea, swell splitting flag  (nproma,nblks_c,ntracer)
      &  ikp(:), ikp1(:),         & ! for discrete approximation of nonlinear transfer
      &  ikm(:), ikm1(:),         & ! --//--
      &  k1w(:,:), k2w(:,:),      & ! --//--
      &  k11w(:,:), k21w(:,:),    & ! --//--
      &  ja1(:,:), ja2(:,:),      & ! --//--
      &  non_lin_tr_ind(:,:,:,:)  & ! tracer index for nonlinear interaction p_diag%non_lin_tr_ind(nfreqs+4,2,ndirs,8)
      &  => NULL()

    TYPE(t_ptr_2d3d), ALLOCATABLE :: freq_ptr(:)       !< pointer array: one pointer for each frequence
    TYPE(t_ptr_2d3d), ALLOCATABLE :: wave_num_c_ptr(:) !< pointer array: one pointer for each frequence
    TYPE(t_ptr_2d3d), ALLOCATABLE :: wave_num_e_ptr(:) !< pointer array: one pointer for each frequence
    TYPE(t_ptr_2d3d), ALLOCATABLE :: tr_ptr(:)         !< pointer array: one pointer for each tracer
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
