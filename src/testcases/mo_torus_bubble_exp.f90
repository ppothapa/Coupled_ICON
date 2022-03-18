!!  Subroutine to initialize a bubble on a torus
!!
!!
!! @par Revision History
!! - first version by Sebastian Rast, MPIM, (2012-01-14)
!! - built on predecessor versions by M. Bergemann  
!! - start from moist profile that contains a temperature and moisture "bubble" on the bottom
!! @par Literature
!! -
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
MODULE mo_torus_bubble_exp
  USE mo_kind,                ONLY: wp
  USE mo_model_domain,        ONLY: t_patch
  USE mo_physical_constants,  ONLY: rd, grav, p0ref,rd_o_cpd, o_m_rdv
  USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag, t_nh_metrics, t_nh_ref
  USE mo_parallel_config,     ONLY: nproma
  USE mo_torus_bubble_exp_config, ONLY: bubble_config
!  USE mo_nh_testcases_nml,    ONLY: rce_tprescr_noise
  USE mo_run_config,          ONLY: iqv
  USE mo_satad,               ONLY: spec_humi, sat_pres_water
  USE mo_hydro_adjust,        ONLY: hydro_adjust_iterative
!  USE mo_random_util,         ONLY: add_random_noise_global

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_nh_state_rce_bubble_glb

  CONTAINS
  
  SUBROUTINE init_nh_state_rce_bubble_glb( ptr_patch, ptr_nh_prog,  ptr_nh_ref, ptr_nh_diag,  &
  &                                 ptr_metrics)

    ! INPUT PARAMETERS:
    TYPE(t_patch),TARGET,  INTENT(IN)   :: &  !< patch on which computation is performed
      &  ptr_patch
    TYPE(t_nh_prog),       INTENT(INOUT):: &  !< prognostic state vector
      &  ptr_nh_prog
    TYPE(t_nh_diag),       INTENT(INOUT):: &  !< diagnostic state vector
      &  ptr_nh_diag
    TYPE(t_nh_metrics),    INTENT(IN)   :: &
      &  ptr_metrics                          !< NH metrics state
    TYPE(t_nh_ref),        INTENT(INOUT):: &  !< reference state vector
      &  ptr_nh_ref

    REAL(wp), ALLOCATABLE :: temp(:,:,:), rh(:,:,:), zeta_xy(:,:)
    REAL(wp)              :: pres(nproma), tt(nproma)
    INTEGER  :: jb,jk,jl  ! loop indices
    INTEGER  :: nblks_c,npromz_c
    INTEGER  :: nlen,nlev
    REAL(wp) :: qv, x, y, z, dz, sat_pres, wat_pres, tv, fwhm, sigma_x, sigma_z, zeta
    REAL(wp), POINTER :: psfc,        t_am,        t0,        gamma0,      &
                       & gamma1,      z0,          t_perturb, relhum_bg,   &
                       & relhum_mx,   hw_frac_x,   hw_z,      x_center,    &
                       & torus_domain_length
    LOGICAL, POINTER  :: lgaussxy

    psfc        => bubble_config%psfc
    t_am        => bubble_config%t_am
    t0          => bubble_config%t0
    gamma0      => bubble_config%gamma0
    gamma1      => bubble_config%gamma1
    z0          => bubble_config%z0
    t_perturb   => bubble_config%t_perturb
    relhum_bg   => bubble_config%relhum_bg
    relhum_mx   => bubble_config%relhum_mx
    hw_frac_x   => bubble_config%hw_frac_x
    hw_z        => bubble_config%hw_z
    x_center    => bubble_config%x_center
    lgaussxy    => bubble_config%lgaussxy
    torus_domain_length => bubble_config%torus_domain_length
    
    ! values for the blocking
    nblks_c  = ptr_patch%nblks_c
    npromz_c = ptr_patch%npromz_c

    ! number of vertical levels
    nlev   = ptr_patch%nlev
    
    ALLOCATE(temp(nproma,nlev,nblks_c))
    ALLOCATE(rh(nproma,nlev,nblks_c))
    ALLOCATE(zeta_xy(nproma,nblks_c))
    
    ! init surface pressure
    ptr_nh_diag%pres_sfc(:,:) = psfc
  
    ! Tracers: all zero by default
    ptr_nh_prog%tracer(:,:,:,:) = 0._wp

    ! Parameters for horizontal Gaussian profile
    fwhm = hw_frac_x * torus_domain_length ! full width half max
    sigma_x = fwhm /  (2._wp * SQRT(2._wp * LOG(2._wp)))! full width half max to sigma_x
    sigma_z = hw_z / (2._wp * SQRT(2._wp * LOG(2._wp)))
    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      END IF
      ! Calculate gaussians in horizontal direction first
      ! Either in x direction only (lgaussxy=.FALSE.) or
      ! in x and y direction with the same sigma_x and x_center
      IF (lgaussxy) THEN
         DO jl = 1, nlen
            x = ptr_patch%cells%cartesian_center(jl, jb)%x(1)
            y = ptr_patch%cells%cartesian_center(jl, jb)%x(2)
            zeta_xy(jl,jb) = gaussian(x_center, sigma_x, 1._wp, x) * &
                           & gaussian(x_center, sigma_x, 1._wp, y)
         END DO
      ELSE
         DO jl = 1, nlen
            x = ptr_patch%cells%cartesian_center(jl, jb)%x(1)
            zeta_xy(jl,jb) = gaussian(x_center, sigma_x, 1._wp, x)
         END DO
      END IF
    END DO
    
    DO jb = 1, nblks_c

      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      END IF

      tt(1:nlen) = t0 + t_perturb * zeta_xy(1:nlen,jb)
      pres(1:nlen) = psfc
      
      DO jk = nlev, 1, -1
        ! Calculation of temperature profile, pressure, Exner pressure, and potential virtual temperature 
        ! We perform a simple integration of the hydrostatic equation from bottom to top of atm.
          DO jl = 1, nlen
            z = ptr_metrics%z_mc(jl, jk, jb)
            IF ( jk == nlev ) THEN
              dz = z
            ELSE
              dz = z-ptr_metrics%z_mc(jl, jk+1, jb)
            END IF
            ! Set a 2D patch of moisture at height 'bub_ver_width' [m] above the
            ! centre of the experiment domain. The max moisture is set to
            ! saturation pressure and decreases exponentially towards 0.
            ! Calculate the temperature profile, (instable in the PBL)
            zeta=zeta_xy(jl,jb)*gaussian(0._wp, sigma_z, 1._wp, z)
            IF ( z <= z0 ) THEN ! we are within the PBL
              temp(jl, jk, jb) = t0 - gamma0 * z
            ELSE
              temp(jl, jk, jb) = t0 - gamma0 * z0 - gamma1 * (z-z0)
            END IF
            temp(jl, jk, jb) = MAX(temp(jl, jk, jb), t_am) + t_perturb * zeta
            sat_pres = sat_pres_water(temp(jl, jk, jb))
            wat_pres = sat_pres * ((relhum_mx-relhum_bg) * zeta + relhum_bg) 
              ! The gaussian is defined such that the third argument being =1 means gaussian(0,sigma_x,1,x)=1
            tt(jl) = 0.5_wp * (tt(jl) + temp(jl, jk, jb))
            pres(jl) = pres(jl) * EXP(-grav/rd/tt(jl) * dz)
            tt(jl) = temp(jl,jk,jb)
            tv                              = tt(jl) / (1._wp - (wat_pres / pres(jl)) * o_m_rdv)
            qv                              = spec_humi(wat_pres, pres(jl))
            rh(jl, jk, jb)                  = wat_pres/sat_pres
            ptr_nh_prog%exner(jl, jk, jb)   = (pres(jl) / p0ref)**rd_o_cpd
            ptr_nh_diag%pres(jl, jk, jb)    = pres(jl)
            ptr_nh_prog%rho(jl, jk, jb)     = pres(jl) / rd / tv
            ptr_nh_prog%theta_v(jl, jk, jb) = tv / ptr_nh_prog%exner(jl, jk, jb) 
            ptr_nh_prog%tracer(jl, jk, jb, iqv) = qv
          END DO !jl

      END DO !jk
!!$       write(0,*) 'exner', MAXVAL(ptr_nh_prog%exner(1:nproma,1:nlev,jb)), MINVAL(ptr_nh_prog%exner(1:nproma,1:nlev,jb)) 
!!$       write(0,*) 'pres', MAXVAL(ptr_nh_diag%pres(1:nproma,1:nlev,jb)), MINVAL(ptr_nh_diag%pres(1:nproma,1:nlev,jb))
    END DO ! jb       

    ! Adjust preliminary profiles to numerics of ICON dynamical core
    ! Relative humidity is constant throughout domain and levels, but needs to be stored in a 3d var.
    CALL hydro_adjust_iterative(                                                                  &
         & ptr_patch,                             ptr_metrics,                                    &
         & temp,                                  rh, ptr_nh_prog%exner,                          &
         & ptr_nh_prog%theta_v,                   ptr_nh_prog%rho,                                &
         & ptr_nh_prog%tracer(:,:,:,iqv),         luse_exner_fg=.TRUE.,                           &
         & opt_exner_lbc=ptr_nh_prog%exner(:,nlev,:)                                              )

    
  !meridional and zonal wind
  ptr_nh_prog%vn = 0._wp
  ptr_nh_ref%vn_ref = ptr_nh_prog%vn

  !vertical wind
  ptr_nh_prog%w = 0._wp
  ptr_nh_ref%w_ref = ptr_nh_prog%w

  END SUBROUTINE init_nh_state_rce_bubble_glb
  FUNCTION gaussian(mu, sigma, delta, x)
    !> @brief
    !> calculate Gaussian profile
    !> @param[in] mu
    !> @param[in] sigma
    !> @param[in] delta
    !> @param[in] x
    !
    REAL(wp), INTENT(in) :: mu    !< expectation of the Gaussian
    REAL(wp), INTENT(in) :: sigma !< std-dev. of Gaussian
    REAL(wp), INTENT(in) :: delta !< max. diviation from 0
    REAL(wp), INTENT(in) :: x     !< x value
    REAL(wp) :: gaussian, xx

    xx = ((x - mu) / sigma)
    gaussian = EXP(-.5_wp * xx * xx) * delta
  END FUNCTION gaussian
  
END MODULE mo_torus_bubble_exp
