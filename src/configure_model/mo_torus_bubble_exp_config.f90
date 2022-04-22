! This code is subject to the DWD and MPI-M-Software-License-Agreement in           
!! its most recent form.                                                             
!! Please see the file LICENSE in the root of the source tree for this code.         
!! Where software is supplied by third parties, it is indicated in the               
!! headers of the routines.                                                          
!!                                                                                   
MODULE mo_torus_bubble_exp_config                                                           
                                                                                     
  USE mo_exception            ,ONLY: message, message_text, print_value, warning, finish                                                                                 
  USE mo_kind                 ,ONLY: wp                                              
  USE mo_parallel_config      ,ONLY: nproma                                          
                                                                                     
  IMPLICIT NONE                                                                      
                                                                                     
  PRIVATE                                                                            
                                                                                     
  ! configuration                                                                    
  PUBLIC ::         bubble_config             !< user specified configuration parameters    
  PUBLIC ::    init_torus_bubble_exp_config   !< allocate and initialize bubble_config   
  PUBLIC ::    eval_torus_bubble_exp_config   !< evaluate bubble_config                  
  PUBLIC ::   print_torus_bubble_exp_config   !< print out                                  
 ! This code is subject to the DWD and MPI-M-Software-License-Agreement in           
!! its most recent form.                                                             
!! Please see the file LICENSE in the root of the source tree for this code.         
!! Where software is supplied by third parties, it is indicated in the               
!! headers of the routines.                                                          
!!                                                                                   

  TYPE t_bubble_config
     ! For a detailed description, see cr2021_08_03_jsr
     REAL(wp) :: psfc      !surface temperature
     REAL(wp) :: t0        !temperature at bottom of atmosphere (maybe different from surface temp)
     REAL(wp) :: t_am      !absolute minimum of temperature (default: 180K)
     REAL(wp) :: gamma0    !lapse rate in lowest atmospheric part
     REAL(wp) :: z0        !below z0 lapse rate gamma0 is applied, above gamma1
     REAL(wp) :: gamma1    !lapse rate above z0. However, temperature cannot fall below t_am
     REAL(wp) :: t_perturb ! maximum temperature perturbation in center of Gaussians
     REAL(wp) :: relhum_bg ! background relative humidity
     REAL(wp) :: relhum_mx ! maximum relative humidity
     REAL(wp) :: hw_frac_x ! half width in x-direction as a fraction of domain length of torus
     REAL(wp) :: hw_z      ! half width in z-direction in metres
     REAL(wp) :: x_center  ! deviation of maximum of Gaussian in x direction
     LOGICAL  :: lgaussxy  ! .TRUE., if half width calculated for x-direction and x_center is applied also to y direction
     REAL(wp) :: torus_domain_length ! domain length of torus (diameter of region with doubly periodic boundary conditions)
  END TYPE t_bubble_config

  TYPE(t_bubble_config), TARGET :: bubble_config

CONTAINS

  SUBROUTINE init_torus_bubble_exp_config
    bubble_config%psfc       = 101325._wp
    bubble_config%t_am       = 180._wp
    bubble_config%t0         = 298.15_wp
    bubble_config%gamma0     = 0.006_wp
    bubble_config%z0         = 6000._wp
    bubble_config%gamma1     = 0.009_wp
    bubble_config%t_perturb  = 10._wp
    bubble_config%relhum_bg  = 0.7_wp
    bubble_config%relhum_mx  = 0.95_wp
    bubble_config%hw_frac_x  = 0.2_wp
    bubble_config%hw_z       = 200._wp
    bubble_config%x_center   = 0._wp
    bubble_config%lgaussxy   = .FALSE.
    bubble_config%torus_domain_length = 0._wp
  END SUBROUTINE init_torus_bubble_exp_config

  SUBROUTINE eval_torus_bubble_exp_config
    CHARACTER(LEN=*), PARAMETER :: routine='eval_torus_bubble_exp_config'
    REAL(wp), POINTER :: t_am,        psfc,      t0,           gamma0,      &
                       & gamma1,      z0,        t_perturb,    relhum_bg,   &
                       & relhum_mx,   hw_frac_x, hw_z,          x_center,   &
                       & torus_dl
    t_am        => bubble_config%t_am
    psfc        => bubble_config%psfc
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
    torus_dl    => bubble_config%torus_domain_length

    CALL message    ('','')                                                          
    CALL message    ('','------------------------------------------------------------------------')
    CALL message    ('','')                                                          
    CALL message    ('','Effective input to the torus bubble experiment')                          
    CALL message    ('','================================')                          
    CALL message    ('','')                                                          
 
    IF (t_am <= 0._wp) THEN
       WRITE (message_text,*) 'Absolute minimum temperature in atmosphere t_am=',t_am, &
            & ' but cannot be below or equal to 0K'
       CALL finish(routine,message_text)
    ELSE
       WRITE (message_text,*) 'Absolute minimum temperature in atmosphere t_am=',t_am
       CALL message('',message_text)
    END IF
    ! t0 can be arbitrary, because the formula sets it to values >= t_am
    IF (psfc <= 0._wp) THEN
       WRITE (message_text,*) 'Surface pressure psfc=',psfc,' but cannot be below or equal to zero'
       CALL finish(routine,message_text)
    ELSE
       WRITE (message_text,*) 'Surface pressure psfc=',psfc
       CALL message('',message_text)
    END IF
    IF (gamma0 < 0._wp) THEN
       WRITE (message_text,*) 'Lapse rate gamma0=',gamma0,' but cannot be below 0'
       CALL finish(routine,message_text)
    ELSE
       WRITE (message_text,*) 'Lapse rate gamm0=',gamma0
       CALL message('',message_text)
    END IF
    IF (gamma1 < 0._wp) THEN
       WRITE (message_text,*) 'Lapse rate gamma1=',gamma1,' but cannot be below 0'
       CALL finish(routine,message_text)
    ELSE
       WRITE (message_text,*) 'Lapse rate gamm0=',gamma1
       CALL message('',message_text)
    END IF
    IF (z0 <= 0._wp) THEN
       WRITE (message_text,*) 'Altitude z0=',z0,' but cannot be below 0'
       CALL finish(routine,message_text)
    ELSE
       WRITE (message_text,*) 'Altitude z0=',z0
       CALL message('',message_text)
    END IF
    IF (t_perturb < 0._wp) THEN
       WRITE (message_text,*) 'Temperature perturbation in atmosphere t_perturb=',t_perturb, &
            & ' but cannot be below 0K'
       CALL finish(routine,message_text)
    ELSE
       WRITE (message_text,*) 'Temperature perturbation in atmosphere t_perturb=',t_perturb
       CALL message('',message_text)
    END IF
    IF (relhum_bg < 0._wp .or. relhum_bg > 1._wp) THEN
       WRITE (message_text,*) 'Background relative humidity relhum_bg=',relhum_bg, &
            ' but must be in [0,1]'
       CALL finish(routine,message_text)
    ELSE
       WRITE (message_text,*) 'Background relative humidity relhum_bg=',relhum_bg
       CALL message('',message_text)
    END IF
    IF (relhum_mx < 0._wp .or. relhum_mx > 1._wp) THEN
       WRITE (message_text,*) 'Relative humidity at maximum of Gaussians relhum_mx=', &
            & relhum_mx,' but must be in [0,1]'
       CALL finish(routine,message_text)
    ELSE
       WRITE (message_text,*) 'Relative humidity at maximum of Gaussians relhum_mx=',relhum_mx
       CALL message('',message_text)
    END IF
    IF (hw_frac_x < 0._wp .or. hw_frac_x > 1._wp) THEN
       WRITE (message_text,*) 'Fraction of domain covered by half width of Gaussian hw_frac_x=', &
            & hw_frac_x,' but must be in [0,1]'
       CALL finish(routine,message_text)
    ELSE
       WRITE (message_text,*) 'Fraction of domain covered by half width of Gaussian hw_frac_x=',hw_frac_x
       CALL message('',message_text)
    END IF
    IF (hw_z < 0._wp) THEN
       WRITE (message_text,*) 'Half width in vertical direction hw_z=',hw_z,' but cannot be below 0'
       CALL finish(routine,message_text)
    ELSE
       WRITE (message_text,*) 'Half width in vertical direction hw_z=',hw_z
       CALL message('',message_text)
    END IF
    IF (torus_dl <= 0._wp) THEN
       WRITE (message_text,*) 'Torus domain length torus_domain_length=',torus_dl,' but it must be larger than 0'
       CALL finish(routine,message_text)
    ELSE
       WRITE (message_text,*) 'Torus domain length torus_domain_length=',torus_dl
       CALL message('',message_text)
    END IF
    CALL message ('','')
  END SUBROUTINE eval_torus_bubble_exp_config
  SUBROUTINE print_torus_bubble_exp_config
    CALL message    ('','========================================================================')                                                                      
    CALL message    ('','')                                                          
    CALL message    ('','Torus bubble experiment configuration')                             
    CALL message    ('','=============================')                             
    CALL message    ('','')                                                          

    CALL print_value ('    bubble_config%t_am       ', bubble_config%t_am      )
    CALL print_value ('    bubble_config%psfc       ', bubble_config%psfc      )
    CALL print_value ('    bubble_config%t0         ', bubble_config%t0        )
    CALL print_value ('    bubble_config%gamma0     ', bubble_config%gamma0    )
    CALL print_value ('    bubble_config%gamma1     ', bubble_config%gamma1    )
    CALL print_value ('    bubble_config%z0         ', bubble_config%z0        )
    CALL print_value ('    bubble_config%t_perturb  ', bubble_config%t_perturb )
    CALL print_value ('    bubble_config%relhum_bg  ', bubble_config%relhum_bg )
    CALL print_value ('    bubble_config%relhum_mx  ', bubble_config%relhum_mx )
    CALL print_value ('    bubble_config%hw_frac_x  ', bubble_config%hw_frac_x )
    CALL print_value ('    bubble_config%hw_z       ', bubble_config%hw_z      )
    CALL print_value ('    bubble_config%x_center   ', bubble_config%x_center  )
    CALL print_value ('    bubble_config%lgaussxz   ', bubble_config%lgaussxy  )
    CALL print_value ('    bubble_config%torus_domain_length',bubble_config%torus_domain_length)
    
  END SUBROUTINE print_torus_bubble_exp_config
END MODULE mo_torus_bubble_exp_config
