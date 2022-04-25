!>
!! Configuration of the parameterization for NWP graupel microphysics,
!! that is used in the AES physics package.
!!
!! @author Monika Esch, MPI-M
!!
!! @par Revision History
!! First version by Monika Esch, MPI-M (2018-06)
!!
!! Based on earlier codes of:
!!     ...
!!
!! References: 
!!     ...
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_aes_mig_config

  USE mo_exception            ,ONLY: message, print_value
  USE mo_kind                 ,ONLY: wp
  USE mo_impl_constants       ,ONLY: max_dom
  USE mo_grid_config          ,ONLY: n_dom

  IMPLICIT NONE
  PRIVATE
  PUBLIC ::                     name   !< name for this unit

  ! configuration
  PUBLIC ::         aes_mig_config   !< user specified configuration parameters
  PUBLIC ::    init_aes_mig_config   !< allocate and initialize aes_mig_config
!!$  PUBLIC ::    eval_aes_mig_config   !< evaluate aes_mig_config
  PUBLIC ::   print_aes_mig_config   !< print out

  !>
  !! Name of this unit
  !!
  CHARACTER(LEN=*), PARAMETER :: name = 'aes_mig'
  
  !>
  !! Configuration type containing parameters and switches for the configuration of the AES microphysics
  !!
  TYPE t_aes_mig_config
     !
     ! configuration parameters
     ! ------------------------
     !
     ! thresholds
     REAL(wp) :: qi0            ! cloud ice threshold for autoconversion
     REAL(wp) :: qc0            ! cloud water threshold for autoconversion
     !
     ! grid scale microphysics
     REAL(wp) :: zceff_min      ! minimum value for sticking efficiency
     REAL(wp) :: v0snow         ! factor in the terminal velocity for snow
     REAL(wp) :: zvz0i          ! terminal fall velocity of ice
     REAL(wp) :: icesedi_exp    ! exponent for density correction of cloud ice sedimentation
     REAL(wp) :: mu_rain        ! parameter in gamma distribution for rain
     REAL(wp) :: rain_n0_factor ! tuning factor for intercept parameter of raindrop size distribution
     !
     LOGICAL  :: ldiag_ttend
     LOGICAL  :: ldiag_qtend
     !
  END TYPE t_aes_mig_config

  !>
  !! Configuration state vectors, for multiple domains/grids.
  !!
  TYPE(t_aes_mig_config), TARGET :: aes_mig_config(max_dom)
  
CONTAINS

  !----

  !>
  !! Initialize the configuration state vector
  !!
  SUBROUTINE init_aes_mig_config
    !
    ! Graupel microphyiscs configuration
    ! --------------------------------------
    !
    ! general thresholds
    aes_mig_config(:)% qi0            = 0.0_wp
    aes_mig_config(:)% qc0            = 0.0_wp
    !
    ! grid scale microphysics
    aes_mig_config(:)% zceff_min      = 0.01_wp
    aes_mig_config(:)% v0snow         = 25.0_wp
    aes_mig_config(:)% zvz0i          = 1.25_wp      ! original value of Heymsfield+Donner 1990: 3.29
    aes_mig_config(:)% icesedi_exp    = 0.33_wp
    aes_mig_config(:)% mu_rain        = 0.0_wp
    aes_mig_config(:)% rain_n0_factor = 1.0_wp
    !
    aes_mig_config(:)% ldiag_ttend    = .FALSE.
    aes_mig_config(:)% ldiag_qtend    = .FALSE.
    !
  END SUBROUTINE init_aes_mig_config

  !----

!!$  !>
!!$  !! Evaluate additional derived parameters
!!$  !!
!!$  SUBROUTINE eval_aes_mig_config
!!$    !
!!$    ...
!!$    !
!!$  END SUBROUTINE eval_aes_mig_config

  !----

  !>
  !! Print out the user controlled configuration state
  !!
  SUBROUTINE print_aes_mig_config
    !
    INTEGER           :: jg
    CHARACTER(LEN=2)  :: cg
    !
    CALL message    ('','')
    CALL message    ('','========================================================================')
    CALL message    ('','')
    CALL message    ('','Graupel microphyiscs configuration')
    CALL message    ('','==================================')
    CALL message    ('','')
    !
    DO jg = 1,n_dom
       !
       WRITE(cg,'(i0)') jg
       !
       CALL message    ('','For domain '//cg)
       CALL message    ('','------------')
       CALL message    ('','')
       CALL print_value('    aes_mig_config('//TRIM(cg)//')% qi0            ',aes_mig_config(jg)% qi0            )
       CALL print_value('    aes_mig_config('//TRIM(cg)//')% qc0            ',aes_mig_config(jg)% qc0            )
       CALL message    ('','')
       CALL print_value('    aes_mig_config('//TRIM(cg)//')% zceff_min      ',aes_mig_config(jg)% zceff_min      )
       CALL print_value('    aes_mig_config('//TRIM(cg)//')% v0snow         ',aes_mig_config(jg)% v0snow         )
       CALL print_value('    aes_mig_config('//TRIM(cg)//')% zvz0i          ',aes_mig_config(jg)% zvz0i          )
       CALL print_value('    aes_mig_config('//TRIM(cg)//')% icesedi_exp    ',aes_mig_config(jg)% icesedi_exp    )
       CALL print_value('    aes_mig_config('//TRIM(cg)//')% mu_rain        ',aes_mig_config(jg)% mu_rain        )
       CALL print_value('    aes_mig_config('//TRIM(cg)//')% rain_n0_factor ',aes_mig_config(jg)% rain_n0_factor )
       CALL message    ('','')
       CALL print_value('    aes_mig_config('//TRIM(cg)//')% ldiag_ttend    ',aes_mig_config(jg)% ldiag_ttend    )
       CALL print_value('    aes_mig_config('//TRIM(cg)//')% ldiag_qtend    ',aes_mig_config(jg)% ldiag_qtend    )
       CALL message    ('','')
       !
    END DO
    !
  END SUBROUTINE print_aes_mig_config

  !----

END MODULE mo_aes_mig_config
