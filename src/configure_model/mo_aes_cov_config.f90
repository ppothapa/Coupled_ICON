!>
!! Configuration of the parameterization for cloud cover,
!! that is used in the AES physics package.
!!
!! @author Marco Giorgetta, MPI-M
!!
!! @par Revision History
!! First version by Marco Giorgetta, MPI-M (2019-06)
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
MODULE mo_aes_cov_config

  USE mo_exception            ,ONLY: finish, message, print_value
  USE mo_kind                 ,ONLY: wp
  USE mo_impl_constants       ,ONLY: max_dom

  IMPLICIT NONE
  PRIVATE
  PUBLIC ::                     name   !< name for this unit

  ! configuration
  PUBLIC ::         aes_cov_config   !< user specified configuration parameters
  PUBLIC ::    init_aes_cov_config   !< allocate and initialize aes_cov_config
  PUBLIC ::    eval_aes_cov_config   !< evaluate aes_cov_config
  PUBLIC ::   print_aes_cov_config   !< print out

  !>
  !! Name of this unit
  !!
  CHARACTER(LEN=*), PARAMETER :: name = 'aes_cov'
  
  !>
  !! Configuration type containing parameters and switches for the configuration of the AES microphysics
  !!
  TYPE t_aes_cov_config
     !
     ! configuration parameters
     ! ------------------------
     !
     REAL(wp) :: cqx      !          critical mass fraction of cloud water + cloud ice in air, in kg/kg
     !
  END TYPE t_aes_cov_config

  !>
  !! Configuration state vectors, for multiple domains/grids.
  !!
  TYPE(t_aes_cov_config), TARGET :: aes_cov_config(max_dom)
  
CONTAINS

  !----

  !>
  !! Initialize the configuration state vector
  !!
  SUBROUTINE init_aes_cov_config
    !
    ! AES cloud cover configuration
    ! -------------------------------
    !
    aes_cov_config(:)% cqx      = 1.0e-8_wp
    !
  END SUBROUTINE init_aes_cov_config

  !----

  !>
  !! Evaluate additional derived parameters
  !!
  SUBROUTINE eval_aes_cov_config(ng)
    !
    INTEGER, INTENT(in) :: ng
    !
    INTEGER             :: jg
    CHARACTER(LEN=2)    :: cg
    !
    DO jg = 1,ng
       !
       IF (aes_cov_config(jg)% cqx < 0.0_wp .OR. 1.0_wp < aes_cov_config(jg)% cqx) THEN
          WRITE(cg,'(i0)') jg
          CALL finish('eval_aes_cov_config', &
               &      'aes_cov_config('//TRIM(cg)//')% cqx <0 or >1 is not allowed')
       END IF
       !
    END DO
    !
  END SUBROUTINE eval_aes_cov_config

  !----

  !>
  !! Print out the user controlled configuration state
  !!
  SUBROUTINE print_aes_cov_config(ng)
    !
    INTEGER, INTENT(in) :: ng
    !
    INTEGER             :: jg
    CHARACTER(LEN=2)    :: cg
    !
    CALL message    ('','')
    CALL message    ('','========================================================================')
    CALL message    ('','')
    CALL message    ('','AES cloud cover configuration')
    CALL message    ('','=============================')
    CALL message    ('','')
    !
    DO jg = 1,ng
       !
       WRITE(cg,'(i0)') jg
       !
       CALL message    ('','For domain '//cg)
       CALL message    ('','------------')
       CALL message    ('','')
       CALL message    ('','---      0/1 cloud cover based on mass fraction of cloud water + cloud ice in air')
       CALL print_value('    aes_cov_config('//TRIM(cg)//')% cqx      ',aes_cov_config(jg)% cqx     )
       CALL message    ('','')
       !
    END DO
    !
  END SUBROUTINE print_aes_cov_config

  !----

END MODULE mo_aes_cov_config
