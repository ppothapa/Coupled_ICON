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

  USE mo_vertical_coord_table ,ONLY: vct_a
  USE mo_physical_constants   ,ONLY: tmelt

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
     ! cloud cover
     INTEGER  :: icov     !          cloud cover scheme
     !                               0:  constant  cloud cover
     !                               1: fractional cloud cover dependent on relative humidity
     !                               2:     0/1    cloud cover dependent on relative humidity
     !                               3:     0/1    cloud cover dependent on cloud condensate
     ! icov=0:
     REAL(wp) :: clcon    !          cloud cover constant, in m2/m2, [0,1]
     !
     ! icov=1 and icov=2:
     REAL(wp) :: csat     !          relative humidity for 100% cloud cover
     !
     ! icov=1:
     REAL(wp) :: crs      !          critical relative humidity at surface
     REAL(wp) :: crt      !          critical relative humidity aloft
     INTEGER  :: nex      !          transition parameter for critical relative humidity profile
     REAL(wp) :: zmaxinv  !          maximum height (m) above sea level for search of inversion layer
     INTEGER  :: jksinv   !          vertical start index for search of inversion layer
     !                               diagnosed in eval_aes_cov_config
     REAL(wp) :: zmininv  !          minimum height (m) above sea level for search of inversion layer
     INTEGER  :: jkeinv   !          vertical end index for search of inversion layer
     !                               diagnosed in eval_aes_cov_config
     REAL(wp) :: cinv     !          fraction of dry adiabatic lapse rate for search of top level of inversion layer over sea
     REAL(wp) :: csatsc   !          lower limit of scaling factor for saturation mixing ratio in layer below inversion
     !                               (csatsc=1 defaults to the standard scheme without accounting of inversion layers)
     ! icov=3:
     REAL(wp) :: cqx      !          critical mass mixing ratio in kg/kg of cloud water + cloud ice
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
    aes_cov_config(:)% icov     = 1
    aes_cov_config(:)% clcon    = 0.0_wp
    aes_cov_config(:)% csat     = 1.0_wp
    aes_cov_config(:)% crs      = 0.968_wp
    aes_cov_config(:)% crt      = 0.8_wp
    aes_cov_config(:)% nex      = 2
    aes_cov_config(:)% zmaxinv  = 2000.0_wp
    aes_cov_config(:)% zmininv  =  200.0_wp
    aes_cov_config(:)% cinv     = 0.25_wp
    aes_cov_config(:)% csatsc   = 0.7_wp
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
    INTEGER             :: jg, jk, klev
    CHARACTER(LEN=2)    :: cg
    !
    klev = SIZE(vct_a)-1
    !
    DO jg = 1,ng
       !
       WRITE(cg,'(i0)') jg
       !
       SELECT CASE (aes_cov_config(jg)% icov)
       CASE (0)
          !
          IF (aes_cov_config(jg)% clcon < 0.0_wp .OR. 1.0_wp < aes_cov_config(jg)% clcon) THEN
             CALL finish('eval_aes_cov_config', &
                  &      'aes_cov_config('//TRIM(cg)//')% clcon <0 or >1 is not allowed')
          END IF
          !
       CASE (1)
          !
          ! diagnose jksinv
          aes_cov_config(jg)% jksinv = 1
          DO jk = 1,klev
             IF ((vct_a(jk)+vct_a(jk+1))*0.5_wp > aes_cov_config(jg)% zmaxinv) THEN
                aes_cov_config(jg)% jksinv = aes_cov_config(jg)% jksinv + 1
             ELSE
                EXIT
             END IF
          END DO
          !
          ! diagnose jkeinv
          aes_cov_config(jg)% jkeinv = klev
          DO jk = klev,1,-1
             IF ((vct_a(jk)+vct_a(jk+1))*0.5_wp < aes_cov_config(jg)% zmininv) THEN
                aes_cov_config(jg)% jkeinv = aes_cov_config(jg)% jkeinv - 1
             ELSE
                EXIT
             END IF
          END DO
          !
          ! check that crs and crt are smaller than csat 
          IF (aes_cov_config(jg)% crs >= aes_cov_config(jg)% csat) THEN
             CALL finish('eval_aes_cov_config', &
                  &      'aes_cov_config('//TRIM(cg)//')% crs >= aes_cov_config('//TRIM(cg)//')% csat is not allowed')
          END IF
          !
          IF (aes_cov_config(jg)% crt >= aes_cov_config(jg)% csat) THEN
             CALL finish('eval_aes_cov_config', &
                  &      'aes_cov_config('//TRIM(cg)//')% crt >= aes_cov_config('//TRIM(cg)//')% csat is not allowed')
          END IF
          !
       CASE (2)
          !
       CASE (3)
          !
          IF (aes_cov_config(jg)% cqx < 0.0_wp .OR. 1.0_wp < aes_cov_config(jg)% cqx) THEN
             CALL finish('eval_aes_cov_config', &
                  &      'aes_cov_config('//TRIM(cg)//')% cqx <0 or >1 is not allowed')
          END IF
          !
       CASE DEFAULT
          !
          CALL finish('eval_aes_cov_config', &
                  &   'aes_cov_config('//TRIM(cg)//')% icov /= 0:3 is not allowed')
       END SELECT
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
       CALL print_value('    aes_cov_config('//TRIM(cg)//')% icov     ',aes_cov_config(jg)% icov    )
       SELECT CASE (aes_cov_config(jg)% icov)
       CASE(0)
          CALL message    ('','---      --> use constant cloud cover')
          CALL print_value('    aes_cov_config('//TRIM(cg)//')% clcon    ',aes_cov_config(jg)% clcon   )
       CASE(1)
          CALL message    ('','---      --> use the fractional cloud cover scheme')
          CALL print_value('    aes_cov_config('//TRIM(cg)//')% csat     ',aes_cov_config(jg)% csat    )
          CALL print_value('    aes_cov_config('//TRIM(cg)//')% crs      ',aes_cov_config(jg)% crs     )
          CALL print_value('    aes_cov_config('//TRIM(cg)//')% crt      ',aes_cov_config(jg)% crt     )
          CALL print_value('    aes_cov_config('//TRIM(cg)//')% nex      ',aes_cov_config(jg)% nex     )
          CALL print_value('    aes_cov_config('//TRIM(cg)//')% zmaxinv  ',aes_cov_config(jg)% zmaxinv )
          CALL print_value('    aes_cov_config('//TRIM(cg)//')% jksinv   ',aes_cov_config(jg)% jksinv  )
          CALL print_value('    aes_cov_config('//TRIM(cg)//')% zmininv  ',aes_cov_config(jg)% zmininv )
          CALL print_value('    aes_cov_config('//TRIM(cg)//')% jkeinv   ',aes_cov_config(jg)% jkeinv  )
          CALL print_value('    aes_cov_config('//TRIM(cg)//')% cinv     ',aes_cov_config(jg)% cinv    )
          CALL print_value('    aes_cov_config('//TRIM(cg)//')% csatsc   ',aes_cov_config(jg)% csatsc  )
       CASE(2)
          CALL message    ('','---      --> use the 0/1 cloud cover scheme with rel. humidity')
          CALL print_value('    aes_cov_config('//TRIM(cg)//')% csat     ',aes_cov_config(jg)% csat    )
       CASE(3)
          CALL message    ('','---      --> use the 0/1 cloud cover scheme with cloud condensate')
          CALL print_value('    aes_cov_config('//TRIM(cg)//')% cqx      ',aes_cov_config(jg)% cqx     )
       END SELECT
       CALL message    ('','')
       !
    END DO
    !
  END SUBROUTINE print_aes_cov_config

  !----

END MODULE mo_aes_cov_config
