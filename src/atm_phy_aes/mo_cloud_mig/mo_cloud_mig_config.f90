!>
!! Configuration of the parameterization for cloud microphysics "graupel",
!! from DWD that is used in the sapphire physics package.
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
MODULE mo_cloud_mig_config

  USE mo_exception            ,ONLY: message, print_value
  USE mo_kind                 ,ONLY: wp
  USE mo_impl_constants       ,ONLY: max_dom

  USE mo_cloud_mig_types      ,ONLY: t_cloud_mig_config

  IMPLICIT NONE
  PRIVATE

  ! configuration
  PUBLIC ::         cloud_mig_config   !< user specified configuration parameters
  PUBLIC ::    init_cloud_mig_config   !< allocate and initialize cloud_mig_config
!!$  PUBLIC ::    eval_cloud_mig_config   !< evaluate cloud_mig_config
  PUBLIC ::   print_cloud_mig_config   !< print out

  !>
  !! Configuration state vectors, for multiple domains/grids.
  !!
  TYPE(t_cloud_mig_config), TARGET :: cloud_mig_config(max_dom)
  
CONTAINS

  !----

  !>
  !! Initialize the configuration state vector
  !!
  SUBROUTINE init_cloud_mig_config
    !
    ! Graupel microphysics configuration
    ! --------------------------------------
    !
    ! no parameter settings available
    !
    ! grid scale microphysics
    !
  END SUBROUTINE init_cloud_mig_config

  !----

!!$  !>
!!$  !! Evaluate additional derived parameters
!!$  !!
!!$  SUBROUTINE eval_cloud_mig_config
!!$    !
!!$    ...
!!$    !
!!$  END SUBROUTINE eval_cloud_mig_config

  !----

  !>
  !! Print out the user controlled configuration state
  !!
  SUBROUTINE print_cloud_mig_config(ng)
    !
    INTEGER, INTENT(in) :: ng
    !
    INTEGER             :: jg
    CHARACTER(LEN=2)    :: cg
    !
    CALL message    ('','')
    CALL message    ('','========================================================================')
    CALL message    ('','')
    CALL message    ('','Cloud microphysics "graupel" configuration')
    CALL message    ('','==========================================')
    CALL message    ('','')
    !
    DO jg = 1,ng
       !
       WRITE(cg,'(i0)') jg
       !
       !
       ! kept as example... nothing yet available...
       !
       CALL message    ('','For domain '//cg)
       CALL message    ('','------------')
       CALL message    ('','')
       !CALL print_value('    cloud_mig_config('//TRIM(cg)//')% zceff_min      ',cloud_mig_config(jg)% zceff_min      )
       !CALL print_value('    cloud_mig_config('//TRIM(cg)//')% v0snow         ',cloud_mig_config(jg)% v0snow         )
       !CALL print_value('    cloud_mig_config('//TRIM(cg)//')% zvz0i          ',cloud_mig_config(jg)% zvz0i          )
       !CALL print_value('    cloud_mig_config('//TRIM(cg)//')% icesedi_exp    ',cloud_mig_config(jg)% icesedi_exp    )
       !CALL print_value('    cloud_mig_config('//TRIM(cg)//')% mu_rain        ',cloud_mig_config(jg)% mu_rain        )
       !CALL print_value('    cloud_mig_config('//TRIM(cg)//')% rain_n0_factor ',cloud_mig_config(jg)% rain_n0_factor )
       CALL message    ('','')
       !
    END DO
    !
  END SUBROUTINE print_cloud_mig_config

  !----

END MODULE mo_cloud_mig_config
