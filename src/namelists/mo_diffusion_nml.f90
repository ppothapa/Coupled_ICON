!>
!! Contains the setup of variables related to horizontal diffusion
!!        
!! @par Revision History
!!   Revision History in mo_global_variables.f90 (r3611)
!!   Modification by Constantin Junk (2011-02-24)
!!     - added new module mo_diffusion_nml
!!     - separated declaration of namelist dynamics_nml from 
!!       mo_global_variables and moved it mo_io_nml
!!     - separated reading of diffusion_nml from subroutine
!!       setup_dynamics and moved it to the new subroutine
!!       setup_diffusion
!!
!! @par Copyright
!! 2002-2006 by DWD and MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!!    violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!    copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!    an according license agreement with DWD and MPI-M.
!! </ol>
!!
!! @par Warranty
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!!
!!
MODULE mo_diffusion_nml

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, message_text, finish, print_value
  USE mo_impl_constants,      ONLY: max_char_length, max_dom, SUCCESS
  USE mo_physical_constants,  ONLY: grav
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_namelist,            ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_master_nml,          ONLY: lrestart
  USE mo_mpi,                 ONLY: p_pe, p_io
  USE mo_diffusion_config,    ONLY: diffusion_config
  USE mo_grid_configuration,  ONLY: global_cell_type
  USE mo_vertical_coord_table,ONLY: vct_a, vct_b, apzero
  USE mo_io_restart_namelist, ONLY: open_tmpfile, store_and_close_namelist,  &
    &                               open_and_restore_namelist, close_tmpfile

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: k2, k4, k6, k2s, k2e, k4s, k4e ! will be removed later KF
  PUBLIC :: read_diffusion_namelist, diffusion_nml_setup

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  !-------------------------------------------------------------------------
  ! 1.0 Namelist variables and auxiliary parameters setting up the
  !     configuration of the dynamical core
  !-------------------------------------------------------------------------
  INTEGER  :: &       ! order of horizontal diffusion
    &  hdiff_order    ! 2: 2nd order linear diffusion on all vertical levels 
                      ! 3: Smagorinsky diffusion for hexagonal model
                      ! 4: 4th order linear diffusion on all vertical levels 
                      ! 5: Smagorinsky diffusion for triangular model
                      ! 24 or 42: 2nd order linear diffusion for upper levels,
                      !           4th order for lower levels
                          

  REAL(wp) :: &       ! (relevant only when hdiff_order = 24 or 42)
    &  k2_pres_max    ! pressure (in Pa) specified by the user
                      ! to determine the lowest vertical level 
                      ! to which 2nd order linear diffusion is applied.
                      ! For the levels with pressure > k2_pres_max, 
                      ! 4th order linear diffusion is applied. 

  INTEGER  :: &       ! (relevant only when hdiff_order = 24 or 42)
    &  k2_klev_max    ! vertical level index specified by the user
                      ! to determine the lowest vertical level 
                      ! to which 2nd order linear diffusion is applied.
                      ! For the levels with k > k2_klev_max, 
                      ! 4th order linear diffusion is applied. 

  REAL(wp) ::                &
    &  hdiff_efdt_ratio,     & ! ratio of e-folding time to (2*)time step
    &  hdiff_min_efdt_ratio, & ! minimum value of hdiff_efdt_ratio (for upper sponge layer)
    &  hdiff_tv_ratio,       & ! the ratio of diffusion coefficient: temp:mom
    &  hdiff_smag_fac,       & ! scaling factor for Smagorinsky diffusion
    &  hdiff_multfac           ! multiplication factor of normalized diffusion coefficient
                               ! for nested domains

  LOGICAL  ::       &
    &  lhdiff_temp, & ! if .TRUE., apply horizontal diffusion to temp.
    &  lhdiff_vn      ! if .TRUE., apply horizontal diffusion to momentum.

  REAL(wp), ALLOCATABLE, DIMENSION(:):: &
    &  k2, &          ! numerical diffusion coefficients
    &  k4, &          ! Values for these parameters are not directly
    &  k6             ! specified by the user, but derived from the ratio 
                      ! between the e-folding time and the model time step
                      ! (hdiff_efdt_ratio above), and the horizontal 
                      ! resolution of the model

  INTEGER  ::    &    ! indices defining to which vertical levels
    &  k2s, k2e, &    ! 2nd and 4th linear diffusion are applied.
    &  k4s, k4e       ! The values are not specified by the user via namelist,
                      ! but determined from k2_klev_max, k2_pres_max
                      ! and the configuration of the vertical coordinate

  NAMELIST /diffusion_nml/  hdiff_order,  hdiff_efdt_ratio,               &
    &                       hdiff_smag_fac,  lhdiff_temp,  lhdiff_vn,     &
    &                       hdiff_tv_ratio,  hdiff_multfac,  k2_klev_max, &
    &                       k2_pres_max,  hdiff_min_efdt_ratio



  CONTAINS
  !>
  !! Initialization of variables for the horizontal diffusion.
  !!
  SUBROUTINE diffusion_nml_setup(i_ndom,parent_id,nlev)

   INTEGER, INTENT(IN) :: parent_id(max_dom-1) !< list of parent ID's
   INTEGER, INTENT(IN) :: i_ndom !< dimension for time level variables
   INTEGER, INTENT(IN) :: nlev

   INTEGER  :: jg, ist, jk, funit
   REAL(wp) :: zpres(nlev+1)

   CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'mo_diffusion_nml/diffusion_nml_setup:'


   !------------------------------------------------------------
   ! 4.0 check the consistency of the parameters
   !------------------------------------------------------------

 !  SELECT CASE(hdiff_order)
 !  CASE(-1)
 !    CALL message(TRIM(routine),'Horizontal diffusion switched off.')
 !  CASE(2,4)
 !    CONTINUE

 !  CASE(3)
 !    IF (global_cell_type==3) CALL finish(TRIM(routine), &
 !    'hdiff_order = 3 invalid for triangular model.')

 !  CASE(5)
 !    IF (global_cell_type==6) CALL finish(TRIM(routine), &
 !    'hdiff_order = 5 invalid for hexagonal model.')

 !  CASE(24,42)
 !    IF (.NOT.(dynamics_config(1)%iequations==IHS_ATM_TEMP)) CALL finish(TRIM(routine), &
 !    'hdiff_order = 24 or 42 only implemented for the hydrostatic atm model')

 !  CASE DEFAULT
 !    CALL finish(TRIM(routine),                     &
 !      & 'Error: Invalid choice of hdiff_order. '// &                
 !      & 'Choose from -1, 2, 3, 4, 5, 24, and 42.')
 !  END SELECT

 !  IF (hdiff_efdt_ratio<=0._wp) THEN
 !    CALL message(TRIM(routine),'No horizontal background diffusion is used')
 !  ENDIF

 !  IF ( lshallow_water ) lhdiff_temp=.FALSE.

    !-----------------------------------------------------------
    ! If using hybrid linear diffusion, set the starting and 
    ! ending vertical level indices for each diffusion order. 
    !-----------------------------------------------------------
    IF ( hdiff_order==24.OR. hdiff_order==42) THEN                                 
                                                                                
      CALL message('','')
      CALL message('----- horizontal diffusion','')

      IF ( k2_pres_max >0._wp) THEN  ! User has specified a pressure value

        CALL print_value('hdiff:  k2_pres_max (Pa) = ', k2_pres_max)

        ! Calculate the pressure values at layer interfaces
        ! assuming surface pressure is apzero.

        zpres(:) = vct_a(:) + vct_b(:)*apzero

        IF ( k2_pres_max <= zpres(1)) THEN
        ! Model does not include mass of the whole atmosphere; User
        ! specified a pressure value located above the model top.
        ! 2nd order diffusion will not be applied. Only 4th order.
        
          k2_klev_max = 0
          CALL print_value('hdiff: ptop (Pa)        = ',vct_a(1))
          CALL message('--- hdiff',' k2_pres_max <= ptop')

        ELSE IF ( k2_pres_max >= zpres(nlev+1)) THEN
        ! User specified a very high pressure. 2nd order diffusion 
        ! will be applied to all vertical levels.

          k2_klev_max = nlev
          CALL print_value('hdiff: pres_sfc (Pa)    = ',zpres(nlev+1))
          CALL message('--- hdiff',' k2_pres_max >= pres_sfc')

        ELSE ! Search for the layer in which k2_pres_max is located.

          DO jk = 1,nlev
            IF (( k2_pres_max > zpres(jk)).AND.( k2_pres_max <= zpres(jk+1))) THEN
               k2_klev_max = jk 
              EXIT
            END IF
          END DO
          CALL print_value('hdiff: half level pressure (-) = ',zpres( k2_klev_max  ))
          CALL print_value('hdiff: half level pressure (+) = ',zpres( k2_klev_max+1))

        END IF ! k2_pres_max
      END IF   ! k2_pres_max >0
      ! If the user didn't specifiy a pressure value, then use the 
      ! default or user-specified level index k2_klev_max.

      k2s = 1
      k2e = k2_klev_max
      k4s = k2_klev_max +1
      k4e = nlev

      ! Inform the user about the configuration

      CALL print_value('hdiff: k2s =',k2s)
      CALL print_value('hdiff: k2e =',k2e)
      CALL print_value('hdiff: k4s =',k4s)
      CALL print_value('hdiff: k4e =',k4e)
      IF (k2e >= k2s) THEN
        WRITE(message_text,'(2(a,i4.4))') '2nd order from level ',k2s,' to ',k2e
        CALL message('--- hdiff',TRIM(message_text))
      END IF
      IF (k4e >= k4s) THEN
        WRITE(message_text,'(2(a,i4.4))') '4nd order from level ',k4s,' to ',k4e
        CALL message('--- hdiff',TRIM(message_text))
      END IF
      CALL message('--- hdiff','------')
      CALL message('','')

    END IF !  hdiff_order==24.OR. hdiff_order==42

    !-----------------------------------------------------
    ! Store the namelist for restart
    !-----------------------------------------------------
    funit = open_tmpfile()
    WRITE(funit,NML=diffusion_nml)
    CALL store_and_close_namelist(funit, 'diffusion_nml')

    ! Write the contents of the namelist to an ASCII file
    IF(p_pe == p_io) WRITE(nnml_output,nml=diffusion_nml)

    !-----------------------------------------------------------------------
    ! calculate diffusion coefficient
    !-----------------------------------------------------------------------
    ALLOCATE( k2(i_ndom),k4(i_ndom),k6(i_ndom), STAT=ist )
    IF (ist/=SUCCESS) CALL finish(TRIM(routine),&                             
      & 'Allocation of k2, k4, k6 failed')                            

    IF ( hdiff_efdt_ratio <= 0._wp) THEN
      k2(:) = 0._wp
      k4(:) = 0._wp
      k6(:) = 0._wp
    ELSE
      k2(1)=1._wp/( hdiff_efdt_ratio*8._wp)
      k4(1)=1._wp/( hdiff_efdt_ratio*64._wp)
      k6(1)=1._wp/( hdiff_efdt_ratio*512._wp)

      DO jg = 2, i_ndom
        k2(jg) = k2(parent_id(jg-1))* hdiff_multfac
        k4(jg) = k4(parent_id(jg-1))* hdiff_multfac
        k6(jg) = k6(parent_id(jg-1))* hdiff_multfac
      ENDDO
    ENDIF

  END SUBROUTINE diffusion_nml_setup
  !-------------

  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Read Namelist for diffusion. 
  !!
  !! This subroutine 
  !! - reads the Namelist for diffusion
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a 
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state (partly)  
  !!
  !! @par Revision History
  !!  by Daniel Reinert, DWD (2011-06-07)
  !!
  SUBROUTINE read_diffusion_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename 
    INTEGER :: istat, funit
    INTEGER :: jg           ! loop index

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_diffusion_nml: read_diffusion_namelist'

    !-----------------------
    ! 1. default settings
    !-----------------------
      lhdiff_temp          = .TRUE.
      lhdiff_vn            = .TRUE.

      hdiff_order          = 4
      hdiff_efdt_ratio     = 1.0_wp
      hdiff_min_efdt_ratio = 1.0_wp
      hdiff_multfac        = 1.0_wp
      hdiff_smag_fac       = 0.15_wp
      hdiff_tv_ratio       = 1.0_wp

      k2_pres_max          = -99.0_wp                                                    
      k2_klev_max          = 0



    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (lrestart) THEN
      funit = open_and_restore_namelist('diffusion_nml')
      READ(funit,NML=diffusion_nml)
      CALL close_tmpfile(funit)
    END IF

    !------------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processors)
    !------------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('diffusion_nml', status=istat)
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, diffusion_nml)
    END SELECT
    CALL close_nml

    !----------------------------------------------------
    ! 4. Sanity check
    !----------------------------------------------------

    SELECT CASE( hdiff_order)
    CASE(-1)
      CALL message(TRIM(routine),'Horizontal diffusion switched off.')
    CASE(2,3,4,5,24,42)
      CONTINUE
    CASE DEFAULT
      CALL finish(TRIM(routine),                         &
        & 'Error: Invalid choice of  hdiff_order. '// &                
        & 'Choose from -1, 2, 3, 4, 5, 24, and 42.')
    END SELECT


    IF ( hdiff_efdt_ratio<=0._wp) THEN
      CALL message(TRIM(routine),'No horizontal background diffusion is used')
    ENDIF

    !----------------------------------------------------
    ! 5. Fill the configuration state
    !----------------------------------------------------

    DO jg= 1,max_dom
      diffusion_config(jg)% hdiff_order          =  hdiff_order
      diffusion_config(jg)% hdiff_efdt_ratio     =  hdiff_efdt_ratio
      diffusion_config(jg)% hdiff_smag_fac       =  hdiff_smag_fac
      diffusion_config(jg)% lhdiff_temp          =  lhdiff_temp
      diffusion_config(jg)% lhdiff_vn            =  lhdiff_vn
      diffusion_config(jg)% hdiff_tv_ratio       =  hdiff_tv_ratio
      diffusion_config(jg)% hdiff_multfac        =  hdiff_multfac
      diffusion_config(jg)% k2_klev_max          =  k2_klev_max
      diffusion_config(jg)% k2_pres_max          =  k2_pres_max
      diffusion_config(jg)% hdiff_min_efdt_ratio =  hdiff_min_efdt_ratio
    ENDDO

    !-----------------------------------------------------
    ! 6. Store the namelist for restart
    !-----------------------------------------------------
    funit = open_tmpfile()
    WRITE(funit,NML=diffusion_nml)                    
    CALL store_and_close_namelist(funit,'diffusion_nml') 

    ! 7. write the contents of the namelist to an ASCII file
    !
    IF(p_pe == p_io) WRITE(nnml_output,nml=diffusion_nml)


  END SUBROUTINE read_diffusion_namelist

END MODULE mo_diffusion_nml
