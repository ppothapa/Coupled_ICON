!>
!! Contains output routines for CDI output
!! This module acts merely as a wrapper for either calling the direct
!! or the asynchronous output routines.
!!
!!
!! @par Revision History
!! Initial implementation by Rainer Johanni (2010-12-02)
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
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
MODULE mo_output

  USE mo_exception,           ONLY: message_text, get_filename_noext !, finish
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH
  USE mo_kind,                ONLY: wp
  USE mo_mpi,                 ONLY: p_pe, p_io, process_mpi_io_size
  USE mo_io_units,            ONLY: filename_max
  USE mo_model_domain_import, ONLY: n_dom, &
    &                               n_dom_start !, nroot, lplane
  USE mo_io_config,           ONLY: out_expname, no_output
  USE mo_impl_constants,      ONLY: ihs_ocean !,            &
!     &                              ihs_atm_temp,         &
!     &                              ihs_atm_theta,        &
!     &                              inh_atmosphere,       &
!     &                              ishallow_water
  USE mo_dynamics_config,     ONLY: iequations, nold, nnow, nnew, nnew_rcf, nnow_rcf 
  USE mo_io_vlist,            ONLY: setup_vlist, destruct_vlist,           &
     &                              open_output_vlist, close_output_vlist, &
     &                              write_vlist
  USE mo_io_async,            ONLY: setup_io_procs, shutdown_io_procs, &
     &                              output_async, set_output_file
  USE mo_datetime,            ONLY: t_datetime,iso8601
!--------------------------------------------------------------------------------------------------
  USE mo_io_output,           ONLY: set_output_time_vl    => set_output_time,     &
                                    init_output_vl        => init_output,         &
     &                              open_output_files_vl  => open_output_files,   &
     &                              close_output_files_vl => close_output_files,  &
     &                              write_output_vl       => write_output
!--------------------------------------------------------------------------------------------------
  USE mo_io_restart,          ONLY: set_restart_time, set_restart_vct,         &
                                  & init_restart, open_writing_restart_files,  &
                                  & write_restart, close_writing_restart_files,&
                                  & finish_restart, set_restart_depth,         &
                                  & set_restart_height
  USE mo_io_restart_attributes,ONLY: set_restart_attribute
  USE mo_model_domain,        ONLY: t_patch, p_patch
  USE mo_interpolation,       ONLY: t_lon_lat_intp
  USE mo_run_config,          ONLY: ltimer
  USE mo_timer,               ONLY: timer_start, timer_stop,&
    &                               timer_write_output !, timer_write_restart_file
  USE mo_meteogram_output,    ONLY: meteogram_flush_file
  USE mo_meteogram_config,    ONLY: meteogram_output_config

  USE mo_oce_state,           ONLY: set_zlev
  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_output'

  LOGICAL :: l_omit_dom   ! flag if "'_DOM',jg" should be omitted in the filenames

  ! RJ: The following parameter is for experimentally enabling output through var_lists
  ! This is not fully functional yet and thus not implemented in a namelist variable
  LOGICAL, PARAMETER :: use_var_lists = .FALSE.

  !------------------------------------------------------------------------------------------------
  !
  ! Public routines:

  PUBLIC :: init_output_files, close_output_files, write_output
  PUBLIC :: create_restart_file

CONTAINS

  !------------------------------------------------------------------------------------------------
  !>
  !! Initialize output file(s)
  !! 
  !! This can be done for the first  or any other time:
  !! * Set the file name
  !! * Add current namelist values as global attributes (netcdf)
  !! * Create data variable definitions
  !! For jfile=1 some additional refinement initializations are done.

  SUBROUTINE init_output_files(jfile, lclose)

    INTEGER, INTENT(IN) :: jfile !> Number of fileset to open
    LOGICAL, INTENT(IN) :: lclose !> lclose old file

    INTEGER :: jg, jlev
    INTEGER :: nlev              !< number of full levels
    CHARACTER(LEN=filename_max) :: outputfile
    CHARACTER(LEN=filename_max) :: grid_filename

    IF ( no_output ) RETURN
    
IF(use_var_lists) THEN
    IF(.NOT.lclose) CALL init_output_vl
    IF(lclose) CALL close_output_files_vl
    CALL open_output_files_vl(jfile)
ELSE
    IF(.NOT.lclose) THEN

      ! This is the first call - initialize

      ! If n_dom=1, i.e. if no grid refinement is used on the sphere, then
      ! do not USE "'_DOM',jg" in output file name, if it is not present
      ! in the input file name. Use l_omit_dom to check if this is the case.
      l_omit_dom = .FALSE.

!       IF (lplane) THEN
!         gridtype='plan'
!       ELSE
!         gridtype='icon'
!       END IF

      DO jg = 1, n_dom

!         jlev = p_patch(jg)%level

        ! Grid file name(s) for input
        !
        ! Allow file names without "DOM" specifier if n_dom=1.
!         IF (n_dom == 1) THEN
!           ! Check if file name without "DOM" specifier exists.
!           WRITE (gridfile(jg),'(a,a,i0,a,i2.2,a)') &
!             &    TRIM(gridtype),'R',nroot,'B',jlev,'-grid.nc'
!           INQUIRE (FILE=gridfile(jg), EXIST=l_omit_dom)
!           ! Otherwise use file name with "DOM" specifier
!           IF (.NOT. l_omit_dom)                                            &
!             &    WRITE (gridfile(jg),'(a,a,i0,2(a,i2.2),a)')                &
!             &    TRIM(gridtype),'R',nroot,'B',jlev,'_DOM',jg,'-grid.nc'
!         ELSE
!           ! n_dom >1 --> "'_DOM',jg" required in file name
!           WRITE (gridfile(jg),'(a,a,i0,2(a,i2.2),a)') &
!             &    TRIM(gridtype),'R',nroot,'B',jlev,'_DOM',jg,'-grid.nc'
!         ENDIF

        ! Set up vlist for this grid level
        ! Please note: setup_vlist only sets up the vlist, it does not open any output file!
        ! The third parameter has to be set to .TRUE. if the current
        ! task actually does I/O for the patch in question.
        ! Compare to the call of open_output_vlist below!

        CALL setup_vlist( TRIM(p_patch(jg)%grid_filename), jg, &
                        & process_mpi_io_size==0 .AND. p_pe==p_io)

      ENDDO

    ELSE

      ! If not called for the first time, close previous output files
      ! (only if we are actually doing output!)
      IF(process_mpi_io_size == 0 .AND. p_pe == p_io) THEN
        DO jg = n_dom, 1, -1
          CALL close_output_vlist(jg)
        ENDDO
      ENDIF

    ENDIF

    DO jg = 1, n_dom

      jlev = p_patch(jg)%level
      nlev = p_patch(jg)%nlev
      grid_filename = get_filename_noext(p_patch(jg)%grid_filename)
      ! Raw data file name(s) for output
      !
      WRITE (outputfile,'(a,a,a,a,i4.4,a)')  &
        &  TRIM(out_expname),"_",TRIM(grid_filename),'_', jfile, '.nc'
      
!       SELECT CASE (iequations)
!         !
!       CASE (ishallow_water)
!         IF (l_omit_dom) THEN
!           WRITE (outputfile,'(a,a,i0,a,i2.2,a,i4.4,a)')  &
!             &  TRIM(out_expname), '_R', nroot, 'B', jlev, '_', jfile, '.nc'
!         ELSE
!           WRITE (outputfile,'(a,a,i2.2,a,i0,a,i2.2,a,i4.4,a)')  &
!             &  TRIM(out_expname), '_DOM', jg, '_R', nroot, 'B', jlev, '_', jfile, '.nc'
!         END IF
!         !
!       CASE (ihs_atm_temp, ihs_atm_theta, inh_atmosphere)
!         IF (l_omit_dom) THEN
!           WRITE (outputfile,'(a,a,i0,a,i2.2,a,i0,a,i4.4,a)')  &
!             &  TRIM(out_expname), '_R', nroot, 'B', jlev, 'L', nlev, '_', jfile, '.nc'
!         ELSE
!           WRITE (outputfile,'(a,a,i2.2,a,i0,a,i2.2,a,i0,a,i4.4,a)')  &
!             &  TRIM(out_expname), '_DOM', jg, '_R', nroot, 'B', jlev, 'L', nlev, '_', jfile, '.nc'
!         END IF
!       CASE (ihs_ocean)
!         WRITE (outputfile,'(a,a,i0,a,i2.2,a,i0,a,i4.4,a)')  &
!           &  TRIM(out_expname), '_O.R', nroot, 'B', jlev, 'L', n_zlev, '_', jfile, '.nc'
!         WRITE(*,'(a,a)') ' control_model: Initial output file for setup_vlist_oce is ', &
!           &              TRIM(outputfile)
!         !
!       CASE DEFAULT
!         CALL finish(modname,'Unsupported value of iequations in init_output_files')
!         !
!       END SELECT

    ! WRITE(0,'(a,a)') ' Initial output file for setup_vlist is ', &
    !     &              TRIM(outputfile)
     !IF (iequations == ihs_ocean) THEN
     !  ! #slo# must be aligned with general output
     !  !CALL setup_vlist_oce( p_patch(1:), TRIM(p_patch(jg)%grid_filename), TRIM(outputfile), jg )
     !ELSE
        IF(process_mpi_io_size == 0) THEN
          IF(p_pe == p_io) CALL open_output_vlist(TRIM(outputfile), jg)
        ELSE
          CALL set_output_file(outputfile, jg)
        ENDIF
     !ENDIF

    ENDDO

    ! Setup I/O PEs if this is the initial call and I/O PEs are enabled
    ! Note that this has to be done AFTER the output files are set!

    IF(jfile == 1 .AND. process_mpi_io_size>0) CALL setup_io_procs()
ENDIF
    
  END SUBROUTINE init_output_files

  !------------------------------------------------------------------------------------------------
  !>
  !! Closes output files and finalizes I/O setup
  !! Note: This routine must only be called for the final close, not when the output files
  !! are switched during the run!

  SUBROUTINE close_output_files

    INTEGER jg

    IF ( no_output ) RETURN

IF(use_var_lists) THEN
    CALL close_output_files_vl
ELSE
    DO jg = n_dom, 1, -1
      IF(process_mpi_io_size == 0 .AND. p_pe == p_io) CALL close_output_vlist(jg)
      CALL destruct_vlist( jg )
    ENDDO

    IF(process_mpi_io_size>0) CALL shutdown_io_procs
ENDIF

  END SUBROUTINE close_output_files

  !------------------------------------------------------------------------------------------------
  !>
  SUBROUTINE write_output(datetime, z_sim_time)

    TYPE(t_datetime),      INTENT(in) :: datetime
    REAL(wp), OPTIONAL, INTENT(in) :: z_sim_time(n_dom)

    ! Local variables
    INTEGER :: jg
    CHARACTER*32 date_iso8601

!    Proposal by Matthias Raschendorfer for correct output
!
!    INTEGER sec
!
!    sec=NINT(datetime%second)-INT(datetime%second)
!    outptime=datetime
!    IF (sec.NE.0) THEN
!       CALL add_time(REAL(sec,wp),0,0,0,outptime)
!    END IF
!
!    IF(process_mpi_io_size == 0) THEN
!      CALL write_vlist(outptime)
!    ELSE
!      CALL output_async(outptime)
!    ENDIF
!
    IF ( no_output ) RETURN

    IF (ltimer) CALL timer_start(timer_write_output)

IF(use_var_lists) THEN
    WRITE(date_iso8601,'(i4.4,i2.2,i2.2,"T",i2.2,i2.2,i2.2,"Z")') &
     & datetime%year, datetime%month, datetime%day,               &
     & datetime%hour, datetime%minute, NINT(datetime%second)
    CALL set_output_time_vl(date_iso8601)
    CALL write_output_vl
ELSE
    IF ( PRESENT(z_sim_time) ) THEN  
      IF(process_mpi_io_size == 0) THEN
        CALL write_vlist(datetime, z_sim_time(1))
        ! write recent samples of meteogram output
        DO jg = 1, n_dom
          IF (meteogram_output_config(jg)%lenabled) THEN
            CALL meteogram_flush_file(jg)
          END IF
        END DO
      ELSE
        CALL output_async(datetime,z_sim_time(1))
      ENDIF
    ELSE
      IF(process_mpi_io_size == 0) THEN
        CALL write_vlist(datetime)
        ! write recent samples of meteogram output
        DO jg = 1, n_dom
          IF (meteogram_output_config(jg)%lenabled) THEN
            CALL meteogram_flush_file(jg)
          END IF
        END DO
      ELSE
        CALL output_async(datetime)
      ENDIF
    ENDIF

ENDIF

    IF (ltimer) CALL timer_stop(timer_write_output)
  END SUBROUTINE write_output

  !-------------
  !>
  !! 
  !! Hui Wan (MPI-M, 2011-05)
  !!
  SUBROUTINE create_restart_file( patch, datetime,            &
                                & jfile, l_have_output,       &
                                & opt_pvct,                   &
                                & opt_t_elapsed_phy,          &
                                & opt_lcall_phy, opt_sim_time,&
                                & opt_jstep_adv_ntsteps,      &
                                & opt_jstep_adv_marchuk_order,&
                                & opt_depth, opt_zheight,     &
                                & opt_zheight_mc ,opt_zheight_ifc )

    TYPE(t_patch),   INTENT(IN) :: patch
    TYPE(t_datetime),INTENT(IN) :: datetime
    INTEGER, INTENT(IN) :: jfile  ! current output file index
    LOGICAL, INTENT(IN) :: l_have_output

    REAL(wp), INTENT(IN), OPTIONAL :: opt_pvct(:)
    INTEGER,  INTENT(IN), OPTIONAL :: opt_depth

    REAL(wp), INTENT(IN), OPTIONAL :: opt_t_elapsed_phy(:,:)
    LOGICAL , INTENT(IN), OPTIONAL :: opt_lcall_phy(:,:)
    REAL(wp), INTENT(IN), OPTIONAL :: opt_sim_time
    INTEGER,  INTENT(IN), OPTIONAL :: opt_jstep_adv_ntsteps
    INTEGER,  INTENT(IN), OPTIONAL :: opt_jstep_adv_marchuk_order
    INTEGER,  INTENT(IN), OPTIONAL :: opt_zheight
    REAL(wp), INTENT(IN), OPTIONAL :: opt_zheight_mc (:,:,:) 
    REAL(wp), INTENT(IN), OPTIONAL :: opt_zheight_ifc(:,:,:)

    INTEGER :: klev, jg, kcell, kvert, kedge, icelltype, izlev, i
    REAL(wp), ALLOCATABLE :: zlevels_full(:), zlevels_half(:)


    CHARACTER(LEN=132) :: string
    CHARACTER(len=MAX_CHAR_LENGTH) :: attname   ! attribute name
    INTEGER :: jp, jp_end   ! loop index and array size

    !IF (ltimer) CALL timer_start(timer_write_restart_file)
    !----------------
    ! Initialization
    klev      = patch%nlev
    jg        = patch%id
    kcell     = patch%n_patch_cells_g
    kvert     = patch%n_patch_verts_g
    kedge     = patch%n_patch_edges_g
    icelltype = patch%cell_type

    CALL set_restart_attribute( 'current_caltime', datetime%caltime )
    CALL set_restart_attribute( 'current_calday' , datetime%calday )

    CALL set_restart_attribute( 'nold'    , nold    (jg))
    CALL set_restart_attribute( 'nnow'    , nnow    (jg))
    CALL set_restart_attribute( 'nnew'    , nnew    (jg))
    CALL set_restart_attribute( 'nnow_rcf', nnow_rcf(jg))
    CALL set_restart_attribute( 'nnew_rcf', nnew_rcf(jg))

    !----------------
    ! additional restart-output for nonhydrostatic model
    IF (PRESENT(opt_sim_time)) THEN
      WRITE(attname,'(a,i2.2)') 'sim_time_DOM',jg
      CALL set_restart_attribute( TRIM(attname), opt_sim_time )
    ENDIF

    !-------------------------------------------------------------
    ! DR
    ! WORKAROUND FOR FIELDS WHICH NEED TO GO INTO THE RESTART FILE,
    ! BUT SO FAR CANNOT BE HANDELED CORRECTLY BY ADD_VAR OR
    ! SET_RESTART_ATTRIBUTE
    !-------------------------------------------------------------

    IF (PRESENT(opt_jstep_adv_ntsteps)) THEN
        WRITE(attname,'(a,i2.2)') 'jstep_adv_ntsteps_DOM',jg
        CALL set_restart_attribute( TRIM(attname), opt_jstep_adv_ntsteps )
    ENDIF

    IF (PRESENT(opt_jstep_adv_marchuk_order)) THEN
        WRITE(attname,'(a,i2.2)') 'jstep_adv_marchuk_order_DOM',jg
        CALL set_restart_attribute( TRIM(attname), opt_jstep_adv_marchuk_order )
    ENDIF

    IF (PRESENT(opt_t_elapsed_phy) .AND. PRESENT(opt_lcall_phy)) THEN
      ! Inquire array size
      jp_end = SIZE(opt_t_elapsed_phy,2)
      DO jp = 1, jp_end
        WRITE(attname,'(a,i2.2,a,i2.2)') 't_elapsed_phy_DOM',jg,'_PHY',jp
        CALL set_restart_attribute( TRIM(attname), opt_t_elapsed_phy(jg,jp) )

        WRITE(attname,'(a,i2.2,a,i2.2)') 'lcall_phy_DOM',jg,'_PHY',jp
        CALL set_restart_attribute( TRIM(attname), opt_lcall_phy(jg,jp) )
      ENDDO
    ENDIF


    IF (l_have_output) THEN
      CALL set_restart_attribute( 'next_output_file', jfile+1 )
    ELSE
      CALL set_restart_attribute( 'next_output_file', jfile   )
    END IF

    IF (PRESENT(opt_pvct)) CALL set_restart_vct( opt_pvct )  ! Vertical coordinate (A's and B's)
    IF (PRESENT(opt_zheight)) THEN                           ! geometrical height for NH 
!      CALL set_restart_height(opt_zheight_ifc ,opt_zheight_mc)
!DR start preliminary fix
      ALLOCATE(zlevels_full(opt_zheight))
      ALLOCATE(zlevels_half(opt_zheight+1))
      DO i = 1, opt_zheight
        zlevels_full(i) = REAL(i,wp)
      END DO
      DO i = 1, opt_zheight+1
        zlevels_half(i) = REAL(i,wp)
      END DO
      CALL set_restart_height(zlevels_half, zlevels_full)
      DEALLOCATE(zlevels_full)
      DEALLOCATE(zlevels_half)
!DR end preliminary fix
    ENDIF
    IF (PRESENT(opt_depth)) THEN                              ! Ocean depth
      izlev = opt_depth
      ALLOCATE(zlevels_full(izlev))
      ALLOCATE(zlevels_half(izlev+1))
      CALL set_zlev(zlevels_half, zlevels_full)
      CALL set_restart_depth(zlevels_half, zlevels_full)
      DEALLOCATE(zlevels_full)
      DEALLOCATE(zlevels_half)
    ELSE
      izlev = 0
    END IF

    CALL init_restart( TRIM(out_expname), &! exp name
                     & '1.2.2',           &! model version
                     & kcell, icelltype,  &! total # of cells, # of vertices per cell
                     & kvert, 9-icelltype,&! total # of vertices, # of vertices per dual cell
                     & kedge, 4,          &! total # of cells, shape of control volume for edge 
                     & klev,              &! total # of vertical layers
                     & izlev,             &! total # of depths below sea
                     & 0)                  ! total numbers of geometric height above ground

    CALL set_restart_time( iso8601(datetime) )  ! Time tag

    ! Open new file, write data, close and then clean-up.
    message_text = get_filename_noext(patch%grid_filename)
    WRITE(string,'(a,a)') 'restart.',TRIM(message_text)

    CALL open_writing_restart_files( TRIM(string) )

#ifdef NOMPI
    CALL write_restart
#else
    CALL write_restart( patch )
#endif

    CALL close_writing_restart_files
    CALL finish_restart

    !IF (ltimer) CALL timer_stop(timer_write_restart_file)

  END SUBROUTINE create_restart_file

END MODULE mo_output
