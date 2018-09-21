!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!>
!! Preliminary read and time interpolation routine for monthly SST and sea ice data
!! 
!! This is  a clone of the respective ECHAM routine
!!
!! U. Schlese, DKRZ,  May 1993, original version
!! U. Schulzweida, MPI, May 1999, netCDF version
!! L. Kornblueh, MPI, November 2001, cleanup for parallel environment
!! U. Schulzweida, MPI, May 2002, blocking (nproma)
!! L. Kornblueh, MPI, February 2013, adapted as temporary reader in ICON using cdi
!!
!! TODO: ctfreez in echam = 271.38, this is 271.45 K
!
MODULE mo_bc_sst_sic
  
  USE mo_kind,               ONLY: dp, i8
  USE mo_exception,          ONLY: finish, message, message_text
  USE mo_mpi,                ONLY: my_process_is_mpi_workroot
  USE mo_scatter,            ONLY: scatter_time_array
  USE mo_model_domain,       ONLY: t_patch
  USE mo_grid_config,        ONLY: n_dom
  USE mo_parallel_config,    ONLY: nproma
  USE mo_physical_constants, ONLY: tf_salt !, tmelt
  USE mo_impl_constants,     ONLY: MAX_CHAR_LENGTH
  USE mo_cdi,                ONLY: streamOpenRead, streamInqVlist, gridInqSize,      &
    &                              vlistInqTaxis, streamInqTimestep, taxisInqVdate,  &
    &                              vlistInqVarGrid, streamClose, streamReadVarSlice
  USE mo_util_cdi,           ONLY: cdiGetStringError
  USE mo_bcs_time_interpolation, ONLY: t_time_interpolation_weights

  IMPLICIT NONE

  PRIVATE

  TYPE t_ext_sea
    REAL(dp), POINTER :: sst(:,:,:) => NULL()
    REAL(dp), POINTER :: sic(:,:,:) => NULL()
  END TYPE t_ext_sea

  TYPE(t_ext_sea), ALLOCATABLE, TARGET :: ext_sea(:)

  CHARACTER(len=MAX_CHAR_LENGTH) :: sst_fn
  CHARACTER(len=MAX_CHAR_LENGTH) :: sic_fn

  PUBLIC :: read_bc_sst_sic
  PUBLIC :: bc_sst_sic_time_interpolation
  PUBLIC :: get_current_bc_sst_sic_year

  INTEGER(i8), SAVE :: current_year = -1

CONTAINS
  
  SUBROUTINE read_bc_sst_sic(year, p_patch)

    INTEGER(i8),   INTENT(in) :: year
    TYPE(t_patch), INTENT(in) :: p_patch

    INTEGER :: jg

    REAL(dp), POINTER :: zin(:,:) => NULL()

    LOGICAL :: lexist

    jg = p_patch%id

    ! allocate once only module memory for the external SST and SIC data
    IF (.NOT. ALLOCATED (ext_sea        )) ALLOCATE (ext_sea(n_dom))
    IF (.NOT. ASSOCIATED(ext_sea(jg)%sst)) ALLOCATE (ext_sea(jg)%sst(nproma, p_patch%nblks_c, 0:13))
    IF (.NOT. ASSOCIATED(ext_sea(jg)%sic)) ALLOCATE (ext_sea(jg)%sic(nproma, p_patch%nblks_c, 0:13))

    ! allocate temporary field for reading SST and SIC on the whole grid jg
    IF (.NOT. ASSOCIATED(zin)) ALLOCATE(zin(p_patch%n_patch_cells_g, 0:13))

    IF (my_process_is_mpi_workroot()) THEN
   
      IF (n_dom > 1) THEN
        WRITE(sst_fn,'(a,i2.2,a)') 'bc_sst_DOM', jg, '.nc'
      ELSE
        sst_fn = 'bc_sst.nc'
      ENDIF

      WRITE(message_text,'(a,a,a,i4)') &
           'Read SST from ', TRIM(sst_fn), ' for ', year
      CALL message('',message_text)
      
      INQUIRE (file=TRIM(sst_fn), exist=lexist)
      IF (lexist) THEN
        CALL read_sst_sic_data(TRIM(sst_fn), year, zin)
      ELSE
        WRITE (message_text,*) 'Could not open file ',TRIM(sst_fn)
        CALL message('',message_text)
        CALL finish ('mo_bc_sst_sic:read_bc_sst_sic', 'run terminated.')
      ENDIF

    ENDIF

    ! local (i.e. cut SST-data into slices of length nproma)
    CALL scatter_time_array(zin, ext_sea(jg)%sst, p_patch%cells%decomp_info%glb_index)

    ! and the same for SIC (don't allocate zin again)
    IF (my_process_is_mpi_workroot()) THEN

      IF (n_dom > 1) THEN
        WRITE(sic_fn,'(a,i2.2,a)') 'bc_sic_DOM', jg, '.nc'
      ELSE
        sic_fn = 'bc_sic.nc'
      ENDIF

      WRITE(message_text,'(a,a,a,i4)') &
           'Read sea ice from ', TRIM(sic_fn), ' for ', year
      CALL message('',message_text)
      
      INQUIRE (file=TRIM(sic_fn), exist=lexist)
      IF (lexist) THEN
        CALL read_sst_sic_data(TRIM(sic_fn), year, zin)
      ELSE
        WRITE (message_text,*) 'Could not open file ', TRIM(sic_fn)
        CALL message('',message_text)
        CALL finish ('mo_bc_sst_sic:read_bc_sst_sic', 'run terminated.')
      ENDIF
      
    ENDIF

    ! local (i.e. cut SIC-data into slices of length nproma)
    CALL scatter_time_array(zin, ext_sea(jg)%sic, p_patch%cells%decomp_info%glb_index)
    
    DEALLOCATE(zin)
    
    IF (jg==n_dom) current_year = year

  END SUBROUTINE read_bc_sst_sic
  
  SUBROUTINE read_sst_sic_data(fn, y, zin)
    
    CHARACTER(len=*), INTENT(in) :: fn
    INTEGER(i8), INTENT(in) :: y
    REAL(dp), POINTER       :: zin(:,:)
    
    INTEGER :: ngridsize
    INTEGER(i8) :: ym1, yp1
    
    INTEGER :: taxisID
    INTEGER :: vlistID, varID, streamID, tsID
    INTEGER :: nmiss, status, vdate, vyear, vmonth

    REAL(dp), ALLOCATABLE :: buffer(:)

    CHARACTER(LEN=MAX_CHAR_LENGTH) :: cdiErrorText

    ym1 = y-1
    yp1 = y+1

    streamID = streamOpenRead(fn)
    IF ( streamID < 0 ) THEN
      CALL cdiGetStringError(streamID, cdiErrorText)
      WRITE(message_text,*) TRIM(cdiErrorText)
      CALL finish('mo_bc_sst_sic:read_sst_sic_data', message_text)
    END IF
    
    vlistID = streamInqVlist(streamID)
    varID = 0
    ngridsize = gridInqSize(vlistInqVarGrid(vlistID, varID))

    ALLOCATE(buffer(ngridsize))

    taxisID = vlistInqTaxis(vlistID)
    tsID = 0

    DO
      status = streamInqTimestep(streamID, tsID)
      IF ( status == 0 ) EXIT
      vdate = taxisInqVdate(taxisID)
      vyear = vdate/10000
      vmonth = (vdate/100)-vyear*100

      IF (INT(vyear,i8) == ym1 .AND. vmonth == 12) THEN
        CALL streamReadVarSlice(streamID, varID, 0, buffer, nmiss)
        zin(:,0) = buffer(:)
      ELSE IF (INT(vyear,i8) == y) THEN
        CALL streamReadVarSlice(streamID, varID, 0, buffer, nmiss)
        zin(:,vmonth) = buffer(:)
      ELSE IF (INT(vyear,i8) == yp1 .AND. vmonth == 1) THEN
        CALL streamReadVarSlice(streamID, varID, 0, buffer, nmiss)
        zin(:,13) = buffer(:)
        EXIT
      ENDIF
      tsID = tsID+1
    END DO

    CALL streamClose(streamID)
    
    DEALLOCATE(buffer)

  END SUBROUTINE read_sst_sic_data

  SUBROUTINE bc_sst_sic_time_interpolation(tiw, mask_lnd, tsw, seaice, siced, p_patch)
    
    TYPE( t_time_interpolation_weights), INTENT(in) :: tiw
    LOGICAL        , INTENT(in)  :: mask_lnd(:,:)  !< logical land-sea mask, .TRUE. means there is no fraction of ocean/sea-ice 
    REAL(dp)       , INTENT(out) :: tsw(:,:) 
    REAL(dp)       , INTENT(out) :: seaice(:,:) 
    REAL(dp)       , INTENT(out) :: siced(:,:) 
    TYPE(t_patch)  , INTENT(in)  :: p_patch

    REAL(dp) :: zts(SIZE(tsw,1),SIZE(tsw,2))
    REAL(dp) :: zic(SIZE(tsw,1),SIZE(tsw,2))

    INTEGER  :: jg

    jg = p_patch%id

    zts(:,:) = tiw%weight1 * ext_sea(jg)%sst(:,:,tiw%month1_index) + tiw%weight2 * ext_sea(jg)%sst(:,:,tiw%month2_index)
    zic(:,:) = tiw%weight1 * ext_sea(jg)%sic(:,:,tiw%month1_index) + tiw%weight2 * ext_sea(jg)%sic(:,:,tiw%month2_index)

    !TODO: missing siced needs to be added

    ! mask_lnd = .TRUE. means grid box is completely covered by land or lake.
    ! Note that lakes and ocean/sea ice are mutually exclusive, i.e. a cell cannot contain both lake and ocean/sea ice.
    WHERE (mask_lnd(:,:))
      seaice(:,:) = 0._dp
      !TODO: check tsw/i/l sequence,dummy setting to some reasonable value for land and ice
      tsw(:,:) = zts(:,:)
    ELSE WHERE
      seaice(:,:) = zic(:,:)*0.01_dp               ! assuming input data is in percent
      ! seaice(:,:) = MAX(0.0_dp, MIN(0.99_dp, zic(:,:)))
      seaice(:,:) = MERGE(0.99_dp, seaice(:,:), seaice(:,:) > 0.99_dp)
      ! IF (seaice(:,:) <= 0.01_dp) seaice(:,:) = 0.0_dp
      seaice(:,:) = MERGE(0.0_dp, seaice(:,:), seaice(:,:) <= 0.01_dp)
      ! IF (seaice(:,:) > 0.0_dp) THEN           ! ice
      !   tsw(:,:)=tf_salt              
      ! ELSE                                     ! water
      !   tsw(:,:)=MAX(zts(:,:), tf_salt)
      ! END IF
      tsw(:,:) = MERGE(tf_salt, MAX(zts(:,:), tf_salt), seaice(:,:) > 0.0_dp) 
    END WHERE
    
    WHERE (seaice(:,:) > 0.0_dp)
      siced(:,:) = MERGE(2._dp, 1._dp, p_patch%cells%center(:,:)%lat > 0.0_dp)
    ELSEWHERE
      siced(:,:) = 0._dp
    ENDWHERE

    !CALL message('','Interpolated sea surface temperature and sea ice cover.')

  END SUBROUTINE bc_sst_sic_time_interpolation

  FUNCTION get_current_bc_sst_sic_year() RESULT(this_year)
    INTEGER(i8) :: this_year
    this_year = current_year
  END FUNCTION get_current_bc_sst_sic_year

END MODULE mo_bc_sst_sic
