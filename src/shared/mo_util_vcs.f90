!!> Auxiliary module for status information of the version control system (VCS).
!!
!!  Most of the routines which are defined within this module are
!!  ISO-C bindings to a related C program "build/version.c". This
!!  program is automatically generated during ICON's make process.
!!
!!  @par Copyright and License
!!
!!  This code is subject to the DWD and MPI-M-Software-License-Agreement in
!!  its most recent form.
!!  Please see the file LICENSE in the root of the source tree for this code.
!!  Where software is supplied by third parties, it is indicated in the
!!  headers of the routines.
!!
MODULE mo_util_vcs

  USE, INTRINSIC :: iso_c_binding,   ONLY: c_int, c_char, c_null_char
  USE, INTRINSIC :: iso_fortran_env, ONLY: compiler_version
  USE mo_util_sysinfo,               ONLY: util_user_name, util_os_system, util_node_name
  USE mo_cdi,                        ONLY: gribapiLibraryVersion
  USE mo_cf_convention,              ONLY: set_cf_global, cf_global_info
  USE mo_exception,                  ONLY: message_text, message
  USE mo_mpi,                        ONLY: my_process_is_global_root

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: show_version
  PUBLIC :: util_repository_url
  PUBLIC :: util_branch_name
  PUBLIC :: util_revision_key

  INTERFACE

    SUBROUTINE private_util_repository_url(name, actual_len) BIND(c,name='repository_url')
      IMPORT :: c_int, c_char
      CHARACTER(c_char), DIMENSION(*), INTENT(inout) :: name
      INTEGER(c_int), INTENT(inout) :: actual_len
    END SUBROUTINE private_util_repository_url

    SUBROUTINE private_util_branch_name(name, actual_len) BIND(c,name='branch_name')
      IMPORT :: c_int, c_char
      CHARACTER(c_char), DIMENSION(*), INTENT(inout) :: name
      INTEGER(c_int), INTENT(inout) :: actual_len
    END SUBROUTINE private_util_branch_name

    SUBROUTINE private_util_revision_key(name, actual_len) BIND(c,name='revision_key')
      IMPORT :: c_int, c_char
      CHARACTER(c_char), DIMENSION(*), INTENT(inout) :: name
      INTEGER(c_int), INTENT(inout) :: actual_len
    END SUBROUTINE private_util_revision_key

    SUBROUTINE private_util_git_tag_name(name, actual_len) BIND(c,name='git_tag_name')
      IMPORT :: c_int, c_char
      CHARACTER(c_char), DIMENSION(*), INTENT(inout) :: name
      INTEGER(c_int), INTENT(inout) :: actual_len
    END SUBROUTINE private_util_git_tag_name

  END INTERFACE

CONTAINS

  SUBROUTINE util_repository_url(name, actual_len)
    CHARACTER(len=*), INTENT(out) :: name
    INTEGER, INTENT(inout) :: actual_len
    INTEGER :: i
    CALL private_util_repository_url(name, actual_len)
    char_loop: DO i = 1 , LEN(name)
      IF (name(i:i) == c_null_char) EXIT char_loop
    END DO char_loop
    name(i:LEN(name)) = ' '
  END SUBROUTINE util_repository_url

  SUBROUTINE util_branch_name(name, actual_len)
    CHARACTER(len=*), INTENT(out) :: name
    INTEGER, INTENT(inout) :: actual_len
    INTEGER :: i
    CALL private_util_branch_name(name, actual_len)
    char_loop: DO i = 1 , LEN(name)
      IF (name(i:i) == c_null_char) EXIT char_loop
    END DO char_loop
    name(i:LEN(name)) = ' '
  END SUBROUTINE util_branch_name

  SUBROUTINE util_git_tag_name(name, actual_len)
    CHARACTER(len=*), INTENT(out) :: name
    INTEGER, INTENT(inout) :: actual_len
    INTEGER :: i
    CALL private_util_git_tag_name(name, actual_len)
    char_loop: DO i = 1 , LEN(name)
      IF (name(i:i) == c_null_char) EXIT char_loop
    END DO char_loop
    name(i:LEN(name)) = ' '
  END SUBROUTINE util_git_tag_name

  SUBROUTINE util_revision_key(name, actual_len)
    CHARACTER(len=*), INTENT(out) :: name
    INTEGER, INTENT(inout) :: actual_len
    INTEGER :: i
    CALL private_util_revision_key(name, actual_len)
    char_loop: DO i = 1 , LEN(name)
      IF (name(i:i) == c_null_char) EXIT char_loop
    END DO char_loop
    name(i:LEN(name)) = ' '
  END SUBROUTINE util_revision_key


  ! --- VCS version
  !
  ! When using SVN, the revision number is printed as it is
  ! generated by the "config/pvcs.pl" command-line tool during the
  ! make process. This revision number is stored in the file
  ! build/<config.guess triplet>/src/version.c.
  !
  ! Comments on the revision number format: Consider, for example
  !
  !     21000:21099M
  !
  ! "M"           : The working copy has local modifications.
  !
  ! 21000:21099 :  This is a mixed-revision working copy. Not all
  !                files have the same revision number and revision
  !                numbers are within the given range. This usually
  !                happens when a developer makes a checkout of
  !                revision 21000, then changes a single file and
  !                commits, resulting in revision 21099 for this
  !                file. As long as the user does not "svn update",
  !                his working copy has mixed revisions.
  SUBROUTINE show_version()
    CHARACTER(len=256) :: repository  = ''
    CHARACTER(len=256) :: branch      = ''
    CHARACTER(len=256) :: git_tag     = ''
    CHARACTER(len=256) :: revision    = ''
    CHARACTER(len=256) :: executable  = ''
    CHARACTER(len=256) :: user_name   = ''
    CHARACTER(len=256) :: os_name     = ''
    CHARACTER(len=256) :: host_name   = ''
    CHARACTER(len=256) :: tmp_string  = ''
    CHARACTER(len=  8) :: date_string = ''
    CHARACTER(len= 10) :: time_string = ''
    INTEGER :: nlen, nlena, nlenb, nlenc, nlend, &
      &        grb_major_version, grb_minor_version, grb_revision_version

    nlen = 256
    CALL util_repository_url(repository, nlen)
    nlen = 256
    CALL util_branch_name(branch, nlen)
    nlen = 256
    CALL util_revision_key(revision, nlen)
    nlen = 256
    CALL util_git_tag_name(git_tag, nlen)
    
    CALL get_command_argument(0, executable, nlend)
    CALL DATE_AND_TIME(date_string, time_string)
    
    tmp_string = ''
    CALL util_os_system (tmp_string, nlena)
    os_name = tmp_string(1:nlena)
    
    tmp_string = ''
    CALL util_user_name (tmp_string, nlenb)
    user_name = tmp_string(1:nlenb)
    
    tmp_string = ''
    CALL util_node_name (tmp_string, nlenc)
    host_name = tmp_string(1:nlenc)
    
    CALL set_cf_global(title       = 'ICON simulation',                                                    &
      &                institution = 'Max Planck Institute for Meteorology/Deutscher Wetterdienst',        &
      &                source      = TRIM(repository)//'@'//TRIM(revision),                                &
      &                history     = executable(1:nlend)//' at '//date_string(1:8)//' '//time_string(1:6), &
      &                references  = 'see MPIM/DWD publications',                                          &
      &                comment     = TRIM(user_name)//' on '//TRIM(host_name)//' ('//TRIM(os_name)//')')
    
    !-------------------------------------------------------------------
    IF (my_process_is_global_root()) THEN
      
      WRITE(message_text,'(a,a)') 'Repository: ', TRIM(repository)
      CALL message('',message_text)
      WRITE(message_text,'(a,a)') 'Branch    : ', TRIM(branch)
      CALL message('',message_text)
      IF (LEN_TRIM(git_tag) > 0) THEN
         WRITE(message_text,'(a,a)') 'Git tag   : ', TRIM(git_tag)
         CALL message('',message_text)
      END IF
      WRITE(message_text,'(a,a)') 'Revision  : ', TRIM(revision)
      CALL message('',message_text)
      WRITE(message_text,'(a,a)') 'Executable: ', TRIM(cf_global_info%history)
      CALL message('',message_text)
      WRITE(message_text,'(a,a)') 'User      : ', TRIM(cf_global_info%comment)
      CALL message('',message_text)
      
      ! --- Ask CDI for the GRIB version in use
      CALL gribapiLibraryVersion(grb_major_version, grb_minor_version, grb_revision_version)
      IF (grb_major_version > 0) THEN
        WRITE(message_text,'(a,i0,a,i0,a,i0)') 'GRIB-API  : ', grb_major_version, ".", grb_minor_version, ".",&
          &                                    grb_revision_version
        CALL message('',message_text)
      END IF
      
      ! --- Cray compiler: print compiler version in use:
#ifdef _CRAYFTN
      WRITE(message_text,'(a,a)') 'Compiler  : ', _RELEASE_STRING
      ! --- Nvidia compiler: use predefined macros because compiler_version() not available
#elif __NVCOMPILER
      WRITE(message_text,'(a,I2,a,I1)') 'Compiler  : ',__NVCOMPILER_MAJOR__,'.',__NVCOMPILER_MINOR__
#else
      WRITE(message_text,'(a)') compiler_version()
#endif
      CALL message('', message_text)
      
      CALL message('','')
    END IF
  END SUBROUTINE show_version

END MODULE mo_util_vcs
