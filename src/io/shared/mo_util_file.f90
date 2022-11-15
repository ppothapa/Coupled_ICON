!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
MODULE mo_util_file

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_CHAR, C_NULL_CHAR, C_LONG
  USE mo_exception, ONLY: finish
  USE mo_impl_constants, ONLY: SUCCESS
  USE mo_kind, ONLY: i8

  IMPLICIT NONE

  PRIVATE

  INTERFACE 
    FUNCTION private_symlink(file, link) RESULT(iret) BIND(C,NAME='symlink')
#if defined (__SUNPRO_F95)
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_CHAR
#else
      IMPORT :: C_INT, C_CHAR
#endif
      INTEGER(C_INT) :: iret
#if defined (__SUNPRO_F95)
      CHARACTER(kind=C_CHAR,len=*), INTENT(in) :: file
      CHARACTER(kind=C_CHAR,len=*), INTENT(in) :: link
#else
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: file
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: link
#endif
    END FUNCTION private_symlink
  END INTERFACE
  
  INTERFACE
    FUNCTION private_unlink(filename) RESULT(iret) BIND(C,NAME='unlink')
#if defined (__SUNPRO_F95)
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_CHAR
#else
      IMPORT :: C_INT, C_CHAR
#endif
      INTEGER(C_INT) :: iret
#if defined (__SUNPRO_F95)
     CHARACTER(kind=C_CHAR,len=*), INTENT(in) :: filename
#else
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: filename
#endif
    END FUNCTION private_unlink
  END INTERFACE
    
  INTERFACE
    FUNCTION private_islink(filename) RESULT(iret) BIND(C,NAME='util_islink')
#if defined (__SUNPRO_F95)
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_CHAR
#else
      IMPORT :: C_INT, C_CHAR
#endif
      INTEGER(C_INT) :: iret
#if defined (__SUNPRO_F95)
      CHARACTER(kind=C_CHAR,len=*), INTENT(in) :: filename
#else
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: filename
#endif
    END FUNCTION private_islink
  END INTERFACE

  INTERFACE 
    FUNCTION private_rename(old_filename, new_filename) RESULT(iret) BIND(C,NAME='rename')
#if defined (__SUNPRO_F95)
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_CHAR
#else
      IMPORT :: C_INT, C_CHAR
#endif
      INTEGER(C_INT) :: iret
#if defined (__SUNPRO_F95)
      CHARACTER(kind=C_CHAR,len=*), INTENT(in) :: old_filename
      CHARACTER(kind=C_CHAR,len=*), INTENT(in) :: new_filename
#else
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: old_filename
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: new_filename
#endif
    END FUNCTION private_rename
  END INTERFACE

  INTERFACE
    FUNCTION private_create_tmpfile(filename, max_len) RESULT(flen) BIND(C,NAME='util_create_tmpfile')
#if defined (__SUNPRO_F95)
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_CHAR, C_INT
#else
      IMPORT :: C_CHAR, C_INT
#endif
      INTEGER(C_INT) :: flen
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(inout) :: filename
      INTEGER(C_INT), VALUE, INTENT(in) :: max_len
    END FUNCTION private_create_tmpfile
  END INTERFACE

  INTERFACE
    FUNCTION private_filesize(filename) RESULT(flen) BIND(C,NAME='util_filesize')
#if defined (__SUNPRO_F95)
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_LONG, C_CHAR
#else
      IMPORT :: C_LONG, C_CHAR
#endif
      INTEGER(C_LONG) :: flen
#if defined (__SUNPRO_F95)
      CHARACTER(kind=C_CHAR,len=*), INTENT(in) :: filename
#else
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: filename
#endif
    END FUNCTION private_filesize
  END INTERFACE

  INTERFACE
    FUNCTION private_file_is_writable(filename) RESULT(iwritable) BIND(C,NAME='util_file_is_writable')
#if defined (__SUNPRO_F95)
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_CHAR
#else
      IMPORT :: C_INT, C_CHAR
#endif
      INTEGER(C_INT) :: iwritable
#if defined (__SUNPRO_F95)
      CHARACTER(kind=C_CHAR,len=*), INTENT(in) :: filename
#else
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: filename
#endif
    END FUNCTION private_file_is_writable
  END INTERFACE

  PUBLIC :: util_symlink
  PUBLIC :: util_unlink
  PUBLIC :: util_islink
  PUBLIC :: util_rename
  PUBLIC :: util_tmpnam
  PUBLIC :: util_filesize
  PUBLIC :: util_file_is_writable
  PUBLIC :: createSymlink
  PUBLIC :: get_filename
  PUBLIC :: get_path

  CHARACTER(*), PARAMETER :: modname = "mo_util_file"

CONTAINS

  FUNCTION util_symlink(file, link) RESULT(iret)
    INTEGER :: iret
    CHARACTER(len=*), INTENT(in) :: file
    CHARACTER(len=*), INTENT(in) :: link
    iret = private_symlink(TRIM(file)//C_NULL_CHAR, TRIM(link)//C_NULL_CHAR)
  END FUNCTION util_symlink

  FUNCTION util_unlink(filename) RESULT(iret)
    INTEGER :: iret
    CHARACTER(len=*), INTENT(in) :: filename
    iret = private_unlink(TRIM(filename)//C_NULL_CHAR)
  END FUNCTION util_unlink

  FUNCTION util_islink(filename) RESULT(islink)
    LOGICAL :: islink
    CHARACTER(len=*), INTENT(in) :: filename
    INTEGER :: iret
    iret = private_islink(TRIM(filename)//C_NULL_CHAR)
    islink = .FALSE.
    IF (iret == 1) islink = .TRUE.
  END FUNCTION util_islink

  FUNCTION util_rename(old_filename, new_filename) RESULT(iret)
    INTEGER :: iret
    CHARACTER(len=*), INTENT(in) :: old_filename
    CHARACTER(len=*), INTENT(in) :: new_filename
    iret = private_rename(TRIM(old_filename)//C_NULL_CHAR, TRIM(new_filename)//C_NULL_CHAR)
  END FUNCTION util_rename

  FUNCTION create_tmpfile(filename, max_len) RESULT(flen)
    INTEGER :: flen
    CHARACTER(len=*), INTENT(out) :: filename
    INTEGER,          INTENT(in)  :: max_len
    ! local variables
    INTEGER :: i
    !
    CHARACTER(C_CHAR), ALLOCATABLE :: c_filename(:)
    !
    ALLOCATE(c_filename(max_len + 1))
    flen = private_create_tmpfile(c_filename, max_len + 1) - 1
    DO i = 1, flen
      filename(i:i) = c_filename(i)
    ENDDO
    DEALLOCATE(c_filename)
  END FUNCTION create_tmpfile

  FUNCTION util_tmpnam(filename, klen) RESULT(flen)
    INTEGER :: flen
    CHARACTER(len=*), INTENT(out) :: filename
    INTEGER,          INTENT(in)  :: klen

    flen = create_tmpfile(filename, klen)
    IF (flen < 0) THEN
      CALL finish('mo_util_file::util_tmpnam', 'Failed to find a tmp filename!')
    END IF
  END FUNCTION util_tmpnam

  FUNCTION util_filesize(filename) RESULT(flen)
    INTEGER(KIND=i8) :: flen
    CHARACTER(len=*), INTENT(in) :: filename
    flen = private_filesize(TRIM(filename)//C_NULL_CHAR)
  END FUNCTION util_filesize

  FUNCTION util_file_is_writable(filename) RESULT(lwritable)
    LOGICAL :: lwritable
    CHARACTER(len=*), INTENT(in) :: filename
    lwritable = (private_file_is_writable(TRIM(filename)//C_NULL_CHAR) == 1)
  END FUNCTION util_file_is_writable

  INTEGER FUNCTION createSymlink(targetPath, linkName) RESULT(error)
    CHARACTER(*), INTENT(IN) :: targetPath, linkName
    INTEGER :: i
    CHARACTER(KIND = C_CHAR) :: linkNameCopy(LEN(linkName) + 1), targetPathCopy(LEN(targetPath) + 1)
    INTERFACE
     INTEGER(C_INT) FUNCTION c_createSymlink(c_targetPath, c_linkName) BIND(C, NAME = "createSymlink")
       IMPORT C_INT, C_CHAR
       CHARACTER(KIND = C_CHAR) :: c_targetPath(*), c_linkName(*)
     END FUNCTION c_createSymlink
    END INTERFACE

    DO i = 1, LEN(targetPath)
      targetPathCopy(i) = targetPath(i:i)
    END DO
    targetPathCopy(i) = C_NULL_CHAR
    DO i = 1, LEN(linkName)
      linkNameCopy(i) = linkName(i:i)
    END DO
    linkNameCopy(i) = C_NULL_CHAR
    error = c_createSymlink(targetPathCopy, linkNameCopy)
  END FUNCTION createSymlink

  ! @return subtring from from the first "/" to the end of @p in_str.
  FUNCTION get_filename(in_str) RESULT(out_str)
    ! Parameters
    CHARACTER(len=*), INTENT(IN) :: in_str
    CHARACTER(len=LEN(in_str))   :: out_str
    ! Local parameters
    INTEGER :: pos

    ! start after last occurrence of '/':
    pos = INDEX(in_str,'/',.TRUE.)
    pos = pos + 1
    out_str = in_str(pos:LEN(in_str))
  END FUNCTION get_filename

  ! @return subtring from the begin of @p in_str to the last "/".
  FUNCTION get_path(in_str) RESULT(out_str)
    ! Parameters
    CHARACTER(len=*), INTENT(IN) :: in_str
    CHARACTER(len=LEN(in_str))   :: out_str
    ! Local parameters
    INTEGER :: pos

    ! start after last occurrence of '/':
    out_str = ""
    pos = INDEX(in_str,'/',.TRUE.)
    IF (pos > 0) THEN
      out_str = in_str(1:pos)
    END IF
  END FUNCTION get_path

END MODULE mo_util_file


