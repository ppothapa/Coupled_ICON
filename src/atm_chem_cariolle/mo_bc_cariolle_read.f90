MODULE mo_bc_cariolle_read
USE mo_kind,               ONLY: wp
USE mo_read_interface,     ONLY: read_bcast_REAL_3D, read_1D,  &
                               & closeFile, openInputFile
IMPLICIT NONE
PRIVATE
PUBLIC    :: read_bcast_real_3d_wrap, read_bcast_real_1d_wrap, & 
           & closeFile_wrap, openInputFile_wrap
CONTAINS
SUBROUTINE read_bcast_real_3d_wrap(file_id, variable_name, n1,n2,n3,a_temp)
!, var)

  INTEGER, INTENT(IN)               :: file_id
  CHARACTER(LEN=*), INTENT(IN)      :: variable_name
  INTEGER, INTENT(IN)               :: n1,n2,n3
  REAL(wp), INTENT(INOUT)           :: a_temp(n1,n2,n3)
  REAL(wp), POINTER                 :: return_pointer(:,:,:)
  write(*,*) 'file_id=',file_id, 'variable_name=',TRIM(variable_name), ' varname end'
  write(*,*) '----------','size=',size(a_temp)
  CALL read_bcast_REAL_3D(file_id=file_id, variable_name=variable_name, &
                          return_pointer=return_pointer)
  write(0,*) 'after read: file_id=',file_id,'size=',size(return_pointer)
  write(0,*) 'SIZE return pointer=',SIZE(return_pointer,1),SIZE(return_pointer,2),SIZE(return_pointer,3), return_pointer(3,2,1)
  a_temp=return_pointer
END SUBROUTINE read_bcast_real_3d_wrap
SUBROUTINE read_bcast_real_1d_wrap(file_id, variable_name, n1, a_temp)
  INTEGER, INTENT(IN)               :: file_id
  CHARACTER(LEN=*), INTENT(IN)      :: variable_name
  INTEGER, INTENT(IN)               :: n1
  REAL(wp), INTENT(INOUT)           :: a_temp(n1)
  REAL(wp), POINTER                 :: return_pointer(:)
  CALL read_1D(file_id, variable_name, return_pointer=return_pointer)
  a_temp=return_pointer
END SUBROUTINE read_bcast_real_1d_wrap
SUBROUTINE closeFile_wrap(file_id)
  INTEGER, INTENT(IN)               :: file_id
  CALL closeFile(file_id)
END SUBROUTINE closeFile_wrap
INTEGER FUNCTION openInputFile_wrap(filename)
  CHARACTER(LEN=*), INTENT(IN)      :: filename
  openInputFile_wrap=openInputFile(filename)
END FUNCTION openInputFile_wrap
END MODULE mo_bc_cariolle_read
