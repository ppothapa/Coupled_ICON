MODULE mo_util_sysinfo

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_CHAR

  IMPLICIT NONE

  PRIVATE

  INTERFACE
    SUBROUTINE private_util_user_name(name, name_len) BIND(C,NAME='util_user_name') 
      IMPORT :: C_INT, C_CHAR
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(out) :: name
      INTEGER(C_INT), INTENT(out) :: name_len
    END SUBROUTINE private_util_user_name
  END INTERFACE

  INTERFACE
    SUBROUTINE private_util_os_system(name, name_len) BIND(C,NAME='util_os_system') 
      IMPORT :: C_INT, C_CHAR
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(out) :: name
      INTEGER(C_INT), INTENT(out) :: name_len
    END SUBROUTINE private_util_os_system
  END INTERFACE

  INTERFACE
    SUBROUTINE private_util_node_name(name, name_len) BIND(C,NAME='util_node_name') 
      IMPORT :: C_INT, C_CHAR
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(out) :: name
      INTEGER(C_INT), INTENT(out) :: name_len
    END SUBROUTINE private_util_node_name
  END INTERFACE

  PUBLIC :: util_user_name
  PUBLIC :: util_os_system
  PUBLIC :: util_node_name

CONTAINS

  SUBROUTINE util_user_name(name, name_len)
    CHARACTER(len=*), INTENT(out) :: name
    INTEGER, INTENT(out) :: name_len
    CALL private_util_user_name(name, name_len)
  END SUBROUTINE util_user_name

  SUBROUTINE util_os_system(name, name_len)
    CHARACTER(len=*), INTENT(out) :: name
    INTEGER, INTENT(out) :: name_len
    CALL private_util_os_system(name, name_len)
  END SUBROUTINE util_os_system

  SUBROUTINE util_node_name(name, name_len)
    CHARACTER(len=*), INTENT(out) :: name
    INTEGER, INTENT(out) :: name_len
    CALL private_util_node_name(name, name_len)
  END SUBROUTINE util_node_name

END MODULE mo_util_sysinfo
