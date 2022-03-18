!>                                                                                   
!! Read configuration parameters as Fortran namelist from an external file.          
!!                                                                                   
!! @author Marco, Giorgetta, MPI-M, 2017-12                                          
!!                                                                                   
!! @par Revision History                                                             
!!                                                                                   
!! @par Copyright and License                                                        
!!                                                                                   
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in           
!! its most recent form.                                                             
!! Please see the file LICENSE in the root of the source tree for this code.         
!! Where software is supplied by third parties, it is indicated in the               
!! headers of the routines.                                                          
!!                                                                                   
MODULE mo_torus_bubble_exp_nml

  USE mo_torus_bubble_exp_config, ONLY: bubble_config, init_torus_bubble_exp_config, &
       & eval_torus_bubble_exp_config, print_torus_bubble_exp_config
  USE mo_process_nml,             ONLY: process_nml
  
  IMPLICIT NONE                                                                      
  PRIVATE                                                                            
  PUBLIC :: process_torus_bubble_exp_nml                                                    
                                                                                     
  NAMELIST /torus_bubble_exp_nml/ bubble_config                             
                                                                                     
CONTAINS                                                                             
                                                                                     
  SUBROUTINE process_torus_bubble_exp_nml(filename)                                         
    !                                                                                
    CHARACTER(LEN=*), INTENT(in) :: filename                                         
    !                                                                                
    CALL init_torus_bubble_exp_config                                                       
    !                                                                                
    CALL process_nml(filename, 'torus_bubble_exp_nml', nml_read, nml_write)                 
    !                                                                                
    CALL eval_torus_bubble_exp_config
    CALL print_torus_bubble_exp_config
  CONTAINS                                                                           
    !                                                                                
    SUBROUTINE nml_read(funit)                                                       
      INTEGER, INTENT(in) :: funit                                                   
      READ(funit, NML=torus_bubble_exp_nml)                                                 
    END SUBROUTINE nml_read                                                          
    !                                                                                
    SUBROUTINE nml_write(funit)                                                      
      INTEGER, INTENT(in) :: funit                                                   
      WRITE(funit, NML=torus_bubble_exp_nml)                                                
    END SUBROUTINE nml_write                                                         
    !                                                                                
  END SUBROUTINE process_torus_bubble_exp_nml
                                              
END MODULE mo_torus_bubble_exp_nml
