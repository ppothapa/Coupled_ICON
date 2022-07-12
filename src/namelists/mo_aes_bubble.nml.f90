!>                                                                                   
!! @brief Read configuration parameters as Fortran namelist from an external file.
!!
!! @contact S. Rast (MPI-M)
!!
!! @par Copyright and License
!!
!! Copyright 2022 Max Planck Institute for Meteorology
!! This program is free software: you can restribute it and/or modify it under the terms of the
!! BSD-3-clause license
!!                                                                                   
MODULE mo_aes_bubble_nml

  USE mo_aes_bubble_config, ONLY: aes_bubble_config, init_aes_bubble_config, &
       &                          eval_aes_bubble_config, print_aes_bubble_config
  USE mo_process_nml,       ONLY: process_nml
  
  IMPLICIT NONE                                                                      
  PRIVATE                                                                            
  PUBLIC :: process_aes_bubble_nml                                                    
                                                                                     
  NAMELIST /aes_bubble_nml/ aes_bubble_config                             
                                                                                     
CONTAINS                                                                             
                                                                                     
  SUBROUTINE process_aes_bubble_nml(filename)                                         
    !                                                                                
    CHARACTER(LEN=*), INTENT(in) :: filename                                         
    !                                                                                
    CALL init_aes_bubble_config                                                       
    !                                                                                
    CALL process_nml(filename, 'aes_bubble_nml', nml_read, nml_write)                 
    !                                                                                
    CALL eval_aes_bubble_config
    CALL print_aes_bubble_config
  CONTAINS                                                                           
    !                                                                                
    SUBROUTINE nml_read(funit)                                                       
      INTEGER, INTENT(in) :: funit                                                   
      READ(funit, NML=aes_bubble_nml)                                                 
    END SUBROUTINE nml_read                                                          
    !                                                                                
    SUBROUTINE nml_write(funit)                                                      
      INTEGER, INTENT(in) :: funit                                                   
      WRITE(funit, NML=aes_bubble_nml)                                                
    END SUBROUTINE nml_write                                                         
    !                                                                                
  END SUBROUTINE process_aes_bubble_nml
                                              
END MODULE mo_aes_bubble_nml
