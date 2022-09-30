!>                                                                                   
!! @brief Read configuration parameters as Fortran namelist from an external file.
!!
!! @contact S. Rast (MPI-M)
!!
!! @par Copyright and License
!!
!! Copyright 2022 Max Planck Institute for Meteorology
!! 
!! Redistribution and use in source and binary forms, with or
!! without modification, are permitted provided that the following
!! conditions are met:
!! 
!! 1. Redistributions of source code must retain the above copyright
!!    notice, this list of conditions and the following disclaimer.
!! 
!! 2. Redistributions in binary form must reproduce the above
!!    copyright notice, this list of conditions and the following
!!    disclaimer in the documentation and/or other materials
!!    provided with the distribution.
!! 
!! 3. Neither the name of the copyright holder nor the names of its
!!    contributors may be used to endorse or promote products
!!    derived from this software without specific prior written
!!    permission.
!! 
!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
!! A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
!! HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
!! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
!! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
!! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
!! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
