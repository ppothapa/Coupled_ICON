!>
!! @brief Contains surface type (tile) indices used by the 
!! turbulent mixing parameterization.
!!
!! @brief Contains subroutines for initializing the AES physics package
!!
!! First version by Hui Wan, 2010-09-21
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
MODULE mo_aes_sfc_indices

  USE mo_exception, ONLY: message, message_text

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: nsfc_type, iwtr, iice, ilnd, igbm   !< index variables
  PUBLIC :: csfc                                !< sfc names
  PUBLIC :: init_sfc_indices                    !< subroutine

  INTEGER :: nsfc_type   !< total number of surface types
  INTEGER :: iwtr = 1    !< index for water-covered surface
  INTEGER :: iice = 2    !< index for ice-covered   surface
  INTEGER :: ilnd = 3    !< index for land          surface
  INTEGER :: igbm        !< index for grid-box mean

  CHARACTER(LEN=3) :: csfc(3) = (/'wtr','ice','lnd'/)

CONTAINS
  !>
  !! Set surface indices according to the simulation setup
  !! (e.g., dynamical core test, aqua-planet, or real
  !! climate simulation).
  !!
  SUBROUTINE init_sfc_indices( ctest_name )

    CHARACTER(len=*),INTENT(IN) :: ctest_name

    SELECT CASE(TRIM(ctest_name))
    CASE('APE','APE_aes','RCE','RCE_glb','RCE_Tconst','RCE_Tprescr','aes_bubble','RCEhydro','CBL_flxconst','RCEMIP_analytical')
      ! Aqua-planet simulation, no land, no ice;
      ! No needed to distinguish the aggregated grid-box mean
      ! and the value on different types of surface

      iwtr      = 1
      nsfc_type = 1
      igbm      = 0
      iice      = 999
      ilnd      = 999

    CASE('APEi','APEc','APEc_nh')
      ! Aqua-planet simulation with ice, but no land;

      iwtr      = 1
      iice      = 2
      nsfc_type = 2
      igbm      = 0
      ilnd      = 999

    CASE('JWw-Moist','LDF-Moist','jabw_m')
      ! Baroclinic wave test, no land, no ice.

      iwtr      = 1
      nsfc_type = 1
      igbm      = 0
      iice      = 999
      ilnd      = 999

      ! maybe worth trying later:
      ! iwtr      = 1
      ! iice      = 2
      ! nsfc_type = 2
      ! igbm      = 0
      ! ilnd      = 999

      ! A wild idea: a dry or completely frozen planet
      ! iice      = 1
      ! ilnd      = 2
      ! nsfc_type = 2
      ! igbm      = 0
      ! iwtr      = 999

    CASE('TPEo','TPEc')
      ! Terra-planet simulation
      ! no ocean, no ice but lakes and ice on lakes ... therefore have to use iice and iwtr
      ! No need to distinguish the aggregated grid-box mean
      ! and the value on different types of surface

      nsfc_type = 3
      iwtr      = 1
      iice      = 2
      ilnd      = 3
      igbm      = 0

    CASE DEFAULT
      ! Standard setup for real-world climate simulation.
      ! Three surface types are considered.

      iwtr      = 1
      iice      = 2
      ilnd      = 3
      nsfc_type = 3
      igbm      = 0

    END SELECT


    WRITE(message_text,*) " "
    CALL message("mo_aes_sfc_indices/init_sfc_indices",TRIM(message_text))

    WRITE(message_text,'(i3,a)')    &
      & nsfc_type, " surface type(s) activated."
    CALL message("mo_aes_sfc_indices/init_sfc_indices",TRIM(message_text))

    WRITE(message_text,'(a,4i4,a)') &
      & "Indices for water, ice, land, and grid-box mean are ", &
      & iwtr, iice, ilnd, igbm, ", respectively."
    CALL message("mo_aes_sfc_indices/init_sfc_indices",TRIM(message_text))

    WRITE(message_text,*) " "
    CALL message("mo_aes_sfc_indices/init_sfc_indices",TRIM(message_text))

  END SUBROUTINE init_sfc_indices
  !-------------

END MODULE mo_aes_sfc_indices
