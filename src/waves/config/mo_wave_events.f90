!>
!! Creation and destruction of mtime events for the wave model
!!
!! Creation and destruction of mtime events for the wave model
!!
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2023-02-09)
!!
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_wave_events

  USE mtime,                       ONLY: datetime, timedelta, newTimedelta, &
    &                                    deallocateTimedelta,               &
    &                                    event, eventGroup, newEvent, addEventToEventGroup
  USE mo_event_manager,            ONLY: addEventGroup, getEventGroup, printEventGroup
  USE mo_time_config,              ONLY: t_time_config

  IMPLICIT NONE

  PRIVATE

  ! subroutine
  PUBLIC :: create_wave_events

  ! events
  PUBLIC :: dummyWaveEvent
  PUBLIC :: waveEventGroup

  TYPE(event),      POINTER :: dummyWaveEvent  => NULL()
  TYPE(eventGroup), POINTER :: waveEventGroup  => NULL()

CONTAINS

  !>
  !! Create mtime events for wave model
  !!
  !! This routine creates mtime events for the wave model and
  !! puts them into suitable event groups.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2023-02-09)
  !!
  SUBROUTINE create_wave_events (time_config)

    TYPE(t_time_config), INTENT(IN) :: time_config  !< information for time control

    TYPE(datetime), POINTER         :: eventStartDate    => NULL(), &
      &                                eventEndDate      => NULL(), &
      &                                eventRefDate      => NULL()
    TYPE(timedelta), POINTER        :: eventInterval     => NULL()
    INTEGER                         :: waveEventsGroupID   !< might be changed to PUBLIC module variable,
                                                           !< if needed outside of this module
    INTEGER                         :: ierr
    LOGICAL                         :: lret


    ! create an event group for the wave model
    waveEventsGroupID = addEventGroup('waveEventGroup')
    waveEventGroup    => getEventGroup(waveEventsGroupID)


    ! create an example dummy event
    eventRefDate   => time_config%tc_exp_startdate
    eventStartDate => time_config%tc_exp_startdate
    eventEndDate   => time_config%tc_exp_stopdate
    eventInterval  => newTimedelta("PT01h")

    dummyWaveEvent => newEvent('dummyWave', eventRefDate, eventStartDate, &
      &                           eventEndDate, eventInterval, errno=ierr)

    lret = addEventToEventGroup(dummyWaveEvent, waveEventGroup)

    CALL printEventGroup(waveEventsGroupID)
    ! cleanup
    CALL deallocateTimedelta(eventInterval)

  END SUBROUTINE create_wave_events

END MODULE mo_wave_events

