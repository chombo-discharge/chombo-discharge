/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_Timer.H
  @brief  Implementation of CD_Timer.H
  @author Robert Marskar
  @todo   Make sure that events that were started but not stopped are reports as unfinished (or somesuch)
*/

#ifndef CD_Timer_H
#define CD_Timer_H

// Std includes
#include <chrono>
#include <map>

// Chombo includes
#include <REAL.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Class which is used for run-time monitoring of events
*/
class Timer {
public:

  /*!
    @brief Clock implements
  */
  using Clock = std::chrono::steady_clock;

  /*!
    @brief Point in time
  */
  using TimePoint = std::chrono::steady_clock::time_point;

  /*!
    @brief Static function which returns the time (in seconds) between now and an arbitrary time in the past
  */
  inline
  static Real wallClock();

  /*!
    @brief Default constructor. This creates a timer without any timer events and without a name. 
  */
  inline
  Timer();

  /*!
    @brief Destructor
  */
  inline
  ~Timer();

  /*!
    @brief Disallowed copy construction
    @param[in] a_other Other Timer
  */
  Timer(const Timer& a_other) = delete;

  /*!
    @brief Disallowed move construction
    @param[in] a_other Other Timer
  */
  Timer(const Timer&& a_other) = delete;

  /*!
    @brief Disallowed copy assignement
    @param[in] a_other Other Timer
  */
  Timer& operator=(const Timer& a_other) = delete;

  /*!
    @brief Disallowed move assignement
    @param[in] a_other Other Timer
  */
  Timer& operator=(const Timer&& a_other) = delete;

  /*!
    @brief Start an event
    @param[in] a_event Event name
    @note This will give a run-time error if the event has already been started (i.e. startEvent(a_event) has not been called prior). 
  */
  inline
  void startEvent(const std::string a_event);

  /*!
    @brief Stop an event
    @param[in] a_event Event name
    @note This will give a run-time error if the event has not already been started (i.e. startEvent(a_event) has not been called prior). 
  */
  inline
  void stopEvent(const std::string a_event);

  /*!
    @brief Reset all events. 
    @note This clears out m_events. 
  */
  inline
  void resetEvents();

protected:

  /*!
    @brief Enum for indexing tuple. A rare case of where unscoped enums are useful. 
  */
  enum EventFields {
    StoppedEvent = 0,
    StartClock   = 1,
    StopClock    = 2
  };

  /*!
    @brief Timer events. First entry is the name of the event, the second is the time.
  */
  std::map<std::string, std::tuple<bool, TimePoint, TimePoint> > m_events;
};

#include <CD_NamespaceFooter.H>

#include <CD_TimerImplem.H>

#endif
