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
#include <tuple>

// Chombo includes
#include <REAL.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Class which is used for run-time monitoring of events
*/
class Timer
{
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
    @brief Duration in seconds
  */
  using Duration = std::chrono::duration<Real>;

  /*!
    @brief Static function which returns the time (in seconds) between now and an arbitrary time in the past
  */
  inline static Real
  wallClock();

  /*!
    @brief Default constructor. This creates a timer without any timer events and without a name. 
  */
  inline Timer() = default;

  /*!
    @brief Default constructor. This creates a timer without any timer events and without a name. 
  */
  inline Timer(const std::string a_process);

  /*!
    @brief Destructor
  */
  inline ~Timer();

  /*!
    @brief Copy construction
    @param[in] a_other Other Timer
  */
  Timer(const Timer& a_other) = default;

  /*!
    @brief Copy assignement
    @param[in] a_other Other Timer
  */
  Timer&
  operator=(const Timer& a_other) = default;

  /*!
    @brief Start an event
    @param[in] a_event Event name
    @note This will give a run-time error if the event has already been started (i.e. startEvent(a_event) has not been called prior). 
  */
  inline void
  startEvent(const std::string a_event) noexcept;

  /*!
    @brief Stop an event
    @param[in] a_event Event name
    @note This will give a run-time error if the event has not already been started (i.e. startEvent(a_event) has not been called prior). 
  */
  inline void
  stopEvent(const std::string a_event) noexcept;

  /*!
    @brief Print all timed events to cout.
    @details This routine prints a header, the timing report for the various event, and a tail. All events are included but unfinished ones
    are not counted towards the elapsed time. 
    @param[in] a_localReportOnly If true, no reduction over mpi
  */
  inline void
  eventReport(std::ostream& a_outputStream, const bool a_localReportOnly = false) const noexcept;

  /*!
    @brief Print all timed events to a file.
    @details This routine prints a header row consisting of the event names. The remaining rows are the timings for the various ranks.
    @param[in] a_localReportOnly If true, no reduction over mpi
  */
  inline void
  writeReportToFile(const std::string a_fileName) const noexcept;

  /*!
    @brief Clear all events
  */
  inline void
  clear() noexcept;

protected:
  /*!
    @brief Enum for indexing tuple. A rare case of where unscoped enums are useful. 
  */
  enum EventFields
  {
    StoppedEvent = 0,
    StartClock   = 1,
    ElapsedTime  = 2
  };

  /*!
    @brief Process name. Used for input/output.
  */
  std::string m_processName;

  /*!
    @brief Timer events. First entry is the name of the event. 
    @details The second is the identifier for the event
  */
  std::map<std::string, std::tuple<bool, TimePoint, Duration>> m_events;

  /*!
    @brief Print report header
  */
  void
  printReportHeader(std::ostream& a_outputStream) const noexcept;

  /*!
    @brief Print report tail
  */
  void
  printReportTail(std::ostream& a_outputStream, const std::pair<Real, Real> a_elapsedTime) const noexcept;

  /*!
    @brief Compute the total time for all finished events. 
    @return First entry is the local time. Second entry is the global time. 
  */
  std::pair<Real, Real>
  computeTotalElapsedTime(const bool a_localReportOnly) const noexcept;
};

#include <CD_NamespaceFooter.H>

#include <CD_TimerImplem.H>

#endif
