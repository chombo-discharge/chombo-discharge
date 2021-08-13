/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_TimerImplem.H
  @brief  Implementation of CD_Timer.H
  @author Robert Marskar
*/

#ifndef CD_TimerImplem_H
#define CD_TimerImplem_H

#include <CD_Timer.H>
#include <CD_NamespaceHeader.H>

inline
Real Timer::wallClock(){

  // TLDR: This just returns the current time in seconds (since an arbitrary time in the past). 
  const auto currentTime       = std::chrono::steady_clock::now();
  const auto durationInSeconds = std::chrono::duration<Real>(currentTime.time_since_epoch());

  return durationInSeconds.count();
}

inline
Timer::~Timer() {
  m_events.clear();
}

inline
void Timer::startEvent(const std::string a_event) {

  // Only start an event if it has not already been started. 
  if (m_events.find(a_event) == m_events.end() ) {
    const TimePoint start = Clock::now();
    
    m_events.emplace(a_event, std::make_tuple(false, start, start));
  }
}

inline
void Timer::stopEvent(const std::string a_event) {
  if (m_events.find(a_event) == m_events.end()) {
    tuple<bool, TimePoint, TimePoint>& event = m_events.at(a_event);

    const TimePoint start = std::get<StartClock>(event);
    const TimePoint stop  = Clock::now();

    event = std::make_tuple(true, start, stop);
  }
  else{
    std::cerr << "Timer::stopEvent -- event '" + a_event + "' has not been started\n";
  }
}
 
#include <CD_NamespaceFooter.H>

#endif
