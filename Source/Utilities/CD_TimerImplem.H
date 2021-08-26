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

// Std includes
#include <sstream>

// Our includes
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
Timer::Timer(const std::string a_processName) {
  m_processName = a_processName;
  m_events.clear();
}

inline
Timer::~Timer() {
  m_events.clear();
}

inline
void Timer::startEvent(const std::string a_event) noexcept {

  // Only start an event if it has not already been started. 
  if (m_events.find(a_event) == m_events.end() ) {
    const TimePoint start = Clock::now();
    
    m_events.emplace(a_event, std::make_tuple(false, start, start));
  }
}

inline
void Timer::stopEvent(const std::string a_event) noexcept {
  if (m_events.find(a_event) != m_events.end()) {
    tuple<bool, TimePoint, TimePoint>& event = m_events.at(a_event);

    const TimePoint start = std::get<StartClock>(event);
    const TimePoint stop  = Clock::now();

    event = std::make_tuple(true, start, stop);
  }
  else{
    std::cerr << "Timer::stopEvent -- event '" + a_event + "' has not been started\n";
  }
}

inline
void Timer::printReportHeader() const noexcept {

  std::stringstream header;
  header << "| ---------------------------------------------------------------------------------------------------------|" << "\n"
	 << "| " + m_processName + " kernel report: "                                                                      << "\n"
         << "| ---------------------------------------------------------------------------------------------------------|" << "\n"
	 << "| " << std::left << std::setw(25) << "Event"
	 << "| " << std::right << std::setw(8) << "Loc."
	 << "| " << std::right << std::setw(8) << "%"    
#ifdef CH_MPI
	 << "| " << std::right << std::setw(8) << "Min."
	 << "| " << std::right << std::setw(8) << "Max."
	 << "| " << std::right << std::setw(8) << "Avg." 
	 << "| " << std::right << std::setw(8) << "Dev."
	 << "| " << std::right << std::setw(8) << "Min rank"
	 << "| " << std::right << std::setw(8) << "Max rank"
	 << "| " << "\n"
#endif
         << "| ---------------------------------------------------------------------------------------------------------|" << "\n";

  pout() << header.str();
}

inline
void Timer::printReportTail(const Real a_elapsedTime) const noexcept {
  std::stringstream tail;
  
  tail << "| ---------------------------------------------------------------------------------------------------------|" << "\n"
       << "| Total elapsed time = " << std::fixed << std::setprecision(3) << a_elapsedTime << " seconds"                 << "\n"
       << "| ---------------------------------------------------------------------------------------------------------|" << "\n";

  pout() << tail.str();
}

inline
void Timer::eventReport() const noexcept{
  // This routine prints a header, the timing report for the various event, and a tail. All events are included
  // but unfinished ones are not counted towards the elapsed time. 
  
  this->printReportHeader();

  const Real totalTime = this->computeTotalElapsedTime();
    
  for (const auto& e : m_events){

    const std::string& eventName     = e.first;
    const bool         finishedEvent = std::get<StoppedEvent> (e.second);

    const int srcRank = 0;
    const int numRank = numProc();

    std::stringstream ssLocalDuration; 
    std::stringstream ssPercentage;    
    std::stringstream ssMinDuration;
    std::stringstream ssMaxDuration;
    std::stringstream ssAvgDuration;
    std::stringstream ssDevDuration;
    std::stringstream ssMinRank;
    std::stringstream ssMaxRank;    

    if(finishedEvent){
      const std::chrono::duration<Real> localDuration = std::get<StopClock>(e.second) - std::get<StartClock>(e.second);
      const Real localDurationInSeconds = localDuration.count();

      ssLocalDuration << std::fixed << std::setprecision(3) << localDurationInSeconds;
      ssPercentage    << std::fixed << std::setprecision(3) << 100.*localDurationInSeconds/totalTime;

#ifdef CH_MPI
      Vector<Real> allDurations(numRank);
      gather(allDurations, localDurationInSeconds, srcRank);

      // Compute average time and standard deviation
      Real average;
      Real sigma;
      
      Real maxTime = std::numeric_limits<Real>::min();
      Real minTime = std::numeric_limits<Real>::max();

      int minRank;
      int maxRank;
      
      if(procID() == 0){
	
	average = 0.0;
	for (int i = 0; i < allDurations.size(); i++){
	  average += allDurations[i];

	  if(allDurations[i] < minTime){
	    minRank = i;
	    minTime = allDurations[i];
	  }
	  if(allDurations[i] > maxTime){
	    maxRank = i;
	    maxTime = allDurations[i];
	  }	  	  
	}
	average *= 1./numRank;

	sigma = 0.0;
	for (int i = 0; i < allDurations.size(); i++){
	  sigma += std::pow(allDurations[i]-average, 2);
	}
	sigma *= 1./numRank;
      }

      broadcast(average, srcRank);
      broadcast(sigma,   srcRank);
      broadcast(maxTime, srcRank);
      broadcast(minTime, srcRank);
      broadcast(minRank, srcRank);
      broadcast(maxRank, srcRank);      

      ssMinDuration << std::fixed << std::setprecision(3) << minTime;
      ssMaxDuration << std::fixed << std::setprecision(3) << maxTime;            
      ssAvgDuration << std::fixed << std::setprecision(3) << average;
      ssDevDuration << std::fixed << std::setprecision(3) << sigma;
      ssMinRank     << std::fixed << std::setprecision(3) << minRank;
      ssMaxRank     << std::fixed << std::setprecision(3) << maxRank;
#endif
    }
    else{
      ssLocalDuration  << " - ";
      ssPercentage     << " - ";
      ssMinDuration    << " - ";
      ssMaxDuration    << " - ";
      ssAvgDuration    << " - ";
      ssDevDuration    << " - ";
      ssMinRank        << " - ";
      ssMaxRank        << " - ";    
    }

    std::stringstream outputString;

    // Now print to pout
    outputString << "| " << std::left << std::setw(25) << eventName
		 << "| " << std::right << std::setw(8) << ssLocalDuration.str()
		 << "| " << std::right << std::setw(8) << ssPercentage.str()
#ifdef CH_MPI
		 << "| " << std::right << std::setw(8) << ssMinDuration.str()
		 << "| " << std::right << std::setw(8) << ssMaxDuration.str()
		 << "| " << std::right << std::setw(8) << ssAvgDuration.str()
		 << "| " << std::right << std::setw(8) << ssDevDuration.str()
		 << "| " << std::right << std::setw(8) << ssMinRank.str()
		 << "| " << std::right << std::setw(8) << ssMaxRank.str()
#endif
		 << "| "
		 << "\n";

    pout() << outputString.str();
  }

  this->printReportTail(totalTime);  
}

inline
Real Timer::computeTotalElapsedTime() const noexcept {

  Real totalTime = 0.0;

  for (const auto& e : m_events){
    const bool finishedEvent = std::get<StoppedEvent>(e.second);

    if(finishedEvent){
      const std::chrono::duration<Real> localDuration = std::get<StopClock>(e.second) - std::get<StartClock>(e.second);
      const Real localDurationInSeconds = localDuration.count();

      totalTime += localDurationInSeconds;
    }
  }

  return totalTime;
}
 
#include <CD_NamespaceFooter.H>

#endif
