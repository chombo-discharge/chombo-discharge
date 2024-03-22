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
#include <iostream>
#include <fstream>
#include <iomanip>

// Chombo includes
#ifdef CH_MPI
#include <SPMD.H>
#endif

// Our includes
#include <CD_Timer.H>
#include <CD_NamespaceHeader.H>

inline Real
Timer::wallClock()
{

  // TLDR: This just returns the current time in seconds (since an arbitrary time in the past).
  const auto currentTime       = std::chrono::steady_clock::now();
  const auto durationInSeconds = std::chrono::duration<Real>(currentTime.time_since_epoch());

  return durationInSeconds.count();
}

inline Timer::Timer(const std::string a_processName)
{
  m_processName = a_processName;
  m_events.clear();
}

inline Timer::~Timer()
{
  m_events.clear();
}

inline void
Timer::clear() noexcept
{
  m_events.clear();
}

inline void
Timer::startEvent(const std::string a_event) noexcept
{
  if (m_events.find(a_event) == m_events.end()) { // Create a new event.
    const TimePoint startTime   = Clock::now();
    const Duration  elapsedTime = Duration(0.0);

    m_events.emplace(a_event, std::make_tuple(false, startTime, elapsedTime));
  }
  else { // If the event is not new, we just leave the elapsed time intact.
    //    const TimePoint startTime = Clock::now();
    std::tuple<bool, TimePoint, Duration>& event        = m_events.at(a_event);
    const bool                             stoppedEvent = std::get<StoppedEvent>(event);

    if (stoppedEvent) { // Only start the event if it is not already running.
      const TimePoint startTime           = Clock::now();
      const Duration& previousElapsedTime = std::get<ElapsedTime>(event);

      event = std::make_tuple(false, startTime, previousElapsedTime);
    }
  }
}

inline void
Timer::stopEvent(const std::string a_event) noexcept
{
  if (m_events.find(a_event) != m_events.end()) {
    std::tuple<bool, TimePoint, Duration>& event = m_events.at(a_event);

    const TimePoint startTime           = std::get<StartClock>(event);
    const TimePoint stopTime            = Clock::now();
    const Duration  curElapsedTime      = stopTime - startTime;
    const Duration  previousElapsedTime = std::get<ElapsedTime>(event);
    const Duration  totalElapsedTime    = previousElapsedTime + curElapsedTime;

    event = std::make_tuple(true, startTime, totalElapsedTime);
  }
  else {
    std::cerr << "Timer::stopEvent -- event '" + a_event + "' has not been started\n";
  }
}

inline void
Timer::printReportHeader(std::ostream& a_outputStream) const noexcept
{

  std::stringstream header;
  header
    << "| ---------------------------------------------------------------------------------------------------------|"
    << "\n"
    << "| " + m_processName + " kernel report: "
    << "\n"
    << "| ---------------------------------------------------------------------------------------------------------|"
    << "\n"
    << "| " << std::left << std::setw(25) << "Event"
    << "| " << std::right << std::setw(8) << "Loc. (s)"
    << "| " << std::right << std::setw(8) << "Loc. (%)"
#ifdef CH_MPI
    << "| " << std::right << std::setw(8) << "Min. (s)"
    << "| " << std::right << std::setw(8) << "Max. (s)"
    << "| " << std::right << std::setw(8) << "Avg. (s)"
    << "| " << std::right << std::setw(8) << "Dev. (s)"
    << "| " << std::right << std::setw(8) << "Min rank"
    << "| " << std::right << std::setw(8) << "Max rank"
    << "| "
    << "\n"
#endif
    << "| ---------------------------------------------------------------------------------------------------------|"
    << "\n";

  a_outputStream << header.str();
}

inline void
Timer::printReportTail(std::ostream& a_outputStream, const std::pair<Real, Real> a_elapsedTime) const noexcept
{
  std::stringstream tail;

  tail << "| ---------------------------------------------------------------------------------------------------------|"
       << "\n"
       << "| Local elapsed time  = " << std::fixed << std::setprecision(4) << a_elapsedTime.first << " seconds"
       << "\n"
       << "| Global elapsed time = " << std::fixed << std::setprecision(4) << a_elapsedTime.second << " seconds"
       << "\n"
       << "| ---------------------------------------------------------------------------------------------------------|"
       << "\n";

  a_outputStream << tail.str();
}

inline void
Timer::eventReport(std::ostream& a_outputStream, const bool a_localReportOnly) const noexcept
{
  // This routine prints a header, the timing report for the various event, and a tail. All events are included
  // but unfinished ones are not counted towards the elapsed time.

  this->printReportHeader(a_outputStream);

  const std::pair<Real, Real> totalTime      = this->computeTotalElapsedTime(a_localReportOnly);
  const Real                  totalTimeLocal = totalTime.first;

  for (const auto& e : m_events) {

    const std::string& eventName     = e.first;
    const bool         finishedEvent = std::get<StoppedEvent>(e.second);

    std::stringstream ssLocalDuration;
    std::stringstream ssPercentage;
    std::stringstream ssMinDuration;
    std::stringstream ssMaxDuration;
    std::stringstream ssAvgDuration;
    std::stringstream ssDevDuration;
    std::stringstream ssMinRank;
    std::stringstream ssMaxRank;

    if (finishedEvent) {
      const Duration localDuration          = std::get<ElapsedTime>(e.second);
      const Real     localDurationInSeconds = localDuration.count();

      ssLocalDuration << std::fixed << std::setprecision(4) << localDurationInSeconds;
      ssPercentage << std::fixed << std::setprecision(4) << 100. * localDurationInSeconds / totalTimeLocal;

#ifdef CH_MPI
      if (!a_localReportOnly) {

        const int srcRank = 0;
        const int numRank = numProc();

        Vector<Real> allDurations(numRank);
        gather(allDurations, localDurationInSeconds, srcRank);

        // Compute average time and standard deviation
        Real average;
        Real sigma;

        Real maxTime = std::numeric_limits<Real>::min();
        Real minTime = std::numeric_limits<Real>::max();

        int minRank;
        int maxRank;

        if (procID() == 0) {

          average = 0.0;
          for (int i = 0; i < allDurations.size(); i++) {
            average += allDurations[i];

            if (allDurations[i] < minTime) {
              minRank = i;
              minTime = allDurations[i];
            }
            if (allDurations[i] > maxTime) {
              maxRank = i;
              maxTime = allDurations[i];
            }
          }
          average *= 1. / numRank;

          sigma = 0.0;
          for (int i = 0; i < allDurations.size(); i++) {
            sigma += sqrt(std::pow(allDurations[i] - average, 2));
          }
          sigma = sqrt(sigma / numRank);
        }

        broadcast(average, srcRank);
        broadcast(sigma, srcRank);
        broadcast(maxTime, srcRank);
        broadcast(minTime, srcRank);
        broadcast(minRank, srcRank);
        broadcast(maxRank, srcRank);

        ssMinDuration << std::fixed << std::setprecision(4) << minTime;
        ssMaxDuration << std::fixed << std::setprecision(4) << maxTime;
        ssAvgDuration << std::fixed << std::setprecision(4) << average;
        ssDevDuration << std::fixed << std::setprecision(4) << sigma;
        ssMinRank << std::fixed << std::setprecision(4) << minRank;
        ssMaxRank << std::fixed << std::setprecision(4) << maxRank;
      }
      else {
        ssMinDuration << " N/A ";
        ssMaxDuration << " N/A ";
        ssAvgDuration << " N/A ";
        ssDevDuration << " N/A ";
        ssMinRank << " N/A ";
        ssMaxRank << " N/A ";
      }
#endif
    }
    else {
      ssLocalDuration << " - ";
      ssPercentage << " - ";
      ssMinDuration << " - ";
      ssMaxDuration << " - ";
      ssAvgDuration << " - ";
      ssDevDuration << " - ";
      ssMinRank << " - ";
      ssMaxRank << " - ";
    }

    std::stringstream outputString;

    // Now print to pout
    outputString << "| " << std::left << std::setw(25) << eventName << "| " << std::right << std::setw(8)
                 << ssLocalDuration.str() << "| " << std::right << std::setw(8) << ssPercentage.str()
#ifdef CH_MPI
                 << "| " << std::right << std::setw(8) << ssMinDuration.str() << "| " << std::right << std::setw(8)
                 << ssMaxDuration.str() << "| " << std::right << std::setw(8) << ssAvgDuration.str() << "| "
                 << std::right << std::setw(8) << ssDevDuration.str() << "| " << std::right << std::setw(8)
                 << ssMinRank.str() << "| " << std::right << std::setw(8) << ssMaxRank.str()
#endif
                 << "| "
                 << "\n";

    a_outputStream << outputString.str();
  }

  this->printReportTail(a_outputStream, totalTime);
}

inline std::pair<Real, Real>
Timer::computeTotalElapsedTime(const bool a_localReportOnly) const noexcept
{

  Real elapsedTimeLocal  = 0.0;
  Real elapsedTimeGlobal = 0.0;

  for (const auto& e : m_events) {
    const bool finishedEvent = std::get<StoppedEvent>(e.second);

    if (finishedEvent) {
      const Duration elapsedTime            = std::get<ElapsedTime>(e.second);
      const Real     localDurationInSeconds = elapsedTime.count();

      elapsedTimeLocal += localDurationInSeconds;
    }
  }

#ifdef CH_MPI
  if (!a_localReportOnly) {
    MPI_Allreduce(&elapsedTimeLocal, &elapsedTimeGlobal, 1, MPI_CH_REAL, MPI_MAX, Chombo_MPI::comm);
  }
  else {
    elapsedTimeGlobal = elapsedTimeLocal;
  }
#else
  elapsedTimeGlobal  = elapsedTimeLocal;
#endif

  return std::make_pair(elapsedTimeLocal, elapsedTimeGlobal);
}

inline void
Timer::writeReportToFile(const std::string a_fileName) const noexcept
{

  // First, go through all events and gather the times on the master rank.
  std::vector<std::string> headerRow;
  headerRow.emplace_back("# MPI rank");

  // These are all the times for the finished events. We organize this is a 2D vector which is M x numRank long. Here,
  // finishedEvents[0] contains the elapsed times for one of the events. The length of finishedEvents[0] is exactly
  // numRank long
  std::vector<std::vector<Real>> finishedEvents;

#if CH_MPI
  constexpr int srcRank  = 0;
  const int     numRanks = numProc();
#else
  const int numRanks = 1;
#endif

  for (const auto& e : m_events) {
    const std::string&                           eventName = e.first;
    const std::tuple<bool, TimePoint, Duration>& eventData = e.second;

    const bool finishedEvent = std::get<StoppedEvent>(eventData);

    // If we had an event that finished, we gather all the data on the master rank.
    if (finishedEvent) {

      const Real   localEventDuration = std::get<ElapsedTime>(eventData).count();
      Vector<Real> allTimes(numRanks, 0.0);
#ifdef CH_MPI

      gather(allTimes, localEventDuration, srcRank);
#else
      allTimes[0] = localEventDuration;
#endif

      finishedEvents.emplace_back(allTimes.stdVector());
      headerRow.emplace_back(eventName);
    }
  }

  // Open a file and start writing.
#if CH_MPI
  if (procID() == srcRank) {
#endif
    std::ofstream f;
    f.open(a_fileName, std::ios_base::trunc);

    const int width = 12;
    for (int col = 0; col < headerRow.size(); col++) {
      f << std::left << std::setw(width) << headerRow[col] << "\t";
    }
    f << "\n";

    // Write the data for rank = 0 etc
    for (int irank = 0; irank < numRanks; irank++) {
      f << std::left << std::setw(width) << irank << "\t";
      for (int event = 0; event < finishedEvents.size(); event++) {
        f << std::left << std::setw(width) << finishedEvents[event][irank] << "\t";
      }
      f << "\n";
    }

    f.close();
#if CH_MPI
  }
#endif
}

#include <CD_NamespaceFooter.H>

#endif
