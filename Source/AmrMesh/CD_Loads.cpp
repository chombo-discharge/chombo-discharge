/* chombo-discharge
 * Copyright © 2024 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_Loads.cpp
  @brief  Implementation of CD_Loads.H
  @author Robert Marskar
*/

// Std includes
#include <algorithm>

// Chombo includes
#include <SPMD.H>
#include <CH_Timer.H>

// Our includes
#include <CD_Loads.H>
#include <CD_NamespaceHeader.H>

Loads::Loads() noexcept
{
  CH_TIME("Loads::Loads");

  this->resetLoads();
}

Loads::~Loads() noexcept
{
  CH_TIME("Loads::~Loads");

  m_loads.clear();
}

std::map<int, Real>&
Loads::getLoads() noexcept
{
  CH_TIME("Loads::getLoads");

  return (m_loads);
}

const std::map<int, Real>&
Loads::getLoads() const noexcept
{
  CH_TIME("Loads::getLoads");

  return (m_loads);
}

void
Loads::resetLoads() noexcept
{
  CH_TIME("Loads::resetLoads");

  m_loads.clear();

  for (int irank = 0; irank < numProc(); irank++) {
    m_loads[irank] = 0;
  }
}

void
Loads::assignLoads(const std::map<int, Real>& a_assignedLoads) noexcept
{
  CH_TIME("Loads::assignLoads(std::map)");

  const int numRanks = numProc();

  if (a_assignedLoads.size() != numRanks) {
    MayDay::Abort("Loads::assignLoads(std::map) -- a_assignedLoads.size() != numProc()");
  }
  else {
    for (int irank = 0; irank < numRanks; irank++) {
      m_loads[irank] = a_assignedLoads.at(irank);
    }
  }
}

void
Loads::assignLoads(const std::vector<Real>& a_assignedLoads) noexcept
{
  CH_TIME("Loads::assignLoads(std::vector)");

  const int numRanks = numProc();

  if (a_assignedLoads.size() != numRanks) {
    MayDay::Abort("Loads::assignLoads(std::vector) -- a_assignedLoads.size() != numProc()");
  }
  else {
    for (int irank = 0; irank < numRanks; irank++) {
      m_loads[irank] = a_assignedLoads[irank];
    }
  }
}

void
Loads::assignLoads(const Vector<Real>& a_assignedLoads) noexcept
{
  CH_TIME("Loads::assignLoads(Vector)");

  const int numRanks = numProc();

  if (a_assignedLoads.size() != numRanks) {
    MayDay::Abort("Loads::assignLoads(Vector) -- a_assignedLoads.size() != numProc()");
  }
  else {
    for (int irank = 0; irank < numRanks; irank++) {
      m_loads[irank] = a_assignedLoads[irank];
    }
  }
}

void
Loads::incrementLoads(const std::map<int, Real>& a_increments) noexcept
{
  CH_TIME("Loads::incrementLoads(std::map)");

  const int numRanks = numProc();

  if (a_increments.size() != numRanks) {
    MayDay::Abort("Loads::incrementLoads(std::map) -- a_increments.size() != numProc()");
  }
  else {
    for (int irank = 0; irank < numRanks; irank++) {
      m_loads[irank] += a_increments.at(irank);
    }
  }
}

void
Loads::incrementLoads(const std::vector<Real>& a_increments) noexcept
{
  CH_TIME("Loads::incrementLoads(std::vector)");

  const int numRanks = numProc();

  if (a_increments.size() != numRanks) {
    MayDay::Abort("Loads::incrementLoads(std::vector) -- a_increments.size() != numProc()");
  }
  else {
    for (int irank = 0; irank < numRanks; irank++) {
      m_loads[irank] += a_increments[irank];
    }
  }
}

void
Loads::incrementLoads(const Vector<Real>& a_increments) noexcept
{
  CH_TIME("Loads::incrementLoads(Vector)");

  const int numRanks = numProc();

  if (a_increments.size() != numRanks) {
    MayDay::Abort("Loads::incrementLoads(Vector) -- a_increments.size() != numProc()");
  }
  else {
    for (int irank = 0; irank < numRanks; irank++) {
      m_loads[irank] += a_increments[irank];
    }
  }
}

void
Loads::incrementLoad(const int a_rank, const Real a_increment) noexcept
{
  CH_TIME("Loads::incrementLoad");

  const int numRanks = numProc();

  if (a_rank >= numRanks) {
    MayDay::Abort("Loads::incrementLoad -- 'a_rank > numProc()'");
  }
  else {
    m_loads[a_rank] += a_increment;
  }
}

std::vector<std::pair<int, Real>>
Loads::getSortedLoads() const noexcept
{
  CH_TIME("Loads::sortedLoads");

  std::vector<std::pair<int, Real>> sortedLoads;

  // Insert loads into a vector so we can sort
  for (const auto& curLoad : m_loads) {
    sortedLoads.emplace_back(curLoad.first, curLoad.second);
  }

  auto loadSort = [](const std::pair<int, Real>& A, const std::pair<int, Real>& B) -> bool {
    return A.second < B.second;
  };

  std::sort(sortedLoads.begin(), sortedLoads.end(), loadSort);

#if 1 // Remove later

  for (int i = 0; i < sortedLoads.size() - 1; i++) {
    if (sortedLoads[i].second > sortedLoads[i + 1].second) {
      MayDay::Abort("Loads::getSortedLoads - crap on sorting method");
    }
  }
#endif

  return sortedLoads;
}

#include <CD_NamespaceFooter.H>
