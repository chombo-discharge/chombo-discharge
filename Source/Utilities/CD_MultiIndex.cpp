/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MultiIndex.cpp
  @brief  Implementation of CD_MultiIndex.H
  @author Robert Marskar
*/

// Our includes
#include <CD_MultiIndex.H>
#include <CD_NamespaceHeader.H>

MultiIndex::MultiIndex(const int a_order)
{
  this->define(a_order);
}

void
MultiIndex::define(const int a_order)
{
  m_order = a_order;

  // Make indices and define maps.
  makeIndices();
  makeMaps();
  reset();
}

void
MultiIndex::reset()
{
  m_iter = m_indices.begin();
}

IntVect
MultiIndex::getCurrentIndex() const
{
  return *m_iter;
}

int
MultiIndex::getOrder() const
{
  return m_order;
}

int
MultiIndex::getNumIndices() const
{
  return m_indices.size();
}

int
MultiIndex::getLinearIndex(const IntVect a_multiIndex) const
{
  int ret = -1;

  if (m_mapToLinearIndex.find(a_multiIndex) != m_mapToLinearIndex.end()) {
    ret = m_mapToLinearIndex.at(a_multiIndex);
  }
  else {
    MayDay::Abort("MultiIndex::getLinearIndex - index out of range!");
  }

  return ret;
}

IntVect
MultiIndex::getMultiIndex(const int a_linearIndex) const
{
  IntVect ret;

  if (m_mapToMultiIndex.find(a_linearIndex) != m_mapToMultiIndex.end()) {
    ret = m_mapToMultiIndex.at(a_linearIndex);
  }
  else {
    MayDay::Abort("MultiIndex::getMultiIndex - index out of range!");
  }

  return ret;
}

int
MultiIndex::operator[](const int a_dir) const
{
  return (*m_iter)[a_dir];
}

bool
MultiIndex::ok() const
{
  return (m_iter != m_indices.end());
}

int
MultiIndex::factorial(const int a_n) const
{
  if (a_n > 1) {
    return a_n * factorial(a_n - 1);
  }
  else {
    return 1;
  };
}

int
MultiIndex::factorial() const
{
  int fact = 1;
  for (int dir = 0; dir < SpaceDim; dir++) {
    fact *= factorial((*m_iter)[dir]);
  }

  return fact;
}

int
MultiIndex::norm() const
{
  return norm((*m_iter));
}

int
MultiIndex::norm(const IntVect a_iv) const
{
  int retval = 0;
  for (int dir = 0; dir < SpaceDim; dir++) {
    retval += a_iv[dir];
  }

  return retval;
}

Real
MultiIndex::pow(const RealVect& a_vec)
{

  Real retval = 1.;

  for (int dir = 0; dir < SpaceDim; dir++) {
    retval *= std::pow(a_vec[dir], (*m_iter)[dir]);
  }

  return retval;
}

void
MultiIndex::operator++()
{
  m_iter++;
}

void
MultiIndex::makeIndices()
{

  IntVect cur = IntVect::Zero;

  while (MultiIndex::norm(cur) <= m_order) {
    m_indices.emplace_back(cur);

    // Now go to the next index in lexigraphical order
    if (norm(cur) < m_order) { // Can raise first index.
      cur[0]++;
    }
    else if (norm(cur) == m_order) { // Can't raise first order, check next index.

      IntVect next = IntVect(D_DECL(0, cur[1] + 1, cur[2]));
      if (norm(next) <= m_order) { // Ok, this is a valid index.
        cur = next;
      }
#if CH_SPACEDIM == 3
      else if (norm(next) > m_order) {
        IntVect next = IntVect(D_DECL(0, 0, cur[2] + 1));
        if (norm(next) <= m_order) {
          cur = next;
        }
        else {
          cur[0]++; // Designed to break out of loop
        }
      }
#endif
      else {
        cur[0]++; // Designed to break out of loop
      }
    }
  }
}

void
MultiIndex::makeMaps()
{

  // Create the map from multi-index to linear index (i.e. column entry).
  for (int k = 0; k < m_indices.size(); k++) {
    m_mapToLinearIndex.emplace(m_indices[k], k);
  }

  for (const auto& lm : m_mapToLinearIndex) {
    m_mapToMultiIndex.emplace(lm.second, lm.first);
  }
}

#include <CD_NamespaceFooter.H>
