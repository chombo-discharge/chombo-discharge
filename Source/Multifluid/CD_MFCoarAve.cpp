/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFCoarAve.cpp
  @brief  Implementation of CD_MFCoarAve.H
  @author Robert Marskar
*/

// Our includes
#include <CD_MFCoarAve.H>
#include <CD_NamespaceHeader.H>

MFCoarAve::MFCoarAve()
{}

MFCoarAve::~MFCoarAve()
{}

MFCoarAve::MFCoarAve(const Vector<RefCountedPtr<EBCoarAve>>& a_aveOps)
{
  this->define(a_aveOps);
}

void
MFCoarAve::define(const Vector<RefCountedPtr<EBCoarAve>>& a_aveOps)
{
  m_aveOps = a_aveOps;
}

const RefCountedPtr<EBCoarAve>&
MFCoarAve::getAveOp(const int a_phase) const
{
  CH_assert(a_phase >= 0 && a_phase < m_aveOps.size());
  return m_aveOps[a_phase];
}

#include <CD_NamespaceFooter.H>
