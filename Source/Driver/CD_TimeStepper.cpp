/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_TimeStepper.cpp
  @brief  Implementation of CD_TimeStepper.H
  @author Robert Marskar
*/

// Our includes
#include <CD_TimeStepper.H>
#include <CD_LoadBalancing.H>
#include <CD_NamespaceHeader.H>

TimeStepper::TimeStepper() {}

TimeStepper::~TimeStepper() {}

void
TimeStepper::setAmr(const RefCountedPtr<AmrMesh>& a_amr)
{
  CH_TIME("TimeStepper::setAmr(RefCountedPtr<AmrMesh>)");
  if (m_verbosity > 5) {
    pout() << "TimeStepper::setAmr(RefCountedPtr<AmrMesh>)" << endl;
  }

  CH_assert(!a_amr.isNull());

  m_amr = a_amr;
}

void
TimeStepper::setComputationalGeometry(const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry)
{
  CH_TIME("TimeStepper::setComputationalGeometry(RefCountedPtr<ComputationalGeometry>)");
  if (m_verbosity > 5) {
    pout() << "TimeStepper::setComputationalGeometry(RefCountedPtr<ComputationalGeometry>)" << endl;
  }

  CH_assert(!a_computationalGeometry.isNull());

  m_computationalGeometry = a_computationalGeometry;
}

bool
TimeStepper::needToRegrid()
{
  CH_TIME("TimeStepper::needToRegrid()");
  if (m_verbosity > 5) {
    pout() << "TimeStepper::needToRegrid()" << endl;
  }

  return false;
}

void
TimeStepper::parseRuntimeOptions()
{
  CH_TIME("TimeStepper::parseRuntimeOptions()");
  if (m_verbosity > 5) {
    pout() << "TimeStepper::parseRuntimeOptions()" << endl;
  }
}

Vector<long int>
TimeStepper::getCheckpointLoads(const std::string a_realm, const int a_level) const
{
  CH_TIME("TimeStepper::getCheckpointLoads(string, int)");
  if (m_verbosity > 5) {
    pout() << "TimeStepper::getCheckpointLoads(string, int)" << endl;
  }

  const DisjointBoxLayout& dbl      = m_amr->getGrids(a_realm)[a_level];
  const Vector<Box>&       boxArray = dbl.boxArray();

  Vector<long int> loads(boxArray.size(), 0L);
  for (int i = 0; i < boxArray.size(); i++) {
    loads[i] = boxArray[i].numPts();
  }

  return loads;
}

bool
TimeStepper::loadBalanceThisRealm(const std::string a_realm) const
{
  CH_TIME("TimeStepper::loadBalanceThisRealm(string)");
  if (m_verbosity > 5) {
    pout() << "TimeStepper::loadBalanceThisRealm(string)" << endl;
  }

  return false;
}

void
TimeStepper::loadBalanceBoxes(Vector<Vector<int>>&             a_procs,
                              Vector<Vector<Box>>&             a_boxes,
                              const std::string                a_realm,
                              const Vector<DisjointBoxLayout>& a_grids,
                              const int                        a_lmin,
                              const int                        a_finestLevel)
{
  CH_TIME(
    "TimeStepper::loadBalanceBoxes(Vector<Vector<int> >, Vector<Vector<Box>, string, Vector<DisjointBoxLayout>, int, int)");
  if (m_verbosity > 5) {
    pout()
      << "TimeStepper::loadBalanceBoxes(Vector<Vector<int> >, Vector<Vector<Box>, string, Vector<DisjointBoxLayout>, int, int)"
      << endl;
  }

  a_procs.resize(1 + a_finestLevel);
  a_boxes.resize(1 + a_finestLevel);

  for (int lvl = 0; lvl <= a_finestLevel; lvl++) {
    a_boxes[lvl] = a_grids[lvl].boxArray();

    LoadBalancing::makeBalance(a_procs[lvl], a_boxes[lvl]);
  }
}

#ifdef CH_USE_HDF5
void
TimeStepper::readCheckpointHeader(HDF5HeaderData& a_header)
{}
#endif

#ifdef CH_USE_HDF5
void
TimeStepper::writeCheckpointHeader(HDF5HeaderData& a_header) const
{}
#endif

void
TimeStepper::prePlot()
{}

#include <CD_NamespaceFooter.H>
