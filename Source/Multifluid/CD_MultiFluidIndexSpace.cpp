/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_MultiFluidIndexSpace.cpp
 * @brief  Implementation of CD_MultiFluidIndexSpace.H
 * @author Robert Marskar
 */

// Chombo includes
#include <AllRegularService.H>

// Our includes
#include <CD_MultiFluidIndexSpace.H>
#include <CD_MemoryReport.H>
#include <CD_NamespaceHeader.H>

MultiFluidIndexSpace::MultiFluidIndexSpace()
{
  m_ebis.resize(phase::numPhases);
  for (int i = 0; i < m_ebis.size(); i++) {
    m_ebis[i] = RefCountedPtr<EBIndexSpace>(new EBIndexSpace());
  }
}

MultiFluidIndexSpace::~MultiFluidIndexSpace() = default;

void
MultiFluidIndexSpace::define(const Box&                      a_domain,
                             const RealVect&                 a_origin,
                             const Real&                     a_dx,
                             const Vector<GeometryService*>& a_geoservers,
                             const bool                      a_distributedData,
                             int                             a_nCellMax,
                             int                             a_max_coar,
                             bool /*a_fix_phase*/)
{

  // Define the gas geoserver
  if (a_distributedData) {
    m_ebis[phase::gas]->setDistributedData();
  }
  m_ebis[phase::gas]->define(a_domain, a_origin, a_dx, *a_geoservers[phase::gas], a_nCellMax, a_max_coar);

  MemoryReport::getMaxMinMemoryUsage();

  // Define the solid state geoserver. This EBIS might not exist.
  if (a_geoservers[phase::solid] == NULL) {
    m_ebis[phase::solid] = RefCountedPtr<EBIndexSpace>(nullptr);
  }
  else {
    if (a_distributedData) {
      m_ebis[phase::solid]->setDistributedData();
    }
    m_ebis[phase::solid]->define(a_domain, a_origin, a_dx, *a_geoservers[phase::solid], a_nCellMax, a_max_coar);

    MemoryReport::getMaxMinMemoryUsage();
  }
}

const RefCountedPtr<EBIndexSpace>&
MultiFluidIndexSpace::getEBIndexSpace(const phase::which_phase a_phase) const
{
  return m_ebis[a_phase];
}

const RefCountedPtr<EBIndexSpace>&
MultiFluidIndexSpace::getEBIndexSpace(const int a_phase) const
{
  return m_ebis[a_phase];
}

RefCountedPtr<EBIndexSpace>&
MultiFluidIndexSpace::getEBIndexSpace(const phase::which_phase a_phase)
{
  return m_ebis[a_phase];
}

RefCountedPtr<EBIndexSpace>&
MultiFluidIndexSpace::getEBIndexSpace(const int a_phase)
{
  return m_ebis[a_phase];
}

int
MultiFluidIndexSpace::numPhases() const
{
  int phases = 0;
  for (int i = 0; i < m_ebis.size(); i++) {
    if (!m_ebis[i].isNull()) {
      phases++;
    }
  }

  return phases;
}

IntVectSet
MultiFluidIndexSpace::interfaceRegion(const IntVectSet& a_gasCutCells, const IntVectSet& a_solidCutCells) const
{
  CH_TIME("MultiFluidIndexSpace::interfaceRegion");

  Vector<IntVectSet> allGasCutCells(numProc());
  Vector<IntVectSet> allSolidCutCells(numProc());

  gather(allGasCutCells, a_gasCutCells, 0);
  gather(allSolidCutCells, a_solidCutCells, 0);

  IntVectSet ret;

  if (procID() == 0) {
    IntVectSet gasCutCells;
    IntVectSet solidCutCells;

    for (int i = 0; i < allGasCutCells.size(); i++) {
      gasCutCells |= allGasCutCells[i];
    }

    for (int i = 0; i < allSolidCutCells.size(); i++) {
      solidCutCells |= allSolidCutCells[i];
    }

    ret = gasCutCells & solidCutCells;
  }

  broadcast(ret, 0);

  return ret;
}

#include <CD_NamespaceFooter.H>
