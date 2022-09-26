/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzEBBC.cpp
  @brief  Implementation of CD_EBHelmholtzEBBC.H
  @author Robert Marskar
*/

// Chombo includes
#include <EBCellFactory.H>
#include <CH_Timer.H>

// Our includes
#include <CD_EBHelmholtzEBBC.H>
#include <CD_NamespaceHeader.H>

constexpr int EBHelmholtzEBBC::m_comp;
constexpr int EBHelmholtzEBBC::m_nComp;

EBHelmholtzEBBC::EBHelmholtzEBBC() { CH_TIME("EBHelmholtzEBBC::EBHelmholtzEBBC()"); }

EBHelmholtzEBBC::~EBHelmholtzEBBC() { CH_TIME("EBHelmholtzEBBC::~EBHelmholtzEBBC()"); }

void
EBHelmholtzEBBC::define(const Location::Cell a_dataLocation,
                        const EBLevelGrid&   a_eblg,
                        const RealVect&      a_probLo,
                        const Real&          a_dx,
                        const int            a_ghostCF)
{
  CH_TIME(
    "EBHelmholtzEBBC::define(Location::Cell, EBLevelGrid, RefCountedPtr<LD<BaseIVFAB<Real> > >, RealVect, Real, int)");

  CH_assert(!a_Bcoef.isNull());
  CH_assert(a_dx > 0.0);
  CH_assert(a_ghostCF >= 0);

  m_dataLocation = a_dataLocation;
  m_eblg         = a_eblg;
  m_probLo       = a_probLo;
  m_dx           = a_dx;
  m_ghostCF      = a_ghostCF;

  this->define();
}

const LayoutData<BaseIVFAB<VoFStencil>>&
EBHelmholtzEBBC::getGradPhiStencils() const
{
  CH_TIME("EBHelmholtzEBBC::getGradPhiStencils()");

  return m_gradPhiStencils;
}

#include <CD_NamespaceFooter.H>
