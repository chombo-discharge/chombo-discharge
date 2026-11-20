/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_ItoKMCPhysics.cpp
 * @brief  Implementation of CD_ItoKMCPhysics.H
 * @author Robert Marskar
 */

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_ItoKMCPhysics.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::ItoKMC;

thread_local bool                                            ItoKMCPhysics::m_hasKMCSolver;
thread_local KMCSolverType                                   ItoKMCPhysics::m_kmcSolver;
thread_local KMCState                                        ItoKMCPhysics::m_kmcState;
thread_local std::vector<std::shared_ptr<const KMCReaction>> ItoKMCPhysics::m_kmcReactionsThreadLocal;
thread_local KMCState                                        ItoKMCPhysics::m_kmcStateScratch;
thread_local std::vector<Real>                               ItoKMCPhysics::m_kmcPropensityScratch;
thread_local std::vector<long long>                          ItoKMCPhysics::m_weightScratch;
thread_local std::vector<Real>                               ItoKMCPhysics::m_parentWeightScratch;
thread_local std::vector<std::shared_ptr<const KMCReaction>> ItoKMCPhysics::m_kmcReactionsDt;
thread_local std::vector<Real>                               ItoKMCPhysics::m_reactiveDtFactorsDt;
thread_local std::vector<Real>                               ItoKMCPhysics::m_kmcPropensityScratchDt;

Vector<std::string>
ItoKMCPhysics::getPlotVariableNames() const noexcept
{
  CH_TIME("ItoKMCPhysics::getPlotVariableNames");

  return Vector<std::string>(0);
}

Vector<Real>
ItoKMCPhysics::getPlotVariables(const RealVect& /*a_E*/,
                                const RealVect& /*a_pos*/,
                                const Vector<Real>& /*a_phi*/,
                                const Vector<RealVect>& /*a_gradPhi*/,
                                const Real /*a_dx*/,
                                const Real /*a_kappa*/) const noexcept
{
  CH_TIME("ItoKMCPhysics::getPlotVariables");

  return Vector<Real>(0);
}

int
ItoKMCPhysics::getNumberOfPlotVariables() const noexcept
{
  CH_TIME("ItoKMCPhysics::getNumberOfPlotVariables");

  return 0;
}

bool
ItoKMCPhysics::needGradients() const noexcept
{
  CH_TIME("ItoKMCPhysics::needGradients");

  return false;
}

void
ItoKMCPhysics::surfaceEmissionEB(Vector<ParticleSoA<ItoParticle>>&       a_secondaryParticles,
                                 Vector<Real>&                           a_cdrFluxes,
                                 Vector<ParticleSoA<Photon>>&            a_secondaryPhotons,
                                 Vector<Real>&                           a_surfaceDensityChange,
                                 const Vector<ParticleSoA<ItoParticle>>& a_primaryParticles,
                                 const Vector<Real>&                     a_cdrFluxesExtrap,
                                 const Vector<ParticleSoA<Photon>>&      a_primaryPhotons,
                                 const Vector<Real>& /*a_surfaceDensities*/,
                                 const RealVect& a_E,
                                 const RealVect& a_physicalCellCenter,
                                 const RealVect& a_cellCentroid,
                                 const RealVect& a_bndryCentroid,
                                 const RealVect& a_bndryNormal,
                                 const Real      a_bndryArea,
                                 const Real      a_dx,
                                 const Real      a_dt,
                                 const bool      a_isDielectric,
                                 const int       a_matIndex) const noexcept
{
  for (int i = 0; i < a_surfaceDensityChange.size(); i++) {
    a_surfaceDensityChange[i] = 0.0;
  }

  this->secondaryEmissionEB(a_secondaryParticles,
                            a_cdrFluxes,
                            a_secondaryPhotons,
                            a_primaryParticles,
                            a_cdrFluxesExtrap,
                            a_primaryPhotons,
                            a_E,
                            a_physicalCellCenter,
                            a_cellCentroid,
                            a_bndryCentroid,
                            a_bndryNormal,
                            a_bndryArea,
                            a_dx,
                            a_dt,
                            a_isDielectric,
                            a_matIndex);
}

bool
ItoKMCPhysics::emitsWithoutPrimariesEB() const noexcept
{
  CH_TIME("ItoKMCPhysics::emitsWithoutPrimariesEB");

  return true;
}

bool
ItoKMCPhysics::emitsWithoutPrimariesDomain() const noexcept
{
  CH_TIME("ItoKMCPhysics::emitsWithoutPrimariesDomain");

  return true;
}

#include <CD_NamespaceFooter.H>
