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
EBHelmholtzEBBC::define(const Location::Cell                           a_dataLocation,
                        const EBLevelGrid&                             a_eblg,
                        const EBLevelGrid&                             a_eblgFiCo,
                        const RefCountedPtr<LevelData<BaseFab<bool>>>& a_amrValidCells,
                        const RealVect&                                a_probLo,
                        const Real&                                    a_dx,
                        const bool                                     a_hasFineAMRLevel,
                        const bool                                     a_isMGLevel,
                        const int                                      a_refToFine,
                        const int                                      a_ghostCF)
{
  CH_TIME("EBHelmholtzEBBC::define");

  MayDay::Warning("EBHelmholtzEBBC::define -- need routine for computing the EB flux through the EB now!");
  MayDay::Warning("EBHelmholtzEBBC::define -- check if AMRMultiGrid can sync so we don't interpolate twice");

  CH_assert(a_dx > 0.0);
  CH_assert(a_ghostCF >= 0);

  m_dataLocation    = a_dataLocation;
  m_eblg            = a_eblg;
  m_eblgFiCo        = a_eblgFiCo;
  m_probLo          = a_probLo;
  m_dx              = a_dx;
  m_amrValidCells   = a_amrValidCells;
  m_hasFineAMRLevel = a_hasFineAMRLevel;
  m_isMGLevel       = a_isMGLevel;
  m_refToFine       = a_refToFine;
  m_ghostCF         = a_ghostCF;

  if (a_hasFineAMRLevel) {
    if (!(a_eblgFiCo.isDefined())) {
      MayDay::Error("EBHelmholtzEBBC::define - logic bust. Should have defined the refined grids");
    }
    if (a_isMGLevel) {
      MayDay::Error("EBHelmholtzEBBC::define - logic bust. Should not be an MG level (1)");
    }
  }
  if (a_eblgFiCo.isDefined()) {
    if (m_isMGLevel) {
      MayDay::Error("EBHelmholtzEBBC::define - logic bust. Should be an MG level (2)");
    }
    if (!m_hasFineAMRLevel) {
      MayDay::Error("EBHelmholtzEBBC::define - logic bust. Should have an AMR level");
    }
  }

  this->define();
}

const LayoutData<BaseIVFAB<VoFStencil>>&
EBHelmholtzEBBC::getGradPhiRelaxStencils() const noexcept
{
  CH_TIME("EBHelmholtzEBBC::getGradPhiRelaxStencils()");

  return m_gradPhiRelaxStencils;
}

void
EBHelmholtzEBBC::applyEBFluxResid(VoFIterator&           a_vofit,
                                  EBCellFAB&             a_Lphi,
                                  const EBCellFAB&       a_phi,
                                  const EBCellFAB* const a_phiFine,
                                  const BaseIVFAB<Real>& a_Bcoef,
                                  const DataIndex&       a_dit,
                                  const Real&            a_beta,
                                  const bool&            a_homogeneousPhysBC) const
{
  CH_TIME("EBHelmholtzEBBC::applyEBFluxResid");
}

#include <CD_NamespaceFooter.H>
