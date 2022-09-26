/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_MFHelmholtzNeumannEBBC.cpp
  @brief  Implementation of CD_MFHelmholtzNeumannEBBC.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_MFHelmholtzNeumannEBBC.H>
#include <CD_BoxLoops.H>
#include <CD_NamespaceHeader.H>

MFHelmholtzNeumannEBBC::MFHelmholtzNeumannEBBC(const int a_phase, const RefCountedPtr<MFHelmholtzJumpBC>& a_jumpBC)
  : MFHelmholtzEBBC(a_phase, a_jumpBC)
{
  CH_TIME("MFHelmholtzNeumannEBBC::MFHelmholtzNeumannEBBC(int, RefCountedPtr<MFHelmholtzJumpBC>)");

  m_useConstant = false;
  m_useFunction = false;
}

MFHelmholtzNeumannEBBC::~MFHelmholtzNeumannEBBC() { CH_TIME("MFHelmholtzNeumannEBBC::~MFHelmholtzNeumannEBBC()"); }

void
MFHelmholtzNeumannEBBC::setDphiDn(const Real a_DphiDn)
{
  CH_TIME("MFHelmholtzNeumannEBBC::setDphiDn(Real)");

  m_multByBco = true;

  m_useConstant = true;
  m_useFunction = false;

  m_constantDphiDn = a_DphiDn;
}

void
MFHelmholtzNeumannEBBC::setDphiDn(const std::function<Real(const RealVect& a_pos)>& a_DphiDn)
{
  CH_TIME("MFHelmholtzNeumannEBBC::setDphiDn(std::function<Real(RealVect)>)");

  m_multByBco = true;

  m_useConstant = false;
  m_useFunction = true;

  m_functionDphiDn = a_DphiDn;
}

void
MFHelmholtzNeumannEBBC::setBxDphiDn(const Real a_BxDphiDn)
{
  CH_TIME("MFHelmholtzNeumannEBBC::setBxDphiDn(Real)");

  this->setDphiDn(a_BxDphiDn);

  m_multByBco = false;
}

void
MFHelmholtzNeumannEBBC::setBxDphiDn(const std::function<Real(const RealVect& a_pos)>& a_BxDphiDn)
{
  CH_TIME("MFHelmholtzNeumannEBBC::setBxDphiDn(std::function<Real(RealVect)>)");

  this->setDphiDn(a_BxDphiDn);

  m_multByBco = false;
}

void
MFHelmholtzNeumannEBBC::defineSinglePhase()
{
  CH_TIME("MFHelmholtzNeumannEBBC::defineSinglePhase()");

  // No stencils to define here, but we must use constant or function-based BCs.
  if (!(m_useConstant || m_useFunction)) {
    MayDay::Error("MFHelmholtzNeumannEBBC::defineSinglePhase - not using constant or function!");
  }
}

void
MFHelmholtzNeumannEBBC::applyEBFluxSinglePhase(VoFIterator&           a_singlePhaseVofs,
                                               EBCellFAB&             a_Lphi,
                                               const EBCellFAB&       a_phi,
                                               const BaseIVFAB<Real>& a_Bcoef,
                                               const DataIndex&       a_dit,
                                               const Real&            a_beta,
                                               const bool&            a_homogeneousPhysBC) const
{
  CH_TIME("MFHelmholtzNeumannEBBC::applyEBFluxSinglePhase(VoFIterator, EBCellFAB, EBCellFAB, DataIndex, Real, bool)");

  // TLDR: For Neumann, we want to add the flux beta*bco*area*(dphi/dn)/dx where the
  //       dx comes from the fact that the term we are computing will be added to kappa*div(F)
  if (!a_homogeneousPhysBC) {

    auto kernel = [&](const VolIndex& vof) -> void {
      Real value = 0.0;
      if (m_useConstant) {
        value = m_constantDphiDn;
      }
      else if (m_useFunction) {
        const RealVect pos = this->getBoundaryPosition(vof, a_dit);
        value              = m_functionDphiDn(pos);
      }

      // B-coefficient, area fraction, and division by dx (from Div(F)) already a part of the boundary weights, but
      // beta is not.
      const EBISBox& ebisbox   = m_eblg.getEBISL()[a_dit];
      const Real     areaFrac  = ebisbox.bndryArea(vof);
      const Real     B         = m_multByBco ? (*m_Bcoef)[a_dit](vof, m_comp) : 1;
      const Real     kappaDivF = a_beta * B * value * areaFrac / m_dx;

      a_Lphi(vof, m_comp) += kappaDivF;
    };

    BoxLoops::loop(a_singlePhaseVofs, kernel);
  }

  return;
}

#include <CD_NamespaceFooter.H>
