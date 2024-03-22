/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzNeumannEBBC.cpp
  @brief  Implementation of CD_EBHelmholtzNeumannEBBC.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_EBHelmholtzNeumannEBBC.H>
#include <CD_BoxLoops.H>
#include <CD_NamespaceHeader.H>

EBHelmholtzNeumannEBBC::EBHelmholtzNeumannEBBC()
{
  CH_TIME("EBHelmholtzNeumannEBBC::EBHelmholtzNeumannEBBC()");

  m_multByBco   = true;
  m_useConstant = false;
  m_useFunction = false;
}

EBHelmholtzNeumannEBBC::~EBHelmholtzNeumannEBBC()
{
  CH_TIME("EBHelmholtzNeumannEBBC::~EBHelmholtzNeumannEBBC()");
}

void
EBHelmholtzNeumannEBBC::setDphiDn(const Real a_DphiDn)
{
  CH_TIME("EBHelmholtzNeumannEBBC::setDphiDn(Real)");

  m_multByBco = true;

  m_useConstant = true;
  m_useFunction = false;

  m_constantDphiDn = a_DphiDn;
}

void
EBHelmholtzNeumannEBBC::setDphiDn(const std::function<Real(const RealVect& a_pos)>& a_DphiDn)
{
  CH_TIME("EBHelmholtzNeumannEBBC::setDphiDn(std::function<Real(RealVect)>)");

  m_multByBco = true;

  m_useConstant = false;
  m_useFunction = true;

  m_functionDphiDn = a_DphiDn;
}

void
EBHelmholtzNeumannEBBC::setBxDphiDn(const Real a_BxDphiDn)
{
  CH_TIME("EBHelmholtzNeumannEBBC::setBxDphiDn(Real)");

  this->setDphiDn(a_BxDphiDn);

  m_multByBco = false;
}

void
EBHelmholtzNeumannEBBC::setBxDphiDn(const std::function<Real(const RealVect& a_pos)>& a_BxDphiDn)
{
  CH_TIME("EBHelmholtzNeumannEBBC::setBxDphiDn(std::function<Real(RealVect)>)");

  this->setDphiDn(a_BxDphiDn);

  m_multByBco = false;
}

void
EBHelmholtzNeumannEBBC::define()
{
  CH_TIME("EBHelmholtzNeumannEBBC::define()");

  CH_assert(m_useConstant || m_useFunction);

  // Also issue a run-time error outside of debugging mode.
  if (!(m_useConstant || m_useFunction)) {
    MayDay::Error("EBHelmholtzNeumannEBBC::define - logic bust, not using constant or function.");
  }

  const DisjointBoxLayout& dbl = m_eblg.getDBL();
  const DataIterator&      dit = dbl.dataIterator();

  const int nbox = dit.size();

  m_gradPhiStencils.define(dbl);

  // Reset the stencil everywhere.
#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    const Box         box     = dbl[din];
    const EBISBox&    ebisbox = m_eblg.getEBISL()[din];
    const EBGraph&    ebgraph = ebisbox.getEBGraph();
    const IntVectSet& ivs     = ebisbox.getIrregIVS(box);

    BaseIVFAB<VoFStencil>& stencils = m_gradPhiStencils[din];

    stencils.define(ivs, ebgraph, m_nComp);

    // Iteration space for kernel
    VoFIterator vofit(ivs, ebgraph);

    auto kernel = [&](const VolIndex& vof) -> void {
      stencils(vof, m_comp).clear();
    };

    BoxLoops::loop(vofit, kernel);
  }
}

void
EBHelmholtzNeumannEBBC::applyEBFlux(VoFIterator&           a_vofit,
                                    EBCellFAB&             a_Lphi,
                                    const EBCellFAB&       a_phi,
                                    const BaseIVFAB<Real>& a_Bcoef,
                                    const DataIndex&       a_dit,
                                    const Real&            a_beta,
                                    const bool&            a_homogeneousPhysBC) const
{
  CH_TIME("EBHelmholtzNeumannEBBC::applyEBFlux(VoFIterator, EBCellFAB, EBCellFAB, DataIndex, Real, bool)");

  CH_assert(m_useConstant || m_useFunction);

  // TLDR: For Neumann, we want to add the flux beta*bco*area*(dphi/dn)/dx where the
  //       dx comes from the fact that the term we are computing will be added to kappa*div(F)
  if (!a_homogeneousPhysBC) {
    auto kernel = [&](const VolIndex& vof) -> void {
      Real value;
      if (m_useConstant) {
        value = m_constantDphiDn;
      }
      else if (m_useFunction) {
        const RealVect pos = this->getBoundaryPosition(vof, a_dit);
        value              = m_functionDphiDn(pos);
      }
      else {
        value = 0.0;
        MayDay::Error("EBHelmholtzNeumannEBBC::applyEBFlux - logic bust");
      }

      // B-coefficient, area fraction, and division by dx (from Div(F)) already a part of the boundary weights, but
      // beta is not.
      const EBISBox& ebisbox   = m_eblg.getEBISL()[a_dit];
      const Real     areaFrac  = ebisbox.bndryArea(vof);
      const Real     B         = m_multByBco ? a_Bcoef(vof, m_comp) : 1;
      const Real     kappaDivF = a_beta * B * value * areaFrac / m_dx;

      a_Lphi(vof, m_comp) += kappaDivF;
    };

    BoxLoops::loop(a_vofit, kernel);
  }

  return;
}

#include <CD_NamespaceFooter.H>
