/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFHelmholtzDirichletEBBC.cpp
  @brief  Implementation of CD_MFHelmholtzDirichletEBBC.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_MFHelmholtzDirichletEBBC.H>
#include <CD_LeastSquares.H>
#include <CD_BoxLoops.H>
#include <CD_NamespaceHeader.H>

MFHelmholtzDirichletEBBC::MFHelmholtzDirichletEBBC(const int a_phase, const RefCountedPtr<MFHelmholtzJumpBC>& a_jumpBC)
  : MFHelmholtzEBBC(a_phase, a_jumpBC)
{
  CH_TIME("MFHelmholtzDirichletEBBC::MFHelmholtzDirichletEBBC(int, RefCountedPtr<MFHelmholtzJumpBC>)");

  // Default settings.
  m_order  = -1;
  m_weight = -1;

  m_useConstant = false;
  m_useFunction = false;
}

MFHelmholtzDirichletEBBC::~MFHelmholtzDirichletEBBC()
{
  CH_TIME("MFHelmholtzDirichletEBBC::~MFHelmholtzDirichletEBBC()");
}

void
MFHelmholtzDirichletEBBC::setValue(const Real a_value)
{
  CH_TIME("MFHelmholtzDirichletEBBC::setValue(Real)");

  m_useConstant = true;
  m_useFunction = false;

  m_constantValue = a_value;
}

void
MFHelmholtzDirichletEBBC::setValue(const std::function<Real(const RealVect& a_pos)>& a_value)
{
  CH_TIME("MFHelmholtzDirichletEBBC::setValue(std:function<Real(RealVect)>)");

  m_useConstant = false;
  m_useFunction = true;

  m_functionValue = a_value;
}

void
MFHelmholtzDirichletEBBC::setOrder(const int a_order)
{
  CH_TIME("MFHelmholtzDirichletEBBC::setOrder(int)");

  CH_assert(a_order > 0);

  m_order = a_order;
}

void
MFHelmholtzDirichletEBBC::setWeight(const int a_weight)
{
  CH_TIME("MFHelmholtzDirichletEBBC::setWeight(int)");

  CH_assert(a_weight >= 0);

  m_weight = a_weight;
}

void
MFHelmholtzDirichletEBBC::setDomainDropOrder(const int a_domainSize) {
  CH_TIME("MFHelmholtzDirichletEBBC::setDomainDropOrder()");

  m_domainDropOrder = a_domainSize;
}

void
MFHelmholtzDirichletEBBC::defineSinglePhase()
{
  CH_TIME("MFHelmholtzDirichletEBBC::defineSinglePhase()");

  CH_assert(m_order > 0);
  CH_assert(m_weight >= 0);
  CH_assert(m_useConstant || m_useFunction);

  // Also issue run-time errors because those errors can break everything. But I don't see how you would get in a position where that happens.
  if (m_order <= 0 || m_weight < 0) {
    MayDay::Error("MFHelmholtzDirichletEBBC - must have order > 0 and weight >= 0");
  }
  if (!(m_useConstant || m_useFunction)) {
    MayDay::Error("MFHelmholtzDirichletEBBC - not using constant or function!");
  }

  // TLDR: We compute the stencil for reconstructing dphi/dn on the boundary. This is done with least squares reconstruction.
  const DisjointBoxLayout& dbl = m_eblg.getDBL();
  const ProblemDomain& domain = m_eblg.getDomain();
  
  // Drop order if we must
  for (int dir = 0; dir < SpaceDim; dir++) {
    if(domain.size()[dir] <= m_domainDropOrder) {
      m_order = 1;
    }
  }  

  for (DataIterator dit(dbl); dit.ok(); ++dit) {
    const Box         box     = dbl[dit()];
    const EBISBox&    ebisbox = m_eblg.getEBISL()[dit()];
    const IntVectSet& ivs     = ebisbox.getIrregIVS(box);

    BaseIVFAB<Real>&       weights  = m_boundaryWeights[dit()];
    BaseIVFAB<VoFStencil>& stencils = m_kappaDivFStencils[dit()];

    const BaseIVFAB<Real>& Bcoef = (*m_Bcoef)[dit()];

    VoFIterator& singlePhaseVofs = m_jumpBC->getSinglePhaseVofs(m_phase, dit());

    auto kernel = [&](const VolIndex& vof) -> void {
      const Real areaFrac = ebisbox.bndryArea(vof);
      const Real B        = Bcoef(vof, m_comp);

      int                         order;
      bool                        foundStencil = false;
      std::pair<Real, VoFStencil> pairSten;

      // Try quadrants first.
      order = m_order;
      while (!foundStencil && order > 0) {
        foundStencil = this->getLeastSquaresBoundaryGradStencil(pairSten,
                                                                vof,
                                                                VofUtils::Neighborhood::Quadrant,
                                                                dit(),
                                                                order,
                                                                m_weight);
        order--;

        // Check if stencil reaches too far across CF
        if (foundStencil) {
          foundStencil = this->isStencilValidCF(pairSten.second, dit());
        }
      }

      // If we couldn't find in a quadrant, try a larger neighborhood
      order = m_order;
      while (!foundStencil && order > 0) {
        foundStencil = this->getLeastSquaresBoundaryGradStencil(pairSten,
                                                                vof,
                                                                VofUtils::Neighborhood::Radius,
                                                                dit(),
                                                                order,
                                                                m_weight);
        order--;

        // Check if stencil reaches too far across CF
        if (foundStencil) {
          foundStencil = this->isStencilValidCF(pairSten.second, dit());
        }
      }

      if (foundStencil) {
        weights(vof, m_comp)  = pairSten.first;
        stencils(vof, m_comp) = pairSten.second;

        // Stencil and weight must also be scaled by the B-coefficient, dx (because it's used in kappa*Div(F)) and the area fraction.
        weights(vof, m_comp) *= B * areaFrac / m_dx;
        stencils(vof, m_comp) *= B * areaFrac / m_dx;
      }
      else {
        // Dead cell. No flux.
        weights(vof, m_comp) = 0.0;
        stencils(vof, m_comp).clear();
      }
    };

    BoxLoops::loop(singlePhaseVofs, kernel);
  }
}

void
MFHelmholtzDirichletEBBC::applyEBFluxSinglePhase(VoFIterator&     a_singlePhaseVofs,
                                                 EBCellFAB&       a_Lphi,
                                                 const EBCellFAB& a_phi,
                                                 const DataIndex& a_dit,
                                                 const Real&      a_beta,
                                                 const bool&      a_homogeneousPhysBC) const
{
  CH_TIME("MFHelmholtzDirichletEBBC::applyEBFluxSinglePhase(VoFIterator, EBCellFAB, EBCellFAB, DataIndex, Real, bool)");

  // Apply the stencil for computing the contribution to kappaDivF. Note divF is sum(faces) B*grad(Phi)/dx and that this
  // is the contribution from the EB face. B/dx is already included in the stencils and boundary weights, but beta is not.

  // Do single phase cells. Only inhomogeneous contribution here.
  if (!a_homogeneousPhysBC) {

    auto kernel = [&](const VolIndex& vof) -> void {
      Real value;

      if (m_useConstant) {
        value = m_constantValue;
      }
      else if (m_useFunction) {
        const RealVect pos = this->getBoundaryPosition(vof, a_dit);
        value              = m_functionValue(pos);
      }
      else {
        value = 0.0;
        MayDay::Error("MFHelmholtzDirichletEBBC::applyEBFluxSinglePhase - logic bust");
      }

      a_Lphi(vof, m_comp) += a_beta * value * m_boundaryWeights[a_dit](vof, m_comp);
    };

    BoxLoops::loop(a_singlePhaseVofs, kernel);
  }

  return;
}

#include <CD_NamespaceFooter.H>
