/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzRobinEBBC.cpp
  @brief  Implementation of CD_EBHelmholtzRobinEBBC.H
  @author Robert Marskar
  @todo   Add least squares implementation of extrapolation stuff. 
*/

// Chombo includes
#include <EBArith.H>

// Our includes
#include <CD_LeastSquares.H>
#include <CD_EBHelmholtzRobinEBBC.H>
#include <CD_BoxLoops.H>
#include <CD_NamespaceHeader.H>

EBHelmholtzRobinEBBC::EBHelmholtzRobinEBBC()
{
  CH_TIME("EBHelmholtzRobinEBBC::EBHelmholtzRobinEBBC()");

  m_order           = -1;
  m_weight          = -1;
  m_domainDropOrder = -1;

  m_useConstant = false;
  m_useFunction = false;
}

EBHelmholtzRobinEBBC::EBHelmholtzRobinEBBC(const int  a_order,
                                           const int  a_weight,
                                           const Real a_A,
                                           const Real a_B,
                                           const Real a_C)
  : EBHelmholtzRobinEBBC()
{
  CH_TIME("EBHelmholtzRobinEBBC::EBHelmholtzRobinEBBC(int, int, Real, Real, Real)");

  CH_assert(a_order > 0);
  CH_assert(a_weight >= 0);

  this->setOrder(a_order);
  this->setWeight(a_weight);
  this->setCoefficients(a_A, a_B, a_C);
}

EBHelmholtzRobinEBBC::EBHelmholtzRobinEBBC(const int                                         a_order,
                                           const int                                         a_weight,
                                           const std::function<Real(const RealVect& a_pos)>& a_A,
                                           const std::function<Real(const RealVect& a_pos)>& a_B,
                                           const std::function<Real(const RealVect& a_pos)>& a_C)
  : EBHelmholtzRobinEBBC()
{
  CH_TIME("EBHelmholtzRobinEBBC::setOrder(int, int, 3xstd::function<Real(RealVect)>)");

  CH_assert(a_order > 0);
  CH_assert(a_weight >= 0);

  this->setOrder(a_order);
  this->setWeight(a_weight);
  this->setCoefficients(a_A, a_B, a_C);
}

EBHelmholtzRobinEBBC::~EBHelmholtzRobinEBBC() { CH_TIME("EBHelmholtzRobinEBBC::~EBHelmholtzRobinEBBC()"); }

void
EBHelmholtzRobinEBBC::setOrder(const int a_order)
{
  CH_TIME("EBHelmholtzRobinEBBC::setOrder(int)");

  CH_assert(a_order > 0);

  m_order = a_order;
}

void
EBHelmholtzRobinEBBC::setWeight(const int a_weight)
{
  CH_TIME("EBHelmholtzRobinEBBC::setWeight(int)");

  CH_assert(a_weight >= 0);

  m_weight = a_weight;
}

void
EBHelmholtzRobinEBBC::setDomainDropOrder(const int a_domainSize)
{
  CH_TIME("EBHelmholtzRobinEBBC::setDomainDropOrder()");

  m_domainDropOrder = a_domainSize;
}

void
EBHelmholtzRobinEBBC::setCoefficients(const Real a_A, const Real a_B, const Real a_C)
{
  CH_TIME("EBHelmholtzRobinEBBC::setCoefficients(Real, Real, Real)");

  m_constantA = a_A;
  m_constantB = a_B;
  m_constantC = a_C;

  m_useConstant = true;
  m_useFunction = false;
}

void
EBHelmholtzRobinEBBC::setCoefficients(const std::function<Real(const RealVect& a_pos)>& a_A,
                                      const std::function<Real(const RealVect& a_pos)>& a_B,
                                      const std::function<Real(const RealVect& a_pos)>& a_C)
{
  CH_TIME("EBHelmholtzRobinEBBC::setCoefficients(3xstd::function<Real(RealVect)>)");

  m_functionA = a_A;
  m_functionB = a_B;
  m_functionC = a_C;

  m_useConstant = false;
  m_useFunction = true;
}

void
EBHelmholtzRobinEBBC::define()
{
  CH_TIME("EBHelmholtzRobinEBBC::define()");

  CH_assert(m_order > 0);
  CH_assert(m_weight >= 0);
  CH_assert(m_useConstant || m_useFunction);

  // Also issue a run-time error because these errors will break everything.
  if (!(m_useConstant || m_useFunction)) {
    MayDay::Error("EBHelmholtzRobinEBBC::define() - not using constant or function!");
  }

  const DisjointBoxLayout& dbl    = m_eblg.getDBL();
  const ProblemDomain&     domain = m_eblg.getDomain();

  // Drop order if we must
  for (int dir = 0; dir < SpaceDim; dir++) {
    if (domain.size()[dir] <= m_domainDropOrder) {
      m_order = 1;
    }
  }

  m_gradPhiStencils.define(dbl);

  for (DataIterator dit(dbl); dit.ok(); ++dit) {
    const Box         box     = dbl[dit()];
    const EBISBox&    ebisbox = m_eblg.getEBISL()[dit()];
    const EBGraph&    ebgraph = ebisbox.getEBGraph();
    const IntVectSet& ivs     = ebisbox.getIrregIVS(box);

    BaseIVFAB<VoFStencil>& stencils = m_gradPhiStencils[dit()];

    stencils.define(ivs, ebgraph, m_nComp);

    // Kernel space
    VoFIterator vofit(ivs, ebgraph);

    auto kernel = [&](const VolIndex& vof) -> void {
      const Real areaFrac = ebisbox.bndryArea(vof);
      const Real helmBco  = (*m_Bcoef)[dit()](vof, m_comp);

      VoFStencil& fluxStencil = stencils(vof, m_comp);

      int  order;
      bool foundStencil = false;

      // Try to find a stencil which uses quadrant-based interpolation, using only the vofs
      // that fall within the quadrant that the cut-cell normal points into.
      order = m_order;
      while (!foundStencil && order > 0) {
        fluxStencil = this->getInterpolationStencil(vof, dit(), VofUtils::Neighborhood::Quadrant, order);
        order--;

        // Check that the stencil doesn't reach into ghost cells it shouldn't!
        if (foundStencil) {
          foundStencil = this->isStencilValidCF(fluxStencil, dit());
        }
      }

      // If the above failed we try a larger neighborhood
      order = m_order;
      while (!foundStencil && order > 0) {
        fluxStencil = this->getInterpolationStencil(vof, dit(), VofUtils::Neighborhood::Radius, order);
        order--;

        // Check that the stencil doesn't reach into ghost cells it shouldn't!
        if (foundStencil) {
          foundStencil = this->isStencilValidCF(fluxStencil, dit());
        }
      }

      // The above stencil is a stencil for interpolating to the cut-cell centroid. We must have the stencil in flux form, using
      // the expression A*phi + B*dphi/dn = C. Set Robin BC constants and scale stencil accordingly.
      if (foundStencil) {
        Real A;
        Real B;
        if (m_useConstant) {
          A = m_constantA;
          B = m_constantB;
        }
        else if (m_useFunction) {
          const RealVect pos = this->getBoundaryPosition(vof, dit());
          A                  = m_functionA(pos);
          B                  = m_functionB(pos);
        }
        else {
          A = 0.0;
          B = 0.0;
          MayDay::Error("EBHelmholtzRobinEBBC::define - logic bust");
        }

        // The normal derivative is dphi/dn = (A*phi - C)/B and the (stencil) flux is
        // kappaDivF = area*b*dphidn/Delta x. Scale accordingly.
        if (std::abs(B) > 0.0) {
          fluxStencil *= A * areaFrac * helmBco / (B * m_dx);
        }
        else {
          fluxStencil.clear();
        }
      }
      else { // Dead cell
        fluxStencil.clear();
      }
    };

    // Run kernel
    BoxLoops::loop(vofit, kernel);
  }
}

void
EBHelmholtzRobinEBBC::applyEBFlux(VoFIterator&           a_vofit,
                                  EBCellFAB&             a_Lphi,
                                  const EBCellFAB&       a_phi,
                                  const BaseIVFAB<Real>& a_Bcoef,
                                  const DataIndex&       a_dit,
                                  const Real&            a_beta,
                                  const bool&            a_homogeneousPhysBC) const
{
  CH_TIME("EBHelmholtzRobinEBBC::applyEBFlux(VoFIterator, EBCellFAB, EBCellFAB, DataIndex, Real, bool)");

  CH_assert(m_useFunction || m_useConstant);

  // Recall that the "flux" is kappaDivF = area*dphi/dn/DeltaX where dphi/dn = (A*phi - C)/B. We already have the phi
  // term in the stencil so only need to add -C/B.
  if (!a_homogeneousPhysBC) {

    // Kernel
    auto kernel = [&](const VolIndex& vof) -> void {
      Real B;
      Real C;
      if (m_useConstant) {
        B = m_constantB;
        C = m_constantC;
      }
      else if (m_useFunction) {
        const RealVect pos = this->getBoundaryPosition(vof, a_dit);
        B                  = m_functionB(pos);
        C                  = m_functionC(pos);
      }
      else {
        B = 0.0;
        C = 0.0;
        MayDay::Error("EBHelmholtzRobinEBBC::applyEBFlux - logic bust");
      }

      const EBISBox& ebisbox   = m_eblg.getEBISL()[a_dit];
      const Real     areaFrac  = ebisbox.bndryArea(vof);
      const Real     helmBco   = (*m_Bcoef)[a_dit](vof, m_comp);
      const Real     kappaDivF = -a_beta * helmBco * areaFrac * C / (m_dx * B);

      if (std::abs(B) > 0.0) {
        a_Lphi(vof, m_comp) += kappaDivF;
      }
    };

    BoxLoops::loop(a_vofit, kernel);
  }
}

VoFStencil
EBHelmholtzRobinEBBC::getInterpolationStencil(const VolIndex&              a_vof,
                                              const DataIndex&             a_dit,
                                              const VofUtils::Neighborhood a_neighborhood,
                                              const int                    a_order) const
{
  CH_TIME("EBHelmholtzRobinEBBC::getInterpolationStencil(VolIndex, DataIndex, VofUtils::Neighborhood, int)");

  CH_assert(a_order > 0);

  // TLDR: This routine will compute a stencil for interpolating the mesh data to the embedded boundary centroid using least squares reconstruction
  //       of the solution. The user will input the desired neighborhood and order of that interpolation. By default, the radius of the stencil is
  //       the same as the order.

  const EBISBox& ebisbox     = m_eblg.getEBISL()[a_dit];
  const bool     useStartVof = !(
    m_weight >
    0); // If we use unweighted least squares we can, in fact, include the cut-cell itself in the interpolation.
  const int radius = a_order;

  // Get the vofs around the cut-cell. Note that if m_weight = 0 we enable the cut-cell itself in the interpolation.
  Vector<VolIndex> vofs;
  switch (a_neighborhood) {
  case VofUtils::Neighborhood::Quadrant:
    vofs = VofUtils::getVofsInQuadrant(a_vof,
                                       ebisbox,
                                       ebisbox.normal(a_vof),
                                       radius,
                                       VofUtils::Connectivity::MonotonePath,
                                       useStartVof);
    break;
  case VofUtils::Neighborhood::Radius:
    vofs = VofUtils::getVofsInRadius(a_vof, ebisbox, radius, VofUtils::Connectivity::MonotonePath, useStartVof);
    break;
  default:
    MayDay::Error(
      "EBHelmholtzRobinEBBC::getInterpolationStencil(VolIndex, DataIndex, VofUtils::Neighborhood) -- logic bust");
  }

  // Build displacements vector, i.e. distances from cell centers/centroids to the cut-cell EB centroid.
  const Vector<RealVect> displacements =
    LeastSquares::getDisplacements(Location::Cell::Boundary, m_dataLocation, a_vof, vofs, ebisbox, m_dx);

  // M = Number of unknowns in Taylor expansion of order a_order.
  // K = Number of equations (displacements)
  const int M = LeastSquares::getTaylorExpansionSize(a_order);
  const int K = displacements.size();

  // If we have enough equations we can get an interpolation stencil.
  VoFStencil interpStencil;
  if (K > M) {
    interpStencil = LeastSquares::computeInterpolationStencil(vofs, displacements, m_weight, a_order);
  }

  return interpStencil;
}

#include <CD_NamespaceFooter.H>
