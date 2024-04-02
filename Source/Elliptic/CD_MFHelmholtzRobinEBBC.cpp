/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_MFHelmholtzRobinEBBC.cpp
  @brief  Implementation of CD_MFHelmholtzRobinEBBC.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_MFHelmholtzRobinEBBC.H>
#include <CD_VofUtils.H>
#include <CD_LeastSquares.H>
#include <CD_BoxLoops.H>
#include <CD_NamespaceHeader.H>

MFHelmholtzRobinEBBC::MFHelmholtzRobinEBBC(const int a_phase, const RefCountedPtr<MFHelmholtzJumpBC>& a_jumpBC)
  : MFHelmholtzEBBC(a_phase, a_jumpBC)
{
  CH_TIME("MFHelmholtzRobinEBBC::MFHelmholtzRobinEBBC(int, RefCountedPtr<MFHelmholtzJumpBC>)");

  m_order  = -1;
  m_weight = -1;

  m_useConstant = false;
  m_useFunction = false;
}

MFHelmholtzRobinEBBC::~MFHelmholtzRobinEBBC()
{
  CH_TIME("MFHelmholtzRobinEBBC::~MFHelmholtzRobinEBBC()");
}

void
MFHelmholtzRobinEBBC::setOrder(const int a_order)
{
  CH_TIME("MFHelmholtzRobinEBBC::setOrder(int)");

  CH_assert(a_order > 0);

  m_order = a_order;
}

void
MFHelmholtzRobinEBBC::setWeight(const int a_weight)
{
  CH_TIME("MFHelmholtzRobinEBBC::setWeight(int)");

  CH_assert(a_weight > 0);

  m_weight = a_weight;
}

void
MFHelmholtzRobinEBBC::setDomainDropOrder(const int a_domainSize)
{
  CH_TIME("MFHelmholtzRobinEBBC::setDomainDropOrder()");

  m_domainDropOrder = a_domainSize;
}

void
MFHelmholtzRobinEBBC::setCoefficients(const Real a_A, const Real a_B, const Real a_C)
{
  CH_TIME("MFHelmholtzRobinEBBC::setCoefficients(Real, Real, Real)");

  m_constantA = a_A;
  m_constantB = a_B;
  m_constantC = a_C;

  m_useConstant = true;
  m_useFunction = false;
}

void
MFHelmholtzRobinEBBC::setCoefficients(const std::function<Real(const RealVect& a_pos)>& a_A,
                                      const std::function<Real(const RealVect& a_pos)>& a_B,
                                      const std::function<Real(const RealVect& a_pos)>& a_C)
{
  CH_TIME("MFHelmholtzRobinEBBC::setCoefficients(3xstd::function<Real(RealVect)>)");

  m_functionA = a_A;
  m_functionB = a_B;
  m_functionC = a_C;

  m_useConstant = false;
  m_useFunction = true;
}

void
MFHelmholtzRobinEBBC::defineSinglePhase()
{
  CH_TIME("MFHelmholtzRobinEBBC::define()");

  CH_assert(m_order > 0);
  CH_assert(m_weight >= 0);
  CH_assert(m_useConstant || m_useFunction);

  // Also issue a run-time error because these errors will break everything.
  if (!(m_useConstant || m_useFunction)) {
    MayDay::Error("MFHelmholtzRobinEBBC::define() - not using constant or function!");
  }
  if (m_order <= 0) {
    MayDay::Error("MFHelmholtzRobinEBBC - must have order > 0");
  }
  if (m_weight < 0) {
    MayDay::Error("MFHelmholtzRobinEBBC - must have weight >= 0");
  }

  const DisjointBoxLayout& dbl    = m_eblg.getDBL();
  const ProblemDomain&     domain = m_eblg.getDomain();
  const DataIterator&      dit    = dbl.dataIterator();

  const int nbox = dit.size();

  // Drop order if we must
  for (int dir = 0; dir < SpaceDim; dir++) {
    if (domain.size()[dir] <= m_domainDropOrder) {
      m_order = 1;
    }
  }

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    const Box      box     = dbl[din];
    const EBISBox& ebisbox = m_eblg.getEBISL()[din];

    BaseIVFAB<Real>&       weights  = m_boundaryWeights[din];
    BaseIVFAB<VoFStencil>& stencils = m_gradPhiStencils[din];

    // Build interpolation stencils.
    VoFIterator& singlePhaseVofs = m_jumpBC->getSinglePhaseVofs(m_phase, din);

    auto kernel = [&](const VolIndex& vof) -> void {
      const Real areaFrac = ebisbox.bndryArea(vof);

      weights(vof, m_comp) = 0.0;

      VoFStencil& fluxStencil = stencils(vof, m_comp);

      int  order;
      bool foundStencil = false;

      // Try to find a stencil which uses quadrant-based interpolation, using only the vofs
      // that fall within the quadrant that the cut-cell normal points into.
      order = m_order;
      while (!foundStencil && order > 0) {
        fluxStencil = this->getInterpolationStencil(vof, din, VofUtils::Neighborhood::Quadrant, order);
        order--;

        // Check that the stencil doesn't reach into ghost cells it shouldn't!
        if (foundStencil) {
          foundStencil = this->isStencilValidCF(fluxStencil, din);
        }
      }

      // If the above failed we try a larger neighborhood
      order = m_order;
      while (!foundStencil && order > 0) {
        fluxStencil = this->getInterpolationStencil(vof, din, VofUtils::Neighborhood::Radius, order);
        order--;

        // Check that the stencil doesn't reach into ghost cells it shouldn't!
        if (foundStencil) {
          foundStencil = this->isStencilValidCF(fluxStencil, din);
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
          const RealVect pos = this->getBoundaryPosition(vof, din);
          A                  = m_functionA(pos);
          B                  = m_functionB(pos);
        }
        else {
          A = 0.0;
          B = 0.0;

          MayDay::Error("MFHelmholtzRobinEBBC::defineSinglePhase - logic bust");
        }

        // The normal derivative is dphi/dn = (A*phi - C)/B and the (stencil) flux is
        // kappaDivF = area*b*dphidn/Delta x. Scale accordingly.
        if (std::abs(B) > 0.0) {
          fluxStencil *= A * areaFrac / (B * m_dx);
        }
        else {
          fluxStencil.clear();
        }
      }
      else {
        // Dead cell. No flux.
        // const std::string baseErr = "MFHelmholtzRobinEBBC::defineSinglePhase - dead cell on domain = ";
        // const std::string vofErr  = " on vof = ";
        // const std::string impErr  = " (this may cause multigrid divergence)";

        // std::cout << baseErr << m_eblg.getDomain() << vofErr << vof << impErr << std::endl;

        fluxStencil.clear();
      }
    };

    BoxLoops::loop(singlePhaseVofs, kernel);
  }
}

void
MFHelmholtzRobinEBBC::applyEBFluxSinglePhase(VoFIterator&           a_singlePhaseVofs,
                                             EBCellFAB&             a_Lphi,
                                             const EBCellFAB&       a_phi,
                                             const BaseIVFAB<Real>& a_Bcoef,
                                             const DataIndex&       a_dit,
                                             const Real&            a_beta,
                                             const bool&            a_homogeneousPhysBC) const
{

  // TLDR: For Robin, the flux is b*dphi/dn = beta*b*A*phi/B - beta*b*C/B and we have stored
  //       the term b*A*phi/B in the interpolation stencil and return it to the operator. The other term we compute below.

  if (!a_homogeneousPhysBC) {
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

        MayDay::Error("MFHelmholtzRobinEBBC::applyEBFluxSinglePhase - logic bust");
      }

      const EBISBox& ebisbox   = m_eblg.getEBISL()[a_dit];
      const Real     areaFrac  = ebisbox.bndryArea(vof);
      const Real     helmBco   = a_Bcoef(vof, m_comp);
      const Real     kappaDivF = -a_beta * helmBco * areaFrac * C / (m_dx * B);

      if (std::abs(B) > 0.0) {
        a_Lphi(vof, m_comp) += kappaDivF;
      }
    };

    BoxLoops::loop(a_singlePhaseVofs, kernel);
  }

  return;
}

VoFStencil
MFHelmholtzRobinEBBC::getInterpolationStencil(const VolIndex&              a_vof,
                                              const DataIndex&             a_dit,
                                              const VofUtils::Neighborhood a_neighborhood,
                                              const int                    a_order) const
{
  CH_TIME("MFHelmholtzRobinEBBC::getInterpolationStencil(VolIndex, DataIndex, VofUtils::Neighborhood, int)");

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
  case VofUtils::Neighborhood::Quadrant: {
    vofs = VofUtils::getVofsInQuadrant(a_vof,
                                       ebisbox,
                                       ebisbox.normal(a_vof),
                                       radius,
                                       VofUtils::Connectivity::MonotonePath,
                                       useStartVof);

    break;
  }
  case VofUtils::Neighborhood::Radius: {
    vofs = VofUtils::getVofsInRadius(a_vof, ebisbox, radius, VofUtils::Connectivity::MonotonePath, useStartVof);

    break;
  }
  default: {
    MayDay::Error(
      "MFHelmholtzRobinEBBC::getInterpolationStencil(VolIndex, DataIndex, VofUtils::Neighborhood) -- logic bust");

    break;
  }
  }

  // Build displacements vector, i.e. distances from cell centers/centroids to the cut-cell EB centroid.
  const Vector<RealVect> displacements = LeastSquares::getDisplacements(Location::Cell::Boundary,
                                                                        m_dataLocation,
                                                                        a_vof,
                                                                        vofs,
                                                                        ebisbox,
                                                                        m_dx);

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
