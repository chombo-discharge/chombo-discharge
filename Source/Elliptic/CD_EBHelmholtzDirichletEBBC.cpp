/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzDirichletEBBC.cpp
  @brief  Implementation of CD_EBHelmholtzDirichletEBBC.H
  @author Robert Marskar
*/

// Chombo includes
#include <EBArith.H>
#include <CH_Timer.H>
#include <NeighborIterator.H>

// Our includes
#include <CD_LeastSquares.H>
#include <CD_EBHelmholtzDirichletEBBC.H>
#include <CD_BoxLoops.H>
#include <CD_NamespaceHeader.H>

EBHelmholtzDirichletEBBC::EBHelmholtzDirichletEBBC()
{
  CH_TIME("EBHelmholtzDirichletEBBC::EBHelmholtzDirichletEBBC()");

  m_order           = -1;
  m_weight          = -1;
  m_domainDropOrder = 0;

  m_useConstant = false;
  m_useFunction = false;
}

EBHelmholtzDirichletEBBC::~EBHelmholtzDirichletEBBC()
{
  CH_TIME("EBHelmholtzDirichletEBBC::~EBHelmholtzDirichletEBBC()");
}

void
EBHelmholtzDirichletEBBC::setOrder(const int a_order)
{
  CH_TIME("EBHelmholtzDirichletEBBC::setOrder(int)");

  CH_assert(a_order > 0);
  m_order = a_order;
}

void
EBHelmholtzDirichletEBBC::setWeight(const int a_weight)
{
  CH_TIME("EBHelmholtzDirichletEBBC::setWeight(int)");

  CH_assert(a_weight >= 0);

  m_weight = a_weight;
}

void
EBHelmholtzDirichletEBBC::setValue(const int a_value)
{
  CH_TIME("EBHelmholtzDirichletEBBC::setValue(int)");

  m_useConstant = true;
  m_useFunction = false;

  m_constantValue = a_value;
}

void
EBHelmholtzDirichletEBBC::setValue(const std::function<Real(const RealVect& a_pos)>& a_value)
{
  CH_TIME("EBHelmholtzDirichletEBBC::setValue(std::function<Real(RealVect)>)");

  m_useConstant = false;
  m_useFunction = true;

  m_functionValue = a_value;
}

void
EBHelmholtzDirichletEBBC::setDomainDropOrder(const int a_domainSize)
{
  CH_TIME("EBHelmholtzDirichletEBBC::setDomainDropOrder()");

  m_domainDropOrder = a_domainSize;
}

void
EBHelmholtzDirichletEBBC::define()
{
  CH_TIME("EBHelmholtzDirichletEBBC::define()");

  CH_assert(m_order > 0);
  CH_assert(m_weight >= 0);

  // Also issue runtime error.
  if (m_order <= 0 || m_weight < 0) {
    MayDay::Error("EBHelmholtzDirichletEBBC - must have order > 0 and weight >= 0");
  }
  if (!(m_useConstant || m_useFunction)) {
    MayDay::Error("EBHelmholtzDirichletEBBC - not using constant or function!");
  }

  const DisjointBoxLayout& dbl    = m_eblg.getDBL();
  const ProblemDomain&     domain = m_eblg.getDomain();

  // Drop order if we must
  for (int dir = 0; dir < SpaceDim; dir++) {
    if (domain.size()[dir] <= m_domainDropOrder) {
      m_order = 1;
    }
  }

  m_boundaryWeightsRelax.define(dbl);
  m_gradPhiRelaxStencils.define(dbl);

  for (DataIterator dit(dbl); dit.ok(); ++dit) {
    const Box         box     = dbl[dit()];
    const EBISBox&    ebisbox = m_eblg.getEBISL()[dit()];
    const EBGraph&    ebgraph = ebisbox.getEBGraph();
    const IntVectSet& ivs     = ebisbox.getIrregIVS(box);

    BaseIVFAB<Real>&       weights  = m_boundaryWeightsRelax[dit()];
    BaseIVFAB<VoFStencil>& stencils = m_gradPhiRelaxStencils[dit()];

    weights.define(ivs, ebgraph, m_nComp);
    stencils.define(ivs, ebgraph, m_nComp);

    // Iteration space for kernel
    VoFIterator vofit(ivs, ebgraph);

    // Kernel
    auto kernel = [&](const VolIndex& vof) -> void {
      const Real areaFrac = ebisbox.bndryArea(vof);

      int                         order;
      bool                        foundStencil = false;
      std::pair<Real, VoFStencil> pairSten;

      // Try quadrants first.
      order = m_order;
      while (!foundStencil && order > 0) {
        foundStencil = this->getGradStencilSingleLevel(pairSten, vof, VofUtils::Neighborhood::Quadrant, dit(), order);
        order--;

        // Check if stencil reaches too far across CF
        if (foundStencil) {
          foundStencil = this->isStencilValidCF(pairSten.second, dit());
        }
      }

      // If we couldn't find in a quadrant, try a larger neighborhood
      order = m_order;
      while (!foundStencil && order > 0) {
        foundStencil = this->getGradStencilSingleLevel(pairSten, vof, VofUtils::Neighborhood::Radius, dit(), order);
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
        weights(vof, m_comp) *= areaFrac / m_dx;
        stencils(vof, m_comp) *= areaFrac / m_dx;
      }
      else {
        // Dead cell. No flux.
        weights(vof, m_comp) = 0.0;
        stencils(vof, m_comp).clear();
      }
    };

    BoxLoops::loop(vofit, kernel);
  }

  this->defineResidualStencils();
}

void
EBHelmholtzDirichletEBBC::defineResidualStencils() noexcept
{
  CH_TIME("EBHelmholtzDirichletEBBC::defineResidualStencils");

  const DisjointBoxLayout& dbl     = m_eblg.getDBL();
  const DisjointBoxLayout& dblFine = m_eblgFiCo.getDBL();

  const EBISLayout& ebisl     = m_eblg.getEBISL();
  const EBISLayout& ebislFine = m_eblgFiCo.getEBISL();

  const ProblemDomain& domain     = m_eblg.getDomain();
  const ProblemDomain& domainFine = m_eblgFiCo.getDomain();

  m_boundaryWeightsResid.define(dbl);
  m_gradPhiResidStencils.define(dbl);
  m_gradPhiResidStencilsFine.define(dbl);

  for (DataIterator dit(dbl); dit.ok(); ++dit) {

    if (m_hasFineAMRLevel && !m_isMGLevel) {
      BaseIVFAB<Real>&             gradPhiResidWeight       = m_boundaryWeightsResid[dit()];
      BaseIVFAB<VoFStencil>&       gradPhiResidStencils     = m_gradPhiResidStencils[dit()];
      BaseIVFAB<VoFStencil>&       gradPhiResidStencilsFine = m_gradPhiResidStencilsFine[dit()];
      const BaseIVFAB<Real>&       gradPhiRelaxWeight       = m_boundaryWeightsRelax[dit()];
      const BaseIVFAB<VoFStencil>& gradPhiRelaxStencil      = m_gradPhiRelaxStencils[dit()];
      const BaseFab<bool>&         validCells               = (*m_amrValidCells)[dit()];

      // Where to define our residual stencils.
      const IntVectSet& ivs     = gradPhiRelaxStencil.getIVS();
      const EBGraph&    ebgraph = ebisl[dit()].getEBGraph();

      gradPhiResidWeight.define(ivs, ebgraph, m_nComp);
      gradPhiResidStencils.define(ivs, ebgraph, m_nComp);
      gradPhiResidStencilsFine.define(ivs, ebgraph, m_nComp);

      // Go through the relaxation stencils. If it reaches into an invalid cell (which can happen near EBCF),
      // recompute a better stencil for the AMR residual.
      auto kernel = [&](const VolIndex& vof) -> void {
        const IntVect iv = vof.gridIndex();

        Vector<VolIndex> validVoFs;
        Vector<VolIndex> invalidVoFs;

        // Check if we need a new stencil.
        if (validCells(iv)) {
          const VoFStencil& relaxStencil = gradPhiRelaxStencil(vof, 0);

          // Go through the stencil and find the VoFs that are involved.
          for (int i = 0; i < relaxStencil.size(); i++) {
            const VolIndex& ivof = relaxStencil.vof(i);
            if (validCells(ivof.gridIndex())) {
              validVoFs.push_back(ivof);
            }
            else {
              invalidVoFs.push_back(ivof);
            }
          }
        }

        if (invalidVoFs.size() > 0) {
          // We need to compute a new stencil for the finite volume flux, which also uses cells on
          // the fine level.

          // Recompute a two-level stencil, replacing the invalid vofs with corresponding fine-level vofs
          Vector<VolIndex> fineVoFs;
          for (int i = 0; i < invalidVoFs.size(); i++) {
            fineVoFs.append(ebisl.refine(invalidVoFs[i], m_refToFine, dit()));
          }

          // Try to find a new stencil which does not reach into invalid cells.
          bool foundStencil = false;
          int  order        = m_order;

          std::tuple<Real, VoFStencil, VoFStencil> stencil;
          while (!foundStencil && order > 0) {
            foundStencil = this->getGradStencilDualLevel(stencil, vof, validVoFs, fineVoFs, dit(), order);
            order--;

            // Check if stencil reaches too far across CF
            if (foundStencil) {
              foundStencil = this->isStencilValidCF(std::get<1>(stencil), dit());
            }
          }

          // If we couldn't find a stencil, just default to the relaxation stencil.
          if (foundStencil) {
            gradPhiResidWeight(vof, m_comp)       = std::get<0>(stencil);
            gradPhiResidStencils(vof, m_comp)     = std::get<1>(stencil);
            gradPhiResidStencilsFine(vof, m_comp) = std::get<2>(stencil);
          }
          else {
            pout() << "EBHelmholtzDirichletEBBC - could not compute AMR stencil for vof = " << vof << endl;

            gradPhiResidWeight(vof, m_comp)   = gradPhiRelaxWeight(vof, m_comp);
            gradPhiResidStencils(vof, m_comp) = gradPhiRelaxStencil(vof, m_comp);
          }
        }
        else {
          // In this case we can just use the relaxation stencil for the residual.
          gradPhiResidWeight(vof, m_comp)   = gradPhiRelaxWeight(vof, m_comp);
          gradPhiResidStencils(vof, m_comp) = gradPhiRelaxStencil(vof, m_comp);
        }
      };

      VoFIterator vofit(ivs, ebgraph);
      BoxLoops::loop(vofit, kernel);
    }
  }
}

void
EBHelmholtzDirichletEBBC::applyEBFluxRelax(VoFIterator&           a_vofit,
                                           EBCellFAB&             a_Lphi,
                                           const EBCellFAB&       a_phi,
                                           const BaseIVFAB<Real>& a_Bcoef,
                                           const DataIndex&       a_dit,
                                           const Real&            a_beta,
                                           const bool&            a_homogeneousPhysBC) const
{
  CH_TIME("EBHelmholtzDirichletEBBC::applyEBFluxRelax(VoFIterator, EBCellFAB, EBCellFAB, DataIndex, Real, bool)");

  CH_assert(m_useConstant || m_useFunction);

  if (!a_homogeneousPhysBC) {

    auto kernel = [&](const VolIndex& vof) -> void {
      Real       value = 0.0;
      const Real B     = a_Bcoef(vof, m_comp);

      if (m_useConstant) {
        value = m_constantValue;
      }
      else if (m_useFunction) {
        const RealVect pos = this->getBoundaryPosition(vof, a_dit);
        value              = m_functionValue(pos);
      }
      else {
        value = 0.0;

        MayDay::Error("EBHelmholtzDirichletEBBC::applyEBFluxRelax -- logic bust");
      }

      // Area fraction, and division by dx (from Div(F)) already a part of the boundary weights, but
      // beta and Bcoef are not.
      a_Lphi(vof, m_comp) += a_beta * B * value * m_boundaryWeightsRelax[a_dit](vof, m_comp);
    };

    BoxLoops::loop(a_vofit, kernel);
  }

  return;
}

void
EBHelmholtzDirichletEBBC::applyEBFluxResid(VoFIterator&           a_vofit,
                                           EBCellFAB&             a_Lphi,
                                           const EBCellFAB&       a_phi,
                                           const EBCellFAB* const a_phiFine,
                                           const BaseIVFAB<Real>& a_Bcoef,
                                           const DataIndex&       a_dit,
                                           const Real&            a_beta) const
{
  CH_TIME("EBHelmholtzDirichletEBBC::getGradStencilSingleLevel");

  auto kernel = [&](const VolIndex& vof) -> void {
    // Value of phi on the boundary
    Real value = 0.0;

    if (m_useConstant) {
      value = m_constantValue;
    }
    else if (m_useFunction) {
      const RealVect pos = this->getBoundaryPosition(vof, a_dit);
      value              = m_functionValue(pos);
    }
    else {
      value = 0.0;

      MayDay::Error("EBHelmholtzDirichletEBBC::applyEBFluxRelax -- logic bust");
    }

    Real DphiDn = value * m_boundaryWeightsResid[a_dit](vof, m_comp);

    // Apply stencils on this level.
    const VoFStencil& stencil = m_gradPhiResidStencils[a_dit](vof, m_comp);
    for (int i = 0; i < stencil.size(); i++) {
      const VolIndex& ivof    = stencil.vof(i);
      const Real&     iweight = stencil.weight(i);

      DphiDn += iweight * a_phi(ivof, m_comp);
    }

    // Apply stencil on this level.
    if (m_hasFineAMRLevel) {
      const EBCellFAB& phiFine = *a_phiFine;

      const VoFStencil& stencil = m_gradPhiResidStencilsFine[a_dit](vof, m_comp);
      for (int i = 0; i < stencil.size(); i++) {
        const VolIndex& ivof    = stencil.vof(i);
        const Real&     iweight = stencil.weight(i);

        DphiDn += iweight * phiFine(ivof, m_comp);
      }
    }

    a_Lphi += a_beta * a_Bcoef(vof, m_comp) * DphiDn;
  };
}

bool
EBHelmholtzDirichletEBBC::getGradStencilSingleLevel(std::pair<Real, VoFStencil>& a_stencil,
                                                    const VolIndex&              a_vof,
                                                    const VofUtils::Neighborhood a_neighborhood,
                                                    const DataIndex&             a_dit,
                                                    const int                    a_order) const
{
  CH_TIME("EBHelmholtzDirichletEBBC::getGradStencilSingleLevel");

  CH_assert(m_order > 0);
  CH_assert(m_weight >= 0);

  bool foundStencil = false;

  const bool addStartVof = false;

  const EBISBox& ebisbox = m_eblg.getEBISL()[a_dit];
  const RealVect normal  = ebisbox.normal(a_vof);

  const VoFStencil gradientStencil = LeastSquares::getBndryGradSten(a_vof,
                                                                    a_neighborhood,
                                                                    m_dataLocation,
                                                                    ebisbox,
                                                                    m_dx,
                                                                    a_order,
                                                                    m_weight,
                                                                    a_order,
                                                                    addStartVof);

  if (gradientStencil.size() > 0 && normal != RealVect::Zero) {

    const VoFStencil DphiDnStencil  = LeastSquares::projectGradSten(gradientStencil, -normal);
    const Real       boundaryWeight = -LeastSquares::sumAllWeights(DphiDnStencil);

    a_stencil = std::make_pair(boundaryWeight, DphiDnStencil);

    foundStencil = true;
  }

  return foundStencil;
}

bool
EBHelmholtzDirichletEBBC::getGradStencilDualLevel(std::tuple<Real, VoFStencil, VoFStencil>& a_stencil,
                                                  const VolIndex&                           a_vof,
                                                  const Vector<VolIndex>&                   a_coarVoFs,
                                                  const Vector<VolIndex>&                   a_fineVoFs,
                                                  const DataIndex&                          a_dit,
                                                  const int                                 a_order) const
{
  CH_TIME("EBHelmholtz::getGradStencilDualLevel");

  bool foundStencil = false;

  // Compute the displacement vectors from the input vof to the coar vofs and fine vofs.
  const EBISBox& ebisBox     = m_eblg.getEBISL()[a_dit];
  const EBISBox& ebisBoxFine = m_eblgFiCo.getEBISL()[a_dit];

  const Real dx     = m_dx;
  const Real dxFine = m_dx / m_refToFine;

  const RealVect normal = ebisBox.normal(a_vof);

  if (normal != RealVect::Zero) {

    // Compute the least squares approximation to the gradient. This version of getGradSten returns separate fine and coarse
    // stencils. We also prune the boundary value from the equations, so we need to reconstruct the weight for that.
    IntVectSet knownTerms;
    knownTerms |= IntVect::Zero;

    const std::pair<VoFStencil, VoFStencil> stencils = LeastSquares::getGradSten(a_vof,
                                                                                 Location::Cell::Boundary,
                                                                                 m_dataLocation,
                                                                                 a_fineVoFs,
                                                                                 a_coarVoFs,
                                                                                 ebisBoxFine,
                                                                                 ebisBox,
                                                                                 dxFine,
                                                                                 dx,
                                                                                 m_weight,
                                                                                 a_order,
                                                                                 knownTerms);

    const VoFStencil& fineSten = stencils.first;
    const VoFStencil& coarSten = stencils.second;

    if (fineSten.size() == 0 && coarSten.size() == 0) {
      foundStencil = false;
    }
    else {
      foundStencil = true;

      Real boundaryWeight = 0.0;

      const VoFStencil fineDphiDnStencil = LeastSquares::projectGradSten(fineSten, -normal);
      const VoFStencil coarDphiDnStencil = LeastSquares::projectGradSten(coarSten, -normal);

      boundaryWeight += LeastSquares::sumAllWeights(fineDphiDnStencil);
      boundaryWeight += LeastSquares::sumAllWeights(coarDphiDnStencil);

      a_stencil = std::make_tuple(boundaryWeight, coarDphiDnStencil, fineDphiDnStencil);
    }
  }

  return foundStencil;
}

bool
EBHelmholtzDirichletEBBC::getJohansenStencil(std::pair<Real, VoFStencil>& a_stencil,
                                             const VolIndex&              a_vof,
                                             const DataIndex&             a_dit,
                                             const int                    a_order) const
{
  CH_TIME("EBHelmholtzDirichletEBBC::getJohansenStencil(std::pair<Real, VoFStencil>, VolIndex, DataIndex, int)");

  CH_assert(a_order == 1 || m_order == 2);

  // Also issue run-time error if the assertions break.
  if (!(a_order == 1 || a_order == 2)) {
    MayDay::Error("EBHelmholtzDirichletEBBC::getJohansenStencil - only order 1 and 2 is supported!");
  }

  bool foundStencil = false;

  const EBISBox& ebisbox  = m_eblg.getEBISL()[a_dit];
  const RealVect normal   = ebisbox.normal(a_vof);
  const RealVect centroid = ebisbox.bndryCentroid(a_vof);

  bool               dropOrder;
  Vector<VoFStencil> pointStencils;
  Vector<Real>       distanceAlongLine;
  EBArith::dataRayCast(dropOrder,
                       pointStencils,
                       distanceAlongLine,
                       normal,
                       centroid,
                       a_vof,
                       ebisbox,
                       m_dx * RealVect::Unit,
                       IntVectSet(),
                       m_comp,
                       a_order);

  if (!dropOrder) {
    Real       weight;
    VoFStencil stencil;

    if (a_order == 1) {
      const Real x1 = distanceAlongLine[0];

      weight  = 1. / x1;
      stencil = pointStencils[0];
      stencil *= -1. / x1;
    }
    else if (a_order == 2) {
      const Real x1    = distanceAlongLine[0];
      const Real x2    = distanceAlongLine[1];
      const Real denom = x2 * x2 * x1 - x1 * x1 * x2;

      VoFStencil phi1Sten = pointStencils[0];
      VoFStencil phi2Sten = pointStencils[1];

      phi1Sten *= -x2 * x2 / denom;
      phi2Sten *= x1 * x1 / denom;

      weight = -x1 * x1 / denom + x2 * x2 / denom;
      stencil += phi1Sten;
      stencil += phi2Sten;
    }

    a_stencil = std::make_pair(weight, stencil);

    foundStencil = true;
  }

  return foundStencil;
}

bool
EBHelmholtzDirichletEBBC::getChomboLsqStencil(std::pair<Real, VoFStencil>& a_stencil,
                                              const VolIndex&              a_vof,
                                              const DataIndex&             a_dit) const
{
  VoFStencil sten;
  Real       bndryWeight;

  EBArith::getLeastSquaresGradSten(sten,
                                   bndryWeight,
                                   a_vof,
                                   m_eblg.getEBISL()[a_dit],
                                   m_dx * RealVect::Unit,
                                   m_eblg.getDomain(),
                                   m_comp);

  a_stencil = std::make_pair(bndryWeight, sten);

  return (sten.size() > 0);
}

#include <CD_NamespaceFooter.H>
