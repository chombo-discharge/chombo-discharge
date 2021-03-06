/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFHelmholtzJumpBC.cpp
  @brief  Implementation of CD_MFHelmholtzJumpBC.H
  @author Robert Marskar
*/

// Our includes
#include <CD_MFHelmholtzJumpBC.H>
#include <CD_VofUtils.H>
#include <CD_LeastSquares.H>
#include <CD_BoxLoops.H>
#include <CD_NamespaceHeader.H>

constexpr int MFHelmholtzJumpBC::m_comp;
constexpr int MFHelmholtzJumpBC::m_nComp;

MFHelmholtzJumpBC::MFHelmholtzJumpBC(const Location::Cell a_dataLocation,
                                     const MFLevelGrid&   a_mflg,
                                     const BcoefPtr&      a_Bcoef,
                                     const Real           a_dx,
                                     const int            a_order,
                                     const int            a_weight,
                                     const int            a_radius,
                                     const int            a_ghostCF)
{
  CH_TIME("MFHelmholtzJumpBC::MFHelmholtzJumpBC");

  CH_assert(a_order > 0);
  CH_assert(a_weight >= 0);
  CH_assert(a_ghostCF > 0);

  m_dataLocation = a_dataLocation;
  m_mflg         = a_mflg;
  m_Bcoef        = a_Bcoef;
  m_dx           = a_dx;
  m_weight       = a_weight;
  m_order        = a_order;
  m_radius       = a_radius;
  m_ghostCF      = a_ghostCF;
  m_numPhases    = m_mflg.numPhases();
  m_multiPhase   = m_numPhases > 1;

  this->defineIterators();
  this->defineStencils();
}

MFHelmholtzJumpBC::~MFHelmholtzJumpBC() { CH_TIME("MFHelmholtzJumpBC::~MFHelmholtzJumpBC()"); }

int
MFHelmholtzJumpBC::getOrder() const
{
  return m_order;
}

int
MFHelmholtzJumpBC::getWeight() const
{
  return m_weight;
}

int
MFHelmholtzJumpBC::getRadius() const
{
  return m_radius;
}

const BaseIVFAB<Real>&
MFHelmholtzJumpBC::getBndryPhi(const int a_phase, const DataIndex& a_dit) const
{
  return m_boundaryPhi[a_dit].getIVFAB(a_phase);
}

VoFIterator&
MFHelmholtzJumpBC::getSinglePhaseVofs(const int a_phase, const DataIndex& a_dit) const
{
  return (*m_singlePhaseVofs.at(a_phase))[a_dit];
}

VoFIterator&
MFHelmholtzJumpBC::getMultiPhaseVofs(const int a_phase, const DataIndex& a_dit) const
{
  return (*m_multiPhaseVofs.at(a_phase))[a_dit];
}

void
MFHelmholtzJumpBC::defineStencils()
{
  CH_TIME("MFHelmholtzJumpBC::defineStencils()");

  CH_assert(m_order > 0);
  CH_assert(m_weight >= 0);

  // TLDR: This routine computes the stencils for approximating dphi/dn on each side of the boundary. If we have multi-valued cells we use an average formulation. These
  //       stencils can later be used to compute the bounadry value on the interface.

  if (m_multiPhase) { // MFHelmholtzJumpBC internals should never be called unless it's a multiphase problem.
    const DisjointBoxLayout& dbl = m_mflg.getGrids();

    m_boundaryPhi.define(dbl);
    m_avgStencils.define(dbl);
    m_avgWeights.define(dbl);
    m_denom.define(dbl);

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      const Box box = dbl[dit()];

      m_boundaryPhi[dit()].define(m_mflg, dit());
      m_avgStencils[dit()].define(m_mflg, dit());
      m_avgWeights[dit()].define(m_mflg, dit());
      m_denom[dit()].define(m_mflg, dit());
    }

    // Build stencils and weights for each phase.
    for (int iphase = 0; iphase < m_numPhases; iphase++) {
      const EBLevelGrid& eblg  = m_mflg.getEBLevelGrid(iphase);
      const EBISLayout&  ebisl = eblg.getEBISL();

      for (DataIterator dit(dbl); dit.ok(); ++dit) {
        const Box         box     = dbl[dit()];
        const EBISBox&    ebisbox = ebisl[dit()];
        const EBGraph&    ebgraph = ebisbox.getEBGraph();
        const IntVectSet& ivs     = m_ivs[dit()];

        // Build stencils like we always do
        BaseIVFAB<VoFStencil> gradStencils(ivs, ebgraph, m_nComp);
        BaseIVFAB<Real>       bndryWeights(ivs, ebgraph, m_nComp);

        // Iteration space for kernel
        VoFIterator vofit(ivs, ebgraph);

        // Kernel
        auto kernel = [&](const VolIndex& vof) -> void {
          int                         order;
          bool                        foundStencil = false;
          std::pair<Real, VoFStencil> pairSten;

          // Try quadrants first.
          order = m_order;
          while (!foundStencil && order > 0) {
            foundStencil =
              this->getLeastSquaresBoundaryGradStencil(pairSten, vof, ebisbox, VofUtils::Neighborhood::Quadrant, order);
            order--;

            // Check if stencil reaches too far across CF
            if (foundStencil) {
              foundStencil = this->isStencilValidCF(pairSten.second, dit());
            }
          }

          // If we couldn't find in a quadrant, try a larger neighborhood
          order = m_order;
          while (!foundStencil && order > 0) {
            foundStencil =
              this->getLeastSquaresBoundaryGradStencil(pairSten, vof, ebisbox, VofUtils::Neighborhood::Radius, order);
            order--;

            // Check if stencil reaches too far across CF
            if (foundStencil) {
              foundStencil = this->isStencilValidCF(pairSten.second, dit());
            }
          }

          if (foundStencil) {
            bndryWeights(vof, m_comp) = pairSten.first;
            gradStencils(vof, m_comp) = pairSten.second;
          }
          else {
            bndryWeights(vof, m_comp) = 0.0;
            gradStencils(vof, m_comp).clear();
          }
        };

        // Execute kernel and build stencils.
        BoxLoops::loop(vofit, kernel);

        // Build the average stencils. Only matters if the cell is a multi-valued cell.
        BaseIVFAB<VoFStencil>& avgStencils = m_avgStencils[dit()].getIVFAB(iphase);
        BaseIVFAB<Real>&       avgWeights  = m_avgWeights[dit()].getIVFAB(iphase);
        const BaseIVFAB<Real>& Bcoef       = (*m_Bcoef)[dit()].getIVFAB(iphase);

        for (IVSIterator ivsIt(ivs); ivsIt.ok(); ++ivsIt) {
          const IntVect iv = ivsIt();

          const VolIndex         curVof(iv, 0);
          const Vector<VolIndex> allVofs = ebisbox.getVoFs(iv);

          Real&       curWeight  = avgWeights(curVof, m_comp);
          VoFStencil& curStencil = avgStencils(curVof, m_comp);

          curWeight = 0.0;
          curStencil.clear();

          Real avgBco = 0.0;
          for (int ivof = 0; ivof < allVofs.size(); ivof++) {
            const VolIndex& vof = allVofs[ivof];

            avgBco += Bcoef(vof, m_comp);
            curWeight += bndryWeights(vof, m_comp);
            curStencil += gradStencils(vof, m_comp);
          }

          const Real invNum = 1. / allVofs.size();

          avgBco *= invNum;
          curWeight *= invNum;
          curStencil *= invNum;

          curStencil *= avgBco;
          curWeight *= avgBco;
        }
      }
    }

    // For efficiency reasons, store 1/(bp*wp + bq*wq). Scale the average stencils by this value as well since we apply it anyways.
    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      constexpr int vofComp     = 0;
      constexpr int firstPhase  = 0;
      constexpr int secondPhase = 1;

      for (IVSIterator ivsIt(m_ivs[dit()]); ivsIt.ok(); ++ivsIt) {
        const IntVect  iv = ivsIt();
        const VolIndex vof0(iv, vofComp);

        VoFStencil& derivStenPhase0 = m_avgStencils[dit()].getIVFAB(firstPhase)(vof0, vofComp);
        VoFStencil& derivStenPhase1 = m_avgStencils[dit()].getIVFAB(secondPhase)(vof0, vofComp);

        const Real& weightPhase0 = m_avgWeights[dit()].getIVFAB(firstPhase)(vof0, vofComp);
        const Real& weightPhase1 = m_avgWeights[dit()].getIVFAB(secondPhase)(vof0, vofComp);

        Real& denomPhase0 = m_denom[dit()].getIVFAB(firstPhase)(vof0, vofComp);
        Real& denomPhase1 = m_denom[dit()].getIVFAB(secondPhase)(vof0, vofComp);

        const Real denom = 1. / (weightPhase0 + weightPhase1);

        denomPhase0 = denom;
        denomPhase1 = denom;

        derivStenPhase0 *= denom;
        derivStenPhase1 *= denom;
      }
    }

    this->resetBC();
  }
}

void
MFHelmholtzJumpBC::defineIterators()
{
  CH_TIME("MFHelmholtzJumpBC::defineIterators()");

  // TLDR: This function defines iterators for iterating over regular cut-cells and over multi-fluid cut-cells. The iterators
  //       must exist for both single-phase and multi-phase vofs. This is true even if we're not actually solving a multiphase
  //       problem because the boundary conditions classes will still need the iterators.

  const DisjointBoxLayout& dbl = m_mflg.getGrids();

  m_ivs.define(dbl);

  for (int iphase = 0; iphase < m_numPhases; iphase++) {
    m_singlePhaseVofs.emplace(iphase, std::make_shared<LayoutData<VoFIterator>>());
    m_multiPhaseVofs.emplace(iphase, std::make_shared<LayoutData<VoFIterator>>());

    LayoutData<VoFIterator>& singlePhaseVofs = *m_singlePhaseVofs.at(iphase);
    LayoutData<VoFIterator>& multiPhaseVofs  = *m_multiPhaseVofs.at(iphase);

    singlePhaseVofs.define(dbl);
    multiPhaseVofs.define(dbl);

    const EBLevelGrid& eblg  = m_mflg.getEBLevelGrid(iphase);
    const EBISLayout&  ebisl = eblg.getEBISL();

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      const Box        box      = dbl[dit()];
      const EBISBox&   ebisbox  = ebisl[dit()];
      const EBGraph&   ebgraph  = ebisbox.getEBGraph();
      const IntVectSet allIrreg = ebisbox.getIrregIVS(box);

      m_ivs[dit()] = m_mflg.interfaceRegion(dbl[dit()], dit());

      IntVectSet singlePhaseCells = allIrreg;
      IntVectSet multiPhaseCells  = IntVectSet();

      singlePhaseCells -= m_ivs[dit()];
      multiPhaseCells |= m_ivs[dit()];

      VoFIterator& singlePhaseVofIt = singlePhaseVofs[dit()];
      VoFIterator& multiPhaseVofIt  = multiPhaseVofs[dit()];

      singlePhaseVofIt.define(singlePhaseCells, ebgraph);
      multiPhaseVofIt.define(multiPhaseCells, ebgraph);
    }
  }
}

bool
MFHelmholtzJumpBC::getLeastSquaresBoundaryGradStencil(std::pair<Real, VoFStencil>& a_stencil,
                                                      const VolIndex&              a_vof,
                                                      const EBISBox&               a_ebisbox,
                                                      const VofUtils::Neighborhood a_neighborhood,
                                                      const int                    a_order) const
{
  CH_TIME("MFHelmholtzJumpBC::getLeastSquarseBoundaryGradStencil(...)");

  // TLDR: This routine computes a stencil for approximating dphi/dn using least squares gradient reconstruction on the EB centroid. We assume that the value
  //       is a "known term" in the expansion.

  bool       foundStencil = false;
  const bool addStartVof  = false;

  const RealVect normal = a_ebisbox.normal(a_vof);

  const VoFStencil gradStencil = LeastSquares::getBndryGradSten(a_vof,
                                                                a_neighborhood,
                                                                m_dataLocation,
                                                                a_ebisbox,
                                                                m_dx,
                                                                a_order,
                                                                m_weight,
                                                                a_order,
                                                                addStartVof);

  if (gradStencil.size() > 0 && normal != RealVect::Zero) {

    const VoFStencil DphiDnStencil  = LeastSquares::projectGradSten(gradStencil, -normal);
    const Real       boundaryWeight = -LeastSquares::sumAllWeights(DphiDnStencil);

    a_stencil = std::make_pair(boundaryWeight, DphiDnStencil);

    foundStencil = true;
  }

  return foundStencil;
}

void
MFHelmholtzJumpBC::resetBC() const
{
  CH_TIME("MFHelmholtzJumpBC::resetBC()");

  for (DataIterator dit = m_boundaryPhi.dataIterator(); dit.ok(); ++dit) {
    for (int iphase = 0; iphase < m_numPhases; iphase++) {
      BaseIVFAB<Real>& bndryPhi = m_boundaryPhi[dit()].getIVFAB(iphase);
      bndryPhi.setVal(0.0);
    }
  }
}

void
MFHelmholtzJumpBC::matchBC(LevelData<BaseIVFAB<Real>>& a_jump,
                           const LevelData<MFCellFAB>& a_phi,
                           const bool                  a_homogeneousPhysBC) const
{
  CH_TIME("MFHelmholtzJumpBC::matchBC(LD<MFCellFAB>, LD<BaseIVFAB<Real>, bool)");

  if (m_multiPhase) {
    for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit) {
      this->matchBC(a_jump[dit()], a_phi[dit()], a_homogeneousPhysBC, dit());
    }
  }
}

inline void
MFHelmholtzJumpBC::matchBC(BaseIVFAB<Real>& a_jump,
                           const MFCellFAB& a_phi,
                           const bool       a_homogeneousPhysBC,
                           const DataIndex& a_dit) const
{
  CH_assert(m_multiPhase);

  // TLDR: This routine computes the boundary value of phi from an expression b1*dphi/dn1 + b2*dphi/dn2 = sigma, where dphi/dn can be represented as
  //
  //          dphi/dn = wB*phiB + sum[w(i) * phi(i)]
  //
  //       This yields an equation which can be solved for phiB. The solution to that is
  //
  //           phiB = sigma/(b1*w1 + b2*w2) - b1*sum[w1(i) * phi1(i)]/(b1*w1 + b2*w2) - b2*sum[w2(i) * phi2(i)]/(b1*w1 + b2*w2).
  //
  //       Because I'm not crazy, I have stored the term 1/(b1*w1 + b2*w2) in m_denom so we can just multiply it in when we need it. Moreover,
  //       this term has already been multiplied into the stencil weights, which is the reason why we only do applyStencil below (without dividing
  //       by the above factor).

  constexpr int vofComp     = 0;
  constexpr int firstPhase  = 0;
  constexpr int secondPhase = 1;

  const EBCellFAB& phiPhase0 = a_phi.getPhase(firstPhase);
  const EBCellFAB& phiPhase1 = a_phi.getPhase(secondPhase);

  const EBISBox& ebisBoxPhase0 = phiPhase0.getEBISBox();
  const EBISBox& ebisBoxPhase1 = phiPhase1.getEBISBox();

  BaseIVFAB<Real>& bndryPhiPhase0 = m_boundaryPhi[a_dit].getIVFAB(firstPhase);
  BaseIVFAB<Real>& bndryPhiPhase1 = m_boundaryPhi[a_dit].getIVFAB(secondPhase);

  for (IVSIterator ivsIt(m_ivs[a_dit]); ivsIt.ok(); ++ivsIt) {
    const IntVect iv = ivsIt();

    const VolIndex vof0 = VolIndex(iv, vofComp);

    const VoFStencil& derivStenPhase0 = m_avgStencils[a_dit].getIVFAB(firstPhase)(vof0, vofComp);
    const VoFStencil& derivStenPhase1 = m_avgStencils[a_dit].getIVFAB(secondPhase)(vof0, vofComp);

    const Real& denomPhase0 = m_denom[a_dit].getIVFAB(firstPhase)(vof0, vofComp);
    //    const Real& denomPhase1 = m_denom[a_dit].getIVFAB(secondPhase)(vof0, vofComp);

    // Compute the average jump.
    Real jump = 0.0;
    if (!a_homogeneousPhysBC) {
      Vector<VolIndex> vofs = ebisBoxPhase0.getVoFs(iv);
      for (const auto& v : vofs.stdVector())
        jump += a_jump(v, m_comp);
      jump *= 1. / vofs.size();
    }

    // Do the matching
    const Real contribPhase0 = this->applyStencil(derivStenPhase0, phiPhase0);
    const Real contribPhase1 = this->applyStencil(derivStenPhase1, phiPhase1);

    const Real phiBndry = denomPhase0 * jump - (contribPhase0 + contribPhase1);

    // Copy the result to the individual phase data holders
    Vector<VolIndex> vofsPhase0 = ebisBoxPhase0.getVoFs(iv);
    Vector<VolIndex> vofsPhase1 = ebisBoxPhase1.getVoFs(iv);

    for (const auto& v : vofsPhase0.stdVector())
      bndryPhiPhase0(v, m_comp) = phiBndry;
    for (const auto& v : vofsPhase1.stdVector())
      bndryPhiPhase1(v, m_comp) = phiBndry;
  }
}

#include <CD_NamespaceFooter.H>
