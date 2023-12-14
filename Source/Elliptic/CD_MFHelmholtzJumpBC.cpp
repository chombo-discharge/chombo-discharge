/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFHelmholtzJumpBC.cpp
  @brief  Implementation of CD_MFHelmholtzJumpBC.H
  @author Robert Marskar
*/

// Chombo includes
#include <EBCellFactory.H>
#include <BaseIVFactory.H>

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
                                     const int            a_ghostCF,
                                     const IntVect        a_ghostPhi)
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
  m_ghostPhi     = a_ghostPhi;

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

const LayoutData<MFInterfaceFAB<VoFStencil>>&
MFHelmholtzJumpBC::getGradPhiStencils() const noexcept
{
  return m_gradPhiStencils;
}

const LayoutData<MFInterfaceFAB<Real>>&
MFHelmholtzJumpBC::getGradPhiWeights() const noexcept
{
  return m_gradPhiWeights;
}

void
MFHelmholtzJumpBC::setBco(const RefCountedPtr<LevelData<MFBaseIVFAB>>& a_Bcoef)
{
  CH_TIME("MFHelmholtzJumpBC::defineStencils()");

  m_Bcoef = a_Bcoef;

  if (m_multiPhase) {
    this->buildAverageStencils();
  }
}
bool
MFHelmholtzJumpBC::isMultiPhase() const noexcept
{
  CH_TIME("MFHelmholtzJumpBC::isMultiPhase");

  return m_multiPhase;
}

void
MFHelmholtzJumpBC::defineStencils()
{
  CH_TIME("MFHelmholtzJumpBC::defineStencils()");

  CH_assert(m_order > 0);
  CH_assert(m_weight >= 0);

  // TLDR: This routine computes the stencils for approximating dphi/dn on each side of the boundary. If we have multi-valued cells we use an average formulation. These
  //       stencils can later be used to compute the bounadry value on the interface.

  // MFHelmholtzJumpBC internals should never be called unless it's a multiphase problem.
  if (m_multiPhase) {
    const DisjointBoxLayout& dbl = m_mflg.getGrids();

    m_gradPhiStencils.define(dbl);
    m_gradPhiWeights.define(dbl);
    m_boundaryPhi.define(dbl);
    m_avgStencils.define(dbl);
    m_avgWeights.define(dbl);
    m_avgVoFs.define(dbl);
    m_denom.define(dbl);

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      const Box box = dbl[dit()];

      m_gradPhiStencils[dit()].define(m_mflg, dit());
      m_gradPhiWeights[dit()].define(m_mflg, dit());
      m_boundaryPhi[dit()].define(m_mflg, dit());
      m_avgStencils[dit()].define(m_mflg, dit());
      m_avgWeights[dit()].define(m_mflg, dit());
      m_denom[dit()].define(m_mflg, dit());
      m_avgVoFs[dit()].define(m_mflg, dit());
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
        BaseIVFAB<VoFStencil>& gradStencils = m_gradPhiStencils[dit()].getIVFAB(iphase);
        BaseIVFAB<Real>&       bndryWeights = m_gradPhiWeights[dit()].getIVFAB(iphase);

        // Iteration space for kernel
        VoFIterator vofit(ivs, ebgraph);

        // Kernel
        auto kernel = [&](const VolIndex& vof) -> void {
          int                         order;
          bool                        foundStencil = false;
          std::pair<Real, VoFStencil> pairSten;

          // Try semi-circle first.
          order = m_order;
          while (!foundStencil && order > 0) {
            foundStencil = this->getLeastSquaresBoundaryGradStencil(pairSten,
                                                                    vof,
                                                                    ebisbox,
                                                                    VofUtils::Neighborhood::SemiCircle,
                                                                    order);
            order--;

            // Check if stencil reaches too far across CF
            if (foundStencil) {
              foundStencil = this->isStencilValidCF(pairSten.second, dit());
            }
          }

          // Try quadrant if that didn't work.
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

          // Last ditch effort: Try a full radius
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
            // const std::string baseErr = "MFHelmholtzJumpBC::defineStencils - dead cell on domain = ";
            // const std::string vofErr  = " on vof = ";
            // const std::string impErr  = " (this may cause multigrid divergence)";

            // std::cout << baseErr << m_mflg.getDomain() << vofErr << vof << impErr << std::endl;

            bndryWeights(vof, m_comp) = 0.0;
            gradStencils(vof, m_comp).clear();
          }
        };

        // Execute kernel and build stencils.
        BoxLoops::loop(vofit, kernel);
      }
    }

    // Build the average stencils.
    this->buildAverageStencils();
    this->resetBC();
  }
}

void
MFHelmholtzJumpBC::buildAverageStencils()
{
  CH_TIME("MFHelmholtzJumpBC::buildAverageStencils()");

  CH_assert(m_multiPhase);

  const DisjointBoxLayout& dbl = m_mflg.getGrids();

  // Compute the average stencils and weights.
  for (int iphase = 0; iphase < m_numPhases; iphase++) {
    const EBLevelGrid& eblg  = m_mflg.getEBLevelGrid(iphase);
    const EBISLayout&  ebisl = eblg.getEBISL();

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      const Box         box     = dbl[dit()];
      const EBISBox&    ebisbox = ebisl[dit()];
      const IntVectSet& ivs     = m_ivs[dit()];

      const BaseIVFAB<VoFStencil>& gradStencils = m_gradPhiStencils[dit()].getIVFAB(iphase);
      const BaseIVFAB<Real>&       bndryWeights = m_gradPhiWeights[dit()].getIVFAB(iphase);

      // Build the average stencils. Only matters if the cell is a multi-valued cell.
      BaseIVFAB<VoFStencil>&       avgStencils = m_avgStencils[dit()].getIVFAB(iphase);
      BaseIVFAB<Real>&             avgWeights  = m_avgWeights[dit()].getIVFAB(iphase);
      BaseIVFAB<Vector<VolIndex>>& avgVoFs     = m_avgVoFs[dit()].getIVFAB(iphase);
      const BaseIVFAB<Real>&       Bcoef       = (*m_Bcoef)[dit()].getIVFAB(iphase);

      for (IVSIterator ivsIt(ivs); ivsIt.ok(); ++ivsIt) {
        const IntVect iv = ivsIt();

        const VolIndex         curVof(iv, 0);
        const Vector<VolIndex> allVofs = ebisbox.getVoFs(iv);

        Real&             curWeight  = avgWeights(curVof, m_comp);
        VoFStencil&       curStencil = avgStencils(curVof, m_comp);
        Vector<VolIndex>& curVoFs    = avgVoFs(curVof, m_comp);

        curWeight = 0.0;
        curStencil.clear();
        curVoFs = allVofs;

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

  // Build the aggregate stencils for performance.
  for (int iphase = 0; iphase < m_numPhases; iphase++) {
    const EBLevelGrid& eblg  = m_mflg.getEBLevelGrid(iphase);
    const EBISLayout&  ebisl = eblg.getEBISL();

    LayoutData<RefCountedPtr<AggStencil<EBCellFAB, BaseIVFAB<Real>>>>& aggStencils = m_aggStencils[iphase];
    aggStencils.define(dbl);

    // Some dummy data for the stencils.
    LevelData<EBCellFAB> phiProxy(dbl, 1, m_ghostPhi, EBCellFactory(ebisl));

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      const EBISBox& ebisbox = ebisl[dit()];
      const EBGraph& ebgraph = ebisbox.getEBGraph();

      Vector<RefCountedPtr<BaseIndex>>   dstBaseIndex;
      Vector<RefCountedPtr<BaseStencil>> dstBaseStencil;

      auto kernel = [&](const VolIndex& vof) -> void {
        const VoFStencil& stencil = (m_avgStencils[dit()].getIVFAB(iphase))(VolIndex(vof.gridIndex(), 0), 0);

        dstBaseIndex.push_back(RefCountedPtr<BaseIndex>(new VolIndex(vof)));
        dstBaseStencil.push_back(RefCountedPtr<BaseStencil>(new VoFStencil(stencil)));
      };

      VoFIterator vofit(m_ivs[dit()], ebgraph);

      BoxLoops::loop(vofit, kernel);

      aggStencils[dit()] = RefCountedPtr<AggStencil<EBCellFAB, BaseIVFAB<Real>>>(
        new AggStencil<EBCellFAB, BaseIVFAB<Real>>(dstBaseIndex,
                                                   dstBaseStencil,
                                                   phiProxy[dit()],
                                                   m_boundaryPhi[dit()].getIVFAB(iphase)));
    }
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

  CH_TIMERS("MFHelmholtzJumpBC::matchBC(patch)");
  CH_TIMER("Setup", t1);
  CH_TIMER("Apply stencils", t2);
  CH_TIMER("Get vofs", t3);
  CH_TIMER("Compute jump", t4);
  CH_TIMER("Compute phiBndry", t5);

  // TLDR: This routine computes the boundary value of phi from an expression b1*dphi/dn1 + b2*dphi/dn2 = sigma, where dphi/dn can be represented as
  //
  //          dphi/dn = wB*phiB + sum[w(i) * phi(i)]
  //
  //       This yields an equation which can be solved for phiB. The solution to that is
  //
  //           phiB = sigma/(b1*w1 + b2*w2) - b1*sum[w1(i) * phi1(i)]/(b1*w1 + b2*w2) - b2*sum[w2(i) * phi2(i)]/(b1*w1 + b2*w2).
  //
  //       Because I'm not crazy, I have stored the term 1/(b1*w1 + b2*w2) in m_denom so we can just multiply it in when we need it. Moreover,
  //       this term has already been multiplied into the stencil weights, which is the reason why we only do the stencil apply below (without dividing
  //       by the above factor).

  CH_START(t1);
  constexpr int vofComp     = 0;
  constexpr int firstPhase  = 0;
  constexpr int secondPhase = 1;

  const EBCellFAB& phiPhase0 = a_phi.getPhase(firstPhase);
  const EBCellFAB& phiPhase1 = a_phi.getPhase(secondPhase);

  const EBISBox& ebisBoxPhase0 = phiPhase0.getEBISBox();
  const EBISBox& ebisBoxPhase1 = phiPhase1.getEBISBox();

  BaseIVFAB<Real>& bndryPhiPhase0 = m_boundaryPhi[a_dit].getIVFAB(firstPhase);
  BaseIVFAB<Real>& bndryPhiPhase1 = m_boundaryPhi[a_dit].getIVFAB(secondPhase);

  const BaseIVFAB<VoFStencil>& avgStencilsPhase0 = m_avgStencils[a_dit].getIVFAB(firstPhase);
  const BaseIVFAB<VoFStencil>& avgStencilsPhase1 = m_avgStencils[a_dit].getIVFAB(secondPhase);

  const BaseIVFAB<Real>& denomFactorPhase0 = m_denom[a_dit].getIVFAB(firstPhase);
  const BaseIVFAB<Real>& denomFactorPhase1 = m_denom[a_dit].getIVFAB(secondPhase);

  const BaseIVFAB<Vector<VolIndex>>& avgVoFsPhase0 = m_avgVoFs[a_dit].getIVFAB(firstPhase);
  const BaseIVFAB<Vector<VolIndex>>& avgVoFsPhase1 = m_avgVoFs[a_dit].getIVFAB(secondPhase);

  const AggStencil<EBCellFAB, BaseIVFAB<Real>>& aggSten0 = *m_aggStencils[firstPhase][a_dit];
  const AggStencil<EBCellFAB, BaseIVFAB<Real>>& aggSten1 = *m_aggStencils[secondPhase][a_dit];
  CH_STOP(t1);

  CH_START(t2);
  aggSten0.apply(bndryPhiPhase0, phiPhase0, 0, false);
  aggSten1.apply(bndryPhiPhase1, phiPhase1, 0, false);
  CH_STOP(t2);

  for (IVSIterator ivsIt(m_ivs[a_dit]); ivsIt.ok(); ++ivsIt) {
    CH_START(t3);
    const IntVect  iv   = ivsIt();
    const VolIndex vof0 = VolIndex(iv, vofComp);

    const Vector<VolIndex>& vofsPhase0 = avgVoFsPhase0(vof0, vofComp);
    const Vector<VolIndex>& vofsPhase1 = avgVoFsPhase1(vof0, vofComp);
    CH_STOP(t3);

    CH_START(t4);
    const Real& denomPhase0 = denomFactorPhase0(vof0, vofComp);

    Real jump = 0.0;
    if (!a_homogeneousPhysBC) {
      for (int i = 0; i < vofsPhase0.size(); i++) {
        jump += a_jump(vofsPhase0[i], m_comp);
      }
      jump *= 1. / vofsPhase0.size();
    }
    CH_STOP(t4);

    CH_START(t5);
    for (int i = 0; i < vofsPhase0.size(); i++) {
      bndryPhiPhase0(vofsPhase0[i], m_comp) += bndryPhiPhase1(vof0, m_comp);
      bndryPhiPhase0(vofsPhase0[i], m_comp) *= -1.0;
    }
    for (int i = 0; i < vofsPhase1.size(); i++) {
      bndryPhiPhase1(vofsPhase1[i], m_comp) = bndryPhiPhase0(vof0, m_comp);
    }

    if (!a_homogeneousPhysBC) {
      for (int i = 0; i < vofsPhase0.size(); i++) {
        bndryPhiPhase0(vofsPhase0[i], m_comp) += denomPhase0 * jump;
      }
      for (int i = 0; i < vofsPhase1.size(); i++) {
        bndryPhiPhase1(vofsPhase1[i], m_comp) += denomPhase0 * jump;
      }
    }

    CH_STOP(t5);
  }
}

#include <CD_NamespaceFooter.H>
