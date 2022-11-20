/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzOp.cpp
  @brief  Implementation of CD_EBHelmholtzOp.H
  @author Robert Marskar
  @todo   Once performance and stability has settled down, remove the debug code in applyOpIrregular
*/

// Chombo includes
#include <ParmParse.H>
#include <EBLevelDataOps.H>
#include <EBCellFactory.H>
#include <CH_Timer.H>

// Our includes
#include <CD_EBHelmholtzOp.H>
#include <CD_LeastSquares.H>
#include <CD_BoxLoops.H>
#include <CD_ParallelOps.H>
#include <CD_NamespaceHeader.H>

constexpr int EBHelmholtzOp::m_nComp;
constexpr int EBHelmholtzOp::m_comp;

EBHelmholtzOp::EBHelmholtzOp(const Location::Cell                             a_dataLocation,
                             const EBLevelGrid&                               a_eblgFine,
                             const EBLevelGrid&                               a_eblg,
                             const EBLevelGrid&                               a_eblgCoFi,
                             const EBLevelGrid&                               a_eblgCoar,
                             const EBLevelGrid&                               a_eblgCoarMG,
                             const RefCountedPtr<EBMultigridInterpolator>&    a_interpolator,
                             const RefCountedPtr<EBFluxRegister>&             a_fluxReg,
                             const RefCountedPtr<EBCoarAve>&                  a_coarAve,
                             const RefCountedPtr<EBHelmholtzDomainBC>&        a_domainBc,
                             const RefCountedPtr<EBHelmholtzEBBC>&            a_ebBc,
                             const RealVect&                                  a_probLo,
                             const Real&                                      a_dx,
                             const int&                                       a_refToFine,
                             const int&                                       a_refToCoar,
                             const bool&                                      a_hasFine,
                             const bool&                                      a_hasCoar,
                             const bool&                                      a_hasMGObjects,
                             const Real&                                      a_alpha,
                             const Real&                                      a_beta,
                             const RefCountedPtr<LevelData<EBCellFAB>>&       a_Acoef,
                             const RefCountedPtr<LevelData<EBFluxFAB>>&       a_Bcoef,
                             const RefCountedPtr<LevelData<BaseIVFAB<Real>>>& a_BcoefIrreg,
                             const IntVect&                                   a_ghostPhi,
                             const IntVect&                                   a_ghostRhs,
                             const Smoother&                                  a_smoother)
  : LevelTGAHelmOp<LevelData<EBCellFAB>, EBFluxFAB>(false), // Time-independent
    m_smoother(a_smoother),
    m_dataLocation(a_dataLocation),
    m_hasMGObjects(a_hasMGObjects),
    m_hasFine(a_hasFine),
    m_hasCoar(a_hasCoar),
    m_refToCoar(a_hasCoar ? a_refToCoar : 1),
    m_refToFine(a_hasFine ? a_refToFine : 1),
    m_ghostPhi(a_ghostPhi),
    m_ghostRhs(a_ghostRhs),
    m_alpha(a_alpha),
    m_beta(a_beta),
    m_dx(a_dx),
    m_probLo(a_probLo),
    m_eblgFine(),
    m_eblg(a_eblg),
    m_eblgCoFi(),
    m_eblgCoar(),
    m_eblgCoarMG(),
    m_domainBc(a_domainBc),
    m_ebBc(a_ebBc),
    m_interpolator(a_interpolator),
    m_fluxReg(a_fluxReg),
    m_coarAve(a_coarAve),
    m_Acoef(a_Acoef),
    m_Bcoef(a_Bcoef),
    m_BcoefIrreg(a_BcoefIrreg)
{

  CH_TIME("EBHelmholtzOp::EBHelmholtzOp(...)");

  CH_assert(!a_Acoef.isNull());
  CH_assert(!a_Bcoef.isNull());
  CH_assert(!a_BcoefIrreg.isNull());

  CH_assert(a_Acoef->nComp() == 1);
  CH_assert(a_Bcoef->nComp() == 1);
  CH_assert(a_BcoefIrreg->nComp() == 1);

  // Default settings. Always solve for comp = 0. If you want something different, copy your
  // input two different data holders before you use AMRMultiGrid.
  m_doInterpCF = true;
  m_doCoarsen  = true;
  m_doExchange = true;
  m_refluxFree = false;
  m_profile    = false;
  m_interval   = Interval(m_comp, m_comp);

  ParmParse pp("EBHelmholtzOp");
  pp.query("reflux_free", m_refluxFree);
  pp.query("profile", m_profile);

  m_timer = Timer("EBHelmholtzOp");

  if (m_hasFine) {
    m_eblgFine = a_eblgFine;
  }

  // If we are using a centroid discretization we must interpolate three ghost cells in the general case. Issue a warning if the
  // interpolator doesn't fill enough ghost cells.
  if (m_dataLocation == Location::Cell::Centroid && m_hasCoar) {
    const int numFill = m_interpolator->getGhostCF();
    if (numFill < 3) {
      MayDay::Warning(
        "EBHelmholtzOp::EBHelmholtzOp() - interpolator is not filling enough ghost cells. Your discretization may be incorrect!");
    }
  }

  // Define restriction and prolongation operators.
  if (m_hasCoar) {
    m_eblgCoFi = a_eblgCoFi;
    m_eblgCoar = a_eblgCoar;

    m_restrictOp.define(m_eblg, m_eblgCoar, m_refToCoar);
    m_prolongOp.define(m_eblg, m_eblgCoar, m_refToCoar);
  }

  if (m_hasMGObjects) {
    constexpr int mgRef = 2;

    m_eblgCoarMG = a_eblgCoarMG;

    m_restrictOpMG.define(m_eblg, m_eblgCoarMG, mgRef);
    m_prolongOpMG.define(m_eblg, m_eblgCoarMG, mgRef);
  }

  // Define BC objects.
  const int ghostCF = m_hasCoar ? m_interpolator->getGhostCF() : 99;
  m_domainBc->define(m_dataLocation, m_eblg, m_probLo, m_dx);
  m_ebBc->define(m_dataLocation, m_eblg, m_probLo, m_dx, ghostCF);

  // Define stencils and compute relaxation terms.
  this->defineStencils();
}

void
EBHelmholtzOp::setAcoAndBco(const RefCountedPtr<LevelData<EBCellFAB>>&       a_Acoef,
                            const RefCountedPtr<LevelData<EBFluxFAB>>&       a_Bcoef,
                            const RefCountedPtr<LevelData<BaseIVFAB<Real>>>& a_BcoefIrreg)
{
  CH_TIME("EBHelmholtzOp::setAcoAndBco()");

  // Set new coefficients and then run defineStencils which will update all the stencils
  // and relaxation weights.
  m_Acoef      = a_Acoef;
  m_Bcoef      = a_Bcoef;
  m_BcoefIrreg = a_BcoefIrreg;

  this->defineStencils();
}

const RefCountedPtr<LevelData<EBCellFAB>>&
EBHelmholtzOp::getAcoef()
{
  return m_Acoef;
}

const RefCountedPtr<LevelData<EBFluxFAB>>&
EBHelmholtzOp::getBcoef()
{
  return m_Bcoef;
}

const RefCountedPtr<LevelData<BaseIVFAB<Real>>>&
EBHelmholtzOp::getBcoefIrreg()
{
  return m_BcoefIrreg;
}

EBHelmholtzOp::~EBHelmholtzOp() { CH_TIME("EBHelmholtzOp::~EBHelmholtzOp()"); }

void
EBHelmholtzOp::turnOffCFInterp()
{
  CH_TIME("EBHelmholtzOp::turnOffCFInterp()");

  m_doInterpCF = false;
}

void
EBHelmholtzOp::turnOnCFInterp()
{
  CH_TIME("EBHelmholtzOp::turnOnCFInterp()");

  m_doInterpCF = true;
}

void
EBHelmholtzOp::turnOffExchange()
{
  CH_TIME("EBHelmholtzOp::turnOffExchange()");

  m_doExchange = false;
}

void
EBHelmholtzOp::turnOnExchange()
{
  CH_TIME("EBHelmholtzOp::turnOnExchange()");

  m_doExchange = true;
}

void
EBHelmholtzOp::turnOffCoarsening()
{
  CH_TIME("EBHelmholtzOp::turnOffCoarsening()");

  m_doCoarsen = false;
}

void
EBHelmholtzOp::turnOnCoarsening()
{
  CH_TIME("EBHelmholtzOp::turnOnCoarsening()");

  m_doCoarsen = true;
}

LevelData<EBFluxFAB>&
EBHelmholtzOp::getFlux() const
{
  CH_TIME("EBHelmholtzOp::getFlux()");

  return m_flux;
}

void
EBHelmholtzOp::defineStencils()
{
  CH_TIME("EBHelmholtzOp::getFlux()");

  // Basic defines.
  EBCellFactory cellFact(m_eblg.getEBISL());
  EBFluxFactory fluxFact(m_eblg.getEBISL());

  m_relCoef.define(m_eblg.getDBL(), m_nComp, IntVect::Zero, cellFact);
  m_flux.define(m_eblg.getDBL(), m_nComp, IntVect::Zero, fluxFact);

  m_vofIterIrreg.define(m_eblg.getDBL());
  m_vofIterMulti.define(m_eblg.getDBL());
  m_vofIterStenc.define(m_eblg.getDBL());
  m_alphaDiagWeight.define(m_eblg.getDBL());
  m_betaDiagWeight.define(m_eblg.getDBL());
  m_relaxStencils.define(m_eblg.getDBL());

  for (int dir = 0; dir < SpaceDim; dir++) {
    m_vofIterDomLo[dir].define(m_eblg.getDBL());
    m_vofIterDomHi[dir].define(m_eblg.getDBL());
    m_centroidFluxStencil[dir].define(m_eblg.getDBL());
  }

  // Get the "colors" for multi-colored relaxation.
  EBArith::getMultiColors(m_colors);

  // First strip of cells on the inside of the computational domain. I.e. the "domain wall" cells where we need BCs.
  for (int dir = 0; dir < SpaceDim; dir++) {
    for (SideIterator sit; sit.ok(); ++sit) {
      Box domainBox = m_eblg.getDomain().domainBox();
      Box sidebox   = adjCellBox(domainBox, dir, sit(), 1);
      sidebox.shift(dir, -sign(sit()));
      m_sideBox.emplace(std::make_pair(dir, sit()), sidebox);
    }
  }

  // This contains the part of the eb flux that contains interior cells.
  const LayoutData<BaseIVFAB<VoFStencil>>& ebFluxStencil = m_ebBc->getGradPhiStencils();

  // Define everything
  for (DataIterator dit(m_eblg.getDBL()); dit.ok(); ++dit) {
    const Box      cellBox = m_eblg.getDBL()[dit()];
    const EBISBox& ebisbox = m_eblg.getEBISL()[dit()];
    const EBGraph& ebgraph = ebisbox.getEBGraph();

    const IntVectSet irregIVS = ebisbox.getIrregIVS(cellBox);
    const IntVectSet multiIVS = ebisbox.getMultiCells(cellBox);

    // Define the cells where we explicitly store stencils. If we use cell-centered data we only
    // need explicit stencils for kappa*L(phi) on the cut-cells. If this is a centroid-based discretization
    // we also need those stencils for cells that share a grid face with a cut-cell. Those cells will have at least one
    // grid face where we can't use centered differencing.
    IntVectSet stencIVS = irregIVS;
    if (m_dataLocation == Location::Cell::Centroid) {

      // Iteration space
      VoFIterator vofit(irregIVS, ebgraph);

      auto kernel = [&](const VolIndex& vof) -> void {
        for (int dir = 0; dir < SpaceDim; dir++) {
          for (SideIterator sit; sit.ok(); ++sit) {
            Vector<VolIndex> neighborVoFs = ebisbox.getVoFs(vof, dir, sit(), 1);

            for (const auto& curVoF : neighborVoFs.stdVector()) {
              if (ebisbox.isIrregular(curVoF.gridIndex())) {
                stencIVS |= curVoF.gridIndex();
              }
            }
          }
        }
      };

      BoxLoops::loop(vofit, kernel);

      stencIVS &= cellBox;
    }

    // Define iterators. These iterators run over the irregular cells, multi-valued cells only, and cells where we have explicit
    // stencils for kappa*L(phi). The domain iterators are iterators for cut-cells that neighbor a domain edge (face)
    m_vofIterIrreg[dit()].define(irregIVS, ebgraph);
    m_vofIterMulti[dit()].define(multiIVS, ebgraph);
    m_vofIterStenc[dit()].define(stencIVS, ebgraph);

    for (int dir = 0; dir < SpaceDim; dir++) {
      const IntVectSet loIrreg = irregIVS & m_sideBox.at(std::make_pair(dir, Side::Lo));
      const IntVectSet hiIrreg = irregIVS & m_sideBox.at(std::make_pair(dir, Side::Hi));

      m_vofIterDomLo[dir][dit()].define(loIrreg, ebgraph);
      m_vofIterDomHi[dir][dit()].define(hiIrreg, ebgraph);
    }

    // Define data holders.
    m_relaxStencils[dit()].define(stencIVS, ebgraph, m_nComp);
    m_alphaDiagWeight[dit()].define(stencIVS, ebgraph, m_nComp);
    m_betaDiagWeight[dit()].define(stencIVS, ebgraph, m_nComp);

    // The below code may seem intimidating at first. What happens is that we explicitly store stencils for all cells that is either a cut-cell
    // or shares a face with a cut-cell. Now, we have to compute stencils explicitly for this subset of cells, and we need representations both
    // of centroid fluxes, i.e. b*grad(phi) (because of refluxing), and also kappa*div(F). The latter is obviously found by summing the finite
    // volume contributions.
    //
    //                    EB
    //                   /
    //    |------------|/ ----------|------------|
    //    |A           /           B|           C|
    //    |           /|            |            |
    //    |          / |            |     x      |
    //    |         /  |       x    |            |
    //    |        /  x|            |            |
    //    |-------/----|------------|------------|
    //    |D     /     |           E|           F|
    //    |     /      |            |            |
    //    |    /       |     x      |     x      |
    //    |   /    x   |            |            |
    //    |  /         |            |            |
    //    |-/----------|------------|------------|
    //     /
    //    EB
    //
    // The sketch above represents the principle behind the discretization. All the faces A-B, B-C, D-E, A-D, B-E require
    // explicit stencils for the fluxes through the faces. This is true even though E is a regular cell and B-C is a regular face.
    //
    // The code computes the required stencils using (essentially) the following procedure:
    //
    //   1. Iterate through all faces in cut-cells, and for a regular cell which shares at least one edge (face in 3D) with a cut-cell.
    //   2. Compute an approximation to the flux on the face centroid
    //   3. Store this flux stencil explicitly (it might be needed for refluxing operations).
    //   4. The face obviously connects two vofs, and if one of those vofs is a part of the subset of cells discussed above, add the flux
    //      contribution to that vof. For example, we need to store a flux stencil for
    //   5. Irregular cells have an extra face, the EB face. Add the flux contribution from that face into storage for kappa*div(F).
    //
    // Referring to the sketch above, we iterate through cells A, B, D, E and compute the face centroid fluxes for all these cells. This includes
    // e.g. the face connecting E and F. However, since F only has regular faces we don't store the stencil explicitly for that cell.
    //
    BaseIVFAB<VoFStencil>& opStencil  = m_relaxStencils[dit()];
    VoFIterator&           vofitStenc = m_vofIterStenc[dit()];
    VoFIterator&           vofitIrreg = m_vofIterIrreg[dit()];

    BoxLoops::loop(vofitStenc, [&](const VolIndex& vof) -> void {
      opStencil(vof, m_comp).clear();
    });

    for (int dir = 0; dir < SpaceDim; dir++) {
      m_centroidFluxStencil[dir][dit()].define(stencIVS, ebgraph, dir, m_nComp);

      BaseIFFAB<VoFStencil>& fluxStencils = m_centroidFluxStencil[dir][dit()];

      // 1.
      FaceIterator faceIt(stencIVS, ebgraph, dir, FaceStop::SurroundingNoBoundary);

      auto kernel = [&](const FaceIndex& face) -> void {
        if (!face.isBoundary()) {
          const VolIndex vofLo = face.getVoF(Side::Lo);
          const VolIndex vofHi = face.getVoF(Side::Hi);

          // 2.
          const VoFStencil fluxSten = this->getFaceCentroidFluxStencil(face, dit());

          // 3.
          fluxStencils(face, m_comp) = fluxSten;

          VoFStencil loKappaDivFSten = fluxSten;
          VoFStencil hiKappaDivFSten = fluxSten;

          // Sign explanation. For vofLo, faceIt() is the face on the "high" side and
          // vice versa for vofHi
          loKappaDivFSten *= ebisbox.areaFrac(face) / m_dx;
          hiKappaDivFSten *= -ebisbox.areaFrac(face) / m_dx;

          // 4. Note that we prune storage here.
          if (stencIVS.contains(vofLo.gridIndex())) {
            opStencil(vofLo, m_comp) += loKappaDivFSten;
          }
          if (stencIVS.contains(vofHi.gridIndex())) {
            opStencil(vofHi, m_comp) += hiKappaDivFSten;
          }
        }
      };

      BoxLoops::loop(faceIt, kernel);
    }

    // 5. Add contributions to the operator from the EB faces.
    BoxLoops::loop(vofitIrreg, [&](const VolIndex& vof) -> void {
      VoFStencil ebSten = ebFluxStencil[dit()](vof, m_comp);

      ebSten *= (*m_BcoefIrreg)[dit()](vof, m_comp);

      opStencil(vof, m_comp) += ebSten;
    });
  }

  // Compute relaxation weights.
  this->computeDiagWeight();
  this->computeRelaxationCoefficient();
  this->makeAggStencil();
}

void
EBHelmholtzOp::setAlphaAndBeta(const Real& a_alpha, const Real& a_beta)
{
  CH_TIME("EBHelmholtzOp::setAlphaAndBeta(Real, Real)");

  m_alpha = a_alpha;
  m_beta  = a_beta;

  // When we change alpha and beta we need to recompute relaxation coefficients...
  this->computeDiagWeight();
  this->computeRelaxationCoefficient();
  this->makeAggStencil();
}

void
EBHelmholtzOp::residual(LevelData<EBCellFAB>&       a_residual,
                        const LevelData<EBCellFAB>& a_phi,
                        const LevelData<EBCellFAB>& a_rhs,
                        bool                        a_homogeneousPhysBC)
{
  CH_TIME("EBHelmholtzOp::residual(LD<EBCellFAB>, LD<EBCellFAB>, LD<EBCellFAB>, bool)");

  // Compute the residual = rhs - L(phi).
  this->applyOp(a_residual, a_phi, nullptr, a_homogeneousPhysBC, true); // residual = L(phi)
  this->axby(a_residual, a_residual, a_rhs, -1.0, 1.0);                 // residual = rhs - L(phi).
}

void
EBHelmholtzOp::preCond(LevelData<EBCellFAB>& a_corr, const LevelData<EBCellFAB>& a_residual)
{
  CH_TIME("EBHelmholtzOp::preCond(LD<EBCellFAB>, LD<EBCellFAB>)");

  EBLevelDataOps::assign(a_corr, a_residual);
  EBLevelDataOps::scale(a_corr, m_relCoef);

  this->relax(a_corr, a_residual, 40);
}

void
EBHelmholtzOp::create(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs)
{
  CH_TIME("EBHelmholtzOp::create(LD<EBCellFAB>, LD<EBCellFAB>)");

  a_lhs.define(m_eblg.getDBL(), a_rhs.nComp(), a_rhs.ghostVect(), EBCellFactory(m_eblg.getEBISL()));
}

void
EBHelmholtzOp::assign(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs)
{
  CH_TIME("EBHelmholtzOp::assign(LD<EBCellFAB>, LD<EBCellFAB>)");

  a_rhs.copyTo(a_lhs);
}

Real
EBHelmholtzOp::dotProduct(const LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs)
{
  CH_TIME("EBHelmholtzOp::dotProduct(LD<EBCellFAB>, LD<EBCellFAB>)");

  ProblemDomain domain;
  Real          volume;

  return EBLevelDataOps::kappaDotProduct(volume, a_lhs, a_rhs, EBLEVELDATAOPS_ALLVOFS, domain);
}

void
EBHelmholtzOp::incr(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs, const Real a_scale)
{
  CH_TIME("EBHelmholtzOp::incr(LD<EBCellFAB>, LD<EBCellFAB>, Real)");

  EBLevelDataOps::incr(a_lhs, a_rhs, a_scale);
}

void
EBHelmholtzOp::axby(LevelData<EBCellFAB>&       a_lhs,
                    const LevelData<EBCellFAB>& a_x,
                    const LevelData<EBCellFAB>& a_y,
                    const Real                  a_a,
                    const Real                  a_b)
{
  CH_TIME("EBHelmholtzOp::axby(LD<EBCellFAB>, LD<EBCellFAB>, LD<EBCellFAB>, Real, Real)");

  EBLevelDataOps::axby(a_lhs, a_x, a_y, a_a, a_b);
}

void
EBHelmholtzOp::scale(LevelData<EBCellFAB>& a_lhs, const Real& a_scale)
{
  CH_TIME("EBHelmholtzOp::axby(LD<EBCellFAB>, Real)");

  EBLevelDataOps::scale(a_lhs, a_scale);
}

Real
EBHelmholtzOp::norm(const LevelData<EBCellFAB>& a_rhs, const int a_order)
{
  CH_TIME("EBHelmholtzOp::axby(LD<EBCellFAB>, int)");

  // TLDR: This computes the Linf norm.
  Real maxNorm = 0.0;

  for (DataIterator dit = a_rhs.dataIterator(); dit.ok(); ++dit) {
    const EBCellFAB& rhs     = a_rhs[dit()];
    const FArrayBox& regRhs  = rhs.getFArrayBox();
    const EBISBox&   ebisbox = rhs.getEBISBox();
    const Box        box     = a_rhs.disjointBoxLayout()[dit()];

    auto regularKernel = [&](const IntVect& iv) -> void {
      if (ebisbox.isRegular(iv)) {
        maxNorm = std::max(maxNorm, std::abs(regRhs(iv, m_comp)));
      }
    };

    auto irregularKernel = [&](const VolIndex& vof) -> void {
      maxNorm = std::max(maxNorm, std::abs(rhs(vof, m_comp)));
    };

    // Launch the kernels.
    BoxLoops::loop(box, regularKernel);
    BoxLoops::loop(m_vofIterIrreg[dit()], irregularKernel);
  }

  maxNorm = ParallelOps::max(maxNorm);

  return maxNorm;
}

void
EBHelmholtzOp::setToZero(LevelData<EBCellFAB>& a_lhs)
{
  CH_TIME("EBHelmholtzOp::setToZero(LD<EBCellFAB>)");

  EBLevelDataOps::setToZero(a_lhs);
}

void
EBHelmholtzOp::createCoarser(LevelData<EBCellFAB>& a_coarse, const LevelData<EBCellFAB>& a_fine, bool a_ghosted)
{
  CH_TIME("EBHelmholtzOp::createCoarser(LD<EBCellFAB>, LD<EBCellFAB>, bool)");

  EBCellFactory factCoar(m_eblgCoarMG.getEBISL());

  a_coarse.define(m_eblgCoarMG.getDBL(), a_fine.nComp(), a_fine.ghostVect(), factCoar);
}

void
EBHelmholtzOp::createCoarsened(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs, const int& a_refRat)
{
  CH_TIME("EBHelmholtzOp::createCoarsened(LD<EBCellFAB>, LD<EBCellFAB>, bool)");

  CH_assert(m_hasCoar);

  EBCellFactory factCoFi(m_eblgCoFi.getEBISL());

  a_lhs.define(m_eblgCoFi.getDBL(), a_rhs.nComp(), a_rhs.ghostVect(), factCoFi);
}

void
EBHelmholtzOp::restrictResidual(LevelData<EBCellFAB>&       a_resCoar,
                                LevelData<EBCellFAB>&       a_phi,
                                const LevelData<EBCellFAB>& a_rhs)
{
  CH_TIME("EBHelmholtzOp::restrictResidual(LD<EBCellFAB>, LD<EBCellFAB>, LD<EBCellFAB>)");

  // Compute the residual on this level first. Make a temporary for that.
  LevelData<EBCellFAB> res;
  this->create(res, a_phi);
  this->setToZero(res);
  this->residual(res, a_phi, a_rhs, true);

  // Restrict it onto the coarse level.
  m_restrictOpMG.restrictResidual(a_resCoar, res, m_interval);
}

void
EBHelmholtzOp::prolongIncrement(LevelData<EBCellFAB>& a_phi, const LevelData<EBCellFAB>& a_correctCoarse)
{
  CH_TIME("EBHelmholtzOp::prolongIncrement(LD<EBCellFAB>, LD<EBCellFAB>)");

  m_prolongOpMG.prolongResidual(a_phi, a_correctCoarse, m_interval);
}

int
EBHelmholtzOp::refToCoarser()
{
  return m_refToCoar;
}

void
EBHelmholtzOp::AMROperator(LevelData<EBCellFAB>&             a_Lphi,
                           const LevelData<EBCellFAB>&       a_phiFine,
                           const LevelData<EBCellFAB>&       a_phi,
                           const LevelData<EBCellFAB>&       a_phiCoar,
                           const bool                        a_homogeneousPhysBC,
                           AMRLevelOp<LevelData<EBCellFAB>>* a_finerOp)
{
  CH_TIME("EBHelmholtzOp::AMROperator(4xLD<EBCellFAB>, bool, AMRLevelOp<LD<EBCellFAB> >*)");

  if (m_profile) {
    m_timer.startEvent("AMROperator");
  }

  if (m_refluxFree) {
    this->refluxFreeAMROperator(a_Lphi, a_phiFine, a_phi, a_phiCoar, a_homogeneousPhysBC, a_finerOp);
  }
  else {
    // Standard approach:
    //
    // If we have a finer level, there is a chance that the stencils from the EB will reach under it, in which case
    // it might fetch potentially bogus data. A clunky way of handling this is to coarsen the data on the fine level
    // first. The best solution would probably be to have the EB flux stencil reach directly into the fine level.
    if (m_hasFine && m_doCoarsen) {
      EBHelmholtzOp* fineOp = (EBHelmholtzOp*)a_finerOp;
      fineOp->coarsen((LevelData<EBCellFAB>&)a_phi, a_phiFine);
    }

    this->applyOp(a_Lphi, a_phi, &a_phiCoar, a_homogeneousPhysBC, false);

    if (m_hasFine) {
      this->reflux(a_Lphi, a_phiFine, a_phi, *a_finerOp);
    }
  }

  if (m_profile) {
    m_timer.stopEvent("AMROperator");
    m_timer.eventReport(pout(), false);
  }
}

void
EBHelmholtzOp::refluxFreeAMROperator(LevelData<EBCellFAB>&             a_Lphi,
                                     const LevelData<EBCellFAB>&       a_phiFine,
                                     const LevelData<EBCellFAB>&       a_phi,
                                     const LevelData<EBCellFAB>&       a_phiCoar,
                                     const bool                        a_homogeneousPhysBC,
                                     AMRLevelOp<LevelData<EBCellFAB>>* a_finerOp)
{
  CH_TIME("EBHelmholtzOp::refluxFreeAMROperator(4xLD<EBCellFAB>, bool, AMRLevelOp<LD<EBCellFAB> >*)");

  // TLDR: This is the "reflux-free" version of the AMROperator which gets rid of the refluxing step and
  //       rather replaces the whole coarse-fine choreopgraphy by conservative flux coarsening. This routine is
  //       a bit more complex than the other version, but that is mostly because the other version reuses the
  //       things I wrote for the relaxation step.

  // If we have a finer level, there is a chance that the stencils from the EB will reach under it, in which case
  // it might fetch potentially bogus data. A clunky way of handling this is to coarsen the data on the fine level
  // first. The best solution would probably be to have the EB flux stencil reach directly into the fine level.
  if (m_hasFine && m_doCoarsen) {
    EBHelmholtzOp* fineOp = (EBHelmholtzOp*)a_finerOp;
    fineOp->coarsen((LevelData<EBCellFAB>&)a_phi, a_phiFine);
  }

  // Make sure this level has updated ghost cells. The guards around these calls
  // are there because MFHelmholtzOp might decide to update ghost cells for us, in
  // which case we won't redo them.
  if (m_doExchange) {
    LevelData<EBCellFAB>& phi = (LevelData<EBCellFAB>&)a_phi;

    phi.exchange();
  }
  if (m_hasCoar && m_doInterpCF) {
    this->inhomogeneousCFInterp((LevelData<EBCellFAB>&)a_phi, a_phiCoar);
  }

  // Compute the centroid-centered fluxes on this level.
  this->computeFlux(a_phi);

  // Finer level must compute fluxes and coarsens them onto this level. We also update the ghost cells,
  // although this MIGHT have been done already we have no guarantee of that. Then we coarsen those fluxes
  // onto this level.
  if (m_hasFine) {
    EBHelmholtzOp* fineOp = (EBHelmholtzOp*)a_finerOp;

    LevelData<EBCellFAB>& phiFine = (LevelData<EBCellFAB>&)a_phiFine;

    phiFine.exchange();
    fineOp->inhomogeneousCFInterp(phiFine, a_phi);

    fineOp->computeFlux(a_phiFine);
    fineOp->coarsen(m_flux, fineOp->getFlux());
  }

  // Fill the domain fluxes.
  for (DataIterator dit(m_eblg.getDBL()); dit.ok(); ++dit) {
    this->fillDomainFlux(m_flux[dit()], a_phi[dit()], m_eblg.getDBL()[dit()], dit());
  }

  // The above calls replaced the fluxes on this level by (conservative) averages of the fluxes on
  // the finer level. We can proceed as usual now, ignoring the fine level.
  const Real inverseDx = 1.0 / m_dx;
  for (DataIterator dit(m_eblg.getDBL()); dit.ok(); ++dit) {
    const Box      cellBox = m_eblg.getDBL()[dit()];
    const EBISBox& ebisbox = m_eblg.getEBISL()[dit()];

    EBCellFAB& Lphi    = a_Lphi[dit()];
    FArrayBox& LphiReg = Lphi.getFArrayBox();

    const EBCellFAB& phi  = a_phi[dit()];
    const EBCellFAB& Aco  = (*m_Acoef)[dit()];
    const EBFluxFAB& Bco  = (*m_Bcoef)[dit()];
    const EBFluxFAB& flux = m_flux[dit()];

    const BaseIVFAB<Real>&       BcoIrreg      = (*m_BcoefIrreg)[dit()];
    const BaseIVFAB<VoFStencil>& ebFluxStencil = m_ebBc->getGradPhiStencils()[dit()];
    const BaseIVFAB<Real>&       alphaDiag     = m_alphaDiagWeight[dit()];

    // Add alpha-term.
    Lphi.setVal(0.0);
    Lphi += Aco;
    Lphi *= phi;
    Lphi *= m_alpha;

    // Add in the beta-term -- b-coefficient and m_beta is already a part of the flux.
    for (int dir = 0; dir < SpaceDim; dir++) {
      const FArrayBox& fluxReg = flux[dir].getFArrayBox();

      auto regularKernel = [&](const IntVect& iv) -> void {
        LphiReg(iv, m_comp) += inverseDx * (fluxReg(iv + BASISV(dir), m_comp) - fluxReg(iv, m_comp));
      };

      BoxLoops::loop(cellBox, regularKernel);
    }

    // Kernel for cut-cells. This will add everything except the inhomogeneous flux through the EB
    // and the domain fluxes in cut-cells.
    auto irregularKernel = [&](const VolIndex& vof) -> void {
      // Stencil that describes the EB flux into the cell. This is the part of the flux due to
      // beta * b * dphi/dn on the EB face.
      const VoFStencil& ebFluxSten = ebFluxStencil(vof, m_comp);

      // Set L(phi) = kappa * alpha * Acoef * phi
      Lphi(vof, m_comp) = alphaDiag(vof, m_comp) * phi(vof, m_comp);

      // Add in the fluxes from the non-EB/non-domain cut-cell faces. When doing this we set
      // L(phi) = kappa * alpha * Acoef * phi + beta * sum_f(b*grad(phi)) where the sum runs
      // over internal faces.
      for (int dir = 0; dir < SpaceDim; dir++) {
        const EBFaceFAB& faceFlux = flux[dir];

        for (SideIterator sit; sit.ok(); ++sit) {
          const Side::LoHiSide side = sit();

          const int               isign = sign(side);
          const Vector<FaceIndex> faces = ebisbox.getFaces(vof, dir, side);

          for (int iface = 0; iface < faces.size(); iface++) {
            const FaceIndex face     = faces[iface];
            const Real      faceArea = ebisbox.areaFrac(face);

            // Distinguish between domain and internal faces; the domain bc object does not multiply the
            // flux by beta.
            Real curFlux = 0.0;
            if (!face.isBoundary()) {
              curFlux = faceFlux(face, m_comp);
            }
            else {
              curFlux = m_domainBc->getFaceFlux(vof, phi, Bco[dir], dir, sit(), dit(), a_homogeneousPhysBC);
              curFlux *= m_beta;
            }

            Lphi(vof, m_comp) += isign * inverseDx * faceArea * curFlux;
          }
        }
      }

      // Add in the flux from the EB face. After this the only term missing the inhomogeneous
      // flux from the EB. I.e. when dphi/dn = w_b * phi_b + sum(w_i * phi(i)) then we are
      // missing only the term w_b * phi_b. This is added by the EBBC object.
      for (int i = 0; i < ebFluxSten.size(); i++) {
        const Real      weight = ebFluxSten.weight(i);
        const VolIndex& ivof   = ebFluxSten.vof(i);

        Lphi(vof, m_comp) += m_beta * BcoIrreg(vof, m_comp) * weight * phi(ivof, m_comp);
      }
    };

    BoxLoops::loop(m_vofIterIrreg[dit()], irregularKernel);

    // Finally, add in the inhomogeneous EB flux.
    m_ebBc->applyEBFlux(m_vofIterIrreg[dit()], Lphi, phi, BcoIrreg, dit(), m_beta, a_homogeneousPhysBC);
  }
}

void
EBHelmholtzOp::AMROperatorNF(LevelData<EBCellFAB>&       a_Lphi,
                             const LevelData<EBCellFAB>& a_phi,
                             const LevelData<EBCellFAB>& a_phiCoar,
                             bool                        a_homogeneousPhysBC)
{
  CH_TIME("EBHelmholtzOp::AMROperatorNF(3xLD<EBCellFAB>, bool)");

  this->applyOp(a_Lphi, a_phi, &a_phiCoar, a_homogeneousPhysBC, false);
}

void
EBHelmholtzOp::AMROperatorNC(LevelData<EBCellFAB>&             a_Lphi,
                             const LevelData<EBCellFAB>&       a_phiFine,
                             const LevelData<EBCellFAB>&       a_phi,
                             bool                              a_homogeneousPhysBC,
                             AMRLevelOp<LevelData<EBCellFAB>>* a_finerOp)
{
  CH_TIME("EBHelmholtzOp::AMROperatorNC(3xLD<EBCellFAB>, bool, AMRLevelOp<LD<EBCellFAB> >)");

  LevelData<EBCellFAB> phiCoar; // Should be safe on the bottom AMR level because only multigrid levels exist below.
  this->AMROperator(a_Lphi, a_phiFine, a_phi, phiCoar, a_homogeneousPhysBC, a_finerOp);
}

void
EBHelmholtzOp::AMRResidual(LevelData<EBCellFAB>&             a_residual,
                           const LevelData<EBCellFAB>&       a_phiFine,
                           const LevelData<EBCellFAB>&       a_phi,
                           const LevelData<EBCellFAB>&       a_phiCoar,
                           const LevelData<EBCellFAB>&       a_rhs,
                           bool                              a_homogeneousPhysBC,
                           AMRLevelOp<LevelData<EBCellFAB>>* a_finerOp)
{
  CH_TIME("EBHelmholtzOp::AMRResidual(5xLD<EBCellFAB>, bool, AMRLevelOp<LD<EBCellFAB> >)");

  // We are computing the residual res = rhs - L(phi). This should include refluxing and all that, and possibly
  // also accounting for EB flux stencils that reach under the fine level (handled in AMROperator).

  this->AMROperator(a_residual, a_phiFine, a_phi, a_phiCoar, a_homogeneousPhysBC, a_finerOp); // a_residual =  L(phi)
  this->axby(a_residual, a_rhs, a_residual, 1., -1.); // a_residual = rhs - -L(phi)
}

void
EBHelmholtzOp::AMRResidualNF(LevelData<EBCellFAB>&       a_residual,
                             const LevelData<EBCellFAB>& a_phi,
                             const LevelData<EBCellFAB>& a_phiCoar,
                             const LevelData<EBCellFAB>& a_rhs,
                             bool                        a_homogeneousPhysBC)
{
  CH_TIME("EBHelmholtzOp::AMRResidualNF(3xLD<EBCellFAB>, bool)");

  // In this case we want to compute the residual on the finest level, knowing that there are no finer levels with which
  // we have to reflux.
  this->AMROperatorNF(a_residual, a_phi, a_phiCoar, a_homogeneousPhysBC); // a_residual = L(phi)
  this->axby(a_residual, a_rhs, a_residual, 1., -1.);                     // a_residual = rhs - L(phi)
}

void
EBHelmholtzOp::AMRResidualNC(LevelData<EBCellFAB>&             a_residual,
                             const LevelData<EBCellFAB>&       a_phiFine,
                             const LevelData<EBCellFAB>&       a_phi,
                             const LevelData<EBCellFAB>&       a_rhs,
                             bool                              a_homogeneousPhysBC,
                             AMRLevelOp<LevelData<EBCellFAB>>* a_finerOp)
{
  CH_TIME("EBHelmholtzOp::AMRResidual(4xLD<EBCellFAB>, bool, AMRLevelOp<LD<EBCellFAB> >)");

  // We are computing the residual res = rhs - L(phi). We know that there are no coarser levels so we don't have to worry
  // about coarse-fine interpolation, but we still have to reflux.

  this->AMROperatorNC(a_residual, a_phiFine, a_phi, a_homogeneousPhysBC, a_finerOp); // a_residual = L(phi)
  this->axby(a_residual, a_rhs, a_residual, 1., -1.);                                // a_residual = rhs - L(phi)
}

void
EBHelmholtzOp::AMRRestrict(LevelData<EBCellFAB>&       a_residualCoarse,
                           const LevelData<EBCellFAB>& a_residual,
                           const LevelData<EBCellFAB>& a_correction,
                           const LevelData<EBCellFAB>& a_coarseCorrection,
                           bool                        a_skip_res)
{
  CH_TIME("EBHelmholtzOp::AMRRestrict(4xLD<EBCellFAB>, bool)");

  // TLDR: We want to restrict the residual onto the coarser level. So we just compute it on this level
  //       and have our restriction operator do the rest. Note that we DO have to worry about coarse-fine
  //       interpolation here (done in applyOp).

  LevelData<EBCellFAB> resThisLevel;
  this->create(resThisLevel, a_residual);

  constexpr bool homogeneousPhysBC = true;
  constexpr bool homogeneousCFBC   = false;

  // We should average a_residual - L(correction, coarCorrection).
  this->applyOp(resThisLevel, a_correction, &a_coarseCorrection, homogeneousPhysBC, homogeneousCFBC);
  this->incr(resThisLevel, a_residual, -1.0);
  this->scale(resThisLevel, -1.0);

  // Restrict residual.
  m_restrictOp.restrictResidual(a_residualCoarse, resThisLevel, m_interval);
}

void
EBHelmholtzOp::AMRProlong(LevelData<EBCellFAB>& a_correction, const LevelData<EBCellFAB>& a_coarseCorrection)
{
  CH_TIME("EBHelmholtzOp::AMRProlong(LD<EBCellFAB>, LD<EBCellFAB>)");

  m_prolongOp.prolongResidual(a_correction, a_coarseCorrection, m_interval);
}

void
EBHelmholtzOp::AMRUpdateResidual(LevelData<EBCellFAB>&       a_residual,
                                 const LevelData<EBCellFAB>& a_correction,
                                 const LevelData<EBCellFAB>& a_coarseCorrection)
{
  CH_TIME("EBHelmholtzOp::AMRUpdateResidual(LD<EBCellFAB>, LD<EBCellFAB>, LD<EBCellFAB>)");

  constexpr bool homogeneousPhysBC = true;
  constexpr bool homogeneousCFBC   = false;

  LevelData<EBCellFAB> lcorr;
  this->create(lcorr, a_correction);

  // lcorr = L(phi, phiCoar)
  this->applyOp(lcorr, a_correction, &a_coarseCorrection, homogeneousPhysBC, homogeneousCFBC);

  // a_residual = a_residual - L(phi)
  this->incr(a_residual, lcorr, -1.0);
}

void
EBHelmholtzOp::applyOp(LevelData<EBCellFAB>& a_Lphi, const LevelData<EBCellFAB>& a_phi, bool a_homogeneousPhysBC)
{
  CH_TIME("EBHelmholtzOp::applyOp(LD<EBCellFAB>, LD<EBCellFAB>, bool)");

  this->applyOp(a_Lphi, a_phi, nullptr, a_homogeneousPhysBC, true);
}

void
EBHelmholtzOp::applyOp(LevelData<EBCellFAB>&             a_Lphi,
                       const LevelData<EBCellFAB>&       a_phi,
                       const LevelData<EBCellFAB>* const a_phiCoar,
                       const bool                        a_homogeneousPhysBC,
                       const bool                        a_homogeneousCFBC)
{
  CH_TIME("EBHelmholtzOp::applyOp(LD<EBCellFAB>, LD<EBCellFAB>, LD<EBCellFAB>, bool, bool)");

  // TLDR: We are computing L(phi) on a level, and there's possibly a coarser level here as well. If there is,
  //       we need to recompute the ghost cells.

  // Nasty cast, but there's no guarantee that a_phi came in with valid ghost cells. We could
  // do a local copy, but that can end up being expensive since this is called on every relaxation.
  LevelData<EBCellFAB>& phi = (LevelData<EBCellFAB>&)a_phi;

  if (m_doExchange) {
    phi.exchange();
  }

  if (m_hasCoar && m_doInterpCF) {
    this->interpolateCF(phi, a_phiCoar, a_homogeneousCFBC);
  }

  // Apply operator in each kernel.
  const DisjointBoxLayout& dbl = a_Lphi.disjointBoxLayout();
  for (DataIterator dit(dbl); dit.ok(); ++dit) {
    this->applyOp(a_Lphi[dit()], phi[dit()], dbl[dit()], dit(), a_homogeneousPhysBC);
  }
}

void
EBHelmholtzOp::applyOp(EBCellFAB&       a_Lphi,
                       EBCellFAB&       a_phi,
                       const Box&       a_cellBox,
                       const DataIndex& a_dit,
                       const bool       a_homogeneousPhysBC)
{
  CH_TIME("EBHelmholtzOp::applyOp(EBCellFAB, EBCellFAB, Box, DataIndex, bool)");

  // TLDR: This is the kernel version of applyOp. It first applies the operator in the regular cells using the standard
  //       5/7 point stencil. Once those are done we apply the operator in the cut-cells.
  const EBISBox& ebisbox = m_eblg.getEBISL()[a_dit];

  if (!ebisbox.isAllCovered()) {
    this->applyOpRegular(a_Lphi, a_phi, a_cellBox, a_dit, a_homogeneousPhysBC);
    this->applyOpIrregular(a_Lphi, a_phi, a_cellBox, a_dit, a_homogeneousPhysBC);
  }
}

void
EBHelmholtzOp::applyOpRegular(EBCellFAB&       a_Lphi,
                              EBCellFAB&       a_phi,
                              const Box&       a_cellBox,
                              const DataIndex& a_dit,
                              const bool       a_homogeneousPhysBC)
{
  CH_TIME("EBHelmholtzOp::applyOpRegular(EBCellFAB, EBCellFAB, Box, DataIndex, bool)");

  // TLDR: This is the regular kernel which computes L(phi) using the standard 5/7 point stencil. Since we want a simple kernel,
  //       we first fill the ghost cells on the domain boundary so that the flux through the boundary is consistent with the flux
  //       from the BC object.

  // Fill a_phi such that centered differences pushes in the domain flux.
  this->applyDomainFlux(a_phi, a_cellBox, a_dit, a_homogeneousPhysBC);

  BaseFab<Real>&       Lphi = a_Lphi.getSingleValuedFAB();
  const BaseFab<Real>& phi  = a_phi.getSingleValuedFAB();
  const BaseFab<Real>& aco  = (*m_Acoef)[a_dit].getSingleValuedFAB();

  // Need a handle to the regular b-coefficient which also exposes it in the kernel. This is
  // the best that I came up with.
  const BaseFab<Real>* bco[SpaceDim];
  for (int dir = 0; dir < SpaceDim; dir++) {
    bco[dir] = &(*m_Bcoef)[a_dit][dir].getSingleValuedFAB();
  }

  // This is the C++ kernel. It adds the diagonal and the Laplacian part.
  const Real factor = m_beta / (m_dx * m_dx);

  auto kernel = [&](const IntVect& iv) -> void {
    Lphi(iv, m_comp) = m_alpha * aco(iv, m_comp) * phi(iv, m_comp) +
                       factor * (+(*bco[0])(iv + BASISV(0), m_comp) * (phi(iv + BASISV(0), m_comp) - phi(iv, m_comp)) -
                                 (*bco[0])(iv, m_comp) * (phi(iv, m_comp) - phi(iv - BASISV(0), m_comp)) +
                                 (*bco[1])(iv + BASISV(1), m_comp) * (phi(iv + BASISV(1), m_comp) - phi(iv, m_comp)) -
                                 (*bco[1])(iv, m_comp) * (phi(iv, m_comp) - phi(iv - BASISV(1), m_comp))
#if CH_SPACEDIM == 3
                                 + (*bco[2])(iv + BASISV(2), m_comp) * (phi(iv + BASISV(2), m_comp) - phi(iv, m_comp)) -
                                 (*bco[2])(iv, m_comp) * (phi(iv, m_comp) - phi(iv - BASISV(2), m_comp))
#endif
                                );
  };

  // Launch the kernel.
  BoxLoops::loop(a_cellBox, kernel);
}

void
EBHelmholtzOp::applyDomainFlux(EBCellFAB&       a_phi,
                               const Box&       a_cellBox,
                               const DataIndex& a_dit,
                               const bool       a_homogeneousPhysBC)
{
  CH_TIME("EBHelmholtzOp::applyDomainFlux(EBCellFAB, Box, DataIndex, bool)");

  // TLDR: We compute the flux on the domain edges and store it in a cell-centered box. We then monkey with the ghost cells outside
  //       the domain so that centered differences on the edge cells inject said flux. This is a simple trick for enforcing the flux
  //       on the domain edges when we later compute the finite volume Laplacian.

  constexpr Real tol = 1.E-15;

  for (int dir = 0; dir < SpaceDim; dir++) {

    BaseFab<Real>& phiFAB = a_phi.getSingleValuedFAB();
    BaseFab<Real>& bco    = (*m_Bcoef)[a_dit][dir].getSingleValuedFAB();

    Box loBox;
    Box hiBox;
    int hasLo;
    int hasHi;
    EBArith::loHi(loBox, hasLo, hiBox, hasHi, m_eblg.getDomain(), a_cellBox, dir);

    if (hasLo == 1) {

      // Fill the domain flux. This might look weird, and we are actually putting the flux in a cell-centered data holder. By this, we implicitly
      // understand that the flux that is stored in the box is the flux that comes in through the lo side in the coordinate direction we are looking.
      FArrayBox faceFlux(loBox, m_nComp);
      m_domainBc->getFaceFlux(faceFlux, phiFAB, bco, dir, Side::Lo, a_dit, a_homogeneousPhysBC);

      // The EBArith loBox is cell-centered interior box abutting the domain side. We want the box immediately outside the domain.
      Box ghostBox = loBox;
      ghostBox.shift(dir, -1);

      // This kernel might look weird, but we have designed our BC classes in such a weird way -- they fill boundary fluxes but the fluxes
      // are stored in a cell-centered box abutting the domain. So, this is just like a "regular" kernel, with the exception of that pesky flux
      // which physically lives on the face but is computationally stored on the cell.
      auto kernel = [&](const IntVect& iv) -> void {
        const Real& B = bco(iv + BASISV(dir), m_comp);

        Real scaledFlux;
        if (std::abs(B) > tol) {
          scaledFlux = faceFlux(iv + BASISV(dir), m_comp) / B;
        }
        else {
          scaledFlux = 0.0;
        }
        phiFAB(iv, m_comp) = phiFAB(iv + BASISV(dir), m_comp) - scaledFlux * m_dx;
      };

      BoxLoops::loop(ghostBox, kernel);
    }

    if (hasHi == 1) {
      // Fill the domain flux. This might look weird, and we are actually putting the flux in a cell-centered data holder. By this, we implicitly
      // understand that the flux that is stored in the box is the flux that comes in through the hi side in the coordinate direction we are looking.
      FArrayBox faceFlux(hiBox, m_nComp);
      m_domainBc->getFaceFlux(faceFlux, phiFAB, bco, dir, Side::Hi, a_dit, a_homogeneousPhysBC);

      // The EBArith hiBox is cell-centered interior box abutting the domain side. We want the box immediately outside the domain.
      Box ghostBox = hiBox;
      ghostBox.shift(dir, 1);

      // This kernel might look weird, but we have designed our BC classes in such a weird way -- they fill boundary fluxes but the fluxes
      // are stored in a cell-centered box abutting the domain. So, this is just like a "regular" kernel, with the exception of that pesky flux
      // which physically lives on the face but is computationally stored on the cell.
      auto kernel = [&](const IntVect& iv) -> void {
        const Real& B = bco(iv - BASISV(dir), m_comp);

        Real scaledFlux;
        if (std::abs(B) > tol) {
          scaledFlux = faceFlux(iv - BASISV(dir), m_comp) / B;
        }
        else {
          scaledFlux = 0.0;
        }
        phiFAB(iv, m_comp) = phiFAB(iv - BASISV(dir), m_comp) + scaledFlux * m_dx;
      };

      BoxLoops::loop(ghostBox, kernel);
    }
  }
}

void
EBHelmholtzOp::fillDomainFlux(EBFluxFAB& a_flux, const EBCellFAB& a_phi, const Box& a_cellBox, const DataIndex& a_dit)
{
  CH_TIME("EBHelmholtzOp::fillDomainFlux(EBFluxFAB, EBCellFAB, Box, DataIndex)");

  // TLDR: We compute the flux on the domain edges and store it in a cell-centered box. We then monkey with the ghost cells outside
  //       the domain so that centered differences on the edge cells inject said flux. This is a simple trick for enforcing the flux
  //       on the domain edges when we later compute the finite volume Laplacian.

  constexpr Real tol = 1.E-15;

  for (int dir = 0; dir < SpaceDim; dir++) {

    BaseFab<Real>&       fluxReg = a_flux[dir].getSingleValuedFAB();
    const BaseFab<Real>& phiReg  = a_phi.getSingleValuedFAB();
    const BaseFab<Real>& bco     = (*m_Bcoef)[a_dit][dir].getSingleValuedFAB();

    Box loBox;
    Box hiBox;
    int hasLo;
    int hasHi;
    EBArith::loHi(loBox, hasLo, hiBox, hasHi, m_eblg.getDomain(), a_cellBox, dir);

    if (hasLo == 1) {

      // Fill the domain flux. This might look weird because we are putting the flux in a cell-centered box. By this, we implicitly
      // understand that the flux that is stored in the box is the flux that comes in through the lo side in the coordinate direction we are looking.
      FArrayBox faceFlux(loBox, m_nComp);
      m_domainBc->getFaceFlux(faceFlux, phiReg, bco, dir, Side::Lo, a_dit, false);

      // Copy flux over to the input data holder and multiply by beta.
      auto kernel = [&](const IntVect& iv) -> void {
        fluxReg(iv, m_comp) = m_beta * faceFlux(iv, m_comp);
      };

      BoxLoops::loop(loBox, kernel);
    }

    if (hasHi == 1) {
      // Fill the domain flux. This might look weird because we are putting the flux in a cell-centered box. By this, we implicitly
      // understand that the flux that is stored in the box is the flux that comes in through the lo side in the coordinate direction we are looking.
      FArrayBox faceFlux(hiBox, m_nComp);
      m_domainBc->getFaceFlux(faceFlux, phiReg, bco, dir, Side::Hi, a_dit, false);

      // Copy flux over to the input data holder and multiply by beta.
      auto kernel = [&](const IntVect& iv) -> void {
        fluxReg(iv + BASISV(dir), m_comp) = m_beta * faceFlux(iv, m_comp);
      };

      BoxLoops::loop(hiBox, kernel);
    }
  }
}

void
EBHelmholtzOp::applyOpIrregular(EBCellFAB&       a_Lphi,
                                const EBCellFAB& a_phi,
                                const Box&       a_cellBox,
                                const DataIndex& a_dit,
                                const bool       a_homogeneousPhysBC)
{
  CH_TIME("EBHelmholtzOp::applyOpIrregular(EBCellFAB, EBCellFAB, Box, DataIndex, bool)");

  // This routine computes L(phi) in a grid patch. This includes the cut-cell itself. Note that such cells can include cells that have both
  // an EB face and a domain face. This routine handles both.

  // Apply the operator in all cells where we needed an explicit stencil. Note that the operator stencils do NOT include stencils for fluxes
  // through the domain faces. That is handled below.
#if 0 // Original code that does not use AggStencil. Leaving this in place for backwards compatibility and debugging, in case this should ever break.
  VoFIterator& vofit = m_vofIterStenc[a_dit];
  for (vofit.reset(); vofit.ok(); ++vofit) {
    const VolIndex& vof = vofit();

    // Finite volume stencil representation of Int[B*grad(phi)]dA
    const VoFStencil& stenc = m_relaxStencils[a_dit](vof, m_comp);

    // kappa * alpha * aco (m_alphaDiagWeight holds kappa* aco)
    const Real& alphaDiag = m_alpha * m_alphaDiagWeight[a_dit](vof, m_comp);

    a_Lphi(vof, m_comp) = alphaDiag * a_phi(vof, m_comp);

    // Add the term corresponding to kappa*div(b*grad(phi)), excluding contributions from domain faces.
    for (int i = 0; i < stenc.size(); i++) {
      const VolIndex& ivof    = stenc.vof(i);
      const Real&     iweight = stenc.weight(i);

      a_Lphi(vof, m_comp) += m_beta * iweight * a_phi(ivof, m_comp); // Note that bco is a part of the stencil weight.
    }
  }
  m_ebBc->applyEBFlux(m_vofIterIrreg[a_dit], a_Lphi, a_phi, (*m_BcoefIrreg)[a_dit], a_dit, m_beta, a_homogeneousPhysBC);
#else // New code that uses AggStencil.
  constexpr bool incrementOnly = false;

  m_aggRelaxStencil[a_dit]->apply(a_Lphi, a_phi, m_alphaDiagWeight[a_dit], m_alpha, m_beta, m_comp, incrementOnly);
  m_ebBc->applyEBFlux(m_vofIterIrreg[a_dit], a_Lphi, a_phi, (*m_BcoefIrreg)[a_dit], a_dit, m_beta, a_homogeneousPhysBC);
#endif

  // Do irregular faces on domain sides. This was not included in the stencils above. m_domainBc should give the centroid-centered flux so we don't do interpolations here.
  for (int dir = 0; dir < SpaceDim; dir++) {

    // Iterators for high and low sides.
    VoFIterator& vofitLo = m_vofIterDomLo[dir][a_dit];
    VoFIterator& vofitHi = m_vofIterDomHi[dir][a_dit];

    const EBFaceFAB& Bcoef = (*m_Bcoef)[a_dit][dir];

    // Kernels for high and low sides.
    auto kernelLo = [&](const VolIndex& vof) -> void {
      const Real flux = m_domainBc->getFaceFlux(vof, a_phi, Bcoef, dir, Side::Lo, a_dit, a_homogeneousPhysBC);
      a_Lphi(vof, m_comp) -= flux * m_beta / m_dx;
    };

    auto kernelHi = [&](const VolIndex& vof) -> void {
      const Real flux = m_domainBc->getFaceFlux(vof, a_phi, Bcoef, dir, Side::Hi, a_dit, a_homogeneousPhysBC);
      a_Lphi(vof, m_comp) += flux * m_beta / m_dx;
    };

    // Execute kernels on lo and hi side
    BoxLoops::loop(vofitLo, kernelLo);
    BoxLoops::loop(vofitHi, kernelHi);
  }
}

void
EBHelmholtzOp::diagonalScale(LevelData<EBCellFAB>& a_rhs, bool a_kappaWeighted)
{
  CH_TIME("EBHelmholtzOp::diagonalScale(LD<EBCellFAB>, bool)");

  // Scale by volume fraction if asked.
  if (a_kappaWeighted) {
    EBLevelDataOps::kappaWeight(a_rhs);
  }

  // Scale by a-coefficient and alpha, too.
  for (DataIterator dit(a_rhs.dataIterator()); dit.ok(); ++dit) {
    a_rhs[dit()] *= (*m_Acoef)[dit()];
  }
}

void
EBHelmholtzOp::divideByIdentityCoef(LevelData<EBCellFAB>& a_rhs)
{
  CH_TIME("EBHelmholtzOp::divideByIdentityCoef(LD<EBCellFAB>)");

  for (DataIterator dit(a_rhs.disjointBoxLayout()); dit.ok(); ++dit) {
    a_rhs[dit()] /= (*m_Acoef)[dit()];
  }
}

void
EBHelmholtzOp::applyOpNoBoundary(LevelData<EBCellFAB>& a_Lphi, const LevelData<EBCellFAB>& a_phi)
{
  CH_TIME("EBHelmholtzOp::applyOpNoBounary(LD<EBCellFAB>, LD<EBCellFAB>)");

  m_doInterpCF = false;
  this->applyOp(a_Lphi, a_phi, true);
  m_doInterpCF = true;
}

void
EBHelmholtzOp::fillGrad(const LevelData<EBCellFAB>& a_phi)
{
  CH_TIME("EBHelmholtzOp::fillGrad(LD<EBCellFAB>, LD)");

  MayDay::Warning("EBHelmholtzOp::fillGrad - not implemented (yet)");
}

void
EBHelmholtzOp::getFlux(EBFluxFAB&                  a_flux,
                       const LevelData<EBCellFAB>& a_data,
                       const Box&                  a_grid,
                       const DataIndex&            a_dit,
                       Real                        a_scale)
{
  CH_TIME("EBHelmholtzOp::getFlux(EBFluxFAB, LD<EBCellFAB>, Box, DataIndex, Real)");

  MayDay::Warning("EBHelmholtzOp::getFlux - not implemented (yet)");
}

void
EBHelmholtzOp::homogeneousCFInterp(LevelData<EBCellFAB>& a_phi)
{
  CH_TIME("EBHelmholtzOp::homogeneousCFInterp(LD<EBCellFAB>)");

  // This routine interpolates ghost cells in a_phi, assuming that the coarse-level data is zero.

  if (m_hasCoar) {
    m_interpolator->coarseFineInterpH(a_phi, m_interval);
  }
}

void
EBHelmholtzOp::inhomogeneousCFInterp(LevelData<EBCellFAB>& a_phiFine, const LevelData<EBCellFAB>& a_phiCoar)
{
  CH_TIME("EBHelmholtzOp::inhomogeneousCFInterp(LD<EBCellFAB>, LD<EBCellFAB>)");

  // This routine interpolates ghost cells in a_phi. Responsibility for this is passed to the input interpolator.

  if (m_hasCoar) {
    m_interpolator->coarseFineInterp(a_phiFine, a_phiCoar, m_interval);
  }
}

void
EBHelmholtzOp::interpolateCF(LevelData<EBCellFAB>&       a_phiFine,
                             const LevelData<EBCellFAB>* a_phiCoar,
                             const bool                  a_homogeneousCFBC)
{
  CH_TIME("EBHelmholtzOp::interpolateCF(LD<EBCellFAB>, LD<EBCellFAB>, bool)");

  // This routine is just a wrapper for the two types of ghost cell interpolation. If a_phiCoar is nullptr then we are always
  // doing homogeneous interpolation. The flag a_homogeneousCFBC should specify this.

  if (m_hasCoar) {
    if (a_homogeneousCFBC) {
      this->homogeneousCFInterp(a_phiFine);
    }
    else {

      // I will call this an error because the user has probably unintentionally called for inhomogeneous interpolation when he/she meant
      // homogeneous interpolation.
      if (a_phiCoar == nullptr) {
        MayDay::Error("EBHelmholtzOp::interpolateCF -- calling inhomogeneousCFInterp with nullptr coarse is an error.");
      }

      this->inhomogeneousCFInterp(a_phiFine, *a_phiCoar);
    }
  }
}

void
EBHelmholtzOp::relax(LevelData<EBCellFAB>& a_correction, const LevelData<EBCellFAB>& a_residual, int a_iterations)
{
  CH_TIME("EBHelmholtzOp::relax(LD<EBCellFAB>, LD<EBCellFAB>, int)");

  // This function performs relaxation steps on the correction. User can switch between different kernels.

  switch (m_smoother) {
  case Smoother::NoRelax: {
    break;
  }
  case Smoother::PointJacobi: {
    this->relaxPointJacobi(a_correction, a_residual, a_iterations);

    break;
  }
  case Smoother::GauSaiRedBlack: {
    this->relaxGSRedBlack(a_correction, a_residual, a_iterations);

    break;
  }
  case Smoother::GauSaiMultiColor: {
    this->relaxGSMultiColor(a_correction, a_residual, a_iterations);

    break;
  }
  default: {
    MayDay::Error("EBHelmholtzOp::relax - bogus relaxation method requested");

    break;
  }
  }
}

void
EBHelmholtzOp::relaxPointJacobi(LevelData<EBCellFAB>&       a_correction,
                                const LevelData<EBCellFAB>& a_residual,
                                const int                   a_iterations)
{
  CH_TIME("EBHelmholtzOp::relaxPointJacobi(LD<EBCellFAB>, LD<EBCellFAB>, int)");

  // TLDR: This function performs point Jacobi relaxation in the form phi^(k+1) = phi^k - (res - L(phi))/|diag(L)|. Here, diag(L) is captured
  //       in m_relCoef. For performance integration with MFHelmholtzOp this routine passed that routine to a separate kernel call.

  LevelData<EBCellFAB> Lcorr;
  this->create(Lcorr, a_residual);

  const DisjointBoxLayout& dbl = m_eblg.getDBL();

  for (int iter = 0; iter < a_iterations; iter++) {
    if (m_doExchange) {
      a_correction.exchange();
    }

    this->homogeneousCFInterp(a_correction);

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      this->pointJacobiKernel(Lcorr[dit()], a_correction[dit()], a_residual[dit()], dbl[dit()], dit());
    }
  }
}

void
EBHelmholtzOp::pointJacobiKernel(EBCellFAB&       a_Lcorr,
                                 EBCellFAB&       a_correction,
                                 const EBCellFAB& a_residual,
                                 const Box&       a_cellBox,
                                 const DataIndex& a_dit)
{
  CH_TIME("EBHelmholtzOp::pointJacobiKernel(EBCellFAB, EBCellFAB, EBCellFAB, Box, DataIndex)");

  // TLDR: This just computes correction = correction - (res - L(phi))/diag(L). Note that we use under-relaxation (with a factor of 0.5) for stability.

  const EBISBox& ebisbox = m_eblg.getEBISL()[a_dit];

  if (!ebisbox.isAllCovered()) {
    this->applyOp(a_Lcorr, a_correction, a_cellBox, a_dit, true);

    a_Lcorr -= a_residual;
    a_Lcorr *= m_relCoef[a_dit];
    a_Lcorr *= 0.5;
    a_correction -= a_Lcorr;
  }
}

void
EBHelmholtzOp::relaxGSRedBlack(LevelData<EBCellFAB>&       a_correction,
                               const LevelData<EBCellFAB>& a_residual,
                               const int                   a_iterations)
{
  CH_TIME("EBHelmholtzOp::relaxGSRedBlack(LD<EBCellFAB>, LD<EBCellFAB>, int)");

  // TLDR: This function performs red-black Gauss-Seidel relaxation. As always, this occurs in the form phi^(k+1) = phi^k - (res - L(phi))/|diag(L)| but
  //       for a red-black update pattern:
  //
  //       For performance integration with MFHelmholtzOp this routine passed that routine to a separate kernel call.

  LevelData<EBCellFAB> Lcorr;
  this->create(Lcorr, a_residual);

  const DisjointBoxLayout& dbl = m_eblg.getDBL();

  for (int iter = 0; iter < a_iterations; iter++) {

    // First do "red" cells, then "black" cells. Note that ghost cell interpolation and exchanges are required between the colors.
    for (int redBlack = 0; redBlack <= 1; redBlack++) {
      if (m_doExchange) {
        a_correction.exchange();
      }

      this->homogeneousCFInterp(a_correction);

      for (DataIterator dit(dbl); dit.ok(); ++dit) {
        this->gauSaiRedBlackKernel(Lcorr[dit()], a_correction[dit()], a_residual[dit()], dbl[dit()], dit(), redBlack);
      }
    }
  }
}

void
EBHelmholtzOp::gauSaiRedBlackKernel(EBCellFAB&       a_Lcorr,
                                    EBCellFAB&       a_corr,
                                    const EBCellFAB& a_resid,
                                    const Box&       a_cellBox,
                                    const DataIndex& a_dit,
                                    const int&       a_redBlack)
{
  CH_TIME("EBHelmholtzOp::gauSaiRedBlackkernel(EBCellFAB, EBCellFAB, EBCellFAB, Box, DataIndex, int)");

  // This is the kernel for computing phi^(k+1) = phi^k - (res - L(phi))/|diag(L)| with a red-black pattern. Here, "red" cells are encoded by a_redBlack=0.

  const EBISBox&   ebisbox = m_eblg.getEBISL()[a_dit];
  const EBCellFAB& relCoef = m_relCoef[a_dit];

  if (!ebisbox.isAllCovered()) {
    this->applyOp(a_Lcorr, a_corr, a_cellBox, a_dit, true);

    BaseFab<Real>&       phiReg  = a_corr.getSingleValuedFAB();
    const BaseFab<Real>& LphiReg = a_Lcorr.getSingleValuedFAB();
    const BaseFab<Real>& rhsReg  = a_resid.getSingleValuedFAB();
    const BaseFab<Real>& relReg  = relCoef.getSingleValuedFAB();

    // Regular kernel. Several ways we can do this -- we can either check if the cell is red/black like we do here, which is the easiest. This should
    // not come at a performance cost, I think. An alternative is to compute an offset on the starting index on the innermost loop, like we used to do
    // with Fortran.
    auto regularKernel = [&](const IntVect& iv) -> void {
      const bool doThisCell = std::abs((iv.sum() + a_redBlack) % 2) == 0;

      if (doThisCell) {
        phiReg(iv, m_comp) += relReg(iv, m_comp) * (rhsReg(iv, m_comp) - LphiReg(iv, m_comp));
      }
    };

    // Irregular red-black kernel.
    auto irregularKernel = [&](const VolIndex& vof) -> void {
      const IntVect& iv = vof.gridIndex();

      const bool doThisCell = std::abs((iv.sum() + a_redBlack) % 2) == 0;

      if (doThisCell) {
        a_corr(vof, m_comp) += relCoef(vof, m_comp) * (a_resid(vof, m_comp) - a_Lcorr(vof, m_comp));
      }
    };

    // Launch the kernels over their respective domains.
    BoxLoops::loop(a_cellBox, regularKernel);
    BoxLoops::loop(m_vofIterMulti[a_dit], irregularKernel);
  }
}

void
EBHelmholtzOp::relaxGSMultiColor(LevelData<EBCellFAB>&       a_correction,
                                 const LevelData<EBCellFAB>& a_residual,
                                 const int                   a_iterations)
{
  CH_TIME("EBHelmholtzOp::relaxGSMultiColor(LD<EBCellFAB>, LD<EBCellFAB>, int)");

  // TLDR: This function performs multi-colored Gauss-Seidel relaxation. As always, this occurs in the form phi^(k+1) = phi^k - (res - L(phi))/|diag(L)| but
  //       using more colors than just red-black. The update pattern here cycles through quadrants/octants in 2D/3D. This is just like red-black except that
  //       we have four/eight colors in 2D/3D.
  //
  //       For performance integration with MFHelmholtzOp this routine passed that routine to a separate kernel call.

  LevelData<EBCellFAB> Lcorr;
  this->create(Lcorr, a_residual);

  const DisjointBoxLayout& dbl = m_eblg.getDBL();

  for (int iter = 0; iter < a_iterations; iter++) {
    for (int icolor = 0; icolor < m_colors.size(); icolor++) {
      if (m_doExchange) {
        a_correction.exchange();
      }

      this->homogeneousCFInterp(a_correction);

      for (DataIterator dit(dbl); dit.ok(); ++dit) {
        this->gauSaiMultiColorKernel(Lcorr[dit()],
                                     a_correction[dit()],
                                     a_residual[dit()],
                                     dbl[dit()],
                                     dit(),
                                     m_colors[icolor]);
      }
    }
  }
}

void
EBHelmholtzOp::gauSaiMultiColorKernel(EBCellFAB&       a_Lcorr,
                                      EBCellFAB&       a_corr,
                                      const EBCellFAB& a_resid,
                                      const Box&       a_cellBox,
                                      const DataIndex& a_dit,
                                      const IntVect&   a_color)
{
  CH_TIME("EBHelmholtzOp::gauSaiRedBlackkernel(EBCellFAB, EBCellFAB, EBCellFAB, Box, DataIndex, int)");

  // This is the kernel for computing phi^(k+1) = phi^k - (res - L(phi))/|diag(L)| with a "multi-colored" pattern. As with red-black, we follow a pattern, but
  // the pattern in this case uses more "colors".

  const EBISBox&   ebisbox = m_eblg.getEBISL()[a_dit];
  const EBCellFAB& relCoef = m_relCoef[a_dit];

  if (!ebisbox.isAllCovered()) {
    this->applyOp(a_Lcorr, a_corr, a_cellBox, a_dit, true);

    BaseFab<Real>&       phiReg  = a_corr.getSingleValuedFAB();
    const BaseFab<Real>& LphiReg = a_Lcorr.getSingleValuedFAB();
    const BaseFab<Real>& rhsReg  = a_resid.getSingleValuedFAB();
    const BaseFab<Real>& relReg  = relCoef.getSingleValuedFAB();

    // Regular cells (well, plus whatever is not multi-valued)
    IntVect loIV = a_cellBox.smallEnd();
    IntVect hiIV = a_cellBox.bigEnd();
    for (int dir = 0; dir < SpaceDim; dir++) {
      if (loIV[dir] % 2 != a_color[dir])
        loIV[dir]++;
    }

    if (loIV <= hiIV) {
      const Box colorBox(loIV, hiIV);

      // Regular kernel. This is just the point Jacobi -- the magic happens in the striding below.
      auto regularKernel = [&](const IntVect& iv) -> void {
        phiReg(iv, m_comp) += relReg(iv, m_comp) * (rhsReg(iv, m_comp) - LphiReg(iv, m_comp));
      };

      // Irregular kernel Does the same as above -- a stride of two.
      auto irregularKernel = [&](const VolIndex& vof) -> void {
        const IntVect& iv = vof.gridIndex();

        // Do stride check.
        bool doThisCell = true;
        for (int dir = 0; dir < SpaceDim; dir++) {
          if (iv[dir] % 2 != a_color[dir]) {
            doThisCell = false;
          }
        }

        if (doThisCell) {
          a_corr(vof, m_comp) += relCoef(vof, m_comp) * (a_resid(vof, m_comp) - a_Lcorr(vof, m_comp));
        }
      };

      // Launch the kernels.
      BoxLoops::loop(colorBox, regularKernel, 2 * IntVect::Unit);
      BoxLoops::loop(m_vofIterMulti[a_dit], irregularKernel);
    }
  }
}

void
EBHelmholtzOp::computeDiagWeight()
{
  CH_TIME("EBHelmholtzOp::computeDiagWeight()");

  // Compute the diagonal term in the kappa-weighted operator. We put the result into m_alphaDiagWeight = kappa * A;
  for (DataIterator dit(m_eblg.getDBL().dataIterator()); dit.ok(); ++dit) {
    const EBISBox& ebisbox = m_eblg.getEBISL()[dit()];

    VoFIterator& vofit = m_vofIterStenc[dit()];

    auto alphaKernel = [&](const VolIndex& vof) -> void {
      const Real volFrac = m_eblg.getEBISL()[dit()].volFrac(vof);
      const Real Aco     = (*m_Acoef)[dit()](vof, m_comp);

      m_alphaDiagWeight[dit()](vof, m_comp) = volFrac * Aco;
    };

    // Compute relaxation factor. Adjust the weight with domain boundary faces.
    auto betaKernel = [&](const VolIndex& vof) -> void {
      const IntVect iv = vof.gridIndex();

      VoFStencil& curStencil = m_relaxStencils[dit()](vof, m_comp);

      Real betaWeight = EBArith::getDiagWeight(curStencil, vof);
      for (int dir = 0; dir < SpaceDim; dir++) {
        for (SideIterator sit; sit.ok(); ++sit) {
          const Box sidebox = m_sideBox.at(std::make_pair(dir, sit()));

          if (sidebox.contains(iv)) {
            Real              weightedAreaFrac = 0.0;
            Vector<FaceIndex> faces            = ebisbox.getFaces(vof, dir, sit());
            for (auto& f : faces.stdVector()) {
              weightedAreaFrac += ebisbox.areaFrac(f) * (*m_Bcoef)[dit()][dir](f, m_comp) / (m_dx * m_dx);
            }
            betaWeight += -weightedAreaFrac;
          }
        }
      }

      m_betaDiagWeight[dit()](vof, m_comp) = betaWeight;
    };

    BoxLoops::loop(vofit, alphaKernel);
    BoxLoops::loop(vofit, betaKernel);
  }
}

void
EBHelmholtzOp::computeRelaxationCoefficient()
{
  CH_TIME("EBHelmholtzOp::computeRelaxationCoefficient()");

  // TLDR: Compute the relaxation coefficient in the operator. This is just the inverted diagonal of kappa*L(phi). It is inverted
  //       for performance reasons (because we divide by the diagonal in the relaxation steps).

  for (DataIterator dit(m_eblg.getDBL()); dit.ok(); ++dit) {
    const Box cellBox = m_eblg.getDBL()[dit()];

    // Add the diagonal term alpha * A
    m_relCoef[dit()].setVal(0.0);
    m_relCoef[dit()] += (*m_Acoef)[dit()];
    m_relCoef[dit()] *= m_alpha;

    // Add in the diagonal term for the variable-coefficient Laplacian
    BaseFab<Real>& regRel = m_relCoef[dit()].getSingleValuedFAB();
    for (int dir = 0; dir < SpaceDim; dir++) {

      // Regular kernel for adding -beta*(bcoef(loFace) + bcoef(hiFace))/dx^2 to the relaxation term.
      BaseFab<Real>& regBcoDir = (*m_Bcoef)[dit()][dir].getSingleValuedFAB();

      const Real factor        = m_beta / (m_dx * m_dx);
      auto       regularKernel = [&](const IntVect& iv) -> void {
        regRel(iv, m_comp) -= factor * (regBcoDir(iv + BASISV(dir), m_comp) + regBcoDir(iv, m_comp));
      };

      BoxLoops::loop(cellBox, regularKernel);
    }

    // At this point we have computed kappa*diag(L), but we need to invert it.
    auto inversionKernel = [&](const IntVect& iv) -> void {
      regRel(iv, m_comp) = 1. / regRel(iv, m_comp);
    };
    BoxLoops::loop(cellBox, inversionKernel);

    // Do the same for the irregular cells
    auto irregularKernel = [&](const VolIndex& vof) -> void {
      // m_alphaDiagWeight holds kappa * A
      const Real alphaWeight = m_alpha * m_alphaDiagWeight[dit()](vof, m_comp);

      // m_betaWeight holds the diagonal part of kappa*div(b*grad(phi)) in the cut-cells.
      const Real betaWeight = m_beta * m_betaDiagWeight[dit()](vof, m_comp);

      m_relCoef[dit()](vof, m_comp) = 1. / (alphaWeight + betaWeight);
    };

    BoxLoops::loop(m_vofIterStenc[dit()], irregularKernel);
  }
}

void
EBHelmholtzOp::makeAggStencil()
{
  CH_TIME("EBHelmholtzOp::makeAggStencil()");

  // These are proxies for when we compute L(phi). The number of ghost cells is important because
  // VCAggStencil will use these proxies for computing offsets in the data.
  LevelData<EBCellFAB> phiProxy(m_eblg.getDBL(), m_nComp, m_ghostPhi, EBCellFactory(m_eblg.getEBISL()));
  LevelData<EBCellFAB> rhsProxy(m_eblg.getDBL(), m_nComp, m_ghostRhs, EBCellFactory(m_eblg.getEBISL()));

  m_aggRelaxStencil.define(m_eblg.getDBL());

  for (DataIterator dit(m_eblg.getDBL()); dit.ok(); ++dit) {

    Vector<RefCountedPtr<BaseIndex>>   dstBaseIndex;
    Vector<RefCountedPtr<BaseStencil>> dstBaseStencil;

    VoFIterator& vofit = m_vofIterStenc[dit()];

    auto kernel = [&](const VolIndex& vof) -> void {
      const VoFStencil opSten = m_relaxStencils[dit()](vof, m_comp);

      dstBaseIndex.push_back(RefCountedPtr<BaseIndex>(new VolIndex(vof)));
      dstBaseStencil.push_back(RefCountedPtr<BaseStencil>(new VoFStencil(opSten)));
    };

    BoxLoops::loop(vofit, kernel);

    m_aggRelaxStencil[dit()] = RefCountedPtr<VCAggStencil>(new VCAggStencil(dstBaseIndex,
                                                                            dstBaseStencil,
                                                                            phiProxy[dit()],
                                                                            rhsProxy[dit()],
                                                                            m_relCoef[dit()],
                                                                            m_alphaDiagWeight[dit()],
                                                                            m_nComp));
  }
}

VoFStencil
EBHelmholtzOp::getFaceCenterFluxStencil(const FaceIndex& a_face, const DataIndex& a_dit) const
{
  CH_TIME("EBHelmholtzOp::getFaceCenterFluxStencil(FaceIndex, DataIndex)");

  // TLDR: This routine computes a regular finite difference stencil for getting a second-order accurate approximation to the face-centered flux (presuming
  //       that the data is cell-centered).
  VoFStencil fluxStencil;

  // BC handles the boundary fluxes.
  if (!a_face.isBoundary()) {
    fluxStencil.add(a_face.getVoF(Side::Hi), 1.0 / m_dx);
    fluxStencil.add(a_face.getVoF(Side::Lo), -1.0 / m_dx);
    fluxStencil *= (*m_Bcoef)[a_dit][a_face.direction()](a_face, m_comp);
  }

  return fluxStencil;
}

VoFStencil
EBHelmholtzOp::getFaceCentroidFluxStencil(const FaceIndex& a_face, const DataIndex& a_dit) const
{
  CH_TIME("EBHelmholtzOp::getFaceCentroidFluxStencil(FaceIndex, DataIndex)");

  // TLDR: This routine computes a second order accurate approximation to the flux on a cut-cell centroid. How this is done differs between discretizations. For
  //       cell-centered discretizations we get the face-centered fluxes and interpolate them to the face centroid. For centroid-based discretizations we have to
  //       compute the stencil directly, using least squares reconstruction. This is much more involved.

  VoFStencil fluxStencil;

  if (!a_face.isBoundary()) { // Domain BC classes handle domain faces.
    const EBISBox& ebisbox = m_eblg.getEBISL()[a_dit];

    const VolIndex loVoF = a_face.getVoF(Side::Lo);
    const VolIndex hiVoF = a_face.getVoF(Side::Hi);

    const bool irregFace = ebisbox.isIrregular(loVoF.gridIndex()) || ebisbox.isIrregular(hiVoF.gridIndex());

    // Centered differencing for regular faces.
    if (!irregFace) {
      fluxStencil = this->getFaceCenterFluxStencil(a_face, a_dit);
    }
    else {

      // Irregular face, but using cell-centered data. Get face-centered stencils and interpolate them to centroids.
      if (m_dataLocation == Location::Cell::Center) {
        const FaceStencil interpolationStencil =
          EBArith::getInterpStencil(a_face, IntVectSet(), ebisbox, m_eblg.getDomain());

        for (int i = 0; i < interpolationStencil.size(); i++) {
          const FaceIndex& iface   = interpolationStencil.face(i);
          const Real&      iweight = interpolationStencil.weight(i);

          // Get the face-centered stencil.
          VoFStencil fluxCenterStencil = this->getFaceCenterFluxStencil(iface, a_dit);

          fluxCenterStencil *= iweight;

          fluxStencil += fluxCenterStencil;
        }
      }

      // Irregular face, but with centroid-centered data. In this case we reconstruct the gradient using a least squares reconstruction of the solution.
      if (m_dataLocation == Location::Cell::Centroid) {
        constexpr int stencilWeight = 2;
        constexpr int stencilRadius = 2;
        constexpr int stencilOrder  = 2;

        const VoFStencil gradSten = LeastSquares::getGradSten(a_face,
                                                              LeastSquares::FaceLocation::Centroid,
                                                              LeastSquares::CellLocation::Centroid,
                                                              ebisbox,
                                                              m_dx,
                                                              stencilRadius,
                                                              stencilWeight,
                                                              stencilOrder,
                                                              IntVectSet());

        if (gradSten.size() > 0) {
          fluxStencil = LeastSquares::projectGradSten(gradSten, BASISREALV(a_face.direction()));
          fluxStencil *= (*m_Bcoef)[a_dit][a_face.direction()](a_face, m_comp);
        }
        else {
          MayDay::Warning(
            "EBHelmholtz::getFaceCentroidFluxStencil -- could not find centroid stencil. Maybe your multigrid runs too deep?");
        }
      }
    }
  }

  return fluxStencil;
}

void
EBHelmholtzOp::getFaceCentroidFlux(EBFaceFAB&       a_fluxCentroid,
                                   const EBCellFAB& a_phi,
                                   const Box&       a_cellBox,
                                   const DataIndex& a_dit,
                                   const int        a_dir)
{
  CH_TIME("EBHelmholtzOp::getFaceCentroidFlux(EBFaceFAB, EBCellFAB, Box, DataIndex, int)");

  // TLDR: This routine computes the face centroid flux on all regular and irregular faces in the input box.

  // First do the regular faces and then the irregular faces. Not that when a centroid discretization is used, some faces can be irregular but nonetheless require
  // an explicit stencil.
  this->computeFaceCenteredFlux(a_fluxCentroid, a_phi, a_cellBox, a_dit, a_dir);
  this->computeFaceCentroidFlux(a_fluxCentroid, a_phi, a_cellBox, a_dit, a_dir);

  a_fluxCentroid *= m_beta;
}

void
EBHelmholtzOp::computeFaceCenteredFlux(EBFaceFAB&       a_fluxCenter,
                                       const EBCellFAB& a_phi,
                                       const Box&       a_cellBox,
                                       const DataIndex& a_dit,
                                       const int        a_dir)
{
  CH_TIME("EBHelmholtzOp::computeFaceCenteredFlux(EBFaceFAB, EBCellFAB, Box, DataIndex, int)");

  // Since a_fluxCenter is used as a basis for interpolation near the cut-cells, we need to fill fluxes on faces in direction that extends one cell
  // outside of the cellbox in directions that are tangential to a_dir. We do not want to include domain faces here since that
  // is taken care of elsewhere.

  // Discard domain faces
  Box compCellBox = a_cellBox;
  compCellBox.grow(a_dir, 1);
  compCellBox &= m_eblg.getDomain();
  compCellBox.grow(a_dir, -1);

  BaseFab<Real>&       regFlux = a_fluxCenter.getSingleValuedFAB();
  const BaseFab<Real>& regPhi  = a_phi.getSingleValuedFAB();
  const BaseFab<Real>& regBco  = (*m_Bcoef)[a_dit][a_dir].getSingleValuedFAB();

  // This kernel does centered differencing, setting the face flux to flux = B*(phi(high) - phi(lo))/dx. Recall that regFlux
  // is face centered but phi lives on in the cell.
  const Real inverseDx     = 1. / m_dx;
  auto       regularKernel = [&](const IntVect& iv) -> void {
    regFlux(iv, m_comp) = regBco(iv, m_comp) * inverseDx * (regPhi(iv, m_comp) - regPhi(iv - BASISV(a_dir), m_comp));
  };

  // Kernel region. All interior faces in compCellBox
  const Box faceBox = surroundingNodes(compCellBox, a_dir);

  // Launch kernel.
  BoxLoops::loop(faceBox, regularKernel);
}

void
EBHelmholtzOp::computeFaceCentroidFlux(EBFaceFAB&       a_flux,
                                       const EBCellFAB& a_phi,
                                       const Box&       a_cellBox,
                                       const DataIndex& a_dit,
                                       const int        a_dir)
{
  CH_TIME("EBHelmholtzOp::computeFaceCentroidFlux(EBFaceFAB, EBCellFAB, Box, DataIndex, int)");

  // This routine computes the face centroid fluxes using precomputed stencils (in defineStencils). This is needed because
  // cut-cell faces require more than centered differencing.

  const BaseIFFAB<VoFStencil>& fluxStencils = m_centroidFluxStencil[a_dir][a_dit];
  const EBGraph&               ebgraph      = fluxStencils.getEBGraph();
  IntVectSet                   ivs          = fluxStencils.getIVS();

  ivs &= a_cellBox;

  FaceIterator faceIt(ivs, ebgraph, a_dir, FaceStop::SurroundingNoBoundary);

  auto kernel = [&](const FaceIndex& face) -> void {
    const VoFStencil& sten = fluxStencils(face, m_comp);

    a_flux(face, m_comp) = 0.0;
    for (int i = 0; i < sten.size(); i++) {
      const VolIndex& ivof    = sten.vof(i);
      const Real&     iweight = sten.weight(i);

      a_flux(face, m_comp) += iweight * a_phi(ivof, m_comp);
    }
  };

  BoxLoops::loop(faceIt, kernel);
}

void
EBHelmholtzOp::computeFlux(const LevelData<EBCellFAB>& a_phi)
{
  CH_TIME("EBHelmholtzOp::computeFlux(LD<EBCellFAB>)");

  const DisjointBoxLayout& dbl = m_eblg.getDBL();

  for (DataIterator dit(dbl); dit.ok(); ++dit) {
    for (int dir = 0; dir < SpaceDim; dir++) {
      this->getFaceCentroidFlux(m_flux[dit()][dir], a_phi[dit()], dbl[dit()], dit(), dir);
    }
  }
}

void
EBHelmholtzOp::incrementFRCoar(const LevelData<EBCellFAB>& a_phi)
{
  CH_TIME("EBHelmholtzOp::incrementFRCoar(LD<EBCellFAB>)");

  CH_assert(m_hasFine);

  // TLDR: We want to compute the flux on this level so we can later set it to the sum of the fluxes through the fine-grid faces.

  const Real scale = 1.0;

  if (m_doExchange) {
    LevelData<EBCellFAB>& phiCoar = (LevelData<EBCellFAB>&)a_phi;
    phiCoar.exchange();
  }

  for (DataIterator dit(m_eblg.getDBL()); dit.ok(); ++dit) {
    for (int dir = 0; dir < SpaceDim; dir++) {
      EBFaceFAB& flux = m_flux[dit()][dir];

      const Box cellBox = m_eblg.getDBL()[dit()];
      for (SideIterator sit; sit.ok(); ++sit) {

        // Computes fluxes for all faces oriented in +/- dir, but which are not boundary faces. Must do this for all faces because
        // the CF interface does not align with a coarse box.
        this->getFaceCentroidFlux(flux, a_phi[dit()], cellBox, dit(), dir);

        // Increment flux register, recall that beta and the b-coefficient are included in getFaceCentroidFlux.
        m_fluxReg->incrementCoarseBoth(flux, scale, dit(), m_interval, dir, sit());
      }
    }
  }
}

void
EBHelmholtzOp::incrementFRFine(const LevelData<EBCellFAB>&       a_phiFine,
                               const LevelData<EBCellFAB>&       a_phi,
                               AMRLevelOp<LevelData<EBCellFAB>>& a_finerOp)
{
  CH_TIME("EBHelmholtzOp::incrementFRFine(LD<EBCellFAB>, LD<EBCellFAB>, AMRLevelOp<LevelData<EBCellFAB>)");

  CH_assert(m_hasFine);

  // TLDR: We want to compute the fluxes on the fine grid so that we can later use these fluxes in the reflux step. This just means
  //       that we must compute the fine-grid fluxes along the refinement bounary and add those fluxes to the flux register.

  // Doing the nasty here -- there's no guarantee that a_phiFine has the ghost cells it needs.
  LevelData<EBCellFAB>& phiFine = (LevelData<EBCellFAB>&)a_phiFine;
  phiFine.exchange();

  EBHelmholtzOp& finerOp = (EBHelmholtzOp&)(a_finerOp);
  finerOp.inhomogeneousCFInterp(phiFine, a_phi);

  LevelData<EBFluxFAB>& fineFlux = finerOp.getFlux();
  const Real            scale    = 1.0;

  // Compute
  for (DataIterator dit(m_eblgFine.getDBL()); dit.ok(); ++dit) {
    for (int dir = 0; dir < SpaceDim; dir++) {
      EBFaceFAB& flux = fineFlux[dit()][dir];

      const Box cellBox = m_eblgFine.getDBL()[dit()];
      for (SideIterator sit; sit.ok(); ++sit) {

        // The CF interface always ends at the boundary of a fine-level grid box, so we can run the computation for
        // the subset of cells that align with the CF.
        Box stripBox = adjCellBox(cellBox, dir, sit(), 1); // This box is outside the cellBox patch
        stripBox.shift(dir, -sign(sit()));                 // This box is inside the cellBox patch

        // Computes fluxes for all faces oriented in +/- dir, but which are not boundary faces.
        finerOp.getFaceCentroidFlux(flux, phiFine[dit()], stripBox, dit(), dir);

        // Increment flux register, recall that beta and the b-coefficient are included in getFaceCentroidFlux.
        m_fluxReg->incrementFineBoth(flux, scale, dit(), m_interval, dir, sit());
      }
    }
  }
}

void
EBHelmholtzOp::reflux(LevelData<EBCellFAB>&             a_Lphi,
                      const LevelData<EBCellFAB>&       a_phiFine,
                      const LevelData<EBCellFAB>&       a_phi,
                      AMRLevelOp<LevelData<EBCellFAB>>& a_finerOp)
{
  CH_TIME("EBHelmholtzOp::reflux(LD<EBCellFAB>, LD<EBCellFAB>, LD<EBCellFAB>, AMRLevelOp<LevelData<EBCellFAB>)");

  // This routine computes the fluxes on the coarse and fine-side of the boundary and does a refluxing operation where
  // we subtract the contribution from the coarse grid fluxes in a_Lphi and add in the contribution from the fine grid fluxes.

  m_fluxReg->setToZero();

  this->incrementFRCoar(a_phi);
  this->incrementFRFine(a_phiFine, a_phi, a_finerOp);

  m_fluxReg->reflux(a_Lphi, m_interval, 1. / m_dx);
}

void
EBHelmholtzOp::coarsen(LevelData<EBCellFAB>& a_phi, const LevelData<EBCellFAB>& a_phiFine)
{
  CH_TIME("EBHelmholtzOp::coarsen(LD<EBCellFAB>, LD<EBCellFAB>)");

  m_coarAve->averageData(a_phi, a_phiFine, m_interval, Average::Conservative);
}

void
EBHelmholtzOp::coarsen(LevelData<EBFluxFAB>& a_flux, const LevelData<EBFluxFAB>& a_fineFlux)
{
  m_coarAve->averageData(a_flux, a_fineFlux, m_interval, Average::Conservative);
}

#include <CD_NamespaceFooter.H>
