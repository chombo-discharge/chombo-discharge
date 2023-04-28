/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBGradient.cpp
  @brief  Implementation of CD_EBGradient.H
  @author Robert Marskar
*/

// Chombo includes
#include <ParmParse.H>
#include <EBArith.H>
#include <EBCellFactory.H>
#include <BaseIVFactory.H>
#include <NeighborIterator.H>
#include <CH_Timer.H>

// Our includes
#include <CD_Timer.H>
#include <CD_EBGradient.H>
#include <CD_ParallelOps.H>
#include <CD_LeastSquares.H>
#include <CD_VofUtils.H>
#include <CD_LoadBalancing.H>
#include <CD_BoxLoops.H>
#include <CD_NamespaceHeader.H>

constexpr int EBGradient::m_comp;
constexpr int EBGradient::m_nComp;

EBGradient::EBGradient(const EBLevelGrid& a_eblg,
                       const EBLevelGrid& a_eblgFine,
                       const EBLevelGrid& a_eblgFiCo,
                       const bool         a_hasFine,
                       const Real         a_dx,
                       const int          a_refRat,
                       const int          a_order,
                       const int          a_weighting,
                       const IntVect      a_ghostVector) noexcept
{
  CH_TIME("EBGradient::EBGradient");

  this->define(a_eblg, a_eblgFine, a_eblgFiCo, a_hasFine, a_dx, a_refRat, a_order, a_weighting, a_ghostVector);
}

void
EBGradient::define(const EBLevelGrid& a_eblg,
                   const EBLevelGrid& a_eblgFine,
                   const EBLevelGrid& a_eblgFiCo,
                   const bool         a_hasFine,
                   const Real         a_dx,
                   const int          a_refRat,
                   const int          a_order,
                   const int          a_weighting,
                   const IntVect      a_ghostVector) noexcept
{
  CH_TIME("EBGradient::define");

  CH_assert(a_order > 0);
  CH_assert(a_weighting >= 0);

  m_eblg         = a_eblg;
  m_eblgFine     = a_eblgFine;
  m_eblgFiCo     = a_eblgFiCo;
  m_dataLocation = Location::Cell::Center;
  m_dx           = a_dx;
  m_refRat       = a_refRat;
  m_order        = a_order;
  m_weighting    = a_weighting;
  m_ghostVector  = a_ghostVector;
  m_hasFine      = a_hasFine;
  m_hasEBCF      = false;

  if (m_hasFine) {
    m_dxFine = m_dx / m_refRat;
  }
  else {
    m_dxFine = 1;
  }

  // Define the level stencils. These are regular finite-difference stencils.
  this->defineLevelStencils();

  bool      makeEBCF = false;
  ParmParse pp("EBGradient");
  pp.query("ebcf", makeEBCF);

  // If there's finer level there might also be an EBCF crossing. Those can be tricky to deal with because
  // the stencil on this level (i.e., the coarse level) might reach underneath the finer level and obtain bogus data.
  // This loop defines a few masks for figuring out which cells lie on the CF region, and which cells can't be used
  // for finite differencing. Once those masks have been defined, we iterate through them and see if there are cells
  // where we can't find good stencils. If we do find such cells, we trigger m_hasEBCF which will define
  // a new set of grids where we run a least squares procedure. THOSE stencils are "dual-level" stencils
  // that ONLY reach into valid regions (i.e., not covered by a finer grid level).
  if (m_hasFine && makeEBCF) {

    // Masks with transient lifetimes.
    LevelData<FArrayBox> coarMaskCF;
    LevelData<FArrayBox> coarMaskInvalid;

    // Define masks on the input grids.
    this->defineMasks(coarMaskCF, coarMaskInvalid);

    // Define iterators for the input grid AND define the simplified buffer grids.
    this->defineIteratorsEBCF(coarMaskCF, coarMaskInvalid);

    if (m_hasEBCF) {
      this->defineStencilsEBCF(coarMaskInvalid);
    }
  }

  // Define optimized stencils.
  this->makeAggStencils();

  m_isDefined = true;
}

EBGradient::~EBGradient() noexcept { CH_TIME("EBGradient::~EBGradient"); }

void
EBGradient::computeLevelGradient(LevelData<EBCellFAB>& a_gradient, const LevelData<EBCellFAB>& a_phi) const noexcept
{
  CH_TIMERS("EBGradient::computeLevelGradient");
  CH_TIMER("EBGradient::regular_cells", t1);
  CH_TIMER("EBGradient::irregular_cells", t2);
  CH_TIMER("EBGradient::set_covered_cells", t3);

  CH_assert(m_isDefined);
  CH_assert(a_gradient.nComp() == SpaceDim);
  CH_assert(a_phi.nComp() == 1);
  CH_assert(a_phi.ghostVect() == m_ghostVector);
  CH_assert(a_gradient.ghostVect() == m_ghostVector);

  // TLDR: This routine computes the level gradient, i.e. using finite difference stencils isolated to this level.
  const DisjointBoxLayout& dbl   = m_eblg.getDBL();
  const EBISLayout&        ebisl = m_eblg.getEBISL();

  for (DataIterator dit(dbl); dit.ok(); ++dit) {
    EBCellFAB&       grad    = a_gradient[dit()];
    const EBCellFAB& phi     = a_phi[dit()];
    const EBISBox&   ebisBox = ebisl[dit()];
    const Box        cellBox = dbl[dit()];

    if (!ebisBox.isAllCovered()) {

      // Regular kernel -- this is the kernel used in regular cells. It just uses 2nd order centered differencing
      // for the first derivative.
      BaseFab<Real>&       gradFAB = grad.getSingleValuedFAB();
      const BaseFab<Real>& phiFAB  = phi.getSingleValuedFAB();
      const Real           idx     = 1. / (2.0 * m_dx);

      auto regularKernel = [&](const IntVect& iv) -> void {
        for (int dir = 0; dir < SpaceDim; dir++) {
          gradFAB(iv, dir) = idx * (phiFAB(iv + BASISV(dir), m_comp) - phiFAB(iv - BASISV(dir), m_comp));
        }
      };

      // Irregular kernel -- this is the kernel used in the cut-cells. It simply applies the finite difference based stencil.
      auto irregularKernel = [&, this](const VolIndex& vof) -> void {
        for (int dir = 0; dir < SpaceDim; dir++) {
          grad(vof, dir) = 0.0;
        }

        // Apply stencil. Note that the stencil "variable" is the gradient component (i.e., direction)
        auto& sten = m_levelStencils[dit()](vof, m_comp);
        for (int i = 0; i < sten.size(); i++) {
          const VolIndex& ivof    = sten.vof(i);
          const Real&     iweight = sten.weight(i);
          const Real&     ivar    = sten.variable(i);

          grad(vof, ivar) += iweight * phi(ivof, m_comp);
        }
      };

      // Now apply the kernels using our nifty BoxLoops.
      CH_START(t1);
      BoxLoops::loop(dbl[dit()], regularKernel);
      CH_STOP(t1);

      // Irregular cells done using AggStencil. Which is faster.
      CH_START(t2);
      for (int dir = 0; dir < SpaceDim; dir++) {
        EBCellFAB alias(Interval(dir, dir), grad);
        m_aggLevelStencils[dir][dit()]->apply(alias, phi, 0, 0, 1, false);
      }
      CH_STOP(t2);
    }

    // Covered data is always bogus.
    CH_START(t3);
    for (int dir = 0; dir < SpaceDim; dir++) {
      a_gradient[dit()].setCoveredCellVal(0.0, dir);
    }
    CH_STOP(t3);
  }
}

void
EBGradient::computeNormalDerivative(LevelData<EBFluxFAB>& a_gradient, const LevelData<EBCellFAB>& a_phi) const noexcept
{
  CH_TIMERS("EBGradient::computeNormalDerivative");
  CH_TIMER("EBGradient::regular_faces", t1);
  CH_TIMER("EBGradient::irregular_faces", t2);

  CH_assert(a_gradient.nComp() == SpaceDim);
  CH_assert(a_phi.nComp() == 1);

  // TLDR: This routine computes the level gradient, i.e. using finite difference stencils isolated to this level. It only does this for interior
  // faces, and it only computes the component of the gradient that is normal to the face.
  const DisjointBoxLayout& dbl    = m_eblg.getDBL();
  const EBISLayout&        ebisl  = m_eblg.getEBISL();
  const ProblemDomain&     domain = m_eblg.getDomain();
  const Real               idx    = 1. / m_dx; // Needed for the kernels.

  for (DataIterator dit(dbl); dit.ok(); ++dit) {
    const Box        cellBox = dbl[dit()];
    const EBCellFAB& phi     = a_phi[dit()];
    const EBISBox&   ebisBox = ebisl[dit()];
    const EBGraph&   ebgraph = ebisBox.getEBGraph();

    for (int dir = 0; dir < SpaceDim; dir++) {
      EBFaceFAB&           grad        = a_gradient[dit()][dir];
      BaseFab<Real>&       regGradient = grad.getSingleValuedFAB();
      const BaseFab<Real>& regPhi      = phi.getSingleValuedFAB();

      // Get all cells that do abut the domain edge.
      Box interiorCells = cellBox;
      interiorCells.grow(dir, 1);
      interiorCells &= domain;
      interiorCells.grow(dir, -1);

      // Turn the interiorCells Box to a face-centered box in the dir-direction. When we do this we get all faces that are oriented along the direction
      // dir, but which are not boundary faces. Also define an iterator for going through the subset of those cells that shared a cell edge/face with a cut-cell.
      const Box         interiorFaces = surroundingNodes(interiorCells, dir);
      const IntVectSet& irregIVS      = ebisBox.getIrregIVS(interiorCells);

      FaceIterator faceit(irregIVS, ebgraph, dir, FaceStop::SurroundingNoBoundary);

      // C++ kernel for regular grid faces. Note that we iterate over the IntVects in the face-centered boxes, so
      // the cell on the high side of the face has index iv and on the low side it has index iv - BASISV(dir);
      auto regularFaceDerivative = [&](const IntVect& iv) -> void {
        regGradient(iv, dir) = idx * regPhi(iv, m_comp) - regPhi(iv - BASISV(dir), m_comp);
      };

      // Cut-cell version of the above.
      auto irregularFaceDerivative = [&](const FaceIndex& face) -> void {
        const VolIndex& vofHi = face.getVoF(Side::Hi);
        const VolIndex& vofLo = face.getVoF(Side::Lo);

        grad(face, dir) = idx * (phi(vofHi, m_comp) - phi(vofLo, m_comp));
      };

      // Launch our C++ kernels.
      CH_START(t1);
      BoxLoops::loop(interiorFaces, regularFaceDerivative);
      CH_STOP(t1);
      CH_START(t2);
      BoxLoops::loop(faceit, irregularFaceDerivative);
      CH_STOP(t2);
    }
  }
}

void
EBGradient::computeAMRGradient(LevelData<EBCellFAB>&       a_gradient,
                               const LevelData<EBCellFAB>& a_phi,
                               const LevelData<EBCellFAB>& a_phiFine) const noexcept
{
  CH_TIMERS("EBGradient::computeAMRGradient");
  CH_TIMER("EBGradient::buffer_def", t1);
  CH_TIMER("EBGradient::copy_and_exchange", t2);
  CH_TIMER("EBGradient::ebcf_calculate", t3);

  CH_assert(m_isDefined);
  CH_assert(a_gradient.nComp() == SpaceDim);
  CH_assert(a_phi.nComp() == 1);
  CH_assert(a_phiFine.nComp() == 1);

  // TLDR: This routine computes the two-level gradient. It first computes the level gradient and then iterates through
  //       all cells that require a "two-level" view of the gradient.
  this->computeLevelGradient(a_gradient, a_phi);

  // Do corrections near the EBCF.
  if (m_hasFine && m_hasEBCF) {
    const DisjointBoxLayout& dbl     = m_eblg.getDBL();
    const DisjointBoxLayout& dblFiCo = m_eblgFiCo.getDBL();

    const EBISLayout& ebisl     = m_eblg.getEBISL();
    const EBISLayout& ebislFiCo = m_eblgFiCo.getEBISL();

    CH_START(t1)
    LevelData<EBCellFAB> phiFiCo(dblFiCo, 1, m_ghostVector, EBCellFactory(ebislFiCo));
    CH_STOP(t1);

    CH_START(t2);
    a_phiFine.copyTo(phiFiCo);
    CH_STOP(t2);

    CH_START(t3);
    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      EBCellFAB&       gradient = a_gradient[dit()];
      const EBCellFAB& phi      = a_phi[dit()];
      const EBCellFAB& phiFine  = phiFiCo[dit()];

      const BaseIVFAB<VoFStencil>& stencilsCoar = m_ebcfStencilsCoar[dit()];
      const BaseIVFAB<VoFStencil>& stencilsFine = m_ebcfStencilsFine[dit()];

      auto kernel = [&](const VolIndex& vof) -> void {
        const VoFStencil& stencilCoar = stencilsCoar(vof, m_comp);
        const VoFStencil& stencilFine = stencilsFine(vof, m_comp);

        for (int dir = 0; dir < SpaceDim; dir++) {
          gradient(vof, dir) = 0.0;
        }

        for (int i = 0; i < stencilCoar.size(); i++) {
          const VolIndex& coarVoF    = stencilCoar.vof(i);
          const Real&     coarWeight = stencilCoar.weight(i);
          const int&      coarVar    = stencilCoar.variable(i);

          gradient(vof, coarVar) += coarWeight * phi(coarVoF, m_comp);
        }

        for (int i = 0; i < stencilFine.size(); i++) {
          const VolIndex& fineVoF    = stencilFine.vof(i);
          const Real&     fineWeight = stencilFine.weight(i);
          const int&      fineVar    = stencilFine.variable(i);

          gradient(vof, fineVar) += fineWeight * phiFine(fineVoF, m_comp);
        }
      };

      BoxLoops::loop(m_ebcfIterator[dit()], kernel);
    }
    CH_STOP(t3);
  }
}

void
EBGradient::defineLevelStencils() noexcept
{
  CH_TIMERS("EBGradient::defineLevelStencils");
  CH_TIMER("EBGradient::define_dbl", t1);
  CH_TIMER("EBGradient::define_patch", t2);
  CH_TIMER("EBGradient::irreg_loop", t3);

  // TLDR: This just defines finite-difference stencils on this level. We need explicit stencils at the domain boundary AND at
  //       the EB. There are other cases, too, but those are handled by defineStencilsEBCF.

  const DisjointBoxLayout& dbl    = m_eblg.getDBL();
  const EBISLayout&        ebisl  = m_eblg.getEBISL();
  const ProblemDomain&     domain = m_eblg.getDomain();

  CH_START(t1);
  m_levelStencils.define(dbl);
  m_levelIterator.define(dbl);
  CH_STOP(t1);

  for (DataIterator dit(dbl); dit.ok(); ++dit) {
    CH_START(t2);
    const Box      cellBox = dbl[dit()];
    const EBISBox& ebisbox = ebisl[dit()];
    const EBGraph& ebgraph = ebisbox.getEBGraph();

    // Find cells where we need explicit stencils.
    IntVectSet bndryIVS = ebisbox.getIrregIVS(cellBox);
    for (int dir = 0; dir < SpaceDim; dir++) {
      Box loBox;
      Box hiBox;
      int hasLo;
      int hasHi;

      EBArith::loHi(loBox, hasLo, hiBox, hasHi, domain, cellBox, dir);

      if (hasLo) {
        bndryIVS |= IntVectSet(loBox);
      }
      if (hasHi) {
        bndryIVS |= IntVectSet(hiBox);
      }
    }

    VoFIterator&           bndryIterator = m_levelIterator[dit()];
    BaseIVFAB<VoFStencil>& bndryStencils = m_levelStencils[dit()];

    bndryIterator.define(bndryIVS, ebgraph);
    bndryStencils.define(bndryIVS, ebgraph, m_nComp);
    CH_STOP(t2);

    auto kernel = [&](const VolIndex& vof) -> void {
      VoFStencil& stencil = bndryStencils(vof, m_comp);

      stencil.clear();

      // Get derivative in each direction. Note that EBArith might decide to drop order.
      for (int dir = 0; dir < SpaceDim; dir++) {
        VoFStencil derivDirStencil;
        EBArith::getFirstDerivStencilWidthOne(derivDirStencil, vof, ebisbox, dir, m_dx, nullptr, dir);
        stencil += derivDirStencil;
      }
    };

    CH_START(t3);
    BoxLoops::loop(bndryIterator, kernel);
    CH_STOP(t3);
  }
}

void
EBGradient::defineMasks(LevelData<FArrayBox>& a_coarMaskCF, LevelData<FArrayBox>& a_coarMaskInvalid) noexcept
{
  CH_TIME("EBGradient::defineMasks");

  CH_assert(m_hasFine);

  // TLDR: This routine defines the masks in a_coarMaskCF and a_coarMaskInvalid. These masks live on the coarse level
  //       and they hold a value of 1 in coarse cells that abut the refinement boundary (a_coarMaskCF) and a value of 1
  //       in cells that are covered by a finer level (a_coarMaskInvalid). These are just constructed by creating some
  //       buffer storage on the (coarsened) fine grid with a value of 1 in the appropriate regions, and then copying
  //       the result to the coarse grid.

  constexpr Real zero   = 0.0;
  constexpr Real one    = 1.0;
  const Interval interv = Interval(m_comp, m_comp);

  // Handle to computational grids.
  const DisjointBoxLayout& dblCoar = m_eblg.getDBL();
  const DisjointBoxLayout& dblFine = m_eblgFine.getDBL();

  const ProblemDomain& domainCoar = m_eblg.getDomain();

  // Create some temporary storage. We coarsen the fine grid and create a Level<FArrayBox> object where each
  // Box is grown by one. We do this because we want to put some data in the ghost cells outside the valid region
  // of the fine level.
  DisjointBoxLayout dblCoFi;
  coarsen(dblCoFi, dblFine, m_refRat);

  LevelData<FArrayBox> coFiMaskCF(dblCoFi, m_nComp, IntVect::Unit);
  LevelData<FArrayBox> coFiMaskInvalid(dblCoFi, m_nComp, IntVect::Zero);

  // Build the mask regions.
  for (DataIterator dit(dblCoFi); dit.ok(); ++dit) {
    FArrayBox& maskCF      = coFiMaskCF[dit()];
    FArrayBox& maskInvalid = coFiMaskInvalid[dit()];
    const Box  cellBox     = dblCoFi[dit()];
    const Box  ghostedBox  = maskCF.box();

    // Set the "coarse-fine" mask. It is set to 1 in cells that are not valid cells. We do this by setting the FArrayBox data to one
    // in the entire patch, and then we iterate through the parts of the patch that overlap with valid grid boxes. That part of the
    // data is set to zero. We end up with a bunch of grid patches which hold a value of 1 in ghost cells that do not overlap with
    // the valid region of the grid.
    maskCF.setVal(one);
    maskCF.setVal(zero, cellBox, m_comp, m_nComp);

    NeighborIterator nit(dblCoFi);
    for (nit.begin(dit()); nit.ok(); ++nit) {
      const Box overlapBox = dblCoFi[nit()] & ghostedBox;

      if (!overlapBox.isEmpty()) {
        maskCF.setVal(zero, overlapBox, m_comp, m_nComp);
      }
    }

    // Set the invalid cell mask on the coarsened fine grids. We set it to 1 everywhere on the fine grid and 0 elsewhere. When we add this to the coarse
    // level we will add 1 into every cell that is covered by a finer level.
    maskInvalid.setVal(one, cellBox, m_comp, m_nComp);
  }

  // Define the input masks. After this, coarMaskCF will have a value of 1 in all cells that abut the fine level. Likewise, coarMaskInvalid
  // will have a value of 1 in all cells that are covered by a finer level, including one layer of ghost cells.
  a_coarMaskCF.define(dblCoar, m_nComp, IntVect::Zero);      // Does not need ghost cells
  a_coarMaskInvalid.define(dblCoar, m_nComp, IntVect::Unit); // Needs one ghost cell.

  for (DataIterator dit(dblCoar); dit.ok(); ++dit) {
    a_coarMaskCF[dit()].setVal(zero);
    a_coarMaskInvalid[dit()].setVal(zero);
  }

  // Ned to copy from the buffer masks to this mask.
  Copier copierCF;
  Copier copierInvalid;

  copierCF.ghostDefine(dblCoFi, dblCoar, domainCoar, IntVect::Unit);
  copierInvalid.ghostDefine(dblCoFi, dblCoar, domainCoar, IntVect::Zero, IntVect::Unit);

  coFiMaskCF.copyTo(interv, a_coarMaskCF, interv, copierCF, LDaddOp<FArrayBox>());
  coFiMaskInvalid.copyTo(interv, a_coarMaskInvalid, interv, copierInvalid, LDaddOp<FArrayBox>());
}

bool
EBGradient::isFiniteDifferenceStencilValid(const IntVect&   a_ivCoar,
                                           const EBISBox&   a_ebisBox,
                                           const FArrayBox& a_invalidRegion) const noexcept
{
  CH_TIME("EBGradient::isFiniteDifferenceStencilValid");

  // TLDR: Routine which checks if a finite difference stencil is "valid" in the sense that
  //       it does not reach into cut-cells in invalid regions. That is, the stencil is not allowed to reach
  //       into coarse-grid cut-cells that lie underneath the fine grid.

  constexpr Real zero = 0.0;

  bool validStencil = true;

  // Get all VoFs in the input cell
  Vector<VolIndex> coarVoFs = a_ebisBox.getVoFs(a_ivCoar);

  for (const auto& coarVoF : coarVoFs.stdVector()) {

    // Build a finite difference stencil.
    VoFStencil gradSten;
    for (int dir = 0; dir < SpaceDim; dir++) {
      VoFStencil derivSten;
      EBArith::getFirstDerivStencilWidthOne(derivSten, coarVoF, a_ebisBox, dir, 1.0, nullptr, dir);
      gradSten += derivSten;
    }

    // Check if the stencil is OK, i.e. that it does not reach into cut-cells.
    for (int i = 0; i < gradSten.size(); i++) {
      const VolIndex& ivof = gradSten.vof(i);
      const IntVect   iv   = ivof.gridIndex();

      const bool isCoveredByFinerCell = a_invalidRegion(iv, m_comp) > zero;
      const bool isIrregularCell      = a_ebisBox.isIrregular(iv);

      if (isCoveredByFinerCell && isIrregularCell) {
        validStencil = false;
      }
    }
  }

  return validStencil;
}

void
EBGradient::defineIteratorsEBCF(const LevelData<FArrayBox>& a_coarMaskCF,
                                const LevelData<FArrayBox>& a_coarMaskInvalid) noexcept
{
  CH_TIME("EBGradient::defineIteratorsEBCF");

  // TLDR: This defines which cells need explicit two-level stencils. We move along the refinement boundary on the coarse
  //       level and check if we can use finite-differencing. If not, we flag the cell and the box. We then generate the
  //       buffer grids, which consist of all boxes that contain at least one cell that requires a least squares stencil.

  constexpr Real zero = 0.0;

  const DisjointBoxLayout& dbl     = m_eblg.getDBL();
  const DisjointBoxLayout& dblFine = m_eblgFine.getDBL();

  const EBISLayout& ebisl     = m_eblg.getEBISL();
  const EBISLayout& ebislFine = m_eblgFine.getEBISL();

  const ProblemDomain& domain     = m_eblg.getDomain();
  const ProblemDomain& domainFine = m_eblgFine.getDomain();

  // Define iterators and stencils for EBCF
  m_ebcfIterator.define(dbl);
  m_ebcfStencilsCoar.define(dbl);
  m_ebcfStencilsFine.define(dbl);

  Vector<int> ebcfRanks;
  Vector<Box> ebcfBoxes;

  LayoutData<IntVectSet> sets(dbl);

  int localHasEBCF = 0;

  for (DataIterator dit(dbl); dit.ok(); ++dit) {
    const Box      cellBox = dbl[dit()];
    const EBISBox& ebisBox = ebisl[dit()];
    const EBGraph& ebgraph = ebisBox.getEBGraph();

    const bool isAllRegular = ebisBox.isAllRegular();
    const bool isAllCovered = ebisBox.isAllCovered();
    const bool isIrregular  = !isAllRegular && !isAllCovered;

    const FArrayBox& coarseFineRegion = a_coarMaskCF[dit()];
    const FArrayBox& invalidRegion    = a_coarMaskInvalid[dit()];

    // Determine cells where we need to drop order.
    IntVectSet& ebcfIVS = sets[dit()];
    ebcfIVS.makeEmpty();

    if (isIrregular) {

      // Iterate through the coarse-fine region and check if the finite difference stencil reaches into a cut-cell.
      auto kernel = [&](const IntVect& ivCoar) -> void {
        if (coarseFineRegion(ivCoar, m_comp) > zero) {
          const bool hasStencil = this->isFiniteDifferenceStencilValid(ivCoar, ebisBox, invalidRegion);

          // Ok, found a stencil that won't work.
          if (!hasStencil) {
            ebcfIVS |= ivCoar;

            localHasEBCF = 1;
          }
        }
      };

      // Launch kernel.
      BoxLoops::loop(cellBox, kernel);
    }

    // Define the iterator.
    m_ebcfIterator[dit()].define(ebcfIVS, ebgraph);
    m_ebcfStencilsCoar[dit()].define(ebcfIVS, ebgraph, m_nComp);
    m_ebcfStencilsFine[dit()].define(ebcfIVS, ebgraph, m_nComp);
  }

  // All ranks bware: An EBCF situation has gotten out of hand.
  m_hasEBCF = ParallelOps::sum(localHasEBCF) > 0;
}

void
EBGradient::defineStencilsEBCF(const LevelData<FArrayBox>& a_coarMaskInvalid) noexcept
{
  CH_TIMERS("EBGradient::defineStencilsEBCF");
  CH_TIMER("EBGradient::build_mask", t2);
  CH_TIMER("EBGradient::compute_stencils", t3);

  CH_assert(m_hasFine);
  CH_assert(m_hasEBCF);
  CH_assert(a_coarMaskInvalid.ghostVect() == IntVect::Unit);

  // TLDR: This code iterates through the buffer grids and defines least squares stencils in the
  //       various cells that need it. The input mask is used to determine which cells we can use
  //       on the fine/coarse levels.

  const DisjointBoxLayout& dbl     = m_eblg.getDBL();
  const DisjointBoxLayout& dblFine = m_eblgFiCo.getDBL();

  const EBISLayout& ebisl     = m_eblg.getEBISL();
  const EBISLayout& ebislFine = m_eblgFiCo.getEBISL();

  const ProblemDomain& domain = m_eblg.getDomain();

  for (DataIterator dit(dbl); dit.ok(); ++dit) {
    const Box cellBox      = dbl[dit()];
    const Box grownBox     = grow(cellBox, 1) & domain;
    const Box grownBoxFine = refine(grownBox, m_refRat);

    const EBISBox& ebisBox     = ebisl[dit()];
    const EBISBox& ebisBoxFine = ebislFine[dit()];

    const FArrayBox& coarMaskInvalid = a_coarMaskInvalid[dit()];

    VoFIterator&           ebcfIterator = m_ebcfIterator[dit()];
    BaseIVFAB<VoFStencil>& coarStencils = m_ebcfStencilsCoar[dit()];
    BaseIVFAB<VoFStencil>& fineStencils = m_ebcfStencilsFine[dit()];

    const BaseIVFAB<VoFStencil>& levelStencils = m_levelStencils[dit()];

    // Make the invalid region mask into a DenseIntVectSet because that's what LeastSquares wants.
    DenseIntVectSet validRegionCoar(grownBox, true);
    DenseIntVectSet validRegionFine(grownBoxFine, false);

    auto regularKernel = [&](const IntVect& iv) -> void {
      if (coarMaskInvalid(iv, m_comp) > 0.0) {
        validRegionCoar -= iv;
        validRegionFine |= refine(Box(iv, iv), m_refRat);
      }
    };

    CH_START(t2);
    BoxLoops::loop(grownBox, regularKernel);
    CH_STOP(t2);

    // Compute the stencils. Try to achieve the highest order we can.
    auto irregularKernel = [&](const VolIndex& vof) -> void {
      VoFStencil& coarStencil = coarStencils(vof, m_comp);
      VoFStencil& fineStencil = fineStencils(vof, m_comp);

      coarStencil.clear();
      fineStencil.clear();

      int  order        = m_order;
      bool foundStencil = false;

      while (order > 0 && !foundStencil) {
        foundStencil = this->getLeastSquaresStencil(coarStencil,
                                                    fineStencil,
                                                    vof,
                                                    m_dataLocation,
                                                    ebisl,
                                                    ebislFine,
                                                    dit(),
                                                    validRegionCoar,
                                                    validRegionFine,
                                                    m_dx,
                                                    m_dxFine,
                                                    order,
                                                    m_weighting);
        order--;
      }

      // As a last-ditch effort, get a finite difference stencil after all.
      if (!foundStencil) {
        coarStencil.clear();
        fineStencil.clear();

        for (int dir = 0; dir < SpaceDim; dir++) {
          VoFStencil derivDirStencil;
          EBArith::getFirstDerivStencilWidthOne(derivDirStencil, vof, ebisBox, dir, m_dx, nullptr, dir);
          coarStencil += derivDirStencil;
        }

        MayDay::Warning("CD_EBGradient::defineStencilsEBCF -- could not find stencil!");
      }
    };

    CH_START(t3);
    BoxLoops::loop(ebcfIterator, irregularKernel);
    CH_STOP(t3);
  }
}

bool
EBGradient::getFiniteDifferenceStencil(VoFStencil&            a_stencil,
                                       const VolIndex&        a_vof,
                                       const EBISBox&         a_ebisBox,
                                       const DenseIntVectSet& a_validRegion,
                                       const Real             a_dx) const noexcept
{
  CH_TIME("EBGradient::getFiniteDifferenceStencil");

  // TLDR: This routine computes a finite difference stencil. If it does not have any valid cells to use (defined by a_validRegion)
  //       it will not return a stencil.

  a_stencil.clear();

  bool foundStencil = false;

  const IntVect iv = a_vof.gridIndex();

  for (int dir = 0; dir < SpaceDim; dir++) {
    VoFStencil derivStencil;

    // Check for cells to the low/high side of this one.
    const IntVect ivLo = iv - BASISV(dir);
    const IntVect ivHi = iv + BASISV(dir);

    const bool isLoCellCovered = a_validRegion[ivLo];
    const bool isHiCellCovered = a_validRegion[ivHi];

    const bool isLoCellIrregular = a_ebisBox.isIrregular(ivLo);
    const bool isHiCellIrregular = a_ebisBox.isIrregular(ivHi);

    const Vector<VolIndex> loVoFs = a_ebisBox.getVoFs(a_vof, dir, Side::Lo, 1);
    const Vector<VolIndex> hiVoFs = a_ebisBox.getVoFs(a_vof, dir, Side::Hi, 1);

    const int numLoVoFs = loVoFs.size();
    const int numHiVoFs = hiVoFs.size();

    const bool useLoCell = (numLoVoFs > 0) && !(isLoCellIrregular && isLoCellCovered);
    const bool useHiCell = (numHiVoFs > 0) && !(isHiCellIrregular && isHiCellCovered);

    // Use centered differencing in dir-direction if we can. Otherwise, drop order to forward/backward differences.
    if (useLoCell && useHiCell) {
      for (int i = 0; i < numHiVoFs; i++) {
        derivStencil.add(hiVoFs[i], 1. / numHiVoFs, dir);
      }

      for (int i = 0; i < numLoVoFs; i++) {
        derivStencil.add(loVoFs[i], -1. / numLoVoFs, dir);
      }

      derivStencil *= 1. / (2.0 * a_dx);
    }
    else if (useLoCell) {
      derivStencil.add(a_vof, 1.0, dir);

      for (int i = 0; i < numLoVoFs; i++) {
        derivStencil.add(loVoFs[i], -1. / numLoVoFs, dir);
      }

      derivStencil *= 1. / a_dx;
    }
    else if (useHiCell) {
      derivStencil.add(a_vof, -1.0, dir);

      for (int i = 0; i < numHiVoFs; i++) {
        derivStencil.add(hiVoFs[i], 1. / numHiVoFs, dir);
      }

      derivStencil *= 1. / a_dx;
    }

    if (useLoCell || useHiCell) {
      foundStencil = true;

      a_stencil += derivStencil;
    }
    else {
      foundStencil = false;
    }
  }

  return foundStencil;
}

bool
EBGradient::getLeastSquaresStencil(VoFStencil&            a_stencilCoar,
                                   VoFStencil&            a_stencilFine,
                                   const VolIndex&        a_vofCoar,
                                   const CellLocation&    a_dataLocation,
                                   const EBISLayout&      a_ebislCoar,
                                   const EBISLayout&      a_ebislFine,
                                   const DataIndex&       a_dit,
                                   const DenseIntVectSet& a_validCellsCoar,
                                   const DenseIntVectSet& a_validCellsFine,
                                   const Real&            a_dxCoar,
                                   const Real&            a_dxFine,
                                   const int&             a_order,
                                   const int&             a_weight) const noexcept
{
  CH_TIME("EBGradient::getLeastSquaresStencil");

  // TLDR: This routine computes a two-level least squares stencil. We include cells that are within a range of 1 on the coarse level, but only
  //       those in a monotone path from the input vof. The fine cells are defined as the cells that are available through a coarsening of the
  //       coarse cells. Since the resulting system might sometimes be too big, we trim the system down to a specified size. Once the least squares
  //       system has been "solved", i.e. a minimization of the truncation order for the various expansions has been achieved, we formulate the stencil
  //       using the output of LeastSquares::computeDualLevelStencils. Note that the solution on the input vof is known, so it is pruned from the
  //       system of equations.

  bool foundStencil;

  a_stencilCoar.clear();
  a_stencilFine.clear();

  const EBISBox& ebisBoxCoar = a_ebislCoar[a_dit];
  const EBISBox& ebisBoxFine = a_ebislFine[a_dit];

  // Get all coarse cells within radius 1.
  const int coarRadius = 1;

  // Get coar vofs in range. The fine VoFs are defined as the VoFs that are available through refinement of the coarse vofs.
  Vector<VolIndex> coarVoFs =
    VofUtils::getVofsInRadius(a_vofCoar, ebisBoxCoar, coarRadius, VofUtils::Connectivity::MonotonePath, false);
  Vector<VolIndex> fineVoFs;
  for (int i = 0; i < coarVoFs.size(); i++) {
    fineVoFs.append(a_ebislCoar.refine(coarVoFs[i], m_refRat, a_dit));
  }

  // Only unique, in case something went wrong.
  //  VofUtils::onlyUnique(coarVoFs);
  VofUtils::onlyUnique(fineVoFs);

  // Only include fine VoFs that lie with a_validCellsFine.
  VofUtils::includeCells(coarVoFs, a_validCellsCoar);
  VofUtils::includeCells(fineVoFs, a_validCellsFine);

  // See if we have enough cells to solve the system of equations. The "-1" is because we interpolate to the cell center/centroid, but we already know
  // phi at this point.
  const int numEquations = fineVoFs.size() + coarVoFs.size();
  const int numUnknowns  = LeastSquares::getTaylorExpansionSize(a_order) - 1;

  if (numEquations > numUnknowns) {

    // In many cases we will have WAY too many equations for the specified order. This is particularly true in 3D
    // because the number of coar vofs included in a radius r from the ghost vof can be 3^3 = 27. In addition, if we
    // use a refinement factor of 4 and just half of those cells are covered by finer cells, we end up with a system
    // size of = 13 + 14*4^3 = 97. That's way more than we need for order 2 which requires 10 equations. Since the SVD
    // decomposition scales as O(n^3), we trim the system size to bring the cost down, discarding the cells that are
    // furthest away.
    std::vector<VolIndex>& fineVofsTrimmedSize = fineVoFs.stdVector();
    std::vector<VolIndex>& coarVofsTrimmedSize = coarVoFs.stdVector();

    // Coordinates of the vof that we will interpolate to (excluding lower-left corner because of the subtraction
    // in the comparators).
    const RealVect x0 = Location::position(a_dataLocation, a_vofCoar, ebisBoxCoar, a_dxCoar);

    // Things we need to capture
    const auto& loc         = a_dataLocation;
    const auto& p           = x0;
    const auto& ebisboxFine = ebisBoxFine;
    const auto& ebisboxCoar = ebisBoxCoar;
    const auto& dxFine      = a_dxFine;
    const auto& dxCoar      = a_dxCoar;

    // For sorting fine vofs, based on distance to the ghost vof. Shortest distance goes first.
    auto comparatorFine = [&loc, &p, &ebisboxFine, &dxFine](const VolIndex& v1, const VolIndex& v2) -> bool {
      const RealVect d1 = Location::position(loc, v1, ebisboxFine, dxFine) - p;
      const RealVect d2 = Location::position(loc, v2, ebisboxFine, dxFine) - p;

      const Real l1 = d1.vectorLength();
      const Real l2 = d2.vectorLength();

      return l1 < l2;
    };

    // For sorting coar vofs, based on distance to the ghost vof. Shortest distance goes first.
    auto comparatorCoar = [&loc, &p, &ebisboxCoar, &dxCoar](const VolIndex& v1, const VolIndex& v2) -> bool {
      const RealVect d1 = Location::position(loc, v1, ebisboxCoar, dxCoar) - p;
      const RealVect d2 = Location::position(loc, v2, ebisboxCoar, dxCoar) - p;

      const Real l1 = d1.vectorLength();
      const Real l2 = d2.vectorLength();

      return l1 < l2;
    };

    // Sort and trim the system size.
    std::sort(fineVofsTrimmedSize.begin(), fineVofsTrimmedSize.end(), comparatorFine);
    std::sort(coarVofsTrimmedSize.begin(), coarVofsTrimmedSize.end(), comparatorCoar);

    const int curFineSize = fineVofsTrimmedSize.size();
    const int curCoarSize = coarVofsTrimmedSize.size();

    fineVofsTrimmedSize.resize(std::min(2 * numUnknowns, curFineSize));
    coarVofsTrimmedSize.resize(std::min(2 * numUnknowns, curCoarSize));

    // Build the displacement vectors
    Vector<RealVect> fineDisplacements;
    Vector<RealVect> coarDisplacements;

    for (const auto& coarVoF : coarVoFs.stdVector()) {
      coarDisplacements.push_back(
        LeastSquares::displacement(a_dataLocation, a_dataLocation, a_vofCoar, coarVoF, ebisBoxCoar, a_dxCoar));
    }

    for (const auto& fineVoF : fineVoFs.stdVector()) {
      fineDisplacements.push_back(LeastSquares::displacement(a_dataLocation,
                                                             a_dataLocation,
                                                             a_vofCoar,
                                                             fineVoF,
                                                             ebisBoxCoar,
                                                             ebisBoxFine,
                                                             a_dxCoar,
                                                             a_dxFine));
    }

    // Now solve the least squares system for the gradient at the input VoF. After this we have found an expression which minimizes the
    // truncation error in the neighborhood of the VoF. The result is is a two-level stencil consisting of cells on both the fine and
    // the coarse levels.
    // These are the unknown terms that we want stencils for.
    IntVectSet derivTerms;
    IntVectSet knownTerms;

    for (int dir = 0; dir < SpaceDim; dir++) {
      derivTerms |= BASISV(dir);
    }
    knownTerms |= IntVect::Zero;

    const std::map<IntVect, std::pair<VoFStencil, VoFStencil>> stencils = LeastSquares::computeDualLevelStencils<
      float>(derivTerms, knownTerms, fineVoFs, coarVoFs, fineDisplacements, coarDisplacements, a_weight, a_order);
    // LeastSquares returns a map over all derivatives (unknowns) in the Taylor series. These are stored
    // as IntVects so that IntVect(1,1) = d^2/(dxdy) and so on. We are after IntVect(1,0) = d/dx, IntVect(0,1) = d/dy. We fetch
    // those and place them in a_stencilFine and a_stencilCoar. We encode the direction in the stencil variable (in a_stencilFine) and
    // a_stencilCoar.
    for (int dir = 0; dir < SpaceDim; dir++) {
      const IntVect deriv = BASISV(dir);

      const VoFStencil& fineDerivStencil = stencils.at(deriv).first;
      const VoFStencil& coarDerivStencil = stencils.at(deriv).second;

      for (int i = 0; i < fineDerivStencil.size(); i++) {
        const VolIndex& vof    = fineDerivStencil.vof(i);
        const Real&     weight = fineDerivStencil.weight(i);

        a_stencilFine.add(vof, weight, dir);
      }

      for (int i = 0; i < coarDerivStencil.size(); i++) {
        const VolIndex& vof    = coarDerivStencil.vof(i);
        const Real&     weight = coarDerivStencil.weight(i);

        a_stencilCoar.add(vof, weight, dir);
      }
    }

    // We have modified the right-hand side of the least squares system by pruning a_vofCoar from the system of equations (it's value is known). So,
    // our least squares solution is actually something like
    //
    //    [df/dx]     [A11 A12 ... ] [f(x0) - f(x)]
    //    [df/dy]  =  [A21 A22 ... ] [f(x1) - f(x)]
    //    [df/dz]     [A31 A32 ... ] [f(x2) - f(x)]
    //
    // where f(x) = phi(a_vofCoar). LeastSquares::computeGradSten does not case about the modified right-hand side and the stencils it returns
    // only account for f(x0), f(x1) etc. We need to add in the contribution from a_vofCoar, which is fortunately just the sum of the weights for
    // each derivative.
    for (int dir = 0; dir < SpaceDim; dir++) {
      Real coarVofWeight = 0.0;

      for (int i = 0; i < a_stencilFine.size(); i++) {
        const int  curVar    = a_stencilFine.variable(i);
        const Real curWeight = a_stencilFine.weight(i);

        if (curVar == dir) {
          coarVofWeight += curWeight;
        }
      }

      for (int i = 0; i < a_stencilCoar.size(); i++) {
        const int  curVar    = a_stencilCoar.variable(i);
        const Real curWeight = a_stencilCoar.weight(i);

        if (curVar == dir) {
          coarVofWeight += curWeight;
        }
      }

      a_stencilCoar.add(a_vofCoar, -1.0 * coarVofWeight, dir);
    }

    foundStencil = true;
  }
  else {
    foundStencil = false;
  }

  return foundStencil;
}

void
EBGradient::makeAggStencils() noexcept
{
  CH_TIMERS("EBGradient::makeAggStencils");
  CH_TIMER("EBGradient::makeAggStencils::define_buffers", t1);
  CH_TIMER("EBGradient::makeAggStencils::dbl_define", t2);
  CH_TIMER("EBGradient::makeAggStencils::populate_stencil", t3);
  CH_TIMER("EBGradient::makeAggStencils::define_aggstencil", t4);

  // Make some proxies for the input/output data holders
  const DisjointBoxLayout& dbl   = m_eblg.getDBL();
  const EBISLayout&        ebisl = m_eblg.getEBISL();

  CH_START(t1);
  LevelData<EBCellFAB> proxy(dbl, 1, m_ghostVector, EBCellFactory(ebisl));
  CH_STOP(t1);

  for (int dir = 0; dir < SpaceDim; dir++) {
    CH_START(t2);
    m_aggLevelStencils[dir].define(dbl);
    CH_STOP(t2);

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      const BaseIVFAB<VoFStencil>& rawStencils = m_levelStencils[dit()];

      // Populate stencils
      CH_START(t3);
      Vector<RefCountedPtr<BaseIndex>>   dstBaseIndex;
      Vector<RefCountedPtr<BaseStencil>> dstBaseStencil;

      auto kernel = [&](const VolIndex& vof) -> void {
        const VoFStencil& gradStencil = rawStencils(vof, m_comp);

        // Extract the stencil in the dir-direction only.
        VoFStencil derivDirStencil;
        for (int i = 0; i < gradStencil.size(); i++) {
          if (gradStencil.variable(i) == dir) {
            derivDirStencil.add(gradStencil.vof(i), gradStencil.weight(i));
          }
        }

        dstBaseIndex.push_back(RefCountedPtr<BaseIndex>(new VolIndex(vof)));
        dstBaseStencil.push_back(RefCountedPtr<BaseStencil>(new VoFStencil(derivDirStencil)));
      };

      BoxLoops::loop(m_levelIterator[dit()], kernel);
      CH_STOP(t3);

      CH_START(t4);
      m_aggLevelStencils[dir][dit()] = RefCountedPtr<AggStencil<EBCellFAB, EBCellFAB>>(
        new AggStencil<EBCellFAB, EBCellFAB>(dstBaseIndex, dstBaseStencil, proxy[dit()], proxy[dit()]));
      CH_STOP(t4);
    }
  }
}

#include <CD_NamespaceFooter.H>
