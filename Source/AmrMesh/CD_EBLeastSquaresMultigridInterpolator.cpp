/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBLeastSquaresMultigridInterpolator.cpp
  @brief  Implementation of CD_EBLeastSquaresMultigridInterpolator.H
  @author Robert Marskar
*/

// Std includes
#include <sstream>

// Chombo includes
#include <CH_Timer.H>
#include <EBCellFactory.H>
#include <NeighborIterator.H>
#include <EBAlias.H>
#include <ParmParse.H>

// Our includes
#include <CD_EBLeastSquaresMultigridInterpolator.H>
#include <CD_VofUtils.H>
#include <CD_ParallelOps.H>
#include <CD_LeastSquares.H>
#include <CD_BoxLoops.H>
#include <CD_NamespaceHeader.H>

constexpr int EBLeastSquaresMultigridInterpolator::m_stenComp;
constexpr int EBLeastSquaresMultigridInterpolator::m_numStenComp;
constexpr int EBLeastSquaresMultigridInterpolator::m_comp;

EBLeastSquaresMultigridInterpolator::EBLeastSquaresMultigridInterpolator(const EBLevelGrid& a_eblgFine,
                                                                         const EBLevelGrid& a_eblgCoFi,
                                                                         const EBLevelGrid& a_eblgCoar,
                                                                         const CellLocation a_dataLocation,
                                                                         const IntVect&     a_ghostVector,
                                                                         const int          a_refRat,
                                                                         const int          a_ghostCF,
                                                                         const int          a_order,
                                                                         const int          a_weighting) noexcept
{
  CH_TIMERS("EBLeastSquaresMultigridInterpolator::EBLeastSquaresMultigridInterpolator");
  CH_TIMER("EBLeastSquaresMultigridInterpolator::define_regular", t1);
  CH_TIMER("EBLeastSquaresMultigridInterpolator::define_irregular", t2);

  CH_assert(a_ghostCF > 0);
  CH_assert(a_refRat >= 2);
  CH_assert(a_refRat % 2 == 0);
  CH_assert(a_eblgFine.getGhost() >= a_ghostVector.max());

  // Not enough ghost cells in input eblg. I don't know why anyone would do that, but make sure we abort if it happens. This is a more
  // transparent error message than the assertion.
  if (a_eblgFine.getGhost() < a_ghostVector.max()) {
    MayDay::Error("EBLeastSquaresMultigridInterpolator::EBLeastSquaresMultigridInterpolator - not enough ghost cells!");
  }

  // Build the regular stencil objects for regular-grid interpolation. I also leave a timer in place until performance and
  // scalability has been investigated. To check performance, add EBLeastSquaresMultigridInterpolator.profile=true to your input script.
  CH_START(t1);
  m_refRat          = a_refRat;
  m_ghostCF         = a_ghostCF;
  m_order           = a_order;
  m_ghostVectorFine = a_ghostVector;
  m_dataLocation    = a_dataLocation;
  m_weight          = a_weighting;
  m_eblgFine        = a_eblgFine;
  m_eblgCoFi        = a_eblgCoFi;
  m_eblgCoar        = a_eblgCoar;

  const int fineRadius = std::max(2, m_order);
  const int coarRadius = std::max(2, fineRadius / m_refRat);

  m_ghostVectorCoFi = coarRadius * IntVect::Unit;
  CH_STOP(t1);

  CH_START(t2);
  this->defineGhostRegions();
  this->defineBuffers();
  this->defineCoarseInterp();
  this->defineStencilsEBCF();
  CH_STOP(t2);
}

EBLeastSquaresMultigridInterpolator::~EBLeastSquaresMultigridInterpolator() noexcept
{
  CH_TIME("EBLeastSquaresMultigridInterpolator::~EBLeastSquaresMultigridInterpolator");
}

int
EBLeastSquaresMultigridInterpolator::getGhostCF() const noexcept
{
  return m_ghostCF;
}

std::pair<DerivStencil, DerivStencil>
EBLeastSquaresMultigridInterpolator::getInterpolationStencilRegular(const IntVect&       a_fineGhost,
                                                                    const DataIndex&     a_din,
                                                                    const int            a_dir,
                                                                    const Side::LoHiSide a_side) const noexcept
{
  CH_TIME("EBLeastSquaresMultigridInterpolator::getInterpolationStencilRegular");

  // TLDR: This version does the same type of interpolation as in regularCoarseFineInterp. I.e., it
  // constructs a Lagrange polynomial for fine-ghost interpolation, using the same centering where the
  // ghost cell is at the origin of a coordinate system.
  const IntVect fineIV = a_fineGhost;
  const IntVect coarIV = coarsen(a_fineGhost, m_refRat);

  const Real x0 = -0.5 * (m_refRat - 1);
  const Real xi = 0.0;
  const Real x1 = 1.0;
  const Real x2 = 2.0;

  const Real L0 = (xi - x1) * (xi - x2) / ((x0 - x1) * (x0 - x2));
  const Real L1 = (xi - x0) * (xi - x2) / ((x1 - x0) * (x1 - x2));
  const Real L2 = (xi - x0) * (xi - x1) / ((x2 - x0) * (x2 - x1));

  const int iHiLo = sign(a_side);

  DerivStencil fineStencil;
  DerivStencil coarStencil;

  const CoarseInterpQuadCF& coarseStencils = (a_side == Side::Lo) ? m_loCoarseInterpCF[a_dir][a_din]
                                                                  : m_hiCoarseInterpCF[a_dir][a_din];

  // Interpolation of the coarse-grid data to the line connecting the fine-grid data with the ghost cell. delta is the
  // displacement vector from coarse-grid cell to fine-grid ghost cell. Note that this is normalized by the coarse grid
  // cell size
  const RealVect delta = (RealVect(fineIV) - m_refRat * RealVect(coarIV) + 0.5 * (1.0 - m_refRat)) / m_refRat;

  coarStencil.accumulate(coarIV, 1.0);
  for (int d = 0; d < SpaceDim; d++) {
    if (d != a_dir) {
      const DerivStencil firstDeriv  = coarseStencils.getFirstDerivStencil(coarIV, d);
      const DerivStencil secondDeriv = coarseStencils.getSecondDerivStencil(coarIV, d);

      for (int i = 0; i < firstDeriv.size(); i++) {
        coarStencil.accumulate(firstDeriv.getIndex(i), firstDeriv.getWeight(i) * delta[d]);
      }
      for (int i = 0; i < secondDeriv.size(); i++) {
        coarStencil.accumulate(secondDeriv.getIndex(i), secondDeriv.getWeight(i) * 0.5 * delta[d] * delta[d]);
      }
    }
  }

#if CH_SPACEDIM == 3
  const DerivStencil mixedDeriv = coarseStencils.getMixedDerivStencil(coarIV);
  for (int d = 0; d < SpaceDim; d++) {
    if (d != a_dir) {
      for (int i = 0; i < mixedDeriv.size(); i++) {
        coarStencil.accumulate(mixedDeriv.getIndex(i), mixedDeriv.getWeight(i) * delta[d]);
      }
    }
  }
#endif

  // Final interpolant uses the coarse-side value and two fine-grid values inside the grid patch.
  coarStencil *= L0;
  fineStencil.accumulate(fineIV - iHiLo * BASISV(a_dir), L1);
  fineStencil.accumulate(fineIV - 2 * iHiLo * BASISV(a_dir), L2);

  return std::make_pair(fineStencil, coarStencil);
}

std::pair<VoFStencil, VoFStencil>
EBLeastSquaresMultigridInterpolator::getInterpolationStencilEB(const VolIndex&  a_fineGhost,
                                                               const DataIndex& a_din) const noexcept
{
  CH_TIME("EBLeastSquaresMultigridInterpolator::getInterpolationStencilEB");

  VoFStencil fineStencil;
  VoFStencil coarStencil;

  if (m_ghostCells[a_din].contains(a_fineGhost.gridIndex())) {
    fineStencil = m_fineStencils[a_din](a_fineGhost, m_comp);
    coarStencil = m_coarStencils[a_din](a_fineGhost, m_comp);
  }
  else {
    MayDay::Abort("EBLeastSquaresMultigridInterpolation::getInterpolationStencilEB - logic bust!");
  }

  return std::make_pair(fineStencil, coarStencil);
}

void
EBLeastSquaresMultigridInterpolator::coarseFineInterp(LevelData<EBCellFAB>&       a_phiFine,
                                                      const LevelData<EBCellFAB>& a_phiCoar,
                                                      const Interval              a_variables) const noexcept
{
  CH_TIMERS("EBLeastSquaresMultigridInterpolator::coarseFineInterp");
  CH_TIMER("EBLeastSquaresMultigridInterpolator::coarseFineInterp::coar_irreg", t1);
  CH_TIMER("EBLeastSquaresMultigridInterpolator::coarseFineInterp::fine_irreg", t2);

  CH_assert(m_ghostCF * IntVect::Unit <= a_phiFine.ghostVect());
  CH_assert(a_phiFine.nComp() > a_variables.end());
  CH_assert(a_phiCoar.nComp() > a_variables.end());

  if (a_phiFine.ghostVect() != m_ghostVectorFine) {
    MayDay::Error("EBLeastSquaresMultigridInterpolator::coarseFineInterp -- number of ghost cells do not match!");
  }

  const DataIterator dit  = a_phiFine.dataIterator();
  const int          nbox = dit.size();

  LevelData<EBCellFAB> phiCoFi(m_eblgCoFi.getDBL(), 1, m_ghostVectorCoFi, EBCellFactory(m_eblgCoFi.getEBISL()));

  // Interpolate all variables near the EB. We will copy a_phiCoar to phiCoFi which holds the data on the coarse grid cells around each fine-grid
  // patch. Note that phiCoFi provides a LOCAL view of the coarse grid around each fine-level patch, so we can apply the stencils directly.
  for (int icomp = a_variables.begin(); icomp <= a_variables.end(); icomp++) {
    const Interval srcInterv = Interval(icomp, icomp);
    const Interval dstInterv = Interval(m_comp, m_comp);

    // Copy to buffer
    a_phiCoar.copyTo(srcInterv, phiCoFi, dstInterv, m_copier);

    // Do regular interpolation as if the EB is not there.
    this->regularCoarseFineInterp(a_phiFine, phiCoFi, icomp, 0);

    // Go through each grid patch and the to-be-interpolated ghost cells across the refinement boundary. We simply
    // apply the stencils here.
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      EBCellFAB&       dstFine = a_phiFine[din];
      const EBCellFAB& srcFine = a_phiFine[din];
      const EBCellFAB& srcCoar = phiCoFi[din];

      // Apply the coarse and fine stencils
      constexpr int numComp = 1;
      CH_START(t1);
      m_aggCoarStencils[din]->apply(dstFine, srcCoar, m_comp, icomp, numComp, false);
      CH_STOP(t1);

      CH_START(t2);
      m_aggFineStencils[din]->apply(dstFine, srcFine, icomp, icomp, numComp, true);
      CH_STOP(t2);
    }
  }
}

void
EBLeastSquaresMultigridInterpolator::coarseFineInterpH(LevelData<EBCellFAB>& a_phiFine,
                                                       const Interval        a_variables) const noexcept
{
  CH_TIME("EBLeastSquaresMultigridInterpolator::coarseFineInterpH(LD<EBCellFAB>, Interval)");

  CH_assert(m_ghostCF * IntVect::Unit <= a_phiFine.ghostVect());
  CH_assert(a_phiFine.nComp() > a_variables.end());

  if (a_phiFine.ghostVect() != m_ghostVectorFine) {
    MayDay::Error("EBLeastSquaresMultigridInterpolator::coarseFineInterp -- number of ghost cells do not match!");
  }

  const DataIterator dit  = m_eblgFine.getDBL().dataIterator();
  const int          nbox = dit.size();

  // TLDR: This routine does the coarse-fine interpolation with the coarse-grid data set to zero.
#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    this->coarseFineInterpH(a_phiFine[din], a_variables, din);
  }
}

void
EBLeastSquaresMultigridInterpolator::coarseFineInterpH(EBCellFAB&       a_phi,
                                                       const Interval   a_variables,
                                                       const DataIndex& a_din) const noexcept
{
  CH_TIMERS("EBLeastSquaresMultigridInterpolator::coarseFineInterpH(LD<EBCellFAB>)");
  CH_TIMER("EBLeastSquaresMultigridInterpolator::regular_interp", t1);
  CH_TIMER("EBLeastSquaresMultigridInterpolator::irregular_interp", t2);

  CH_assert(a_phi.nComp() > a_variables.end());

  // TLDR: This routine does the coarse-fine interpolation with the coarse-grid data set to zero. This is the kernel version,
  //       operating on a grid patch. It first does a direct kernel for regular data, and then does the interpolation near the
  //       EB after that.
  const Real dxFine = 1.0;
  const Real dxCoar = 1.0 * m_refRat;
  const Real c1     = 2 * (dxCoar - dxFine) / (dxCoar + dxFine);
  const Real c2     = -(dxCoar - dxFine) / (dxCoar + 3 * dxFine);

  for (int ivar = a_variables.begin(); ivar <= a_variables.end(); ivar++) {

    // Do the regular interp on all sides of the current patch.
    CH_START(t1);
    for (int dir = 0; dir < SpaceDim; dir++) {
      for (SideIterator sit; sit.ok(); ++sit) {
        const IntVect shift = sign(sit()) * BASISV(dir);

        // Regular homogeneous interpolation on side sit() in direction dir
        const Box ghostBox = m_cfivs[a_din].at(std::make_pair(dir, sit()));

        if (!ghostBox.isEmpty()) {
          BaseFab<Real>& phiReg = a_phi.getSingleValuedFAB();

          // C++ kernel for homogeneous interpolation along a line, assumning that coarse-grid
          // data is zero.
          auto interpHomo = [&](const IntVect& iv) -> void {
            phiReg(iv, ivar) = c1 * phiReg(iv - shift, ivar) + c2 * phiReg(iv - 2 * shift, ivar);
          };

          // Apply the kernel.
          BoxLoops::loop(ghostBox, interpHomo);
        }
      }
    }
    CH_STOP(t1);

    // Apply fine stencil near the EB. It might look weird that we apply the stencil to a_phi AND put the result in a_phi. This
    // is because the stencils are defined in the ghost cells we will fill, but the stencils only reach into valid data. So, we
    // are, in fact, not writing to data that is used by the other stencils
    CH_START(t2);
    constexpr int numComp = 1;
    m_aggFineStencils[a_din]->apply(a_phi, a_phi, ivar, ivar, numComp, false);
    CH_STOP(t2);
  }
}

void
EBLeastSquaresMultigridInterpolator::defineGhostRegions() noexcept
{
  CH_TIME("EBLeastSquaresMultigridInterpolator::defineGhostRegions");

  const DisjointBoxLayout& dbl    = m_eblgFine.getDBL();
  const DataIterator&      dit    = dbl.dataIterator();
  const ProblemDomain&     domain = m_eblgFine.getDomain();
  const EBISLayout&        ebisl  = m_eblgFine.getEBISL();
  const int                nbox   = dit.size();

  // Define the "regular" ghost interpolation regions. This is just one cell wide since the operator stencil
  // has a width of 1 in regular cells.
  m_cfivs.define(dbl);

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    const Box cellBox = dbl[din];

    std::map<std::pair<int, Side::LoHiSide>, Box>& cfivsBoxes = m_cfivs[din];

    for (int dir = 0; dir < SpaceDim; dir++) {
      for (SideIterator sit; sit.ok(); ++sit) {

        IntVectSet cfivs = IntVectSet(adjCellBox(cellBox, dir, sit(), 1));

        NeighborIterator nit(dbl); // Subtract the other boxes if they intersect this box.
        for (nit.begin(din); nit.ok(); ++nit) {
          cfivs -= dbl[nit()];
        }

        cfivs &= domain;
        cfivs.recalcMinBox();

        cfivsBoxes.emplace(std::make_pair(dir, sit()), cfivs.minBox());
      }
    }
  }

  // This hook is for the interpolation over the coarse-fine boundary near the EB. In this case we may
  // require more ghost cells to be interpolated (defined by m_ghostCF). This routine computes those
  // cells, including all ghost cells that are within range m_ghostCF from the cut-cell.
  m_ghostCells.define(dbl);

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    const Box      cellBox = dbl[din];
    const EBISBox& ebisbox = ebisl[din];

    if (ebisbox.isAllRegular() || ebisbox.isAllCovered()) {
      m_ghostCells[din] = IntVectSet();
    }
    else {
      // 1. Define the width of the ghost layer region around current (irregular grid) patch
      Box grownBox = grow(cellBox, m_ghostCF);
      grownBox &= domain;

      m_ghostCells[din] = IntVectSet(grownBox);

      NeighborIterator nit(dbl);
      for (nit.begin(din); nit.ok(); ++nit) {
        m_ghostCells[din] -= dbl[nit()];
      }
      m_ghostCells[din] -= cellBox;

      // 2. Only include ghost cells that are within range m_ghostCF of an irregular grid cell
      IntVectSet irreg = ebisbox.getIrregIVS(cellBox);
      irreg.grow(m_ghostCF);
      m_ghostCells[din] &= irreg;
    }
  }
}

void
EBLeastSquaresMultigridInterpolator::defineBuffers() noexcept
{
  CH_TIME("EBLeastSquaresMultigridInterpolator::defineBuffers()");

  m_copier.define(m_eblgCoar.getDBL(), m_eblgCoFi.getDBL(), m_ghostVectorCoFi);
}

void
EBLeastSquaresMultigridInterpolator::defineCoarseInterp() noexcept
{
  CH_TIME("EBLeastSquaresMultigridInterpolator::defineCoarseInterp");

  const DisjointBoxLayout& dblFine    = m_eblgFine.getDBL();
  const DataIterator&      dit        = dblFine.dataIterator();
  const ProblemDomain&     domainCoar = m_eblgCoar.getDomain();
  const int                nbox       = dit.size();

  for (int dir = 0; dir < SpaceDim; dir++) {
    m_loCoarseInterpCF[dir].define(dblFine);
    m_hiCoarseInterpCF[dir].define(dblFine);

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      for (SideIterator sit; sit.ok(); ++sit) {
        const Side::LoHiSide side      = sit();
        const Box            interpBox = m_cfivs[din].at(std::make_pair(dir, side));

        CoarseInterpQuadCF& interp = (side == Side::Lo) ? m_loCoarseInterpCF[dir][din] : m_hiCoarseInterpCF[dir][din];

        interp.define(dblFine, domainCoar, din, interpBox, m_refRat, dir);
      }
    }
  }
}

void
EBLeastSquaresMultigridInterpolator::defineStencilsEBCF() noexcept
{
  CH_TIMERS("EBLeastSquaresMultigridInterpolator::defineStencilsEBCF");
  CH_TIMER("EBLeastSquaresMultigridInterpolator::define_grids", t1);
  CH_TIMER("EBLeastSquaresMultigridInterpolator::define_intvectsets", t2);
  CH_TIMER("EBLeastSquaresMultigridInterpolator::neighbor_loop", t3);
  CH_TIMER("EBLeastSquaresMultigridInterpolator::define_structures", t4);
  CH_TIMER("EBLeastSquaresMultigridInterpolator::compute_stencils", t5);

  // This routine defines stencils for all the ghost cells we need to fill across the EBCF boundary.
  const int comp = 0;

  const Real dxFine = 1.0;
  const Real dxCoar = dxFine * m_refRat;

  const ProblemDomain& domFine = m_eblgFine.getDomain();
  const ProblemDomain& domCoar = m_eblgCoFi.getDomain();

  const DisjointBoxLayout& dblFine = m_eblgFine.getDBL();
  const DisjointBoxLayout& dblCoar = m_eblgCoFi.getDBL();

  const DataIterator& ditFine = dblFine.dataIterator();
  const DataIterator& ditCoar = dblCoar.dataIterator();

  const EBISLayout& ebislFine = m_eblgFine.getEBISL();
  const EBISLayout& ebislCoar = m_eblgCoFi.getEBISL();

  const int nboxFine = ditFine.size();
  const int nboxCoar = ditCoar.size();

  CH_START(t1);
  m_fineStencils.define(dblFine);
  m_coarStencils.define(dblFine);
  m_ghostIterFine.define(dblFine);
  CH_STOP(t1);

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nboxFine; mybox++) {
    const DataIndex& din = ditFine[mybox];

    const Box origFineBox    = dblFine[din];
    const Box origCoarBox    = dblCoar[din];
    const Box ghostedFineBox = grow(origFineBox, m_ghostVectorFine);
    const Box grownCoarBox   = grow(origCoarBox, m_ghostVectorCoFi);

    // Define the valid regions such that the interpolation does not include coarse grid cells that fall beneath the fine level,
    // and no fine cells outside the CF.
    CH_START(t2);
    DenseIntVectSet validFineCells(origFineBox, true);
    DenseIntVectSet validCoarCells(grownCoarBox, true);

    // Subtract the coarse cells under the fine box.
    Box coarsenedFineBox = origFineBox;
    coarsenedFineBox.coarsen(m_refRat);
    validCoarCells -= coarsenedFineBox;
    CH_STOP(t2);

    // Same for parts of the current (grown) patch that overlaps with neighboring boxes.
    CH_START(t3);
    NeighborIterator nit(dblFine);
    for (nit.begin(din); nit.ok(); ++nit) {

      // = neighboring grid patch on the fine level. Use it's cells if they are also ghost cells in the current patch
      const Box neighborBoxFine = dblFine[nit()];

      // Overlapping region between the ghosted box and a neighboring patch (with only valid cells). We are
      // allowed to use these cells.
      const Box fineOverlap = ghostedFineBox & neighborBoxFine;
      validFineCells |= DenseIntVectSet(fineOverlap, true);

      // validCoarCells are NOT allowed to use the region the fine box.
      Box neighborBoxCoar = neighborBoxFine;
      neighborBoxCoar.coarsen(m_refRat);
      neighborBoxCoar &= grownCoarBox;

      //      validCoarCells -= DenseIntVectSet(neighborBoxCoar, true);
      validCoarCells -= neighborBoxCoar;
    }

    // Restrict to domain
    validFineCells &= domFine;
    validCoarCells &= domCoar;
    CH_STOP(t3);

    // Now go through each ghost cell and get an interpolation stencil to specified order.
    CH_START(t4);
    const EBISBox& ebisboxFine = m_eblgFine.getEBISL()[din];
    const EBISBox& ebisboxCoar = m_eblgCoFi.getEBISL()[din];

    const EBGraph& fineGraph = ebisboxFine.getEBGraph();

    m_fineStencils[din].define(m_ghostCells[din], fineGraph, m_numStenComp);
    m_coarStencils[din].define(m_ghostCells[din], fineGraph, m_numStenComp);
    m_ghostIterFine[din].define(m_ghostCells[din], fineGraph);
    CH_STOP(t4);

    // Build stencils.
    CH_START(t5);
    VoFIterator& vofit = m_ghostIterFine[din];

    auto kernel = [&](const VolIndex& ghostVofFine) -> void {
      const VolIndex& ghostVofCoar = ebislFine.coarsen(ghostVofFine, m_refRat, din);

      VoFStencil& fineSten = m_fineStencils[din](ghostVofFine, comp);
      VoFStencil& coarSten = m_coarStencils[din](ghostVofFine, comp);

      int  order        = m_order;
      bool foundStencil = false;

      // Try to fine a two-level stencil for this ghost cell. Drop order if we can't find it (and print a warning).
      while (order > 0 && !foundStencil) {
        foundStencil = this->getStencil(fineSten,
                                        coarSten,
                                        m_dataLocation,
                                        ghostVofFine,
                                        ghostVofCoar,
                                        ebisboxFine,
                                        ebisboxCoar,
                                        validFineCells,
                                        validCoarCells,
                                        dxFine,
                                        dxCoar,
                                        order,
                                        m_weight);

        order--;

        if (!foundStencil) {
          pout() << "EBLeastSquaresMultigridInterpolator -- on domain = " << m_eblgFine.getDomain()
                 << ", dropping order for vof = " << ghostVofFine << endl;
        }
      }

      // Drop to order 0 if we never found a stencil, and issue an error code.
      if (!foundStencil) {
        fineSten.clear();
        coarSten.clear();

        coarSten.add(ghostVofCoar, 1.0);

        pout()
          << "EBLeastSquaresMultigridInterpolator::defineStencilsEBCF -- could not find stencil and dropping to order 0"
          << endl;
      }
    };

    BoxLoops::loop(vofit, kernel);
    CH_STOP(t5);
  }

  // We now have all the stencils we need. Make them into an AggStencil for performance
  // optimization.
  this->makeAggStencils();
}

bool
EBLeastSquaresMultigridInterpolator::getStencil(VoFStencil&            a_stencilFine,
                                                VoFStencil&            a_stencilCoar,
                                                const CellLocation&    a_dataLocation,
                                                const VolIndex&        a_ghostVofFine,
                                                const VolIndex&        a_ghostVofCoar,
                                                const EBISBox&         a_ebisboxFine,
                                                const EBISBox&         a_ebisboxCoar,
                                                const DenseIntVectSet& a_validFineCells,
                                                const DenseIntVectSet& a_validCoarCells,
                                                const Real&            a_dxFine,
                                                const Real&            a_dxCoar,
                                                const int&             a_order,
                                                const int&             a_weight) const noexcept
{
  CH_TIMERS("EBLeastSquaresMultigridInterpolator::getStencil");
  CH_TIMER("get_vofs", t1);
  CH_TIMER("trim_vofs", t2);
  CH_TIMER("compute_displacements", t3);
  CH_TIMER("compute_stencil", t4);

  // On input, we know which ghost cell we want to interpolate to, and we happen to have a map of valid cells in a_validFineCells and a_validCoarCells. We use
  // that information to build an overdetermined linear system of equations that interpolate to a_ghostVofFine to order a_order. If we don't have enough equations,
  // this routine will return false, and will not create stencils.

  bool foundStencil = true;

  // I think these radii are good -- but there's no hard limit here. Increase the radii if you
  // see that the stencil drops order.
  CH_START(t1);
  const int fineRadius = std::max(2, a_order);
  const int coarRadius = std::max(2, fineRadius / m_refRat);

  Vector<VolIndex> fineVofs;
  Vector<VolIndex> coarVofs;

  // Get all Vofs in specified radii. Don't use cells that are not in a_validFineCells or in a_validCoarCells.
  fineVofs = VofUtils::getVofsInRadius(a_ghostVofFine,
                                       a_ebisboxFine,
                                       fineRadius,
                                       VofUtils::Connectivity::MonotonePath,
                                       false);
  coarVofs = VofUtils::getVofsInRadius(a_ghostVofCoar,
                                       a_ebisboxCoar,
                                       coarRadius,
                                       VofUtils::Connectivity::MonotonePath,
                                       true);

  VofUtils::includeCells(fineVofs, a_validFineCells);
  VofUtils::includeCells(coarVofs, a_validCoarCells);

  const int numEquations = coarVofs.size() + fineVofs.size();
  const int numUnknowns  = LeastSquares::getTaylorExpansionSize(a_order);
  CH_STOP(t1);

  if (numEquations >= numUnknowns) { // We have enough equations to get a stencil.
    // In many cases we will have WAY too many equations for the specified order. This is particularly true in 3D
    // because the number of coar vofs included in a radius r from the ghost vof can be (1 + 2*r)^3. So for r = 2
    // this = 125 cells, not counting the fine cells. Since singular value decomposition scales like O(n^3), the
    // penalty for including too may equations can be quite severe.
    CH_START(t2);
    std::vector<VolIndex>& fineVofsTrimmedSize = fineVofs.stdVector();
    std::vector<VolIndex>& coarVofsTrimmedSize = coarVofs.stdVector();

    // Coordinates of the ghost vof that we will interpolate to (excluding lower-left corner because of the subtraction
    // in the comparators).
    const RealVect x0 = Location::position(a_dataLocation, a_ghostVofFine, a_ebisboxFine, a_dxFine);

    // Things we need to capture
    const auto& loc         = a_dataLocation;
    const auto& p           = x0;
    const auto& ebisBoxFine = a_ebisboxFine;
    const auto& ebisBoxCoar = a_ebisboxCoar;
    const auto& dxFine      = a_dxFine;
    const auto& dxCoar      = a_dxCoar;

    // For sorting fine vofs, based on distance to the ghost vof. Shortest distance goes first.
    auto comparatorFine = [&loc, &p, &ebisBoxFine, &dxFine](const VolIndex& v1, const VolIndex& v2) -> bool {
      const RealVect d1 = Location::position(loc, v1, ebisBoxFine, dxFine) - p;
      const RealVect d2 = Location::position(loc, v2, ebisBoxFine, dxFine) - p;

      const Real l1 = d1.vectorLength();
      const Real l2 = d2.vectorLength();

      return l1 < l2;
    };

    // For sorting coar vofs, based on distance to the ghost vof. Shortest distance goes first.
    auto comparatorCoar = [&loc, &p, &ebisBoxCoar, &dxCoar](const VolIndex& v1, const VolIndex& v2) -> bool {
      const RealVect d1 = Location::position(loc, v1, ebisBoxCoar, dxCoar) - p;
      const RealVect d2 = Location::position(loc, v2, ebisBoxCoar, dxCoar) - p;

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
    CH_STOP(t2);

    // Build displacement vectors.
    CH_START(t3);
    Vector<RealVect> fineDisplacements;
    Vector<RealVect> coarDisplacements;

    for (const auto& fineVof : fineVofs.stdVector()) {
      fineDisplacements.push_back(
        LeastSquares::displacement(a_dataLocation, a_dataLocation, a_ghostVofFine, fineVof, a_ebisboxFine, a_dxFine));
    }

    for (const auto& coarVof : coarVofs.stdVector()) {
      coarDisplacements.push_back(LeastSquares::displacement(a_dataLocation,
                                                             a_dataLocation,
                                                             a_ghostVofFine,
                                                             coarVof,
                                                             a_ebisboxFine,
                                                             a_ebisboxCoar,
                                                             a_dxFine,
                                                             a_dxCoar));
    }
    CH_STOP(t3);

    // LeastSquares computes all unknown terms in a Taylor expansion up to specified order. We want the 0th order term, i.e. the interpolated value,
    // which in multi-index notation is the term (0,0), i.e. IntVect::Zero. The format of the two-level least squares routine is such that the
    // fine stencil lies on the first index. This can be confusing, but the LeastSquares uses a very compact notation.
    CH_START(t4);
    IntVect    interpStenIndex = IntVect::Zero;
    IntVectSet derivs          = IntVectSet(interpStenIndex);
    IntVectSet knownTerms      = IntVectSet();

    std::map<IntVect, std::pair<VoFStencil, VoFStencil>> stencils = LeastSquares::computeDualLevelStencils<
      Real>(derivs, knownTerms, fineVofs, coarVofs, fineDisplacements, coarDisplacements, a_weight, a_order);

    a_stencilFine = stencils.at(interpStenIndex).first;
    a_stencilCoar = stencils.at(interpStenIndex).second;

    foundStencil = true;
    CH_STOP(t4);
  }
  else {
    foundStencil = false;
  }

  return foundStencil;
}

void
EBLeastSquaresMultigridInterpolator::makeAggStencils() noexcept
{
  CH_TIME("EBLeastSquaresMultigridInterpolator::makeAggStencils");

  const DisjointBoxLayout& dblFine = m_eblgFine.getDBL();
  const DisjointBoxLayout& dblCoFi = m_eblgCoFi.getDBL();

  const DataIterator ditFine  = dblFine.dataIterator();
  const int          nboxFine = ditFine.size();

  // Make some proxies.
  LevelData<EBCellFAB> phiProxyFine(dblFine, 1, m_ghostVectorFine, EBCellFactory(m_eblgFine.getEBISL()));
  LevelData<EBCellFAB> phiProxyCoar(dblCoFi, 1, m_ghostVectorCoFi, EBCellFactory(m_eblgCoFi.getEBISL()));

  // Define the fine-grid stencils.
  m_aggFineStencils.define(dblFine);
  m_aggCoarStencils.define(dblCoFi);

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nboxFine; mybox++) {
    const DataIndex& din = ditFine[mybox];

    Vector<RefCountedPtr<BaseIndex>>   dstBaseIndex;
    Vector<RefCountedPtr<BaseStencil>> dstBaseStencilFine;
    Vector<RefCountedPtr<BaseStencil>> dstBaseStencilCoar;

    VoFIterator& vofit = m_ghostIterFine[din];

    auto kernel = [&](const VolIndex& vofFine) -> void {
      const VoFStencil& stencilFine = m_fineStencils[din](vofFine, m_stenComp);
      const VoFStencil& stencilCoar = m_coarStencils[din](vofFine, m_stenComp);

      dstBaseIndex.push_back(RefCountedPtr<BaseIndex>(new VolIndex(vofFine)));
      dstBaseStencilFine.push_back(RefCountedPtr<BaseStencil>(new VoFStencil(stencilFine)));
      dstBaseStencilCoar.push_back(RefCountedPtr<BaseStencil>(new VoFStencil(stencilCoar)));
    };

    BoxLoops::loop(vofit, kernel);

    m_aggFineStencils[din] = RefCountedPtr<AggStencil<EBCellFAB, EBCellFAB>>(
      new AggStencil<EBCellFAB, EBCellFAB>(dstBaseIndex, dstBaseStencilFine, phiProxyFine[din], phiProxyFine[din]));

    m_aggCoarStencils[din] = RefCountedPtr<AggStencil<EBCellFAB, EBCellFAB>>(
      new AggStencil<EBCellFAB, EBCellFAB>(dstBaseIndex, dstBaseStencilCoar, phiProxyCoar[din], phiProxyFine[din]));
  }
}

void
EBLeastSquaresMultigridInterpolator::regularCoarseFineInterp(LevelData<EBCellFAB>&       a_finePhi,
                                                             const LevelData<EBCellFAB>& a_coarPhi,
                                                             const int                   a_fineVar,
                                                             const int                   a_coarVar) const noexcept

{
  CH_TIMERS("EBLeastSquaresMultigridInterpolator::regularCoarseFineInterp");
  CH_TIMER("EBLeastSquaresMultigridInterpolator::regularCoarseFineInterp::coarse_interp", t1);
  CH_TIMER("EBLeastSquaresMultigridInterpolator::regularCoarseFineInterp::fine_interp", t2);

  CH_assert(a_finePhi.nComp() > a_fineVar);
  CH_assert(a_coarPhi.nComp() > a_coarVar);

  // Compute Lagrangian polynomial interpolant for fine interpolation. For the fine interpolation we assume that
  // the ghost cell has "coordinate" x = 0, and we compute the Lagrangian interpolant L0(x), L1(x), etc. The other
  // cells that are involved have positions as follows:
  const Real x0 = -0.5 * (m_refRat - 1);
  const Real xi = 0.0;
  const Real x1 = 1.0;
  const Real x2 = 2.0;

  const Real L0 = (xi - x1) * (xi - x2) / ((x0 - x1) * (x0 - x2));
  const Real L1 = (xi - x0) * (xi - x2) / ((x1 - x0) * (x1 - x2));
  const Real L2 = (xi - x0) * (xi - x1) / ((x2 - x0) * (x2 - x1));

  const DisjointBoxLayout& dblFine    = m_eblgFine.getDBL();
  const ProblemDomain&     domainCoar = m_eblgCoFi.getDomain();
  const DataIterator&      ditFine    = dblFine.dataIterator();
  const int                nboxFine   = ditFine.size();

  // We are interpolating the first layer of ghost cells to O(h^3). To do this, we must first do an interpolation on the
  // coarse grid, and then cubic interpolation on the fine grid.
#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nboxFine; mybox++) {
    const DataIndex& din = ditFine[mybox];

    const Box fineBox = dblFine[din];

    FArrayBox&       finePhi = a_finePhi[din].getFArrayBox();
    const FArrayBox& coarPhi = a_coarPhi[din].getFArrayBox();

    for (int dir = 0; dir < SpaceDim; dir++) {
      for (SideIterator sit; sit.ok(); ++sit) {
        const int iHiLo     = sign(sit());
        const Box interpBox = m_cfivs[din].at(std::make_pair(dir, sit()));

        // Coarse-side interpolation stencil. This does interpolation orthogonal to direction 'dir'
        const CoarseInterpQuadCF& coarseStencils = (sit() == Side::Lo) ? m_loCoarseInterpCF[dir][din]
                                                                       : m_hiCoarseInterpCF[dir][din];

        // Adds first derivative to the Taylor expansion.
        auto applyDerivs = [&](const IntVect& fineIV) -> void {
          const IntVect coarIV = coarsen(fineIV, m_refRat);

          finePhi(fineIV, a_fineVar) = coarPhi(coarIV, a_coarVar);

          // Displacement vector from coarse-grid cell to fine-grid ghost cell. Note that this is normalized by
          // the coarse grid cell size.
          const RealVect delta = (RealVect(fineIV) - m_refRat * RealVect(coarIV) + 0.5 * (1.0 - m_refRat)) / m_refRat;

          // Add contributions from first and second derivatives.
          for (int d = 0; d < SpaceDim; d++) {
            if (d != dir) {
              const Real firstDeriv  = coarseStencils.computeFirstDeriv(coarPhi, coarIV, d, a_coarVar);
              const Real secondDeriv = coarseStencils.computeSecondDeriv(coarPhi, coarIV, d, a_coarVar);

              finePhi(fineIV, a_fineVar) += firstDeriv * delta[d] + 0.5 * secondDeriv * delta[d] * delta[d];
            }
          }

#if CH_SPACEDIM == 3
          // Add contribution from mixed derivative. Only in 3D.
          Real mixedDeriv = coarseStencils.computeMixedDeriv(coarPhi, coarIV, a_coarVar);
          for (int d = 0; d < SpaceDim; d++) {
            if (d != dir) {
              mixedDeriv *= delta[d];
            }
          }
          finePhi(fineIV, a_fineVar) += mixedDeriv;
#endif
        };

        // We've put the coarse-grid interpolation into finePhi(fineIV, a_fineVar). Now use that value
        // when doing quadratic interpolation with the additional fine-grid data.
        auto interpOnFine = [&](const IntVect& fineIV) -> void {
          const Real phi0 = finePhi(fineIV, a_fineVar);
          const Real phi1 = finePhi(fineIV - iHiLo * BASISV(dir), a_fineVar);
          const Real phi2 = finePhi(fineIV - 2 * iHiLo * BASISV(dir), a_fineVar);

          finePhi(fineIV, a_fineVar) = phi0 * L0 + phi1 * L1 + phi2 * L2;
        };

        CH_START(t1);
        BoxLoops::loop(interpBox, applyDerivs);
        CH_STOP(t1);

        CH_START(t2);
        BoxLoops::loop(interpBox, interpOnFine);
        CH_STOP(t2);
      }
    }
  }
}

#include <CD_NamespaceFooter.H>
