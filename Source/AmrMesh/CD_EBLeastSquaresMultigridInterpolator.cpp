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
#include <EBCellFactory.H>
#include <EBLevelDataOps.H>
#include <NeighborIterator.H>
#include <EBAlias.H>
#include <CH_Timer.H>
#include <ParmParse.H>

// Our includes
#include <CD_Timer.H>
#include <CD_EBLeastSquaresMultigridInterpolator.H>
#include <CD_VofUtils.H>
#include <CD_LeastSquares.H>
#include <CD_BoxLoops.H>
#include <CD_NamespaceHeader.H>

constexpr int EBLeastSquaresMultigridInterpolator::m_stenComp;
constexpr int EBLeastSquaresMultigridInterpolator::m_numStenComp;
constexpr int EBLeastSquaresMultigridInterpolator::m_comp;

EBLeastSquaresMultigridInterpolator::EBLeastSquaresMultigridInterpolator(const EBLevelGrid& a_eblgFine,
                                                                         const EBLevelGrid& a_eblgCoar,
                                                                         const CellLocation a_dataLocation,
                                                                         const IntVect&     a_ghostVector,
                                                                         const int          a_refRat,
                                                                         const int          a_ghostCF,
                                                                         const int          a_order,
                                                                         const int          a_weighting)
{
  CH_TIME("EBLeastSquaresMultigridInterpolator::EBLeastSquaresMultigridInterpolator");

  CH_assert(a_ghostCF > 0);
  CH_assert(a_refRat % 2 == 0);

  const DisjointBoxLayout& gridsFine = a_eblgFine.getDBL();
  const DisjointBoxLayout& gridsCoar = a_eblgCoar.getDBL();

  // Build the regular stencil objects for regular-grid interpolation. For now, we use QuadCFInterp for this. I also leave
  // a timer in place until performance and scalability has been investigated. To check performance, add
  // EBLeastSquaresMultigridInterpolator.profile=true to your input script.
  Timer timer("EBLeastSquaresMultigridInterpolator::EBLeastSquaresMultigridInterpolator");

  m_refRat       = a_refRat;
  m_ghostCF      = a_ghostCF;
  m_order        = a_order;
  m_ghostVector  = a_ghostVector;
  m_dataLocation = a_dataLocation;
  m_weight       = a_weighting;

  timer.startEvent("Define QuadCFInterp");
  //  QuadCFInterp::define(gridsFine, &gridsCoar, 1.0, a_refRat, 1, a_eblgFine.getDomain());
  m_scalarInterpolator =
    RefCountedPtr<QuadCFInterp>(new QuadCFInterp(gridsFine, &gridsCoar, 1.0, a_refRat, 1, a_eblgFine.getDomain()));
  m_vectorInterpolator = RefCountedPtr<QuadCFInterp>(
    new QuadCFInterp(gridsFine, &gridsCoar, 1.0, a_refRat, SpaceDim, a_eblgFine.getDomain()));
  timer.stopEvent("Define QuadCFInterp");

  timer.startEvent("Define grids");
  this->defineGrids(a_eblgFine, a_eblgCoar);
  timer.stopEvent("Define grids");

  timer.startEvent("Define ghost regions");
  this->defineGhostRegions();
  timer.stopEvent("Define ghost regions");

  timer.startEvent("Define buffers");
  this->defineBuffers();
  timer.stopEvent("Define buffers");

  timer.startEvent("Define stencils");
  this->defineStencilsEBCF();
  timer.stopEvent("Define stencils");

  timer.startEvent("Set coarsening ratio");
  if (m_eblgFine.getMaxCoarseningRatio() < m_refRat) {
    m_eblgFine.setMaxCoarseningRatio(m_refRat, m_eblgFine.getEBIS());
  }
  timer.stopEvent("Set coarsening ratio");

  ParmParse pp("EBLeastSquaresMultigridInterpolator");
  bool      profile = false;
  pp.query("profile", profile);
  if (profile) {
    timer.eventReport(pout());
  }
}

int
EBLeastSquaresMultigridInterpolator::getGhostCF() const
{
  return m_ghostCF;
}

EBLeastSquaresMultigridInterpolator::~EBLeastSquaresMultigridInterpolator()
{
  CH_TIME("EBLeastSquaresMultigridInterpolator::~EBLeastSquaresMultigridInterpolator");
}

void
EBLeastSquaresMultigridInterpolator::coarseFineInterp(LevelData<EBCellFAB>&       a_phiFine,
                                                      const LevelData<EBCellFAB>& a_phiCoar,
                                                      const Interval              a_variables)
{
  CH_TIME("EBLeastSquaresMultigridInterpolator::coarseFineInterp");

  CH_assert(m_ghostCF * IntVect::Unit <= a_phiFine.ghostVect());

  if (a_phiFine.ghostVect() != m_ghostVector) {
    MayDay::Error("EBLeastSquaresMultigridInterpolator::coarseFineInterp -- number of ghost cells do not match!");
  }

  // TLDR: This routine does the inhomogeneous coarse-fine interpolation, i.e. the coarse data is not set to zero.

  // Do the regular interpolation -- let QuadCFInterp take care of that.
  LevelData<FArrayBox> fineAlias;
  LevelData<FArrayBox> coarAlias;

  aliasEB(fineAlias, (LevelData<EBCellFAB>&)a_phiFine);
  aliasEB(coarAlias, (LevelData<EBCellFAB>&)a_phiCoar);

  //  QuadCFInterp::coarseFineInterp(fineAlias, coarAlias);
  RefCountedPtr<QuadCFInterp> interpolator;
  if (a_variables == Interval(0, 0)) {
    interpolator = m_scalarInterpolator;
  }
  else if (a_variables == Interval(0, SpaceDim - 1)) {
    interpolator = m_vectorInterpolator;
  }
  interpolator->coarseFineInterp(fineAlias, coarAlias);

  // Interpolate all variables near the EB. We will copy a_phiCoar to m_grownCoarData which holds the data on the coarse grid cells around each fine-grid
  // patch. Note that m_grownCoarData provides a LOCAL view of the coarse grid around each fine-level patch, so we can apply the stencils directly.
  for (int icomp = a_variables.begin(); icomp <= a_variables.end(); icomp++) {

    // Copy data to scratch data holders. We need the coarse-data to be accessible by the fine data (through the same iterators), so we can't
    // get away from a copy here (I think).
    const Interval srcInterv = Interval(icomp, icomp);
    const Interval dstInterv = Interval(m_comp, m_comp);
    a_phiCoar.copyTo(srcInterv, m_grownCoarData, dstInterv);

    // Go through each grid patch and the to-be-interpolated ghost cells across the refinement boundary. We simply
    // apply the stencils here.
    for (DataIterator dit = a_phiFine.dataIterator(); dit.ok(); ++dit) {
      EBCellFAB&       dstFine = a_phiFine[dit()];
      const EBCellFAB& srcFine = a_phiFine[dit()];
      const EBCellFAB& srcCoar = m_grownCoarData[dit()];

      // Apply the coarse stencil
      constexpr int numComp = 1;
      m_aggCoarStencils[dit()]->apply(dstFine,
                                      srcCoar,
                                      m_comp,
                                      icomp,
                                      numComp,
                                      false); // true/false => increment/not increment
      m_aggFineStencils[dit()]->apply(dstFine,
                                      srcFine,
                                      icomp,
                                      icomp,
                                      numComp,
                                      true); // true/false => increment/not increment
    }
  }
}

void
EBLeastSquaresMultigridInterpolator::coarseFineInterpH(LevelData<EBCellFAB>& a_phiFine,
                                                       const Interval        a_variables) const
{
  CH_TIME("EBLeastSquaresMultigridInterpolator::coarseFineInterpH(LD<EBCellFAB>, Interval)");

  CH_assert(m_ghostCF * IntVect::Unit <= a_phiFine.ghostVect());

  if (a_phiFine.ghostVect() != m_ghostVector) {
    MayDay::Error("EBLeastSquaresMultigridInterpolator::coarseFineInterp -- number of ghost cells do not match!");
  }

  // TLDR: This routine does the coarse-fine interpolation with the coarse-grid data set to zero.
  for (DataIterator dit(m_eblgFine.getDBL()); dit.ok(); ++dit) {
    this->coarseFineInterpH(a_phiFine[dit()], a_variables, dit());
  }
}

void
EBLeastSquaresMultigridInterpolator::coarseFineInterpH(EBCellFAB&       a_phi,
                                                       const Interval   a_variables,
                                                       const DataIndex& a_dit) const
{
  CH_TIME("EBLeastSquaresMultigridInterpolator::coarseFineInterpH(LD<EBCellFAB>, Interval, DataIndex)");

  // TLDR: This routine does the coarse-fine interpolation with the coarse-grid data set to zero. This is the kernel version,
  //       operating on a grid patch. It first does a direct kernel for regular data, and then does the interpolation near the
  //       EB after that.
  const Real dxFine = 1.0;
  const Real dxCoar = 1.0 * m_refRat;

  for (int ivar = a_variables.begin(); ivar <= a_variables.end(); ivar++) {

    // Do the regular interp on all sides of the current patch.
    for (int dir = 0; dir < SpaceDim; dir++) {
      for (SideIterator sit; sit.ok(); ++sit) {

        // Regular homogeneous interpolation on side sit() in direction dir
        const Box ghostBox = m_cfivs[a_dit].at(std::make_pair(dir, sit()));

        if (!ghostBox.isEmpty()) {
          BaseFab<Real>& phiReg = a_phi.getSingleValuedFAB();

          // C++ kernel for homogeneous interpolation along a line, assumning that coarse-grid
          // data is zero.
          const Real    c1    = 2 * (dxCoar - dxFine) / (dxCoar + dxFine);
          const Real    c2    = -(dxCoar - dxFine) / (dxCoar + 3 * dxFine);
          const IntVect shift = sign(sit()) * BASISV(dir);

          auto interpHomo = [&](const IntVect& iv) -> void {
            phiReg(iv, ivar) = c1 * phiReg(iv - shift, ivar) + c2 * phiReg(iv - 2 * shift, ivar);
          };

          // Apply the kernel.
          BoxLoops::loop(ghostBox, interpHomo);
        }
      }
    }

    // Apply fine stencil near the EB. It might look weird that we apply the stencil to a_phi AND put the result in a_phi. This
    // is because the stencils are defined in the ghost cells we will fill, but the stencils only reach into valid data. So, we
    // are, in fact, not writing to data that is used by the other stencils
    constexpr int numComp = 1;
    m_aggFineStencils[a_dit]->apply(a_phi, a_phi, ivar, ivar, numComp, false);
  }
}

void
EBLeastSquaresMultigridInterpolator::defineGhostRegions()
{
  CH_TIME("EBLeastSquaresMultigridInterpolator::defineGhostRegions");

  const DisjointBoxLayout& dbl    = m_eblgFine.getDBL();
  const ProblemDomain&     domain = m_eblgFine.getDomain();
  const EBISLayout&        ebisl  = m_eblgFine.getEBISL();

  // Define the "regular" ghost interpolation regions. This is just one cell wide since the operator stencil
  // has a width of 1 in regular cells.
  m_cfivs.define(dbl);
  for (DataIterator dit(dbl); dit.ok(); ++dit) {
    const Box cellBox = dbl[dit()];

    std::map<std::pair<int, Side::LoHiSide>, Box>& cfivsBoxes = m_cfivs[dit()];

    for (int dir = 0; dir < SpaceDim; dir++) {
      for (SideIterator sit; sit.ok(); ++sit) {

        IntVectSet cfivs = IntVectSet(adjCellBox(cellBox, dir, sit(), 1));

        NeighborIterator nit(dbl); // Subtract the other boxes if they intersect this box.
        for (nit.begin(dit()); nit.ok(); ++nit) {
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
  for (DataIterator dit(dbl); dit.ok(); ++dit) {
    const Box      cellBox = dbl[dit()];
    const EBISBox& ebisbox = ebisl[dit()];

    if (ebisbox.isAllRegular() || ebisbox.isAllCovered()) {
      m_ghostCells[dit()] = IntVectSet();
    }
    else {
      // 1. Define the width of the ghost layer region around current (irregular grid) patch
      Box grownBox = grow(cellBox, m_ghostCF);
      grownBox &= domain;

      m_ghostCells[dit()] = IntVectSet(grownBox);

      NeighborIterator nit(dbl);
      for (nit.begin(dit()); nit.ok(); ++nit) {
        m_ghostCells[dit()] -= dbl[nit()];
      }
      m_ghostCells[dit()] -= cellBox;

      // 2. Only include ghost cells that are within range m_ghostCF of an irregular grid cell
      IntVectSet irreg = ebisbox.getIrregIVS(cellBox);
      irreg.grow(2 * m_refRat);
      m_ghostCells[dit()] &= irreg;
    }
  }
}

void
EBLeastSquaresMultigridInterpolator::defineGrids(const EBLevelGrid& a_eblgFine, const EBLevelGrid& a_eblgCoar)
{
  CH_TIME("EBLeastSquaresMultigridInterpolator::defineGrids");

  // Define the fine level grids. Use a shallow copy if we can. Otherwise, if the user has asked for too many ghost cells to be filled
  // then we need to fill in the data again.
  const int ebisGhostFine = a_eblgFine.getGhost();

  if (
    ebisGhostFine >=
    m_ghostVector
      .max()) { // Not enough ghost cells in input eblg. I don't know why anyone would do that, but here's a hook for safety.
    m_eblgFine = a_eblgFine;
  }
  else { // Need to define from scratch
    m_eblgFine = EBLevelGrid(a_eblgFine.getDBL(), a_eblgFine.getDomain(), m_ghostVector.max(), a_eblgFine.getEBIS());
  }

  // Coarse is just a copy.
  m_eblgCoar = a_eblgCoar;

  // Define the coarsened fine grids -- needs same number of ghost cells as m_eblgFine.
  coarsen(m_eblgCoFi, m_eblgFine, m_refRat);
  m_eblgCoFi.setMaxRefinementRatio(m_refRat);
}

void
EBLeastSquaresMultigridInterpolator::defineBuffers()
{
  CH_TIME("EBLeastSquaresMultigridInterpolator::defineBuffers()");

  // TLDR: On both the fine and the coarse level we want to have a BoxLayout which is just like the DisjointBoxLayout except that
  //       the boxes should be grown on each level.

  // This is the number of ghost cells that we will place on the coarse grid. We use it to ensure that we have enough equations.
  const int fineRadius = std::max(2, m_order);
  const int coarRadius = std::max(2, fineRadius / m_refRat);

  const DisjointBoxLayout& coFiGrids = m_eblgCoFi.getDBL();
  const EBISLayout&        coFiEBISL = m_eblgCoFi.getEBISL();

  m_grownCoarData.define(coFiGrids, 1, coarRadius * IntVect::Unit, EBCellFactory(coFiEBISL));
}

void
EBLeastSquaresMultigridInterpolator::defineStencilsEBCF()
{
  CH_TIME("EBLeastSquaresMultigridInterpolator::defineStencilsEBCF");

  // This routine defines stencils for all the ghost cells we need to fill across the EBCF boundary.
  const int comp = 0;

  const Real dxFine = 1.0;
  const Real dxCoar = dxFine * m_refRat;

  const ProblemDomain& domFine = m_eblgFine.getDomain();
  const ProblemDomain& domCoar = m_eblgCoFi.getDomain();

  const DisjointBoxLayout& dblFine = m_eblgFine.getDBL();
  const DisjointBoxLayout& dblCoar = m_eblgCoFi.getDBL();

  const EBISLayout& ebislFine = m_eblgFine.getEBISL();
  const EBISLayout& ebislCoar = m_eblgCoFi.getEBISL();

  m_fineStencils.define(dblFine);
  m_coarStencils.define(dblFine);
  m_vofIterFine.define(dblFine);

  for (DataIterator dit(dblFine); dit.ok(); ++dit) {
    const Box origFineBox    = dblFine[dit()];
    const Box ghostedFineBox = grow(origFineBox, m_ghostVector);
    const Box grownCoarBox   = m_grownCoarData[dit()].box();

    // Define the valid regions such that the interpolation does not include coarse grid cells that fall beneath the fine level,
    // and no fine cells outside the CF.
    DenseIntVectSet validFineCells(
      origFineBox,
      true); // Original patch. We will this patch's ghost cells if those ghost cells overlap with other patches on this level.
    DenseIntVectSet validCoarCells(grownCoarBox, true); //

    // Subtract the coarse cells under the fine box.
    Box coarsenedFineBox = origFineBox;
    coarsenedFineBox.coarsen(m_refRat);
    validCoarCells -= coarsenedFineBox;

    // Same for parts of the current (grown) patch that overlaps with neighboring boxes.
    NeighborIterator nit(dblFine);
    for (nit.begin(dit()); nit.ok(); ++nit) {
      const Box neighborBoxFine = dblFine
        [nit()]; // = neighboring grid patch on the fine level. Use it's cells if they are also ghost cells in the current patch

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

    // Now go through each ghost cell and get an interpolation stencil to specified order.
    const EBISBox& ebisboxFine = m_eblgFine.getEBISL()[dit()];
    const EBISBox& ebisboxCoar = m_eblgCoFi.getEBISL()[dit()];

    const EBGraph& fineGraph = ebisboxFine.getEBGraph();

    m_fineStencils[dit()].define(m_ghostCells[dit()], fineGraph, m_numStenComp);
    m_coarStencils[dit()].define(m_ghostCells[dit()], fineGraph, m_numStenComp);
    m_vofIterFine[dit()].define(m_ghostCells[dit()], fineGraph);

    // Build stencils.
    VoFIterator& vofit = m_vofIterFine[dit()];

    auto kernel = [&](const VolIndex& ghostVofFine) -> void {
      const VolIndex& ghostVofCoar = ebislFine.coarsen(ghostVofFine, m_refRat, dit());

      VoFStencil& fineSten = m_fineStencils[dit()](ghostVofFine, comp);
      VoFStencil& coarSten = m_coarStencils[dit()](ghostVofFine, comp);

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
                                                const int&             a_weight)
{
  CH_TIME("EBLeastSquaresMultigridInterpolator::getStencil");

  // On input, we know which ghost cell we want to interpolate to, and we happen to have a map of valid cells in a_validFineCells and a_validCoarCells. We use
  // that information to build an overdetermined linear system of equations that interpolate to a_ghostVofFine to order a_order. If we don't have enough equations,
  // this routine will return false, and will not create stencils.

  //Timer timer("EBLeastSquaresMultigridInterpolator::getStencil");

  bool foundStencil = true;

  // I think these radii are good -- but there's no hard limit here. Increase the radii if you
  // see that the stencil drops order.
  const int fineRadius = std::max(2, a_order);
  const int coarRadius = std::max(2, fineRadius / m_refRat);

  Vector<VolIndex> fineVofs;
  Vector<VolIndex> coarVofs;

  // Get all Vofs in specified radii. Don't use cells that are not in a_validFineCells or in a_validCoarCells.
  //timer.startEvent("Get Vofs");
  fineVofs =
    VofUtils::getVofsInRadius(a_ghostVofFine, a_ebisboxFine, fineRadius, VofUtils::Connectivity::MonotonePath, false);
  coarVofs =
    VofUtils::getVofsInRadius(a_ghostVofCoar, a_ebisboxCoar, coarRadius, VofUtils::Connectivity::MonotonePath, true);
  //timer.stopEvent("Get Vofs");

  //  const int sizeBefore = fineVofs.size();
  //timer.startEvent("include cells");
  VofUtils::includeCells(fineVofs, a_validFineCells);
  VofUtils::includeCells(coarVofs, a_validCoarCells);
  //timer.stopEvent("include cells");

  //  const int sizeAfter = fineVofs.size();

  //  std::cout << sizeBefore << "\t" << sizeAfter << std::endl;

  const int numEquations = coarVofs.size() + fineVofs.size();
  const int numUnknowns  = LeastSquares::getTaylorExpansionSize(a_order);

  if (numEquations >= numUnknowns) { // We have enough equations to get a stencil.
    //timer.startEvent("sort");
    // In many cases we will have WAY too many equations for the specified order. This is particularly true in 3D
    // because the number of coar vofs included in a radius r from the ghost vof can be (1 + 2*r)^3. So for r = 2
    // this = 125 cells, not counting the fine cells. Since singular value decomposition scales like O(n^3), the
    // penalty for including too may equations can be quite severe.
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
    //timer.stopEvent("sort");

    // Build displacement vectors
    Vector<RealVect> fineDisplacements;
    Vector<RealVect> coarDisplacements;

    //timer.startEvent("dist");
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
    //timer.stopEvent("dist");

    // LeastSquares computes all unknown terms in a Taylor expansion up to specified order. We want the 0th order term, i.e. the interpolated value,
    // which in multi-index notation is the term (0,0), i.e. IntVect::Zero. The format of the two-level least squares routine is such that the
    // fine stencil lies on the first index. This can be confusing, but the LeastSquares uses a very compact notation.
    IntVect    interpStenIndex = IntVect::Zero;
    IntVectSet derivs          = IntVectSet(interpStenIndex);
    IntVectSet knownTerms      = IntVectSet();

    //timer.startEvent("stencil");
    std::map<IntVect, std::pair<VoFStencil, VoFStencil>> stencils =
      LeastSquares::computeDualLevelStencils<float>(derivs,
                                                    knownTerms,
                                                    fineVofs,
                                                    coarVofs,
                                                    fineDisplacements,
                                                    coarDisplacements,
                                                    a_weight,
                                                    a_order);

    a_stencilFine = stencils.at(interpStenIndex).first;
    a_stencilCoar = stencils.at(interpStenIndex).second;
    //timer.stopEvent("stencil");

    foundStencil = true;
  }
  else {
    foundStencil = false;
  }

  //timer.eventReport(true);

  return foundStencil;
}

void
EBLeastSquaresMultigridInterpolator::makeAggStencils()
{
  CH_TIME("EBLeastSquaresMultigridInterpolator::makeAggStencils");

  // Make some proxies.
  LevelData<EBCellFAB> phiProxy(m_eblgFine.getDBL(), 1, m_ghostVector, EBCellFactory(m_eblgFine.getEBISL()));

  // Define the fine-grid stencils.
  m_aggFineStencils.define(m_eblgFine.getDBL());
  m_aggCoarStencils.define(m_eblgFine.getDBL());

  for (DataIterator dit(m_eblgFine.getDBL()); dit.ok(); ++dit) {

    Vector<RefCountedPtr<BaseIndex>>   dstBaseIndex;
    Vector<RefCountedPtr<BaseStencil>> dstBaseStencilFine;
    Vector<RefCountedPtr<BaseStencil>> dstBaseStencilCoar;

    VoFIterator& vofit = m_vofIterFine[dit()];

    auto kernel = [&](const VolIndex& vofFine) -> void {
      const VoFStencil& stencilFine = m_fineStencils[dit()](vofFine, m_stenComp);
      const VoFStencil& stencilCoar = m_coarStencils[dit()](vofFine, m_stenComp);

      dstBaseIndex.push_back(RefCountedPtr<BaseIndex>(new VolIndex(vofFine)));
      dstBaseStencilFine.push_back(RefCountedPtr<BaseStencil>(new VoFStencil(stencilFine)));
      dstBaseStencilCoar.push_back(RefCountedPtr<BaseStencil>(new VoFStencil(stencilCoar)));
    };

    BoxLoops::loop(vofit, kernel);

    m_aggFineStencils[dit()] = RefCountedPtr<AggStencil<EBCellFAB, EBCellFAB>>(
      new AggStencil<EBCellFAB, EBCellFAB>(dstBaseIndex, dstBaseStencilFine, phiProxy[dit()], phiProxy[dit()]));

    m_aggCoarStencils[dit()] = RefCountedPtr<AggStencil<EBCellFAB, EBCellFAB>>(
      new AggStencil<EBCellFAB, EBCellFAB>(dstBaseIndex, dstBaseStencilCoar, m_grownCoarData[dit()], phiProxy[dit()]));
  }
}

#include <CD_NamespaceFooter.H>
