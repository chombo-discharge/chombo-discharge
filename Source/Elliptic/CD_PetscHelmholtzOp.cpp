/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_PetscHelmholtzOp.cpp
 * @brief  Implementation of CD_PetscHelmholtzOp.H
 * @author Robert Marskar
 */

#ifdef CH_USE_PETSC

// Std includes
#include <cmath>

// Chombo includes
#include <CH_Timer.H>
#include <MayDay.H>
#include <ParmParse.H>
#include <BoxIterator.H>
#include <EBArith.H>

// Our includes
#include <CD_PetscHelmholtzOp.H>
#include <CD_DomainFluxIFFABFactory.H>
#include <CD_BoxLoops.H>
#include <CD_ParallelOps.H>
#include <CD_NamespaceHeader.H>

PetscHelmholtzOp::PetscHelmholtzOp() noexcept
{
  CH_TIME("PetscHelmholtzOp::PetscHelmholtzOp");

  m_isDefined = false;
  m_verbose   = false;
  m_matrix    = nullptr;
  m_alpha     = 0.0;
  m_beta      = 0.0;

  ParmParse pp("PetscHelmholtzOp");
  pp.query("verbose", m_verbose);
}

PetscHelmholtzOp::~PetscHelmholtzOp() noexcept
{
  CH_TIME("PetscHelmholtzOp::~PetscHelmholtzOp");
  if (m_verbose) {
    pout() << "PetscHelmholtzOp::~PetscHelmholtzOp" << endl;
  }

  this->clear();
}

void
PetscHelmholtzOp::clear() noexcept
{
  CH_TIME("PetscHelmholtzOp::clear");
  if (m_verbose) {
    pout() << "PetscHelmholtzOp::clear" << endl;
  }

  if (m_matrix != nullptr) {
    PetscCallVoid(MatDestroy(&m_matrix));
  }

  m_petscGrid = RefCountedPtr<PetscGrid>();
  m_domainBc  = RefCountedPtr<PetscHelmholtzDomainBC>();
  m_ebBc      = RefCountedPtr<PetscHelmholtzEBBC>();
  m_jumpBc    = RefCountedPtr<PetscHelmholtzJumpBC>();
  m_probLo    = RealVect::Zero;
  m_isDefined = false;

  m_interpolators.resize(0);
  m_dx.resize(0);
  m_aCoef.clear();
  m_bCoef.clear();
  m_bCoefIrreg.clear();
  m_boundaryWeights[0].clear();
  m_boundaryWeights[1].clear();
  m_bCoefFiCo.clear();

  for (int dir = 0; dir < SpaceDim; dir++) {
    for (int iside = 0; iside < 2; iside++) {
      m_domainFluxStencils[dir][iside].resize(0);
      m_fineCFFluxStencils[dir][iside].resize(0);
      m_coarCFFluxStencils[dir][iside].resize(0);
    }
  }
}

void
PetscHelmholtzOp::define(const RefCountedPtr<PetscGrid>&              a_petscGrid,
                         const Vector<MFMultigridInterpolator>&       a_interpolators,
                         const RefCountedPtr<PetscHelmholtzDomainBC>& a_domainBc,
                         const RefCountedPtr<PetscHelmholtzEBBC>&     a_ebBc,
                         const RefCountedPtr<PetscHelmholtzJumpBC>&   a_jumpBc,
                         const MFAMRCellData&                         a_aCoef,
                         const MFAMRFluxData&                         a_bCoef,
                         const MFAMRIVData&                           a_bCoefIrreg,
                         const RealVect&                              a_probLo,
                         const Vector<Real>&                          a_dx,
                         const Real                                   a_alpha,
                         const Real                                   a_beta) noexcept
{
  CH_TIME("PetscHelmholtzOp::define");
  if (m_verbose) {
    pout() << "PetscHelmholtzOp::define" << endl;
  }

  CH_assert(!a_petscGrid.isNull());
  CH_assert(a_petscGrid->isDefined());
  CH_assert(!a_domainBc.isNull());

  const int numLevels = a_petscGrid->getMFLevelGrids().size();

  // Checked in every build: the operator indexes these per level. The coefficients are checked in fillMatrix.
  if (a_dx.size() < numLevels || a_interpolators.size() < numLevels) {
    MayDay::Abort("PetscHelmholtzOp::define -- dx or the interpolators do not cover every level");
  }

  this->clear();

  m_petscGrid     = a_petscGrid;
  m_interpolators = a_interpolators;
  m_domainBc      = a_domainBc;
  m_ebBc          = a_ebBc;
  m_jumpBc        = a_jumpBc;
  m_aCoef         = a_aCoef;
  m_bCoef         = a_bCoef;
  m_bCoefIrreg    = a_bCoefIrreg;
  m_probLo        = a_probLo;
  m_dx            = a_dx;
  m_alpha         = a_alpha;
  m_beta          = a_beta;

  this->checkSupported();

  m_domainBc->define(m_petscGrid, m_probLo, m_dx);

  if (!m_ebBc.isNull()) {
    m_ebBc->define(m_petscGrid, m_interpolators, m_probLo, m_dx);
  }

  if (!m_jumpBc.isNull()) {
    m_jumpBc->define(m_petscGrid, m_interpolators, m_probLo, m_dx);
  }

  this->createMatrix();

  m_isDefined = true;

  this->fillMatrix();
}

void
PetscHelmholtzOp::setAlphaAndBeta(const Real a_alpha, const Real a_beta) noexcept
{
  CH_TIME("PetscHelmholtzOp::setAlphaAndBeta");
  if (m_verbose) {
    pout() << "PetscHelmholtzOp::setAlphaAndBeta" << endl;
  }

  m_alpha = a_alpha;
  m_beta  = a_beta;
}

void
PetscHelmholtzOp::setCoefficients(const MFAMRCellData& a_aCoef,
                                  const MFAMRFluxData& a_bCoef,
                                  const MFAMRIVData&   a_bCoefIrreg) noexcept
{
  CH_TIME("PetscHelmholtzOp::setCoefficients");
  if (m_verbose) {
    pout() << "PetscHelmholtzOp::setCoefficients" << endl;
  }

  m_aCoef      = a_aCoef;
  m_bCoef      = a_bCoef;
  m_bCoefIrreg = a_bCoefIrreg;
}

const Mat&
PetscHelmholtzOp::getMatrix() const noexcept
{
  CH_TIME("PetscHelmholtzOp::getMatrix");
  if (m_verbose) {
    pout() << "PetscHelmholtzOp::getMatrix" << endl;
  }

  CH_assert(m_isDefined);

  return m_matrix;
}

void
PetscHelmholtzOp::checkSupported() const noexcept
{
  CH_TIME("PetscHelmholtzOp::checkSupported");
  if (m_verbose) {
    pout() << "PetscHelmholtzOp::checkSupported" << endl;
  }

  const Vector<RefCountedPtr<MFLevelGrid>>& mflg = m_petscGrid->getMFLevelGrids();

  // ProblemDomain::contains is true across a periodic boundary, so fillMatrix would take such a face for an interior
  // face whose neighbor is a periodic image outside the domain box, which has no PETSc row. Periodicity is a property
  // of the domain, so checking the coarsest level covers all levels.
  const ProblemDomain& domain = mflg[0]->getDomain();

  for (int dir = 0; dir < SpaceDim; dir++) {
    if (domain.isPeriodic(dir)) {
      MayDay::Abort("PetscHelmholtzOp -- periodic domains are not supported for now");
    }
  }

  // Every level above the coarsest needs an interpolator for its coarse-fine ghost cells.
  for (int lvl = 1; lvl < mflg.size(); lvl++) {
    for (int iphase = 0; iphase < m_petscGrid->getNumPhases(); iphase++) {
      if (this->getInterpolator(lvl, iphase).isNull()) {
        MayDay::Abort("PetscHelmholtzOp -- missing coarse-fine interpolator");
      }
    }
  }
}

Real
PetscHelmholtzOp::getBcoef(const IntVect&       a_iv,
                           const int            a_dir,
                           const Side::LoHiSide a_side,
                           const int            a_level,
                           const DataIndex&     a_din,
                           const int            a_phase) const noexcept
{
  CH_TIME("PetscHelmholtzOp::getBcoef");
  if (m_verbose) {
    pout() << "PetscHelmholtzOp::getBcoef" << endl;
  }

  CH_assert(a_dir >= 0 && a_dir < SpaceDim);

  // The face on the low side of cell iv has index iv, the face on the high side has index iv + BASISV(dir).
  const IntVect    face = (a_side == Side::Lo) ? a_iv : a_iv + BASISV(a_dir);
  const FArrayBox& bco  = (*m_bCoef[a_level])[a_din].getPhase(a_phase)[a_dir].getFArrayBox();

  CH_assert(bco.box().contains(face));

  return bco(face, 0);
}

RefCountedPtr<EBMultigridInterpolator>
PetscHelmholtzOp::getInterpolator(const int a_level, const int a_phase) const noexcept
{
  CH_TIME("PetscHelmholtzOp::getInterpolator");
  if (m_verbose) {
    pout() << "PetscHelmholtzOp::getInterpolator" << endl;
  }

  return (a_level > 0) ? m_interpolators[a_level].getInterpolator(a_phase) : RefCountedPtr<EBMultigridInterpolator>();
}

IntVectSet
PetscHelmholtzOp::getCoveredCells(const int a_level, const DataIndex& a_din, const int a_phase) const noexcept
{
  CH_TIME("PetscHelmholtzOp::getCoveredCells");
  if (m_verbose) {
    pout() << "PetscHelmholtzOp::getCoveredCells" << endl;
  }

  IntVectSet coveredCells;

  const Vector<RefCountedPtr<MFLevelGrid>>& mflg    = m_petscGrid->getMFLevelGrids();
  const BaseFab<PetscAMRCell>&              cells   = (*m_petscGrid->getAMRToPetsc()[a_level])[a_din];
  const EBISBox&                            ebisbox = mflg[a_level]->getEBLevelGrid(a_phase).getEBISL()[a_din];
  const Box                                 box     = mflg[a_level]->getGrids()[a_din];

  // Only the cut-face interpolations use these cells, and there are no covered cells on the finest level.
  if (a_level < static_cast<int>(mflg.size()) - 1 && !(ebisbox.isAllRegular())) {
    const Box region = grow(box, 2) & cells.box();

    for (BoxIterator bit(region); bit.ok(); ++bit) {
      if (cells(bit()).isCoveredByFinerGrid()) {
        coveredCells |= bit();
      }
    }
  }

  return coveredCells;
}

bool
PetscHelmholtzOp::isInterfaceCell(const IntVect& a_iv, const int a_level, const DataIndex& a_din) const noexcept
{
  CH_TIME("PetscHelmholtzOp::isInterfaceCell");
  if (m_verbose) {
    pout() << "PetscHelmholtzOp::isInterfaceCell" << endl;
  }

  bool isInterface = false;

  if (m_petscGrid->getNumPhases() == 2) {
    const MFLevelGrid& mflg = *m_petscGrid->getMFLevelGrids()[a_level];

    isInterface = mflg.getEBLevelGrid(0).getEBISL()[a_din].isIrregular(a_iv) &&
                  mflg.getEBLevelGrid(1).getEBISL()[a_din].isIrregular(a_iv);
  }

  return isInterface;
}

VoFStencil
PetscHelmholtzOp::getCutFaceFluxStencil(const FaceIndex&  a_face,
                                        const IntVectSet& a_coveredCells,
                                        const int         a_level,
                                        const DataIndex&  a_din,
                                        const int         a_phase) const noexcept
{
  CH_TIME("PetscHelmholtzOp::getCutFaceFluxStencil");
  if (m_verbose) {
    pout() << "PetscHelmholtzOp::getCutFaceFluxStencil" << endl;
  }

  const EBLevelGrid& eblg    = m_petscGrid->getMFLevelGrids()[a_level]->getEBLevelGrid(a_phase);
  const EBISBox&     ebisbox = eblg.getEBISL()[a_din];
  const EBFaceFAB&   bco     = (*m_bCoef[a_level])[a_din].getPhase(a_phase)[a_face.direction()];
  const Real         dx      = m_dx[a_level];

  const FaceStencil interp = EBArith::getInterpStencil(a_face, a_coveredCells, ebisbox, eblg.getDomain());

  VoFStencil flux;

  for (int i = 0; i < interp.size(); i++) {
    const FaceIndex& iface = interp.face(i);

    // As in EBHelmholtzOp, a domain face in the interpolation contributes nothing.
    if (!(iface.isBoundary())) {
      VoFStencil faceFlux;

      faceFlux.add(iface.getVoF(Side::Hi), 1.0 / dx);
      faceFlux.add(iface.getVoF(Side::Lo), -1.0 / dx);

      faceFlux *= bco(iface, 0) * interp.weight(i);

      flux += faceFlux;
    }
  }

  return flux;
}

PetscInt
PetscHelmholtzOp::getCoarseRow(const IntVect&   a_ivGhost,
                               const int        a_refRat,
                               const int        a_level,
                               const DataIndex& a_din,
                               const int        a_phase) const noexcept
{
  CH_TIME("PetscHelmholtzOp::getCoarseRow");
  if (m_verbose) {
    pout() << "PetscHelmholtzOp::getCoarseRow" << endl;
  }

  CH_assert(a_level > 0);

  // The coarse row is found through the coarsened fine grids, which share this grid index.
  const IntVect                ivCoar    = coarsen(a_ivGhost, a_refRat);
  const BaseFab<PetscAMRCell>& cellsCoFi = (*m_petscGrid->getAMRToPetscCoFi()[a_level - 1])[a_din];

  if (!cellsCoFi.box().contains(ivCoar)) {
    MayDay::Abort("PetscHelmholtzOp::getCoarseRow -- the coarse cell across a refinement boundary is outside the "
                  "coarsened grid map");
  }

  const PetscInt rowCoar = cellsCoFi(ivCoar).getPetscRow(a_phase);

  if (rowCoar < 0) {
    MayDay::Abort("PetscHelmholtzOp::getCoarseRow -- the coarse cell across a refinement boundary is not a row");
  }

  return rowCoar;
}

void
PetscHelmholtzOp::createMatrix() noexcept
{
  CH_TIME("PetscHelmholtzOp::createMatrix");
  if (m_verbose) {
    pout() << "PetscHelmholtzOp::createMatrix" << endl;
  }

  const Vector<RefCountedPtr<MFLevelGrid>>& mflg = m_petscGrid->getMFLevelGrids();

  const int numLevels = mflg.size();
  const int numPhases = m_petscGrid->getNumPhases();

  // On one level a row has the cell and at most one neighbor per face; the Robin BC references the interior neighbor,
  // which is a neighbor already. Rows next to a refinement boundary also hold the coarse-fine interpolation stencils,
  // whose size is not known in advance, so on several levels the preallocation is an estimate and PETSc may have to
  // allocate more. Neighbors in other boxes can be on other ranks.
  //
  // The rows of cut cells, and of the cells next to them, also hold the EB face stencils and the interpolations of the
  // cut-face fluxes to their centroids, whose sizes are not known in advance either. With cut cells the preallocation
  // is the same estimate as on several levels.
  int hasCutCells = 0;

  for (int phase = 0; phase < numPhases; phase++) {
    for (int lvl = 0; lvl < numLevels; lvl++) {
      const EBLevelGrid&       eblgLvl = mflg[lvl]->getEBLevelGrid(phase);
      const DisjointBoxLayout& dbl     = eblgLvl.getDBL();
      const DataIterator&      dit     = dbl.dataIterator();
      const int                nbox    = dit.size();

#pragma omp parallel for schedule(runtime) reduction(max : hasCutCells)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din = dit[mybox];

        if (!(eblgLvl.getEBISL()[din].isAllRegular())) {
          hasCutCells = 1;
        }
      }
    }
  }

  const bool isEstimate = (numLevels > 1) || (ParallelOps::max(hasCutCells) > 0);

  const PetscInt numLocalRows   = m_petscGrid->getNumLocalRows();
  const PetscInt numLevelPoints = 2 * SpaceDim + 1;
  const PetscInt numDiagonal    = isEstimate ? 4 * numLevelPoints : numLevelPoints;
  const PetscInt numOffDiagonal = numDiagonal - 1;

  if (m_matrix != nullptr) {
    PetscCallVoid(MatDestroy(&m_matrix));
  }

  PetscCallVoid(MatCreateAIJ(Chombo_MPI::comm,
                             numLocalRows,
                             numLocalRows,
                             PETSC_DETERMINE,
                             PETSC_DETERMINE,
                             numDiagonal,
                             nullptr,
                             numOffDiagonal,
                             nullptr,
                             &m_matrix));

  // On one level without cut cells the preallocation is exact, so fail loudly if it is exceeded.
  PetscCallVoid(MatSetOption(m_matrix, MAT_NEW_NONZERO_ALLOCATION_ERR, isEstimate ? PETSC_FALSE : PETSC_TRUE));

  // Storage for the weights of the boundary data on the domain faces, per phase.
  for (int phase = 0; phase < numPhases; phase++) {
    m_boundaryWeights[phase].resize(numLevels);

    for (int lvl = 0; lvl < numLevels; lvl++) {
      const EBLevelGrid& eblgLvl = mflg[lvl]->getEBLevelGrid(phase);

      m_boundaryWeights[phase][lvl] = RefCountedPtr<LevelData<DomainFluxIFFAB>>(
        new LevelData<DomainFluxIFFAB>(eblgLvl.getDBL(),
                                       1,
                                       IntVect::Zero,
                                       DomainFluxIFFABFactory(eblgLvl.getEBISL(), eblgLvl.getDomain())));
    }
  }
}

void
PetscHelmholtzOp::addToMatrix(const PetscInt a_row, const PetscStencil& a_stencil) noexcept
{
  CH_TIME("PetscHelmholtzOp::addToMatrix");
  if (m_verbose) {
    pout() << "PetscHelmholtzOp::addToMatrix" << endl;
  }

  CH_assert(a_row >= 0);

  const PetscInt numEntries = a_stencil.size();

  if (numEntries > 0) {
    Vector<PetscInt>    cols(numEntries);
    Vector<PetscScalar> vals(numEntries);

    for (int i = 0; i < numEntries; i++) {
      cols[i] = a_stencil.row(i);
      vals[i] = a_stencil.weight(i);
    }

    // MatSetValues is not thread safe. PetscCallAbort, not PetscCallVoid, since a return cannot leave the critical
    // section.
#pragma omp critical(PetscHelmholtzOp_MatSetValues)
    {
      PetscCallAbort(Chombo_MPI::comm, MatSetValues(m_matrix, 1, &a_row, numEntries, &cols[0], &vals[0], ADD_VALUES));
    }
  }
}

void
PetscHelmholtzOp::fillMatrix() noexcept
{
  CH_TIME("PetscHelmholtzOp::fillMatrix");
  if (m_verbose) {
    pout() << "PetscHelmholtzOp::fillMatrix" << endl;
  }

  CH_assert(m_isDefined);

  const Vector<RefCountedPtr<MFLevelGrid>>& mflg = m_petscGrid->getMFLevelGrids();

  const int numLevels = mflg.size();
  const int numPhases = m_petscGrid->getNumPhases();

  // The coefficients are indexed with the grid indices of the PetscGrid, so they must be defined on its grids, e.g. not
  // on another realm's grids or on the grids from before a regrid. Checked in every build.
  if (m_aCoef.size() < numLevels || m_bCoef.size() < numLevels) {
    MayDay::Abort("PetscHelmholtzOp::fillMatrix -- the coefficients do not cover every level");
  }

  for (int lvl = 0; lvl < numLevels; lvl++) {
    const DisjointBoxLayout& dbl = mflg[lvl]->getGrids();

    if (m_aCoef[lvl].isNull() || m_bCoef[lvl].isNull() || !(m_aCoef[lvl]->disjointBoxLayout() == dbl) ||
        !(m_bCoef[lvl]->disjointBoxLayout() == dbl)) {
      MayDay::Abort("PetscHelmholtzOp::fillMatrix -- the coefficients are not defined on the grids of the PetscGrid");
    }
  }

  // Zero the entries of a matrix that has been filled before. A refill keeps the sparsity pattern. The entries are
  // added (ADD_VALUES), since a coarse row also gets the refluxing contributions of the fine faces across it.
  PetscBool isAssembled = PETSC_FALSE;

  PetscCallVoid(MatAssembled(m_matrix, &isAssembled));

  if (isAssembled) {
    PetscCallVoid(MatZeroEntries(m_matrix));
  }

  // The rows of each phase. The interface faces couple them.
  for (int phase = 0; phase < numPhases; phase++) {
    for (int lvl = 0; lvl < numLevels; lvl++) {
      const Real dx = m_dx[lvl];

      const EBLevelGrid&       eblg   = mflg[lvl]->getEBLevelGrid(phase);
      const DisjointBoxLayout& dbl    = eblg.getDBL();
      const ProblemDomain&     domain = eblg.getDomain();

      const LevelData<BaseFab<PetscAMRCell>>& amrToPetsc = *m_petscGrid->getAMRToPetsc()[lvl];

      // Refinement ratio to the coarser level, for the coarse side of the refinement boundary. The coarse cell sees a
      // fine face through an area fraction 1/refRat^(SpaceDim-1) and with its own dx = refRat * dx, hence
      // 1/refRat^SpaceDim.
      const int  refRat      = (lvl > 0) ? m_petscGrid->getRefinementRatios()[lvl - 1] : 1;
      const Real refluxScale = 1.0 / std::pow(refRat, SpaceDim);

      const DataIterator& dit  = dbl.dataIterator();
      const int           nbox = dit.size();

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex&             din             = dit[mybox];
        const Box                    box             = dbl[din];
        const BaseFab<PetscAMRCell>& cells           = amrToPetsc[din];
        const FArrayBox&             aco             = (*m_aCoef[lvl])[din].getPhase(phase).getFArrayBox();
        const EBISBox&               ebisbox         = eblg.getEBISL()[din];
        DomainFluxIFFAB&             boundaryWeights = (*m_boundaryWeights[phase][lvl])[din];

        const RefCountedPtr<EBMultigridInterpolator> interpolator = this->getInterpolator(lvl, phase);

        // Cells that the interpolation of cut-face fluxes to the face centroids must not use.
        const IntVectSet coveredCells = this->getCoveredCells(lvl, din, phase);

        auto kernel = [&](const IntVect& iv) -> void {
          const PetscInt row = cells(iv).getPetscRow(phase);

          // Cells covered by the finer level, or by the EB, are not rows.
          if (row < 0) {
            return;
          }

          const VolIndex vof(iv, 0);
          const bool     isIrregular = ebisbox.isIrregular(iv);

          // Stencils over cells on the next coarser level, this level, and the next finer level, mapped once at the
          // end, plus the PETSc stencils of the faces that touch coarse-fine ghost cells and of the EB face.
          VoFStencil   stencilCoar;
          VoFStencil   stencil;
          VoFStencil   stencilFine;
          PetscStencil rowStencil;

          // kappa * alpha * a. Regular cells have kappa = 1.
          if (isIrregular) {
            stencil.add(vof, m_alpha * ebisbox.volFrac(vof) * aco(iv, 0));
          }
          else {
            stencil.add(vof, m_alpha * aco(iv, 0));
          }

          // beta * div(b * grad(phi)): each face contributes beta/dx * s * area * b * F, where F is the derivative of
          // phi in the +dir direction on the face and s = -1 on the low side and +1 on the high side.
          for (int dir = 0; dir < SpaceDim; dir++) {
            for (SideIterator sit; sit.ok(); ++sit) {
              const Side::LoHiSide side    = sit();
              const Real           s       = sign(side);
              const IntVect        ivNeigh = iv + s * BASISV(dir);

              // A face covered by the EB carries no flux. Multi-valued cells are not rows.
              const Vector<FaceIndex> faces = ebisbox.getFaces(vof, dir, side);

              if (faces.size() == 0) {
                continue;
              }

              if (faces.size() != 1) {
                MayDay::Abort("PetscHelmholtzOp::fillMatrix -- a cell must have at most one face in each direction");
              }

              const FaceIndex& face = faces[0];

              if (!domain.contains(ivNeigh) && isIrregular) {
                // Domain face of a cut cell: the domain BC fluxes of the face and its neighbouring domain face,
                // interpolated to the face centroid, as EBHelmholtzDomainBC does.
                const Real factor = m_beta * s * this->getBcoef(iv, dir, side, lvl, din, phase) *
                                    ebisbox.areaFrac(face) / dx;

                const FaceStencil interp = EBArith::getInterpStencil(face, coveredCells, ebisbox, domain);

                for (int i = 0; i < interp.size(); i++) {
                  const VolIndex ivof = interp.face(i).getVoF(flip(side));

                  VoFStencil bcStencilCoar;
                  VoFStencil bcStencil;
                  VoFStencil bcStencilFine;

                  Real bcWeight = 0.0;

                  m_domainBc->getFaceFluxStencil(bcWeight,
                                                 bcStencilCoar,
                                                 bcStencil,
                                                 bcStencilFine,
                                                 phase,
                                                 lvl,
                                                 din,
                                                 ivof,
                                                 dir,
                                                 side);

                  m_petscGrid->replaceGhostCells(bcStencil, bcStencilCoar, interpolator, lvl, din);

                  bcStencilCoar *= factor * interp.weight(i);
                  bcStencil *= factor * interp.weight(i);
                  bcStencilFine *= factor * interp.weight(i);

                  stencilCoar += bcStencilCoar;
                  stencil += bcStencil;
                  stencilFine += bcStencilFine;
                }

                // The boundary data enter the right-hand side through the same interpolation, see
                // addBoundaryConditionsToRHS.
                boundaryWeights(dir, side)(face, 0) = 0.0;
              }
              else if (!domain.contains(ivNeigh)) {
                const Real factor = m_beta * s * this->getBcoef(iv, dir, side, lvl, din, phase) / dx;

                VoFStencil bcStencilCoar;
                VoFStencil bcStencil;
                VoFStencil bcStencilFine;

                Real bcWeight = 0.0;

                m_domainBc->getFaceFluxStencil(bcWeight,
                                               bcStencilCoar,
                                               bcStencil,
                                               bcStencilFine,
                                               phase,
                                               lvl,
                                               din,
                                               vof,
                                               dir,
                                               side);

                boundaryWeights(dir, side)(face, 0) = bcWeight;

                m_petscGrid->replaceGhostCells(bcStencil, bcStencilCoar, interpolator, lvl, din);

                bcStencilCoar *= factor;
                bcStencil *= factor;
                bcStencilFine *= factor;

                stencilCoar += bcStencilCoar;
                stencil += bcStencil;
                stencilFine += bcStencilFine;
              }
              else if (cells(ivNeigh).isCoveredByFinerGrid()) {
                // The neighbor is covered by the finer level: this face is replaced by the fine faces under it, which
                // add their fluxes to this row from the finer level.
              }
              else if (isIrregular) {
                // Face of a cut cell: the flux at the face centroid, as in EBHelmholtzOp. b is in the flux stencil. The
                // stencil can hold coarse-fine ghost cells, also when the neighbor is a row. The faces of regular cells
                // are whole, with the centroid at the center, so they take the two-point flux below.
                const Real factor = m_beta * s * ebisbox.areaFrac(face) / dx;

                VoFStencil flux = this->getCutFaceFluxStencil(face, coveredCells, lvl, din, phase);
                VoFStencil fluxCoar;

                if (cells(ivNeigh).getPetscRow(phase) >= 0) {
                  m_petscGrid->replaceGhostCells(flux, fluxCoar, interpolator, lvl, din);

                  flux *= factor;
                  fluxCoar *= factor;

                  stencil += flux;
                  stencilCoar += fluxCoar;
                }
                else if (cells(ivNeigh).isGhostCF()) {
                  // Fine side of a refinement boundary: this row gets the flux, and the coarse cell across the face
                  // gets it with the refluxing weight, as for the regular faces below.
                  const PetscStencil petscFlux = m_petscGrid->getPetscStencilWithGhostCells(VoFStencil(),
                                                                                            flux,
                                                                                            VoFStencil(),
                                                                                            interpolator,
                                                                                            lvl,
                                                                                            din,
                                                                                            phase);

                  PetscStencil fineContribution = petscFlux;
                  fineContribution *= factor;
                  rowStencil += fineContribution;

                  PetscStencil coarContribution = petscFlux;
                  coarContribution *= -factor * refluxScale;

                  this->addToMatrix(this->getCoarseRow(ivNeigh, refRat, lvl, din, phase), coarContribution);
                }
                else {
                  // Checked in every build: skipping the face would silently drop it from the operator.
                  MayDay::Abort("PetscHelmholtzOp::fillMatrix -- neighbor cell is not a row, a coarse-fine ghost cell, "
                                "or covered by the finer level");
                }
              }
              else if (cells(ivNeigh).getPetscRow(phase) >= 0) {
                const Real factor = m_beta * s * this->getBcoef(iv, dir, side, lvl, din, phase) / dx;

                // F = s * (phi(ivNeigh) - phi(iv))/dx, with both cells rows on this level.
                stencil.add(VolIndex(ivNeigh, 0), factor * s / dx);
                stencil.add(vof, -factor * s / dx);
              }
              else if (cells(ivNeigh).isGhostCF()) {
                const Real factor = m_beta * s * this->getBcoef(iv, dir, side, lvl, din, phase) / dx;

                // The neighbor is a coarse-fine ghost cell, filled from this level and the coarser level. F is then
                // over both levels, and the same fine face flux also enters the coarse cell on the other side of the
                // face.
                CH_assert(lvl > 0);

                const std::pair<VoFStencil, VoFStencil> ghost = interpolator->getGhostStencil(VolIndex(ivNeigh, 0),
                                                                                              din);

                VoFStencil fluxFine = ghost.first;
                VoFStencil fluxCoar = ghost.second;

                fluxFine.add(vof, -1.0);
                fluxFine *= s / dx;
                fluxCoar *= s / dx;

                const PetscStencil flux = m_petscGrid->getPetscStencilWithCoarserLevel(fluxFine,
                                                                                       fluxCoar,
                                                                                       lvl,
                                                                                       din,
                                                                                       phase);

                PetscStencil fineContribution = flux;
                fineContribution *= factor;
                rowStencil += fineContribution;

                // Refluxing: the coarse cell sees this face on its opposite side, through an area fraction
                // 1/refRat^(SpaceDim-1) and with its own dx, so its weight is -factor/refRat^SpaceDim.
                PetscStencil coarContribution = flux;
                coarContribution *= -factor * refluxScale;

                this->addToMatrix(this->getCoarseRow(ivNeigh, refRat, lvl, din, phase), coarContribution);
              }
              else {
                // Checked in every build: skipping the face would silently drop it from the operator.
                MayDay::Abort(
                  "PetscHelmholtzOp::fillMatrix -- neighbor cell is not a row, a coarse-fine ghost cell, or "
                  "covered by the finer level");
              }
            }
          }

          // The EB face: beta/dx * boundary area * b * F, with F = dphi/dn out of the cell. On an interface cell, which
          // is cut in both phases, the jump condition gives b * F over the rows of both phases.
          if (isIrregular && ebisbox.bndryArea(vof) > 0.0) {
            Real         ebWeight = 0.0;
            PetscStencil ebStencil;

            if (this->isInterfaceCell(iv, lvl, din)) {
              m_jumpBc->getInterfaceFluxStencil(ebWeight,
                                                ebStencil,
                                                phase,
                                                lvl,
                                                din,
                                                iv,
                                                (*m_bCoefIrreg[lvl])[din].getIVFAB(0)(vof, 0),
                                                (*m_bCoefIrreg[lvl])[din].getIVFAB(1)(vof, 0));

              ebStencil *= m_beta * ebisbox.bndryArea(vof) / dx;
            }
            else {
              m_ebBc->getEBFluxStencil(ebWeight, ebStencil, phase, lvl, din, vof);

              ebStencil *= m_beta * (*m_bCoefIrreg[lvl])[din].getIVFAB(phase)(vof, 0) * ebisbox.bndryArea(vof) / dx;
            }

            rowStencil += ebStencil;
          }

          rowStencil += m_petscGrid->getPetscStencil(stencil, lvl, din, phase);

          // Checked in every build: there is no coarser level below level 0, and no finer level above the finest.
          if (stencilCoar.size() > 0) {
            if (lvl == 0) {
              MayDay::Abort("PetscHelmholtzOp::fillMatrix -- a stencil reaches a coarser level from level 0");
            }

            rowStencil += m_petscGrid->getPetscStencilWithCoarserLevel(VoFStencil(), stencilCoar, lvl, din, phase);
          }

          if (stencilFine.size() > 0) {
            if (lvl == numLevels - 1) {
              MayDay::Abort("PetscHelmholtzOp::fillMatrix -- a stencil reaches a finer level from the finest level");
            }

            rowStencil += m_petscGrid->getPetscStencilWithFinerLevel(VoFStencil(), stencilFine, lvl, din, phase);
          }

          rowStencil.compress();

          this->addToMatrix(row, rowStencil);
        };

        BoxLoops::loop(box, kernel);
      }
    }
  }

  PetscCallVoid(MatAssemblyBegin(m_matrix, MAT_FINAL_ASSEMBLY));
  PetscCallVoid(MatAssemblyEnd(m_matrix, MAT_FINAL_ASSEMBLY));
}

void
PetscHelmholtzOp::addBoundaryConditionsToRHS(Vec& a_rhs, const EBAMRIVData& a_sigma) const noexcept
{
  CH_TIME("PetscHelmholtzOp::addBoundaryConditionsToRHS");
  if (m_verbose) {
    pout() << "PetscHelmholtzOp::addBoundaryConditionsToRHS" << endl;
  }

  CH_assert(m_isDefined);

  const Vector<RefCountedPtr<MFLevelGrid>>& mflg = m_petscGrid->getMFLevelGrids();

  const int numLevels = mflg.size();
  const int numPhases = m_petscGrid->getNumPhases();

  // Write through the local array, which threads can share: every row this adds to is a cell of this rank, and a cell
  // is only touched by the thread that has its box.
  PetscInt     rowBegin = -1;
  PetscInt     rowEnd   = -1;
  PetscScalar* rhs      = nullptr;

  PetscCallVoid(VecGetOwnershipRange(a_rhs, &rowBegin, &rowEnd));
  PetscCallVoid(VecGetArray(a_rhs, &rhs));

  for (int phase = 0; phase < numPhases; phase++) {
    for (int lvl = 0; lvl < numLevels; lvl++) {
      const Real dx = m_dx[lvl];

      const EBLevelGrid&       eblg   = mflg[lvl]->getEBLevelGrid(phase);
      const DisjointBoxLayout& dbl    = eblg.getDBL();
      const ProblemDomain&     domain = eblg.getDomain();

      const LevelData<BaseFab<PetscAMRCell>>& amrToPetsc = *m_petscGrid->getAMRToPetsc()[lvl];

      const DataIterator& dit  = dbl.dataIterator();
      const int           nbox = dit.size();

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex&             din             = dit[mybox];
        const Box                    box             = dbl[din];
        const BaseFab<PetscAMRCell>& cells           = amrToPetsc[din];
        const EBISBox&               ebisbox         = eblg.getEBISL()[din];
        const DomainFluxIFFAB&       boundaryWeights = (*m_boundaryWeights[phase][lvl])[din];

        // Cells that the interpolation of cut-face fluxes to the face centroids must not use, as in fillMatrix.
        const IntVectSet coveredCells = this->getCoveredCells(lvl, din, phase);

        for (int dir = 0; dir < SpaceDim; dir++) {
          for (SideIterator sit; sit.ok(); ++sit) {
            const Side::LoHiSide side = sit();
            const Real           s    = sign(side);

            // Skip unless the cell layer just outside the box is outside the domain. If so, shifting that layer one
            // cell back gives the cells in this box that abut the domain side.
            Box bndryBox = adjCellBox(box, dir, side, 1);

            if (domain.contains(bndryBox.smallEnd())) {
              continue;
            }

            bndryBox.shift(dir, -sign(side));

            auto kernel = [&](const IntVect& iv) -> void {
              const PetscInt row = cells(iv).getPetscRow(phase);

              // Cells covered by the finer level are not rows, and have no boundary weights.
              if (row < 0) {
                return;
              }

              CH_assert(row >= rowBegin && row < rowEnd);

              const VolIndex          vof(iv, 0);
              const Vector<FaceIndex> faces = ebisbox.getFaces(vof, dir, side);

              // A domain face covered by the EB carries no flux.
              if (faces.size() == 0) {
                return;
              }

              if (faces.size() != 1) {
                MayDay::Abort(
                  "PetscHelmholtzOp::addBoundaryConditionsToRHS -- a cell must have at most one face on the "
                  "domain side");
              }

              if (ebisbox.isIrregular(iv)) {
                // Domain face of a cut cell: the boundary data of the face and its neighbouring domain face, each with
                // its own w_B, interpolated to the face centroid as in fillMatrix.
                const FaceStencil interp = EBArith::getInterpStencil(faces[0], coveredCells, ebisbox, domain);

                Real weightedData = 0.0;

                for (int i = 0; i < interp.size(); i++) {
                  const VolIndex ivof = interp.face(i).getVoF(flip(side));

                  VoFStencil bcStencilCoar;
                  VoFStencil bcStencil;
                  VoFStencil bcStencilFine;

                  Real bcWeight = 0.0;

                  m_domainBc->getFaceFluxStencil(bcWeight,
                                                 bcStencilCoar,
                                                 bcStencil,
                                                 bcStencilFine,
                                                 phase,
                                                 lvl,
                                                 din,
                                                 ivof,
                                                 dir,
                                                 side);

                  weightedData += interp.weight(i) * bcWeight *
                                  m_domainBc->getBoundaryData(phase, lvl, din, ivof, dir, side);
                }

                const Real area = ebisbox.areaFrac(faces[0]);
                const Real c = m_beta * s * this->getBcoef(iv, dir, side, lvl, din, phase) * area / dx * weightedData;

                rhs[row - rowBegin] -= c;
              }
              else {
                const Real weight = boundaryWeights(dir, side)(faces[0], 0);
                const Real data   = m_domainBc->getBoundaryData(phase, lvl, din, vof, dir, side);
                const Real c      = m_beta * s * this->getBcoef(iv, dir, side, lvl, din, phase) / dx * weight * data;

                rhs[row - rowBegin] -= c;
              }
            };

            BoxLoops::loop(bndryBox, kernel);
          }
        }

        // The EB faces: beta/dx * boundary area * b * w_B * g, as the EB flux in fillMatrix.
        const IntVectSet irregCells = ebisbox.getIrregIVS(box);

        for (VoFIterator vofit(irregCells, ebisbox.getEBGraph()); vofit.ok(); ++vofit) {
          const VolIndex& vof = vofit();
          const PetscInt  row = cells(vof.gridIndex()).getPetscRow(phase);

          if (row >= 0 && ebisbox.bndryArea(vof) > 0.0) {
            CH_assert(row >= rowBegin && row < rowEnd);

            Real         ebWeight = 0.0;
            PetscStencil ebStencil;

            const IntVect& iv = vof.gridIndex();

            if (this->isInterfaceCell(iv, lvl, din)) {
              // The jump flux includes b, see fillMatrix.
              m_jumpBc->getInterfaceFluxStencil(ebWeight,
                                                ebStencil,
                                                phase,
                                                lvl,
                                                din,
                                                iv,
                                                (*m_bCoefIrreg[lvl])[din].getIVFAB(0)(VolIndex(iv, 0), 0),
                                                (*m_bCoefIrreg[lvl])[din].getIVFAB(1)(VolIndex(iv, 0), 0));

              const Real data = m_jumpBc->getJumpData(lvl, din, iv, a_sigma);

              rhs[row - rowBegin] -= m_beta * ebisbox.bndryArea(vof) / dx * ebWeight * data;
            }
            else {
              m_ebBc->getEBFluxStencil(ebWeight, ebStencil, phase, lvl, din, vof);

              const Real bco  = (*m_bCoefIrreg[lvl])[din].getIVFAB(phase)(vof, 0);
              const Real data = m_ebBc->getBoundaryData(phase, lvl, din, vof);

              rhs[row - rowBegin] -= m_beta * bco * ebisbox.bndryArea(vof) / dx * ebWeight * data;
            }
          }
        }
      }
    }
  }

  PetscCallVoid(VecRestoreArray(a_rhs, &rhs));
}

void
PetscHelmholtzOp::defineDomainFluxStencils() noexcept
{
  CH_TIME("PetscHelmholtzOp::defineDomainFluxStencils");
  if (m_verbose) {
    pout() << "PetscHelmholtzOp::defineDomainFluxStencils" << endl;
  }

  MayDay::Abort("PetscHelmholtzOp::defineDomainFluxStencils -- not implemented yet");
}

void
PetscHelmholtzOp::defineCoarseFineFluxStencils() noexcept
{
  CH_TIME("PetscHelmholtzOp::defineCoarseFineFluxStencils");
  if (m_verbose) {
    pout() << "PetscHelmholtzOp::defineCoarseFineFluxStencils" << endl;
  }

  MayDay::Abort("PetscHelmholtzOp::defineCoarseFineFluxStencils -- not implemented yet");
}

void
PetscHelmholtzOp::copyFineBcoefToCoarse() noexcept
{
  CH_TIME("PetscHelmholtzOp::copyFineBcoefToCoarse");
  if (m_verbose) {
    pout() << "PetscHelmholtzOp::copyFineBcoefToCoarse" << endl;
  }

  MayDay::Abort("PetscHelmholtzOp::copyFineBcoefToCoarse -- not implemented yet");
}

#include <CD_NamespaceFooter.H>

#endif
