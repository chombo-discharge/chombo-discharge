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

  if (m_matrix != nullptr) {
    PetscCallVoid(MatDestroy(&m_matrix));
  }
}

void
PetscHelmholtzOp::define(const RefCountedPtr<PetscGrid>&                       a_petscGrid,
                         const Vector<RefCountedPtr<EBMultigridInterpolator>>& a_interpolators,
                         const RefCountedPtr<PetscHelmholtzDomainBC>&          a_domainBc,
                         const MFAMRCellData&                                  a_aCoef,
                         const MFAMRFluxData&                                  a_bCoef,
                         const RealVect&                                       a_probLo,
                         const Vector<Real>&                                   a_dx,
                         const Real                                            a_alpha,
                         const Real                                            a_beta) noexcept
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

  m_petscGrid     = a_petscGrid;
  m_interpolators = a_interpolators;
  m_domainBc      = a_domainBc;
  m_aCoef         = a_aCoef;
  m_bCoef         = a_bCoef;
  m_probLo        = a_probLo;
  m_dx            = a_dx;
  m_alpha         = a_alpha;
  m_beta          = a_beta;

  this->checkSupported();

  m_domainBc->define(m_petscGrid, m_probLo, m_dx);

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
PetscHelmholtzOp::setCoefficients(const MFAMRCellData& a_aCoef, const MFAMRFluxData& a_bCoef) noexcept
{
  CH_TIME("PetscHelmholtzOp::setCoefficients");
  if (m_verbose) {
    pout() << "PetscHelmholtzOp::setCoefficients" << endl;
  }

  m_aCoef = a_aCoef;
  m_bCoef = a_bCoef;
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
    if (m_interpolators[lvl].isNull()) {
      MayDay::Abort("PetscHelmholtzOp -- missing coarse-fine interpolator");
    }
  }

  // No embedded boundaries: the gas phase must be regular everywhere, and any other phase must have no cells (which is
  // what a geometry without dielectrics gives).
  int hasEB = 0;

  for (int lvl = 0; lvl < mflg.size(); lvl++) {
    const DisjointBoxLayout& dbl = mflg[lvl]->getGrids();

    for (int iphase = 0; iphase < m_petscGrid->getNumPhases(); iphase++) {
      const EBLevelGrid&  eblg  = mflg[lvl]->getEBLevelGrid(iphase);
      const EBISLayout&   ebisl = eblg.getEBISL();
      const DataIterator& dit   = dbl.dataIterator();
      const int           nbox  = dit.size();

#pragma omp parallel for schedule(runtime) reduction(max : hasEB)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din     = dit[mybox];
        const Box&       box     = dbl[din];
        const EBISBox&   ebisbox = ebisl[din];

        const bool isSupported = (iphase == 0) ? ebisbox.isAllRegular() : ebisbox.isCovered(box);

        if (!isSupported) {
          hasEB = 1;
        }
      }
    }
  }

  if (ParallelOps::max(hasEB) > 0) {
    MayDay::Abort("PetscHelmholtzOp -- embedded boundaries and multiple phases are not supported for now");
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
  CH_assert(a_dir >= 0 && a_dir < SpaceDim);

  // The face on the low side of cell iv has index iv, the face on the high side has index iv + BASISV(dir).
  const IntVect    face = (a_side == Side::Lo) ? a_iv : a_iv + BASISV(a_dir);
  const FArrayBox& bco  = (*m_bCoef[a_level])[a_din].getPhase(a_phase)[a_dir].getFArrayBox();

  CH_assert(bco.box().contains(face));

  return bco(face, 0);
}

void
PetscHelmholtzOp::createMatrix() noexcept
{
  CH_TIME("PetscHelmholtzOp::createMatrix");
  if (m_verbose) {
    pout() << "PetscHelmholtzOp::createMatrix" << endl;
  }

  constexpr int phase = 0;

  const Vector<RefCountedPtr<MFLevelGrid>>& mflg = m_petscGrid->getMFLevelGrids();

  const int numLevels = mflg.size();

  // On one level a row has the cell and at most one neighbor per face; the Robin BC references the interior neighbor,
  // which is a neighbor already. Rows next to a refinement boundary also hold the coarse-fine interpolation stencils,
  // whose size is not known in advance, so on several levels the preallocation is an estimate and PETSc may have to
  // allocate more. Neighbors in other boxes can be on other ranks.
  const PetscInt numLocalRows   = m_petscGrid->getNumLocalRows();
  const PetscInt numLevelPoints = 2 * SpaceDim + 1;
  const PetscInt numDiagonal    = (numLevels == 1) ? numLevelPoints : 4 * numLevelPoints;
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

  // On one level the preallocation is exact, so fail loudly if it is exceeded.
  PetscCallVoid(MatSetOption(m_matrix, MAT_NEW_NONZERO_ALLOCATION_ERR, (numLevels == 1) ? PETSC_TRUE : PETSC_FALSE));

  // Storage for the weights of the boundary data on the domain faces.
  m_boundaryWeights.resize(numLevels);

  for (int lvl = 0; lvl < numLevels; lvl++) {
    const EBLevelGrid& eblgLvl = mflg[lvl]->getEBLevelGrid(phase);

    m_boundaryWeights[lvl] = RefCountedPtr<LevelData<DomainFluxIFFAB>>(
      new LevelData<DomainFluxIFFAB>(eblgLvl.getDBL(),
                                     1,
                                     IntVect::Zero,
                                     DomainFluxIFFABFactory(eblgLvl.getEBISL(), eblgLvl.getDomain())));
  }
}

void
PetscHelmholtzOp::addToMatrix(const PetscInt a_row, const PetscStencil& a_stencil) noexcept
{
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
PetscHelmholtzOp::replaceGhostCells(VoFStencil&      a_stencil,
                                    VoFStencil&      a_stencilCoar,
                                    const int        a_level,
                                    const DataIndex& a_din) const noexcept
{
  const BaseFab<PetscAMRCell>& cells = (*m_petscGrid->getAMRToPetsc()[a_level])[a_din];

  VoFStencil kept;

  for (int i = 0; i < a_stencil.size(); i++) {
    const VolIndex& vof    = a_stencil.vof(i);
    const Real      weight = a_stencil.weight(i);
    const IntVect&  iv     = vof.gridIndex();

    if (cells.box().contains(iv) && cells(iv).isGhostCF()) {
      CH_assert(a_level > 0);

      const std::pair<VoFStencil, VoFStencil> ghost = m_interpolators[a_level]->getGhostStencil(vof, a_din);

      VoFStencil ghostFine = ghost.first;
      VoFStencil ghostCoar = ghost.second;

      ghostFine *= weight;
      ghostCoar *= weight;

      kept += ghostFine;
      a_stencilCoar += ghostCoar;
    }
    else {
      kept.add(vof, weight);
    }
  }

  a_stencil = kept;
}

void
PetscHelmholtzOp::fillMatrix() noexcept
{
  CH_TIME("PetscHelmholtzOp::fillMatrix");
  if (m_verbose) {
    pout() << "PetscHelmholtzOp::fillMatrix" << endl;
  }

  CH_assert(m_isDefined);

  constexpr int phase = 0;

  const Vector<RefCountedPtr<MFLevelGrid>>& mflg = m_petscGrid->getMFLevelGrids();

  const int numLevels = mflg.size();

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

  for (int lvl = 0; lvl < numLevels; lvl++) {
    const Real dx = m_dx[lvl];

    const EBLevelGrid&       eblg   = mflg[lvl]->getEBLevelGrid(phase);
    const DisjointBoxLayout& dbl    = eblg.getDBL();
    const ProblemDomain&     domain = eblg.getDomain();

    const LevelData<BaseFab<PetscAMRCell>>& amrToPetsc = *m_petscGrid->getAMRToPetsc()[lvl];

    // Refinement ratio to the coarser level, for the coarse side of the refinement boundary. The coarse cell sees a
    // fine face through an area fraction 1/refRat^(SpaceDim-1) and with its own dx = refRat * dx, hence
    // 1/refRat^SpaceDim.
    const int  refRat      = (lvl > 0) ? std::lround(m_dx[lvl - 1] / dx) : 1;
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
      DomainFluxIFFAB&             boundaryWeights = (*m_boundaryWeights[lvl])[din];

      auto kernel = [&](const IntVect& iv) -> void {
        const PetscInt row = cells(iv).getPetscRow(phase);

        // Cells covered by the finer level are not rows.
        if (row < 0) {
          return;
        }

        const VolIndex vof(iv, 0);

        // Stencils over cells on the next coarser level, this level, and the next finer level, mapped once at the
        // end, plus the PETSc stencils of the faces that touch coarse-fine ghost cells.
        VoFStencil   stencilCoar;
        VoFStencil   stencil;
        VoFStencil   stencilFine;
        PetscStencil rowStencil;

        stencil.add(vof, m_alpha * aco(iv, 0));

        // beta * div(b * grad(phi)): each face contributes beta/dx * s * b * F, where F is the derivative of phi in
        // the +dir direction on the face and s = -1 on the low side and +1 on the high side.
        for (int dir = 0; dir < SpaceDim; dir++) {
          for (SideIterator sit; sit.ok(); ++sit) {
            const Side::LoHiSide side    = sit();
            const Real           s       = sign(side);
            const Real           factor  = m_beta * s * this->getBcoef(iv, dir, side, lvl, din, phase) / dx;
            const IntVect        ivNeigh = iv + s * BASISV(dir);

            if (!domain.contains(ivNeigh)) {
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

              // Without embedded boundaries the cell has exactly one face on the domain side. Checked in every build:
              // using only the first of several faces would silently drop flux.
              const Vector<FaceIndex> faces = ebisbox.getFaces(vof, dir, side);

              if (faces.size() != 1) {
                MayDay::Abort("PetscHelmholtzOp::fillMatrix -- a cell must have exactly one face on the domain side");
              }

              boundaryWeights(dir, side)(faces[0], 0) = bcWeight;

              this->replaceGhostCells(bcStencil, bcStencilCoar, lvl, din);

              bcStencilCoar *= factor;
              bcStencil *= factor;
              bcStencilFine *= factor;

              stencilCoar += bcStencilCoar;
              stencil += bcStencil;
              stencilFine += bcStencilFine;
            }
            else if (cells(ivNeigh).getPetscRow(phase) >= 0) {
              // F = s * (phi(ivNeigh) - phi(iv))/dx, with both cells rows on this level.
              stencil.add(VolIndex(ivNeigh, 0), factor * s / dx);
              stencil.add(vof, -factor * s / dx);
            }
            else if (cells(ivNeigh).isGhostCF()) {
              // The neighbor is a coarse-fine ghost cell, filled from this level and the coarser level. F is then over
              // both levels, and the same fine face flux also enters the coarse cell on the other side of the face.
              CH_assert(lvl > 0);

              const std::pair<VoFStencil, VoFStencil> ghost = m_interpolators[lvl]->getGhostStencil(
                VolIndex(ivNeigh, 0),
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
              // 1/refRat^(SpaceDim-1) and with its own dx, so its weight is -factor/refRat^SpaceDim. The coarse row is
              // found through the coarsened fine grids, which share this grid index.
              const IntVect                ivCoar    = coarsen(ivNeigh, refRat);
              const BaseFab<PetscAMRCell>& cellsCoFi = (*m_petscGrid->getAMRToPetscCoFi()[lvl - 1])[din];

              if (!cellsCoFi.box().contains(ivCoar)) {
                MayDay::Abort("PetscHelmholtzOp::fillMatrix -- the coarse cell across a refinement boundary is outside "
                              "the coarsened grid map");
              }

              const PetscInt rowCoar = cellsCoFi(ivCoar).getPetscRow(phase);

              if (rowCoar < 0) {
                MayDay::Abort("PetscHelmholtzOp::fillMatrix -- the coarse cell across a refinement boundary is not a "
                              "row");
              }

              PetscStencil coarContribution = flux;
              coarContribution *= -factor * refluxScale;

              this->addToMatrix(rowCoar, coarContribution);
            }
            else if (cells(ivNeigh).isCoveredByFinerGrid()) {
              // The neighbor is covered by the finer level: this face is replaced by the fine faces under it, which
              // add their fluxes to this row from the finer level.
            }
            else {
              // Checked in every build: skipping the face would silently drop it from the operator.
              MayDay::Abort("PetscHelmholtzOp::fillMatrix -- neighbor cell is not a row, a coarse-fine ghost cell, or "
                            "covered by the finer level");
            }
          }
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

  PetscCallVoid(MatAssemblyBegin(m_matrix, MAT_FINAL_ASSEMBLY));
  PetscCallVoid(MatAssemblyEnd(m_matrix, MAT_FINAL_ASSEMBLY));
}

void
PetscHelmholtzOp::addBoundaryConditionsToRHS(Vec& a_rhs) const noexcept
{
  CH_TIME("PetscHelmholtzOp::addBoundaryConditionsToRHS");
  if (m_verbose) {
    pout() << "PetscHelmholtzOp::addBoundaryConditionsToRHS" << endl;
  }

  CH_assert(m_isDefined);

  constexpr int phase = 0;

  const Vector<RefCountedPtr<MFLevelGrid>>& mflg = m_petscGrid->getMFLevelGrids();

  const int numLevels = mflg.size();

  // Write through the local array, which threads can share: every row this adds to is a cell of this rank, and a cell
  // is only touched by the thread that has its box.
  PetscInt     rowBegin = -1;
  PetscInt     rowEnd   = -1;
  PetscScalar* rhs      = nullptr;

  PetscCallVoid(VecGetOwnershipRange(a_rhs, &rowBegin, &rowEnd));
  PetscCallVoid(VecGetArray(a_rhs, &rhs));

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
      const DomainFluxIFFAB&       boundaryWeights = (*m_boundaryWeights[lvl])[din];

      for (int dir = 0; dir < SpaceDim; dir++) {
        for (SideIterator sit; sit.ok(); ++sit) {
          const Side::LoHiSide side = sit();
          const Real           s    = sign(side);

          // Skip unless the cell layer just outside the box is outside the domain. If so, shifting that layer one cell
          // back gives the cells in this box that abut the domain side.
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

            if (faces.size() != 1) {
              MayDay::Abort("PetscHelmholtzOp::addBoundaryConditionsToRHS -- a cell must have exactly one face on the "
                            "domain side");
            }

            const Real weight = boundaryWeights(dir, side)(faces[0], 0);
            const Real data   = m_domainBc->getBoundaryData(phase, lvl, din, vof, dir, side);
            const Real c      = m_beta * s * this->getBcoef(iv, dir, side, lvl, din, phase) / dx * weight * data;

            rhs[row - rowBegin] -= c;
          };

          BoxLoops::loop(bndryBox, kernel);
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
PetscHelmholtzOp::defineEBFluxStencils() noexcept
{
  CH_TIME("PetscHelmholtzOp::defineEBFluxStencils");
  if (m_verbose) {
    pout() << "PetscHelmholtzOp::defineEBFluxStencils" << endl;
  }

  MayDay::Abort("PetscHelmholtzOp::defineEBFluxStencils -- not implemented yet");
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
