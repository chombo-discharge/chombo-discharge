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

  CH_assert(a_dx.size() >= numLevels);
  CH_assert(a_interpolators.size() >= numLevels);
  CH_assert(a_aCoef.size() >= numLevels);
  CH_assert(a_bCoef.size() >= numLevels);

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

  // Single phase for now, see checkSupported.
  constexpr int phase = 0;

  Vector<EBLevelGrid> eblg(numLevels);

  for (int lvl = 0; lvl < numLevels; lvl++) {
    eblg[lvl] = m_petscGrid->getMFLevelGrids()[lvl]->getEBLevelGrid(phase);
  }

  m_domainBc->define(Location::Cell::Center, eblg, m_probLo, m_dx);

  this->assembleMatrix();

  m_isDefined = true;
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

  // Faces across a periodic boundary would be taken for domain faces, see assembleMatrix. Periodicity is a property of
  // the domain, so checking the coarsest level covers all levels.
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
      const EBISLayout& ebisl = mflg[lvl]->getEBLevelGrid(iphase).getEBISL();

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
        const Box&     box     = dbl[dit()];
        const EBISBox& ebisbox = ebisl[dit()];

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

    PetscCallVoid(MatSetValues(m_matrix, 1, &a_row, numEntries, &cols[0], &vals[0], ADD_VALUES));
  }
}

void
PetscHelmholtzOp::assembleMatrix() noexcept
{
  CH_TIME("PetscHelmholtzOp::assembleMatrix");
  if (m_verbose) {
    pout() << "PetscHelmholtzOp::assembleMatrix" << endl;
  }

  constexpr int phase = 0;

  const int numLevels = m_petscGrid->getMFLevelGrids().size();

  // On one level a row has the cell and at most one neighbor per face; the Robin BC references the interior neighbor,
  // which is a neighbor already. Rows next to a refinement boundary also hold the coarse-fine interpolation stencils,
  // whose size is not known in advance, so on several levels the preallocation is an estimate and PETSc may have to
  // allocate more.
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

  for (int lvl = 0; lvl < numLevels; lvl++) {
    const Real dx = m_dx[lvl];

    const EBLevelGrid&       eblg   = m_petscGrid->getMFLevelGrids()[lvl]->getEBLevelGrid(phase);
    const DisjointBoxLayout& dbl    = eblg.getDBL();
    const ProblemDomain&     domain = eblg.getDomain();

    const LevelData<BaseFab<PetscAMRCell>>& amrToPetsc = *m_petscGrid->getAMRToPetsc()[lvl];

    // Refinement ratio to the coarser level, for the coarse side of the refinement boundary.
    const int refRat = (lvl > 0) ? std::lround(m_dx[lvl - 1] / dx) : 1;

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
      const DataIndex&             din   = dit();
      const Box                    box   = dbl[din];
      const BaseFab<PetscAMRCell>& cells = amrToPetsc[din];
      const FArrayBox&             aco   = (*m_aCoef[lvl])[din].getPhase(phase).getFArrayBox();

      auto kernel = [&](const IntVect& iv) -> void {
        const PetscInt row = cells(iv).getPetscRow(phase);

        // Cells covered by the finer level are not rows.
        if (row < 0) {
          return;
        }

        const VolIndex vof(iv, 0);

        // Stencil over level lvl cells that are rows, mapped once at the end, plus the PETSc stencils of the faces
        // that touch coarse-fine ghost cells.
        VoFStencil   levelStencil;
        PetscStencil rowStencil;

        levelStencil.add(vof, m_alpha * aco(iv, 0));

        // beta * div(b * grad(phi)): each face contributes beta/dx * s * b * F, where F is the derivative of phi in the
        // +dir direction on the face and s = -1 on the low side and +1 on the high side.
        for (int dir = 0; dir < SpaceDim; dir++) {
          for (SideIterator sit; sit.ok(); ++sit) {
            const Side::LoHiSide side    = sit();
            const Real           s       = sign(side);
            const Real           bco     = this->getBcoef(iv, dir, side, lvl, din, phase);
            const Real           factor  = m_beta * s * bco / dx;
            const IntVect        ivNeigh = iv + s * BASISV(dir);

            if (!domain.contains(ivNeigh)) {
              const std::pair<VoFStencil, VoFStencil> faceStencil = m_domainBc->getFaceFluxStencil(vof,
                                                                                                   dir,
                                                                                                   side,
                                                                                                   lvl,
                                                                                                   din);

              // Boundary stencils that reach under the finer level are not supported yet. Checked in every build:
              // dropping the second part would silently change the operator.
              if (faceStencil.second.size() != 0) {
                MayDay::Abort("PetscHelmholtzOp::assembleMatrix -- domain BC stencils on the finer level are not "
                              "supported yet");
              }

              VoFStencil bcStencil = faceStencil.first;
              bcStencil *= factor;

              levelStencil += bcStencil;
            }
            else if (cells(ivNeigh).getPetscRow(phase) >= 0) {
              // F = s * (phi(ivNeigh) - phi(iv))/dx, with both cells rows on this level.
              levelStencil.add(VolIndex(ivNeigh, 0), factor * s / dx);
              levelStencil.add(vof, -factor * s / dx);
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

              CH_assert(cellsCoFi.box().contains(ivCoar));

              const PetscInt rowCoar = cellsCoFi(ivCoar).getPetscRow(phase);

              if (rowCoar < 0) {
                MayDay::Abort("PetscHelmholtzOp::assembleMatrix -- the coarse cell across a refinement boundary is not "
                              "a row");
              }

              PetscStencil coarContribution = flux;
              coarContribution *= -factor / std::pow(refRat, SpaceDim);

              this->addToMatrix(rowCoar, coarContribution);
            }
            else if (cells(ivNeigh).isCoveredByFinerGrid()) {
              // The neighbor is covered by the finer level: this face is replaced by the fine faces under it, which
              // add their fluxes to this row from the finer level.
            }
            else {
              // Checked in every build: skipping the face would silently drop it from the operator.
              MayDay::Abort("PetscHelmholtzOp::assembleMatrix -- neighbor cell is not a row, a coarse-fine ghost "
                            "cell, or covered by the finer level");
            }
          }
        }

        rowStencil += m_petscGrid->getPetscStencil(levelStencil, lvl, din, phase);
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

  const int numLevels = m_petscGrid->getMFLevelGrids().size();

  for (int lvl = 0; lvl < numLevels; lvl++) {
    const Real dx = m_dx[lvl];

    const EBLevelGrid&       eblg   = m_petscGrid->getMFLevelGrids()[lvl]->getEBLevelGrid(phase);
    const DisjointBoxLayout& dbl    = eblg.getDBL();
    const ProblemDomain&     domain = eblg.getDomain();

    const LevelData<BaseFab<PetscAMRCell>>& amrToPetsc = *m_petscGrid->getAMRToPetsc()[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
      const DataIndex&             din   = dit();
      const Box                    box   = dbl[din];
      const BaseFab<PetscAMRCell>& cells = amrToPetsc[din];

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

            // Cells covered by the finer level are not rows.
            if (row < 0) {
              return;
            }

            const VolIndex vof(iv, 0);
            const Real     g = m_beta * s * this->getBcoef(iv, dir, side, lvl, din, phase) / dx *
                           m_domainBc->getFaceFluxConstant(vof, dir, side, lvl, din);
            const PetscScalar val = -g;

            PetscCallVoid(VecSetValue(a_rhs, row, val, ADD_VALUES));
          };

          BoxLoops::loop(bndryBox, kernel);
        }
      }
    }
  }

  PetscCallVoid(VecAssemblyBegin(a_rhs));
  PetscCallVoid(VecAssemblyEnd(a_rhs));
}

#include <CD_NamespaceFooter.H>

#endif
