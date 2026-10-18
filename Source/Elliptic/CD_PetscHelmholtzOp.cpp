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
PetscHelmholtzOp::define(const RefCountedPtr<PetscGrid>&              a_petscGrid,
                         const RefCountedPtr<PetscHelmholtzDomainBC>& a_domainBc,
                         const MFAMRCellData&                         a_aCoef,
                         const MFAMRFluxData&                         a_bCoef,
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
  CH_assert(a_dx.size() >= a_petscGrid->getMFLevelGrids().size());

  m_petscGrid = a_petscGrid;
  m_domainBc  = a_domainBc;
  m_aCoef     = a_aCoef;
  m_bCoef     = a_bCoef;
  m_probLo    = a_probLo;
  m_dx        = a_dx;
  m_alpha     = a_alpha;
  m_beta      = a_beta;

  this->checkSupported();

  // Single level and single phase for now, see checkSupported.
  constexpr int level = 0;
  constexpr int phase = 0;

  m_domainBc->define(Location::Cell::Center,
                     m_petscGrid->getMFLevelGrids()[level]->getEBLevelGrid(phase),
                     m_probLo,
                     m_dx[level]);

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

  if (mflg.size() != 1) {
    MayDay::Abort("PetscHelmholtzOp -- only a single level is supported for now");
  }

  // Faces across a periodic boundary would be taken for domain faces, see assembleMatrix.
  const ProblemDomain& domain = mflg[0]->getDomain();

  for (int dir = 0; dir < SpaceDim; dir++) {
    if (domain.isPeriodic(dir)) {
      MayDay::Abort("PetscHelmholtzOp -- periodic domains are not supported for now");
    }
  }

  // No embedded boundaries: the gas phase must be regular everywhere, and any other phase must have no cells (which is
  // what a geometry without dielectrics gives).
  const DisjointBoxLayout& dbl = mflg[0]->getGrids();

  int hasEB = 0;

  for (int iphase = 0; iphase < m_petscGrid->getNumPhases(); iphase++) {
    const EBISLayout& ebisl = mflg[0]->getEBLevelGrid(iphase).getEBISL();

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
      const Box&     box     = dbl[dit()];
      const EBISBox& ebisbox = ebisl[dit()];

      const bool isSupported = (iphase == 0) ? ebisbox.isAllRegular() : ebisbox.isCovered(box);

      if (!isSupported) {
        hasEB = 1;
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
PetscHelmholtzOp::assembleMatrix() noexcept
{
  CH_TIME("PetscHelmholtzOp::assembleMatrix");
  if (m_verbose) {
    pout() << "PetscHelmholtzOp::assembleMatrix" << endl;
  }

  constexpr int level = 0;
  constexpr int phase = 0;

  const Real dx = m_dx[level];

  const EBLevelGrid&       eblg   = m_petscGrid->getMFLevelGrids()[level]->getEBLevelGrid(phase);
  const DisjointBoxLayout& dbl    = eblg.getDBL();
  const ProblemDomain&     domain = eblg.getDomain();

  const LevelData<BaseFab<PetscAMRCell>>& amrToPetsc = *m_petscGrid->getAMRToPetsc()[level];

  // A row has the cell and at most one neighbor per face. The Robin BC references the interior neighbor, which is a
  // neighbor already. Neighbors in other boxes can be on other ranks.
  const PetscInt numLocalRows = m_petscGrid->getNumLocalRows();
  const PetscInt maxEntries   = 2 * SpaceDim + 1;

  if (m_matrix != nullptr) {
    PetscCallVoid(MatDestroy(&m_matrix));
  }

  PetscCallVoid(MatCreateAIJ(Chombo_MPI::comm,
                             numLocalRows,
                             numLocalRows,
                             PETSC_DETERMINE,
                             PETSC_DETERMINE,
                             maxEntries,
                             nullptr,
                             maxEntries - 1,
                             nullptr,
                             &m_matrix));

  // Fail loudly if the preallocation is too small, rather than silently reallocating.
  PetscCallVoid(MatSetOption(m_matrix, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE));

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
    const DataIndex&             din   = dit();
    const Box                    box   = dbl[din];
    const BaseFab<PetscAMRCell>& cells = amrToPetsc[din];
    const FArrayBox&             aco   = (*m_aCoef[level])[din].getPhase(phase).getFArrayBox();

    auto kernel = [&](const IntVect& iv) -> void {
      const VolIndex vof(iv, 0);

      VoFStencil stencil;

      stencil.add(vof, m_alpha * aco(iv, 0));

      // beta * div(b * grad(phi)): each face contributes beta/dx * s * b * F, where F is the derivative of phi in the
      // +dir direction on the face and s = -1 on the low side and +1 on the high side.
      for (int dir = 0; dir < SpaceDim; dir++) {
        for (SideIterator sit; sit.ok(); ++sit) {
          const Side::LoHiSide side    = sit();
          const Real           s       = sign(side);
          const Real           factor  = m_beta * s * this->getBcoef(iv, dir, side, level, din, phase) / dx;
          const IntVect        ivNeigh = iv + s * BASISV(dir);

          if (domain.contains(ivNeigh)) {
            // F = s * (phi(ivNeigh) - phi(iv))/dx.
            stencil.add(VolIndex(ivNeigh, 0), factor * s / dx);
            stencil.add(vof, -factor * s / dx);
          }
          else {
            const std::pair<VoFStencil, VoFStencil> faceStencil = m_domainBc->getFaceFluxStencil(vof, dir, side, din);

            // Single level: the finer-level part must be empty.
            CH_assert(faceStencil.second.size() == 0);

            VoFStencil bcStencil = faceStencil.first;
            bcStencil *= factor;

            stencil += bcStencil;
          }
        }
      }

      const PetscStencil petscStencil = m_petscGrid->getPetscStencil(stencil, level, din, phase);
      const PetscInt     row          = cells(iv).getPetscRow(phase);

      CH_assert(row >= 0);

      Vector<PetscInt>    cols(petscStencil.size());
      Vector<PetscScalar> vals(petscStencil.size());

      for (int i = 0; i < petscStencil.size(); i++) {
        cols[i] = petscStencil.row(i);
        vals[i] = petscStencil.weight(i);
      }

      PetscCallVoid(MatSetValues(m_matrix, 1, &row, petscStencil.size(), &cols[0], &vals[0], INSERT_VALUES));
    };

    BoxLoops::loop(box, kernel);
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

  constexpr int level = 0;
  constexpr int phase = 0;

  const Real dx = m_dx[level];

  const EBLevelGrid&       eblg   = m_petscGrid->getMFLevelGrids()[level]->getEBLevelGrid(phase);
  const DisjointBoxLayout& dbl    = eblg.getDBL();
  const ProblemDomain&     domain = eblg.getDomain();

  const LevelData<BaseFab<PetscAMRCell>>& amrToPetsc = *m_petscGrid->getAMRToPetsc()[level];

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
          const VolIndex vof(iv, 0);
          const Real     g = m_beta * s * this->getBcoef(iv, dir, side, level, din, phase) / dx *
                         m_domainBc->getFaceFluxConstant(vof, dir, side, din);
          const PetscInt    row = cells(iv).getPetscRow(phase);
          const PetscScalar val = -g;

          CH_assert(row >= 0);

          PetscCallVoid(VecSetValue(a_rhs, row, val, ADD_VALUES));
        };

        BoxLoops::loop(bndryBox, kernel);
      }
    }
  }

  PetscCallVoid(VecAssemblyBegin(a_rhs));
  PetscCallVoid(VecAssemblyEnd(a_rhs));
}

#include <CD_NamespaceFooter.H>

#endif
