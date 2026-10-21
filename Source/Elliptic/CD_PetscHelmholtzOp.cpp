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

  if (mflg.size() != 1) {
    MayDay::Abort("PetscHelmholtzOp -- only a single level is supported for now");
  }

  // ProblemDomain::contains is true across a periodic boundary, so fillMatrix would take such a face for an interior
  // face whose neighbor is a periodic image outside the domain box, which has no PETSc row.
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
    const EBLevelGrid& eblg  = mflg[0]->getEBLevelGrid(iphase);
    const EBISLayout&  ebisl = eblg.getEBISL();

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
PetscHelmholtzOp::createMatrix() noexcept
{
  CH_TIME("PetscHelmholtzOp::createMatrix");
  if (m_verbose) {
    pout() << "PetscHelmholtzOp::createMatrix" << endl;
  }

  constexpr int phase = 0;

  const Vector<RefCountedPtr<MFLevelGrid>>& mflg = m_petscGrid->getMFLevelGrids();

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

  // Storage for the weights of the boundary data on the domain faces.
  m_boundaryWeights.resize(mflg.size());

  for (int lvl = 0; lvl < mflg.size(); lvl++) {
    const EBLevelGrid& eblgLvl = mflg[lvl]->getEBLevelGrid(phase);

    m_boundaryWeights[lvl] = RefCountedPtr<LevelData<DomainFluxIFFAB>>(
      new LevelData<DomainFluxIFFAB>(eblgLvl.getDBL(),
                                     1,
                                     IntVect::Zero,
                                     DomainFluxIFFABFactory(eblgLvl.getEBISL(), eblgLvl.getDomain())));
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

  constexpr int level = 0;
  constexpr int phase = 0;

  const Real dx = m_dx[level];

  const Vector<RefCountedPtr<MFLevelGrid>>& mflg = m_petscGrid->getMFLevelGrids();

  const EBLevelGrid&       eblg   = mflg[level]->getEBLevelGrid(phase);
  const DisjointBoxLayout& dbl    = eblg.getDBL();
  const ProblemDomain&     domain = eblg.getDomain();

  const LevelData<BaseFab<PetscAMRCell>>& amrToPetsc = *m_petscGrid->getAMRToPetsc()[level];

  // Zero the entries of a matrix that has been filled before. A refill keeps the sparsity pattern.
  PetscBool isAssembled = PETSC_FALSE;

  PetscCallVoid(MatAssembled(m_matrix, &isAssembled));

  if (isAssembled) {
    PetscCallVoid(MatZeroEntries(m_matrix));
  }

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
    const DataIndex&             din             = dit();
    const Box                    box             = dbl[din];
    const BaseFab<PetscAMRCell>& cells           = amrToPetsc[din];
    const FArrayBox&             aco             = (*m_aCoef[level])[din].getPhase(phase).getFArrayBox();
    const EBISBox&               ebisbox         = eblg.getEBISL()[din];
    DomainFluxIFFAB&             boundaryWeights = (*m_boundaryWeights[level])[din];

    auto kernel = [&](const IntVect& iv) -> void {
      const VolIndex vof(iv, 0);

      // Stencils over cells on the next coarser level, this level, and the next finer level.
      VoFStencil stencilCoar;
      VoFStencil stencil;
      VoFStencil stencilFine;

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
            VoFStencil bcStencilCoar;
            VoFStencil bcStencil;
            VoFStencil bcStencilFine;

            Real bcWeight = 0.0;

            m_domainBc->getFaceFluxStencil(bcWeight,
                                           bcStencilCoar,
                                           bcStencil,
                                           bcStencilFine,
                                           phase,
                                           level,
                                           din,
                                           vof,
                                           dir,
                                           side);

            // Without embedded boundaries the cell has exactly one face on the domain side.
            const Vector<FaceIndex> faces = ebisbox.getFaces(vof, dir, side);

            CH_assert(faces.size() == 1);

            boundaryWeights(dir, side)(faces[0], 0) = bcWeight;

            bcStencilCoar *= factor;
            bcStencil *= factor;
            bcStencilFine *= factor;

            stencilCoar += bcStencilCoar;
            stencil += bcStencil;
            stencilFine += bcStencilFine;
          }
        }
      }

      PetscStencil petscStencil = m_petscGrid->getPetscStencil(stencil, level, din, phase);

      if (stencilCoar.size() > 0) {
        petscStencil += m_petscGrid->getPetscStencilWithCoarserLevel(VoFStencil(), stencilCoar, level, din, phase);
      }

      if (stencilFine.size() > 0) {
        petscStencil += m_petscGrid->getPetscStencilWithFinerLevel(VoFStencil(), stencilFine, level, din, phase);
      }

      petscStencil.compress();

      const PetscInt row = cells(iv).getPetscRow(phase);

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

  const Vector<RefCountedPtr<MFLevelGrid>>& mflg = m_petscGrid->getMFLevelGrids();

  const EBLevelGrid&       eblg   = mflg[level]->getEBLevelGrid(phase);
  const DisjointBoxLayout& dbl    = eblg.getDBL();
  const ProblemDomain&     domain = eblg.getDomain();

  const LevelData<BaseFab<PetscAMRCell>>& amrToPetsc = *m_petscGrid->getAMRToPetsc()[level];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
    const DataIndex&             din             = dit();
    const Box                    box             = dbl[din];
    const BaseFab<PetscAMRCell>& cells           = amrToPetsc[din];
    const EBISBox&               ebisbox         = eblg.getEBISL()[din];
    const DomainFluxIFFAB&       boundaryWeights = (*m_boundaryWeights[level])[din];

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
          const VolIndex          vof(iv, 0);
          const Vector<FaceIndex> faces = ebisbox.getFaces(vof, dir, side);

          CH_assert(faces.size() == 1);

          const Real weight = boundaryWeights(dir, side)(faces[0], 0);
          const Real data   = m_domainBc->getBoundaryData(phase, level, din, vof, dir, side);
          const Real c      = m_beta * s * this->getBcoef(iv, dir, side, level, din, phase) / dx * weight * data;

          const PetscInt    row = cells(iv).getPetscRow(phase);
          const PetscScalar val = -c;

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
