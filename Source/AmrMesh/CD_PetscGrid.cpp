/*!
  chombo-discharge
  Copyright © 2026 SINTEF Energy Research.
  Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
*/

/*!
  @file   CD_PetscGrid.cpp
  @brief  Implementation of CD_PetscGrid.H
  @author Robert Marskar
*/

#ifdef CH_USE_PETSC

// Chombo includes
#include <CH_Timer.H>
#include <ParmParse.H>

// Our includes
#include <CD_PetscGrid.H>
#include <CD_BoxLoops.H>
#include <CD_ParallelOps.H>
#include <CD_NamespaceHeader.H>

#if 1
#warning "Not sure how I want to handle memory management within PetscGrid. Maybe pass this off to the outside world??"
#warning "We should REALLY time how fast transfers between Chombo and PETSc really are"
#warning "I really want a debug function for writing the PETSc grid to a file, showing all DOFs, etc."
#endif

PetscGrid::PetscGrid() noexcept
{
  CH_TIME("PetscGrid::PetscGrid");

  m_isDefined = false;
  m_verbose   = false;
  m_debug     = false;
  m_profile   = false;

  m_finestLevel = -1;
  m_numGhost    = -1;
  m_numPhases   = -1;

  m_numLocalRows  = -1;
  m_numGlobalRows = -1;
  m_localRowBegin = -1;
}

PetscGrid::~PetscGrid() noexcept
{
  CH_TIME("PetscGrid::~PetscGrid");
}

void
PetscGrid::define(const Vector<RefCountedPtr<MFLevelGrid>>&              a_levelGrids,
                  const Vector<RefCountedPtr<MFLevelGrid>>&              a_levelGridsCoFi,
                  const Vector<RefCountedPtr<MFLevelGrid>>&              a_levelGridsFiCo,
                  const Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>& a_validCells,
                  const Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>& a_coarHaloCF,
                  const Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>& a_fineHaloCF,
                  const int                                              a_finestLevel,
                  const int                                              a_numGhost) noexcept
{
  CH_TIME("PetscGrid::define");
  if (m_verbose) {
    pout() << "PetscGrid::define" << endl;
  }

  ParmParse pp("PetscGrid");

  pp.query("verbose", m_verbose);
  pp.query("debug", m_debug);
  pp.query("profile", m_profile);

  m_levelGrids     = a_levelGrids;
  m_levelGridsCoFi = a_levelGridsCoFi;
  m_levelGridsFiCo = a_levelGridsFiCo;
  m_validCells     = a_validCells;
  m_coarHaloCF     = a_coarHaloCF;
  m_fineHaloCF     = a_fineHaloCF;

  m_finestLevel = a_finestLevel;
  m_numGhost    = a_numGhost;
  m_numPhases   = m_levelGrids[0]->numPhases();

  CH_assert(m_numGhost >= 0);
  CH_assert(m_finestLevel >= 0);

  this->defineAMRCells();
  this->definePetscRows();
  //  this->defineCoFiBuffers();

  m_isDefined = true;
}

void
PetscGrid::clear() noexcept
{
  CH_TIME("PetscGrid::clear");
  if (m_verbose) {
    pout() << "PetscGrid::clear" << endl;
  }

  if (m_isDefined) {
    PetscCallVoid(ISLocalToGlobalMappingDestroy(&m_localToGlobalIS));
  }
}

void
PetscGrid::defineAMRCells() noexcept
{
  CH_TIME("PetscGrid::defineAMRCells");
  if (m_verbose) {
    pout() << "PetscGrid::defineAMRCells" << endl;
  }

  const int curComp = 0;
  const int numComp = 1;

  m_amrToPetsc.resize(1 + m_finestLevel);

  for (int lvl = m_finestLevel; lvl >= 0; lvl--) {
    const ProblemDomain&     domain = m_levelGrids[lvl]->getDomain();
    const DisjointBoxLayout& dbl    = m_levelGrids[lvl]->getGrids();
    const DataIterator&      dit    = dbl.dataIterator();
    const int                nbox   = dit.size();

    m_amrToPetsc[lvl] = RefCountedPtr<LevelData<BaseFab<PetscAMRCell>>>(
      new LevelData<BaseFab<PetscAMRCell>>(dbl, numComp, m_numGhost * IntVect::Unit));

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din            = dit[mybox];
      const Box        cellBox        = dbl[din];
      const Box        ghostedCellBox = grow(cellBox, m_numGhost);

      BaseFab<PetscAMRCell>& amrToPetsc = (*m_amrToPetsc[lvl])[din];
      const BaseFab<bool>&   validCells = (*m_validCells[lvl])[din];

      auto setValidCells = [&](const IntVect& iv) -> void {
        amrToPetsc(iv).setPetscRow(0, -1);
        amrToPetsc(iv).setPetscRow(1, -1);
        amrToPetsc(iv).setDomainBoundaryCell(false);
        amrToPetsc(iv).setGhostCF(true);

        if (cellBox.contains(iv)) {
          amrToPetsc(iv).setGhostCF(false);

          if (validCells(iv)) {
            amrToPetsc(iv).setCoveredByFinerGrid(false);
          }
          else {
            amrToPetsc(iv).setCoveredByFinerGrid(true);
          }
        }
      };

      BoxLoops::loop(ghostedCellBox, setValidCells);

      // Fix CF region flags
      if (lvl < m_finestLevel) {
        const BaseFab<bool>& haloCells = (*m_coarHaloCF[lvl])[din];

        auto setOuterCFRegion = [&](const IntVect& iv) -> void {
          amrToPetsc(iv).setCoarCF(haloCells(iv));
        };

        BoxLoops::loop(cellBox, setOuterCFRegion);
      }

      if (lvl > 0) {
        const BaseFab<bool>& haloCells = (*m_fineHaloCF[lvl])[din];

        auto setInnerCFRegion = [&](const IntVect& iv) -> void {
          amrToPetsc(iv).setFineCF(haloCells(iv));
        };

        BoxLoops::loop(cellBox, setInnerCFRegion);
      }

      // Fix up domain side flags
      for (int dir = 0; dir < SpaceDim; dir++) {
        const Box domainBoxLo = adjCellLo(domain.domainBox(), dir, -1);
        const Box domainBoxHi = adjCellHi(domain.domainBox(), dir, -1);

        auto setDomainFlag = [&](const IntVect& iv) -> void {
          if (domainBoxLo.contains(iv) || domainBoxHi.contains(iv)) {
            amrToPetsc(iv).setDomainBoundaryCell(true);
          }
        };

        BoxLoops::loop(cellBox, setDomainFlag);
      }
    }

    m_amrToPetsc[lvl]->exchange();
  }
}

void
PetscGrid::definePetscRows() noexcept
{
  CH_TIME("PetscGrid::definePetscRows");
  if (m_verbose) {
    pout() << "PetscGrid::definePetscRows" << endl;
  }

  CH_assert(!(m_levelGrids[0].isNull()));

  m_numLocalRows  = 0;
  m_numGlobalRows = 0;

  m_numRowsPerRank.resize(numProc(), 0);

  for (int iphase = 0; iphase < m_numPhases; iphase++) {
    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      const EBLevelGrid&       eblg  = m_levelGrids[lvl]->getEBLevelGrid(iphase);
      const EBISLayout&        ebisl = eblg.getEBISL();
      const DisjointBoxLayout& dbl   = eblg.getDBL();
      const DataIterator&      dit   = dbl.dataIterator();
      const int                nbox  = dit.size();

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex din = dit[mybox];

        const Box&        cellBox = dbl[din];
        const EBISBox&    ebisBox = ebisl[din];
        const IntVectSet& ivs     = ebisBox.getIrregIVS(cellBox);
        const EBGraph&    ebgraph = ebisBox.getEBGraph();

        BaseFab<PetscAMRCell>& amrToPetsc = (*m_amrToPetsc[lvl])[din];

        auto regularKernel = [&](const IntVect& iv) -> void {
          if (!(amrToPetsc(iv).isCoveredByFinerGrid()) && !ebisBox.isCovered(iv)) {
            if (ebisBox.isMultiValued(iv)) {
              MayDay::Abort("PetscGrid::definePetscRows -- multi-valued cells are not permitted");
            }

            PetscDOF dof;

            dof.gridLevel = lvl;
            dof.gridIndex = din;
            dof.gridCell  = iv;
            dof.phase     = iphase;

            m_petscToAMR.push_back(dof);

            amrToPetsc(iv).setPetscRow(iphase, m_numLocalRows);

            m_numLocalRows = m_numLocalRows + 1;
          }
        };

        BoxLoops::loop(cellBox, regularKernel);
      }
    }
  }

  // Gather the number of rows per MPI rank. And figure out the total number of rows.
  m_numRowsPerRank = ParallelOps::gather(m_numLocalRows);
  m_numGlobalRows  = ParallelOps::sum(m_numLocalRows);

  // Compute the global row offset for this rank.
  m_localRowBegin = 0;
  for (PetscInt i = 0; i < procID(); i++) {
    m_localRowBegin += m_numRowsPerRank[i];
  }

  // Build local to global mapping for PETSc, so that we can extract local subvectors.
  PetscInt* gidx;
  PetscCallVoid(PetscMalloc1(m_numLocalRows, &gidx));
  for (PetscInt i = 0; i < m_numLocalRows; i++) {
    gidx[i] = m_localRowBegin + i;
  }

  PetscCallVoid(
    ISLocalToGlobalMappingCreate(Chombo_MPI::comm, 1, m_numLocalRows, gidx, PETSC_OWN_POINTER, &m_localToGlobalIS));

  if (m_debug) {
    pout() << "PetscGrid::definePetscRows - numLocalUnknowns = " << m_numLocalRows << endl;
    pout() << "PetscGrid::definePetscRows - numGlobalUnknowns = " << m_numGlobalRows << endl;
  }
}

void
PetscGrid::defineCoFiBuffers() noexcept
{
  CH_TIME("PetscGrid::defineCoFiBuffers");
  if (m_verbose) {
    pout() << "PetscGrid::defineCoFiBuffers" << endl;
  }

#warning "Need to begin here next time"

  const int numComp = 1;
  const int comp    = 0;

  m_amrToPetscCoFi.resize(1 + m_finestLevel);
  m_amrToPetscFiCo.resize(1 + m_finestLevel);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const bool hasCoar = lvl > 0;
    const bool hasFine = lvl < m_finestLevel;

    if (hasFine) {
      const DisjointBoxLayout& dbl  = m_levelGridsCoFi[lvl]->getGrids();
      const DataIterator       dit  = dbl.dataIterator();
      const int                nbox = dit.size();

      m_amrToPetscCoFi[lvl] = RefCountedPtr<LevelData<BaseFab<PetscAMRCell>>>(
        new LevelData<BaseFab<PetscAMRCell>>(dbl, numComp, m_numGhost * IntVect::Unit));

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din = dit[mybox];

        BaseFab<PetscAMRCell>& amrToPetsc = (*m_amrToPetscCoFi[lvl])[din];

        amrToPetsc.setVal(PetscAMRCell());
      }

      m_amrToPetsc[lvl]->copyTo(*m_amrToPetscCoFi[lvl]);
      m_amrToPetscCoFi[lvl]->exchange();
    }

    if (hasCoar) {
      const DisjointBoxLayout& dbl  = m_levelGridsFiCo[lvl]->getGrids();
      const DataIterator       dit  = dbl.dataIterator();
      const int                nbox = dit.size();

      m_amrToPetscFiCo[lvl] = RefCountedPtr<LevelData<BaseFab<PetscAMRCell>>>(
        new LevelData<BaseFab<PetscAMRCell>>(dbl, numComp, m_numGhost * IntVect::Unit));

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din = dit[mybox];

        BaseFab<PetscAMRCell>& amrToPetsc = (*m_amrToPetscFiCo[lvl])[din];

        amrToPetsc.setVal(PetscAMRCell());
      }

      m_amrToPetsc[lvl]->copyTo(*m_amrToPetscFiCo[lvl]);
      m_amrToPetscFiCo[lvl]->exchange();
    }
  }
}

void
PetscGrid::create(Vec& x) noexcept
{
  CH_TIME("PetscGrid::create");
  if (m_verbose) {
    pout() << "PetscGrid::create" << endl;
  }

  CH_assert(m_isDefined);

  PetscCallVoid(VecCreate(Chombo_MPI::comm, &x));
  PetscCallVoid(VecSetSizes(x, m_numLocalRows, m_numGlobalRows));
  PetscCallVoid(VecSetLocalToGlobalMapping(x, m_localToGlobalIS));
  PetscCallVoid(VecSetFromOptions(x));
}

void
PetscGrid::destroy(Vec& x) noexcept
{
  CH_TIME("PetscGrid::destroy");
  if (m_verbose) {
    pout() << "PetscGrid::destroy" << endl;
  }

  PetscCallVoid(VecDestroy(&x));
}

void
PetscGrid::setValue(Vec& a_x, const PetscScalar a_value) const noexcept
{
  CH_TIME("PetscGrid::setValue");
  if (m_verbose) {
    pout() << "PetscGrid::setValue" << endl;
  }

  PetscInt     startIdx = -1;
  PetscInt     endIdx   = -1;
  PetscScalar* arr      = nullptr;

  PetscCallVoid(VecGetOwnershipRange(a_x, &startIdx, &endIdx));
  PetscCallVoid(VecGetArray(a_x, &arr));

  CH_assert(endIdx <= startIdx);
  CH_assert(endIdx >= 0);
  CH_assert(startIdx >= 0);
  CH_assert(endIdx - startIdx == m_numLocalRows - 1);

  for (PetscInt i = 0; i < m_numLocalRows; i++) {
    arr[i] = a_value;
  }

  PetscCallVoid(VecRestoreArray(a_x, &arr));
}

void
PetscGrid::putChomboInPetsc(Vec& a_x, const MFAMRCellData& a_y) const noexcept
{
  CH_TIME("PetscGrid::putChomboInPetsc");
  if (m_verbose) {
    pout() << "PetscGrid::putChomboInPetsc" << endl;
  }

  CH_assert(m_isDefined);

  // Get the local PETSc vector and do basic error checks.
  PetscInt     startIdx = -1;
  PetscInt     endIdx   = -1;
  PetscScalar* arr      = nullptr;

  PetscCallVoid(VecGetOwnershipRange(a_x, &startIdx, &endIdx));
  PetscCallVoid(VecGetArray(a_x, &arr));

  CH_assert(endIdx <= startIdx);
  CH_assert(endIdx >= 0);
  CH_assert(startIdx >= 0);
  CH_assert(endIdx - startIdx == m_numLocalRows - 1);

  for (int iphase = 0; iphase < m_numPhases; iphase++) {
    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      const EBLevelGrid&       eblg  = m_levelGrids[lvl]->getEBLevelGrid(iphase);
      const EBISLayout&        ebisl = eblg.getEBISL();
      const DisjointBoxLayout& dbl   = eblg.getDBL();
      const DataIterator&      dit   = dbl.dataIterator();
      const int                nbox  = dit.size();

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex din     = dit[mybox];
        const Box&      cellBox = dbl[din];

        const FArrayBox&             data       = (*a_y[lvl])[din].getPhase(iphase).getFArrayBox();
        const BaseFab<PetscAMRCell>& amrToPetsc = (*m_amrToPetsc[lvl])[din];

        auto kernel = [&](const IntVect iv) -> void {
          const PetscInt& k = amrToPetsc(iv).getPetscRow(iphase);

          if (k >= 0) {
            arr[k] = data(iv, 0);
          }
        };

        BoxLoops::loop(cellBox, kernel);
      }
    }
  }

  PetscCallVoid(VecRestoreArray(a_x, &arr));
}

void
PetscGrid::putPetscInChombo(MFAMRCellData& a_y, const Vec& a_x) const noexcept
{
  CH_TIME("PetscGrid::putPetscInChombo");
  if (m_verbose) {
    pout() << "PetscGrid::putPetscInChombo" << endl;
  }

  CH_assert(m_isDefined);

  // Get the local PETSc vector and do basic error checks.
  PetscInt     startIdx = -1;
  PetscInt     endIdx   = -1;
  PetscScalar* arr      = nullptr;

  PetscCallVoid(VecGetOwnershipRange(a_x, &startIdx, &endIdx));
  PetscCallVoid(VecGetArray(a_x, &arr));

  CH_assert(endIdx <= startIdx);
  CH_assert(endIdx >= 0);
  CH_assert(startIdx >= 0);
  CH_assert(endIdx - startIdx == m_numLocalRows - 1);

  for (PetscInt i = 0; i < m_numLocalRows; i++) {
    const PetscDOF dof = m_petscToAMR[i];

    const int       phase     = dof.phase;
    const int       gridLevel = dof.gridLevel;
    const DataIndex gridIndex = dof.gridIndex;
    const IntVect   gridCell  = dof.gridCell;

    FArrayBox& data = (*a_y[gridLevel])[gridIndex].getPhase(phase).getFArrayBox();

    data(gridCell, 0) = arr[i];

#warning "Debug code here"
#if 0
    const PetscAMRCell& amrCell = (*m_amrToPetsc[gridLevel])[gridIndex](gridCell);

    data(gridCell, 0) = 0.0;

    // Ensure that this cell is NOT covered by a finer grid.
    if (amrCell.isCoveredByFinerGrid()) {
      MayDay::Abort("logic bust -- cannot be covered by another cell");
    }

    // Ensure that this cell does not lie across the CF boundary
    if (amrCell.isGhostCF()) {
      MayDay::Abort("logic bust -- cannot be a ghost cell");
    }

    // Fill with values = 1 if the cell lies on the FINE side of the CF boundary.
    Box bx = Box(gridCell, gridCell);
    bx.grow(1);
    bx &= m_levelGrids[gridLevel]->getDomain();

    for (BoxIterator bit(bx); bit.ok(); ++bit) {
      const PetscAMRCell& cell2 = (*m_amrToPetsc[gridLevel])[gridIndex](bit());

      if (cell2.isGhostCF()) {
        data(gridCell, 0) = 1.0;
      }
    }

    // Fill boundary cells with values += 2
    if (amrCell.isDomainBoundaryCell()) {
      data(gridCell, 0) += 2.0;
    }

    data(gridCell, 0) = m_localRowBegin + (*m_amrToPetsc[gridLevel])[gridIndex](gridCell).getPetscRow(phase);
#endif
  }

  PetscCallVoid(VecRestoreArray(a_x, &arr));
}

void
PetscGrid::dumpPetscGrid(const std::string a_filename) const noexcept
{
  CH_TIME("PetscGrid::dumpPetscGrid");
  if (m_verbose) {
    pout() << "PetscGrid::dumpPetscGrid" << endl;
  }

#ifdef CH_USE_HDF5
#warning "Not implemented (yet)"
#endif
}

#include <CD_NamespaceFooter.H>

#endif
