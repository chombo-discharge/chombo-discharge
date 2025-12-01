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
#warning "At the end -- let's see if we can trim some memory from this class"
#warning "I must probably build the Chombo->Petsc maps including ghost cells (stencils will reach out of patches)"
#warning "Not sure how I want to handle memory management within PetscGrid. Maybe pass this off to the outside world??"
#warning "We should REALLY time how fast transfers between Chombo and PETSc really are"
#warning "I need to figure out which cells are ghost cells, covered by the geometry, or covered by a finer grid."
#warning "Priority 1: Build maps over which cells are which. Must be viewable from each patch".
#warning "Priority 2: Map global row numbers to cells on each phase. Must be viewable from each patch".
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
PetscGrid::define(const Vector<RefCountedPtr<MFLevelGrid>>&                 a_levelGrids,
                  const Vector<RefCountedPtr<MFLevelGrid>>&                 a_levelGridsCoFi,
                  const Vector<RefCountedPtr<MFLevelGrid>>&                 a_levelGridsFiCo,
                  const Vector<RefCountedPtr<LevelData<BaseFab<AMRCell>>>>& a_amrCells,
                  const int                                                 a_finestLevel,
                  const int                                                 a_numGhost) noexcept
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
  m_amrCells       = a_amrCells;
  m_finestLevel    = a_finestLevel;
  m_numGhost       = a_numGhost;
  m_numPhases      = m_levelGrids[0]->numPhases();

  CH_assert(m_numGhost >= 0);
  CH_assert(m_finestLevel >= 0);

  this->definePetscRows();
  this->defineCellFlags();
  this->defineRowViews();

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
    Vector<RefCountedPtr<LevelData<BaseFab<PetscInt>>>>& LIDs = m_localRows[iphase];
    Vector<RefCountedPtr<LevelData<BaseFab<PetscInt>>>>& GIDs = m_globalRows[iphase];

    LIDs.resize(1 + m_finestLevel);
    GIDs.resize(1 + m_finestLevel);

    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      const EBLevelGrid&       eblg  = m_levelGrids[lvl]->getEBLevelGrid(iphase);
      const EBISLayout&        ebisl = eblg.getEBISL();
      const DisjointBoxLayout& dbl   = eblg.getDBL();
      const DataIterator&      dit   = dbl.dataIterator();
      const int                nbox  = dit.size();

      LIDs[lvl] = RefCountedPtr<LevelData<BaseFab<PetscInt>>>(new LevelData<BaseFab<PetscInt>>(dbl, 1, IntVect::Zero));
      GIDs[lvl] = RefCountedPtr<LevelData<BaseFab<PetscInt>>>(new LevelData<BaseFab<PetscInt>>(dbl, 1, IntVect::Zero));

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex         din      = dit[mybox];
        const BaseFab<AMRCell>& amrCells = (*m_amrCells[lvl])[din];
        const Box&              cellBox  = dbl[din];
        const EBISBox&          ebisBox  = ebisl[din];
        const IntVectSet&       ivs      = ebisBox.getIrregIVS(cellBox);
        const EBGraph&          ebgraph  = ebisBox.getEBGraph();

        BaseFab<PetscInt>& lid = (*LIDs[lvl])[din];
        BaseFab<PetscInt>& gid = (*GIDs[lvl])[din];

        lid.setVal(-1);
        gid.setVal(-1);

        auto regularKernel = [&](const IntVect& iv) -> void {
          if (!(amrCells(iv).isCoveredByFinerGrid()) && !ebisBox.isCovered(iv)) {
            if (ebisBox.isMultiValued(iv)) {
              MayDay::Abort("PetscGrid::definePetscRows -- multi-valued cells are not permitted");
            }

            PetscDOF dof;

            dof.gridLevel = lvl;
            dof.gridIndex = din;
            dof.gridCell  = iv;
            dof.phase     = iphase;

            m_petscToAMR.push_back(dof);

            lid(iv) = m_numLocalRows;
            gid(iv) = m_numLocalRows;

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

  // Above, LIDs = GIDs, but we now know the PETSc row where we begin indexing for this rank, so we now
  // increment GIDs by m_localRowBegin.
  for (int iphase = 0; iphase < m_numPhases; iphase++) {
    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      const DisjointBoxLayout& dbl  = m_levelGrids[lvl]->getGrids();
      const DataIterator&      dit  = dbl.dataIterator();
      const int                nbox = dit.size();

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex din     = dit[mybox];
        const Box       cellBox = dbl[din];

        BaseFab<PetscInt>& gid = (*m_globalRows[iphase][lvl])[din];

        auto regularKernel = [&](const IntVect& iv) -> void {
          if (gid(iv) >= 0) {
            gid(iv) += m_localRowBegin;
          }
        };

        BoxLoops::loop(cellBox, regularKernel);
      }
    }
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
PetscGrid::defineRowViews() noexcept
{
  CH_TIME("PetscGrid::defineRowViews");
  if (m_verbose) {
    pout() << "PetscGrid::defineRowViews" << endl;
  }
#if 0
  const int     numComps = 1;
  const IntVect ghostVec = m_numGhost * IntVect::Unit;

  // Allocate, but do not assign anything yet.
  for (int iphase = 0; iphase < m_numPhases; iphase++) {
    Vector<RefCountedPtr<LevelData<BaseFab<PetscInt>>>>& rows     = m_rows[iphase];
    Vector<RefCountedPtr<LevelData<BaseFab<PetscInt>>>>& rowsCoFi = m_rowsCoFi[iphase];
    Vector<RefCountedPtr<LevelData<BaseFab<PetscInt>>>>& rowsFiCo = m_rowsFiCo[iphase];

    rows.resize(1 + m_finestLevel);
    rowsCoFi.resize(1 + m_finestLevel);
    rowsFiCo.resize(1 + m_finestLevel);

    // Define all row data -- set everything to an invalid row for now.
    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      const bool hasCoar = (lvl > 0);
      const bool hasFine = (lvl < m_finestLevel);

      {
        const DisjointBoxLayout& dbl  = m_levelGrids[lvl]->getGrids();
        const DataIterator&      dit  = dbl.dataIterator();
        const int                nbox = dit.size();

        rows[lvl] = RefCountedPtr<LevelData<BaseFab<PetscInt>>>(
          new LevelData<BaseFab<PetscInt>>(dbl, numComps, ghostVec));

#pragma omp parallel for schedule(runtime)
        for (int mybox = 0; mybox < nbox; mybox++) {
          (*rows[lvl])[dit[mybox]].setVal(-1);
        }
      }

      if (hasFine) {
        const DisjointBoxLayout& dbl  = m_levelGridsCoFi[lvl]->getGrids();
        const DataIterator&      dit  = dbl.dataIterator();
        const int                nbox = dit.size();

        rowsCoFi[lvl] = RefCountedPtr<LevelData<BaseFab<PetscInt>>>(
          new LevelData<BaseFab<PetscInt>>(dbl, numComps, ghostVec));

#pragma omp parallel for schedule(runtime)
        for (int mybox = 0; mybox < nbox; mybox++) {
          (*rowsCoFi[lvl])[dit[mybox]].setVal(-1);
        }
      }

      if (hasCoar) {
        const DisjointBoxLayout& dbl  = m_levelGridsFiCo[lvl]->getGrids();
        const DataIterator&      dit  = dbl.dataIterator();
        const int                nbox = dit.size();

        rowsFiCo[lvl] = RefCountedPtr<LevelData<BaseFab<PetscInt>>>(
          new LevelData<BaseFab<PetscInt>>(dbl, numComps, ghostVec));

#pragma omp parallel for schedule(runtime)
        for (int mybox = 0; mybox < nbox; mybox++) {
          (*rowsFiCo[lvl])[dit[mybox]].setVal(-1);
        }
      }
    }
  }
#endif
}

void
PetscGrid::defineCellFlags() noexcept
{
  CH_TIME("PetscGrid::defineCellFlags");
  if (m_verbose) {
    pout() << "PetscGrid::defineCellFlags" << endl;
  }
#if 0
  const int     numComps = 1;
  const IntVect ghostVec = m_numGhost * IntVect::Unit;

  m_cellFlags.resize(1 + m_finestLevel);
  m_cellFlagsCoFi.resize(1 + m_finestLevel);
  m_cellFlagsFiCo.resize(1 + m_finestLevel);

  // Allocate storages -- set everything to an invalid cell.
  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const bool hasCoar = (lvl > 0);
    const bool hasFine = (lvl < m_finestLevel);

    {
      const DisjointBoxLayout& dbl  = m_levelGrids[lvl]->getGrids();
      const DataIterator&      dit  = dbl.dataIterator();
      const int                nbox = dit.size();

      m_cellFlags[lvl] = RefCountedPtr<LevelData<BaseFab<int>>>(new LevelData<BaseFab<int>>(dbl, numComps, ghostVec));

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        //        (*m_cellFlags[lvl])[dit[mybox]].setVal(PetscGrid::InvalidCell);
      }
    }

    if (hasFine) {
      const DisjointBoxLayout& dbl  = m_levelGridsCoFi[lvl]->getGrids();
      const DataIterator&      dit  = dbl.dataIterator();
      const int                nbox = dit.size();

      m_cellFlagsCoFi[lvl] = RefCountedPtr<LevelData<BaseFab<int>>>(
        new LevelData<BaseFab<int>>(dbl, numComps, ghostVec));

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        //        (*m_cellFlagsCoFi[lvl])[dit[mybox]].setVal(PetscGrid::InvalidCell);
      }
    }

    if (hasCoar) {
      const DisjointBoxLayout& dbl  = m_levelGridsFiCo[lvl]->getGrids();
      const DataIterator&      dit  = dbl.dataIterator();
      const int                nbox = dit.size();

      m_cellFlagsFiCo[lvl] = RefCountedPtr<LevelData<BaseFab<int>>>(
        new LevelData<BaseFab<int>>(dbl, numComps, ghostVec));

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        //        (*m_cellFlagsFiCo[lvl])[dit[mybox]].setVal(PetscGrid::InvalidCell);
      }
    }
  }

  // Fill cell flags
  for (int lvl = m_finestLevel; lvl >= 0; lvl--) {
    const DisjointBoxLayout& dbl      = m_levelGrids[lvl]->getGrids();
    const DataIterator&      dit      = dbl.dataIterator();
    const int                numBoxes = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < numBoxes; mybox++) {
      const DataIndex& din     = dit[mybox];
      const Box        cellBox = dbl[din];

      BaseFab<int>&        cellFlags  = (*m_cellFlags[lvl])[din];
      const BaseFab<bool>& amrCells = (*m_amrCells[lvl])[din];

      auto regularKernel = [&](const IntVect& iv) -> void {

      };
    }
  }
#endif
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

#warning "Debug code in PetscGrid::create -- scheduled for removal"
#if 0
  PetscInt    indices[m_numLocalRows];
  PetscScalar values[m_numLocalRows];

  for (PetscInt i = 0; i < m_numLocalRows; i++) {
    indices[i] = i;
    values[i]  = 1.0 * procID();
  }

  PetscCallVoid(VecSetValuesLocal(x, m_numLocalRows, indices, values, INSERT_VALUES));
  PetscCallVoid(VecAssemblyBegin(x));
  PetscCallVoid(VecAssemblyEnd(x));
  PetscCallVoid(VecView(x, PETSC_VIEWER_STDOUT_WORLD));
#endif
#endif
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

        const FArrayBox&         data       = (*a_y[lvl])[din].getPhase(iphase).getFArrayBox();
        const BaseFab<PetscInt>& localIndex = (*m_localRows[iphase][lvl])[din];

        auto kernel = [&](const IntVect iv) -> void {
          const PetscInt& k = localIndex(iv);
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

#warning "Debug code enabled"
#if 1
    const AMRCell& amrCell = (*m_amrCells[gridLevel])[gridIndex](gridCell);

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
      const AMRCell& cell2 = (*m_amrCells[gridLevel])[gridIndex](bit());

      if (cell2.isGhostCF()) {
        data(gridCell, 0) = 1.0;
      }
    }

    // Fill boundary cells with values += 2
    if(amrCell.isDomainBoundaryCell()){
      data(gridCell,0) += 2.0;
    }

    data(gridCell,0) = (*m_amrCells[gridLevel])[gridIndex](gridCell).getNumPhases();
#endif
  }

  PetscCallVoid(VecRestoreArray(a_x, &arr));
}

#include <CD_NamespaceFooter.H>
