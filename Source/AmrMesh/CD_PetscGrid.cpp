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

#warning "I must probably build the Chombo->Petsc maps including ghost cells (stencils will reach out of patches)"
#warning "Not sure how I want to handle memory management within PetscGrid. Maybe pass this off to the outside world??"
#warning "We should REALLY time how fast transfers between Chombo and PETSc really are"
#warning "PetscGrid should take the coarsened and refined MFLevelGrid as input"
#warning "I really want the MFBaseFab<PetscInt>, no?"
#warning "How much PETSc math should we add to this class?"
#warning "I need to figure out which cells are ghost cells, covered by the geometry, or covered by a finer grid."

// These are all negative because we fill the DOFs on each level by their global
// row #, which is always positive. Cells that require special handling are assigned
// a negative value so we can flag them correctly later on.
const int PetscGrid::InvalidCell  = -4;
const int PetscGrid::CoveredCell  = -3;
const int PetscGrid::GhostCell    = -2;
const int PetscGrid::BoundaryCell = -1;

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

  m_numLocalDOFs  = -1;
  m_numGlobalDOFs = -1;
  m_localDOFBegin = -1;
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
  m_finestLevel    = a_finestLevel;
  m_numGhost       = a_numGhost;
  m_numPhases      = m_levelGrids[0]->numPhases();

  CH_assert(m_numGhost >= 0);
  CH_assert(m_finestLevel >= 0);

  this->definePetscDOFs();
  this->allocateLocalViews();
  this->fillLocalViews();

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
PetscGrid::definePetscDOFs() noexcept
{
  CH_TIME("PetscGrid::definePetscDOFs");
  if (m_verbose) {
    pout() << "PetscGrid::definePetscDOFs" << endl;
  }

  CH_assert(!(m_levelGrids[0].isNull()));

  m_numLocalDOFs  = 0;
  m_numGlobalDOFs = 0;

  m_numDOFsPerRank.resize(numProc(), 0);

  for (int iphase = 0; iphase < m_numPhases; iphase++) {
    Vector<RefCountedPtr<LevelData<BaseFab<PetscInt>>>>& LIDs = m_localIndices[iphase];
    Vector<RefCountedPtr<LevelData<BaseFab<PetscInt>>>>& GIDs = m_globalIndices[iphase];

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
        const DataIndex      din        = dit[mybox];
        const BaseFab<bool>& validCells = (*m_validCells[lvl])[din];
        const Box&           cellBox    = dbl[din];
        const EBISBox&       ebisBox    = ebisl[din];
        const IntVectSet&    ivs        = ebisBox.getIrregIVS(cellBox);
        const EBGraph&       ebgraph    = ebisBox.getEBGraph();

        BaseFab<PetscInt>& lid = (*LIDs[lvl])[din];
        BaseFab<PetscInt>& gid = (*GIDs[lvl])[din];

        lid.setVal(-1);
        gid.setVal(-1);

        auto regularKernel = [&](const IntVect& iv) -> void {
          if (validCells(iv) && !ebisBox.isCovered(iv)) {
            if (ebisBox.isMultiValued(iv)) {
              MayDay::Abort("PetscGrid::definePetscDOFs -- multi-valued cells are not permitted");
            }

            PetscDOF dof;

            dof.gridLevel = lvl;
            dof.gridIndex = din;
            dof.gridCell  = iv;
            dof.phase     = iphase;

            m_petscToAMR.push_back(dof);

            lid(iv) = m_numLocalDOFs;
            gid(iv) = m_numLocalDOFs;

            m_numLocalDOFs = m_numLocalDOFs + 1;
          }
        };

        BoxLoops::loop(cellBox, regularKernel);
      }
    }
  }

  // Gather the number of DOFs per MPI rank. And figure out the total number of DOFs.
  m_numDOFsPerRank = ParallelOps::gather(m_numLocalDOFs);
  m_numGlobalDOFs  = ParallelOps::sum(m_numLocalDOFs);

  // Compute the global row offset for this rank.
  m_localDOFBegin = 0;
  for (PetscInt i = 0; i < procID(); i++) {
    m_localDOFBegin += m_numDOFsPerRank[i];
  }

  // Above, LIDs = GIDs, but we now know the PETSc row where we begin indexing for this rank, so we now
  // increment GIDs by m_localDOFBegin.
  for (int iphase = 0; iphase < m_numPhases; iphase++) {
    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      const DisjointBoxLayout& dbl  = m_levelGrids[lvl]->getGrids();
      const DataIterator&      dit  = dbl.dataIterator();
      const int                nbox = dit.size();

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex din     = dit[mybox];
        const Box       cellBox = dbl[din];

        BaseFab<PetscInt>& gid = (*m_globalIndices[iphase][lvl])[din];

        auto regularKernel = [&](const IntVect& iv) -> void {
          if (gid(iv) >= 0) {
            gid(iv) += m_localDOFBegin;
          }
        };

        BoxLoops::loop(cellBox, regularKernel);
      }
    }
  }

  // Build local to global mapping for PETSc, so that we can extract local subvectors.
  PetscInt* gidx;
  PetscCallVoid(PetscMalloc1(m_numLocalDOFs, &gidx));
  for (PetscInt i = 0; i < m_numLocalDOFs; i++) {
    gidx[i] = m_localDOFBegin + i;
  }

  PetscCallVoid(
    ISLocalToGlobalMappingCreate(Chombo_MPI::comm, 1, m_numLocalDOFs, gidx, PETSC_OWN_POINTER, &m_localToGlobalIS));

  if (m_debug) {
    pout() << "PetscGrid::definePetscDOFs - numLocalUnknowns = " << m_numLocalDOFs << endl;
    pout() << "PetscGrid::definePetscDOFs - numGlobalUnknowns = " << m_numGlobalDOFs << endl;
  }
}

void
PetscGrid::allocateLocalViews() noexcept
{
  CH_TIME("PetscGrid::allocateLocalViews");
  if (m_verbose) {
    pout() << "PetscGrid::allocateLocalViews" << endl;
  }

  const int     numComps = 1;
  const IntVect ghostVec = m_numGhost * IntVect::Unit;

  m_levelDOFs.resize(1 + m_finestLevel);
  m_levelDOFsCoFi.resize(1 + m_finestLevel);
  m_levelDOFsFiCo.resize(1 + m_finestLevel);

  // Define all DOF data -- set everything to invalid cell for now.
  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const bool hasCoar = (lvl > 0);
    const bool hasFine = (lvl < m_finestLevel);

    {
      const DisjointBoxLayout& dbl  = m_levelGrids[lvl]->getGrids();
      const DataIterator&      dit  = dbl.dataIterator();
      const int                nbox = dit.size();

      m_levelDOFs[lvl] = RefCountedPtr<LevelData<BaseFab<PetscInt>>>(
        new LevelData<BaseFab<PetscInt>>(dbl, numComps, ghostVec));

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        (*m_levelDOFs[lvl])[dit[mybox]].setVal(PetscGrid::InvalidCell);
      }
    }

    if (hasFine) {
      const DisjointBoxLayout& dbl  = m_levelGridsCoFi[lvl]->getGrids();
      const DataIterator&      dit  = dbl.dataIterator();
      const int                nbox = dit.size();

      m_levelDOFsCoFi[lvl] = RefCountedPtr<LevelData<BaseFab<PetscInt>>>(
        new LevelData<BaseFab<PetscInt>>(dbl, numComps, ghostVec));

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        (*m_levelDOFsCoFi[lvl])[dit[mybox]].setVal(PetscGrid::InvalidCell);
      }
    }

    if (hasCoar) {
      const DisjointBoxLayout& dbl  = m_levelGridsFiCo[lvl]->getGrids();
      const DataIterator&      dit  = dbl.dataIterator();
      const int                nbox = dit.size();

      m_levelDOFsFiCo[lvl] = RefCountedPtr<LevelData<BaseFab<PetscInt>>>(
        new LevelData<BaseFab<PetscInt>>(dbl, numComps, ghostVec));

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        (*m_levelDOFsFiCo[lvl])[dit[mybox]].setVal(PetscGrid::InvalidCell);
      }
    }
  }
}

void
PetscGrid::fillLocalViews() noexcept
{
  CH_TIME("PetscGrid::fillLocalViews");
  if (m_verbose) {
    pout() << "PetscGrid::fillLocalViews" << endl;
  }

  const int     numComps = 1;
  const IntVect ghostVec = m_numGhost * IntVect::Unit;
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
  PetscCallVoid(VecSetSizes(x, m_numLocalDOFs, m_numGlobalDOFs));
  PetscCallVoid(VecSetLocalToGlobalMapping(x, m_localToGlobalIS));
  PetscCallVoid(VecSetFromOptions(x));

#warning "Debug code in PetscGrid::create -- scheduled for removal"
#if 0
  PetscInt    indices[m_numLocalDOFs];
  PetscScalar values[m_numLocalDOFs];

  for (PetscInt i = 0; i < m_numLocalDOFs; i++) {
    indices[i] = i;
    values[i]  = 1.0 * procID();
  }

  PetscCallVoid(VecSetValuesLocal(x, m_numLocalDOFs, indices, values, INSERT_VALUES));
  PetscCallVoid(VecAssemblyBegin(x));
  PetscCallVoid(VecAssemblyEnd(x));
  PetscCallVoid(VecView(x, PETSC_VIEWER_STDOUT_WORLD));
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
  CH_assert(endIdx - startIdx == m_numLocalDOFs - 1);

  for (PetscInt i = 0; i < m_numLocalDOFs; i++) {
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
  CH_assert(endIdx - startIdx == m_numLocalDOFs - 1);

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
        const BaseFab<PetscInt>& localIndex = (*m_localIndices[iphase][lvl])[din];

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
  CH_assert(endIdx - startIdx == m_numLocalDOFs - 1);

  for (PetscInt i = 0; i < m_numLocalDOFs; i++) {
    const PetscDOF dof = m_petscToAMR[i];

    const int       phase     = dof.phase;
    const int       gridLevel = dof.gridLevel;
    const DataIndex gridIndex = dof.gridIndex;
    const IntVect   gridCell  = dof.gridCell;

    FArrayBox& data = (*a_y[gridLevel])[gridIndex].getPhase(phase).getFArrayBox();

    data(gridCell, 0) = arr[i];
  }

  PetscCallVoid(VecRestoreArray(a_x, &arr));
}

#include <CD_NamespaceFooter.H>

#endif
