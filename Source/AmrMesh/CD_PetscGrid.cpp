/* chombo-discharge
 * Copyright © 2026 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
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

PetscGrid::PetscGrid() noexcept
{
  CH_TIME("PetscGrid::PetscGrid");

  m_isDefined = false;
  m_verbose   = false;
  m_debug     = false;
  m_profile   = false;

  m_numLocalDOFs  = -1;
  m_numGlobalDOFs = -1;
  m_finestLevel   = -1;
  m_numGhost      = -1;
}

PetscGrid::~PetscGrid() noexcept
{
  CH_TIME("PetscGrid::~PetscGrid");
}

void
PetscGrid::define(const Vector<RefCountedPtr<MFLevelGrid>>&              a_amrGrids,
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

  m_amrGrids    = a_amrGrids;
  m_validCells  = a_validCells;
  m_finestLevel = a_finestLevel;
  m_numGhost    = a_numGhost;

  CH_assert(m_numGhost >= 0);
  CH_assert(m_finestLevel >= 0);

  this->buildPetscMaping();

  m_isDefined = true;
}

void
PetscGrid::clear() noexcept
{
  CH_TIME("PetscGrid::clear");
  if (m_verbose) {
    pout() << "PetscGrid::clear" << endl;
  }
}

void
PetscGrid::buildPetscMaping() noexcept
{
  CH_TIME("PetscGrid::buildPetscMaping");
  if (m_verbose) {
    pout() << "PetscGrid::buildPetscMaping" << endl;
  }

  m_numLocalDOFs  = 0;
  m_numGlobalDOFs = 0;
  m_numDOFsPerRank.resize(numProc(), 0);

  const int myRank    = procID();
  const int numPhases = m_amrGrids[0]->numPhases();

  for (int iphase = 0; iphase < numPhases; iphase++) {
    Vector<RefCountedPtr<LayoutData<BaseFab<PetscInt>>>>& GIDs = m_globalIndices[iphase];
    Vector<RefCountedPtr<LayoutData<BaseFab<PetscInt>>>>& LIDs = m_localIndices[iphase];

    GIDs.resize(1 + m_finestLevel);
    LIDs.resize(1 + m_finestLevel);

    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      const EBLevelGrid&       eblg  = m_amrGrids[lvl]->getEBLevelGrid(iphase);
      const EBISLayout&        ebisl = eblg.getEBISL();
      const DisjointBoxLayout& dbl   = eblg.getDBL();
      const DataIterator&      dit   = dbl.dataIterator();
      const int                nbox  = dit.size();

      GIDs[lvl] = RefCountedPtr<LayoutData<BaseFab<PetscInt>>>(new LayoutData<BaseFab<PetscInt>>(dbl));
      LIDs[lvl] = RefCountedPtr<LayoutData<BaseFab<PetscInt>>>(new LayoutData<BaseFab<PetscInt>>(dbl));

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex      din        = dit[mybox];
        const BaseFab<bool>& validCells = (*m_validCells[lvl])[din];
        const Box&           cellBox    = dbl[din];
        const EBISBox&       ebisBox    = ebisl[din];
        const IntVectSet&    ivs        = ebisBox.getIrregIVS(cellBox);
        const EBGraph&       ebgraph    = ebisBox.getEBGraph();

        BaseFab<PetscInt>& gid = (*GIDs[lvl])[din];
        BaseFab<PetscInt>& lid = (*LIDs[lvl])[din];

        gid.setVal(-1);
        lid.setVal(-1);

        auto regularKernel = [&](const IntVect& iv) -> void {
          if (validCells(iv) && !ebisBox.isCovered(iv)) {
            if (ebisBox.isMultiValued(iv)) {
              MayDay::Abort("PetscGrid::buildPetscMapping -- multi-valued cells are not permitted");
            }

            PetscDOF dof;

            dof.gridLevel = lvl;
            dof.gridIndex = mybox;
            dof.gridCell  = iv;
            dof.phase     = iphase;

            m_petscToAMR.push_back(dof);

            lid(iv)        = m_numLocalDOFs;
            m_numLocalDOFs = m_numLocalDOFs + 1;
          }
        };

        BoxLoops::loop(cellBox, regularKernel);
      }
    }
  }

  m_numDOFsPerRank = ParallelOps::gather(m_numLocalDOFs);
  m_numGlobalDOFs  = ParallelOps::sum(m_numLocalDOFs);

  m_localDOFBegin = 0;
  for (int i = 0; i < procID(); i++) {
    m_localDOFBegin += m_numDOFsPerRank[i];
  }

  if (m_debug) {
    pout() << "PetscGrid::buildPetscMaping - numLocalUnknowns = " << m_numLocalDOFs << endl;
    pout() << "PetscGrid::buildPetscMaping - numGlobalUnknowns = " << m_numGlobalDOFs << endl;
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

  MayDay::Abort("PetscGrid::create - not implemented");
}

void
PetscGrid::destroy(Vec& x) noexcept
{
  CH_TIME("PetscGrid::destroy");
  if (m_verbose) {
    pout() << "PetscGrid::destroy" << endl;
  }

  CH_assert(m_isDefined);

  MayDay::Abort("PetscGrid::destroy - not implemented");
}

void
PetscGrid::putChomboInPetsc(Vec& a_x, const MFAMRCellData& a_y) const noexcept
{
  CH_TIME("PetscGrid::putChomboInPetsc");
  if (m_verbose) {
    pout() << "PetscGrid::putChomboInPetsc" << endl;
  }

  MayDay::Abort("PetscGrid::putChomboInPetsc -- not implemented");
}

void
PetscGrid::putPetscInChombo(MFAMRCellData& a_y, const Vec& a_x) const noexcept
{
  CH_TIME("PetscGrid::putPetscInChombo");
  if (m_verbose) {
    pout() << "PetscGrid::putPetscInChombo" << endl;
  }

  MayDay::Abort("PetscGrid::putPetscInChombo -- not implemented");
}

#include <CD_NamespaceFooter.H>

#endif
