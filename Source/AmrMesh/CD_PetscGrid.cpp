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

  this->buildLocalToGlobalMapping();

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
PetscGrid::buildLocalToGlobalMapping() noexcept
{
  CH_TIME("PetscGrid::buildLocalToGlobalMapping");
  if (m_verbose) {
    pout() << "PetscGrid::buildLocalToGlobalMapping" << endl;
  }

  m_numLocalDOFs  = 0;
  m_numGlobalDOFs = 0;
  m_numDOFsPerRank.resize(numProc(), 0);

  const int myRank = procID();

  const int numPhases = m_amrGrids[0]->numPhases();

  // 1. Figure out the number of unknowns per rank.
  for (int iphase = 0; iphase < numPhases; iphase++) {
    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      const EBLevelGrid&       eblg  = m_amrGrids[lvl]->getEBLevelGrid(iphase);
      const EBISLayout&        ebisl = eblg.getEBISL();
      const DisjointBoxLayout& dbl   = eblg.getDBL();
      const DataIterator&      dit   = dbl.dataIterator();
      const int                nbox  = dit.size();

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex      din        = dit[mybox];
        const BaseFab<bool>& validCells = (*m_validCells[lvl])[din];
        const Box&           cellBox    = dbl[din];
        const EBISBox&       ebisBox    = ebisl[din];
        const IntVectSet&    ivs        = ebisBox.getIrregIVS(cellBox);
        const EBGraph&       ebgraph    = ebisBox.getEBGraph();

        auto regularKernel = [&](const IntVect& iv) -> void {
          if (validCells(iv) && ebisBox.isRegular(iv)) {
            m_numLocalDOFs += 1;
          }
        };

        auto irregularKernel = [&](const VolIndex& vof) -> void {
          const IntVect iv = vof.gridIndex();
          if (validCells(iv) && ebisBox.isIrregular(iv)) {
            m_numLocalDOFs += 1;
          }
        };

        VoFIterator vofit(ivs, ebgraph);

        BoxLoops::loop(cellBox, regularKernel);
        BoxLoops::loop(vofit, irregularKernel);
      }
    }
  }

  m_numDOFsPerRank = ParallelOps::gather(m_numLocalDOFs);
  m_numGlobalDOFs  = ParallelOps::sum(m_numLocalDOFs);

  if (m_debug) {
    pout() << "PetscGrid::buildLocalToGlobalMapping - numLocalUnknowns = " << m_numLocalDOFs << endl;
    pout() << "PetscGrid::buildLocalToGlobalMapping - numGlobalUnknowns = " << m_numGlobalDOFs << endl;
  }
}

void
PetscGrid::putChomboInPetsc(Vec& a_x, const MFAMRCellData& a_y) const noexcept
{
  CH_TIME("PetscGrid::putChomboInPetsc");
  if (m_verbose) {
    pout() << "PetscGrid::putChomboInPetsc" << endl;
  }
}

void
PetscGrid::putPetscInChombo(MFAMRCellData& a_y, const Vec& a_x) const noexcept
{
  CH_TIME("PetscGrid::putPetscInChombo");
  if (m_verbose) {
    pout() << "PetscGrid::putPetscInChombo" << endl;
  }
}

#include <CD_NamespaceFooter.H>

#endif
