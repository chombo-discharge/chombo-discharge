/**
 * chombo-discharge
 * Copyright © 2026 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/**
 * @file   CD_PetscGrid.cpp
 * @brief  Implementation of CD_PetscGrid.H
 * @author Robert Marskar
 */

#ifdef CH_USE_PETSC

// Std includes
#include <sstream>

// Chombo includes
#include <CH_Timer.H>
#include <AMRIO.H>
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
  if (m_verbose) {
    pout() << "PetscGrid::~PetscGrid" << endl;
  }
}

void
PetscGrid::define(const Vector<RefCountedPtr<MFLevelGrid>>&              a_levelGrids,
                  const Vector<RefCountedPtr<MFLevelGrid>>&              a_levelGridsCoFi,
                  const Vector<RefCountedPtr<MFLevelGrid>>&              a_levelGridsFiCo,
                  const Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>& a_validCells,
                  const Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>& a_coarHaloCF,
                  const Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>& a_fineHaloCF,
                  const Vector<int>&                                     a_refinementRatios,
                  const Vector<Real>&                                    a_dx,
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

  CH_assert(a_finestLevel >= 0);
  CH_assert(a_numGhost >= 0);
  CH_assert(a_levelGrids.size() > a_finestLevel);
  CH_assert(a_levelGridsCoFi.size() > a_finestLevel);
  CH_assert(a_levelGridsFiCo.size() > a_finestLevel);
  CH_assert(a_validCells.size() > a_finestLevel);
  CH_assert(a_coarHaloCF.size() > a_finestLevel);
  CH_assert(a_fineHaloCF.size() > a_finestLevel);
  CH_assert(a_refinementRatios.size() >= a_finestLevel);
  CH_assert(a_dx.size() > a_finestLevel);
  CH_assert(!(a_levelGrids[0].isNull()));

  m_levelGrids     = a_levelGrids;
  m_levelGridsCoFi = a_levelGridsCoFi;
  m_levelGridsFiCo = a_levelGridsFiCo;
  m_validCells     = a_validCells;
  m_coarHaloCF     = a_coarHaloCF;
  m_fineHaloCF     = a_fineHaloCF;
  m_refRat         = a_refinementRatios;
  m_dx             = a_dx;

  m_finestLevel = a_finestLevel;
  m_numGhost    = a_numGhost;
  m_numPhases   = m_levelGrids[0]->numPhases();

  // m_petscToAMR and m_vofIterators hold at most two phases. Checked in every build: a third phase would write past
  // the end of both arrays.
  if (m_numPhases != 1 && m_numPhases != 2) {
    MayDay::Abort("PetscGrid::define -- only one or two phases are supported");
  }

  this->defineAMRCells();
  this->definePetscRows();
  this->defineCoFiBuffers();

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

  // The index set is gone, so a second clear() must not destroy it again.
  m_isDefined = false;
}

bool
PetscGrid::isDefined() const noexcept
{
  CH_TIME("PetscGrid::isDefined");
  if (m_verbose) {
    pout() << "PetscGrid::isDefined" << endl;
  }

  return m_isDefined;
}

int
PetscGrid::getNumPhases() const noexcept
{
  CH_TIME("PetscGrid::getNumPhases");
  if (m_verbose) {
    pout() << "PetscGrid::getNumPhases" << endl;
  }

  CH_assert(m_isDefined);

  return m_numPhases;
}

const Vector<RefCountedPtr<MFLevelGrid>>&
PetscGrid::getMFLevelGrids() const noexcept
{
  CH_TIME("PetscGrid::getMFLevelGrids");
  if (m_verbose) {
    pout() << "PetscGrid::getMFLevelGrids" << endl;
  }

  CH_assert(m_isDefined);

  return m_levelGrids;
}

const Vector<RefCountedPtr<MFLevelGrid>>&
PetscGrid::getMFLevelGridsCoFi() const noexcept
{
  CH_TIME("PetscGrid::getMFLevelGridsCoFi");
  if (m_verbose) {
    pout() << "PetscGrid::getMFLevelGridsCoFi" << endl;
  }

  CH_assert(m_isDefined);

  return m_levelGridsCoFi;
}

const Vector<RefCountedPtr<MFLevelGrid>>&
PetscGrid::getMFLevelGridsFiCo() const noexcept
{
  CH_TIME("PetscGrid::getMFLevelGridsFiCo");
  if (m_verbose) {
    pout() << "PetscGrid::getMFLevelGridsFiCo" << endl;
  }

  CH_assert(m_isDefined);

  return m_levelGridsFiCo;
}

const Vector<RefCountedPtr<LevelData<BaseFab<PetscAMRCell>>>>&
PetscGrid::getAMRToPetsc() const noexcept
{
  CH_TIME("PetscGrid::getAMRToPetsc");
  if (m_verbose) {
    pout() << "PetscGrid::getAMRToPetsc" << endl;
  }

  CH_assert(m_isDefined);

  return m_amrToPetsc;
}

const Vector<RefCountedPtr<LevelData<BaseFab<PetscAMRCell>>>>&
PetscGrid::getAMRToPetscCoFi() const noexcept
{
  CH_TIME("PetscGrid::getAMRToPetscCoFi");
  if (m_verbose) {
    pout() << "PetscGrid::getAMRToPetscCoFi" << endl;
  }

  CH_assert(m_isDefined);

  return m_amrToPetscCoFi;
}

const Vector<RefCountedPtr<LevelData<BaseFab<PetscAMRCell>>>>&
PetscGrid::getAMRToPetscFiCo() const noexcept
{
  CH_TIME("PetscGrid::getAMRToPetscFiCo");
  if (m_verbose) {
    pout() << "PetscGrid::getAMRToPetscFiCo" << endl;
  }

  CH_assert(m_isDefined);

  return m_amrToPetscFiCo;
}

Vector<RefCountedPtr<LayoutData<VoFIterator>>>&
PetscGrid::getVoFIterator(const int a_iphase) const noexcept
{
  CH_TIME("PetscGrid::getVoFIterator");
  if (m_verbose) {
    pout() << "PetscGrid::getVoFIterator" << endl;
  }

  CH_assert(m_isDefined);
  CH_assert(a_iphase >= 0 && a_iphase < m_numPhases);

  return m_vofIterators[a_iphase];
}

PetscStencil
PetscGrid::getPetscStencil(const VoFStencil& a_stencil,
                           const int         a_level,
                           const DataIndex&  a_din,
                           const int         a_phase) const noexcept
{
  CH_TIME("PetscGrid::getPetscStencil");
  if (m_verbose) {
    pout() << "PetscGrid::getPetscStencil" << endl;
  }

  CH_assert(m_isDefined);
  CH_assert(a_level >= 0 && a_level <= m_finestLevel);
  CH_assert(a_phase >= 0 && a_phase < m_numPhases);

  PetscStencil petscStencil;

  this->addRows(petscStencil, a_stencil, (*m_amrToPetsc[a_level])[a_din], a_phase, "stencil");

  petscStencil.compress();

  return petscStencil;
}

PetscStencil
PetscGrid::getPetscStencilWithFinerLevel(const VoFStencil& a_stencil,
                                         const VoFStencil& a_fineStencil,
                                         const int         a_level,
                                         const DataIndex&  a_din,
                                         const int         a_phase) const noexcept
{
  CH_TIME("PetscGrid::getPetscStencilWithFinerLevel");
  if (m_verbose) {
    pout() << "PetscGrid::getPetscStencilWithFinerLevel" << endl;
  }

  CH_assert(m_isDefined);
  CH_assert(a_level >= 0 && a_level < m_finestLevel);
  CH_assert(a_phase >= 0 && a_phase < m_numPhases);

  PetscStencil petscStencil;

  // The level a_level+1 rows on the refinement of grid a_din share its DataIndex.
  this->addRows(petscStencil, a_stencil, (*m_amrToPetsc[a_level])[a_din], a_phase, "stencil");
  this->addRows(petscStencil, a_fineStencil, (*m_amrToPetscFiCo[a_level + 1])[a_din], a_phase, "fine stencil");

  petscStencil.compress();

  return petscStencil;
}

PetscStencil
PetscGrid::getPetscStencilWithCoarserLevel(const VoFStencil& a_stencil,
                                           const VoFStencil& a_coarStencil,
                                           const int         a_level,
                                           const DataIndex&  a_din,
                                           const int         a_phase) const noexcept
{
  CH_TIME("PetscGrid::getPetscStencilWithCoarserLevel");
  if (m_verbose) {
    pout() << "PetscGrid::getPetscStencilWithCoarserLevel" << endl;
  }

  CH_assert(m_isDefined);
  CH_assert(a_level > 0 && a_level <= m_finestLevel);
  CH_assert(a_phase >= 0 && a_phase < m_numPhases);

  PetscStencil petscStencil;

  // The level a_level-1 rows on the coarsening of grid a_din share its DataIndex.
  this->addRows(petscStencil, a_stencil, (*m_amrToPetsc[a_level])[a_din], a_phase, "stencil");
  this->addRows(petscStencil, a_coarStencil, (*m_amrToPetscCoFi[a_level - 1])[a_din], a_phase, "coarse stencil");

  petscStencil.compress();

  return petscStencil;
}

void
PetscGrid::addRows(PetscStencil&                a_petscStencil,
                   const VoFStencil&            a_stencil,
                   const BaseFab<PetscAMRCell>& a_map,
                   const int                    a_phase,
                   const std::string&           a_what) const noexcept
{
  CH_TIME("PetscGrid::addRows");
  if (m_verbose) {
    pout() << "PetscGrid::addRows" << endl;
  }

  CH_assert(a_phase >= 0 && a_phase < m_numPhases);

  // Both checks abort in every build rather than being assertions. An entry outside the map would read out of bounds,
  // and PETSc silently ignores negative indices in MatSetValues, so a missing row would drop the entry from the matrix
  // without an error. Each check is a single compare on the normal path; the diagnosis is in abortInvalidEntry.
  for (int i = 0; i < a_stencil.size(); i++) {
    const IntVect& iv = a_stencil.vof(i).gridIndex();

    if (!(a_map.box().contains(iv))) {
      this->abortInvalidEntry(iv, a_map, a_phase, a_what);
    }

    const PetscInt row = a_map(iv).getPetscRow(a_phase);

    if (row < 0) {
      this->abortInvalidEntry(iv, a_map, a_phase, a_what);
    }

    a_petscStencil.add(row, a_stencil.weight(i));
  }
}

void
PetscGrid::abortInvalidEntry(const IntVect&               a_iv,
                             const BaseFab<PetscAMRCell>& a_map,
                             const int                    a_phase,
                             const std::string&           a_what) const noexcept
{
  CH_TIME("PetscGrid::abortInvalidEntry");
  if (m_verbose) {
    pout() << "PetscGrid::abortInvalidEntry" << endl;
  }

  std::stringstream msg;
  msg << "PetscGrid::addRows -- " << a_what << " entry " << a_iv << " on phase " << a_phase;

  if (!(a_map.box().contains(a_iv))) {
    msg << " lies outside the PETSc map, whose box with ghost cells is " << a_map.box();
  }
  else if (a_map(a_iv).isCoveredByFinerGrid()) {
    msg << " is not a PETSc row: it is covered by a finer level";
  }
  else if (a_map(a_iv).isGhostCF()) {
    msg << " is not a PETSc row: it is a ghost cell outside the level (across a refinement boundary or the domain "
           "boundary)";
  }
  else {
    msg << " is not a PETSc row: it is covered by the EB on this phase";
  }

  MayDay::Abort(msg.str().c_str());
}

void
PetscGrid::defineAMRCells() noexcept
{
  CH_TIME("PetscGrid::defineAMRCells");
  if (m_verbose) {
    pout() << "PetscGrid::defineAMRCells" << endl;
  }

  CH_assert(m_finestLevel >= 0);
  CH_assert(m_levelGrids.size() > m_finestLevel);
  CH_assert(m_validCells.size() > m_finestLevel);
  CH_assert(m_coarHaloCF.size() > m_finestLevel);
  CH_assert(m_fineHaloCF.size() > m_finestLevel);

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

      BoxLoops::loop<D_DECL(1, 1, 1)>(ghostedCellBox, setValidCells);

      // Fix CF region flags
      if (lvl < m_finestLevel) {
        const BaseFab<bool>& haloCells = (*m_coarHaloCF[lvl])[din];

        auto setOuterCFRegion = [&](const IntVect& iv) -> void {
          amrToPetsc(iv).setCoarCF(haloCells(iv));
        };

        BoxLoops::loop<D_DECL(1, 1, 1)>(cellBox, setOuterCFRegion);
      }

      if (lvl > 0) {
        const BaseFab<bool>& haloCells = (*m_fineHaloCF[lvl])[din];

        auto setInnerCFRegion = [&](const IntVect& iv) -> void {
          amrToPetsc(iv).setFineCF(haloCells(iv));
        };

        BoxLoops::loop<D_DECL(1, 1, 1)>(cellBox, setInnerCFRegion);
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

        BoxLoops::loop<D_DECL(1, 1, 1)>(cellBox, setDomainFlag);
      }
    }
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
  CH_assert(m_amrToPetsc.size() > m_finestLevel);

  m_numLocalRows  = 0;
  m_numGlobalRows = 0;

  m_numRowsPerRank.resize(numProc(), 0);

  for (int iphase = 0; iphase < m_numPhases; iphase++) {
    m_petscToAMR[iphase].resize(1 + m_finestLevel);
    m_vofIterators[iphase].resize(1 + m_finestLevel);

    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      const EBLevelGrid&       eblg  = m_levelGrids[lvl]->getEBLevelGrid(iphase);
      const EBISLayout&        ebisl = eblg.getEBISL();
      const DisjointBoxLayout& dbl   = eblg.getDBL();
      const DataIterator&      dit   = dbl.dataIterator();
      const int                nbox  = dit.size();

      m_petscToAMR[iphase][lvl]   = RefCountedPtr<LayoutData<Vector<PetscDOF>>>(new LayoutData<Vector<PetscDOF>>(dbl));
      m_vofIterators[iphase][lvl] = RefCountedPtr<LayoutData<VoFIterator>>(new LayoutData<VoFIterator>(dbl));

      // PS: This loop is not thread safe!
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex din     = dit[mybox];
        const Box&      cellBox = dbl[din];
        const EBISBox&  ebisBox = ebisl[din];
        const EBGraph&  ebgraph = ebisBox.getEBGraph();

        // Cut cells that are rows. The kernel aborts on multi-valued cells, so each of these is a single VoF.
        IntVectSet irregIVS;

        BaseFab<PetscAMRCell>& amrToPetsc = (*m_amrToPetsc[lvl])[din];
        Vector<PetscDOF>&      petscToAMR = (*m_petscToAMR[iphase][lvl])[din];

        petscToAMR.resize(0);

        auto regularKernel = [&](const IntVect& iv) -> void {
          if (!(amrToPetsc(iv).isCoveredByFinerGrid()) && !ebisBox.isCovered(iv)) {
            if (ebisBox.isMultiValued(iv)) {
              MayDay::Abort("PetscGrid::definePetscRows -- multi-valued cells are not permitted!");
            }
            if (ebisBox.isIrregular(iv)) {
              irregIVS |= iv;
            }

            PetscDOF dof;

            dof.gridCell = iv;
            dof.localRow = m_numLocalRows;

            petscToAMR.push_back(dof);
            amrToPetsc(iv).setPetscRow(iphase, m_numLocalRows);

            m_numLocalRows = m_numLocalRows + 1;
          }
        };

        BoxLoops::loop<D_DECL(1, 1, 1)>(cellBox, regularKernel);

        (*m_vofIterators[iphase][lvl])[din].define(irregIVS, ebgraph);
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

  // The rows installed above are local to this rank. Offset them to global rows, then fill ghost cells so that they
  // hold the global rows of neighbouring boxes, including boxes owned by other ranks.
  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl  = m_levelGrids[lvl]->getGrids();
    const DataIterator&      dit  = dbl.dataIterator();
    const int                nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex din = dit[mybox];

      BaseFab<PetscAMRCell>& amrToPetsc = (*m_amrToPetsc[lvl])[din];

      auto offsetKernel = [&](const IntVect& iv) -> void {
        for (int iphase = 0; iphase < m_numPhases; iphase++) {
          const PetscInt localRow = amrToPetsc(iv).getPetscRow(iphase);

          if (localRow >= 0) {
            amrToPetsc(iv).setPetscRow(iphase, m_localRowBegin + localRow);
          }
        }
      };

      BoxLoops::loop<D_DECL(1, 1, 1)>(dbl[din], offsetKernel);
    }

    m_amrToPetsc[lvl]->exchange();
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

  CH_assert(m_amrToPetsc.size() > m_finestLevel);
  CH_assert(m_levelGridsCoFi.size() > m_finestLevel);
  CH_assert(m_levelGridsFiCo.size() > m_finestLevel);

  const int numComp = 1;

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
PetscGrid::create(Vec& a_x) noexcept
{
  CH_TIME("PetscGrid::create");
  if (m_verbose) {
    pout() << "PetscGrid::create" << endl;
  }

  CH_assert(m_isDefined);

  PetscCallVoid(VecCreate(Chombo_MPI::comm, &a_x));
  PetscCallVoid(VecSetSizes(a_x, m_numLocalRows, m_numGlobalRows));
  PetscCallVoid(VecSetLocalToGlobalMapping(a_x, m_localToGlobalIS));
  PetscCallVoid(VecSetFromOptions(a_x));
}

void
PetscGrid::destroy(Vec& a_x) noexcept
{
  CH_TIME("PetscGrid::destroy");
  if (m_verbose) {
    pout() << "PetscGrid::destroy" << endl;
  }

  CH_assert(m_isDefined);

  PetscCallVoid(VecDestroy(&a_x));
}

void
PetscGrid::setValue(Vec& a_x, const PetscScalar a_value) const noexcept
{
  CH_TIME("PetscGrid::setValue");
  if (m_verbose) {
    pout() << "PetscGrid::setValue" << endl;
  }

  CH_assert(m_isDefined);

  PetscInt     startIdx = -1;
  PetscInt     endIdx   = -1;
  PetscScalar* arr      = nullptr;

  PetscCallVoid(VecGetOwnershipRange(a_x, &startIdx, &endIdx));
  PetscCallVoid(VecGetArray(a_x, &arr));

  CH_assert(endIdx >= startIdx);
  CH_assert(endIdx >= 0);
  CH_assert(startIdx >= 0);
  CH_assert(endIdx - startIdx == m_numLocalRows);

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
  CH_assert(a_y.size() > m_finestLevel);

  // Get the local PETSc vector and do basic error checks.
  PetscInt     startIdx = -1;
  PetscInt     endIdx   = -1;
  PetscScalar* arr      = nullptr;

  PetscCallVoid(VecGetOwnershipRange(a_x, &startIdx, &endIdx));
  PetscCallVoid(VecGetArray(a_x, &arr));

  CH_assert(endIdx >= startIdx);
  CH_assert(endIdx >= 0);
  CH_assert(startIdx >= 0);
  CH_assert(endIdx - startIdx == m_numLocalRows);

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
          const PetscInt globalRow = amrToPetsc(iv).getPetscRow(iphase);

          if (globalRow >= 0) {
            CH_assert(globalRow >= m_localRowBegin && globalRow < m_localRowBegin + m_numLocalRows);

            arr[globalRow - m_localRowBegin] = data(iv, 0);
          }
        };

        BoxLoops::loop<D_DECL(1, 1, 1)>(cellBox, kernel);
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
  CH_assert(a_y.size() > m_finestLevel);

  // Get the local PETSc vector and do basic error checks.
  PetscInt     startIdx = -1;
  PetscInt     endIdx   = -1;
  PetscScalar* arr      = nullptr;

  PetscCallVoid(VecGetOwnershipRange(a_x, &startIdx, &endIdx));
  PetscCallVoid(VecGetArray(a_x, &arr));

  CH_assert(endIdx >= startIdx);
  CH_assert(endIdx >= 0);
  CH_assert(startIdx >= 0);
  CH_assert(endIdx - startIdx == m_numLocalRows);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl  = m_levelGrids[lvl]->getGrids();
    const DataIterator&      dit  = dbl.dataIterator();
    const int                nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      for (int iphase = 0; iphase < m_numPhases; iphase++) {
        const Vector<PetscDOF>& petscToAMR = (*m_petscToAMR[iphase][lvl])[din];

        FArrayBox& data = (*a_y[lvl])[din].getPhase(iphase).getFArrayBox();

        for (int i = 0; i < petscToAMR.size(); i++) {
          CH_assert(petscToAMR[i].localRow >= 0 && petscToAMR[i].localRow < m_numLocalRows);

          data(petscToAMR[i].gridCell, 0) = arr[petscToAMR[i].localRow];
        }
      }
    }
  }

  PetscCallVoid(VecRestoreArray(a_x, &arr));
}

void
PetscGrid::dumpPetscGrid(const std::string& a_filename) const noexcept
{
  CH_TIME("PetscGrid::dumpPetscGrid");
  if (m_verbose) {
    pout() << "PetscGrid::dumpPetscGrid" << endl;
  }

  CH_assert(m_isDefined);

  const int numComp = 8;

  Vector<LevelData<FArrayBox>*> amrData(1 + m_finestLevel);
  Vector<DisjointBoxLayout>     amrGrids(1 + m_finestLevel);
  Vector<std::string>           varNames(numComp);

  varNames[0] = "valid_cell";
  varNames[1] = "coar_cf";
  varNames[2] = "fine_cf";
  varNames[3] = "ghost_cf";
  varNames[4] = "domain_cell";
  varNames[5] = "num_rows";
  varNames[6] = "petsc_row_phase0";
  varNames[7] = "petsc_row_phase1";

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl  = m_levelGrids[lvl]->getGrids();
    const DataIterator&      dit  = dbl.dataIterator();
    const int                nbox = dit.size();

    amrData[lvl]  = new LevelData<FArrayBox>(dbl, numComp, m_numGhost * IntVect::Unit);
    amrGrids[lvl] = dbl;

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      FArrayBox& data = (*amrData[lvl])[din];

      data.setVal(-1.0);

      const BaseFab<PetscAMRCell>& amrToPetsc = (*m_amrToPetsc[lvl])[din];

      auto kernel = [&](const IntVect& iv) -> void {
        const PetscAMRCell& amrCell = amrToPetsc(iv, 0);

        data(iv, 0) = amrCell.isCoveredByFinerGrid() ? 0.0 : 1.0;
        data(iv, 1) = amrCell.isCoarCF() ? 1.0 : 0.0;
        data(iv, 2) = amrCell.isFineCF() ? 1.0 : 0.0;
        data(iv, 3) = amrCell.isGhostCF() ? 1.0 : 0.0;
        data(iv, 4) = amrCell.isDomainBoundaryCell() ? 1.0 : 0.0;

        data(iv, 5) = 0.0;
        if (amrCell.getPetscRow(0) >= 0) {
          data(iv, 5) += 1.0;
          data(iv, 6) = amrCell.getPetscRow(0);
        }
        if (amrCell.getPetscRow(1) >= 0) {
          data(iv, 5) += 1.0;
          data(iv, 7) = amrCell.getPetscRow(1);
        }
      };

      BoxLoops::loop<D_DECL(1, 1, 1)>(data.box(), kernel);
    }
  }

#ifdef CH_USE_HDF5
  WriteAMRHierarchyHDF5(a_filename,
                        amrGrids,
                        amrData,
                        varNames,
                        m_levelGrids[0]->getDomain().domainBox(),
                        m_dx[0],
                        0.0,
                        0.0,
                        m_refRat,
                        1 + m_finestLevel);
#endif

  // Free memory
  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    delete amrData[lvl];
  }
}

#include <CD_NamespaceFooter.H>

#endif
