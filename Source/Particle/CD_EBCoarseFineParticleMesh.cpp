/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_EBCoarseFineParticleMesh.cpp
 * @brief  Implementation of CD_EBCoarseFineParticleMesh.H
 * @author Robert Marskar
 */

// Chombo includes
#include <CH_Timer.H>
#include <EBCellFactory.H>
#include <ParmParse.H>

// Our includes
#include <CD_EBCoarseFineParticleMesh.H>
#include <CD_EBAddOp.H>
#include <CD_BoxLoops.H>
#include <CD_NamespaceHeader.H>

constexpr int EBCoarseFineParticleMesh::m_comp;
constexpr int EBCoarseFineParticleMesh::m_nComp;

EBCoarseFineParticleMesh::EBCoarseFineParticleMesh() noexcept : m_isDefined(false), m_verbose(false)
{
  CH_TIME("EBCoarseFineParticleMesh::EBCoarseFineParticleMesh");
}

EBCoarseFineParticleMesh::EBCoarseFineParticleMesh(const EBLevelGrid& a_eblgCoar,
                                                   const EBLevelGrid& a_eblgFine,
                                                   const int          a_refRat,
                                                   const IntVect      a_ghost) noexcept
  : m_verbose(false)
{
  CH_TIME("EBCoarseFineParticleMesh::EBCoarseFineParticleMesh");

  this->define(a_eblgCoar, a_eblgFine, a_refRat, a_ghost);
}

EBCoarseFineParticleMesh::~EBCoarseFineParticleMesh() noexcept
{
  CH_TIME("EBCoarseFineParticleMesh::~EBCoarseFineParticleMesh");
}

void
EBCoarseFineParticleMesh::define(const EBLevelGrid& a_eblgCoar,
                                 const EBLevelGrid& a_eblgFine,
                                 const int          a_refRat,
                                 const IntVect      a_ghost) noexcept
{
  CH_TIME("EBCoarseFineParticleMesh::define");

  ParmParse pp("EBCoarseFineParticleMesh");
  pp.query("verbose", m_verbose);

  if (m_verbose) {
    pout() << "EBCoarseFineParticleMesh::define" << endl;
  }

  CH_assert(a_refRat % 2 == 0);
  CH_assert(a_refRat >= 2);

  m_eblgCoar = a_eblgCoar;
  m_eblgFine = a_eblgFine;
  m_refRat   = a_refRat;
  m_ghost    = a_ghost;

  // Make the coarsened fine grids.
  coarsen(m_eblgCoFi, m_eblgFine, m_refRat);

  // valid -> valid
  m_copierCoarToCoFi.define(m_eblgCoar.getDBL(), m_eblgCoFi.getDBL(), m_eblgCoar.getDomain(), false, IntVect::Zero);

  // valid+ghost -> valid
  m_copierCoFiToCoarIncludeGhosts.ghostDefine(m_eblgCoFi.getDBL(),
                                              m_eblgCoar.getDBL(),
                                              m_eblgCoar.getDomain(),
                                              m_ghost);
  this->defineVoFIterators();

  m_isDefined = true;
}

void
EBCoarseFineParticleMesh::defineVoFIterators() noexcept
{
  CH_TIME("EBCoarseFineParticleMesh::defineVoFIterators");
  if (m_verbose) {
    pout() << "EBCoarseFineParticleMesh::defineVoFIterators" << endl;
  }

  const DisjointBoxLayout& dblFine    = m_eblgFine.getDBL();
  const ProblemDomain&     domainFine = m_eblgFine.getDomain();
  const EBISLayout&        ebislFine  = m_eblgFine.getEBISL();

  const DisjointBoxLayout& dblCoFi   = m_eblgCoFi.getDBL();
  const EBISLayout&        ebislCoFi = m_eblgCoFi.getEBISL();

  m_vofIterFineGhosts.define(dblFine);
  m_vofIterCoFiGhosts.define(dblCoFi);
  m_vofIterFineCut.define(dblFine);

  const DataIterator& ditFine  = dblFine.dataIterator();
  const int           nboxFine = ditFine.size();

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nboxFine; mybox++) {
    const DataIndex& din = ditFine[mybox];

    const Box&     cellBoxFine = dblFine[din];
    const EBISBox& ebisBoxFine = ebislFine[din];
    const EBISBox& ebisBoxCoFi = ebislCoFi[din];

    const Box grownBoxFine = grow(cellBoxFine, m_ghost) & domainFine;
    const Box grownBoxCoFi = coarsen(grownBoxFine, m_refRat);

    // Coarse cells that the regular-cell sum in addFineGhostsToCoarse cannot handle: multi-valued coarse cells, and
    // coarse cells with a multi-valued cell among the fine ghost cells.
    IntVectSet multiCellsFine = ebisBoxFine.getMultiCells(grownBoxFine);
    multiCellsFine -= cellBoxFine;

    IntVectSet coarCells = ebisBoxCoFi.getMultiCells(grownBoxCoFi);
    coarCells |= coarsen(multiCellsFine, m_refRat);

    // The fine ghost cells in those coarse cells.
    IntVectSet fineCells = coarCells;
    fineCells.refine(m_refRat);
    fineCells &= grownBoxFine;
    fineCells -= cellBoxFine;

    m_vofIterFineGhosts[din].define(fineCells, ebisBoxFine.getEBGraph());
    m_vofIterCoFiGhosts[din].define(coarCells, ebisBoxCoFi.getEBGraph());

    // Valid fine cells in cut coarse cells.
    IntVectSet cutCells = ebisBoxCoFi.getIrregIVS(coarsen(cellBoxFine, m_refRat));
    cutCells.refine(m_refRat);
    cutCells &= cellBoxFine;

    m_vofIterFineCut[din].define(cutCells, ebisBoxFine.getEBGraph());
  }
}

void
EBCoarseFineParticleMesh::addFineGhostsToCoarse(LevelData<EBCellFAB>&       a_coarData,
                                                const LevelData<EBCellFAB>& a_fineData) const noexcept
{
  CH_TIME("EBCoarseFineParticleMesh::addFineGhostsToCoarse");
  if (m_verbose) {
    pout() << "EBCoarseFineParticleMesh::addFineGhostsToCoarse" << endl;
  }

  CH_assert(m_isDefined);
  CH_assert(a_coarData.nComp() == 1);
  CH_assert(a_fineData.nComp() == 1);
  CH_assert(m_refRat % 2 == 0);
  CH_assert(a_coarData.ghostVect() == m_ghost);
  CH_assert(a_fineData.ghostVect() == m_ghost);

  // clang-format off
  // TLDR: This routine will take the ghost cells that are in a_fineData and add them to the coarse level. We use buffers for this,
  //       copying the data from a_fineData to our nifty m_bufferFine buffer holder. The ghost cells in that scratch data are
  //       coarsened onto yet another buffer. Finally, we add the contents in that buffer to the coarse data.
  // clang-format on

  const DisjointBoxLayout& dblFine    = m_eblgFine.getDBL();
  const ProblemDomain&     domainFine = m_eblgFine.getDomain();
  const EBISLayout&        ebislFine  = m_eblgFine.getEBISL();

  const DisjointBoxLayout& dblCoFi = m_eblgCoFi.getDBL();

  const DataIterator& ditFine = dblFine.dataIterator();
  const DataIterator& ditCoFi = dblCoFi.dataIterator();

  const int nboxFine = ditFine.size();

  const Real factor = 1. / pow(m_refRat, SpaceDim);

  LevelData<EBCellFAB> bufferFine(m_eblgFine.getDBL(), m_nComp, m_ghost, EBCellFactory(m_eblgFine.getEBISL()));
  LevelData<EBCellFAB> bufferCoFi(m_eblgCoFi.getDBL(), m_nComp, m_ghost, EBCellFactory(m_eblgCoFi.getEBISL()));

  // Copy the fine data to scratch and reset the interior cells. We do this by copying everything to the fine scratch
  // data, and then set all the valid data in each box to zero. The exchange() operation will then take care of ghost
  // cells that. overlap with valid regions in different boxes (we could use a NeighborIterator to achieve the same).
  // Note that this is critical for the result because we essentially end up setting all valid data to zero so we don't
  // accidentally add mass to invalid region of the coarse grid (i.e., the region underneath the fine grid).
  a_fineData.localCopyTo(bufferFine);
#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nboxFine; mybox++) {
    const DataIndex& din = ditFine[mybox];

    const Box& box = dblFine[din];

    bufferFine[din].setVal(0.0, box, 0, m_nComp);
  }
  bufferFine.exchange();

  // Coarsen the fine grid data. We do this by AVERAGING the fine-grid data onto the coarse grid. This actually involves
  // entire patches and not just individual ghost cells regions, but that's ok because the rest of the data (in the
  // valid fine regions) are set to zero above, so the coarse data underneath those regions will be zero, also.
#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nboxFine; mybox++) {
    const DataIndex& din = ditFine[mybox];

    EBCellFAB&       coFiData = bufferCoFi[din];
    const EBCellFAB& fineData = bufferFine[din];

    coFiData.setVal(0.0);

    // Do all the regular cells.
    FArrayBox&       coFiDataReg = coFiData.getFArrayBox();
    const FArrayBox& fineDataReg = fineData.getFArrayBox();

    const Box fineBox = fineDataReg.box() & domainFine;
    for (int comp = 0; comp < m_nComp; comp++) {

      auto regularKernel = [&](const IntVect& ivFine) -> void {
        const IntVect ivCoar = coarsen(ivFine, m_refRat);

        coFiDataReg(ivCoar, comp) += fineDataReg(ivFine, comp) * factor;
      };

      // Not vectorizable: fine->coarse scatter via coarsen(ivFine) is a non-contiguous write
      // (m_refRat^D fine cells alias the same coarse cell). Covered fine cells hold no mass, so this sum is right for
      // every coarse cell that has one VoF and no multi-valued finer cell; the loops below redo the others.
      BoxLoops::loop<D_DECL(1, 1, 1)>(fineBox, regularKernel);
    }

    // Multi-valued cells: reset each coarse VoF and add the fine ghost VoFs that coarsen to it, with the same weight as
    // the regular cells.
    VoFIterator& vofitCoar = m_vofIterCoFiGhosts[din];
    VoFIterator& vofitFine = m_vofIterFineGhosts[din];

    auto resetCoar = [&](const VolIndex& coarVoF) -> void {
      for (int comp = 0; comp < m_nComp; comp++) {
        coFiData(coarVoF, comp) = 0.0;
      }
    };

    auto addFineToCoar = [&](const VolIndex& fineVoF) -> void {
      const VolIndex coarVoF = ebislFine.coarsen(fineVoF, m_refRat, din);

      for (int comp = 0; comp < m_nComp; comp++) {
        coFiData(coarVoF, comp) += fineData(fineVoF, comp) * factor;
      }
    };

    BoxLoops::loop(vofitCoar, resetCoar);
    BoxLoops::loop(vofitFine, addFineToCoar);
  }

  const Interval interv(0, m_nComp - 1);

  bufferCoFi.copyTo(interv, a_coarData, interv, m_copierCoFiToCoarIncludeGhosts, EBAddOp());
}

void
EBCoarseFineParticleMesh::addInvalidCoarseToFine(LevelData<EBCellFAB>&       a_fineData,
                                                 const LevelData<EBCellFAB>& a_coarData,
                                                 const bool                  a_cutCellsHoldDensity) const noexcept
{
  CH_TIME("EBCoarseFineParticleMesh::addInvalidCoarseToFine");
  if (m_verbose) {
    pout() << "EBCoarseFineParticleMesh::addInvalidCoarseToFine" << endl;
  }

  CH_assert(m_isDefined);
  CH_assert(a_fineData.nComp() == 1);
  CH_assert(a_coarData.nComp() == 1);
  CH_assert(a_fineData.ghostVect() == m_ghost);
  CH_assert(a_coarData.ghostVect() == m_ghost);

  // Piecewise constant interpolation of the density: every valid fine VoF gets the density of the coarse VoF it
  // coarsens to. Mass per cell volume, i.e. volume fraction times density, is scaled by kappaFine/kappaCoar in cut
  // cells; the fine volume fractions of a coarse VoF sum to refRat^SpaceDim times its own, so this conserves mass.
  // Cut cells that already hold the density take the coarse value as is. The coarse data is first copied onto the
  // coarsened fine grids, so each fine box reads its coarse values locally.
  const DisjointBoxLayout& dblFine   = m_eblgFine.getDBL();
  const EBISLayout&        ebislFine = m_eblgFine.getEBISL();
  const EBISLayout&        ebislCoFi = m_eblgCoFi.getEBISL();

  const Interval interv(0, m_nComp - 1);

  LevelData<EBCellFAB> bufferCoFi(m_eblgCoFi.getDBL(), m_nComp, IntVect::Zero, EBCellFactory(m_eblgCoFi.getEBISL()));

  a_coarData.copyTo(interv, bufferCoFi, interv, m_copierCoarToCoFi);

  const DataIterator& dit  = dblFine.dataIterator();
  const int           nbox = dit.size();

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    const Box        fineBox     = dblFine[din];
    const EBCellFAB& coFiData    = bufferCoFi[din];
    const EBISBox&   ebisBoxFine = ebislFine[din];
    const EBISBox&   ebisBoxCoFi = ebislCoFi[din];

    EBCellFAB interpData(ebislFine[din], fineBox, m_nComp);

    BaseFab<Real>&       interpDataReg = interpData.getSingleValuedFAB();
    const BaseFab<Real>& coFiDataReg   = coFiData.getSingleValuedFAB();

    // Not vectorizable: coarse->fine broadcast via coarsen(ivFine) is a strided read. Cells in cut coarse cells are
    // overwritten per VoF below.
    auto regularKernel = [&](const IntVect& ivFine) -> void {
      interpDataReg(ivFine, m_comp) = coFiDataReg(coarsen(ivFine, m_refRat), m_comp);
    };

    auto cutKernel = [&](const VolIndex& fineVoF) -> void {
      const VolIndex coarVoF = ebislFine.coarsen(fineVoF, m_refRat, din);

      // A coarse VoF without volume has children without volume, which take nothing.
      const Real coarVolFrac = ebisBoxCoFi.volFrac(coarVoF);

      Real scale = 1.0;

      if (!a_cutCellsHoldDensity) {
        scale = (coarVolFrac > 0.0) ? ebisBoxFine.volFrac(fineVoF) / coarVolFrac : 0.0;
      }

      interpData(fineVoF, m_comp) = coFiData(coarVoF, m_comp) * scale;
    };

    BoxLoops::loop<D_DECL(1, 1, 1)>(fineBox, regularKernel);
    BoxLoops::loop(m_vofIterFineCut[din], cutKernel);

    a_fineData[din].plus(interpData, fineBox, 0, 0, m_nComp);
  }
}

#include <CD_NamespaceFooter.H>
