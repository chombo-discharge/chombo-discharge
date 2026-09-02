/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_PhaseRealm.cpp
 * @brief  Implementation of CD_PhaseRealm.H
 * @author Robert Marskar
 */

// Std includes
#include <cmath>
#include <limits>

// Chombo includes
#include <EBArith.H>
#include <ParmParse.H>
#include <BoxIterator.H>
#include <BaseIVFactory.H>
#include <EBCellFactory.H>
#include <EBFluxFactory.H>

// Our includes
#include <CD_PhaseRealm.H>
#include <CD_Timer.H>
#include <CD_LoadBalancing.H>
#include <CD_EBLeastSquaresMultigridInterpolator.H>
#include <CD_MemoryReport.H>
#include <CD_BoxLoops.H>
#include <CD_DataOps.H>
#include <CD_Location.H>
#include <CD_ParallelOps.H>
#include <CD_MirrorDeposition.H>
#include <CD_NamespaceHeader.H>

PhaseRealm::PhaseRealm() : m_isDefined(false), m_profile(false), m_verbose(false)
{
  CH_TIME("PhaseRealm::PhaseRealm");

  // Default settings

  this->registerOperator(s_eb_gradient);
  this->registerOperator(s_eb_irreg_interp);

  // Adding this for debugging purposes.
  ParmParse pp("PhaseRealm");
  pp.query("profile", m_profile);
  pp.query("verbosity", m_verbose);
}

PhaseRealm::~PhaseRealm() = default;

void
PhaseRealm::define(const Vector<DisjointBoxLayout>&      a_grids,
                   const Vector<ProblemDomain>&          a_domains,
                   const Vector<int>&                    a_refRat,
                   const Vector<Real>&                   a_dx,
                   const RealVect&                       a_probLo,
                   const int                             a_finestLevel,
                   const int                             a_ebGhost,
                   const int                             a_numGhost,
                   const int                             a_lsfGhost,
                   const int                             a_redistRad,
                   const int                             a_mgInterpOrder,
                   const int                             a_mgInterpRadius,
                   const int                             a_mgInterpWeight,
                   const int                             a_mirrorBandRadius,
                   const Real                            a_mirrorCondTol,
                   const Real                            a_mirrorCreaseTol,
                   const CellCentroidInterpolation::Type a_centroidStencil,
                   const EBCentroidInterpolation::Type   a_ebStencil,
                   const RefCountedPtr<BaseIF>&          a_baseif,
                   const RefCountedPtr<EBIndexSpace>&    a_ebis)
{
  CH_TIME("PhaseRealm::define");

  m_grids                         = a_grids;
  m_domains                       = a_domains;
  m_refinementRatios              = a_refRat;
  m_dx                            = a_dx;
  m_probLo                        = a_probLo;
  m_finestLevel                   = a_finestLevel;
  m_numEbGhostsCells              = a_ebGhost;
  m_numGhostCells                 = a_numGhost;
  m_numLsfGhostCells              = a_lsfGhost;
  m_redistributionRadius          = a_redistRad;
  m_mirrorBandRadius              = a_mirrorBandRadius;
  m_mirrorConditionTolerance      = a_mirrorCondTol;
  m_mirrorCreaseTolerance         = a_mirrorCreaseTol;
  m_multigridInterpolationOrder   = a_mgInterpOrder;
  m_multigridInterpolationRadius  = a_mgInterpRadius;
  m_multigridInterpolationWeight  = a_mgInterpWeight;
  m_cellCentroidInterpolationType = a_centroidStencil;
  m_ebCentroidInterpolationType   = a_ebStencil;
  m_baseif                        = a_baseif;
  m_ebis                          = a_ebis;

  if (!m_ebis.isNull()) {
    m_isDefined = true;
  }

  m_hasEbCf = true;

  // Adding this for debugging purposes.
  ParmParse pp("PhaseRealm");
  pp.query("profile", m_profile);
  pp.query("verbosity", m_verbose);
}

void
PhaseRealm::setGrids(const Vector<DisjointBoxLayout>& a_grids, const int a_finestLevel)
{
  CH_TIME("PhaseRealm::setGrids");
  if (m_verbose) {
    pout() << "PhaseRealm::setGrids" << endl;
  }

  if (m_isDefined) {
    m_grids       = a_grids;
    m_finestLevel = a_finestLevel;
  }
}

void
PhaseRealm::preRegrid()
{
  CH_TIME("PhaseRealm::preRegrid");
  if (m_verbose) {
    pout() << "PhaseRealm::preRegrid" << endl;
  }

  m_grids.resize(0);
  m_ebisl.resize(0);
  m_eblg.resize(0);
  m_eblgCoFi.resize(0);
  m_eblgFiCo.resize(0);
  m_vofIter.resize(0);
  m_multiCutVofIter.resize(0);
  m_faceIter.resize(0);
  m_faceIterNoBoundary.resize(0);
  m_faceIterTanGhost.resize(0);
  m_multiCutFaceIter.resize(0);
  m_coarAve.resize(0);
  m_multigridInterpolator.resize(0);
  m_ebFineInterp.resize(0);
  m_ebReflux.resize(0);
  m_redistributionOp.resize(0);
  m_gradientOp.resize(0);
  m_levelset.resize(0);
  m_regularCells.resize(0);
  m_coveredCells.resize(0);
  m_notCoveredCells.resize(0);
  m_irregularCells.resize(0);
  m_cellCentroidInterpolation.resize(0);
  m_ebCentroidInterpolation.resize(0);
  m_nonConservativeDivergence.resize(0);
}

void
PhaseRealm::regridBase(const int a_lmin)
{
  CH_TIME("PhaseRealm::regridBase");
  if (m_verbose) {
    pout() << "PhaseRealm::regridBase" << endl;
  }

  if (m_isDefined) {

    Timer timer("PhaseRealm::regridBase(int)");

    if (m_profile) {
      pout() << "before/after levelgrid define" << endl;
      MemoryReport::getMaxMinMemoryUsage();
    }
    timer.startEvent("Define EBLevelGrid");
    this->defineEBLevelGrid(a_lmin);
    timer.stopEvent("Define EBLevelGrid");
    if (m_profile) {
      MemoryReport::getMaxMinMemoryUsage();
      pout() << endl;
    }

    if (m_profile) {
      pout() << "before/after vofit define" << endl;
      MemoryReport::getMaxMinMemoryUsage();
    }
    timer.startEvent("Define VoFIterators");
    this->defineVofIterator(a_lmin);
    timer.stopEvent("Define VoFIterators");
    if (m_profile) {
      MemoryReport::getMaxMinMemoryUsage();
      pout() << endl;
    }

    // Built here (rather than in regridOperators) so the cell masks are available to load-balancing
    // routines, which run after regridBase but before regridOperators.
    if (m_profile) {
      pout() << "before/after masks define" << endl;
      MemoryReport::getMaxMinMemoryUsage();
    }
    timer.startEvent("Cell masks");
    this->defineMasks(a_lmin, m_numGhostCells);
    timer.stopEvent("Cell masks");
    if (m_profile) {
      MemoryReport::getMaxMinMemoryUsage();
      pout() << endl;
    }

    if (m_profile) {
      timer.eventReport(pout());
    }
  }
}

void
PhaseRealm::regridOperators(const int a_lmin)
{
  CH_TIME("PhaseRealm::regridOperators");
  if (m_verbose) {
    pout() << "PhaseRealm::regridOperators" << endl;
  }

  if (m_isDefined) {

    Timer timer("PhaseRealm::regridOperators(int)");

    if (m_profile) {
      pout() << "before/after coarave define" << endl;
      MemoryReport::getMaxMinMemoryUsage();
    }
    timer.startEvent("EBCoarAve");
    this->defineEBCoarAve(a_lmin);
    timer.stopEvent("EBCoarAve");
    if (m_profile) {
      MemoryReport::getMaxMinMemoryUsage();
      pout() << endl;
    }

    if (m_profile) {
      pout() << "before/after mg define" << endl;
      MemoryReport::getMaxMinMemoryUsage();
    }
    timer.startEvent("Multigrid interpolator");
    this->defineEBMultigrid(a_lmin);
    timer.stopEvent("Multigrid interpolator");
    if (m_profile) {
      MemoryReport::getMaxMinMemoryUsage();
      pout() << endl;
    }

    if (m_profile) {
      pout() << "before/after fillpatch define" << endl;
      MemoryReport::getMaxMinMemoryUsage();
    }
    timer.startEvent("Ghost interp");
    this->defineFillPatch(a_lmin);
    timer.stopEvent("Ghost interp");
    if (m_profile) {
      MemoryReport::getMaxMinMemoryUsage();
      pout() << endl;
    }

    if (m_profile) {
      pout() << "before/after pwlinterp define" << endl;
      MemoryReport::getMaxMinMemoryUsage();
    }
    timer.startEvent("PWL interp");
    this->defineEBCoarseToFineInterp(a_lmin);
    timer.stopEvent("PWL interp");
    if (m_profile) {
      MemoryReport::getMaxMinMemoryUsage();
      pout() << endl;
    }

    if (m_profile) {
      pout() << "before/after fluxreg define" << endl;
      MemoryReport::getMaxMinMemoryUsage();
    }
    timer.startEvent("Flux register");
    this->defineFluxReg(a_lmin, 1);
    timer.stopEvent("Flux register");
    if (m_profile) {
      MemoryReport::getMaxMinMemoryUsage();
      pout() << endl;
    }

    if (m_profile) {
      pout() << "before/after redist define" << endl;
      MemoryReport::getMaxMinMemoryUsage();
    }
    timer.startEvent("EB redist");
    this->defineRedistOper(a_lmin, 1);
    timer.stopEvent("EB redist");
    if (m_profile) {
      MemoryReport::getMaxMinMemoryUsage();
      pout() << endl;
    }

    if (m_profile) {
      pout() << "before/after gradsten define" << endl;
      MemoryReport::getMaxMinMemoryUsage();
    }
    timer.startEvent("Gradient stencil");
    this->defineGradSten(a_lmin);
    timer.stopEvent("Gradient stencil");
    if (m_profile) {
      MemoryReport::getMaxMinMemoryUsage();
      pout() << endl;
    }

    if (m_profile) {
      pout() << "before/after irregsten define" << endl;
      MemoryReport::getMaxMinMemoryUsage();
    }
    timer.startEvent("Irreg stencil");
    this->defineIrregSten();
    timer.stopEvent("Irreg stencil");
    if (m_profile) {
      MemoryReport::getMaxMinMemoryUsage();
      pout() << endl;
    }

    if (m_profile) {
      pout() << "before/after nonconsdivsten define" << endl;
      MemoryReport::getMaxMinMemoryUsage();
    }
    timer.startEvent("Non-conservative stencil");
    this->defineNonConservativeDivergence(a_lmin);
    timer.stopEvent("Non-conservative stencil");
    if (m_profile) {
      MemoryReport::getMaxMinMemoryUsage();
      pout() << endl;
    }

    if (m_profile) {
      pout() << "before/after particlemesh define" << endl;
      MemoryReport::getMaxMinMemoryUsage();
    }
    timer.startEvent("Particle-mesh");
    this->defineParticleMesh();
    timer.stopEvent("Particle-mesh");
    if (m_profile) {
      MemoryReport::getMaxMinMemoryUsage();
      pout() << endl;
    }

    if (m_profile) {
      pout() << "before/after mirror surface data define" << endl;
      MemoryReport::getMaxMinMemoryUsage();
    }

    timer.startEvent("Mirror surface data");
    this->defineMirrorSurfaceData(a_lmin);
    timer.stopEvent("Mirror surface data");

    if (m_profile) {
      MemoryReport::getMaxMinMemoryUsage();
      pout() << endl;
    }

    if (m_profile) {
      pout() << "before/after levelset define" << endl;
      MemoryReport::getMaxMinMemoryUsage();
    }
    timer.startEvent("Levelset");
    this->defineLevelSet(a_lmin, m_numLsfGhostCells);
    timer.stopEvent("Levelset");
    if (m_profile) {
      MemoryReport::getMaxMinMemoryUsage();
      pout() << endl;
    }

    if (m_profile) {
      timer.eventReport(pout());
    }
  }
}

void
PhaseRealm::registerOperator(const std::string& a_operator)
{
  CH_TIME("PhaseRealm::registerOperator");
  if (m_verbose) {
    pout() << "PhaseRealm::registerOperator" << endl;
  }

  // These are the supported operators - issue an error if we ask for something that is not supported.
  if (!(a_operator == s_eb_coar_ave || a_operator == s_eb_fill_patch || a_operator == s_eb_fine_interp ||
        a_operator == s_eb_flux_reg || a_operator == s_eb_redist || a_operator == s_noncons_div ||
        a_operator == s_eb_gradient || a_operator == s_particle_mesh || a_operator == s_eb_irreg_interp ||
        a_operator == s_eb_multigrid || a_operator == s_levelset || a_operator == s_mirror_deposition)) {

    const std::string str = "PhaseRealm::registerOperator - unknown operator '" + a_operator + "' requested";
    MayDay::Error(str.c_str());
  }

  if (!this->queryOperator(a_operator)) {
    m_operatorMap.emplace(a_operator, true);
  }
}

bool
PhaseRealm::queryOperator(const std::string& a_operator) const
{
  CH_TIME("PhaseRealm::queryOperator");
  if (m_verbose) {
    pout() << "PhaseRealm::queryOperator" << endl;
  }

  bool ret = false;

  if (m_isDefined) {
    ret = true;

    if (m_operatorMap.find(a_operator) == m_operatorMap.end()) {
      ret = false;
    }
  }

  return ret;
}

void
PhaseRealm::defineEBLevelGrid(const int a_lmin)
{
  CH_TIME("PhaseRealm::defineEBLevelGrid");
  if (m_verbose) {
    pout() << "PhaseRealm::defineEBLevelGrid" << endl;
  }

  m_eblg.resize(1 + m_finestLevel);
  m_eblgCoFi.resize(1 + m_finestLevel);
  m_eblgFiCo.resize(1 + m_finestLevel);
  m_ebisl.resize(1 + m_finestLevel);

  for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++) {
    m_eblg[lvl] = RefCountedPtr<EBLevelGrid>(
      new EBLevelGrid(m_grids[lvl], m_domains[lvl], m_numEbGhostsCells, &(*m_ebis)));

    const bool hasCoar = lvl > 0;
    const bool hasFine = lvl < m_finestLevel;

    if (hasCoar) {
      m_eblg[lvl]->setMaxCoarseningRatio(m_refinementRatios[lvl - 1], &(*m_ebis));
    }

    if (hasFine) {
      m_eblg[lvl]->setMaxRefinementRatio(m_refinementRatios[lvl]);
    }

    m_ebisl[lvl] = m_eblg[lvl]->getEBISL();

    // Define the coarsened grids.
    if (lvl > 0) {
      m_eblgCoFi[lvl - 1] = RefCountedPtr<EBLevelGrid>(new EBLevelGrid());

      coarsen(*m_eblgCoFi[lvl - 1], *m_eblg[lvl], m_refinementRatios[lvl - 1]);
      m_eblgCoFi[lvl - 1]->getEBISL().setMaxRefinementRatio(m_refinementRatios[lvl - 1], m_eblg[lvl]->getEBIS());
    }

    // Define the refined grids. Here m_eblgFiCo contains grids on level lvl
    if (lvl < m_finestLevel) {
      m_eblgFiCo[lvl + 1] = RefCountedPtr<EBLevelGrid>(new EBLevelGrid());

      refine(*m_eblgFiCo[lvl + 1], *m_eblg[lvl], m_refinementRatios[lvl]);
      m_eblgFiCo[lvl + 1]->getEBISL().setMaxCoarseningRatio(m_refinementRatios[lvl], m_eblg[lvl]->getEBIS());
    }
  }
}

void
PhaseRealm::defineVofIterator(const int a_lmin)
{
  CH_TIME("PhaseRealm::defineVofIterator");
  if (m_verbose) {
    pout() << "PhaseRealm::defineVofIterator" << endl;
  }

  m_vofIter.resize(1 + m_finestLevel);
  m_multiCutVofIter.resize(1 + m_finestLevel);
  m_faceIter.resize(1 + m_finestLevel);
  m_faceIterNoBoundary.resize(1 + m_finestLevel);
  m_faceIterTanGhost.resize(1 + m_finestLevel);
  m_multiCutFaceIter.resize(1 + m_finestLevel);

  for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++) {

    m_vofIter[lvl]         = RefCountedPtr<LayoutData<VoFIterator>>(new LayoutData<VoFIterator>(m_grids[lvl]));
    m_multiCutVofIter[lvl] = RefCountedPtr<LayoutData<VoFIterator>>(new LayoutData<VoFIterator>(m_grids[lvl]));
    m_faceIter[lvl]        = RefCountedPtr<LayoutData<std::array<FaceIterator, SpaceDim>>>(
      new LayoutData<std::array<FaceIterator, SpaceDim>>(m_grids[lvl]));
    m_faceIterNoBoundary[lvl] = RefCountedPtr<LayoutData<std::array<FaceIterator, SpaceDim>>>(
      new LayoutData<std::array<FaceIterator, SpaceDim>>(m_grids[lvl]));
    m_faceIterTanGhost[lvl] = RefCountedPtr<LayoutData<std::array<FaceIterator, SpaceDim>>>(
      new LayoutData<std::array<FaceIterator, SpaceDim>>(m_grids[lvl]));
    m_multiCutFaceIter[lvl] = RefCountedPtr<LayoutData<std::array<FaceIterator, SpaceDim>>>(
      new LayoutData<std::array<FaceIterator, SpaceDim>>(m_grids[lvl]));

    const DisjointBoxLayout& dbl = m_grids[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const int nbox = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      VoFIterator& vofit         = (*m_vofIter[lvl])[din];
      VoFIterator& multiCutVofit = (*m_multiCutVofIter[lvl])[din];

      const Box&        cellBox  = m_grids[lvl].get(din);
      const EBISBox&    ebisbox  = m_ebisl[lvl][din];
      const EBGraph&    ebgraph  = ebisbox.getEBGraph();
      const IntVectSet& irreg    = ebisbox.getIrregIVS(cellBox);
      const IntVectSet& multiCut = ebisbox.getMultiCells(cellBox);

      vofit.define(irreg, ebgraph);
      multiCutVofit.define(multiCut, ebgraph);

      // Face iterators over cut-cell faces in the valid box, with and without domain boundary faces.
      std::array<FaceIterator, SpaceDim>& faceIter           = (*m_faceIter[lvl])[din];
      std::array<FaceIterator, SpaceDim>& faceIterNoBoundary = (*m_faceIterNoBoundary[lvl])[din];
      for (int dir = 0; dir < SpaceDim; dir++) {
        faceIter[dir].define(irreg, ebgraph, dir, FaceStop::SurroundingWithBoundary);
        faceIterNoBoundary[dir].define(irreg, ebgraph, dir, FaceStop::SurroundingNoBoundary);
      }

      // Per-direction tangential-ghost face iterators. For faces in direction dir the IVS covers
      // irregular cells in a box grown by 1 in the two tangential directions only; dir itself is
      // kept at the valid-box extent. With FaceStop::SurroundingNoBoundary, this gives exactly
      // the cut-cell faces that reach into the first tangential ghost layer while staying within
      // 1 cell of the valid box in the normal (dir) direction.
      std::array<FaceIterator, SpaceDim>& faceIterTanGhost = (*m_faceIterTanGhost[lvl])[din];
      for (int dir = 0; dir < SpaceDim; dir++) {
        Box tanBox = cellBox;
        for (int tanDir = 0; tanDir < SpaceDim; tanDir++) {
          if (tanDir != dir) {
            tanBox.grow(tanDir, 1);
          }
        }
        tanBox &= m_domains[lvl];
        const IntVectSet irregTan = ebisbox.getIrregIVS(tanBox);
        faceIterTanGhost[dir].define(irregTan, ebgraph, dir, FaceStop::SurroundingNoBoundary);
      }

      // Multi-cut face iterators over the valid box (SurroundingNoBoundary). Used as a second pass
      // after box loops over getSingleValuedFAB to avoid double-processing singly-cut faces.
      std::array<FaceIterator, SpaceDim>& multiCutFaceIter = (*m_multiCutFaceIter[lvl])[din];
      for (int dir = 0; dir < SpaceDim; dir++) {
        multiCutFaceIter[dir].define(multiCut, ebgraph, dir, FaceStop::SurroundingNoBoundary);
      }
    }
  }
}

void
PhaseRealm::defineLevelSet(const int a_lmin, const int a_numGhost)
{
  CH_TIME("PhaseRealm::defineLevelSet");
  if (m_verbose) {
    pout() << "PhaseRealm::defineLevelset" << endl;
  }

  constexpr Real minVal = std::numeric_limits<Real>::min();

  const bool doThisOperator = this->queryOperator(s_levelset);

  m_levelset.resize(1 + m_finestLevel);

  if (doThisOperator) {

    const int comp  = 0;
    const int ncomp = 1;

    for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++) {
      const Real dx = m_dx[lvl];

      const DisjointBoxLayout& dbl = m_grids[lvl];
      const DataIterator&      dit = dbl.dataIterator();

      m_levelset[lvl] = RefCountedPtr<LevelData<FArrayBox>>(
        new LevelData<FArrayBox>(dbl, ncomp, a_numGhost * IntVect::Unit));

      const int nbox = dit.size();
#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din = dit[mybox];

        FArrayBox& fab = (*m_levelset[lvl])[din];
        const Box  bx  = fab.box();

        if (!m_baseif.isNull()) {
          // Not auto-vectorizable: m_baseif->value(pos) is a virtual call on the polymorphic implicit
          // function, evaluated per cell. This is a one-time (regrid) setup loop.
          auto kernel = [&](const IntVect& iv) -> void {
            const RealVect pos = m_probLo + (0.5 * RealVect::Unit + RealVect(iv)) * dx;

            fab(iv, comp) = m_baseif->value(pos);
          };

          BoxLoops::loop<D_DECL(1, 1, 1)>(bx, kernel);
        }
        else {
          fab.setVal(minVal, comp);
        }
      }
    }
  }
}

void
PhaseRealm::defineMasks(const int a_lmin, const int a_numGhost)
{
  CH_TIME("PhaseRealm::defineMasks");
  if (m_verbose) {
    pout() << "PhaseRealm::defineMasks" << endl;
  }

  constexpr int comp  = 0;
  constexpr int ncomp = 1;

  const IntVect ghost = a_numGhost * IntVect::Unit;

  m_regularCells.resize(1 + m_finestLevel);
  m_coveredCells.resize(1 + m_finestLevel);
  m_notCoveredCells.resize(1 + m_finestLevel);
  m_irregularCells.resize(1 + m_finestLevel);

  for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl   = m_grids[lvl];
    const EBISLayout&        ebisl = m_ebisl[lvl];
    const DataIterator&      dit   = dbl.dataIterator();

    m_regularCells[lvl] = RefCountedPtr<LevelData<EBCellFAB>>(
      new LevelData<EBCellFAB>(dbl, ncomp, ghost, EBCellFactory(ebisl)));
    m_coveredCells[lvl] = RefCountedPtr<LevelData<EBCellFAB>>(
      new LevelData<EBCellFAB>(dbl, ncomp, ghost, EBCellFactory(ebisl)));
    m_notCoveredCells[lvl] = RefCountedPtr<LevelData<EBCellFAB>>(
      new LevelData<EBCellFAB>(dbl, ncomp, ghost, EBCellFactory(ebisl)));
    m_irregularCells[lvl] = RefCountedPtr<LevelData<EBCellFAB>>(
      new LevelData<EBCellFAB>(dbl, ncomp, ghost, EBCellFactory(ebisl)));

    const int nbox = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      EBCellFAB& regular    = (*m_regularCells[lvl])[din];
      EBCellFAB& covered    = (*m_coveredCells[lvl])[din];
      EBCellFAB& notCovered = (*m_notCoveredCells[lvl])[din];
      EBCellFAB& irregular  = (*m_irregularCells[lvl])[din];

      const EBISBox& ebisbox = regular.getEBISBox();

      // Regular mask = 1 in regular cells, 0 in covered/irregular cells.
      regular.setVal(1.0);
      regular.setCoveredCellVal(0.0, comp);

      // Covered mask = 1 in covered cells, 0 elsewhere.
      covered.setVal(0.0);
      covered.setCoveredCellVal(1.0, comp);

      // Non-covered mask = 1 in regular/irregular cells, 0 in covered cells.
      notCovered.setVal(1.0);
      notCovered.setCoveredCellVal(0.0, comp);

      // Irregular mask = 1 in irregular cells, 0 elsewhere (raised below).
      irregular.setVal(0.0);

      // Lower the regular mask and raise the irregular mask on the irregular cells of the valid box.
      // Write the single-valued FABs directly since the masks are read as regular-grid data.
      BaseFab<Real>&   regularReg   = regular.getSingleValuedFAB();
      BaseFab<Real>&   irregularReg = irregular.getSingleValuedFAB();
      const IntVectSet irregIVS     = ebisbox.getIrregIVS(dbl[din]);

      for (IVSIterator ivsIt(irregIVS); ivsIt.ok(); ++ivsIt) {
        const IntVect& iv = ivsIt();

        regularReg(iv, comp)   = 0.0;
        irregularReg(iv, comp) = 1.0;
      }
    }

    // Fill ghost cells so the masks are consistent across grid patches.
    m_regularCells[lvl]->exchange();
    m_coveredCells[lvl]->exchange();
    m_notCoveredCells[lvl]->exchange();
    m_irregularCells[lvl]->exchange();
  }
}

void
PhaseRealm::defineEBCoarAve(const int a_lmin)
{
  CH_TIME("PhaseRealm::defineEBCoarAve");
  if (m_verbose) {
    pout() << "PhaseRealm::defineEBCoarAve" << endl;
  }

  const bool doThisOperator = this->queryOperator(s_eb_coar_ave);

  m_coarAve.resize(1 + m_finestLevel);

  if (doThisOperator) {
    for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++) {

      const bool hasCoar = lvl > 0;

      if (hasCoar) {
        m_coarAve[lvl] = RefCountedPtr<EBCoarAve>(
          new EBCoarAve(*m_eblg[lvl], *m_eblg[lvl - 1], *m_eblgCoFi[lvl - 1], m_refinementRatios[lvl - 1]));
      }
    }
  }
}

void
PhaseRealm::defineEBMultigrid(const int a_lmin)
{
  CH_TIME("PhaseRealm::defineEBMultigrid");
  if (m_verbose) {
    pout() << "PhaseRealm::defineEBMultigrid" << endl;
  }

  const bool doThisOperator = this->queryOperator(s_eb_multigrid);

  m_multigridInterpolator.resize(1 + m_finestLevel);

  if (doThisOperator) {

    for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++) {

      const bool hasCoar = lvl > 0;

      // Interpolator for ghost cells on level l is stored on level l.
      if (hasCoar) {
        const EBLevelGrid& eblgFine = *m_eblg[lvl];
        const EBLevelGrid& eblgCoFi = *m_eblgCoFi[lvl - 1];
        const EBLevelGrid& eblgCoar = *m_eblg[lvl - 1];

        m_multigridInterpolator[lvl] = RefCountedPtr<EBMultigridInterpolator>(
          new EBLeastSquaresMultigridInterpolator(eblgFine,
                                                  eblgCoFi,
                                                  eblgCoar,
                                                  Location::Cell::Center,
                                                  m_numGhostCells * IntVect::Unit,
                                                  m_refinementRatios[lvl - 1],
                                                  m_multigridInterpolationRadius,
                                                  m_multigridInterpolationOrder,
                                                  m_multigridInterpolationWeight));
      }
    }
  }
}

void
PhaseRealm::defineFillPatch(const int a_lmin)
{
  CH_TIME("PhaseRealm::defineFillPatch");
  if (m_verbose) {
    pout() << "PhaseRealm::defineFillPatch" << endl;
  }

  const bool doThisOperator = this->queryOperator(s_eb_fill_patch);

  m_ghostCellInterpolator.resize(1 + m_finestLevel);

  if (doThisOperator) {

    const int     radius = m_numGhostCells;
    const IntVect ghost  = m_numGhostCells * IntVect::Unit;

    for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++) {
      const bool hasCoar = lvl > 0;

      // Filling ghost cells on level l from coarse data on level l-1 is stored on level l
      if (hasCoar) {
        m_ghostCellInterpolator[lvl] = RefCountedPtr<EBGhostCellInterpolator>(
          new EBGhostCellInterpolator(*m_eblg[lvl],
                                      *m_eblgCoFi[lvl - 1],
                                      *m_eblg[lvl - 1],
                                      ghost,
                                      m_refinementRatios[lvl - 1],
                                      radius));
      }
    }
  }
}

void
PhaseRealm::defineEBCoarseToFineInterp(const int a_lmin)
{
  CH_TIME("PhaseRealm::defineEBCoarseToFineInterp");
  if (m_verbose) {
    pout() << "PhaseRealm::defineEBCoarseToFineInterp" << endl;
  }

  const bool doThisOperator = this->queryOperator(s_eb_fine_interp);

  m_ebFineInterp.resize(1 + m_finestLevel);

  if (doThisOperator) {
    for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++) {

      const bool hasCoar = lvl > 0;

      // Interpolator for filling data on level l from level l-1 lives on level l
      if (hasCoar) {
        m_ebFineInterp[lvl] = RefCountedPtr<EBCoarseToFineInterp>(
          new EBCoarseToFineInterp(*m_eblg[lvl], *m_eblgCoFi[lvl - 1], *m_eblg[lvl - 1], m_refinementRatios[lvl - 1]));
      }
    }
  }
}

void
PhaseRealm::defineFluxReg(const int a_lmin, const int /*a_regsize*/)
{
  CH_TIME("PhaseRealm::defineFluxReg");
  if (m_verbose) {
    pout() << "PhaseRealm::defineFluxReg" << endl;
  }

  ParmParse pp("PhaseRealm");
  bool      useChomboFluxRegister = false;
  pp.query("use_chombo_fluxreg", useChomboFluxRegister);

  const bool doThisOperator = this->queryOperator(s_eb_flux_reg);

  m_ebReflux.resize(1 + m_finestLevel);

  if (doThisOperator) {

    for (int lvl = std::max(0, a_lmin - 1); lvl <= m_finestLevel; lvl++) {

      const bool hasFine = lvl < m_finestLevel;

      // Flux register for matching fluxes between level l and level l+1 (the fine level) lives on level l
      if (hasFine) {
        m_ebReflux[lvl] = RefCountedPtr<EBReflux>(
          new EBReflux(*m_eblg[lvl], *m_eblg[lvl + 1], *m_eblgCoFi[lvl], m_refinementRatios[lvl]));
      }
    }
  }
}

void
PhaseRealm::defineRedistOper(const int a_lmin, const int /*a_regsize*/)
{
  CH_TIME("PhaseRealm::defineRedistOper");
  if (m_verbose) {
    pout() << "PhaseRealm::defineRedistOper" << endl;
  }

  // TLDR: The operator for redistributing on level l lives on level l.

  const bool doThisOperator = this->queryOperator(s_eb_redist);

  m_redistributionOp.resize(1 + m_finestLevel);

  if (doThisOperator) {

    for (int lvl = std::max(0, a_lmin - 2); lvl <= m_finestLevel; lvl++) {
      const bool hasCoar = lvl > 0;
      const bool hasFine = lvl < m_finestLevel;

      int refToCoar = -1;
      int refToFine = -1;

      EBLevelGrid eblgCoar;
      EBLevelGrid eblgCoarsened;
      EBLevelGrid eblg;
      EBLevelGrid eblgRefined;
      EBLevelGrid eblgFine;

      if (hasCoar) {
        eblgCoar      = *m_eblg[lvl - 1];
        eblgCoarsened = *m_eblgCoFi[lvl - 1];

        refToCoar = m_refinementRatios[lvl - 1];
      }

      eblg = *m_eblg[lvl];

      if (hasFine) {
        eblgRefined = *m_eblgFiCo[lvl + 1];
        eblgFine    = *m_eblg[lvl + 1];

        refToFine = m_refinementRatios[lvl];
      }

      const bool redistributeOutside = true;

      m_redistributionOp[lvl] = RefCountedPtr<EBFluxRedistribution>(new EBFluxRedistribution(eblgCoar,
                                                                                             eblgCoarsened,
                                                                                             eblg,
                                                                                             eblgRefined,
                                                                                             eblgFine,
                                                                                             refToCoar,
                                                                                             refToFine,
                                                                                             redistributeOutside));
    }
  }
}

void
PhaseRealm::defineParticleMesh()
{
  CH_TIME("PhaseRealm::defineParticleMesh");
  if (m_verbose) {
    pout() << "PhaseRealm::defineParticleMesh" << endl;
  }

  m_particleMesh.define(m_eblg, m_refinementRatios, m_dx, m_probLo, m_numGhostCells, m_finestLevel);
  m_surfaceDeposition.define(m_eblg, m_eblgCoFi, m_eblgFiCo, m_refinementRatios, m_dx, m_probLo, m_finestLevel, 1);
}

void
PhaseRealm::defineMirrorSurfaceData(const int a_lmin)
{
  CH_TIME("PhaseRealm::defineMirrorSurfaceData");

  if (m_verbose) {
    pout() << "PhaseRealm::defineMirrorSurfaceData" << endl;
  }

  const bool doThisOperator = this->queryOperator(s_mirror_deposition);

  m_mirrorSurfaceData.resize(1 + m_finestLevel);
  m_mirrorHasCutCells.resize(1 + m_finestLevel, 0);

  if (!doThisOperator) {
    return;
  }

  // Two independent requirements, and they must stay independent: eb_ghost sizes the region over which EBISBox
  // geometric data is valid, while num_ghost sizes the DATA and is therefore what the exchanges below can deliver.
  // Deriving one from the other has been a recurring mistake; do not.
  if (m_numEbGhostsCells < 2) {
    MayDay::Error("PhaseRealm::defineMirrorSurfaceData - mirrored deposition needs 'AmrMesh.eb_ghost' >= 2 so that "
                  "EBISBox is valid over the 5^D curvature-fit stencil");
  }

  if (m_numGhostCells < m_mirrorBandRadius) {
    MayDay::Error("PhaseRealm::defineMirrorSurfaceData - mirrored deposition needs 'AmrMesh.num_ghost' >= "
                  "'AmrMesh.mirror_band_radius' so that the exchange delivers the cut-cell data the band reads");
  }

  constexpr int fitStencilRadius = 2;

  // Take the cut-cell centroid and normal from the implicit function rather than from EBISBox, which is a large
  // accuracy win because the fit is S ~ dn/dx and so amplifies any position error by 1/dx. Measured mean
  // |d(2H)|/|2H| against the analytic surface:
  //
  //                            2-D sphere R=4dx    3-D sphere    3-D torus
  //   EBISBox normal+centroid       6.97E-3          1.12E-3       3.549E-2
  //   implicit function             9.97E-7          9.53E-8       3.564E-2
  //
  // The torus is moved by neither, and its 3.5% is a THIRD mechanism -- most likely that a torus is not umbilic, so
  // the shape operator genuinely varies across the 5^D stencil and a fit assuming it constant carries model error.
  // Do not read the two sphere columns as evidence about that case.
  //
  // OFF by default, and the trade is deliberate. Taking the geometry from the implicit function costs roughly
  // 2*SpaceDim + 1 evaluations of BaseIF::value per cut cell, doubled when the centroid moves and the gradient is
  // re-evaluated at the projected point -- about thirteen in 3-D, at every regrid. That is nothing for an analytic
  // sphere and potentially a great deal for a tessellated geometry, where each call is a BVH traversal. EBISBox
  // already holds the same information discretely and hands it over for free.
  //
  // What it costs in accuracy, measured against the analytic surface as mean |d(2H)|/|2H|:
  //
  //                                2-D sphere R=4dx    3-D sphere    3-D torus
  //   EBISBox (this default)            6.97E-3          1.12E-3       3.549E-2
  //   implicit function                 9.97E-7          9.53E-8       3.564E-2
  //
  // The whole difference is the CENTROID, not the normal: EBISBox::normal already comes from the implicit function
  // and agrees with it to ~1E-5 degrees in 2-D, whereas the boundary centroid comes from Chombo's moment
  // reconstruction and sits 0.015-0.027 dx off the true surface. The fit is S ~ dn/dx, so that position error is
  // amplified by 1/dx. In 3-D the normals carry error too, up to 1.6 degrees on a torus, and both matter.
  //
  // Turn it on for a geometry whose implicit function is cheap and whose curvature has to be right. The torus is
  // moved by neither, so do not reach for this to fix that case.
  //
  // One switch rather than one per quantity: two independent toggles is four states to reason about, and the
  // diagnostic below would then have to be read very carefully to know which was in play.
  constexpr bool useImplicitFunctionGeometry = false;

  // Finite-difference step for the gradient, in cells. Two-sided failure: too small cancels, too large truncates.
  // 1E-2 sits well inside both for double precision on an analytic implicit function. Note this truncation is why
  // the implicit-function normal is slightly WORSE than EBISBox's in 2-D, where Chombo's normals are already exact
  // to ~1E-5 degrees; in 3-D they carry up to 1.6 degrees of error and the implicit function wins decisively.
  const Real gradientStep = 1.E-2;

  // Curvature correction on or off.  With it OFF, pass B is skipped entirely: the shape operator stays zero, every
  // cut cell keeps the Status::Planar that pass A gave it, and the whole scheme degrades to a PLANE mirror with no
  // special-casing anywhere downstream -- MirrorDeposition::reflect already does the right thing for S = 0, since
  // Sw vanishes, nhat collapses to n_c, and the Jacobian becomes (1 - 0 + 0)/(1 + 0 + 0) = 1.
  //
  // The point of the switch is to find out whether the correction earns its cost in a production run. What it buys,
  // measured against analytic surfaces: the fitted curvature is accurate to a fraction of a percent on a resolved
  // sphere, and the exact Jacobian matters most where the surface is tightly curved relative to dx -- a concave
  // cavity at R = 3*dx runs at J = 28.9 under CIC and J = 483.8 under TSC, so on that geometry a plane mirror
  // under-weights the image by two to three orders of magnitude. On a gently curved surface it is a small effect.
  //
  // What it costs is the least-squares fit over the 5^D neighbourhood of every cut cell, once per regrid: a
  // symmetric solve plus a residual pass per cell, with the conditioning and crease guards on top.
  constexpr bool useCurvatureCorrection = true;

  // Largest projection, in cells, that will be believed. Not a tuning knob for accuracy -- a tripwire. Composite
  // geometries are min/max of several implicit functions and are therefore NOT differentiable at a seam; a one-sided
  // gradient there can throw the Newton step onto the wrong facet, moving the centroid by order a cell rather than a
  // fortieth of one. Refuse those and keep Chombo's centroid, which is at least on the right piece of surface. The
  // measured shifts are 0.015-0.027 dx mean and 0.045 dx worst, so this carries about a factor of ten of headroom
  // and should fire only on a seam, never on ordinary discretization.
  const Real centroidShiftLimit = 0.5;

  // Counters, reduced at the end. fitOk/fallback are about the CURVATURE; bandNoCutCell is about the GRIDS, which is
  // why they are counted apart -- a run that trips the first has under-resolved geometry, a run that trips the second
  // has a level whose grids do not cover the embedded-boundary band.
  long long fitOk               = 0;
  long long fitTooFewNeighbours = 0;
  long long fitIllConditioned   = 0;
  long long fitCrease           = 0;
  long long bandNoCutCell       = 0;

  // How far the implicit function's geometry is from EBISBox's. Reported so a silent change of source, or a
  // geometry with no implicit function at all, cannot go unnoticed.
  Real      sumNormalAngle     = 0.0;
  Real      maxNormalAngle     = 0.0;
  long long numNormalCmp       = 0;
  long long numNormalFlipped   = 0;
  Real      sumCentroidShift   = 0.0;
  Real      maxCentroidShift   = 0.0;
  long long numCentroidRefused = 0;

  for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl   = m_grids[lvl];
    const EBISLayout&        ebisl = m_ebisl[lvl];
    const DataIterator&      dit   = dbl.dataIterator();
    const Real               dx    = m_dx[lvl];

    m_mirrorSurfaceData[lvl] = RefCountedPtr<LevelData<EBCellFAB>>(
      new LevelData<EBCellFAB>(dbl, MirrorDeposition::numComp, m_numGhostCells * IntVect::Unit, EBCellFactory(ebisl)));

    // Zero first, so that a cell nothing writes carries Status::None by default rather than by accident.
    DataOps::setValue(*m_mirrorSurfaceData[lvl], 0.0);

    const int nbox = dit.size();

    // ---------------------------------------------------------------------------------------------------------
    // Pass A -- seed the cut cells of the VALID box with their centroid and normal.
    // ---------------------------------------------------------------------------------------------------------
    long long numCutCells = 0;

#pragma omp parallel for schedule(runtime)                                                                         \
  reduction(+ : numCutCells, sumNormalAngle, numNormalCmp, numNormalFlipped, sumCentroidShift, numCentroidRefused) \
  reduction(max : maxNormalAngle, maxCentroidShift)

    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      EBCellFAB&     surf    = (*m_mirrorSurfaceData[lvl])[din];
      BaseFab<Real>& surfReg = surf.getSingleValuedFAB();
      const EBISBox& ebisbox = surf.getEBISBox();

      // Use the realm's own cut-cell iterator rather than re-deriving the irregular set: it is built in regridBase
      // from the same EBISBox and is the iteration every other operator here uses, so the two cannot drift apart.
      VoFIterator& vofit = (*m_vofIter[lvl])[din];

      for (vofit.reset(); vofit.ok(); ++vofit) {
        const VolIndex& vof = vofit();
        const IntVect&  iv  = vof.gridIndex();

        // Multiply-cut cells are refined away by this project's workflows, but nothing in the library enforces that,
        // so skip them rather than pick one of the two VoFs arbitrarily -- a position alone cannot select a VoF.
        if (ebisbox.numVoFs(iv) != 1) {
          continue;
        }

        // Location::position multiplies by dx (see its 'ret *= a_dx'), so this is a physical position. Do NOT
        // difference raw bndryCentroid values from different cells: that drops the inter-cell offset and inflates the
        // fitted curvature by about an order of magnitude.
        RealVect       xc   = m_probLo + Location::position(Location::Cell::Boundary, vof, ebisbox, dx);
        const RealVect ebNc = ebisbox.normal(vof);
        RealVect       nc   = ebNc;

        // Implicit-function geometry, see useImplicitFunctionGeometry above. Gated on the switch as well as on the
        // implicit function existing, because the comparison below evaluates BaseIF::value whether or not the result
        // is used -- leaving it running when the switch is off would pay the whole cost for a diagnostic.
        if (useImplicitFunctionGeometry && !m_baseif.isNull()) {
          const Real h = gradientStep * dx;

          const auto gradientAt = [&](const RealVect& a_x) -> RealVect {
            RealVect g = RealVect::Zero;

            for (int dir = 0; dir < SpaceDim; dir++) {
              RealVect lo = a_x;
              RealVect hi = a_x;

              lo[dir] -= h;
              hi[dir] += h;

              g[dir] = (m_baseif->value(hi) - m_baseif->value(lo)) / (2.0 * h);
            }

            return g;
          };

          RealVect grad    = gradientAt(xc);
          Real     gradLen = grad.vectorLength();

          // One Newton step onto phi = 0. EBISBox::normal comes from the implicit function and was measured to be
          // exact, but the boundary CENTROID comes from Chombo's moment reconstruction, which is genuinely discrete.
          // Since the fit is S ~ dn/dx, a centroid error eps enters the curvature as eps/dx -- so this is the other
          // of the fit's only two inputs, and the one still untested.
          if (gradLen > 0.0) {
            const RealVect xProj = xc - (m_baseif->value(xc) / (gradLen * gradLen)) * grad;
            const Real     shift = (xProj - xc).vectorLength() / dx;

            sumCentroidShift += shift;
            maxCentroidShift = std::max(maxCentroidShift, shift);

            if (shift > centroidShiftLimit) {
              numCentroidRefused++;
            }
            else if (useImplicitFunctionGeometry) {
              xc      = xProj;
              grad    = gradientAt(xc);
              gradLen = grad.vectorLength();
            }
          }

          if (gradLen > 0.0) {
            RealVect lsfNc = grad / gradLen;

            // Sign the level-set normal to the EB normal's half-space rather than assuming a convention, then report
            // the angle BEFORE that alignment so a genuinely inverted normal shows up as ~180 degrees rather than
            // being silently folded to ~0.
            const Real cosRaw = std::max(-1.0, std::min(1.0, lsfNc.dotProduct(ebNc)));

            // Report the angle AFTER folding out the sign convention, so a perfect match reads 0 rather than 180 and
            // the residual disagreement is legible. The flip count below says whether a convention difference exists
            // at all, so nothing is hidden by folding.
            const Real angle = std::acos(std::min(1.0, std::abs(cosRaw))) * 180.0 / M_PI;

            sumNormalAngle += angle;
            maxNormalAngle = std::max(maxNormalAngle, angle);
            numNormalCmp++;

            if (cosRaw < 0.0) {
              lsfNc = -lsfNc;
              numNormalFlipped++;
            }

            if (useImplicitFunctionGeometry) {
              nc = lsfNc;
            }
          }
        }

        // Provisional: a valid plane, upgraded to a fitted patch by pass B if the fit succeeds.
        surfReg(iv, MirrorDeposition::compStatus) = static_cast<Real>(MirrorDeposition::Status::Planar);

        for (int dir = 0; dir < SpaceDim; dir++) {
          surfReg(iv, MirrorDeposition::compCentroid + dir) = xc[dir];
          surfReg(iv, MirrorDeposition::compNormal + dir)   = nc[dir];
        }

        numCutCells++;
      }
    }

    m_mirrorSurfaceData[lvl]->exchange();

    // ---------------------------------------------------------------------------------------------------------
    // Pass B -- fit the shape operator over the 5^D neighbourhood of each cut cell. Skipped entirely when the
    // curvature correction is off, which leaves S = 0 and Status::Planar and gives a plane mirror.
    // ---------------------------------------------------------------------------------------------------------
    if (useCurvatureCorrection) {

#pragma omp parallel for schedule(runtime) reduction(+ : fitOk, fitTooFewNeighbours, fitIllConditioned, fitCrease)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din = dit[mybox];

        EBCellFAB&     surf    = (*m_mirrorSurfaceData[lvl])[din];
        BaseFab<Real>& surfReg = surf.getSingleValuedFAB();
        const EBISBox& ebisbox = surf.getEBISBox();

        const Box validBox = dbl[din];

        // Clip to the EBISBox region as well as to the data FAB. A LevelData ghost box is not clipped to the problem
        // domain, but EBISBox only holds data over its own region, so a box on the domain edge would otherwise ask
        // EBISBox about cells it has nothing for.
        const Box readBox = grow(validBox, fitStencilRadius) & surfReg.box() & ebisbox.getRegion();

        for (BoxIterator bit(validBox); bit.ok(); ++bit) {
          const IntVect& iv = bit();

          if (surfReg(iv, MirrorDeposition::compStatus) < 0.5) {
            continue;
          }

          RealVect xi;
          RealVect ni;

          for (int dir = 0; dir < SpaceDim; dir++) {
            xi[dir] = surfReg(iv, MirrorDeposition::compCentroid + dir);
            ni[dir] = surfReg(iv, MirrorDeposition::compNormal + dir);
          }

          // Orthonormal tangent basis. Cross the normal with the axis it is least aligned with, so the cross product is
          // never near-degenerate. This basis is LOCAL TO THE FIT and is discarded at the end of it -- the shape
          // operator is stored in world coordinates precisely so that no later reader has to agree with this choice.
          RealVect tangents[MirrorDeposition::numTangents];
          {
            int minDir = 0;

            for (int dir = 1; dir < SpaceDim; dir++) {
              if (std::abs(ni[dir]) < std::abs(ni[minDir])) {
                minDir = dir;
              }
            }

            RealVect axis = RealVect::Zero;
            axis[minDir]  = 1.0;

#if CH_SPACEDIM == 2
            tangents[0] = RealVect(-ni[1], ni[0]);
#else
            RealVect t1 = RealVect(ni[1] * axis[2] - ni[2] * axis[1],
                                   ni[2] * axis[0] - ni[0] * axis[2],
                                   ni[0] * axis[1] - ni[1] * axis[0]);
            RealVect t2 = RealVect::Zero;

            const Real t1len = t1.vectorLength();

            if (t1len > 0.0) {
              t1 /= t1len;
              t2 = RealVect(ni[1] * t1[2] - ni[2] * t1[1],
                            ni[2] * t1[0] - ni[0] * t1[2],
                            ni[0] * t1[1] - ni[1] * t1[0]);

              const Real t2len = t2.vectorLength();

              if (t2len > 0.0) {
                t2 /= t2len;
              }
            }

            tangents[0] = t1;
            tangents[1] = t2;
#endif
          }

          // Least-squares fit of the shape operator: minimize sum_j |v_j - S u_j|^2 over the tangential displacement
          // u_j and tangential normal-difference v_j of the neighbouring cut cells.
          constexpr int T = MirrorDeposition::numTangents;

          Real A[T][T] = {};
          Real B[T][T] = {};
          Real uNormSq = 0.0;

          int numUsable = 0;

          // Neighbour data is read FROM THE FIELD, never from EBISBox. That is what enforces the delivered-not-merely-
          // existing rule: a cut cell that no box on this level owns was never seeded and never exchanged, so it reads
          // back as Status::None here even though EBISBox can see it. Reading EBISBox instead would fit a patch to
          // geometry this level has no data for.
          Box window(iv - fitStencilRadius * IntVect::Unit, iv + fitStencilRadius * IntVect::Unit);
          window &= readBox;

          for (BoxIterator nbrIt(window); nbrIt.ok(); ++nbrIt) {
            const IntVect& jv = nbrIt();

            if (jv == iv) {
              continue;
            }

            if (surfReg(jv, MirrorDeposition::compStatus) < 0.5) {
              continue;
            }

            if (ebisbox.numVoFs(jv) != 1) {
              continue;
            }

            RealVect xj;
            RealVect nj;

            for (int dir = 0; dir < SpaceDim; dir++) {
              xj[dir] = surfReg(jv, MirrorDeposition::compCentroid + dir);
              nj[dir] = surfReg(jv, MirrorDeposition::compNormal + dir);
            }

            const RealVect dxj = xj - xi;
            const RealVect dnj = nj - ni;

            Real u[T];
            Real v[T];

            for (int a = 0; a < T; a++) {
              u[a] = tangents[a].dotProduct(dxj);
              v[a] = tangents[a].dotProduct(dnj);
            }

            for (int a = 0; a < T; a++) {
              uNormSq += u[a] * u[a];

              for (int b = 0; b < T; b++) {
                A[a][b] += u[a] * u[b];
                B[a][b] += u[a] * v[b];
              }
            }

            numUsable++;
          }

          bool accepted = false;
          Real S[T][T]  = {};

          if (numUsable < 4) {
            fitTooFewNeighbours++;
          }
          else {
            // Condition test on A, normalized by the natural scale (numUsable * dx^2) so the threshold is a pure
            // number in both dimensions. A count of usable neighbours is necessary but NOT sufficient: four cut cells
            // can be nearly collinear in the tangent plane -- a ridge, a fin, a cylinder whose cut cells run along its
            // axis -- and then A is near-singular while the count looks healthy.
            const Real scale = 1.0 / (numUsable * dx * dx);

#if CH_SPACEDIM == 2
            const Real lambdaMin = A[0][0] * scale;
#else
            const Real tr   = (A[0][0] + A[1][1]) * scale;
            const Real det  = (A[0][0] * A[1][1] - A[0][1] * A[1][0]) * scale * scale;
            const Real disc = std::max(0.0, 0.25 * tr * tr - det);

            const Real lambdaMin = 0.5 * tr - std::sqrt(disc);
#endif

            if (!(lambdaMin > m_mirrorConditionTolerance)) {
              fitIllConditioned++;
            }
            else {
              // Solve S^T = A^{-1} B in closed form. A 2x2 symmetric inverse and a 2x2 eigenproblem are about fifteen
              // allocation-free lines; LaPackUtils::computeSVD would give the conditioning number directly but
              // allocates std::vectors per call, inside what is an OpenMP loop over every cut cell.
#if CH_SPACEDIM == 2
              S[0][0] = B[0][0] / A[0][0];
#else
              const Real detA = A[0][0] * A[1][1] - A[0][1] * A[1][0];

              const Real inv00 = A[1][1] / detA;
              const Real inv01 = -A[0][1] / detA;
              const Real inv11 = A[0][0] / detA;

              // S^T = A^{-1} B, then transpose into S.
              const Real st00 = inv00 * B[0][0] + inv01 * B[1][0];
              const Real st01 = inv00 * B[0][1] + inv01 * B[1][1];
              const Real st10 = inv01 * B[0][0] + inv11 * B[1][0];
              const Real st11 = inv01 * B[0][1] + inv11 * B[1][1];

              S[0][0] = st00;
              S[0][1] = st10;
              S[1][0] = st01;
              S[1][1] = st11;
#endif

              // Symmetrize. The shape operator is symmetric for a smooth surface; the fitted one is not exactly, and
              // the antisymmetric part is fit noise.
              for (int a = 0; a < T; a++) {
                for (int b = a + 1; b < T; b++) {
                  const Real avg = 0.5 * (S[a][b] + S[b][a]);

                  S[a][b] = avg;
                  S[b][a] = avg;
                }
              }

              // Crease detector. The residual is made dimensionless by dx, so the threshold does not move with
              // refinement. A sphere union, a facet edge or a dielectric corner is not a quadratic patch, and this is
              // what says so.
              Real residualSq = 0.0;

              for (BoxIterator nbrIt(window); nbrIt.ok(); ++nbrIt) {
                const IntVect& jv = nbrIt();

                if (jv == iv) {
                  continue;
                }

                if (surfReg(jv, MirrorDeposition::compStatus) < 0.5) {
                  continue;
                }

                if (ebisbox.numVoFs(jv) != 1) {
                  continue;
                }

                RealVect xj;
                RealVect nj;

                for (int dir = 0; dir < SpaceDim; dir++) {
                  xj[dir] = surfReg(jv, MirrorDeposition::compCentroid + dir);
                  nj[dir] = surfReg(jv, MirrorDeposition::compNormal + dir);
                }

                const RealVect dxj = xj - xi;
                const RealVect dnj = nj - ni;

                for (int a = 0; a < T; a++) {
                  Real pred = 0.0;

                  for (int b = 0; b < T; b++) {
                    pred += S[a][b] * tangents[b].dotProduct(dxj);
                  }

                  const Real r = tangents[a].dotProduct(dnj) - pred;

                  residualSq += r * r;
                }
              }

              const Real residual = (uNormSq > 0.0) ? dx * std::sqrt(residualSq / uNormSq) : 0.0;

              if (residual > m_mirrorCreaseTolerance) {
                fitCrease++;
              }
              else {
                accepted = true;
              }
            }
          }

          if (accepted) {
            Real shape[MirrorDeposition::numShapeComp];

            MirrorDeposition::liftShapeOperator(tangents, S, shape);

            for (int c = 0; c < MirrorDeposition::numShapeComp; c++) {
              surfReg(iv, MirrorDeposition::compShape + c) = shape[c];
            }

            surfReg(iv, MirrorDeposition::compStatus) = static_cast<Real>(MirrorDeposition::Status::Fitted);

            // S_c annihilates the normal by construction; if it does not, the lift or the tangent basis is wrong, and
            // the failure is otherwise silent because tr(S_c) and the Jacobian stay plausible.
            CH_assert(MirrorDeposition::applyShapeOperator(shape, ni).vectorLength() < 1.E-12);

            fitOk++;
          }
          else {
            // Status stays Planar and the shape components stay zero, which gives J = 1 exactly. A refused fit still
            // reflects -- against a plane -- because dropping it removes the correction in the very cells that need it.
            for (int c = 0; c < MirrorDeposition::numShapeComp; c++) {
              surfReg(iv, MirrorDeposition::compShape + c) = 0.0;
            }
          }
        }
      }

      // Publishes the fitted shape operator to the ghosts that pass C reads. Inside the guard because with no fit
      // there is nothing new to publish -- pass A's exchange already delivered the centroids and normals.
      m_mirrorSurfaceData[lvl]->exchange();
    }

    // ---------------------------------------------------------------------------------------------------------
    // Pass C -- extend the patch outwards to the band.
    //
    // The ORDER of this pass relative to the exchange above is load-bearing. Extend first and exchange afterwards,
    // and a band cell at a patch edge whose nearest cut cell lies one patch over gets assigned a too-far cut cell --
    // and the exchange cannot repair it, because both patches computed from the same incomplete data and agree on
    // the wrong answer. This is the one place where the valid-then-exchange pattern of defineMasks does not carry
    // over.
    // ---------------------------------------------------------------------------------------------------------
#pragma omp parallel for schedule(runtime) reduction(+ : bandNoCutCell)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      EBCellFAB&     surf    = (*m_mirrorSurfaceData[lvl])[din];
      BaseFab<Real>& surfReg = surf.getSingleValuedFAB();
      const EBISBox& ebisbox = surf.getEBISBox();

      const Box validBox = dbl[din];
      const Box readBox  = grow(validBox, m_mirrorBandRadius) & surfReg.box() & ebisbox.getRegion();

      // Pass C reads the seeded/fitted values and writes new ones into the same holder, so it needs its own copy of
      // the band cells' results -- otherwise a band cell written early in the iteration becomes a candidate source
      // for one written later, and the patch propagates outwards instead of being assigned.
      BaseFab<Real> extReg(validBox, MirrorDeposition::numComp);
      extReg.setVal(0.0);

      for (BoxIterator bit(validBox); bit.ok(); ++bit) {
        const IntVect& iv = bit();

        // Cut cells keep what pass B gave them.
        if (surfReg(iv, MirrorDeposition::compStatus) >= 0.5) {
          for (int c = 0; c < MirrorDeposition::numComp; c++) {
            extReg(iv, c) = surfReg(iv, c);
          }

          continue;
        }

        if (ebisbox.isCovered(iv)) {
          continue;
        }

        // Euclidean-nearest over the neighbourhood, NOT a Chebyshev sweep. The error ladder that justifies the
        // quadratic patch over a plane was measured with argmin |x_band - x_cut|^2 over all cut cells; sweeps compute
        // a Chebyshev-nearest assignment with an arbitrary tie-break, and the two differ exactly where more than one
        // piece of surface is within the band -- which is the geometry the patch exists to get right.
        const RealVect cellCentre = m_probLo + (RealVect(iv) + 0.5 * RealVect::Unit) * dx;

        Box window(iv - m_mirrorBandRadius * IntVect::Unit, iv + m_mirrorBandRadius * IntVect::Unit);
        window &= readBox;

        Real    bestDistSq = std::numeric_limits<Real>::max();
        IntVect bestIv     = IntVect::Zero;
        bool    found      = false;

        for (BoxIterator nbrIt(window); nbrIt.ok(); ++nbrIt) {
          const IntVect& jv = nbrIt();

          if (surfReg(jv, MirrorDeposition::compStatus) < 0.5) {
            continue;
          }

          RealVect xj;

          for (int dir = 0; dir < SpaceDim; dir++) {
            xj[dir] = surfReg(jv, MirrorDeposition::compCentroid + dir);
          }

          const RealVect delta  = xj - cellCentre;
          const Real     distSq = delta.dotProduct(delta);

          // Fixed lexicographic tie-break, so the assignment does not depend on iteration order or on the box
          // decomposition.
          bool better = distSq < bestDistSq;

          if (!better && found && distSq <= bestDistSq) {
            for (int dir = SpaceDim - 1; dir >= 0; dir--) {
              if (jv[dir] != bestIv[dir]) {
                better = jv[dir] < bestIv[dir];

                break;
              }
            }
          }

          if (better) {
            bestDistSq = distSq;
            bestIv     = jv;
            found      = true;
          }
        }

        if (found) {
          for (int c = 0; c < MirrorDeposition::numComp; c++) {
            extReg(iv, c) = surfReg(bestIv, c);
          }
        }
        else {
          // Reachable where a level's grids do not cover the embedded-boundary band. A diagnostic, not an abort:
          // particles here simply do not reflect.
          bandNoCutCell++;
        }
      }

      // Commit.
      for (BoxIterator bit(validBox); bit.ok(); ++bit) {
        const IntVect& iv = bit();

        for (int c = 0; c < MirrorDeposition::numComp; c++) {
          surfReg(iv, c) = extReg(iv, c);
        }
      }
    }

    // No exchange after pass C. Surface data is read only at a particle's OWN cell, and particles live in the valid
    // region of the patch that owns them.

    // All-reduced, and that is required rather than tidy -- see the member's documentation.
    m_mirrorHasCutCells[lvl] = (ParallelOps::sum(numCutCells) > 0) ? 1 : 0;
  }

  if (m_verbose) {
    fitOk               = ParallelOps::sum(fitOk);
    fitTooFewNeighbours = ParallelOps::sum(fitTooFewNeighbours);
    fitIllConditioned   = ParallelOps::sum(fitIllConditioned);
    fitCrease           = ParallelOps::sum(fitCrease);
    bandNoCutCell       = ParallelOps::sum(bandNoCutCell);

    // The reductions run unconditionally: they are collectives, and "this rank saw a cut
    // cell" is NOT a collective condition -- gating them on it deadlocks the moment one rank's patch has no cut
    // cells. Only the printing is conditional, on the reduced count.
    const Real      sumAll      = ParallelOps::sum(sumNormalAngle);
    const Real      maxAll      = ParallelOps::max(maxNormalAngle);
    const long long cmpAll      = ParallelOps::sum(numNormalCmp);
    const long long flipAll     = ParallelOps::sum(numNormalFlipped);
    const Real      sumShiftAll = ParallelOps::sum(sumCentroidShift);
    const Real      maxShiftAll = ParallelOps::max(maxCentroidShift);
    const long long refusedAll  = ParallelOps::sum(numCentroidRefused);

    // Say which source built the surface data, and distinguish the two ways it can be EBISBox. Chosen and
    // unavailable are very different situations -- one is a deliberate trade, the other is a geometry that cannot
    // offer the alternative -- and both are otherwise silent.
    if (cmpAll == 0) {
      if (!useImplicitFunctionGeometry) {
        pout() << "PhaseRealm::defineMirrorSurfaceData - using EBISBox centroid and normal by configuration; the "
                  "implicit function is the more accurate source but costs evaluations per cut cell (see "
                  "useImplicitFunctionGeometry)"
               << endl;
      }
      else {
        pout() << "PhaseRealm::defineMirrorSurfaceData - no implicit function on this realm; using EBISBox centroid "
                  "and normal, which are the less accurate source (see defineMirrorSurfaceData)"
               << endl;
      }
    }
    else {
      pout() << "PhaseRealm::defineMirrorSurfaceData - LSF vs EBISBox normal over " << cmpAll << " cut cells: mean "
             << sumAll / cmpAll << " deg, max " << maxAll << " deg, " << flipAll << " inverted; using the "
             << (useImplicitFunctionGeometry ? "implicit-function" : "EBISBox") << " normal" << endl;

      pout() << "PhaseRealm::defineMirrorSurfaceData - centroid off the true surface by mean " << sumShiftAll / cmpAll
             << " dx, max " << maxShiftAll << " dx; using the "
             << (useImplicitFunctionGeometry ? "projected" : "EBISBox") << " centroid, " << refusedAll
             << " refused as beyond " << centroidShiftLimit << " dx" << endl;
    }

    if (!useCurvatureCorrection) {
      pout() << "PhaseRealm::defineMirrorSurfaceData - curvature correction DISABLED; the shape operator is zero and "
                "the scheme is a plane mirror (see useCurvatureCorrection)"
             << endl;
    }

    pout() << "PhaseRealm::defineMirrorSurfaceData - curvature fit: " << fitOk << " ok, " << fitTooFewNeighbours
           << " too few neighbours, " << fitIllConditioned << " ill-conditioned, " << fitCrease << " crease; "
           << bandNoCutCell << " band cells found no cut cell" << endl;
  }
}

void
PhaseRealm::defineGradSten(const int /*a_lmin*/)
{
  CH_TIME("PhaseRealm::defineGradSten");
  if (m_verbose) {
    pout() << "PhaseRealm::defineGradSten" << endl;
  }

  const bool doThisOperator = this->queryOperator(s_eb_gradient);

  m_gradientOp.resize(1 + m_finestLevel);

  if (doThisOperator) {

    // Define gradient operator.
    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {

      const bool hasFine = lvl < m_finestLevel;

      EBLevelGrid eblg;
      EBLevelGrid eblgFine;
      EBLevelGrid eblgFiCo;

      int refRat;

      eblg = *m_eblg[lvl];
      if (hasFine) {
        refRat   = m_refinementRatios[lvl];
        eblgFine = *m_eblg[lvl + 1];
        eblgFiCo = *m_eblgFiCo[lvl + 1];
      }
      else {
        refRat   = 1;
        eblgFine = EBLevelGrid();
        eblgFiCo = EBLevelGrid();
      }

      const int order  = 2;
      const int weight = 1;

      m_gradientOp[lvl] = RefCountedPtr<EBGradient>(new EBGradient(eblg,
                                                                   eblgFine,
                                                                   eblgFiCo,
                                                                   hasFine,
                                                                   m_dx[lvl],
                                                                   refRat,
                                                                   order,
                                                                   weight,
                                                                   m_numGhostCells * IntVect::Unit));
    }
  }
}

void
PhaseRealm::defineIrregSten()
{
  CH_TIME("PhaseRealm::defineIrregSten");
  if (m_verbose) {
    pout() << "PhaseRealm::defineIrregSten" << endl;
  }

  const bool doThisOperator = this->queryOperator(s_eb_irreg_interp);

  m_cellCentroidInterpolation.resize(1 + m_finestLevel);
  m_ebCentroidInterpolation.resize(1 + m_finestLevel);

  if (doThisOperator) {
    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      m_cellCentroidInterpolation[lvl] = RefCountedPtr<CellCentroidInterpolation>(
        new CellCentroidInterpolation(*m_eblg[lvl], m_dx[lvl], m_cellCentroidInterpolationType));

      m_ebCentroidInterpolation[lvl] = RefCountedPtr<EBCentroidInterpolation>(
        new EBCentroidInterpolation(*m_eblg[lvl], m_dx[lvl], m_ebCentroidInterpolationType));
    }
  }
}

void
PhaseRealm::defineNonConservativeDivergence(const int /*a_lmin*/)
{
  CH_TIME("PhaseRealm::defineNonConservativeDivergence");
  if (m_verbose) {
    pout() << "PhaseRealm::defineNonConservativeDivergence" << endl;
  }

  const bool doThisOperator = this->queryOperator(s_noncons_div);

  m_nonConservativeDivergence.resize(1 + m_finestLevel);

  if (doThisOperator) {
    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      m_nonConservativeDivergence[lvl] = RefCountedPtr<EBNonConservativeDivergence>(
        new EBNonConservativeDivergence(*m_eblg[lvl], m_redistributionRadius));
    }
  }
}

const RefCountedPtr<EBIndexSpace>&
PhaseRealm::getEBIndexSpace() const
{
  return m_ebis;
}

const Vector<int>&
PhaseRealm::getRefinementRatios() const
{
  return m_refinementRatios;
}

const Vector<Real>&
PhaseRealm::getDx() const
{
  return m_dx;
}

const Vector<DisjointBoxLayout>&
PhaseRealm::getGrids() const
{
  return m_grids;
}

const Vector<ProblemDomain>&
PhaseRealm::getDomains() const
{
  return m_domains;
}

const Vector<EBISLayout>&
PhaseRealm::getEBISLayout() const
{
  return m_ebisl;
}

const Vector<RefCountedPtr<EBLevelGrid>>&
PhaseRealm::getEBLevelGrid() const
{
  return m_eblg;
}

const Vector<RefCountedPtr<EBLevelGrid>>&
PhaseRealm::getEBLevelGridCoFi() const
{
  return m_eblgCoFi;
}

const Vector<RefCountedPtr<EBLevelGrid>>&
PhaseRealm::getEBLevelGridFiCo() const
{
  return m_eblgFiCo;
}

Vector<RefCountedPtr<LayoutData<VoFIterator>>>&
PhaseRealm::getVofIterator() const
{
  return m_vofIter;
}

Vector<RefCountedPtr<LayoutData<VoFIterator>>>&
PhaseRealm::getMultiCutVofIterator() const
{
  return m_multiCutVofIter;
}

Vector<RefCountedPtr<LayoutData<std::array<FaceIterator, SpaceDim>>>>&
PhaseRealm::getFaceIterator() const
{
  return m_faceIter;
}

Vector<RefCountedPtr<LayoutData<std::array<FaceIterator, SpaceDim>>>>&
PhaseRealm::getFaceIteratorNoBoundary() const
{
  return m_faceIterNoBoundary;
}

Vector<RefCountedPtr<LayoutData<std::array<FaceIterator, SpaceDim>>>>&
PhaseRealm::getFaceIteratorWithTangentialGhosts() const
{
  return m_faceIterTanGhost;
}

Vector<RefCountedPtr<LayoutData<std::array<FaceIterator, SpaceDim>>>>&
PhaseRealm::getMultiCutFaceIterator() const
{
  return m_multiCutFaceIter;
}

const Vector<RefCountedPtr<EBGradient>>&
PhaseRealm::getGradientOp() const
{
  if (!this->queryOperator(s_eb_gradient)) {
    MayDay::Error("PhaseRealm::getGradientOp - operator not registered!");
  }

  return m_gradientOp;
}

const Vector<RefCountedPtr<CellCentroidInterpolation>>&
PhaseRealm::getCellCentroidInterpolation() const
{
  if (!this->queryOperator(s_eb_irreg_interp)) {
    MayDay::Error("PhaseRealm::getCellCentroidInterpolation - operator not registered!");
  }

  return m_cellCentroidInterpolation;
}

const Vector<RefCountedPtr<EBCentroidInterpolation>>&
PhaseRealm::getEBCentroidInterpolation() const
{
  if (!this->queryOperator(s_eb_irreg_interp)) {
    MayDay::Error("PhaseRealm::getEBCentroidInterpolation - operator not registered!");
  }

  return m_ebCentroidInterpolation;
}

const Vector<RefCountedPtr<EBNonConservativeDivergence>>&
PhaseRealm::getNonConservativeDivergence() const
{
  if (!this->queryOperator(s_noncons_div)) {
    MayDay::Error("PhaseRealm::getNonConservativeDivergence - operator not registered!");
  }

  return m_nonConservativeDivergence;
}

Vector<RefCountedPtr<EBCoarAve>>&
PhaseRealm::getCoarseAverage() const
{
  if (!this->queryOperator(s_eb_coar_ave)) {
    MayDay::Error("PhaseRealm::getCoarseAverage - operator not registered!");
  }

  return m_coarAve;
}

Vector<RefCountedPtr<EBMultigridInterpolator>>&
PhaseRealm::getMultigridInterpolator() const
{
  if (!this->queryOperator(s_eb_multigrid)) {
    MayDay::Error("PhaseRealm::getEBMultigridInterpolator - operator not registered!");
  }

  return m_multigridInterpolator;
}

Vector<RefCountedPtr<EBGhostCellInterpolator>>&
PhaseRealm::getGhostCellInterpolator() const
{
  if (!this->queryOperator(s_eb_fill_patch)) {
    MayDay::Error("PhaseRealm::getFillPatch - operator not registered!");
  }

  return m_ghostCellInterpolator;
}

Vector<RefCountedPtr<EBCoarseToFineInterp>>&
PhaseRealm::getFineInterp() const
{
  if (!this->queryOperator(s_eb_fine_interp)) {
    MayDay::Error("PhaseRealm::getFineInterp - operator not registered!");
  }

  return m_ebFineInterp;
}

Vector<RefCountedPtr<EBReflux>>&
PhaseRealm::getFluxRegister() const
{
  if (!this->queryOperator(s_eb_flux_reg)) {
    MayDay::Error("PhaseRealm::getFluxRegister - operator not registered!");
  }

  return m_ebReflux;
}

Vector<RefCountedPtr<EBFluxRedistribution>>&
PhaseRealm::getRedistributionOp() const
{
  if (!this->queryOperator(s_eb_redist)) {
    MayDay::Error("PhaseRealm::getRedistributionOp - operator not registered!");
  }

  return m_redistributionOp;
}

EBAMRParticleMesh&
PhaseRealm::getParticleMesh() const
{
  if (!this->queryOperator(s_particle_mesh)) {
    MayDay::Error("PhaseRealm::getParticleMesh - operator not registered!");
  }

  return m_particleMesh;
}

EBAMRSurfaceDeposition&
PhaseRealm::getSurfaceDeposition() const
{
  if (!this->queryOperator(s_particle_mesh)) {
    MayDay::Error("PhaseRealm::getSurfaceDepostion - operator not registered!");
  }

  return m_surfaceDeposition;
}

const EBAMRFAB&
PhaseRealm::getLevelset() const
{
  if (!this->queryOperator(s_levelset)) {
    MayDay::Error("PhaseRealm::getLevelset - operator not registered!");
  }

  return m_levelset;
}

const EBAMRCellData&
PhaseRealm::getRegularCells() const
{
  return m_regularCells;
}

const EBAMRCellData&
PhaseRealm::getCoveredCells() const
{
  return m_coveredCells;
}

const EBAMRCellData&
PhaseRealm::getNotCoveredCells() const
{
  return m_notCoveredCells;
}

const EBAMRCellData&
PhaseRealm::getIrregularCells() const
{
  return m_irregularCells;
}

const EBAMRCellData&
PhaseRealm::getMirrorSurfaceData() const
{
  return m_mirrorSurfaceData;
}

const Vector<int>&
PhaseRealm::getMirrorHasCutCells() const
{
  return m_mirrorHasCutCells;
}

#include <CD_NamespaceFooter.H>
