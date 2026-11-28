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
#include <EBIndexSpace.H>
#include <EBGraphFactory.H>
#include <EBDataFactory.H>

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
#include <CD_VofUtils.H>
#include <CD_NamespaceHeader.H>

constexpr int PhaseRealm::s_redistributionRadius;

PhaseRealm::PhaseRealm() : m_isDefined(false), m_profile(false), m_verbose(false), m_mirrorCurvatureCorrection(true)
{
  CH_TIME("PhaseRealm::PhaseRealm");

  // Default settings

  this->registerOperator(s_eb_gradient);
  this->registerOperator(s_eb_irreg_interp);

  // Adding this for debugging purposes. These are read with query rather than get -- they are optional debug knobs
  // with working defaults, not physics parameters, so an input that says nothing about them is not an error.
  ParmParse pp("PhaseRealm");
  pp.query("profile", m_profile);
  pp.query("verbosity", m_verbose);
  pp.query("mirror_curvature_correction", m_mirrorCurvatureCorrection);
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
                   const int                             a_srdOrder,
                   const Real                            a_srdTargetVolFrac,
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
  m_srdOrder                      = a_srdOrder;
  m_srdTargetVolFrac              = a_srdTargetVolFrac;
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
  m_eblgFineHalf.resize(0);
  m_fiCoFineCells.resize(0);
  m_ebisFiCo = RefCountedPtr<EBIndexSpace>(nullptr);
  m_mgLevelGrids.resize(0);
  m_mgLevelGridsCoFi.resize(0);
  m_mgCoarAve.resize(0);
  m_finerVoFs.resize(0);
  m_mgFinerVoFs.resize(0);
  m_mgRestrict.resize(0);
  m_mgProlong.resize(0);
  m_mgLevelRestrict.resize(0);
  m_mgLevelProlong.resize(0);
  m_vofIter.resize(0);
  m_multiCutVofIter.resize(0);
  m_faceIter.resize(0);
  m_faceIterNoBoundary.resize(0);
  m_faceIterTanGhost.resize(0);
  m_multiCutFaceIter.resize(0);
  m_domainFaceIterLo.resize(0);
  m_domainFaceIterHi.resize(0);
  m_domainFaceIterTanGhostLo.resize(0);
  m_domainFaceIterTanGhostHi.resize(0);
  m_coarAve.resize(0);
  m_multigridInterpolator.resize(0);
  m_ebFineInterp.resize(0);
  m_ebReflux.resize(0);
  m_massRedistributionOp.resize(0);
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

    // Built with the coarsened grids, so that load-balancing routines, which run after regridBase but before
    // regridOperators, can define coarse-fine operators.
    if (m_profile) {
      pout() << "before/after finer VoF lists" << endl;
      MemoryReport::getMaxMinMemoryUsage();
    }
    timer.startEvent("Finer VoF lists");
    this->defineFinerVoFs(a_lmin);
    timer.stopEvent("Finer VoF lists");
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
      pout() << "before/after multigrid finer VoF lists" << endl;
      MemoryReport::getMaxMinMemoryUsage();
    }
    timer.startEvent("Multigrid finer VoF lists");
    this->defineMultigridFinerVoFs();
    timer.stopEvent("Multigrid finer VoF lists");
    if (m_profile) {
      MemoryReport::getMaxMinMemoryUsage();
      pout() << endl;
    }

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
      pout() << "before/after multigrid coarave define" << endl;
      MemoryReport::getMaxMinMemoryUsage();
    }
    timer.startEvent("Multigrid EBCoarAve");
    this->defineMultigridCoarseAverage();
    timer.stopEvent("Multigrid EBCoarAve");
    timer.startEvent("Multigrid restrict/prolong");
    this->defineMultigridLevelTransfer();
    timer.stopEvent("Multigrid restrict/prolong");
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
    this->defineMassRedistOper(a_lmin, 1);
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
        a_operator == s_eb_flux_reg || a_operator == s_eb_mass_redist || a_operator == s_eb_state_redist ||
        a_operator == s_noncons_div || a_operator == s_eb_gradient || a_operator == s_particle_mesh ||
        a_operator == s_eb_irreg_interp || a_operator == s_eb_multigrid || a_operator == s_levelset ||
        a_operator == s_mirror_deposition)) {

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

    if (hasCoar) {
      m_eblg[lvl]->setMaxCoarseningRatio(m_refinementRatios[lvl - 1], &(*m_ebis));
    }

    m_ebisl[lvl] = m_eblg[lvl]->getEBISL();

    // Define the coarsened grids.
    if (lvl > 0) {
      m_eblgCoFi[lvl - 1] = RefCountedPtr<EBLevelGrid>(new EBLevelGrid());

      coarsen(*m_eblgCoFi[lvl - 1], *m_eblg[lvl], m_refinementRatios[lvl - 1]);
    }
  }

  this->defineEBLevelGridFiCo();
}

void
PhaseRealm::defineEBLevelGridFiCo()
{
  CH_TIME("PhaseRealm::defineEBLevelGridFiCo");
  if (m_verbose) {
    pout() << "PhaseRealm::defineEBLevelGridFiCo" << endl;
  }

  m_eblgFiCo.resize(1 + m_finestLevel);
  m_eblgFineHalf.resize(1 + m_finestLevel);
  m_fiCoFineCells.resize(1 + m_finestLevel);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    m_eblgFiCo[lvl]      = RefCountedPtr<EBLevelGrid>(nullptr);
    m_eblgFineHalf[lvl]  = RefCountedPtr<EBLevelGrid>(nullptr);
    m_fiCoFineCells[lvl] = RefCountedPtr<LevelData<BaseFab<bool>>>(nullptr);
  }

  m_ebisFiCo = RefCountedPtr<EBIndexSpace>(nullptr);

  if (m_finestLevel == 0) {
    return;
  }

  // Sources of the half-resolution layouts that coarsening the copied geometry by 4 steps through.
  for (int lvl = 1; lvl <= m_finestLevel; lvl++) {
    const int refRat = m_refinementRatios[lvl - 1];

    if (refRat != 2 && refRat != 4) {
      MayDay::Abort("PhaseRealm::defineEBLevelGridFiCo - refinement ratios must be 2 or 4");
    }

    if (refRat == 4) {
      m_eblgFineHalf[lvl] = RefCountedPtr<EBLevelGrid>(new EBLevelGrid());

      coarsen(*m_eblgFineHalf[lvl], *m_eblg[lvl], 2);
    }
  }

  // The copying index space spans the realm's domains, finest first, each a coarsening by two of the one before.
  Vector<ProblemDomain> domains;

  ProblemDomain domain = m_domains[m_finestLevel];
  domains.push_back(domain);

  while (domain.domainBox() != m_domains[0].domainBox()) {
    domain.coarsen(2);
    domains.push_back(domain);
  }

  m_ebisFiCo = RefCountedPtr<EBIndexSpace>(new EBIndexSpace());

  m_ebisFiCo->define(domains,
                     m_probLo,
                     m_dx[m_finestLevel],
                     [this](LevelData<EBGraph>&      a_graph,
                            LevelData<EBData>&       a_data,
                            const DisjointBoxLayout& a_grids,
                            const ProblemDomain&     a_domain,
                            const int                a_graphGhost,
                            const int                a_dataGhost) -> void {
                       this->fillCopiedLayout(a_graph, a_data, a_grids, a_domain, a_graphGhost, a_dataGhost);
                     });

  m_fiCoFineCells.resize(1 + m_finestLevel);
  m_fiCoFineCells[0] = RefCountedPtr<LevelData<BaseFab<bool>>>(nullptr);

  for (int lvl = 0; lvl < m_finestLevel; lvl++) {
    const int refRat = m_refinementRatios[lvl];

    // The gradient stencils refine coarse cells up to two cells outside the coarse box, and the mass redistribution
    // up to the redistribution radius. A copy costs no geometry work, so the ghost cells reach both.
    constexpr int gradientStencilWidth = 2;

    const int ghostFiCo = std::max(std::max(m_numEbGhostsCells, m_numGhostCells),
                                   refRat * std::max(gradientStencilWidth, s_redistributionRadius));

    DisjointBoxLayout dblFiCo;
    refine(dblFiCo, m_grids[lvl], refRat);

    m_eblgFiCo[lvl + 1] = RefCountedPtr<EBLevelGrid>(
      new EBLevelGrid(dblFiCo, m_domains[lvl + 1], ghostFiCo, &(*m_ebisFiCo)));

    m_eblgFiCo[lvl + 1]->getEBISL().setMaxCoarseningRatio(refRat, &(*m_ebisFiCo));

#ifndef NDEBUG
    // Mark the cells of level lvl+1, which only the consumers' assertions read.
    m_fiCoFineCells[lvl + 1] = RefCountedPtr<LevelData<BaseFab<bool>>>(
      new LevelData<BaseFab<bool>>(dblFiCo, 1, ghostFiCo * IntVect::Unit));

    LevelData<BaseFab<bool>> fineCells(m_grids[lvl + 1], 1, IntVect::Zero);

    const DataIterator& ditFine  = m_grids[lvl + 1].dataIterator();
    const DataIterator& ditFiCo  = dblFiCo.dataIterator();
    const int           nboxFine = ditFine.size();
    const int           nboxFiCo = ditFiCo.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nboxFine; mybox++) {
      fineCells[ditFine[mybox]].setVal(true);
    }

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nboxFiCo; mybox++) {
      (*m_fiCoFineCells[lvl + 1])[ditFiCo[mybox]].setVal(false);
    }

    fineCells.copyTo(Interval(0, 0), *m_fiCoFineCells[lvl + 1], Interval(0, 0));
#endif
  }
}

void
PhaseRealm::fillCopiedLayout(LevelData<EBGraph>&      a_graph,
                             LevelData<EBData>&       a_data,
                             const DisjointBoxLayout& a_grids,
                             const ProblemDomain&     a_domain,
                             const int                a_graphGhost,
                             const int                a_dataGhost) const
{
  CH_TIME("PhaseRealm::fillCopiedLayout");
  if (m_verbose) {
    pout() << "PhaseRealm::fillCopiedLayout" << endl;
  }

  // Every data cell lies in the region the boxes are classified over.
  CH_assert(a_dataGhost <= a_graphGhost);

  const EBLevelGrid* source = nullptr;

  // Any of the realm's levels can be the source, the coarsest included, although the copied layouts are only ever
  // filled on the domains of the finer levels and of the half steps between ratio-4 levels.
  for (int lvl = 0; lvl <= m_finestLevel && source == nullptr; lvl++) {
    if (m_domains[lvl].domainBox() == a_domain.domainBox()) {
      source = &(*m_eblg[lvl]);
    }
    else if (!m_eblgFineHalf[lvl].isNull() && m_eblgFineHalf[lvl]->getDomain().domainBox() == a_domain.domainBox()) {
      source = &(*m_eblgFineHalf[lvl]);
    }
  }

  if (source == nullptr) {
    MayDay::Abort("PhaseRealm::fillCopiedLayout - no level to copy from on the requested domain");
  }

  const DisjointBoxLayout& dblSrc   = source->getDBL();
  const EBISLayout&        ebislSrc = source->getEBISL();
  const Interval           interv(0, 0);
  const IntVect            graphGhost = a_graphGhost * IntVect::Unit;
  const IntVect            dataGhost  = a_dataGhost * IntVect::Unit;

  const DataIterator& ditSrc  = dblSrc.dataIterator();
  const DataIterator& ditDst  = a_grids.dataIterator();
  const int           nboxSrc = ditSrc.size();
  const int           nboxDst = ditDst.size();

  // Classify each destination cell from the source boxes: 0 where no source cell lands, 1 under all-regular source
  // boxes, 3 under all-covered ones, 2 under the others.
  LevelData<BaseFab<int>> classSrc(dblSrc, 1, IntVect::Zero);
  LevelData<BaseFab<int>> classDst(a_grids, 1, graphGhost);

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nboxSrc; mybox++) {
    const DataIndex& din = ditSrc[mybox];

    const EBISBox& ebisBox = ebislSrc[din];

    classSrc[din].setVal(ebisBox.isAllRegular() ? 1 : (ebisBox.isAllCovered() ? 3 : 2));
  }

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nboxDst; mybox++) {
    classDst[ditDst[mybox]].setVal(0);
  }

  classSrc.copyTo(interv, classDst, interv);

  // The source graph, aliased box by box.
  EBGraphFactory graphFact(a_domain);

  LevelData<EBGraph> graphSrc(dblSrc, 1, IntVect::Zero, graphFact);
  LevelData<EBGraph> graphCopy(a_grids, 1, graphGhost, graphFact);

  a_graph.define(a_grids, 1, graphGhost, graphFact);

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nboxSrc; mybox++) {
    const DataIndex& din = ditSrc[mybox];

    graphSrc[din] = ebislSrc[din].getEBGraph();
  }

  // Each destination box is whole and covered when no regular or cut source cell lands in it, whole and regular when
  // regular source cells fill it, and copied cell by cell otherwise. The cells no source reaches lie outside the finer
  // level, where the copy is not read.
  constexpr int wholeCovered = 0;
  constexpr int wholeRegular = 1;
  constexpr int copied       = 2;

  std::vector<int> boxKind(nboxDst, copied);

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nboxDst; mybox++) {
    const DataIndex& din = ditDst[mybox];

    const Box           region = grow(a_grids[din], a_graphGhost) & a_domain;
    const BaseFab<int>& cls    = classDst[din];

    bool anyOpen    = false;
    bool allRegular = true;

    for (BoxIterator bit(region); bit.ok(); ++bit) {
      const int c = cls(bit(), 0);

      anyOpen    = anyOpen || (c == 1) || (c == 2);
      allRegular = allRegular && (c == 1);
    }

    boxKind[mybox] = !anyOpen ? wholeCovered : (allRegular ? wholeRegular : copied);
  }

  // A box copied cell by cell starts all regular. EBGraph::copy skips a covered source while the destination is all
  // covered, and the per-cell graph that a later cut source makes it allocate starts out regular, which would turn the
  // skipped cells regular. From an all-regular start every copy order gives the source's states. A whole covered box
  // starts covered, so that the covered sources copied into it change nothing and store nothing a cell at a time.
#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nboxDst; mybox++) {
    if (boxKind[mybox] == wholeCovered) {
      graphCopy[ditDst[mybox]].setToAllCovered();
    }
    else {
      graphCopy[ditDst[mybox]].setToAllRegular();
    }
  }

  graphSrc.copyTo(interv, graphCopy, interv);

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nboxDst; mybox++) {
    const DataIndex& din = ditDst[mybox];

    if (boxKind[mybox] == wholeCovered) {
      a_graph[din].setToAllCovered();
    }
    else if (boxKind[mybox] == wholeRegular) {
      a_graph[din].setToAllRegular();
    }
    else {
      a_graph[din] = graphCopy[din];
    }
  }

  // The data, built over the graph and copied.
  EBDataFactory dataFact;

  LevelData<EBData> dataSrc(dblSrc, 1, IntVect::Zero, dataFact);

  a_data.define(a_grids, 1, dataGhost, dataFact);

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nboxSrc; mybox++) {
    const DataIndex& din = ditSrc[mybox];

    dataSrc[din] = ebislSrc[din].getEBData();
  }

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nboxDst; mybox++) {
    const DataIndex& din = ditDst[mybox];

    const Box localBox = grow(a_grids[din], a_dataGhost) & a_domain;

    a_data[din].defineVoFData(a_graph[din], localBox);
    a_data[din].defineFaceData(a_graph[din], localBox);
  }

  dataSrc.copyTo(interv, a_data, interv);
}

void
PhaseRealm::defineMultigridLevels(const Vector<Vector<DisjointBoxLayout>>& a_mgGrids,
                                  const Vector<Vector<ProblemDomain>>&     a_mgDomains)
{
  CH_TIME("PhaseRealm::defineMultigridLevels");
  if (m_verbose) {
    pout() << "PhaseRealm::defineMultigridLevels" << endl;
  }

  CH_assert(a_mgGrids.size() == a_mgDomains.size());

  constexpr int mgRefRat = 2;

  const int numAmrLevels = static_cast<int>(a_mgGrids.size());

  m_mgLevelGrids.resize(0);
  m_mgLevelGridsCoFi.resize(0);
  m_mgCoarAve.resize(0);
  m_mgFinerVoFs.resize(0);
  m_mgLevelRestrict.resize(0);
  m_mgLevelProlong.resize(0);

  if (m_isDefined && numAmrLevels > 0) {
    CH_assert(numAmrLevels == 1 + m_finestLevel);

    m_mgLevelGrids.resize(numAmrLevels);
    m_mgLevelGridsCoFi.resize(numAmrLevels);

    for (int lvl = 0; lvl < numAmrLevels; lvl++) {
      const int numMgLevels = static_cast<int>(a_mgGrids[lvl].size());

      CH_assert(numMgLevels >= 1);
      CH_assert(static_cast<int>(a_mgDomains[lvl].size()) == numMgLevels);
      CH_assert(a_mgDomains[lvl][0] == m_domains[lvl]);

      m_mgLevelGrids[lvl].resize(numMgLevels);
      m_mgLevelGridsCoFi[lvl].resize(numMgLevels - 1);

      m_mgLevelGrids[lvl][0] = m_eblg[lvl];

      for (int img = 1; img < numMgLevels; img++) {
        m_mgLevelGrids[lvl][img] = RefCountedPtr<EBLevelGrid>(
          new EBLevelGrid(a_mgGrids[lvl][img], a_mgDomains[lvl][img], m_numGhostCells, &(*m_ebis)));
      }

      for (int img = 0; img < numMgLevels - 1; img++) {
        m_mgLevelGridsCoFi[lvl][img] = RefCountedPtr<EBLevelGrid>(new EBLevelGrid());

        coarsen(*m_mgLevelGridsCoFi[lvl][img], *m_mgLevelGrids[lvl][img], mgRefRat);
      }
    }
  }
}

void
PhaseRealm::defineMultigridCoarseAverage()
{
  CH_TIME("PhaseRealm::defineMultigridCoarseAverage");
  if (m_verbose) {
    pout() << "PhaseRealm::defineMultigridCoarseAverage" << endl;
  }

  constexpr int mgRefRat = 2;

  // Only multigrid solvers use this, so it stays empty unless one registered the multigrid operator.
  const int numAmrLevels = this->queryOperator(s_eb_multigrid) ? static_cast<int>(m_mgLevelGrids.size()) : 0;

  m_mgCoarAve.resize(numAmrLevels);

  for (int lvl = 0; lvl < numAmrLevels; lvl++) {
    const int numMgLevels = static_cast<int>(m_mgLevelGrids[lvl].size());

    m_mgCoarAve[lvl].resize(numMgLevels - 1);

    for (int img = 0; img < numMgLevels - 1; img++) {
      m_mgCoarAve[lvl][img] = RefCountedPtr<EBCoarAve>(new EBCoarAve(*m_mgLevelGrids[lvl][img],
                                                                     *m_mgLevelGrids[lvl][img + 1],
                                                                     *m_mgLevelGridsCoFi[lvl][img],
                                                                     *m_mgFinerVoFs[lvl][img],
                                                                     mgRefRat));
    }
  }
}

void
PhaseRealm::defineMultigridLevelTransfer()
{
  CH_TIME("PhaseRealm::defineMultigridLevelTransfer");
  if (m_verbose) {
    pout() << "PhaseRealm::defineMultigridLevelTransfer" << endl;
  }

  constexpr int mgRefRat = 2;

  // Only multigrid solvers use this, so it stays empty unless one registered the multigrid operator.
  const int numAmrLevels = this->queryOperator(s_eb_multigrid) ? static_cast<int>(m_mgLevelGrids.size()) : 0;

  m_mgLevelRestrict.resize(numAmrLevels);
  m_mgLevelProlong.resize(numAmrLevels);

  for (int lvl = 0; lvl < numAmrLevels; lvl++) {
    const int numMgLevels = static_cast<int>(m_mgLevelGrids[lvl].size());

    m_mgLevelRestrict[lvl].resize(numMgLevels - 1);
    m_mgLevelProlong[lvl].resize(numMgLevels - 1);

    for (int img = 0; img < numMgLevels - 1; img++) {
      const EBLevelGrid& eblgFine = *m_mgLevelGrids[lvl][img];
      const EBLevelGrid& eblgCoFi = *m_mgLevelGridsCoFi[lvl][img];
      const EBLevelGrid& eblgCoar = *m_mgLevelGrids[lvl][img + 1];

      m_mgLevelRestrict[lvl][img] = RefCountedPtr<EBMGRestrict>(
        new EBMGRestrict(eblgFine, eblgCoFi, eblgCoar, *m_mgFinerVoFs[lvl][img], mgRefRat));
      m_mgLevelProlong[lvl][img] = RefCountedPtr<EBMGProlong>(new EBMGProlong(eblgFine, eblgCoFi, eblgCoar, mgRefRat));
    }
  }
}

void
PhaseRealm::defineFinerVoFs(const int a_lmin)
{
  CH_TIME("PhaseRealm::defineFinerVoFs");
  if (m_verbose) {
    pout() << "PhaseRealm::defineFinerVoFs" << endl;
  }

  m_finerVoFs.resize(1 + m_finestLevel);

  for (int lvl = std::max(1, a_lmin); lvl <= m_finestLevel; lvl++) {
    m_finerVoFs[lvl - 1] = RefCountedPtr<LayoutData<BaseIVFAB<Vector<VolIndex>>>>(
      new LayoutData<BaseIVFAB<Vector<VolIndex>>>());

    this->fillFinerVoFs(*m_finerVoFs[lvl - 1], *m_eblgCoFi[lvl - 1], *m_eblg[lvl], m_refinementRatios[lvl - 1]);
  }
}

void
PhaseRealm::defineMultigridFinerVoFs()
{
  CH_TIME("PhaseRealm::defineMultigridFinerVoFs");
  if (m_verbose) {
    pout() << "PhaseRealm::defineMultigridFinerVoFs" << endl;
  }

  constexpr int mgRefRat = 2;

  // Only multigrid solvers use this, so it stays empty unless one registered the multigrid operator.
  const int numAmrLevels = this->queryOperator(s_eb_multigrid) ? static_cast<int>(m_mgLevelGrids.size()) : 0;

  m_mgFinerVoFs.resize(numAmrLevels);

  for (int lvl = 0; lvl < numAmrLevels; lvl++) {
    const int numMgLevels = static_cast<int>(m_mgLevelGrids[lvl].size());

    m_mgFinerVoFs[lvl].resize(numMgLevels - 1);

    for (int img = 0; img < numMgLevels - 1; img++) {
      m_mgFinerVoFs[lvl][img] = RefCountedPtr<LayoutData<BaseIVFAB<Vector<VolIndex>>>>(
        new LayoutData<BaseIVFAB<Vector<VolIndex>>>());

      this->fillFinerVoFs(*m_mgFinerVoFs[lvl][img], *m_mgLevelGridsCoFi[lvl][img], *m_mgLevelGrids[lvl][img], mgRefRat);
    }
  }
}

void
PhaseRealm::fillFinerVoFs(LayoutData<BaseIVFAB<Vector<VolIndex>>>& a_finerVoFs,
                          const EBLevelGrid&                       a_eblgCoFi,
                          const EBLevelGrid&                       a_eblgFine,
                          const int                                a_refRat) const
{
  CH_TIME("PhaseRealm::fillFinerVoFs");
  if (m_verbose) {
    pout() << "PhaseRealm::fillFinerVoFs" << endl;
  }

  const DisjointBoxLayout& dblCoFi   = a_eblgCoFi.getDBL();
  const EBISLayout&        ebislCoFi = a_eblgCoFi.getEBISL();
  const EBISLayout&        ebislFine = a_eblgFine.getEBISL();

  a_finerVoFs.define(dblCoFi);

  const DataIterator& dit  = dblCoFi.dataIterator();
  const int           nbox = dit.size();

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    const EBISBox& ebisBoxCoFi = ebislCoFi[din];
    const Box      region      = VofUtils::getCoarseCellsWithAllFinerCells(ebislFine[din].getRegion(), a_refRat) &
                       ebisBoxCoFi.getRegion();
    const IntVectSet irregCells = ebisBoxCoFi.getIrregIVS(region);

    BaseIVFAB<Vector<VolIndex>>& finerVoFs = a_finerVoFs[din];

    finerVoFs.define(irregCells, ebisBoxCoFi.getEBGraph(), 1);

    for (VoFIterator vofit(irregCells, ebisBoxCoFi.getEBGraph()); vofit.ok(); ++vofit) {
      finerVoFs(vofit(), 0) = VofUtils::getFinerVoFs(vofit(), ebislFine, din, a_refRat);
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
  m_domainFaceIterLo.resize(1 + m_finestLevel);
  m_domainFaceIterHi.resize(1 + m_finestLevel);
  m_domainFaceIterTanGhostLo.resize(1 + m_finestLevel);
  m_domainFaceIterTanGhostHi.resize(1 + m_finestLevel);

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
    m_domainFaceIterLo[lvl] = RefCountedPtr<LayoutData<std::array<FaceIterator, SpaceDim>>>(
      new LayoutData<std::array<FaceIterator, SpaceDim>>(m_grids[lvl]));
    m_domainFaceIterHi[lvl] = RefCountedPtr<LayoutData<std::array<FaceIterator, SpaceDim>>>(
      new LayoutData<std::array<FaceIterator, SpaceDim>>(m_grids[lvl]));
    m_domainFaceIterTanGhostLo[lvl] = RefCountedPtr<LayoutData<std::array<FaceIterator, SpaceDim>>>(
      new LayoutData<std::array<FaceIterator, SpaceDim>>(m_grids[lvl]));
    m_domainFaceIterTanGhostHi[lvl] = RefCountedPtr<LayoutData<std::array<FaceIterator, SpaceDim>>>(
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

      // Domain-face iterators over the cells that lie next to each domain side, regular cells included. The first
      // pair uses the cells of the valid box, the second the valid box grown by one cell in the tangential directions.
      // The cell strips are empty for boxes that do not touch that side.
      std::array<FaceIterator, SpaceDim>& domainFaceIterLo         = (*m_domainFaceIterLo[lvl])[din];
      std::array<FaceIterator, SpaceDim>& domainFaceIterHi         = (*m_domainFaceIterHi[lvl])[din];
      std::array<FaceIterator, SpaceDim>& domainFaceIterTanGhostLo = (*m_domainFaceIterTanGhostLo[lvl])[din];
      std::array<FaceIterator, SpaceDim>& domainFaceIterTanGhostHi = (*m_domainFaceIterTanGhostHi[lvl])[din];

      const Box& domainBox = m_domains[lvl].domainBox();

      for (int dir = 0; dir < SpaceDim; dir++) {
        Box tanBox = cellBox;
        for (int tanDir = 0; tanDir < SpaceDim; tanDir++) {
          if (tanDir != dir) {
            tanBox.grow(tanDir, 1);
          }
        }
        tanBox &= domainBox;

        const Box stripLo = adjCellLo(domainBox, dir, -1);
        const Box stripHi = adjCellHi(domainBox, dir, -1);

        domainFaceIterLo[dir].define(IntVectSet(stripLo & cellBox), ebgraph, dir, FaceStop::AllBoundaryOnly);
        domainFaceIterHi[dir].define(IntVectSet(stripHi & cellBox), ebgraph, dir, FaceStop::AllBoundaryOnly);

        domainFaceIterTanGhostLo[dir].define(IntVectSet(stripLo & tanBox), ebgraph, dir, FaceStop::AllBoundaryOnly);
        domainFaceIterTanGhostHi[dir].define(IntVectSet(stripHi & tanBox), ebgraph, dir, FaceStop::AllBoundaryOnly);
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
        m_coarAve[lvl] = RefCountedPtr<EBCoarAve>(new EBCoarAve(*m_eblg[lvl],
                                                                *m_eblg[lvl - 1],
                                                                *m_eblgCoFi[lvl - 1],
                                                                *m_finerVoFs[lvl - 1],
                                                                m_refinementRatios[lvl - 1]));
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
  m_mgRestrict.resize(1 + m_finestLevel);
  m_mgProlong.resize(1 + m_finestLevel);

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

        // Restriction onto the coarsened fine grids, and prolongation from the coarse level.
        m_mgRestrict[lvl] = RefCountedPtr<EBMGRestrict>(
          new EBMGRestrict(eblgFine, eblgCoFi, eblgCoFi, *m_finerVoFs[lvl - 1], m_refinementRatios[lvl - 1]));
        m_mgProlong[lvl] = RefCountedPtr<EBMGProlong>(
          new EBMGProlong(eblgFine, eblgCoFi, eblgCoar, m_refinementRatios[lvl - 1]));
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
          new EBCoarseToFineInterp(*m_eblg[lvl],
                                   *m_eblgCoFi[lvl - 1],
                                   *m_eblg[lvl - 1],
                                   m_finerVoFs[lvl - 1],
                                   m_refinementRatios[lvl - 1]));
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
          new EBReflux(*m_eblg[lvl], *m_eblg[lvl + 1], *m_eblgCoFi[lvl], *m_finerVoFs[lvl], m_refinementRatios[lvl]));
      }
    }
  }
}

void
PhaseRealm::defineStateRedistribution(const Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>& a_coveredByFiner,
                                      const Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>& a_coarseFineGhost)
{
  CH_TIME("PhaseRealm::defineStateRedistribution");
  if (m_verbose) {
    pout() << "PhaseRealm::defineStateRedistribution" << endl;
  }

  // The neighbourhoods depend on all levels, so the operator is rebuilt on every regrid.
  m_stateRedistribution = RefCountedPtr<EBStateRedistribution>();

  if (this->queryOperator(s_eb_state_redist)) {
    if (m_finestLevel > 0 && !this->queryOperator(s_eb_mass_redist)) {
      MayDay::Abort("PhaseRealm::defineStateRedistribution -- state redistribution needs the mass redistribution "
                    "operator on AMR hierarchies");
    }

    m_stateRedistribution = RefCountedPtr<EBStateRedistribution>(new EBStateRedistribution(m_eblg,
                                                                                           m_eblgCoFi,
                                                                                           m_refinementRatios,
                                                                                           m_finestLevel,
                                                                                           m_srdOrder,
                                                                                           m_srdTargetVolFrac,
                                                                                           m_massRedistributionOp,
                                                                                           a_coveredByFiner,
                                                                                           a_coarseFineGhost));
  }
}

void
PhaseRealm::defineMassRedistOper(const int a_lmin, const int /*a_regsize*/)
{
  CH_TIME("PhaseRealm::defineMassRedistOper");
  if (m_verbose) {
    pout() << "PhaseRealm::defineMassRedistOper" << endl;
  }

  // TLDR: The operator for redistributing on level l lives on level l.

  const bool doThisOperator = this->queryOperator(s_eb_mass_redist);

  m_massRedistributionOp.resize(1 + m_finestLevel);

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

      m_massRedistributionOp[lvl] = RefCountedPtr<EBMassRedistribution>(
        new EBMassRedistribution(eblgCoar,
                                 eblgCoarsened,
                                 eblg,
                                 eblgRefined,
                                 hasFine ? m_fiCoFineCells[lvl + 1] : RefCountedPtr<LevelData<BaseFab<bool>>>(nullptr),
                                 eblgFine,
                                 refToCoar,
                                 refToFine,
                                 s_redistributionRadius));
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
  m_surfaceDeposition.define(m_eblg, m_eblgCoFi, m_finerVoFs, m_refinementRatios, m_dx, m_probLo, m_finestLevel, 1);
}

void
PhaseRealm::defineMirrorSurfaceData(const int a_lmin)
{
  CH_TIME("PhaseRealm::defineMirrorSurfaceData");

  if (m_verbose) {
    pout() << "PhaseRealm::defineMirrorSurfaceData" << endl;
  }

  // Builds the surface data that mirrored deposition reflects particles across. Every cut cell, and every cell within
  // mirror_band_radius of one, gets a status, a boundary centroid x_c, a normal n_c and a shape operator S_c stored
  // in world coordinates. Three passes per level, with an exchange between them:
  //
  //   Pass A -- seed each cut cell of the valid box with its centroid and normal, marked Planar.
  //   Pass B -- fit S_c by least squares over the 5^D neighbourhood, upgrading the cell to Fitted. Skipped when
  //             PhaseRealm.mirror_curvature_correction is false, which leaves S_c zero and a plane mirror.
  //   Pass C -- extend the record outwards to the band, taking it from the Euclidean-nearest cut cell.
  //
  // Two things are load-bearing and easy to undo by accident. The pass order relative to the exchanges: pass C must
  // extend before the data is exchanged, or a band cell at a patch edge is assigned a cut cell that is not its
  // nearest and no later exchange can repair it. And passes B and C read their neighbours from the FIELD rather than
  // from EBISBox, so a cut cell this level does not own reads back as Status::None instead of contributing geometry
  // the level has no data for.

  const bool doThisOperator = this->queryOperator(s_mirror_deposition);

  m_mirrorSurfaceData.resize(1 + m_finestLevel);
  m_mirrorHasCutCells.resize(1 + m_finestLevel, 0);

  if (!doThisOperator) {
    return;
  }

  // Independent requirements: eb_ghost sizes the region EBISBox is valid over, num_ghost sizes the data the
  // exchanges can deliver. Do not derive one from the other.
  if (m_numEbGhostsCells < 2) {
    MayDay::Error("PhaseRealm::defineMirrorSurfaceData - mirrored deposition needs 'AmrMesh.eb_ghost' >= 2 so that "
                  "EBISBox is valid over the 5^D curvature-fit stencil");
  }

  if (m_numGhostCells < m_mirrorBandRadius) {
    MayDay::Error("PhaseRealm::defineMirrorSurfaceData - mirrored deposition needs 'AmrMesh.num_ghost' >= "
                  "'AmrMesh.mirror_band_radius' so that the exchange delivers the cut-cell data the band reads");
  }

  constexpr int fitStencilRadius = 2;

  // Cut-cell geometry source. EBISBox by default: the implicit function is far more accurate -- almost all of it
  // from the centroid rather than the normal -- but costs O(2*SpaceDim) BaseIF::value calls per cut cell per regrid,
  // which is expensive for a tessellated geometry. See the declaration for the measured trade.
  constexpr bool useImplicitFunctionGeometry = false;

  // Finite-difference step for the gradient, in cells. Large enough not to cancel, small enough not to truncate.
  const Real gradientStep = 1.E-2;

  // A tripwire, not a tuning knob. Composite geometries are not differentiable at a seam, where a one-sided gradient
  // can throw the Newton step onto the wrong facet. Ordinary shifts are ~0.02 dx, so this only fires on a seam.
  const Real centroidShiftLimit = 0.5;

  // Counted apart because they mean different things: the fit counters say the geometry is under-resolved,
  // bandNoCutCell says a level's grids do not cover the band.
  long long fitOk               = 0;
  long long fitTooFewNeighbours = 0;
  long long fitIllConditioned   = 0;
  long long fitCrease           = 0;
  long long bandNoCutCell       = 0;

  // How far the implicit function's geometry is from EBISBox's, so a silent change of source cannot go unnoticed.
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

      // The realm's own cut-cell iterator, so this cannot drift from what every other operator iterates.
      VoFIterator& vofit = (*m_vofIter[lvl])[din];

      for (vofit.reset(); vofit.ok(); ++vofit) {
        const VolIndex& vof = vofit();
        const IntVect&  iv  = vof.gridIndex();

        // A position alone cannot select one of two VoFs, so skip multi-cut cells.
        if (ebisbox.numVoFs(iv) != 1) {
          continue;
        }

        // A cut cell the surface only touches, at a corner or an edge, has no boundary area and a zero normal. It holds
        // no surface to reflect against, and as a source it would hand its zero normal to the cells around it.
        if (!(ebisbox.bndryArea(vof) > 0.0)) {
          continue;
        }

        // Physical position. Do NOT difference raw bndryCentroid values between cells -- that drops the inter-cell
        // offset and inflates the fitted curvature by an order of magnitude.
        RealVect       xc   = m_probLo + Location::position(Location::Cell::Boundary, vof, ebisbox, dx);
        const RealVect ebNc = ebisbox.normal(vof);
        RealVect       nc   = ebNc;

        // Gated on the switch as well as on the implicit function existing: the comparison below evaluates
        // BaseIF::value whether or not the result is used.
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

          // One Newton step onto phi = 0.
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

            // Sign to the EB normal's half-space rather than assuming a convention.
            const Real cosRaw = std::max(-1.0, std::min(1.0, lsfNc.dotProduct(ebNc)));

            // Angle with the sign convention folded out, so a perfect match reads 0. The flip count below says
            // whether a convention difference exists at all.
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
    // Pass B -- fit the shape operator over the 5^D neighbourhood of each cut cell.
    // ---------------------------------------------------------------------------------------------------------
    if (m_mirrorCurvatureCorrection) {

#pragma omp parallel for schedule(runtime) reduction(+ : fitOk, fitTooFewNeighbours, fitIllConditioned, fitCrease)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din = dit[mybox];

        EBCellFAB&     surf    = (*m_mirrorSurfaceData[lvl])[din];
        BaseFab<Real>& surfReg = surf.getSingleValuedFAB();
        const EBISBox& ebisbox = surf.getEBISBox();

        const Box validBox = dbl[din];

        // Clipped to the EBISBox region as well as the FAB: a LevelData ghost box is not clipped to the domain.
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

          // Orthonormal tangent basis, crossing the normal with the axis it is least aligned with so the cross
          // product is never near-degenerate. Local to the fit, and discarded at the end of it.
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

          // Least-squares fit: minimize sum_j |v_j - S u_j|^2 over the neighbours' tangential displacement u_j and
          // tangential normal-difference v_j.
          constexpr int T = MirrorDeposition::numTangents;

          Real A[T][T] = {};
          Real B[T][T] = {};
          Real uNormSq = 0.0;

          int numUsable = 0;

          // Read from the FIELD, never from EBISBox: a cut cell no box on this level owns reads back as
          // Status::None here even though EBISBox can see it.
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
            // Condition test on A, normalized so the threshold is dimensionless. A healthy neighbour count is not
            // sufficient -- collinear cut cells, as on a ridge or a fin, leave A near-singular.
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
              // Solved in closed form rather than by SVD, which would allocate per call inside this loop.
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

              // Symmetrize; the antisymmetric part of the fitted operator is noise.
              for (int a = 0; a < T; a++) {
                for (int b = a + 1; b < T; b++) {
                  const Real avg = 0.5 * (S[a][b] + S[b][a]);

                  S[a][b] = avg;
                  S[b][a] = avg;
                }
              }

              // Crease detector: a facet edge or a corner is not a quadratic patch. Residual made dimensionless by
              // dx so the threshold does not move with refinement.
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

            // S_c annihilates the normal by construction. A wrong lift or tangent basis breaks this and nothing
            // else -- tr(S_c) and the Jacobian stay plausible.
            CH_assert(MirrorDeposition::applyShapeOperator(shape, ni).vectorLength() < 1.E-12);

            fitOk++;
          }
          else {
            // Status stays Planar and S_c stays zero, giving J = 1. A refused fit still reflects, against a plane.
            for (int c = 0; c < MirrorDeposition::numShapeComp; c++) {
              surfReg(iv, MirrorDeposition::compShape + c) = 0.0;
            }
          }
        }
      }

      // Publish the fit to the ghosts pass C reads. Inside the guard: with no fit there is nothing new to publish.
      m_mirrorSurfaceData[lvl]->exchange();
    }

    // ---------------------------------------------------------------------------------------------------------
    // Pass C -- extend the patch outwards to the band. Runs AFTER the exchange above; see the summary.
    // ---------------------------------------------------------------------------------------------------------
#pragma omp parallel for schedule(runtime) reduction(+ : bandNoCutCell)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      EBCellFAB&     surf    = (*m_mirrorSurfaceData[lvl])[din];
      BaseFab<Real>& surfReg = surf.getSingleValuedFAB();
      const EBISBox& ebisbox = surf.getEBISBox();

      const Box validBox = dbl[din];
      const Box readBox  = grow(validBox, m_mirrorBandRadius) & surfReg.box() & ebisbox.getRegion();

      // Own copy of the results, or a band cell written early becomes a source for one written later and the patch
      // propagates outwards instead of being assigned.
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

        // Euclidean-nearest, not a Chebyshev sweep. The two differ exactly where more than one piece of surface is
        // within the band, which is the geometry the patch exists to get right.
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

          // Fixed tie-break, so the assignment does not depend on iteration order or box decomposition.
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
          // Reachable where a level's grids do not cover the band. Particles there simply do not reflect.
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

    // No exchange after pass C: surface data is read only at a particle's own cell.

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

    if (!m_mirrorCurvatureCorrection) {
      pout() << "PhaseRealm::defineMirrorSurfaceData - curvature correction DISABLED; the shape operator is zero and "
                "the scheme is a plane mirror (PhaseRealm.mirror_curvature_correction)"
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

      m_gradientOp[lvl] = RefCountedPtr<EBGradient>(
        new EBGradient(eblg,
                       eblgFine,
                       eblgFiCo,
                       hasFine ? m_fiCoFineCells[lvl + 1] : RefCountedPtr<LevelData<BaseFab<bool>>>(nullptr),
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
        new EBNonConservativeDivergence(*m_eblg[lvl], s_redistributionRadius));
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

const Vector<Vector<RefCountedPtr<EBLevelGrid>>>&
PhaseRealm::getMultigridLevelGrids() const
{
  return m_mgLevelGrids;
}

const Vector<Vector<RefCountedPtr<EBCoarAve>>>&
PhaseRealm::getMultigridCoarseAverage() const
{
  if (!this->queryOperator(s_eb_multigrid)) {
    MayDay::Error("PhaseRealm::getMultigridCoarseAverage - operator not registered!");
  }

  return m_mgCoarAve;
}

const Vector<RefCountedPtr<LayoutData<BaseIVFAB<Vector<VolIndex>>>>>&
PhaseRealm::getFinerVoFs() const
{
  return m_finerVoFs;
}

const Vector<Vector<RefCountedPtr<LayoutData<BaseIVFAB<Vector<VolIndex>>>>>>&
PhaseRealm::getMultigridFinerVoFs() const
{
  if (!this->queryOperator(s_eb_multigrid)) {
    MayDay::Error("PhaseRealm::getMultigridFinerVoFs - operator not registered!");
  }

  return m_mgFinerVoFs;
}

const Vector<RefCountedPtr<EBMGRestrict>>&
PhaseRealm::getMultigridRestriction() const
{
  if (!this->queryOperator(s_eb_multigrid)) {
    MayDay::Error("PhaseRealm::getMultigridRestriction - operator not registered!");
  }

  return m_mgRestrict;
}

const Vector<RefCountedPtr<EBMGProlong>>&
PhaseRealm::getMultigridProlongation() const
{
  if (!this->queryOperator(s_eb_multigrid)) {
    MayDay::Error("PhaseRealm::getMultigridProlongation - operator not registered!");
  }

  return m_mgProlong;
}

const Vector<Vector<RefCountedPtr<EBMGRestrict>>>&
PhaseRealm::getMultigridLevelRestriction() const
{
  if (!this->queryOperator(s_eb_multigrid)) {
    MayDay::Error("PhaseRealm::getMultigridLevelRestriction - operator not registered!");
  }

  return m_mgLevelRestrict;
}

const Vector<Vector<RefCountedPtr<EBMGProlong>>>&
PhaseRealm::getMultigridLevelProlongation() const
{
  if (!this->queryOperator(s_eb_multigrid)) {
    MayDay::Error("PhaseRealm::getMultigridLevelProlongation - operator not registered!");
  }

  return m_mgLevelProlong;
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

Vector<RefCountedPtr<LayoutData<std::array<FaceIterator, SpaceDim>>>>&
PhaseRealm::getDomainFaceIterator(const Side::LoHiSide a_side) const
{
  return (a_side == Side::Lo) ? m_domainFaceIterLo : m_domainFaceIterHi;
}

Vector<RefCountedPtr<LayoutData<std::array<FaceIterator, SpaceDim>>>>&
PhaseRealm::getDomainFaceIteratorWithTangentialGhosts(const Side::LoHiSide a_side) const
{
  return (a_side == Side::Lo) ? m_domainFaceIterTanGhostLo : m_domainFaceIterTanGhostHi;
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

const RefCountedPtr<EBStateRedistribution>&
PhaseRealm::getStateRedistribution() const
{
  if (!this->queryOperator(s_eb_state_redist)) {
    MayDay::Error("PhaseRealm::getStateRedistribution - operator not registered!");
  }

  return m_stateRedistribution;
}

Vector<RefCountedPtr<EBMassRedistribution>>&
PhaseRealm::getMassRedistributionOp() const
{
  if (!this->queryOperator(s_eb_mass_redist)) {
    MayDay::Error("PhaseRealm::getMassRedistributionOp - operator not registered!");
  }

  return m_massRedistributionOp;
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
