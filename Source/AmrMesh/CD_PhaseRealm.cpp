/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_PhaseRealm.cpp
  @brief  Implementation of CD_PhaseRealm.H
  @author Robert Marskar
*/

// Chombo includes
#include <EBArith.H>
#include <ParmParse.H>
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
#include <CD_NamespaceHeader.H>

PhaseRealm::PhaseRealm()
{
  CH_TIME("PhaseRealm::PhaseRealm");

  // Default settings
  m_isDefined = false;
  m_profile   = false;
  m_verbose   = false;

  this->registerOperator(s_eb_gradient);
  this->registerOperator(s_eb_irreg_interp);

  // Adding this for debugging purposes.
  ParmParse pp("PhaseRealm");
  pp.query("profile", m_profile);
  pp.query("verbosity", m_verbose);
}

PhaseRealm::~PhaseRealm()
{}

void
PhaseRealm::define(const Vector<DisjointBoxLayout>&      a_grids,
                   const Vector<ProblemDomain>&          a_domains,
                   const Vector<int>&                    a_refRat,
                   const Vector<Real>&                   a_dx,
                   const RealVect                        a_probLo,
                   const int                             a_finestLevel,
                   const int                             a_ebGhost,
                   const int                             a_numGhost,
                   const int                             a_lsfGhost,
                   const int                             a_redistRad,
                   const int                             a_mgInterpOrder,
                   const int                             a_mgInterpRadius,
                   const int                             a_mgInterpWeight,
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
  m_coarAve.resize(0);
  m_multigridInterpolator.resize(0);
  m_ebFineInterp.resize(0);
  m_ebReflux.resize(0);
  m_redistributionOp.resize(0);
  m_gradientOp.resize(0);
  m_levelset.resize(0);
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
PhaseRealm::registerOperator(const std::string a_operator)
{
  CH_TIME("PhaseRealm::registerOperator");
  if (m_verbose) {
    pout() << "PhaseRealm::registerOperator" << endl;
  }

  // These are the supported operators - issue an error if we ask for something that is not supported.
  if (!(a_operator.compare(s_eb_coar_ave) == 0 || a_operator.compare(s_eb_fill_patch) == 0 ||
        a_operator.compare(s_eb_fine_interp) == 0 || a_operator.compare(s_eb_flux_reg) == 0 ||
        a_operator.compare(s_eb_redist) == 0 || a_operator.compare(s_noncons_div) == 0 ||
        a_operator.compare(s_eb_gradient) == 0 || a_operator.compare(s_particle_mesh) == 0 ||
        a_operator.compare(s_eb_irreg_interp) == 0 || a_operator.compare(s_eb_multigrid) == 0 ||
        a_operator.compare(s_levelset) == 0)) {

    const std::string str = "PhaseRealm::registerOperator - unknown operator '" + a_operator + "' requested";
    MayDay::Error(str.c_str());
  }

  if (!this->queryOperator(a_operator)) {
    m_operatorMap.emplace(a_operator, true);
  }
}

bool
PhaseRealm::queryOperator(const std::string a_operator) const
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

  for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++) {

    m_vofIter[lvl] = RefCountedPtr<LayoutData<VoFIterator>>(new LayoutData<VoFIterator>(m_grids[lvl]));

    const DisjointBoxLayout& dbl = m_grids[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const int nbox = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      VoFIterator& vofit = (*m_vofIter[lvl])[din];

      const Box&        cellBox = m_grids[lvl].get(din);
      const EBISBox&    ebisbox = m_ebisl[lvl][din];
      const EBGraph&    ebgraph = ebisbox.getEBGraph();
      const IntVectSet& irreg   = ebisbox.getIrregIVS(cellBox);

      vofit.define(irreg, ebgraph);
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
          auto kernel = [&](const IntVect& iv) -> void {
            const RealVect pos = m_probLo + (0.5 * RealVect::Unit + RealVect(iv)) * dx;

            fab(iv, comp) = m_baseif->value(pos);
          };

          BoxLoops::loop(bx, kernel);
        }
        else {
          fab.setVal(minVal, comp);
        }
      }
    }
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

    const int     comps  = SpaceDim;
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
PhaseRealm::defineFluxReg(const int a_lmin, const int a_regsize)
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

    const int comps = a_regsize;

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
PhaseRealm::defineRedistOper(const int a_lmin, const int a_regsize)
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
PhaseRealm::defineGradSten(const int a_lmin)
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
PhaseRealm::defineNonConservativeDivergence(const int a_lmin)
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

Vector<RefCountedPtr<LayoutData<VoFIterator>>>&
PhaseRealm::getVofIterator() const
{
  return m_vofIter;
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

#include <CD_NamespaceFooter.H>
