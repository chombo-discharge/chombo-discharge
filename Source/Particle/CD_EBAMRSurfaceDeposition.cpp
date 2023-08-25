/* chombo-discharge
 * Copyright © 2023 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBAMRSurfaceDeposition.cpp
  @brief  Implementation of CD_EBAMRSurfaceDeposition.H
  @author Robert Marskar
*/

// Chombo includes
#include <ParmParse.H>
#include <BaseIVFactory.H>
#include <NeighborIterator.H>
#include <CH_Timer.H>

// Our includes
#include <CD_IrregAddOp.H>
#include <CD_EBAMRSurfaceDeposition.H>
#include <CD_NamespaceHeader.H>

EBAMRSurfaceDeposition::EBAMRSurfaceDeposition() noexcept
{
  CH_TIME("EBAMRSurfaceDeposition::EBAMRSurfaceDeposition(weak)");

  m_debug   = false;
  m_verbose = false;
}

EBAMRSurfaceDeposition::EBAMRSurfaceDeposition(const Vector<RefCountedPtr<EBLevelGrid>>& a_ebGrids,
                                               const Vector<RefCountedPtr<EBLevelGrid>>& a_ebGridsCoarsenedFine,
                                               const Vector<RefCountedPtr<EBLevelGrid>>& a_ebGridsRefinedCoar,
                                               const Vector<int>&                        a_refRat,
                                               const Vector<Real>&                       a_dx,
                                               const RealVect&                           a_probLo,
                                               const int                                 a_finestLevel,
                                               const int                                 a_radius) noexcept
{
  CH_TIME("EBAMRSurfaceDeposition::EBAMRSurfaceDeposition(full)");

  this->define(a_ebGrids,
               a_ebGridsCoarsenedFine,
               a_ebGridsRefinedCoar,
               a_refRat,
               a_dx,
               a_probLo,
               a_finestLevel,
               a_radius);
}

EBAMRSurfaceDeposition::~EBAMRSurfaceDeposition() noexcept
{
  CH_TIME("EBAMRSurfaceDeposition::~EBAMRSurfaceDeposition");
}

void
EBAMRSurfaceDeposition::define(const Vector<RefCountedPtr<EBLevelGrid>>& a_ebGrids,
                               const Vector<RefCountedPtr<EBLevelGrid>>& a_ebGridsCoarsenedFine,
                               const Vector<RefCountedPtr<EBLevelGrid>>& a_ebGridsRefinedCoar,
                               const Vector<int>&                        a_refRat,
                               const Vector<Real>&                       a_dx,
                               const RealVect&                           a_probLo,
                               const int                                 a_finestLevel,
                               const int                                 a_radius) noexcept
{
  CH_TIME("EBAMRSurfaceDeposition::define");
  if (m_verbose) {
    pout() << "EBAMRSurfaceDeposition::define" << endl;
  }

  m_debug                = false;
  m_verbose              = false;
  m_ebGrids              = a_ebGrids;
  m_ebGridsCoarsenedFine = a_ebGridsCoarsenedFine;
  m_ebGridsRefinedCoar   = a_ebGridsRefinedCoar;
  m_refRat               = a_refRat;
  m_dx                   = a_dx;
  m_probLo               = a_probLo;
  m_finestLevel          = a_finestLevel;
  m_radius               = a_radius;

  CH_assert(m_finestLevel >= 0);
  CH_assert(m_radius >= 0);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const bool hasFine = lvl < m_finestLevel;
    const bool hasCoar = lvl > 0;

    CH_assert(!(a_ebGrids[lvl].isNull()));

    if (lvl < m_finestLevel) {
      CH_assert(a_refRat[lvl] >= 2);
      CH_assert(a_refRat[lvl] % 2 == 0);
    }

    if (hasCoar) {
      CH_assert(!(a_ebGridsRefinedCoar[lvl].isNull()));
    }
    if (hasFine) {
      CH_assert(!(a_ebGridsCoarsenedFine[lvl].isNull()));
    }
  }

  // Put in debug mode or not.
  ParmParse pp("EBAMRSurfaceDeposition");
  pp.query("debug", m_debug);
  pp.query("verbose", m_verbose);

  this->defineBuffers();
  this->defineDataMotion();
  this->defineDepositionStencils();
  this->defineCoarseToFineStencils();
  this->defineFineToCoarseStencils();

  m_isDefined = true;
}

void
EBAMRSurfaceDeposition::defineBuffers() noexcept
{
  CH_TIMERS("EBAMRSurfaceDeposition::defineBuffers");
  CH_TIMER("EBAMRSurfaceDeposition::defineBuffers::define_level", t1);
  CH_TIMER("EBAMRSurfaceDeposition::defineBuffers::define_refined_coar", t2);
  CH_TIMER("EBAMRSurfaceDeposition::defineBuffers::define_coarsened_fine", t3);
  if (m_verbose) {
    pout() << "EBAMRSurfaceDeposition::defineBuffers" << endl;
  }

  // TLDR: This defines the buffer data that is required for
  //       1. Depositing particles on each AMR level
  //       2. Moving mass deposited beneath an AMR level into the next finer level
  //       3. Moving mass deposited across a refinement boundary into the coarse level.

  constexpr int nComp = 1;

  CH_START(t1);
  m_data.resize(1 + m_finestLevel);
  m_refinedCoarData.resize(1 + m_finestLevel);
  m_coarsenedFineData.resize(1 + m_finestLevel);

  // Define data on this level
  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl    = m_ebGrids[lvl]->getDBL();
    const DataIterator&      dit    = dbl.dataIterator();
    const EBISLayout&        ebisl  = m_ebGrids[lvl]->getEBISL();
    const ProblemDomain&     domain = m_ebGrids[lvl]->getDomain();

    LayoutData<IntVectSet> irregCells(dbl);

    const int nbox = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      const Box      box     = dbl[din];
      const EBISBox& ebisbox = ebisl[din];

      irregCells[din] = ebisbox.getIrregIVS(grow(box, m_radius) & domain);
    }

    m_data[lvl] = RefCountedPtr<LevelData<BaseIVFAB<Real>>>(
      new LevelData<BaseIVFAB<Real>>(dbl, nComp, m_radius * IntVect::Unit, BaseIVFactory<Real>(ebisl, irregCells)));
  }
  CH_STOP(t1);

  // Define refined coarse data
  CH_START(t2);
  for (int lvl = 1; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& refinedCoarDBL   = m_ebGridsRefinedCoar[lvl]->getDBL();
    const DataIterator&      dit              = refinedCoarDBL.dataIterator();
    const EBISLayout&        refinedCoarEBISL = m_ebGridsRefinedCoar[lvl]->getEBISL();

    LayoutData<IntVectSet> irregCells(refinedCoarDBL);

    const int nbox = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      const Box      box     = refinedCoarDBL[din];
      const EBISBox& ebisbox = refinedCoarEBISL[din];

      irregCells[din] = ebisbox.getIrregIVS(box);
    }

    m_refinedCoarData[lvl] = RefCountedPtr<LevelData<BaseIVFAB<Real>>>(
      new LevelData<BaseIVFAB<Real>>(refinedCoarDBL,
                                     nComp,
                                     IntVect::Zero,
                                     BaseIVFactory<Real>(refinedCoarEBISL, irregCells)));
  }
  CH_STOP(t2);

  // Define coarsened fine data
  CH_START(t3);
  for (int lvl = 0; lvl <= m_finestLevel - 1; lvl++) {
    const DisjointBoxLayout& coarsenedFineDBL    = m_ebGridsCoarsenedFine[lvl]->getDBL();
    const EBISLayout&        coarsenedFineEBISL  = m_ebGridsCoarsenedFine[lvl]->getEBISL();
    const ProblemDomain&     coarsenedFineDomain = m_ebGridsCoarsenedFine[lvl]->getDomain();
    const DataIterator&      dit                 = coarsenedFineDBL.dataIterator();

    LayoutData<IntVectSet> irregCells(coarsenedFineDBL);

    const int nbox = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      const Box      box     = coarsenedFineDBL[din];
      const EBISBox& ebisbox = coarsenedFineEBISL[din];

      irregCells[din] = ebisbox.getIrregIVS(grow(box, m_radius) & coarsenedFineDomain);
    }

    m_coarsenedFineData[lvl] = RefCountedPtr<LevelData<BaseIVFAB<Real>>>(
      new LevelData<BaseIVFAB<Real>>(coarsenedFineDBL,
                                     nComp,
                                     m_radius * IntVect::Unit,
                                     BaseIVFactory<Real>(coarsenedFineEBISL, irregCells)));
  }
  CH_STOP(t3);
}

void
EBAMRSurfaceDeposition::defineDataMotion() noexcept
{
  CH_TIME("EBAMRSurfaceDeposition::defineDataMotion");
  if (m_verbose) {
    pout() << "EBAMRSurfaceDeposition::defineDataMotion" << endl;
  }

  // TLDR: Define copiers for moving data between AMR levels and refinement/coarsenings of AMR levels.

  m_validToValidCopiers.resize(1 + m_finestLevel);
  m_copierLevel.resize(1 + m_finestLevel);
  m_copierRefinedCoarToFineNoGhosts.resize(1 + m_finestLevel);
  m_copierCoarsenedFineToCoar.resize(1 + m_finestLevel);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const bool hasCoar = lvl > 0;
    const bool hasFine = lvl < m_finestLevel;

    const EBLevelGrid&       eblg   = *m_ebGrids[lvl];
    const ProblemDomain&     domain = eblg.getDomain();
    const DisjointBoxLayout& dbl    = eblg.getDBL();

    const bool doExchange = true;

    // Valid to valid copier.
    m_validToValidCopiers[lvl].define(dbl, dbl, domain);

    // Define Copier as going from valid -> valid+ghost.
    m_copierLevel[lvl].define(dbl, dbl, domain, m_radius * IntVect::Unit, doExchange);

    // Reverse so that this goes from valid+ghost -> valid.
    m_copierLevel[lvl].reverse();

    if (hasCoar) {
      const DisjointBoxLayout& refinedCoarDBL = m_ebGridsRefinedCoar[lvl]->getDBL();

      m_copierRefinedCoarToFineNoGhosts[lvl].define(refinedCoarDBL, dbl, domain);
    }

    if (hasFine) {
      m_copierCoarsenedFineToCoar[lvl].ghostDefine(m_ebGridsCoarsenedFine[lvl]->getDBL(),
                                                   m_ebGrids[lvl]->getDBL(),
                                                   m_ebGrids[lvl]->getDomain(),
                                                   m_radius * IntVect::Unit);
    }
  }
}

void
EBAMRSurfaceDeposition::defineDepositionStencils() noexcept
{
  CH_TIMERS("EBAMRSurfaceDeposition::defineDepositionStencils");
  CH_TIMER("EBAMRSurfaceDeposition::defineDepositionStencils::stencil_define", t1);
  CH_TIMER("EBAMRSurfaceDeposition::defineDepositionStencils::calc_weights", t2);
  if (m_verbose) {
    pout() << "EBAMRSurfaceDeposition::defineDepositionStencils" << endl;
  }

  constexpr int nComp = 1;

  CH_START(t1);
  m_depositionStencils.resize(1 + m_finestLevel);

  // Define over valid cut-cells (i.e., cut-cells not covered by a finer grid)
  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl   = m_ebGrids[lvl]->getDBL();
    const EBISLayout&        ebisl = m_ebGrids[lvl]->getEBISL();
    const DataIterator&      dit   = dbl.dataIterator();

    m_depositionStencils[lvl] = RefCountedPtr<LayoutData<BaseIVFAB<VoFStencil>>>(
      new LayoutData<BaseIVFAB<VoFStencil>>(dbl));

    const int nbox = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      const Box        box     = dbl[din];
      const EBISBox&   ebisbox = ebisl[din];
      const EBGraph&   ebgraph = ebisbox.getEBGraph();
      const IntVectSet ivs     = ebisbox.getIrregIVS(box);

      (*m_depositionStencils[lvl])[din].define(ivs, ebgraph, 1);
    }
  }
  CH_STOP(t1);

  // Compute deposition weights for current cells.
  CH_START(t2);
  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl    = m_ebGrids[lvl]->getDBL();
    const EBISLayout&        ebisl  = m_ebGrids[lvl]->getEBISL();
    const ProblemDomain&     domain = m_ebGrids[lvl]->getDomain();
    const DataIterator&      dit    = dbl.dataIterator();
    const Real               dx     = m_dx[lvl];
    const int                nbox   = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      const Box      box     = dbl[din];
      const EBISBox& ebisBox = ebisl[din];
      const EBGraph& ebGraph = ebisBox.getEBGraph();

      BaseIVFAB<VoFStencil>& stencils = (*m_depositionStencils[lvl])[din];

      // Build stencils.
      const IntVectSet& ivs     = stencils.getIVS();
      const EBGraph&    ebgraph = stencils.getEBGraph();

      for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit) {
        const VolIndex curVoF = vofit();
        const IntVect  curIV  = curVoF.gridIndex();

        VoFStencil& stencil = stencils(curVoF, 0);
        stencil.clear();

        Real             totalArea = 0.0;
        Vector<VolIndex> stencilVoFs;

        // These are the cells that we are going to deposit into.
        const IntVectSet neighborhood = ebisBox.getIrregIVS(grow(Box(curIV, curIV), m_radius));

        for (IVSIterator ivsIt(neighborhood); ivsIt.ok(); ++ivsIt) {
          const IntVect iv = ivsIt();

          const Vector<VolIndex> vofs = ebisBox.getVoFs(iv);

          for (int i = 0; i < vofs.size(); i++) {
            const Real areaEB = ebisBox.bndryArea(vofs[i]);
            if (areaEB > std::numeric_limits<Real>::min()) {
              totalArea += areaEB * std::pow(dx, SpaceDim - 1);

              stencilVoFs.push_back(vofs[i]);
            }
          }
        }

        if (totalArea > std::numeric_limits<Real>::min()) {
          totalArea = 1.0 / totalArea;
        }
        else {
          totalArea = 0.0;
        }

        for (int i = 0; i < stencilVoFs.size(); i++) {
          if (domain.contains(stencilVoFs[i].gridIndex())) {
            stencil.add(stencilVoFs[i], totalArea);
          }
        }
      };
    }
  }
  CH_START(t2);
}

void
EBAMRSurfaceDeposition::defineCoarseToFineStencils() noexcept
{
  CH_TIME("EBAMRSurfaceDeposition::defineCoarseToFineStencils");
  if (m_verbose) {
    pout() << "EBAMRSurfaceDeposition::defineCoarseToFineStencils" << endl;
  }

  m_coarseToFineStencils.resize(1 + m_finestLevel);

  for (int lvl = 1; lvl <= m_finestLevel; lvl++) {
    const Real dxCoar   = m_dx[lvl - 1];
    const Real dxFine   = m_dx[lvl];
    const Real dxFactor = std::pow(dxCoar / dxFine, SpaceDim - 1);

    const DisjointBoxLayout& dblCoar = m_ebGrids[lvl - 1]->getDBL();
    const DisjointBoxLayout& dblFine = m_ebGridsRefinedCoar[lvl]->getDBL();

    const EBISLayout& ebislCoar = m_ebGrids[lvl - 1]->getEBISL();
    const EBISLayout& ebislFine = m_ebGridsRefinedCoar[lvl]->getEBISL();

    m_coarseToFineStencils[lvl] = RefCountedPtr<LayoutData<BaseIVFAB<VoFStencil>>>(
      new LayoutData<BaseIVFAB<VoFStencil>>(dblCoar));

    const DataIterator& dit  = dblCoar.dataIterator();
    const int           nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      const Box boxFine = dblFine[din];
      const Box boxCoar = dblCoar[din];

      CH_assert(refine(boxCoar, m_refRat[lvl - 1]) == boxFine);

      const EBISBox& ebisBoxCoar = ebislCoar[din];
      const EBISBox& ebisBoxFine = ebislFine[din];

      const EBGraph&   ebGraphCoar = ebisBoxCoar.getEBGraph();
      const IntVectSet irregCoar   = ebisBoxCoar.getIrregIVS(dblCoar[din]);

      BaseIVFAB<VoFStencil>& interpStencils = (*m_coarseToFineStencils[lvl])[din];

      interpStencils.define(irregCoar, ebGraphCoar, 1);

      for (VoFIterator vofit(irregCoar, ebGraphCoar); vofit.ok(); ++vofit) {
        const VolIndex& coarVoF = vofit();

        VoFStencil& stencil = interpStencils(coarVoF, 0);

        // Figure out which cut-cells on the fine grid comes from refining the current one on the coarse grid
        const Vector<VolIndex> refinedCoarVoFs = ebislCoar.refine(coarVoF, m_refRat[lvl - 1], din);

        Vector<VolIndex> fineIrregVoFs;
        for (int i = 0; i < refinedCoarVoFs.size(); i++) {
          if (ebisBoxFine.isIrregular(refinedCoarVoFs[i].gridIndex())) {
            fineIrregVoFs.push_back(refinedCoarVoFs[i]);
          }
        }

        Real fineArea = 0.0;
        for (int i = 0; i < fineIrregVoFs.size(); i++) {
          fineArea += ebisBoxFine.bndryArea(fineIrregVoFs[i]);
        }

        if (fineArea < std::numeric_limits<Real>::min()) {
          for (int i = 0; i < fineIrregVoFs.size(); i++) {
            stencil.add(fineIrregVoFs[i], 0.0);
          }
        }
        else {
          const Real coarArea   = ebisBoxCoar.bndryArea(coarVoF);
          const Real stenWeight = coarArea * dxFactor / fineArea;

          for (int i = 0; i < fineIrregVoFs.size(); i++) {
            stencil.add(fineIrregVoFs[i], stenWeight);
          }
        }
      }
    }
  }
}

void
EBAMRSurfaceDeposition::defineFineToCoarseStencils() noexcept
{
  CH_TIME("EBAMRSurfaceDeposition::defineFineToCoarseStencils");
  if (m_verbose) {
    pout() << "EBAMRSurfaceDeposition::defineFineToCoarseStencils" << endl;
  }

  m_fineToCoarseStencils.resize(1 + m_finestLevel);

  for (int lvl = 1; lvl <= m_finestLevel; lvl++) {
    const Real dxCoar   = m_dx[lvl - 1];
    const Real dxFine   = m_dx[lvl];
    const Real dxFactor = std::pow(dxFine / dxCoar, SpaceDim - 1);

    const DisjointBoxLayout& dblCoar = m_ebGridsCoarsenedFine[lvl - 1]->getDBL();
    const DisjointBoxLayout& dblFine = m_ebGrids[lvl]->getDBL();

    const EBISLayout& ebislCoar = m_ebGridsCoarsenedFine[lvl - 1]->getEBISL();
    const EBISLayout& ebislFine = m_ebGrids[lvl]->getEBISL();

    const ProblemDomain& domainFine = m_ebGrids[lvl]->getDomain();

    m_fineToCoarseStencils[lvl] = RefCountedPtr<LayoutData<BaseIVFAB<VoFStencil>>>(
      new LayoutData<BaseIVFAB<VoFStencil>>(dblFine));

    const DataIterator& dit = dblFine.dataIterator();

    const int nbox = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      const Box boxFine = dblFine[din];
      const Box boxCoar = dblCoar[din];

      const EBISBox& ebisBoxCoar = ebislCoar[din];
      const EBISBox& ebisBoxFine = ebislFine[din];
      const EBGraph& ebGraphFine = ebisBoxFine.getEBGraph();

      CH_assert(refine(boxCoar, m_refRat[lvl - 1]) == boxFine);

      // Build the CFIVS around this patch; the stencils are defined along this patch.
      IntVectSet irregCFIVS = ebisBoxFine.getIrregIVS((grow(boxFine, m_radius) & domainFine));

      NeighborIterator nit(dblFine);
      for (nit.begin(din); nit.ok(); ++nit) {
        irregCFIVS -= dblFine[nit()];
      }
      irregCFIVS -= boxFine;

      // Define data holder over all cut-cells that are ghost cells on the coarse-fine interface.
      BaseIVFAB<VoFStencil>& coarsenStencils = (*m_fineToCoarseStencils[lvl])[din];

      coarsenStencils.define(irregCFIVS, ebGraphFine, 1);

      // Define the stencils -- this ensures that we do conservative coarsening of the data.
      for (VoFIterator vofit(irregCFIVS, ebGraphFine); vofit.ok(); ++vofit) {
        const VolIndex& fineGhostVoF = vofit();
        const VolIndex& coarVoF      = ebislFine.coarsen(fineGhostVoF, m_refRat[lvl - 1], din);

        CH_assert(ebisBoxFine.isIrregular(fineGhostVoF.gridIndex()));
        CH_assert(ebisBoxCoar.isIrregular(coarVoF.gridIndex()));

        VoFStencil& stencil = coarsenStencils(fineGhostVoF, 0);
        stencil.clear();

        const Real fineArea = ebisBoxFine.bndryArea(fineGhostVoF);
        const Real coarArea = ebisBoxCoar.bndryArea(coarVoF);

        if (fineArea < std::numeric_limits<Real>::min() || coarArea < std::numeric_limits<Real>::min()) {
          stencil.add(coarVoF, 0.0);
        }
        else {
          const Real weight = fineArea * dxFactor / coarArea;

          stencil.add(coarVoF, weight);
        }
      }
    }
  }
}

void
EBAMRSurfaceDeposition::addInvalidCoarseDataToFineData() const noexcept
{
  CH_TIME("EBAMRSurfaceDeposition::addInvalidCoarseDataToFineData");
  if (m_verbose) {
    pout() << "EBAMRSurfaceDeposition::addInvalidCoarseDataToFineData" << endl;
  }

  for (int lvl = 1; lvl <= m_finestLevel; lvl++) {

    const DisjointBoxLayout& dblCoar = m_ebGrids[lvl - 1]->getDBL();
    const DataIterator&      dit     = dblCoar.dataIterator();
    const int                nbox    = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      BaseIVFAB<Real>&       fineData = (*m_refinedCoarData[lvl])[din];
      const BaseIVFAB<Real>& coarData = (*m_data[lvl - 1])[din];

      const BaseIVFAB<VoFStencil>& stencils = (*m_coarseToFineStencils[lvl])[din];

      const IntVectSet irregCoar   = stencils.getIVS();
      const EBGraph&   ebGraphCoar = stencils.getEBGraph();

      fineData.setVal(0.0);

      // Piecewise constant conservative interpolation of coarse-grid data to the fine grid.
      for (VoFIterator vofit(irregCoar, ebGraphCoar); vofit.ok(); ++vofit) {
        const VolIndex&   coarVoF = vofit();
        const VoFStencil& stencil = stencils(coarVoF, 0);

        const Real& coarVal = coarData(coarVoF, 0);

        for (int i = 0; i < stencil.size(); i++) {
          const VolIndex& fineVoF    = stencil.vof(i);
          const Real&     fineWeight = stencil.weight(i);

          fineData(fineVoF, 0) += coarVal * fineWeight;
        }
      }
    }

    // Add the interpolated data to the fine-grid data.
    const Interval variables = Interval(0, 0);
    m_refinedCoarData[lvl]->copyTo(variables,
                                   *m_data[lvl],
                                   variables,
                                   m_copierRefinedCoarToFineNoGhosts[lvl],
                                   IrregAddOp());
  }
}

void
EBAMRSurfaceDeposition::addFineGhostDataToValidCoarData() const noexcept
{
  CH_TIME("EBAMRSurfaceDeposition::addFineGhostDataToValidCoarData");
  if (m_verbose) {
    pout() << "EBAMRSurfaceDeposition::addFineGhostDataToValidCoarData" << endl;
  }

  for (int lvl = 0; lvl < m_finestLevel; lvl++) {
    const DisjointBoxLayout& dblFine = m_ebGrids[lvl + 1]->getDBL();
    const DataIterator&      dit     = dblFine.dataIterator();
    const int                nbox    = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      const BaseIVFAB<VoFStencil>& stencils = (*m_fineToCoarseStencils[lvl + 1])[din];

      const IntVectSet& ghostsFine = stencils.getIVS();
      const EBGraph&    graphFine  = stencils.getEBGraph();

      BaseIVFAB<Real>&       coarData = (*m_coarsenedFineData[lvl])[din];
      const BaseIVFAB<Real>& fineData = (*m_data[lvl + 1])[din];

      coarData.setVal(0.0);

      for (VoFIterator vofit(ghostsFine, graphFine); vofit.ok(); ++vofit) {
        const VolIndex&   fineGhost = vofit();
        const VoFStencil& stencil   = stencils(fineGhost, 0);

        const Real fineVal = fineData(fineGhost, 0);

        for (int i = 0; i < stencil.size(); i++) {
          const VolIndex& coarVoF    = stencil.vof(i);
          const Real&     coarWeight = stencil.weight(i);

          coarData(coarVoF, 0) += fineVal * coarWeight;
        }
      }
    }

    const Interval interv(0, 0);

    m_coarsenedFineData[lvl]->copyTo(interv, *m_data[lvl], interv, m_copierCoarsenedFineToCoar[lvl], IrregAddOp());
  }
}

#include <CD_NamespaceFooter.H>
