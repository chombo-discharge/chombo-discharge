/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBAMRParticleMesh.cpp
  @brief  Implementation of CD_EBAMRParticleMesh.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>
#include <NeighborIterator.H>

// Our includes
#include <CD_EBAMRParticleMesh.H>
#include <CD_NamespaceHeader.H>

// These should be removed when we merge
#warning "Must update all solvers that use deposition and set transition to default"
#warning "EBParticleMesh infrastructure uses temporary buffers, but this should be prealloacted memory."
#warning "The deposition routines would probably do well to update using non-nested box loops"
#warning "The deposition routines should be single-signature only -- with a return function on the component"
#warning "CIC and refinement factor 4 using the transition method. This should be a mask width of 2!"
#warning "EBAMRSurfaceDeposition should also be updated with the new signatures"
#warning "McPhoto, ItoSolver, and CdrSolver should not be allowed to use EBParticleMesh directly"
#warning "CheckDocs.py has triggered"

EBAMRParticleMesh::EBAMRParticleMesh()
{
  CH_TIME("EBAMRParticleMesh::EBAMRParticleMesh()");

  m_isDefined = false;
}

EBAMRParticleMesh::EBAMRParticleMesh(const Vector<RefCountedPtr<EBLevelGrid>>& a_eblgs,
                                     const Vector<int>&                        a_refRat,
                                     const Vector<Real>&                       a_dx,
                                     const RealVect&                           a_probLo,
                                     const int&                                a_ghost,
                                     const int                                 a_finestLevel)
{
  CH_TIME("EBAMRParticleMesh::EBAMRParticleMesh(full)");

  this->define(a_eblgs, a_refRat, a_dx, a_probLo, a_ghost, a_finestLevel);
}

EBAMRParticleMesh::~EBAMRParticleMesh()
{
  CH_TIME("EBAMRParticleMesh::~EBAMRParticleMesh()");
}

void
EBAMRParticleMesh::define(const Vector<RefCountedPtr<EBLevelGrid>>& a_eblgs,
                          const Vector<int>&                        a_refRat,
                          const Vector<Real>&                       a_dx,
                          const RealVect&                           a_probLo,
                          const int&                                a_ghost,
                          const int                                 a_finestLevel)
{
  CH_TIME("EBAMRParticleMesh::define");

  m_eblgs       = a_eblgs;
  m_refRat      = a_refRat;
  m_dx          = a_dx;
  m_probLo      = a_probLo;
  m_ghost       = a_ghost;
  m_finestLevel = a_finestLevel;
  m_verbose     = false;

  ParmParse pp("EBAMRParticleMesh");
  pp.query("verbose", m_verbose);

  this->defineLevelMotion();
  this->defineCoarseFineMotion();
  this->defineEBParticleMesh();
  this->defineOuterHaloMasks();
  this->defineTransitionMasks();

  m_isDefined = true;
}

void
EBAMRParticleMesh::defineLevelMotion()
{
  CH_TIME("EBAMRParticleMesh::defineLevelMotion");
  if (m_verbose) {
    pout() << "EBAMRParticleMesh::defineLevelMotion" << endl;
  }

  // TLDR: Define level Copiers. These are defined such that we can move data from valid+ghost -> valid. We need this because when we deposit particles
  //       we will also deposit into ghost cells that overlap with patches on the same level. The data in those ghost cells needs to find its way into
  //       the neighboring patches, which is what these Copiers do.
  m_levelCopiers.resize(1 + m_finestLevel);
  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const EBLevelGrid&       eblg   = *m_eblgs[lvl];
    const ProblemDomain&     domain = eblg.getDomain();
    const DisjointBoxLayout& dbl    = eblg.getDBL();

    // Define Copier. Note that the below code is basically the same as ghostDefine().
    const bool doExchange = true;

    // Define Copier as going from valid -> valid+ghost.
    m_levelCopiers[lvl].define(dbl, dbl, domain, m_ghost * IntVect::Unit, doExchange);

    // Define Copier as going from valid+ghost -> valid.
    m_levelCopiers[lvl].reverse();
  }
}

void
EBAMRParticleMesh::defineCoarseFineMotion()
{
  CH_TIME("EBAMRParticleMesh::defineCoarseFineMotion");
  if (m_verbose) {
    pout() << "EBAMRParticleMesh::defineCoarseFineMotion" << endl;
  }

  m_coarseFinePM.resize(1 + m_finestLevel);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {

    const bool hasCoar = (lvl > 0);

    if (hasCoar) {
      m_coarseFinePM[lvl] = RefCountedPtr<EBCoarseFineParticleMesh>(
        new EBCoarseFineParticleMesh(*m_eblgs[lvl - 1], *m_eblgs[lvl], m_refRat[lvl - 1], m_ghost * IntVect::Unit));
    }
    else {
      m_coarseFinePM[lvl] = RefCountedPtr<EBCoarseFineParticleMesh>(nullptr);
    }
  }
}

void
EBAMRParticleMesh::defineEBParticleMesh()
{
  CH_TIMERS("EBAMRParticleMesh::defineEBParticleMesh");
  CH_TIMER("EBAMRParticleMesh::defineEBParticleMesh::basic_defines", t1);
  CH_TIMER("EBAMRParticleMesh::defineEBParticleMesh::define_level", t2);
  CH_TIMER("EBAMRParticleMesh::defineEBParticleMesh::define_fico", t3);
  if (m_verbose) {
    pout() << "EBAMRParticleMesh::defineEBParticleMesh" << endl;
  }

  m_ebParticleMesh.resize(1 + m_finestLevel);
  m_ebParticleMeshFiCo.resize(1 + m_finestLevel);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    CH_START(t1);
    const ProblemDomain&     domain = m_eblgs[lvl]->getDomain();
    const DisjointBoxLayout& dbl    = m_eblgs[lvl]->getDBL();
    const DataIterator&      dit    = dbl.dataIterator();
    const EBISLayout&        ebisl  = m_eblgs[lvl]->getEBISL();

    const bool hasCoar = lvl > 0;

    m_ebParticleMesh[lvl] = RefCountedPtr<LayoutData<EBParticleMesh>>(new LayoutData<EBParticleMesh>(dbl));
    CH_START(t1);

    // Define the "regular" particle-mesh interpolation objects. These are defined on the input grids.
    CH_START(t2);
    const int nbox = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      const Box      cellBox = dbl[din];
      const EBISBox& ebisBox = ebisl[din];

      EBParticleMesh& particleMesh = (*m_ebParticleMesh[lvl])[din];

      particleMesh.define(domain, cellBox, ebisBox, m_dx[lvl] * RealVect::Unit, m_probLo);
    }
    CH_STOP(t2);

    // These are "special" particle-mesh interpolation objects for when we need to deposit coarse-level particles on a refined grid.
    CH_START(t3);
    if (hasCoar) {
      const EBLevelGrid&       eblgFiCo  = m_coarseFinePM[lvl]->getEblgFiCo();
      const ProblemDomain&     domain    = eblgFiCo.getDomain();
      const DisjointBoxLayout& dblFiCo   = eblgFiCo.getDBL();
      const DataIterator&      ditFiCo   = dblFiCo.dataIterator();
      const EBISLayout&        ebislFiCo = eblgFiCo.getEBISL();

      m_ebParticleMeshFiCo[lvl] = RefCountedPtr<LayoutData<EBParticleMesh>>(new LayoutData<EBParticleMesh>(dblFiCo));

      const int nboxFiCo = ditFiCo.size();
#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nboxFiCo; mybox++) {
        const DataIndex& din = ditFiCo[mybox];

        const Box      cellBox = dblFiCo[din];
        const EBISBox& ebisBox = ebislFiCo[din];

        EBParticleMesh& particleMesh = (*m_ebParticleMeshFiCo[lvl])[din];

        particleMesh.define(domain, cellBox, ebisBox, m_dx[lvl] * RealVect::Unit, m_probLo);
      }
    }
    else {
      m_ebParticleMeshFiCo[lvl] = RefCountedPtr<LayoutData<EBParticleMesh>>(nullptr);
    }
    CH_START(t3);
  }
}

const EBParticleMesh&
EBAMRParticleMesh::getEBParticleMesh(const int a_lvl, const DataIndex& a_dit) const
{
  CH_TIME("EBAMRParticleMesh::getEBParticleMesh");
  if (m_verbose) {
    pout() << "EBAMRParticleMesh::getEBParticleMesh" << endl;
  }

  CH_assert(a_lvl >= 0);
  CH_assert(a_lvl <= m_finestLevel);

  return (*m_ebParticleMesh[a_lvl])[a_dit];
}

void
EBAMRParticleMesh::defineOuterHaloMasks()
{
  CH_TIME("EBAMRParticleMesh::defineOuterHaloMasks");
  if (m_verbose) {
    pout() << "EBAMRParticleMesh::defineOuterHaloMasks" << endl;
  }

  constexpr int comp    = 0;
  constexpr int numComp = 1;

  m_outerHaloMasks.clear();

  for (int ighost = 1; ighost <= m_ghost; ighost++) {

    Vector<RefCountedPtr<LevelData<BaseFab<bool>>>> mask(1 + m_finestLevel);

    for (int lvl = 0; lvl < m_finestLevel; lvl++) {
      const DisjointBoxLayout& grids     = m_eblgs[lvl]->getDBL();
      const DisjointBoxLayout& gridsFine = m_eblgs[lvl + 1]->getDBL();

      const ProblemDomain& domain     = m_eblgs[lvl]->getDomain();
      const ProblemDomain& domainFine = m_eblgs[lvl + 1]->getDomain();

      // Create the coarsened fine grid.
      DisjointBoxLayout gridsCoFi;
      coarsen(gridsCoFi, gridsFine, m_refRat[lvl]);

      const DataIterator& dit     = grids.dataIterator();
      const DataIterator& ditFine = gridsFine.dataIterator();
      const DataIterator& ditCoFi = gridsCoFi.dataIterator();

      const int numBoxes     = dit.size();
      const int numBoxesFine = ditFine.size();
      const int numBoxesCoFi = ditCoFi.size();

      // Allocate data on this level, and reset the mask.
      mask[lvl] = RefCountedPtr<LevelData<BaseFab<bool>>>(new LevelData<BaseFab<bool>>(grids, numComp, IntVect::Zero));

      LevelData<BaseFab<bool>>& levelMask = *mask[lvl];

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < numBoxes; mybox++) {
        const DataIndex& din = dit[mybox];

        levelMask[din].setVal(false);
      }

      IntVectSet halo;

      // Go through the coarsened fine grid and set the halo to true
#pragma omp parallel for schedule(runtime) reduction(+ : halo)
      for (int mybox = 0; mybox < numBoxesCoFi; mybox++) {
        const DataIndex& din      = ditCoFi[mybox];
        const Box&       coFiBox  = gridsCoFi[din];
        const Box        grownBox = grow(coFiBox, ighost) & domain;

        // Subtract non-ghosted box and neighbor boxes
        IntVectSet myHalo(grownBox);
        myHalo -= coFiBox;

        NeighborIterator nit(gridsCoFi);
        for (nit.begin(din); nit.ok(); ++nit) {
          myHalo -= gridsCoFi[nit()];
        }

        halo |= myHalo;
      }

      // TLDR: In the above, we found the coarse-grid cells surrounding the fine level, viewed from the fine grids.
      //       Below, we create that view from the coarse grid. We use a LevelData<FArrayBox> on the coarsened fine grid,
      //       whose "ghost cells" can added to the _actual_ coarse grid. We then loop through those cells and set the
      //       mask directly.
      LevelData<FArrayBox> coFiMask(gridsCoFi, numComp, ighost * IntVect::Unit);
      LevelData<FArrayBox> coarMask(grids, numComp, IntVect::Zero);

#pragma omp parallel
      {
#pragma omp for schedule(runtime)
        for (int mybox = 0; mybox < numBoxesFine; mybox++) {
          coFiMask[ditFine[mybox]].setVal(0.0);
        }

#pragma omp for schedule(runtime)
        for (int mybox = 0; mybox < numBoxes; mybox++) {
          coarMask[dit[mybox]].setVal(0.0);
        }

        // Run through the halo and set the halo cells to 1 on the coarsened fine grids. Since dblCoFi was a coarsening of the fine
        // grid, the mask value in the valid region (i.e., not including ghosts) is always zero.
#pragma omp for schedule(runtime)
        for (int mybox = 0; mybox < numBoxesFine; mybox++) {
          const DataIndex& din = ditFine[mybox];

          const Box        region  = coFiMask[din].box();
          const IntVectSet curHalo = halo & region;

          for (IVSIterator ivsit(curHalo); ivsit.ok(); ++ivsit) {
            coFiMask[din](ivsit(), comp) = 1.0;
          }
        }
      }

      // Add the result to the coarse-grid data holder.
      Copier copier;
      copier.ghostDefine(gridsCoFi, grids, domain, ighost * IntVect::Unit);
      coFiMask.copyTo(Interval(comp, comp), coarMask, Interval(comp, comp), copier, LDaddOp<FArrayBox>());

      // Go through the grids and make the boolean mask. If there are no valid cells
      // we undefine the BaseFab.
#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < numBoxes; mybox++) {
        const DataIndex& din = dit[mybox];
        const Box        box = grids[din];

        BaseFab<bool>&   boolMask = levelMask[din];
        const FArrayBox& realMask = coarMask[din];

        bool emptyMask = true;

        auto kernel = [&](const IntVect& iv) -> void {
          if (realMask(iv, comp) > 0.0) {
            boolMask(iv, comp) = true;
            emptyMask          = false;
          }
        };

        BoxLoops::loop(box, kernel);

        // Undefine the BaseFab if the mask is empty. This means we can never do a copy.
        if (emptyMask) {
          boolMask.clear();
        }
      }
    }

    // Explicitly set mask on finest level to be nullptr, as there is no finer level.
    mask[m_finestLevel] = RefCountedPtr<LevelData<BaseFab<bool>>>(nullptr);

    m_outerHaloMasks.emplace(ighost, mask);
  }
}

void
EBAMRParticleMesh::defineTransitionMasks()
{
  CH_TIME("EBAMRParticleMesh::defineTransitionMasks");
  if (m_verbose) {
    pout() << "EBAMRParticleMesh::defineTransitionMasks" << endl;
  }

  constexpr int  comp    = 0;
  constexpr int  numComp = 1;
  const Interval interv  = Interval(comp, comp);

  m_transitionMasks.clear();

  for (int ighost = 1; ighost <= m_ghost; ighost++) {

    Vector<RefCountedPtr<LevelData<BaseFab<bool>>>> mask(1 + m_finestLevel);

    for (int lvl = 0; lvl < m_finestLevel; lvl++) {
      const EBLevelGrid& eblgCoar = *m_eblgs[lvl];
      const EBLevelGrid& eblgFine = *m_eblgs[lvl + 1];
      const EBLevelGrid& eblgFiCo = m_coarseFinePM[lvl + 1]->getEblgFiCo();

      const DisjointBoxLayout& gridsCoar = eblgCoar.getDBL();
      const DisjointBoxLayout& gridsFine = eblgFine.getDBL();
      const DisjointBoxLayout& gridsFiCo = eblgFiCo.getDBL();

      const ProblemDomain& domainCoar = eblgCoar.getDomain();
      const ProblemDomain& domainFine = eblgFine.getDomain();

      const DataIterator& ditCoar = gridsCoar.dataIterator();
      const DataIterator& ditFine = gridsFine.dataIterator();
      const DataIterator& ditFiCo = gridsFiCo.dataIterator();

      const int numBoxesCoar = ditCoar.size();
      const int numBoxesFine = ditFine.size();
      const int numBoxesFiCo = ditFiCo.size();

      // Allocate data on this level, and set the mask to false everywhere. We should not need ghost cells because
      // the mask is true only on the valid region.
      mask[lvl] = RefCountedPtr<LevelData<BaseFab<bool>>>(
        new LevelData<BaseFab<bool>>(gridsFiCo, numComp, IntVect::Zero));

      LevelData<BaseFab<bool>>& levelMask = *mask[lvl];
      LevelData<FArrayBox>      cfivsFine(gridsFine, numComp, ighost * IntVect::Unit);
      LevelData<FArrayBox>      cfivsFiCo(gridsFiCo, numComp, IntVect::Zero);

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < numBoxesFiCo; mybox++) {
        const DataIndex& din = ditFiCo[mybox];

        levelMask[din].setVal(false);
        cfivsFiCo[din].setVal(0.0);
      }

      // Build the CFIVS on the fine level.
#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < numBoxesFine; mybox++) {
        const DataIndex& din      = ditFine[mybox];
        const Box&       cellBox  = gridsFine[din];
        const Box&       ghostBox = grow(cellBox, ighost) & domainFine;

        cfivsFine[din].setVal(1.0, ghostBox, comp);
        cfivsFine[din].setVal(0.0, cellBox, comp);

        NeighborIterator nit(gridsFine);
        for (nit.begin(din); nit.ok(); ++nit) {
          const Box neighborBox = gridsFine[nit()];
          const Box mutualBox   = ghostBox & neighborBox;

          cfivsFine[din].setVal(0.0, mutualBox, comp);
        }
      }

      // Copy the result over to the refined coarse grids.
      Copier copier;
      copier.ghostDefine(gridsFine, gridsFiCo, domainFine, ighost * IntVect::Unit);
      cfivsFine.copyTo(interv, cfivsFiCo, interv, copier, LDaddOp<FArrayBox>());

      // Iterate through the coarse grid and set the mask.
#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < numBoxesFiCo; mybox++) {
        const DataIndex& din     = ditFiCo[mybox];
        const Box&       cellBox = gridsFiCo[din];

        BaseFab<bool>&   boolMask = levelMask[din];
        const FArrayBox& realMask = cfivsFiCo[din];

        bool emptyMask = true;

        auto kernel = [&](const IntVect& iv) -> void {
          if (realMask(iv, comp) > 0.5) {
            boolMask(iv, comp) = true;
            emptyMask          = false;
          }
        };

        BoxLoops::loop(cellBox, kernel);

        // Undefine the BaseFab if the mask is empty. This means we can never do a copy.
        if (emptyMask) {
          boolMask.clear();
        }
      }
    }

    // Explicitly set mask on finest level to be nullptr, as there is no finer level.
    mask[m_finestLevel] = RefCountedPtr<LevelData<BaseFab<bool>>>(nullptr);

    m_transitionMasks.emplace(ighost, mask);
  }
}

int
EBAMRParticleMesh::getTransitionMaskWidth(const DepositionType a_depositionType, const int a_refRat) const
{
  CH_TIME("EBAMRParticleMesh::getTransitionMaskWidth");
  if (m_verbose) {
    pout() << "EBAMRParticleMesh::getTransitionMaskWidth" << endl;
  }

  int maskWidth = 0;

  if (a_depositionType == DepositionType::CIC) {
    maskWidth = a_refRat / 2;
  }
  else if (a_depositionType == DepositionType::TSC) {
    maskWidth = a_refRat;
  }
  else {
    MayDay::Abort("EBAMRParticleMesh::getTransitionMaskWidth - logic bust");
  }

  if (maskWidth > m_ghost) {
    MayDay::Abort("EBAMRParticleMesh::transferMaskParticlesTransition -- not enough ghost cells!");
  }

  return maskWidth;
}

#include <CD_NamespaceFooter.H>
