#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// InfiniteMLCSolver.cpp
// petermc, 11 Jun 2004
// Compile with XTRACPPFLAGS=-DCH_USE_TIMER if you want to use the Timer.

#include <iostream>
#include <iomanip>
using std::ifstream;
using std::ios;

#ifdef CH_USE_TIMER
#ifdef HALEM_PROC_SPEED
#include <cstdio>
#include <sys/sysinfo.h>
#include <machine/hal_sysinfo.h>
#endif
#include "CH_Timer.H"

Timer TimeMLCdefine_all     ("MLCdefine all", 10);
Timer TimeInitBarrier       ("Initial Barrier",       TimeMLCdefine_all);
Timer TimeBegin             ("Begin",                 TimeMLCdefine_all);
Timer TimeBarrierPreLocal   ("Pre-local Barrier",     TimeMLCdefine_all);
Timer TimeLocalInit         ("Local Init",            TimeMLCdefine_all);
Timer TimeBarrierPostLocal  ("Post-local Barrier",    TimeMLCdefine_all);
Timer TimeBarrierPreGlobal  ("Pre-global Barrier",    TimeMLCdefine_all);
Timer TimeGlobalInit        ("Global Init",           TimeMLCdefine_all);
Timer TimeBarrierPostGlobal ("Post-global Barrier",   TimeMLCdefine_all);
Timer TimeFinalBarrier      ("Final Barrier",         TimeMLCdefine_all);

Timer TimeMLCsolve_all      ("MLCsolve all", 12);
Timer TimeBarrierPreInit    ("Pre-Barrier Initial",   TimeMLCsolve_all);
Timer TimeSetupArrays       ("Set up arrays",         TimeMLCsolve_all);
Timer TimeInitialLocal      ("Initial Local Solve",   TimeMLCsolve_all);
Timer TimeLocalResidual     ("Initial Local Residual",TimeMLCsolve_all);
Timer TimeBarrierPreResid   ("Pre-Barrier Resid",     TimeMLCsolve_all);
Timer TimeCopyResidual      ("Copy Residual",         TimeMLCsolve_all);
Timer TimeBarrierPostResid  ("Post-Barrier Resid",    TimeMLCsolve_all);
Timer TimeReduce            ("Reduce",                TimeMLCsolve_all);
Timer TimeGlobalSolve       ("Global Solve",          TimeMLCsolve_all);
Timer TimeBarrierPreFine    ("Pre-Barrier Fine",      TimeMLCsolve_all);
Timer TimeCopyOverlapFine   ("Copy Fine Overlap",     TimeMLCsolve_all);
Timer TimeBarrierPostFine   ("Post-Barrier Fine",     TimeMLCsolve_all);
Timer TimeBarrierPreCrse    ("Pre-Barrier Crse",      TimeMLCsolve_all);
Timer TimeCopyOverlapCrse   ("Copy Coarse Overlap",   TimeMLCsolve_all);
Timer TimeBarrierPostCrse   ("Post-Barrier Crse",     TimeMLCsolve_all);
Timer TimeBarrierPreCopy    ("Pre-Barrier Copy",      TimeMLCsolve_all);
Timer TimeCopyRhs           ("Copy rhs",              TimeMLCsolve_all);
Timer TimeBarrierPostCopy   ("Post-Barrier Copy",     TimeMLCsolve_all);
Timer TimeFinalLocal        ("Final Local Solve",     TimeMLCsolve_all);
#endif

#include "parstream.H"
#include "InfiniteMLCSolver.H"
#include "InfiniteDomainF.H"
#include "MayDay.H"
#include "LayoutIterator.H"
#include "EvalOperatorF_F.H"
//#include "AMRIO.H"
using std::cout;
using std::endl;
using std::cerr;
// not used here, but needed for library
#include "HalveF_F.H"
#include "PrecondF_F.H"
#include "TotalChargeF_F.H"
#include "GradientsF_F.H"

#include "NamespaceHeader.H"

// ---------------------------------------------------------
// default constructor
InfiniteMLCSolver::InfiniteMLCSolver()
{
  setDefaultValues();
}


// ---------------------------------------------------------
void
InfiniteMLCSolver::define(const Box&                 a_fineBaseBox,
                          const Box&                 a_fineDomain,
                          int                        a_refToCoarse,
                          const RealVect&            a_dx,
                          int                        a_s1Local,
                          int                        a_s2Local,
                          int                        a_patchSizeLocal,
                          int                        a_dstFaceCoarseningLocal,
                          PoissonDirichlet::OperatorType a_opLocalInitial,
                          int                        a_multipoleOrderLocal,
                          int                        a_degreeNormalDerivativeLocal,
                          int                        a_interpBorderLocal,
                          int                        a_s1Global,
                          int                        a_s2Global,
                          int                        a_patchSizeGlobal,
                          int                        a_dstFaceCoarseningGlobal,
                          PoissonDirichlet::OperatorType a_opGlobal,
                          int                        a_multipoleOrderGlobal,
                          int                        a_degreeNormalDerivativeGlobal,
                          int                        a_interpBorderGlobal,
                          int                        a_interpBorderCF,
                          int                        a_minBufferCoarse,
                          PoissonDirichlet::OperatorType a_opLocalFinal)
{
#ifdef CH_USE_TIMER
  TimeMLCdefine_all.start();
  TimeInitBarrier.start();
  barrier();
  TimeInitBarrier.stop();
  TimeBegin.start();
#endif

  clearMemory();
  m_verbose = 0;
  m_fineBaseBox = a_fineBaseBox - a_fineBaseBox.smallEnd();
  m_fineDomain = a_fineDomain;
  m_refToCoarse = a_refToCoarse;
  m_dx = a_dx;
  m_s1Local = a_s1Local;
  m_s2Local = a_s2Local;
  m_patchSizeLocal = a_patchSizeLocal;
  m_dstFaceCoarseningLocal = a_dstFaceCoarseningLocal;
  m_opLocalInitial = a_opLocalInitial;
  m_multipoleOrderLocal = a_multipoleOrderLocal;
  m_degreeNormalDerivativeLocal = a_degreeNormalDerivativeLocal;
  m_interpBorderLocal = a_interpBorderLocal;
  m_s1Global = a_s1Global;
  m_s2Global = a_s2Global;
  m_patchSizeGlobal = a_patchSizeGlobal;
  m_dstFaceCoarseningGlobal = a_dstFaceCoarseningGlobal;
  m_opGlobal = a_opGlobal;
  m_multipoleOrderGlobal = a_multipoleOrderGlobal;
  m_degreeNormalDerivativeGlobal = a_degreeNormalDerivativeGlobal;
  m_interpBorderGlobal = a_interpBorderGlobal;
  m_interpLayer = a_interpBorderCF;
  m_minBufferCoarse = a_minBufferCoarse;
  m_opLocalFinal = a_opLocalFinal;

  // m_bufferLocal = m_s1Local + m_s2Local;  // FIX:  this can be a bit too high
  m_bufferLocal = m_minBufferCoarse * m_refToCoarse;
  m_bufferLocalSoln = m_s1Local + m_s2Local;

  IntVect baseSize = m_fineBaseBox.size();

  /*
    Define initial local infinite-domain solver for fine level.
  */

  // We use the coarse solution for interpolating the
  // Dirichlet boundary condition in the final local solve.
  m_coarseAddRadiusLocal = m_interpLayer + m_minBufferCoarse;

#ifdef CH_USE_TIMER
  TimeBegin.stop();
  TimeBarrierPreLocal.start();
  barrier();
  TimeBarrierPreLocal.stop();
  TimeLocalInit.start();
#endif

  // This solver will give us initial local solution on
  // grow(fineBoxBaseNodes, m_bufferLocalSoln).
  Box fineBoxBaseNodes = surroundingNodes(m_fineBaseBox);
  m_solverLocalInitial.define(fineBoxBaseNodes,
                              m_s1Local, m_s2Local,
                              m_patchSizeLocal, m_patchSizeLocal,
                              m_dx,
                              m_opLocalInitial,
                              m_multipoleOrderLocal,
                              m_degreeNormalDerivativeLocal,
                              m_interpBorderLocal,
                              false,  // not parallel
                              true,  // yes, get coarse
                              m_refToCoarse,
                              m_coarseAddRadiusLocal
                              );
  m_solverLocalInitial.setVerbose(m_verbose);

  /*
    Define final local Dirichlet solver for fine level.
  */

  m_poisson.define(fineBoxBaseNodes, m_dx, m_opLocalFinal);

#ifdef CH_USE_TIMER
  TimeLocalInit.stop();
  TimeBarrierPostLocal.start();
  barrier();
  TimeBarrierPostLocal.stop();
#endif

  /*
    Coefficients of interpolation from global solution to local.
  */

  int ncoeffs = (m_refToCoarse - 1) * (2 * (m_interpLayer + 1));
  m_interpCoeffs = new Real[ncoeffs];
  FORT_GETPOLYINTERPC(m_interpCoeffs, &m_refToCoarse, &m_interpLayer);

  /*
    Now define solver for global coarse level.

    residGlobalCoarse lives on NODEs of m_coarseGlobalDomain.
  */

  m_dxCoarse = Real(m_refToCoarse) * m_dx;

  // NEW, 25 Jan 2005
  // int bufferCoarseForGlobal = 2 * ceil((m_minBufferCoarse - 1) / 2.);
  int bufferCoarseForGlobal = m_minBufferCoarse - 1;
  if (bufferCoarseForGlobal % 2 != 0) bufferCoarseForGlobal++;
  m_coarseGlobalDomain = grow(coarsen(m_fineDomain, m_refToCoarse),
                              bufferCoarseForGlobal);

  int coarseLength = m_coarseGlobalDomain.longside();
  // lengths of m_coarseGlobalDomain are powers of 2

  // If coarseLength < 8 then get0buffersize complains.
  CH_assert(coarseLength >= 8);

  m_coarseGlobalDomainPhi = grow(m_coarseGlobalDomain, m_s1Global + m_s2Global);

#ifdef TIMER
  TimeBarrierPreGlobal.start();
  barrier();
  TimeBarrierPreGlobal.stop();
  TimeGlobalInit.start();
#endif

  Box coarseGlobalNodes = surroundingNodes(m_coarseGlobalDomain);
  m_solverGlobal.define(coarseGlobalNodes, m_s1Global, m_s2Global,
                        m_patchSizeGlobal, m_patchSizeGlobal,
                        m_dxCoarse,
                        m_opGlobal,
                        m_multipoleOrderGlobal,
                        m_degreeNormalDerivativeGlobal,
                        m_interpBorderGlobal,
                        true // parallel
                        );
  m_solverGlobal.setVerbose(m_verbose);

  m_isDefined = true;

#ifdef CH_USE_TIMER
  TimeGlobalInit.stop();
  TimeBarrierPostGlobal.start();
  barrier();
  TimeBarrierPostGlobal.stop();
#ifdef MEMORY_USAGE
  Real end_memory = get_memory_usage_from_OS();
#ifdef CH_MPI
  // This is for finding avg/min/max of end_memory over all procs.
  // Writes to pout.0.
  Real avg_memory, min_memory, max_memory;
  gather_memory_from_procs(end_memory, avg_memory, min_memory, max_memory);
#endif
#endif
  TimeFinalBarrier.start();
  barrier();
  TimeFinalBarrier.stop();

  TimeMLCdefine_all.stop();

  Real pctTime = 100. / TimeMLCdefine_all.wc_time();

  pout() << "InfiniteMLCSolver define done. "
#ifdef MEMORY_USAGE
         << " mem= " << end_memory
#endif
         << "MB wall-clock time: "    << setiosflags(ios::fixed)
         << TimeMLCdefine_all.wc_time()      << " sec" << endl;
  pout() << "Initial barrier:      " << TimeInitBarrier.wc_time()
         << " seconds"
         << " (" << (TimeInitBarrier.wc_time() * pctTime)
         << "%)" << endl;
  pout() << "Define begin          " << TimeBegin.wc_time()
         << " seconds"
         << " (" << (TimeBegin.wc_time() * pctTime)
         << "%)" << endl;
  pout() << "Pre-local Barrier:    " << TimeBarrierPreLocal.wc_time()
         << " seconds"
         << " (" << (TimeBarrierPreLocal.wc_time() * pctTime)
         << "%)" << endl;
  pout() << "Local Solver Init:    " << TimeLocalInit.wc_time()
         << " seconds"
         << " (" << (TimeLocalInit.wc_time() * pctTime)
         << "%)" << endl;
  pout() << "Post-local Barrier:   " << TimeBarrierPostLocal.wc_time()
         << " seconds"
         << " (" << (TimeBarrierPostLocal.wc_time() * pctTime)
         << "%)" << endl;
  pout() << "Pre-global Barrier:   " << TimeBarrierPreGlobal.wc_time()
         << " seconds"
         << " (" << (TimeBarrierPreGlobal.wc_time() * pctTime)
         << "%)" << endl;
  pout() << "Global Solver Init:   " << TimeGlobalInit.wc_time()
         << " seconds"
         << " (" << (TimeGlobalInit.wc_time() * pctTime)
         << "%)" << endl;
  pout() << "Post-global Barrier:  " << TimeBarrierPostGlobal.wc_time()
         << " seconds"
         << " (" << (TimeBarrierPostGlobal.wc_time() * pctTime)
         << "%)" << endl;
  pout() << "Final define barrier: " << TimeFinalBarrier.wc_time()
         << " seconds"
         << " (" << (TimeFinalBarrier.wc_time() * pctTime)
         << "%)" << endl;
#endif

  // Timer::TimerSummary(); // causes segfault
}


// ---------------------------------------------------------
// abbreviated define function
void
InfiniteMLCSolver::define(
                          const Box&                 a_fineBaseBox,
                          const Box&                 a_fineDomain,
                          int                        a_refToCoarse,
                          const RealVect&            a_dx)
{
  IntVect domainSize = a_fineDomain.size();

  Box fineBoxBase = a_fineBaseBox - a_fineBaseBox.smallEnd();
  IntVect baseSize = fineBoxBase.size();

  // global solve should be O(H^4) accurate, with H = global mesh spacing
  int accuracyGlobal = 4;
  int multipoleOrderGlobal = 2*accuracyGlobal - 1;
  int degreeNormalDerivativeGlobal = accuracyGlobal;
  int interpBorderGlobal = accuracyGlobal - 1;

  PoissonDirichlet::OperatorType opGlobal = PoissonDirichlet::POINT27;

  PoissonDirichlet::OperatorType opLocalFinal = PoissonDirichlet::POINT7;

  int accuracyCF = 4;
  int interpBorderCF = accuracyCF/2 - 1;

  int minBufferCoarse = 2;
  int minBufferLocal = minBufferCoarse * a_refToCoarse;

  if (minBufferLocal > baseSize[0] ||
      minBufferLocal > baseSize[1] ||
      minBufferLocal > baseSize[2])
    {
      cerr << "buffer width " << minBufferLocal
           << " exceeds patch size " << baseSize
           << endl;
      MayDay::Error("returning");
      return;
    }

  int bufferLocalSoln;
  int patchSizeLocal;
  int s1Local, s2Local;
  bool alphaZeroLocal = true;
  if (alphaZeroLocal)
    {
      int errcode = 0;
      FORT_GET0BUFFERSIZE3D(baseSize.dataPtr(), &minBufferLocal,
                            &bufferLocalSoln, &patchSizeLocal, &errcode);
      if (errcode != 0)
        {
          cerr << "get0buffersize3d returned error code "
               << errcode << endl;
          MayDay::Error("returning");
          return;
        }
      s1Local = 0;
      s2Local = bufferLocalSoln;
    }
  else
    {
      int errcode = 0;
      FORT_GETBUFFERSIZE3D(baseSize.dataPtr(), &minBufferLocal,
                           &bufferLocalSoln, &patchSizeLocal, &errcode);
      if (errcode != 0)
        {
          cerr << "getbuffersize3d returned error code "
               << errcode << endl;
          MayDay::Error("returning");
          return;
        }
      s1Local = bufferLocalSoln / 2;
      s2Local = bufferLocalSoln / 2;
    }

  PoissonDirichlet::OperatorType opLocalInitial = PoissonDirichlet::POINT27;

  // local solve should be O(h^4) accurate, with h = local mesh spacing
  int accuracyLocal = 4;
  int multipoleOrderLocal = 2*accuracyLocal - 1;
  int degreeNormalDerivativeLocal = accuracyLocal;
  int interpBorderLocal = accuracyLocal - 1;

  // We use the coarse solution for interpolating the
  // Dirichlet boundary condition in the final local solve.
  int coarseAddRadiusLocal = interpBorderCF + minBufferCoarse;

//   if (verbose > 0)
//     cout << "from nbox = " << nbox
//          << " get minBufferLocal = " << minBufferLocal
//          << " and use bufferLocalSoln = " << bufferLocalSoln
//          << endl;
  // Set coarseLength to length of coarse global domain.
  // The amount by which to grow is set to
  // 2*ceil((minBufferCoarse-1)/2) because it must be at least
  // minBufferCoarse-1 in order to get all of the residual,
  // and it must be even so that the domain length is divisible by 4
  // because in the single-grid solver, patch lengths must be
  // divisible by 4.
  // int bufferCoarseForGlobal = 2 * ceil((bufferCoarse - 1) / 2.);
  int bufferCoarseForGlobal = minBufferCoarse - 1;
  if (bufferCoarseForGlobal % 2 != 0) bufferCoarseForGlobal++;
  // int coarseLength = ndom/a_refToCoarse + 2*bufferCoarseForGlobal;
  IntVect coarse3d =
    domainSize / a_refToCoarse +
    (2*bufferCoarseForGlobal) * IntVect::Unit;

  // If coarseLength < 8 then get0buffersize complains.
  CH_assert(coarse3d >= 8 * IntVect::Unit);

  // Now we find the buffer size bufferGlobal to use for global solution:
  // we want bufferGlobal >= coarseAddRadiusLocal.
  int bufferGlobal;
  int s1Global, s2Global;
  int patchSizeGlobal;
  bool alphaZeroGlobal = true;
  if (alphaZeroGlobal)
    {
      int errcode = 0;
      FORT_GET0BUFFERSIZE3D(coarse3d.dataPtr(), &coarseAddRadiusLocal,
                            &bufferGlobal, &patchSizeGlobal, &errcode);
      if (errcode != 0)
        {
          cerr << "get0buffersize3d returned error code " << errcode << endl;
          MayDay::Error("returning");
          return;
        }
      s1Global = 0;
      s2Global = bufferGlobal;
    }
  else
    {
      int errcode = 0;
      FORT_GETBUFFERSIZE3D(coarse3d.dataPtr(), &coarseAddRadiusLocal,
                           &bufferGlobal, &patchSizeGlobal, &errcode);
      if (errcode != 0)
        {
          cerr << "getbuffersize3d returned error code " << errcode << endl;
          MayDay::Error("returning");
          return;
        }
      s1Global = bufferGlobal / 2;
      s2Global = bufferGlobal / 2;
    }
//   if (verbose > 0)
//     cout << "from coarseLength = " << coarseLength
//          << " we get bufferGlobal = " << bufferGlobal
//          << endl;

  define(a_fineBaseBox, a_fineDomain, a_refToCoarse, a_dx,
         s1Local, s2Local, patchSizeLocal, patchSizeLocal,
         opLocalInitial, multipoleOrderLocal,
         degreeNormalDerivativeLocal, interpBorderLocal,
         s1Global, s2Global, patchSizeGlobal, patchSizeGlobal,
         opGlobal, multipoleOrderGlobal,
         degreeNormalDerivativeGlobal, interpBorderGlobal,
         interpBorderCF, minBufferCoarse,
         opLocalFinal);
}


// ---------------------------------------------------------
// complete constructor
InfiniteMLCSolver::
InfiniteMLCSolver(
                  const Box&                 a_fineBaseBox,
                  const Box&                 a_fineDomain,
                  int                        a_refToCoarse,
                  const RealVect&            a_dx,
                  int                        a_s1Local,
                  int                        a_s2Local,
                  int                        a_patchSizeLocal,
                  int                        a_dstFaceCoarseningLocal,
                  PoissonDirichlet::OperatorType a_opLocalInitial,
                  int                        a_multipoleOrderLocal,
                  int                        a_degreeNormalDerivativeLocal,
                  int                        a_interpBorderLocal,
                  int                        a_s1Global,
                  int                        a_s2Global,
                  int                        a_patchSizeGlobal,
                  int                        a_dstFaceCoarseningGlobal,
                  PoissonDirichlet::OperatorType a_opGlobal,
                  int                        a_multipoleOrderGlobal,
                  int                        a_degreeNormalDerivativeGlobal,
                  int                        a_interpBorderGlobal,
                  int                        a_interpBorderCF,
                  int                        a_minBufferCoarse,
                  PoissonDirichlet::OperatorType a_opLocalFinal)
{
  setDefaultValues();
  define(a_fineBaseBox, a_fineDomain, a_refToCoarse, a_dx,
         a_s1Local, a_s2Local, a_patchSizeLocal, a_dstFaceCoarseningLocal,
         a_opLocalInitial, a_multipoleOrderLocal,
         a_degreeNormalDerivativeLocal, a_interpBorderLocal,
         a_s1Global, a_s2Global, a_patchSizeGlobal, a_dstFaceCoarseningGlobal,
         a_opGlobal, a_multipoleOrderGlobal,
         a_degreeNormalDerivativeGlobal, a_interpBorderGlobal,
         a_interpBorderCF, a_minBufferCoarse, a_opLocalFinal);
}


// ---------------------------------------------------------
// abbreviated constructor
InfiniteMLCSolver::
InfiniteMLCSolver(
                  const Box&                 a_fineBaseBox,
                  const Box&                 a_fineDomain,
                  int                        a_refToCoarse,
                  const RealVect&            a_dx)
{
  setDefaultValues();
  define(a_fineBaseBox, a_fineDomain, a_refToCoarse, a_dx);
}


// ---------------------------------------------------------
// destructor
InfiniteMLCSolver::~InfiniteMLCSolver()
{
  clearMemory();
}


// ---------------------------------------------------------
bool
InfiniteMLCSolver::isDefined() const
{
  return m_isDefined;
}


// ---------------------------------------------------------
void
InfiniteMLCSolver::setDefaultValues()
{
  m_isDefined = false;
  m_verbose = 0;
}


// ---------------------------------------------------------
void
InfiniteMLCSolver::setVerbose(int a_verbose)
{
  m_verbose = a_verbose;
}


// ---------------------------------------------------------
// returns InfiniteMLCSolver to a basically undefined state
void
InfiniteMLCSolver::clearMemory()
{
  if (m_isDefined)
    {
      delete m_interpCoeffs;
    }
}


// ---------------------------------------------------------
// solve
void
InfiniteMLCSolver::solve(LevelData<NodeFArrayBox>&    a_phi,
                         LevelData<NodeFArrayBox>&    a_rhs)
{
#ifdef CH_USE_TIMER
  TimeMLCsolve_all.start();
  TimeBarrierPreInit.start();
  barrier();
  TimeBarrierPreInit.stop();
  TimeSetupArrays.start();
#endif

  if (m_verbose > 0)
    pout() << "in InfiniteMLCSolver::solve" << endl;

  CH_assert(isDefined());

  const Interval intvl0(0, 0);

  const DisjointBoxLayout& phiGrids = a_phi.getBoxes();
  const DisjointBoxLayout& rhsGrids = a_rhs.getBoxes();

#ifndef NDEBUG
  // Check that all boxes:
  // - are contained in m_fineDomain;
  // - have the same dimensions as m_fineBaseBox;
  // - are all translations from the same point by multiples of m_fineBaseBox.
  IntVect baseSize = m_fineBaseBox.size();
  LayoutIterator litFirst = rhsGrids.layoutIterator();
  const Box& rhsFirstBox = rhsGrids.get(litFirst());
  IntVect rhsFirstOrigin = rhsFirstBox.smallEnd();
  for (LayoutIterator lit = phiGrids.layoutIterator(); lit.ok(); ++lit)
    {
      const Box& bx = phiGrids.get(lit());
      CH_assert(m_fineDomain.contains(bx));
      CH_assert(bx.size() == baseSize);
      IntVect bxOff = bx.smallEnd() - rhsFirstOrigin;
      CH_assert(baseSize * (bxOff / baseSize) == bxOff);
    }
  for (LayoutIterator lit = rhsGrids.layoutIterator(); lit.ok(); ++lit)
    {
      const Box& bx = rhsGrids.get(lit());
      CH_assert(m_fineDomain.contains(bx));
      CH_assert(bx.size() == baseSize);
      IntVect bxOff = bx.smallEnd() - rhsFirstOrigin;
      CH_assert(baseSize * (bxOff / baseSize) == bxOff);
    }
#endif

  /*
    Coarsened and/or expanded versions of rhsGrids
  */

  DisjointBoxLayout rhsCoarseLocalGrids;

  coarsen(rhsCoarseLocalGrids, rhsGrids, m_refToCoarse);

  BoxLayout rhsGridsBuffer;
  rhsGridsBuffer.deepCopy(rhsGrids);
  rhsGridsBuffer.grow(m_bufferLocal);
  rhsGridsBuffer.close();

  BoxLayout rhsCoarseLocalGridsBufferm1;
  rhsCoarseLocalGridsBufferm1.deepCopy(rhsCoarseLocalGrids);
  rhsCoarseLocalGridsBufferm1.grow(m_minBufferCoarse-1);
  rhsCoarseLocalGridsBufferm1.close();

  BoxLayout rhsCoarseLocalGridsBufferInterp;
  rhsCoarseLocalGridsBufferInterp.deepCopy(rhsCoarseLocalGrids);
  rhsCoarseLocalGridsBufferInterp.grow(m_coarseAddRadiusLocal);
  rhsCoarseLocalGridsBufferInterp.close();

  /*
    Coarsened and/or expanded versions of phiGrids
  */

  DisjointBoxLayout phiCoarseLocalGrids;
  coarsen(phiCoarseLocalGrids, phiGrids, m_refToCoarse);

  BoxLayout phiCoarseLocalGridsInterp;
  phiCoarseLocalGridsInterp.deepCopy(phiCoarseLocalGrids);
  phiCoarseLocalGridsInterp.grow(m_interpLayer);
  phiCoarseLocalGridsInterp.close();

  // We use phiLocalCoarse for two purposes:
  // (1) evaluating the Mehrstellen operator on phiLocalCoarse,
  //     summing up the results over all boxes, and using the sum
  //     as right-hand side for a global solve;
  // (2) interpolating from phiLocalCoarse to Dirichlet
  //     boundary conditions for the final local solve.
  //     (Here we take its redistributed form, intCoarseData.)

  // In (2), we interpolate
  // to:  fine points in buffer layer, width bufferLocal
  // from:  coarse points in coarse buffer layer,
  //        width m_minBufferCoarse + interpLayer
  //
  // Recall rhsCoarseLocalGridsBufferInterp has buffer size
  // m_coarseAddRadiusLocal = m_interpLayer + m_minBufferCoarse;
  // m_minBufferCoarse = bufferLocal / refToCoarse.
  // We need the extra m_interpLayer in order to do interpolation.
  // LevelData<NodeFArrayBox> phiLocalCoarse(coarseLocalGrids, 1,
  // m_coarseAddRadiusLocal * IntVect::Unit);

  // Define phiLocalCoarse on
  // grow(refine(rhsGrids, m_refToCoarse), m_coarseAddRadiusLocal),
  // where m_coarseAddRadiusLocal = m_interpLayer + m_minBufferCoarse.
  BoxLayoutData<NodeFArrayBox> phiLocalCoarse(rhsCoarseLocalGridsBufferInterp, 1);
  unsigned long int sizeLayout = 0;
  for (DataIterator dit = phiLocalCoarse.dataIterator(); dit.ok(); ++dit)
    {
      const Box& bx = rhsCoarseLocalGridsBufferInterp.get(dit());
      sizeLayout += phiLocalCoarse[dit()].size(bx, intvl0);
    }
  pout() << "InfiniteMLCSolver:  phiLocalCoarse "
         << (sizeLayout / sizeof(Real)) << " scalars, "
         << (Real(sizeLayout) / 1048576.) << " MB" << endl;

  // residLocalCoarse is obtained from evaluating the
  // Mehrstellen operator on phiLocalCoarse with
  // buffer layer of width m_minBufferCoarse.
  // So residLocalCoarse has buffer layer of width m_minBufferCoarse - 1.
  //
  // Since phiLocalCoarse has a buffer layer of width
  // coarseAddRadiusLocal, we could evaluate the operator up
  // to buffer layer of width coarseAddRadiusLocal - 1, but
  // for the purpose of finding the global right-hand side,
  // we evaluate it only up to width m_minBufferCoarse - 1.

  // LevelData<NodeFArrayBox> residLocalCoarse(coarseLocalGrids, 1,
  // (m_minBufferCoarse-1) * IntVect::Unit);

  // Define residLocalCoarse on
  // grow(refine(rhsGrids, m_refToCoarse), m_minBufferCoarse-1).
  BoxLayoutData<NodeFArrayBox> residLocalCoarse(rhsCoarseLocalGridsBufferm1, 1);
  sizeLayout = 0;
  for (DataIterator dit = residLocalCoarse.dataIterator(); dit.ok(); ++dit)
    {
      const Box& bx = rhsCoarseLocalGridsBufferm1.get(dit());
      sizeLayout += residLocalCoarse[dit()].size(bx, intvl0);
    }
  pout() << "InfiniteMLCSolver:  residLocalCoarse "
         << (sizeLayout / sizeof(Real)) << " scalars, "
         << (Real(sizeLayout) / 1048576.) << " MB" << endl;

  ProblemDomain pdFine(m_fineDomain);
  if (m_verbose > 0)
    pout() << "pdFine == " << pdFine << endl;

  // From grow(rhsGrids, m_bufferLocal) = rhsGridsBuffer to phiGrids.
  LayoutData< Vector<RefCountedPtr<NodeFArrayBox> > > intFineData;

#ifdef CH_USE_TIMER
  TimeSetupArrays.stop();
  TimeInitialLocal.start();
  m_solverLocalInitial.beginTimers();
#endif

  { // scope for phiLocalFine, which is a big data holder

    // Define phiLocalFine on
    // grow(rhsGrids, m_bufferLocal),
    // where m_bufferLocal = m_refToCoarse * m_minBufferCoarse.
    BoxLayoutData<NodeFArrayBox> phiLocalFine(rhsGridsBuffer, 1);
    for (DataIterator dit = phiLocalFine.dataIterator(); dit.ok(); ++dit)
      {
        const Box& bx = rhsGridsBuffer.get(dit());
        sizeLayout += phiLocalFine[dit()].size(bx, intvl0);
      }
    pout() << "InfiniteMLCSolver:  phiLocalFine "
           << (sizeLayout / sizeof(Real)) << " scalars, "
           << (Real(sizeLayout) / 1048576.) << " MB" << endl;

    /*****************************************************************/
    // Initial local solve for phiLocalFine and phiLocalCoarse,
    // and evaluation residLocalCoarse = operator(phiLocalCoarse).
    /*****************************************************************/

    for (DataIterator dit = rhsGrids.dataIterator(); dit.ok(); ++dit)
      {
        NodeFArrayBox& rhsNFAB = a_rhs[dit()];
        NodeFArrayBox& phiLocalFineNFAB = phiLocalFine[dit()];
        NodeFArrayBox& phiLocalCoarseNFAB = phiLocalCoarse[dit()];

        // phiLocalFineFab is defined on NODEs of a box of
        // grow(rhsGrids, m_bufferLocal)
        // where m_bufferLocal = m_refToCoarse * m_minBufferCoarse.
        FArrayBox& phiLocalFineFab = phiLocalFineNFAB.getFab();

        // phiLocalCoarseFab is defined on NODEs of a box of
        // grow(refine(rhsGrids, m_refToCoarse), m_coarseAddRadiusLocal),
        // where m_coarseAddRadiusLocal = m_interpLayer + m_minBufferCoarse.
        FArrayBox& phiLocalCoarseFab = phiLocalCoarseNFAB.getFab();

        // rhsFab is defined on NODEs of a box of rhsGrids.
        FArrayBox& rhsFab = rhsNFAB.getFab();

        // Do solve.
        // Store fine solution in phiLocalFineFab,
        // coarse solution in phiLocalCoarseFab.

        if (m_verbose > 0)
          pout() << "calling solveWithCoarse on local" << endl;

        IntVect shiftForSolver = rhsFab.smallEnd();
        rhsFab.shift(-shiftForSolver);
        IntVect shiftCoarseForSolver = shiftForSolver / m_refToCoarse;
        phiLocalCoarseFab.shift(-shiftCoarseForSolver);

        if (m_bufferLocalSoln == m_bufferLocal)
          {
            phiLocalFineFab.shift(-shiftForSolver);
            m_solverLocalInitial.solveWithCoarse(phiLocalFineFab,
                                                 phiLocalCoarseFab,
                                                 rhsFab);
            phiLocalFineFab.shift(shiftForSolver);
          }
        else // m_bufferLocalSoln > m_bufferLocal
          {
            Box bxExt = grow(phiLocalFineFab.box(),
                             m_bufferLocalSoln - m_bufferLocal);

            // phiLocalFineFabExt is defined on NODEs of a box of
            // grow(rhsGrids, m_bufferLocalSoln)
            // where m_bufferLocalSoln > m_bufferLocal.
            FArrayBox phiLocalFineFabExt(bxExt, 1);

            phiLocalFineFabExt.shift(-shiftForSolver);
            m_solverLocalInitial.solveWithCoarse(phiLocalFineFabExt,
                                                 phiLocalCoarseFab,
                                                 rhsFab);
            phiLocalFineFabExt.shift(shiftForSolver);
            phiLocalFineFab.copy(phiLocalFineFabExt);
          }
        rhsFab.shift(shiftForSolver);
        phiLocalCoarseFab.shift(shiftCoarseForSolver);

#ifndef NDEBUG
        // dummy statement in order to get around gdb bug
        int dummy_unused = 0; dummy_unused = 0;
#endif
      }

#ifdef CH_USE_TIMER
    m_solverLocalInitial.endTimers();
    TimeInitialLocal.stop();
#endif

    if (m_verbose > 0)
      pout() << "getting generalCopyTo for fine" << endl;

#ifdef CH_USE_TIMER
    TimeBarrierPreFine.start();
    barrier();
    TimeBarrierPreFine.stop();
    TimeCopyOverlapFine.start();
#endif

    phiLocalFine.generalCopyTo(phiGrids, // destination layout
                               intFineData, // destination
                               intvl0, // copy component 0
                               pdFine); // domain of source and destination

    // After this, we no longer need phiLocalFine.
    // Everything we'll use is in intFineData.

    sizeLayout = 0;
    for (DataIterator dit = phiGrids.dataIterator(); dit.ok(); ++dit)
      {
        Vector< RefCountedPtr<NodeFArrayBox> > intFineDataPatch =
          intFineData[dit()];
        for (int indFine = 0; indFine < intFineDataPatch.size(); indFine++)
          {
            const NodeFArrayBox& intFineNFAB = *intFineDataPatch[indFine];
            sizeLayout += intFineNFAB.size(intFineNFAB.box(), intvl0);
          }
      }
    pout() << "InfiniteMLCSolver:  intFineData "
         << (sizeLayout / sizeof(Real)) << " scalars, "
         << (Real(sizeLayout) / 1048576.) << " MB" << endl;

#ifdef CH_USE_TIMER
    TimeCopyOverlapFine.stop();
    TimeBarrierPostFine.start();
#endif
    if (m_verbose > 0)
      pout() << "got generalCopyTo for fine" << endl;
  } // end of scope for phiLocalFine
  barrier();
#ifdef CH_USE_TIMER
  TimeBarrierPostFine.stop();
  TimeLocalResidual.start();
#endif

  for (DataIterator dit = rhsGrids.dataIterator(); dit.ok(); ++dit)
    {
      if (m_verbose > 0)
        pout() << "evaluating local residual" << endl;

      // Evaluate operator on coarse local solution phiLocalCoarse,
      // and store in residLocalCoarse.

      NodeFArrayBox& phiLocalCoarseNFAB = phiLocalCoarse[dit()];

      // phiLocalCoarseFab is defined on NODEs of a box of
      // grow(refine(rhsGrids, m_refToCoarse), m_coarseAddRadiusLocal),
      // where m_coarseAddRadiusLocal = m_interpLayer + m_minBufferCoarse.
      FArrayBox& phiLocalCoarseFab = phiLocalCoarseNFAB.getFab();

      // residLocalCoarseFab is defined on NODEs of a box of
      // grow(refine(rhsGrids, m_refToCoarse), m_minBufferCoarse-1).
      FArrayBox& residLocalCoarseFab = residLocalCoarse[dit()].getFab();

      const Box& bxResid = residLocalCoarseFab.box();
      int intOpLocalInitial = (int) m_opLocalInitial;
      CH_assert(phiLocalCoarseFab.box().contains(grow(bxResid, 1)));
      FORT_EVALOP(CHF_FRA1(residLocalCoarseFab, 0),
                  CHF_CONST_INT(intOpLocalInitial),
                  CHF_CONST_REALVECT(m_dxCoarse),
                  CHF_BOX(bxResid),
                  CHF_CONST_FRA1(phiLocalCoarseFab, 0));
#ifndef NDEBUG
      // dummy statement in order to get around gdb bug
      int dummy_unused = 0; dummy_unused = 0;
#endif
    }
#ifdef CH_USE_TIMER
  TimeLocalResidual.stop();
#endif

  /*****************************************************************/
  // Add up residLocalCoarse to get residGlobalCoarse, which lives on
  // NODEs of
  // m_coarseGlobalDomain = grow(coarsen(m_fineDomain, m_refToCoarse),
  //                             2*ceil((m_minBufferCoarse-1)/2))
  /*****************************************************************/

  if (m_verbose > 0)
    pout() << "doing reduce" << endl;

  // Need pdCoarse, pdFine for generalCopyTo.
  ProblemDomain pdCoarse(m_coarseGlobalDomain);
  if (m_verbose > 0)
    pout() << "pdCoarse == " << pdCoarse << endl;

  // reducePlus(residGlobalCoarse, residLocalCoarse);
  int nproc = numProc();
  Vector<int> procs(nproc);
  for (int iproc = 0; iproc < nproc; iproc++) procs[iproc] = iproc;

  Vector<Box> domainBoxes(nproc);
  domainBoxes.assign(m_coarseGlobalDomain);
  BoxLayout coarseGlobalLayout(domainBoxes, procs);

  LayoutData< Vector<RefCountedPtr<NodeFArrayBox> > > intResidGlobalCoarseData;

  if (m_verbose > 0)
    pout() << "getting generalCopyTo for global" << endl;

#ifdef CH_USE_TIMER
  TimeBarrierPreResid.start();
  barrier();
  TimeBarrierPreResid.stop();
  TimeCopyResidual.start();
#endif

  residLocalCoarse.generalCopyTo(coarseGlobalLayout, // destination layout
                                 intResidGlobalCoarseData, // destination
                                 intvl0, // copy component 0
                                 pdCoarse); // domain of source and destination

  // After this, we no longer need residLocalCoarse.
  // Everything we'll use is in intResidGlobalCoarseData.

#ifdef CH_USE_TIMER
  TimeCopyResidual.stop();
  TimeBarrierPostResid.start();
  barrier();
  TimeBarrierPostResid.stop();
#endif

  if (m_verbose > 0)
    pout() << "got generalCopyTo for global" << endl;

#ifdef CH_USE_TIMER
  TimeReduce.start();
#endif

  // phiGlobalCoarse lives on NODEs of
  // m_coarseGlobalDomainPhi = grow(coarsen(m_fineDomain, m_refToCoarse),
  //                                2*ceil((m_minBufferCoarse-1)/2) + bufferGlobal)
  NodeFArrayBox phiGlobalCoarse(m_coarseGlobalDomainPhi, 1);
  FArrayBox& phiGlobalCoarseFab = phiGlobalCoarse.getFab();
  if (m_verbose > 0)
    pout() << "phiGlobalCoarseFab.box() == "
         << phiGlobalCoarseFab.box() << endl;

  { // scope for residGlobalCoarse
    NodeFArrayBox residGlobalCoarse(m_coarseGlobalDomain, 1);
    FArrayBox& residGlobalCoarseFab = residGlobalCoarse.getFab();
    if (m_verbose > 0)
      pout() << "residGlobalCoarseFab.box() == "
           << residGlobalCoarseFab.box() << endl;
    residGlobalCoarseFab.setVal(0.);
    // only one element in coarseGlobalLayout
    for (DataIterator dcgl = coarseGlobalLayout.dataIterator(); dcgl.ok(); ++dcgl)
      {
        // Everything in this loop can be replaced by:
        // plusReduce(residGlobalCoarse, intResidGlobalCoarseData[dcgl()]);
        Vector< RefCountedPtr<NodeFArrayBox> > intResidGlobalCoarseDataHere =
          intResidGlobalCoarseData[dcgl()];
        plusReduce(residGlobalCoarse, intResidGlobalCoarseDataHere);
        // No need to delete RefCountedPtr.
        // for (int ind = 0; ind < intResidGlobalCoarseDataHere.size(); ind++)
        // delete intResidGlobalCoarseDataHere[ind];
      }
#ifndef NDEBUG
    // viewFAB(&residGlobalCoarseFab);
#endif

    /*****************************************************************/
    // Do global solve.
    /*****************************************************************/

    if (m_verbose > 0)
      pout() << "calling solve on global" << endl;

#ifdef CH_USE_TIMER
    TimeReduce.stop();
    TimeGlobalSolve.start();
    m_solverGlobal.beginTimers();
#endif
    m_solverGlobal.solve(phiGlobalCoarseFab,
                         residGlobalCoarseFab);

#ifdef CH_USE_TIMER
    m_solverGlobal.endTimers();
    TimeGlobalSolve.stop();
#endif

  } // end of scope for residGlobalCoarseFab

#ifndef NDEBUG
  // viewFAB(&phiGlobalCoarseFab);
#endif

  /*****************************************************************/
  // Final local solve.
  /*****************************************************************/

  /*
    Interpolate from
    psi := phiGlobalCoarse - sum(phiLocalCoarse)
    on the boundaries of the boxes
    to get part of the Dirichlet boundary conditions.
  */

  // phiFine will hold the final solution.
  LevelData<NodeFArrayBox> phiFine(phiGrids, 1); // no ghosts

  // From grow(rhsGrids / m_refToCoarse, m_coarseAddRadiusLocal)
  //      = m_coarseLocalGridsBufferInterp
  // to grow(rhsGrids / m_refToCoarse, m_interpLayer)
  //    = m_coarseLocalGridsInterp
  // Recall:
  // m_coarseAddRadiusLocal = m_interpLayer + m_bufferLocal / m_refToCoarse.
  LayoutData< Vector<RefCountedPtr<NodeFArrayBox> > > intCoarseData;
  if (m_verbose > 0)
    pout() << "getting generalCopyTo for coarse" << endl;

#ifdef CH_USE_TIMER
  TimeBarrierPreCrse.start();
  barrier();
  TimeBarrierPreCrse.stop();
  TimeCopyOverlapCrse.start();
#endif

  phiLocalCoarse.generalCopyTo(phiCoarseLocalGridsInterp, // destination layout
                               intCoarseData, // destination
                               intvl0, // copy component 0
                               pdCoarse); // domain of source and destination

  // After this, we no longer need phiLocalCoarse.
  // Everything we'll use is in intCoarseData.

#ifdef CH_USE_TIMER
  TimeCopyOverlapCrse.stop();
  TimeBarrierPostCrse.start();
  barrier();
  TimeBarrierPostCrse.stop();
#endif

  if (m_verbose > 0)
    pout() << "got generalCopyTo for coarse" << endl;

#ifdef CH_USE_TIMER
  TimeBarrierPreCopy.start();
  barrier();
  TimeBarrierPreCopy.stop();
  TimeCopyRhs.start();
#endif

  LevelData<NodeFArrayBox>* rhsFinalPtr;

  if (rhsGrids == phiGrids)
    {
      rhsFinalPtr = &a_rhs;
    }
  else
    {
      // need this to copy rhs from rhsGrids to phiGrids
      rhsFinalPtr = new LevelData<NodeFArrayBox>(phiGrids, 1);
      LevelData<NodeFArrayBox>& rhsFinal = *rhsFinalPtr;
      // LevelData<NodeFArrayBox> rhsFinal(phiGrids, 1);
      sizeLayout = 0;
      for (DataIterator dit = phiGrids.dataIterator(); dit.ok(); ++dit)
        {
          const Box& bx = phiGrids.get(dit());
          sizeLayout += rhsFinal[dit()].size(bx, intvl0);
        }
      pout() << "InfiniteMLCSolver:  rhsFinal "
             << (sizeLayout / sizeof(Real)) << " scalars, "
             << (Real(sizeLayout) / 1048576.) << " MB" << endl;

      // rhs is zero, for boxes not in rhsGrids
      for (DataIterator dit = phiGrids.dataIterator(); dit.ok(); ++dit)
        rhsFinal[dit()].getFab().setVal(0.);
      // Recall that copyTo copies where and only where
      // the CELL-centered boxes intersect.
      a_rhs.copyTo(intvl0, rhsFinal, intvl0);
    }
  LevelData<NodeFArrayBox>& rhsFinal = *rhsFinalPtr;

#ifdef CH_USE_TIMER
  TimeCopyRhs.stop();
  TimeBarrierPostCopy.start();
  barrier();
  TimeBarrierPostCopy.stop();
  TimeFinalLocal.start();
#endif

  for (DataIterator dit = phiGrids.dataIterator(); dit.ok(); ++dit)
    {
      // we'll solve for phiFineFab, which lives on NODEs of
      // a box of m_fineGrids.
      FArrayBox& phiFineFab = a_phi[dit()].getFab();
      const Box& bxCell = phiGrids.get(dit());
      // int nbox = bxCell.size(0);
      // int nboxCoarse = nbox / m_refToCoarse;
      IntVect baseSize = bxCell.size();
      IntVect baseCoarseSize = baseSize / m_refToCoarse;
      // bx is the NODE-centered basic box shifted from [0:nbox]^3
      const Box& bx = phiFineFab.box();
      Box bxCoarse(coarsen(bx, m_refToCoarse));

      // pout() << "Final local solve on " << bxCell << " solved on " << bx << endl;

      Vector< RefCountedPtr<NodeFArrayBox> > intFineDataPatch =
        intFineData[dit()];
      Vector< RefCountedPtr<NodeFArrayBox> > intCoarseDataPatch =
        intCoarseData[dit()];

      /*
        Get intFineBoxes, the list of distinct boxes in intFineDataPatch,
        and intCoarseBoxes, the same boxes (in the same order) coarsened
        by m_refToCoarse and grown by m_interpLayer.
      */
      Vector<Box> intFineBoxes, intCoarseBoxes;
      intFineBoxes.clear(); intCoarseBoxes.clear();
      for (int indFine = 0; indFine < intFineDataPatch.size(); indFine++)
        {
          const Box& srcFineBox = intFineDataPatch[indFine]->box();
          // Check whether srcFineBox is in intFineBoxes.
          // If it isn't, then we'll add it.
          int numIntBoxesFound = intFineBoxes.size();
          bool foundBox = false;
          for (int indInt = 0; indInt < numIntBoxesFound; indInt++)
            {
              if (srcFineBox == intFineBoxes[indInt])
                {
                  foundBox = true;
                  break;
                }
            }
          if (!foundBox)
            {
              // srcFineBox not found in intFineBoxes, so add it.
              intFineBoxes.push_back(srcFineBox);
              Box srcFineBoxCoarsenedInterp =
                grow(coarsen(srcFineBox, m_refToCoarse),
                     m_interpLayer);
              intCoarseBoxes.push_back(srcFineBoxCoarsenedInterp);
            }
        } // end loop over components of intFineDataPatch

      /*
        Find:
        - indicesFine[ind], the list of indices of intFineDataPatch
          that lie on the box intFineBoxes[ind];
        - indicesCoarse[ind], the list of indices of intCoarseDataPatch
          that lie on the box intCoarseBoxes[ind].
       */
      int numIntBoxes = intFineBoxes.size(); // == intCoarseBoxes.size()
      Vector< Vector<int> > indicesFine(numIntBoxes);
      Vector< Vector<int> > indicesCoarse(numIntBoxes);
      for (int indInt = 0; indInt < numIntBoxes; indInt++)
        {
          // pout() << "intFineBoxes[" << indInt << "] == " << intFineBoxes[indInt];
          // pout() << ", intCoarseBoxes[" << indInt << "] == " << intCoarseBoxes[indInt];
          // pout() << endl;
          indicesFine[indInt].clear();
          for (int indFine = 0; indFine < intFineDataPatch.size(); indFine++)
            {
              const Box& srcFineBox = intFineDataPatch[indFine]->box();
              if (srcFineBox == intFineBoxes[indInt])
                indicesFine[indInt].push_back(indFine);
            }
          indicesCoarse[indInt].clear();
          for (int indCoarse = 0; indCoarse < intCoarseDataPatch.size(); indCoarse++)
            {
              const Box& srcCoarseBox = intCoarseDataPatch[indCoarse]->box();
              if (srcCoarseBox == intCoarseBoxes[indInt])
                indicesCoarse[indInt].push_back(indCoarse);
            }
        } // end loop over components of indicesFine and indicesCoarse

      Vector<FArrayBox*> bcFaceFabLo(SpaceDim);
      Vector<FArrayBox*> bcFaceFabHi(SpaceDim);

      // Boundary conditions:
      // for each of the six faces of the box,
      // - get projection of psi on the face plus interpLayer

      SideIterator sit;
      for (sit.begin(); sit.ok(); sit.next())
        {
          Side::LoHiSide hiorloBface = sit();
          IntVect bcEnd = bx.sideEnd(hiorloBface);
          IntVect bxCoarseEnd = bxCoarse.sideEnd(hiorloBface);
          for (int idirBface = 0; idirBface < SpaceDim; idirBface++)
            {

              // the two directions parallel to this face are
              // ipar0 and ipar1.
              int ipar0 = (idirBface == 0) ? 1 : 0;
              int ipar1 = (idirBface == 2) ? 1 : 2;

              // NODE-centered bxbc contains NODEs of boundary face.
              Box bxbc(bx);
              bxbc.setRange(idirBface, bcEnd[idirBface]);

              // FACE-centered bcFace on this face.
              Box bcFace(bxbc);
              bcFace.enclosedCells(ipar0);
              bcFace.enclosedCells(ipar1);

              FArrayBox* bcFaceFabPtr = new FArrayBox(bxbc, 1); // 1 component
              FArrayBox& bcFaceFab = *bcFaceFabPtr;

              // Now add interpolation of phiGlobalCoarseFab.

              Box bxbcCoarseLayered(coarsen(bxbc, m_refToCoarse));
              bxbcCoarseLayered.grow(ipar0, m_interpLayer);
              bxbcCoarseLayered.grow(ipar1, m_interpLayer);

              FArrayBox bcGlobalCoarse(bxbcCoarseLayered, 1);
              bcGlobalCoarse.copy(phiGlobalCoarseFab);

              int nbox0 = baseSize[ipar0];
              int nbox1 = baseSize[ipar1];

              int nbox0coarse = baseCoarseSize[ipar0];
              int nbox1coarse = baseCoarseSize[ipar1];

              FORT_INTERPFACE(bcFaceFab.dataPtr(),
                              bcGlobalCoarse.dataPtr(),
                              m_interpCoeffs,
                              &nbox0, &nbox1,
                              &m_refToCoarse,
                              &nbox0coarse, &nbox1coarse,
                              &m_interpLayer);

              // NODE-centered psibcFaceFab on coarse level
              // from which we interpolate

              // psi := phiGlobalCoarse - sum(phiLocalCoarse)
              // where the sum is over phiLocalCoarse patches that
              // contain bxCoarseLayeredbc.

              FArrayBox bcLocalFab(bxbc, 1); // 1 component
              bcLocalFab.setVal(0.);

              /*
                Define psiFab on bxbcCoarseLayered,
                which consists of the NODEs of coarsened bxbc,
                plus a layer of width interpLayer.

                Then set
                psisubFab = phiGlobalCoarseFab -
                            sum_(ind)(phiLocalCoarseFab[ind])
                where ind loops through all boxes with
                NODEs intersecting this face.
              */

              // FArrayBox psiFab(bxbcCoarseLayered, 1);
              // psiFab.setVal(0.);

              for (int indInt = 0; indInt < numIntBoxes; indInt++)
                {
                  // Box of intFineDataPatch has no buffer layer.
                  const Box& boxIF = intFineBoxes[indInt];
                  // FACE-centered intersection.
                  Box intersectFace(boxIF);
                  intersectFace.surroundingNodes(idirBface);
                  intersectFace &= bcFace;

                  if (! intersectFace.isEmpty() )
                    {
                      // subface is NODE-centered
                      Box subface(intersectFace);
                      subface.surroundingNodes(ipar0);
                      subface.surroundingNodes(ipar1);
                      // check that edges are on lines divisible by
                      // refToCoarse.
                      CH_assert(subface ==
                             refine(coarsen(subface, m_refToCoarse),
                                    m_refToCoarse));

                      FArrayBox bcSubfaceLocalFab(subface, 1);

                      const Box& boxIC = intCoarseBoxes[indInt];

                      Box subfaceCoarse(coarsen(subface, m_refToCoarse));
                      Box subfaceCoarseLayered(subfaceCoarse);
                      subfaceCoarseLayered.grow(ipar0, m_interpLayer);
                      subfaceCoarseLayered.grow(ipar1, m_interpLayer);
                      FArrayBox psisubFab(subfaceCoarseLayered, 1);
                      psisubFab.setVal(0.);

                      /*
                        Set psisubFab = - sum(phiLocalCoarse).
                       */

                      Vector<int>& indicesCoarseHere = indicesCoarse[indInt];
                      for (int indCoarseInd = 0;
                           indCoarseInd < indicesCoarseHere.size();
                           indCoarseInd++)
                        {
                          // Everything in this loop can be replaced by:
                          // psisubFab.minus(intCoarseDataPatch[indicesCoarseHere[indCoarseInd]]->getFab());
                          int indCoarse = indicesCoarseHere[indCoarseInd];

                          CH_assert(intCoarseDataPatch[indCoarse]->box() == boxIC);

                          // phiLocalCoarse has buffer of size m_interpLayer.
                          const NodeFArrayBox& phiLocalCoarseNFAB =
                            *intCoarseDataPatch[indCoarse];
                          CH_assert(phiLocalCoarseNFAB.box() == boxIC);

                          const FArrayBox& phiLocalCoarseFab =
                            phiLocalCoarseNFAB.getFab();

                          CH_assert( phiLocalCoarseFab.box()
                                  .contains(subfaceCoarseLayered) );

                          psisubFab.minus(phiLocalCoarseFab);
                        }

                      /*
                        Now interpolate from psisubFab
                        on subfaceCoarseLayered
                        to bcSubfaceLocalFab on subface.
                      */

                      int nsub0 = subface.size(ipar0) - 1;
                      int nsub1 = subface.size(ipar1) - 1;
                      int nsub0c = nsub0 / m_refToCoarse;
                      int nsub1c = nsub1 / m_refToCoarse;

                      FORT_INTERPFACE(bcSubfaceLocalFab.dataPtr(),
                                      psisubFab.dataPtr(),
                                      m_interpCoeffs,
                                      &nsub0, &nsub1,
                                      &m_refToCoarse,
                                      &nsub0c, &nsub1c,
                                      &m_interpLayer);

                      /*
                        Add sum(phiLocal) to bcSubfaceLocalFab.
                       */

                      Vector<int>& indicesFineHere = indicesFine[indInt];
                      for (int indFineInd = 0;
                           indFineInd < indicesFineHere.size();
                           indFineInd++)
                        {
                          // Everything in this loop can be replaced by:
                          // bcSubfaceLocalFab.plus(intFineDataPatch[indicesFineHere[indFineInd]]->getFab());

                          int indFine = indicesFineHere[indFineInd];

                          const NodeFArrayBox& phiLocalFineNFAB =
                            *intFineDataPatch[indFine];
                          CH_assert(phiLocalFineNFAB.box() == boxIF);

                          // Need to take the CELL-centered intersection
                          // because if you take the NODE-centered intersection
                          // then you may get NODEs on the face.

                          const FArrayBox& phiLocalFineFab =
                            phiLocalFineNFAB.getFab();

                          // bcSubfaceLocalFab lives on subface.
                          CH_assert(phiLocalFineFab.box().contains(subface));
                          bcSubfaceLocalFab.plus(phiLocalFineFab);
                        } // end loop over components of indicesFineHere

                      // Modify edges of bcSubfaceLocalFab that are not
                      // edges of bxbc.
                      halveInternalEdges(bcSubfaceLocalFab, bxbc);

                      bcLocalFab.plus(bcSubfaceLocalFab);
                    } // if (! intersectFace.isEmpty() )
                } // end loop over components of intFineBoxes

              bcFaceFab.plus(bcLocalFab);

              // Now we have filled in all of bcFaceFab.
              // It should be smooth.

              // Copy bcFaceFab pointer to bcFaceFabLo or bcFaceFabHi.

              switch (hiorloBface)
                {
                case Side::Lo:
                  {
                    bcFaceFabLo[idirBface] = bcFaceFabPtr;
                    break;
                  }
                case Side::Hi:
                  {
                    bcFaceFabHi[idirBface] = bcFaceFabPtr;
                    break;
                  }
                default:
                  {
                    MayDay::Error("undefined side");
                  }
                }

            } // end iteration of idirBface over dimensions
        } // end iteration of sit() over sides

      /*
        Solve Poisson problem with inhomogeneous Dirichlet b.c.
      */

      const NodeFArrayBox& rhsNFAB = rhsFinal[dit()];
      const FArrayBox& rhsFab = rhsNFAB.getFab();
      phiFineFab.copy(rhsFab);
      phiFineFab.copy(*(bcFaceFabLo[0]));
      phiFineFab.copy(*(bcFaceFabHi[0]));
      phiFineFab.copy(*(bcFaceFabLo[1]));
      phiFineFab.copy(*(bcFaceFabHi[1]));
      phiFineFab.copy(*(bcFaceFabLo[2]));
      phiFineFab.copy(*(bcFaceFabHi[2]));
      IntVect shiftBox = phiFineFab.box().smallEnd();
      phiFineFab.shift(-shiftBox);
      m_poisson.solveInhomogeneousInPlace(phiFineFab);
      phiFineFab.shift(shiftBox);

      // free up storage
      for (int idirBface = 0; idirBface < SpaceDim; idirBface++)
        {
          delete bcFaceFabLo[idirBface];
          delete bcFaceFabHi[idirBface];
        }

      // Don't need intFineData for this patch anymore.
      // No need to delete RefCountedPtr.
      // for (int indFine = 0; indFine < intFineDataPatch.size(); indFine++)
      // delete intFineDataPatch[indFine];

      // Don't need intCoarseData for this patch anymore.
      // No need to delete RefCountedPtr.
      // for (int ind = 0; ind < intCoarseDataPatch.size(); ind++)
      // delete intCoarseDataPatch[ind];
    } // end loop over patches

  if (! (rhsGrids == phiGrids) )
    delete rhsFinalPtr;

#ifdef CH_USE_TIMER
  TimeFinalLocal.stop();
  TimeMLCsolve_all.stop();

#ifdef MEMORY_USAGE
  Real end_memory = get_memory_usage_from_OS();
#endif
  // Find size of each fine box.
  int baseEachFine = m_fineBaseBox.size(0);
  int actualEachFine = baseEachFine + 2*m_bufferLocalSoln;

  int numRhsHere = rhsGrids.numBoxes(procID());
  int numPhiHere = phiGrids.numBoxes(procID());

  int baseGlobal = m_coarseGlobalDomain.size(0);
  int actualGlobal = m_coarseGlobalDomainPhi.size(0);

  // Find size of box for coarse residual
  int residCoarseBox = baseEachFine / m_refToCoarse + 2*(m_minBufferCoarse-1);

#ifdef MEMORY_USAGE
#ifdef CH_MPI
  // This is for finding avg/min/max of end_memory over all procs.
  // Writes to pout.0.
  Real avg_memory, min_memory, max_memory;
  gather_memory_from_procs(end_memory, avg_memory, min_memory, max_memory);
#endif
#endif

  Real pctTime = 100. / TimeMLCsolve_all.wc_time();

  pout() << "MLC solve done. "
#ifdef MEMORY_USAGE
         << " mem= " << end_memory
#endif
         << "MB wall-clock time: "    << setiosflags(ios::fixed)
         << TimeMLCsolve_all.wc_time() << " sec" << endl;
  pout() << "Pre-barrier initial:   " << TimeBarrierPreInit.wc_time()
         << " seconds"
         << " (" << (TimeBarrierPreInit.wc_time() * pctTime)
         << "%)" << endl;
  pout() << "Set up arrays:         " << TimeSetupArrays.wc_time()
         << " seconds"
         << " (" << (TimeSetupArrays.wc_time() * pctTime)
         << "%)" << endl;
  if (m_bufferLocalSoln == m_bufferLocal)
    pout() << "Init local solve:      " << TimeInitialLocal.wc_time()
           << " seconds"
           << " (" << (TimeInitialLocal.wc_time() * pctTime)
           << "%)"
           << " on " << numRhsHere << " size " << actualEachFine
           << "^3 (base " << baseEachFine << "^3)" << endl;
  else
    pout() << "Init local solve+copy: " << TimeInitialLocal.wc_time()
           << " seconds"
           << " (" << (TimeInitialLocal.wc_time() * pctTime)
           << "%)"
           << " on " << numRhsHere << " size " << actualEachFine
           << "^3 (base " << baseEachFine << "^3)" << endl;
  pout() << "Pre-barrier fine:      " << TimeBarrierPreFine.wc_time()
         << " seconds"
         << " (" << (TimeBarrierPreFine.wc_time() * pctTime)
         << "%)" << endl;
  pout() << "Copying fine overlap:  " << TimeCopyOverlapFine.wc_time()
         << " seconds"
         << " (" << (TimeCopyOverlapFine.wc_time() * pctTime)
         << "%)" << endl;
  pout() << "Post-barrier fine:     " << TimeBarrierPostFine.wc_time()
         << " seconds"
         << " (" << (TimeBarrierPostFine.wc_time() * pctTime)
         << "%)" << endl;
  pout() << "Get local residual:    " << TimeLocalResidual.wc_time()
         << " seconds"
         << " (" << (TimeLocalResidual.wc_time() * pctTime)
         << "%)"
         << " on " << numRhsHere << " size " << residCoarseBox
         << "^3" << endl;
  pout() << "Pre-barrier residual:  " << TimeBarrierPreResid.wc_time()
         << " seconds"
         << " (" << (TimeBarrierPreResid.wc_time() * pctTime)
         << "%)" << endl;
  pout() << "Residual copying:      " << TimeCopyResidual.wc_time()
         << " seconds"
         << " (" << (TimeCopyResidual.wc_time() * pctTime)
         << "%)" << endl;
  pout() << "Post-barrier residual: " << TimeBarrierPostResid.wc_time()
         << " seconds"
         << " (" << (TimeBarrierPostResid.wc_time() * pctTime)
         << "%)" << endl;
  pout() << "Reduce:                " << TimeReduce.wc_time()
         << " seconds"
         << " (" << (TimeReduce.wc_time() * pctTime)
         << "%)"
         << " on 1 size " << actualGlobal
         << "^3 (base " << baseGlobal << "^3)" << endl;
  pout() << "Global solve:          " << TimeGlobalSolve.wc_time()
         << " seconds"
         << " (" << (TimeGlobalSolve.wc_time() * pctTime)
         << "%)"
         << " on 1 size " << actualGlobal
         << "^3 (base " << baseGlobal << "^3)" << endl;
  pout() << "Pre-barrier crse:      " << TimeBarrierPreCrse.wc_time()
         << " seconds"
         << " (" << (TimeBarrierPreCrse.wc_time() * pctTime)
         << "%)" << endl;
  pout() << "Copying crse overlap:  " << TimeCopyOverlapCrse.wc_time()
         << " seconds"
         << " (" << (TimeCopyOverlapCrse.wc_time() * pctTime)
         << "%)" << endl;
  pout() << "Post-barrier crse:     " << TimeBarrierPostCrse.wc_time()
         << " seconds"
         << " (" << (TimeBarrierPostCrse.wc_time() * pctTime)
         << "%)" << endl;
  pout() << "Pre-barrier copy:      " << TimeBarrierPreCopy.wc_time()
         << " seconds"
         << " (" << (TimeBarrierPreCopy.wc_time() * pctTime)
         << "%)" << endl;
  pout() << "Copy rhs:              " << TimeCopyRhs.wc_time()
         << " seconds"
         << " (" << (TimeCopyRhs.wc_time() * pctTime)
         << "%)" << endl;
  pout() << "Post-barrier copy:     " << TimeBarrierPostCopy.wc_time()
         << " seconds"
         << " (" << (TimeBarrierPostCopy.wc_time() * pctTime)
         << "%)" << endl;
  pout() << "Final local solve:     " << TimeFinalLocal.wc_time()
         << " seconds"
         << " (" << (TimeFinalLocal.wc_time() * pctTime)
         << "%)"
         << " on " << numPhiHere << " size " << baseEachFine << "^3" << endl;
#endif

  // Timer::TimerSummary(); // causes segfault
}


void InfiniteMLCSolver::plusReduce(NodeFArrayBox&                a_sum,
                                   const Vector< RefCountedPtr<NodeFArrayBox> >& a_data)
{
  FArrayBox& sumFab = a_sum.getFab();
  int ncomp = sumFab.nComp();
  // CH_assert(a_data.nComp() == ncomp);
  // const Box& sumFabNodes = sumFab.box();
  sumFab.setVal(0.);
  for (int ind = 0; ind < a_data.size(); ind++)
    {
      const FArrayBox& dataFab = a_data[ind]->getFab();
      CH_assert(dataFab.nComp() == ncomp);
      // Add over domain of intersection of sumFab and dataFab.
      sumFab.plus(dataFab, 0, 0, ncomp);
    }
}


void InfiniteMLCSolver::halveInternalEdges(FArrayBox&  a_fab,
                                           const Box&  a_bxFace)
{
  IntVect faceSizes(a_bxFace.size());
  int* ipar = new int[2];
  // ipar[0], ipar[1] are the directions with sides of length > 1
  ipar[0] = (faceSizes[0] == 1) ? 1 : 0;
  ipar[1] = (faceSizes[2] == 1) ? 1 : 2;

  const Box& fabBox(a_fab.box());
  CH_assert(a_bxFace.contains(fabBox));
  for (int ipardir = 0; ipardir < 2; ipardir++)
    {
      int idimEdge = ipar[ipardir];
      // If left edge of fabBox is NOT at left edge of a_bxFace,
      // then halve values of a_fab at left edge.
      int thisSmallEnd = fabBox.smallEnd(idimEdge);
      if (thisSmallEnd != a_bxFace.smallEnd(idimEdge))
        {
          Box edgeNodes(fabBox);
          edgeNodes.setBig(idimEdge, thisSmallEnd);
          a_fab.mult(0.5, edgeNodes);
          // pout() << "correction on " << edgeNodes << endl;
        }
      // If right edge of fabBox is NOT at right edge of a_bxFace,
      // then halve values of a_fab at left edge.
      int thisBigEnd = fabBox.bigEnd(idimEdge);
      if (thisBigEnd != a_bxFace.bigEnd(idimEdge))
        {
          Box edgeNodes(fabBox);
          edgeNodes.setSmall(idimEdge, thisBigEnd);
          a_fab.mult(0.5, edgeNodes);
          // pout() << "correction on " << edgeNodes << endl;
        }
    }
  delete ipar;
}

#include "NamespaceFooter.H"
