#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// InfiniteNodeSolver.cpp
// petermc, 26 Jan 2004
// Compile with XTRACPPFLAGS=-DTIMER if you want to use the Timer.

#ifdef TIMER
#include <iostream>
#include <iomanip> // not in AMRGodunov
using std::ifstream;
using std::ios;
#endif

#ifdef HALEM_PROC_SPEED
#include <cstdio>
#include <sys/sysinfo.h>
#include <machine/hal_sysinfo.h>
#endif

#ifdef TIMER
#include "memusage.H"
#include "Timer.H"
#endif

#include "InfiniteNodeSolver.H"
#include "InfiniteDomainF.H"
#include "MayDay.H"
#include "Projections.H"
// #include "TotalChargeF_F.H"
using std::cout;
using std::endl;
using std::cerr;

#ifdef CH_MPI
#include "SPMD.H"
#endif

#include "NamespaceHeader.H"

#ifdef TIMER
Timer TimeSingledefine_all   ("NodeSolver define all", 20);
Timer TimeDefineArrays       ("Define Arrays",       TimeSingledefine_all);
Timer TimeDefineNormals      ("Define Normals",      TimeSingledefine_all);
Timer TimePoisson            ("Poisson",             TimeSingledefine_all);
Timer TimeMultipoles         ("Multipoles",          TimeSingledefine_all);

Timer TimeSinglesolve_all    ("NodeSolver solve all", 22);
Timer TimeArrays             ("Arrays",              TimeSinglesolve_all);
Timer TimeHomoDirichlet      ("Homo Dirichlet",      TimeSinglesolve_all);
Timer TimeNormals            ("Normals",             TimeSinglesolve_all);
Timer TimeMultipoleCoeffs    ("MultipoleCoeffs",     TimeSinglesolve_all);
Timer TimeMultipoleEval      ("MultipoleEval",       TimeSinglesolve_all);
Timer TimeMultipoleEvalOuter ("MultipoleEvalOuter",  TimeSinglesolve_all);
Timer TimeMultipoleClear     ("MultipoleClear",      TimeSinglesolve_all);
Timer TimeInhomoDirichlet    ("Inhomo Dirichlet",    TimeSinglesolve_all);
Timer TimeCopying            ("Copying",             TimeSinglesolve_all);
#endif

// ---------------------------------------------------------
// default constructor
InfiniteNodeSolver::InfiniteNodeSolver()
{
  setDefaultValues();
}


// ---------------------------------------------------------
void
InfiniteNodeSolver::define(const Box&        a_domain,
                           int               a_s1,
                           int               a_s2,
                           int               a_patchSize,
                           int               a_dstFaceCoarsening,
                           const RealVect&   a_dx,
                           PoissonDirichlet::OperatorType      a_op,
                           int               a_multipoleOrder,
                           int               a_degreeNormalDerivative,
                           int               a_interpBorder,
                           bool              a_parallel,
                           bool              a_getOuterCoarse,
                           int               a_refToCoarse,
                           int               a_coarseAddRadius)
{
  // The current model is to set up just one object of the
  // InfiniteNodeSolver class and use it for all boxes of the
  // same size and shape.  That way, no repetition of calculation
  // of the matrix, etc.

#ifdef TIMER
  TimeSingledefine_all.start();
  TimeDefineArrays.start();
#endif

  CH_assert(a_domain.ixType() == IndexType::TheNodeType());
  clearMemory();
  m_domain = a_domain; // NODE-centered now
  IntVect domainLength = m_domain.bigEnd() - m_domain.smallEnd();
//   if (domainLength[0] != domainLength[1] ||
//       domainLength[0] != domainLength[2])
//     {
//       cerr << "InfiniteNodeSolver requires cubic domain"
//            << endl;
//       MayDay::Error("returning");
//       return;
//     }
  m_s1 = a_s1;
  m_s2 = a_s2;
  m_patchSize = a_patchSize;
  // m_buffer = a_buffer;
  m_dx = a_dx;
  m_op = a_op;
  m_multipoleOrder = a_multipoleOrder;
  m_degreeNormalDerivative = a_degreeNormalDerivative;
  m_interpBorder = a_interpBorder;

  m_srcBox = grow(m_domain, m_s1); // == D1, NODE-centered
  m_dstBox = grow(m_srcBox, m_s2); // == D2, NODE-centered

  m_phirhsSrcFab.define(m_srcBox, 1);
  m_phirhsDstFab.define(m_dstBox, 1);

  m_dstFaceCoarsening = m_patchSize;

  // Vector<FArrayBox*> m_derivs;
  m_derivs.resize(2*SpaceDim);
  int faceID = 0;
  for (SideIterator sit; sit.ok(); sit.next())
    for (int idir = 0; idir < SpaceDim; idir++)
      {
        Box bxFace = bdryBox(m_srcBox, idir, sit(), 1);
        m_derivs[faceID] = new FArrayBox(bxFace, 1);
        faceID++;
      }

#ifdef TIMER
  TimeDefineArrays.stop();
  TimeDefineNormals.start();
#endif

  // Normals m_normals;
  m_normals.define(m_srcBox, m_dx, m_degreeNormalDerivative);

#ifdef TIMER
  TimeDefineNormals.stop();
  TimePoisson.start();
#endif

  // PoissonDirichlet m_solverPoissonSrc, m_solverPoissonDst;
  m_solverPoissonSrc.define(m_srcBox, m_dx, m_op);
  m_solverPoissonDst.define(m_dstBox, m_dx, m_op);

#ifdef TIMER
  TimePoisson.stop();
  TimeMultipoles.start();
#endif

  m_parallel = a_parallel;
  if (a_getOuterCoarse)
    {
      m_refToCoarse = a_refToCoarse;
      m_getOuterCoarse = (a_coarseAddRadius >= (m_s1 + m_s2) / a_refToCoarse);
    }
  else
    {
      m_refToCoarse = 1;
      m_getOuterCoarse = false;
    }

  // Multipoles m_multipoles;
  if (! m_getOuterCoarse)
    {
      m_coarseAddRadius = 0;
      m_multipoles.define(m_srcBox,
                          m_dstBox,
                          m_dx,
                          m_patchSize,
                          m_multipoleOrder,
                          m_dstFaceCoarsening,
                          m_interpBorder,
                          m_parallel);
    }
  else // m_getOuterCoarse
    {
      // petermc, 6 oct 2004, reduced coarseAddRadius
      // by floor((m_s1 + m_s2) / m_refToCoarse)
      m_coarseAddRadius = a_coarseAddRadius - ((m_s1 + m_s2) / m_refToCoarse);
      m_multipoles.define(m_srcBox,
                          m_dstBox,
                          m_dx,
                          m_patchSize,
                          m_multipoleOrder,
                          m_dstFaceCoarsening,
                          m_interpBorder,
                          m_parallel,
                          true,
                          m_refToCoarse,
                          m_coarseAddRadius);
    }

#ifdef TIMER
  TimeMultipoles.stop();
  TimeSingledefine_all.stop();
  Real end_memory = get_memory_usage_from_OS();

#ifdef CH_MPI
  Real avg_memory, min_memory, max_memory;
  gather_memory_from_procs(end_memory, avg_memory, min_memory, max_memory);
#endif

  Real pctTime = 100. / TimeSingledefine_all.wc_time();

  pout() << "InfiniteNodeSolver define done.  mem= " << end_memory
         << "MB wall-clock time: "    << setiosflags(ios::fixed)
         << TimeSingledefine_all.wc_time()      << " sec" << endl;
  pout() << "Define arrays:        " << TimeDefineArrays.wc_time()
         << " seconds"
         << " (" << (TimeDefineArrays.wc_time() * pctTime)
         << "%)" << endl;
  pout() << "Define normals:       " << TimeDefineNormals.wc_time()
         << " seconds"
         << " (" << (TimeDefineNormals.wc_time() * pctTime)
         << "%)" << endl;
  pout() << "Define poisson:       " << TimePoisson.wc_time()
         << " seconds"
         << " (" << (TimePoisson.wc_time() * pctTime)
         << "%)" << endl;
  pout() << "Define multipoles:    " << TimeMultipoles.wc_time()
         << " seconds"
         << " (" << (TimeMultipoles.wc_time() * pctTime)
         << "%)" << endl;
  pout() << resetiosflags(ios::fixed);

  TimeSingledefine_all.clear();
  TimeDefineArrays.clear();
  TimeDefineNormals.clear();
  TimePoisson.clear();
  TimeMultipoles.clear();

#endif

  m_isDefined = true;
}


// ---------------------------------------------------------
// complete constructor
InfiniteNodeSolver::InfiniteNodeSolver(const Box&        a_domain,
                                       int               a_s1,
                                       int               a_s2,
                                       int               a_patchSize,
                                       int               a_dstFaceCoarsening,
                                       const RealVect&   a_dx,
                                       PoissonDirichlet::OperatorType   a_op,
                                       int               a_multipoleOrder,
                                       int               a_degreeNormalDerivative,
                                       int               a_interpBorder,
                                       bool              a_parallel,
                                       bool              a_getOuterCoarse,
                                       int               a_refToCoarse,
                                       int               a_coarseAddRadius)
{
  setDefaultValues();
  define(a_domain, a_s1, a_s2, a_patchSize, a_dstFaceCoarsening,
         a_dx, a_op, a_multipoleOrder,
         a_degreeNormalDerivative, a_interpBorder, a_parallel,
         a_getOuterCoarse, a_refToCoarse, a_coarseAddRadius);
}


// ---------------------------------------------------------
void
InfiniteNodeSolver::define(const Box&        a_domain,
                           const RealVect&   a_dx,
                           int               a_order,
                           bool              a_parallel,
                           bool              a_getOuterCoarse,
                           int               a_refToCoarse,
                           int               a_coarseAddRadius)
{
  IntVect domainLength = a_domain.bigEnd() - a_domain.smallEnd();
  int s1, s2, patchSize;
  int bufferMin = 0;
  bool alphaZero = true;
  if (alphaZero)
    {
      s1 = 0;
      int errcode = 0;
      FORT_GET0BUFFERSIZE3D(domainLength.dataPtr(), &bufferMin,
                            &s2, &patchSize, &errcode);
      if (errcode != 0)
        {
          cerr << "get0buffersize returned error code " << errcode << endl;
          MayDay::Error("returning");
          return;
        }
    }
  else
    {
      int s1s2;
      int errcode = 0;
      FORT_GETBUFFERSIZE3D(domainLength.dataPtr(), &bufferMin,
                           &s1s2, &patchSize, &errcode);
      if (errcode != 0)
        {
          cerr << "getbuffersize returned error code " << errcode << endl;
          MayDay::Error("returning");
          return;
        }
      s1 = s1s2 / 2;
      s2 = s1s2 / 2;
    }
  PoissonDirichlet::OperatorType op;
  if (a_order == 2)
    {
      op = PoissonDirichlet::POINT7;
    }
  else if (a_order == 4)
    {
      op = PoissonDirichlet::POINT19;
    }
  else
    {
      cout << "InfiniteNodeSolver must have a_order either 2 or 4, not "
           << a_order << endl;
      MayDay::Error("returning");
      return;
    }
  int multipoleOrder = 2 * a_order - 1;
  int degreeNormalDerivative = a_order;
  int interpBorder = a_order - 1;
  define(a_domain, s1, s2, patchSize, patchSize,
         a_dx, op, multipoleOrder,
         degreeNormalDerivative, interpBorder, a_parallel,
         a_getOuterCoarse, a_refToCoarse, a_coarseAddRadius);
}


// ---------------------------------------------------------
// abbreviated constructor
InfiniteNodeSolver::InfiniteNodeSolver(const Box&        a_domain,
                                       const RealVect&   a_dx,
                                       int               a_order,
                                       bool              a_parallel,
                                       bool              a_getOuterCoarse,
                                       int               a_refToCoarse,
                                       int               a_coarseAddRadius)
{
  setDefaultValues();
  define(a_domain, a_dx, a_order,
         a_getOuterCoarse, a_refToCoarse, a_coarseAddRadius);
}


// ---------------------------------------------------------
// destructor
InfiniteNodeSolver::~InfiniteNodeSolver()
{
  clearMemory();
}


// ---------------------------------------------------------
bool
InfiniteNodeSolver::isDefined() const
{
  return m_isDefined;
}


// ---------------------------------------------------------
void
InfiniteNodeSolver::setDefaultValues()
{
  m_isDefined = false;
  m_refToCoarse = -1;
  // m_dx = -1.0;
  m_verbose = 0;
}


// ---------------------------------------------------------
void
InfiniteNodeSolver::setVerbose(int a_verbose)
{
  m_verbose = a_verbose;
}


// ---------------------------------------------------------
// returns InfiniteNodeSolver to a basically undefined state
void
InfiniteNodeSolver::clearMemory()
{
  // if (m_levelOpPtr != NULL)
  // delete m_levelOpPtr;
  // m_levelOpPtr = NULL;
  if (m_isDefined)
    {
      // delete stuff
      int faceID = 0;
      for (SideIterator sit; sit.ok(); sit.next())
        for (int idir = 0; idir < SpaceDim; idir++)
          {
            delete m_derivs[faceID];
            faceID++;
          }
    }
}


// ---------------------------------------------------------
// solve on just this level, inhomogeneous bcs.  arguments live on m_domain.
void
InfiniteNodeSolver::solve(FArrayBox& a_phi,
                          const FArrayBox& a_rhs)
{
  solveUnified(a_phi, NULL, a_rhs);
}


// ---------------------------------------------------------
// solve on just this level, inhomogeneous bcs
void
InfiniteNodeSolver::solveWithCoarse(FArrayBox& a_phi,
                                    FArrayBox& a_phiCoarse,
                                    const FArrayBox& a_rhs)
{
  solveUnified(a_phi, &a_phiCoarse, a_rhs);
}


// ---------------------------------------------------------
// solve on just this level, inhomogeneous bcs
void
InfiniteNodeSolver::solveUnified(FArrayBox& a_phi,
                                 FArrayBox* a_phiCoarsePtr,
                                 const FArrayBox& a_rhs)
{
#ifdef TIMER
  m_uses++;
  TimeSinglesolve_all.start();
  TimeArrays.start();
#endif

  if (m_verbose > 0)
    cout << "in InfiniteNodeSolver::solve" << endl;
  CH_assert(isDefined());
  CH_assert(a_phi.nComp() == 1);
  if (a_phiCoarsePtr != NULL)
    CH_assert(a_phiCoarsePtr->nComp() == 1);
  CH_assert(a_rhs.nComp() == 1);

  // petermc, 25 Jan 2005:  set rhs box to be same as domain, no shifting
  CH_assert(a_rhs.box() == m_domain);

  /*
    First step in James's algorithm:
    solve op(phi) = rhs for phi on m_srcBox.
  */

  // petermc, 25 Jan 2005:  no shifting of base
  // m_phirhsSrcFab has a domain at least as large as a_rhs.
  m_phirhsSrcFab.setVal(0.);
  CH_assert(m_phirhsSrcFab.box().contains(a_rhs.box()));
  // petermc, 1 Dec 2004, made copy over innerBox only,
  // because if a_rhs is nonzero outside innerBox, then
  // solveHomogeneousInPlace returns funny results on the faces.
  Box innerBox = grow(m_phirhsSrcFab.box(), -1);
  m_phirhsSrcFab.copy(a_rhs, innerBox);

#ifdef TIMER
  TimeArrays.stop();
  TimeHomoDirichlet.start();
#endif
  // input:  m_phirhsSrcFab contains rhs
  m_solverPoissonSrc.solveHomogeneousInPlace(m_phirhsSrcFab);
  // output:  m_phirhsSrcFab contains phi
#ifdef TIMER
  TimeHomoDirichlet.stop();
  TimeNormals.start();
#endif

  /*
    Second step in James's algorithm:
    get normal derivatives of phi on faces of m_srcBox.
  */

  // Vector<FArrayBox*> m_derivs;
  m_normals.evalNormal(m_derivs, m_phirhsSrcFab);
#ifdef TIMER
  TimeNormals.stop();
  TimeMultipoleCoeffs.start();
#endif

  /*
    Third step in James's algorithm:
    get potential on faces of m_dstBox
    due to normal derivatives on faces of m_srcBox.
  */

  m_multipoles.setCoeffs(m_derivs);

  m_phirhsDstFab.setVal(0.);
  CH_assert(m_phirhsDstFab.box().contains(a_rhs.box()));
  m_phirhsDstFab.copy(a_rhs);

#ifdef TIMER
  TimeMultipoleCoeffs.stop();
  TimeMultipoleEval.start();
#endif
  // Get result on faces of m_dstBox.  We write over phirhsDstFab on faces.
  m_multipoles.eval(m_phirhsDstFab);

#ifdef TIMER
  TimeMultipoleEval.stop();
  TimeMultipoleEvalOuter.start();
#endif
  if (a_phiCoarsePtr != NULL && m_getOuterCoarse)
    // Note that this leaves a_phiCoarse undefined in the interior
    // where it overlaps with D2 ...
    m_multipoles.evalOuterCoarse(*a_phiCoarsePtr);

#ifdef TIMER
  TimeMultipoleEvalOuter.stop();
  TimeMultipoleClear.start();
#endif
  m_multipoles.clearCoeffs();

#ifdef TIMER
  TimeMultipoleClear.stop();
  TimeInhomoDirichlet.start();
#endif
  /*
    Fourth step in James's algorithm:
    solve op(phi) = rhs for phi on m_dstBox.
  */

  m_solverPoissonDst.solveInhomogeneousInPlace(m_phirhsDstFab);

#ifdef TIMER
  TimeInhomoDirichlet.stop();
  TimeCopying.start();
#endif
  // petermc, 25 Jan 2005:  no shifting
  CH_assert(m_phirhsDstFab.box().contains(a_phi.box()));
  a_phi.copy(m_phirhsDstFab);

  // ... Now we fill in D2.
  // petermc, 21 Jan 2005:  remove condition && m_getOuterCoarse
  if (a_phiCoarsePtr != NULL) // && m_getOuterCoarse)
    projectToCoarse(*a_phiCoarsePtr, m_phirhsDstFab, m_refToCoarse);

#ifdef TIMER
  TimeCopying.stop();
  TimeSinglesolve_all.stop();
#endif
}


// ---------------------------------------------------------
void
InfiniteNodeSolver::beginTimers()
{
#ifdef TIMER
  m_uses = 0;
#endif
}


// ---------------------------------------------------------
void
InfiniteNodeSolver::endTimers()
{
#ifdef TIMER

  Real end_memory = get_memory_usage_from_OS();

#ifdef CH_MPI
  Real avg_memory, min_memory, max_memory;
  gather_memory_from_procs(end_memory, avg_memory, min_memory, max_memory);
#endif

  int dstLength = m_dstBox.bigEnd(0) - m_dstBox.smallEnd(0);
  int srcLength = m_srcBox.bigEnd(0) - m_srcBox.smallEnd(0);

  Real pctTime = 100. / TimeSinglesolve_all.wc_time();

  pout() << "InfiniteNodeSolver solve done.  mem= " << end_memory
         << "MB wall-clock time: "    << setiosflags(ios::fixed)
         << TimeSinglesolve_all.wc_time()      << " sec" << endl;
  pout() << "Set up arrays:            " << TimeArrays.wc_time()
         << " seconds"
         << " (" << (TimeArrays.wc_time() * pctTime)
         << "%) on " << m_uses << endl;
  pout() << "Homo Dirichlet:           " << TimeHomoDirichlet.wc_time()
         << " seconds"
         << " (" << (TimeHomoDirichlet.wc_time() * pctTime)
         << "%) on " << m_uses
         << " size " << srcLength << endl;
  pout() << "Calculate normals:        " << TimeNormals.wc_time()
         << " seconds"
         << " (" << (TimeNormals.wc_time() * pctTime)
         << "%) on " << m_uses
         << " size " << srcLength << endl;
  pout() << "Get multipole coeffs:     " << TimeMultipoleCoeffs.wc_time()
         << " seconds"
         << " (" << (TimeMultipoleCoeffs.wc_time() * pctTime)
         << "%) on " << m_uses
         << " from " << srcLength << " to " << dstLength << endl;
  pout() << "Evaluate multipole:       " << TimeMultipoleEval.wc_time()
         << " seconds"
         << " (" << (TimeMultipoleEval.wc_time() * pctTime)
         << "%) on " << m_uses
         << " from " << srcLength << " to " << dstLength << endl;
  pout() << "Evaluate outer multipole: " << TimeMultipoleEvalOuter.wc_time()
         << " seconds"
         << " (" << (TimeMultipoleEvalOuter.wc_time() * pctTime)
         << "%) on " << m_uses
         << " layer " << m_coarseAddRadius << endl;
  pout() << "Clear multipole:          " << TimeMultipoleClear.wc_time()
         << " seconds"
         << " (" << (TimeMultipoleClear.wc_time() * pctTime)
         << "%) on " << m_uses
         << " from " << srcLength << " to " << dstLength << endl;
  pout() << "Inhomo Dirichlet:         " << TimeInhomoDirichlet.wc_time()
         << " seconds"
         << " (" << (TimeInhomoDirichlet.wc_time() * pctTime)
         << "%) on " << m_uses
         << " size " << dstLength << endl;
  pout() << resetiosflags(ios::fixed);

  TimeSinglesolve_all.clear();
  TimeArrays.clear();
  TimeHomoDirichlet.clear();
  TimeNormals.clear();
  TimeMultipoleCoeffs.clear();
  TimeMultipoleEval.clear();
  TimeMultipoleEvalOuter.clear();
  TimeMultipoleClear.clear();
  TimeInhomoDirichlet.clear();
  TimeCopying.clear();

#endif
}

#include "NamespaceFooter.H"
