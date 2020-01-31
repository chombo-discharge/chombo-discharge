#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// DTGraves, Tues, July 6, 1999

#include <cmath>
#include <iostream>
#include <iomanip>
using std::cout;
using std::endl;
using std::setprecision;
using std::ios;
using std::setiosflags;

#include "parstream.H"
#include "AMRSolver.H"
#include "LayoutIterator.H"
#include "NamespaceHeader.H"

/*****************/
/*****************/
void AMRSolver::setDefaultValues()
{
  m_numLevels = -1;
  m_finestLevel = -1;

  m_refRatio.resize(0);

  m_tolerance = 1.0e-10;
  m_numVCyclesBottom = 1;
#ifdef CH_USE_FLOAT
  m_tolerance = 10.*sqrt(m_tolerance);
#endif
  m_errorTolerance = m_tolerance;
  m_operatorTolerance = 1.0e-4;

  m_maxIter = 42;
  m_minIter = 4;

  m_numSmoothUp = 4;
  m_numSmoothDown = 4;

  // default is to use maxNorm
  m_normType = 0;

  m_ncomp = 1;
  m_verbose = false;

  m_isDefined = false;
}

/*****************/
/*****************/
AMRSolver::AMRSolver():m_amrmgLevel()
{
  setDefaultValues();
}

/*****************/
/*****************/
void AMRSolver::clear()
{
  for (int ilev = 0; ilev < m_amrmgLevel.size(); ilev++)
    {
      if (m_amrmgLevel[ilev] != NULL)
        {
          delete m_amrmgLevel[ilev];
        }
    }
}

/*****************/
/*****************/
AMRSolver::~AMRSolver()
{
  clear();
}

/*****************/
/*****************/
AMRSolver::AMRSolver(const Vector<DisjointBoxLayout>& a_gridsLevel,
                     const Vector<Box>&               a_domainLevel,
                     const Vector<Real>&              a_dxLevel,
                     const Vector<int>&               a_refRatio,
                     int                              a_numLevels,
                     int                              a_lBase,
                     const LevelOp* const             a_opin,
                     int                              a_ncomp,
                     bool                             a_limitCoarsening)
  :m_amrmgLevel()
{
  setDefaultValues();

  Vector<ProblemDomain> physdomains(a_domainLevel.size());
  for (int lev = 0; lev < physdomains.size(); lev++)
    {
      physdomains[lev] = ProblemDomain(a_domainLevel[lev]);
    }

  define(a_gridsLevel, physdomains, a_dxLevel,
         a_refRatio, a_numLevels, a_lBase, a_opin, a_ncomp,
         a_limitCoarsening);
}

/*****************/
/*****************/
AMRSolver::AMRSolver(const Vector<DisjointBoxLayout>& a_gridsLevel,
                     const Vector<ProblemDomain>&     a_domainLevel,
                     const Vector<Real>&              a_dxLevel,
                     const Vector<int>&               a_refRatio,
                     int                              a_numLevels,
                     int                              a_lBase,
                     const LevelOp* const             a_opin,
                     int                              a_ncomp,
                     bool                             a_limitCoarsening)
  :m_amrmgLevel()
{
  setDefaultValues();

  define(a_gridsLevel, a_domainLevel, a_dxLevel,
         a_refRatio, a_numLevels, a_lBase, a_opin, a_ncomp,
         a_limitCoarsening);
}

/*****************/
/*****************/
void AMRSolver::define(const Vector<DisjointBoxLayout>& a_gridsLevel,
                       const Vector<Box>&               a_domainLevel,
                       const Vector<Real>&              a_dxLevel,
                       const Vector<int>&               a_refRatio,
                       int                              a_numLevels,
                       int                              a_lBase,
                       const LevelOp* const             a_opin,
                       int                              a_ncomp,
                       bool                             a_limitCoarsening)
{
  setDefaultValues();

  Vector<ProblemDomain> physdomains(a_domainLevel.size());
  for (int lev = 0; lev < physdomains.size(); lev++)
    {
      physdomains[lev] = ProblemDomain(a_domainLevel[lev]);
    }

  define(a_gridsLevel, physdomains, a_dxLevel, a_refRatio, a_numLevels,
         a_lBase, a_opin, a_ncomp, a_limitCoarsening);
}

/*****************/
/*****************/
void AMRSolver::define(const Vector<DisjointBoxLayout>& a_gridsLevel,
                       const Vector<ProblemDomain>&     a_domainLevel,
                       const Vector<Real>&              a_dxLevel,
                       const Vector<int>&               a_refRatio,
                       int                              a_numLevels,
                       int                              a_lBase,
                       const LevelOp* const             a_opin,
                       int                              a_ncomp,
                       bool                             a_limitCoarsening)
{
  CH_TIME("AMRSolver::define");
  clear();

  m_isDefined = true;

  m_gridsLevel = a_gridsLevel;
  m_domainLevel = a_domainLevel;
  m_numLevels = a_numLevels;
  m_finestLevel = a_numLevels - 1;

  m_refRatio = a_refRatio;
  m_dxLevel = a_dxLevel;

  m_lBase = a_lBase;

  m_ncomp = a_ncomp;

  // define amrlevelmgs
  m_amrmgLevel.resize(m_numLevels, NULL);

  // (dfm -- 12/10/01) if lBase > 1, then we only need to define levels
  // from lBase-1 (needed for c-f boundary conditions)
  int startLev= a_lBase;
  // if (startLev>0) startLev--;

  for (int ilev = startLev; ilev < m_numLevels; ilev++)
    {
      m_amrmgLevel[ilev] = new AMRLevelMG(this, ilev, a_opin, m_ncomp);
    }

  CH_assert(m_lBase >= 0);

  // define levelsolver stuff
  const DisjointBoxLayout& levsolvGrids= m_gridsLevel[m_lBase];
  const ProblemDomain&      levsolvDom  = m_domainLevel[m_lBase];
  const Real&     levsolvDx   = m_dxLevel[m_lBase];
  const DisjointBoxLayout* levsolvBase = NULL;

  // int levsolvRef = 2;
  int levsolvRef = -1;

  if (m_lBase > 0)
    {
      levsolvBase = &m_gridsLevel[m_lBase-1];
      levsolvRef = m_refRatio[m_lBase-1];
    }

  m_levelSolver.define(levsolvGrids, levsolvBase,
                       levsolvDom,   levsolvDx,
                       levsolvRef, a_opin, m_ncomp,
                       a_limitCoarsening);

  m_levelSolver.setMaxIter(m_numVCyclesBottom-1);
  m_levelSolver.setNumSmoothUp(m_numSmoothUp);
  m_levelSolver.setNumSmoothDown(m_numSmoothDown);
}

/*****************/
///has full define function been called?
/*****************/
bool AMRSolver::isDefined() const
{
  return m_isDefined;
}

/**
   Set number of multigrid smoothings on way up v-cycle.
   Default == 4
*/
void AMRSolver::setNumSmoothUp(int a_numSmoothUp)
{
  CH_assert(a_numSmoothUp >= 0);

  m_numSmoothUp = a_numSmoothUp;

  for (int ilev = m_lBase; ilev < m_numLevels; ilev++)
    {
      m_amrmgLevel[ilev]->setNumSmoothUp(m_numSmoothUp);
    }
  if (m_levelSolver.isDefined())
  {
    m_levelSolver.setNumSmoothUp(m_numSmoothUp);
  }

}

void AMRSolver::setNumVCyclesBottom(int a_numVCyclesBottom)
{
  CH_assert(a_numVCyclesBottom >= 0);

  m_numVCyclesBottom = a_numVCyclesBottom;

  if (m_levelSolver.isDefined())
    {
      m_levelSolver.setMaxIter(m_numVCyclesBottom);
    }
}

/**
   Set number of multigrid smoothings on way down v-cycle;
   Default == 4.
*/
void AMRSolver::setNumSmoothDown(int a_numSmoothDown)
{
  CH_assert(a_numSmoothDown >= 0);

  m_numSmoothDown = a_numSmoothDown;

  for (int ilev = m_lBase; ilev < m_numLevels; ilev++)
    {
      m_amrmgLevel[ilev]->setNumSmoothDown(m_numSmoothDown);
    }
  if (m_levelSolver.isDefined())
    {
      m_levelSolver.setNumSmoothDown(m_numSmoothDown);
    }
}

/**
   Set tolerance of iterative solution.  Default is 1.0e-10.
**/
void AMRSolver::setTolerance(Real a_tolerance)
{
  CH_assert(a_tolerance >= 0);

  m_tolerance = a_tolerance;
  m_errorTolerance = a_tolerance;
}

/**
   Set tolerance for stopping due to hung convergence
   Default is 1.0e-5
**/
void AMRSolver::setOperatorTolerance(Real a_tolerance)
{
  CH_assert(a_tolerance >= 0);

  m_operatorTolerance = a_tolerance;
}

/**
   Set max number of iterations.  default is 42.
**/
void AMRSolver::setMaxIter(int a_maxIter)
{
  CH_assert(a_maxIter >= 0);

  m_maxIter = a_maxIter;
}

/**
   Set min number of iterations.  default is 4.
**/
void AMRSolver::setMinIter(int a_minIter)
{
  CH_assert(a_minIter >= 0);

  m_minIter = a_minIter;
}

void
AMRSolver::setConvergenceMetric(Real a_metric, int a_comp)
{
    // may need to resize this vector
  if (a_comp >= m_convergenceMetrics.size())
  {
    Vector<Real> tempVect = m_convergenceMetrics;
    m_convergenceMetrics.resize(a_comp+1);
    for (int i=0; i<tempVect.size(); i++)
      {
        m_convergenceMetrics[i] = tempVect[i];
      }
    for (int i=tempVect.size(); i<a_comp; i++)
      {
        m_convergenceMetrics[i] = 1.0;
      }
  } // end if we need to resize vector

  m_convergenceMetrics[a_comp] = a_metric;

  // pass this on through to the LevelSolver
  if (m_levelSolver.isDefined())
    {
      m_levelSolver.setConvergenceMetric(a_metric, a_comp);
    }
}

/// sets norm type
void
AMRSolver::setNormType(int a_normType)
{
  m_normType = a_normType;
  m_levelSolver.setNormType(a_normType);
}

/// Solves on hierarchy to tolerance m_tolerance
void AMRSolver::solveAMR(Vector<LevelData<FArrayBox> *>&       a_phiLevel,
                         const Vector<LevelData<FArrayBox> *>& a_rhsLevel,
                         bool a_initializePhiToZero)
{
  CH_TIME("AMRSolver::solveAMR");
  CH_assert(isDefined());
  CH_assert(a_phiLevel.size() > m_finestLevel);
  CH_assert(a_rhsLevel.size() > m_finestLevel);

  int ncomp = a_phiLevel[m_lBase]->nComp();
  if (ncomp > m_ncomp)
    {
      MayDay::Warning("AMRSolver::solveAMR -- phi.nComp > solver.ncomp");
      ncomp = m_ncomp;
    }

  Vector<LevelData<FArrayBox> *> residVect(m_finestLevel+1,NULL);
  Vector<LevelData<FArrayBox> *> corrVect (m_finestLevel+1,NULL);

  // allocate residual and correction arrays
  // initial guess at the solution is zero
  for (int ilev = m_lBase; ilev <= m_finestLevel; ilev++)
    {
      const DisjointBoxLayout& levelGrids = a_phiLevel[ilev]->getBoxes();

      residVect[ilev] = new LevelData<FArrayBox>(levelGrids, ncomp);
      corrVect[ilev]  = new LevelData<FArrayBox>(levelGrids, ncomp, IntVect::Unit);

      LevelData<FArrayBox>& phiLev = *a_phiLevel[ilev];
      LevelData<FArrayBox>& corrLev = *corrVect[ilev];

      DataIterator ditLev = phiLev.dataIterator();
      for (ditLev.reset(); ditLev.ok(); ++ditLev)
        {
          if (a_initializePhiToZero)
            {
              phiLev [ditLev()].setVal(0.);
            }
          corrLev[ditLev()].setVal(0.0);
        }

    }

  // if m_lBase > 0, need to also create lBase-1 level correction for
  // cf BC's (set to 0)
  if (m_lBase > 0)
    {
      const DisjointBoxLayout& crseGrids = a_phiLevel[m_lBase-1]->getBoxes();

      corrVect[m_lBase-1] = new LevelData<FArrayBox>(crseGrids, ncomp, IntVect::Unit);

      // set crse correction to 0 as well
      LevelData<FArrayBox>& crseCorr = *corrVect[m_lBase-1];

      DataIterator crseDit = crseCorr.dataIterator();
      for (crseDit.begin(); crseDit.ok(); ++crseDit)
        {
          crseCorr[crseDit()].setVal(0.0);
        }
    }

  // compute initial residual
  Vector<Real> currentRes = computeResidualNorm(residVect, a_phiLevel,
                                                a_rhsLevel, m_normType);

  // if convergence metrics haven't been set, set them here
  if (m_convergenceMetrics.size() == 0)
    {
      Vector<Real> normRHS = computeNorm(a_rhsLevel, m_normType);

      // if normRHS = 0, replace with initial residual...
      for (int comp=0; comp<normRHS.size(); comp++)
        {
          if (normRHS[comp] == 0.0) normRHS[comp] = currentRes[comp];
        }

      m_convergenceMetrics = normRHS;
      // also should pass this on through to levelsolver
      if (m_levelSolver.isDefined())
        {
          for (int comp=0; comp<m_convergenceMetrics.size(); comp++)
            {
              m_levelSolver.setConvergenceMetric(m_convergenceMetrics[comp],
                                                 comp);
            }
        }
    }

  Vector<Real> oldRes = currentRes;

  if (m_verbose)
    {
      pout() << "AMRSolver, lBase = " << m_lBase  << endl;
      for (int comp = 0; comp < ncomp; comp++)
        {
          std::ios::fmtflags origFlags = pout().flags();
          int origPrecision = pout().precision();
          pout() << "initial max(residual[" << comp << "]) = "
                 << setiosflags(ios::showpoint)
                 << setiosflags(ios::scientific)
                 << setprecision(12)
                 << currentRes[comp]
                 << endl;
          pout().flags(origFlags);
          pout().precision(origPrecision);
        }
    }

  // if initial Residual = 0, don't bother solving...
  bool done = true;
  for (int comp = 0; comp < ncomp; comp++)
    {
      if (currentRes[comp] != 0) done=false;
    }

  if (done)
    {
      if (m_verbose)
        {
          pout() << "AMRSolver::solveAMR--initial residual == 0\n returning... \n";
        }
    }

  int iter = 0;

  while (!done)
    {
      AMRVCycleMG(corrVect, residVect);

      // this is a clumsy way to do this (increment
      // solution with correction, reset corr to 0,
      // then recompute residual), but it should at least
      // do the right thing
      for (int ilev = m_lBase; ilev <= m_finestLevel; ilev++)
        {
          LevelData<FArrayBox>& levelCorr = *corrVect[ilev];
          LevelData<FArrayBox>& levelPhi = *a_phiLevel[ilev];

          DataIterator levelDit = levelPhi.dataIterator();
          for (levelDit.begin(); levelDit.ok(); ++levelDit)
            {
              levelPhi[levelDit()].plus(levelCorr[levelDit()],
                                        levelCorr[levelDit()].box(),
                                        0, 0, ncomp);
              levelCorr[levelDit()].setVal(0.0);
            }
        }

      oldRes = currentRes;
      currentRes = computeResidualNorm(residVect, a_phiLevel, a_rhsLevel,
                                       m_normType);

      iter++;

      done = true;

      for (int comp = 0; comp < ncomp; comp++)
        {
          // DFM (3/24/2000) second test is in case where convergence
          // hangs due to machine precision (or solvability) issues...
          //done = done &&
          //((iter > m_maxIter)
          //|| ((iter > m_minIter) && (currentRes[comp] >
          // oldRes[comp]*(1.0-m_operatorTolerance)))             || (currentRes[comp] <= m_tolerance*initRes[comp]));

          // DFM (1/16/04) remove hang test, since it's not clear we
          // really want to do it this way
          done = done &&
            ((iter > m_maxIter)
             || (currentRes[comp] <= m_tolerance*m_convergenceMetrics[comp]));

          if (m_verbose)
            {
              std::ios::fmtflags origFlags = pout().flags();
              int origPrecision = pout().precision();
              Real rate;
              if (currentRes[comp] != 0.0) rate = oldRes[comp]/currentRes[comp]; else rate = 0;
              pout() << " AMRSolver iteration #  " << iter
                     << ": Max(res[" << comp << "]) = "
                     << setiosflags(ios::showpoint)
                     << setiosflags(ios::scientific)
                     << setprecision(12)
                     << currentRes[comp]
                     << "   rate: "
                     << setiosflags(ios::fixed)
                     << setprecision(3)
                     << rate
                     << endl;
              pout().flags(origFlags);
              pout().precision(origPrecision);
            }
        }
    } // end while (!done)

  for (int comp = 0; comp < ncomp; comp++)
    {
      if (currentRes[comp] <= m_tolerance*m_convergenceMetrics[comp])
        {
          if (m_verbose)
            {
              pout() << "AMRSolver: " << iter << " iterations, final max(res["
                     << comp << "]) = "
                     << currentRes[comp] << endl;
            }
        }
      else if (iter > m_maxIter)
        {
          if (m_verbose)
            {
              pout() << "AMRSolver: reached maximum number of "
                     << iter << " iterations, final max(res)[" << comp
                     << "] = " << currentRes[comp] << endl;
            }
        }
      else if (currentRes[comp] > oldRes[comp]*(1.0-m_operatorTolerance))
        {
          if (m_verbose)
            {
              pout() << "AMRSolver: reached solver hang point; " << iter
                     << " iterations, final max(res[" << comp
                     << "]) = " << currentRes[comp] << endl;
            }
        }
      else
        {
          if (m_verbose)
            {
              pout() << "AMRSolver NOT CONVERGED! - final max(res) = "
                     << currentRes[comp] << " after "
                     << iter << " iterations" << endl;
            }
        }
    }

  // clean up storage
  for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
      if (residVect[ilev] != NULL)
        {
          delete residVect[ilev];
          residVect[ilev] = NULL;
        }

      if (corrVect[ilev] != NULL)
        {
          delete corrVect[ilev];
          corrVect[ilev] = NULL;
        }
    }
}

/**
     set whether the solver does i/o.  default is true
   */
void
AMRSolver::setVerbose(bool a_verbose)
{
  CH_assert(isDefined());

  m_verbose = a_verbose;
}

/// Does one relaxation V-cycle using a MG solver
/// implicit assumption that we're in residual-correction form
/// uses homogeneous physical boundary conditions, but inhomogenous
/// CF interp
void AMRSolver::AMRVCycleMG(Vector<LevelData<FArrayBox> *>&       a_corrLevel,
                            const Vector<LevelData<FArrayBox> *>& a_residLevel)
{
  CH_TIME("AMRSolver::AMRVCycleMG");
  CH_assert(isDefined());
  CH_assert(a_corrLevel.size() > m_finestLevel);
  CH_assert(a_residLevel.size() > m_finestLevel);

  // this is a kluge to make this work in residual-correction form
  // copy a_residLevel into each level's m_resid
  Interval compInterval = a_residLevel[m_lBase]->interval();

  for (int ilev= m_lBase; ilev <= m_finestLevel; ilev++)
    {
      LevelData<FArrayBox>& levelResid = m_amrmgLevel[ilev]->m_resid;

      a_residLevel[ilev]->copyTo(compInterval, levelResid, compInterval);
    }

  // sweep down vcycle
  for (int ilev = m_finestLevel; ilev > m_lBase; ilev--)
    {
      m_amrmgLevel[ilev]->downSweep(a_corrLevel,a_residLevel);
    }

  // solve at level lBase
  LevelData<FArrayBox> & bottomCorr = m_amrmgLevel[m_lBase]->m_corr;
  const LevelData<FArrayBox> & bottomRes = m_amrmgLevel[m_lBase]->m_resid;

  // this only does m_numVCyclesBottom iterations
  if (m_finestLevel == m_lBase)
    {
      // if only one level is being solved, it makes sense to
      // do bottom smoothing on resid/phi directly
      m_levelSolver.levelSolveH(*a_corrLevel[m_lBase],
                                *a_residLevel[m_lBase]);
    }
  else
    {
      m_levelSolver.levelSolveH(bottomCorr, bottomRes);

      // add correction to lBase's phi
      LevelData<FArrayBox> & bottomPhi = *a_corrLevel[m_lBase];

      DataIterator dit = bottomPhi.dataIterator();
      for (dit.reset(); dit.ok(); ++dit)
        {
          bottomPhi[dit()] += bottomCorr[dit()];
        }
    }

  // sweep back up vcycle
  for (int ilev = m_lBase+1; ilev <= m_finestLevel; ilev++)
    {
      m_amrmgLevel[ilev]->upSweep(a_corrLevel,a_residLevel);
    }
}

/**
    Calculate norm of multilevel residual on levels lBase to lmax.
    Does not include data covered by finer levels.
    */
Vector<Real> AMRSolver::computeResidualNorm(
                            Vector<LevelData<FArrayBox> *>&       a_phiLevel,
                            const Vector<LevelData<FArrayBox> *>& a_rhsLevel,
                            int                                   a_normType)
{
  CH_TIME("AMRSolver::computeResidualNorm_2");
  CH_assert(isDefined());
  CH_assert(a_phiLevel.size() > m_finestLevel);
  CH_assert(a_rhsLevel.size() > m_finestLevel);
  CH_assert((a_normType >= 0) && (a_normType <= 2));

  int ncomp = a_phiLevel[m_lBase]->nComp();

  if (ncomp > m_ncomp)
    {
      ncomp = m_ncomp;
    }

  Vector<Real> normTot(ncomp,0);
  Vector<Real> normLevel(ncomp,0.0);

  for (int ilev = m_finestLevel; ilev >= m_lBase; ilev--)
    {
      m_amrmgLevel[ilev]->computeAMRResidual(a_phiLevel, a_rhsLevel);

      normLevel = m_amrmgLevel[ilev]->computeResidualNorm(a_normType);

      for (int comp=0; comp<ncomp; comp++)
        {
          if (a_normType == 0)
            {
              normTot[comp] = Max(normTot[comp], normLevel[comp]);
            }
          else
            {
              normTot[comp] += normLevel[comp];
            }
        }
    } // end loop over levels

  if (a_normType == 2)
    {
      for (int comp=0; comp<ncomp; comp++)
        {
          normTot[comp] = sqrt(normTot[comp]);
        }
    }

  return normTot;
}

/**
    Calculate norm of multilevel residual on levels lBase to lmax.
    Does not include data covered by finer levels.  Also returns
    residual in argument a_residLevel
    */
Vector<Real> AMRSolver::computeResidualNorm(
                            Vector<LevelData<FArrayBox> *>&       a_residLevel,
                            Vector<LevelData<FArrayBox> *>&       a_phiLevel,
                            const Vector<LevelData<FArrayBox> *>& a_rhsLevel,
                            int                                   a_normType)
{
  CH_TIME("AMRSolver::computeResidualNorm_1");
  CH_assert(isDefined());
  CH_assert(a_phiLevel.size() > m_finestLevel);
  CH_assert(a_rhsLevel.size() > m_finestLevel);
  CH_assert((a_normType >= 0) && (a_normType <= 2));

  int ncomp = a_phiLevel[m_lBase]->nComp();

  if (ncomp > m_ncomp)
    {
      ncomp = m_ncomp;
    }

  Vector<Real> normTot(ncomp,0);
  Vector<Real> normLevel(ncomp,0.0);

  for (int ilev = m_finestLevel; ilev >= m_lBase; ilev--)
    {
      LevelData<FArrayBox>& levelResid = *a_residLevel[ilev];
      m_amrmgLevel[ilev]->computeAMRResidual(levelResid, a_phiLevel, a_rhsLevel);

      normLevel = m_amrmgLevel[ilev]->computeNorm(levelResid,a_normType);

      for (int comp=0; comp<ncomp; comp++)
        {
          if (a_normType == 0)
            {
              normTot[comp] = Max(normTot[comp], normLevel[comp]);
            }
          else
            {
              normTot[comp] += normLevel[comp];
            }
        }
    } // end loop over levels

  if (a_normType == 2)
    {
      for (int comp=0; comp<ncomp; comp++)
        {
          normTot[comp] = sqrt(normTot[comp]);
        }
    }

  return normTot;
}

/**
   Calculate multilevel residual on level ilev.
    */
void AMRSolver::computeAMRResidual(
                    Vector<LevelData<FArrayBox> *>&       a_phiLevel,
                    const Vector<LevelData<FArrayBox> *>& a_rhsLevel,
                    LevelData<FArrayBox>&                 a_res,
                    int                                   a_ilev)
{
  CH_TIME("AMRSolver::computeAMRResidual");
  CH_assert(isDefined());
  CH_assert(a_ilev <= m_finestLevel);
  CH_assert(a_ilev >= 0);
  CH_assert(a_phiLevel.size() > m_finestLevel);
  CH_assert(a_rhsLevel.size() > m_finestLevel);
  CH_assert(a_rhsLevel[a_ilev]->getBoxes() == m_gridsLevel[a_ilev]);
  CH_assert(a_phiLevel[a_ilev]->getBoxes() == m_gridsLevel[a_ilev]);
  CH_assert(a_res.getBoxes() == m_gridsLevel[a_ilev]);

  // compute residual internal to amrmgLevel
  m_amrmgLevel[a_ilev]->computeAMRResidual(a_phiLevel,a_rhsLevel);

  const LevelData<FArrayBox> & amrmgRes = m_amrmgLevel[a_ilev]->m_resid;

  CH_assert(amrmgRes.getBoxes() == m_gridsLevel[a_ilev]);

  // copy residual into output array
  DataIterator dit = amrmgRes.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
    {
      a_res[dit()].copy(amrmgRes[dit()]);
    }
}

/**
    Calculate multilevel L(phi).  includes refluxing and all that
    This is the three-level operator.
*/
void AMRSolver::applyAMROperator(Vector<LevelData<FArrayBox> *>& a_phiLevel,
                                 LevelData<FArrayBox>&           a_LofPhi,
                                 int                             a_ilev)
{
  CH_assert(isDefined());
  CH_assert(a_ilev <= m_finestLevel);
  CH_assert(a_ilev >= 0);
  CH_assert(a_phiLevel.size() > m_finestLevel);

  m_amrmgLevel[a_ilev]->applyAMROperator(a_phiLevel,a_LofPhi);
}

/**
   Calculate multilevel L(phi) with homogeneous physical boundary conditions
   (but with inhomogeneous C/F BCs)
   This is the three level operator.
*/
void AMRSolver::applyAMROperatorHphys(
                    Vector<LevelData<FArrayBox>* >& a_phiLevel,
                    LevelData<FArrayBox>&           a_LofPhi,
                    int                             a_ilev)
{
  CH_assert(isDefined());
  CH_assert(a_ilev <= m_finestLevel);
  CH_assert(a_ilev >= 0);
  CH_assert(a_phiLevel.size() > m_finestLevel);

  m_amrmgLevel[a_ilev]->applyAMROperatorHphys(a_phiLevel,a_LofPhi);
}

/* computes multilevel norm on valid regions over levels lbase->lmax for all components
   in a_phiVect.
 */
Vector<Real>
AMRSolver::computeNorm(const Vector<LevelData<FArrayBox>* >& a_phiVect,
                       int a_normType)
{
  int ncomp = a_phiVect[m_lBase]->nComp();
  Vector<Real> normTot(ncomp, 0);
  Vector<Real> normLevel(ncomp,0);

  for (int ilev = m_lBase; ilev<=m_finestLevel; ilev++)
    {
      LevelData<FArrayBox>& levelPhi = *a_phiVect[ilev];
      normLevel = m_amrmgLevel[ilev]->computeNorm(levelPhi, a_normType);

      for (int comp=0; comp<ncomp; comp++)
        {
          if (a_normType == 0)
            {
              normTot[comp] = Max(normTot[comp], normLevel[comp]);
            }
          else
            {
              normTot[comp] += normLevel[comp];
            }
        } // end loop over components
    } // end loop over levels

  if (a_normType == 2)
  {
    for (int comp=0; comp<ncomp; comp++)
      {
        normTot[comp] = sqrt(normTot[comp]);
      }
  }

  return normTot;

}

// accessor for the residual in AMRMGLevel
LevelData<FArrayBox>&
AMRSolver::getResid(int a_level)
{
  CH_assert (a_level>= 0);
  CH_assert (a_level < m_amrmgLevel.size());

  LevelData<FArrayBox>& levelResid = m_amrmgLevel[a_level]->m_resid;

  return levelResid;

}

// accessor for the residual in AMRMGLevel
LevelData<FArrayBox>&
AMRSolver::getCorr(int a_level)
{
  CH_assert (a_level>= 0);
  CH_assert (a_level < m_amrmgLevel.size());

  LevelData<FArrayBox>& levelCorr = m_amrmgLevel[a_level]->m_corr;

  return levelCorr;

}
#include "NamespaceFooter.H"
