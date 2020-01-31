#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// DTGraves, Fri, July 23, 1999

#include <cmath>

#include "LevelSolver.H"
#include "IntVectSet.H"
#include "parstream.H"
#include "LayoutIterator.H"
#include "NamespaceHeader.H"

using std::cout;
using std::endl;

bool LevelSolver::isDefined() const
{
  return m_isDefined;
}

void LevelSolver::setDefaultValues()
{
  m_levelOpPtr = NULL;

  // set this to an invalid value
  m_nRefCrse = -1;
  m_dxLevel = -1.0;

  m_maxIter = 33;
  m_minIter = 4;

  m_tolerance = 1.0e-10;
  m_operatorTolerance = 1.0e-4;
#ifdef CH_USE_FLOAT
  m_tolerance = sqrt(m_tolerance);
#endif

  // default is maxNorm
  m_normType = 0;

  m_verbose = false;

  m_isDefined = false;
}

void LevelSolver::setVerbose(bool a_verbose)
{
  m_verbose = a_verbose;
}

void LevelSolver::setNumSmoothUp(int a_numSmoothUp)
{
  CH_assert(isDefined());

  m_numSmoothUp = a_numSmoothUp;
  m_levelMG.setNumSmoothUp(a_numSmoothUp);
}

void LevelSolver::setNumSmoothDown(int a_numSmoothDown)
{
  CH_assert(isDefined());

  m_numSmoothDown = a_numSmoothDown;
  m_levelMG.setNumSmoothDown(a_numSmoothDown);
}

void LevelSolver::setTolerance(Real a_tolerance)
{
  m_tolerance = a_tolerance;
}

// set tolerance for stopping due to hung convergence
void LevelSolver::setOperatorTolerance(Real a_tolerance)
{
  CH_assert(a_tolerance >= 0.0);

  m_operatorTolerance = a_tolerance;
}

void LevelSolver::setMaxIter(int a_maxIter)
{
  m_maxIter = a_maxIter;
}

void LevelSolver::setMinIter(int a_minIter)
{
  m_minIter = a_minIter;
}

void
LevelSolver::setConvergenceMetric(Real a_metric, int a_comp)
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
  // pass this on through to the LevelMG (and from there to the LevelOps)
  m_levelMG.setConvergenceMetric(a_metric, a_comp);
}

void
LevelSolver::resetConvergenceMetrics()
{
  m_convergenceMetrics.clear();
  // don't reset levelMG convergence metrics,
  // since when the're re-set in the LevelSolver, the new
  // metrics will be passed through.
}

void
LevelSolver::setNormType(int a_normType)
{
  m_normType = a_normType;
}

LevelSolver::LevelSolver()
{
  setDefaultValues();
}

LevelSolver::~LevelSolver()
{
  clearMemory();
}

void LevelSolver::define(const DisjointBoxLayout& a_grids,
                         const DisjointBoxLayout* a_baseGrids,
                         const Box&               a_domain,
                         Real                     a_dxLevel,
                         int                      a_nRefCrse,
                         const LevelOp* const     a_opin,
                         int                      a_ncomp,
                         bool                     a_limitCoarsening)
{
  ProblemDomain physdomain(a_domain);

  define(a_grids, a_baseGrids, physdomain, a_dxLevel, a_nRefCrse,
         a_opin, a_ncomp, a_limitCoarsening);
}

void LevelSolver::define(const DisjointBoxLayout& a_grids,
                         const DisjointBoxLayout* a_baseGrids,
                         const ProblemDomain&     a_domain,
                         Real                     a_dxLevel,
                         int                      a_nRefCrse,
                         const LevelOp* const     a_opin,
                         int                      a_ncomp,
                         bool                     a_limitCoarsening)
{
  clearMemory();

  m_grids = a_grids;
  m_domain = a_domain;
  m_nRefCrse = a_nRefCrse;
  m_dxLevel = a_dxLevel;

  m_resid.  define(m_grids, a_ncomp, IntVect::Zero);
  m_scratch.define(m_grids, a_ncomp, IntVect::Zero);
  m_corr.   define(m_grids, a_ncomp, IntVect::Unit);

  m_isDefined = true;

  m_levelOpPtr = a_opin->new_levelop();
  m_levelOpPtr->define(m_grids, a_baseGrids, m_dxLevel,
                       m_nRefCrse, m_domain, false, a_ncomp);

  // coarsen boxes once to start, because coarsest level boxes
  // must also be coarsenable (for proper coarse-fine interpolation)
  // however, if grid merely covers entire physical domain, then
  // want to coarsen all the way down
  int nCoarserLevels = -1;

  // if limiting coarsening, only coarsen down to next AMR level
  // only relevant if there's a coarser AMR level defined
  if (a_limitCoarsening && (a_baseGrids != NULL) )
    {
      nCoarserLevels = 1;
      int temp = m_nRefCrse/2;
      // funny way to do a log base 2
      while (temp > 1)
        {
          nCoarserLevels++;
          temp /= 2;
        } // end log base 2
    }  // end if we're limiting the coarsening
  else
    {
      // begin if we're coarsening as far as possible

      // in this case, if there is a coarse-fine interface, we can
      // only coarsen down to grids with a dimension of 2, since we
      // need 2 fine cells to do coarse-fine interpolation.
      // However, if the entire domain is covered, then we can coarsen
      // as much as possible (down to a 1x1 grid), since we don't need
      // to do any coarse-fine interpolation.

      // (DFM 10/21/02) don't bother doing _any_ of this if m_grids
      // is an empty set
      if (m_grids.size() > 0)
        {
//           // why do it this way rather than as a dense IVS?
//           IntVectSet ivsCheck(m_domain.domainBox());

//           LayoutIterator lit = m_grids.layoutIterator();
//           for (lit.reset(); lit.ok(); ++lit)
//             {
//               ivsCheck -= m_grids.get(lit());
//             }

//           if (ivsCheck.isEmpty())
//             {
//               nCoarserLevels += 1;
//             }
          // There is a much cheaper way to determine if the entire domain is covered.  count points (bvs)
          long long totalPts = m_domain.domainBox().numPts();
          LayoutIterator lit = m_grids.layoutIterator();
          for (lit.reset(); lit.ok(); ++lit)
            {
              totalPts -= m_grids.get(lit()).numPts();
            }
          if (totalPts == 0)
            {
              nCoarserLevels += 1;
            }


          int refrattot = 2;

          bool getOut = false;
          while (!getOut)
            {
              bool addone = true;
              for (lit.reset(); lit.ok(); ++lit)
                {
                  Box blocal = m_grids.get(lit());

                  if (refine(coarsen(blocal,refrattot),refrattot) != blocal)
                    {
                      addone = false;
                    }

                  if (!addone)
                    {
                      break;
                    }
                }

              if (addone)
                {
                  nCoarserLevels += 1;
                  refrattot *= 2;
                }
              else
                {
                  getOut = true;
                }
            } // end while loop
        } // end if there are actually grids

    } // end if we're coarsening as far as possible

  nCoarserLevels = Max(nCoarserLevels, 0);

  m_levelMG.define(m_grids, a_baseGrids, m_dxLevel,
                   m_nRefCrse, m_domain, nCoarserLevels, a_opin, a_ncomp);
}

// complete constructor
LevelSolver::LevelSolver(const DisjointBoxLayout& a_grids,
                         const DisjointBoxLayout* a_baseGrids,
                         const Box&               a_domain,
                         Real                     a_dxLevel,
                         int                      a_nRefCrse,
                         const LevelOp* const     a_opin,
                         int                      a_ncomp,
                         bool                     a_limitCoarsening)
{
  setDefaultValues();

  ProblemDomain physdomain(a_domain);

  define(a_grids, a_baseGrids, physdomain, a_dxLevel,
         a_nRefCrse, a_opin, a_ncomp, a_limitCoarsening);
}

// complete constructor
LevelSolver::LevelSolver(const DisjointBoxLayout& a_grids,
                         const DisjointBoxLayout* a_baseGrids,
                         const ProblemDomain&     a_domain,
                         Real                     a_dxLevel,
                         int                      a_nRefCrse,
                         const LevelOp* const     a_opin,
                         int                      a_ncomp,
                         bool                     a_limitCoarsening)
{
  setDefaultValues();

  define(a_grids, a_baseGrids, a_domain, a_dxLevel,
         a_nRefCrse, a_opin, a_ncomp, a_limitCoarsening);
}

// returns LevelSolver to a basically undefined state
void LevelSolver::clearMemory()
{
  if (m_levelOpPtr != NULL)
    {
      delete m_levelOpPtr;
    }

  m_levelOpPtr = NULL;
}

// solve on just this level, inhomogeneous bcs
void LevelSolver::levelSolve(LevelData<FArrayBox>&       a_phi,
                             const LevelData<FArrayBox>* a_phic,
                             const LevelData<FArrayBox>& a_rhs,
                             bool                        a_initializePhiToZero)
{
  CH_assert(isDefined());
  CH_assert(a_phi.nComp() == a_rhs.nComp());
  CH_assert(a_phi.nComp() == m_resid.nComp());

  /* compute initial residual
     -- can't use residual() fn because
     don't want to zero out higher levels
     -- also note that we use LNf,
     because this is a LEVEL solve
  */

  CH_assert(m_levelOpPtr->isDefined());
  DataIterator dit = a_phi.dataIterator();

  if (a_initializePhiToZero)
    {
      //initial guess of phi is zero
      for (dit.reset(); dit.ok(); ++dit)
        {
          a_phi[dit()].setVal(0.0);
        }
    }

  m_levelOpPtr->applyOpI(a_phi, a_phic, m_resid);

  for (dit.reset(); dit.ok(); ++dit)
    {
      m_resid[dit()] -= a_rhs[dit()];
      m_resid[dit()].negate();

      m_corr[dit()].setVal(0.0);
    }

  Vector<Real> oldRes(m_resid.nComp());

  bool done = true;

  for (int comp = 0; comp < m_resid.nComp(); comp++)
    {
      Interval thisInterval(comp,comp);

      // unscaled norm
      oldRes[comp] = norm(m_resid, thisInterval, m_normType);

      // do scaling
      if (m_normType != 0)
        {
          Real exponent = SpaceDim;
          exponent /= m_normType;
          Real scale = pow(m_dxLevel, exponent);
          oldRes[comp] *= scale;
        }

      // if initRes = 0, don't bother solving
      if (oldRes[comp] > 0.01*m_tolerance)
        {
          done=false;
        }
    }

  if (done)
    {
      return;
    }

  // if this hasn't already been set,
  // compute norm(RHS) for use as convergence metrics
  if (m_convergenceMetrics.size() == 0)
    {
      m_convergenceMetrics = Vector<Real>(a_rhs.nComp(),0);

      for (int comp=0; comp<a_rhs.nComp(); comp++)
        {
          Interval thisInterval(comp,comp);
          // unscaled norm
          m_convergenceMetrics[comp] = norm(a_rhs, thisInterval, m_normType);

          // do scaling if necessary
          if (m_normType != 0)
            {
              Real exponent = SpaceDim;
              exponent /= m_normType;
              Real scale = pow(m_dxLevel, exponent);
              m_convergenceMetrics[comp] *= scale;
            }

          // if norm(rhs) is 0, replace with norm(residual)
          if (m_convergenceMetrics[comp] == 0.0)
            {
              m_convergenceMetrics[comp] = oldRes[comp];
            }
          m_levelMG.setConvergenceMetric(m_convergenceMetrics[comp],
                                         comp);
        }
    } // end if we are computing norm(RHS) as a convergence metric

  Real currentRes = 0.0;
  int iter = 0;
  while (!done)
    {
      m_levelMG.mgRelax(m_corr, m_resid, true);

      m_levelOpPtr->applyOpH(m_corr, m_scratch);

      for (dit.reset(); dit.ok(); ++dit)
        {
          m_scratch[dit()] -= m_resid[dit()];
          m_scratch[dit()].negate();
        }

      iter++;

      done = true;
      for (int comp=0; comp<m_resid.nComp(); comp++)
        {
          Interval thisInterval(comp,comp);

          // unscaled norm
          currentRes = norm(m_scratch, thisInterval, m_normType);

          // do scaling if necessary
          if (m_normType != 0)
            {
              Real exponent = SpaceDim;
              exponent /= m_normType;
              Real scale = pow(m_dxLevel, exponent);
              currentRes *= scale;
            }

          if (m_verbose)
            {
              pout() << "LevelSolve iteration # " << iter
                     << " Max(res) = " << currentRes
                     << endl;
            }

          // DFM (9/20/02) third test is to catch case where
          // convergence hangs due to machine precision (or
          // solvability) issues.  This follows the approach
          // used in AMRSolver
          //done = done && ((iter > m_maxIter)
          //                || (currentRes <= m_tolerance*initRes[comp] )
          //                || ((iter > m_minIter)
          //                    && (currentRes
          //                        > (1.0-m_operatorTolerance)*oldRes[comp])));

          // DFM (1/16/04) remove hang test (same as in AMRSolver),
          // since it's not clear we really want to do it this way
          done = done && ((iter > m_maxIter)
                          || (currentRes <= m_tolerance*m_convergenceMetrics[comp]));

          // save currentRes into oldRes
          oldRes[comp] = currentRes;
        } // end loop over components
    }  // end solver iterations [end while (!done)]

  // update phi
  for (dit.reset(); dit.ok(); ++dit)
    {
      a_phi[dit()] += m_corr[dit()];
    }

  // reset convergence metrics on exit to prevent re-using them
  // if the solver is being re-used
  resetConvergenceMetrics();

}

// solve on just this level, homogeneous bcs
void LevelSolver::levelSolveH(LevelData<FArrayBox>&       a_phi,
                              const LevelData<FArrayBox>& a_rhs,
                              bool                        a_initializePhiToZero)
{
  CH_assert(isDefined());
  CH_assert(a_phi.nComp() == a_rhs.nComp());
  CH_assert(a_phi.nComp() == m_resid.nComp());

  /** compute initial residual
      -- can't use residual() fn because
      don't want to zero out higher levels
      -- also note that we use LNf,
      because this is a LEVEL solve
  **/

  CH_assert(m_levelOpPtr->isDefined());

  DataIterator dit = a_phi.dataIterator();

  if (a_initializePhiToZero)
    {
      //initial guess of phi is zero
      for (dit.reset(); dit.ok(); ++dit)
        {
          a_phi[dit()].setVal(0.0);
        }
    }

  m_levelOpPtr->applyOpH(a_phi, m_resid);

  for (dit.reset(); dit.ok(); ++dit)
    {
      m_resid[dit()] -= a_rhs[dit()];
      m_resid[dit()].negate();
      m_corr[dit()].setVal(0.0);
    }

  Vector<Real> oldRes(m_resid.nComp());

  bool done = true;

  for (int comp = 0; comp < m_resid.nComp(); comp++)
    {
      Interval thisInterval(comp,comp);
      // unscaled norm...
      oldRes[comp] = norm(m_resid, thisInterval, m_normType);

      // do scaling if necessary
      if (m_normType != 0)
        {
          Real exponent = SpaceDim;
          exponent /= m_normType;
          Real scale = pow(m_dxLevel, exponent);
          oldRes[comp] *= scale;
        }

      // if initRes = 0, don't bother solving
      if (oldRes[comp] > 0.01*m_tolerance)
        {
          done = false;
        }
    }

  // if this hasn't already been set,
  // compute norm(RHS) for use as convergence metrics
  if (m_convergenceMetrics.size() == 0)
    {
      m_convergenceMetrics = Vector<Real>(a_rhs.nComp(),0);

      for (int comp=0; comp<a_rhs.nComp(); comp++)
        {
          Interval thisInterval(comp,comp);
          // unscaled norm
          m_convergenceMetrics[comp] = norm(a_rhs, thisInterval, m_normType);

          // do scaling if necessary
          if (m_normType != 0)
            {
              Real exponent = SpaceDim;
              exponent /= m_normType;
              Real scale = pow(m_dxLevel, exponent);
              m_convergenceMetrics[comp] *= scale;
            }

          // if norm(rhs) is 0, replace with norm(residual)
          if (m_convergenceMetrics[comp] == 0.0)
            {
              m_convergenceMetrics[comp] = oldRes[comp];
            }
          m_levelMG.setConvergenceMetric(oldRes[comp], comp);
        }
    } // end if we are computing norm(RHS) as a convergence metric

  if (done)
    {
      return;
    }

  Real currentRes = 0.0;
  int iter = 0;
  while (!done)
    {
      iter++;

      // compute modified residual LofPhi = res - LNf(corr)
      // use homogeneous form of ApplyOp because want to use
      // homogeneous CF BC's for correction
      m_levelMG.mgRelax(m_corr, m_resid, true);

      m_levelOpPtr->applyOpH(m_corr, m_scratch);

      for (dit.reset(); dit.ok(); ++dit)
        {
          m_scratch[dit()] -= m_resid[dit()];
          m_scratch[dit()].negate();
        }

      done = true;
      for (int comp = 0; comp < m_scratch.nComp(); comp++)
        {
          Interval thisInterval(comp,comp);
          // unscaled norm
          currentRes = norm(m_scratch, thisInterval, m_normType);

          // do scaling if necessary
          if (m_normType != 0)
            {
              Real exponent = SpaceDim;
              exponent /= m_normType;
              Real scale = pow(m_dxLevel, exponent);
              currentRes *= scale;
            }

          if (m_verbose)
            {
              pout() << "LevelSolveH iteration # " << iter
                     << ": Max(res) = " << currentRes
                     << endl;
            }

          // DFM (9/20/02) third test is to catch case where
          // convergence hangs due to machine precision (or
          // solvability) issues.  This follows the approach
          // used in AMRSolver
          //done = done && ((iter > m_maxIter)
          //                || (currentRes <= m_tolerance*initRes[comp] )
          //                    || ((iter < m_minIter)
          //                        && (currentRes
          //                           > (1.0-m_operatorTolerance)*oldRes[comp])));

          // DFM (1/16/04) remove hang test because it's not clear that we want
          // to do it this way...
          done = done && ((iter > m_maxIter)
                          || (currentRes <= m_tolerance*m_convergenceMetrics[comp] ));

          oldRes[comp] = currentRes;
        } // end loop over components
    }  // end solver iterations

  // update phi
  for (dit.reset(); dit.ok(); ++dit)
    {
      a_phi[dit()] += m_corr[dit()];
    }
}
#include "NamespaceFooter.H"
