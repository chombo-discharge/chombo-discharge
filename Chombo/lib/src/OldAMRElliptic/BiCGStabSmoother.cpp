#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// DFMartin, Sun, May 5, 2002

#include "DotProduct.H"
#include "LayoutIterator.H"
#include "BoxIterator.H"

#include "BiCGStabSmoother.H"
#include "NamespaceHeader.H"

BiCGStabSmoother::BiCGStabSmoother()
{
  m_maxIter = 40;

#ifdef CH_USE_FLOAT
  m_solverTol = 1.0e-7;
  m_small = 1.0e-30;
#endif
#ifdef CH_USE_DOUBLE
  m_solverTol = 1.0e-10;
  m_small = 1.0e-60;
#endif
  // m_convergeSmall = 1.0e-60;
  m_convergeSmall = 1.0e-2;

  m_enableRestart = true;
  m_verbose = false;
}

BiCGStabSmoother::~BiCGStabSmoother()
{
}

void BiCGStabSmoother::setMaxIter(int a_maxIter)
{
  m_maxIter = a_maxIter;
}

void BiCGStabSmoother::setConvergenceMetric(Real a_metric, int a_comp)
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
}

void BiCGStabSmoother::setSolverTol(Real a_solverTol)
{
  m_solverTol = a_solverTol;
}

void BiCGStabSmoother::setEnableRestart(bool a_enableRestart)
{
  m_enableRestart = a_enableRestart;
}

void BiCGStabSmoother::setVerbose(bool a_verbose)
{
  m_verbose = a_verbose;
}

BaseBottomSmoother* BiCGStabSmoother::new_bottomSmoother() const
{
  BiCGStabSmoother* newsmoother = new BiCGStabSmoother();
  newsmoother->setVerbose(m_verbose);

  if (newsmoother == NULL)
    {
      MayDay::Error("Out of Memory in BiCGStabSmoother::new_bottomSmoother");
    }

  return static_cast<BaseBottomSmoother*>(newsmoother);
}

/***********************/
// This does a GSRB Pre/Conditioned BiCGStab on a level
// for the bottom solver.
/***********************/
void BiCGStabSmoother::doBottomSmooth(LevelData<FArrayBox>&       a_phi,
                                      const LevelData<FArrayBox>& a_rhs,
                                      LevelOp*                    a_levelopPtr)
{
  bool alreadyRestarted = false;

  // if we won't allow resarts, then act like we already have...
  if (!m_enableRestart)
    {
      alreadyRestarted = true;
    }

  Real small2 = 1.0e-8;
  int ncomp = a_rhs.nComp();

  CH_assert(a_rhs.nComp() == a_phi.nComp());

  // in the case where the convergence metrics haven't already been
  // defined, default to norm(rhs) for the bottom smooth
  if (m_convergenceMetrics.size() != ncomp)
    {
      m_convergenceMetrics = Vector<Real>(ncomp, 0);
      for (int comp=0; comp<ncomp; comp++)
        {
          Interval thisInterval(comp,comp);
          int normType = 0;
          m_convergenceMetrics[comp] = norm(a_rhs, thisInterval, normType);
        } // end loop over components
    } // end if we need to resize convergence metrics

  const DisjointBoxLayout grids = a_rhs.getBoxes();

  CH_assert (grids == a_phi.getBoxes());

  LevelData<FArrayBox> corremf(grids, ncomp, IntVect::Unit);
  LevelData<FArrayBox> residmf(grids, ncomp, IntVect::Zero);
  LevelData<FArrayBox> rtwidmf(grids, ncomp, IntVect::Zero);
  LevelData<FArrayBox> shatmf (grids, ncomp, IntVect::Unit);
  LevelData<FArrayBox> phatmf (grids, ncomp, IntVect::Unit);
  LevelData<FArrayBox> pmf    (grids, ncomp, IntVect::Zero);
  LevelData<FArrayBox> vmf    (grids, ncomp, IntVect::Zero);
  LevelData<FArrayBox> tmf    (grids, ncomp, IntVect::Zero);

  DataIterator dit = residmf.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      // this is necessary to prevent uninitialized memory reads
      residmf[dit()].setVal(0.0);
      rtwidmf[dit()].setVal(0.0);
      shatmf[dit()] .setVal(0.0);
      phatmf[dit()] .setVal(0.0);
      pmf[dit()]    .setVal(0.0);
      tmf[dit()]    .setVal(0.0);
      corremf[dit()].setVal(0.0);
      vmf[dit()]    .setVal(0.0);
    }

  // compute initial residual
  a_levelopPtr->applyOpH(a_phi,residmf);

  for (dit.begin(); dit.ok(); ++dit)
    {
      residmf[dit()] -= a_rhs[dit()];
      residmf[dit()].negate();
    }

  Vector<Real> oldNorm (ncomp,0.0);
  Vector<Real> rnorm   (ncomp,0.0);
  Vector<Real> norms   (ncomp,0.0);

  // compute residual and rhs norms
  // finish is norm(residual)/nrom(rhs) < m_solverTol
  // for all components
  for (int comp = 0; comp < ncomp; comp++)
    {
      Interval compInterval(comp,comp);

      rnorm[comp] = norm(residmf, compInterval, 0);

      oldNorm[comp] = rnorm[comp];

      if (m_verbose)
        {
          pout() << "Bottom Smoother norm of initial residual [" << comp
                 << "] = " << rnorm[comp] << endl;
        }
    }

  for (dit.begin(); dit.ok(); ++dit)
    {
      rtwidmf[dit()].copy(residmf[dit()],0,0,ncomp);
    }

  Vector<Vector<Real> > rho(ncomp);
  Vector<Vector<Real> > alpha(ncomp);
  // Vector<Vector<Real> > beta(ncomp);
  Vector<Vector<Real> > omega(ncomp);

  // dummy omega[0], alpha[0] (not used)
  for (int comp=0; comp<ncomp; comp++)
    {
      omega[comp].push_back(0.0);
      alpha[comp].push_back(0.0);
      // beta[comp].push_back(0.0);
    }

  // not sure if this is the right thing to do.
  // option 1: if any rhoDots are bad, add in current correction and exit
  // option 2: keep going with other components?
  bool initialResidualZero = true;

  // if initial residual is already 0, don't do anything here
  for (int comp=0; comp<ncomp; comp++)
    {
      if (rnorm[comp] > m_small) initialResidualZero = false;
    }

  bool finished = initialResidualZero;
  bool badDotProd = false;

  int iter;
  for (iter = 1; (iter <= m_maxIter) && !finished; iter++)
    {
      Vector<Real> rhodot(ncomp);
      badDotProd = false;

      for (int comp = 0; comp < ncomp; comp++)
        {
          Interval compInterval(comp,comp);

          rhodot[comp] = DotProduct(residmf, rtwidmf, grids, compInterval);

          if (Abs(rhodot[comp]) < m_small)
            {
              // cerr << "bicbstab failed in procedure ";
              // cerr << "PoissonOp::bottomSmoother" << endl;
              // cerr << "at iteration " << iter;
              badDotProd = true;
            }

          rho[comp].push_back(rhodot[comp]);
        } // end loop over components

      if (badDotProd)
        {
          for (dit.begin(); dit.ok(); ++dit)
            {
              a_phi[dit()].plus(corremf[dit()],0,0,ncomp);
            }

          if (m_verbose)
            {
              pout() << "BiCGStabSmoother: BadDotProd -- returning!  iter = "
                     << iter << endl;
            }

          return;
        }

      // if we're restarting the computation, act as if it was
      // the first iteration
      if (iter == 1 || alreadyRestarted)
        {
          for (dit.begin(); dit.ok(); ++dit)
            {
              pmf[dit()].copy(residmf[dit()]);
            }
        }
      else
        {
          for (int comp = 0; comp < ncomp; comp++)
            {
              Real betaim1 = (rho[comp][iter-1] / rho[comp][iter-2])
                           * (alpha[comp][iter-1] / omega[comp][iter-1]);

              // beta[comp].push_back(betaim1);

              Real omegaim1 = omega[comp][iter-1];

              // p(i) = r(i-1) + beta(i-1)*(p(i-1) -omega(i-1)*v(i-1))

              for (dit.begin(); dit.ok(); ++dit)
                {
                  FArrayBox& pfab = pmf[dit()];
                  FArrayBox& vfab = vmf[dit()];
                  FArrayBox& rfab = residmf[dit()];

                  Box fabbox = grids.get(dit());
                  FArrayBox temp(fabbox, 1);

                  temp.copy(vfab,comp,0,1);

                  temp *= omegaim1;
                  temp.negate();
                  temp.plus(pfab,comp,0,1);
                  temp *= betaim1;
                  temp.plus(rfab,comp,0,1);

                  pfab.copy(temp,0,comp,1);
                } // end loop over grids
            } // end loop over components
        } // end if i ==  1 's else

      // solve M(phat) = p(i)
      a_levelopPtr->levelPreconditioner(phatmf, pmf);

      // v(i) = A(phat)
      a_levelopPtr->applyOpH(phatmf,vmf);

      // alpha(i) = rho(i-1)/(rtwid^T v(i))
      bool badCorrection = false;
      for (int comp = 0; comp < ncomp; comp++)
        {
          Interval compInt(comp,comp);
          Real denom = DotProduct(vmf, rtwidmf, grids,compInt);
          Real alphai;

          if (Abs(denom) > small2*Abs(rho[comp][iter-1]))
            {
              alphai = rho[comp][iter-1]/denom;
            }
          else
            {
              alphai = 0.0;
              badCorrection = true;
            }

          alpha[comp].push_back(alphai);

          // s = r(i-1) - alpha(i)*v(i)
        } // end loop over components

      if (!badCorrection)
        {
          // increment residual and solution
          for (dit.begin(); dit.ok(); ++dit)
            {
              FArrayBox& rfab    = residmf[dit()];
              FArrayBox& vfab    = vmf[dit()];
              FArrayBox& phatfab = phatmf[dit()];
              FArrayBox& corrfab = corremf[dit()];

              for (int comp = 0; comp < ncomp; comp++)
                {
                  // increment residual: resid = resid - alpha*v
                  rfab.plus(vfab,-1.0*alpha[comp][iter],comp,comp,1);

                  // increment correction: corr = corr + alpha*phat
                  corrfab.plus(phatfab,alpha[comp][iter],comp,comp,1);
                }
            } // end loop over grids
        } // end if the first correction step is OK
      else
        {
          // if correction was bad, don't increment anything,
          // but set residual to zero (forces recomputation of resid)
          for (dit.begin(); dit.ok(); ++dit)
            {
              residmf[dit()].setVal(0.0);
            }
        }

      // loop over components and compute residuals
      for (int comp = 0; comp < ncomp; comp++)
        {
          Interval compInt(comp,comp);

          norms[comp] = norm(residmf, compInt, 0);
        } // end loop over comps

      // check norm of s.
      // if small enough, step to end (forces recomputation of residual, etc)
      bool allDone = true;
      for (int comp = 0; comp < ncomp; comp++)
        {
          if (norms[comp] > m_solverTol*m_convergenceMetrics[comp])
            {
              allDone = false;
            }
        }

      // if not "done" here, do second correction step
      if (!allDone)
        {
          // solve Mshat = s
          a_levelopPtr->levelPreconditioner(shatmf, residmf);

          // t = A(shat)
          a_levelopPtr->applyOpH(shatmf,tmf);

          // w(i) = tTresid/tTt
          for (int comp=0; comp<ncomp; comp++)
            {
              Interval compInt(comp,comp);
              Real numerw = DotProduct(tmf, residmf, grids,compInt);
              Real denomw = DotProduct(tmf, tmf, grids,compInt);
              Real omegai = 0;

              if (Abs(denomw) > m_small)
                {
                  omegai = numerw/denomw;
                }

              omega[comp].push_back(omegai);
            } // end loop over comps

          // corr(i) = corr(i-1) + omega(i)*shat
          // resid(i) = resid(i) - omega(i)*t

          for (dit.begin(); dit.ok(); ++dit)
            {
              FArrayBox& corrfab = corremf[dit()];
              FArrayBox& shatfab = shatmf[dit()];
              FArrayBox& tfab    = tmf[dit()];

              for (int comp = 0; comp < ncomp; comp++)
                {
                  // increment correction: corr = corr+omega*shat
                  corrfab.plus(shatfab,omega[comp][iter],
                               comp,comp,1);

                  // increment residual: resid = resid - omega*t
                  residmf[dit()].plus(tfab,-1.0*omega[comp][iter],
                                      comp,comp,1);
                }
            } // end loop over grids

          // reset this just in case (if we made it this
          // far, then we can assume that we're not caught in a
          // loop
          // (dfm, 5/29/02) actually, don't reset this -- only
          // allow one restart.
          // alreadyRestarted = false;
        } // end second update loop
      else
        {
          // if we skipped second update step, still need to store something
          // in omega.  try using 0
          for (int comp = 0; comp < ncomp; comp++)
            {
              omega[comp].push_back(0.0);
              // beta[comp].push_back(0.0);
            }
        }

      // now check residual.
      allDone = true;
      for (int comp = 0; comp < ncomp; comp++)
        {
          Interval compInt(comp,comp);

          rnorm[comp] = norm(residmf, compInt, 0);

          if (m_verbose && (rnorm[comp] > 0))
            {
              pout() << "BottomSmoother: after " << iter
                     << " iterations, residnorm[" << comp << "] = "
                     << rnorm[comp] << endl;
            }

          if ((rnorm[comp] > m_solverTol*m_convergenceMetrics[comp]) &&
              (Abs(omega[comp][iter]) > m_small))
            {
              allDone = false;
            }
#if 0
          // commented out version with hanging test.  also, I think
          // the second && should be an ||? (DFM 4/8/04)

          // second test is to eliminate solver hanging
          allDone = allDone && ( ((rnorm[comp]
                                   < m_solverTol*m_convergenceMetrics[comp])
                                 && (Abs(omega[comp][iter]) < m_small))
                                 || (rnorm[comp]
                                     > oldNorm[comp]*(1.0-m_convergeSmall)));
#endif

          oldNorm[comp] = rnorm[comp];
        } // end loop over comps

      // if residnorm is small enough (or if correction is
      // small enough), then recompute residual from scratch
      // and check again (to avoid roundoff-drift issues).  if
      // we _still_ think we're done, then exit
      if (allDone)
        {
          // recompute residual from scratch -- first increment
          // a_phi with correction -- also reset correction to zero
          for (dit.begin(); dit.ok(); ++dit)
            {
              a_phi[dit()].plus(corremf[dit()],0,0,ncomp);
              corremf[dit()].setVal(0.0);
            }

          a_levelopPtr->applyOpH(a_phi,residmf);

          for (dit.begin(); dit.ok(); ++dit)
            {
              residmf[dit()] -= a_rhs[dit()];
              residmf[dit()] *= -1.0;
            }

          // now recompute residual norms and check for convergence
          allDone = true;
          for (int comp=0; comp<ncomp; comp++)
            {
              Interval compInterval(comp,comp);
              rnorm[comp] = norm(residmf, compInterval, 0);

              if (rnorm[comp] > m_solverTol*m_convergenceMetrics[comp])
                {
                  allDone = false;
                }
            } // end loop over components for residnorm checking

          // if we've already restarted and haven't gotten
          // anywhere, then we've stalled and should probably
          // exit.  also, if newly-computed residual still indicates
          // convergence, then we're done.
          if (alreadyRestarted || allDone)
            {
              finished = true;
            }
          else
            {
              // otherwise, go back and try again with new residual
              if (m_verbose)
                {
                  pout() << "Bottom Smoother -- restarting iterations"
                         << endl;
                }

              // also reset \tilde{r}
              for (dit.begin(); dit.ok(); ++dit)
                {
                  rtwidmf[dit()].copy(residmf[dit()],0,0,ncomp);
                }

              finished = false;
              alreadyRestarted = true;
            }

          if (finished && m_verbose)
            {
              for (int comp=0; comp<ncomp; comp++)
                {
                  pout() << "Bottom Smoother final Norm[" << comp << "] = "
                         << rnorm[comp] << endl;
                }
            }
        } // end if we're recomputing the residual
    } // end main loop over BiCGStab iterations

  bool converged = true;

  for (int comp = 0; comp < ncomp; comp++)
    {
      if (rnorm[comp] >= m_solverTol*m_convergenceMetrics[comp])
        {
          converged = false;
        }
    }

  if (badDotProd == true || initialResidualZero == true)
    {
      // if bad dot product or initial residual is zero, then don't bother
      // telling me it's not converged.
    }
  else
    {
      if (!converged && m_verbose)
        {
          pout() << "BiCGStabSmoother not Converged!  iter = " << iter << endl;

          for (int comp = 0; comp < ncomp; comp++)
            {
              pout() << "comp = " << comp
                     << " resid/metric = "
                     << rnorm[comp]/m_convergenceMetrics[comp]
                     << endl;
            }
        }
    }
}
#include "NamespaceFooter.H"
