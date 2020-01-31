#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cmath>

#include "TGASolver.H"
#include "FluxBox.H"
#include "HOExtrapBC.H"
#include "ExtrapBC.H"
#include "Extrap0thOrderBC.H"
#include "timeInterp.H"
#include "NamespaceHeader.H"

// small parameter used in defining TGA constants for viscous solves
#define TGA_EPS 1.0e-8


TGASolver::TGASolver()
{
  m_solverTol = 1.0e-7;
  m_helmop = NULL;
  m_noAlphaOp = NULL;
  m_noBetaOp = NULL;
}

TGASolver::TGASolver(const DisjointBoxLayout& a_levelGrids,
                     const DisjointBoxLayout* a_crseGrids,
                     const ProblemDomain& a_domain,
                     Real a_dxLevel,
                     int a_nRefCrse,
                     const BaseHelmholtzOp* a_helmop,
                     int a_ncomp)
{
  m_solverTol = 1.0e-7;
  m_helmop = NULL;
  m_noAlphaOp = NULL;
  m_noBetaOp = NULL;

  define(a_levelGrids, a_crseGrids, a_domain, a_dxLevel,
         a_nRefCrse, a_helmop, a_ncomp);
}

TGASolver::~TGASolver()
{
  if (m_helmop != NULL)
    {
      delete m_helmop;
      m_helmop = NULL;
    }

  if (m_noAlphaOp != NULL)
    {
      delete m_noAlphaOp;
      m_noAlphaOp = NULL;
    }

  if (m_noBetaOp != NULL)
    {
      delete m_noBetaOp;
      m_noBetaOp = NULL;
    }
}

void
TGASolver::define(const DisjointBoxLayout& a_levelGrids,
                  const DisjointBoxLayout* a_crseGrids,
                  const ProblemDomain& a_domain,
                  Real a_dxLevel,
                  int a_nRefCrse,
                  const BaseHelmholtzOp* a_helmop,
                  int a_ncomp)
{
  CH_TIME("TGA Solver::define");
  m_dx = a_dxLevel;
  m_domain = a_domain;
  if (m_helmop != NULL)
    {
      delete m_helmop;
      m_helmop = NULL;
    }
  m_helmop = (BaseHelmholtzOp*) a_helmop->new_levelop();
  if (m_noAlphaOp != NULL)
    {
      delete m_noAlphaOp;
      m_noAlphaOp = NULL;
    }
  m_noAlphaOp = (BaseHelmholtzOp*) a_helmop->new_levelop();
  if (m_noBetaOp != NULL)
    {
      delete m_noBetaOp;
      m_noBetaOp = NULL;
    }
  m_noBetaOp = (BaseHelmholtzOp*) a_helmop->new_levelop();
  // check this because if crse level doesn't exist,
  // then we can get away with homogeneousOnly in levelop define
  bool homogeneousOnly = false;
  m_helmop->define(a_levelGrids, a_crseGrids, a_dxLevel,
                   a_nRefCrse, a_domain, homogeneousOnly,
                   a_ncomp);

  m_noAlphaOp->define(a_levelGrids, a_crseGrids, a_dxLevel,
                      a_nRefCrse, a_domain, homogeneousOnly,
                      a_ncomp);
  m_noAlphaOp->scaleAlpha(0.0);

  m_noBetaOp->define(a_levelGrids, a_crseGrids, a_dxLevel,
                     a_nRefCrse, a_domain, homogeneousOnly,
                     a_ncomp);
  m_noBetaOp->scaleBeta(0.0);

  {
  CH_TIME("TGA Solver::define define levelSolver");
  m_levelSolver.define(a_levelGrids, a_crseGrids, a_domain,
                       a_dxLevel, a_nRefCrse, a_helmop,
                       a_ncomp);
  m_levelSolver.setTolerance(m_solverTol);
  }
}

void
TGASolver::updateSoln(LevelData<FArrayBox>& a_phiNew,
                      LevelData<FArrayBox>& a_phiOld,
                      LevelFluxRegister* a_FineFluxRegPtr,
                      LevelFluxRegister* a_CrseFluxRegPtr,
                      const LevelData<FArrayBox>* a_crsePhiOldPtr,
                      Real a_crseOldTime,
                      const LevelData<FArrayBox>* a_crsePhiNewPtr,
                      Real a_crseNewTime,
                      // note that a_src can't be const because of BCs
                      LevelData<FArrayBox>& a_src,
                      Real a_oldTime,
                      Real a_dt)
{
  // first compute parameters needed for TGA solves
  Real a = 2.0 - sqrt(2.0) - TGA_EPS;
  Real discr = sqrt(a*a - 4.0*a + 2.0);
  Real r1 = (2.0*a - 1.0)/(a + discr);
  Real r2 = (2.0*a - 1.0)/(a - discr);

  int ncomp = a_phiOld.nComp();
  const DisjointBoxLayout levelGrids = a_phiNew.getBoxes();
  // probably easiest if we just keep this dataIterator around
  DataIterator dit = a_phiNew.dataIterator();

  LevelData<FArrayBox> tempStorage(levelGrids, ncomp,
                                   IntVect::Unit);
  LevelData<FArrayBox> intermediateSoln(levelGrids, ncomp,
                                        IntVect::Unit);
  LevelData<FArrayBox> diffusiveSource(levelGrids, ncomp);

  LevelData<FluxBox> diffusiveFlux(levelGrids, ncomp);

  LevelData<FArrayBox>* crseDataPtr = NULL;
  if (a_crsePhiOldPtr != NULL)
    {
      const DisjointBoxLayout& crseGrids = a_crsePhiOldPtr->getBoxes();
      crseDataPtr = new LevelData<FArrayBox>(crseGrids, ncomp);
    }



  // set bogus values for debugging purposes
#ifdef NDEBUG
  for (dit.reset(); dit.ok(); ++dit)
    {
      tempStorage[dit()].setVal(666.666);
      intermediateSoln[dit()].setVal(666.666);
      diffusiveSource[dit()].setVal(666.666);
      diffusiveFlux[dit()].setVal(666.666);
    }
#endif

  // first need to compute diffused source term.  Want to use
  // Extrapolation coarse-fine BC's for this, so the easiest
  // way to do this is to use the PhysBC functionality, to set
  // _all_ ghost cells to extrap BCs, then call an exchange to
  // go back and reset the ghost cells which are in the valid
  // domain on this level

  // only need to set C/F BCs if a coarser level exists
  if (a_crsePhiOldPtr != NULL)
    {
      // first set up DomainGhostBC object
      DomainGhostBC bogusBC;
      Interval comps = a_src.interval();

      for (int dir=0; dir<SpaceDim; dir++)
        {
          SideIterator sit;
          for (sit.begin(); sit.ok(); ++sit)
            {
              ExtrapBC thisBC(dir, sit(), comps);
              bogusBC.setBoxGhostBC(thisBC);
            }
        }

      // now loop over grids on this level
      for (dit.reset(); dit.ok(); ++dit)
        {
          // i think a better idea would be to create a CFExtrap class,
          // eventually.
          Box thisBCbox(levelGrids[dit()]);
          for (int BCdir=0; BCdir<SpaceDim; BCdir++)
            {
              // first do low direction.  idea is to grow box and intersect
              // with domain, then shrink.  if the result isn't equal to
              // the original box, then we have a physical boundary.
              // in that case, grow the BC box in this direction to move
              // the apparent boundary away from the edge of the grid box
              Box testBox(thisBCbox);
              testBox.growLo(BCdir, 1);
              testBox &= m_domain;
              testBox.growLo(BCdir, -1);
              if (testBox != thisBCbox) thisBCbox.growLo(BCdir,1);

              // now do hi side
              testBox = thisBCbox;
              testBox.growHi(BCdir,1);
              testBox &= m_domain;
              testBox.growHi(BCdir, -1);
              if (testBox != thisBCbox) thisBCbox.growHi(BCdir,1);
            } // end loop over directions

          FArrayBox& thisSrc = a_src[dit()];
          bogusBC.applyHomogeneousBCs(thisSrc, thisBCbox, m_dx);
        }
      // at this point, ALL ghost cells which are not physical BCs
      // are set by extrap.

      // reset valid cells on this level to correct values
      a_src.exchange(a_src.interval());

    }

  // now can compute modified source. first apply scaled
  // operator (note that since we've already set extrap
  // BCs, we pass in a null pointer to the operation, indicating
  // that the coarse-fine BCs need not be set.
  // also note that the operator already implicitly includes the
  // diffusion coefficient, so since we want the diffused source
  // term to be (I + (half - a)dt*nu*L), we only have to scale
  // coefficient by (half-a)*dt

  // since the operator in m_helmop may have a non-unity alpha,
  // use the version of the operator with alpha set to 0

  LevelData<FArrayBox>* nullDataPtr = NULL;
  Real scaleMult = -(0.5 - a)*a_dt;
  m_noAlphaOp->scaleBeta(scaleMult);
  m_noAlphaOp->applyOpI(a_src, nullDataPtr, diffusiveSource);
  // now add in the a_src
  for (dit.begin(); dit.ok(); ++dit)
    {
      diffusiveSource[dit()] += a_src[dit()];
    }

  // increment flux
  FArrayBox tempFlux;
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& thisSrc = a_src[dit()];
      FluxBox& thisFlux = diffusiveFlux[dit()];
      for (int dir=0; dir<SpaceDim; ++dir)
        {
          m_noAlphaOp->getFlux(tempFlux, thisSrc, dit(), dir);
          // this one gets negated because it's on the "other side"
          // of the equation
          tempFlux.negate();
          tempFlux *= a_dt;
          thisFlux[dir].copy(tempFlux);
        }
    }

  // return operator to original scaling
  m_noAlphaOp->scaleBeta(1.0/scaleMult);

  // first need to do explicit part --
  // RHS for intermediate solve is
  // (alpha + (1-a)*dt*nu*L)phi + dt*f, (note opposite
  // sign) where f is
  // what is now stored in the diffusiveSrc

  // first need to compute coarse-level BC at this level's old time
  if (crseDataPtr != NULL)
    {
      timeInterp(*crseDataPtr, a_oldTime,
                 *a_crsePhiOldPtr, a_crseOldTime,
                 *a_crsePhiNewPtr, a_crseNewTime,
                 crseDataPtr->interval());

    }


  scaleMult = -(1.0-a)*a_dt;
  m_helmop->scaleBeta(scaleMult);
  m_helmop->applyOpI(a_phiOld, crseDataPtr, tempStorage);
  // increment flux
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& thisOldPhi = a_phiOld[dit()];
      FluxBox& thisFlux = diffusiveFlux[dit()];
      for (int dir=0; dir<SpaceDim; ++dir)
        {
          m_helmop->getFlux(tempFlux, thisOldPhi, dit(), dir);
          // this one also gets negated because it's on the right-hand
          // side of the equation
          tempFlux.negate();
          thisFlux[dir] += tempFlux;
        }
    }

  // return m_helmop to original scaling
  m_helmop->scaleBeta(1.0/scaleMult);

  // now compute RHS for intermediate solve
  for (dit.begin(); dit.ok(); ++dit)
    {
      diffusiveSource[dit()] *= a_dt;
      diffusiveSource[dit()] += tempStorage[dit()];
    }



  // now construct coarse-level BC for intermediate solve
  // intermediate solution will be at time = oldTime + (1-r1)dt
  if (crseDataPtr != NULL)
    {
      Real intermediateTime = a_oldTime + (1-r1)*a_dt;

      timeInterp(*crseDataPtr, intermediateTime,
                 *a_crsePhiOldPtr, a_crseOldTime,
                 *a_crsePhiNewPtr, a_crseNewTime,
                 crseDataPtr->interval());

    }

  // intermediate solve is (alpha - r2*dt*nu*L)phi = RHS
  // where RHS is what is now in diffusiveSrc
  scaleMult = r2*a_dt;
  m_levelSolver.scaleBeta(scaleMult);

  // a good first guess at the solution is the old time value
  a_phiOld.copyTo(a_phiOld.interval(), intermediateSoln,
                  intermediateSoln.interval());

  // in this case, set convergence metrics in levelSolver to be
  // norm of old phi. take norm over all components at once.
  Interval thisInterval(0,a_phiOld.nComp()-1);
  int normType = 0;
  Real normPhi = norm(a_phiOld, thisInterval, normType);
  if (normPhi == 0)
    {
      normPhi = norm(diffusiveSource, thisInterval, normType);
    }

  for (int comp=0; comp<a_phiOld.nComp(); comp++)
    {
      m_levelSolver.setConvergenceMetric(normPhi, comp);
    }

  bool initializeSolnToZero = false;
  m_levelSolver.levelSolve(intermediateSoln, crseDataPtr,
                           diffusiveSource, initializeSolnToZero);
  // return levelSolver to original scaling
  m_levelSolver.scaleBeta(1.0/scaleMult);

  // increment flux with piece from intermediate solution
  m_helmop->scaleBeta(scaleMult);
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& thisIntermediateSoln = intermediateSoln[dit()];
      FluxBox& thisFlux = diffusiveFlux[dit()];
      for (int dir=0; dir<SpaceDim; ++dir)
        {
          m_helmop->getFlux(tempFlux, thisIntermediateSoln, dit(), dir);
          thisFlux[dir] += tempFlux;
        }
    }
  // reset coefficient of helmOp
  m_helmop->scaleBeta(1.0/scaleMult);

  // now construct CF-BC for final solve

  if (crseDataPtr != NULL)
    {
      Real newTime = a_oldTime + a_dt;
      timeInterp(*crseDataPtr, newTime,
                 *a_crsePhiOldPtr, a_crseOldTime,
                 *a_crsePhiNewPtr, a_crseNewTime,
                 crseDataPtr->interval());
    }


  // final solve is (I - r1*dt*nu*L)phiNew = alpha*intermediateSoln
  // multiply intermediate solution by alpha using noBetaOp
  // (boundary conditions are irrelevant, so ensure no C/F BC's are
  // done by using nullDataPtr as crse-level data.
  // we are re-using diffusiveSrc to hold the RHS
  m_noBetaOp->applyOpI(intermediateSoln, nullDataPtr, diffusiveSource);

  scaleMult = r1*a_dt;
  m_levelSolver.scaleBeta(scaleMult);

  for (int comp=0; comp<a_phiOld.nComp(); comp++)
    {
      m_levelSolver.setConvergenceMetric(normPhi, comp);
    }

  // good first guess at soln is intermediate soln
  intermediateSoln.copyTo(intermediateSoln.interval(),
                          a_phiNew, a_phiNew.interval());

  m_levelSolver.levelSolve(a_phiNew, crseDataPtr,
                           diffusiveSource, initializeSolnToZero);

  // rescale levelsolver, just in case we want to reuse this solver
  m_levelSolver.scaleBeta(1.0/scaleMult);

  // increment diffusive flux with final piece (from new soln)
  m_helmop->scaleBeta(scaleMult);
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& thisPhi = a_phiNew[dit()];
      FluxBox& thisFlux = diffusiveFlux[dit()];
      for (int dir=0; dir<SpaceDim; ++dir)
        {
          m_helmop->getFlux(tempFlux, thisPhi, dit(), dir);
          thisFlux[dir] += tempFlux;
        }
    }
  // reset scale for helmop
  m_helmop->scaleBeta(1.0/scaleMult);

  // now increment flux registers -- note that because of the way
  // we defined the fluxes, the dt multiplier is already in the
  // flux
  if (a_FineFluxRegPtr != NULL)
    {
      Real fluxMult = 1.0;
      Interval comps = diffusiveFlux.interval();
      for (dit.begin(); dit.ok(); ++dit)
        {
          FluxBox& thisFlux = diffusiveFlux[dit()];
          for (int dir=0; dir<SpaceDim; ++dir)
            {
              a_FineFluxRegPtr->incrementCoarse(thisFlux[dir],
                                                fluxMult, dit(), comps,
                                                comps, dir);
            }
        }
    } // end if there is a finer-level

  if (a_CrseFluxRegPtr != NULL)
    {
      Real fluxMult = 1.0;
      Interval comps = diffusiveFlux.interval();

      for (dit.begin(); dit.ok(); ++dit)
        {
          FluxBox& thisFlux = diffusiveFlux[dit()];
          for (int dir=0; dir<SpaceDim; ++dir)
            {
              a_CrseFluxRegPtr->incrementFine(thisFlux[dir],
                                              fluxMult, dit(), comps,
                                              comps, dir, Side::Lo);


              a_CrseFluxRegPtr->incrementFine(thisFlux[dir],
                                              fluxMult, dit(), comps,
                                              comps, dir, Side::Hi);
            }
        }
    } // end if there is a coarser level




  // clean up storage, and we should be done...
  if (crseDataPtr != NULL)
    {
      delete crseDataPtr;
      crseDataPtr = NULL;
    }

}


// in this function, we call the updateSoln function, compute
// the update, then subtract off the extra pieces to return the
// diffusive part of the update
void
TGASolver::computeDiffusion(LevelData<FArrayBox>& a_DiffusiveTerm,
                            LevelData<FArrayBox>& a_phiOld,
                            LevelFluxRegister* a_FineFluxRegPtr,
                            LevelFluxRegister* a_crseFluxRegPtr,
                            const LevelData<FArrayBox>* a_crsePhiOldPtr,
                            Real a_crseOldTime,
                            const LevelData<FArrayBox>* a_crsePhiNewPtr,
                            Real a_crseNewTime,
                            LevelData<FArrayBox>& a_src,
                            Real a_oldTime,
                            Real a_dt)
{
  CH_TIME("TGA Solver::computeDiffusion");
  // first compute updated solution
  LevelData<FArrayBox> tempSoln(a_phiOld.getBoxes(),
                                a_phiOld.nComp(),
                                IntVect::Unit);

  DataIterator tsdit = tempSoln.dataIterator();
  for (tsdit.begin(); tsdit.ok(); ++tsdit)
    {
      tempSoln[tsdit()].setVal(0.0);
    }

  {
  CH_TIME("TGA Solver::computeDiffusion updateSoln");
  updateSoln(tempSoln, a_phiOld, a_FineFluxRegPtr,
             a_crseFluxRegPtr, a_crsePhiOldPtr, a_crseOldTime,
             a_crsePhiNewPtr, a_crseNewTime, a_src, a_oldTime, a_dt);
  }

  // now subtract everything off to leave us with diffusive term
  DataIterator dit = a_DiffusiveTerm.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      // first subtract off old-time state to leave us with dphi
      tempSoln[dit()] -= a_phiOld[dit()];
      // now divide by dt to get dphi/dt
      tempSoln[dit()].divide(a_dt);
    }

  // now pass through the no-beta operator in case there is an alpha coeff
  // this will multiply by alpha, if there is one
  m_noBetaOp->applyOpH(tempSoln, a_DiffusiveTerm);

  for (dit.begin(); dit.ok(); ++dit)
    {
      // and finally, subtract off a_src
      a_DiffusiveTerm[dit()].minus(a_src[dit()]);
    }

  // what's left should be the time-centered diffusive part of the update


}


BaseHeatSolver*
TGASolver::new_heatSolver() const
{
  TGASolver* new_ptr = new TGASolver();
  if (m_helmop != NULL)
    {
      new_ptr->m_helmop = static_cast<BaseHelmholtzOp*>(m_helmop->new_levelop());
    }
  else
    {
      new_ptr->m_helmop = NULL;
    }

  // copy solver tolerance
  new_ptr->m_solverTol = m_solverTol;


  return static_cast<BaseHeatSolver*>(new_ptr);
}

void
TGASolver::setSolverTolerance(Real a_solverTol)
{
  m_solverTol = a_solverTol;
  if (m_levelSolver.isDefined())
    {
      m_levelSolver.setTolerance(m_solverTol);
    }
}

Real
TGASolver::solverTolerance() const
{
  return m_solverTol;
}

#include "NamespaceFooter.H"
