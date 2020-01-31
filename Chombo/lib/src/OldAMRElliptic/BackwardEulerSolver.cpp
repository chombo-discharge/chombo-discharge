#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "BackwardEulerSolver.H"
#include "FluxBox.H"
#include "HOExtrapBC.H"
#include "timeInterp.H"
#include "NamespaceHeader.H"


BackwardEulerSolver::BackwardEulerSolver()
{
  m_solverTol = 1.0e-7;
  m_helmop = NULL;
  m_noBetaOp = NULL;
}

BackwardEulerSolver::BackwardEulerSolver(const DisjointBoxLayout& a_levelGrids,
                                         const DisjointBoxLayout* a_crseGrids,
                                         const ProblemDomain& a_domain,
                                         Real a_dxLevel,
                                         int a_nRefCrse,
                                         const BaseHelmholtzOp* a_helmop,
                                         int a_ncomp)
{
  m_solverTol = 1.0e-7;
  m_helmop = NULL;
  m_noBetaOp = NULL;
  define(a_levelGrids, a_crseGrids, a_domain, a_dxLevel,
         a_nRefCrse, a_helmop, a_ncomp);
}

BackwardEulerSolver::~BackwardEulerSolver()
{
  if (m_helmop != NULL)
    {
      delete m_helmop;
      m_helmop = NULL;
    }
  if (m_noBetaOp != NULL)
    {
      delete m_noBetaOp;
      m_noBetaOp = NULL;
    }

}

void
BackwardEulerSolver::define(const DisjointBoxLayout& a_levelGrids,
                            const DisjointBoxLayout* a_crseGrids,
                            const ProblemDomain& a_domain,
                            Real a_dxLevel,
                            int a_nRefCrse,
                            const BaseHelmholtzOp* a_helmop,
                            int a_ncomp)
{
  m_dx = a_dxLevel;
  m_domain = a_domain;
  if (m_helmop!=NULL) delete m_helmop;
  m_helmop = (BaseHelmholtzOp*) a_helmop->new_levelop();
  if (m_noBetaOp!=NULL) delete m_noBetaOp;
  m_noBetaOp = (BaseHelmholtzOp*) a_helmop->new_levelop();
  // check this because if crse level doesn't exist,
  // then we can get away with homogeneousOnly in levelop define
  bool homogeneousOnly = false;
  m_helmop->define(a_levelGrids, a_crseGrids, a_dxLevel,
                   a_nRefCrse, a_domain, homogeneousOnly,
                   a_ncomp);
  // don't need to worry about C/F BCs since beta will be set to 0
  homogeneousOnly = true;
  m_noBetaOp->define(a_levelGrids, a_crseGrids, a_dxLevel,
                   a_nRefCrse, a_domain, homogeneousOnly,
                   a_ncomp);
  m_noBetaOp->scaleBeta(0.0);

  m_levelSolver.define(a_levelGrids, a_crseGrids, a_domain,
                       a_dxLevel, a_nRefCrse, a_helmop,
                       a_ncomp);
  m_levelSolver.setTolerance(m_solverTol);

}

void
BackwardEulerSolver::updateSoln(LevelData<FArrayBox>& a_phiNew,
                                LevelData<FArrayBox>& a_phiOld,
                                LevelFluxRegister* a_FineFluxRegPtr,
                                LevelFluxRegister* a_CrseFluxRegPtr,
                                const LevelData<FArrayBox>* a_crsePhiOldPtr,
                                Real a_crseOldTime,
                                const LevelData<FArrayBox>* a_crsePhiNewPtr,
                                Real a_crseNewTime,
                                // note that a_src isn't const because of BCs
                                LevelData<FArrayBox>& a_src,
                                Real a_oldTime,
                                Real a_dt)
{


  int ncomp = a_phiOld.nComp();
  const DisjointBoxLayout levelGrids = a_phiNew.getBoxes();
  // probably easiest if we just keep this dataIterator around
  DataIterator dit = a_phiNew.dataIterator();

  LevelData<FArrayBox> RHS(levelGrids, ncomp);

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
      diffusiveFlux[dit()].setVal(666.666);
    }
#endif


  // first construct RHS for solve
  // RHS = phiOld + dt*Src

  // note that we can't simply do a copy because of the possibility of
  // there being a coefficient for phi -- operator is (alpha*I + beta*lap)
  // to ensure that this case works as well, use clone of original op
  // with the coefficient for beta set to 0.  doing an applyOp will then
  // give you alpha*u_old

   m_noBetaOp->applyOpH(a_phiOld, RHS);

   // now add in src
   for (dit.begin(); dit.ok(); ++dit)
     {
       RHS[dit()].plus(a_src[dit()], a_dt);
     }



   // now construct coarse-level BC for solve
   // intermediate solution will be at newTime = oldTime + dt
   if (crseDataPtr != NULL)
     {
       Real newTime = a_oldTime + a_dt;

       timeInterp(*crseDataPtr, newTime,
                  *a_crsePhiOldPtr, a_crseOldTime,
                  *a_crsePhiNewPtr, a_crseNewTime,
                  crseDataPtr->interval());

     }
   // solve is (I - dt*nu*L)phi = RHS
   Real scaleMult = a_dt;
   m_levelSolver.scaleBeta(scaleMult);

   // in this case, set convergence metrics in levelSolver to be
   // norm of old phi. take norm over all components at once.
   Interval thisInterval(0,a_phiOld.nComp()-1);
   int normType = 0;
   Real normPhi = norm(a_phiOld, thisInterval, normType);
   if (normPhi == 0)
     {
       normPhi = norm(RHS, thisInterval, normType);
     }

   for (int comp=0; comp<a_phiOld.nComp(); comp++)
     {
       m_levelSolver.setConvergenceMetric(normPhi, comp);
     }


   //a good first guess at the solution  is the old time value
   a_phiOld.copyTo(a_phiOld.interval(), a_phiNew,
                   a_phiNew.interval());

   bool initializeSolnToZero = false;
   m_levelSolver.levelSolve(a_phiNew, crseDataPtr,
                            RHS, initializeSolnToZero);
   // return levelSolver to original scaling
   m_levelSolver.scaleBeta(1.0/scaleMult);

   // increment flux with diffusive piece
   m_helmop->scaleBeta(scaleMult);
   FArrayBox tempFlux;
   for (dit.begin(); dit.ok(); ++dit)
     {
       FArrayBox& thisSoln = a_phiNew[dit()];
       FluxBox& thisFlux = diffusiveFlux[dit()];
       for (int dir=0; dir<SpaceDim; ++dir)
         {
           m_helmop->getFlux(tempFlux, thisSoln, dit(), dir);
           thisFlux[dir].copy(tempFlux);
         }
     }
   // reset coefficient of helmOp
   m_helmop->scaleBeta(1.0/scaleMult);


   // now increment flux registers -- note that because of the way
   // we defined the fluxes, the dt multiplier is already in the
   // flux
   if (a_FineFluxRegPtr != NULL)
     {
       Real fluxMult = 1.0;
       //Real fluxMult = 0.0;
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
BackwardEulerSolver::computeDiffusion(LevelData<FArrayBox>& a_DiffusiveTerm,
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

  // first compute updated solution
  LevelData<FArrayBox> tempSoln(a_phiOld.getBoxes(),
                                a_phiOld.nComp(),
                                IntVect::Unit);

  updateSoln(tempSoln, a_phiOld, a_FineFluxRegPtr,
             a_crseFluxRegPtr, a_crsePhiOldPtr, a_crseOldTime,
             a_crsePhiNewPtr, a_crseNewTime, a_src, a_oldTime, a_dt);

  // now subtract everything off to leave us with diffusive term
  DataIterator dit = a_DiffusiveTerm.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      // first subtract off old-time state to leave us with dphi
      tempSoln[dit()] -= a_phiOld[dit()];
      // now divide by dt to get dphi/dt
      tempSoln[dit()].divide(a_dt);
    }
  // now pass this through the no-beta operator in case there
  // is an alpha coefficient
  // this will multiply by alpha, if there is one

  m_noBetaOp->applyOpH(tempSoln,a_DiffusiveTerm);

  for (dit.begin(); dit.ok(); ++dit)
    {
      // and finally, subtract off a_src
      a_DiffusiveTerm[dit()].minus(a_src[dit()]);
    }

  // do this just to make it easier to look at
  //a_DiffusiveTerm.exchange(a_DiffusiveTerm.interval());

  // what's left should be the time-centered diffusive part of the update
}


BaseHeatSolver*
BackwardEulerSolver::new_heatSolver() const
{
  BackwardEulerSolver* new_ptr = new BackwardEulerSolver();
  if (m_helmop != NULL)
    {
      new_ptr->m_helmop=static_cast<BaseHelmholtzOp*>(m_helmop->new_levelop());
    }
  else
    {
      new_ptr->m_helmop = NULL;
    }

  if (m_noBetaOp != NULL)
    {
      new_ptr->m_noBetaOp=static_cast<BaseHelmholtzOp*>(m_helmop->new_levelop());
    }
  else
    {
      new_ptr->m_noBetaOp = NULL;
    }

  new_ptr->m_solverTol = m_solverTol;

  return static_cast<BaseHeatSolver*>(new_ptr);
}


void
BackwardEulerSolver::setSolverTolerance(Real a_solverTol)
{
  m_solverTol = a_solverTol;
  if (m_levelSolver.isDefined())
    {
      m_levelSolver.setTolerance(m_solverTol);
    }
}


Real
BackwardEulerSolver::solverTolerance() const
{
  return m_solverTol;
}

#include "NamespaceFooter.H"
