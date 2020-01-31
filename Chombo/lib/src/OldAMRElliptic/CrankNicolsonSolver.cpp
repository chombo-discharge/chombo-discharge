#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "CrankNicolsonSolver.H"
#include "FluxBox.H"
#include "HOExtrapBC.H"
#include "timeInterp.H"
#include "NamespaceHeader.H"

CrankNicolsonSolver::CrankNicolsonSolver()
{
  m_helmop = NULL;
  m_solverTol = 1.0e-7;
}

CrankNicolsonSolver::CrankNicolsonSolver(const DisjointBoxLayout& a_levelGrids,
                                         const DisjointBoxLayout* a_crseGrids,
                                         const ProblemDomain& a_domain,
                                         Real a_dxLevel,
                                         int a_nRefCrse,
                                         const BaseHelmholtzOp* a_helmop,
                                         int a_ncomp)
{
  m_helmop = NULL;
  m_solverTol = 1.0e-7;
  define(a_levelGrids, a_crseGrids, a_domain, a_dxLevel,
         a_nRefCrse, a_helmop, a_ncomp);
}

CrankNicolsonSolver::~CrankNicolsonSolver()
{
  if (m_helmop != NULL)
    {
      delete m_helmop;
      m_helmop = NULL;
    }
}

void
CrankNicolsonSolver::define(const DisjointBoxLayout& a_levelGrids,
                            const DisjointBoxLayout* a_crseGrids,
                            const ProblemDomain& a_domain,
                            Real a_dxLevel,
                            int a_nRefCrse,
                            const BaseHelmholtzOp* a_helmop,
                            int a_ncomp)
{
  m_dx = a_dxLevel;
  m_domain = a_domain;
  if (m_helmop != NULL)
    {
      delete m_helmop;
    }
  m_helmop = (BaseHelmholtzOp*) a_helmop->new_levelop();
  // check this because if crse level doesn't exist,
  // then we can get away with homogeneousOnly in levelop define
  bool homogeneousOnly = false;
  m_helmop->define(a_levelGrids, a_crseGrids, a_dxLevel,
                   a_nRefCrse, a_domain, homogeneousOnly,
                   a_ncomp);
  m_levelSolver.define(a_levelGrids, a_crseGrids, a_domain,
                       a_dxLevel, a_nRefCrse, a_helmop,
                       a_ncomp);
  m_levelSolver.setTolerance(m_solverTol);

}

void
CrankNicolsonSolver::updateSoln(LevelData<FArrayBox>& a_phiNew,
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


  int ncomp = a_phiOld.nComp();
  const DisjointBoxLayout levelGrids = a_phiNew.getBoxes();
  // probably easiest if we just keep this dataIterator around
  DataIterator dit = a_phiNew.dataIterator();

  LevelData<FArrayBox> diffusiveSource(levelGrids, ncomp);
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
      diffusiveSource[dit()].setVal(666.666);
      RHS[dit()].setVal(666.666);
      diffusiveFlux[dit()].setVal(666.666);
    }
#endif

  // now compute RHS for diffusion solve:
  // RHS = (I + 0.5*nu*dt)phi_old + dt*source
  // We can use the existing
  // operator by first scaling it.  Implicit assumption is that
  // helmholtzop currently has coefficient == -nu -- for RHS, need to
  // rescale coefficient by -dt/2, then apply operator to old-time phi.
  // finally, add dt*source to get RHS
  Real scaleMult = -0.5*a_dt;
  m_helmop->scaleBeta(scaleMult);

  // first need to compute coarse-level BC at this level's old time
  if (crseDataPtr != NULL)
    {
      timeInterp(*crseDataPtr, a_oldTime,
                 *a_crsePhiOldPtr, a_crseOldTime,
                 *a_crsePhiNewPtr, a_crseNewTime,
                 crseDataPtr->interval());
    }

  m_helmop->applyOpI(a_phiOld, crseDataPtr, RHS);
  // increment flux
  FArrayBox tempFlux;

  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& thisOldPhi = a_phiOld[dit()];
      FluxBox& thisFlux = diffusiveFlux[dit()];
      for (int dir=0; dir<SpaceDim; ++dir)
        {
          m_helmop->getFlux(tempFlux, thisOldPhi, dit(), dir);
          tempFlux.negate();
          thisFlux[dir].copy(tempFlux);
        }
    }

  // return m_helmop to original scaling
  m_helmop->scaleBeta(1.0/scaleMult);

  // now add in source term
  for (dit.begin(); dit.ok(); ++dit)
    {
      RHS[dit()].plus(a_src[dit()],a_dt);
    }

  // now construct CF-BC for viscous solve

  if (crseDataPtr != NULL)
    {
      Real newTime = a_oldTime + a_dt;
      timeInterp(*crseDataPtr, newTime,
                 *a_crsePhiOldPtr, a_crseOldTime,
                 *a_crsePhiNewPtr, a_crseNewTime,
                 crseDataPtr->interval());
    }


  // solve is (I - 0.5*dt*nu*L)phiNew = RHS
  scaleMult = 0.5*a_dt;
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


  // good first guess at soln is old soln
  a_phiOld.copyTo(a_phiOld.interval(),
                  a_phiNew, a_phiNew.interval());

  bool initializeSolnToZero = false;
  m_levelSolver.levelSolve(a_phiNew, crseDataPtr,
                           RHS, initializeSolnToZero);

  // rescale levelsolver, just in case we want to reuse this solver
  m_levelSolver.scaleBeta(1.0/scaleMult);

  // increment diffusive flux with piece from new soln
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
CrankNicolsonSolver::computeDiffusion(LevelData<FArrayBox>& a_DiffusiveTerm,
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
  updateSoln(a_DiffusiveTerm, a_phiOld, a_FineFluxRegPtr,
             a_crseFluxRegPtr, a_crsePhiOldPtr, a_crseOldTime,
             a_crsePhiNewPtr, a_crseNewTime, a_src, a_oldTime, a_dt);

  // now subtract everything off to leave us with diffusive term
  DataIterator dit = a_DiffusiveTerm.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      // first subtract off old-time state to leave us with dphi
      a_DiffusiveTerm[dit()] -= a_phiOld[dit()];
      // now divide by dt to get dphi/dt
      a_DiffusiveTerm[dit()].divide(a_dt);
      // and finally, subtract off a_src
      a_DiffusiveTerm[dit()].minus(a_src[dit()]);
    }

  // what's left should be the time-centered diffusive part of the update
}


BaseHeatSolver*
CrankNicolsonSolver::new_heatSolver() const
{
  CrankNicolsonSolver* new_ptr = new CrankNicolsonSolver();
  if (m_helmop != NULL)
    {
      new_ptr->m_helmop = static_cast<BaseHelmholtzOp*>(m_helmop->new_levelop());
    }
  else
    {
      new_ptr->m_helmop = NULL;
    }

  new_ptr->m_solverTol = m_solverTol;

  return static_cast<BaseHeatSolver*>(new_ptr);
}


void
CrankNicolsonSolver::setSolverTolerance(Real a_solverTol)
{
  m_solverTol = a_solverTol;
  if (m_levelSolver.isDefined())
    {
      m_levelSolver.setTolerance(m_solverTol);
    }
}


Real
CrankNicolsonSolver::solverTolerance() const
{
  return m_solverTol;
}


#include "NamespaceFooter.H"
