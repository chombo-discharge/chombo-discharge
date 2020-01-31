#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// DTGraves, Thurs, July 15, 1999

#include <cstdlib>
#include <iostream>
#include "SPACE.H"
#include <cmath>

#include "SPMD.H"
#include "AMRSolver.H"
#include "AMRLevelMG.H"
#include "LayoutIterator.H"
#include "ParmParse.H"
#include "NamespaceHeader.H"

void AMRLevelMG::setDefaultValues()
{
  m_isDefined = false;
  m_parent = NULL;
  m_level = -1;
  m_isDefined = false;
  m_levelopPtr = NULL;
}

void AMRLevelMG::clearMemory()
{
  if (m_levelopPtr != NULL)
    {
      delete m_levelopPtr;
    }

  m_levelopPtr = NULL;
  m_residualCopier.clear();
  m_fineExchangeCopier.clear();

}

/// Constructor
AMRLevelMG::AMRLevelMG()
{
  setDefaultValues();
}

/// Destructor
AMRLevelMG::~AMRLevelMG()
{
  clearMemory();
}

/// Define level
AMRLevelMG::AMRLevelMG(const AMRSolver* const a_parent,
                       int                    a_level,
                       const LevelOp* const   a_opin,
                       int                    a_ncomp)
{
  setDefaultValues();
  define(a_parent, a_level, a_opin, a_ncomp);
}

/// Define level
void AMRLevelMG::define(const AMRSolver* const a_parent,
                        int                    a_level,
                        const LevelOp* const   a_opin,
                        int                    a_ncomp)
{
  clearMemory();

  CH_assert(a_parent != NULL);
  CH_assert(a_parent->isDefined());
  CH_assert(a_level >= 0);
  CH_assert(a_level <= a_parent->m_finestLevel);

  m_isDefined = true;

  m_parent = a_parent;
  m_level = a_level;

  int nCoarserLevels;
  const DisjointBoxLayout& grids = a_parent->m_gridsLevel[a_level];
  const DisjointBoxLayout* baseBAPtr = NULL;
  // put reasonable value into refratio
  int reftoCoarse = 2;
  const ProblemDomain& domain = a_parent->m_domainLevel[a_level];

  if (m_level > 0)
    {
      const DisjointBoxLayout& baseGrids = a_parent->m_gridsLevel[a_level-1];

      baseBAPtr = &baseGrids;
      reftoCoarse = a_parent->m_refRatio[a_level-1];

      m_levfluxreg.define(grids, baseGrids, domain, reftoCoarse, a_ncomp);

      // nCoarserLevels = log2(reftoCoarse) - 1
      // bug fixed by petermc, 15 Apr 2002
      nCoarserLevels = -1;
      for (int interRatio = reftoCoarse; interRatio >= 2; interRatio /= 2)
        {
          nCoarserLevels++;
        }

      // (DFM 8/31/04) -- only define mginterp if level > 0
      m_mginterp.define(grids, a_ncomp, reftoCoarse, domain);
    }
  else
    {
      nCoarserLevels = 0;
    }

  Real dx = a_parent->m_dxLevel[a_level];

  m_levelMG.define(grids, baseBAPtr, dx,
                   reftoCoarse, domain,
                   nCoarserLevels, a_opin,
                   a_ncomp);

  //residual,lofphi has no ghost cells, corr and phi have 1 ghost cell
  m_lofPhi .define(grids, a_ncomp, IntVect::Zero);
  m_resid  .define(grids, a_ncomp, IntVect::Zero);
  m_corr   .define(grids, a_ncomp, IntVect::Unit);
  m_dcorr  .define(grids, a_ncomp, IntVect::Unit);
  m_phiSave.define(grids, a_ncomp, IntVect::Unit);

  if (m_level > 0)
    {
      coarsen(m_coarsenedGrids, grids, reftoCoarse);
      m_resC.define(m_coarsenedGrids, a_ncomp,IntVect::Zero);
      m_residualCopier.define(m_coarsenedGrids, a_parent->m_gridsLevel[a_level-1]);
      m_averageOp.define(grids, m_coarsenedGrids, a_ncomp, reftoCoarse);
    }

  if (m_level < a_parent->m_finestLevel)
    {
      m_fineExchangeCopier.define(a_parent->m_gridsLevel[a_level+1],
                                  a_parent->m_gridsLevel[a_level+1],
                                  IntVect::Unit);
    }

  m_levelopPtr = a_opin->new_levelop();
  m_levelopPtr->define(grids, baseBAPtr,
                       dx, reftoCoarse,
                       domain, false,a_ncomp);

  //ParmParse pp("amrlevelmg");
}

void AMRLevelMG::setNumSmoothUp(int a_numSmoothUp)
{
  CH_assert(isDefined());

  m_levelMG.setNumSmoothUp(a_numSmoothUp);
}

void AMRLevelMG::setNumSmoothDown(int a_numSmoothDown)
{
  CH_assert(isDefined());

  m_levelMG.setNumSmoothDown(a_numSmoothDown);
}

/// complete, 3-level operator
void AMRLevelMG::applyAMROperator(Vector<LevelData<FArrayBox>* >& a_phiLevel,
                                  LevelData<FArrayBox>&           a_LofPhi)
{
  CH_TIME("AMRLevelMG::applyAMROperator");
  CH_assert(isDefined());

  LevelData<FArrayBox>& phi = *a_phiLevel[m_level];

  DataIterator dit = a_LofPhi.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
  {
    a_LofPhi[dit()].setVal(0.0);
  }

  LevelData<FArrayBox>* crsePhi = NULL;

  // operator on this level including interpolated bc's from Level-1
  if (m_level > 0)
    {
      crsePhi = a_phiLevel[m_level-1];
    }

  LevelOp& levelop = *m_levelopPtr;
  levelop.applyOpI(phi, crsePhi, a_LofPhi);

  // now reflux to enforce flux-matching from finer levels...
  if (m_level < m_parent->m_finestLevel)
    {
      reflux(a_phiLevel, a_LofPhi);
    }
}

/// complete, 3-level operator with homogeneous physical boundary conditions
void AMRLevelMG::applyAMROperatorHphys(
                     Vector<LevelData<FArrayBox>* >& a_phiLevel,
                     LevelData<FArrayBox>&           a_LofPhi)
{
  CH_assert(isDefined());

  LevelData<FArrayBox>& phi = *a_phiLevel[m_level];

  DataIterator dit = a_LofPhi.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
    {
      a_LofPhi[dit()].setVal(0.0);
    }

  LevelData<FArrayBox>* crsePhi = NULL;

  // operator on this level including interpolated bc's from Level-1
  if (m_level > 0)
    {
      crsePhi = a_phiLevel[m_level-1];
    }

  LevelOp& levelop = *m_levelopPtr;
  levelop.applyOpIcfHphys(phi, crsePhi, a_LofPhi);

  // now reflux to enforce flux-matching from finer levels...
  if (m_level < m_parent->m_finestLevel)
    {
      reflux(a_phiLevel, a_LofPhi);
    }
}

/**
      compute complete, 3-level residual
      and put it into local data m_resid
  */
void AMRLevelMG::computeAMRResidual(
                     Vector<LevelData<FArrayBox>* >&       a_phiLevel,
                     const Vector<LevelData<FArrayBox>* >& a_rhsLevel)
{
  CH_TIME("AMRLevelMG::computeAMRResidual");
  CH_assert(isDefined());

  LevelData<FArrayBox>& rhs = *a_rhsLevel[m_level];

  applyAMROperator(a_phiLevel, m_resid);

  DataIterator dit = rhs.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
    {
      FArrayBox& residmf = m_resid[dit()];
      FArrayBox& rhsmf = rhs[dit()];

      residmf -= rhsmf;
      residmf.negate();
    }
}

/**
      compute complete, 3-level residual with homogeneous physical BCs
      and inhomogeneous C/F BCs) and put it into local data m_resid
  */
void AMRLevelMG::computeAMRResidualHphys(
                     Vector<LevelData<FArrayBox>* >&       a_phiLevel,
                     const Vector<LevelData<FArrayBox>* >& a_rhsLevel)
{
  CH_assert(isDefined());

  LevelData<FArrayBox>& rhs = *a_rhsLevel[m_level];

  applyAMROperatorHphys(a_phiLevel, m_resid);

  DataIterator dit = rhs.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
    {
      FArrayBox& residmf = m_resid[dit()];
      FArrayBox& rhsmf = rhs[dit()];

      residmf -= rhsmf;
      residmf.negate();
    }
}

/**
      compute complete, 3-level residual
      and put it into argument a_resid
  */
void AMRLevelMG::computeAMRResidual(
                     LevelData<FArrayBox>&                 a_resid,
                     Vector<LevelData<FArrayBox>* >&       a_phiLevel,
                     const Vector<LevelData<FArrayBox>* >& a_rhsLevel)
{
  CH_TIME("AMRLevelMG::computeAMRResidual(resid)");
  CH_assert(isDefined());

  LevelData<FArrayBox>& rhs = *a_rhsLevel[m_level];

  applyAMROperator(a_phiLevel, a_resid);

  DataIterator dit = rhs.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
    {
      FArrayBox& residmf = a_resid[dit()];
      FArrayBox& rhsmf = rhs[dit()];

      residmf -= rhsmf;
      residmf.negate();
    }
}

/**
    Sweep down v-cycle -- assumes phi,rhs already in residual-correction
    form (homogenous physical BCs, inhomogeneous C/F BCs)
*/
void AMRLevelMG::downSweep(Vector<LevelData<FArrayBox>* >&       a_phiLevel,
                           const Vector<LevelData<FArrayBox>* >& a_rhsLevel)
{
  CH_assert(isDefined());

  LevelData<FArrayBox>& phi = *a_phiLevel[m_level];
  LevelData<FArrayBox>& rhs = *a_rhsLevel[m_level];

  CH_assert(phi.nComp() == rhs.nComp());

  DataIterator dit = phi.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
    {
      m_phiSave[dit()].copy(phi[dit()]);
      m_corr[dit()].setVal(0.);
    }

  // Apply smoother.
  smooth(m_corr,m_resid);

  // Update phi.
  for (dit.reset(); dit.ok(); ++dit)
    {
      phi[dit()] += m_corr[dit()];
    }

  // Compute residual for the next coarser level.
  // Form residual on next coarser level.
  if (m_level > 0)
    {
      // Initialize correction at the next coarser level.
      // this has to be done before applyOpI
      //int reftoCoarse = m_parent->m_refRatio[m_level-1];
      AMRLevelMG& nextCoarser = *m_parent->m_amrmgLevel[m_level-1];

      DataIterator ditCoar = nextCoarser.m_corr.dataIterator();
      for (ditCoar.reset(); ditCoar.ok(); ++ditCoar)
        {
          nextCoarser.m_corr[ditCoar()].setVal(0.);
        }

      // petermc, 18 Apr 2002, replaced this with applyOpH
      // because coarser level is zero anyway.
      // m_levelopPtr->applyOpIcfHphys(m_corr, &nextCoarser.m_corr, m_lofPhi);
      m_levelopPtr->applyOpH(m_corr, m_lofPhi);

      for (dit.reset(); dit.ok(); ++dit)
        {
          m_lofPhi[dit()] -= m_resid[dit()];
          m_lofPhi[dit()].negate();
        }

      m_averageOp.averageToCoarse(m_resC, m_lofPhi);

      nextCoarser.computeAMRResidualHphys(a_phiLevel, a_rhsLevel);

      // overwrite residual on Level-1 with
      //coarsened residual from this level
      m_resC.copyTo(m_resC.interval(),
                    nextCoarser.m_resid,
                    nextCoarser.m_resid.interval(), m_residualCopier);
    }
}

/**
      Sweep up v-cycle -- assumes phi, rhs already in residual-correction
      form (homogeneous physical BCs, inhomogeneous C/F BC's)
  */
void AMRLevelMG::upSweep(Vector<LevelData<FArrayBox> *>&       a_phiLevel,
                         const Vector<LevelData<FArrayBox> *>& a_rhsLevel)
{

  CH_assert(isDefined());

  LevelData<FArrayBox>&       phi = *a_phiLevel[m_level];
  const LevelData<FArrayBox>& rhs = *a_rhsLevel[m_level];

  CH_assert(phi.nComp() == rhs.nComp());

  // Interpolate corrections from next coarser level to the correction
  // field at this level.
  if (m_level > 0)
    {
      // Initialize correction at the next coarser level.
      //this uses the two-level only operator (Lnf)
      AMRLevelMG& nextCoarser = *m_parent->m_amrmgLevel[m_level-1];
      m_mginterp.interpToFine(m_corr, nextCoarser.m_corr);

      m_levelopPtr->applyOpIcfHphys(m_corr,
                                    &nextCoarser.m_corr,
                                    m_lofPhi);
    }

  // Finish resid calc,
  // initialize the correction to the correction,
  // and apply smoother.
  DataIterator dit = m_corr.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
    {
      m_resid[dit()] -= m_lofPhi[dit()];
      m_dcorr[dit()].setVal(0.);
    }

  smooth(m_dcorr,m_resid);

  // Increment correction.
  // Increment saved value of phi with corr, and overwrite phi.
  for (dit.reset(); dit.ok(); ++dit)
    {
      m_corr[dit()] += m_dcorr[dit()];
      m_phiSave[dit()] += m_corr[dit()];
      phi[dit()].copy(m_phiSave[dit()]);
    }
}

/** set convergence metric -- this is simply a pass-through to LevelOps
 */
void
AMRLevelMG::setConvergenceMetric(Real a_metric, int a_comp)
{

}

/**  smooth a_soln */
void AMRLevelMG::smooth(LevelData<FArrayBox>&       a_soln,
                        const LevelData<FArrayBox>& a_rhs)
{
  CH_assert(isDefined());
  CH_assert(a_soln.ghostVect() >= IntVect::Unit);
  CH_assert(a_soln.getBoxes() == m_parent->m_gridsLevel[m_level]);
  CH_assert(a_rhs.getBoxes() == m_parent->m_gridsLevel[m_level]);

  // no bottom solve here
  m_levelMG.mgRelax(a_soln, a_rhs, false);
}

/** compute and return norm of internal data m_resid.  */
Vector<Real> AMRLevelMG::computeResidualNorm(int a_normType) const
{
  CH_assert(isDefined());

  Vector<Real> normPeterson(m_resid.nComp(),0.0);
  normPeterson = computeNorm(m_resid, a_normType);

  return normPeterson;
}

LevelOp* AMRLevelMG::levelOpPtr() const
{
  return m_levelopPtr;
}

/// has define function been called?  if not, most fcns won't work
bool AMRLevelMG::isDefined() const
{
  return m_isDefined;
}

// returns normType norm of mfab
Vector<Real>
AMRLevelMG::computeNorm(const LevelData<FArrayBox>& a_mfinput,
                        int                         a_normType) const
{
  CH_TIME("AMRLevelMG::computeNorm");
  CH_assert(isDefined());
  CH_assert((a_normType >= 0)&&(a_normType <= 2));

  // create temps so inputs do not get changed
  LevelData<FArrayBox>  mfab;
  mfab.define(a_mfinput, a_mfinput.interval());

  // compute the bloody norm
  Vector<Real> normTot(mfab.nComp(),0.0);
  DataIterator dit = mfab.dataIterator();
  int ncomp = mfab.nComp();

  if (m_level < m_parent->m_finestLevel)
    {
      const DisjointBoxLayout& fineGrids = m_parent->m_gridsLevel[m_level+1];
      int nrefFine =  m_parent->m_refRatio[m_level];

      for (dit.ok(); dit.ok(); ++dit)
        {
          LayoutIterator litFine = fineGrids.layoutIterator();
          for (litFine.reset(); litFine.ok(); ++litFine)
            {
              Box overlayBox = mfab[dit()].box();
              Box coarsenedGrid = coarsen(fineGrids[litFine()], nrefFine);

              overlayBox &= coarsenedGrid;

              if (!overlayBox.isEmpty())
                {
                  mfab[dit()].setVal(0.0,overlayBox,0, ncomp);
                }
            }
        }
    }

  for (int comp=0; comp<mfab.nComp(); comp++)
    {
      Interval thisInterval(comp,comp);

      // this function does the wacky boxlib-type norm
      normTot[comp] = norm(mfab, thisInterval, a_normType);

      // so we normalize it for sanity's sake
      if (a_normType != 0)
        {
          Real dx = m_parent->m_dxLevel[m_level];
          Real expon = SpaceDim/a_normType;
          Real normfac = pow(dx, expon);
          normTot[comp] *= normfac;
        }
    } // end loop over components

  return normTot;
}

/**
   reflux enforces flux-matching from finer levels.
   Steps involved:
   a)initialize flux register with coarse flux
   b)increment flux register with sum of fine fluxes
   b)modify solution with flux difference
*/
void AMRLevelMG::reflux(Vector<LevelData<FArrayBox> *>& a_phiLevel,
                       LevelData<FArrayBox>&            a_Lofphi)
{
  CH_TIME("AMRLevelMG::reflux");
  CH_assert(isDefined());
  CH_assert(m_level < m_parent->m_finestLevel);

  AMRLevelMG& higher = *(m_parent->m_amrmgLevel[m_level + 1]);

  CH_assert(higher.isDefined());

  LevelFluxRegister& finelevfr = higher.m_levfluxreg;

  // finelevfr.setToZero(); called in initFRCoarse
  initFRCoarse(a_phiLevel);

  incrementFRFine(a_phiLevel);
  Real dxlevel = m_parent->m_dxLevel[m_level];

  CH_assert(dxlevel > 0);

  Real scale = 1.0/dxlevel;
  finelevfr.reflux(a_Lofphi, scale);
}

///internally useful functions
/**
   Fill fluxregisters with coarse flux
*/
void AMRLevelMG::initFRCoarse(Vector<LevelData<FArrayBox> *>& a_phiLevel)
{
  CH_TIME("AMRLevelMG::initFRCoarse");
  CH_assert(isDefined());
  CH_assert(m_level < m_parent->m_finestLevel);

  AMRLevelMG& higher = *(m_parent->m_amrmgLevel[m_level + 1]);

  CH_assert(higher.isDefined());

  LevelFluxRegister& finelevfr = higher.m_levfluxreg;

  // looping through every box in this grid
  // and fill fluxregisters with coarse flux
  finelevfr.setToZero();
  LevelData<FArrayBox>& phic = *a_phiLevel[m_level];

  Interval interv(0,a_phiLevel[m_level]->nComp()-1);

  DataIterator dit = phic.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
    {
      const FArrayBox& coarfab = phic[dit()];

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          FArrayBox coarflux;
          m_levelopPtr->getFlux(coarflux, coarfab, dit(), idir);

          Real scale = 1.0;
          finelevfr.incrementCoarse(coarflux,
                                    scale,dit(),
                                    interv,interv,idir);
        }
    }
}

/**
    increment fluxregisters with fine flux
*/
void AMRLevelMG::incrementFRFine(Vector<LevelData<FArrayBox> *>& a_phiLevel)
{
  CH_TIME("AMRLevelMG::incrementFRFine");
  CH_assert(isDefined());
  CH_assert(m_level < m_parent->m_finestLevel);

  AMRLevelMG& higher = *(m_parent->m_amrmgLevel[m_level+1]);

  CH_assert(higher.isDefined());

  LevelFluxRegister& finelevfr = higher.m_levfluxreg;

  // flux register stuff
  // looping through every box in finer grid
  // and increment fluxregisters with fine flux
  LevelData<FArrayBox>& phif = *a_phiLevel[m_level + 1];
  LevelData<FArrayBox>& phic = *a_phiLevel[m_level];
  Interval interv = phif.interval();

  // need to do cfinterpolation in aggregate before
  // any getFineData calls are made
  higher.m_levelopPtr->CFInterp(phif, phic);
  phif.exchange(phif.interval(),  m_fineExchangeCopier);
  const DisjointBoxLayout& grids = phif.getBoxes();
  IntVect phiGhost = phif.ghostVect();

  DataIterator ditf = phif.dataIterator();
  for (ditf.reset(); ditf.ok(); ++ditf)
    {
      const FArrayBox& phifFab = phif[ditf()];
      const Box& gridbox = grids[ditf()];

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          int normalGhost = phiGhost[idir];
          SideIterator sit;
          for (sit.begin(); sit.ok(); sit.next())
            {
              Side::LoHiSide hiorlo = sit();
              Box fabbox;

              if (sit() == Side::Lo)
                {
                  // assumption here that the stencil required
                  // to compute the flux in the normal direction
                  // is 2* the number of ghost cells for phi
                  // (which is a reasonable assumption, and probably
                  // better than just assuming you need one cell on
                  // either side of the interface
                  // (dfm 8-4-06)
                  fabbox = adjCellLo(gridbox,idir, 2*normalGhost);
                  fabbox.shift(idir, 1);
                }
              else
                {
                  fabbox = adjCellHi(gridbox,idir, 2*normalGhost);
                  fabbox.shift(idir, -1);
                }

              // just in case we need ghost cells in the transverse direction
              // (dfm 8-4-06)
              for (int otherDir=0; otherDir<SpaceDim; ++otherDir)
                {
                  if (otherDir != idir)
                    {
                      fabbox.grow(otherDir, phiGhost[otherDir]);
                    }
                }

              CH_assert(!fabbox.isEmpty());

              FArrayBox phifab(fabbox, phif.nComp());
              phifab.copy(phifFab);

              FArrayBox fineflux;
              higher.m_levelopPtr->getFlux(fineflux, phifab, ditf(), idir);

              Real scale = 1.0;
              finelevfr.incrementFine(fineflux, scale, ditf(),
                                      interv, interv, idir, hiorlo);
            }
        }
    }
}
#include "NamespaceFooter.H"
