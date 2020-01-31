#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// DTGraves, Weds, July 21, 1999

#include "PoissonOp.H"
#include "PoissonOpF_F.H"
#include "LayoutIterator.H"
#include "BoxIterator.H"
#include "BiCGStabSmoother.H"
#include "CFStencil.H"
#include "CH_Timer.H"
#include "NamespaceHeader.H"

void PoissonOp::setDomainGhostBC(const DomainGhostBC& a_dombcIn)
{
  m_domghostbc = a_dombcIn;
  m_isBCDefined = true;
}

void PoissonOp::setBottomSmoother(const BaseBottomSmoother& a_bottomSmoother)
{
  if (m_bottomSmootherPtr != NULL)
    {
      delete m_bottomSmootherPtr;
    }

  m_bottomSmootherPtr = a_bottomSmoother.new_bottomSmoother();
}

bool PoissonOp::isDefined() const
{
  return (m_isDefined && m_isBCDefined);
}

LevelOp* PoissonOp::new_levelop() const
{
  // only check to see if the boundary conditions are
  // defined.  not necessary to have whole thing defined.
  // the boundary conditions need to be defined because
  // the interface does not know about the boundary
  // condtions.  the solvers call the other define functions
  // in our inimitable two-stage construction
  CH_assert(m_isBCDefined);

  PoissonOp* newop = new PoissonOp();

  if (newop == NULL)
    {
      MayDay::Error("Out of Memory in PoissonOp::new_levelop");
    }

  newop->setDomainGhostBC(m_domghostbc);

  if (m_bottomSmootherPtr != NULL)
    {
      newop->setBottomSmoother(*m_bottomSmootherPtr);
    }

  newop->m_numPrecondSmooth = m_numPrecondSmooth;

  return static_cast<LevelOp*>(newop);
}

/***********************/
// return to undefined state
/***********************/
void PoissonOp::setDefaultValues()
{
  m_dxLevel = -1.0;
  m_refRatio = -1;

  m_ihcfiEnabled = false;
  m_isBCDefined = false;

  // default is to use BiCGStab as bottom smoother.
  m_bottomSmootherPtr = new BiCGStabSmoother;

  m_numPrecondSmooth = 2;

  m_isDefined = false;
}

/***********************/
// return to undefined state
/***********************/
void PoissonOp::clearMemory()
{
  m_dxLevel = -1.0;
  m_refRatio = -1;

  m_ihcfiEnabled = false;

  // don't touch bottom smoother here, since we want it to be able
  // to pass through the define statement unaltered

  m_isDefined = false;
}

/***********************/
// default constructor
/***********************/
PoissonOp::PoissonOp()
  :LevelOp()
{
  setDefaultValues();
}

/***********************/
/***********************/
void PoissonOp::define(const DisjointBoxLayout& a_grids,
                       const DisjointBoxLayout* a_baseBAPtr,
                       Real                     a_dxLevel,
                       int                      a_refRatio,
                       const Box&               a_domain,
                       bool                     a_homogeneousOnly,
                       int                      a_ncomp)
{
  ProblemDomain probdomain(a_domain);

  define(a_grids, a_baseBAPtr, a_dxLevel, a_refRatio, probdomain,
         a_homogeneousOnly, a_ncomp);
}

void PoissonOp::define(const DisjointBoxLayout& a_grids,
                       const DisjointBoxLayout* a_baseBAPtr,
                       Real                     a_dxLevel,
                       int                      a_refRatio,
                       const ProblemDomain&     a_domain,
                       bool                     a_homogeneousOnly,
                       int                      a_ncomp)
{
  CH_TIME("PoissonOp::define");
  clearMemory();

  m_isDefined = true;

  CH_assert(!a_domain.isEmpty());
  CH_assert(a_grids.checkPeriodic(a_domain));

  m_dxLevel = a_dxLevel;
  m_domain = a_domain;
  m_grids = a_grids;

  m_exchangeCopier.define(a_grids, a_grids, IntVect::Unit, true);

  m_refRatio = a_refRatio;
  m_dxCrse = m_refRatio*m_dxLevel;

  m_ihcfiEnabled = (!a_homogeneousOnly);

  m_ncomp = a_ncomp;

  if (m_ihcfiEnabled)
    {
      m_quadCFI.define(a_grids,a_baseBAPtr, a_dxLevel,
                       a_refRatio, m_ncomp, m_domain);
    }

  if (a_baseBAPtr != NULL)
    {
      m_baseBA = *a_baseBAPtr;
    }

  Vector<Box> periodicBoxes;
  CFStencil::buildPeriodicVector(periodicBoxes, m_domain, m_grids);
  DataIterator lit = m_grids.dataIterator();
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      m_loCFIVS[idir].define(m_grids);
      m_hiCFIVS[idir].define(m_grids);

      for (lit.begin(); lit.ok(); ++lit)
        {
          m_loCFIVS[idir][lit].define(m_domain, m_grids.get(lit),
                                      periodicBoxes, idir,Side::Lo);
          m_hiCFIVS[idir][lit].define(m_domain, m_grids.get(lit),
                                      periodicBoxes, idir,Side::Hi);
        }
    }
}

// a_reftofine is the refinement between this operator and a_opfine
void PoissonOp::define(const LevelOp* a_opfine,
                       int            a_reftoFine)
{
  clearMemory();

  const PoissonOp* opfineptr = dynamic_cast<const PoissonOp*>(a_opfine);

  if (opfineptr == NULL)
    {
      MayDay::Error("PoissonOp::define: casting failed");
    }

  const PoissonOp& opfine = *opfineptr;

  CH_assert(opfine.isDefined());
  CH_assert(a_reftoFine > 0);

  m_isDefined = true;
  m_ihcfiEnabled = false;

  setDomainGhostBC(opfine.m_domghostbc);

  m_dxLevel = (opfine.m_dxLevel)*a_reftoFine;
  m_domain = coarsen(opfine.m_domain,a_reftoFine);

  coarsen(m_grids, opfine.m_grids, a_reftoFine);

  m_exchangeCopier.define(m_grids, m_grids, IntVect::Unit, true);

  m_dxCrse = opfine.m_dxCrse;
  m_ncomp = opfine.m_ncomp;

  // this refratio is the ratio between this level and the
  // next level in the amr hierarchy
  m_refRatio = -1;

  // don't need to define the cfinterp since inhomogeneous cfinterpolaton
  // is never done here.  We shall leave m_quadCFI and base_ba undefined and
  // and just leave the refratio -1.  We DO need fineinterpIVS for doing
  // homogeneous cfinterpolation
  DataIterator lit = m_grids.dataIterator();
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      m_loCFIVS[idir].define(m_grids);
      m_hiCFIVS[idir].define(m_grids);

      for (lit.begin(); lit.ok(); ++lit)
        {
          m_loCFIVS[idir][lit()].define(m_domain, m_grids.get(lit()),
                                        m_grids, idir,Side::Lo);
          m_hiCFIVS[idir][lit()].define(m_domain, m_grids.get(lit()),
                                        m_grids, idir,Side::Hi);
        }
    }

  setBottomSmoother(*(opfine.m_bottomSmootherPtr));
}

/***********************/
//  apply coarse-fine boundary conditions -- assume that phi grids
//  are grown by one
/***********************/
void PoissonOp::CFInterp(LevelData<FArrayBox>&       a_phif,
                         const LevelData<FArrayBox>& a_phic)
{
  CH_assert(isDefined());
  CH_assert(m_ihcfiEnabled);
  CH_assert( a_phif.ghostVect() >= IntVect::Unit);
  CH_assert(m_quadCFI.isDefined());

  m_quadCFI.coarseFineInterp(a_phif, a_phic);
}

/***********************/
// does homogeneous coarse/fine interpolation
/***********************/
void PoissonOp::homogeneousCFInterp(LevelData<FArrayBox>& a_phif)
{
  CH_assert(isDefined());
  CH_assert(a_phif.ghostVect() >= IntVect::Unit);

  DataIterator dit = a_phif.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const DataIndex& datInd = dit();

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          SideIterator sit;
          for (sit.begin(); sit.ok(); sit.next())
            {
              homogeneousCFInterp(a_phif,datInd,idir,sit());
            }
        }
    }
}

/***********************/
// does homogeneous coarse/fine interpolation
/***********************/
void PoissonOp::homogeneousCFInterp(LevelData<FArrayBox>& a_phif,
                                    const DataIndex&      a_datInd,
                                    int                   a_idir,
                                    Side::LoHiSide        a_hiorlo)
{
  CH_TIME("PoissonOp::homogeneousCFInterp");
  CH_assert(isDefined());
  CH_assert((a_idir >= 0) && (a_idir  < SpaceDim));
  CH_assert((a_hiorlo == Side::Lo )||(a_hiorlo == Side::Hi ));
  CH_assert( a_phif.ghostVect() >= IntVect::Unit);
  //  CH_assert (m_ncomp == a_phif.nComp());

  const CFIVS* cfivsPtr = NULL;

  if (a_hiorlo == Side::Lo)
    {
      cfivsPtr = &m_loCFIVS[a_idir][a_datInd];
    }
  else
    {
      cfivsPtr = &m_hiCFIVS[a_idir][a_datInd];
    }

  if (cfivsPtr->isPacked())
    {
      int ihiorlo = sign(a_hiorlo);
      FArrayBox& phiFab = a_phif[a_datInd];
      FORT_OLDINTERPHOMO(CHF_FRA(phiFab),
                         CHF_BOX(cfivsPtr->packedBox()),
                         CHF_CONST_REAL(m_dxLevel),
                         CHF_CONST_REAL(m_dxCrse),
                         CHF_CONST_INT(a_idir),
                         CHF_CONST_INT(ihiorlo));

    }
  else
    {
      const IntVectSet& interpIVS = cfivsPtr->getFineIVS();

      if (!interpIVS.isEmpty())
        {
          // Assuming homogenous, interpolate on fine ivs
          interpOnIVSHomo(a_phif, a_datInd, a_idir,
                          a_hiorlo, interpIVS);
        }
    }
}

void PoissonOp::interpOnIVSHomo(LevelData<FArrayBox>& a_phif,
                                const DataIndex&      a_datInd,
                                const int             a_idir,
                                const Side::LoHiSide  a_hiorlo,
                                const IntVectSet&     a_interpIVS)
{
  CH_assert(isDefined());
  CH_assert((a_idir >= 0) && (a_idir  < SpaceDim));
  CH_assert((a_hiorlo == Side::Lo )||(a_hiorlo == Side::Hi ));
  CH_assert(a_phif.ghostVect() >= IntVect::Unit);

  IVSIterator fineIVSit(a_interpIVS);
  FArrayBox& a_phi = a_phif[a_datInd];

  CH_assert(m_ncomp = a_phi.nComp());

  // much of these scalar values can be precomputed and stored if
  // we ever need to speed-up this function (ndk)
  Real x1 = m_dxLevel;
  Real x2 = 0.5*(3.0*m_dxLevel+m_dxCrse);

  Real denom = 1.0-((x1+x2)/x1);
  Real idenom = 1/(denom); // divide is more expensive usually

  Real x = 2.0*m_dxLevel;
  Real xsquared = x*x;

  Real m1 = 1.0/(x1*x1);
  Real m2 = 1.0/(x1*(x1-x2));

  Real q1 = 1.0/(x1-x2);
  Real q2 = x1+x2;

  int ihilo = sign(a_hiorlo);
  Real pa,pb,a,b;

  IntVect ivf;
  for (fineIVSit.begin(); fineIVSit.ok(); ++fineIVSit)
    {
      ivf = fineIVSit();

      // quadratic interpolation
      for (int ivar = 0; ivar < m_ncomp; ivar++)
        {
          ivf[a_idir]-=2*ihilo;
          pa = a_phi(ivf, ivar);

          ivf[a_idir]+=ihilo;
          pb = a_phi(ivf, ivar);

          a = ((pb-pa)*m1 - (pb)*m2)*idenom;
          b = (pb)*q1 - a*q2;

          a_phi(fineIVSit(), ivar) = a*xsquared + b*x + pa;

          // this returns the ivf to the original point
          ivf[a_idir]+=ihilo;
        } //end loop over components
    } //end loop over fine intvects
}

PoissonOp::~PoissonOp()
{
  clearMemory();

  if (m_bottomSmootherPtr != NULL)
    {
      delete m_bottomSmootherPtr;
      m_bottomSmootherPtr = NULL;
    }
}

/***********************/
// this smoother assumes problem has
// already been put into residual-correction
//  form, so that crse/fine BC's are homogeneous
/***********************/
void
PoissonOp::smooth(LevelData<FArrayBox>&       a_phi,
                  const LevelData<FArrayBox>& a_rhs)
{
  CH_assert(isDefined());

  levelGSRB(a_phi,a_rhs);
}

/**************************/
// this preconditioner first initializes phihat to (IA)phihat = rhshat
// (diagonization of L -- A is the matrix version of L)
// then smooths with a couple of passes of levelGSMC
/**************************/
void
PoissonOp::levelPreconditioner(LevelData<FArrayBox>&       a_phihat,
                               const LevelData<FArrayBox>& a_rhshat)
{
  // diagonal term of this operator is 4/h/h in 2D, 6/h/h in 3D,
  // so inverse of this is our initial multiplier
  Real mult = -m_dxLevel*m_dxLevel/(2.0*SpaceDim);
  Interval comps = a_phihat.interval();

  CH_assert(a_phihat.nComp() == a_rhshat.nComp());

  // don't need to use a Copier -- plain copy will do
  DataIterator dit = a_phihat.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      a_phihat[dit()].copy(a_rhshat[dit()].box(), comps,
                           a_rhshat[dit()].box(), a_rhshat[dit()],
                           comps);

      a_phihat[dit()] *= mult;
    }

  // don't need to fill in ghost cells, so we construct
  // our own copier here.
  //
  // Copier thisCopier(a_rhshat.getBoxes(), a_phihat.getBoxes(),
  //                   IntVect::Zero);
  // a_rhshat.copyTo(comps, a_phihat, comps, thisCopier);
  //
  // DataIterator dit = a_phihat.dataIterator();
  // for (dit.begin(); dit.ok(); ++dit)
  //   {
  //     a_phihat[dit()] *= mult;
  //   }


  for (int i=0; i<m_numPrecondSmooth; i++)
    {
      levelGSRB(a_phihat, a_rhshat);
      levelGSRB(a_phihat, a_rhshat);
    }
}

/***********************/
//  levelGSRB does GSRB on a level.  it has no knowledge of overlying
//  finer grids, and assumes coarse/fine BC's are homogeneous
/***********************/
void PoissonOp::levelGSRB(LevelData<FArrayBox>&       a_phi,
                          const LevelData<FArrayBox>& a_rhs)
{
  CH_assert(isDefined());
  CH_assert(a_phi.isDefined());
  CH_assert(a_rhs.isDefined());
  CH_assert(a_phi.ghostVect() >= IntVect::Unit);
  CH_assert(a_phi.nComp() == a_rhs.nComp());

  // do first red, then black passes
  for (int whichPass = 0; whichPass <= 1; whichPass++)
    {
      // do coarse/fine and copy bc's
      //should be done patch by patch.
      homogeneousCFInterp(a_phi);

      // now step through grids...
      Real dx = m_dxLevel;
      DataIterator dit = a_phi.dataIterator();

      //fill in intersection of ghostcells and a_phi's boxes
      a_phi.exchange(a_phi.interval(), m_exchangeCopier);

      for (dit.begin(); dit.ok(); ++dit)
        {
          // invoke physical BC's where necessary
          m_domghostbc.applyHomogeneousBCs(a_phi[dit()],
                                           m_domain,
                                           m_dxLevel);
        }

      // fill in intersection of ghostcells and a_phi's boxes
      // dfm -- for a 5 point stencil, this should not be necessary
      // a_phi.exchange(a_phi.interval(), m_exchangeCopier);
      for (dit.begin(); dit.ok(); ++dit)
        {
#ifndef NDEBUG
          FArrayBox& thisPhi = a_phi[dit()];
          const FArrayBox& thisRhs = a_rhs[dit()];
          const Box& thisBox = m_grids.get(dit());

          CH_assert(thisRhs.box().contains(thisBox));
          CH_assert(thisPhi.box().contains(grow(thisBox,1)));
#endif

          FORT_OLDGSRBLEVELLAP(CHF_FRA(a_phi[dit()]),
                               CHF_CONST_FRA(a_rhs[dit()]),
                               CHF_BOX(m_grids.get(dit())),
                               CHF_CONST_REAL(dx),
                               CHF_CONST_INT(whichPass));

        } // end loop through grids
    } // end loop through red-black
}
/***************/
/***************/
void PoissonOp::levelGSMC(LevelData<FArrayBox>&       a_phi,
                          const LevelData<FArrayBox>& a_rhs)
{
  int comp = 0;
  int ncomp = 1;
  LevelData<FArrayBox> residLD(m_grids, ncomp, IntVect::Zero);
  LevelData<FArrayBox> lofphi(m_grids, ncomp, IntVect::Zero);

  IntVect color = IntVect::Zero;
  IntVect limit = IntVect::Unit;

  Real weight = 0.0;

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      weight += -2.0 / (m_dxLevel*m_dxLevel);
    }
  weight = 1.0/weight;

  // Loop over all possibilities (in all dimensions)
  while (color[SpaceDim-1] <= limit[SpaceDim-1])
    {
      // Get the residual everywhere
      applyOpH(a_phi, lofphi);
      for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
        {
          residLD[dit()].copy(a_rhs[dit()]);
          residLD[dit()] -= lofphi[dit()];
        }

      int ibox = 0;
      for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
        {
          FArrayBox& phi = a_phi[dit()];
          const Box& box = m_grids.get(dit());

          FArrayBox phiSave(box, ncomp);
          phiSave.copy(phi);

          const FArrayBox& resid = residLD[dit()];

          IntVect loIV = box.smallEnd();
          IntVect hiIV = box.bigEnd();

          for (int idir = 0; idir < SpaceDim; idir++)
            {
              if (loIV[idir] % 2 != color[idir])
                {
                  loIV[idir]++;
                }
            }

          if (loIV <= hiIV)
            {
              Box coloredBox(loIV,hiIV);

              FORT_GSMULTICOLOR(CHF_FRA1(phi,comp),
                                CHF_CONST_FRA1(phiSave,comp),
                                CHF_CONST_REAL(weight),
                                CHF_CONST_FRA1(resid,comp),
                                CHF_BOX(coloredBox));
            }

          ibox++;
        }

      color[0]++;

      // If greater that the maximum it can be...
      if (color[0] > limit[0])
        {
          // For all dimensions...
          for (int i = 0; i < SpaceDim-1; i++)
            {
              color[i] = 0;
              color[i+1]++;

              // If the next dimension's power is okay go on
              if (color[i+1] <= limit[i+1])
                {
                  break;
                }
            }
        }
    }
}

/***********************/
// evaluate operator, inhomogeneous C/F boundary conditions
// note that operator applies C/F interpolation and copy bc's
/***********************/
void PoissonOp::applyOpI(LevelData<FArrayBox>&       a_phi,
                         const LevelData<FArrayBox>* a_phicPtr,
                         LevelData<FArrayBox>&       a_lofPhi)
{
  CH_TIME("PoissonOp::applyOpI");
  CH_assert(isDefined());
  CH_assert(m_ihcfiEnabled);
  CH_assert(a_phi.nComp() == a_lofPhi.nComp());
  CH_assert(a_phi.nComp() == m_ncomp);

  if (a_phicPtr != NULL)
    {
      // apply C/F boundary conditions...
      CFInterp(a_phi,*a_phicPtr);
    }

  //fill in intersection of ghostcells and phi's boxes
  a_phi.exchange(a_phi.interval(), m_exchangeCopier);

  DataIterator dit = a_phi.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      m_domghostbc.applyInhomogeneousBCs(a_phi[dit()],
                                         m_domain,
                                         m_dxLevel);
    }

  //fill in intersection of ghostcells and phi's boxes
  // for a 5-point operator, this should not be necessary
  //a_phi.exchange(a_phi.interval(), m_exchangeCopier);

  Real dx = m_dxLevel;

  for (dit.begin(); dit.ok(); ++dit)
    {
#ifndef NDEBUG
      FArrayBox& thisPhi = a_phi[dit()];
      FArrayBox& thisLPhi = a_lofPhi[dit()];
      const Box& thisBox = m_grids.get(dit());

      CH_assert(thisPhi.box().contains(grow(thisBox,1)));
      CH_assert(thisLPhi.box().contains(thisBox));
#endif

      FORT_OLDOPERATORLAP(CHF_FRA(a_lofPhi[dit()]),
                          CHF_CONST_FRA(a_phi[dit()]),
                          CHF_BOX(m_grids.get(dit())),
                          CHF_CONST_REAL(dx));
    }
}

/***********************/
// evaluate operator, inhomogeneous C/F boundary conditions
// note that operator applies C/F interpolation and copy bc's
// homogeneous physical bcs
/***********************/
void PoissonOp::applyOpIcfHphys(LevelData<FArrayBox>&       a_phi,
                                const LevelData<FArrayBox>* a_phicPtr,
                                LevelData<FArrayBox>&       a_lofPhi)
{
  CH_assert(isDefined());
  CH_assert(m_ihcfiEnabled);
  CH_assert(a_phi.nComp() == a_lofPhi.nComp());
  CH_assert(a_phi.nComp() == m_ncomp);

  if (a_phicPtr != NULL)
    {
      // apply C/F boundary conditions...
      CFInterp(a_phi,*a_phicPtr);
    }

  // fill in intersection of ghostcells and phi's boxes
  a_phi.exchange(a_phi.interval(), m_exchangeCopier);

  DataIterator dit = a_phi.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      m_domghostbc.applyHomogeneousBCs(a_phi[dit()],
                                       m_domain,
                                       m_dxLevel);
    }

  // fill in intersection of ghostcells and phi's boxes
  // for a 5 point operator, this shouldn't be necessary
  // a_phi.exchange(a_phi.interval(), m_exchangeCopier);

  Real dx = m_dxLevel;

  for (dit.begin(); dit.ok(); ++dit)
    {
#ifndef NDEBUG
      FArrayBox& thisPhi = a_phi[dit()];
      FArrayBox& thisLPhi = a_lofPhi[dit()];
      const Box& thisBox = m_grids.get(dit());

      CH_assert(thisPhi.box().contains(grow(thisBox,1)));
      CH_assert(thisLPhi.box().contains(thisBox));
#endif

      FORT_OLDOPERATORLAP(CHF_FRA(a_lofPhi[dit()]),
                          CHF_CONST_FRA(a_phi[dit()]),
                          CHF_BOX(m_grids.get(dit())),
                          CHF_CONST_REAL(dx));
    }
}

/***********************/
// evaluate operator, inhomogeneous C/F boundary conditions
// note that operator applies C/F interpolation and copy bc's
// homogeneous physical bcs
/***********************/
void PoissonOp::applyOpH(LevelData<FArrayBox>& a_phi,
                         LevelData<FArrayBox>& a_lofPhi)
{
  CH_assert(isDefined());
  CH_assert (a_phi.nComp() == a_lofPhi.nComp());
  // CH_assert (a_phi.nComp() == m_ncomp);

  homogeneousCFInterp(a_phi);

  DataIterator dit = a_phi.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      m_domghostbc.applyHomogeneousBCs(a_phi[dit()],
                                       m_domain,
                                       m_dxLevel);
    }

  // fill in intersection of ghostcells and phi's boxes
  a_phi.exchange(a_phi.interval(), m_exchangeCopier);

  Real dx = m_dxLevel;

  for (dit.begin(); dit.ok(); ++dit)
    {

#ifndef NDEBUG
      FArrayBox& thisPhi = a_phi[dit()];
      FArrayBox& thisLPhi = a_lofPhi[dit()];
      const Box& thisBox = m_grids.get(dit());

      CH_assert(thisPhi.box().contains(grow(thisBox,1)));
      CH_assert(thisLPhi.box().contains(thisBox));
#endif

      FORT_OLDOPERATORLAP(CHF_FRA(a_lofPhi[dit()]),
                          CHF_CONST_FRA(a_phi[dit()]),
                          CHF_BOX(m_grids.get(dit())),
                          CHF_CONST_REAL(dx));
    }
}

/***********************/
// evaluate operator, homogeneous C/F boundary conditions
// note -- grids for phi need not correspond to bx
// inhomogeneous physical boundary conditions
/***********************/
void PoissonOp::applyOpHcfIphys(LevelData<FArrayBox>& a_phi,
                                LevelData<FArrayBox>& a_lofPhi)
{
  CH_assert(isDefined());
  CH_assert (a_phi.nComp() == a_lofPhi.nComp());
  CH_assert (a_phi.nComp() == m_ncomp);

  homogeneousCFInterp(a_phi);

  DataIterator dit = a_phi.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      m_domghostbc.applyInhomogeneousBCs(a_phi[dit()],
                                         m_domain,
                                         m_dxLevel);
    }

  // fill in intersection of ghostcells and phi's boxes
  a_phi.exchange(a_phi.interval(), m_exchangeCopier);

  Real dx = m_dxLevel;

  for (dit.begin(); dit.ok(); ++dit)
    {

#ifndef NDEBUG
      FArrayBox& thisPhi = a_phi[dit()];
      FArrayBox& thisLPhi = a_lofPhi[dit()];
      const Box& thisBox = m_grids.get(dit());

      CH_assert(thisPhi.box().contains(grow(thisBox,1)));
      CH_assert(thisLPhi.box().contains(thisBox));
#endif

      FORT_OLDOPERATORLAP(CHF_FRA(a_lofPhi[dit()]),
                          CHF_CONST_FRA(a_phi[dit()]),
                          CHF_BOX(m_grids.get(dit())),
                          CHF_CONST_REAL(dx));
    }
}

/***********************/
// This calls whatever has been set as the bottom solver.
/***********************/
void PoissonOp::bottomSmoother(LevelData<FArrayBox>&       a_phi,
                               const LevelData<FArrayBox>& a_rhs)
{
  CH_assert(isDefined());
  CH_assert(m_bottomSmootherPtr != NULL);

  m_bottomSmootherPtr->doBottomSmooth(a_phi, a_rhs, this);
}

/**
   get flux( == flux at THIS level)
   The fluxes live on the cell faces with direction dir.
   Fluxes are computed for all interior edges of data.
   The flux fab is resized inside the routine.
*/
void PoissonOp::getFlux(FArrayBox&       a_fineFlux,
                        const FArrayBox& a_data,
                        const DataIndex& a_datInd,
                        int              a_dir)
{
  CH_TIME("PoissonOp::getFlux");
  CH_assert(isDefined());
  CH_assert(a_dir >= 0);
  CH_assert(a_dir <  SpaceDim);
  CH_assert(!a_data.box().isEmpty());

  Box edgebox = surroundingNodes(a_data.box(),a_dir);
  edgebox.grow(a_dir, -1);

  // if this fails, the data box was too small (one cell wide, in fact)
  CH_assert(!edgebox.isEmpty());

  a_fineFlux.resize(edgebox, a_data.nComp());

  FORT_OLDGETFLUX(CHF_FRA(a_fineFlux),
                  CHF_CONST_FRA(a_data),
                  CHF_BOX(edgebox),
                  CHF_CONST_REAL(m_dxLevel),
                  CHF_CONST_INT(a_dir));
  /*
  BoxIterator bit(edgebox);

  for ( bit.begin(); bit.ok(); bit.next())
    {
      IntVect iv = bit();
      IntVect shiftiv = BASISV(a_dir);
      IntVect ivlo = iv - shiftiv;
      IntVect ivhi = iv;

      CH_assert(a_data.box().contains(ivlo));
      CH_assert(a_data.box().contains(ivhi));

      for (int ivar = 0; ivar < a_data.nComp(); ivar++)
        {
          Real phihi = a_data(ivhi,ivar);
          Real philo = a_data(ivlo,ivar);
          Real gradphi =(phihi-philo)/m_dxLevel;

          a_fineFlux(iv,ivar) = gradphi;
        }
    }
  */
}

void
PoissonOp::setConvergenceMetric(Real a_metric, int a_comp)
{
  m_bottomSmootherPtr->setConvergenceMetric(a_metric, a_comp);
}

/// number of GSRB's to do in the LevelPreconditioner call
void
PoissonOp::numPrecondSmooth(int a_numGSRB)
{
  m_numPrecondSmooth = a_numGSRB;
}

#include "NamespaceFooter.H"
