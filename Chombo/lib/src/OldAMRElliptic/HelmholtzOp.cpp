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

#include "LayoutIterator.H"
#include "DotProduct.H"
#include "BoxIterator.H"

#include "NoOpSmoother.H"
#include "BiCGStabSmoother.H"

#include "PoissonOpF_F.H"
#include "HelmholtzOp.H"
#include "HelmholtzOpF_F.H"
#include "NamespaceHeader.H"

// has full define function been called?
bool HelmholtzOp::isDefined() const
{
  return ((m_hcoeffDefined)&&
          (m_isDefined && m_isBCDefined));
}

/***********************/
// default constructor
/***********************/
HelmholtzOp::HelmholtzOp()
{
  setDefaultValues();
}

/***********************/
/***********************/
void HelmholtzOp::setDomainGhostBC(const DomainGhostBC& a_dombcin)
{
  m_domghostbc = a_dombcin;
  m_isBCDefined = true;
}

void HelmholtzOp::setBottomSmoother(const BaseBottomSmoother& a_bottomSmoother)
{
  if (m_bottomSmootherPtr != NULL)
    {
      delete m_bottomSmootherPtr;
    }

  m_bottomSmootherPtr = a_bottomSmoother.new_bottomSmoother();
}

void HelmholtzOp::define(const DisjointBoxLayout& a_grids,
                         const DisjointBoxLayout* a_baseBAPtr,
                         Real                     a_dxLevel,
                         int                      a_refRatio,
                         const Box&               a_domain,
                         bool                     a_homogeneousOnly,
                         int                      a_ncomp)
{
  ProblemDomain physdomain(a_domain);

  define(a_grids, a_baseBAPtr, a_dxLevel, a_refRatio,
         physdomain, a_homogeneousOnly, a_ncomp);
}

void HelmholtzOp::define(const DisjointBoxLayout& a_grids,
                         const DisjointBoxLayout* a_baseBAPtr,
                         Real                     a_dxLevel,
                         int                      a_refRatio,
                         const ProblemDomain&     a_domain,
                         bool                     a_homogeneousOnly,
                         int                      a_ncomp)
{
  clearMemory();

  m_isDefined = true;

  CH_assert(!a_domain.isEmpty());

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
                       a_refRatio, m_ncomp,  m_domain);
    }

  if (a_baseBAPtr != NULL)
    {
      m_baseBA = *a_baseBAPtr;
    }

  DataIterator lit = m_grids.dataIterator();
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      m_loCFIVS[idir].define(m_grids);
      m_hiCFIVS[idir].define(m_grids);

      for (lit.reset(); lit.ok(); ++lit)
        {
          m_loCFIVS[idir][lit()].define(m_domain, m_grids.get(lit()),
                                        m_grids, idir,Side::Lo);
          m_hiCFIVS[idir][lit()].define(m_domain, m_grids.get(lit()),
                                        m_grids, idir,Side::Hi);
        }
    }
}

// a_reftoFine is the refinement between this operator and a_opfine
void HelmholtzOp::define(const LevelOp* a_opfine,
                         int            a_reftoFine)
{
  clearMemory();

  const HelmholtzOp* opfineptr = dynamic_cast<const HelmholtzOp*>(a_opfine);
  if (opfineptr == NULL)
    {
      MayDay::Error("HelmholtzOp::define: casting failed");
    }

  const HelmholtzOp& opfine = *opfineptr;

  CH_assert(opfine.isDefined());
  CH_assert(a_reftoFine > 0);

  m_isDefined = true;
  m_ihcfiEnabled = false;

  setDomainGhostBC(opfine.m_domghostbc);

  m_dxLevel = (opfine.m_dxLevel)*a_reftoFine;

  m_domain = coarsen(opfine.m_domain, a_reftoFine);
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

      for (lit.reset(); lit.ok(); ++lit)
        {
          m_loCFIVS[idir][lit()].define(m_domain, m_grids.get(lit()),
                                        m_grids, idir,Side::Lo);
          m_hiCFIVS[idir][lit()].define(m_domain, m_grids.get(lit()),
                                        m_grids, idir,Side::Hi);
        }
    }
}

// virtual constructor workaround
LevelOp* HelmholtzOp::new_levelop() const
{
  CH_assert(m_hcoeffDefined);
  CH_assert(m_isBCDefined);

  HelmholtzOp* oshPtr = new HelmholtzOp();
  if (oshPtr == NULL)
    {
      MayDay::Error("Out of Memory in HelmholtzOp::new_levelop");
    }

  oshPtr->setAlpha(m_alphaCoeff);
  oshPtr->setHelmCoeff(m_helmCoeff);
  oshPtr->setDomainGhostBC(m_domghostbc);

  if (m_bottomSmootherPtr != NULL)
    {
      oshPtr->setBottomSmoother(*m_bottomSmootherPtr);
    }

  return static_cast<LevelOp*>(oshPtr);
}

/***********************/
/***********************/
void HelmholtzOp::clearMemory()
{
  m_exchangeCopier.clear();
  // leave bottom smoother untouched here because we want it to
  // persist through define statement if necessary
}

/***********************/
/***********************/
void HelmholtzOp::setDefaultValues()
{
  m_dxLevel = -1.0;
  m_refRatio = -1;

  m_isBCDefined = false;
  m_ihcfiEnabled = false;
  m_hcoeffDefined = false;

  // default is to use BiCGStab as bottom smoother
  m_bottomSmootherPtr = new BiCGStabSmoother;
  // m_bottomSmootherPtr = new NoOpSmoother;

  // alpha defaults to 1 (basic Helmholtz eqn)
  m_alphaCoeff = 1.0;

  m_isDefined = false;
}

// set helmholtz coefficient
void HelmholtzOp::setHelmCoeff(Real a_helmcoeff)
{
  m_helmCoeff = a_helmcoeff;
  m_hcoeffDefined = true;
}

// scale helmholtz coefficient
void HelmholtzOp::scaleHelmCoeff(Real a_scale)
{
  scaleBeta(a_scale);
}

void HelmholtzOp::scaleBeta(Real a_scale)
{
  // this doesn't make any sense unless the coefficient has
  // already been defined.
  CH_assert (m_hcoeffDefined);

  m_helmCoeff *= a_scale;
}

// set alpha coefficient
void HelmholtzOp::setAlpha(Real a_alpha)
{
  m_alphaCoeff = a_alpha;
}

// rescale alpha in (alpha*I + beta*laplacian)
void HelmholtzOp::scaleAlpha(Real a_scale)
{
  m_alphaCoeff = a_scale*m_alphaCoeff;
}

/***********************/
// apply coarse-fine boundary conditions -- assume that phi grids
// are grown by one
/***********************/
void HelmholtzOp::CFInterp(LevelData<FArrayBox>&       a_phif,
                           const LevelData<FArrayBox>& a_phic)
{
  CH_assert(isDefined());
  CH_assert(m_ihcfiEnabled);
  CH_assert(m_quadCFI.isDefined());
  CH_assert( a_phif.ghostVect() >= IntVect::Unit);

  m_quadCFI.coarseFineInterp(a_phif, a_phic);
}

/***********************/
// does homogeneous coarse/fine interpolation
/***********************/
void HelmholtzOp::homogeneousCFInterp(LevelData<FArrayBox>& a_phif)
{
  CH_assert(isDefined());
  CH_assert( a_phif.ghostVect() >= IntVect::Unit);

  DataIterator dit = a_phif.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
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
void HelmholtzOp::homogeneousCFInterp(LevelData<FArrayBox>& a_phif,
                                      const DataIndex&      a_datInd,
                                      int                   a_idir,
                                      Side::LoHiSide        a_hiorlo)
{
  CH_assert(isDefined());
  CH_assert((a_idir >= 0) && (a_idir  < SpaceDim));
  CH_assert((a_hiorlo == Side::Lo )||(a_hiorlo == Side::Hi ));
  CH_assert( a_phif.ghostVect() >= IntVect::Unit);
  CH_assert( m_ncomp == a_phif.nComp());

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

      FORT_OLDINTERPHOMO(CHF_FRA(a_phif[a_datInd]),
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
          // int ihilo = sign(a_hiorlo);
          // Box phistarbox = interpIVS.minBox();
          // phistarbox.shift(a_idir, ihilo);
          // FArrayBox phistar(phistarbox, m_ncomp);
          // //hence the homogeneous...
          // phistar.setVal(0.);

          // given phistar, interpolate on fine ivs
          interpOnIVSHomo(a_phif, // phistar,
                          a_datInd, a_idir, a_hiorlo,
                          interpIVS);
        }
    }
}

// void HelmholtzOp::interpOnIVS(LevelData<FArrayBox>& a_phif,
//                               const FArrayBox&      a_phistar,
//                               const DataIndex&      a_datInd,
//                               const int             a_idir,
//                               const Side::LoHiSide  a_hiorlo,
//                               const IntVectSet&     a_interpIVS)
// {
//   CH_assert(isDefined());
//   CH_assert((a_idir >= 0) && (a_idir  < SpaceDim));
//   CH_assert((a_hiorlo == Side::Lo )||(a_hiorlo == Side::Hi ));
//   CH_assert( a_phif.ghostVect() >= IntVect::Unit);
//
//   IVSIterator fineIVSit(a_interpIVS);
//   FArrayBox& a_phi = a_phif[a_datInd];
//
//   CH_assert (m_ncomp = a_phi.nComp());
//
//   int ihilo = sign(a_hiorlo);
//   for (fineIVSit.begin(); fineIVSit.ok(); ++fineIVSit)
//     {
//       IntVect ivf = fineIVSit();
//
//       // quadratic interpolation
//       for (int ivar = 0; ivar < m_ncomp; ivar++)
//         {
//           Real pa =      a_phi(ivf -2*ihilo*BASISV(a_idir), ivar);
//           Real pb =      a_phi(ivf -  ihilo*BASISV(a_idir), ivar);
//           Real pc =  a_phistar(ivf +  ihilo*BASISV(a_idir), ivar);
//
//           // phi = ax**2 + bx + c, x = 0 at pa
//           Real hf = m_dxLevel;
//           Real hc = m_dxCrse;
//
//           Real x1 = hf;
//           Real x2 = 0.5*(3.0*hf+hc);
//
//           Real denom = 1.0-((x1+x2)/x1);
//           Real a = ((pb-pa)/(x1*x1)) - ((pb-pc)/(x1*(x1-x2)));
//           a /= denom;
//
//           Real b = (pb-pc)/(x1-x2) - a*(x1+x2);
//           Real c = pa;
//
//           Real x = 2.0*hf;
//
//           a_phi(ivf,ivar) = a*x*x + b*x + c;
//         } //end loop over components
//     } //end loop over fine intvects
// }

void HelmholtzOp::interpOnIVSHomo(LevelData<FArrayBox>& a_phif,
                                  const DataIndex&      a_datInd,
                                  const int             a_idir,
                                  const Side::LoHiSide  a_hiorlo,
                                  const IntVectSet&     a_interpIVS)
{
  CH_assert(isDefined());
  CH_assert((a_idir >= 0) && (a_idir  < SpaceDim));
  CH_assert((a_hiorlo == Side::Lo ) || (a_hiorlo == Side::Hi ));
  CH_assert(a_phif.ghostVect() >= IntVect::Unit);

  IVSIterator fineIVSit(a_interpIVS);
  FArrayBox& a_phi = a_phif[a_datInd];

  CH_assert(m_ncomp = a_phi.nComp());

  Real x1 = m_dxLevel;
  Real x2 = 0.5*(3.0*m_dxLevel+m_dxCrse);

  Real denom = 1.0-((x1+x2)/x1);
  Real idenom = 1/(denom); // divide is more expensive usually

  Real x = 2.0*m_dxLevel;
  Real xsquared = x*x;

  Real m1 = 1/(x1*x1);
  Real m2 = 1/(x1*(x1-x2));

  Real q1 = 1/(x1-x2);
  Real q2 = x1+x2;

  int ihilo = sign(a_hiorlo);

  for (fineIVSit.begin(); fineIVSit.ok(); ++fineIVSit)
    {
      IntVect ivf = fineIVSit();

      // quadratic interpolation
      for (int ivar = 0; ivar < m_ncomp; ivar++)
        {
          ivf[a_idir] -= 2*ihilo;
          Real pa = a_phi(ivf, ivar);

          ivf[a_idir] += ihilo;
          Real pb = a_phi(ivf, ivar);

          Real a = (pb-pa)*m1 - (pb)*m2;
          a *= idenom;
          Real b = (pb)*q1 - a*q2;

          a_phi(fineIVSit(),ivar) = a*xsquared + b*x + pa;

          // this returns ivf to the original location for use when nvar > 1
          ivf[a_idir]+=ihilo;
        } //end loop over components
    } //end loop over fine intvects
}

/***********************/
/***********************/
HelmholtzOp::~HelmholtzOp()
{
  clearMemory();

  // since clearMemory doesn't delete bottom smoother ptr, need to do it
  // ourselves.
  if (m_bottomSmootherPtr != NULL)
    {
      delete m_bottomSmootherPtr;

      m_bottomSmootherPtr = NULL;
    }
}

/***********************/
// this smoother assumes problem has
// already been put into residual-correction
// form, so that crse/fine BC's are homogeneous
/***********************/
void HelmholtzOp::smooth(LevelData<FArrayBox>&       a_phi,
                         const LevelData<FArrayBox>& a_rhs)
{
  CH_assert(isDefined());
  CH_assert(a_phi.nComp() == a_rhs.nComp());
  CH_assert(a_phi.nComp() == m_ncomp);

  levelGSRB(a_phi,a_rhs);
}

/**************************/
// this preconditioner first initializes phihat to (IA)phihat = rhshat
// (diagonization of L -- A is the matrix version of L)
// then smooths with a couple of passes of levelGSRB
/**************************/
void HelmholtzOp::levelPreconditioner(LevelData<FArrayBox>&       a_phihat,
                                      const LevelData<FArrayBox>& a_rhshat)
{
  // diagonal intialization
  Real mult = 1.0*m_alphaCoeff - m_helmCoeff*m_dxLevel*m_dxLevel/(2*SpaceDim);
  Interval comps = a_phihat.interval();

  CH_assert(a_phihat.nComp() == a_rhshat.nComp());

  // don't need to fill in ghost cells, so we construct
  // our own copier here.
  Copier thisCopier(a_rhshat.getBoxes(), a_phihat.getBoxes(),
                    IntVect::Zero);

  a_rhshat.copyTo(comps, a_phihat, comps, thisCopier);

  DataIterator dit = a_phihat.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      a_phihat[dit()] *= mult;
    }

  levelGSRB(a_phihat, a_rhshat);
  levelGSRB(a_phihat, a_rhshat);
}

/***********************/
//  levelGSRB does GSRB on a level.  it has no knowledge of overlying
//  finer grids, and assumes coarse/fine BC's are homogeneous
/***********************/
void HelmholtzOp::levelGSRB(LevelData<FArrayBox>&       a_phi,
                            const LevelData<FArrayBox>& a_rhs)
{
  CH_assert(isDefined());

  // do first red, then black passes
  for (int whichPass =0; whichPass <= 1; whichPass++)
    {
      // do coarse/fine and copy bc's
      //should be done patch by patch.
      homogeneousCFInterp(a_phi);

      //fill in intersection of ghostcells and a_phi's boxes
      a_phi.exchange(a_phi.interval(), m_exchangeCopier);

      // now step through grids...
      Real dx = m_dxLevel;
      DataIterator dit = a_phi.dataIterator();
      for (dit.reset(); dit.ok(); ++dit)
        {
#ifndef NDEBUG
          FArrayBox& thisPhi = a_phi[dit()];
          const FArrayBox& thisRhs = a_rhs[dit()];
          const Box& thisBox = m_grids.get(dit());
          CH_assert(thisRhs.box().contains(thisBox));
          CH_assert(thisPhi.box().contains(grow(thisBox,1)));
#endif

          // invoke physical BC's where necessary
          // can do this here (w/o a following exchange)
          // because we use a 5-pt stencil
          m_domghostbc.applyHomogeneousBCs(a_phi[dit()],
                                           m_domain,
                                           m_dxLevel);

          FORT_GSRBLEVELHELM(CHF_FRA(a_phi[dit()]),
                             CHF_CONST_FRA(a_rhs[dit()]),
                             CHF_BOX(m_grids.get(dit())),
                             CHF_CONST_REAL(dx),
                             CHF_CONST_REAL(m_alphaCoeff),
                             CHF_CONST_REAL(m_helmCoeff),
                             CHF_CONST_INT(whichPass));
        } // end loop through grids
    } // end loop through red-black
}

/***********************/
// evaluate operator, inhomogeneous C/F boundary conditions
// note that operator applies C/F interpolation and copy bc's
/***********************/
void HelmholtzOp::applyOpI(LevelData<FArrayBox>&       a_phi,
                           const LevelData<FArrayBox>* a_phicPtr,
                           LevelData<FArrayBox>&       a_lofPhi)
{
  CH_assert(isDefined());
  CH_assert(m_ihcfiEnabled);
  CH_assert (a_phi.nComp() == a_lofPhi.nComp());
  CH_assert (a_phi.nComp() == m_ncomp);

  if (a_phicPtr != NULL)
    {
      // apply C/F boundary conditions...
      CFInterp(a_phi,*a_phicPtr);
    }

  // fill in intersection of ghostcells and phi's boxes
  a_phi.exchange(a_phi.interval(), m_exchangeCopier);

  Real dx = m_dxLevel;

  DataIterator dit = a_phi.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
    {
#ifndef NDEBUG
      FArrayBox& thisPhi = a_phi[dit()];
      FArrayBox& thisLPhi = a_lofPhi[dit()];
      const Box& thisBox = m_grids.get(dit());
      CH_assert(thisPhi.box().contains(grow(thisBox,1)));
      CH_assert(thisLPhi.box().contains(thisBox));
#endif

      // can do this here (w/o a following exchange)
      // because we use a 5-pt stencil
      m_domghostbc.applyInhomogeneousBCs(a_phi[dit()],
                                         m_domain,
                                         m_dxLevel);

      FORT_OPERATORHELM(CHF_FRA(a_lofPhi[dit()]),
                        CHF_CONST_FRA(a_phi[dit()]),
                        CHF_BOX(m_grids.get(dit())),
                        CHF_CONST_REAL(dx),
                        CHF_CONST_REAL(m_alphaCoeff),
                        CHF_CONST_REAL(m_helmCoeff));
    }
}

/***********************/
// evaluate operator, inhomogeneous C/F boundary conditions
// note that operator applies C/F interpolation and copy bc's
// homogeneous physical bcs
/***********************/
void HelmholtzOp::applyOpIcfHphys(LevelData<FArrayBox>&       a_phi,
                                  const LevelData<FArrayBox>* a_phicPtr,
                                  LevelData<FArrayBox>&       a_lofPhi)
{
  CH_assert(isDefined());
  CH_assert(m_ihcfiEnabled);
  CH_assert (a_phi.nComp() == a_lofPhi.nComp());
  CH_assert (a_phi.nComp() == m_ncomp);

  if (a_phicPtr != NULL)
    {
      // apply C/F boundary conditions...
      CFInterp(a_phi,*a_phicPtr);
    }

  // fill in intersection of ghostcells and phi's boxes
  a_phi.exchange(a_phi.interval(), m_exchangeCopier);

  Real dx = m_dxLevel;

  DataIterator dit = a_phi.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
    {
#ifndef NDEBUG
      FArrayBox& thisPhi = a_phi[dit()];
      FArrayBox& thisLPhi = a_lofPhi[dit()];
      const Box& thisBox = m_grids.get(dit());
      CH_assert(thisPhi.box().contains(grow(thisBox,1)));
      CH_assert(thisLPhi.box().contains(thisBox));
#endif

      // can do this here (w/o a following exchange)
      // because we use a 5-pt stencil
      m_domghostbc.applyHomogeneousBCs(a_phi[dit()],
                                       m_domain,
                                       m_dxLevel);

      FORT_OPERATORHELM(CHF_FRA(a_lofPhi[dit()]),
                        CHF_CONST_FRA(a_phi[dit()]),
                        CHF_BOX(m_grids.get(dit())),
                        CHF_CONST_REAL(dx),
                        CHF_CONST_REAL(m_alphaCoeff),
                        CHF_CONST_REAL(m_helmCoeff));
    }
}

/***********************/
// evaluate operator, inhomogeneous C/F boundary conditions
// note that operator applies C/F interpolation and copy bc's
// homogeneous physical bcs
/***********************/
void HelmholtzOp::applyOpH(LevelData<FArrayBox>& a_phi,
                           LevelData<FArrayBox>& a_lofPhi)
{
  CH_assert(isDefined());
  CH_assert (a_phi.nComp() == a_lofPhi.nComp());
  CH_assert (a_phi.nComp() == m_ncomp);

  homogeneousCFInterp(a_phi);

  // fill in intersection of ghostcells and phi's boxes
  a_phi.exchange(a_phi.interval(), m_exchangeCopier);

  Real dx = m_dxLevel;

  DataIterator dit = a_phi.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
    {
#ifndef NDEBUG
      FArrayBox& thisPhi = a_phi[dit()];
      FArrayBox& thisLPhi = a_lofPhi[dit()];
      const Box& thisBox = m_grids.get(dit());
      CH_assert(thisPhi.box().contains(grow(thisBox,1)));
      CH_assert(thisLPhi.box().contains(thisBox));
#endif

      // can do this here (w/o a following exchange)
      // because we use a 5-pt stencil
      m_domghostbc.applyHomogeneousBCs(a_phi[dit()],
                                       m_domain,
                                       m_dxLevel);

      FORT_OPERATORHELM(CHF_FRA(a_lofPhi[dit()]),
                        CHF_CONST_FRA(a_phi[dit()]),
                        CHF_BOX(m_grids.get(dit())),
                        CHF_CONST_REAL(dx),
                        CHF_CONST_REAL(m_alphaCoeff),
                        CHF_CONST_REAL(m_helmCoeff));
    }
}

/***********************/
// evaluate operator, homogeneous C/F boundary conditions
// note -- grids for phi need not correspond to bx
// inhomogeneous physical boundary conditions
/***********************/
void HelmholtzOp::applyOpHcfIphys(LevelData<FArrayBox>& a_phi,
                                  LevelData<FArrayBox>& a_lofPhi)
{
  CH_assert(isDefined());
  CH_assert (a_phi.nComp() == a_lofPhi.nComp());
  CH_assert (a_phi.nComp() == m_ncomp);

  homogeneousCFInterp(a_phi);

  // fill in intersection of ghostcells and phi's boxes
  a_phi.exchange(a_phi.interval(), m_exchangeCopier);

  Real dx = m_dxLevel;

  DataIterator dit = a_phi.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
    {
#ifndef NDEBUG
      FArrayBox& thisPhi = a_phi[dit()];
      FArrayBox& thisLPhi = a_lofPhi[dit()];
      const Box& thisBox = m_grids.get(dit());
      CH_assert(thisPhi.box().contains(grow(thisBox,1)));
      CH_assert(thisLPhi.box().contains(thisBox));
#endif

      // can do this here (w/o a following exchange)
      // because we use a 5-pt stencil
      m_domghostbc.applyInhomogeneousBCs(a_phi[dit()],
                                         m_domain,
                                         m_dxLevel);

      FORT_OPERATORHELM(CHF_FRA(a_lofPhi[dit()]),
                        CHF_CONST_FRA(a_phi[dit()]),
                        CHF_BOX(m_grids.get(dit())),
                        CHF_CONST_REAL(dx),
                        CHF_CONST_REAL(m_alphaCoeff),
                        CHF_CONST_REAL(m_helmCoeff));
    }
}

/***********************/
// This calls whatever has been set as the bottom solver.
/***********************/
void HelmholtzOp::bottomSmoother(LevelData<FArrayBox>&       a_phi,
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
void HelmholtzOp::getFlux(FArrayBox&       a_fineFlux,
                          const FArrayBox& a_data,
                          const DataIndex& a_datInd,
                          int              a_dir)
{
  CH_assert(isDefined());
  CH_assert(a_dir >= 0);
  CH_assert(a_dir  < SpaceDim);
  CH_assert(!a_data.box().isEmpty());

  Box edgebox = surroundingNodes(a_data.box(),a_dir);
  edgebox.grow(a_dir, -1);

  CH_assert(!edgebox.isEmpty());

  a_fineFlux.resize(edgebox, a_data.nComp());

  BoxIterator bit(edgebox);
  for (bit.begin(); bit.ok(); bit.next())
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

          a_fineFlux(iv,ivar) = m_helmCoeff*gradphi;
        }
    }
}

void
HelmholtzOp::setConvergenceMetric(Real a_metric, int a_comp)
{
  m_bottomSmootherPtr->setConvergenceMetric(a_metric, a_comp);
}
#include "NamespaceFooter.H"
