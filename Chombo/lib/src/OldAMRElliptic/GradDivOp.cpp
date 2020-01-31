#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "GradDivOp.H"
#include "GradDivOpF_F.H"
#include "DotProduct.H"
#include "LayoutIterator.H"
#include "BoxIterator.H"
#include "BiCGStabSmoother.H"
#include "NamespaceHeader.H"

void
GradDivOp::setDomainGhostBC(const DomainGhostBC& a_dombcIn)
{
  m_domghostbc = a_dombcIn;
  m_isBCDefined = true;
}

void
GradDivOp::setTanGradBC(const DomainGhostBC& a_dombcIn)
{
  m_tangradbc = a_dombcIn;
  m_isGradBCDefined = true;
}

bool
GradDivOp::isDefined() const
{
  return (m_isDefined && m_isBCDefined && m_isGradBCDefined);
}

/***********************/
// return to undefined state
/***********************/
void
GradDivOp::setDefaultValues()
{
  m_dxLevel = -1.0;
  m_refRatio = -1;
  m_isDefined = false;
  m_ihcfiEnabled = false;
  m_isBCDefined = false;
  m_isGradBCDefined = false;

}

/***********************/
// return to undefined state
/***********************/
void
GradDivOp::clearMemory()
{
  m_dxLevel = -1.0;
  m_refRatio = -1;
  m_isDefined = false;
  m_ihcfiEnabled = false;
  m_exchangeCopier.clear();
  // don't touch bottom smoother here, since we want it to be able to
  // pass through the define statement unaltered.
}
/***********************/
// default constructor
/***********************/
GradDivOp::GradDivOp()
{
  setDefaultValues();
}

/***********************/
/***********************/

void
GradDivOp::define(
                  const DisjointBoxLayout& a_grids,
                  const DisjointBoxLayout* a_baseBAPtr,
                  Real  a_dxLevel,
                  int a_refRatio,
                  const Box& a_domain,
                  bool a_homogeneousOnly,
                  int a_ncomp)
{
  ProblemDomain probdomain(a_domain);
  define(a_grids, a_baseBAPtr, a_dxLevel, a_refRatio, probdomain,
         a_homogeneousOnly, a_ncomp);
}



void
GradDivOp::define(
                  const DisjointBoxLayout& a_grids,
                  const DisjointBoxLayout* a_baseBAPtr,
                  Real  a_dxLevel,
                  int a_refRatio,
                  const ProblemDomain& a_domain,
                  bool a_homogeneousOnly,
                  int a_ncomp)
{
  clearMemory();
  m_isDefined = true;
  CH_assert(!a_domain.isEmpty());
  CH_assert(a_grids.checkPeriodic(a_domain));

  m_dxLevel = a_dxLevel;
  m_domain = a_domain;
  m_grids = a_grids;
  m_exchangeCopier.define(a_grids, a_grids, IntVect::Unit);
  m_refRatio = a_refRatio;
  m_dxCrse = m_refRatio*m_dxLevel;
  m_ihcfiEnabled = (!a_homogeneousOnly);
  m_ncomp = a_ncomp;
  if (m_ihcfiEnabled)
    m_quadCFI.define(a_grids,a_baseBAPtr, a_dxLevel,
                     a_refRatio, m_ncomp, m_domain);
  if (a_baseBAPtr != NULL)
    {
      m_baseBA = *a_baseBAPtr;
    }

  // allocate storage for tangential gradients
  m_tanGrad.define(a_grids, m_ncomp*SpaceDim, IntVect::Unit);


  DataIterator lit = m_grids.dataIterator();
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      m_loCFIVS[idir].define(m_grids);
      m_hiCFIVS[idir].define(m_grids);
      m_loTanStencilSets[idir].define(m_grids);
      m_hiTanStencilSets[idir].define(m_grids);

      for (lit.begin(); lit.ok(); ++lit)
        {
          m_loCFIVS[idir][lit()].define(m_domain, m_grids.get(lit()),
                                        m_grids, idir,Side::Lo);
          m_hiCFIVS[idir][lit()].define(m_domain, m_grids.get(lit()),
                                        m_grids, idir,Side::Hi);

          const IntVectSet& fineIVSlo = m_loCFIVS[idir][lit()].getFineIVS();
          m_loTanStencilSets[idir][lit()].define(fineIVSlo,
                                                 m_domain, idir);

          const IntVectSet& fineIVShi = m_hiCFIVS[idir][lit()].getFineIVS();
          m_hiTanStencilSets[idir][lit()].define(fineIVShi,
                                                 m_domain, idir);

          if (idir==0)
            m_tanGrad[lit()].setVal(666.66);
        }
    }

}

// a_reftofine is the refinement between this operator and a_opfine
void GradDivOp::define(const GradDivOp* a_opfinePtr,
                       int a_reftoFine)
{
  clearMemory();
  const GradDivOp& opfine = *a_opfinePtr;

  CH_assert(opfine.isDefined());
  CH_assert(a_reftoFine > 0);
  m_isDefined = true;
  m_ihcfiEnabled = false;
  setDomainGhostBC(opfine.m_domghostbc);

  m_dxLevel = (opfine.m_dxLevel)*a_reftoFine;
  m_domain = coarsen(opfine.m_domain,a_reftoFine);
  coarsen(m_grids, opfine.m_grids, a_reftoFine);
  m_exchangeCopier.define(m_grids, m_grids, IntVect::Unit);
  m_dxCrse = opfine.m_dxCrse;
  m_ncomp = opfine.m_ncomp;

  // define storage for tangential gradients
  m_tanGrad.define(m_grids, m_ncomp*SpaceDim, IntVect::Unit);

  //this refratio is the ratio between this level and the
  //next level in the amr hierarchy
  m_refRatio = -1;
  //don't need to define the cfinterp since inhomogeneous cfinterpolaton
  //is never done here.  We shall leave m_quadCFI and base_ba undefined and
  //and just leave the refratio -1.  We DO need fineinterp_ivs for doing
  //homogeneous cfinterpolation
  DataIterator lit = m_grids.dataIterator();
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      m_loCFIVS[idir].define(m_grids);
      m_hiCFIVS[idir].define(m_grids);
      m_loTanStencilSets[idir].define(m_grids);
      m_hiTanStencilSets[idir].define(m_grids);
      for (lit.begin(); lit.ok(); ++lit)
        {
          m_loCFIVS[idir][lit()].define(m_domain, m_grids.get(lit()),
                                        m_grids, idir,Side::Lo);
          m_hiCFIVS[idir][lit()].define(m_domain, m_grids.get(lit()),
                                        m_grids, idir,Side::Hi);

          const IntVectSet& fineIVSlo = m_loCFIVS[idir][lit()].getFineIVS();
          m_loTanStencilSets[idir][lit()].define(fineIVSlo,
                                                 m_domain, idir);

          const IntVectSet& fineIVShi = m_hiCFIVS[idir][lit()].getFineIVS();
          m_hiTanStencilSets[idir][lit()].define(fineIVShi,
                                                 m_domain, idir);
        }
    }
}


/***********************/
//  apply coarse-fine boundary conditions -- assume that U grids
//  are grown by one
/***********************/
void
GradDivOp::CFInterp(LevelData<FArrayBox>& a_Uf,
                    const LevelData<FArrayBox>& a_Uc)
{
  CH_assert(isDefined());
  CH_assert(m_ihcfiEnabled);
  CH_assert( a_Uf.ghostVect() >= IntVect::Unit);
  CH_assert(m_quadCFI.isDefined());

  // need to do an exchange to ensure that tangential derivatives
  // are computed correctly
  a_Uf.exchange(a_Uf.interval());

  // store the tangential gradients for later use...
  m_quadCFI.coarseFineInterp(a_Uf, m_tanGrad, a_Uc);
}


/***********************/
// does homogeneous coarse/fine interpolation
/***********************/
void
GradDivOp::homogeneousCFInterp(LevelData<FArrayBox>& a_Uf)
{
  CH_assert(isDefined());
  CH_assert( a_Uf.ghostVect() >= IntVect::Unit);

  // need to do this to be sure that tangential derivatives are computed
  // correctly
  a_Uf.exchange(a_Uf.interval());
  DataIterator dit = a_Uf.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const DataIndex& datInd = dit();

      // first fill in cells for U
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          SideIterator sit;
          for (sit.begin(); sit.ok(); sit.next())
            {
              homogeneousCFInterpU(a_Uf,datInd,idir,sit());
            }
        }

      // now fill in tangential gradient cells
      for (int idir = 0; idir<SpaceDim; idir++)
        {
          SideIterator sit;
          for (sit.begin(); sit.ok(); sit.next())
            {
              homogeneousCFInterpTanGrad(m_tanGrad, a_Uf,
                                         datInd,idir,sit());
            }
        }
    }
}

/***********************/
// does homogeneous coarse/fine interpolation for U
/***********************/
void
GradDivOp::homogeneousCFInterpU(LevelData<FArrayBox>& a_Uf,
                                const DataIndex& a_datInd,
                                int a_idir,
                                Side::LoHiSide a_hiorlo)
{
  CH_assert(isDefined());
  CH_assert((a_idir >= 0) && (a_idir  < SpaceDim));
  CH_assert((a_hiorlo == Side::Lo )||(a_hiorlo == Side::Hi ));
  CH_assert( a_Uf.ghostVect() >= IntVect::Unit);

  const CFIVS* cfivs_ptr = NULL;
  if (a_hiorlo == Side::Lo)
    cfivs_ptr = &m_loCFIVS[a_idir][a_datInd];
  else
    cfivs_ptr = &m_hiCFIVS[a_idir][a_datInd];

  const IntVectSet& interp_ivs = cfivs_ptr->getFineIVS();
  if (!interp_ivs.isEmpty())
    {
      int ihilo = sign(a_hiorlo);
      Box Ustarbox = interp_ivs.minBox();
      Ustarbox.shift(a_idir, ihilo);
      FArrayBox Ustar(Ustarbox, m_ncomp);
      //hence the homogeneous...
      Ustar.setVal(0.);

      //given Ustar, interpolate on fine ivs to fill ghost cells for U
      interpUOnIVS(a_Uf, Ustar, a_datInd, a_idir, a_hiorlo,
                     interp_ivs);
    }
}


void
GradDivOp::interpUOnIVS(LevelData<FArrayBox>& a_Uf,
                        const FArrayBox& a_Ustar,
                        const DataIndex& a_datInd,
                        const int a_idir,
                        const Side::LoHiSide a_hiorlo,
                        const IntVectSet& a_interpIVS)
{
  CH_assert(isDefined());
  CH_assert((a_idir >= 0) && (a_idir  < SpaceDim));
  CH_assert((a_hiorlo == Side::Lo )||(a_hiorlo == Side::Hi ));
  CH_assert( a_Uf.ghostVect() >= IntVect::Unit);
  IVSIterator fine_ivsit(a_interpIVS);
  FArrayBox& a_U = a_Uf[a_datInd];
  CH_assert (m_ncomp = a_U.nComp());

  Real x1 = m_dxLevel;
  Real x2 = 0.5*(3.*m_dxLevel+m_dxCrse);
  Real denom = 1.0-((x1+x2)/x1);
  Real idenom = 1.0/(denom); // divide is more expensive usually
  Real x = 2.*m_dxLevel;
  Real xsquared = x*x;

  Real m1 = 1/(x1*x1);
  Real m2 = 1/(x1*(x1-x2));

  Real q1 = 1/(x1-x2);
  Real q2 = x1+x2;

  int ihilo = sign(a_hiorlo);
  IntVect ai = -2*ihilo*BASISV(a_idir);
  IntVect bi = -  ihilo*BASISV(a_idir);
  IntVect ci =    ihilo*BASISV(a_idir);


  for (fine_ivsit.begin(); fine_ivsit.ok(); ++fine_ivsit)
    {
      IntVect ivf = fine_ivsit();
      // quadratic interpolation
      for (int ivar = 0; ivar < m_ncomp; ivar++)
        {
          Real pa =      a_U(ivf + ai, ivar);
          Real pb =      a_U(ivf + bi, ivar);
          Real pc =  a_Ustar(ivf + ci, ivar);
          //U = ax**2 + bx + c, x = 0 at pa
          Real a = (pb-pa)*m1 - (pb-pc)*m2;
          a *= idenom;
          Real b = (pb-pc)*q1 - a*q2;
          Real c = pa;
          a_U(ivf,ivar) = a*xsquared + b*x + c;
        } //end loop over components
    } //end loop over fine intvects
}

/***********************/
// does homogeneous coarse/fine interpolation for tangential gradients
/***********************/
void
GradDivOp::homogeneousCFInterpTanGrad(LevelData<FArrayBox>& a_tanGrad,
                                      const LevelData<FArrayBox>& a_U,
                                      const DataIndex& a_DatInd,
                                      int a_idir,
                                      Side::LoHiSide a_hiorlo)
{

  // sanity checks
  CH_assert (isDefined());
  CH_assert ((a_idir >= 0) && (a_idir < SpaceDim));
  CH_assert ((a_hiorlo == Side::Lo) || (a_hiorlo == Side::Hi));
  CH_assert ( a_tanGrad.ghostVect() >= IntVect::Unit);

  const TensorFineStencilSet* cfstencil_ptr = NULL;
  if (a_hiorlo == Side::Lo)
    cfstencil_ptr = &m_loTanStencilSets[a_idir][a_DatInd];
  else
    cfstencil_ptr = &m_hiTanStencilSets[a_idir][a_DatInd];


  const FArrayBox& U = a_U[a_DatInd];
  FArrayBox& tanGrad = a_tanGrad[a_DatInd];

  // loop over gradient directions
  for (int gradDir = 0; gradDir<SpaceDim; gradDir++)
    {
      if (gradDir != a_idir)
        {

          // first do centered stencil
          const IntVectSet& centeredIVS =
            cfstencil_ptr->getCenteredStencilSet(gradDir);

          int ihilo = sign(a_hiorlo);

          if (!centeredIVS.isEmpty())
            {
              // do centered computation
              IVSIterator cntrd_ivs(centeredIVS);
              // want to average fine-grid gradient with coarse
              // grid gradient, which is 0 (which is where the
              // extra factor of one-half comes from)

              Real gradMult = (0.5/m_dxLevel);
              for (cntrd_ivs.begin(); cntrd_ivs.ok(); ++cntrd_ivs)
                {
                  IntVect ivf = cntrd_ivs();
                  IntVect fineULoc = ivf - ihilo*BASISV(a_idir);
                  IntVect fineULoc2 = fineULoc - ihilo*BASISV(a_idir);
                  // loop over variables
                  for (int ivar = 0; ivar<a_U.nComp(); ivar++)
                    {
                      Real fineHi = U(fineULoc2+BASISV(gradDir),ivar);
                      Real fineLo = U(fineULoc2-BASISV(gradDir),ivar);
                      Real fineGrada = gradMult*(fineHi-fineLo);

                      fineHi = U(fineULoc+BASISV(gradDir),ivar);
                      fineLo = U(fineULoc-BASISV(gradDir),ivar);
                      Real fineGradb = gradMult*(fineHi-fineLo);

                      int gradComp = SpaceDim*ivar + gradDir;

                      Real h = m_dxLevel;
                      Real a = (2./h/h)*(fineGrada*(m_refRatio+1.0)
                                         -fineGradb*(m_refRatio+3.0))/
                        (m_refRatio*m_refRatio + 4*m_refRatio + 3.0);
                      Real b = (fineGradb-fineGrada)/h - a*h;
                      Real c = fineGrada;
                      Real x = 2.*h;
                      tanGrad(ivf,gradComp) = a*x*x + b*x + c;

                    }
                } // end loop over centered difference cells
            } // end if there are centered cells

          // now do forward-difference cells
          const IntVectSet& forwardIVS =
            cfstencil_ptr->getForwardStencilSet(gradDir);

          if (!forwardIVS.isEmpty())
            {
              // do forward-difference computations
              IVSIterator fwd_ivs(forwardIVS);
              // set up multipliers for gradient; since we want to average
              // fine-grid gradient with coarse-grid gradient (which is 0),
              // include an extra factor of one-half here.
              Real mult0 = -1.5/m_dxLevel;
              Real mult1 = 2.0/m_dxLevel;
              Real mult2 = -0.5/m_dxLevel;

              for (fwd_ivs.begin(); fwd_ivs.ok(); ++fwd_ivs)
                {
                  IntVect ivf = fwd_ivs();
                  IntVect fineULoc = ivf - ihilo*BASISV(a_idir);
                  IntVect fineULoc2 = fineULoc - ihilo*BASISV(a_idir);
                  //now loop overvariables
                  for (int var= 0; var<a_U.nComp(); var++)
                    {
                      Real fine0 = U(fineULoc2,var);
                      Real fine1 = U(fineULoc2+BASISV(gradDir),var);
                      Real fine2 = U(fineULoc2+2*BASISV(gradDir),var);
                      Real fineGrada = mult0*fine0 +mult1*fine1 +mult2*fine2;

                      fine0 = U(fineULoc,var);
                      fine1 = U(fineULoc+BASISV(gradDir),var);
                      fine2 = U(fineULoc+2*BASISV(gradDir),var);
                      Real fineGradb = mult0*fine0 +mult1*fine1 +mult2*fine2;

                      int gradComp = var*SpaceDim + gradDir;
                      // now compute gradient

                      Real h = m_dxLevel;
                      Real a = (2./h/h)*(fineGrada*(m_refRatio+1.0)
                                         -fineGradb*(m_refRatio+3.0))/
                        (m_refRatio*m_refRatio + 4*m_refRatio + 3.0);
                      Real b = (fineGradb-fineGrada)/h - a*h;
                      Real c = fineGrada;
                      Real x = 2.*h;

                      tanGrad(ivf,gradComp) = a*x*x + b*x + c;

#if 0
                      tanGrad(ivf,gradComp)=0.5*(mult0*fine0 + mult1*fine1
                                                 + mult2*fine2);
#endif
                    } // end loop over variables
                } // end loop over forward-difference locations
            } // end if there are forward-difference cells

          // now do backward-difference cells
          const IntVectSet& backwardIVS =
            cfstencil_ptr->getBackwardStencilSet(gradDir);

          if (!backwardIVS.isEmpty())
            {
              IVSIterator back_ivs(backwardIVS);
              // set up multipliers for gradient -- since we want to average
              // fine-grid gradient with coarse-grid gradient (which is 0),
              // include an extra factor of one-half here.
              Real mult0 = -1.5/m_dxLevel;
              Real mult1 = 2.0/m_dxLevel;
              Real mult2 = -0.5/m_dxLevel;

              for (back_ivs.begin(); back_ivs.ok(); ++back_ivs)
                {
                  IntVect ivf = back_ivs();
                  IntVect fineULoc = ivf - ihilo*BASISV(a_idir);
                  IntVect fineULoc2 = fineULoc - ihilo*BASISV(a_idir);
                  // now loop over variables
                  for (int var=0; var<a_U.nComp(); var++)
                    {
                      Real fine0 = U(fineULoc2,var);
                      Real fine1 = U(fineULoc2-BASISV(gradDir),var);
                      Real fine2 = U(fineULoc2-2*BASISV(gradDir),var);
                      Real fineGrada = mult0*fine0 +mult1*fine1 +mult2*fine2;

                      fine0 = U(fineULoc,var);
                      fine1 = U(fineULoc-BASISV(gradDir),var);
                      fine2 = U(fineULoc-2*BASISV(gradDir),var);
                      Real fineGradb = mult0*fine0 +mult1*fine1 +mult2*fine2;

                      int gradComp = var*SpaceDim + gradDir;

                      Real h = m_dxLevel;
                      Real a = (2./h/h)*(fineGrada*(m_refRatio+1.0)
                                         -fineGradb*(m_refRatio+3.0))/
                        (m_refRatio*m_refRatio + 4*m_refRatio + 3.0);
                      Real b = (fineGradb-fineGrada)/h - a*h;
                      Real c = fineGrada;
                      Real x = 2.*h;
                      tanGrad(ivf,gradComp) = a*x*x + b*x + c;

#if 0
                      tanGrad(ivf,gradComp) = 0.5(mult0*fine0 +mult1*fine1
                                                  +mult2*fine2);
#endif

                    } // end loop over variables
                } // end loop over backward-difference cells
            } // end if there are backward-difference cells


        } // end if gradDir is a tangential direction
    } // end loop over gradient directions

}

/***********************/
/***********************/
GradDivOp::~GradDivOp()
{
  clearMemory();
}


/***********************/
// evaluate operator, inhomogeneous C/F boundary conditions
// note that operator applies C/F interpolation and copy bc's
/***********************/
void
GradDivOp::applyOpI(LevelData<FArrayBox>& a_LofU,
                    LevelData<FArrayBox>& a_U,
                    const LevelData<FArrayBox>* a_UcPtr)

{
  CH_assert(isDefined());
  CH_assert(m_ihcfiEnabled);
  CH_assert (a_U.nComp() == a_LofU.nComp());
  CH_assert (a_U.nComp() == m_ncomp);

  if (a_UcPtr != NULL)
    {
      // apply C/F boundary conditions...
      CFInterp(a_U,*a_UcPtr);
    }

  //fill in intersection of ghostcells and U's boxes
  a_U.exchange(a_U.interval(), m_exchangeCopier);

  DataIterator dit = a_U.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      m_domghostbc.applyInhomogeneousBCs(a_U[dit()],
                                         m_domain,
                                         m_dxLevel);
    }

  //fill in intersection of ghostcells and U's boxes
  // for a 5-point operator, this should not be necessary
  //a_U.exchange(a_U.interval(), m_exchangeCopier);

  LevelData<FluxBox> faceDiv(m_grids, 1);

  // compute grad(U) over interior cells
  computeTanGradInterior(a_U);

  // fill in physical BC's for tangential gradients
  for (dit.begin(); dit.ok(); ++dit)
    {
      m_tangradbc.applyInhomogeneousBCs(m_tanGrad[dit()],
                                        m_domain,
                                        m_dxLevel);
    }

  // now compute div(U) on faces
  computeFaceDiv(faceDiv, a_U, m_tanGrad);


  for (dit.begin(); dit.ok(); ++dit)
    {

      const FArrayBox& thisU = a_U[dit()];
      FArrayBox& thisLU = a_LofU[dit()];
      const Box& thisBox = m_grids.get(dit());
#ifndef NDEBUG
      CH_assert(thisU.box().contains(grow(thisBox,1)));
      CH_assert(thisLU.box().contains(thisBox));
#endif

      thisLU.setVal(0.0);

      // now loop over directions, incrementing LofU
      // with directional component of Laplacian(U) + grad(div(U))
      for (int dir=0; dir<SpaceDim; dir++)
        {
          const FArrayBox& dirDiv = faceDiv[dit()][dir];
          FORT_INCREMENTGRADDIVOP(CHF_FRA(thisLU),
                                  CHF_CONST_FRA(thisU),
                                  CHF_CONST_FRA(dirDiv),
                                  CHF_BOX(thisBox),
                                  CHF_REAL(m_dxLevel),
                                  CHF_INT(dir));
        }
    }
}


/***********************/
// evaluate operator, inhomogeneous C/F boundary conditions
// note that operator applies C/F interpolation and copy bc's
// homogeneous physical bcs
/***********************/
void
GradDivOp::applyOpIcfHphys(LevelData<FArrayBox>& a_LofU,
                           LevelData<FArrayBox>& a_U,
                           const LevelData<FArrayBox>* a_UcPtr)

{
  CH_assert(isDefined());
  CH_assert(m_ihcfiEnabled);
  CH_assert (a_U.nComp() == a_LofU.nComp());
  CH_assert (a_U.nComp() == m_ncomp);

  if (a_UcPtr != NULL)
    {
      // apply C/F boundary conditions...
      CFInterp(a_U,*a_UcPtr);

    }

  //fill in intersection of ghostcells and U's boxes
  a_U.exchange(a_U.interval(), m_exchangeCopier);

  DataIterator dit = a_U.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      m_domghostbc.applyHomogeneousBCs(a_U[dit()],
                                       m_domain,
                                       m_dxLevel);

    }

  //fill in intersection of ghostcells and U's boxes
  // for a 5 point operator, this shouldn't be necessary
  //a_U.exchange(a_U.interval(), m_exchangeCopier);

  LevelData<FluxBox> faceDiv(m_grids, 1);
  // compute tangential grad(U) over interior cells
  computeTanGradInterior(a_U);

  // fill in physical BC's for tangential gradients
  for (dit.begin(); dit.ok(); ++dit)
    {
      m_tangradbc.applyHomogeneousBCs(m_tanGrad[dit()],
                                      m_domain,
                                      m_dxLevel);
    }

  // now compute div(U) on faces
  computeFaceDiv(faceDiv, a_U, m_tanGrad);


  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& thisU = a_U[dit()];
      FArrayBox& thisLU = a_LofU[dit()];
      const Box& thisBox = m_grids.get(dit());
#ifndef NDEBUG
      CH_assert(thisU.box().contains(grow(thisBox,1)));
      CH_assert(thisLU.box().contains(thisBox));
#endif
      thisLU.setVal(0.0);

      // now loop over directions, incrementing LofU
      // with directional component of Laplaician(U) + grad(div(U))
      for (int dir=0; dir<SpaceDim; dir++)
        {
          const FArrayBox& dirDiv = faceDiv[dit()][dir];
          FORT_INCREMENTGRADDIVOP(CHF_FRA(thisLU),
                                  CHF_CONST_FRA(thisU),
                                  CHF_CONST_FRA(dirDiv),
                                  CHF_BOX(thisBox),
                                  CHF_CONST_REAL(m_dxLevel),
                                  CHF_INT(dir));
        }
    }
}


/***********************/
// evaluate operator, inhomogeneous C/F boundary conditions
// note that operator applies C/F interpolation and copy bc's
// homogeneous physical bcs
/***********************/
void
GradDivOp::applyOpH(LevelData<FArrayBox>& a_LofU,
                    LevelData<FArrayBox>& a_U)

{
  CH_assert(isDefined());
  CH_assert (a_U.nComp() == a_LofU.nComp());

  homogeneousCFInterp(a_U);

  DataIterator dit = a_U.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      m_domghostbc.applyHomogeneousBCs(a_U[dit()],
                                       m_domain,
                                       m_dxLevel);
    }

  //fill in intersection of ghostcells and U's boxes
  a_U.exchange(a_U.interval(), m_exchangeCopier);

  LevelData<FluxBox> faceDiv(m_grids, 1);

  // compute grad(U) over interior cells
  computeTanGradInterior(a_U);

  // fill in physical BC's for tangential gradients
  for (dit.begin(); dit.ok(); ++dit)
    {
      m_tangradbc.applyHomogeneousBCs(m_tanGrad[dit()],
                                      m_domain,
                                      m_dxLevel);
    }


  // now compute div(U) on faces
  computeFaceDiv(faceDiv, a_U, m_tanGrad);


  for (dit.begin(); dit.ok(); ++dit)
    {

      FArrayBox& thisU = a_U[dit()];
      FArrayBox& thisLU = a_LofU[dit()];
      const Box& thisBox = m_grids.get(dit());
#ifndef NDEBUG
      CH_assert(thisU.box().contains(grow(thisBox,1)));
      CH_assert(thisLU.box().contains(thisBox));
#endif

      thisLU.setVal(0.0);

      // now loop over directions, incrementing LofU
      // with directional component of Laplacian(U) + grad(div(U))
      for (int dir=0; dir<SpaceDim; dir++)
        {
          const FArrayBox& dirDiv = faceDiv[dit()][dir];
          FORT_INCREMENTGRADDIVOP(CHF_FRA(thisLU),
                                  CHF_CONST_FRA(thisU),
                                  CHF_CONST_FRA(dirDiv),
                                  CHF_BOX(thisBox),
                                  CHF_REAL(m_dxLevel),
                                  CHF_INT(dir));
        }
    }
}


/***********************/
// evaluate operator, homogeneous C/F boundary conditions
// note -- grids for U need not correspond to bx
// inhomogeneous physical boundary conditions
/***********************/
void
GradDivOp::applyOpHcfIphys(LevelData<FArrayBox>& a_LofU,
                           LevelData<FArrayBox>& a_U)

{
  CH_assert(isDefined());
  CH_assert (a_U.nComp() == a_LofU.nComp());
  CH_assert (a_U.nComp() == m_ncomp);

  homogeneousCFInterp(a_U);

  DataIterator dit = a_U.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      m_domghostbc.applyInhomogeneousBCs(a_U[dit()],
                                         m_domain,
                                         m_dxLevel);
    }

  //fill in intersection of ghostcells and U's boxes
  a_U.exchange(a_U.interval(), m_exchangeCopier);

  LevelData<FluxBox> faceDiv(m_grids, 1);
  // compute grad(U) over interior cells
  computeTanGradInterior(a_U);

  // fill in physical BC's for tangential gradients
  for (dit.begin(); dit.ok(); ++dit)
    {
      m_tangradbc.applyInhomogeneousBCs(m_tanGrad[dit()],
                                        m_domain,
                                        m_dxLevel);
    }

  // now compute div(U) on faces
  computeFaceDiv(faceDiv, a_U, m_tanGrad);



  for (dit.begin(); dit.ok(); ++dit)
    {

      FArrayBox& thisU = a_U[dit()];
      FArrayBox& thisLU = a_LofU[dit()];
      const Box& thisBox = m_grids.get(dit());
#ifndef NDEBUG
      CH_assert(thisU.box().contains(grow(thisBox,1)));
      CH_assert(thisLU.box().contains(thisBox));
#endif

      thisLU.setVal(0.0);

      // now loop over directions, incrementing LofU
      // with directional component of Laplacian(U) + grad(div(U))
      for (int dir=0; dir<SpaceDim; dir++)
        {
          const FArrayBox& dirDiv = faceDiv[dit()][dir];
          FORT_INCREMENTGRADDIVOP(CHF_FRA(thisLU),
                                  CHF_CONST_FRA(thisU),
                                  CHF_CONST_FRA(dirDiv),
                                  CHF_BOX(thisBox),
                                  CHF_REAL(m_dxLevel),
                                  CHF_INT(dir));
        }
    }
}

/***********************/
// evaluate operator, inhomogeneous C/F boundary conditions
// note that operator applies C/F interpolation and copy bc's
/***********************/
void
GradDivOp::computeFaceDiv(LevelData<FluxBox>& a_divU,
                          LevelData<FArrayBox>& a_U,
                          const LevelData<FArrayBox>* a_UcPtr)

{
  CH_assert(isDefined());
  CH_assert(m_ihcfiEnabled);
  CH_assert (a_U.nComp() == SpaceDim*a_divU.nComp());
  CH_assert (a_U.nComp() == m_ncomp);

  if (a_UcPtr != NULL)
    {
      // apply C/F boundary conditions...
      CFInterp(a_U,*a_UcPtr);
    }

  //fill in intersection of ghostcells and U's boxes
  a_U.exchange(a_U.interval(), m_exchangeCopier);

  DataIterator dit = a_U.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      m_domghostbc.applyInhomogeneousBCs(a_U[dit()],
                                         m_domain,
                                         m_dxLevel);
    }

  //fill in intersection of ghostcells and U's boxes
  // for a 5-point operator, this should not be necessary
  //a_U.exchange(a_U.interval(), m_exchangeCopier);

  LevelData<FluxBox> faceDiv(m_grids, 1);

  // compute grad(U) over interior cells
  computeTanGradInterior(a_U);

  // fill in physical BC's for tangential gradients
  for (dit.begin(); dit.ok(); ++dit)
    {
      m_tangradbc.applyInhomogeneousBCs(m_tanGrad[dit()],
                                        m_domain,
                                        m_dxLevel);
    }

  // now compute div(U) on faces
  computeFaceDiv(a_divU, a_U, m_tanGrad);


}


/***********************/
// evaluate operator, inhomogeneous C/F boundary conditions
// note that operator applies C/F interpolation and copy bc's
// homogeneous physical bcs
/***********************/
void
GradDivOp::computeFaceDivH(LevelData<FluxBox>& a_divU,
                           LevelData<FArrayBox>& a_U)

{
  CH_assert(isDefined());
  CH_assert (a_U.nComp() == SpaceDim*a_divU.nComp());

  homogeneousCFInterp(a_U);

  DataIterator dit = a_U.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      m_domghostbc.applyHomogeneousBCs(a_U[dit()],
                                       m_domain,
                                       m_dxLevel);
    }

  //fill in intersection of ghostcells and U's boxes
  a_U.exchange(a_U.interval(), m_exchangeCopier);

  LevelData<FluxBox> faceDiv(m_grids, 1);

  // compute grad(U) over interior cells
  computeTanGradInterior(a_U);

  // fill in physical BC's for tangential gradients
  for (dit.begin(); dit.ok(); ++dit)
    {
      m_tangradbc.applyHomogeneousBCs(m_tanGrad[dit()],
                                      m_domain,
                                      m_dxLevel);
    }


  // now compute div(U) on faces
  computeFaceDiv(faceDiv, a_U, m_tanGrad);

}




/**
    get flux( == flux at THIS level)
    The fluxes live on the cell faces with direction dir.
    Fluxes are computed for all interior edges of data.
    The flux fab is resized inside the routine.
*/
void
GradDivOp::getFlux(FArrayBox& a_fineFlux,
                   const FArrayBox& a_data,
                   const DataIndex& a_datInd,
                   int a_dir)
{
  CH_assert(isDefined());
  CH_assert(a_dir >= 0);
  CH_assert(a_dir  < SpaceDim);
  CH_assert(!a_data.box().isEmpty());

  Box edgebox = surroundingNodes(a_data.box(),a_dir);
  edgebox.grow(a_dir, -1);

  // need to be sure we have an updated div(U)
  // only do tangential directions here
  for (int gradDir=0; gradDir<SpaceDim; gradDir++)
    if (gradDir != a_dir)
      computeTanGradInterior(a_data, a_datInd, gradDir);

  // should we apply physical boundary conditions here?
  // guessing no, since we shouldn't ever need fluxes
  // at the physical boundary anyway.  also, don't know
  // a priori whether to do homogeneous or inhomogeneous BC's

  FArrayBox faceDiv(edgebox,1);

  computeFaceDiv(faceDiv, a_data, m_tanGrad[a_datInd], a_datInd,a_dir);


  //if this fails, the data box was too small (one cell wide, in fact)
  CH_assert(!edgebox.isEmpty());

  a_fineFlux.resize(edgebox, a_data.nComp());
  BoxIterator bit(edgebox);
  for (bit.begin(); bit.ok(); bit.next())
    {
      IntVect iv = bit();
      a_fineFlux(iv,a_dir) += faceDiv(iv,0);
    }
}


// note that this assumes that the tangential gradients have already
// been computed (passed in to make dependency explicit).  Also
// assumes that all boundary conditions have been set, including
// appropriate BC's on a_tanGrad!
void
GradDivOp::computeFaceDiv(LevelData<FluxBox>& a_div,
                          const LevelData<FArrayBox>& a_vel,
                          const LevelData<FArrayBox>& a_tanGrad)
{
  // first make certain that components all match up.
  int ncomp = a_div.nComp();
  CH_assert (SpaceDim*ncomp <= a_vel.nComp());
  CH_assert (SpaceDim*SpaceDim*ncomp <= a_tanGrad.nComp());

  // now check grids
  const DisjointBoxLayout& levelGrids = a_div.getBoxes();
  // now check grids (actually, can't do this because it's likely in
  // levelsolver, etc, that the the grids will be equivalent, but
  // generated in different places, so that the == operator will
  // return false even though the boxes are the same)

  // CH_assert (levelGrids == a_vel.getBoxes());
  // CH_assert (levelGrids == a_tanGrad.getBoxes());


  // loop over grids and sum over directions
  DataIterator dit = levelGrids.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
    {
      FluxBox& thisDiv = a_div[dit()];
      const FArrayBox& thisVel = a_vel[dit()];
      const FArrayBox& thisTanGrad = a_tanGrad[dit()];

      // now loop over face directions
      for (int faceDir=0; faceDir<SpaceDim; faceDir++)
        {
          FArrayBox& thisDivDir = thisDiv[faceDir];

          // now call function for single face direction and dataIndex
          computeFaceDiv(thisDivDir, thisVel, thisTanGrad, dit(),faceDir);

        } // end loop over face directions
    } // end loop over grids
}


// computes for a single face direction and dataIndex
void
GradDivOp::computeFaceDiv(FArrayBox& a_div,
                          const FArrayBox& a_vel,
                          const FArrayBox& a_tanGrad,
                          const DataIndex& a_dataInd,
                          int a_faceDir)
{
  int ncomp = a_div.nComp();
  a_div.setVal(0.0);
  IntVect shiftiv = BASISV(a_faceDir);
  Real oneOnDx = 1.0/m_dxLevel;

  const Box& faceDivBox = a_div.box();
  BoxIterator divIt(faceDivBox);

  for (divIt.begin(); divIt.ok(); divIt.next())
    {

      IntVect iv(divIt());
      IntVect ivlo = iv - shiftiv;

      // now loop over directions for divergence computation
      // this could probably stand to go into fortran...
      for (int divDir=0; divDir<SpaceDim; divDir++)
        {
          if (divDir == a_faceDir)
            {
              CH_assert (a_vel.box().contains(ivlo));
              CH_assert (a_vel.box().contains(iv));
              for (int comp=0; comp<ncomp; comp++)
                {
                  int velcomp = SpaceDim*comp + divDir;
                  a_div(iv,comp) += oneOnDx*(a_vel(iv,velcomp)
                                             - a_vel(ivlo,velcomp));
                } // end loop over components
            } // end if a normal direction
          else
            {
              // take average of cell-centered tangential derivatives
              // and increment div
              CH_assert (a_tanGrad.box().contains(ivlo));
              CH_assert (a_tanGrad.box().contains(iv));
              for (int comp=0; comp<ncomp; comp++)
                {
                  int velcomp = SpaceDim*comp + divDir;
                  int gradcomp = SpaceDim*velcomp + divDir;
                  a_div(iv,comp) += 0.5*(a_tanGrad(iv,gradcomp)
                                         + a_tanGrad(ivlo, gradcomp));
                } // end loop over components
            } // end if this is a tangential direction
        }  // end loop over divergence direction components
    } // end loop over box
}


// this computes tangential gradients of a_vel in interior cells
// and also puts extrapolated values in cells outside physical domain
// (is this what we want?)
void
GradDivOp::computeTanGradInterior(const LevelData<FArrayBox>& a_vel)
{
  // first make certain all components match up
  int ncomp = a_vel.nComp();
  CH_assert (SpaceDim*ncomp <= m_tanGrad.nComp());

  const DisjointBoxLayout& levelGrids = a_vel.getBoxes();
  // now check grids (actually, can't do this because it's likely in
  // levelsolver, etc, that the the grids will be equivalent, but
  // generated in different places, so that the == operator will
  // return false even though the boxes are the same)
  // CH_assert (levelGrids == m_tanGrad.getBoxes());

  DataIterator dit = levelGrids.dataIterator();

  for (dit.reset(); dit.ok(); ++dit)
    {
      computeTanGradInterior(a_vel[dit()], dit());
    }

  // end with exchange (don't think I need to do a cornercopier exchange)
  m_tanGrad.exchange(m_tanGrad.interval());

}

void
GradDivOp::computeTanGradInterior(const FArrayBox& a_vel,
                                  const DataIndex& a_datInd)
{
  for (int gradDir=0; gradDir<SpaceDim; gradDir++)
    {
      computeTanGradInterior(a_vel, a_datInd, gradDir);
    }
}

void
GradDivOp::computeTanGradInterior(const FArrayBox& a_vel,
                                  const DataIndex& a_datInd,
                                  int a_gradDir)
{
  // this box wil tell us which cells can use normal stencil and
  // which will use shifted stencil
  Box domainBox(m_domain.domainBox());
  domainBox.grow(1);
  domainBox &= m_domain;

  const Box& gridBox = m_tanGrad.getBoxes()[a_datInd];
  FArrayBox& thisGrad = m_tanGrad[a_datInd];

  IntVect shiftVect = BASISV(a_gradDir);
  Box centeredBox(gridBox);
  centeredBox &= a_vel.box();
  centeredBox.grow(a_gradDir,1);
  centeredBox &= domainBox;
  centeredBox.grow(a_gradDir,-1);

  Real centeredMult = 0.5/m_dxLevel;

  // first go ahead and do centered-difference cells
  BoxIterator centeredIt(centeredBox);
  for (centeredIt.begin(); centeredIt.ok(); centeredIt.next())
    {
      IntVect iv(centeredIt());
      IntVect ivhi = iv + shiftVect;
      IntVect ivlo = iv - shiftVect;
      CH_assert (a_vel.box().contains(ivlo));
      CH_assert (a_vel.box().contains(ivhi));

      for (int comp=0; comp<a_vel.nComp(); comp++)
        {
          // this is for component of divergence, so it's
          // really a normal derivative
          if (comp == a_gradDir)
            {
              int gradComp = comp*SpaceDim + a_gradDir;
              thisGrad(iv,gradComp) = centeredMult*(a_vel(ivhi,comp)
                                                    -a_vel(ivlo,comp));
            } // end if component is a tangential one
        } // end loop over vel components
    } // end loop over centered-difference cells


  // create domain check boxes
  Box domCheckBoxHi = adjCellHi(domainBox,a_gradDir,1);
  domCheckBoxHi.shift(a_gradDir, -1);

  Box domCheckBoxLo = adjCellLo(domainBox,a_gradDir,1);
  domCheckBoxLo.shift(a_gradDir,1);

  if (!(centeredBox == gridBox))
    {
      // we have some one-sided differences to do
      // first check hi
      Box hiBox = adjCellHi(centeredBox,a_gradDir,1);
      hiBox &= domCheckBoxHi;
      if (!hiBox.isEmpty())
        {
          Real mult0 = 1.5/m_dxLevel;
          Real mult1 = -2.0/m_dxLevel;
          Real mult2 = 0.5/m_dxLevel;

          BoxIterator hiSideIt(hiBox);
          for (hiSideIt.begin(); hiSideIt.ok(); hiSideIt.next())
            {
              IntVect iv(hiSideIt());
              IntVect ivlo = iv - shiftVect;
              IntVect ivlower = ivlo - shiftVect;
              // basic assumption here that grid interior
              // is at least 2 cells wide (otherwise i don't
              // know what to do)
              CH_assert (a_vel.box().contains(iv));
              CH_assert (a_vel.box().contains(ivlo));
              CH_assert (a_vel.box().contains(ivlower));

              for (int comp=0; comp<a_vel.nComp(); comp++)
                {
                  // remember, only compute _tangential_ gradients
                  if (comp == a_gradDir)
                    {
                      int gradComp = comp*SpaceDim + a_gradDir;
                      Real vel0 = a_vel(iv,comp);
                      Real vel1 = a_vel(ivlo,comp);
                      Real vel2 = a_vel(ivlower,comp);

                      thisGrad(iv,gradComp) = mult0*vel0
                        + mult1*vel1 + mult2*vel2;

                    } // end if component is tangential
                } // end loop over tangential velocity comps
            } // end loop over hi-side cells
        } // end if there is a hi-side box

      //  now check lo-side
      Box loBox = adjCellLo(centeredBox, a_gradDir,1);
      loBox &= domCheckBoxLo;
      if (!loBox.isEmpty())
        {
          Real mult0 = -1.5/m_dxLevel;
          Real mult1 = 2.0/m_dxLevel;
          Real mult2 = -0.5/m_dxLevel;

          BoxIterator loSideIt(loBox);
          for (loSideIt.begin(); loSideIt.ok(); ++loSideIt)
            {
              IntVect iv(loSideIt());
              IntVect ivhi = iv + shiftVect;
              IntVect ivhigher = ivhi +shiftVect;
              // basic assumption here that grid interior
              // is at least 2 cells wide (otherwise, i really
              // don't know what to do)
              CH_assert (a_vel.box().contains(iv));
              CH_assert (a_vel.box().contains(ivhi));
              CH_assert (a_vel.box().contains(ivhigher));

              for (int comp=0; comp<a_vel.nComp(); comp++)
                {
                  // only compute tangential gradients
                  if (comp == a_gradDir)
                    {
                      int gradComp = comp*SpaceDim + a_gradDir;
                      Real vel0 = a_vel(iv,comp);
                      Real vel1 = a_vel(ivhi,comp);
                      Real vel2 = a_vel(ivhigher, comp);

                      thisGrad(iv,gradComp) = mult0*vel0
                        + mult1*vel1 + mult2*vel2;
                    } // end if component is tangential
                } // end loop over velocity components
            } // end loop over lo-side box
        } // end if there is a lo-side box
    } // end if need to do any one-sided differences

}


const LevelData<FArrayBox>&
GradDivOp::getTanGrad() const
{
  return m_tanGrad;
}


#include "NamespaceFooter.H"
