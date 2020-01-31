#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "TensorOp.H"
#include "QuadCFInterp.H"
#include "TensorOpF_F.H"
#include "DotProduct.H"
#include "LayoutIterator.H"
#include "BoxIterator.H"
#include "BiCGStabSmoother.H"
#include "NamespaceHeader.H"

void
TensorOp::setDomainGhostBC(const DomainGhostBC& a_dombcIn)
{
  m_domghostbc = a_dombcIn;
  m_isBCDefined = true;
}

void
TensorOp::setTanGradBC(const DomainGhostBC& a_dombcIn)
{
  m_tangradbc = a_dombcIn;
  m_isGradBCDefined = true;
}

void
TensorOp::setBottomSmoother(const BaseBottomSmoother& a_bottomSmoother)
{
  if (m_bottom_smoother_ptr != NULL)
    {
      delete m_bottom_smoother_ptr;
    }
  m_bottom_smoother_ptr = a_bottomSmoother.new_bottomSmoother();
}

bool
TensorOp::isDefined() const
{
  return (m_isDefined && m_isBCDefined && m_isGradBCDefined);
}
/***********************/
/***********************/
LevelOp*
TensorOp::new_levelop() const
{
  //only check to see if the boundary conditions are
  //defined.  not necessary to have whole thing defined.
  //the boundary conditions need to be defined because
  //the interface does not know about the boundary
  //condtions.  the solvers call the other define functions
  //in our inimitable two-stage construction
  CH_assert(m_isBCDefined);
  CH_assert(m_isGradBCDefined);
  TensorOp* newop = new TensorOp();
  if (newop == NULL)
    {
      MayDay::Error("Out of Memory in TensorOp::new_levelop");
    }
  newop->setDomainGhostBC(m_domghostbc);
  newop->setTanGradBC(m_tangradbc);
  if (m_bottom_smoother_ptr != NULL)
    {
      newop->setBottomSmoother(*m_bottom_smoother_ptr);
    }

  return static_cast<LevelOp*>(newop);
}
/***********************/
// return to undefined state
/***********************/
void
TensorOp::setDefaultValues()
{
  m_dxLevel = -1.0;
  m_refRatio = -1;
  m_isDefined = false;
  m_ihcfiEnabled = false;
  m_isBCDefined = false;
  m_isGradBCDefined = false;
  // default is to use BiCGStab as bottom smoother
  m_bottom_smoother_ptr = new BiCGStabSmoother;

}

/***********************/
// return to undefined state
/***********************/
void
TensorOp::clearMemory()
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
TensorOp::TensorOp():LevelOp()
{
  setDefaultValues();
}

/***********************/
/***********************/

void
TensorOp::define(
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
TensorOp::define(
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

        }
    }

}

// a_reftofine is the refinement between this operator and a_opfine
void TensorOp::define(
                      const LevelOp* a_opfine,
                      int a_reftoFine)
{
  clearMemory();
  const TensorOp* opfineptr =
    dynamic_cast<const TensorOp*>(a_opfine);
  if (opfineptr == NULL)
    MayDay::Error("TensorOp::define: casting failed");
  const TensorOp& opfine = *opfineptr;
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
  setBottomSmoother(*(opfine.m_bottom_smoother_ptr));

}


/***********************/
//  apply coarse-fine boundary conditions -- assume that phi grids
//  are grown by one
/***********************/
void
TensorOp::CFInterp(LevelData<FArrayBox>& a_phif,
                   const LevelData<FArrayBox>& a_phic)
{
  CH_assert(isDefined());
  CH_assert(m_ihcfiEnabled);
  CH_assert( a_phif.ghostVect() >= IntVect::Unit);
  CH_assert(m_quadCFI.isDefined());

  // need to do an exchange to ensure that tangential derivatives
  // are computed correctly
  a_phif.exchange(a_phif.interval());

  // store the tangential gradients for later use...
  m_quadCFI.coarseFineInterp(a_phif, m_tanGrad, a_phic);
}


/***********************/
// does homogeneous coarse/fine interpolation
/***********************/
void
TensorOp::homogeneousCFInterp(LevelData<FArrayBox>& a_phif)
{
  CH_assert(isDefined());
  CH_assert( a_phif.ghostVect() >= IntVect::Unit);

  QuadCFInterp::homogeneousCFInterp(a_phif, m_tanGrad,
                                    m_loCFIVS, m_hiCFIVS,
                                    m_dxLevel, m_dxCrse, m_ncomp,
                                    m_loTanStencilSets,m_hiTanStencilSets);
}

/***********************/
/***********************/
TensorOp::~TensorOp()
{
  clearMemory();
  if (m_bottom_smoother_ptr != NULL)
    {
      delete m_bottom_smoother_ptr;
      m_bottom_smoother_ptr = NULL;
    }
}

/***********************/
//this smoother assumes problem has
//already been put into residual-correction
//  form, so that crse/fine BC's are homogeneous
/***********************/
void
TensorOp::smooth(LevelData<FArrayBox>& a_phi,
                 const LevelData<FArrayBox>& a_rhs)
{
  CH_assert(isDefined());
  levelGSRB(a_phi,a_rhs);
}


/**************************/
// this preconditioner first initializes phihat to (IA)phihat = rhshat
// (diagonization of L -- A is the matrix version of L)
// then smooths with a couple of passes of levelGSRB
void
TensorOp::levelPreconditioner(LevelData<FArrayBox>& a_phihat,
                              const LevelData<FArrayBox>& a_rhshat)
{
  // for first stab, use same initialization as Poisson's equation
  Real mult = 0.25*m_dxLevel*m_dxLevel;
  Interval comps = a_phihat.interval();
  CH_assert(a_phihat.nComp() == a_rhshat.nComp());

  a_rhshat.copyTo(comps, a_phihat, comps);

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
void
TensorOp::levelGSRB(LevelData<FArrayBox>& a_phi,
                    const LevelData<FArrayBox>& a_rhs)
{
  CH_assert(isDefined());
  CH_assert(a_phi.isDefined());
  CH_assert(a_rhs.isDefined());
  CH_assert(a_phi.ghostVect() >= IntVect::Unit);
  CH_assert(a_phi.nComp() == a_rhs.nComp());

  // local temp storage
  LevelData<FluxBox> faceDiv(m_grids, 1);


  // do first red, then black passes
  for (int whichPass =0; whichPass <= 1; whichPass++)
    {
      // do coarse/fine and copy bc's
      //should be done patch by patch.
      homogeneousCFInterp(a_phi);

      //fill in intersection of ghostcells and a_phi's boxes
      a_phi.exchange(a_phi.interval(), m_exchangeCopier);

      // now step through grids...
      DataIterator dit = a_phi.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          // invoke physical BC's where necessary
          m_domghostbc.applyHomogeneousBCs(a_phi[dit()],
                                           m_domain,
                                           m_dxLevel);
        }

      //fill in intersection of ghostcells and phi's boxes
      // for a 5-point operator, this should not be necessary
      //a_phi.exchange(a_phi.interval(), m_exchangeCopier);

      computeTanGradInterior(a_phi);

      // fill in homogeneous physical BCs for tangential gradients
      for (dit.begin(); dit.ok(); ++dit)
        {
          m_tangradbc.applyHomogeneousBCs(m_tanGrad[dit()],
                                          m_domain,
                                          m_dxLevel);
        }

      computeFaceDiv(faceDiv, a_phi, m_tanGrad);



      for (dit.begin(); dit.ok(); ++dit)
        {

#ifndef NDEBUG
          FArrayBox& thisPhi = a_phi[dit()];
          const FArrayBox& thisRhs = a_rhs[dit()];
          const Box& thisBox = m_grids.get(dit());
          CH_assert(thisRhs.box().contains(thisBox));
          CH_assert(thisPhi.box().contains(grow(thisBox,1)));
#endif


          // first increment with laplacian part
          FORT_GSRBLEVELTENSORLAP(CHF_FRA(a_phi[dit()]),
                                  CHF_CONST_FRA(a_rhs[dit()]),
                                  CHF_BOX(m_grids.get(dit())),
                                  CHF_REAL(m_dxLevel),
                                  CHF_CONST_INT(whichPass));

          // need to loop over direction doing GSRB for
          // grad(div) part
          for (int dir=0; dir<SpaceDim; dir++)
            {
              const FArrayBox& dirDiv = faceDiv[dit()][dir];

              FORT_GSRBLEVELTENSORDIR(CHF_FRA(a_phi[dit()]),
                                      CHF_CONST_FRA(dirDiv),
                                      CHF_BOX(m_grids.get(dit())),
                                      CHF_REAL(m_dxLevel),
                                      CHF_CONST_INT(whichPass),
                                      CHF_CONST_INT(dir));
            } // end loop over directions

        } // end loop through grids

    } // end loop through red-black
}

/***********************/
// evaluate operator, inhomogeneous C/F boundary conditions
// note that operator applies C/F interpolation and copy bc's
/***********************/
void
TensorOp::applyOpI(LevelData<FArrayBox>& a_phi,
                   const LevelData<FArrayBox>* a_phicPtr,
                   LevelData<FArrayBox>& a_lofPhi)
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

  LevelData<FluxBox> faceDiv(m_grids, 1);

  // compute grad(phi) over interior cells
  computeTanGradInterior(a_phi);

  // fill in physical BC's for tangential gradients
  for (dit.begin(); dit.ok(); ++dit)
    {
      m_tangradbc.applyInhomogeneousBCs(m_tanGrad[dit()],
                                        m_domain,
                                        m_dxLevel);
    }

  // now compute div(phi) on faces
  computeFaceDiv(faceDiv, a_phi, m_tanGrad);


  for (dit.begin(); dit.ok(); ++dit)
    {

      const FArrayBox& thisPhi = a_phi[dit()];
      FArrayBox& thisLPhi = a_lofPhi[dit()];
      const Box& thisBox = m_grids.get(dit());
#ifndef NDEBUG
      CH_assert(thisPhi.box().contains(grow(thisBox,1)));
      CH_assert(thisLPhi.box().contains(thisBox));
#endif

      thisLPhi.setVal(0.0);

      // now loop over directions, incrementing LofPhi
      // with directional component of Laplacian(phi) + grad(div(phi))
      for (int dir=0; dir<SpaceDim; dir++)
        {
          const FArrayBox& dirDiv = faceDiv[dit()][dir];
          FORT_INCREMENTTENSOROP(CHF_FRA(thisLPhi),
                                 CHF_CONST_FRA(thisPhi),
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
TensorOp::applyOpIcfHphys(
                          LevelData<FArrayBox>& a_phi,
                          const LevelData<FArrayBox>* a_phicPtr,
                          LevelData<FArrayBox>& a_lofPhi)
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

  //fill in intersection of ghostcells and phi's boxes
  a_phi.exchange(a_phi.interval(), m_exchangeCopier);

  DataIterator dit = a_phi.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      m_domghostbc.applyHomogeneousBCs(a_phi[dit()],
                                       m_domain,
                                       m_dxLevel);

    }

  //fill in intersection of ghostcells and phi's boxes
  // for a 5 point operator, this shouldn't be necessary
  //a_phi.exchange(a_phi.interval(), m_exchangeCopier);

  LevelData<FluxBox> faceDiv(m_grids, 1);
  // compute tangential grad(phi) over interior cells
  computeTanGradInterior(a_phi);

  // fill in physical BC's for tangential gradients
  for (dit.begin(); dit.ok(); ++dit)
    {
      m_tangradbc.applyHomogeneousBCs(m_tanGrad[dit()],
                                      m_domain,
                                      m_dxLevel);
    }

  // now compute div(phi) on faces
  computeFaceDiv(faceDiv, a_phi, m_tanGrad);


  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& thisPhi = a_phi[dit()];
      FArrayBox& thisLPhi = a_lofPhi[dit()];
      const Box& thisBox = m_grids.get(dit());
#ifndef NDEBUG
      CH_assert(thisPhi.box().contains(grow(thisBox,1)));
      CH_assert(thisLPhi.box().contains(thisBox));
#endif
      thisLPhi.setVal(0.0);

      // now loop over directions, incrementing LofPhi
      // with directional component of Laplaician(phi) + grad(div(phi))
      for (int dir=0; dir<SpaceDim; dir++)
        {
          const FArrayBox& dirDiv = faceDiv[dit()][dir];
          FORT_INCREMENTTENSOROP(CHF_FRA(thisLPhi),
                                 CHF_CONST_FRA(thisPhi),
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
TensorOp::applyOpH(LevelData<FArrayBox>& a_phi,
                   LevelData<FArrayBox>& a_lofPhi)
{
  CH_assert(isDefined());
  CH_assert (a_phi.nComp() == a_lofPhi.nComp());

  homogeneousCFInterp(a_phi);

  DataIterator dit = a_phi.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      m_domghostbc.applyHomogeneousBCs(a_phi[dit()],
                                       m_domain,
                                       m_dxLevel);
    }

  //fill in intersection of ghostcells and phi's boxes
  a_phi.exchange(a_phi.interval(), m_exchangeCopier);

  LevelData<FluxBox> faceDiv(m_grids, 1);

  // compute grad(phi) over interior cells
  computeTanGradInterior(a_phi);

  // fill in physical BC's for tangential gradients
  for (dit.begin(); dit.ok(); ++dit)
    {
      m_tangradbc.applyHomogeneousBCs(m_tanGrad[dit()],
                                      m_domain,
                                      m_dxLevel);
    }


  // now compute div(phi) on faces
  computeFaceDiv(faceDiv, a_phi, m_tanGrad);


  for (dit.begin(); dit.ok(); ++dit)
    {

      FArrayBox& thisPhi = a_phi[dit()];
      FArrayBox& thisLPhi = a_lofPhi[dit()];
      const Box& thisBox = m_grids.get(dit());
#ifndef NDEBUG
      CH_assert(thisPhi.box().contains(grow(thisBox,1)));
      CH_assert(thisLPhi.box().contains(thisBox));
#endif

      thisLPhi.setVal(0.0);

      // now loop over directions, incrementing LofPhi
      // with directional component of Laplacian(phi) + grad(div(phi))
      for (int dir=0; dir<SpaceDim; dir++)
        {
          const FArrayBox& dirDiv = faceDiv[dit()][dir];
          FORT_INCREMENTTENSOROP(CHF_FRA(thisLPhi),
                                 CHF_CONST_FRA(thisPhi),
                                 CHF_CONST_FRA(dirDiv),
                                 CHF_BOX(thisBox),
                                 CHF_REAL(m_dxLevel),
                                 CHF_INT(dir));
        }
    }
}


/***********************/
// evaluate operator, homogeneous C/F boundary conditions
// note -- grids for phi need not correspond to bx
// inhomogeneous physical boundary conditions
/***********************/
void
TensorOp::applyOpHcfIphys(LevelData<FArrayBox>& a_phi,
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

  //fill in intersection of ghostcells and phi's boxes
  a_phi.exchange(a_phi.interval(), m_exchangeCopier);

  LevelData<FluxBox> faceDiv(m_grids, 1);
  // compute grad(phi) over interior cells
  computeTanGradInterior(a_phi);

  // fill in physical BC's for tangential gradients
  for (dit.begin(); dit.ok(); ++dit)
    {
      m_tangradbc.applyInhomogeneousBCs(m_tanGrad[dit()],
                                        m_domain,
                                        m_dxLevel);
    }

  // now compute div(phi) on faces
  computeFaceDiv(faceDiv, a_phi, m_tanGrad);



  for (dit.begin(); dit.ok(); ++dit)
    {

      FArrayBox& thisPhi = a_phi[dit()];
      FArrayBox& thisLPhi = a_lofPhi[dit()];
      const Box& thisBox = m_grids.get(dit());
#ifndef NDEBUG
      CH_assert(thisPhi.box().contains(grow(thisBox,1)));
      CH_assert(thisLPhi.box().contains(thisBox));
#endif

      thisLPhi.setVal(0.0);

      // now loop over directions, incrementing LofPhi
      // with directional component of Laplacian(phi) + grad(div(phi))
      for (int dir=0; dir<SpaceDim; dir++)
        {
          const FArrayBox& dirDiv = faceDiv[dit()][dir];
          FORT_INCREMENTTENSOROP(CHF_FRA(thisLPhi),
                                 CHF_CONST_FRA(thisPhi),
                                 CHF_CONST_FRA(dirDiv),
                                 CHF_BOX(thisBox),
                                 CHF_REAL(m_dxLevel),
                                 CHF_INT(dir));
        }
    }
}



/***********************/
// This calls whatever has been set as the bottom solver.
/***********************/


void
TensorOp::bottomSmoother(LevelData<FArrayBox>& a_phi,
                         const LevelData<FArrayBox>& a_rhs)
{
  CH_assert(isDefined());
  CH_assert(m_bottom_smoother_ptr != NULL);
  m_bottom_smoother_ptr->doBottomSmooth(a_phi, a_rhs, this);
}


/**
    get flux( == flux at THIS level)
    The fluxes live on the cell faces with direction dir.
    Fluxes are computed for all interior edges of data.
    The flux fab is resized inside the routine.
*/
void
TensorOp::getFlux(FArrayBox& a_fineFlux,
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

  // need to be sure we have an updated div(phi)
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

          if (ivar == a_dir)
            a_fineFlux(iv,ivar) += faceDiv(iv,0);

        }
    }
}


void
TensorOp::setConvergenceMetric(Real a_metric, int a_comp)
{
  m_bottom_smoother_ptr->setConvergenceMetric(a_metric, a_comp);
}


// note that this assumes that the tangential gradients have already
// been computed (passed in to make dependency explicit).  Also
// assumes that all boundary conditions have been set, including
// appropriate BC's on a_tanGrad!
void
TensorOp::computeFaceDiv(LevelData<FluxBox>& a_div,
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
TensorOp::computeFaceDiv(FArrayBox& a_div,
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
      for (int divDir=0; divDir<SpaceDim; divDir++)
        {
          if (divDir == a_faceDir)
            {
              CH_assert (a_vel.box().contains(ivlo));
              CH_assert (a_vel.box().contains(iv));
              for (int comp=0; comp<ncomp; comp++)
                {
                  int velcomp = TensorCFInterp::gradIndex(comp,divDir);
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
                  int velcomp = TensorCFInterp::gradIndex(comp,divDir);
                  int gradcomp = TensorCFInterp::gradIndex(velcomp,divDir);
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
TensorOp::computeTanGradInterior(const LevelData<FArrayBox>& a_vel)
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
TensorOp::computeTanGradInterior(const FArrayBox& a_vel,
                                 const DataIndex& a_datInd)
{
  for (int gradDir=0; gradDir<SpaceDim; gradDir++)
    {
      computeTanGradInterior(a_vel, a_datInd, gradDir);
    }
}

void
TensorOp::computeTanGradInterior(const FArrayBox& a_vel,
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
              int gradComp = TensorCFInterp::gradIndex(comp,a_gradDir);
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
                      int gradComp = TensorCFInterp::gradIndex(comp,a_gradDir);
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
                      int gradComp = TensorCFInterp::gradIndex(comp,a_gradDir);
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
TensorOp::getTanGrad() const
{
  return m_tanGrad;
}


#include "NamespaceFooter.H"
