#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// GaussianAdvectFun.cpp
// petermc, 4 Feb 2011

// #include <math.h>
#include "GaussianAdvectFun.H"
#include "BoxIterator.H"
#include "FunCylindricalF_F.H"

#include "NamespaceHeader.H"

// ---------------------------------------------------------
// default constructor
GaussianAdvectFun::GaussianAdvectFun()
{
  //  setDefaultValues();
}


// ---------------------------------------------------------
void GaussianAdvectFun::define(Real a_r0,
                               RealVect a_center,
                               RealVect a_x0,
                               Real a_initialTime,
                               RealVect a_rotationCenter,
                               Real a_omega)
{
  //  m_ibc = new GaussianAdvectMultiMappedIBC;
  m_r0 = a_r0;
  m_rr0powInv2 = 1. / (m_r0 * m_r0);
  m_rr0powInvD = 1.;
  for (int idir = 0; idir < SpaceDim; idir++) m_rr0powInvD /= m_r0;
  m_center = a_center;
  m_x0 = a_x0;
  m_initialTime = a_initialTime;
  m_rotationCenter = a_rotationCenter;
  m_omega = a_omega;
  //  m_ibc->setSolidBodyRotation(m_rotationCenter, m_omega);
  //  m_ibc->setParams(m_r0, m_center, m_x0, m_initialTime);
  //  m_ibc->setTime(0.);
  m_isDefined = true;
}


// ---------------------------------------------------------
// destructor
GaussianAdvectFun::~GaussianAdvectFun()
{
  //  delete m_ibc;
}


// ---------------------------------------------------------
void GaussianAdvectFun::setCoordSys(MultiBlockCoordSys* a_coordSysPtr)
{
  m_coordSysPtr = a_coordSysPtr;
  //  m_ibc->setCoordSys(a_coordSysPtr);
}

// ---------------------------------------------------------
void GaussianAdvectFun::setFunctionFromPhysical(FArrayBox&         a_funFab,
                                                const Box&         a_bx,
                                                const FArrayBox&   a_coordsFab)
{
  CH_assert(isDefined());
  int ncomp = 1;
  // m_ibc->pointVal(a_funFab, a_bx, true);
  // Box of current grid
  //Box intersectBox(a_phi.box());
  Box bxShrunken = grow(a_bx, -4);
  int blockNumber = m_coordSysPtr->whichBlock(bxShrunken);
  const NewCoordSys* coordSysCurrentPtr =
    m_coordSysPtr->getCoordSys(blockNumber);

  // compute initial values

  // compute current center
  RealVect center = m_center;
  // amount of rotation
  Real twoPi = 8.0*atan(1.0);
  Real thetaOffset = twoPi*m_omega*(m_initialTime);
  RealVect localCenter = m_center - m_rotationCenter;
  Real radCenter = sqrt(D_TERM(localCenter[0]*localCenter[0],
                               +localCenter[1]*localCenter[1],
                               +0));
  Real oldTheta = acos(localCenter[0]/radCenter);
  Real newTheta = oldTheta + thetaOffset;
  D_TERM(center[0] = radCenter*cos(newTheta);,
         center[1] = radCenter*sin(newTheta);,
         center[2] = center[2];)

    center += m_rotationCenter;

  // 0.5 for cell-centered, 0 for node-centering
  RealVect meshOffset(0.5*RealVect::Unit);
  for (int dir=0; dir<SpaceDim; dir++)
    {
      if (a_bx.type(dir) == IndexType::NODE) meshOffset[dir] = 0.0;
    }

  FArrayBox XiFab(a_bx, SpaceDim);
  FArrayBox XFab(a_bx, SpaceDim);
  coordSysCurrentPtr->getCenterMappedCoordinates(XiFab, a_bx);
  coordSysCurrentPtr->realCoord(XFab, XiFab, a_bx);

  BoxIterator bit(a_bx);
  for (bit.begin();bit.ok();++bit)
    {
      IntVect iv = bit();
      // You could replace XFab with a_coordsFab, because they are the same.
      RealVect X = RealVect(D_DECL(XFab(iv, 0), XFab(iv, 1), XFab(iv, 2)));
      Real dist = 0.;
      for (int idir = 0;idir < SpaceDim; idir++)
        {
          Real loc = m_x0[idir] + X[idir] - center[idir];
          //          if (m_domain.isPeriodic(idir)  )
          //            {
          //              loc = min(abs(loc), abs(loc + 1.0));
          //              loc = min(abs(loc), abs(loc - 1.0));
          //            }
          dist = dist + loc*loc;
        }
      for (int icomp = 0;icomp < ncomp;icomp++)
        {
          a_funFab(iv, icomp) = exp(-dist * m_rr0powInv2) * m_rr0powInvD;
        }
    }
  bool includeJ = false;
  if (includeJ)
    {
      FArrayBox JFab(a_bx, 1);
      coordSysCurrentPtr->pointwiseJ(JFab, XiFab, a_bx);
      for (int icomp = 0;icomp < ncomp;icomp++)
        {
          a_funFab.mult(JFab, a_bx, 0, icomp);
        }
    }
}


// ---------------------------------------------------------
Real GaussianAdvectFun::setFunctionMax(Real   a_bxWidth,
                                       Real   a_outerRadius)
{
  CH_assert(isDefined());
  Real maxfun = 1.0 / (m_r0 * m_r0);
  return maxfun;
}

#include "NamespaceFooter.H"
