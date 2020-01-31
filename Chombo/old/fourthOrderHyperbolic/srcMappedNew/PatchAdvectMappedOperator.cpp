#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cstdio>

#include "PatchAdvectMappedOperator.H"
#include "MOLAdvectionPhysics.H"

#include "NamespaceHeader.H"

//////////////////////////////////////////////////////////////////////////////
// Constructor - set up some defaults
PatchAdvectMappedOperator::PatchAdvectMappedOperator() : PatchMappedConsOperator()
{
  m_advVelAvgPtr = NULL;
  m_advVelCenPtr = NULL;
  m_advVelTransformedPtr = NULL;
}

//////////////////////////////////////////////////////////////////////////////
// Destructor - free up storage
PatchAdvectMappedOperator::~PatchAdvectMappedOperator()
{
}


//////////////////////////////////////////////////////////////////////////////
void
PatchAdvectMappedOperator::getAllFluxes(FluxBox&        a_FfaceAvg,
                                        FluxBox&        a_FfaceCen,
                                        const FluxBox&  a_WfaceAvg,
                                        const FluxBox&  a_WfaceCen)
{
  // Need to cast non-const in order to alias
  FluxBox& WfaceAvg = (FluxBox&) a_WfaceAvg;
  FluxBox& WfaceCen = (FluxBox&) a_WfaceCen;
  CH_assert(m_isCurrentBoxSet);
  MOLAdvectionPhysics* advectionPhysics = (MOLAdvectionPhysics*) m_molPhysics;
  Box bx1inDomain = grow(m_currentBox, 1);
  bx1inDomain &= m_domain;
  for (int idir = 0; idir < SpaceDim; idir++)
    { // idir specifies which faces:  either x or y or z
      // Call setVelocityFab instead of setVelocityFlux
      // because we set all flux components, but all on idir-faces.
      Box faceBox1(bx1inDomain);
      faceBox1.surroundingNodes(idir);
      advectionPhysics->setVelocityFab(&(*m_advVelAvgPtr)[idir]);
      getDirFluxes(a_FfaceAvg[idir], WfaceAvg[idir], faceBox1);

      Box faceBox0(m_currentBox);
      faceBox0.surroundingNodes(idir);
      advectionPhysics->setVelocityFab(&(*m_advVelCenPtr)[idir]);
      getDirFluxes(a_FfaceCen[idir], WfaceCen[idir], faceBox0);
    }
}


//////////////////////////////////////////////////////////////////////////////
void
PatchAdvectMappedOperator::preRiemann(FArrayBox&  a_WLeft,
                                      FArrayBox&  a_WRight,
                                      int         a_dir,
                                      const Box&  a_box)
{
  // This function doesn't change a_WLeft or a_WRight;
  // instead, it transforms the advection velocity.
  FArrayBox& advVelFab = (*m_advVelAvgPtr)[a_dir];
  CH_assert(advVelFab.box().contains(a_box));
  setBasisVectors(a_box, a_dir);
  // If we transformed m_advVelAvgPtr (instead of making a copy)
  // and sent it to m_molPhysics, then this class would still have
  // the transformed m_advVelAvgPtr, and that might mess things up later.
  m_advVelTransformedPtr = new FArrayBox(a_box, SpaceDim);
  m_advVelTransformedPtr->copy(advVelFab);
  forwardBasisTransform(*m_advVelTransformedPtr);
  MOLAdvectionPhysics* advectionPhysics =
    (MOLAdvectionPhysics*) m_molPhysics;
  advectionPhysics->setVelocityFab(m_advVelTransformedPtr);
}


//////////////////////////////////////////////////////////////////////////////
void
PatchAdvectMappedOperator::postRiemann(FArrayBox&  a_Wface,
                                       int         a_dir,
                                       const Box&  a_box)
{
  unsetBasisVectors();
  delete m_advVelTransformedPtr;
}


//////////////////////////////////////////////////////////////////////////////
void
PatchAdvectMappedOperator::setAdvVelAvg(const FluxBox* a_advVelAvgPtr)
{
  m_advVelAvgPtr = (FluxBox*) a_advVelAvgPtr;
}

//////////////////////////////////////////////////////////////////////////////
void
PatchAdvectMappedOperator::setAdvVelCen(const FluxBox* a_advVelCenPtr)
{
  m_advVelCenPtr = (FluxBox*) a_advVelCenPtr;
}

//////////////////////////////////////////////////////////////////////////////
void
PatchAdvectMappedOperator::deconvolveAdvVel()
{
  m_advVelAvgPtr->copy(*m_advVelCenPtr);
  m_util.deconvolveFace(*m_advVelAvgPtr, *m_advVelCenPtr);
}

#include "NamespaceFooter.H"
