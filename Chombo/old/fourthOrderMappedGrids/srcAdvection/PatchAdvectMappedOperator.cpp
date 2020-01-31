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
  m_advVelFacePtr = NULL;
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
  CH_TIME("PatchAdvectMappedOperator::getAllFluxes");
  CH_assert(m_isCurrentBoxSet);
  MOLAdvectionPhysics* advectionPhysics = (MOLAdvectionPhysics*) m_molPhysics;
  Box bx1inDomain = grow(m_currentBox, 1);
  bx1inDomain &= m_domain;
  for (int idir = 0; idir < SpaceDim; idir++)
    { // idir specifies which faces:  either x or y or z
      // Call setVelocityFab instead of setVelocityFlux
      // because we set all flux components, but all on idir-faces.
      advectionPhysics->setVelocityFab(&(*m_advVelFacePtr)[idir]);

      // a_FfaceAvg[idir] need not be high-order, because we are
      // using it only to take its laplacian,
      // in PatchConsOperator::getNormalFlux().
      Box faceBox1(bx1inDomain);
      faceBox1.surroundingNodes(idir);
      getDirFluxes(a_FfaceAvg[idir], a_WfaceAvg[idir], faceBox1);

      Box faceBox0(m_currentBox);
      faceBox0.surroundingNodes(idir);
      getDirFluxes(a_FfaceCen[idir], a_WfaceCen[idir], faceBox0);
    }
}


//////////////////////////////////////////////////////////////////////////////
void
PatchAdvectMappedOperator::preRiemann(FArrayBox&  a_WLeft,
                                      FArrayBox&  a_WRight,
                                      int         a_dir,
                                      const Box&  a_box)
{
  // No need to transform the advection velocity.
  FArrayBox& advVelFab = (*m_advVelFacePtr)[a_dir];
  CH_assert(advVelFab.box().contains(a_box));
  MOLAdvectionPhysics* advectionPhysics =
    (MOLAdvectionPhysics*) m_molPhysics;
  advectionPhysics->setVelocityFab(&advVelFab);
}


//////////////////////////////////////////////////////////////////////////////
void
PatchAdvectMappedOperator::postRiemann(FArrayBox&  a_Wface,
                                       int         a_dir,
                                       const Box&  a_box)
{
}


//////////////////////////////////////////////////////////////////////////////
void
PatchAdvectMappedOperator::setAdvVelFace(const FluxBox* a_advVelFacePtr)
{
  m_advVelFacePtr = (FluxBox*) a_advVelFacePtr;
}

#include "NamespaceFooter.H"
