#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "SingleBlockCSAdaptor.H"
#include "CartesianBlockCS.H"

#include "NamespaceHeader.H"

//-----------------------------------------------------------------------
SingleBlockCSAdaptor::
SingleBlockCSAdaptor(const NewCoordSys* const a_coordSys,
                     const ProblemDomain& a_problemDomain):
  MultiBlockCoordSys()
{
  CH_assert(a_coordSys != NULL);

  m_coordSysVect.resize(1);
  m_coordSysVect[0] = const_cast<NewCoordSys*>(a_coordSys);
  m_gotCoordSysVect = true;

  m_mappingBlocks.resize(1);
  m_mappingBlocks[0] = a_problemDomain.domainBox();
  m_gotMappingBlocks = true;

  m_boundaries.resize(1);
  Tuple<BlockBoundary, 2*SpaceDim>& blockBoundaries = m_boundaries[0];
  for (int idir = 0; idir < SpaceDim; idir++)
  {
    blockBoundaries[idir].define(0);
    blockBoundaries[idir + SpaceDim].define(0);
  }
  m_gotBoundaries = true;

  // Initialize block transformations
  initializeBlockTransformations();
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
SingleBlockCSAdaptor::
~SingleBlockCSAdaptor()
{
  delete m_coordSysVect[0];
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
SingleBlockCSAdaptor::
blockRemapping(RealVect& a_xi_valid,
               int& a_n_valid,
               const RealVect& a_xiSrc,
               int a_nSrc) const
{
  CH_assert(a_nSrc == 0);
  a_n_valid = 0;
  RealVect X = m_coordSysVect[0]->realCoord(a_xiSrc);
  for (int d = 0; d < SpaceDim; ++d)
    a_xi_valid[d] = X[d];
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
SingleBlockCSAdaptorFactory::
SingleBlockCSAdaptorFactory(const NewCoordSysFactory* const a_coordSysFactory):
  m_coordSysFactory(a_coordSysFactory)
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
SingleBlockCSAdaptorFactory::
~SingleBlockCSAdaptorFactory()
{
  delete m_coordSysFactory;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
MultiBlockCoordSys*
SingleBlockCSAdaptorFactory::
getCoordSys(const ProblemDomain& a_levelDomain,
            const RealVect& a_dx) const
{
  SingleBlockCSAdaptor* coordSysPtr = new SingleBlockCSAdaptor(m_coordSysFactory->getCoordSys(a_levelDomain, a_dx), a_levelDomain);
  return (static_cast <MultiBlockCoordSys*> (coordSysPtr));
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"
