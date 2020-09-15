#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "NewCoordSys.H"
#include "AnalyticCS.H"

#include "NamespaceHeader.H"

//-----------------------------------------------------------------------
AnalyticCS::
AnalyticCS(const RealVect& a_dx,
           RefCountedPtr<VectorFunction> a_X,
           RefCountedPtr<VectorFunction> a_Xi,
           RefCountedPtr<TensorFunction> a_J):
  NewFourthOrderCoordSys(),
  m_X(a_X),
  m_Xi(a_Xi),
  m_J(a_J),
  m_J0(),
  m_Xi0()
{
  CH_assert(!a_X.isNull());
  CH_assert(!a_Xi.isNull());
  CH_assert(!a_J.isNull());
  m_dx = a_dx;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
AnalyticCS::
~AnalyticCS()
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
RealVect
AnalyticCS::
realCoord(const RealVect& a_Xi) const
{
  return (*m_X)(a_Xi);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
RealVect
AnalyticCS::
mappedCoord(const RealVect& a_x) const
{
  return (*m_Xi)(a_x);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
Real
AnalyticCS::
dXdXi(const RealVect& a_Xi, int a_dirX, int a_dirXi) const
{
  // Have we already computed this Jacobian?
  if (a_Xi != m_Xi0)
  {
    m_Xi0 = a_Xi;
    m_J0 = (*m_J)(a_Xi);
  }
  return m_J0(a_dirX, a_dirXi);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
AnalyticCSFactory::
AnalyticCSFactory(RefCountedPtr<VectorFunction> a_X,
                  RefCountedPtr<VectorFunction> a_Xi,
                  RefCountedPtr<TensorFunction> a_J):
  NewCoordSysFactory(),
  m_X(a_X),
  m_Xi(a_Xi),
  m_J(a_J)
{
  CH_assert(!a_X.isNull());
  CH_assert(!a_Xi.isNull());
  CH_assert(!a_J.isNull());
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
NewCoordSys*
AnalyticCSFactory::
getCoordSys(const ProblemDomain& a_domain,
            const RealVect& a_dx) const
{
  AnalyticCS* newCSPtr = new AnalyticCS(a_dx, m_X, m_Xi, m_J);
  return dynamic_cast<NewCoordSys*>(newCSPtr);
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"
