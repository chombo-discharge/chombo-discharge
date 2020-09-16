#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "ParabolaIF.H"

#include "NamespaceHeader.H"

ParabolaIF::ParabolaIF(const Real    &     a_a,
                       const RealVect&     a_point,
                       const bool    &     a_inside)
  :BaseIF()
{
  m_a      = a_a;     
  m_point  = a_point; 
  m_inside = a_inside;

}

Real 
ParabolaIF::
value(const IndexTM<int,SpaceDim> & a_partialDerivative,
      const IndexTM<Real,SpaceDim>& a_point) const
{
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      if(a_partialDerivative[idir] < 0)
        {
          MayDay::Error("invalid partial derivative");
        }
    }
  int order = a_partialDerivative.sum();

  Real retval = LARGEREALVAL;

  if (order == 0)
    {
      retval = value(a_point);
    }
  else if (order == 1)
    {
      if(a_partialDerivative[0] == 1)
        retval = 2*m_a*(a_point[0] - m_point[0]);
      else if (a_partialDerivative[1] == 1)
        retval = -1.0;
      else
        retval = 0.0;
    }
  else if (order == 2)
    {
      if(a_partialDerivative[0] == 2)
        retval = 2*m_a;
      else
        retval = 0.0;
    }
  else
    {
      retval = 0.0;
    }

  if(!m_inside)
    retval = -retval;

  //  pout() << "retval2  = " << retval << endl;
  return retval;

}

Real 
ParabolaIF::
value(const RealVect& a_point) const
{
  IndexTM<Real,GLOBALDIM> pt;

  //check does SpaceDim = GLOBALDIM
  if (GLOBALDIM==3 && SpaceDim==2)
    {
      MayDay::Abort("HyperPlaneIF should be wrapped in ReferenceHeightIF when GLOBALDIM==3 and SpaceDim==2");
    }
  else
    {
      for (int idir = 0; idir < SpaceDim; ++idir)
        {
          pt[idir] = a_point[idir];
        }
    }

  //  pout() << "value 1  = " << value(pt) << endl;
  return value(pt);
}

Real 
ParabolaIF::
value(const IndexTM<Real,SpaceDim>& a_point) const
{
  Real x = a_point[0];
  Real y = a_point[1];
  Real x0 = m_point[0];
  Real y0 = m_point[1];
  Real retval = m_a*(x-x0)*(x-x0) + y0 - y;
  if(!m_inside)
    retval = -retval;
  //  pout() << "retval  = " << retval << endl;
  return retval;
}

BaseIF* 
ParabolaIF::
newImplicitFunction() const
{
  ParabolaIF* parabolaPtr = new ParabolaIF(m_a,m_point,m_inside);
  return static_cast<BaseIF*>(parabolaPtr);
}

#include "NamespaceFooter.H"
