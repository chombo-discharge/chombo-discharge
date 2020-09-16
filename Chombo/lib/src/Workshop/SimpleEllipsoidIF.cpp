#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "SimpleEllipsoidIF.H"

#include "NamespaceHeader.H"
/***/
Real 
SimpleEllipsoidIF::
value(const IndexTM<int ,SpaceDim> & a_partialDerivative,
      const IndexTM<Real,SpaceDim> & a_point) const
{
  int order = a_partialDerivative.sum();
  Real retval = LARGEREALVAL;

  RealVect X;
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      X[idir] = a_point[idir] - m_X0[idir];
    }

  if (order == 0)
    {
      retval = value(a_point);
    }
  else if (order == 1)
    {
      bool found = false;
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          if(a_partialDerivative[idir] == 1)
            {
              found = true;
              retval = 2*X[idir]/m_A2[idir];
            }
        }
      if(!found) MayDay::Error("logic error");
    }
  else if (order == 2)
    {
      
      bool found = false;
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          if(a_partialDerivative[idir] == 2)
            {
              found = true;
              retval = 2/m_A2[idir];
            }
        }
      if(!found) //mixed deriv
        {
          retval = 0;
        }
    }
  else
    {
      retval = 0.0;
    }

  if(!m_inside)
    retval = -retval;
  return retval;
}
/***/
Real 
SimpleEllipsoidIF::
value(const RealVect& a_point) const
{
  IndexTM<Real,GLOBALDIM> pt;

  for (int idir = 0; idir < SpaceDim; ++idir)
    {
      pt[idir] = a_point[idir];
    }

  return value(pt);
}
/***/
Real 
SimpleEllipsoidIF::
value(const IndexTM<Real,SpaceDim>& a_point) const
{
  
  RealVect X;
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      X[idir] = a_point[idir] - m_X0[idir];
    }

  Real retval = 0;
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      retval += (X[idir]*X[idir])/m_A2[idir];
    }

  retval -= m_R*m_R;
  if(!m_inside)
    retval = -retval;
  return retval;
}
/***/
#include "NamespaceFooter.H"
