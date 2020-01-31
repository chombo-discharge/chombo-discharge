#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "SineIF.H"
#include "CONSTANTS.H"
#include "NamespaceHeader.H"

SineIF::
SineIF(const RealVect & a_A,
       const RealVect & a_point,
       const RealVect & a_F,
       const bool     & a_inside)
  :BaseIF()
{
  m_A      = a_A     ;
  m_point  = a_point ;
  m_F      = a_F     ;
  m_inside = a_inside;
  m_factor = PI*m_A*m_F;
  m_piF    = PI*m_F;
}
///
Real 
SineIF::
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
  int zdir = SpaceDim-1;
  bool mixedDeriv = false;
  int  ideriv = -1;
  bool foundANonZero = false;
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      if(a_partialDerivative[idir] > 0)
        {
          if(!foundANonZero)
            {
              foundANonZero = true;
              ideriv = idir;
            }
          else
            {
              mixedDeriv  = true;
            }
        }
    }
      
  if(!foundANonZero)
    {
      if(order != 0) 
        {
          MayDay::Error("Logic Error in SineIF");
        }
      //no derivatives
      //just evaluate the function
      retval = value(a_point);
      return retval; //need to return here because otherwise the inside/outside thing screws up
    }
  else if(mixedDeriv)
    {
      //all mixed derivs are zero since this is a sum of f(x) + g(y) + h(z)
      retval = 0;
    }
  else if ((a_partialDerivative[zdir] == 1))
    {
      //df/dy = -1
      retval = -1.0;
    }
  else if (a_partialDerivative[zdir] > 1)
    {
      //all higher derivs in y = 0
      retval = 0.0;
    }
  else
    {
      if(ideriv < 0) 
        {
          MayDay::Error("Logic Error 2 in SineIF");
        }
      //derivs in the sine directions are more complicated
      Real magni =  m_A[ideriv]*pow(m_piF[ideriv], a_partialDerivative[ideriv]);
      //this is my wacky way of changing signs every couple  of derivs
      int imodfour = a_partialDerivative[ideriv]%4;
      int imodtwo  = a_partialDerivative[ideriv]%2;
      Real rsign;
      if((imodfour == 0) || (imodfour == 1))
        {
          //deriv of sin is cos 
          rsign = 1;
        }
      else
        {
          CH_assert((imodfour == 2) || (imodfour == 3));
          //deriv of cos is  -sin
          rsign  = -1;
        }
      Real argum = m_piF[ideriv]*(a_point[ideriv] - m_point[ideriv]);
      if(imodtwo == 0)
        {
          retval = rsign*magni*sin(argum);
        }
      else
        {
          CH_assert(imodtwo == 1);
          retval = rsign*magni*cos(argum);
        }
    }

  if(!m_inside)
    retval = -retval;

  //  pout() << "retval2  = " << retval << endl;
  return retval;

}
///
Real 
SineIF::
value(const RealVect& a_point) const
{
  IndexTM<Real, SpaceDim> pt;

  for (int idir = 0; idir < SpaceDim; ++idir)
    {
      pt[idir] = a_point[idir];
    }

  //  pout() << "value 1  = " << value(pt) << endl;
  return value(pt);
}
///
Real 
SineIF::
value(const IndexTM<Real,SpaceDim>& a_point) const
{
  int zdir = SpaceDim-1;
  Real retval = -(a_point[zdir] - m_point[zdir]);
  //loop over non-z directions and add in the sine functions
  for(int idir = 0; idir < zdir; idir++)
    {
      Real magni = m_A[idir];
      Real argum = m_piF[idir]*(a_point[idir] - m_point[idir]);
      retval += magni*sin(argum);
    }

  if(!m_inside)
    retval = -retval;

  return retval;
}
////
BaseIF* 
SineIF::
newImplicitFunction() const
{
  SineIF* sinePtr = new SineIF(m_A,m_point,m_F,m_inside);
  return static_cast<BaseIF*>(sinePtr);
}

#include "NamespaceFooter.H"
