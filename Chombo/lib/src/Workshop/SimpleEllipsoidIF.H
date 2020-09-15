#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _SIMPLEELLIPSOIDIF_H_
#define _SIMPLEELLIPSOIDIF_H_

#include "MayDay.H"
#include "RealVect.H"

#include "BaseIF.H"

#include "NamespaceHeader.H"

///
/**
   This implicit function specifies a SimpleEllipsoid.
   pout() << "SimpleEllipsoid geometry" << endl;
   pout() << "(x-xo/a)^2 + (y-y0/b)^2 + (z-z/c)^ - r^2 = f(x,y,z)" << endl;

*/
class SimpleEllipsoidIF: public BaseIF
{
public:
  ///
  /**
     "(x-xo/A)^2 + (y-y0/B)^2 + (z-z0/C)^ - r^2 = f(x,y,z)" << endl;
   */
  SimpleEllipsoidIF(const RealVect& a_ABE,
              const RealVect& a_X0,
              const Real    & a_R,
              const bool    & a_inside)
  {
    m_ABE=      a_ABE;   
    m_X0=       a_X0;    
    m_R=        a_R;     
    m_inside=   a_inside;

    for(int idir = 0; idir < SpaceDim; idir++)
      {
        CH_assert(a_ABE[idir] > 1.e-10);
        m_A2[idir] = a_ABE[idir]*a_ABE[idir];
      }
  }

  /// Destructor
  virtual ~SimpleEllipsoidIF()
  {
  }

  virtual Real value(const IndexTM<int,SpaceDim> & a_partialDerivative,
                     const IndexTM<Real,SpaceDim>& a_point) const;

  virtual Real value(const RealVect& a_point) const;

  virtual Real value(const IndexTM<Real,SpaceDim>& a_point) const;

  virtual BaseIF* newImplicitFunction() const
  {
    SimpleEllipsoidIF* SimpleEllipsoidPtr = new SimpleEllipsoidIF(m_ABE,m_X0, m_R, m_inside);
    return static_cast<BaseIF*>(SimpleEllipsoidPtr);
  }

protected:
  RealVect m_ABE;
  RealVect m_X0;
  Real     m_R;
  bool     m_inside;
  RealVect m_A2;

private:
  SimpleEllipsoidIF();

};

#include "NamespaceFooter.H"
#endif
