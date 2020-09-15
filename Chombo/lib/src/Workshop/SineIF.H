#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _SINEIF_H_
#define _SINEIF_H_

#include "MayDay.H"
#include "RealVect.H"

#include "BaseIF.H"

#include "NamespaceHeader.H"

///
/**
    This implicit function specifies a 
    sine function
    given constants 
    x0, y0 
    A  (magnitude) 
    F  (frequency (first two comps))
    2d:
    f(x,y) = y - y0 + A*sin(pi*F(x-x0))
    3d:
    f(x,y,z) = z - z0 + A0*sin(pi*F0(x-x0))+ A1*sin(pi*F1(y-y0))
    
 */
class SineIF: public BaseIF
{
public:
  ///
  /**
     last comp of frequency ignored
     2d:
     point = (x0, y0)
     freq = (f0, f1)
    f(x,y) = y - y0 + A0*sin(pi*f0(x-x0))
    3d:
    point = (x0, y0, z0)
    freq = (f0, f1, f2)
    f(x,y,z) = z - z0 + A0*sin(pi*(x-x0)) + A1*sin(pi*(y*y0))
    
   */
  SineIF(const RealVect & a_A,
         const RealVect & a_point,
         const RealVect & a_F,
         const bool     & a_inside);

  /// Destructor
  virtual ~SineIF()
  {
  }

  virtual Real value(const IndexTM<int,SpaceDim> & a_partialDerivative,
                     const IndexTM<Real,SpaceDim>& a_point) const;

  virtual Real value(const RealVect& a_point) const;

  virtual Real value(const IndexTM<Real,SpaceDim>& a_point) const;

  virtual BaseIF* newImplicitFunction() const;

protected:

  RealVect m_A;      
  RealVect m_point;
  RealVect m_F;      
  bool     m_inside;
  RealVect m_factor;
  RealVect m_piF;
private:
  SineIF();

};

#include "NamespaceFooter.H"
#endif
