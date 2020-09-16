#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PARABOLAIF_H_
#define _PARABOLAIF_H_

#include "MayDay.H"
#include "RealVect.H"

#include "BaseIF.H"

#include "NamespaceHeader.H"

///
/**
    This implicit function specifies a parabola.
 */
class ParabolaIF: public BaseIF
{
public:
  ///
  /**
      y = a*(x-x0)^2 + y0, (x0,y0) = a_point
   */
  ParabolaIF(const Real& a_a,
             const RealVect& a_point,
             const bool&     a_inside);

  /// Destructor
  virtual ~ParabolaIF()
  {
  }

  virtual Real value(const IndexTM<int,SpaceDim> & a_partialDerivative,
                     const IndexTM<Real,SpaceDim>& a_point) const;

  virtual Real value(const RealVect& a_point) const;

  virtual Real value(const IndexTM<Real,SpaceDim>& a_point) const;

  virtual BaseIF* newImplicitFunction() const;

protected:
   Real          m_a;
   RealVect      m_point;
   bool          m_inside;
private:
  ParabolaIF();

};

#include "NamespaceFooter.H"
#endif
