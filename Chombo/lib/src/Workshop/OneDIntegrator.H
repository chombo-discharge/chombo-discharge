#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _ONEDINTEGRATOR_H_
#define _ONEDINTEGRATOR_H_

#include "RealVect.H"
#include "ProblemDomain.H"
#include "IndexTM.H"
#include "EB_TYPEDEFS.H"

#include "Notation.H"
#include "GeometryService.H"
#include "IndexedMoments.H"
#include "NamespaceHeader.H"

///base class for integration
class  BaseOneDFunc
{
public:
  BaseOneDFunc()
  {;  }

  virtual ~BaseOneDFunc()
  {; }

  virtual Real value(const Real& location) const = 0;
};

///
/**
   Gives an approximation to the numerical integral of a function over a region
 */
class OneDIntegrator
{

public:

  /// Default constructor
  OneDIntegrator()
  {
  }

  /// Default destructor
  virtual ~OneDIntegrator()
  {
  }

  ///
  /**
     integrate the function over the region dividing the 
     by dividing up the region of integration into 
     bits of fixed size = size/num_bits (midpoint rule)
  */
  Real integral(const RefCountedPtr<BaseOneDFunc>&  a_func,
                const Real                       &  a_xstart,
                const Real                       &  a_xend,
                const int                        &  a_num_bits
                ) const;

};

#include "NamespaceFooter.H"

#endif
