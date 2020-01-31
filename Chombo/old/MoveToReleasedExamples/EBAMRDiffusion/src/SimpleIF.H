#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _SIMPLEIF_H_
#define _SIMPLEIF_H_

#include "MayDay.H"
#include "RealVect.H"

#include "BaseIF.H"

#include "NamespaceHeader.H"

///
/**
    This implicit function specifies a simple geometry for a cell.
 */
class SimpleIF: public BaseIF
{
public:
  /// Null constructor
  SimpleIF();

  /// Copy constructor
  SimpleIF(const SimpleIF& a_inputIF);

  /// Destructor
  virtual ~SimpleIF();

  /// Return the value of the function at a_point.
  virtual Real value(const RealVect& a_point) const;

  /// A factory to make more of this object
  virtual BaseIF* newImplicitFunction() const;

protected:

private:
  void operator=(const SimpleIF& a_inputIF)
  {
    MayDay::Abort("SimpleIF doesn't allow assignment");
  }
};

#include "NamespaceFooter.H"
#endif
