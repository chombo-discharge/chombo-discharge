#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "SimpleIF.H"

#include "NamespaceHeader.H"

SimpleIF::SimpleIF()
{
}

SimpleIF::SimpleIF(const SimpleIF& a_inputIF)
{
}

SimpleIF::~SimpleIF()
{
}

Real SimpleIF::value(const RealVect& a_point) const
{
  // The overall function value
  Real retval;

  // Set z in case this is a 2D problem
  Real x;
  Real y;
  Real z = 0.5;

  // Set coordinates depending on the number of dimensions
  D_TERM(
  x = a_point[0];,
  y = a_point[1];,
  z = a_point[2];)

  // Never let x be less than -2.25 due to the sqrt(2.25+x) below
  if (x < -2.25)
  {
    x = -2.25;
  }

  // Compute teh overall function
  Real R;
  R = (-1.50 + sqrt(2.25+x))*(-1.50 + sqrt(2.25+x)) + (y/2.00)*(y/2.00) + z*z;

  retval = R*R - R*sqrt(R) - 0.80*z*R + 0.70*z*z*z;

  // Return the function value
  return retval;
}

BaseIF* SimpleIF::newImplicitFunction() const
{
  // Create and return a copy of this object
  SimpleIF* simplePtr = new SimpleIF();
  return static_cast<BaseIF*>(simplePtr);
}

#include "NamespaceFooter.H"
