#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "KeratocyteIF.H"

#include "NamespaceHeader.H"

KeratocyteIF::KeratocyteIF()
{
}

KeratocyteIF::KeratocyteIF(const KeratocyteIF& a_inputIF)
{
}

KeratocyteIF::~KeratocyteIF()
{
}

Real KeratocyteIF::value(const RealVect& a_point) const
{
  // The overall function value
  Real retval;

  // Set z in case this is a 2D problem (a slight perturbation to original
  // problem
  Real x;
  Real y;
  Real z = 0.2002;

  // Set coordinates depending on the number of dimensions
  D_TERM(
  x = a_point[0];,
  y = a_point[1];,
  z = a_point[2];)

  // Compute the four functions used to define the cell
  Real cell1 = (x/5.00)*(x/5.00) + (y/20.00)*(y/20.00) - 1.0001;
  Real cell2 = 0.06*x + z - 0.60;
  Real cell3 = ((x+3.50)/3.50)*((x+3.50)/3.50)
             + ( y      /5.00)*( y      /5.00)
             + ((z-1.00)/3.50)*((z-1.00)/3.50) - 1.0001;
  Real cell4 = 0.2001 - z;

  // Compute (((cell1 intersect cell2) union cell3) intersect cell4
  Real cell = Max(cell1,cell2);
  cell      = Min(cell ,cell3);
  cell      = Max(cell ,cell4);

  // Compute the two functions used to define the nucleus
  Real nucleus1 = (x+3.50)*(x+3.50) + y*y + (z-2.00)*(z-2.00) - 5.7601;
  Real nucleus2 = 0.3001 - z;

  // Take there intersection
  Real nucleus = Max(nucleus1,nucleus2);

  // Take cell minus nucleus = cell intersect (complement nucleus)
  retval = Max(cell,-nucleus);

  // Return the function value
  return retval;
}

BaseIF* KeratocyteIF::newImplicitFunction() const
{
  // Create and return a copy of this object
  KeratocyteIF* keratocytePtr = new KeratocyteIF();
  return static_cast<BaseIF*>(keratocytePtr);
}

#include "NamespaceFooter.H"
