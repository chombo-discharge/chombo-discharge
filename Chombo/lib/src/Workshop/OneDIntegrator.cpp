#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


#include "OneDIntegrator.H"
#include "NamespaceHeader.H"

//integrate the function over the region dividing the 
//by dividing up the region of integration into 
//bits of fixed size = size/num_bits
Real 
OneDIntegrator::
integral(const RefCountedPtr<BaseOneDFunc>&  a_func,
         const Real                       &  a_xstart,
         const Real                       &  a_xend,
         const int                        &  a_num_bits
         ) const
{
  CH_assert(a_xstart < a_xend);
  CH_assert(a_num_bits > 0);

  Real retval = 0;
  Real dx = (a_xend - a_xstart)/a_num_bits;
  for(int ibit = 0; ibit < a_num_bits; ibit++)
    {
      Real xloc   = a_xstart + (ibit+0.5)*dx;
      Real valloc = a_func->value(xloc);
      retval += dx*valloc;
    }

  return retval;
}


#include "NamespaceFooter.H"

