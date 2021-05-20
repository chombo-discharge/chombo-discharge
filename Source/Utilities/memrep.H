/*!
  @file DataOps.H
  @brief Agglomeration of useful data operations
  @author Robert Marskar
  @date Nov 2017
*/

#ifndef _MEMREP_
#define _MEMREP_

#include "RealVect.H"

#include "CD_NamespaceHeader.H"

class memrep {
public:
  
  static void getMaxMin_memory();

  static void getMaxMin_memory(Real& a_max_peak, Real& a_min_peak, Real& a_min_unfreed, Real& a_max_unfreed);

  static void get_memory(Vector<Real>& a_peak, Vector<Real>& a_unfreed);
};

#include "CD_NamespaceFooter.H"

#endif
