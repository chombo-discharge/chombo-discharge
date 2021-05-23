/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MemoryReport.H
  @brief  Utility functions for fetching memory usage
  @author Robert Marskar
*/

#ifndef CD_MemoryReport_H
#define CD_MemoryReport_H

// Chombo includes
#include <RealVect.H>

// Our includes
#include <CD_NamespaceHeader.H>

namespace MemoryReport {
  
  void getMaxMinMemoryUsage();

  void getMaxMinMemoryUsage(Real& a_max_peak, Real& a_min_peak, Real& a_min_unfreed, Real& a_max_unfreed);

  void getMemoryUsage(Vector<Real>& a_peak, Vector<Real>& a_unfreed);
};

#include <CD_NamespaceFooter.H>

#endif
