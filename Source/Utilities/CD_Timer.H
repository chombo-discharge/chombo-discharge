/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_Timer.H
  @brief  Implementation of CD_Timer.H
  @author Robert Marskar
*/

// Std includes
#include <chrono>

// Chombo includes
#include <REAL.H>

// Our includes
#include <CD_NamespaceHeader.H>

class Timer {
public:
  
  inline
  static Real wallClock(){
#ifdef CH_MPI
    return Timer::wallClock();
#else
    return 0.0;
#endif
  }

};

#include <CD_NamespaceFooter.H>
