#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CH_COUNTERS_H_
#define _CH_COUNTERS_H_

#include <sys/time.h>
#include <iostream>

#ifdef CH_PAPI
#include "papi.h"
#endif

#include "FortranNameMacro.H"

// These objects are outside all namespaces.  ch_timer is a COMMON BLOCK with the Fortran code
// and is defined in CH_CounterF.f.  In the longer run I will need to make ch_flops 
// atomic, or threadprivate, and do the same for the m_flop data member in CH_Timer
#include "BaseNamespaceHeader.H"

extern "C"{
  extern struct{
      long long int ch_flops;
  } FORTRAN_BASENAME(CH_TIMER,ch_timer);
}


inline double TimerGetTimeStampWC()
{
  struct timeval tv;   //  Values from call to gettimeofday
  struct timezone tz;
  gettimeofday(&tv, &tz);
  return((double)tv.tv_sec + 0.000001 * (double)tv.tv_usec);
}


#define NCOUNTERS 3
extern long long int ch_counters[NCOUNTERS];
extern int ch_eventset;

inline long long int& ch_flops(){ return  FORTRAN_BASENAME(CH_TIMER,ch_timer).ch_flops;}


void CountersInit();

inline void readCounters()
{
#ifdef CH_PAPI
  PAPI_read(ch_eventset, ch_counters);
#endif
}

void streamDump(std::ostream& os);

#include "BaseNamespaceFooter.H"
#endif
