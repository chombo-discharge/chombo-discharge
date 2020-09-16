#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef _CH_THREAD_H_
#define _CH_THREAD_H_

#include "BaseNamespaceHeader.H"

//returns true if we are on thread 0 (or if not threaded)
extern bool onThread0();

//returns the value of OMP_NUM_THREADS (or 1 if not threaded)
extern int  getMaxThreads();

#include "BaseNamespaceFooter.H"
#endif
