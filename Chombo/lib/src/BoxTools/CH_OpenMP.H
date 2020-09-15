#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CH_OPENMP_H_
#define _CH_OPENMP_H_

#include "DataIterator.H"
//#include "FABSTACKTEMP.H"
//#include "FABSTACKSPACEDIM.H"
//#include "FABSTACKSDSQUARED.H"
#include "NamespaceHeader.H"



#if _OPENMP==1
#include <omp.h>
#define CH_SCHEDULE dynamic

#endif
                
#include "NamespaceFooter.H"
#endif
