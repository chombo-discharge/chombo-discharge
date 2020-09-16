#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PYFUNCTIONS_H__
#define _PYFUNCTIONS_H__

#include "Python.h"
#include "RefCountedPtr.H"

// This handles the initialization of the Functions module.
PyMODINIT_FUNC initFunctions();

#endif
