#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CHOMBOPLOT_H_
#define _CHOMBOPLOT_H_

#include "datatype.h"
/* This file contains the various things related to the main body of the
 * program.  It is pretty sparse, and really shouldn't be too cluttered
 * up.
 *
 */

/* protos */
extern void init_display(int argc, char **argv, datatype *me);

#ifndef TRUE
#define TRUE  1
#define FALSE 0
#endif

#endif
