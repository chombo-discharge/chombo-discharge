#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//
// These are the inevitable functions that people
// can't avoid using when writing a SPMD program.
// It is a minimal set from which more elaborate
// functionality can be generated.  As always, a
// user is free to utilize the entire MPI programming
// on their own platform.  The functions are
// assured to work on all machines supported.
//

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include "datatype.h"

#include "libsx.h"             /* defines libsx stuff */
#include "chomboPlot.h"        /* where program specific stuff is defined */
#include "callbacks.h"         /* prototypes for callback functions       */
#ifdef USE_ARRAYVIEW
#include "ArrayView.H"
#endif
#include "SPMD.H"

int main(int argc, char **argv)
{
  if (numProc() > 1)
      MayDay::Error("chombo plot not designed for parallel runs");

  datatype data;

  init_display(argc, argv, &data);  /* setup the display */
  MainLoop();                         /* go right into the main loop */
  return(0);
}
