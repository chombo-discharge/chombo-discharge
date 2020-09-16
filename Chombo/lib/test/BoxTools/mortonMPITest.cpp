#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// dtg 2018

#include <cmath>
#include "BoxLayout.H"
#include "Misc.H"
#include "BRMeshRefine.H"
#include "DebugDump.H"
#include "BoxLayout.H"
#include "UsingNamespace.H"

/**
   mortonTest returns:
   0: all tests passed
*/
int
mortonTest()
{
  int retval = 0;
#ifdef CH_MPI
  int ncellmax = 16;
  int domsize = 8192;
  Box domain(IntVect::Zero, (domsize-1)*IntVect::Unit);
  Vector<Box> vbox;
  domainSplit(domain, vbox, ncellmax, 1);

  Vector<Box> serbox = vbox;
  Vector<Box> mpibox = vbox;
  
  mortonOrdering(mpibox);
  serialMortonOrdering(serbox);
 
  if(mpibox.size() != vbox.size())
  {
    retval = -2;
    return retval;
  }
  if(serbox.size() != vbox.size())
  {
    retval = -3;
    return retval;
  }
  for(int ibox = 0; ibox < vbox.size(); ibox++)
  {
    if(serbox[ibox] != mpibox[ibox])
    {
      retval = -4;
      return retval;
    }
  }
  
#endif
  return retval;
}

/// Code:
int
main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  int icode = mortonTest();
  if (icode != 0)
    {
      pout() << " mortonOrdering failed with error code " << icode   << endl;
    }
  else
    {
      pout() << "mortonOrdering passed all tests"   << endl;
    }
#ifdef CH_MPI
  MPI_Finalize();
#endif
  return icode;
}


