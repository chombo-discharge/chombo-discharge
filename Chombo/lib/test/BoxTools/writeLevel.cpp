#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "REAL.H"
#include "parstream.H"

#include "LevelData.H"
#include "LoadBalance.H"
#include "FArrayBox.H"
#include "IntVectSet.H"
#include "AMRIO.H"
#include "BoxIterator.H"
#include "BRMeshRefine.H" // contains domainSplit function
#include "CH_Attach.H"
#include "UsingNamespace.H"

using std::endl;

static bool verbose = true;


/***************/
/***************/
void dumpmemoryatexit();
/***************/
/***************/


void
parseTestOptions( int argc ,char* argv[] ) ;

int writeLevelTest();

int
main(int argc,char **argv)
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
  //  registerDebugger();
#endif
  parseTestOptions( argc ,argv ) ;

  int retval = 0;

  retval = writeLevelTest();

  if (retval == 0)
    {
      pout() << "writeLevel passed all tests!" << endl;
    }
  else
    {
      pout() << "writeLevel failed at least one test, return value = "
             << retval << endl;
    }

#ifdef CH_MPI
  MPI_Finalize();
#endif

  return retval;
}


int writeLevelTest()
{
  int domLength =  128;
  int maxBoxSize =   8;


  Box domainBox(IntVect::Zero, (domLength-1)*IntVect::Unit);
  
  // start with non-periodic case
  ProblemDomain domain(domainBox);

  pout() << "calling domainSplit" << std::endl;
  Vector<Box> gridBoxes;
  domainSplit(domainBox, gridBoxes, maxBoxSize);
  int numBoxes = gridBoxes.size();
  Vector<int> procAssign(numBoxes,0);

  pout() << "calling loadbalance" << std::endl;
  // for mpi case, distribute boxes among processors
  LoadBalance(procAssign,gridBoxes);
  pout() << "calling dbl define" << std::endl;
  DisjointBoxLayout grids(gridBoxes, procAssign, domain);

  const int nComp = 1;
  IntVect ghostVect = IntVect::Unit;
  LevelData<FArrayBox> data(grids, nComp, ghostVect);

  writeLevel(&data);

  return 0;
}


///
// Parse the standard test options (-v -q) out of the command line.
// Stop parsing when a non-option argument is found.
///
void
parseTestOptions( int argc ,char* argv[] )
{
  for ( int i = 1 ; i < argc ; ++i )
  {
    if ( argv[i][0] == '-' ) //if it is an option
    {
      // compare 3 chars to differentiate -x from -xx
      if ( strncmp( argv[i] ,"-v" ,3 ) == 0 )
      {
        verbose = true ;
        // argv[i] = "" ;
      }
      else if ( strncmp( argv[i] ,"-q" ,3 ) == 0 )
      {
        verbose = false ;
        // argv[i] = "" ;
      }
      else
      {
        break ;
      }
    }
  }
  return ;
}
