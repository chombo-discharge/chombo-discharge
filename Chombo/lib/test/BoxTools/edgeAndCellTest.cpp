#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cmath>
#include "Misc.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "BoxIterator.H"
#include "LayoutIterator.H"
#include "LoadBalance.H"
#include "CellToEdge.H"
#include "EdgeToCell.H"
#include "LoadBalance.H"
#include "DebugDump.H"
#include "UsingNamespace.H"

/// Global variables for test output:

static const char *pgmname = "edgeandcelltest" ;
static const char *indent = "   " ,*indent2 = "      " ;
static bool verbose = true ;

int g_nx = 32;

int edgeToCellTest()
{
  int nx = g_nx;
  int retflag = 0;
  Real rightval = 4.0;
  Real wrongval = -6.0;
  Box cellDom(IntVect::Zero, (nx-1)*IntVect::Unit);
  FArrayBox cellStuff(cellDom, 1);
  FluxBox   faceStuff(cellDom, 1);
  
  cellStuff.setVal(wrongval);
  faceStuff.setVal(rightval);

  for(int idir = 0; idir < SpaceDim; idir++)
  {
    EdgeToCell(faceStuff, 0, cellStuff, 0, cellDom, idir);
    Real maxval = cellStuff.max(0);
    if(std::abs(maxval - rightval) > 1.0e-6)
    {
      retflag = -(idir+3);
      return retflag;
    }
  }
                                                     
  return retflag;
}
int cellToEdgeTest()
{
  int nx = g_nx;
  int retflag = 0;
  Real rightval = 4.0;
  Real wrongval = -6.0;
  Box cellDom(IntVect::Zero, (nx-1)*IntVect::Unit);
  FArrayBox cellStuff(grow(cellDom, 1), 1);
  FluxBox   faceStuff(cellDom, 1);
  
  cellStuff.setVal(rightval);
  faceStuff.setVal(wrongval);

  CellToEdge(cellStuff, faceStuff);
  for(int idir = 0; idir < SpaceDim; idir++)
  {
    Real maxval = faceStuff[idir].max(0);
    if(std::abs(maxval - rightval) > 1.0e-6)
    {
      retflag = -(idir+3);
      return retflag;
    }
  }
                                                     
  return retflag;
}
/*****************/
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
int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  bool passed = true;
  int icode = 0;

  parseTestOptions( argc ,argv ) ;
  if ( verbose )
    pout() << indent2 << "Beginning " << pgmname << " ..." << endl ;

  icode = cellToEdgeTest();

  if (icode != 0)
    {
      pout() << "celltoedgetest: failed with error code " << icode       << endl;
      passed = false;
    }
  icode = edgeToCellTest();

  if (icode != 0)
    {
      pout() << "edgetocelltest: failed with error code " << icode       << endl;
      passed = false;
    }

  if (passed)
    pout() << indent << pgmname
           << ": passed all tests"
           << endl;
#ifdef CH_MPI
  MPI_Finalize();
#endif

  return icode;
}
