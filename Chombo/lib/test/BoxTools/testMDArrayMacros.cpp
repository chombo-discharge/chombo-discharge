#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>
#include <cstring>
using std::endl;
#include "parstream.H"
#include "FArrayBox.H"
#include "CHArray.H"
#include "BaseFabMacros.H"
#ifdef CH_MPI
#include <mpi.h>
#endif
#include "UsingNamespace.H"

/// Prototypes:

void
parseTestOptions( int argc ,char* argv[] ) ;

int
testMDArrayMacros();

/// Global variables for handling output:
static const char *pgmname = "testMDArrayMacros" ;
static const char *indent = "   ", *indent2 = "      " ;
static bool verbose = true ;

/// Code:

int
main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  parseTestOptions( argc ,argv ) ;

  if (verbose)
    pout() << indent2 << "Beginning " << pgmname << " ..." << endl ;

  ///
  // Run the tests
  ///
  int ret = testMDArrayMacros() ;

  if (ret == 0)
    {
      pout() << indent << pgmname << " passed all tests" << endl ;
    }
  else
    {
      pout() << indent << pgmname << " failed " << ret << " test(s)" << endl ;
    }

#ifdef CH_MPI
  MPI_Finalize();
#endif
  return ret;
}


/*******************************************************************************
 *
 * Routine testMDArrayMacros
 *
 ******************************************************************************/

int
testMDArrayMacros()
{
  int errors = 0;

/*--------------------------------------------------------------------*
 * Parameters
 *--------------------------------------------------------------------*/

  const int n = 4;

/*--------------------------------------------------------------------*
 * Test MDArray access
 *--------------------------------------------------------------------*/

#if CXXSTD>=14

//--BaseFab

  {  // +ve
    Box box(IntVect::Unit, n*IntVect::Unit);
    FArrayBox fab(box, 2);

    for (int iC = 0; iC != 2; ++iC)
      {
        MD_BOXLOOP(box, i)
          {
            if (&(fab[MD_IX(i, iC)]) !=
                &(fab(IntVect(D_DECL6(i0, i1, i2, i3, i4, i5)), iC)))
              {
                ++errors;
              }
          }
      }
  }

  {  // -ve
    Box box(-n*IntVect::Unit, -IntVect::Unit);
    FArrayBox fab(box, 2);

    for (int iC = 0; iC != 2; ++iC)
      {
        MD_BOXLOOP(box, i)
          {
            if (&(fab[MD_IX(i, iC)]) !=
                &(fab(IntVect(D_DECL6(i0, i1, i2, i3, i4, i5)), iC)))
              {
                ++errors;
              }
          }
      }
  }

  {  // spanning 0
    Box box(-IntVect::Unit, (n-2)*IntVect::Unit);
    FArrayBox fabA(box, 2);

    for (int iC = 0; iC != 2; ++iC)
      {
        MD_BOXLOOP(box, i)
          {
            if (&(fabA[MD_IX(i, iC)]) !=
                &(fabA(IntVect(D_DECL6(i0, i1, i2, i3, i4, i5)), iC)))
              {
                ++errors;
              }
          }
      }

    // Test move
    FArrayBox fabB(std::move(fabA));
    for (int iC = 0; iC != 2; ++iC)
      {
        MD_BOXLOOP(box, i)
          {
            if (&(fabB[MD_IX(i, iC)]) !=
                &(fabB(IntVect(D_DECL6(i0, i1, i2, i3, i4, i5)), iC)))
              {
                ++errors;
              }
          }
      }

    // Test move assignment
    FArrayBox fabC(Box(-IntVect::Unit, IntVect::Unit), 1);
    fabC = std::move(fabB);
    for (int iC = 0; iC != 2; ++iC)
      {
        MD_BOXLOOP(box, i)
          {
            if (&(fabC[MD_IX(i, iC)]) !=
                &(fabC(IntVect(D_DECL6(i0, i1, i2, i3, i4, i5)), iC)))
              {
                ++errors;
              }
          }
      }

    // Test shift
    fabC.shift(0, 2);
    for (int iC = 0; iC != 2; ++iC)
      {
        MD_BOXLOOP(fabC.box(), i)
          {
            if (&(fabC[MD_IX(i, iC)]) !=
                &(fabC(IntVect(D_DECL6(i0, i1, i2, i3, i4, i5)), iC)))
              {
                ++errors;
              }
          }
      }

    // Test shift
    fabC.shift(0, -4);
    for (int iC = 0; iC != 2; ++iC)
      {
        MD_BOXLOOP(fabC.box(), i)
          {
            if (&(fabC[MD_IX(i, iC)]) !=
                &(fabC(IntVect(D_DECL6(i0, i1, i2, i3, i4, i5)), iC)))
              {
                ++errors;
              }
          }
      }
  }

#endif

  return errors;
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
      if (argv[i][0] == '-') //if it is an option
        {
          // compare 3 chars to differentiate -x from -xx
          if (strncmp( argv[i] ,"-v" ,3 ) == 0)
            {
              verbose = true ;
              // argv[i] = "" ;
            }
          else if (strncmp( argv[i] ,"-q" ,3 ) == 0)
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
