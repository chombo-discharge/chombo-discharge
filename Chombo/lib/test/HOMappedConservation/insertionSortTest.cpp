#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cstring>

#include "SPMD.H"
#include "parstream.H"
#include "InsertionSort.H"

#include "UsingNamespace.H"

//--Prototypes

int
insertionSortTest();

void
parseTestOptions(int argc ,char* argv[]);

//--Global variables for handling output:

static const char *pgmname = "insertionSortTest";
static const char *indent = "   ", *indent2 = "      ";
static bool verbose = true;

/*--------------------------------------------------------------------*
 *  Entry
 *--------------------------------------------------------------------*/

int
main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  parseTestOptions( argc ,argv ) ;
  if ( verbose )
    pout() << indent2 << "Beginning " << pgmname << " ..." << std::endl ;

  int eekflag = insertionSortTest();
  if (eekflag == 0)
    {
      pout() << indent << pgmname
             << ": test passed." << std::endl;
    }
  else
    {
      pout() << indent << pgmname
             << ": test FAILED with error code "
             << eekflag << std::endl;
    }
#ifdef CH_MPI
  MPI_Finalize();
#endif
  return eekflag;
}

/*--------------------------------------------------------------------*
 *  Insertion sort testing
 *--------------------------------------------------------------------*/

int insertionSortTest()
{
  int status = 0;
  srand(0);
  const int n = 1000;
  int *a = new int[n];
  int *b = new int[n];
  int *ix = new int[n];
  for (int i = 0; i != n; ++i)
    {
      a[i] = (int)((((double)rand())/RAND_MAX)*1000);
      b[i] = a[i];
      ix[i] = i;
    }
  Sort::insertion(n, ix, Sort::CmpGreaterIndex<int>(a));
  Sort::Move2Array<int, int> mover(a, b);
  Sort::rearrangeToIndex(n, ix, mover);
  for (int i = 1; i != n; ++i)
    {
      if (a[i-1] < a[i] || b[i-1] < b[i] || ix[i-1] > ix[i])
        {
          status |= 1;
        }
    }
  if (status & 1)
    {
      pout() << "Insertion index sort and/or rearrange for 2 arrays failed\n";
    }

  // Test the mask
  int c[10] =
  {
    0, 8, 6, 3, 4, 5, 7, 1, 2, 9
  };
  int m[10] =
  {
    0, 1, 1, 0, 0, 0, 1, 1, 1, 0
  };
  Sort::insertion(10, c, std::less<int>(), m);
  for (int i = 1; i != 10; ++i)
    {
      if (c[i-1] > c[i])
        {
          status |= 2;
        }
    }
  if (status & 2)
    {
      pout() << "Insertion sort with mask failed\n";
    }
  return status;
}

/*--------------------------------------------------------------------*
 *  Parse the standard test options (-v -q) out of the command line.
 *--------------------------------------------------------------------*/

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
        }
    }
  return ;
}
