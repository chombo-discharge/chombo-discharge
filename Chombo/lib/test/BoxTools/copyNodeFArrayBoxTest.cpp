#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// petermc, 15 Jun 2016

#include <cmath>
#include <limits> // contains numeric_limits<Real>
#include "Misc.H"
#include "NodeFArrayBox.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "BoxIterator.H"
#include "LayoutIterator.H"
#include "CH_Attach.H"
#include "UsingNamespace.H"

/// Global variables for test output:

static const char *pgmname = "copyNodeFArrayBoxTest" ;
static const char *indent = "   " ,*indent2 = "      " ;
static bool verbose = true ;
static Real precision = 2 * numeric_limits<Real>::epsilon();
static Real srcVal = 3.14;
static Real destVal = 3.89;

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


void setValLevelDataNodeFArrayBox(LevelData<NodeFArrayBox>& a_data,
                                  Real a_val)
{
  const DisjointBoxLayout& layout = a_data.disjointBoxLayout();
  for (DataIterator dit = a_data.dataIterator(); dit.ok(); ++dit)
    {
      const Box& cells = layout[dit];
      FArrayBox& fab = a_data[dit].getFab();
      // Don't fill ghost nodes of a_data.
      Box nodes = surroundingNodes(cells);
      fab.setVal(a_val, nodes, 0);
    }
}


bool copied(const LevelData<NodeFArrayBox>& a_dest,
            const Box& a_bx)
{
  int numFailed = 0;
  for (DataIterator dit = a_dest.dataIterator(); dit.ok(); ++dit)
    {
      const FArrayBox& destFab = a_dest[dit].getFab();
      FArrayBox intersectFab(a_bx, 1);
      intersectFab.copy(destFab);
      intersectFab -= srcVal;
      if (intersectFab.norm(0) > precision)
        {
          numFailed++;
        }
    }
#ifdef CH_MPI
  int recv;
  int resultMPI = MPI_Allreduce(&numFailed, &recv,
                                1, MPI_INT,
                                MPI_MAX, Chombo_MPI::comm);
  if (resultMPI != MPI_SUCCESS)
    { //bark!!!
      MayDay::Error("communication error on MPI_Allreduce");
    }
  numFailed = recv;
#endif
  return (numFailed == 0);
}


int copyNodeFArrayBoxTest()
{
  int icode = 0;

  int boxlen = 8;
  IntVect lo0 = IntVect::Zero;
  IntVect hi0 = (boxlen - 1) * IntVect::Unit;
  Box b0cells(lo0, hi0, IndexType::TheCellType());
  Box b1cells(b0cells);
  b1cells.shift(boxlen * BASISV(0));

  Box b0nodes = surroundingNodes(b0cells);
  Box b1nodes = surroundingNodes(b1cells);
  Box intersection = b0nodes & b1nodes;

  int proc0 = 0;
  int proc1 = numProc() - 1;

  /*
    Test: source and destination layouts are the same 2-box array,
    on the same processors.
  */
  {
    Vector<Box> twoBoxArray(2);
    twoBoxArray[0] = b0cells;
    twoBoxArray[1] = b1cells;
    Vector<int> twoBoxProcs(2);
    twoBoxProcs[0] = proc0;
    twoBoxProcs[1] = proc1;
    DisjointBoxLayout twoBoxLayout(twoBoxArray, twoBoxProcs);
    LevelData<NodeFArrayBox> src(twoBoxLayout, 1);
    LevelData<NodeFArrayBox> dest(twoBoxLayout, 1);
    setValLevelDataNodeFArrayBox(src, srcVal);
    setValLevelDataNodeFArrayBox(dest, destVal);
    src.copyTo(dest);
    if (copied(dest, intersection))
      {
        pout() << indent << pgmname
               << ": two-box same-procs test passed" << endl;
      }
    else
      {
        pout() << indent << pgmname
               << ": two-box same-procs test FAILED" << endl;
        icode = -1;
      }
  }

  /* Test: from box 0 to box 1.  THIS TEST FAILS.  So don't even try it.
   */
  /*
  {
    Vector<Box> srcBoxArray(1, b0cells);
    Vector<Box> destBoxArray(1, b1cells);
    Vector<int> procs(1, 0);
    DisjointBoxLayout srcLayout(srcBoxArray, procs);
    DisjointBoxLayout destLayout(destBoxArray, procs);
    LevelData<NodeFArrayBox> src(srcLayout, 1);
    LevelData<NodeFArrayBox> dest(destLayout, 1);
    setValLevelDataNodeFArrayBox(src, srcVal);
    setValLevelDataNodeFArrayBox(dest, destVal);
    src.copyTo(dest);
    if (copied(dest, intersection))
      {
        pout() << indent << pgmname
               << ": single adjacent test passed" << endl;
      }
    else
      {
        pout() << indent << pgmname
               << ": single adjacent test FAILED" << endl;
        icode = -2;
      }
  }
  */

  /* Test: from box 0 to box 1 via expanded box
   */
  {
    Vector<Box> srcBoxArray(1, b0cells);
    Vector<Box> destBoxArray(1, b1cells);
    Vector<int> procs(1, 0);
    DisjointBoxLayout srcLayout(srcBoxArray, procs);
    DisjointBoxLayout destLayout(destBoxArray, procs);
    LevelData<NodeFArrayBox> src(srcLayout, 1);
    LevelData<NodeFArrayBox> dest(destLayout, 1);
    setValLevelDataNodeFArrayBox(src, srcVal);
    setValLevelDataNodeFArrayBox(dest, destVal);
    // src.copyTo(dest);
    BoxLayout destExpandedLayout;
    destExpandedLayout.deepCopy(destLayout);
    destExpandedLayout.grow(1);
    destExpandedLayout.closeNoSort();
    BoxLayoutData<NodeFArrayBox> destExpanded(destExpandedLayout, 1);
    src.copyTo(destExpanded);
    for (DataIterator dit = destLayout.dataIterator(); dit.ok(); ++dit)
      {
        dest[dit].copy(destExpanded[dit]);
      }
    if (copied(dest, intersection))
      {
        pout() << indent << pgmname
               << ": single adjacent test to expanded dest passed" << endl;
      }
    else
      {
        pout() << indent << pgmname
               << ": single adjacent test to expanded dest FAILED" << endl;
        icode = -3;
      }
  }

  /* Test: from box 0 to box 1 via Copier
   */
  {
    Vector<Box> srcBoxArray(1, b0cells);
    Vector<Box> destBoxArray(1, b1cells);
    Vector<int> procs(1, 0);
    DisjointBoxLayout srcLayout(srcBoxArray, procs);
    DisjointBoxLayout destLayout(destBoxArray, procs);
    LevelData<NodeFArrayBox> src(srcLayout, 1, IntVect::Unit);
    LevelData<NodeFArrayBox> dest(destLayout, 1);
    setValLevelDataNodeFArrayBox(src, srcVal);
    setValLevelDataNodeFArrayBox(dest, destVal);
    Copier copier;
    copier.ghostDefine(srcLayout, destLayout,
                       srcLayout.physDomain(), IntVect::Unit);
    src.copyTo(dest, copier);
    if (copied(dest, intersection))
      {
        pout() << indent << pgmname
               << ": single adjacent test with Copier passed" << endl;
      }
    else
      {
        pout() << indent << pgmname
               << ": single adjacent test with Copier FAILED" << endl;
        icode = -4;
      }
  }

#ifdef CH_MPI
  /*
    Test: source and destination layouts are the same 2-box array,
    on different processors.
  */
  if (numProc() > 1)
  {
    Vector<Box> twoBoxArray(2);
    twoBoxArray[0] = b0cells;
    twoBoxArray[1] = b1cells;
    Vector<int> twoBoxSrcProcs(2);
    twoBoxSrcProcs[0] = proc0;
    twoBoxSrcProcs[1] = proc1;
    Vector<int> twoBoxDestProcs(2);
    twoBoxDestProcs[0] = proc1;
    twoBoxDestProcs[1] = proc0;
    DisjointBoxLayout twoBoxSrcLayout(twoBoxArray, twoBoxSrcProcs);
    DisjointBoxLayout twoBoxDestLayout(twoBoxArray, twoBoxDestProcs);
    LevelData<NodeFArrayBox> src(twoBoxSrcLayout, 1);
    LevelData<NodeFArrayBox> dest(twoBoxDestLayout, 1);
    setValLevelDataNodeFArrayBox(src, srcVal);
    setValLevelDataNodeFArrayBox(dest, destVal);
    src.copyTo(dest);
    if (copied(dest, intersection))
      {
        pout() << indent << pgmname
               << ": two-box opposite-procs test passed" << endl;
      }
    else
      {
        pout() << indent << pgmname
               << ": two-box opposite-procs test FAILED" << endl;
        icode = -5;
      }
  }
#endif

  return icode;
}


/// Code:
int
main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  bool passed = true;
  int icode = 0;
  //registerDebugger();
  parseTestOptions( argc ,argv ) ;
  if ( verbose )
    pout() << indent2 << "Beginning " << pgmname << " ..." << endl ;

  icode = copyNodeFArrayBoxTest();
  if (icode != 0)
    {
      pout() << indent << pgmname
           << " failed with error code " << icode
           << endl;
      passed = false;
    }

  if (passed)
    {
      pout() << indent << pgmname
           << " passed all tests"
           << endl;
    }
#ifdef CH_MPI
  MPI_Finalize();
#endif
  return icode;
}
