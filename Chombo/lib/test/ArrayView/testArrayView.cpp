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
// testArrayView.cpp
//
// test the ArrayView functions
//
#include <cstring>

#include "IntVect.H"
#include "Box.H"
#include "BoxIterator.H"
#include "BaseFab.H"
#include "REAL.H"
#include "FArrayBox.H"
#include "LayoutData.H"
#include "BoxLayout.H"
#include "DisjointBoxLayout.H"
#include "LayoutIterator.H"
#include "LevelData.H"
#include "LoadBalance.H"
#include "ArrayView.H"
#include "UsingNamespace.H"

/// Prototypes:

void
parseTestOptions(int argc ,char* argv[]) ;

void
parseInput(int argc, char* argv[]);

void
testArrayView();

/// Global variables for handling output:

static const char *pgmname = "testArrayView" ;
static const char *indent = "   ", *indent2 = "      " ;
static bool verbose = true ;

/// Code:

int
main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
  if (procID() == 0)
  {
    pout() << "Sorry, Arrayview not implemented for parallel neither pass nor fail"<<endl;
  }
  MPI_Finalize();
  return 0;
#endif
  parseTestOptions( argc ,argv ) ;

  if ( verbose )
    pout() << indent2 << "Beginning " << pgmname << " ..." << endl ;

  parseInput(argc, argv);

  testArrayView();

  if ( verbose )
    pout() << indent2 << "... " << pgmname << " finished" << endl ;
}

void
parseInput(int argc, char* argv[])
{
  // set port offset for socket creation
  arrayview_user_port_offset = 0;
  // assume there is a number following the `-p' so dont
  // look at the last argument
  for (int count = 1; count < argc-1; ++count)
  {
    const char* arg = argv[count];
    if (strncmp(arg, "-p" ,3) == 0)  //compare 3 to differentiate -p from -px
    {
      arrayview_user_port_offset = atoi(argv[++count]);
    }
  }
}

void
testArrayView()
{
// single grid tests
  Box b(IntVect::Zero, 3*IntVect::Unit);
  BoxIterator bit(b);

// FArrayBox test
  FArrayBox fab(b,SpaceDim);
  for (bit.begin(); bit.ok(); ++bit)
  {
    IntVect iv = bit();
    for (int d = 0; d < SpaceDim; ++d)
    {
      fab(iv,d) = iv[d];
    }
  }
  pout() << indent << pgmname << ": testing FAB" << endl ;
  ArrayView(&fab);

// BaseFab<Real> test
  BaseFab<Real> basefab(b,SpaceDim);
  for (bit.begin(); bit.ok(); ++bit)
  {
    IntVect iv = bit();
    for (int d = 0; d < SpaceDim; ++d)
    {
      basefab(iv,d) = iv[d];
    }
  }
  pout() << indent << pgmname << ": testing BaseFAB<Real>" << endl ;
  ArrayView(&basefab);

// multiple grids tests
    Vector<Box> grids;
#if (CH_SPACEDIM == 2)
    grids.push_back(Box(IntVect(0,0), IntVect(3,3)));
    grids.push_back(Box(IntVect(4,0), IntVect(7,3)));
    grids.push_back(Box(IntVect(0,4), IntVect(3,7)));
    grids.push_back(Box(IntVect(4,4), IntVect(7,7)));
#elif (CH_SPACEDIM == 3)
    grids.push_back(Box(IntVect(0,0,0), IntVect(3,3,3)));
    grids.push_back(Box(IntVect(4,0,0), IntVect(7,3,3)));
    grids.push_back(Box(IntVect(0,4,0), IntVect(3,7,3)));
    grids.push_back(Box(IntVect(0,0,4), IntVect(3,3,7)));
#endif
    Vector<int> assign;
    LoadBalance(assign, grids);
    BoxLayout boxlayout (grids, assign);
    boxlayout.close();
    DataIterator lit = boxlayout.dataIterator();

// LayoutData<BaseFab<Real> > test
    LayoutData<BaseFab<Real> > layoutdata(boxlayout);
    for (lit.begin(); lit.ok(); ++lit)
    {
      Box b = boxlayout.get(lit());
      BaseFab<Real>& fab = layoutdata[lit()];
      fab.define(b,SpaceDim);
      BoxIterator bit(b);
      for (bit.begin(); bit.ok(); ++bit)
      {
        IntVect iv = bit();
        for (int d = 0; d < SpaceDim; ++d)
        {
          fab(iv,d) = iv[d];
        }
      }
    }
    pout() << indent << pgmname << ": testing MultiBaseFAB<Real>" << endl ;
    MultiArrayViewRealBaseFab(&layoutdata);

//  BoxLayoutData<FArrayBox> test
    BoxLayoutData<FArrayBox> boxlayoutdata(boxlayout,SpaceDim);
    for (lit.begin(); lit.ok(); ++lit)
    {
      FArrayBox& fab = boxlayoutdata[lit()];
      Box b = fab.box();
      BoxIterator bit(b);
      for (bit.begin(); bit.ok(); ++bit)
      {
        IntVect iv = bit();
        for (int d = 0; d < SpaceDim; ++d)
        {
          fab(iv,d) = iv[d];
        }
      }
    }
    pout() << indent << pgmname << ": testing MultiFAB(boxlayoutdata)" << endl ;
    MultiArrayViewFab(&boxlayoutdata);

    LoadBalance(assign, grids);
    DisjointBoxLayout dboxlayout(grids, assign);
    dboxlayout.close();
    DataIterator dlit = dboxlayout.dataIterator();

//  LevelData<FArrayBox> test
    LevelData<FArrayBox> leveldata(dboxlayout,SpaceDim);
    for (dlit.begin(); dlit.ok(); ++dlit)
    {
      FArrayBox& fab = leveldata[dlit()];
      Box b = fab.box();
      BoxIterator bit(b);
      for (bit.begin(); bit.ok(); ++bit)
      {
        IntVect iv = bit();
        for (int d = 0; d < SpaceDim; ++d)
        {
          fab(iv,d) = iv[d];
        }
      }
    }
    pout() << indent << pgmname << ": testing MultiFAB(leveldata)" << endl ;
    MultiArrayViewFab(&leveldata);

//  LevelData<BaseFab<Real> > test
    LevelData<BaseFab<Real> > leveldatabasefab(dboxlayout,SpaceDim);
    for (dlit.begin(); dlit.ok(); ++dlit)
    {
      BaseFab<Real>& fab = leveldatabasefab[dlit()];
      Box b = fab.box();
      BoxIterator bit(b);
      for (bit.begin(); bit.ok(); ++bit)
      {
        IntVect iv = bit();
        for (int d = 0; d < SpaceDim; ++d)
        {
          fab(iv,d) = iv[d];
        }
      }
    }
    pout() << indent << pgmname << ": testing MultiBaseFAB<Real>(leveldatabasefab)" << endl ;
    MultiArrayViewRealBaseFab(&leveldatabasefab);
}

///
// Parse the standard test options (-p -q) out of the command line.
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
