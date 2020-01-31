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
using std::ifstream;
using std::ios;
#include <string>

#include "FABView.H"

#include "ParmParse.H"
#include "CH_HDF5.H"
#include "parstream.H"

#ifdef CH_USE_TIMER
#include "CH_Timer.H"
#endif

#include "IndexedMoments.H"
#include "MomentIterator.H"

#ifdef CH_MPI
#include "CH_Attach.H"
#endif

#include "NamespaceHeader.H"


static const Real s_tol = 1e-13;
static const Real s_dx = 0.5;

void testIndexedMoments1D();
void testIndexedMoments2D();
void testIndexedMoments3D();

void testMomentIterator1D();
void testMomentIterator2D();
void testMomentIterator3D();

// amrAdvect is a function (as opposed to inline in main()) to get
// around MPI scoping problems
// setupFixedGrids allows fixed grids to be read in and used in this AMR
// computation example
void setupFixedGrids(Vector<Vector<Box> >& a_amrGrids,
                     const ProblemDomain&  a_domain,
                     int                   a_maxLevel,
                     int                   a_maxGridSize,
                     int                   a_blockFactor,
                     int                   a_verbosity,
                     std::string           a_gridFile);

// One more function for MPI
void dumpmemoryatexit();

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

int main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  // Start MPI
  MPI_Init(&a_argc,&a_argv);
  setChomboMPIErrorHandler();
#endif

#ifdef USE_ARRAYVIEW
  trapfpe();
#endif

  /*
  // Check for an input file
  char* inFile = NULL;
  if (a_argc > 1)
  {
    inFile = a_argv[1];
  }
  else
  {
    pout() << "Usage: " << a_argv[0] << " <inputfile>" << endl;
    pout() << "No input file specified" << endl;
    return -1;
  }

  // Parse the command line and the input file (if any)
  ParmParse pp(a_argc-2,a_argv+2,NULL,inFile);
  */

  // Run the tests
  testIndexedMoments1D();
  testMomentIterator1D();
  if (SpaceDim > 1)
    {
      testIndexedMoments2D();
      testMomentIterator2D();
    }
  if (SpaceDim > 2)
    {
      testIndexedMoments3D();
      testMomentIterator3D();
    }

#ifdef CH_MPI
  // Exit MPI
  dumpmemoryatexit();
#endif

  CH_TIMER_REPORT();

#ifdef CH_MPI
  MPI_Finalize();
#endif
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
void testIndexedMoments1D()
{
  pout() << "Index moment tests for 1D" << endl;
  const int D1 = 1;

  {
    const int p = 4;
    IndexedMoments<D1,p> moms;
    IndexedMoments<D1,p> regmoms;
    regmoms.setRegular(s_dx);
    int size = (p+D1);
    pout() << "For p=" << p << ", " << size << " moments." << endl;

    int ix=0;
    Vector<Real> expected(size);
    Vector<Real> regexpected(size);
    Vector<IndexTM<int,D1> > indexes(size);
    // Box iterators will have the right ordering
    Box pbox(IntVect::Zero,IntVect(D_DECL(p,0,0)));
    BoxIterator bit(pbox);
    for (bit.begin(); bit.ok(); ++bit)
    {
      if (bit().sum() > p)
        continue; 

      // calculate a fake moment and store it
      IndexTM<int,D1> itm(bit()[0]);
      indexes[ix] = itm;
      expected[ix] = pow(0.5,bit()[0]);
      moms[ix] = expected[ix];

      // calculate a 1D regular moments
      Real moment;
      if (itm[0] % 2)
        moment = 0;
      else
        moment = pow(0.5*s_dx, itm.sum()) * s_dx / (Real) (itm[0]+1);
      regexpected[ix] = moment;

      pout() << "Index " << ix << "=" << itm << endl;
      ++ix;
    }

    // Check the results vs. expected
    Vector<IndexTM<int,D1> > testix = moms.getMonomialPowers();
    for (ix=0; ix < size; ++ix)
    {
      /**/
      pout() << indexes[ix] << " moment = " 
        << moms[indexes[ix]] 
        << ", expected " << expected[ix] << endl;
      /**/
      CH_assert((testix[ix] - indexes[ix]).sum() == 0);
      CH_assert(fabs(moms[indexes[ix]] - expected[ix]) < s_tol);
      CH_assert(fabs(regmoms[indexes[ix]] - regexpected[ix]) < s_tol);
    }
  }
  pout() << "OK!" << endl;
}

void testMomentIterator1D()
{
  pout() << "moment iterator tests for 1D" << endl;
  const int D1 = 1;


  const int p = 4;
  IndexedMoments<D1,p> indmoms;
  MomentIterator<D1,p> iter;

  Vector<Real> vecmoms(indmoms.size());
  for (int imom = 0; imom < indmoms.size(); imom++)
    {
      Real ans = imom + 1 + imom*imom;
      vecmoms[imom] = ans;
      indmoms[imom] = ans;
    }

  int vecindex = 0;
  for (iter.reset(); iter.ok(); ++iter)
    {
      IndexTM<int, D1> index = iter();
      if (Abs(indmoms[index]-vecmoms[vecindex]) > s_tol)
        {
          MayDay::Error("itertest1d: moments do not match");
        }
      vecindex++;
    }

  if (vecmoms.size() != vecindex) 
    {
      MayDay::Error("iteratortest1D sizes do not seem to match");
    }
  pout() << "OK!" << endl;
}

void testMomentIterator2D()
{
  pout() << "moment iterator tests for 2D" << endl;
  const int D1 = 2;

  const int p = 4;
  IndexedMoments<D1,p> indmoms;
  MomentIterator<D1,p> iter;

  Vector<Real> vecmoms(indmoms.size());
  for (int imom = 0; imom < vecmoms.size(); imom++)
    {
      Real ans = imom + 1 + imom*imom;
      vecmoms[imom] = ans;
      indmoms[imom] = ans;
    }

  int vecindex = 0;
  for (iter.reset(); iter.ok(); ++iter)
    {
      IndexTM<int, D1> index = iter();
      if (Abs(indmoms[index]-vecmoms[vecindex]) > s_tol)
        {
          MayDay::Error("itertest1d: moments do not match");
        }
      vecindex++;
    }

  if (vecmoms.size() != vecindex) 
    {
      MayDay::Error("iteratortest1D sizes do not seem to match");
    }
  pout() << "OK!" << endl;
}

void testMomentIterator3D()
{
  pout() << "moment iterator tests for 3D" << endl;
  const int D1 = 3;

  const int p = 4;
  IndexedMoments<D1,p> indmoms;
  MomentIterator<D1,p> iter;

  Vector<Real> vecmoms(indmoms.size());
  for (int imom = 0; imom < vecmoms.size(); imom++)
    {
      Real ans = imom + 1 + imom*imom;
      vecmoms[imom] = ans;
      indmoms[imom] = ans;
    }

  int vecindex = 0;
  for (iter.reset(); iter.ok(); ++iter)
    {
      IndexTM<int, D1> index = iter();
      if (Abs(indmoms[index]-vecmoms[vecindex]) > s_tol)
        {
          MayDay::Error("itertest1d: moments do not match");
        }
      vecindex++;
    }

  if (vecmoms.size() != vecindex) 
    {
      MayDay::Error("iteratortest1D sizes do not seem to match");
    }
  pout() << "OK!" << endl;
}

void testIndexedMoments2D()
{
  pout() << "Index moment tests for 2D" << endl;
  const int D2 = 2;

  {
    const int p = 4;
    IndexedMoments<D2,p> moms;
    IndexedMoments<D2,p> regmoms;
    regmoms.setRegular(s_dx);
    int size = (p+D2)*(p+D2-1)/2;
    pout() << "For p=" << p << ", " << size << " moments." << endl;

    int ix=0;
    Vector<Real> expected(size);
    Vector<Real> regexpected(size);
    Vector<IndexTM<int,D2> > indexes(size);
    // Box iterators will have the right ordering
    Box pbox(IntVect::Zero,IntVect(D_DECL(p,p,0)));
    BoxIterator bit(pbox);
    for (bit.begin(); bit.ok(); ++bit)
    {
      if (bit().sum() > p)
        continue; 

      // calculate a fake moment and store it
      IndexTM<int,D2> itm(bit()[0], bit()[1]);
      indexes[ix] = itm;
      expected[ix] = pow(0.5,bit()[0])*pow(1.5,bit()[1]);
      moms[ix] = expected[ix];

      // calculate a 2D regular moments
      Real moment;
      if ((itm[0] % 2) || (itm[1] % 2))
        moment = 0;
      else
      {
        moment = pow(0.5*s_dx, itm.sum()) * pow(s_dx, D2);
        moment /= (Real) (itm[0]+1)*(itm[1]+1);
      }
      regexpected[ix] = moment;

      pout() << "Index " << ix << "=" << itm << endl;
      ++ix;
    }

    // Check the results vs. expected
    Vector<IndexTM<int,D2> > testix = moms.getMonomialPowers();
    for (ix=0; ix < size; ++ix)
    {
      /**/
      pout() << indexes[ix] << " moment = " 
        << moms[indexes[ix]] 
        << ", expected " << expected[ix] << endl;
      /**/
      CH_assert((testix[ix] - indexes[ix]).sum() == 0);
      CH_assert(fabs(moms[indexes[ix]] - expected[ix]) < s_tol);
      CH_assert(fabs(regmoms[indexes[ix]] - regexpected[ix]) < s_tol);
    }
  }
  pout() << "OK!" << endl;
}

void testIndexedMoments3D()
{
  pout() << "Index moment tests for 3D" << endl;
  const int D3 = 3;


  {
    const int p = 4;
    IndexedMoments<D3,p> moms;
    IndexedMoments<D3,p> regmoms;
    regmoms.setRegular(s_dx);
    int size = (p+D3)*(p+D3-1)*(p+D3-2)/(2*3);
    pout() << "For p=" << p << ", " << size << " moments" << endl;

    int ix=0;
    Vector<Real> expected(size);
    Vector<Real> regexpected(size);
    Vector<IndexTM<int,D3> > indexes(size);
    // Box iterators will have the right ordering
    Box pbox(IntVect::Zero,IntVect(D_DECL(p,p,p)));
    BoxIterator bit(pbox);
    for (bit.begin(); bit.ok(); ++bit)
    {
      if (bit().sum() > p)
        continue; 

      // calculate a fake moment and store it
      IndexTM<int,D3> itm(bit()[0], bit()[1], bit()[2]);
      indexes[ix] = itm;
      expected[ix] = pow(0.5,bit()[0])*pow(1.5,bit()[1])*pow(0.2,bit()[2]);
      moms[ix] = expected[ix];

      // calculate a 2D regular moments
      Real moment;
      if ((itm[0] % 2) || (itm[1] % 2) || (itm[2] % 2))
        moment = 0;
      else
      {
        moment = pow(0.5*s_dx, itm.sum()) * pow(s_dx, D3);
        moment /= (Real) (itm[0]+1)*(itm[1]+1)*(itm[2]+1);
      }
      regexpected[ix] = moment;

      pout() << "Index " << ix << "=" << itm << endl;
      ++ix;
    }
    CH_assert(ix == size);

    // Check the results vs. expected
    Vector<IndexTM<int,D3> > testix = moms.getMonomialPowers();
    for (ix=0; ix < size; ++ix)
    {
      /**/
      pout() << indexes[ix] << " moment = " 
        << moms[indexes[ix]] 
        << ", expected " << expected[ix] << endl;
      /**/
      CH_assert((testix[ix] - indexes[ix]).sum() == 0);
      CH_assert(fabs(moms[indexes[ix]] - expected[ix]) < s_tol);
      CH_assert(fabs(regmoms[indexes[ix]] - regexpected[ix]) < s_tol);
    }
  }
  pout() << "OK!" << endl;
}

#include "NamespaceFooter.H"
