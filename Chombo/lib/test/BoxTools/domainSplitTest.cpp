#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cstdio>
#include <cmath>
#include "Box.H"
#include "SPMD.H"
#include "REAL.H"
#include "MayDay.H"
#include "IntVectSet.H"
#include "BRMeshRefine.H"
#include "parstream.H"
#include <iostream>
#include <fstream>
#include <string>
#ifdef CH_MPI
#include "mpi.h"
#endif
#include "UsingNamespace.H"

using std::endl;

/// Global variables for test output:
static const char *pgmname = "domainSplitTest" ;
static const char *indent = "   " ,*indent2 = "      " ;
static bool verbose = true ;

const int destProc = 0;

/// Prototypes:

void
parseTestOptions( int argc ,char* argv[] ) ;
/**
   domainSplitTest returns:
   0: all tests passed
   1: boxes coming out of domain split are not in
   the original domain.
   2: boxes coming out of domain split have a side
   that is too long
   3: domainsplit returned an error code
   4: boxes do not cover  domain
  */
int domainSplitTest(void);

int main(int argc, char* argv[])
{
  using std::ostream;
  bool passed = true;
  int icode = 0;

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  parseTestOptions( argc ,argv ) ;
  if ( verbose )
    pout () << indent2 << "Beginning " << pgmname << " ..." << endl ;

  icode = domainSplitTest();

  if (procID() == destProc)
  {
      if (icode != 0)
      {
          pout() << indent << pgmname
               << ": failed with error code " << icode
               << endl;
          passed = false;
      }
      if (passed)
      {
          pout() << indent << pgmname
               << ": passed all tests"
               << endl;
      }
  }

#ifdef CH_MPI
  MPI_Finalize();
#endif

  return icode;
}
int domainSplitTest()
{
  int retflag = 0;
  int domlen = 100;
  int maxlen = 25;
  if (SpaceDim > 3)
    {
      domlen = 50;
      maxlen = 13;
    }

  int blockf = 2;

  Box domain(IntVect::Zero, (domlen-1)*IntVect::Unit);

  Vector<Box> vbox;

  domainSplit(domain, vbox, maxlen, blockf);

  // useful to switch between dense and tree IVS's...
  //IntVectSet::setMaxDense(1000000000);
  IntVectSet ivsDom(domain);
  for (int ibox = 0; ibox < vbox.size(); ibox++)
    {
      const Box& boxloc = vbox[ibox];
      if (!domain.contains(boxloc))
        {
          if (verbose)
            pout() << "error: " << boxloc << "not contained in " << domain  << endl;
          retflag += 1;
        }
      if ( verbose )
        pout() << "subtracting box " << boxloc << " from ivs" << endl;
      ivsDom -= boxloc;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          if (boxloc.size(idir) > maxlen)
            {
              if (verbose)
                pout() << "error: " << boxloc << "is too big in " << idir
                     << " direction (maxlen == " << maxlen << ")"  << endl;
              retflag += 2;
            }
        }
    }
  if (!ivsDom.isEmpty())
    {
      if (verbose)
        pout() << "error: domain not covered " << endl;
      retflag += 4;
    }
  if (verbose)
    {
      pout() << "domain = " << domain << endl;
      pout() << "max box size = " << maxlen << endl;
      pout() << "boxes = " << endl;
      for (int ibox = 0; ibox < vbox.size(); ibox++)
        {
          const Box& boxloc = vbox[ibox];
          pout() << boxloc << endl;
        }
      pout() << "box sizes = " << endl;
      for (int ibox = 0; ibox < vbox.size(); ibox++)
        {
          const Box& boxloc = vbox[ibox];
          for (int idir = 0; idir < SpaceDim; idir ++)
            pout() << boxloc.size(idir)  << "   ";
          pout() << endl;
        }
    }

  // test to see if granularity is screwing things up
  {
    Box testDomainBox(IntVect::Zero, IntVect(D_DECL6(1, 63, 15,
                                                     15,15,15)));
    int ncellmax = 16;
    int block_factor = 2;
    Vector<Box> boxes;

    domainSplit(testDomainBox, boxes, ncellmax, block_factor);

    // check results
    bool passedDomainSplitTest = true;
    if (boxes.size() == 0 )
      {
        // test fails
        if (verbose)
          {
            pout() << "no boxes returned" << endl;
          }
        passedDomainSplitTest = false;
      }

    // check to be sure that boxes cover domain and no box is
    // greater than ncellmax
    IntVectSet uncoveredIVS(testDomainBox);
    for (int n=0; n<boxes.size(); n++)
      {
        if (!testDomainBox.contains(boxes[n]))
          { // fail because boxes[n] contains cells outside testDomainBox
            if (verbose)
              {
                pout() << "box " << boxes[n]
                       << " contains cells outside " << testDomainBox << endl;
              }
            passedDomainSplitTest = false;
          }

        if (!uncoveredIVS.contains(boxes[n]))
          { // fail because boxes[n] contains cells already covered
            if (verbose)
              {
                pout() << "box " << boxes[n]
                       << " contains previously covered cells" << endl;
              }
            passedDomainSplitTest = false;
          }
        uncoveredIVS -= boxes[n];

        const IntVect& lengthVect = boxes[n].size();
        for (int idir = 0; idir < SpaceDim; idir++)
          {
            if (lengthVect[idir] > ncellmax)
              { // fail because boxes[n] is too big
                if (verbose)
                  {
                    pout() << "box " << boxes[n]
                           << " has length > " << ncellmax
                           << " in dimension " << idir << endl;
                  }
                passedDomainSplitTest = false;
              }
            if (lengthVect[idir] % block_factor != 0)
              { // fail because length is not a multiple of block_factor
                if (verbose)
                  {
                    pout() << "box " << boxes[n]
                           << " not blocked by " << block_factor
                           << " in dimension " << idir << endl;
                  }
                passedDomainSplitTest = false;
              }
          }
      }
    if (!uncoveredIVS.isEmpty())
      { // fail because not all of testDomainBox is covered by boxes
        if (verbose)
          {
            pout() << "boxes do not cover all of " << testDomainBox << endl;
          }
        passedDomainSplitTest = false;
      }

    if (verbose)
      pout() << indent << pgmname << " domainSplit thin domain test: "
             << (( passedDomainSplitTest ) ? "passed" : "failed")
             << endl ;
    if (!passedDomainSplitTest) retflag += 8;
  }

  // test that nonzero low end is processed correctly
  {
    IntVect lo = IntVect(D_DECL6(-32, 64, 96, 16, 16, 16));
    int block_factor = 4;
    IntVect num_blocks = IntVect(D_DECL6(16, 12, 12, 2, 2, 2));
    IntVect length_total = block_factor * num_blocks;
    IntVect hi = lo + length_total - IntVect::Unit;
    Box testDomainBox(lo, hi);
    IntVect refineDirs = IntVect(D_DECL6(0, 0, 1, 1, 1, 1));
    int ncellmax = 8 * block_factor;

    Vector<Box> boxes;
    domainSplit(testDomainBox, boxes, ncellmax, block_factor);

    // check results
    bool passedDomainSplitTest = true;
    if (boxes.size() == 0)
      {
        // test fails
        if (verbose)
          {
            pout() << "no boxes returned" << endl;
          }
        passedDomainSplitTest = false;
      }

    // check to be sure that boxes cover domain and no box is
    // greater than ncellmax
    IntVectSet uncoveredIVS(testDomainBox);
    for (int n=0; n<boxes.size(); n++)
      {
        if (!testDomainBox.contains(boxes[n]))
          { // fail because boxes[n] contains cells outside testDomainBox
            if (verbose)
              {
                pout() << "box " << boxes[n]
                       << " contains cells outside " << testDomainBox << endl;
              }
            passedDomainSplitTest = false;
          }
        if (!uncoveredIVS.contains(boxes[n]))
          { // fail because boxes[n] contains cells already covered
            if (verbose)
              {
                pout() << "box " << boxes[n]
                       << " contains previously covered cells" << endl;
              }
            passedDomainSplitTest = false;
          }
        uncoveredIVS -= boxes[n];

        const IntVect& lengthVect = boxes[n].size();
        for (int idir = 0; idir < SpaceDim; idir++)
          {
            if (lengthVect[idir] > ncellmax)
              { // fail because boxes[n] is too big
                if (verbose)
                  {
                    pout() << "box " << boxes[n]
                           << " has length > " << ncellmax
                           << " in dimension " << idir << endl;
                  }
                passedDomainSplitTest = false;
              }
            if (lengthVect[idir] % block_factor != 0)
              { // fail because length is not a multiple of block_factor
                if (verbose)
                  {
                    pout() << "box " << boxes[n]
                           << " not blocked by " << block_factor
                           << " in dimension " << idir << endl;
                  }
                passedDomainSplitTest = false;
              }
          }
      }
    if (!uncoveredIVS.isEmpty())
      { // fail because not all of testDomainBox is covered by boxes
        if (verbose)
          {
            pout() << "boxes do not cover all of " << testDomainBox << endl;
          }
        passedDomainSplitTest = false;
      }

    if (verbose)
      {
        pout() << indent << pgmname << " domainSplit off-zero domain test: "
               << (( passedDomainSplitTest ) ? "passed" : "failed")
               << endl ;
      }
    if (!passedDomainSplitTest) retflag += 16;
  }

  // test that this does the right thing when some dimensions held fixed
  {
    IntVect lo = IntVect(D_DECL6(-32, 64, 96, 16, 16, 16));
    int block_factor = 4;
    IntVect num_blocks = IntVect(D_DECL6(16, 12, 12, 2, 2, 2));
    IntVect length_total = block_factor * num_blocks;
    IntVect hi = lo + length_total - IntVect::Unit;
    Box testDomainBox(lo, hi);
    IntVect refineDirs = IntVect(D_DECL6(0, 0, 1, 1, 1, 1));
    int ncellmax = 8 * block_factor;
    
    Vector<Box> boxes;
    domainSplit(testDomainBox, boxes, ncellmax, block_factor, refineDirs);

    // check results
    bool passedDomainSplitTest = true;
    if (boxes.size() == 0)
      {
        // test fails
        if (verbose)
          {
            pout() << "no boxes returned" << endl;
          }
        passedDomainSplitTest = false;
      }

    // check to be sure that boxes cover domain and no box has length
    // greater than ncellmax in dimensions where refinedDirs == 1.
    // In dimensions where refinedDirs == 0, check that boxes have full length.
    IntVectSet uncoveredIVS(testDomainBox);
    for (int n=0; n<boxes.size(); n++)
      {
        if (!testDomainBox.contains(boxes[n]))
          { // fail because boxes[n] contains cells outside testDomainBox
            if (verbose)
              {
                pout() << "box " << boxes[n]
                       << " contains cells outside " << testDomainBox << endl;
              }
            passedDomainSplitTest = false;
          }

        if (!uncoveredIVS.contains(boxes[n]))
          { // fail because boxes[n] contains cells already covered
            if (verbose)
              {
                pout() << "box " << boxes[n]
                       << " contains previously covered cells" << endl;
              }
            passedDomainSplitTest = false;
          }
        uncoveredIVS -= boxes[n];

        const IntVect& lengthVect = boxes[n].size();
        for (int idir = 0; idir < SpaceDim; idir++)
          {
            if (refineDirs[idir] == 0)
              { // should not be splitting in this dimension
                if (lengthVect[idir] != length_total[idir])
                  { // fail because boxes[n] not full length in this dimension
                    if (verbose)
                      {
                        pout() << "box " << boxes[n]
                               << " not full length " << length_total[idir]
                               << " in dimension " << idir << endl;
                      }
                    passedDomainSplitTest = false;
                  }
              }
            else
              {
                if (lengthVect[idir] > ncellmax)
                  { // fail because boxes[n] is too big in this dimension
                    if (verbose)
                      {
                        pout() << "box " << boxes[n]
                               << " has length > " << ncellmax
                               << " in dimension " << idir << endl;
                      }
                    passedDomainSplitTest = false;
                  }
                if (lengthVect[idir] % block_factor != 0)
                  { // fail because length is not a multiple of block_factor
                    if (verbose)
                      {
                        pout() << "box " << boxes[n]
                               << " not blocked by " << block_factor
                               << " in dimension " << idir << endl;
                      }
                    passedDomainSplitTest = false;
                  }
              }
          }
      }
    if (!uncoveredIVS.isEmpty())
      { // fail because not all of testDomainBox is covered by boxes
        if (verbose)
          {
            pout() << "boxes do not cover all of " << testDomainBox << endl;
          }
        passedDomainSplitTest = false;
      }

    if (verbose)
      {
        pout() << indent << pgmname << " domainSplit anisotropic test: "
               << (( passedDomainSplitTest ) ? "passed" : "failed")
               << endl ;
      }
    if (!passedDomainSplitTest) retflag += 32;
  }

  return retflag;
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
