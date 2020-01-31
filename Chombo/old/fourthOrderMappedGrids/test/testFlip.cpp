#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <fstream>
#include <iostream>
using std::cout;
using std::endl;
using std::ofstream;


#include "parstream.H"
#include "flipGrids.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "DataIterator.H"
#include "BoxIterator.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "FABView.H"

/// Global variables for handling output:
static const char* pgmname = "testFlip" ;
static const char* indent = "   ";
static const char* indent2 = "      " ;
//static bool verbose = false ;
static bool verbose = true ;

#ifdef CH_USE_DOUBLE
static Real precision = 1.0e-15;
#else
static Real precision = 1.0e-7;
#endif

///
// Parse the standard test options (-v -q -h) and
// app-specific options (-S <domain_size>) out of the command line
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
              //argv[i] = "" ;
            }
          else if ( strncmp( argv[i] ,"-q" ,3 ) == 0 )
            {
              verbose = false ;
              //argv[i] = "" ;
            }
          else if ( strncmp( argv[i] ,"-h" ,3 ) == 0 )
            {
              pout() << "usage: " << pgmname << " [-hqv]" << std::endl ;
              exit( 99 ) ;
            }
        }
    }
  return;
}


int
testFlip();

int
main(int argc ,char* argv[])
{
#ifdef CH_MPI
  MPI_Init (&argc, &argv);
#endif
  parseTestOptions( argc ,argv ) ;
  if ( verbose )
    pout () << indent2 << "Beginning " << pgmname << " ..." << endl ;

  int status = testFlip();

  if ( status == 0 )
    cout << indent << "Flip test" << " passed." << endl ;
  else
    cout << indent << "Flip test" << " failed with return code "
         << status << endl ;



  if ( status == 0 )
    cout << indent << pgmname << " passed." << endl ;
  else
    cout << indent << pgmname << " failed with return code " << status << endl ;


#ifdef CH_MPI
  MPI_Finalize ();
#endif
  return status ;
}

int
testFlip()
{
  int status = 0;

  // set up boxLayout which straddles the x=0 plane...
  Box domainBox(-16*IntVect::Unit, 15*IntVect::Unit);
  ProblemDomain entireDomain(domainBox);
  entireDomain.setPeriodic(1,true);

  int maxBoxSize = 8;
  Vector<Box> gridBoxes;
  domainSplit(domainBox, gridBoxes, maxBoxSize);

  Vector<int> procAssign(gridBoxes.size(), 0);
  LoadBalance(procAssign, gridBoxes);

  DisjointBoxLayout grids(gridBoxes, procAssign, entireDomain);

  int reflectDir = 0;
  int reflectCoord = 0;
  Vector<Tuple<DataIndex,2> > boxCorrelation;
  DisjointBoxLayout reflectedBoxes;

  // this is the box we want to have reflected data for
  Box dataBox(4*IntVect::Unit, 11*IntVect::Unit);

  getFlippedGrids(reflectedBoxes, boxCorrelation,
                  grids, dataBox, reflectDir, reflectCoord);


  // now test results. first loop over all boxes in grids
  {
    // note that in this test, all cells in the dataBox
    // should be contained somewhere in the reflectedBoxes
    IntVectSet testIVS(dataBox);

    DataIterator dit = grids.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
      {
        Box intersectBox = grids[dit];
        intersectBox &= dataBox;
        if (!intersectBox.isEmpty() )
          {
            DataIndex reflectDit;
            bool foundMatch = false;
            for (int n=0; n<boxCorrelation.size(); n++)
              {
                Tuple<DataIndex,2> thisCorrelation = boxCorrelation[n];
                if (thisCorrelation[0] == dit())
                  {
                    reflectDit = thisCorrelation[1];
                    foundMatch = true;
                  }
              }
            if (foundMatch == false)
              {
                if (verbose)
                  {
                    pout() << "   box correlation not found! " << endl;
                  }
                status += 1;
              }
            else
              {
                const Box reflectBox = reflectedBoxes[reflectDit];

                BoxIterator bit(intersectBox);
                for (bit.begin(); bit.ok(); ++bit)
                  {
                    // location in original box
                    IntVect iv = bit();
                    testIVS -= iv;
                    // corresponding location in reflected box
                    iv[0] = -1*iv[0];

                    if (!reflectBox.contains(iv))
                      {
                        if (verbose)
                          {
                            pout() << "   reflected point " << iv
                                   << " not found in reflected box!"
                                   << endl;
                          }
                        status += 100;
                      } // end if reflected box contains the reflected point
                  } // end loop over points in reflect Box
              } // end if we found a box correlation

          } // end if this grid intersects the dataBox
      } // end loop over grids

    // this test also only works in serial
    if (!testIVS.isEmpty() && (numProc() ==1) )
      {
        if (verbose)
          {
            pout() << "    all reflection points not contained" << endl;
          }
        status += 1000;
      }
  } // end testing results

  {
    // make sure we can do the copyTo
    LevelData<FArrayBox> data(grids, 1);
    LevelData<FArrayBox> reflectedData(reflectedBoxes, 1);

    // initialize data
    DataIterator dit = data.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
      {
        BoxIterator bit(data[dit].box());
        for (bit.begin(); bit.ok(); ++bit)
          {
            IntVect iv = bit();
            data[dit](iv,0) = iv[0];
          }
      }

    data.copyTo(reflectedData);

    // now check to be sure that copyTo went ok
    DataIterator newDit = reflectedData.dataIterator();
    for (newDit.begin(); newDit.ok(); ++newDit)
      {
        BoxIterator bit(reflectedData[newDit].box());
        for (bit.begin(); bit.ok(); ++bit)
          {
            IntVect iv = bit();
            if (reflectedData[newDit](iv,0) != iv[0])
              {
                if (verbose)
                  {
                    pout() << "   reflected data fails data copyTo test at "
                           << iv << endl;
                  }
                status += 5000;
              }
          }
      }
  } // end copyTo testing

  return status;

}
