#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// Purpose:
//
//  Test the ListBox data holder for particles
//
// Usage:
//
//  <program-name> [-q|-v] ...
//
//  where:
//    -q means run quietly (only pass/fail messages printed)
//    -v means run verbosely (all messages printed)
//
//  Default is `-v'
//
//  Reading in the arguments is terminated if a non-recognized option is passed.
//

// Include files:
#include <cstdio>
#include <string.h>

#include "parstream.H"
#include "ListBox.H"
#include "ListBoxFactory.H"
#include "BinItem.H"
#include "RealVect.H"

#ifdef CH_MPI
#include "mpi.h"
#endif

#include "UsingNamespace.H"

//////////////////////////////////////////////////////////////
using std::endl;

void parseTestOptions(int argc, char* argv[]);

void initData(ListBox<BinItem>& a_data, const Real a_dx);

void print( ListBox<BinItem> & a_bf );

/// Global variables for handling output

static const char *pgmname = "testListBox";
static const char *indent = "   ";

static bool verbose = true ;

/// Code:

int
main(int argc ,char *argv[] )
{
  int status = 0;

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  {
    parseTestOptions(argc, argv);

    int boxSize = 10;
    Real dx = 1.0/boxSize;
    RealVect meshSpacing(D_DECL(dx,dx,dx));
    RealVect origin(D_DECL(0,0,0));

    Box baseDomain(IntVect::Zero, (boxSize - 1) * IntVect::Unit);

    { // first try a single box on this level
      ListBox<BinItem> thisListBox;
      thisListBox.define(baseDomain, meshSpacing, origin);

      initData(thisListBox, dx);

      // now test copy function
      ListBox<BinItem> copyListBox(baseDomain, meshSpacing, origin);
      Interval copyInterval(0,0);
      copyListBox.copy(baseDomain, copyInterval, baseDomain, thisListBox,
                      copyInterval);

      if ( verbose ) print( thisListBox );
      if ( verbose ) print( copyListBox ) ;

      {
        // now test linearIn/Out
        static const int NB=10000000 ;
        char *buf = new char[NB];
        char *buf2 = new char[NB];
        for ( int i=0 ; i<NB ; ++i )
        {
          buf[i] = '@' ;
        }
        thisListBox.linearOut( buf ,thisListBox.box()) ;
        if ( verbose ) print( thisListBox ) ;
        ListBox<BinItem> thatListBox( baseDomain, meshSpacing, origin ) ;
        thatListBox.linearIn( buf, thisListBox.box()) ;
        if ( verbose ) print( thatListBox ) ;
        for ( int i=0 ; i<NB ; ++i )
        {
          buf2[i] = '@' ;
        }
        thatListBox.linearOut( buf2 ,thisListBox.box()) ;
        int is_diff = strncmp( buf ,buf2 ,NB ) ;
        if ( is_diff != 0 && verbose )
        {
          pout() << " fail: linearOut()/linearIn() are not symmetric" << endl ;
        }
        status += is_diff ;
        int bfsize  = thisListBox.size( thisListBox.box()) ;
        int bfsize2 = copyListBox.size( thisListBox.box()) ;
        if ( bfsize != bfsize2 )
        {
          ++status ;
          if ( verbose )
          {
            pout() << "fail: linear size()s are different: " << bfsize << "," << bfsize2 << endl ;
          }
        }

	// to count the number of bytes that get written out, we count backwards from the
	// last end of the buffer, stopping at the last non-@ character.
        int at ,at2;
        for ( at=NB-1  ; at>=0  ; --at  )
        {
          if ( buf [at ] != '@' ) break ;
        }
        for ( at2=NB-1 ; at2>=0 ; --at2 )
        {
          if ( buf2[at2] != '@' ) break ;
        }
        if ( ( at == 0 || at2 == 0 ) )
        {
          ++status ;
          if ( verbose )
          {
            pout() << "fail: linear overran buffer" << endl ;
          }
        }

	// That worked, unless the last byte written out is '@' when interpreted as a char.
	// To handle that case we must also use a second char array with a different token.
        char *test = new char[NB];
        char *test2 = new char[NB];

        for ( int i=0 ; i<NB ; ++i )
        {
          test[i] = '?';
	  test2[i] = '?';
        }

        thisListBox.linearOut( test ,thisListBox.box());
        thatListBox.linearOut( test2, thisListBox.box()) ;	

        int bt ,bt2;
        for ( bt=NB-1  ; bt>=0  ; --bt  )
        {
          if ( test [bt ] != '?' ) break ;
        }
        for ( bt2=NB-1 ; bt2>=0 ; --bt2 )
        {
          if ( test2[bt2] != '?' ) break ;
        }
        if ( ( bt == 0 || bt2 == 0 ) )
        {
          ++status ;
          if ( verbose )
          {
            pout() << "fail: linear overran buffer" << endl ;
          }
        }

	// the last character can't be both '@' and '?', so we use the largest.
	int outsize = max(at + 1, bt + 1);
	int outsize2 = max(at2 + 1, bt2 + 1);

        if ( outsize != bfsize )
        {
          ++status ;
          if ( verbose )
          {
            pout() << "fail: first linearOut() disagrees with size(): " << outsize << "," << bfsize << endl ;
          }
        }
        if ( outsize2 != bfsize2 )
        {
          ++status ;
          if ( verbose )
          {
            pout() << "fail: second linearOut() disagrees with size(): " << outsize2 << "," << bfsize2 << endl ;
          }
        }
        delete [] buf;
        delete [] buf2;
      }
    }
    // done
    pout() << indent << pgmname << ": "
           << ( (status == 0) ? "passed all tests" : "failed at least one test,")
           << endl;
#ifdef CH_MPI
  MPI_Finalize();
#endif

  return status ;
  }
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

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
            }
          else if ( strncmp( argv[i] ,"-q" ,3 ) == 0 )
            {
              verbose = false ;
            }
          else
            {
              break ;
            }
        }
    }
  return ;
}
//----------------------------------------------------------------

void
initData(ListBox<BinItem>& a_data, const Real a_dx)
{
  // at the moment, place one particle in a_data
  RealVect position(D_DECL(0.2,0.2,0.2));
  RealVect increment(D_DECL(0.05,0.05,0.0));
  BinItem particle(position);
  CH_XD::List<BinItem> thisList;
  thisList.append(particle);
  // now add additional particles
  int numExtraParticles = 20;
  for (int n=0; n<numExtraParticles; n++)
  {
    particle.position() = increment + particle.position();
    thisList.append(particle);
  }

  a_data.addItemsDestructive(thisList);
}

void
print( ListBox<BinItem>& lb )
{
  pout() << "(ListBox " << lb.box() << ' ';

      // loop through items, and add to rlist
      for (ListIterator<BinItem> li( lb.listItems() ); ++li; )
        {
          pout() << li() << ' ';
        }

  pout() << ")" << endl ;
}
