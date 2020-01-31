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
//  Test the infrastucture for particles
//
// Usage:
//  <program-name> [-q|-v] ...
//
//  where:
//    -q means run quietly (only pass/fail messages printed)
//    -v means run verbosely (all messages printed)
//    -writePlots meanse write out the results in a set of hdf5 plotfiles
//    ... all non-option arguments are ignored (Chombo convention)
//
//  Default is `-v'
//
//  Reading in the arguments is terminated if a non-recognized option is passed.
//

// Include files:
#include <cstdio>
#include <string.h>
#include <assert.h>

#include "parstream.H"
#include "BinFab.H"
#include "BinFabFactory.H"
#include "BinItem.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "RealVect.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "AMRIO.H"

#ifdef CH_MPI
#include "mpi.h"
#endif

#include "UsingNamespace.H"

//////////////////////////////////////////////////////////////
using std::endl;

void parseTestOptions(int argc, char* argv[]);

void initData(BinFab<BinItem>& a_data, const Real a_dx);

void convertParticlesToReals(const BinFab<BinItem>& a_bins,
                             FArrayBox& a_data);

int moveParticles(BinFab<BinItem>& a_particles, const Real a_dt);

void print( BinFab<BinItem> & a_bf );

int countItems( BinFab<BinItem> & a_bf );

/// Global variables for handling output

static const char *pgmname = "testBinFab";
static const char *indent = "   ";

static bool verbose = true ;
static bool writePlotFiles = false;

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

    //int baseDomainSize = 8;
    int baseDomainSize = 10;
    Real dx = 1.0/baseDomainSize;
    RealVect meshSpacing(D_DECL(dx,dx,dx));
    RealVect origin(D_DECL(0,0,0));
    //int numGhost = 0;
    //IntVect ghostVect(numGhost*IntVect::Unit);

    Box baseDomBox(IntVect::Zero, (baseDomainSize-1)*IntVect::Unit);

    ProblemDomain baseDomain(baseDomBox);

    { // first try a single box on this level
      BinFab<BinItem> thisBinFab;
      thisBinFab.define(baseDomBox, meshSpacing, origin);

      initData(thisBinFab, dx);

      FArrayBox dumpFab(baseDomBox,1);
      if (writePlotFiles)
      {
        convertParticlesToReals(thisBinFab, dumpFab);

#ifdef CH_USE_HDF5
        writeFABname(&dumpFab, "firstTest.0.hdf5");
#endif
      }

      // now test copy function
     BinFab<BinItem> copyBinFab(baseDomBox, meshSpacing, origin);
     Interval copyInterval(0,0);
     copyBinFab.copy(baseDomBox, copyInterval, baseDomBox, thisBinFab,
                     copyInterval);

      if (writePlotFiles)
      {
        dumpFab.setVal(0.0);
        //    writeFABname(&dumpFab, "zeroTest.hdf5");
        convertParticlesToReals(copyBinFab, dumpFab);
#ifdef CH_USE_HDF5
        writeFABname(&dumpFab, "copyTest.hdf5");
#endif
      }

      {
        // now test linearIn/Out
        static const int NB=10000000 ;
        char *buf = new char[NB];
        char *buf2 = new char[NB];
        for ( int i=0 ; i<NB ; ++i )
        {
          buf[i] = '@' ;
        }
        thisBinFab.linearOut( buf ,thisBinFab.box() ,thisBinFab.interval() ) ;
        if ( verbose ) print( thisBinFab ) ;
        BinFab<BinItem> thatBinFab( baseDomBox, meshSpacing, origin ) ;
        thatBinFab.linearIn( buf, thisBinFab.box() ,thisBinFab.interval() ) ;
        // add second linearIn to make sure we're not appending instead of
        // over-writing -- will cause test to fail if we're appending,
        // but should have no difference if we're overwriting what's there.
        thatBinFab.linearIn( buf, thisBinFab.box() ,thisBinFab.interval() ) ;
        if ( verbose ) print( thatBinFab ) ;
        for ( int i=0 ; i<NB ; ++i )
        {
          buf2[i] = '@' ;
        }
        thatBinFab.linearOut( buf2 ,thisBinFab.box() ,thisBinFab.interval() ) ;
        int is_diff = strncmp( buf ,buf2 ,NB ) ;
        if ( is_diff != 0 && verbose )
        {
          pout() << " fail: linearOut()/linearIn() are not symmetric" << endl ;
        }
        status += is_diff ;
        int bfsize  = thisBinFab.size( thisBinFab.box() ,thisBinFab.interval() ) ;
        int bfsize2 = thatBinFab.size( thisBinFab.box() ,thisBinFab.interval() ) ;
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

	thisBinFab.linearOut( test ,thisBinFab.box() ,thisBinFab.interval() ) ;
	thatBinFab.linearOut( test2 ,thisBinFab.box() ,thisBinFab.interval() ) ;

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

      int numSteps = 10;
      Real dt = 0.1;
      for (int step=0; step<numSteps; step++)
      {
        int lost_p = moveParticles(thisBinFab, dt);
        if ( lost_p > 0 && verbose )
        {
          pout() << " info: firstTest: moveParticles() lost " << lost_p << "particles." << endl ;
        }
        if (writePlotFiles)
        {
          convertParticlesToReals(thisBinFab, dumpFab);

          char iter_str[80];
          sprintf(iter_str, "%s%d.hdf5", "firstTest.", step+1);
#ifdef CH_USE_HDF5
          writeFABname(&dumpFab, iter_str);
#endif
        }
      }
    }
    {
      // now try with a LevelData of these things.
      int maxBoxSize = 5;
      Vector<Box> boxes;
      domainSplit(baseDomain, boxes, maxBoxSize);

      Vector<int> procAssign;
      int eekflag = LoadBalance(procAssign, boxes);
      assert(eekflag == 0);

      DisjointBoxLayout grids(boxes, procAssign, baseDomain);

      BinFabFactory<BinItem > factory(meshSpacing, origin);

      DataIterator dit = grids.dataIterator();

      LevelData<BinFab<BinItem> > levelBins(grids, 1,
                                            IntVect::Unit,
                                            factory);

      LevelData<FArrayBox> dumpData(grids, 1);

      if (writePlotFiles)
      {
        for (dit.begin(); dit.ok(); ++dit)
        {
          initData(levelBins[dit()], dx);
          convertParticlesToReals(levelBins[dit()], dumpData[dit()]);
        }

        // dump initial data
#ifdef CH_USE_HDF5
        writeLevelname(&dumpData, "binDump.00.hdf5");
#endif
      }

      // copy test
      LevelData<BinFab<BinItem> > copyTester(grids, 1,
                                             IntVect::Unit,
                                             factory);

      levelBins.copyTo(levelBins.interval(), copyTester, levelBins.interval());

      if (writePlotFiles)
      {
        // dump tester
        for (dit.begin(); dit.ok(); ++dit)
        {
          convertParticlesToReals(copyTester[dit()], dumpData[dit()]);
        }

#ifdef CH_USE_HDF5
        writeLevelname(&dumpData, "copyDump.hdf5");
#endif
      }

      // now move particles

      int numSteps = 10;
      // need to be this small to satisfy CFL condition imposed by
      // ghost-layer size
      Real dt = 0.04;
      for (int step = 0; step < numSteps; step++)
      {
        // do exchange
        levelBins.exchange(levelBins.interval());
        for (dit.begin(); dit.ok(); ++dit)
        {
          BinFab<BinItem>& thisBinFab = levelBins[dit()];
          FArrayBox thisDataFab(thisBinFab.box(),1);
          convertParticlesToReals(thisBinFab, thisDataFab);

          int lost_p = moveParticles(levelBins[dit()],dt);
          if ( lost_p > 0 && verbose )
          {
            pout() << " info: LevelData test: moveParticles() lost " << lost_p << "particles." << endl ;
          }

          convertParticlesToReals(thisBinFab, thisDataFab);
          convertParticlesToReals(levelBins[dit()], dumpData[dit()]);
        }

        if (writePlotFiles)
        {
          // now dump data
          char iter_str[80];
          if (step < 9)
            sprintf(iter_str, "%s%d.hdf5", "binDump.0", step+1);
          else
            sprintf(iter_str, "%s%d.hdf5", "binDump.", step+1);
#ifdef CH_USE_HDF5
          writeLevelname(&dumpData, iter_str);
#endif
        }
      } // end loop over advance steps
    } // end LevelData test

    // addItemsDestructive tests:
    {
      { // Create one particle for every cell and add them destructively:
        // the list should be returned empty and the BinFab should have all
        // the particles.
        CH_XD::List<BinItem> particles_l ;
        RealVect halfcell(meshSpacing/2.0) ;
        Box box(IntVect::Zero, (baseDomainSize-1)*IntVect::Unit) ;
        BinItem item ;
        int np = 0 ;
        for ( BoxIterator bi(box) ; bi.ok() ; ++bi )
        {
          item.setPosition( bi() * meshSpacing + halfcell );
          particles_l.append( item );
          ++np ;
        }
        if ( verbose ) pout() << indent << "created " << np << " particles..." ;
        // check for failure
        if ( np != box.volume() )
        {
          ++status ;
          pout() << endl << "error: expecting " << box.volume() << " particles." << endl ;
        }
        BinFab<BinItem> thisBinFab;
        thisBinFab.define(box, meshSpacing, origin);
        thisBinFab.addItemsDestructive( particles_l );
        int nbf = countItems( thisBinFab ) ;
        if ( verbose )
          pout() << " added " << nbf << " particles, " << particles_l.length() << " left over..." ;
        // check for failure
        if ( nbf != np )
        {
          ++status ;
          pout() << endl << "error: expecting " << np << " added particles." << endl ;
        }
        if ( particles_l.isNotEmpty() )
        {
          ++status ;
          if ( nbf == np ) pout() << endl ;
          pout() << "error: expecting 0 particles left over." << endl ;
        }
      }
      { // Create one particle for every cell and add them destructively to a sub-box
        CH_XD::List<BinItem> particles_l ;
        RealVect halfcell(meshSpacing/2.0) ;
        Box box(IntVect::Zero, (baseDomainSize-1)*IntVect::Unit) ;
        Box subbox( grow( box, -2 ) ) ;
        BinItem item ;
        int np = 0 ;
        for ( BoxIterator bi(box) ; bi.ok() ; ++bi )
        {
          item.setPosition( bi() * meshSpacing + halfcell );
          particles_l.append( item );
          ++np ;
        }
        if ( verbose ) pout() << indent << "created " << np << " particles..." ;
        BinFab<BinItem> thisBinFab;
        thisBinFab.define(box, meshSpacing, origin);
        // remove only particles in the given Box
        thisBinFab.addItemsDestructive( particles_l ,subbox );
        int nbf = countItems( thisBinFab ) ;
        if ( verbose ) pout() << " added " << nbf << " particles, " << particles_l.length() << " left over..." ;
        // check for failure
        if ( nbf != np )
        {
          ++status ;
          pout() << "error: expecting " << np << " added particles." << endl ;
        }
        int nsb = subbox.volume() ;
        int nnsb = box.volume() - nsb ;
        if ( particles_l.length() != nnsb )
        {
          ++status ;
          pout() << "error: expecting " << nnsb << " particles left over." << endl ;
        }
      }
      { // Create one particle for every cell and add them destructively outside a sub-box
        CH_XD::List<BinItem> particles_l ;
        RealVect halfcell(meshSpacing/2.0) ;
        Box box(IntVect::Zero, (baseDomainSize-1)*IntVect::Unit) ;
        Box subbox( grow( box, -2 ) ) ;
        BinItem item ;
        int np = 0 ;
        for ( BoxIterator bi(box) ; bi.ok() ; ++bi )
        {
          item.setPosition( bi() * meshSpacing + halfcell );
          particles_l.append( item );
          ++np ;
        }
        if ( verbose ) pout() << indent << "created " << np << " particles..." ;
        BinFab<BinItem> thisBinFab;
        thisBinFab.define(box, meshSpacing, origin);
        // remove only particles NOT in the given Box
        thisBinFab.addItemsDestructive( particles_l ,subbox ,false );
        int nbf = countItems( thisBinFab ) ;
        if ( verbose )
          pout() << " added " << nbf << " particles, " << particles_l.length() << " left over." << endl ;
        // check for failure
        if ( nbf != np )
        {
          ++status ;
          pout() << endl << "error: expecting " << np << " added particles." << endl ;
        }
        int nsb = subbox.volume() ;
        if ( particles_l.length() != nsb )
        {
          ++status ;
          pout() << "error: expecting " << nsb << " particles left over." << endl ;
        }
      }
    }
    // done
    pout() << indent << pgmname << ": "
           << ( (status == 0) ? "passed all tests" : "failed at least one test,")
           << endl;
  }

#ifdef CH_MPI
  MPI_Finalize();
#endif

  return status ;
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
              // argv[i] = "" ;
            }
          else if ( strncmp( argv[i] ,"-q" ,3 ) == 0 )
            {
              verbose = false ;
              // argv[i] = "" ;
            }
          else if ( strncmp( argv[i] ,"-writePlots" ,3 ) == 0 )
            {
#ifdef CH_USE_HDF5
              writePlotFiles = true ;
#else
              pout() << "warning: " << pgmname
                   << ": ignoring '-writePlots' option because HDF5 not enabled" << endl ;
#endif
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

//----------------------------------------------------------------

void
initData(BinFab<BinItem>& a_data, const Real a_dx)
{
  // at the moment, place one particle in a_data
  RealVect position(D_DECL(0.2,0.2,0.2));
  RealVect increment(D_DECL(0.05,0.05,0.0));
  //RealVect increment(D_DECL(0.01,0.01,0.0));
  BinItem particle(position);
  CH_XD::List<BinItem> thisList;
  thisList.append(particle);
  // now add additional particles
  //int numExtraParticles = 10;
  int numExtraParticles = 20;
  for (int n=0; n<numExtraParticles; n++)
  {
    particle.position() = increment + particle.position();
    thisList.append(particle);
  }

  a_data.addItems(thisList);

}

//----------------------------------------------------------------

void
convertParticlesToReals(const BinFab<BinItem>& a_bins,
                        FArrayBox& a_data)
{
  Box intersectBox(a_data.box());
  intersectBox &= a_bins.box();

  BoxIterator bit(intersectBox);

  for (bit.reset(); bit.ok(); ++bit)
    {
      const CH_XD::List<BinItem>& thisList = a_bins(bit(), 0);
      Real numHere = (Real) thisList.length();
      a_data(bit()) = numHere;
    }
}

//----------------------------------------------------------------

int
moveParticles(BinFab<BinItem>& a_particles, const Real a_dt)
{
  BoxIterator bit(a_particles.box());
  RealVect velocity(D_DECL(1.0,0.5,0.0));
  //RealVect velocity(D_DECL(1.0,0.0,0.0));

  RealVect amountToMove = a_dt*velocity;

  for (bit.reset(); bit.ok(); ++bit)
    {
      CH_XD::List<BinItem>& thisList = a_particles(bit(),0);

      if (thisList.length() > 0)
        {
          // now loop over the items in the list and move them
          ListIterator<BinItem> lit(thisList);
          for (lit.rewind(); lit; ++lit)
            {
              BinItem& thisParticle = thisList[lit];
              thisParticle.position() = thisParticle.position() + amountToMove;
            }
        }
    }
  // finally, do rebinning
  a_particles.reBin();
  return 0;

}

//----------------------------------------------------------------

void
print( BinFab<BinItem> & bf )
{
  pout() << "(BinFab " << bf.box() ;
  for ( BoxIterator i(bf.box()) ; i.ok() ; ++i )
  {
    if ( bf(i(),0).length() > 0 )
    {
      IntVect iv = i() ;
      pout() << " (" << D_TERM(iv[0], << "," << iv[1], << "," << iv[2] ) << ") "
           << bf(i(),0).length() ;
      for ( ListIterator<BinItem> p(bf(iv,0)) ; p.ok() ; ++p )
      {
        pout() << p() ;
      }
    }
  }
  pout() << ")" << endl ;
}

//----------------------------------------------------------------

int countItems( BinFab<BinItem> & a_bf )
{
  int count = 0 ;
  for ( BoxIterator i(a_bf.box()) ; i.ok() ; ++i )
  {
    count += a_bf(i(),0).length() ;
  }
  return count ;
}
