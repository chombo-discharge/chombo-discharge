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
//  Test the infrastucture for particles. This tests the fillGhosts
//  function of ParticleData. It creates a ParticleData, calls fillGhosts,
//  and makes sure that you end up with the right number of extra particles
//  on each ListBox.
//   
// Usage:
//  <program-name> [-q|-v] ...
//
//  where:
//    -q means run quietly (only pass/fail messages printed)
//    -v means run verbosely (all messages printed)
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
#include "ListBox.H"
#include "ListBoxFactory.H"
#include "BinItem.H"
#include "RealVect.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "ParticleData.H"

#ifdef CH_MPI
#include "mpi.h"
#endif

#include "UsingNamespace.H"

//////////////////////////////////////////////////////////////
using std::endl;

void parseTestOptions(int argc, char* argv[]);

void initData(ParticleData<BinItem>& a_data, 
              const Real a_dx,
              const int a_numParticles);

/// Global variables for handling output

static const char *pgmname = "testGhostParticles";
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

    int fixedBoxSize = 32;
    int domainDimension = 256;
    int numGhost = 4;

    int partPerBox = pow((double) (fixedBoxSize + 2*numGhost), CH_SPACEDIM);

    Real dx = 1.0 / domainDimension;
    RealVect meshSpacing(D_DECL(dx, dx, dx));
    RealVect origin(D_DECL(0, 0, 0));

    Box domainBox(IntVect::Zero, (domainDimension - 1) * IntVect::Unit);
    const bool isPeriodic[SpaceDim] = {D_DECL(true, true, true)};
    ProblemDomain probDomain(domainBox, isPeriodic);
    
    Vector<Box> boxes;
    domainSplit(domainBox, boxes, fixedBoxSize);
    Vector<int> procAssign;
    int eekflag = LoadBalance(procAssign, boxes);
    assert(eekflag == 0);
    
    DisjointBoxLayout grids(boxes, procAssign, domainBox); 
    DataIterator dit(grids);
    
    ParticleData<BinItem> particles(grids,
				    probDomain,
				    fixedBoxSize,
				    meshSpacing, 
				    origin);

    // put the ghost particles in a new ParticleData, with ghost "cells"  
    BoxLayout ghostGrids;
    ghostGrids.deepCopy(grids);
    ghostGrids.grow(numGhost);
    ghostGrids.closeNoSort();

    ParticleData<BinItem> particlesWithGhosts(ghostGrids, 
					      probDomain,
					      fixedBoxSize,
					      meshSpacing,
					      origin);

    initData(particles, dx, domainDimension);

    particles.fillGhosts(particlesWithGhosts, numGhost);
    
    assert(particlesWithGhosts.isClosed());

    for (DataIterator dit(ghostGrids); dit.ok(); ++dit)
      {
	unsigned idx = grids.index(dit());
	size_t numGhosts = particlesWithGhosts[dit].numItems();
	if (verbose)
	  {
	    pout() << "Box " << idx << " has " << numGhosts;
	    pout() << " particles." << endl;
	  }
	if (numGhosts != partPerBox)
	  {
	    if (verbose)
	      {
		pout() << "Box has wrong number of particles after ghost fill." << endl;
	      }
	    status += 1;
	  }
      }
  }

  // done
  pout() << indent << pgmname << ": "
	 << ( (status == 0) ? "passed all tests" : "failed at least one test,")
	 << endl;

  CH_TIMER_REPORT();

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
initData(ParticleData<BinItem>& a_data, 
         const Real a_dx, 
         const int a_numParticles)
{

  CH_TIMERS("testGhostParticles(initData)");

  int N_hi = a_numParticles - 1;
  IntVect lo(D_DECL(0, 0, 0));
  IntVect hi(D_DECL(N_hi, N_hi, N_hi));
  Box domain(lo, hi);

  BoxLayout BL = a_data.getBoxes();
  DataIterator dit(BL);
  for (dit.begin(); dit.ok(); ++dit)
  {
    const Box& thisBox = BL[dit];
    CH_XD::List<BinItem> thisList;
    BoxIterator bit(thisBox);
    for (bit.begin(); bit.ok(); ++bit)
      {
	IntVect iv = bit();
	RealVect position = ((RealVect)iv + 0.5)*a_dx;
	BinItem particle(position);
	thisList.append(particle);
      }
    a_data[dit].addItemsDestructive(thisList);
  }
}
