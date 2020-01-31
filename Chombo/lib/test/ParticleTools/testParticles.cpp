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
//  Test the infrastucture for particles. This tests the Particle 
// class itself, making sure the linearization, getting / setting,
// and equality operations work as expected.
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

#include "parstream.H"
#include "Particle.H"
#include "RealVect.H"

#ifdef CH_MPI
#include "mpi.h"
#endif

#include "UsingNamespace.H"

void parseTestOptions(int argc, char* argv[]);


static const char *pgmname = "testParticles";
static const char *indent = "   ";

static bool verbose = true ;

//////////////////////////////////////////////////////////////
int
main(int argc ,char *argv[] )
{
  int status = 0;

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  parseTestOptions(argc, argv);

  // define two particles with different values
  Real mass1 = 1.0;
  RealVect position1 = RealVect(D_DECL(1.0, 1.0, 1.0));
  RealVect velocity1 = RealVect(D_DECL(0.5, 0.5, 0.5));

  Particle particle1;
  particle1.define(mass1, position1, velocity1);

  Real mass2 = 2.0;
  RealVect position2 = RealVect(D_DECL(-1.0, -1.0, -1.0));
  RealVect velocity2 = RealVect(D_DECL(-0.5, -0.5, -0.5));

  Particle particle2;
  particle2.define(mass2, position2, velocity2);

  // test that particle1 has the correct values. This tests that 
  // the assignment and access functions are working correctly. 
  if (mass1 != particle1.mass())
    {
      if (verbose) 
	{
	  pout() << " fail: particle1 mass is not correct" << std::endl;
	}
      status = 1;
    }
  
  if (position1 != particle1.position())
    {
      if (verbose) 
	{
	  pout() << " fail: particle1 position is not correct" << std::endl;
	}
      status = 1;
    }

  if (velocity1 != particle1.velocity())
    {
      if (verbose) 
	{
	  pout() << " fail: particle1 velocity is not correct" << std::endl;
	}
      status = 1;
    }

  // acceleration should have been set to zero by default
  if (RealVect(D_DECL(0.0, 0.0, 0.0)) != particle1.acceleration())
    {
      if (verbose) 
	{
	  pout() << " fail: particle1 acceleration is not correct" << std::endl;
	}
      status = 1;
    }

  // test the equality operators
  if (particle1 == particle2)
    {
      if (verbose) 
	{
	  pout() << " fail: particle equality operator is giving the wrong answer" << std::endl;
	}
      status = 1;
    }

  if (&particle1 == &particle2)
    {
      if (verbose) 
	{
	  pout() << " fail: particle equality operator is giving the wrong answer" << std::endl;
	}
      status = 1;
    }

  // now test linearIn/Out
  static const int NB=100;
  char *buf = new char[NB];

  for ( int i=0 ; i<NB ; ++i )
    {
      buf[i] = '@' ;
    }

  
  particle1.linearOut(buf);

  particle2.linearIn(buf);

  // particle1 and particle2 should now have equal values
  if (particle1 != particle2)
    {
      if (verbose) 
	{
	  pout() << " fail: particle linearIn / linearOut operations are not working properly" << std::endl;
	}
      status = 1;
    }

  // but the pointers should still be different
  if (&particle1 == &particle2)
    {
      if (verbose) 
	{
	  pout() << " fail: particle linearIn / linearOut operations are not working properly" << std::endl;
	}
      status = 1;
    }

  // test size
  if (particle1.size() != particle2.size())
    {
      if (verbose) 
	{
	  pout() << " fail: particle size is not correct" << std::endl;
	}
      status = 1;
    }

  // size of the basic particle be should 1 Real plus 3 RealVects
  if (particle1.size() != (1 + 3*SpaceDim) * sizeof(Real))
    {
      if (verbose) 
	{
	  pout() << " fail: particle size is not correct" << std::endl;
	}
      status = 1;
    }

  // done
  pout() << indent << pgmname << ": "
	 << ( (status == 0) ? "passed all tests" : "failed at least one test,")
	 << std::endl;

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
