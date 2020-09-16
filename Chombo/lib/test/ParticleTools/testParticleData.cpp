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
//  Test the infrastucture for particles. This test creates a level's 
//  worth of ListBoxes and populates them with particles. It then shifts
//  the particles so that some of them move out of their original boxes.
//  We then call remap the particles to their new boxes, moving them to 
//  proper processor if needed. 
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

void moveParticles(ParticleData<BinItem>& a_particles, 
                   const Real a_dx);

/// Global variables for handling output

static const char *pgmname = "testParticleData";
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

    int maxBoxSize = 32;
    int domainDimension = 64;

    int numParticles = D_TERM(domainDimension,*domainDimension,*domainDimension);
    int numOutCast;
    if (domainDimension == maxBoxSize)
      {
        numOutCast = 0;
      }
    else
      {
        numOutCast = (domainDimension / maxBoxSize) * pow( (double) domainDimension, SpaceDim-1);
      }

    Real dx = 1.0 / domainDimension;
    RealVect meshSpacing(D_DECL(dx,dx,dx));
    RealVect origin(D_DECL(0,0,0));

    Box domainBox(IntVect::Zero, (domainDimension - 1) * IntVect::Unit);
    ProblemDomain probDomain(domainBox);

    { 
      Vector<Box> boxes;
      domainSplit(domainBox, boxes, maxBoxSize);
      Vector<int> procAssign;
      int eekflag = LoadBalance(procAssign, boxes);
      assert(eekflag == 0);

      DisjointBoxLayout grids(boxes, procAssign, domainBox); 
      DataIterator dit(grids);

      ParticleData<BinItem> levelParticles;
      levelParticles.define(grids,
                            probDomain,
                            maxBoxSize,
                            meshSpacing,
                            origin);

      if (verbose)
      {
        pout() << "Filling initial conditions. Number of particles: ";
      }

      initData(levelParticles, dx, domainDimension);

      if (verbose)
      {      
        pout() << levelParticles.numValid();
        pout() << ". Expected: " << numParticles << "." << endl;
      }

      if (levelParticles.numValid() != numParticles) 
      {
        ++status;
        if (verbose)
        {
          pout() << "Fail. Lost particles during intialization." << endl;
        }
      }

      if (verbose)
      {
        pout() << "Initial number of misplaced particles: ";
      }

      levelParticles.gatherOutcast();
     
      if (verbose)
      {
        pout() << levelParticles.numOutcast();
        pout() << ". Expected: " << 0 << "." << endl;
      }

      if (levelParticles.numOutcast() != 0)
      {
        ++status;
        if (verbose)
        {
          pout() << "Fail. There should be no outcast particles immediately after initialization." << endl;
        }
      }

      if (verbose)
      {
        pout() << "Begin moving particles." << endl;
      }

      for (int i = 0; i < 2; i++)
      {
        moveParticles(levelParticles, dx);

      if (verbose)
      {
        pout() << "Number of outcasts: ";
      }

      levelParticles.gatherOutcast();

      if (verbose)
      {
        pout() << levelParticles.numOutcast();
        pout() << ". Expected: " << numOutCast << "." << endl;
      }

        if (levelParticles.numOutcast() != numOutCast) 
        {
          ++status;
          if (verbose)
          {
            pout() << "Fail. Wrong number of outcast particles after move." << endl;
          }
        }

      if (verbose)
      {
        pout() << "Number that are still in the right place: "; 
        pout() << levelParticles.numValid();
        pout() << ". Expected: " << numParticles - numOutCast << "." << endl;
      }

        if (levelParticles.numValid() != numParticles - numOutCast) 
        {
          ++status;
          if (verbose)
          {
            pout() << "Fail. Wrong number of valid particles after move." << endl;
          }
        }

        levelParticles.remapOutcast();

      if (verbose)
      {
        pout() << "Number of outcasts after rebin: ";
      }

      levelParticles.gatherOutcast();

      if (verbose)
      {
        pout() << levelParticles.numOutcast();
        pout() << ". Expected: " << 0 << "." << endl;
      }

        if (levelParticles.numOutcast() != 0) 
        {
          ++status;
          if (verbose)
          {
            pout() << "Fail. Wrong number of outcast particles after rebin." << endl;
          }
        }

      if (verbose)
      {
        pout() << "Number of correctly sorted particles after rebin: "; 
        pout() << levelParticles.numValid();
        pout() << ". Expected: " << numParticles << "." << endl;
      }

        if (levelParticles.numValid() != numParticles) 
        {
          ++status;
          if (verbose)
          {
            pout() << "Fail. Wrong number of valid particles after rebin." << endl;
          }
        }
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
initData(ParticleData<BinItem>& a_data, 
         const Real a_dx, 
         const int a_numParticles)
{

  int N_hi = a_numParticles - 1;
  IntVect lo(D_DECL(0, 0, 0));
  IntVect hi(D_DECL(N_hi, N_hi, N_hi));
  Box domain(lo, hi);

  CH_XD::List<BinItem> thisList;
  BoxIterator bit(domain);
  for (bit.begin(); bit.ok(); ++bit)
  {
    IntVect iv = bit();
    RealVect position = ((RealVect)iv + 0.5)*a_dx;
    BinItem particle(position);
    thisList.append(particle);
  }

  BoxLayout BL = a_data.getBoxes();
  DataIterator dit(BL);
  for (dit.begin(); dit.ok(); ++dit)
  {
    a_data[dit].addItemsDestructive(thisList);
  }

}

void moveParticles(ParticleData<BinItem>& a_data, const Real dx)
{

  RealVect shift(D_DECL(1.0, 0.0, 0.0));
  shift *= dx;

  BoxLayout BL = a_data.getBoxes();
  DataIterator dit(BL);
  for (dit.begin(); dit.ok(); ++dit)
  {
    for (ListIterator<BinItem> li (a_data[dit].listItems()); li.ok(); ++li)
    {
      RealVect x=li().position();
      RealVect new_position = x + shift;

      // enforce periodic boundary conditions
      for (int d = 0; d < CH_SPACEDIM; d++)
      {
        if (new_position[d] < 0.0)
        {
          new_position[d] += 1.0;
        }

        if (new_position[d]  > 1.0)
        {
          new_position[d] -= 1.0;
        }
      }

      li().setPosition(new_position);
    }
  }
}
