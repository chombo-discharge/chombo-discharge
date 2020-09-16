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
//  Test the IO routines for particles. This test creates a ParticleData 
//  and writes it out to an hdf5 file. It then reads the particles back
//  in and verifies that the correct results are obtained. This process
//  is repeated for each particle type defined in ParticleTools.
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

#include "ParticleIO.H"
#include "CH_HDF5.H"
#include "parstream.H"
#include "ListBox.H"
#include "ListBoxFactory.H"
#include "BinItem.H"
#include "Particle.H"
#include "FArrayBox.H"
#include "RealVect.H"
#include "BRMeshRefine.H"
#include "AMRIO.H"
#include "LoadBalance.H"
#include "ParticleData.H"

#ifdef CH_MPI
#include "mpi.h"
#endif

#include "UsingNamespace.H"

//////////////////////////////////////////////////////////////
using std::endl;

void parseTestOptions(int argc, char* argv[]);

// overload initData functions for each particle type
void initData(ParticleData<BinItem>& a_data, 
              const Real a_dx,
              const int a_numParticles);

void initData(ParticleData<Particle>& a_data, 
              const Real a_dx,
              const int a_numParticles);

template <class P>
int testParticleIO();

/// Global variables for handling output
static const char *pgmname = "testParticleIO";
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

    status += testParticleIO<BinItem>();
    status += testParticleIO<Particle>();
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

template <class P>
int testParticleIO()
{
  int status = 0;
#ifdef CH_USE_HDF5
  int maxBoxSize = 32;
  int domainDimension = 64;
    
  Real dx = 1.0 / domainDimension;
  RealVect meshSpacing(D_DECL(dx,dx,dx));
  RealVect origin(D_DECL(0,0,0));

  Box domainBox(IntVect::Zero, (domainDimension - 1) * IntVect::Unit);
  ProblemDomain probDomain(domainBox);

  Vector<Box> boxes;
  domainSplit(domainBox, boxes, maxBoxSize);
  Vector<int> procAssign;
  int eekflag = LoadBalance(procAssign, boxes);
  assert(eekflag == 0);
  
  DisjointBoxLayout grids(boxes, procAssign, domainBox); 
  DataIterator dit(grids);
  
  ListBoxFactory<P> myFactory(meshSpacing, origin);
  
  ParticleData<P> particlesBefore(grids, 
                                  probDomain,
                                  maxBoxSize,
                                  meshSpacing,
                                  origin);

  ParticleData<P> particlesAfter(grids, 
                                  probDomain,
                                  maxBoxSize,
                                  meshSpacing,
                                  origin);
  
  initData(particlesBefore, dx, domainDimension);
  
  int numParticlesBefore;
  numParticlesBefore = particlesBefore.numParticles();
  if (verbose)
    {
      pout() << "Saving " << numParticlesBefore << " particles." << endl;
    }
  
  // test output filename
  string filename = "testParticleIO.hdf5";    
  
  // create file, write particles, close
  HDF5Handle outputHandle(filename, HDF5Handle::CREATE);
  writeParticlesToHDF(outputHandle, particlesBefore, "particles");
  outputHandle.close();
  
  // read particles
  HDF5Handle inputHandle(filename, HDF5Handle::OPEN_RDONLY);
  readParticlesFromHDF(inputHandle, particlesAfter, "particles");
  inputHandle.close();
  
  int numParticlesAfter;
  numParticlesAfter = particlesAfter.numParticles();
  
  if (verbose)
    {
      pout() << "Read in " << numParticlesAfter << " particles." << endl;
    }
  
  if(numParticlesBefore != numParticlesAfter)
    {
      ++status;
      if (verbose)
        {
          pout() << "Fail. Lost particles during input / output." << endl;
        }
    }

  for (dit.reset(); dit.ok(); ++dit)
    {
      List<P> before = particlesBefore[dit()].listItems();
      List<P> after  = particlesAfter[dit()].listItems();
      
      if (before.length() != after.length())
        {
          ++status;
          if (verbose)
            {
              pout() << "Fail. Box has lost particles during I/O." << endl;
            }
        }
      
      else
        {
          ListIterator<P> litBefore(before);
          ListIterator<P> litAfter(after);
          for (litBefore.begin(); litBefore; ++litBefore)
            {
              P& p1 = litBefore();
              P& p2 = litAfter();
              if (p1 != p2)
                {
                  ++status;
                  if (verbose)
                    {
                      pout() << "Fail. Particle data has changed during IO." << endl;
                    }
                }
              ++litAfter;
            }
        }
      
      
    }

#endif
  return status;
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

//----------------------------------------------------------------

void
initData(ParticleData<Particle>& a_data, 
         const Real a_dx, 
         const int a_numParticles)
{

  int N_hi = a_numParticles - 1;
  IntVect lo(D_DECL(0, 0, 0));
  IntVect hi(D_DECL(N_hi, N_hi, N_hi));
  Box domain(lo, hi);

  CH_XD::List<Particle> thisList;
  BoxIterator bit(domain);
  for (bit.begin(); bit.ok(); ++bit)
  {
    IntVect iv = bit();
    RealVect position = ((RealVect)iv + 0.5)*a_dx;
    Real mass = 1.0;
    Particle particle(mass, position);
    thisList.append(particle);
  }

  BoxLayout BL = a_data.getBoxes();
  DataIterator dit(BL);
  for (dit.begin(); dit.ok(); ++dit)
  {
    a_data[dit].addItemsDestructive(thisList);
  }

}
