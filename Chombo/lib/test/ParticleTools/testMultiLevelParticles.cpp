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
//  Test the infrastucture for particles. This test creates two AMR 
//  levels worth of Particles. It then moves them using a constant velocity
//  field in a way that particles will move across level boundaries. This
//  tests the functions that pass around particles between refinement levels. 
//   
// Usage:
//  <program-name> [-q|-v] ...
//
//  where:
//    -q means run quietly (only pass/fail messages printed)
//    -v means run verbosely (all messages printed)
//    ... all non-option arguments are ignored (Chombo convention)
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
#include "ParticleValidRegion.H"

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

void advance(ParticleData<BinItem>& a_data, const Real a_delta);

/// Global variables for handling output

static const char *pgmname = "testMultiLevelParticles";
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

    // do not change these for this test or the numbers below will not make sense.
    int maxBoxSize = 8;
    int domainDimension = 32;
    
    int numParticles = pow((double) domainDimension, SpaceDim);

    { 

      RealVect origin(D_DECL(0,0,0));

      // create the coarse level
      Real dx_crse = 1.0 / domainDimension;
      RealVect crseMeshSpacing(D_DECL(dx_crse, dx_crse, dx_crse));

      Box crseDomainBox(IntVect::Zero, (domainDimension - 1) * IntVect::Unit);
      ProblemDomain crseDomain(crseDomainBox);

      Vector<Box> crseBoxes;
      domainSplit(crseDomainBox, crseBoxes, maxBoxSize);
      Vector<int> procAssign;
      int eekflag = LoadBalance(procAssign, crseBoxes);
      assert(eekflag == 0);

      DisjointBoxLayout crseGrids(crseBoxes, procAssign, crseDomain); 

      ParticleData<BinItem> crseParticles(crseGrids, 
					  crseDomain,
					  maxBoxSize,
					  crseMeshSpacing,
					  origin);

      ParticleValidRegion crseValidRegion(crseGrids, NULL, 1, 0);

      // now create the fine level
      int RefRatio = 2;
      int n_buffer_cells = 2;
      Real dx_fine = dx_crse / 2.0;
      RealVect fineMeshSpacing(D_DECL(dx_fine, dx_fine, dx_fine));

      Box fineDomainBox = refine(crseDomainBox, RefRatio);
      ProblemDomain fineDomain = refine(crseDomain, RefRatio);

      int fineBoxDimension = 32;
      Box fineLevelBox(IntVect(D_DECL(16, 16, 16)), IntVect(D_DECL(47, 47, 47)));
      Vector<Box> fineBoxes;
      domainSplit(fineLevelBox, fineBoxes, maxBoxSize);
      eekflag = LoadBalance(procAssign, fineBoxes);
      assert(eekflag == 0);

      DisjointBoxLayout fineGrids(fineBoxes, procAssign, fineDomain);

      ParticleData<BinItem> fineParticles(fineGrids,
					  fineDomain,
					  maxBoxSize,
					  fineMeshSpacing,
					  origin);

      ParticleValidRegion fineValidRegion(fineGrids, crseValidRegion.mask(), RefRatio, n_buffer_cells);

      if (verbose)
      {
        pout() << "Filling initial conditions." << endl;
        pout() << "    Number of particles: ";
      }

      initData(crseParticles, dx_crse, domainDimension);

      if (verbose)
      {      
        pout() << crseParticles.numValid();
        pout() << ". Expected: " << numParticles << "." << endl;
      }

      if (crseParticles.numValid() != numParticles) 
      {
        ++status;
        if (verbose)
        {
          pout() << "    Fail. Lost particles during intialization." << endl;
        }
      }

    if (verbose)
    {
      pout() << "Transfering particles from crse to fine. " << endl; 
    }

    // now find all the particles that should go on the fine level
    collectValidParticles(fineParticles.outcast(), 
                          crseParticles, 
                          fineValidRegion.mask(),
                          fineMeshSpacing, 
                          RefRatio);

    // Here we compute the expected answers. "numBuffer" is the number of 
    // particles that lie on the "buffer region" at coarse-fine level boundaries, 
    // and are thus sent to the crse instead of the fine level.
    int numBuffer = pow((double) fineBoxDimension/RefRatio, SpaceDim) 
      - pow((double) (fineBoxDimension/RefRatio) - n_buffer_cells, SpaceDim);
    int crseExpected = pow((double) domainDimension, SpaceDim) // the number of particles in the entire domain
      - pow((double) (fineBoxDimension/RefRatio), SpaceDim) // the number that lie on the region covered by the fine level
      + numBuffer; // the number that live in the buffer region.
    int fineExpected = pow((double) fineBoxDimension/RefRatio, SpaceDim) // the number that live on the fine level
      - numBuffer; // the number that live in the buffer region.

    int crseActual = crseParticles.numValid();

    if (verbose)
    {
      pout() << "    Number of particles on coarse level after transfer: "; 
      pout() << crseActual << ". ";
      pout() << "Expected: " << crseExpected << endl;
    }

    fineParticles.remapOutcast();

    int fineActual = fineParticles.numValid();

    if (verbose)
    {
      pout() << "    Number of particles on fine level after transfer: "; 
      pout() << fineActual << ". ";
      pout() << "Expected: " << fineExpected << endl;
    }   

    if (verbose)
    {
      pout() << "    Number of outcasts on the fine level after transfer: "; 
      pout() << fineParticles.numOutcast() << ". ";
      pout() << "Expected: " << 0 << endl;
    }       

    if (verbose)
    {
      pout() << "    Number lost: "; 
      pout() << numParticles - crseParticles.numValid() - fineParticles.numValid() << ". ";
      pout() << "Expected: " << 0 << endl;
    }       

    if (crseExpected != crseActual || fineExpected != fineActual)
    {
      status++;
      if (verbose)
      {
        pout() << "    Fail. Lost particles during intial transfer." << endl;
      }
    }

    if (verbose)
    {
      pout() << "Moving particles. " << endl; 
    }

    // now we translate all the particles. Particles on the coarse level will be
    // shifted by one dx_crse.
    advance(crseParticles, dx_crse);

    // Particles on the fine level will be moved *twice* by a distance dx_fine each time.
    advance(fineParticles, dx_fine);
    advance(fineParticles, dx_fine);    

    if (verbose)
    {
      pout() << "Transfering particles to correct level after move. " << endl; 
    }    

    // The two levels have now advanced to the same time. 
    // Some of the particles have changed levels. First, collect all the particles that
    // go from crse to fine:
    collectValidParticles(fineParticles.outcast(), 
                          crseParticles, 
                          fineValidRegion.mask(),
                          fineMeshSpacing, 
                          RefRatio);

    int expectedNumTransfer = pow((double) fineBoxDimension/RefRatio - n_buffer_cells, SpaceDim-1);

    if (verbose)
    {
      pout() << "    Transfering " << fineParticles.numOutcast() << " particles from crse to fine. ";
      pout() << "Expected: " << expectedNumTransfer << endl;
    }
    
    // now take the particles that have left the fine level and put them on the crse
    collectValidParticles(crseParticles.outcast(), fineParticles, 
                          fineValidRegion.mask(), fineMeshSpacing, 1, true);

    if (verbose)
    {
      pout() << "    Transfering " << crseParticles.numOutcast() << " particles from fine to crse. ";
      pout() << "Expected: " << expectedNumTransfer << endl;
    }

    if ( fineParticles.numOutcast() != expectedNumTransfer || 
	 crseParticles.numOutcast() != expectedNumTransfer)
    {
      status++;
      if (verbose)
        {
          pout() << "    Fail. Wrong number of particles transfered after move" << endl;
        }
    }

    // finally, remap outcasts again for both the fine and the crse level
    fineParticles.remapOutcast();
    crseParticles.remapOutcast();

    if (verbose)
    {
      pout() << "Checking for particle conservation. " << endl; 
    }  

    crseActual = crseParticles.numValid();
    fineActual = fineParticles.numValid();

    if (verbose)
    {
      pout() << "    Number of particles on coarse level after transfer: "; 
      pout() << crseActual << ". ";
      pout() << "Expected: " << crseExpected << endl;
    }

    if (verbose)
    {
      pout() << "    Number of particles on fine level after transfer: "; 
      pout() << fineActual << ". ";
      pout() << "Expected: " << fineExpected << endl;
    }   

    if (verbose)
    {
      pout() << "    Number of outcasts on the crse level after transfer: "; 
      pout() << crseParticles.numOutcast() << ". ";
      pout() << "Expected: " << 0 << endl;
    }       

    if (verbose)
    {
      pout() << "    Number of outcasts on the fine level after transfer: "; 
      pout() << fineParticles.numOutcast() << ". ";
      pout() << "Expected: " << 0 << endl;
    }       

    if (verbose)
    {
      pout() << "    Number lost: "; 
      pout() << numParticles - crseParticles.numValid() - fineParticles.numValid() << ". ";
      pout() << "Expected: " << 0 << endl;
    }       

    if (crseExpected != crseActual || fineExpected != fineActual)
    {
      status++;
      if (verbose)
      {
        pout() << "    Fail. Lost particles after transfer." << endl;
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

// advance moves the particles, collects together those that have
// changed boxes, and puts these "outcasts" back in the right places.
void advance(ParticleData<BinItem>& a_data, const Real a_delta)
{
  
  // move the particles by delta
  moveParticles(a_data, a_delta);

  // collect outcast particles
  a_data.gatherOutcast();

  // rebin
  a_data.remapOutcast();

}
