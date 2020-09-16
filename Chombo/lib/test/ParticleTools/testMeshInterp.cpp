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
//  Test the particle-mesh interaction functions. This creates a Particle
//  with a funky position, and deposits it onto the mesh using NGP, CIC, TSC, 
//  and W4 interpolants. It also tests interpolating the field back to the 
//  particle position using the same kernels. In each case, the results is
//  compared to an exact answer.
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
#include "MeshInterp.H"
#include "RealVect.H"
#include "List.H"
#include "BoxIterator.H"

#ifdef CH_MPI
#include "mpi.h"
#endif

#include "UsingNamespace.H"

static const char *pgmname = "testMeshInterp";
static const char *indent = "   ";

//////////////////////////////////////////////////////////////
int
main(int argc ,char *argv[] )
{
  int status = 0;

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  // basic properties of the mesh
  int N = 16;
  Real domainLength = 1.0;
  Real dx = domainLength / N;
  RealVect leftEdge(D_DECL6(0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
  RealVect meshSpacing(D_DECL6(dx, dx, dx, dx, dx, dx));
  IntVect smallEnd(D_DECL6(0, 0, 0, 0, 0, 0));
  IntVect bigEnd(D_DECL6(N, N, N, N, N, N));

  // create our MeshInterp object
  Box domainBox(smallEnd, bigEnd);
  MeshInterp mesh(domainBox, meshSpacing, leftEdge);

  // The FArrayBox that the particle will be deposited on
  FArrayBox rho(domainBox, 1);
  rho.setVal(0.0);

  // Properties of the particle. We give it a somewhat 
  // irregular position - note that if you change this
  // you will also need to change the expected answer below.
  Real mass = 1.0;
  Real x = 4.75*dx;
  RealVect position = RealVect(D_DECL6(x, x, x, x, x, x));

  // Define our Particle
  Particle p;
  p.define(mass, position);

  // make a list with just the one particle
  List<Particle> particleList;
  particleList.append(p);

  // First, test the deposition functions. The correct answers for depositing
  // the above partice are hard-coded below.

  // Deposit CIC
  InterpType interp = CIC;
  mesh.deposit(particleList, rho, interp);

  // Check answer
  BoxIterator bit(domainBox);
  for (bit.reset(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      
      Real cellVolume = pow(dx, SpaceDim);
      Real expectedAnswer = 1.0 / cellVolume;

      for (int idir = 0; idir < SpaceDim; idir++)
        {

          
          if (iv[idir] == 4)
            {
              expectedAnswer *= 0.75;
            }

          else if (iv[idir] == 5)
            {
              expectedAnswer *= 0.25;
            }
            
          else
            {
              expectedAnswer *= 0.0;
            }
        }
         
      if (rho.get(iv, 0) != expectedAnswer)
        {
          pout() << " fail: CIC deposited density is not correct at point " << iv << std::endl;
          pout() << indent << " got: " << rho.get(iv, 0) 
                 << " expected: " << expectedAnswer << std::endl;
          status = 1;         
        }

    }

  // reset rho to zero
  rho *= 0.0;

  // Deposit
  interp = NGP;
  mesh.deposit(particleList, rho, interp);

  // Check answer
  for (bit.reset(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      
      Real cellVolume = pow(dx, SpaceDim);
      Real expectedAnswer = 1.0 / cellVolume;

      for (int idir = 0; idir < SpaceDim; idir++)
        {

          
          if (iv[idir] != 4)
            {
              expectedAnswer *= 0.0;
            }

        }
         
      if (rho.get(iv, 0) != expectedAnswer)
        {
          pout() << " fail: NGP deposited density is not correct at point " << iv << std::endl; 
  pout() << indent << " got: " << rho.get(iv, 0) 
                 << " expected: " << expectedAnswer << std::endl;
          status = 1;         
        }

    }

  // reset rho to zero
  rho *= 0.0;

  // test TSC deposit
  interp = TSC;
  mesh.deposit(particleList, rho, interp);

  // Check answer
  for (bit.reset(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      
      Real cellVolume = pow(dx, SpaceDim);
      Real expectedAnswer = 1.0 / cellVolume;

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          if (iv[idir] == 3)
            {
              expectedAnswer *= 1.0 / 32.0;
            }

          else if (iv[idir] == 4)
            {
              expectedAnswer *= 22.0 / 32.0;
            }

          else if (iv[idir] == 5)
            {
              expectedAnswer *= 9.0 / 32.0;
            }

          else
            {
              expectedAnswer *= 0.0;
            }

        }
         
      if (rho.get(iv, 0) != expectedAnswer)
        {
          pout() << " fail: TSC deposited density is not correct at point " << iv << std::endl;
          pout() << indent << " got: " << rho.get(iv, 0) 
                 << " expected: " << expectedAnswer << std::endl;
          status = 1;         
        }

    }

  // reset rho to zero
  rho *= 0.0;

  // test W4 deposit
  interp = W4;
  mesh.deposit(particleList, rho, interp);

  // Check answer
  for (bit.reset(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      
      Real cellVolume = pow(dx, SpaceDim);
      Real expectedAnswer = 1.0 / cellVolume;

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          if (iv[idir] == 3)
            {
              expectedAnswer *= (-9.0/128.0);
            }

          else if (iv[idir] == 4)
            {
              expectedAnswer *= (111.0/128.0);
            }

          else if (iv[idir] == 5)
            {
              expectedAnswer *= (29.0/128.0);
            }

          else if (iv[idir] == 6)
            {
              expectedAnswer *= (-3.0/128.0);
            }
          
          else
            {
              expectedAnswer *= 0.0;
            }

        }
         
      if (rho.get(iv, 0) != expectedAnswer)
        {
          pout() << " fail: W4 deposited density is not correct at point " << iv << std::endl;
          pout() << indent << " got: " << rho.get(iv, 0) 
                 << " expected: " << expectedAnswer << std::endl;
          status = 1;         
        }

    }

  // Now, we test the interpolation functions. We need a spacedim-component field for this.
  FArrayBox field(domainBox, SpaceDim);
  for (bit.reset(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          field.set(iv, idir, iv[idir]);
        }
    }

  interp = NGP;
  mesh.interpolate(particleList, field, interp);
  for (ListIterator<Particle> lit(particleList); lit; ++lit)
    {
      Particle& this_particle = lit();
      RealVect& this_acceleration = this_particle.acceleration();
      RealVect expectedField = RealVect(D_DECL6(4.0, 4.0, 4.0, 4.0, 4.0, 4.0));
      if (this_acceleration != expectedField)
        {
          pout() << " fail: NGP interpolated field is not correct" << std::endl;
          pout() << indent << " got: " << this_acceleration 
                 << " expected: " << expectedField << std::endl;
          status = 1;         
        }
      this_particle.setAcceleration(RealVect(D_DECL6(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)));
    }

  interp = CIC;
  mesh.interpolate(particleList, field, interp);
  for (ListIterator<Particle> lit(particleList); lit; ++lit)
    {
      Particle& this_particle = lit();
      RealVect& this_acceleration = this_particle.acceleration();
      RealVect expectedField = RealVect(D_DECL6(4.25, 4.25, 4.25, 4.25, 4.25, 4.25));
      if (this_acceleration != expectedField)
        {
          pout() << " fail: NGP interpolated field is not correct" << std::endl;
          pout() << indent << " got: " << this_acceleration 
                 << " expected: " << expectedField << std::endl;
          status = 1;         
        }
      this_particle.setAcceleration(RealVect(D_DECL6(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)));
    }


  interp = TSC;
  mesh.interpolate(particleList, field, interp);
  for (ListIterator<Particle> lit(particleList); lit; ++lit)
    {
      Particle& this_particle = lit();
      RealVect& this_acceleration = this_particle.acceleration();
      RealVect expectedField = RealVect(D_DECL6(4.25, 4.25, 4.25, 4.25, 4.25, 4.25));
      if (this_acceleration != expectedField)
        {
          pout() << " fail: NGP interpolated field is not correct" << std::endl;
          pout() << indent << " got: " << this_acceleration 
                 << " expected: " << expectedField << std::endl;
          status = 1;         
        }
      this_particle.setAcceleration(RealVect(D_DECL6(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)));
    }  

  interp = W4;
  mesh.interpolate(particleList, field, interp);
  for (ListIterator<Particle> lit(particleList); lit; ++lit)
    {
      Particle& this_particle = lit();
      RealVect& this_acceleration = this_particle.acceleration();
      RealVect expectedField = RealVect(D_DECL6(4.25, 4.25, 4.25, 4.25, 4.25, 4.25));
      if (this_acceleration != expectedField)
        {
          pout() << " fail: NGP interpolated field is not correct" << std::endl;
          pout() << indent << " got: " << this_acceleration 
                 << " expected: " << expectedField << std::endl;
          status = 1;         
        }
      this_particle.setAcceleration(RealVect(D_DECL6(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)));
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
