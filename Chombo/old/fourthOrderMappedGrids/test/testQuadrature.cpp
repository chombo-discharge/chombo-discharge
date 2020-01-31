#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>
using std::cout;
using std::endl;

#include "parstream.H"
#include "GaussianQuadrature.H"
#include "NewtonCotesQuadrature.H"

#ifdef CH_MPI
#include <mpi.h>
#endif

/// Global variables for handling output:
static const char* pgmname = "testQuadrature" ;
static const char* indent = "   ";
static const char* indent2 = "      " ;
static bool verbose = true ;

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
  return ;
}

int
testQuadrature();

int
main(int argc ,char* argv[])
{
#ifdef CH_MPI
  MPI_Init (&argc, &argv);
#endif
  parseTestOptions( argc ,argv ) ;
  if ( verbose )
    pout () << indent2 << "Beginning " << pgmname << " ..." << endl ;

  int status = testQuadrature();

  if ( status == 0 )
    cout << indent << "Quadrature test" << " passed." << endl ;
  else
    cout << indent << "Quadrature test" << " failed with return code " << status << endl ;



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
testQuadrature()
{
  int returnCode = 0;

  int numPointsMin = 2;
  int numPointsMax = 5;
  for (int numPoints = numPointsMin; numPoints<=numPointsMax; numPoints++)
    {
      GaussianQuadrature testGaussianQuad(numPoints);
      NewtonCotesQuadrature testNCQuad(numPoints);

      Real avgRateGaussian = 0.0;
      Real avgRateNC = 0.0;
      int numRates = 0;

      pout() << "Testing " << numPoints << "-point Gaussian Quadrature" << endl;
      int dir1 = SpaceDim-1;
      int dir2 = SpaceDim-2;

      int nx = 4;
      //Real length = 2.0;
      Real length = 1.0;
      int numRefine = 5;
      Real oldGaussianErr = 0.0;
      Real oldNCErr = 0.0;
      Real Pi = 3.14159;

      for (int ref=0; ref<numRefine; ref++)
        {
          // integrate from 0->length
          Real dx = length/nx;

          {
            // Gaussian quadrature

            const Vector<QuadratureElement> coeffs= testGaussianQuad.coefficients(dir1, dir2);

            Real integral = 0;
            RealVect cellLength(dx*RealVect::Unit);
            Real scale = testGaussianQuad.weightMult(cellLength, dir1, dir2);


            // start at center of lowest "cell"
            Real xCenter = 0.5*dx;
            for (int i=0; i<nx; i++)
              {

                // do quadrature at this point
                for (int i=0; i<coeffs.size(); i++)
                  {
                    Real x = xCenter + coeffs[i].location[0]*dx/2;
                    //Real f = x*x*x*x;
                    Real f = sin(Pi*x);
                    integral += coeffs[i].weight*f*scale;
                  }
                // index to next point
                xCenter += dx;
              } // end loop over "cells"

            //Real exactVal = pow(length, 5)/5.0;
            Real exactVal =  (1 - cos(Pi*length))/Pi;
            if (SpaceDim < 3) exactVal = 0.0;
            Real error = exactVal - integral;
            pout() << "nPts = " << nx
                   << ",  Gaussian computed = " << integral
                   << ", exact value = " << exactVal
                   << ", error = " << error;

            if (ref > 0)
              {
                Real rate = log(Abs(oldGaussianErr/error))/log(2.0);
                pout () << "   Rate = " << rate;
                avgRateGaussian += rate;
                numRates++;
              }

            pout () << endl;
            oldGaussianErr = error;
          }


          {
            // Newton-Cotes quadrature
            const Vector<QuadratureElement> coeffs= testNCQuad.coefficients(dir1, dir2);

            Real integral = 0;
            RealVect cellLength(dx*RealVect::Unit);
            Real scale = testNCQuad.weightMult(cellLength, dir1, dir2);


            // start at center of lowest "cell"
            Real xCenter = 0.5*dx;
            for (int i=0; i<nx; i++)
              {

                // do quadrature at this point
                for (int i=0; i<coeffs.size(); i++)
                  {
                    Real x = xCenter + coeffs[i].location[0]*dx/2;
                    //Real f = x*x*x*x;
                    Real f = sin(Pi*x);
                    integral += coeffs[i].weight*f*scale;
                  }
                // index to next point
                xCenter += dx;
              } // end loop over "cells"

            //Real exactVal = pow(length, 5)/5.0;
            Real exactVal =  (1 - cos(Pi*length))/Pi;
            if (SpaceDim < 3) exactVal = 0.0;
            Real error = exactVal - integral;
            pout() << "       " << nx
                   << ",  Newton-Cotes computed = " << integral
                   << ", exact value = " << exactVal
                   << ", error = " << error;

            if (ref > 0)
              {
                Real rate = log(Abs(oldNCErr/error))/log(2.0);
                pout () << "   Rate = " << rate;
                avgRateNC += rate;
              }

            pout() << endl;

            oldNCErr = error;
          }


          nx *= 2;
        } // end loop over refinements

      avgRateGaussian /= numRates;
      avgRateNC /= numRates;
      pout() << endl;
      pout()  << numPoints << "-point Gaussian Quadrature avg rate = " << avgRateGaussian << endl;
      pout()  << numPoints << "-point Newton-Cotes Quadrature avg rate = " << avgRateNC << endl;
      pout() << endl;

    } // end loop over numPoints

  return returnCode;

}

