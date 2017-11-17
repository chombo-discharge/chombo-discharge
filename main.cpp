/*!
  @file main.cpp
  @brief Main file for running the chombo-streamer code
  @author Robert Marskar
*/

// Chombo files
#include "ParmParse.H"
#include "EBIndexSpace.H"

#include "ComputationalGeometry.H"
#include "sphere_sphere_geometry.H"

//
int main(int argc, char* argv[]){

#ifdef CH_MPI
  MPI_Init(&argc,&argv);
#endif

  // --------------------------------------------
  // EBIndexSpace for ion domain and field domain
  // --------------------------------------------
  EBIndexSpace* fieldEBIS = new EBIndexSpace();

  sphere_sphere_geometry::dumpScript();

#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif
}
