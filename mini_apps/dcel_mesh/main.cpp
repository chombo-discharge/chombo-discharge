/*!
  @file  main.cpp
  @brief Example main file for running the chombo-streamer code
  @author Robert Marskar
*/

#include "dcel_vert.H"
#include "dcel_edge.H"

#include <ParmParse.H>


int main(int argc, char* argv[]){

#ifdef CH_MPI
  MPI_Init(&argc,&argv);
#endif

#ifdef CH_MPI
  MPI_Finalize();
#endif
}
