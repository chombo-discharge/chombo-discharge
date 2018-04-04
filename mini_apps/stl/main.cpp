/*!
  @file  main.cpp
  @brief Example main file for running the chombo-streamer code
  @author Robert Marskar
*/

#include "triangle.H"

#include <ParmParse.H>

int main(int argc, char* argv[]){

#ifdef CH_MPI
  MPI_Init(&argc,&argv);
#endif

  // Build argument list from input file
  char* inputFile = argv[1];
  ParmParse PP(argc-2,argv+2,NULL,inputFile);


  const RealVect x0 = RealVect::Zero;      // 
  const RealVect x1 = RealVect(BASISV(0));
  const RealVect x2 = RealVect(BASISV(1));
  const RealVect n  = -RealVect(BASISV(2));
  
  triangle tri(x0, x1, x2, n);

  const RealVect test1 = (x0 + x1 + x2)/3.0 + 0.0*RealVect(BASISV(2)); // Right above the centroid. Distance should be 1
  const RealVect test2 = 1.1*RealVect(BASISV(0)) + 0.1*RealVect(BASISV(1)) + RealVect(BASISV(2)); // Outside edge
  const RealVect test3 = 1.1*RealVect(BASISV(0)) + 0.0*RealVect(BASISV(1)) + RealVect(BASISV(2)); // Outside vertex
  

  pout() << "\n Inside test" << endl;
  pout() << tri.signed_distance(test1) << endl;

  pout() << "\n Outside edge test" << endl;
  pout() << tri.signed_distance(test2) << endl;

  pout() << "\n Outside vertex test" << endl;
  pout() << tri.signed_distance(test3) << endl;

  

#ifdef CH_MPI
  MPI_Finalize();
#endif
}
