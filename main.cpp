/*!
  @file main.cpp
  @brief Main file for running the chombo-streamer code
  @author Robert Marskar
*/

// Chombo files
#include "ParmParse.H"
#include "EBIndexSpace.H"

#include "computational_geometry.H"
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

  const int nCells = 512;
  const RealVect& a_probLo = -RealVect::Unit;
  physical_domain physdom(-RealVect::Unit, RealVect::Unit);
  ProblemDomain probdom(IntVect::Zero, (nCells - 1)*IntVect::Unit);
  const Real& finestdx = (physdom.get_prob_lo()[0] - physdom.get_prob_hi()[0])/nCells;

  computational_geometry* compgeom = static_cast<computational_geometry*> (new sphere_sphere_geometry());
  compgeom->build_geometries(physdom, probdom, finestdx, 8);

  // Real r;
  // ParmParse pp("sphere");
  // pp.get("electrode_radius", r);

#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif
}
