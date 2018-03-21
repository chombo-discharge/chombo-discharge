/*!
  @file  main.cpp
  @brief Example main file for running the chombo-streamer code
  @author Robert Marskar
*/

#include "plasma_engine.H"
#include "plasma_kinetics.H"
#include "splitstep_tga.H"
#include "wire_wire_tagger.H"
#include "wire_wire_geometry.H"
#include "air7.H"

#include <ParmParse.H>

/*!
  @brief Potential
*/
Real potential_curve(const Real a_time){
  Real potential = 6.E4;
  
  return potential;
}

int main(int argc, char* argv[]){

#ifdef CH_MPI
  MPI_Init(&argc,&argv);
#endif

  // Build argument list from input file
  char* inputFile = argv[1];
  ParmParse PP(argc-2,argv+2,NULL,inputFile);


  RefCountedPtr<physical_domain> physdom         = RefCountedPtr<physical_domain> (new physical_domain());
  RefCountedPtr<time_stepper> timestepper        = RefCountedPtr<time_stepper>(new splitstep_tga());
  RefCountedPtr<amr_mesh> amr                    = RefCountedPtr<amr_mesh> (new amr_mesh());
  RefCountedPtr<cell_tagger> tagger              = RefCountedPtr<cell_tagger> (new wire_wire_tagger());
  RefCountedPtr<computational_geometry> compgeom = RefCountedPtr<computational_geometry> (new wire_wire_geometry());
  RefCountedPtr<plasma_kinetics> plaskin         = RefCountedPtr<plasma_kinetics> (new air7());
  RefCountedPtr<plasma_engine> engine            = RefCountedPtr<plasma_engine> (new plasma_engine(physdom,
												   compgeom,
												   plaskin,
												   timestepper,
												   amr,
												   tagger));
  engine->set_potential(potential_curve);
  engine->setup_and_run();


#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif
}
