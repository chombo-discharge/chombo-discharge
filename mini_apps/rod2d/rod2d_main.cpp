/*!
  @file  main.cpp
  @brief Example main file for running the chombo-streamer code
  @author Robert Marskar
*/

#include "plasma_engine.H"
#include "plasma_kinetics.H"
#include "rk2.H"
#include "field_tagger.H"


#include "air_bolsig.H"
#include "morrow_lowke.H"
#include "sphere_sphere_geometry.H"
#include "rod_slab.H"

#include <ParmParse.H>

/*!
  @brief Potential
*/
Real potential_curve(const Real a_time){
  Real potential = 4.E3;

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
  RefCountedPtr<time_stepper> timestepper        = RefCountedPtr<time_stepper>(new rk2());
  RefCountedPtr<amr_mesh> amr                    = RefCountedPtr<amr_mesh> (new amr_mesh());
  RefCountedPtr<cell_tagger> tagger              = RefCountedPtr<cell_tagger> (new field_tagger());
  RefCountedPtr<computational_geometry> compgeom = RefCountedPtr<computational_geometry> (new rod_slab());
  RefCountedPtr<plasma_kinetics> plaskin         = RefCountedPtr<plasma_kinetics> (new air_bolsig());
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
