/*!
  @file   mechshaft_main.cpp
  @brief  Main file for the mechanical shaft simulations. 
  @author Robert Marskar
*/

#include "plasma_engine.H"
#include "plasma_kinetics.H"
#include "rk2.H"
#include "field_tagger.H"
#include "morrow_lowke.H"
#include "morrow_lowke.H"
#include "mechanical_shaft.H"

#include <ParmParse.H>

/*!
  @brief Potential
*/
Real potential_curve(const Real a_time){
  return 8.E3;
}

int main(int argc, char* argv[]){

#ifdef CH_MPI
  MPI_Init(&argc,&argv);
#endif

  // Build argument list from input file
  char* inputFile = argv[1];
  ParmParse PP(argc-2,argv+2,NULL,inputFile);


  RefCountedPtr<physical_domain> physdom         = RefCountedPtr<physical_domain> (new physical_domain());
  RefCountedPtr<plasma_kinetics> plaskin         = RefCountedPtr<plasma_kinetics> (new morrow_lowke());
  RefCountedPtr<time_stepper> timestepper        = RefCountedPtr<time_stepper>(new rk2());
  RefCountedPtr<amr_mesh> amr                    = RefCountedPtr<amr_mesh> (new amr_mesh());
  RefCountedPtr<cell_tagger> tagger              = RefCountedPtr<cell_tagger> (new field_tagger());
  RefCountedPtr<computational_geometry> compgeom = RefCountedPtr<computational_geometry> (new mechanical_shaft());
  RefCountedPtr<plasma_engine> engine            = RefCountedPtr<plasma_engine> (new plasma_engine(physdom,
												   compgeom,
												   plaskin,
												   timestepper,
												   amr,
												   tagger));

  // Run plasma engine
  engine->set_potential(potential_curve);
  engine->setup_and_run();


#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif
}
