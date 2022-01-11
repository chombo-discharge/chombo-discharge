#include "CD_Driver.H"
#include <CD_GeoCoarsener.H>
#include <CD_FieldSolverFactory.H>
#include <CD_FieldSolverMultigrid.H>
#include <CD_CdrLayoutImplem.H>
#include <CD_CdrGodunov.H>
#include <CD_RtLayoutImplem.H>
#include <CD_EddingtonSP1.H>
#include <CD_CdrPlasmaGenericModel.H>
#include <CD_RegularGeometry.H>
#include <CD_CdrPlasmaGodunovStepper.H>
#include <CD_CdrPlasmaStreamerTagger.H>
#include "ParmParse.H"

#include <nlohmann/json.hpp>

using json = nlohmann::json;

// This is the potential curve (constant in this case). Modify it if you want to.
Real g_potential;
Real potential_curve(const Real a_time){
  return g_potential;
}

using namespace ChomboDischarge;
using namespace Physics::CdrPlasma;

int main(int argc, char* argv[]){

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  // Build class options from input script and command line options
  const std::string input_file = argv[1];
  ParmParse pp(argc-2, argv+2, NULL, input_file.c_str());  


  CdrPlasmaGenericModel model = CdrPlasmaGenericModel();


#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif
}
