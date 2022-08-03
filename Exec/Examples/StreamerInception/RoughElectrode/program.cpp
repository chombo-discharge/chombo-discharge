#include <CD_Driver.H>
#include <CD_TracerParticleSolver.H>
#include <CD_RoughSphere.H>
#include <CD_FieldSolverMultigrid.H>
#include <CD_StreamerInceptionStepper.H>
#include <CD_LookupTable.H>
#include <CD_DataParser.H>
#include <ParmParse.H>

using namespace ChomboDischarge;
using namespace Physics::StreamerInception;

int main(int argc, char* argv[]){

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  // Build class options from input script and command line options
  const std::string input_file = argv[1];
  ParmParse pp(argc-2, argv+2, NULL, input_file.c_str());

  // Read ionization and attachment coefficients and make them into functions.
  constexpr Real N = 2.45E25;
  
  LookupTable<2> ionizationData = DataParser::fractionalFileReadASCII("transport_data.txt",
								      "E/N (Td)	Townsend ioniz. coef. alpha/N (m2)",
								      "");
  LookupTable<2> attachmentData = DataParser::fractionalFileReadASCII("transport_data.txt",
								      "E/N (Td)	Townsend attach. coef. eta/N (m2)",
								      "");

  ionizationData.sort(0);
  attachmentData.sort(0);
  
  ionizationData.setTableSpacing(TableSpacing::Exponential);
  attachmentData.setTableSpacing(TableSpacing::Exponential);

  ionizationData.scale<0>(N * 1.E-21);
  attachmentData.scale<0>(N * 1.E-21);

  ionizationData.scale<1>(N);
  attachmentData.scale<1>(N);
  
  ionizationData.makeUniform(500);
  attachmentData.makeUniform(500);
  
  auto alphaEff = [&](const Real& E) -> Real {
    const Real alpha = ionizationData.getEntry<1>(E);
    const Real eta = attachmentData.getEntry<1>(E);

    return alpha - eta;
  };  

  // Set geometry and AMR 
  RefCountedPtr<ComputationalGeometry> compgeom = RefCountedPtr<ComputationalGeometry> (new RoughSphere());
  RefCountedPtr<AmrMesh> amr                    = RefCountedPtr<AmrMesh> (new AmrMesh());

  // Set up time stepper 
  auto timestepper = RefCountedPtr<StreamerInceptionStepper<>> (new StreamerInceptionStepper<>());

  // Set ionization coefficient. 
  timestepper->setAlpha(alphaEff);  

  // Set up the Driver and run it
  RefCountedPtr<Driver> engine = RefCountedPtr<Driver> (new Driver(compgeom, timestepper, amr));
  engine->setupAndRun(input_file);

#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif
}
