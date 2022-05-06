#include <CD_Driver.H>
#include <CD_FieldSolverMultigrid.H>
#include <CD_CoaxialCable.H>
#include <CD_FieldStepper.H>
#include <CD_DischargeIO.H>
#include <ParmParse.H>

using namespace ChomboDischarge;
using namespace Physics::Electrostatics;

// Parameters from input script
constexpr Real R0   = 1.0;
constexpr Real R1   = 2.0;
constexpr Real R2   = 3.0;
constexpr Real eps  = 5.0;
constexpr Real phi0 = 1.0;
constexpr Real a    = phi0/( log(R1/R2) - log(R1/R0)/eps );


// Exact solution for "coaxial cable" with relative permittivity eps and potential phi0 on the inner conductor.
Real exactSolution(const RealVect& x) {
  const Real r = sqrt(x[0]*x[0] + x[1]*x[1]);

  Real phi = 0.0;
    
  if(r >= R0 && r <= R1) {
    phi = phi0 + a/eps*log(r/R0);
  }
  else if(r >= R1 && r <= R2) {
    phi = a * log(r/R2);
  }

  return phi;
}

int main(int argc, char* argv[]){

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  // Build class options from input script and command line options
  const std::string input_file = argv[1];
  ParmParse pp(argc-2, argv+2, NULL, input_file.c_str());

  // Set geometry and AMR 
  RefCountedPtr<ComputationalGeometry> compgeom = RefCountedPtr<ComputationalGeometry> (new CoaxialCable());
  RefCountedPtr<AmrMesh> amr                    = RefCountedPtr<AmrMesh> (new AmrMesh());
  RefCountedPtr<GeoCoarsener> geocoarsen        = RefCountedPtr<GeoCoarsener> (new GeoCoarsener());
  RefCountedPtr<CellTagger> tagger              = RefCountedPtr<CellTagger> (NULL);

  // Set up basic Poisson, potential = 1 
  auto timestepper = RefCountedPtr<FieldStepper<FieldSolverMultigrid> >
     (new FieldStepper<FieldSolverMultigrid>());

  // Set up the Driver and run it
  RefCountedPtr<Driver> engine = RefCountedPtr<Driver> (new Driver(compgeom, timestepper, amr, tagger, geocoarsen));
  engine->setupAndRun(input_file);

  // Compute the rror
  MFAMRCellData err;
  amr->allocate(err, "primal", 1);
  DataOps::setValue(err, exactSolution, amr->getProbLo(), amr->getDx(), 0);
  DataOps::incr(err, timestepper->getPotential(), -1.0);

  // Extract error from gas and solid sides
  EBAMRCellData errGas = amr->alias(phase::gas,   err);
  EBAMRCellData errSol = amr->alias(phase::solid, err);  

  // Compute the various norms
  const Real Linf = DataOps::norm(*errGas[0], amr->getDomains()[0], 0, true);
  const Real L1   = DataOps::norm(*errGas[0], amr->getDomains()[0], 1, true);
  const Real L2   = DataOps::norm(*errGas[0], amr->getDomains()[0], 2, true);  

  if(procID() == 0) {
    std::cout << "Linf = " << Linf << std::endl;
    std::cout << "L1   = " << L1   << std::endl;
    std::cout << "L2   = " << L2   << std::endl;    
  }
  
#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif
}
