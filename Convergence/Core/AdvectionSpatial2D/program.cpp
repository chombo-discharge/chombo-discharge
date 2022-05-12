#include <CD_Driver.H>
#include <CD_CdrGodunov.H>
#include <CD_RegularGeometry.H>
#include <CD_AdvectionDiffusionStepper.H>
#include <CD_AdvectionDiffusionTagger.H>
#include <CD_DischargeIO.H>
#include <ParmParse.H>


using namespace ChomboDischarge;
using namespace Physics::AdvectionDiffusion;

int
main(int argc, char* argv[])
{

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  // Build class options from input script and command line options
  const std::string input_file = argv[1];
  ParmParse         pp(argc - 2, argv + 2, NULL, input_file.c_str());

  // Set up the exact solution.
  ParmParse pp2("AdvectionDiffusion");
  ParmParse pp3("Driver");  
  Real omega;
  Real blobAmplitude;
  Real blobRadius;
  Real stopTime;
  RealVect blobCenter;
  Vector<Real> v;
  pp2.get("omega", omega);
  pp2.get("blob_amplitude", blobAmplitude);
  pp2.get("blob_radius", blobRadius);
  pp2.getarr("blob_center", v, 0, SpaceDim);
  pp3.get("stop_time", stopTime);
  blobCenter = RealVect(D_DECL(v[0],v[1],v[2]));

  // This is the initial shape. 
  auto initShape = [c=blobCenter,a=blobAmplitude,r=blobRadius](const RealVect& x) -> Real {
		     const RealVect d = x - c;
		     const Real d2 = d.dotProduct(d);
		     const Real r2 = r*r;

		     return a*exp(-0.5*d2*d2/(r2*r2));
		   };
  auto exactSolution = [w=omega, c=blobCenter, T = stopTime, s=initShape](const RealVect& x) -> Real {

			 // Transform input point to polar coordinates
			 const Real r = x.vectorLength();
			 const Real theta = atan2(x[1],x[0]);

			 // Angle in rotated coordinate system
			 const Real alpha = theta - w*T;

			 const RealVect y = r*RealVect(D_DECL(cos(alpha), sin(alpha), 0.0));


			 return s(RealVect(D_DECL(y[0], y[1], 0.0)));
		       };


  // Initialize norms and grid cells. 
  std::vector<IntVect> nCells{64 * IntVect::Unit,
                              128 * IntVect::Unit,
                              256 * IntVect::Unit,
                              512 * IntVect::Unit};

  std::vector<std::array<Real, 3>> norms;  

  // Set geometry and AMR
  RefCountedPtr<ComputationalGeometry> compgeom   = RefCountedPtr<ComputationalGeometry>(new RegularGeometry());
  RefCountedPtr<AmrMesh>               amr        = RefCountedPtr<AmrMesh>(new AmrMesh());
  RefCountedPtr<GeoCoarsener>          geocoarsen = RefCountedPtr<GeoCoarsener>(new GeoCoarsener());

  // Set up basic AdvectionDiffusion
  RefCountedPtr<CdrSolver>   solver      = RefCountedPtr<CdrSolver>(new CdrGodunov());
  RefCountedPtr<TimeStepper> timestepper = RefCountedPtr<TimeStepper>(new AdvectionDiffusionStepper(solver));
  RefCountedPtr<CellTagger>  tagger      = RefCountedPtr<CellTagger>(new AdvectionDiffusionTagger(solver, amr));

  // Run the various cases
  for (const auto& cells : nCells) {
    amr->setCoarsestGrid(cells);
    amr->buildDomains();  

    // Set up the Driver and run it
    RefCountedPtr<Driver> engine = RefCountedPtr<Driver>(new Driver(compgeom, timestepper, amr, tagger, geocoarsen));
    engine->setupAndRun(input_file);

    // Compute the solution error.
    EBAMRCellData error;
    EBAMRCellData exact;
    EBAMRCellData compu;
    EBAMRCellData outpu;

    amr->allocate(error, "primal", phase::gas, 1);
    amr->allocate(exact, "primal", phase::gas, 1);
    amr->allocate(compu, "primal", phase::gas, 1);
    amr->allocate(outpu, "primal", phase::gas, 3);    

    // Set exact/computed solution and the erro.
    DataOps::setValue(exact, exactSolution, amr->getProbLo(), amr->getDx(),0);
    DataOps::setValue(error, exactSolution, amr->getProbLo(), amr->getDx(),0);    
    
    compu.copy(solver->getPhi());

    DataOps::incr(error, compu, -1.0);

    exact[0]->copyTo(Interval(0,0), *outpu[0], Interval(0,0));
    compu[0]->copyTo(Interval(0,0), *outpu[0], Interval(1,1));
    error[0]->copyTo(Interval(0,0), *outpu[0], Interval(2,2));        
    

    // Compute the solution error norms
    const Real Linf = DataOps::norm(*error[0], amr->getDomains()[0], 0, true);
    const Real L1   = DataOps::norm(*error[0], amr->getDomains()[0], 1, true);
    const Real L2   = DataOps::norm(*error[0], amr->getDomains()[0], 2, true);

    norms.emplace_back(std::array<Real, 3>{Linf, L1, L2});    

    DischargeIO::writeEBHDF5(outpu, "output.hdf5");    
  }

  // Compute convergence rates
#ifdef CH_MPI
  if (procID() == 0) {
#endif
    std::cout << "# cells"
              << "\t"
              << "Linf error"
              << "\t"
              << "L1 error"
              << "\t"
              << "L2 error" << std::endl;
    for (int i = 0; i < norms.size(); i++) {
      std::cout << nCells[i][0] << "\t" << std::get<0>(norms[i]) << "\t" << std::get<1>(norms[i]) << "\t"
                << std::get<2>(norms[i]) << std::endl;
    }      
#ifdef CH_MPI
  }
#endif  

#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif
}
