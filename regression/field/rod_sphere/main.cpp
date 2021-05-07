#include "driver.H"
#include "field_solver_multigrid.H"
#include "rod_dielectric.H"
#include "field_stepper.H"
#include "ParmParse.H"
#include "LeastSquares.H"

using namespace physics::poisson;

int main(int argc, char* argv[]){

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  // Build class options from input script and command line options
  const std::string input_file = argv[1];
  ParmParse pp(argc-2, argv+2, NULL, input_file.c_str());

  // Set geometry and AMR 
  RefCountedPtr<computational_geometry> compgeom = RefCountedPtr<computational_geometry> (new rod_dielectric());
  RefCountedPtr<amr_mesh> amr                    = RefCountedPtr<amr_mesh> (new amr_mesh());
  RefCountedPtr<geo_coarsener> geocoarsen        = RefCountedPtr<geo_coarsener> (new geo_coarsener());
  RefCountedPtr<cell_tagger> tagger              = RefCountedPtr<cell_tagger> (NULL);

  // Set up basic Poisson, potential = 1 
  auto timestepper = RefCountedPtr<field_stepper<field_solver_multigrid> >
     (new field_stepper<field_solver_multigrid>());

  // Set up the driver and run it
  RefCountedPtr<driver> engine = RefCountedPtr<driver> (new driver(compgeom, timestepper, amr, tagger, geocoarsen));
  engine->setup_and_run(input_file);

#if 1 // Try to make a least squares interpolation
  RealVect p0( 0.0,  0.0);
  RealVect p1(-1, -1);
  RealVect p2(-1,  1);
  RealVect p3( 1,  1);
  RealVect p4( 1, -1);

  VolIndex v0(IntVect( 0,  0), 0);
  VolIndex v1(IntVect(-1, -1), 0);
  VolIndex v2(IntVect(-1,  1), 0);
  VolIndex v3(IntVect( 1,  1), 0);
  VolIndex v4(IntVect( 1, -1), 0);

  Vector<VolIndex> vofs;
  vofs.push_back(v1);
  vofs.push_back(v2);
  vofs.push_back(v3);
  vofs.push_back(v4);

  Vector<RealVect> deltas;
  deltas.push_back(p0-p1);
  deltas.push_back(p0-p2);
  deltas.push_back(p0-p3);
  deltas.push_back(p0-p4);

  Real f1 = 1.0;
  Real f2 = 1.0;
  Real f3 = 1.0;
  Real f4 = 1.0;

  VoFStencil sten = LeastSquares::computeInterpolationStencil(vofs, deltas, 0, 1);

  if(procID() == 0){
    Real f0 = 0.0;
    for (int i = 0; i < sten.size(); i++){
      const VolIndex& vof = sten.vof(i);
      const Real w        = sten.weight(i);

      std::cout << vof << "\t" << w << std::endl;
      if(vof == v1) f0 += f1*w;
      if(vof == v2) f0 += f2*w;
      if(vof == v3) f0 += f3*w;
      if(vof == v4) f0 += f4*w;
    }

    std::cout << f0 << std::endl;
  }
  
#endif

#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif
}
