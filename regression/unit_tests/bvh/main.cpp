#include "driver.H"
#include "bvh.H"
#include "point_mass.H"

#include <iostream>
#include <fstream>
#include <chrono>
#include <random>

#include <ParmParse.H>

int main(int argc, char* argv[]){

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  // Build class options from input script and command line options
  char* input_file = argv[1];
  ParmParse parmParse(argc-2, argv+2, NULL, input_file);

  ParmParse pp("bvh");

  int seed, ppc, num_points;
  int min_weight, max_weight;

  pp.get("seed",       seed);
  pp.get("ppc",        ppc);
  pp.get("num_points", num_points);
  pp.get("min_weight", min_weight);
  pp.get("max_weight", max_weight);

  if(seed < 0) seed = std::chrono::system_clock::now().time_since_epoch().count();

  // Make RNG, uniform weights, and uniform positions
  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<Real> ranFloat(0.0, 1.0);
  std::uniform_int_distribution<std::mt19937::result_type> ranInt(min_weight, max_weight);


  // Create initial particles
  std::vector<point_mass> inputParticles(0);
  Real Mass = 0.0;
  while(inputParticles.size() < num_points){
    const RealVect pos = RealVect(D_DECL(ranFloat(rng), ranFloat(rng), ranFloat(rng)));
    const Real mass    = 1.0*ranInt(rng);

    if(pos[1] < pos[0]+0.25){
      inputParticles.push_back(point_mass(pos, mass, 0.0));
      Mass += mass;
    }
  }

#if 0
  const RealVect bigPos1 = RealVect(D_DECL(0.25, 0.25, .5));
  const RealVect bigPos2 = RealVect(D_DECL(0.5, 0.5, .5));
  const RealVect bigPos3 = RealVect(D_DECL(0.75, 0.75, .5));
  const Real bigMass = 1000.;
  Mass += 3*bigMass;
  inputParticles.push_back(point_mass(bigPos1, bigMass));
  inputParticles.push_back(point_mass(bigPos2, bigMass));
  inputParticles.push_back(point_mass(bigPos3, bigMass));
#endif


  // Do particle merging/splitting

  bvh_tree<point_mass> tree(inputParticles, Mass);
  const Real t0 = MPI_Wtime();
  tree.build_tree(ppc);
  const Real t1 = MPI_Wtime();

  // Create output particles
  std::vector<point_mass> outputParticles(0);
  const std::vector<std::shared_ptr<bvh_node<point_mass> > >& leaves = tree.get_leaves();
  for (int i = 0; i < leaves.size(); i++){
    const point_mass newParticle(leaves[i]->get_data());

    outputParticles.push_back(newParticle);
  }

  std::cout << "time per input particle = " << (t1-t0)/num_points << std::endl;

  // Write input and output particles
  ofstream inputPar, outputPar;
  
  inputPar.open ("input_particles.dat");
  outputPar.open ("output_particles.dat");

#if CH_SPACEDIM==2
    inputPar << "x" << "\t" << "y" << "\t" << "color" << "\n";
#else
    inputPar << "x" << "\t" << "y" << "\t" << "z" << "\t" << "color" << "\n";
#endif

  for (const auto& p : inputParticles){
    const Real m        = p.mass();
    const RealVect& pos = p.pos();
#if CH_SPACEDIM==2
    inputPar << pos[0] << "\t" << pos[1] << "\t" << m << "\n";
#else
    inputPar << pos[0] << "\t" << pos[1] << "\t" << pos[2] << "\t" << m << "\n";
#endif
  }

#if CH_SPACEDIM==2
  outputPar << "x" << "\t" << "y" << "\t" << "color" << "\n";
#else
  outputPar << "x" << "\t" << "y" << "\t" << "z" << "\t" << "color" << "\n";
#endif
  for (const auto& p : outputParticles){
    const Real m        = p.mass();
    const RealVect& pos = p.pos();
#if CH_SPACEDIM==2
    outputPar << pos[0] << "\t" << pos[1] << "\t" << m << "\n";
#else
    outputPar << pos[0] << "\t" << pos[1] << "\t" << pos[2] << "\t" << m << "\n";
#endif
  }
  
  inputPar.close();
  outputPar.close();
  

  
#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif
}
