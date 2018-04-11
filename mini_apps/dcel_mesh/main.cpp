/*!
  @file  main.cpp
  @brief Example main file for running the chombo-streamer code
  @author Robert Marskar
*/

#include "dcel_vert.H"
#include "dcel_edge.H"
#include "dcel_poly.H"
#include "dcel_mesh.H"
#include "ply_reader.H" 

#include "plasma_engine.H"
#include "plasma_kinetics.H"
#include "rk2.H"
#include "field_tagger.H"
#include "air7.H"
#include "dcel_geometry.H"

#include <ParmParse.H>

Real potential_curve(const Real a_time){
  return 0.0;
}

int main(int argc, char* argv[]){

#ifdef CH_MPI
  MPI_Init(&argc,&argv);
#endif

  // Build argument list from input file
  char* inputFile = argv[1];
  ParmParse PP(argc-2,argv+2,NULL,inputFile);

  // Build surface mesh
  std::string str;
  int tree_depth;
  int max_elem;
  ParmParse pp("dcel");
  pp.get("which_file", str);
  pp.get("tree_depth", tree_depth);
  pp.get("max_elem", max_elem);
  dcel_mesh* plymesh = new dcel_mesh();
  ply_reader::read_ascii(*plymesh, str);
  plymesh->reconcile_polygons(false);
  plymesh->build_tree(tree_depth, max_elem);


  RefCountedPtr<physical_domain> physdom         = RefCountedPtr<physical_domain> (new physical_domain());
  RefCountedPtr<time_stepper> timestepper        = RefCountedPtr<time_stepper>(new rk2());
  RefCountedPtr<amr_mesh> amr                    = RefCountedPtr<amr_mesh> (new amr_mesh());
  RefCountedPtr<cell_tagger> tagger              = RefCountedPtr<cell_tagger> (new field_tagger());
  RefCountedPtr<computational_geometry> compgeom = RefCountedPtr<computational_geometry> (new dcel_geometry(plymesh));
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
  MPI_Finalize();
#endif
}
