/*!
  @file  main.cpp
  @brief Example main file for running the chombo-streamer code
  @author Robert Marskar
*/

#include "dcel_vert.H"
#include "dcel_edge.H"
#include "dcel_poly.H"
#include "dcel_mesh.H"

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

  const RealVect x0 = RealVect::Zero;
  const RealVect x1 = RealVect(BASISV(0));
  const RealVect x2 = RealVect(BASISV(1));
  const RealVect x3 = RealVect(BASISV(2));

  const RealVect n0 = -RealVect(BASISV(2)); // Normal vector for face 0
  const RealVect n1 = -RealVect(BASISV(1)); // Normal vector for face 1
  const RealVect n2 = -RealVect(BASISV(0)); // Normal vector for face 2
  const RealVect n3 =  RealVect(1,1,1);     // Normal vector for face 3

  dcel_vert* v0 = new dcel_vert();
  dcel_vert* v1 = new dcel_vert();
  dcel_vert* v2 = new dcel_vert();
  dcel_vert* v3 = new dcel_vert();

  dcel_edge* e0  = new dcel_edge(); // Inside p0, points towards v0. next=e1,  prev=e2,  pair=e3
  dcel_edge* e1  = new dcel_edge(); // Inside p0, points towards v2. next=e2,  prev=e0,  pair=
  dcel_edge* e2  = new dcel_edge(); // Inside p0, points towards v1. next=e0,  prev=e1,  pair=
  dcel_edge* e3  = new dcel_edge(); // Inside p1, points towards v1, next=e4,  prev=e5,  pair=e0
  dcel_edge* e4  = new dcel_edge(); // Inside p1, points towards v3, next=e5,  prev=e3,  pair=
  dcel_edge* e5  = new dcel_edge(); // Inside p1, points towards v0, next=e3,  prev=e4,  pair=
  dcel_edge* e6  = new dcel_edge(); // Inside p2, points towards v3, next=e7,  prev=e8,  pair=e5
  dcel_edge* e7  = new dcel_edge(); // Inside p2, points towards v2, next=e8,  prev=e6,  pair=
  dcel_edge* e8  = new dcel_edge(); // Inside p2, points towards v0, next=e6,  prev=e7,  pair=e1
  dcel_edge* e9  = new dcel_edge(); // Inside p3, points towards v1, next=e10, prev=e11, pair=e4
  dcel_edge* e10 = new dcel_edge(); // Inside p3, points towards v2, next=e11, prev=e9,  pair=e2
  dcel_edge* e11 = new dcel_edge(); // Inside p3, points towards v3, next=e9,  prev=e10, pair=e7

  dcel_poly* p0 = new dcel_poly(); // Triangle in xy-plane
  dcel_poly* p1 = new dcel_poly(); // Triangle in xz plane
  dcel_poly* p2 = new dcel_poly(); // Triangle in yz-plane
  dcel_poly* p3 = new dcel_poly(); // "Diagonal" triangle


  // Define vertices. The edge is an outgoing edge. 
  v0->define(x0, e1);
  v1->define(x1, e0);
  v2->define(x2, e2);
  v3->define(x3, e5);

  // Define half edges, forget normal vector for now
  e0->define(v0,  e3,   e1,  e2);
  e1->define(v2,  e8,   e2,  e0);
  e2->define(v1,  e10,  e0,  e1);
  e3->define(v1,  e0,   e4,  e5);
  e4->define(v3,  e9,   e5,  e3);
  e5->define(v0,  e6,   e3,  e4);
  e6->define(v3,  e5,   e7,  e8);
  e7->define(v2,  e11,  e8,  e6);
  e8->define(v0,  e1,   e6,  e7);
  e9->define(v1,  e4,   e10, e11);
  e10->define(v2, e2,   e11, e9);
  e11->define(v3, e7,   e9,  e10);

  // Define polygons
  p0->define(n0, e0);
  p1->define(n1, e3);
  p2->define(n2, e6);
  p3->define(n3, e9);

  // Define mesh from existing polygons and edges
  Vector<dcel_poly*> polygons;
  polygons.push_back(p0);
  polygons.push_back(p1);
  polygons.push_back(p2);
  polygons.push_back(p3);

  Vector<dcel_edge*> edges;
  edges.push_back(e0);
  edges.push_back(e1);
  edges.push_back(e2);
  edges.push_back(e3);
  edges.push_back(e4);
  edges.push_back(e5);
  edges.push_back(e6);
  edges.push_back(e7);
  edges.push_back(e8);
  edges.push_back(e9);
  edges.push_back(e10);
  edges.push_back(e11);

  Vector<dcel_vert*> vertices;
  vertices.push_back(v0);
  vertices.push_back(v1);
  vertices.push_back(v2);
  vertices.push_back(v3);

  // Create poly mesh and find the distance to a point

  dcel_mesh* mesh = new dcel_mesh(polygons, edges, vertices);
  mesh->reconcile_polygons();

  RefCountedPtr<physical_domain> physdom         = RefCountedPtr<physical_domain> (new physical_domain());
  RefCountedPtr<time_stepper> timestepper        = RefCountedPtr<time_stepper>(new rk2());
  RefCountedPtr<amr_mesh> amr                    = RefCountedPtr<amr_mesh> (new amr_mesh());
  RefCountedPtr<cell_tagger> tagger              = RefCountedPtr<cell_tagger> (new field_tagger());
  RefCountedPtr<computational_geometry> compgeom = RefCountedPtr<computational_geometry> (new dcel_geometry(mesh));
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
