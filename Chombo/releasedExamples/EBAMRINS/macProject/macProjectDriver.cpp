#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>
using std::cerr;

#include "ParmParse.H"
#include "LoadBalance.H"
#include "BRMeshRefine.H"

#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "EBCellFactory.H"

#include "EBFABView.H"
#include "EBDebugDump.H"

#include "EBLevelDataOps.H"
#include "BaseIF.H"
#include "SphereIF.H"
#include "PlaneIF.H"
#include "GeometryShop.H"
#include "DebugDump.H"
#include "EBDebugDump.H"
#include "EBLevelDataOps.H"
#include "LoadBalance.H"
#include "EBLevelGrid.H"
#include "CH_Timer.H"
#include "EBPhysIBC.H"
#include "EBArith.H"
#include "ExtrapAdvectBC.H"
#include "DirichletPoissonDomainBC.H"
#include   "NeumannPoissonDomainBC.H"
#include   "NeumannPoissonEBBC.H"
#include   "EBSimpleSolver.H"
#include "WrappedGShop.H"
#include "SetupFunctions.H"
#include "EBLevelMACProjector.H"
#include "memusage.H"
#include "memtrack.H"

#define PI 3.141592653589793

                
/********/
void runProjection()
{
  CH_TIME("run advection");

  Real coveredval = -1;

  int nx      = 32;

  Real max_vel_mag = 1.0;
  Real max_vel_rad = 0.25;
  int nStream    = 8;
  ParmParse pp;

  pp.get("nstream", nStream);
  pp.get("covered_value", coveredval);

  pp.get("max_vel_mag"  , max_vel_mag);
  pp.get("max_vel_rad"  , max_vel_rad);

  pout() << "num_streams     = " << nStream         << endl;
  pout() << "max_vel_mag     = " << max_vel_mag     << endl;
  pout() << "max_vel_rad     = " << max_vel_rad     << endl;

  Real dx;
  DisjointBoxLayout grids;

  pout() << "defining geometry" << endl;


  Real geomCen;
  Real geomRad;
  int whichGeom;
  defineGeometry(grids, dx, geomCen, geomRad, whichGeom, nx);
  int nghost = 4;
  IntVect dataGhostIV =   nghost*IntVect::Unit;


  Box domain = grids.physDomain().domainBox();
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  EBLevelGrid eblg(grids, domain, nghost, ebisPtr);

  pout() << "inititializing data"   << endl;

  EBFluxFactory fluxfact(eblg.getEBISL());
  LevelData<EBFluxFAB>  fcVel(grids,        1, dataGhostIV,  fluxfact);
  LevelData<EBFluxFAB>  fcGph(grids,        1, dataGhostIV,  fluxfact);
  
  InitializeVel(fcVel, eblg, dx, geomCen, geomRad, max_vel_mag, max_vel_rad);
  
  int mgCycle = 1;
  int numSmooth;
  int maxIter;
  Real tolerance;
  pp.get("mgCycle", mgCycle);
  pp.get("numSmooth", numSmooth);
  pp.get("maxIter", maxIter);
  pp.get("tolerance", tolerance);
  int maxDepth = -1;
  Real time = 0;

  pout() << "initializing projection " << endl;
  NeumannPoissonEBBCFactory*       neumEBPtr = new   NeumannPoissonEBBCFactory();
  NeumannPoissonDomainBCFactory*   neumDMPtr = new   NeumannPoissonDomainBCFactory();
  DirichletPoissonDomainBCFactory* diriDMPtr = new DirichletPoissonDomainBCFactory();
  neumEBPtr->setValue(0.);
  neumDMPtr->setValue(0.);
  diriDMPtr->setValue(0.);

  RefCountedPtr<LinearSolver<LevelData<EBCellFAB> > > bottomSolverPtr(new EBSimpleSolver());

  //LinearSolver<LevelData<EBCellFAB> >  baseSolver = static_cast<LinearSolver< LevelData<EBCellFAB> >& >(bottomsolver);

  RefCountedPtr<BaseEBBCFactory>         EBBCPhi(neumEBPtr);
  RefCountedPtr<BaseDomainBCFactory> domainBCPhi(neumDMPtr);
  RefCountedPtr<BaseDomainBCFactory> domainBCVel(diriDMPtr);

  EBLevelMACProjector projector(eblg.getDBL(),
                                eblg.getEBISL(),
                                eblg.getDomain(),
                                dx*RealVect::Unit,
                                RealVect::Zero,
                                *bottomSolverPtr,
                                EBBCPhi,
                                domainBCPhi,
                                domainBCVel,
                                numSmooth, 
                                mgCycle,
                                maxIter,
                                tolerance, maxDepth, time,
                                dataGhostIV, dataGhostIV);

  pout() << "projecting velocity field " << endl;

  projector.project(fcVel, fcGph);

  pout() << "leaving " << endl;

}
/******/
int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  // Scoping trick
  {
#ifdef CH_MPI
    MPI_Barrier(Chombo_MPI::comm);
#endif

    if (argc < 2)
      {
        cerr << " usage " << argv[0] << " <input_file_name> " << endl;
        exit(0);
      }

    char* inFile = argv[1];
    ParmParse pp(argc-2,argv+2,NULL,inFile);

    //solve the stinking problem and output everything
    runProjection();

    EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
    ebisPtr->clear();

  }
  // End scoping trick

#ifdef CH_MPI
  CH_TIMER_REPORT();
  dumpmemoryatexit();
  MPI_Finalize();
#endif
}
