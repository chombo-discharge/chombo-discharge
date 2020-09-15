#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cstdio>
#include <cstring>
#include <cassert>
#include <cmath>

#include <vector>
#include <memory>

#include <iostream>
#include <fstream>
#include <sstream>

#include "Proto.H"
#include "GodunovAdvectionOp.H"
#include "BCG_Integrator.H"
#include "SGMultigrid.H"
#include "LevelData.H"
#include "BaseFab.H"
#include "ParmParse.H"
#include "Proto.H"
#include "ProtoInterface.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "DebugDump.H"
#include "AMRIO.H"

#define PI 3.141592653589793

typedef Proto::Var<double,   1> Scalar;
typedef Proto::Var<double, DIM> PVector;
typedef Proto::Box Bx;
using   Proto::Point;
using   Proto::BoxData;
using     std::cout;
using     std::endl;
class RunParams
{
public:
  RunParams()
  {
    nstepmax = 0;
    nx       = 64;
    outinterv= 10;
    tmax     = 2.0;
    domsize  = 1.0;
    vortrad = 0.125*domsize;
    vortloc = 0.5*domsize;
    cfl      = 0.1;
    viscosity= 0.0002;
    presiter = 1;
    resetDx();
  }

  int nstepmax;
  int maxgrid;
  int presiter;
  int nx;
  int outinterv;
  double viscosity;
  double tmax;
  double domsize;
  double vortrad;
  double vortloc;
  double dx;
  double cfl;
  void resetDx()
  {
    dx = domsize/nx;
    vortrad = 0.125*domsize;
    vortloc = 0.5*domsize;
  }

  void coarsen(int a_refrat)
  {
    dx *= a_refrat;
    nx /= a_refrat;;
    
  }
  void print() const
  {

    cout << "Navier Stokes simulation weird litttle vortex.  Periodic bcs." << endl;
    cout << "parameters: "                           << endl;
    cout << "nx                  =  "   << nx        << endl;
    cout << "max_grid            =  "   << maxgrid   << endl;
    cout << "output interval     =  "   << outinterv << endl;
    cout << "nstepmax            =  "   << nstepmax  << endl;
    cout << "tmax                =  "   << tmax      << endl;
    cout << "domain size         =  "   << domsize   << endl;
    cout << "viscosity           =  "   << viscosity << endl;
    cout << "pressure iterations =  "   << presiter  << endl;
    cout << "CFL number          =  "   << cfl       << endl;
  }
};                  
////////////
void
parseInputs(RunParams& a_params)
{
  ParmParse pp;
  pp.get("nx"                 , a_params.nx);
  pp.get("max_grid"           , a_params.maxgrid);
  pp.get("pressure_iterations", a_params.presiter);
  pp.get("max_step"           , a_params.nstepmax);
  pp.get("output_interval"    , a_params.outinterv);
  pp.get("viscosity"          , a_params.viscosity);
  pp.get("max_time"           , a_params.tmax);
  pp.get("vortex_radius"      , a_params.vortrad);
  pp.get("vortex_center"      , a_params.vortloc);
  pp.get("domain_size"        , a_params.domsize);
  pp.get("cfl"                , a_params.cfl);
  a_params.resetDx();
  a_params.print();
}

PROTO_KERNEL_START unsigned int InitializeVelF(Point&            a_p,
                                               PVector         & a_U,
                                               RunParams         a_params)
{
  
  double xrel[DIM];
  double radsq = 0;
  double rad0sq = a_params.vortrad*a_params.vortrad;
  {
    double  x[DIM];
    double x0[DIM];
    for(int idir = 0; idir < DIM; idir++)
    {    
      x[idir] = (a_p[idir] + 0.5)*a_params.dx;
      x0[idir] = a_params.vortloc;
    }


    for(int idir = 0; idir < DIM; idir++)
    {    
      xrel[idir] = x[idir] - x0[idir];
    }
    for(int idir = 0; idir < DIM; idir++)
    {  
      radsq += xrel[idir]*xrel[idir];
    }
  }
  for(int idir = 0; idir < DIM; idir++)
  {  
    double uval = 0;
    if(radsq < rad0sq)
    {
      double cosval = cos(0.5*PI*(radsq/rad0sq));
      double phival = cosval*cosval;
      if(idir == 0)
      {
        uval = -(xrel[1])*phival;
      }
      else if(idir == 1)
      {
        uval =  xrel[0]*phival;
      }
    }
    else
    {
      uval = 0.0;
    }
    a_U(idir) = uval;
  }
  return 0;
}
PROTO_KERNEL_END(InitializeVelF, InitializeVel)



///
void getDtAndMaxWave(double                     & a_dt,
                     double                     & a_maxwave,
                     const LevelData<FArrayBox> & a_velocity,
                     const RunParams            & a_params)
{
  a_maxwave = norm(a_velocity, Interval(0, DIM-1), 0);
  a_dt = a_params.cfl*a_params.dx/a_maxwave;
}
void writeFile(const LevelData<FArrayBox> & a_data,
               const string               & a_basename,
               const int                  & a_stepnum)
{
  char filename[1024];
  sprintf(filename, "%s.step%d.%d.hdf5", a_basename.c_str(), a_stepnum, DIM);
  writeLevelname(&a_data, filename);
}
///
void navierRun(const RunParams& a_params)
{
  int nghost = 3;
  IntVect domLo = IntVect::Zero;
  IntVect domHi  = (a_params.nx - 1)*IntVect::Unit;
  constexpr bool is_periodic[] = {true, true, true};

  ProblemDomain domain(domLo, domHi, is_periodic);

  Vector<Box> boxes;
  unsigned int blockfactor = 8;
  domainSplit(domain, boxes, a_params.maxgrid, blockfactor);
  Vector<int> procs;
  LoadBalance(procs, boxes);
  DisjointBoxLayout grids(boxes, procs, domain);
  LevelData<FArrayBox > gradpres(grids, DIM, nghost*IntVect::Unit);
  LevelData<FArrayBox > velocity(grids, DIM, nghost*IntVect::Unit);
  LevelData<FArrayBox >  scratch(grids, DIM, IntVect::Zero);
#if DIM==2
  LevelData<FArrayBox > vorticity(grids, 1,   IntVect::Zero);
#else
  LevelData<FArrayBox > vorticity(grids, DIM, IntVect::Zero);
#endif

  DataIterator dit = grids.dataIterator();
#pragma omp parallel
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    Bx patch = grids[dit[ibox]];
    BoxData<double, DIM> velbd;
    ProtoCh::aliasBoxData<double, DIM>(velbd, velocity[dit[ibox]]);
    forallInPlace_p(InitializeVel, patch, velbd, a_params);
  }
  velocity.exchange();
  
  //run the simulation
  int   nstep = 0; 
  double time = 0;
  BCG_Integrator op(grids, a_params.dx, a_params.viscosity, nghost);

  cout << "initial projection" << endl;
  op.ccProject(velocity, gradpres);
#pragma omp parallel
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    gradpres[dit[ibox]].setVal(0.);
  }
  cout << "iterating for initial pressure" << endl;
  double maxwave = 0;
  double dt = 0;
  getDtAndMaxWave(dt, maxwave, velocity, a_params);

  for(int iiter = 0; iiter < a_params.presiter; iiter++)
  {
    velocity.copyTo(scratch);
    op.advanceSolution(velocity, gradpres, dt);
    scratch.copyTo(velocity);
  }
  cout << "starting navier stokes run" << endl;
  if(a_params.outinterv >= 0)
  {
    op.computeVorticity(vorticity, velocity);
    writeFile(velocity, string("velo"), nstep);
    writeFile(gradpres, string("velo"), nstep);
    writeFile(vorticity,string("vort"), nstep);
  }

  while((nstep < a_params.nstepmax) && (time < a_params.tmax))
  {

    cout << "start of step = "     << nstep << ", old time = " << time << ", maxwave = " << maxwave << " dt = "  << dt  << endl; 

    op.advanceSolution(velocity, gradpres, dt);
    
    nstep++;
    time += dt;
    getDtAndMaxWave(dt, maxwave, velocity, a_params);

    if(a_params.outinterv >= 0)
    {
      op.computeVorticity(vorticity, velocity);
      writeFile(velocity, string("velo"), nstep);
      writeFile(gradpres, string("velo"), nstep);
      writeFile(vorticity,string("vort"), nstep);
    }
  }

  cout << "finished navier run" << endl;

}


////////////
int main(int a_argc, char* a_argv[])
{

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  {
    if (a_argc < 2)
    {
      cerr<< " usage " << a_argv[0] << " <input_file_name> " << endl;
      exit(0);
    }
    char* in_file = a_argv[1];
    ParmParse  pp(a_argc-2,a_argv+2,NULL,in_file);

    RunParams params;
    parseInputs(params);
    navierRun(params);
  }
  
#ifdef CH_MPI
  dumpmemoryatexit();
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif

  return 0;
}

