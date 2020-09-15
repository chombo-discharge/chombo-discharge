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
#include "EulerRK4.H"
#include "LevelData.H"
#include "BaseFab.H"
#include "FArrayBox.H"
#include "ParmParse.H"
#include "LoadBalance.H"
#include "BRMeshRefine.H"
#include "ProtoInterface.H"
#include "Box.H"
#include "AMRIO.H"

#define PI 3.141592653589793
#define NUMCOMPS DIM+2

typedef Proto::Var<double,DIM> V;
typedef Proto::Var<double,NUMCOMPS> State;

typedef Proto::Box Bx;
using   Proto::Point;
using   Proto::BoxData;
using   Proto::Stencil;
using   Proto::RK4;
using ProtoCh::getProtoBox;
using ProtoCh::aliasBoxData;
using     std::cout;
using     std::endl;
using     std::shared_ptr;

class RunParams
{
public:
  RunParams()
  {
gamma    = 1.4;
    nstepmax = 0;
    nx       = 64;
    outinterv= 10;
    tmax     = 2.0;
    domsize  = 1.0;
    cfl      = 0.1;
    resetDx();
  }

  int nstepmax;
  int maxgrid;
  int nx;
  int outinterv;
  double tmax;
  double domsize;
  double gamma;
  double dx;
  double cfl;
  void resetDx()
  {
    dx = domsize/nx;
  }

  void coarsen(int a_refrat)
  {
    dx *= a_refrat;
    nx /= a_refrat;;
    
  }
  void print() const
  {

    pout() << "Compressible Euler with periodic bcs." << endl;

    pout() << "parameters: "                           << endl;
    pout() << "nx                  =  "   << nx        << endl;
    pout() << "max_grid            =  "   << maxgrid   << endl;
    pout() << "output interval     =  "   << outinterv << endl;
    pout() << "nstepmax            =  "   << nstepmax  << endl;
    pout() << "tmax                =  "   << tmax      << endl;
    pout() << "domain size         =  "   << domsize   << endl;
    pout() << "CFL number          =  "   << cfl       << endl;
    pout() << "gamma               =  "   << gamma     << endl;
  }
};                  
////////////
void
parseInputs(RunParams& a_params)
{
  ParmParse pp;
  pp.get("nx"                 , a_params.nx);
  pp.get("max_grid"           , a_params.maxgrid);
  pp.get("max_step"           , a_params.nstepmax);
  pp.get("output_interval"    , a_params.outinterv);
  pp.get("max_time"           , a_params.tmax);
  pp.get("domain_size"        , a_params.domsize);
  pp.get("cfl"                , a_params.cfl);
  a_params.resetDx();
  a_params.print();
}

PROTO_KERNEL_START 
unsigned int InitializeStateF(State& a_U,
                              const V& a_x)
{
  double gamma = 1.4;
  double rho0 = gamma;
  double p0 = 1.;
  double umag = .5/sqrt(1.*(DIM));
  double rho = rho0;
  for (int dir = 0; dir < DIM; dir++)
  {
    rho += .1*rho0*sin(2*2*PI*a_x(0));
  }
  double p = p0 + (rho - rho0)*gamma*p0/rho0;
  a_U(0) = rho;
  double ke = 0.;
  for (int dir = 1; dir <= DIM; dir++)
  {
    ke += umag*umag;
    a_U(dir) = rho*umag;
  }
  ke *= .5;
  a_U(NUMCOMPS-1) = p/(gamma-1.0) + rho*ke;
  return 0;
}
PROTO_KERNEL_END(InitializeStateF, InitializeState)

//=================================================================================================
PROTO_KERNEL_START
void iotaFuncF(Point           & a_p,
               V               & a_X,
               double            a_h)
{
  for (int ii = 0; ii < DIM; ii++)
  {
    a_X(ii) = a_p[ii]*a_h + 0.5*a_h;
  }
}
PROTO_KERNEL_END(iotaFuncF,iotaFunc)

/***/
void 
writeData(int step,LevelData<FArrayBox> & a_U)
{
  char filename[1024];
  sprintf(filename, "state.step%d.%d.hdf5",step, DIM);
  writeLevelname(&a_U, filename);
}
/***/
void eulerRun(const RunParams& a_params)
{
  CH_TIME("eulerRun");
  double tstop = a_params.tmax;
  int maxStep  = a_params.nstepmax;
  int nGhost = NGHOST;

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
  //this sneaks in parameters not covered by the RK4 template class
  //not the most elegant solution but it works for single level
  EulerOp::s_dx    = a_params.dx;
  EulerOp::s_gamma = a_params.gamma;

  shared_ptr<LevelData<FArrayBox> > Uptr(new LevelData<FArrayBox>(grids, NUMCOMPS, nGhost*IntVect::Unit));
  LevelData<FArrayBox>&  U = *Uptr;

  EulerState state(Uptr);
  RK4<EulerState,EulerRK4Op,EulerDX> rk4;
  
  Stencil<double> Lap2nd = Stencil<double>::Laplacian();

  pout() << "before initializestate"<< endl;
  DataIterator dit = grids.dataIterator();
#pragma omp parallel for
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    Box grid = grids[dit[ibox]];
    Box grownBox = U[dit[ibox]].box();
    Bx valid = getProtoBox(grid);
    Bx grnbx = getProtoBox(grownBox);
    BoxData<double, DIM> x(grnbx);
    forallInPlace_p(iotaFunc, grnbx, x, EulerOp::s_dx);
    BoxData<double, NUMCOMPS> ubd;
    aliasBoxData<double, NUMCOMPS>(ubd, U[dit[ibox]]);
    
    forallInPlace(InitializeState,grnbx,ubd,x);

    //smooth initial data
    BoxData<double, NUMCOMPS> lapu = Lap2nd(ubd, 1.0/24.0);
    ubd += lapu;
  }

  double maxwave = EulerOp::maxWave(*state.m_U);
  double dt = .25*a_params.cfl*a_params.dx/maxwave;
  pout() << "initial maximum wave speed = " << maxwave << ", dt = "<< dt << endl;

  pout() << "after initializestate"<< endl;
  U.exchange(state.m_exchangeCopier);

  double time = 0.;

  pout() << "starting time loop"<< endl;
  if(a_params.outinterv > 0)
  {
    writeData(0, U);
  }
  for (int k = 0;(k < maxStep) && (time < tstop);k++)
  {
    {
      CH_TIME("rk4_advance");
      rk4.advance(time,dt,state);
    }
    //this was computed during the advance.
    //so the standard trick is to reuse it.
    maxwave = state.m_velSave;

    time += dt;

    double dtnew = a_params.cfl*a_params.dx/maxwave; double dtold = dt;
    dt = min(1.1*dtold, dtnew);

    pout() <<"nstep = " << k << " time = " << time << ", dt = " << dt << endl;
    if((a_params.outinterv > 0) && (k%a_params.outinterv == 0))
    {
      writeData(k+1, U);
    }
  }

}
////////////
int main(int a_argc, char* a_argv[])
{

#ifdef CH_MPI
  MPI_Init(&a_argc, &a_argv);
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
    eulerRun(params);
  }
  
#ifdef CH_MPI
  dumpmemoryatexit();
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif

  return 0;
}
