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



typedef Proto::Var<double,   1> Scalar;
typedef Proto::Box Bx;
using   Proto::Point;
using   Proto::BoxData;
using     std::cout;
using     std::endl;

#define PI 3.141592653589793

class SolveParams
{
public:
  SolveParams()
  {
    relaxOnly = 0;
    maxiter   = 27;
    numsmooth = 4;
    nx        = 64;
    domsize   = 1.0;
    tol       = 1.0e-9;
    alpha = 1.0;
    beta  = -1.0;
    blobrad  = 0.1;
    maxgrid = 16;
    dofileio = 0;
    useDenseStencil = 0;
    resetDx();
  }

  int useDenseStencil;
  int relaxOnly;
  int maxiter;
  int numsmooth;
  int nstepmax;
  int nx;
  int maxgrid;
  int dofileio;
  double domsize;
  double tol;
  double dx;
  double alpha;
  double beta;
  double blobrad;
  double blobloc[DIM];
  
  void resetDx()
  {
    dx = domsize/nx;

    for(int idir = 0; idir < DIM; idir++)
    {
      blobloc[idir] = domsize/2.;
    }
  }

  void print() const
  {
    pout() << "Multigrid Solve of alpha I + beta Lapl phi = rhs  with periodic boundary conditions " << endl;

    pout() << "multigrid solve parameters: " << endl;
    pout() << "max iterations =  "   << maxiter   << endl;
    pout() << "num smooths    =  "   << numsmooth << endl;
    pout() << "nx             =  "   << nx        << endl;
    pout() << "max_grid       =  "   << maxgrid   << endl;
    pout() << "domain size    =  "   << domsize   << endl;
    pout() << "tolerance      =  "   << tol       << endl;
    pout() << "dx             =  "   << dx        << endl;
    pout() << "alpha          =  "   << alpha     << endl;
    pout() << "beta           =  "   << beta      << endl;
    pout() << "blobrad        =  "   << blobrad   << endl;
    pout() << "using multicolor gauss seidel smoothing  " << endl;
    if(relaxOnly == 1)
    {
      pout() << "doing relax only"  << endl;
    }
    else
    {
      pout() << "doing full multigrid solve"  << endl;
    }
    if(dofileio == 0)
    {
      pout() << "file output turned off" << endl;
    }
    else
    {
      pout() << "file output turned on" << endl;
    }
    if(useDenseStencil == 0)
    {
      pout() << "using standard sparse Laplacian stencil"  << endl;
    }
    else if(DIM==2)
    {
      pout() << "using dense (9 pt) Laplacian stencil"  << endl;
    }
    else 
    {
      pout() << "using dense (27 pt) Laplacian stencil"  << endl;
    }
  }
};                  
////////////
void
parseInputs(SolveParams & a_params)
{
  ParmParse pp;
  pp.get("nx"         , a_params.nx);
  pp.get("do_file_output", a_params.dofileio);
  pp.get("use_dense_stencil", a_params.useDenseStencil);
  pp.get("max_iter"   , a_params.maxiter);
  pp.get("relax_only" , a_params.relaxOnly);
  pp.get("num_smooth" , a_params.numsmooth);
  pp.get("domain_size", a_params.domsize);
  pp.get("tolerance"  , a_params.tol);
  pp.get("alpha"      , a_params.alpha);
  pp.get("beta"       , a_params.beta);
  pp.get("max_grid"   , a_params.maxgrid);
  pp.get("blob_rad"   , a_params.blobrad);
  a_params.resetDx();
  a_params.print();
}

/****************/
PROTO_KERNEL_START unsigned int setRHSF(Point&               a_p,
                                        Scalar            & a_rhs,
                                        SolveParams         a_params)
{
  double x[DIM];
  for(int idir = 0; idir < DIM; idir++)
  {    
    x[idir] = (a_p[idir] + 0.5)*a_params.dx;
  }

  double rad0sq = a_params.blobrad*a_params.blobrad;
  const double* const x0 = a_params.blobloc;
  double radsq = 0;
  for(int idir = 0; idir < DIM; idir++)
  {  
    radsq += (x[idir]-x0[idir])*(x[idir]- x0[idir]);
  }
  if(radsq < rad0sq)
  {
    double cosval = cos(0.5*PI*(radsq/rad0sq));
    //an attempt at having a rhs that sums to zero
    double rhsval =  (x[0]-x0[0])*cosval*cosval; 

    a_rhs(0) = rhsval;
  }
  else
  {
    a_rhs(0) = 0.0;
  }
  return 0;
}
PROTO_KERNEL_END(setRHSF, setRHS) 

PROTO_KERNEL_START void initParabolaT(Point& p, Scalar& data)
{
  data(0) = 0;
  for(int idir = 0; idir < DIM; idir ++)
  {
    data(0) += p[idir]*p[idir];
  }
}
PROTO_KERNEL_END(initParabolaT, initParabola)

/***/
double absMax(const LevelData<FArrayBox >& a_resid)
{
  double absmaxval = 0;
  DataIterator dit = a_resid.dataIterator();
#pragma omp parallel for
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    BoxData<double, 1> resbd;
    ProtoCh::aliasBoxData<double, 1>(resbd, a_resid[dit[ibox]]);
    double maxgrid = resbd.absMax();
#pragma omp critical
    absmaxval = std::max(maxgrid, absmaxval);
  }

#ifdef CH_MPI
  double recv;
  int result = MPI_Allreduce(&absmaxval, &recv, 1, MPI_DOUBLE,
                             MPI_MAX, Chombo_MPI::comm);
  if (result != MPI_SUCCESS)
  {
    MayDay::Error("sorry, but I had a communcation error on absmax");
  }
  absmaxval = recv;
#endif

  return absmaxval;
}
/***/
void
multigridSolve(const SolveParams& a_params)
{
  CH_TIME("Multigrid_Solve");
  SolveParams params = a_params;

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
  LevelData<FArrayBox > phi(grids, 1, IntVect::Unit);
  LevelData<FArrayBox > rhs(grids, 1, IntVect::Zero);
  LevelData<FArrayBox > res(grids, 1, IntVect::Zero);

  DataIterator dit = grids.dataIterator();
#pragma omp parallel
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    Bx patch = grids[dit[ibox]];
    BoxData<double, 1> rhsbd;
    ProtoCh::aliasBoxData<double, 1>(rhsbd, rhs[dit[ibox]]);
    forallInPlace_p(setRHS, patch, rhsbd, params);
    phi[dit[ibox]].setVal(0.);
  }

  
  SGMultigrid::s_numSmoothUp     =    a_params.numsmooth;
  SGMultigrid::s_numSmoothDown   =    a_params.numsmooth;
  SGMultigrid::s_numSmoothBottom = 16*a_params.numsmooth;

  int iter = 0;
  double rhsabsmax = absMax(rhs);
  double resStart = rhsabsmax;
  resStart = std::max(resStart, a_params.tol);
  double resIter  = resStart;
  pout() << "iter = " << iter << ", ||resid|| = " << resIter << endl;
  bool useDense = (a_params.useDenseStencil != 0);
  SGMultigrid solver(a_params.alpha, a_params.beta, a_params.dx, grids, useDense);
  if(a_params.relaxOnly == 0)
  {
    CH_TIME("full_multigrid_solve");
    while((resIter > a_params.tol*resStart) && (iter <  a_params.maxiter))
    {
      solver.vCycle(phi, rhs);
      solver.residual(res, phi, rhs);
  
      iter++;
      resIter = absMax(res);
      pout() << "iter = " << iter << ", ||resid|| = " << resIter << endl;
    }
  }
  else
  {
    CH_TIME("relaxation");
    for(int irelax = 0; irelax <a_params.maxiter ; irelax++)
    {
      pout() << "relax iter = " << iter << endl;
      solver.relax(phi, rhs);
      iter++;
    }
  }
#ifdef CH_USE_HDF5
  if(a_params.dofileio != 0)
  {
    CH_TIME("writing_files");
    pout() << "writing files" << endl;
    writeLevelname(&rhs, "rhs.hdf5");
    writeLevelname(&phi, "phi.hdf5");
  }
#endif
  
}
/**/
int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  // scoping...
  {
    if (argc < 2)
      {
        cerr<< " usage " << argv[0] << " <input_file_name> " << endl;
        exit(0);
      }
    char* in_file = argv[1];
    ParmParse  pp(argc-2,argv+2,NULL,in_file);

    SolveParams params;
    parseInputs(params);

    multigridSolve(params);
  }
  //end scoping trick
  
#ifdef CH_MPI
  dumpmemoryatexit();
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif

  return 0;

}  
