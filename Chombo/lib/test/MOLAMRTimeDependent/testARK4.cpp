#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


#include "parstream.H"
using std::endl;
#include "FArrayBox.H"
#include "LevelData.H"
#include "LevelDataOps.H"
#include "LayoutIterator.H"
#include "LoadBalance.H"
#include "BRMeshRefine.H"

// ARK4 test classes
#include "ARK4DenseOutput.H"
#include "TestRhsData.H"
#include "TestSolnData.H"
#include "TestImExOp.H"

#include "UsingNamespace.H"

/// Global variables for handling output:
static const char* pgmname = "testARK4" ;
static const char* indent = "   ";
static const char* indent2 = "      " ;

/// Prototypes:
int
testARK4();

void
parseTestOptions(int argc ,char* argv[]) ;

int
main(int argc ,char* argv[])
{
#ifdef CH_MPI
  MPI_Init (&argc, &argv);
#endif
  pout () << indent2 << "Beginning " << pgmname << " ..." << endl ;

  int status = testARK4();

  if ( status == 0 )
    pout() << indent << pgmname << " passed." << endl ;
  else
    pout() << indent << pgmname << " failed with return code " << status << endl ;

#ifdef CH_MPI
  MPI_Finalize ();
#endif
  return status ;
}

// Check that the data is a constant, and return the constant
Real ldfabVal(LevelData<FArrayBox>& a_data)
{
  Real ldmin;
  Real ldmax;
  ldmin = CH_BADVAL;
  ldmax = -CH_BADVAL;
  DataIterator dit(a_data.getBoxes());
  for (dit.begin(); dit.ok(); ++dit)
  {
    Real min = a_data[dit].min();
    ldmin = (ldmin < min) ? ldmin : min;
    Real max = a_data[dit].max();
    ldmax = (ldmax > max) ? ldmax : max;
  }
  CH_assert(ldmax == ldmin);
  return ldmax;
}

int
testARK4 ()
{
  CH_TIMERS("testARK4");

  int numCells = 32;
  IntVect loVect = IntVect::Zero;
  IntVect hiVect = (numCells-1)*IntVect::Unit;
  Box domainBox(loVect, hiVect);
  ProblemDomain baseDomain(domainBox);

  int maxBoxSize = numCells/2;

  Vector<Box> vectBoxes;
  domainSplit(baseDomain, vectBoxes, maxBoxSize, 1);
  Vector<int> procAssign(vectBoxes.size(), 0);
  LoadBalance(procAssign, vectBoxes);
  DisjointBoxLayout dbl(vectBoxes, procAssign, baseDomain);

  // Set up the data classes
  LevelData<FArrayBox> data(dbl,1,IntVect::Zero);
  LevelData<FArrayBox> accum(dbl,1,IntVect::Zero);
  TestSolnData soln;
  // soln.define(dbl,1,IntVect::Zero);
  soln.aliasData(data, &accum);

  int nDenseCoefs = 4;
  Vector<TestRhsData*> denseCoefs(nDenseCoefs);
  for (int icoef = 0; icoef < nDenseCoefs; ++icoef)
  {
    denseCoefs[icoef] = new TestRhsData();
    denseCoefs[icoef]->define(dbl,1,IntVect::Zero);
  }

  // Do a convergence study, across 3 time step sizes
  Real basedt = 1;
  Real Nres = 4;
  Vector<Real> errors(Nres,0);
  Vector<Real> denseErrs(Nres,0);

  // Exact solution is exp((cE + cI)*t)*(1+t)
  // So that:
  //   explicit op = (1/(1+t) + cE) * phi;
  //   implicit op = cI * phi;
  Real coef = TestImExOp::s_cE + TestImExOp::s_cI;

  bool denseOutput = true;
  ARK4DenseOutput<TestSolnData, TestRhsData, TestImExOp> ark;
  ark.define(soln, basedt, denseOutput); 
  ark.start0end1(true, true);
  LevelDataOps<FArrayBox> ops;
  for (int res=0; res < Nres; ++res)
  {
    Real time = 0;
    int Nstep = pow((Real) 2,(Real) res+2);
    Real dt = basedt / (Real) Nstep;
    pout() << "Time step: " << dt << endl;
    ark.resetDt(dt);
    // Set the initial condition
    Real phi0 = 1.0;
    ops.setVal(data, phi0);
    ops.setVal(accum, 0);
    // advance nstep
    for (int step = 0; step < Nstep; ++step)
    {
      ark.advance(time, soln);
      time += dt;
    }
    Real exact = exp(coef*time)*(1 + time);
    Real val = ldfabVal(data);
    Real error = (exact - val);
    pout() << "Soln at time " << time << " = " << 
      val << ", error = " << error << endl;
    errors[res] = error;

    Real accumDiff = ldfabVal(accum);
    pout() << "Accumulated RHS = " << accumDiff <<
      " , difference from soln change = " << (accumDiff - (val - phi0)) << endl;

    // Test dense output for the last time step
    Real theta = 1 - 1/sqrt(2);
    Real tint = time - (1-theta)*dt;
    exact = exp(coef*tint)*(1 + tint);
    ark.denseOutputCoefs(denseCoefs);
    // Compute time-interpolated value in soln
    soln.zero();
    soln.increment(*denseCoefs[0]);
    for (int icoef = 1; icoef < nDenseCoefs; ++icoef)
    {
      Real factor = pow(theta, icoef);
      soln.increment(*denseCoefs[icoef],factor);
    }
    val = ldfabVal(data);
    error = exact - val;
    denseErrs[res] = error;
    pout() << "Dense output at time " << tint << " = " << 
      val << ", error = " << error << endl;
  }

  pout() << "Orders of convergence: " << endl;
  Real rate = 4;
  for (int res=1; res < Nres; ++res)
  {
    Real ratio = errors[res] / errors[res-1];
    Real solnrate = (-log(ratio) / log(2));
    pout() << "  soln: " << solnrate << endl;
    ratio = denseErrs[res] / denseErrs[res-1];
    Real denserate = (-log(ratio) / log(2));
    pout() << "  dense output: " << denserate << endl;
    rate = min(solnrate, denserate);
  }

  return (rate > 3.8) ? 0 : 1;
}
