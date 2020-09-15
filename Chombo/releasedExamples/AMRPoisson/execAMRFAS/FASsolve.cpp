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
#include <fstream>

#include "LevelData.H"
#include "FArrayBox.H"
#include "ParmParse.H"
#include "Vector.H"
#include "AMRIO.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "ProblemDomain.H"
#include "BCFunc.H"
#include "AMRNonLinearPoissonOp.H"
#include "AMRPoissonOp.H"
#include "AMRFASMultiGrid.H"
#include "AMRMultiGrid.H"
#include "BiCGStabSolver.H"
#include "BoxIterator.H"
#include "CONSTANTS.H"
#include "memusage.H"
#include "computeNorm.H"
#include "FABView.H"

#include "UsingNamespace.H"

int s_verbosity = 1;

enum probTypes {exact,
  inexact,
  numProbTypes};


int s_probtype = exact;

//  -----------------------------------------
// boundary condition stuff
//  -----------------------------------------
///
/**
 */
class GlobalBCRS
{
public:
  static std::vector<bool> s_printedThatLo, s_printedThatHi;
  static std::vector<int> s_bcLo, s_bcHi;
  static RealVect s_trigvec;
  static bool s_areBCsParsed, s_valueParsed, s_trigParsed;
};

std::vector<bool> GlobalBCRS::s_printedThatLo = std::vector<bool>(SpaceDim, false);
std::vector<bool> GlobalBCRS::s_printedThatHi = std::vector<bool>(SpaceDim, false);
std::vector<int>  GlobalBCRS::s_bcLo = std::vector<int>();
std::vector<int>  GlobalBCRS::s_bcHi = std::vector<int>();
RealVect          GlobalBCRS::s_trigvec = RealVect::Zero;
bool              GlobalBCRS::s_areBCsParsed= false;
bool              GlobalBCRS::s_valueParsed= false;
bool              GlobalBCRS::s_trigParsed= false;

void ParseValue(Real* pos,
                int* dir,
                Side::LoHiSide* side,
                Real* a_values)
{
  a_values[0]=0.;
}

void ParseBC(FArrayBox& a_state,
             const Box& a_valid,
             const ProblemDomain& a_domain,
             Real a_dx,
             bool a_homogeneous)
{

  if (!a_domain.domainBox().contains(a_state.box()))
  {
    Box valid = a_valid;
    for (int i=0; i<CH_SPACEDIM; ++i)
    {
      // don't do anything if periodic
      if (!a_domain.isPeriodic(i))
      {
        Box ghostBoxLo = adjCellBox(valid, i, Side::Lo, 1);
        Box ghostBoxHi = adjCellBox(valid, i, Side::Hi, 1);
        if (!a_domain.domainBox().contains(ghostBoxLo))
        {
          DiriBC(a_state,
                 valid,
                 a_dx,
                 true,
                 ParseValue,
                 i,
                 Side::Lo,
                 1);
        }

        if (!a_domain.domainBox().contains(ghostBoxHi))
        {
          DiriBC(a_state,
                 valid,
                 a_dx,
                 true,
                 ParseValue,
                 i,
                 Side::Hi,
                 1);
        }
      } // end if is not periodic in ith direction
    }
  }
}

void setRHS(Vector<LevelData<FArrayBox>* > a_rhs,
            Vector<ProblemDomain>& a_amrDomains,
            Vector<int>& a_refRatios,
            Vector<Real>& a_amrDx,
            int a_finestLevel)
{
  CH_TIME("setRHS");

  Real gamma = 0;
  ParmParse pp("solver");
  pp.query("gamma", gamma);

  for (int lev=0; lev<=a_finestLevel; lev++)
  {
    LevelData<FArrayBox>& levelRhs = *(a_rhs[lev]);
    const DisjointBoxLayout& levelGrids = levelRhs.getBoxes();

    // rhs is cell-centered...
    RealVect ccOffset = 0.5*a_amrDx[lev]*RealVect::Unit;

    DataIterator levelDit = levelGrids.dataIterator();
    for (levelDit.begin(); levelDit.ok(); ++levelDit)
    {
      FArrayBox& thisRhs = levelRhs[levelDit];

      if (s_probtype == exact)
      {


        BoxIterator bit(thisRhs.box());
        for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect iv = bit();
          RealVect loc(iv);
          loc *= a_amrDx[lev];
          loc += ccOffset;

          D_TERM(Real x = loc[0];, Real y = loc[1];, Real z = loc[2];)
          Real mult_arg = D_TERM((x-x*x),*(y-y*y),*(z-z*z));
          Real plus_arg = D_TERM((x-x*x),+(y-y*y),+(z-z*z));

          thisRhs(iv, 0) = ( 2*plus_arg + gamma*mult_arg*exp(mult_arg));

        }
      }
      else if (s_probtype == inexact)
      {

        BoxIterator bit(thisRhs.box());
        for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect iv = bit();
          RealVect loc(iv);
          loc *= a_amrDx[lev];
          loc += ccOffset;

          Real x=1;
          Real y=0;
          Real z=1;
          D_TERM(x = loc[0];, y = loc[1];,  z = loc[2];)
            thisRhs(iv, 0) = ((9*M_PI*M_PI + gamma*exp((x*x-x*x*x)*sin(3*M_PI*(y*z))))*(x*x-x*x*x) +6*x -2)*sin(3*M_PI*y);

        }
      }
      else
      {
        MayDay::Error("undefined problem type");
      }
    } // end loop over grids on this level
  } // end loop over levels
}


void setExact(Vector<LevelData<FArrayBox>* > a_rhs,
              Vector<ProblemDomain>& a_amrDomains,
              Vector<int>& a_refRatios,
              Vector<Real>& a_amrDx,
              int a_finestLevel)
{
  CH_TIME("setExact");

  Real gamma = 0;
  ParmParse pp("solver");
  pp.query("gamma", gamma);

  for (int lev=0; lev<=a_finestLevel; lev++)
  {
    LevelData<FArrayBox>& levelRhs = *(a_rhs[lev]);
    const DisjointBoxLayout& levelGrids = levelRhs.getBoxes();

    // rhs is cell-centered...
    RealVect ccOffset = 0.5*a_amrDx[lev]*RealVect::Unit;

    DataIterator levelDit = levelGrids.dataIterator();
    for (levelDit.begin(); levelDit.ok(); ++levelDit)
    {
      FArrayBox& thisRhs = levelRhs[levelDit];

      if (s_probtype == exact)
      {
        BoxIterator bit(thisRhs.box());
        for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect iv = bit();
          RealVect loc(iv);
          loc *= a_amrDx[lev];
          loc += ccOffset;


          D_TERM(Real x = loc[0];, Real y = loc[1];, Real z = loc[2];)
          Real arg = D_TERM((x-x*x),*(y-y*y),*(z-z*z));
          thisRhs(iv, 0) = arg;

        }
      }
      else if (s_probtype == inexact)
      {
        BoxIterator bit(thisRhs.box());
        for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect iv = bit();
          RealVect loc(iv);
          loc *= a_amrDx[lev];
          loc += ccOffset;

          Real x=1;
          Real y=0;
          Real z=1;
          D_TERM(x = loc[0];, y = loc[1];, z = loc[2];)
          thisRhs(iv, 0) = (x*x - x*x*x)*sin(3*M_PI*y*z);

        }
      }
      else
      {

        //MayDay::Error("undefined problem type");
      }
    } // end loop over grids on this level
  } // end loop over levels
}



void
setupGrids(Vector<DisjointBoxLayout>& a_amrGrids,
           Vector<ProblemDomain>& a_amrDomains,
           Vector<int>& a_refRatios,
           Vector<Real>& a_amrDx,
           int& a_finestLevel)
{
  CH_TIME("setupGrids");

  a_finestLevel = 0;
  ParmParse ppGrids("grids");

  // get grid generation parameters
  int maxLevel, maxBoxSize, blockFactor;
  Real fillRatio;

  ppGrids.get("max_level", maxLevel);
  ppGrids.get("max_box_size",maxBoxSize);
  ppGrids.get("block_factor", blockFactor);
  ppGrids.get("fillRatio", fillRatio);

  // note that there only need to be numLevels-1 refinement ratios
  a_refRatios.resize(maxLevel);
  ppGrids.getarr("ref_ratio", a_refRatios, 0, maxLevel);

  Vector<int>  is_periodic_int;
  bool is_periodic[SpaceDim];
  ppGrids.getarr("is_periodic", is_periodic_int, 0, SpaceDim);
  for (int dir=0; dir<SpaceDim; dir++)
  {
    is_periodic[dir] = (is_periodic_int[dir] == 1);
  }

  IntVect numCells;
  Vector<int> incells(SpaceDim);
  ppGrids.getarr("num_cells", incells, 0, SpaceDim);
  numCells = IntVect(D_DECL6(incells[0],incells[1],incells[2],
                             incells[3],incells[4],incells[5]) );

  RealVect domainSize = RealVect::Unit;
  if (ppGrids.contains("domain_size"))
  {
    Vector<Real> insize(SpaceDim);
    ppGrids.getarr("domain_size", insize, 0, SpaceDim);
    domainSize = RealVect(D_DECL6(insize[0],insize[1],insize[2],
                                  insize[3],insize[4],insize[5]) );
  }

  // resize dataholders
  int maxNumLevels = maxLevel +1;
  a_amrGrids.resize(maxNumLevels);
  a_amrDomains.resize(maxNumLevels);
  a_amrDx.resize(maxNumLevels,-1);
  a_finestLevel = 0;

  // assumes dx=dy=dz
  a_amrDx[0] = domainSize[0]/numCells[0];

  IntVect domLo = IntVect::Zero;
  IntVect domHi  = numCells - IntVect::Unit;

  ProblemDomain baseDomain(domLo, domHi, is_periodic);
  a_amrDomains[0] = baseDomain;

  // set up refined domains, etc
  for (int lev=1; lev<= maxLevel; lev++)
  {
    a_amrDomains[lev] = a_amrDomains[lev-1];
    a_amrDomains[lev].refine(a_refRatios[lev-1]);
    a_amrDx[lev] = a_amrDx[lev-1]/a_refRatios[lev-1];
  }

  Vector<Vector<Box> > vectBoxes(maxLevel+1);

  // local scope. for base-level grid generation
  {
    CH_TIME("BaseGridCreation");
    // generate base level grids

    domainSplit(baseDomain, vectBoxes[0], maxBoxSize, blockFactor);

    Vector<int> procAssign(vectBoxes[0].size(), 0);

    LoadBalance(procAssign, vectBoxes[0]);

    DisjointBoxLayout baseGrids(vectBoxes[0], procAssign, baseDomain);

    a_amrGrids[0] = baseGrids;
  }


  if (maxLevel > 0)
  {
    bool read_grids = false;
    ppGrids.query("read_in_grids", read_grids);
    if (read_grids)
    {
      for (int ilev = 1; ilev <= maxLevel; ilev++)
      {
        const ProblemDomain& levDomain = a_amrDomains[ilev];

        Vector<Box>   boxes;
        char boxCountVar[100];
        int boxCount;
        sprintf(boxCountVar, "level_%d_box_count", ilev);
        ppGrids.get(boxCountVar, boxCount);
        boxes.resize(boxCount);
        for (int ibox = 0; ibox < boxCount; ibox++)
        {
          char boxLoVar[100];
          char boxHiVar[100];
          sprintf(boxLoVar, "level_%d_box_%d_lo", ilev, ibox);
          sprintf(boxHiVar, "level_%d_box_%d_hi", ilev, ibox);
          Vector<int> boxLo, boxHi;
          ppGrids.getarr(boxLoVar, boxLo, 0, SpaceDim);
          ppGrids.getarr(boxHiVar, boxHi, 0, SpaceDim);
          IntVect ivLo(D_DECL(boxLo[0], boxLo[1], boxLo[2]));
          IntVect ivHi(D_DECL(boxHi[0], boxHi[1], boxHi[2]));
          boxes[ibox] = Box(ivLo, ivHi);
          if (!levDomain.contains(boxes[ibox]))
          {
            MayDay::Error("box outside of domain");
          }
        }
        //check to see if level 0 domain is covered
        if (ilev == 0)
        {
          IntVectSet ivDom(levDomain.domainBox());
          for (int ibox = 0; ibox < boxes.size(); ibox++)
          {
            ivDom -= boxes[ibox];
          }
          if (!ivDom.isEmpty())
          {
            MayDay::Error("level 0 boxes must cover the domain");
          }
        }
        Vector<int>  proc(boxes.size());
        LoadBalance(proc,boxes);
        a_amrGrids[ilev] = DisjointBoxLayout(boxes, proc, levDomain);
        a_finestLevel++;
      }

    }
    else
    {
      // tag on grad(rhs)
      int bufferSize = 2;
      ppGrids.query("buffer_size", bufferSize);
      BRMeshRefine meshGen(a_amrDomains[0],
                           a_refRatios,
                           fillRatio,
                           blockFactor,
                           bufferSize,
                           maxBoxSize);

      // to be used by MeshRefine...
      Vector<Vector<Box> > oldMeshes(maxLevel+1);
      oldMeshes[0] = vectBoxes[0];
      for (int lev=1; lev<oldMeshes.size(); lev++)
      {
        oldMeshes[lev].push_back(a_amrDomains[lev].domainBox());
      }

      Real refineThresh;
      ppGrids.get("refine_threshold", refineThresh);

      Real threshSqr = refineThresh*refineThresh;

      bool moreLevels = true;
      while (moreLevels)
      {
        // tag based on grad(rhs)
        // first need to allocate RHS
        Vector<LevelData<FArrayBox>* > tempRHS(a_finestLevel+1, NULL);
        for (int lev=0; lev<= a_finestLevel; lev++)
        {
          // note that we add a ghost cell to simplify gradients
          tempRHS[lev] = new LevelData<FArrayBox>(a_amrGrids[lev],
                                                  1, IntVect::Unit);
        }

        setRHS(tempRHS, a_amrDomains, a_refRatios, a_amrDx,
               a_finestLevel);

        Vector<IntVectSet> tags(a_finestLevel+1);

        for (int lev=0; lev<a_finestLevel+1; lev++)
        {
          const DisjointBoxLayout& levelGrids = a_amrGrids[lev];
          const LevelData<FArrayBox>& levelRHS = *tempRHS[lev];
          IntVectSet& levelTags = tags[lev];

          // compute mag(gradient)
          DataIterator dit = levelGrids.dataIterator();
          for (dit.begin(); dit.ok(); ++dit)
          {
            const FArrayBox& rhsFab = levelRHS[dit];
            // local storage for gradient
            FArrayBox gradFab(levelGrids[dit],1);
            gradFab.setVal(0.0);
            Real thisGrad;

            BoxIterator bit(levelGrids[dit]);
            for (bit.begin(); bit.ok(); ++bit)
            {
              IntVect iv=bit();
              for (int dir=0; dir<SpaceDim; dir++)
              {
                // use mag(undivided gradient)
                IntVect hi = iv + BASISV(dir);
                IntVect lo = iv - BASISV(dir);
                thisGrad = rhsFab(hi,0) - rhsFab(lo,0);
                gradFab(iv,0) += (thisGrad*thisGrad);
              } // end loop over directions
            } // end loop over cells

            //gradFab now has mag(grad*dx)^2

            // tag where mag(gradient) > tolerance^Dim
            for (bit.begin(); bit.ok(); ++bit)
            {
              IntVect iv = bit();
              if (gradFab(iv,0) > threshSqr)
              {
                levelTags |= iv;
              }
            } // end loop over cells
          } // end loop over grids on this level

        } // end loop over levels


        // call meshRefine.
        for (int lev=1; lev<=a_finestLevel; lev++)
        {
          oldMeshes[lev] = vectBoxes[lev];
        }

        int topLevel = a_finestLevel;
        int newFinestLevel =  meshGen.regrid(vectBoxes,
                                             tags,
                                             0,
                                             topLevel,
                                             oldMeshes);


        // define new grids if necessary and test to see if we're done
        if (newFinestLevel > a_finestLevel)
        {
          a_finestLevel = newFinestLevel;

          // setup new grid hierarchy
          for (int lev=1; lev<=a_finestLevel; lev++)
          {
            Vector<int> procAssign(vectBoxes[lev].size(),0);
            LoadBalance(procAssign, vectBoxes[lev]);
            DisjointBoxLayout levelGrids(vectBoxes[lev],
                                         procAssign,
                                         a_amrDomains[lev]);
            a_amrGrids[lev] = levelGrids;
          }
          if (s_verbosity>2) pout() << "setupGrids: "<< a_finestLevel <<") size " << a_amrGrids[a_finestLevel].size() << endl;
        }
        else
        {
          moreLevels = false;
        }

        if (a_finestLevel == maxLevel)
        {
          moreLevels = false;
        }

        // clean up before starting again
        for (int lev=0; lev<tempRHS.size(); lev++)
        {
          delete tempRHS[lev];
        }

      } // end while (moreLevels)

    }

    // fill in remaining levels with empty DisjointBoxLayouts
    for (int lev= a_finestLevel+1; lev<=maxLevel; lev++)
    {
      a_amrGrids[lev] = DisjointBoxLayout();
    }

  }


}


void
setupSolver(AMRMultiGrid<LevelData<FArrayBox> > *a_amrSolver,
            LinearSolver<LevelData<FArrayBox> >& a_bottomSolver,
            const Vector<DisjointBoxLayout>& a_amrGrids,
            const Vector<ProblemDomain>& a_amrDomains,
            const Vector<int>& a_refRatios,
            const Vector<Real>& a_amrDx,
            int a_finestLevel)
{
  CH_TIME("setupSolver");

  ParmParse ppSolver("solver");

  bool nonlinearOp = true;
  ppSolver.query("nonlinearOp", nonlinearOp);

  int numLevels = a_finestLevel+1;

  if (nonlinearOp)
  {
  AMRNonLinearPoissonOpFactory opFactory;

  // solving nonlinear poisson problem here
  Real alpha = 0.0;
  Real beta = 1.0;
  Real gamma = -1.0;

  ppSolver.query("gamma", gamma);

  opFactory.define(a_amrDomains[0],
                   a_amrGrids,
                   a_refRatios,
                   a_amrDx[0],
                   &ParseBC, alpha, beta,gamma);


  AMRLevelOpFactory<LevelData<FArrayBox> >& castFact = (AMRLevelOpFactory<LevelData<FArrayBox> >& ) opFactory;

  a_amrSolver->define(a_amrDomains[0], castFact,
                      &a_bottomSolver, numLevels);

  }
  else
  {
        AMRPoissonOpFactory opFactory;

      // solving poisson problem here
      Real alpha = 0.0;
      Real beta = -1.0;

      opFactory.define(a_amrDomains[0],
                       a_amrGrids,
                       a_refRatios,
                       a_amrDx[0],
                       &ParseBC, alpha, beta);

      AMRLevelOpFactory<LevelData<FArrayBox> >& castFact = (AMRLevelOpFactory<LevelData<FArrayBox> >& ) opFactory;

      a_amrSolver->define(a_amrDomains[0], castFact,
                          &a_bottomSolver, numLevels);
  }

  // multigrid solver parameters
  int numSmooth, numMG, maxIter;
  int numBottom;
  ppSolver.query("num_bottom", numBottom);
  Real eps, hang;
  ppSolver.get("num_smooth", numSmooth);
  ppSolver.get("num_mg",     numMG);
  ppSolver.get("max_iterations", maxIter);
  ppSolver.get("tolerance", eps);
  ppSolver.get("hang",      hang);

  Real normThresh = 1.0e-30;
  a_amrSolver->setSolverParameters(numSmooth, numSmooth, numBottom,
                                   numMG, maxIter, eps, hang, normThresh);
  a_amrSolver->m_verbosity = s_verbosity-1;

  // optional parameters
  ppSolver.query("num_pre", a_amrSolver->m_pre);
  ppSolver.query("num_post", a_amrSolver->m_post);

}

int runSolver()
{
  CH_TIME("runSolver");

  int status = 0;
  ParmParse ppMain("main");
  ParmParse ppSolver("solver");

  ppMain.query("verbosity", s_verbosity);
  ppMain.query("problem", s_probtype);

  // set up grids&
  Vector<DisjointBoxLayout> amrGrids;
  Vector<ProblemDomain> amrDomains;
  Vector<int> refRatios;
  Vector<Real> amrDx;
  int finestLevel;

  setupGrids(amrGrids, amrDomains, refRatios, amrDx, finestLevel);

  // initialize solver
  bool FASmultigrid = true;
  ppSolver.query("FASmultigrid", FASmultigrid);

  AMRMultiGrid<LevelData<FArrayBox> > *amrSolver;

  if (FASmultigrid)
  {
    amrSolver = new AMRFASMultiGrid<LevelData<FArrayBox> >();
  }
  else
  {
    amrSolver = new AMRMultiGrid<LevelData<FArrayBox> >();
  }


  BiCGStabSolver<LevelData<FArrayBox> > bottomSolver;
  bottomSolver.m_verbosity = s_verbosity-2;
  setupSolver(amrSolver, bottomSolver, amrGrids, amrDomains,
              refRatios, amrDx, finestLevel);


  // allocate solution and RHS, initialize RHS
  int numLevels = amrGrids.size();
  Vector<LevelData<FArrayBox>* > phi(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > rhs(numLevels, NULL);
  // this is for convenience
  Vector<LevelData<FArrayBox>* > resid(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > exact(numLevels, NULL);
  Vector<LevelData<FArrayBox>* > error(numLevels, NULL);

  for (int lev=0; lev<=finestLevel; lev++)
  {
    const DisjointBoxLayout& levelGrids = amrGrids[lev];
    phi[lev] = new LevelData<FArrayBox>(levelGrids, 1, IntVect::Unit);
    rhs[lev] = new LevelData<FArrayBox>(levelGrids, 1, IntVect::Zero);
    resid[lev] = new LevelData<FArrayBox>(levelGrids, 1, IntVect::Zero);
    exact[lev] = new LevelData<FArrayBox>(levelGrids, 1, IntVect::Zero);
    error[lev] = new LevelData<FArrayBox>(levelGrids, 1, IntVect::Zero);
  }

  bool zeroInitialGuess = true;
  setRHS(rhs, amrDomains, refRatios, amrDx, finestLevel );
  setExact(exact, amrDomains, refRatios, amrDx, finestLevel );

  // Start with exact solution
  zeroInitialGuess = false;
  setExact(phi, amrDomains, refRatios, amrDx, finestLevel );

  // do solve
  int iterations = 1;
  ppMain.get("iterations", iterations);
  ppSolver.query("zeroInitialGuess", zeroInitialGuess);

  for (int iiter = 0; iiter < iterations; iiter++)
  {
    pout() << "about to go into solve" << endl;
    amrSolver->solve(phi, rhs, finestLevel, 0, zeroInitialGuess);
    pout() << "done solve" << endl;
  }

  // Compute error
  for (int lev=0; lev<=finestLevel; lev++)
  {
    phi[lev]->copyTo(*error[lev]);
    for (DataIterator dit = error[lev]->dataIterator(); dit.ok(); ++dit)
    {
      (*error[lev])[dit].minus((*exact[lev])[dit]);
    }
  }

  // Compute error metrics
  Real max, L1, L2;
  max = computeNorm(error,    refRatios, amrDx[0], Interval(0,0), 0);
  L1 = computeNorm(error,    refRatios, amrDx[0], Interval(0,0), 1);
  L2 = computeNorm(error,    refRatios, amrDx[0], Interval(0,0), 2);

  char errStr[1024];
  sprintf(errStr, "Error = %1.2e (max), %1.2e (L1), %1.2e (L2)", max, L1, L2 );
  pout() << errStr << endl;

  // write results to file

  bool writePlots = true;
  ppMain.query("writePlotFiles", writePlots);

  int numPlotComps = 5;
#ifdef CH_USE_HDF5

  if (writePlots)
  {
    int numLevels = finestLevel +1;
    Vector<LevelData<FArrayBox>* > plotData(numLevels, NULL);

    pout() << "Write Plots. norm=" << amrSolver->computeAMRResidual(resid,phi,rhs,finestLevel,0) << endl;

    for (int lev=0; lev<numLevels; lev++)
    {
      plotData[lev] = new LevelData<FArrayBox>(amrGrids[lev],
                                               numPlotComps, IntVect::Zero);

      Interval phiInterval(0,0);
      phi[lev]->copyTo(phiInterval, *plotData[lev], phiInterval);
      Interval rhsInterval(1,1);
      rhs[lev]->copyTo(phiInterval, *plotData[lev], rhsInterval);
      Interval resInterval(2,2);
      resid[lev]->copyTo(phiInterval, *plotData[lev], resInterval);
      Interval exactInterval(3,3);
      exact[lev]->copyTo(phiInterval, *plotData[lev], exactInterval);
      Interval errInterval(4,4);
      error[lev]->copyTo(phiInterval, *plotData[lev], errInterval);
    }

    string fname = "FASOut.";

    char suffix[30];
    sprintf(suffix, "%dd.hdf5",SpaceDim);
    fname += suffix;

    Vector<string> varNames(numPlotComps);
    varNames[0] = "phi";
    varNames[1] = "rhs";
    varNames[2] = "res";
    varNames[3] = "exact";
    varNames[4] = "error";

    Real bogusVal = 1.0;

    WriteAMRHierarchyHDF5(fname,
                          amrGrids,
                          plotData,
                          varNames,
                          amrDomains[0].domainBox(),
                          amrDx[0],
                          bogusVal,
                          bogusVal,
                          refRatios,
                          numLevels);

    // clean up
    for (int lev=0; lev<plotData.size(); lev++)
    {
      delete plotData[lev];
    }
  } // end if writing plots
#endif // end if HDF5

  // clean up
  for (int lev=0; lev<phi.size(); lev++)
  {
    delete phi[lev];
    delete rhs[lev];
    delete resid[lev];
    delete exact[lev];
    delete error[lev];
  }

  delete amrSolver;

  return status;
}


/*****/
int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  int status = 0;

  // scoping...
  {
    if (argc < 2)
    {
      cerr<< " usage " << argv[0] << " <input_file_name> " << endl;
      exit(0);
    }
    char* in_file = argv[1];
    ParmParse  pp(argc-2,argv+2,NULL,in_file);

    int solverStatus = runSolver();
    status += solverStatus;
  }
  //end scoping trick

#ifdef CH_MPI
  dumpmemoryatexit();
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif

  return(status);
}
