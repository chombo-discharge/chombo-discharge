#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "BRMeshRefine.H"
#include "RefCountedPtr.H"
#include "AMRMultiGrid.H"
#include "AMRPoissonOp.H"
#include "LoadBalance.H"

#include "AMRBoxesAndRanksIO.H"
#include "DiffusionSolver.H"
#include "DiffusionSolverF_F.H"

void DiffValue(Real* pos,
               int* dir,
               Side::LoHiSide* side,
               Real* a_values)
{
  a_values[0]=0.0;
}

void DiffBC(FArrayBox& a_state,
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
          Box ghostBoxLo = adjCellBox(valid, i, Side::Lo, 1);
          Box ghostBoxHi = adjCellBox(valid, i, Side::Hi, 1);
          if (!a_domain.domainBox().contains(ghostBoxLo))
            {
              NeumBC(a_state,
                     valid,
                     a_dx,
                     a_homogeneous,
                     DiffValue,
                     i,
                     Side::Lo);
            }

          if (!a_domain.domainBox().contains(ghostBoxHi))
            {
              NeumBC(a_state,
                     valid,
                     a_dx,
                     a_homogeneous,
                     DiffValue,
                     i,
                     Side::Hi);
            }
        }
    }
}

DiffusionSolver::DiffusionSolver(const DiffusionParams & a_params)
{
  // Save a copy of the solver parameters
  m_params = a_params;

  // The integrator hasn't been defined yet
  m_integrator = NULL;
}

DiffusionSolver::~DiffusionSolver()
{
  // Delete data holder for the solution at the old time (on all AMR levels)
  for (int ilev = 0; ilev < m_solnOld.size(); ilev++)
    {
      delete m_solnOld[ilev];
    }

  // Delete data holder for the solution at the new time (on all AMR levels)
  for (int ilev = 0; ilev < m_solnNew.size(); ilev++)
    {
      delete m_solnNew[ilev];
    }

  // Delete data holder for the source/sink term (on all AMR levels)
  for (int ilev = 0; ilev < m_source.size(); ilev++)
    {
      delete m_source[ilev];
    }

  // Delete the integrator (if defined)
  if (m_integrator != NULL)
    {
      delete m_integrator;
    }
}

void DiffusionSolver::init()
{
  CH_TIME("DiffusionSolver::init");

  // Initialize the index space
  initIndexSpace();

  // Initialize the data
  initData();

  // Initialize the solver
  initSolver();
}

void DiffusionSolver::run()
{
  CH_TIME("DiffusionSolver::run");

  // Compute the number of time steps
  int numSteps = m_params.m_endTime / m_params.m_dt;

  // Iterate until the end time is reached
  int step;
  for (step = 0; step < numSteps; step++)
    {
      // Set and print the current time
      Real time = step * m_params.m_dt;

      pout() << "time = " << time << "\n";

      // Write the solution if it's the right time
      if (m_params.m_outputInterval > 0 &&
          step % m_params.m_outputInterval == 0)
        {
          writeOutput(step,time);
        }

      // Set the source/sink term (combined here)
      setSource();

      // Advance one time step
      {
        CH_TIME("DiffusionSolver::advance");
        m_integrator->oneStep(m_solnNew,
                              m_solnOld,
                              m_source,
                              m_params.m_dt,
                              0,
                              m_params.m_numLevels-1,
                              false);
      }

      // Copy the new solution to the old solution
      {
        LevelDataOps<FArrayBox> ops;
        CH_TIME("DiffusionSolver::copy");
        for (int ilev = 0; ilev < m_params.m_numLevels; ilev++)
          {
            ops.assign(*m_solnOld[ilev],*m_solnNew[ilev]);
          }
      }
    }

  // Write the solution at the end
  if (m_params.m_outputInterval > 0)
    {
      writeOutput(step,m_params.m_endTime);
    }
}


void DiffusionSolver::initIndexSpace()
{
  Vector< Vector<Box> > boxes;
  Vector< Vector<int> > ranks;
  Box baseProbBox;
  readABRfile(boxes, ranks, m_params.m_refRatio,
              baseProbBox,
              m_params.m_abrFile);

  m_params.m_coarsestDomain = ProblemDomain(baseProbBox);
  //sanity
  if (m_params.m_refRatio.size() == 0) m_params.m_refRatio.resize(1,2);

  m_params.m_coarsestDomain = ProblemDomain(baseProbBox);
  m_grids.resize(m_params.m_numLevels);
  for (int ilev = 0; ilev < m_params.m_numLevels; ilev++)
    {
      Vector<int> procs;

      mortonOrdering(boxes[ilev]);

      LoadBalance(procs, boxes[ilev]);

      m_grids[ilev] = DisjointBoxLayout(boxes[ilev], procs);
    }
}

void DiffusionSolver::initData()
{
  CH_TIME("DiffusionSolver::initData");
  LevelDataOps<FArrayBox> ops;

  // Make pointers for data holders at each AMR level
  m_solnOld.resize(m_params.m_numLevels);
  m_solnNew.resize(m_params.m_numLevels);

  m_source.resize(m_params.m_numLevels);

  // The simulation only has one variable
  int oneVariable = 1;

  for (int ilev = 0; ilev < m_params.m_numLevels; ilev++)
    {
      // At each level, allocate and initialize (as needed) space for the old
      // solution, the new solution, and the source/sink terms
      m_solnOld[ilev] = new LevelData<FArrayBox>(m_grids[ilev],
                                                 oneVariable,
                                                 m_params.m_numGhostSoln);

      ops.setVal(*(m_solnOld[ilev]),m_params.m_initialValue);

      m_solnNew[ilev] = new LevelData<FArrayBox>(m_grids[ilev],
                                                 oneVariable,
                                                 m_params.m_numGhostSoln);

      m_source[ilev] = new LevelData<FArrayBox>(m_grids[ilev],
                                                oneVariable,
                                                m_params.m_numGhostSource);

    }

  for (int ilev = 0; ilev < m_params.m_numLevels; ilev++)
    {
      // At each level, count and report the
      // number of computational cells
      Real totalFortran = 0.0;


      DisjointBoxLayout     fineGrids   = m_grids[ilev];

      for (DataIterator dit=fineGrids.dataIterator(); dit.ok(); ++dit)
        {
          const Box&        curBox      = fineGrids[dit()];
          totalFortran += curBox.numPts();
        }

#ifdef CH_MPI
      Real mpiTotal;
      int status = MPI_Allreduce(&totalFortran, &mpiTotal, 1, MPI_CH_REAL,
                                 MPI_SUM, Chombo_MPI::comm);

      if (status != MPI_SUCCESS)
        {
          MayDay::Error("MPI error summing 'totalFortran'");
        }

      totalFortran = mpiTotal;

#endif

      long long totalFortranInt = totalFortran;

      std::ios::fmtflags origFlags = pout().flags();
      int origWidth = pout().width();
      int origPrecision = pout().precision();

      pout() << "Level " << ilev << ":" << "\n";
      pout() << setiosflags(ios::right);
      pout() << "  Total computation cells: " << setw(10) << totalFortranInt << "\n";
      pout() << "\n";

      pout().flags(origFlags);
      pout().width(origWidth);
      pout().precision(origPrecision);
    }
}

void DiffusionSolver::initSolver()
{
  CH_TIME("DiffusionSolver::initSolver");

  // This is the multigrid solver used for backward Euler
  RefCountedPtr<AMRMultiGrid<LevelData<FArrayBox > > >
    solver(new AMRMultiGrid<LevelData<FArrayBox> >);

  // Define the data holders and interpolators
  RealVect vectDx = m_params.m_dx * RealVect::Unit;

  // Minimize the ghost cell filling in the operator relaxation
  if (m_params.m_mgLazyRelax)
    {
      AMRPoissonOp::s_relaxMode = 0;
    }
  else
    {
      AMRPoissonOp::s_relaxMode = 1;
    }

  // Set up a factory to produce Helmholtz operators at various resolutions
  AMRPoissonOpFactory operatorFactory;
  operatorFactory.define(m_params.m_coarsestDomain,
                         m_grids,
                         m_params.m_refRatio,
                         m_params.m_dx,
                         DiffBC,
                         1.0, m_params.m_diffusionConstant);

  // Set the verbosity of the bottom solver for multigrid
  m_bottomSolver.m_verbosity = 0;

  // Define the multigrid solver and set various parameters
  solver->define(m_params.m_coarsestDomain,
                 operatorFactory,
                 &m_bottomSolver,
                 m_params.m_numLevels);

  Real normThresh = 1.0e-30;

  solver->setSolverParameters(m_params.m_mgNumSmooths,
                              m_params.m_mgNumSmooths,
                              m_params.m_mgNumSmooths,
                              m_params.m_mgNumCycles,
                              m_params.m_mgIterMax,
                              m_params.m_mgToler,
                              m_params.m_mgHangToler,
                              normThresh);
  solver->m_verbosity = 3;
  solver->init(m_solnOld,m_source,m_params.m_numLevels-1,0);

  // Create the backward Euler solver based on the multigrid solver
  m_integrator = new BackwardEuler(solver,
                                   operatorFactory,
                                   m_params.m_coarsestDomain,
                                   m_params.m_refRatio,
                                   m_params.m_numLevels);
}

void DiffusionSolver::setSource()
{
  CH_TIME("DiffusionSolver::setSource");

  Real dx = m_params.m_dx;

  // Set the source/sink term on each level
  for (int ilev = 0; ilev < m_params.m_numLevels; ilev++)
    {
      // Grids on this level
      DisjointBoxLayout& curGrids = m_grids[ilev];

      // Data holders for the source/sink and old solution on this level
      LevelData<FArrayBox>* sourceLevel = m_source[ilev];
      LevelData<FArrayBox>* solnLevel   = m_solnOld[ilev];

      // Set the source/sink term for each grid
      for (DataIterator dit=curGrids.dataIterator(); dit.ok(); ++dit)
        {
          // Current grid
          const Box& curBox = curGrids[dit()];

          // Single grid data holders
          FArrayBox& sourceFAB = (*sourceLevel)[dit()];
          FArrayBox& solnFAB   =   (*solnLevel)[dit()];

          // Call Fortran to set the source/sink for all the single valued cells
          FORT_SIMPLESOURCESINK(CHF_FRA1(sourceFAB,0),
                                CHF_CONST_FRA1(solnFAB,0),
                                CHF_CONST_REAL(m_params.m_sourceScaling),
                                CHF_CONST_REAL(m_params.m_sinkScaling),
                                CHF_CONST_REAL(dx),
                                CHF_CONST_REALVECT(m_params.m_loCorner),
                                CHF_BOX(curBox));

        }

      if (ilev < m_params.m_refRatio.size())
        dx /= m_params.m_refRatio[ilev];
    }
}

void DiffusionSolver::writeOutput(const int  & a_step,
                                  const Real & a_time)
{
  CH_TIME("DiffusionSolver::writeOutput");

  // Generate the full output file name
  char suffix[128];
  sprintf(suffix,".%06d.%dd.hdf5",a_step,SpaceDim);

  string filename = m_params.m_outputPrefix + suffix;

  // Name the solution variable
  Vector<string> names(1);
  names[0] = "G";

  // Write the output file
  WriteAMRHierarchyHDF5(filename,
                        m_grids,
                        m_solnOld,
                        m_params.m_coarsestDomain.domainBox(),
                        m_params.m_refRatio,
                        m_params.m_numLevels);
}
