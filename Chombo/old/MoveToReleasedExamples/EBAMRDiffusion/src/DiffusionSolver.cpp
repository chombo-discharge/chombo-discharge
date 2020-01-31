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

#include "EBIndexSpace.H"
#include "EBEllipticLoadBalance.H"
#include "EBLevelDataOps.H"
#include "EBAMRDataOps.H"
#include "EBAMRIO.H"
#include "EBLevelGrid.H"
#include "EBQuadCFInterp.H"
#include "EBAMRPoissonOpFactory.H"
#include "BaseDomainBC.H"
#include "NeumannPoissonDomainBC.H"
#include "NeumannPoissonEBBC.H"
#include "NeumannConductivityDomainBC.H"
#include "NeumannConductivityEBBC.H"
#include "GeometryShop.H"
#include "AllRegularService.H"

#include "DiffusionSolver.H"
#include "SimpleIF.H"
#include "KeratocyteIF.H"
#include "NewKeratocyteIF.H"
#include "AMRBoxesAndRanksIO.H"
#include "ParmParse.H"
#include "DiffusionSolverF_F.H"

DiffusionSolver::DiffusionSolver(const DiffusionParams & a_params)
{
  // Save a copy of the solver parameters
  m_params = a_params;
  m_time = 0;

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

  // Initialize the geometry
  initGeometry();

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
    m_time = time;
    pout() << "time = " << time << "\n";

    // If coefficients are changing, need to reinitialize solver.
    if (m_params.m_useVariableCoeff && (step > 0))
    {
      initSolver();
    }

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
      CH_TIME("DiffusionSolver::copy");
      EBAMRDataOps::assign(m_solnOld,m_solnNew);
    }
  }

  // Write the solution at the end
  if (m_params.m_outputInterval > 0)
  {
    writeOutput(step,m_params.m_endTime);
  }
}

void DiffusionSolver::initGeometry()
{
  CH_TIME("DiffusionSolver::initGeometry");

  // Compute the finest domain and cell size
  ProblemDomain fineDomain = m_params.m_coarsestDomain;
  RealVect fineDx = m_params.m_dx * RealVect::Unit;

  for (int ilev = 0; ilev < m_params.m_numLevels-1; ilev++)
  {
    fineDomain.refine(m_params.m_refRatio[ilev]);
    fineDx /= m_params.m_refRatio[ilev];
  }

  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  int maxCoarsen = -1;

  if (m_params.m_geometry == "noEB")
  {
    // If there is no EB, no implicit function is needed
    AllRegularService allReg;

    ebisPtr->define(fineDomain,
                    m_params.m_loCorner,
                    fineDx[0],
                    allReg,
                    m_params.m_maxBoxSize,
                    maxCoarsen);
  }
  else
  {
    // Create the implicit function for one of the known geometries
    BaseIF* geometry;

    if (m_params.m_geometry == "simple")
    {
      geometry = new SimpleIF();
    }
    else if (m_params.m_geometry == "keratocyte")
    {
      geometry = new KeratocyteIF();
    }
    else if (m_params.m_geometry == "new_keratocyte")
    {
      geometry = new NewKeratocyteIF();
    }
    else
    {
      string errorMessage = "DiffusionSolver::initGeometry:  Unknown geometry: '" + m_params.m_geometry + "'";
      MayDay::Error(errorMessage.c_str());
    }

    // Set up the infrastructure for the geometry
    int verbosity = 0;
    GeometryShop workshop(*geometry,verbosity,fineDx);

    ebisPtr->define(fineDomain,
                    m_params.m_loCorner,
                    fineDx[0],
                    workshop,
                    m_params.m_maxBoxSize,
                    maxCoarsen);

    delete geometry;
  }
}

void DiffusionSolver::initIndexSpace()
{
  CH_TIME("DiffusionSolver::initIndexSpace");

  // Create a set of grids at the coarsest level
  Vector<int> coarsestProcs;
  Vector<Box> coarsestBoxes;

  domainSplit(m_params.m_coarsestDomain,
              coarsestBoxes,
              m_params.m_maxBoxSize,
              m_params.m_blockFactor);

  // Order them using a space filling curve
  mortonOrdering(coarsestBoxes);

  // Do load balancing using timing information
  EBEllipticLoadBalance(coarsestProcs,
                        coarsestBoxes,
                        m_params.m_coarsestDomain);

  // Create more of the infrastructure for geometry
  m_grids.resize(m_params.m_numLevels);
  m_grids[0] = DisjointBoxLayout(coarsestBoxes,
                                 coarsestProcs,
                                 m_params.m_coarsestDomain);

  m_ebisl.resize(m_params.m_numLevels);

  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  ebisPtr->fillEBISLayout(m_ebisl[0],
                          m_grids[0],
                          m_params.m_coarsestDomain,
                          m_params.m_numGhostEBISLayout);

  if (m_params.m_numLevels > 1)
  {
    // If there is more than one level, the finer levels need to created by
    // "BRMeshRefine"
    BRMeshRefine meshRefine;

    meshRefine.define(m_params.m_coarsestDomain,
                      m_params.m_refRatio,
                      m_params.m_fillRatio,
                      m_params.m_blockFactor,
                      m_params.m_nestingRadius,
                      m_params.m_maxBoxSize);

    // Compute the second finest domain
    ProblemDomain secondFinestDomain = m_params.m_coarsestDomain;

    for (int ilev = 0; ilev < m_params.m_numLevels-2; ilev++)
    {
      secondFinestDomain.refine(m_params.m_refRatio[ilev]);
    }

    // Tags for creating the finer levels
    Vector<IntVectSet> tags(m_params.m_numLevels);

    // Get the depth of the second to finest level
    EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
    int depth = ebisPtr->getLevel(secondFinestDomain);

    if (m_params.m_tagType == "boundary")
    {
      // Get the indices of all the irregular cells at the second to finest
      // level
      tags[m_params.m_numLevels-2] = ebisPtr->irregCells(depth);
      tags[m_params.m_numLevels-2].grow(2);
    }
    else if (m_params.m_tagType == "interior")
    {
      // Get the boxes covering the interior and irregular cells on the second
      // to finest level
      DisjointBoxLayout curGrids = ebisPtr->getFlowGrids(depth);

      // Get the indices of all the cells associated with the above boxes
      for (DataIterator dit=curGrids.dataIterator(); dit.ok(); ++dit)
      {
        tags[m_params.m_numLevels-2] |= curGrids[dit()];
      }
    }
    else if (m_params.m_tagType == "sink")
    {
      double fineDx = m_params.m_dx;
      for (int i = 0; i < m_params.m_numLevels - 2; i ++)
      {
        fineDx /= m_params.m_refRatio[i];
      }

      DisjointBoxLayout curGrids = ebisPtr->getFlowGrids(depth);

      for (DataIterator dit = curGrids.dataIterator(); dit.ok(); ++dit)
      {
        const Box& box = curGrids.get(dit());
        for (BoxIterator bit(box); bit.ok(); ++bit)
        {
          const IntVect& iv = bit();

          RealVect x = (iv + 0.5*RealVect::Unit) * fineDx + m_params.m_loCorner;

          if ((x[0] > -2.0) && (x[0] * x[0] / (4.7 * 4.7) + x[1] * x[1] / (19.7 * 19.7) > 1.0))
          {
            tags[m_params.m_numLevels-2] |= iv;
          }
        }
      }
    }
    else if (m_params.m_tagType == "new_sink")
    {
      double fineDx = m_params.m_dx;
      for (int i = 0; i < m_params.m_numLevels - 2; i ++)
      {
        fineDx /= m_params.m_refRatio[i];
      }

      DisjointBoxLayout curGrids = ebisPtr->getFlowGrids(depth);

      for (DataIterator dit = curGrids.dataIterator(); dit.ok(); ++dit)
      {
        const Box& box = curGrids.get(dit());
        for (BoxIterator bit(box); bit.ok(); ++bit)
        {
          const IntVect& iv = bit();

          RealVect x = (iv + 0.5*RealVect::Unit) * fineDx + m_params.m_loCorner;

          if ((x[0] > -2.0) && (x[0] * x[0] / (4.95 * 4.95) + x[1] * x[1] / (19.95 * 19.95) > 1.0))
          {
            tags[m_params.m_numLevels-2] |= iv;
          }
        }
      }
    }
    else
    {
      string errorMessage = "DiffusionSolver::initIndexSpace:  Unknown tag type: '" + m_params.m_tagType + "'";
      MayDay::Error(errorMessage.c_str());
    }

    Vector<Vector<Box> > oldBoxes(m_params.m_numLevels);
    Vector<Vector<Box> > newBoxes;

    // Need the boxes on the coarsest level and the tags on the second to
    // finest level to make all the boxes on all the levels
    oldBoxes[0] = coarsestBoxes;

    // Make all the boxes on all the levels
    meshRefine.regrid(newBoxes,
                      tags,
                      0,
                      m_params.m_numLevels-1,
                      oldBoxes);

    // Go through all the new levels, Morton order the boxes, load balance the
    // result, and create the data structures needed
    ProblemDomain curDomain = m_params.m_coarsestDomain;
    for (int ilev = 1; ilev < m_params.m_numLevels; ilev++)
    {
      Vector<int> curProcs;

      curDomain.refine(m_params.m_refRatio[ilev-1]);

      // Order them using a space filling curve
      mortonOrdering(newBoxes[ilev]);

      // Do load balancing using timing information
      EBEllipticLoadBalance(curProcs,
                            newBoxes[ilev],
                            curDomain);

      // Create more of the infrastructure for geometry
      m_grids[ilev] = DisjointBoxLayout(newBoxes[ilev],
                                        curProcs,
                                        curDomain);

      ebisPtr->fillEBISLayout(m_ebisl[ilev],
                              m_grids[ilev],
                              curDomain,
                              m_params.m_numGhostEBISLayout);
    }
  }

  ParmParse pp;
  if (pp.contains("abr_file_output_name"))
  {
    string abrFile;
    pp.get("abr_file_output_name",  abrFile);

    Vector< Vector<Box> > boxes(m_params.m_numLevels);
    Vector< Vector<int> > ranks(m_params.m_numLevels);

    for (int ilev = 0; ilev < m_params.m_numLevels; ilev++)
    {
      boxes[ilev] = m_grids[ilev].boxArray();
      ranks[ilev] = m_grids[ilev].procIDs();
    }

    Box baseProbBox = m_params.m_coarsestDomain.domainBox();
    writeABRfile(boxes, ranks,
                 m_params.m_refRatio,
                 m_params.m_numLevels,
                 numProc(),
                 baseProbBox,
                 abrFile);
  }
}

void DiffusionSolver::initData()
{
  CH_TIME("DiffusionSolver::initData");

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
    EBCellFactory ebCellFactory(m_ebisl[ilev]);

    m_solnOld[ilev] = new LevelData<EBCellFAB>(m_grids[ilev],
                                               oneVariable,
                                               m_params.m_numGhostSoln,
                                               ebCellFactory);
    EBLevelDataOps::setVal(*(m_solnOld[ilev]),m_params.m_initialValue);

    m_solnNew[ilev] = new LevelData<EBCellFAB>(m_grids[ilev],
                                               oneVariable,
                                               m_params.m_numGhostSoln,
                                               ebCellFactory);

    m_source[ilev] = new LevelData<EBCellFAB>(m_grids[ilev],
                                              oneVariable,
                                              m_params.m_numGhostSource,
                                              ebCellFactory);
  }

  for (int ilev = 0; ilev < m_params.m_numLevels; ilev++)
  {
    // At each level, count and report the number of irregular cells and the
    // number of computational cells
    Real totalFortran = 0.0;
    Real totalIrreg   = 0.0;

    DisjointBoxLayout     fineGrids   = m_grids[ilev];
    LevelData<EBCellFAB>* fineDataPtr = m_solnOld[ilev];

    for (DataIterator dit=fineGrids.dataIterator(); dit.ok(); ++dit)
    {
      const Box&        curBox      = fineGrids[dit()];
      const EBISBox&    curEBISBox  = (*fineDataPtr)[dit()].getEBISBox();
      const IntVectSet& curIrregIVS = curEBISBox.getIrregIVS(curBox);

      totalFortran += curBox.numPts();
      totalIrreg   += curIrregIVS.numPts();
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

    status = MPI_Allreduce(&totalIrreg, &mpiTotal, 1, MPI_CH_REAL,
                           MPI_SUM, Chombo_MPI::comm);

    if (status != MPI_SUCCESS)
    {
      MayDay::Error("MPI error summing 'totalIrreg'");
    }

    totalIrreg = mpiTotal;
#endif

    long long totalFortranInt = totalFortran;
    long long totalIrregInt   = totalIrreg;

    std::ios::fmtflags origFlags = pout().flags();
    int origWidth = pout().width();
    int origPrecision = pout().precision();

    pout() << "Level " << ilev << ":" << "\n";
    pout() << setiosflags(ios::right);
    pout() << "  Total computation cells: " << setw(10) << totalFortranInt << "\n";
    pout() << "  Total irregular cells:   " << setw(10) << totalIrregInt   << "\n";
    pout() << "\n";

    pout().flags(origFlags);
    pout().width(origWidth);
    pout().precision(origPrecision);
  }
}

void DiffusionSolver::initSolver()
{
  CH_TIME("DiffusionSolver::initSolver");
  if (m_params.m_useVariableCoeff)
  {
    initSolverVariableCoeff();
  }
  else
  {
    initSolverConstantCoeff();
  }
}

void DiffusionSolver::getEBLGAndQuadCFI(Vector<EBLevelGrid>                   & a_ebLevelGrids,
                                        Vector<RefCountedPtr<EBQuadCFInterp> >& a_quadCFInterp)
{
  a_ebLevelGrids.resize(m_grids.size());
  a_quadCFInterp.resize(m_grids.size());

  // Define the data holders and interpolators
  ProblemDomain levelDomain = m_params.m_coarsestDomain;
  ProblemDomain coarserDomain;
  for (int ilev = 0; ilev < m_grids.size(); ilev++)
  {
    a_ebLevelGrids[ilev].define(m_grids[ilev],m_ebisl[ilev],levelDomain);

    if (ilev > 0)
    {
      int numVariables = 1;

      a_quadCFInterp[ilev] = RefCountedPtr<EBQuadCFInterp>
        (new EBQuadCFInterp(m_grids[ilev],
                            m_grids[ilev-1],
                            m_ebisl[ilev],
                            m_ebisl[ilev-1],
                            coarserDomain,
                            m_params.m_refRatio[ilev-1],
                            numVariables,
                            *(a_ebLevelGrids[ilev].getCFIVS())));
    }

    coarserDomain = levelDomain;

    if (ilev < m_grids.size()-1)
    {
      levelDomain.refine(m_params.m_refRatio[ilev]);
    }
  }
}

void DiffusionSolver::setConductivityCoefs(LevelData<EBCellFAB>          &    a_aco,
                                           LevelData<EBFluxFAB>          &    a_bco,
                                           LevelData<BaseIVFAB<Real> >   &    a_bcoIrreg,
                                           const Real                    &    a_dx)
{
  CH_TIME("DiffusionSolver::setConductivityCoefs");
  const DisjointBoxLayout& dbl = a_aco.disjointBoxLayout();
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
  {
    a_aco[dit()].setVal(1.0);
    const EBISBox& ebisBox = a_aco[dit()].getEBISBox();
    const Box&         box = dbl.get(dit());
    IntVectSet ivsIrreg = ebisBox.getIrregIVS(box);
    for (VoFIterator vofit(ivsIrreg, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      RealVect bndryCentroid = ebisBox.bndryCentroid(vofit());
      bndryCentroid *= a_dx;
      RealVect bndryLoc = EBArith::getVofLocation(vofit(), a_dx*RealVect::Unit, bndryCentroid);
      Real bcopoint;
      FORT_POINTGETBCODIFF(CHF_REAL(bcopoint),
                           CHF_REALVECT(bndryLoc),
                           CHF_REAL(m_time),
                           CHF_REAL(m_params.m_diffusionConstant),
                           CHF_REAL(m_params.m_diffusionEps)
                           );

      a_bcoIrreg[dit()](vofit(), 0) = bcopoint;
    }

    for (int idir = 0; idir < SpaceDim; idir++)
    {
      EBFaceFAB&       bcodir = a_bco[dit()][idir];
      BaseFab<Real>& svbcodir =   bcodir.getSingleValuedFAB();
      //do all as regular faces first
      FORT_GETBCODIFF(CHF_FRA1(svbcodir,0),
                      CHF_CONST_REAL(m_params.m_dx),
                      CHF_CONST_REAL(m_time),
                      CHF_BOX(svbcodir.box()),
                      CHF_CONST_REAL(m_params.m_diffusionConstant),
                      CHF_CONST_REAL(m_params.m_diffusionEps),
                      CHF_CONST_INT(idir)
                      );

      for (FaceIterator faceit(ivsIrreg, ebisBox.getEBGraph(), idir, FaceStop::SurroundingWithBoundary); faceit.ok(); ++faceit)
      {
        const FaceIndex& face =  faceit();
        if (box.contains(face.gridIndex(Side::Lo)) || box.contains(face.gridIndex(Side::Hi)))
        {
          RealVect bndryCentroid = ebisBox.centroid(face);
          bndryCentroid *= a_dx;
          RealVect faceLoc = EBArith::getFaceLocation(face, a_dx*RealVect::Unit, bndryCentroid);
          Real bcopoint;
          FORT_POINTGETBCODIFF(CHF_REAL(bcopoint),
                               CHF_REALVECT(faceLoc),
                               CHF_REAL(m_time),
                               CHF_REAL(m_params.m_diffusionConstant),
                               CHF_REAL(m_params.m_diffusionEps)
                               );
          bcodir(face, 0) = bcopoint;
        }
      }
    }
  }
}

void DiffusionSolver::defineVariableCoeffs(Vector<RefCountedPtr<LevelData<EBCellFAB> > >&           a_aco,
                                           Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&           a_bco,
                                           Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >&    a_bcoIrreg)
{
  CH_TIME("DiffusionSolver::defineConductivityCoef");
  a_aco.resize(        m_params.m_numLevels);
  a_bco.resize(        m_params.m_numLevels);
  a_bcoIrreg.resize(   m_params.m_numLevels);
  Real dxLev = m_params.m_dx;
  ProblemDomain domLev =   m_params.m_coarsestDomain;
  for (int ilev = 0; ilev < m_params.m_numLevels; ilev++)
  {
    LayoutData<IntVectSet> irregSets(m_grids[ilev]);
    int nghost = 1;
    for (DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
    {
      //has to correspond to number of ghost cells
      Box grownBox = grow(m_grids[ilev].get(dit()), nghost);
      grownBox &= domLev;
      irregSets[dit()] = m_ebisl[ilev][dit()].getIrregIVS(grownBox);
    }

    EBFluxFactory        ebfluxfact(m_ebisl[ilev]);
    EBCellFactory        ebcellfact(m_ebisl[ilev]);
    BaseIVFactory<Real>  baseivfact(m_ebisl[ilev], irregSets);

    a_aco[ilev]         = RefCountedPtr<LevelData<EBCellFAB       > >(new LevelData<EBCellFAB       >(m_grids[ilev], 1, nghost*IntVect::Unit, ebcellfact));
    a_bco[ilev]         = RefCountedPtr<LevelData<EBFluxFAB       > >(new LevelData<EBFluxFAB       >(m_grids[ilev], 1, nghost*IntVect::Unit, ebfluxfact));
    a_bcoIrreg[ilev]    = RefCountedPtr<LevelData<BaseIVFAB<Real> > >(new LevelData<BaseIVFAB<Real> >(m_grids[ilev], 1, nghost*IntVect::Unit, baseivfact));

    setConductivityCoefs(*a_aco[ilev], *a_bco[ilev], *a_bcoIrreg[ilev], dxLev);
    dxLev /=      m_params.m_refRatio[ilev];
    domLev.refine(m_params.m_refRatio[ilev]);
  }
}

void DiffusionSolver::getVariableCoeffOpFactory(RefCountedPtr<EBConductivityOpFactory>& a_factory)
{
  // Set up the no flux domain and embedded boundary conditions
  RefCountedPtr<NeumannConductivityDomainBCFactory> domBC(new NeumannConductivityDomainBCFactory());
  domBC->setValue(0.0);

  RefCountedPtr<NeumannConductivityEBBCFactory>      ebBC(new NeumannConductivityEBBCFactory());
  ebBC->setValue(0.0);

  Vector<EBLevelGrid>  eblg;
  Vector<RefCountedPtr<EBQuadCFInterp> > quadCFI;
  getEBLGAndQuadCFI(eblg, quadCFI);
  Vector<RefCountedPtr<LevelData<EBCellFAB> > >           aco;
  Vector<RefCountedPtr<LevelData<EBFluxFAB> > >           bco;
  Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >    bcoIrreg;

  defineVariableCoeffs(aco, bco, bcoIrreg);

  //coefficients come in through the =coefficients.
  Real unity = 1.0;
//  int relaxType = 0;
//  pout() << "using multicolored gauss seidel" << endl;
  int relaxType= m_params.m_mgRelaxType;
  if (relaxType == 1)
  {
    pout() << "using multi-colored gauss seidel relaxation" << endl;
  }
  else if (relaxType == 0)
  {
    pout() << "using point jacobi relaxation" << endl;
  }
  else if (relaxType == 2)
  {
    pout() << "using gsrb fast relaxation" << endl;
  }
  else
  {
    MayDay::Error("bogus relaxType for variable coefficients");
  }

  a_factory = RefCountedPtr<EBConductivityOpFactory>
    (new EBConductivityOpFactory(eblg, quadCFI, unity, unity, aco, bco, bcoIrreg,
                                 m_params.m_dx,  m_params.m_refRatio, domBC, ebBC,
                                 m_params.m_numGhostSoln, m_params.m_numGhostSource, relaxType));

}

void DiffusionSolver::initSolverVariableCoeff()
{
  CH_TIME("DiffusionSolver::initSolverVariableCoeff");
  // This is the multigrid solver used for backward Euler
  RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB > > >
    solver(new AMRMultiGrid<LevelData<EBCellFAB> >);

  RefCountedPtr<EBConductivityOpFactory> operatorFactory;

  getVariableCoeffOpFactory(operatorFactory);

  // Set the verbosity of the bottom solver for multigrid
  m_bottomSolver.m_verbosity = 0;

  // Define the multigrid solver and set various parameters
  solver->define(m_params.m_coarsestDomain,
                 *operatorFactory,
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
  if (m_integrator  != NULL) delete m_integrator;
  m_integrator = new EBBackwardEuler(solver,
                                     *operatorFactory,
                                     m_params.m_coarsestDomain,
                                     m_params.m_refRatio,
                                     m_params.m_numLevels);

}

void DiffusionSolver::initSolverConstantCoeff()
{
  CH_TIME("DiffusionSolver::initSolverConstantCoeff");

  // This is the multigrid solver used for backward Euler
  RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB > > >
    solver(new AMRMultiGrid<LevelData<EBCellFAB> >);

  // These are data holders and interpolators for the Helmholtz operator
  Vector<EBLevelGrid> ebLevelGrids(m_grids.size());
  Vector<RefCountedPtr<EBQuadCFInterp> > quadCFInterp(m_grids.size());

  // Define the data holders and interpolators
  ProblemDomain levelDomain = m_params.m_coarsestDomain;
  ProblemDomain coarserDomain;
  for (int ilev = 0; ilev < m_grids.size(); ilev++)
    {
      ebLevelGrids[ilev].define(m_grids[ilev],m_ebisl[ilev],levelDomain);

      if (ilev > 0)
        {
          int numVariables = 1;

          quadCFInterp[ilev] = RefCountedPtr<EBQuadCFInterp>
            (new EBQuadCFInterp(m_grids[ilev],
                                m_grids[ilev-1],
                                m_ebisl[ilev],
                                m_ebisl[ilev-1],
                                coarserDomain,
                                m_params.m_refRatio[ilev-1],
                                numVariables,
                                *(ebLevelGrids[ilev].getCFIVS())));
        }

      coarserDomain = levelDomain;

      if (ilev < m_grids.size()-1)
        {
          levelDomain.refine(m_params.m_refRatio[ilev]);
        }
    }

  RealVect vectDx = m_params.m_dx * RealVect::Unit;

  // Set up the no flux domain and embedded boundary conditions
  RefCountedPtr<NeumannPoissonDomainBCFactory> domainBCFactory(new NeumannPoissonDomainBCFactory());
  domainBCFactory->setValue(0.0);

  RefCountedPtr<NeumannPoissonEBBCFactory> ebBCFactory(new NeumannPoissonEBBCFactory());
  ebBCFactory->setValue(0.0);

  // Minimize the ghost cell filling in the operator relaxation
  EBAMRPoissonOp::doLazyRelax(m_params.m_mgLazyRelax);

  // Set up a factory to produce Helmholtz operators at various resolutions
  EBAMRPoissonOpFactory operatorFactory(ebLevelGrids,
                                        m_params.m_refRatio,
                                        quadCFInterp,
                                        vectDx,
                                        m_params.m_loCorner,
                                        m_params.m_mgNumPrecondIter,
                                        m_params.m_mgRelaxType,
                                        domainBCFactory,
                                        ebBCFactory,
                                        1.0,
                                        m_params.m_diffusionConstant,
                                        0.0,
                                        m_params.m_numGhostSoln,
                                        m_params.m_numGhostSource,
                                        m_params.m_numLevels);

  operatorFactory.setWhichReflux(m_params.m_whichReflux);

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
  if (m_integrator != NULL)
  {
    delete m_integrator;
  }
  m_integrator = new EBBackwardEuler(solver,
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
    LevelData<EBCellFAB>* sourceLevel = m_source[ilev];
    LevelData<EBCellFAB>* solnLevel   = m_solnOld[ilev];

    // Set the source/sink term for each grid
    for (DataIterator dit=curGrids.dataIterator(); dit.ok(); ++dit)
    {
      // Current grid
      const Box& curBox = curGrids[dit()];

      // Single grid data holders
      EBCellFAB& sourceEBCellFAB = (*sourceLevel)[dit()];
      EBCellFAB& solnEBCellFAB   = (*solnLevel)[dit()];

      // Single grid data holders for single valued cells
      BaseFab<double>& sourceFAB = sourceEBCellFAB.getSingleValuedFAB();
      BaseFab<double>& solnFAB   = solnEBCellFAB.getSingleValuedFAB();

      // Check the geometry type to determine the source/sink to use
      if (m_params.m_geometry == "noEB" ||
          m_params.m_geometry == "simple")
      {
        // Call Fortran to set the source/sink for all the single valued cells
        FORT_SIMPLESOURCESINK(CHF_FRA1(sourceFAB,0),
                              CHF_CONST_FRA1(solnFAB,0),
                              CHF_CONST_REAL(m_params.m_sourceScaling),
                              CHF_CONST_REAL(m_params.m_sinkScaling),
                              CHF_CONST_REAL(dx),
                              CHF_CONST_REALVECT(m_params.m_loCorner),
                              CHF_BOX(curBox));

        // Get all the machinery need to update multi-valued cells
        const EBISBox& curEBISBox = sourceEBCellFAB.getEBISBox();
        const EBGraph& curEBGraph = curEBISBox.getEBGraph();
        IntVectSet multiCells = curEBISBox.getMultiCells(curBox);

        // Set the source/sink for all the multi-value cells
        for (VoFIterator vofit(multiCells,curEBGraph); vofit.ok(); ++vofit)
        {
          const VolIndex& volIndex = vofit();
          const IntVect&  iv = volIndex.gridIndex();

          Real source;

          // Call Fortran for a single value and assign the result
          FORT_SIMPLESOURCESINKPOINT(CHF_REAL(source),
                                     CHF_CONST_REAL(solnEBCellFAB(volIndex,0,1)),
                                     CHF_CONST_REAL(m_params.m_sourceScaling),
                                     CHF_CONST_REAL(m_params.m_sinkScaling),
                                     CHF_CONST_REAL(dx),
                                     CHF_CONST_REALVECT(m_params.m_loCorner),
                                     CHF_CONST_INTVECT(iv));

          sourceEBCellFAB(volIndex,0,1) = source;
        }
      }
      else if (m_params.m_geometry == "keratocyte")
      {
        // Call Fortran to set the source/sink for all the single valued cells
        FORT_KERASOURCESINK(CHF_FRA1(sourceFAB,0),
                            CHF_CONST_FRA1(solnFAB,0),
                            CHF_CONST_REAL(m_params.m_sourceScaling),
                            CHF_CONST_REAL(m_params.m_sinkScaling),
                            CHF_CONST_REAL(dx),
                            CHF_CONST_REALVECT(m_params.m_loCorner),
                            CHF_BOX(curBox));

        // Get all the machinery need to update multi-valued cells
        const EBISBox& curEBISBox = sourceEBCellFAB.getEBISBox();
        const EBGraph& curEBGraph = curEBISBox.getEBGraph();
        IntVectSet multiCells = curEBISBox.getMultiCells(curBox);

        // Set the source/sink for all the multi-value cells
        for (VoFIterator vofit(multiCells,curEBGraph); vofit.ok(); ++vofit)
        {
          const VolIndex& volIndex = vofit();
          const IntVect&  iv = volIndex.gridIndex();

          Real source;

          // Call Fortran for a single value and assign the result
          FORT_KERASOURCESINKPOINT(CHF_REAL(source),
                                   CHF_CONST_REAL(solnEBCellFAB(volIndex,0,1)),
                                   CHF_CONST_REAL(m_params.m_sourceScaling),
                                   CHF_CONST_REAL(m_params.m_sinkScaling),
                                   CHF_CONST_REAL(dx),
                                   CHF_CONST_REALVECT(m_params.m_loCorner),
                                   CHF_CONST_INTVECT(iv));

          sourceEBCellFAB(volIndex,0,1) = source;
        }
      }
      else if (m_params.m_geometry == "new_keratocyte")
      {
        // Call Fortran to set the source/sink for all the single valued cells
        FORT_NEWKERASOURCESINK(CHF_FRA1(sourceFAB,0),
                               CHF_CONST_FRA1(solnFAB,0),
                               CHF_CONST_REAL(m_params.m_sourceScaling),
                               CHF_CONST_REAL(m_params.m_sinkScaling),
                               CHF_CONST_REAL(dx),
                               CHF_CONST_REALVECT(m_params.m_loCorner),
                               CHF_BOX(curBox));

        // Get all the machinery need to update multi-valued cells
        const EBISBox& curEBISBox = sourceEBCellFAB.getEBISBox();
        const EBGraph& curEBGraph = curEBISBox.getEBGraph();
        IntVectSet multiCells = curEBISBox.getMultiCells(curBox);

        // Set the source/sink for all the multi-value cells
        for (VoFIterator vofit(multiCells,curEBGraph); vofit.ok(); ++vofit)
        {
          const VolIndex& volIndex = vofit();
          const IntVect&  iv = volIndex.gridIndex();

          Real source;

          // Call Fortran for a single value and assign the result
          FORT_NEWKERASOURCESINKPOINT(CHF_REAL(source),
                                      CHF_CONST_REAL(solnEBCellFAB(volIndex,0,1)),
                                      CHF_CONST_REAL(m_params.m_sourceScaling),
                                      CHF_CONST_REAL(m_params.m_sinkScaling),
                                      CHF_CONST_REAL(dx),
                                      CHF_CONST_REALVECT(m_params.m_loCorner),
                                      CHF_CONST_INTVECT(iv));

          sourceEBCellFAB(volIndex,0,1) = source;
        }
      }
      else
      {
        string errorMessage = "DiffusionSolver::setSource:  Unknown geometry: '" + m_params.m_geometry + "'";
        MayDay::Error(errorMessage.c_str());
      }
    }

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

  // Leave the values in the covered cells alone
  bool replaceCovered = false;
  Vector<Real> dummy;

  // Write the output file
  writeEBHDF5(filename,
              m_grids,
              m_solnOld,
              names,
              m_params.m_coarsestDomain,
              m_params.m_dx,
              m_params.m_dt,
              a_time,
              m_params.m_refRatio,
              m_params.m_numLevels,
              replaceCovered,
              dummy);
}
