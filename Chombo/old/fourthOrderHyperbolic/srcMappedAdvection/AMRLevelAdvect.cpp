#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iomanip>

#include <string>
#include "parstream.H"

#include "Box.H"
#include "Vector.H"
#include "IntVectSet.H"
#include "FArrayBox.H"
#include "NodeFArrayBox.H"
#include "DisjointBoxLayout.H"
#include "LayoutIterator.H"
#include "CH_HDF5.H"
#include "SPMD.H"
#include "LoadBalance.H"
#include "LevelFluxRegister.H"
#include "ProblemDomain.H"
#include "BoxIterator.H"
#include "computeSum.H"
#include "PiecewiseLinearFillPatch.H"
#include "computeNorm.H"
#include "EdgeToCell.H"
#include "CellToEdge.H"
//#include "SingleLevelDivergence.H"
#include "computeMappedDt.H"
#include "mappedCTU.H"
#include "mappedDCU.H"
#include "MappedPhysIBC.H"

#include "NodeAMRIO.H"
#include "MappedAdvectionPhysics.H"
#include "AdvectPhysicsF_F.H"
#include "zalesakLimiter.H"
#include "minValRedistribution.H"

#include "LevelRK4.H"

#include "AMRLevel.H"
#include "AMRLevelAdvect.H"

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// Constructor
AMRLevelAdvect::AMRLevelAdvect()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelAdvect default constructor" << endl;
  }

  m_levelStep = 0;
  m_cfl = 0.8;
  m_spaceOrder = 4;
  m_limitFaceValues = true;
  m_enforceMinVal = false;
  m_minVal = -100000.0;
  m_lowOrderFluxScheme = 1;
  m_redistributeNegativeVal = false;
  m_maxRedistributionPasses = 0;
  m_highOrderFluxesOnly = false;
  m_lowOrderFluxesOnly = false;
  m_domainLength = 1.0;
  m_refineThresh = 0.2;
  m_initial_dt_multiplier = 0.1;
  m_coordSysFactPtr = NULL;
}

////////////////////////////////////////////////////////////////////////////////

// Destructor
AMRLevelAdvect::~AMRLevelAdvect()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelAdvect destructor" << endl;
  }
  delete m_coordSysPtr;
}

////////////////////////////////////////////////////////////////////////////////

// Define new AMR level
void AMRLevelAdvect::define(AMRLevel*            a_coarserLevelPtr,
                            const ProblemDomain& a_problemDomain,
                            int                  a_level,
                            int                  a_refRatio)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelAdvect::define " << a_level << endl;
  }

  // Call inherited define
  AMRLevel::define(a_coarserLevelPtr,
                   a_problemDomain,
                   a_level,
                   a_refRatio);

  // Get setup information from the next coarser level
  if (a_coarserLevelPtr != NULL)
  {
    AMRLevelAdvect* amrAdvectPtr = dynamic_cast<AMRLevelAdvect*>(a_coarserLevelPtr);

    if (amrAdvectPtr != NULL)
    {
      m_cfl = amrAdvectPtr->m_cfl;
      m_spaceOrder = amrAdvectPtr->m_spaceOrder;
      m_limitFaceValues = amrAdvectPtr->m_limitFaceValues;
      m_enforceMinVal = amrAdvectPtr->m_enforceMinVal;
      m_minVal = amrAdvectPtr->m_minVal;
      m_lowOrderFluxScheme = amrAdvectPtr->m_lowOrderFluxScheme;
      m_redistributeNegativeVal = amrAdvectPtr->m_redistributeNegativeVal;
      m_maxRedistributionPasses = amrAdvectPtr->m_maxRedistributionPasses;
      m_highOrderFluxesOnly = amrAdvectPtr->m_highOrderFluxesOnly;
      m_lowOrderFluxesOnly = amrAdvectPtr->m_lowOrderFluxesOnly;
      m_useHyperviscosity = amrAdvectPtr->m_useHyperviscosity;
      m_hyperviscosity = amrAdvectPtr->m_hyperviscosity;
      m_domainLength = amrAdvectPtr->m_domainLength;
      m_refineThresh = amrAdvectPtr->m_refineThresh;
      m_tagBufferSize = amrAdvectPtr->m_tagBufferSize;
    }
    else
    {
      MayDay::Error("AMRLevelAdvect::define: a_coarserLevelPtr is not castable to AMRLevelAdvect*");
    }
  }

  // Compute the grid spacing
  m_dx = m_domainLength / a_problemDomain.domainBox().longside();

  // Nominally, one layer of ghost cells is maintained permanently and
  // individual computations may create local data with more

  // Begin application-dependent code - dfm

  // if we're limiting, need 3 ghost cells in phi. otherwise, 2 will do
  m_numGhost = 3;
  m_numStates = 1;
  m_stateNames.push_back("phi_0");
  //m_stateNames.push_back("phi_1");


  // End application-dependent code - dfm

}

////////////////////////////////////////////////////////////////////////////////

// Advance by one timestep
Real AMRLevelAdvect::advance()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelAdvect::advance level " << m_level << " to time " << m_time << endl;
  }

  // Begin application-dependent code - dfm



  // Copy the new to the old
  m_phiNew.copyTo(m_phiNew.interval(),
                  m_phiOld,
                  m_phiOld.interval());

  // End application-dependent code - dfm

  //XXX -- unused
  //XXXReal newDt = 0.0;

  // Set up arguments to LevelGodunov::step based on whether there are
  // coarser and finer levels

  // Undefined flux register in case we need it
  LevelFluxRegister dummyFR;

  // Undefined leveldata in case we need it
  LevelData<FArrayBox> dummyData;

  // Set arguments to dummy values and then fix if real values are available
  LevelFluxRegister* coarserFR = &dummyFR;
  LevelFluxRegister* finerFR   = &dummyFR;

  Real tCoarserOld = 0.0;
  Real tCoarserNew = 0.0;

  LevelData<FArrayBox>* coarserPhiOld = &dummyData;
  LevelData<FArrayBox>* coarserPhiNew = &dummyData;
  // A coarser level exists
  if (m_hasCoarser)
  {
    AMRLevelAdvect* coarserPtr = getCoarserLevel();

    // Recall that my flux register goes between my level and the next
    // finer level
    coarserFR = &coarserPtr->m_fluxRegister;

    coarserPhiOld = &coarserPtr->m_phiOld;
    coarserPhiNew = &coarserPtr->m_phiNew;

    tCoarserNew = coarserPtr->m_time;
    tCoarserOld = tCoarserNew - coarserPtr->m_dt;
  }

  // A finer level exists
  if (m_hasFiner)
  {
    // Recall that my flux register goes between my level and the next
    // finer level
    finerFR = &m_fluxRegister;
    finerFR->setToZero();
  }

  // Advance advect equation by one time step using 4th-order
  // Runge-Kutta.

  m_levelAdvectOperator.setAdvectionVel(&m_advVel);
  m_levelAdvectOperator.setBC(m_basicIBCPtr);
  m_levelAdvectOperator.setCoordSys(m_coordSysPtr);
  m_levelAdvectOperator.resetFluxes();
//  m_levelAdvectOperator.estimateEigenvalues( m_phiNew, m_time, m_dt );

  RK4LevelAdvance<LevelData<FArrayBox>,
    LevelFluxRegister,
    LevelAdvectOperator>(m_phiNew,
                         m_phiOld,
                         *coarserPhiOld, tCoarserOld,
                         *coarserPhiNew, tCoarserNew,
                         *coarserFR, *finerFR,
                         m_time, m_dt, m_levelAdvectOperator);

  // if we're enforcing min or max values, compute 1st-order CTU
  // fluxes here
  LevelData<FluxBox> CTUflux;
  if (m_enforceMinVal)
    {
      LevelData<FluxBox>& fluxes = m_levelAdvectOperator.getFluxes();
      LevelData<FluxBox>& NTvel = m_levelAdvectOperator.getNTvel();

      // note that fluxes and NTvel returned by levelAdvectOperator
      // are multiplied by dt
      DataIterator dit = fluxes.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          fluxes[dit] *= 1.0/m_dt;
          NTvel[dit] *= 1.0/m_dt;
        }

      CTUflux.define(fluxes.getBoxes(), fluxes.nComp(), fluxes.ghostVect());

      // note that this returns dt*flux
      computeLowOrderFlux(CTUflux, m_phiOld, NTvel, m_time, m_dt);

      if (!m_highOrderFluxesOnly)
        {
          // now apply Zalesak limiter to enforce the min value
          applyZalesakLimiter(fluxes, CTUflux, m_phiOld, m_minVal,
                              m_dx, m_dt, m_lowOrderFluxesOnly);
        }

      FourthOrderCoordSys* FOCS = dynamic_cast<FourthOrderCoordSys*>(m_coordSysPtr);
      CH_assert(FOCS);

      // now need to re-do advective update with limited fluxes
      FOCS->simpleDivergence(m_phiNew, fluxes);

      Real hDinv = 1.0/(D_TERM(m_dx,*m_dx,*m_dx));
      for (dit.begin(); dit.ok(); ++dit)
        {
          FArrayBox& thisPhiNew = m_phiNew[dit];
          FArrayBox& thisPhiOld = m_phiOld[dit];

          thisPhiNew *= hDinv;
          thisPhiNew *= -m_dt;
          thisPhiNew += thisPhiOld;
        }
    } // end if we're applying the Zalesak limiter

  if (m_redistributeNegativeVal)
  {
     applyMinValRedistribution( m_phiNew, 0.0, m_maxRedistributionPasses );
  }

  // Update the time and store the new timestep
  ++m_levelStep;
  m_time += m_dt;
  Real returnDt = m_dt;

  m_dtNew = returnDt;

  return returnDt;
}

////////////////////////////////////////////////////////////////////////////////

// Things to do after a timestep
void AMRLevelAdvect::postTimeStep()
{
  // Used for conservation tests
  static Real orig_integral = 0.0;
  static Real last_integral = 0.0;
  static bool first = true;

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelAdvect::postTimeStep " << m_level << endl;
  }

  // Begin application-dependent code - PC.

  if (m_hasFiner)
  {
    // Reflux
    Real scale = -1.0/m_dx;
    m_fluxRegister.reflux(m_phiNew,scale);

    // Average from finer level data
    AMRLevelAdvect* amrAdvectFinerPtr = getFinerLevel();

    amrAdvectFinerPtr->m_coarseAverage.averageToCoarse(m_phiNew,
                                                       amrAdvectFinerPtr->m_phiNew);


/*
    amrAdvectFinerPtr->m_coarseAverage.averageToCoarse(m_phiNew,
                                                    amrAdvectFinerPtr->m_phiNew);
*/
    amrAdvectFinerPtr->m_levelAdvectOperator.avgdown(
                                                     amrAdvectFinerPtr->m_phiNew,
                                                     m_phiNew);
  }
  // End application-dependent code - PC.

  if (s_verbosity >= 2 && m_level == 0)
  {
    int nRefFine = 1;

    pout() << "AMRLevelAdvect::postTimeStep:" << endl;
    pout() << "  Sums:" << endl;
    for (int comp = 0; comp < m_numStates; comp++)
    {
      Interval curComp(comp,comp);
      // Begin application-dependent code - PC.
      Real integral = computeSum(m_phiNew,NULL,nRefFine,m_dx,curComp);
      // End application-dependent code - PC.

      pout() << "    " << setw(23)
                       << setprecision(16)
                       << setiosflags(ios::showpoint)
                       << setiosflags(ios::scientific)
                       << integral
             << " --- " << m_stateNames[comp];

      if (comp == 0 )
      {
        pout() << " (" << setw(23)
                       << setprecision(16)
                       << setiosflags(ios::showpoint)
                       << setiosflags(ios::scientific)
                       << (integral-last_integral)
               << " " << setw(23)
                      << setprecision(16)
                      << setiosflags(ios::showpoint)
                      << setiosflags(ios::scientific)
                      << (integral-orig_integral)
               << ")";
        if (first)
        {
          orig_integral = integral;
          first = false;
        }
        last_integral = integral;
      }
      pout() << endl;
    }
  }
  if (s_verbosity >= 1 && m_level == 0)
    {

      writeErrorNorms();

      //Real energy = computeEnergy();
      //pout() << " integrated energy = " << energy << endl ;

    }
  else if (m_level > 0)
    {
      MayDay::Warning("Error norm computation not correct for refined levels");
    }

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelAdvect::postTimeStep " << m_level << " finished" << endl;
  }
}

////////////////////////////////////////////////////////////////////////////////

// Create tags for regridding
void AMRLevelAdvect::tagCells(IntVectSet& a_tags)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelAdvect::tagCells " << m_level << endl;
  }

  // Since tags are calculated using only current time step data, use
  // the same tagging function for initialization and for regridding.
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelAdvect::tagCellsInit " << m_level << endl;
  }

  // Begin application-dependent code - PC.

  // Create tags based on advection term

  LevelData<FArrayBox> lOfPhi(m_phiNew.getBoxes(),m_numStates);

  // Set up arguments to LevelAdvectOperator::eval based on whether there are
  // coarser and finer levels

  // Undefined flux register in case we need it
  LevelFluxRegister dummyFR;

  // Undefined leveldata in case we need it
  LevelData<FArrayBox> dummyData;

  // Set arguments to dummy values and then fix if real values are available
  LevelFluxRegister* coarserFR = &dummyFR;
  LevelFluxRegister* finerFR   = &dummyFR;

  Real tCoarserOld = 0.0;
  Real tCoarserNew = 0.0;

  LevelData<FArrayBox>* coarserPhiOld = &dummyData;
  LevelData<FArrayBox>* coarserPhiNew = &dummyData;
  // A coarser level exists
  if (m_hasCoarser)
  {
    AMRLevelAdvect* coarserPtr = getCoarserLevel();

    coarserPhiOld = &coarserPtr->m_phiOld;
    coarserPhiNew = &coarserPtr->m_phiNew;

    tCoarserNew = coarserPtr->m_time;
    tCoarserOld = tCoarserNew - coarserPtr->m_dt;
  }
  IntVectSet localTags;
  const DisjointBoxLayout& levelDomain = m_phiNew.disjointBoxLayout();
  LevelAdvectOperator lwo;
  if (m_hasCoarser)
    {
        lwo.define(m_grids,
                   getCoarserLevel()->m_grids,
                   m_problem_domain,
                   m_ref_ratio,
                   m_numStates,
                   RealVect(D_DECL(m_dx, m_dx, m_dx)),
                   m_hasCoarser,
                   m_hasFiner);
    }
  else
    {
        lwo.define(m_grids,
                   DisjointBoxLayout(),
                   m_problem_domain,
                   m_ref_ratio,
                   m_numStates,
                   RealVect(D_DECL(m_dx, m_dx, m_dx)),
                   m_hasCoarser,
                   m_hasFiner);
    }
  lwo.setAdvectionVel(&m_advVel);
  lwo.setBC(m_basicIBCPtr);
  lwo.setCoordSys(m_coordSysPtr);

  lwo.evalRHS(lOfPhi,m_phiNew,
              *finerFR,*coarserFR,
              *coarserPhiOld,tCoarserOld,
              *coarserPhiNew,tCoarserNew,
              m_time,m_dt);

  DataIterator dit = levelDomain.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
  {
    const Box& b = levelDomain[dit()];
    FArrayBox& lPhiLocal = lOfPhi[dit()];

    FArrayBox tagVals(b,1);
    FORT_MAGNITUDEF(CHF_FRA1(tagVals,0),
                    CHF_CONST_FRA(lPhiLocal),
                    CHF_BOX(b));


  // Tag where gradient exceeds threshold

    BoxIterator bit(b);
    for (bit.begin(); bit.ok(); ++bit)
    {
      const IntVect& iv = bit();

      if (tagVals(iv) >= m_refineThresh/m_dx/m_dx)
      {
        localTags |= iv;
      }
    }
  }

  // End application-dependent code - PC.

  localTags.grow(m_tagBufferSize);

  // Need to do this in two steps unless a IntVectSet::operator &=
  // (ProblemDomain) operator is defined
  Box localTagsBox = localTags.minBox();
  localTagsBox &= m_problem_domain;
  localTags &= localTagsBox;
  a_tags = localTags;
  /*
  Vector<IntVectSet> allTags;

  const int destProc = uniqueProc(SerialTask::compute);

  gather(allTags,localTags,destProc);

  if (procID() == uniqueProc(SerialTask::compute))
  {
    for (int i = 0; i < allTags.size(); ++i)
    {
      a_tags |= allTags[i];
    }
  }
  */
}

////////////////////////////////////////////////////////////////////////////////

// Create tags at initialization
void AMRLevelAdvect::tagCellsInit(IntVectSet& a_tags)
{
  tagCells(a_tags);
}

////////////////////////////////////////////////////////////////////////////////

// Set up data on this level after regridding
void AMRLevelAdvect::regrid(const Vector<Box>& a_newGrids)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelAdvect::regrid " << m_level << endl;
  }

  // Save original grids and load balance
  m_level_grids = a_newGrids;
  m_grids = loadBalance(a_newGrids);

  if (s_verbosity >= 4)
  {
    // Indicate/guarantee that the indexing below is only for reading
    // otherwise an error/assertion failure occurs
    const DisjointBoxLayout& constGrids = m_grids;

    pout() << "new grids: " << endl;

    for (LayoutIterator lit = constGrids.layoutIterator(); lit.ok(); ++lit)
    {
      pout() << constGrids[lit()] << endl;
    }
  }

  // Save data for later
  // Begin application-dependent code - PC.

  LevelData<FArrayBox> phiOld;
  phiOld.define(m_phiNew);
  m_phiNew.copyTo(m_phiNew.interval(),
                  phiOld,
                  phiOld.interval());


  // Reshape state with new grids
  IntVect ivGhost = m_numGhost*IntVect::Unit;
  m_phiNew.define(m_grids,m_numStates,ivGhost);
  m_phiOld.define(m_grids,m_numStates,ivGhost);
  m_advVel.define(m_grids,SpaceDim,ivGhost);

  // Set up data structures
  levelSetup();

  // Re-create the coordinate system
  m_coordSysPtr->regrid(m_grids);

  // Interpolate from coarser level
  if (m_hasCoarser)
  {
    AMRLevelAdvect* amrAdvectCoarserPtr = getCoarserLevel();
    m_fineInterp.interpToFine(m_phiNew,amrAdvectCoarserPtr->m_phiNew);

  // Begin application-dependent code - PC.

  }

  // Copy from old state
  phiOld.copyTo(phiOld.interval(),
                m_phiNew,
                m_phiNew.interval());

}

////////////////////////////////////////////////////////////////////////////////

// Initialize grids
void AMRLevelAdvect::initialGrid(const Vector<Box>& a_newGrids)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelAdvect::initialGrid " << m_level << endl;
  }

  // Save original grids and load balance
  m_level_grids = a_newGrids;
  m_grids = loadBalance(a_newGrids);

  if (s_verbosity >= 4)
  {
    // Indicate/guarantee that the indexing below is only for reading
    // otherwise an error/assertion failure occurs
    const DisjointBoxLayout& constGrids = m_grids;

    pout() << "new grids: " << endl;
    for (LayoutIterator lit = constGrids.layoutIterator(); lit.ok(); ++lit)
    {
      pout() << constGrids[lit()] << endl;
    }
  }

  // Define old and new state data structures
  // Begin application-dependent code - PC.

  IntVect ivGhost = m_numGhost*IntVect::Unit;
  m_phiNew.define(m_grids,m_numStates,ivGhost);
  m_phiOld.define(m_grids,m_numStates,ivGhost);
  m_advVel.define(m_grids,SpaceDim,ivGhost);
  // End application-dependent code - PC.

  // Create the coordinate system
  IntVect ghostVect(D_DECL(m_numGhost, m_numGhost, m_numGhost));
  m_coordSysPtr = m_coordSysFactPtr->getCoordSys(m_grids,
                                                 m_problem_domain,
                                                 ghostVect);

  // Set up data structures
  levelSetup();


}

////////////////////////////////////////////////////////////////////////////////

// Initialize data
void AMRLevelAdvect::initialData()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelAdvect::initialData " << m_level << endl;
  }
// Begin application-dependent code - PC.

  m_basicIBCPtr->initialize(m_phiNew,m_problem_domain,*m_coordSysPtr,m_dx);


  m_basicIBCPtr->advVel(m_advVel,m_problem_domain,*m_coordSysPtr,m_dx,m_time);

// End application-dependent code - PC.
}

////////////////////////////////////////////////////////////////////////////////

// Things to do after initialization
void AMRLevelAdvect::postInitialize()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelAdvect::postInitialize " << m_level << endl;
  }

  if (s_verbosity >= 1 && m_level == 0)
    {
      writeErrorNorms();
    }
  else if (m_level > 0)
    {
      MayDay::Warning("Error norm computation not correct for refined levels");
    }

}

////////////////////////////////////////////////////////////////////////////////

#ifdef CH_USE_HDF5

// Write checkpoint header
void AMRLevelAdvect::writeCheckpointHeader(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelAdvect::writeCheckpointHeader" << endl;
  }

  // Setup the number of components
  HDF5HeaderData header;
  header.m_int["num_components"] = m_numStates*2;

  // Setup the component names
  char compStr[30];
  for (int comp = 0; comp < m_numStates*2; ++comp)
  {
    sprintf(compStr,"component_%d",comp);
    header.m_string[compStr] = m_stateNames[comp];
  }

  // Write the header
  header.writeToFile(a_handle);

  if (s_verbosity >= 3)
  {
    pout() << header << endl;
  }
}

// Write checkpoint data for this level
void AMRLevelAdvect::writeCheckpointLevel(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelAdvect::writeCheckpointLevel" << endl;
  }

  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;

  a_handle.setGroup(label);

  // Setup the level header information
  HDF5HeaderData header;

  header.m_int ["ref_ratio"]       = m_ref_ratio;
  header.m_int ["tag_buffer_size"] = m_tagBufferSize;
  header.m_real["dx"]              = m_dx;
  header.m_real["dt"]              = m_dt;
  header.m_real["time"]            = m_time;
  header.m_box ["prob_domain"]     = m_problem_domain.domainBox();

  // Setup the periodicity info
  D_TERM(
         if (m_problem_domain.isPeriodic(0))
           header.m_int ["is_periodic_0"] = 1;
         else
           header.m_int ["is_periodic_0"] = 0; ,

         if (m_problem_domain.isPeriodic(1))
           header.m_int ["is_periodic_1"] = 1;
         else
           header.m_int ["is_periodic_1"] = 0; ,

         if (m_problem_domain.isPeriodic(2))
           header.m_int ["is_periodic_2"] = 1;
         else
           header.m_int ["is_periodic_2"] = 0; );

  // Write the header for this level
  header.writeToFile(a_handle);

  if (s_verbosity >= 3)
  {
    pout() << header << endl;
  }

  // Write the data for this level
  LevelData<FArrayBox> outData(m_phiNew.getBoxes(),2*m_numStates);
  Interval interval0(0,m_numStates-1);
  Interval interval1(m_numStates,2*m_numStates-1);
  m_phiNew.copyTo(interval0,outData,interval0);
  write(a_handle,outData.boxLayout());
  write(a_handle,outData,"data");
}

// Read checkpoint header
void AMRLevelAdvect::readCheckpointHeader(HDF5Handle& a_handle)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelAdvect::readCheckpointHeader" << endl;
  }

  // Reader the header
  HDF5HeaderData header;
  header.readFromFile(a_handle);

  if (s_verbosity >= 3)
  {
    pout() << "hdf5 header data:" << endl;
    pout() << header << endl;
  }

  // Get the number of components
  if (header.m_int.find("num_components") == header.m_int.end())
  {
    MayDay::Error("AMRLevelAdvect::readCheckpointHeader: checkpoint file does not have num_components");
  }

  int numStates = header.m_int["num_components"];
  if (numStates != m_numStates*2)
  {
    MayDay::Error("AMRLevelAdvect::readCheckpointHeader: num_components in checkpoint file does not match solver");
  }

  // Get the component names
  std::string stateName;
  char compStr[60];
  for (int comp = 0; comp < m_numStates*2; ++comp)
  {
    sprintf(compStr,"component_%d",comp);
    if (header.m_string.find(compStr) == header.m_string.end())
    {
      MayDay::Error("AMRLevelAdvect::readCheckpointHeader: checkpoint file does not have enough component names");
    }

    stateName = header.m_string[compStr];
    if (stateName != m_stateNames[comp])
    {
      MayDay::Error("AMRLevelAdvect::readCheckpointHeader: state_name in checkpoint does not match solver");
    }
  }
}

// Read checkpoint data for this level
void AMRLevelAdvect::readCheckpointLevel(HDF5Handle& a_handle)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelAdvect::readCheckpointLevel" << endl;
  }

  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;

  // Read the header for this level
  a_handle.setGroup(label);

  HDF5HeaderData header;
  header.readFromFile(a_handle);

  if (s_verbosity >= 3)
  {
    pout() << "hdf5 header data:" << endl;
    pout() << header << endl;
  }

  // Get the refinement ratio
  if (header.m_int.find("ref_ratio") == header.m_int.end())
  {
    MayDay::Error("AMRLevelAdvect::readCheckpointLevel: file does not contain ref_ratio");
  }
  m_ref_ratio = header.m_int["ref_ratio"];

  if (s_verbosity >= 2)
  {
    pout() << "read ref_ratio = " << m_ref_ratio << endl;
  }

  // Get the tag buffer size
  if (header.m_int.find("tag_buffer_size") == header.m_int.end())
  {
    MayDay::Error("AMRLevelAdvect::readCheckpointLevel: file does not contain tag_buffer_size");
  }
  m_tagBufferSize=  header.m_int["tag_buffer_size"];

  if (s_verbosity >= 2)
  {
    pout() << "read tag_buffer_size = " << m_tagBufferSize << endl;
  }

  // Get dx
  if (header.m_real.find("dx") == header.m_real.end())
  {
    MayDay::Error("AMRLevelAdvect::readCheckpointLevel: file does not contain dx");
  }
  m_dx = header.m_real["dx"];

  if (s_verbosity >= 2)
  {
    pout() << "read dx = " << m_dx << endl;
  }

  // Get dt
  if (header.m_real.find("dt") == header.m_real.end())
  {
    MayDay::Error("AMRLevelAdvect::readCheckpointLevel: file does not contain dt");
  }
  m_dt = header.m_real["dt"];

  if (s_verbosity >= 2)
  {
    pout() << "read dt = " << m_dt << endl;
  }

  // Get time
  if (header.m_real.find("time") == header.m_real.end())
  {
    MayDay::Error("AMRLevelAdvect::readCheckpointLevel: file does not contain time");
  }
  m_time = header.m_real["time"];

  if (s_verbosity >= 2)
  {
    pout() << "read time = " << m_time << endl;
  }

  // Get the problem domain
  if (header.m_box.find("prob_domain") == header.m_box.end())
  {
    MayDay::Error("AMRLevelAdvect::readCheckpointLevel: file does not contain prob_domain");
  }

  Box domainBox = header.m_box["prob_domain"];

  // Get the periodicity info -- this is more complicated than it really
  // needs to be in order to preserve backward compatibility
  bool isPeriodic[SpaceDim];
  D_TERM(if (!(header.m_int.find("is_periodic_0") == header.m_int.end()))
           isPeriodic[0] =  (header.m_int["is_periodic_0"] == 1);
         else
           isPeriodic[0] = false; ,

         if (!(header.m_int.find("is_periodic_1") == header.m_int.end()))
           isPeriodic[1] =  (header.m_int["is_periodic_1"] == 1);
         else
           isPeriodic[1] = false; ,

         if (!(header.m_int.find("is_periodic_2") == header.m_int.end()))
           isPeriodic[2] =  (header.m_int["is_periodic_2"] == 1);
         else
           isPeriodic[2] = false;);

  m_problem_domain = ProblemDomain(domainBox,isPeriodic);

  // Get the grids
  Vector<Box> grids;
  const int gridStatus = read(a_handle,grids);

  if (gridStatus != 0)
  {
    MayDay::Error("AMRLevelAdvect::readCheckpointLevel: file does not contain a Vector<Box>");
  }

  // Create level domain
  m_grids = loadBalance(grids);

  // Indicate/guarantee that the indexing below is only for reading
  // otherwise an error/assertion failure occurs
  const DisjointBoxLayout& constGrids = m_grids;

  LayoutIterator lit = constGrids.layoutIterator();
  for (lit.begin(); lit.ok(); ++lit)
  {
    const Box& b = constGrids[lit()];
    m_level_grids.push_back(b);
  }

  if (s_verbosity >= 4)
  {
    pout() << "read level domain: " << endl;
    LayoutIterator lit = constGrids.layoutIterator();
    for (lit.begin(); lit.ok(); ++lit)
    {
      const Box& b = constGrids[lit()];
      pout() << lit().intCode() << ": " << b << endl;
    }
    pout() << endl;
  }

  // Reshape state with new grids
  IntVect ivGhost = m_numGhost*IntVect::Unit;
  LevelData<FArrayBox> inData ;
  inData.define(m_grids,m_numStates*2,ivGhost);
  int dataStatus = read<FArrayBox>(a_handle, inData, "data", m_grids);
  if (dataStatus != 0)
  {
    MayDay::Error("AMRLevelAdvect::readCheckpointLevel: file does not contain all state data");
  }
  m_phiNew.define(m_grids,m_numStates,ivGhost);
  m_phiOld.define(m_grids,m_numStates,ivGhost);
  m_advVel.define(m_grids,SpaceDim,ivGhost);

  Interval interval0(0,m_numStates-1);
  Interval interval1(m_numStates,2*m_numStates-1);
  inData.copyTo(interval0,m_phiNew,interval0);

  // Create the coordinate system
  IntVect ghostVect(D_DECL(m_numGhost, m_numGhost, m_numGhost));
  m_coordSysPtr = m_coordSysFactPtr->getCoordSys(m_grids,
                                                 m_problem_domain,
                                                 ghostVect);

  // Set up data structures
  levelSetup();


}

// Write plotfile header
void AMRLevelAdvect::writePlotHeader(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelAdvect::writePlotHeader" << endl;
  }

  // Setup the number of components
  HDF5HeaderData header;
  char compStr[30];

  /////////////////////////////////////////////////////////////////////////
  // Begin application-dependent code - PC.
  /////////////////////////////////////////////////////////////////////////

  // Setup the component names
  int count = 0;

  sprintf(compStr,"component_%d",count);
  header.m_string[compStr] = "phi_0J";
  count += m_numStates;

  sprintf(compStr,"component_%d",count);
  header.m_string[compStr] = "phi_0";
  count += m_numStates;

  sprintf(compStr,"component_%d",count);
  header.m_string[compStr] = "LOfPhi_0";
  count += m_numStates;

  sprintf(compStr,"component_%d",count);
  header.m_string[compStr] = "error_0";
  count += m_numStates;

  sprintf(compStr,"component_%d",count);
  header.m_string[compStr] = "J";
  count += 1;

  D_TERM(sprintf(compStr,"component_%d",count);
         header.m_string[compStr] = "xVel";,
         sprintf(compStr,"component_%d",count+1);
         header.m_string[compStr] = "yVel";,
         sprintf(compStr,"component_%d",count+2);
         header.m_string[compStr] = "zVel";)
  count += SpaceDim;

  D_TERM(sprintf(compStr,"component_%d",count);
         header.m_string[compStr] = "xLoc";,
         sprintf(compStr,"component_%d",count+1);
         header.m_string[compStr] = "yLoc";,
         sprintf(compStr,"component_%d",count+2);
         header.m_string[compStr] = "zLoc";)
  count += SpaceDim;

  header.m_int["num_components"] = count;

  /////////////////////////////////////////////////////////////////////////
  // End application-dependent code - PC.
  /////////////////////////////////////////////////////////////////////////

  // Write the header
  header.writeToFile(a_handle);

  if (s_verbosity >= 3)
  {
     pout() << header << endl;
  }
}

// Write plotfile data for this level
void AMRLevelAdvect::writePlotLevel(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelAdvect::writePlotLevel" << endl;
  }

  // if we're on level 0, write out the mapped-grid geometry info
  writeMappedPlotFile();


  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;

  a_handle.setGroup(label);

  // Setup the level header information
  HDF5HeaderData header;

  header.m_int ["ref_ratio"]   = m_ref_ratio;
  header.m_real["dx"]          = m_dx;
  header.m_real["dt"]          = m_dt;
  header.m_real["time"]        = m_time;
  header.m_box ["prob_domain"] = m_problem_domain.domainBox();

  // Write the header for this level
  header.writeToFile(a_handle);

  if (s_verbosity >= 3)
  {
    pout() << header << endl;
  }

  /////////////////////////////////////////////////////////////////////////
  // Begin application-dependent code - PC.
  /////////////////////////////////////////////////////////////////////////

  // Write the data for this level
  IntVect ivGhost = m_numGhost*IntVect::Unit;
  LevelData<FArrayBox> lOfPhi(m_phiNew.getBoxes(),m_numStates);
  LevelData<FArrayBox> phiDummy(m_phiNew.getBoxes(),m_numStates,ivGhost);
  LevelData<FArrayBox> location(m_phiNew.getBoxes(),SpaceDim);

  // Create state intervals
  Interval baseSrcInterval(0,m_numStates-1);

  int count = 0;
  Interval phiJDstInterval(count,count+m_numStates-1);
  count += m_numStates;

  Interval phiDstInterval(count,count+m_numStates-1);
  count += m_numStates;

  Interval lOfPhiDstInterval(count,count+m_numStates-1);
  count += m_numStates;

  Interval errDstInterval(count,count+m_numStates-1);
  count += m_numStates;

  Interval jDstInterval(count,count);
  count += 1;

  Interval velDstInterval(count,count+SpaceDim-1);
  count += SpaceDim;

  Interval locDstInterval(count,count+SpaceDim-1);
  count += SpaceDim;

  // Set up arguments to LevelGodunov::step based on whether there are
  // coarser and finer levels
  // Make copy of phi.
  LevelData<FArrayBox> outData(m_phiNew.getBoxes(),count);

  m_phiNew.copyTo(baseSrcInterval,outData,phiJDstInterval);

  m_phiNew.copyTo(baseSrcInterval,phiDummy,baseSrcInterval);
  const LevelData<FArrayBox>& cellAvgJ = m_coordSysPtr->getJ();
  DataIterator dit = phiDummy.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
  {
     Box intersectBox( phiDummy[dit].box() );
     intersectBox &= cellAvgJ[dit].box();
     for (int n=0; n<m_numStates; n++)
     {
        phiDummy[dit].divide( cellAvgJ[dit], 0, n, 1 );
     }
  }
  phiDummy.copyTo(baseSrcInterval,outData,phiDstInterval);

  // Undefined flux register in case we need it
  LevelFluxRegister dummyFR;

  // Set arguments to dummy values and then fix if real values are available
  LevelFluxRegister* coarserFR = &dummyFR;
  LevelFluxRegister* finerFR   = &dummyFR;

  Real tCoarserOld = 0.0;
  Real tCoarserNew = 0.0;

  LevelData<FArrayBox> dummyData(m_phiNew.getBoxes(),1);
  LevelData<FArrayBox>* coarserPhiOld = &dummyData;
  LevelData<FArrayBox>* coarserPhiNew = &dummyData;
  // A coarser level exists
  if (m_hasCoarser)
  {
    AMRLevelAdvect* coarserPtr = getCoarserLevel();

    coarserPhiOld = &coarserPtr->m_phiOld;
    coarserPhiNew = &coarserPtr->m_phiNew;

    tCoarserNew = coarserPtr->m_time;
    tCoarserOld = tCoarserNew - coarserPtr->m_dt;
  }
  LevelAdvectOperator lwo;
  if (m_hasCoarser)
    {
        lwo.define(m_grids,
                   getCoarserLevel()->m_grids,
                   m_problem_domain,
                   m_ref_ratio,
                   m_numStates,
                   RealVect(D_DECL(m_dx, m_dx, m_dx)),
                   m_hasCoarser,
                   m_hasFiner);
    }
  else
    {
        lwo.define(m_grids,
                   DisjointBoxLayout(),
                   m_problem_domain,
                   m_ref_ratio,
                   m_numStates,
                   RealVect(D_DECL(m_dx, m_dx, m_dx)),
                   m_hasCoarser,
                   m_hasFiner);
    }
  LevelData<FArrayBox>& nonConstAdvVel = const_cast<LevelData<FArrayBox>&>(m_advVel);
  lwo.spaceOrder(m_spaceOrder);
  lwo.limitFaceValues(m_limitFaceValues);
  lwo.useHyperviscosity(m_useHyperviscosity);
  lwo.hyperviscosity(m_hyperviscosity);
  lwo.setAdvectionVel(&nonConstAdvVel);
  lwo.setBC(m_basicIBCPtr);
  lwo.setCoordSys(m_coordSysPtr);

  m_phiNew.copyTo(baseSrcInterval,phiDummy,baseSrcInterval);
  lwo.evalRHS(lOfPhi,phiDummy,
              *finerFR,*coarserFR,
              *coarserPhiOld,tCoarserOld,
              *coarserPhiNew,tCoarserNew,
              m_time,m_dt);

  lOfPhi.copyTo(baseSrcInterval,outData,lOfPhiDstInterval);

  // re-use phiDummy as error
  m_basicIBCPtr->setExactOnPhysicalDomain(true);
  m_basicIBCPtr->exactSoln(phiDummy,m_problem_domain,*m_coordSysPtr,m_dx,m_time);
  dit = phiDummy.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
  {
     FArrayBox phiComputed( m_grids[dit], m_numStates );
     phiComputed.copy( m_phiNew[dit], m_grids[dit] );
     phiComputed.divide( cellAvgJ[dit], m_grids[dit], 0, 0, m_numStates );
     phiDummy[dit].minus( phiComputed, m_grids[dit], 0, 0, m_numStates);
  }
  phiDummy.copyTo(baseSrcInterval,outData,errDstInterval);

  // Jacobian LevelData here!
  cellAvgJ.copyTo(cellAvgJ.interval(), outData, jDstInterval);

  // velocity data
  m_advVel.copyTo(m_advVel.interval(), outData, velDstInterval);

  // create location LevelData here!
  dit=location.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
  {
     FArrayBox& thisLoc = location[dit];
     const Box& thisBox = thisLoc.box();
     BoxIterator bit(thisBox);
     for (bit.begin(); bit.ok(); ++bit)
     {
        IntVect iv = bit();
        RealVect xi( D_DECL( (iv[0]+0.5)*m_dx,
                             (iv[1]+0.5)*m_dx,
                             (iv[2]+0.5)*m_dx) );
        RealVect x( m_coordSysPtr->realCoord( xi ) );
        for (int icomp=0; icomp<SpaceDim; icomp++)
        {
           thisLoc(iv,icomp) = x[icomp];
        }
     }
  }
  location.copyTo(location.interval(), outData, locDstInterval);

  // write the BoxLayout and the data
  write(a_handle,outData.boxLayout());
  write(a_handle,outData,"data");

  /////////////////////////////////////////////////////////////////////////
  // End application-dependent code - PC.
  /////////////////////////////////////////////////////////////////////////
}

void
AMRLevelAdvect::writeMappedPlotFile() const
{
  // only do this on level 0
  if (m_level == 0)
    {
      // gather AMR Levels and create node-centered dataset of
      // node locations of mapped grids


      Vector<AMRLevel*> vectAMRLevels;
      {
        // cast away const for this to call this function
        AMRLevelAdvect* nonConstThis = const_cast<AMRLevelAdvect*>(this);

        vectAMRLevels = nonConstThis->getAMRLevelHierarchy();
      }

      int numLevels = vectAMRLevels.size();

      Vector<int> vectRefRatio(numLevels,0);
      Vector<DisjointBoxLayout> vectGrids(numLevels);
      Vector<LevelData<NodeFArrayBox>* > vectNodeLoc(numLevels, NULL);

      const AMRLevelAdvect*  levelPtr = this;
      // loop over levels and set things up
      for (int lev=0; lev<numLevels; lev++)
        {
          const DisjointBoxLayout& levelGrids = levelPtr->m_grids;
          vectGrids[lev] = levelPtr->m_grids;
          vectRefRatio[lev] = levelPtr->m_ref_ratio;

          vectNodeLoc[lev] = new LevelData<NodeFArrayBox>(levelGrids,
                                                         SpaceDim);

          const CoordSys<FArrayBox,FluxBox>* levelCS = levelPtr->m_coordSysPtr;
          Real levelDx = levelPtr->m_dx;

          DataIterator dit = levelGrids.dataIterator();
          for (dit.begin(); dit.ok(); ++dit)
            {
              NodeFArrayBox& thisNodeFAB = (*vectNodeLoc[lev])[dit];
              // node-centered FAB
              FArrayBox& thisFAB = thisNodeFAB.getFab();
              BoxIterator bit(thisFAB.box());
              for (bit.begin(); bit.ok(); ++bit)
                {
                  IntVect iv = bit();
                  // location in index space
                  RealVect nodeIndexLoc(iv);
                  nodeIndexLoc *= levelDx;
                  // now convert to real space location
                  RealVect nodeRealLoc = levelCS->realCoord(nodeIndexLoc);
                  D_EXPR6(thisFAB(iv,0) = nodeRealLoc[0],
                          thisFAB(iv,1) = nodeRealLoc[1],
                          thisFAB(iv,2) = nodeRealLoc[2],
                          thisFAB(iv,3) = nodeRealLoc[3],
                          thisFAB(iv,4) = nodeRealLoc[4],
                          thisFAB(iv,5) = nodeRealLoc[5]);
                }
            } // end loop over grids

          // advance to next level
          levelPtr = levelPtr->getFinerLevel();
        } // end loop over levels


      // create names
      Vector<string> locationNames(SpaceDim);
      D_EXPR6(locationNames[0] = "x",
              locationNames[1] = "y",
              locationNames[2] = "z",
              locationNames[3] = "u",
              locationNames[4] = "v",
              locationNames[5] = "w");


      // create filename
      char iter_str[80];

      sprintf(iter_str,
              "%s%06d.%dd.map.hdf5",
              m_plotfile_prefix.c_str(), m_levelStep, SpaceDim);

      string fileName(iter_str);

      // now call nodal WriteAMRHierarchy function...
      WriteAMRHierarchyHDF5(fileName,
                            vectGrids,
                            vectNodeLoc,
                            locationNames,
                            m_problem_domain.domainBox(),
                            m_dx,
                            m_dt,
                            m_time,
                            vectRefRatio,
                            numLevels);

      // now clean up our mess
      for (int lev=0; lev<vectNodeLoc.size(); lev++)
        {
          delete vectNodeLoc[lev];
        }


    }
}

#endif

////////////////////////////////////////////////////////////////////////////////

// Returns the dt computed earlier for this level
Real AMRLevelAdvect::computeDt()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelAdvect::computeDt " << m_level << endl;
  }

  FourthOrderCoordSys* FOCS = dynamic_cast<FourthOrderCoordSys*>(m_coordSysPtr);
  CH_assert(FOCS);

  // there's probably a better way to do this, but I think this should
  // work.
  // first, average advVel to faces:
  LevelData<FluxBox> faceAdvVel(m_grids, SpaceDim);
  // this shouldn't need to be 4th-order, so just do 2nd order cellToFace
  CellToEdge(m_advVel, faceAdvVel);


  Real newdt = computeMappedDt(faceAdvVel, FOCS, m_cfl);
  //Real newdt = -1;

  // backup plan if newdt is bogus...
  if (newdt < 0) newdt =  m_dt;

  return newdt;
}

////////////////////////////////////////////////////////////////////////////////

// Compute dt using initial data
Real AMRLevelAdvect::computeInitialDt()
{
    Real newDT = computeDt();
    newDT = min(m_initial_dt_multiplier * m_dx, newDT);

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelAdvect::computeInitialDt on level " << m_level << " = " << newDT << endl;
  }

  return newDT;
}

////////////////////////////////////////////////////////////////////////////////

// Set the CFL number
void AMRLevelAdvect::CFL(Real a_cfl)
{
  m_cfl = a_cfl;
}

// Set the spatial order of accuracy
void AMRLevelAdvect::spaceOrder(int a_spaceOrder)
{
  m_spaceOrder = a_spaceOrder;
}

// sets whether to limit face values in levelOperator
void AMRLevelAdvect::limitFaceValues(bool a_limitFaceValues)
{
  m_limitFaceValues = a_limitFaceValues;
}

// sets whether to enforce a min value for advected quantity
void AMRLevelAdvect::enforceMinVal(bool a_enforceMinVal, Real a_minVal)
{
  m_enforceMinVal = a_enforceMinVal;
  m_minVal = a_minVal;
}


void AMRLevelAdvect::lowOrderFluxScheme(std::string a_lowOrderFluxScheme)
{
   if (a_lowOrderFluxScheme=="CTU")
   {
      m_lowOrderFluxScheme = 1;
   }
   else if (a_lowOrderFluxScheme=="DCU")
   {
      m_lowOrderFluxScheme = 0;
   }
   else
   {
      MayDay::Error("AMRLevelAdvect::lowOrderFluxScheme: unknown scheme");
   }
}


void AMRLevelAdvect::redistributeNegativeVal(bool a_redistributeNegativeVal,
                                             int a_maxRedistributionPasses)
{
  m_redistributeNegativeVal = a_redistributeNegativeVal;
  m_maxRedistributionPasses = a_maxRedistributionPasses;
}


void AMRLevelAdvect::useHyperviscosity(bool a_useHyperviscosity)
{
   m_useHyperviscosity = a_useHyperviscosity;
}

void AMRLevelAdvect::hyperviscosity(Real a_hyperviscosity)
{
   CH_assert(a_hyperviscosity>=0);
   m_hyperviscosity = a_hyperviscosity;
}


////////////////////////////////////////////////////////////////////////////////

// Set the physical dimension of the longest side of the domain
void AMRLevelAdvect::domainLength(Real a_domainLength)
{
  m_domainLength = a_domainLength;
}

////////////////////////////////////////////////////////////////////////////////

// Set the refinement threshold
void AMRLevelAdvect::refinementThreshold(Real a_refineThresh)
{
  m_refineThresh = a_refineThresh;
}

////////////////////////////////////////////////////////////////////////////////

// Set the tag buffer size
void AMRLevelAdvect::tagBufferSize(int a_tagBufferSize)
{
  m_tagBufferSize = a_tagBufferSize;
}

////////////////////////////////////////////////////////////////////////////////
void AMRLevelAdvect::coordinateSystem(CoordSysFactory<FArrayBox,FluxBox>* a_coordSysFact)
{
  m_coordSysFactPtr = a_coordSysFact;
}

////////////////////////////////////////////////////////////////////////////////

// Create a load-balanced DisjointBoxLayout from a collection of Boxes
DisjointBoxLayout AMRLevelAdvect::loadBalance(const Vector<Box>& a_grids)
{
  // Load balance and create boxlayout
  Vector<int> procMap;
  //if (procID() == uniqueProc(SerialTask::compute))
  //{
  //  LoadBalance(procMap,a_grids);
  //}
  //broadcast(procMap,uniqueProc(SerialTask::compute));

  // appears to be faster for all procs to do the loadbalance (ndk)
  LoadBalance(procMap,a_grids);

  if (s_verbosity >= 4)
  {
    pout() << "AMRLevelAdvect::loadBalance: procesor map: " << endl;
    for (int igrid = 0; igrid < a_grids.size(); ++igrid)
    {
      pout() << igrid << ": " << procMap[igrid] << "  " << a_grids[igrid].volume() << endl;
    }
    pout() << endl;
  }

  DisjointBoxLayout dbl(a_grids,procMap,m_problem_domain);
  dbl.close();

  return dbl;
}

////////////////////////////////////////////////////////////////////////////////

// Setup menagerie of data structures
// this should be called after the coordSys has been set.
void AMRLevelAdvect::levelSetup()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelAdvect::levelSetup " << m_level << endl;
  }

  AMRLevelAdvect* amrAdvectCoarserPtr = getCoarserLevel();
  AMRLevelAdvect* amrAdvectFinerPtr   = getFinerLevel();

  m_hasCoarser = (amrAdvectCoarserPtr != NULL);
  m_hasFiner   = (amrAdvectFinerPtr   != NULL);

  if (m_hasCoarser)
  {
    int nRefCrse = m_coarser_level_ptr->refRatio();

    m_coarseAverage.define(m_grids,
                           m_numStates,
                           nRefCrse);

    m_fineInterp.define(m_grids,
                        m_numStates,
                        nRefCrse,
                        m_problem_domain);

    const DisjointBoxLayout& coarserLevelDomain = amrAdvectCoarserPtr->m_grids;

    // Maintain levelAdvectOperator
    m_levelAdvectOperator.define(m_grids,
                                 coarserLevelDomain,
                                 m_problem_domain,
                                 nRefCrse,
                                 m_numStates,
                                 RealVect(D_DECL(m_dx, m_dx, m_dx)),
                                 m_hasCoarser,
                                 m_hasFiner);

    // This may look twisted but you have to do this this way because the
    // coarser levels get setup before the finer levels so, since a flux
    // register lives between this level and the next FINER level, the finer
    // level has to do the setup because it is the only one with the
    // information at the time of construction.

    // Maintain flux registers
    amrAdvectCoarserPtr->m_fluxRegister.define(m_grids,
                                               amrAdvectCoarserPtr->m_grids,
                                               m_problem_domain,
                                               amrAdvectCoarserPtr->m_ref_ratio,
                                               m_numStates);
    amrAdvectCoarserPtr->m_fluxRegister.setToZero();
  }
  else
    {
      m_levelAdvectOperator.define(m_grids,
                                   DisjointBoxLayout(),
                                   m_problem_domain,
                                   m_ref_ratio,
                                   m_numStates,
                                   RealVect(D_DECL(m_dx, m_dx, m_dx)),
                                   m_hasCoarser,
                                   m_hasFiner);
    }
  // this should happen whether or not there's a coarser level
  m_levelAdvectOperator.spaceOrder(m_spaceOrder);
  m_levelAdvectOperator.limitFaceValues(m_limitFaceValues);
  m_levelAdvectOperator.useHyperviscosity(m_useHyperviscosity);
  m_levelAdvectOperator.hyperviscosity(m_hyperviscosity);

  // if we're enforcing a min value, then we need a PatchGodunov
  // for the CTU scheme
  if (m_enforceMinVal && m_lowOrderFluxScheme==1)
    {
      MappedAdvectionPhysics advectPhys;
      advectPhys.define(m_problem_domain, m_dx);
      advectPhys.setNComp(m_numStates);

      MappedPhysIBC physibc;
      physibc.define(m_problem_domain, m_dx);
      physibc.setBasicIBC(m_basicIBCPtr);
      CH_assert(m_coordSysPtr != NULL);
      physibc.setCoordSys(m_coordSysPtr);
      advectPhys.setPhysIBC(dynamic_cast<PhysIBC*>(&physibc));

      // defaults for doing CTU
      int normalPredOrder = 0;
      bool useFourthOrderSlopes = false;
      bool usePrimLimit = false;
      bool useCharLimit = false;
      bool useFlattening = false;
      bool useArtVisc = false;
      Real artVisc = 0.0;

      m_CTUpatchGod.define(m_problem_domain,
                           m_dx,
                           &advectPhys,
                           normalPredOrder,
                           useFourthOrderSlopes,
                           usePrimLimit,
                           useCharLimit,
                           useFlattening,
                           useArtVisc,
                           artVisc);

    }


}

////////////////////////////////////////////////////////////////////////////////

// Get the next coarser level
AMRLevelAdvect* AMRLevelAdvect::getCoarserLevel() const
{
  AMRLevelAdvect* amrAdvectCoarserPtr = NULL;

  if (m_coarser_level_ptr != NULL)
  {
    amrAdvectCoarserPtr = dynamic_cast<AMRLevelAdvect*>(m_coarser_level_ptr);

    if (amrAdvectCoarserPtr == NULL)
    {
      MayDay::Error("AMRLevelAdvect::getCoarserLevel: dynamic cast failed");
    }
  }

  return amrAdvectCoarserPtr;
}

////////////////////////////////////////////////////////////////////////////////

// Get the next finer level
AMRLevelAdvect* AMRLevelAdvect::getFinerLevel() const
{
  AMRLevelAdvect* amrAdvectFinerPtr = NULL;

  if (m_finer_level_ptr != NULL)
  {
    amrAdvectFinerPtr = dynamic_cast<AMRLevelAdvect*>(m_finer_level_ptr);

    if (amrAdvectFinerPtr == NULL)
    {
      MayDay::Error("AMRLevelAdvect::getFinerLevel: dynamic cast failed");
    }
  }

  return amrAdvectFinerPtr;
}

/// compute max- and min-preserving low-order fluxes
/** note that this takes the face-centered N^T*vel (which is the
    area-weighted normal velocity on faces)
*/
void
AMRLevelAdvect::computeLowOrderFlux(LevelData<FluxBox>& a_CTUflux,
                                    const LevelData<FArrayBox>& a_oldPhi,
                                    const LevelData<FluxBox>& a_NTvel,
                                    Real a_time,
                                    Real a_dt)
{
  FourthOrderCoordSys* FOCS = dynamic_cast<FourthOrderCoordSys*>(m_coordSysPtr);
  CH_assert(FOCS);

  // call utility function over in
  // Chombo/example/fourthOrderMappedGrids/advection
  if (m_lowOrderFluxScheme==0)
  {
     computeDCUFlux(a_CTUflux, a_oldPhi,
                    a_NTvel, FOCS,
                    m_time, a_dt);
  }
  else
  {
     computeCTUFlux(a_CTUflux, a_oldPhi,
                    a_NTvel, FOCS,
                    m_CTUpatchGod,
                    m_time, a_dt);
  }
}

#if 0
/// apply Zalesak limiter to enforce minVal
/** a_flux comes in with high-order flux and is modified in place using
    the low-order flux
*/
void
AMRLevelAdvect::applyZalesakLimiter(LevelData<FluxBox>& a_flux,
                                    const LevelData<FluxBox>& a_lowOrderFlux,
                                    const LevelData<FArrayBox>& a_oldPhi,
                                    Real a_minVal)
{
  // first test - just use low-order fluxes
  if (m_lowOrderFluxesOnly)
    {
      a_lowOrderFlux.copyTo(a_flux);
    }
  else
    {
      int ncomp = a_flux.nComp();
      CH_assert(a_oldPhi.nComp() == ncomp);
      CH_assert(a_lowOrderFlux.nComp() == ncomp);

      // actually compute limited values


      LevelData<FluxBox> antiDiffusiveFlux(m_grids,
                                           ncomp,
                                           a_flux.ghostVect());



      // don't think I need ghost cells for this guy...
      LevelData<FArrayBox> lowOrderUpdate(m_grids,
                                          ncomp);

      // implementation is simpler if we include a ghost cell on this guy
      LevelData<FArrayBox> limiterVal(m_grids, ncomp, IntVect::Unit);

      // compute low-order update
      Real fakeDx = 1.0;
      SingleLevelDivergence::levelDivergenceMAC(lowOrderUpdate,
                                                a_lowOrderFlux,
                                                fakeDx);

      Real hDinv = 1.0/(D_TERM(m_dx,*m_dx,*m_dx));
      DataIterator dit = a_oldPhi.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          Box gridBox = m_grids[dit];

          FArrayBox& thisLowOrderPhi = lowOrderUpdate[dit];
          const FArrayBox& thisOldPhi = a_oldPhi[dit];

          thisLowOrderPhi *= -m_dt;
          thisLowOrderPhi *= hDinv;
          thisLowOrderPhi += thisOldPhi;

          // also compute antidiffusive flux while we're here
          const FluxBox& thisLowOrderFlux = a_lowOrderFlux[dit];
          FluxBox& thisHighOrderFlux = a_flux[dit];

          FluxBox& thisAntiDiffusiveFlux = antiDiffusiveFlux[dit];
          thisAntiDiffusiveFlux.copy(thisHighOrderFlux);
          thisAntiDiffusiveFlux -= thisLowOrderFlux;

          // now compute sum of fluxes out of each cell
          FArrayBox sumFluxesOut(m_grids[dit], ncomp);
          sumFluxesOut.setVal(0.0);

          // will eventually want to move this to fortran
          for (int dir=0; dir<SpaceDim; dir++)
            {
              IntVect offset = BASISV(dir);
              FArrayBox& antiDiffusiveFluxDir = thisAntiDiffusiveFlux[dir];
              BoxIterator bit(gridBox);
              for (bit.begin(); bit.ok(); ++bit)
                {
                  IntVect iv = bit();
                  for (int comp=0; comp<ncomp; comp++)
                    {
                      // lowSide
                      if (antiDiffusiveFluxDir(iv,comp) < 0.0)
                        {
                          sumFluxesOut(iv,comp) -= antiDiffusiveFluxDir(iv,comp);
                        }
                      // hiSide
                      if (antiDiffusiveFluxDir(iv+offset,comp) > 0.0)
                        {
                          sumFluxesOut(iv,comp) += antiDiffusiveFluxDir(iv+offset,comp);
                        }
                    } // end loop over components
                } // end loop over cells
            } // end loop over directions

          // sumFluxesOut now contains the sum of outward-directed
          // anti-diffusivefluxes for each cell

          // now compute limiter values

          // set default value to one so that any ghost cells outside
          // valid domain have no effect
          FArrayBox& thisLimiterVal = limiterVal[dit];
          thisLimiterVal.setVal(1.0);

          Real cellVol = pow(m_dx, SpaceDim);
          //Real cellVol = m_dx;
          //#if 0
          BoxIterator bit(gridBox);
          for (bit.begin(); bit.ok(); ++bit)
            {
              IntVect iv = bit();
              for (int comp=0; comp<ncomp; comp++)
                {
                  // limiter value a la Zalesak
                  Real diff = (thisLowOrderPhi(iv,comp) - a_minVal)*cellVol;
                  if (sumFluxesOut(iv,comp) > 0)
                    {
                      thisLimiterVal(iv,comp) =min(1.0,
                                                   diff/(m_dt*sumFluxesOut(iv,comp)));
                    }
                  else
                    {
                      thisLimiterVal(iv,comp) = 0.0;
                    }
                }
            }
          //#endif
        }

      // fill in ghost cells appropriately
      limiterVal.exchange();

      // now limit fluxes
      for (dit.begin(); dit.ok(); ++dit)
        {
          const Box& gridBox = m_grids[dit];
          FArrayBox& thisLimiterVal = limiterVal[dit];
          FluxBox& thisAntiDiffusiveFlux = antiDiffusiveFlux[dit];

          for (int dir=0; dir<SpaceDim; dir++)
            {
              IntVect offset = BASISV(dir);
              FArrayBox& limitedFlux = thisAntiDiffusiveFlux[dir];
              FArrayBox& flux = a_flux[dit][dir];
              const FArrayBox& lowOrderFlux = a_lowOrderFlux[dit][dir];

              Box faceBox(gridBox);
              faceBox.surroundingNodes(dir);
              BoxIterator bit(faceBox);
              for (bit.begin(); bit.ok(); ++bit)
                {
                  for (int comp=0; comp<ncomp; comp++)
                    {
                      // for each face, compute limiter value
                      IntVect iv = bit();
                      Real C = 0;
                      // rely on sign of antidiffusive flux to determine
                      // which limit value to use.
                      Real Aloc = limitedFlux(iv,comp);
                      // If we want to enforce max as well as min, this
                      // will become a min of the two limiter values
                      // A comes from right, use high-side limiter val
                      if (Aloc < 0)
                        {
                          C = thisLimiterVal(iv, comp);
                        }
                      else
                        {
                          // use low-side limiter val
                          C=thisLimiterVal(iv-offset, comp);
                        }

                      limitedFlux(iv,comp) = C*limitedFlux(iv,comp);

                      // finally, compute final limited flux value
                      flux(iv,comp) = lowOrderFlux(iv,comp) + limitedFlux(iv,comp);
                    }
                }
            } // end loop over flux directions
        } // end loop over grid boxes
    } // end if we're doing limiting
}

#endif

/// write out error norms
void
AMRLevelAdvect::writeErrorNorms() const
{
  if (m_level == 0)
    {
      // compute error norms
      // for now, single-level; implement multilevel later
      LevelData<FArrayBox> error(m_grids, m_numStates, m_phiNew.ghostVect());
      LevelData<FArrayBox> phiComputed(m_grids, m_numStates,
                                       IntVect::Zero);
      m_basicIBCPtr->setExactOnPhysicalDomain(true);
//      m_basicIBCPtr->setExactOnPhysicalDomain(false);
      m_basicIBCPtr->exactSoln(error,m_problem_domain,*m_coordSysPtr,m_dx,m_time);
      const LevelData<FArrayBox>& cellAvgJ = m_coordSysPtr->getJ();
      DataIterator dit = m_grids.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
      {
        //FArrayBox phiComputed( m_grids[dit], m_numStates );
        phiComputed[dit].copy( m_phiNew[dit], m_grids[dit] );
        phiComputed[dit].divide( cellAvgJ[dit], m_grids[dit], 0, 0,
                                 m_numStates );
        error[dit].minus( phiComputed[dit], m_grids[dit], 0, 0, m_numStates);
        //error[dit].divide( cellAvgJ[dit], m_grids[dit], 0, 0, m_numStates );
      }

      pout() << "Error norms: " << endl;

      for (int comp=0; comp<m_numStates; comp++)
        {
          Interval compInt(comp,comp);
          DisjointBoxLayout* finerGridsPtr = NULL;
          int nRefFine = 1;
          Real L1Norm = computeNorm(error, finerGridsPtr,
                                    nRefFine, m_dx, compInt,
                                    1);

          Real L2Norm = computeNorm(error, finerGridsPtr,
                                    nRefFine, m_dx, compInt,
                                    2);

          Real maxNorm = computeNorm(error, finerGridsPtr,
                                     nRefFine, m_dx, compInt,
                                     0);

          Real maxPhi = computeMax(phiComputed, finerGridsPtr,
                                   nRefFine, compInt);

          Real minPhi = computeMin(phiComputed, finerGridsPtr,
                                   nRefFine, compInt);

          pout() << " t= " << m_time
                 << "    comp " << comp
                 << "  L1(err) = " << L1Norm
                 << "  L2(err) = " << L2Norm
                 << "  max(err) = " << maxNorm << endl;

          // report max and min vals:
          pout() << " t= " << m_time
                 << "    comp " << comp
                 << "  max = " << maxPhi
                 << "  min = " << minPhi << endl;


        } // end loop over comps
    }
}

