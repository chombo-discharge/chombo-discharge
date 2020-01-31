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

#include "AMRLevelMappedCons.H"

#include "NodeFArrayBox.H"
#include "LayoutIterator.H"
#include "CH_HDF5.H"
#include "SPMD.H"
#include "BoxIterator.H"
#include "computeSum.H"
#include "CellToEdge.H"
#include "AMRIO.H"
#include "NodeAMRIO.H"
#include "WaveVelocityF_F.H"
#include "FourthOrderUtil.H"
#include "computeMappedNewDt.H"
#include "PhysMappedIBC.H"

#include "NamespaceHeader.H"

// initialize plotfile options
bool AMRLevelMappedCons::s_writeMap = true;
bool AMRLevelMappedCons::s_writeJ = true;
bool AMRLevelMappedCons::s_writeVolume = true;
bool AMRLevelMappedCons::s_writeError = false;

//////////////////////////////////////////////////////////////////////////////

// Constructor
AMRLevelMappedCons::AMRLevelMappedCons()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelMappedCons default constructor" << endl;
  }
  m_levelConsOperatorPtr = new LevelMappedConsOperator();
  setDefaultValues();
}

//////////////////////////////////////////////////////////////////////////////

// Destructor
AMRLevelMappedCons::~AMRLevelMappedCons()
{
  if (m_coordSysPtr != NULL)
    {
      delete m_coordSysPtr;
      m_coordSysPtr = NULL;
    }
  if (m_mblexPtr != NULL)
    {
      delete m_mblexPtr;
      m_mblexPtr = NULL;
    }
}

//////////////////////////////////////////////////////////////////////////////

void
AMRLevelMappedCons::setDefaultValues()
{
  AMRLevelCons::setDefaultValues();
  m_coordSysFactPtr = NULL;
  m_dtFromCells = false;
  m_levelStep = 0;
  m_useSourceTerm = false;
  m_sourceTermPtr = NULL;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedCons::useSourceTerm(bool a_useSourceTerm)
{
  m_useSourceTerm = a_useSourceTerm;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedCons::sourceTerm(const LevelSourceTerm* const a_sourceTerm)
{
  if (a_sourceTerm != NULL)
    {
      m_sourceTermPtr = a_sourceTerm->new_sourceTerm();
    }
}

//////////////////////////////////////////////////////////////////////////////

// Advance by one timestep
Real AMRLevelMappedCons::advance()
{
  ((LevelMappedConsOperator*) m_levelConsOperatorPtr)->setCoordSys(m_coordSysPtr);

  AMRLevelCons::advanceU();
  ++m_levelStep;
  Real returnDt = computeNewDt();
  return returnDt;
}


//////////////////////////////////////////////////////////////////////////////

// Set a_U, with J, at time a_time.
void AMRLevelMappedCons::setData(LevelData<FArrayBox>& a_U,
                                 Real a_time) const
{
  setDataMapped(a_U, a_time, true);
}


//////////////////////////////////////////////////////////////////////////////

// Set a_U at time a_time.
void AMRLevelMappedCons::setDataMapped(LevelData<FArrayBox>& a_U,
                                       Real a_time,
                                       bool a_includeJ) const
{
  CH_TIME("AMRLevelMappedCons::setData");
  // Begin application-dependent code - PC.

  PhysMappedIBC* physIBCPtr = (PhysMappedIBC*) m_molPhysics->getPhysIBC();
  physIBCPtr->setCoordSys(m_coordSysPtr);
  physIBCPtr->setTime(a_time);
  if (m_initialAverage)
    { // call to new function, petermc, 19 Dec 2008
      // FOR MAPPED, THIS MAY BE A NEW FUNCTION.
      // fourthOrderAverage(m_UnewOnly, m_problem_domain);
      // fourthOrderAverage(m_Unew, m_problem_domain);
      IntVect ghostExpanded = m_ghostVect + IntVect::Unit;
      LevelData<FArrayBox> Uexpanded(m_grids, m_numStates, ghostExpanded);
      if (a_includeJ)
        physIBCPtr->initializeWithJ(Uexpanded);
      else
        physIBCPtr->initialize(Uexpanded);
      const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& boundaries =
        m_coordSysPtr->boundaries();
      const Vector<Box>& mappingBlocks = m_coordSysPtr->mappingBlocks();
      // WHAT ABOUT IF THERE ARE MULTIPLE BOXES IN A CELL?  DIFFERENT?
      DataIterator dit = m_grids.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          FArrayBox& UexpandedFab = Uexpanded[dit];
          FArrayBox& UFab = a_U[dit];
          const Box& bx = UFab.box();
          // Make domainBox big enough to avoid boundary conditions,
          // EXCEPT on external boundaries.
          Box domainBox = grow(bx, 2);
          const Box& bxBase = m_grids[dit];
          int blockNum = m_coordSysPtr->whichBlock(bxBase);
          const Box& blockBox = mappingBlocks[blockNum];
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              int blockSmallEnd = blockBox.smallEnd(idir);
              if (bxBase.smallEnd(idir) == blockSmallEnd)
                {
                  const BlockBoundary& bbLo = boundaries[blockNum][idir];
                  if (bbLo.isDomainBoundary())
                    domainBox.setSmall(idir, blockSmallEnd);
                }
              int blockBigEnd = blockBox.bigEnd(idir);
              if (bxBase.bigEnd(idir) == blockBigEnd)
                {
                  const BlockBoundary& bbHi = boundaries[blockNum][idir + SpaceDim];
                  if (bbHi.isDomainBoundary())
                    domainBox.setBig(idir, blockBigEnd);
                }
            }
          ProblemDomain domain(domainBox);
          fourthOrderAverageCell(UexpandedFab, domain, bx);
          UFab.copy(UexpandedFab);
        }
    }
  else
    {
      if (a_includeJ)
        physIBCPtr->initializeWithJ(a_U);
      else
        physIBCPtr->initialize(a_U);
    }

  // No need to exchange, because sufficient ghost cells have been filled in.

  //  m_mblexPtr->interpGhosts(a_U);
  // m_levelConsOperatorPtr->exchangeGhosts(a_U);

  // End application-dependent code - PC.
}

//////////////////////////////////////////////////////////////////////////////

#ifdef CH_USE_HDF5

// Read checkpoint header
void AMRLevelMappedCons::readCheckpointHeader(HDF5Handle& a_handle)
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelMappedCons::readCheckpointHeader" << endl;
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
      MayDay::Error("AMRLevelMappedCons::readCheckpointHeader: checkpoint file does not have num_components");
    }

  int numStates = header.m_int["num_components"];
  if (numStates != m_numStates*2)
    {
      MayDay::Error("AMRLevelMappedCons::readCheckpointHeader: num_components in checkpoint file does not match solver");
    }

  // Get the component names
  std::string stateName;
  char compStr[60];
  for (int comp = 0; comp < m_numStates*2; ++comp)
    {
      sprintf(compStr,"component_%d",comp);
      if (header.m_string.find(compStr) == header.m_string.end())
        {
          MayDay::Error("AMRLevelMappedCons::readCheckpointHeader: checkpoint file does not have enough component names");
        }

      stateName = header.m_string[compStr];
      if (stateName != m_stateNames[comp])
        {
          MayDay::Error("AMRLevelMappedCons::readCheckpointHeader: state_name in checkpoint does not match solver");
        }
    }

  if (header.m_int.find("iteration") == header.m_int.end())
    {
      MayDay::Error("AMR::restart: checkpoint file does not contain iteration");
    }
  m_levelStep = header.m_int ["iteration"];
}


//////////////////////////////////////////////////////////////////////////////

// Write plotfile header
void AMRLevelMappedCons::writePlotHeader(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelMappedCons::writePlotHeader" << endl;
  }

  // Setup the number of components
  HDF5HeaderData header;
  header.m_int["num_components"] = numPlotComps();

  // Set up the component names
  int comp = 0;
  char compStr[30];
  for (int state = 0; state < m_numStates; ++state)
    { // <UJ> averaged over cell
      sprintf(compStr,"component_%d", comp);
      char stateNameChars[60];
      sprintf(stateNameChars, "%sJ", m_stateNames[state].c_str());
      header.m_string[compStr] = stateNameChars;
      comp++;
    }
  for (int state = 0; state < m_numStates; ++state)
    { // <U> averaged over cell
      sprintf(compStr,"component_%d", comp);
      header.m_string[compStr] = m_stateNames[state];
      comp++;
    }
  //  for (int comp = 0; comp < m_numStates; ++comp)
  //    { // U at cell center
  //      sprintf(compStr,"component_%d",comp + 2*m_numStates);
  //      char stateNameChars[60];
  //      sprintf(stateNameChars, "%sCen", m_stateNames[comp].c_str());
  //      header.m_string[compStr] = stateNameChars;
  //    }

  if (s_writeJ)
    { // J
      sprintf(compStr,"component_%d", comp);
      header.m_string[compStr] = "J";
      comp++;
    }

  if (s_writeVolume)
    { // volume
      sprintf(compStr,"component_%d", comp);
      header.m_string[compStr] = "volume";
      comp++;
    }

  if (s_writeError)
    {
      for (int state = 0; state < m_numStates; ++state)
        { // <UJ> averaged over cell
          sprintf(compStr,"component_%d", comp);
          char stateNameChars[60];
          sprintf(stateNameChars, "err%sJ", m_stateNames[state].c_str());
          header.m_string[compStr] = stateNameChars;
          comp++;
        }
      for (int state = 0; state < m_numStates; ++state)
        { // <U> averaged over cell
          sprintf(compStr,"component_%d", comp);
          char stateNameChars[60];
          sprintf(stateNameChars, "err%s", m_stateNames[state].c_str());
          header.m_string[compStr] = stateNameChars;
          comp++;
        }
    }

  // Write the header
  header.writeToFile(a_handle);
  a_handle.setGroup("/Expressions");
  HDF5HeaderData expressions;
  m_molPhysics->expressions(expressions);
  expressions.writeToFile(a_handle);

  if (s_verbosity >= 3)
    {
      pout() << header << endl;
    }
}

//////////////////////////////////////////////////////////////////////////////

// Write plotfile data for this level
void
AMRLevelMappedCons::writePlotLevel(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelMappedCons::writePlotLevel" << endl;
  }

  // if we're on level 0, write out the mapped-grid geometry info
  if (s_writeMap) writeMappedPlotFile();

  // Setup the level string
  char levelStr[20];
  sprintf(levelStr, "%d", m_level);
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

  const DisjointBoxLayout& levelGrids = m_Unew.getBoxes();
  int numcomp = numPlotComps();

  // AMRNavierStokes has ghost vector of IntVect::Unit,
  // then does plotData.exchange(plotData.interval()) after getPlotData(),
  // and calls write(a_handle, plotData, "data", IntVect::Unit).
  // Do this in order to avoid mismatched contours in VisIt.
  LevelData<FArrayBox> plotData(levelGrids, numcomp, IntVect::Unit);
  getPlotData(plotData);
  plotData.exchange();

  write(a_handle, levelGrids);
  write(a_handle, plotData, "data", IntVect::Unit);
}

//////////////////////////////////////////////////////////////////////////////

void
AMRLevelMappedCons::getPlotData(LevelData<FArrayBox>& a_plot_data) const
{
  CH_assert (a_plot_data.nComp() == numPlotComps());
  const DisjointBoxLayout& levelGrids = a_plot_data.getBoxes();
  DataIterator dit = levelGrids.dataIterator();

  // Create state intervals
  Interval baseSrcInterval(0,m_numStates-1);

  int plot_data_counter = 0;

  LevelData<FArrayBox>& UnewCopy = (LevelData<FArrayBox>&) m_Unew;
  UnewCopy.exchange();
  // for <UJ> == m_Unew
  for (dit.begin(); dit.ok(); ++dit)
    {
      a_plot_data[dit].copy(m_Unew[dit], 0, plot_data_counter, m_numStates);
    }
  plot_data_counter += m_numStates;

  // const LevelData<FArrayBox>& cellAvgJ = m_coordSysPtr->getJ();
  LevelData<FArrayBox> cellAvgJ(levelGrids, 1, m_ghostVect);
  LevelData<FArrayBox> volume(levelGrids, 1, m_ghostVect);
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& bxBase = levelGrids[dit];
      Box bx = grow(bxBase, m_ghostVect);

      const NewCoordSys* coordSysBlockPtr =
        m_coordSysPtr->getCoordSys(bxBase);

      FArrayBox& cellAvgJFab = cellAvgJ[dit];
      coordSysBlockPtr->getAvgJ(cellAvgJFab, bx);

      FArrayBox& volumeFab = volume[dit];
      // to get volumeFab on bx, need NFlub defined on a bx grown by 1
      Box bx1 = grow(bx, 1);
      FluxBox NFlub(bx1, SpaceDim * SpaceDim);
      coordSysBlockPtr->getN(NFlub, bx1);
      coordSysBlockPtr->cellVol(volumeFab, NFlub, bx);
    }

  // phiDummy = m_Unew / cellAvgJ
  LevelData<FArrayBox> phiDummy(levelGrids, m_numStates, m_ghostVect);
  ((LevelMappedConsOperator*) m_levelConsOperatorPtr)->cellUJToCellU(phiDummy, m_Unew);
  // for <U>
  for (dit.begin(); dit.ok(); ++dit)
    {
      a_plot_data[dit].copy(phiDummy[dit], 0, plot_data_counter, m_numStates);
    }
  plot_data_counter += m_numStates;

  if (s_writeJ)
    {
      // Jacobian LevelData here!
      // for <J>
      for (dit.begin(); dit.ok(); ++dit)
        {
          a_plot_data[dit].copy(cellAvgJ[dit], 0, plot_data_counter, 1);
        }
      plot_data_counter++;
    }

  if (s_writeVolume)
    {
      // Volume LevelData here!
      // for volume of each cell
      for (dit.begin(); dit.ok(); ++dit)
        {
          a_plot_data[dit].copy(volume[dit], 0, plot_data_counter, 1);
        }
      plot_data_counter++;
    }

  if (s_writeError)
    {
      LevelData<FArrayBox> errPhiJ(levelGrids, m_numStates, m_ghostVect);
      // Initialize errPhiJ with exact solution.
      setDataMapped(errPhiJ, m_time, true);
      for (dit.begin(); dit.ok(); ++dit)
        {
          FArrayBox& errPhiJFab = errPhiJ[dit];
          const FArrayBox& UnewFab = m_Unew[dit];
          errPhiJFab -= UnewFab;
          errPhiJFab.negate();
          a_plot_data[dit].copy(errPhiJFab, 0, plot_data_counter, m_numStates);
        }
      plot_data_counter += m_numStates;

      LevelData<FArrayBox> errPhi(levelGrids, m_numStates, m_ghostVect);
      // Initialize errPhi with exact solution.
      setDataMapped(errPhi, m_time, false);
      for (dit.begin(); dit.ok(); ++dit)
        {
          FArrayBox& errPhiFab = errPhi[dit];
          const FArrayBox& phiDummyFab = phiDummy[dit];
          errPhiFab -= phiDummyFab;
          errPhiFab.negate();
          a_plot_data[dit].copy(errPhiFab, 0, plot_data_counter, m_numStates);
        }
      plot_data_counter += m_numStates;
    }

  // Now done with plot_data_counter
  CH_assert(plot_data_counter == numPlotComps());
}

//////////////////////////////////////////////////////////////////////////////

int
AMRLevelMappedCons::numPlotComps() const
{
  int numcomp = 2*m_numStates;
  if (s_writeJ) numcomp++;
  if (s_writeVolume) numcomp++;
  if (s_writeError) numcomp += 2*m_numStates;
  return numcomp;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedCons::writeMappedPlotFile() const
{
  int realDim = SpaceDim;
  // only do this on level 0
  if (m_level == 0)
    {
      // gather AMR Levels and create node-centered dataset of
      // node locations of mapped grids

      Vector<AMRLevel*> vectAMRLevels;
      {
        // cast away const for this to call this function
        AMRLevelMappedCons* nonConstThis = const_cast<AMRLevelMappedCons*>(this);

        vectAMRLevels = nonConstThis->getAMRLevelHierarchy();
      }

      int numLevels = vectAMRLevels.size();

      Vector<int> vectRefRatio(numLevels,0);
      Vector<DisjointBoxLayout> vectGrids(numLevels);
      Vector<LevelData<NodeFArrayBox>* > vectNodeLoc(numLevels, NULL);

      const AMRLevelMappedCons* levelPtr = this;
      // loop over levels and set things up
      for (int lev=0; lev<numLevels; lev++)
        {
          const DisjointBoxLayout& levelGrids = levelPtr->m_grids;
          vectGrids[lev] = levelPtr->m_grids;
          vectRefRatio[lev] = levelPtr->m_ref_ratio;

          const MultiBlockCoordSys* levelCS = levelPtr->m_coordSysPtr;
          realDim = levelCS->realDim();
          // Real levelDx = levelPtr->m_dx;
          vectNodeLoc[lev] = new LevelData<NodeFArrayBox>(levelGrids,
                                                          realDim,
                                                          IntVect::Unit);
          DataIterator dit = levelGrids.dataIterator();
          for (dit.begin(); dit.ok(); ++dit)
            {
              const Box& bxBase = levelGrids[dit];

              const NewCoordSys* levelCSblock =
                levelCS->getCoordSys(bxBase);

              const RealVect& levelDx = levelCSblock->dx();

              NodeFArrayBox& thisNodeFAB = (*vectNodeLoc[lev])[dit];
              // node-centered FAB
              FArrayBox& thisFAB = thisNodeFAB.getFab();
              const Box& thisBox = thisFAB.box();
              FArrayBox xiFab(thisBox, SpaceDim);
              BoxIterator bit(thisBox);
              for (bit.begin(); bit.ok(); ++bit)
                {
                  IntVect iv = bit();
                  // location in index space
                  RealVect nodeIndexLoc(iv);
                  nodeIndexLoc *= levelDx;
                  // now convert to real space location
                  // RealVect nodeRealLoc = levelCSblock->realCoord(nodeIndexLoc);
                  D_EXPR6(xiFab(iv,0) = nodeIndexLoc[0],
                          xiFab(iv,1) = nodeIndexLoc[1],
                          xiFab(iv,2) = nodeIndexLoc[2],
                          xiFab(iv,3) = nodeIndexLoc[3],
                          xiFab(iv,4) = nodeIndexLoc[4],
                          xiFab(iv,5) = nodeIndexLoc[5]);
                }
              levelCSblock->realCoord(thisFAB, xiFab, thisBox);
            } // end loop over grids

          // advance to next level
          levelPtr = (AMRLevelMappedCons*) levelPtr->getFinerLevel();
        } // end loop over levels

      // create names
      Vector<string> locationNames(realDim);
      if (realDim > 0) locationNames[0] = "x";
      if (realDim > 1) locationNames[1] = "y";
      if (realDim > 2) locationNames[2] = "z";
      if (realDim > 3) locationNames[3] = "u";
      if (realDim > 4) locationNames[4] = "v";
      if (realDim > 5) locationNames[5] = "w";

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

//////////////////////////////////////////////////////////////////////////////

// Compute dt using initial data
Real AMRLevelMappedCons::computeInitialDt()
{
  // Real newDT = m_initial_dt_multiplier * computeNewDt();
  // petermc, 16 Oct 2009:
  // It looks funny, but this is what DanM's AMRLevelAdvect does.
  // petermc, 29 Oct 2010:  Hm, newDT now has an extra factor of m_cfl.
  Real newDT = computeNewDt();
  newDT = min(m_initial_dt_multiplier * m_dx, newDT);

  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelMappedCons::computeInitialDt on level " << m_level << " = " << newDT << endl;
    }

  return newDT;
}

//////////////////////////////////////////////////////////////////////////////

// Compute dt using initial data
Real AMRLevelMappedCons::computeNewDt()
{
  LevelData<FArrayBox> cellAvgU(m_grids, m_numStates, m_ghostVect);
  // cellAvgU = m_Unew / cellAvgJ
  ((LevelMappedConsOperator*) m_levelConsOperatorPtr)->cellUJToCellU(cellAvgU, m_Unew);
  Real newDT;
  if (m_dtFromCells)
    { // the old method
      // petermc, 29 Oct 2010, multiplied by m_cfl
      newDT = m_cfl * m_dx / getMaxWaveSpeed(cellAvgU);
    }
  else
    {
      LevelData<FArrayBox> cellAvgW(m_grids, m_numStates, m_ghostVect);
      DataIterator dit = m_grids.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          const FArrayBox& UcellFab = cellAvgU[dit];
          FArrayBox& WcellFab = cellAvgW[dit];
          Box bx = WcellFab.box();
          m_molPhysics->consToPrim(WcellFab, UcellFab, bx);
        }

      // cellVel:  cell-averaged velocities
      LevelData<FArrayBox> cellVel;
      Interval velInt = m_molPhysics->velocityInterval();
      aliasLevelData(cellVel, &cellAvgW, velInt);

      // cellWaveVel:  velocity +/- speed of sound
      LevelData<FArrayBox> cellWaveVel(m_grids, SpaceDim, m_ghostVect);
      for (dit.begin(); dit.ok(); ++dit)
        {
          const FArrayBox& UcellFab = cellAvgU[dit];
          Box bx = UcellFab.box();
          FArrayBox soundSpeedFab(bx, 1);
          m_molPhysics->soundSpeed(soundSpeedFab, UcellFab, bx);
          FArrayBox& cellVelFab = cellVel[dit];
          FArrayBox& cellWaveVelFab = cellWaveVel[dit];
          FORT_WAVEVELOCITY(CHF_FRA(cellWaveVelFab),
                            CHF_CONST_FRA(cellVelFab),
                            CHF_CONST_FRA1(soundSpeedFab, 0),
                            CHF_BOX(bx));
        }

      // faceVel:  face-averaged wave velocities
      LevelData<FluxBox> faceWaveVel(m_grids, SpaceDim, m_ghostVect);
      // second-order is sufficient
      CellToEdge(cellWaveVel, faceWaveVel);
      // seems I need m_dx factor here
      newDT = computeMappedNewDt(faceWaveVel, m_coordSysPtr, m_cfl);
    }
  m_dtNew = newDT;
  return newDT;
}

////////////////////////////////////////////////////////////////////////////////
void AMRLevelMappedCons::coordinateSystem(MultiBlockCoordSysFactory* a_coordSysFact)
{
  m_coordSysFactPtr = a_coordSysFact;
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedCons::transferSettingsToLevelOp()
{
  // This should happen whether or not there's a coarser level.
  AMRLevelCons::transferSettingsToLevelOp();

  // It's kind of inelegant to have the defining of
  // m_coordSysPtr and m_geom and m_mblexPtr
  // to occur in a function called transferSettingsToLevelOp(),
  // but doing so allows me to avoid having a separate function
  // AMRLevelMappedCons::levelSetup().

  // Create the coordinate system
  RealVect dxVect = m_dx * RealVect::Unit;
  // m_coordSysPtr = (NewFourthOrderCoordSys*)
  //   m_coordSysFactPtr->getCoordSys(m_problem_domain,
  //                                  dxVect);
  m_coordSysPtr = m_coordSysFactPtr->getCoordSys(m_problem_domain,
                                                 dxVect);
  // MultiBlockLevelGeom m_geom;
  m_geom.define(m_coordSysPtr, m_grids, m_numGhost + 1);
  // MultiBlockLevelExchangeAverage* m_mblexPtr;
  m_mblexPtr = new MultiBlockLevelExchangeAverage();
  m_mblexPtr->define(&m_geom, m_numGhost, m_spaceOrder);

  // Finally do what this function is supposed to do.
  LevelMappedConsOperator& levelConsOperator =
    (LevelMappedConsOperator&) *m_levelConsOperatorPtr;
  levelConsOperator.setCoordSys(m_coordSysPtr);
  levelConsOperator.setLevelExchange(m_mblexPtr);

  // Is this in the right place?
  if (m_sourceTermPtr != NULL)
    m_sourceTermPtr->define(m_coordSysPtr, m_molPhysics, m_grids);
  levelConsOperator.useSourceTerm(m_useSourceTerm);
  // If using source term, transfer m_sourceTermPtr.
  levelConsOperator.setSourceTerm(m_sourceTermPtr);
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedCons::transferSettingsFromAMRLevel(AMRLevelMappedCons* a_amrConsPtr)
{
  AMRLevelCons::transferSettingsFromAMRLevel(a_amrConsPtr);
  m_dtFromCells = a_amrConsPtr->m_dtFromCells;
}

#include "NamespaceFooter.H"
