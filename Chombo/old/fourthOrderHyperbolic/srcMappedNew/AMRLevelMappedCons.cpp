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
#include "AdvectOpF_F.H"
// #include "Divergence.H"
#include "PolytropicPhysicsF_F.H"
#include "FourthOrderUtil.H"
#include "computeMappedNewDt.H"
// #include "DebugInclude.H"
// for kludge
#include "PhysMappedIBC.H"

#include "NamespaceHeader.H"

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
}

//////////////////////////////////////////////////////////////////////////////

void
AMRLevelMappedCons::setDefaultValues()
{
  AMRLevelCons::setDefaultValues();
  m_coordSysFactPtr = NULL;
  m_dtFromCells = false;
  m_levelStep = 0;
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

// Initialize data
void AMRLevelMappedCons::initialData()
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelMappedCons::initialData " << m_level << endl;
    }
  // Begin application-dependent code - PC.

  PhysMappedIBC* physIBCPtr = (PhysMappedIBC*) m_molPhysics->getPhysIBC();
  // Kludge!  Need this for mapped IBC.
  // ((GaussianMappedIBC*) physIBCPtr)->setCoordSys(m_coordSysPtr);
  physIBCPtr->setCoordSys(m_coordSysPtr);
  physIBCPtr->setTime(m_time);
  physIBCPtr->initialize(m_Unew);

  //  m_Unew.exchange();
  if (m_initialAverage)
    { // call to new function, petermc, 19 Dec 2008
      // FOR MAPPED, THIS MAY BE A NEW FUNCTION.
      fourthOrderAverage(m_Unew, m_problem_domain);
    }
  //  m_mblexPtr->interpGhosts(m_Unew);
  m_levelConsOperatorPtr->exchangeGhosts(m_Unew);

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
  header.m_int["num_components"] = 2*m_numStates + 2; // or 3*m_numStates+2 if Ucen written out

  // Set up the component names
  char compStr[30];
  for (int comp = 0; comp < m_numStates; ++comp)
    { // <UJ> averaged over cell
      sprintf(compStr,"component_%d",comp);
      char stateNameChars[60];
      sprintf(stateNameChars, "%sJ", m_stateNames[comp].c_str());
      header.m_string[compStr] = stateNameChars;
    }
  for (int comp = 0; comp < m_numStates; ++comp)
    { // <U> averaged over cell
      sprintf(compStr,"component_%d",comp + m_numStates);
      header.m_string[compStr] = m_stateNames[comp];
    }
  //  for (int comp = 0; comp < m_numStates; ++comp)
  //    { // U at cell center
  //      sprintf(compStr,"component_%d",comp + 2*m_numStates);
  //      char stateNameChars[60];
  //      sprintf(stateNameChars, "%sCen", m_stateNames[comp].c_str());
  //      header.m_string[compStr] = stateNameChars;
  //    }
  { // J
    int comp = 2*m_numStates; // or 3*m_numStates if including Ucen
    sprintf(compStr,"component_%d",comp);
    header.m_string[compStr] = "J";

    // volume
    comp++;
    sprintf(compStr,"component_%d",comp);
    header.m_string[compStr] = "volume";
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
void AMRLevelMappedCons::writePlotLevel(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelMappedCons::writePlotLevel" << endl;
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

  // Create state intervals
  Interval baseSrcInterval(0,m_numStates-1);

  int count = 0;
  Interval phiJDstInterval(count,count+m_numStates-1);
  count += m_numStates;

  Interval phiDstInterval(count,count+m_numStates-1);
  count += m_numStates;

  //  Interval phiCenInterval(count,count+m_numStates-1);
  //  count += m_numStates;

  Interval jDstInterval(count,count);
  count += 1;

  Interval volDstInterval(count,count);
  count += 1;

  //  Interval velDstInterval(count,count+SpaceDim-1);
  //  count += SpaceDim;
  //
  //  Interval locDstInterval(count,count+SpaceDim-1);
  //  count += SpaceDim;
  LevelData<FArrayBox>& UnewCopy = (LevelData<FArrayBox>&) m_Unew;
  UnewCopy.exchange();
  // Set up arguments to LevelGodunov::step based on whether there are
  // coarser and finer levels
  const DisjointBoxLayout& layout = m_Unew.getBoxes();
  LevelData<FArrayBox> phiDummy(layout, m_numStates, m_ghostVect);
  // Make copy of Unew, which is <UJ>.
  LevelData<FArrayBox> outData(layout, count);
  // from baseSrcInterval to phiJDstInterval:  <UJ>
  m_Unew.copyTo(baseSrcInterval, outData, phiJDstInterval);

  // from baseSrcInterval to baseSrcInterval
  //  m_Unew.copyTo(baseSrcInterval, phiDummy, baseSrcInterval);
  // const LevelData<FArrayBox>& cellAvgJ = m_coordSysPtr->getJ();
  LevelData<FArrayBox> cellAvgJ(layout, 1, m_ghostVect);
  LevelData<FArrayBox> volume(layout, 1, m_ghostVect);
  DataIterator dit = layout.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      Box bx = grow(layout[dit], m_ghostVect);

      FArrayBox& cellAvgJFab = cellAvgJ[dit];
      m_coordSysPtr->getAvgJ(cellAvgJFab, bx);

      FArrayBox& volumeFab = volume[dit];
      // to get volumeFab on bx, need NFlub defined on a bx grown by 1
      Box bx1 = grow(bx, 1);
      FluxBox NFlub(bx1, SpaceDim * SpaceDim);
      m_coordSysPtr->getN(NFlub, bx1);
      m_coordSysPtr->cellVol(volumeFab, NFlub, bx);
    }
  // Jacobian LevelData here!
  cellAvgJ.copyTo(cellAvgJ.interval(), outData, jDstInterval);

  // Volume LevelData here!
  volume.copyTo(volume.interval(), outData, volDstInterval);

  // phiDummy = m_Unew / cellAvgJ
  ((LevelMappedConsOperator*) m_levelConsOperatorPtr)->cellUJToCellU(phiDummy, m_Unew);
  // phiDstInterval to hold <U>
  phiDummy.copyTo(baseSrcInterval, outData, phiDstInterval);

  //  // Now fill phiDummy with cell-centered U.
  //  fourthOrderAverage(phiDummy, m_problem_domain, -1);
  //  // phiCenInterval to hold U at cell centers
  //  phiDummy.copyTo(baseSrcInterval, outData, phiCenInterval);

  // write the BoxLayout and the data
  write(a_handle, layout);
  write(a_handle, outData, "data");
}


//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedCons::writeMappedPlotFile() const
{
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

          vectNodeLoc[lev] = new LevelData<NodeFArrayBox>(levelGrids,
                                                          SpaceDim);

          const NewCoordSys* levelCS = levelPtr->m_coordSysPtr;
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
          levelPtr = (AMRLevelMappedCons*) levelPtr->getFinerLevel();
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

//////////////////////////////////////////////////////////////////////////////

// Compute dt using initial data
Real AMRLevelMappedCons::computeInitialDt()
{
  //  Real newDT = m_initial_dt_multiplier * m_dx / getMaxWaveSpeed(m_Unew);
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
  Real newDT;
  if (m_dtFromCells)
    { // the old method
      newDT = m_dx / getMaxWaveSpeed(m_Unew);
    }
  else
    {
      LevelData<FArrayBox> cellAvgU(m_grids, m_numStates, m_ghostVect);
      // cellAvgU = m_Unew / cellAvgJ
      ((LevelMappedConsOperator*) m_levelConsOperatorPtr)->cellUJToCellU(cellAvgU, m_Unew);

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
void AMRLevelMappedCons::coordinateSystem(NewCoordSysFactory* a_coordSysFact)
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
  m_coordSysPtr = (NewFourthOrderCoordSys*)
    m_coordSysFactPtr->getCoordSys(m_problem_domain,
                                   dxVect);

  // Finally do what this function is supposed to do.
  LevelMappedConsOperator& levelConsOperator =
    (LevelMappedConsOperator&) *m_levelConsOperatorPtr;
  levelConsOperator.setCoordSys(m_coordSysPtr);
}

//////////////////////////////////////////////////////////////////////////////

void AMRLevelMappedCons::transferSettingsFromAMRLevel(AMRLevelMappedCons* a_amrConsPtr)
{
  AMRLevelCons::transferSettingsFromAMRLevel(a_amrConsPtr);
  m_dtFromCells = a_amrConsPtr->m_dtFromCells;
}

#include "NamespaceFooter.H"
