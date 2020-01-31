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

#include "AMRLevelShallowWaterMapped.H"
#include "LevelShallowWaterMappedOperator.H"
#include "SWintegrator.H"
#include "CubedSphere2DPanelCS.H"
#include "CubedSphere2DF_F.H"
#include "FourthOrderUtil.H"
#include "ShallowWaterPhysicsF_F.H"

#include "CH_HDF5.H"
#include "SPMD.H"
#include "AMRIO.H"
#include "computeMappedNewDt.H"
// #include "DebugInclude.H"
// for kludge
#include "PhysShallowWaterMappedIBC.H"

#include "NamespaceHeader.H"

//////////////////////////////////////////////////////////////////////////////

// Constructor
AMRLevelShallowWaterMapped::AMRLevelShallowWaterMapped()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelShallowWaterMapped default constructor" << endl;
  }
  m_levelConsOperatorPtr = new LevelShallowWaterMappedOperator();
  setDefaultValues();
}

//////////////////////////////////////////////////////////////////////////////
// Destructor - free up storage
AMRLevelShallowWaterMapped::~AMRLevelShallowWaterMapped()
{
}

//////////////////////////////////////////////////////////////////////////////

// Initialize data
void
AMRLevelShallowWaterMapped::setDataMapped(LevelData<FArrayBox>& a_U,
                                          Real a_time,
                                          bool a_includeJ) const
{
  AMRLevelMappedCons::setDataMapped(a_U, a_time, a_includeJ);

  // LevelShallowWaterMappedOperator& levelShallowWaterOperator =
  //    (LevelShallowWaterMappedOperator&) *m_levelConsOperatorPtr;
}

//////////////////////////////////////////////////////////////////////////////

/*
// Compute dt using initial data
Real
AMRLevelShallowWaterMapped::computeNewDt()
{
  LevelShallowWaterMappedOperator& levelShallowWaterOperator =
    (LevelShallowWaterMappedOperator&) *m_levelConsOperatorPtr;
  Real newDT;
  if (m_dtFromCells)
    { // the old method
      // petermc, 29 Oct 2010, multiplied by m_cfl
      LevelData<FArrayBox> cellAvgU(m_grids, m_numStates, m_ghostVect);
      // cellAvgU = m_Unew / cellAvgJ
      levelShallowWaterOperator.cellUJToCellU(cellAvgU, m_Unew);
      newDT = m_cfl * m_dx / getMaxWaveSpeed(cellAvgU);
    }
  else
    {
      LevelData<FluxBox>& advVelFace = levelShallowWaterOperator.advVelFace();
      newDT = computeMappedNewDt(advVelFace, m_coordSysPtr, m_cfl);
    }
  m_dtNew = newDT;
  return newDT;
}
*/

//////////////////////////////////////////////////////////////////////////////

int
AMRLevelShallowWaterMapped::indexForTagging()
{
  // FROM ADVECTION:  tag on the advected quantity
  return 0;
}

//////////////////////////////////////////////////////////////////////////////

void
AMRLevelShallowWaterMapped::transferSettingsToLevelOp()
{
  AMRLevelMappedCons::transferSettingsToLevelOp();
  m_mblexPtr->defineVector();
  ((LevelShallowWaterMappedOperator*) m_levelConsOperatorPtr)->defineMetricStuff();
}

////////////////////////////////////////////////////////////////////////////////

#ifdef CH_USE_HDF5

// Write plotfile header
void AMRLevelShallowWaterMapped::writePlotHeader(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelShallowWaterMapped::writePlotHeader" << endl;
  }

  // Setup the number of components
  HDF5HeaderData header;
  char compStr[30];

  /////////////////////////////////////////////////////////////////////////
  // Begin application-dependent code - PC.
  /////////////////////////////////////////////////////////////////////////

  // Setup the component names
  int comp = 0;

  // <UJ> averaged over cell

  sprintf(compStr,"component_%d",comp);
  header.m_string[compStr] = "heightJ";
  comp++;

  sprintf(compStr,"component_%d",comp);
  header.m_string[compStr] = "hvel1J";
  comp++;

  sprintf(compStr,"component_%d",comp);
  header.m_string[compStr] = "hvel2J";
  comp++;

  // <U> averaged over cell

  sprintf(compStr,"component_%d",comp);
  header.m_string[compStr] = "height";
  comp++;

  sprintf(compStr,"component_%d",comp);
  header.m_string[compStr] = "hvel1";
  comp++;

  sprintf(compStr,"component_%d",comp);
  header.m_string[compStr] = "hvel2";
  comp++;

  // Velocities.  Averages or centers?  Still 4th-order?

  sprintf(compStr,"component_%d",comp);
  header.m_string[compStr] = "vel1";
  comp++;

  sprintf(compStr,"component_%d",comp);
  header.m_string[compStr] = "vel2";
  comp++;

  sprintf(compStr,"component_%d",comp);
  header.m_string[compStr] = "velZonal";
  comp++;

  sprintf(compStr,"component_%d",comp);
  header.m_string[compStr] = "velMeridional";
  comp++;

  // relative vorticity
  sprintf(compStr,"component_%d",comp);
  header.m_string[compStr] = "zeta";
  comp++;

  // potential enstrophy
  sprintf(compStr,"component_%d",comp);
  header.m_string[compStr] = "potEnstrophy";
  comp++;

  // total energy
  sprintf(compStr,"component_%d",comp);
  header.m_string[compStr] = "E";
  comp++;

  if (s_writeJ)
    { // J
      sprintf(compStr,"component_%d",comp);
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
      // <UJ> averaged over cell
      sprintf(compStr,"component_%d",comp);
      header.m_string[compStr] = "errheightJ";
      comp++;
      sprintf(compStr,"component_%d",comp);
      header.m_string[compStr] = "errhvel1J";
      comp++;
      sprintf(compStr,"component_%d",comp);
      header.m_string[compStr] = "errhvel2J";
      comp++;

      // <U> averaged over cell
      sprintf(compStr,"component_%d",comp);
      header.m_string[compStr] = "errheight";
      comp++;
      sprintf(compStr,"component_%d",comp);
      header.m_string[compStr] = "errhvel1";
      comp++;
      sprintf(compStr,"component_%d",comp);
      header.m_string[compStr] = "errhvel2";
      comp++;
    }

  header.m_int["num_components"] = comp;

  /////////////////////////////////////////////////////////////////////////
  // End application-dependent code - PC.
  /////////////////////////////////////////////////////////////////////////

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

void
AMRLevelShallowWaterMapped::getPlotData(LevelData<FArrayBox>& a_plot_data) const
{
  CH_assert (a_plot_data.nComp() == numPlotComps());
  const DisjointBoxLayout& levelGrids = a_plot_data.getBoxes();
  DataIterator dit = levelGrids.dataIterator();
  int plot_data_counter = 0;

  /*
    <UJ> and <U>
  */

  LevelData<FArrayBox>& UnewCopy = (LevelData<FArrayBox>&) m_Unew;
  UnewCopy.exchange();
  for (dit.begin(); dit.ok(); ++dit)
    {
      a_plot_data[dit].copy(m_Unew[dit], 0, plot_data_counter, m_numStates);
    }
  plot_data_counter += m_numStates;

  LevelData<FArrayBox> cellAvgJ(levelGrids, 1, m_ghostVect);
  LevelData<FArrayBox> volume(levelGrids, 1, m_ghostVect);
  LevelData<FArrayBox> cellCenters(levelGrids, SpaceDim, m_ghostVect);
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& bxBase = levelGrids[dit];
      Box bx = grow(bxBase, m_ghostVect);

      const NewCoordSys* coordSysBlockPtr =
        m_coordSysPtr->getCoordSys(bxBase);

      FArrayBox& cellAvgJFab = cellAvgJ[dit];
      coordSysBlockPtr->getAvgJ(cellAvgJFab, bx);

      FArrayBox& cellCentersFab = cellCenters[dit];
      coordSysBlockPtr->getCenterMappedCoordinates(cellCentersFab, bx);

      FArrayBox& volumeFab = volume[dit];
      // to get volumeFab on bx, need NFlub defined on a bx grown by 1
      Box bx1 = grow(bx, 1);
      FluxBox NFlub(bx1, SpaceDim * SpaceDim);
      coordSysBlockPtr->getN(NFlub, bx1);
      coordSysBlockPtr->cellVol(volumeFab, NFlub, bx);
    }

  // Uavg = <U>
  LevelData<FArrayBox> Uavg(levelGrids, m_numStates, m_ghostVect);
  ((LevelMappedConsOperator*) m_levelConsOperatorPtr)->cellUJToCellU(Uavg, m_Unew);
  m_levelConsOperatorPtr->exchangeGhosts(Uavg);
  for (dit.begin(); dit.ok(); ++dit)
    {
      a_plot_data[dit].copy(Uavg[dit], 0, plot_data_counter, m_numStates);
    }
  plot_data_counter += m_numStates;

  /*
    Velocities
  */

  IntVect ghostVect1 = m_ghostVect - IntVect::Unit;
  IntVect ghostVect2 = ghostVect1 - IntVect::Unit;
  
  LevelData<FArrayBox> Wcen(levelGrids, m_numStates, ghostVect1);
  LevelData<FArrayBox> Wavg(levelGrids, m_numStates, ghostVect2);
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& bxBase = levelGrids[dit];
      Box bx = grow(bxBase, m_ghostVect);
      Box bx1 = grow(bxBase, ghostVect1);
      Box bx2 = grow(bxBase, ghostVect2);

      const FArrayBox& UavgFab = Uavg[dit];
      FArrayBox UcenFab(bx, m_numStates);
      UcenFab.copy(UavgFab);
      fourthOrderAverageCell(UcenFab, -1);

      FArrayBox& WcenFab = Wcen[dit];
      WcenFab.copy(UcenFab); // on bx1
      WcenFab.divide(UcenFab, UHGT, WVELX);
      WcenFab.divide(UcenFab, UHGT, WVELY);

      FArrayBox& WavgFab = Wavg[dit]; // on bx2
      WavgFab.copy(WcenFab); // on bx1
      fourthOrderAverageCell(WavgFab);
    }

  Interval WheightInterval(WHGT, WHGT);
  LevelData<FArrayBox> heightcen;
  aliasLevelData(heightcen, &Wcen, WheightInterval);

  Interval WvelInterval(WVELX, WVELY);
  LevelData<FArrayBox> velCScen;
  aliasLevelData(velCScen, &Wcen, WvelInterval);
  LevelData<FArrayBox> velCSavg;
  aliasLevelData(velCSavg, &Wavg, WvelInterval);

  for (dit.begin(); dit.ok(); ++dit)
    {
      a_plot_data[dit].copy(velCSavg[dit], 0, plot_data_counter, SpaceDim);
    }
  plot_data_counter += SpaceDim;

  // Convert velocity from equiangular to RLL.
  LevelData<FArrayBox> velRLLcen(levelGrids, SpaceDim, ghostVect1);
  // actually ghostVect2, but keep ghostVect1 to call fourthOrderAverageCell
  LevelData<FArrayBox> velRLLavg(levelGrids, SpaceDim, ghostVect1);
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& bxBase = levelGrids[dit];
      Box bx = grow(bxBase, m_ghostVect);
      Box bx1 = grow(bxBase, ghostVect1);
      Box bx2 = grow(bxBase, ghostVect2);

      const CubedSphere2DPanelCS* coordSysBlockPtr =
        (CubedSphere2DPanelCS*) (m_coordSysPtr->getCoordSys(bxBase));

      FArrayBox xiFab(bx1, SpaceDim);
      xiFab.copy(cellCenters[dit]);

      const FArrayBox& velCScenFab = velCScen[dit]; // on bx1
      FArrayBox& velRLLcenFab = velRLLcen[dit];
      coordSysBlockPtr->fabVectorTransformEquiangularToLatLon(xiFab,
                                                              velCScenFab,
                                                              velRLLcenFab);
      FArrayBox& velRLLavgFab = velRLLavg[dit]; // on bx1
      velRLLavgFab.copy(velRLLcenFab);
      fourthOrderAverageCell(velRLLavgFab);
      a_plot_data[dit].copy(velRLLavgFab, 0, plot_data_counter, SpaceDim);
    }
  plot_data_counter += SpaceDim;

  /*
    Other quantities
  */

  // zeta = relative vorticity = curl of relative wind velocity.
  // We'll take the velocity in velCScen, which has ghostVect1.
  LevelData<FArrayBox> zeta(levelGrids, 1, ghostVect2);
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& bxBase = levelGrids[dit];
      Box bx = grow(bxBase, m_ghostVect);
      Box bx1 = grow(bxBase, ghostVect1);
      Box bx2 = grow(bxBase, ghostVect2);

      const CubedSphere2DPanelCS* coordSysBlockPtr =
        (CubedSphere2DPanelCS*) (m_coordSysPtr->getCoordSys(bxBase));
      FArrayBox& zetaFab = zeta[dit];
      // const FArrayBox& velCScenFab = velCScen[dit];
      // coordSysBlockPtr->curl(velCScenFab, zetaFab, bx2);
      const FArrayBox& velRLLcenFab = velRLLcen[dit];
      coordSysBlockPtr->curlSpherical(velRLLcenFab, zetaFab, bx2);
      a_plot_data[dit].copy(zetaFab, 0, plot_data_counter, 1);
    }
  plot_data_counter++;

  Real potEnstrophyTotal = 0.;

  // potential enstrophy = (zeta + f)^2/(2*height),
  // f=2*Omega*sin(theta) is the Coriolis parameter
  // where theta is latitude (0 at equator, +pi/2 at N pole, -pi/2 at S pole)
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& bxBase = levelGrids[dit];
      Box bx = grow(bxBase, m_ghostVect);
      Box bx1 = grow(bxBase, ghostVect1);
      Box bx2 = grow(bxBase, ghostVect2);

      const CubedSphere2DPanelCS* coordSysBlockPtr =
        (CubedSphere2DPanelCS*) (m_coordSysPtr->getCoordSys(bxBase));

      FArrayBox xiFab(bx2, SpaceDim);
      xiFab.copy(cellCenters[dit]);

      FArrayBox lonlatFab(bx2, SpaceDim);
      coordSysBlockPtr->fabTransformEquiangularToLonLat(xiFab, lonlatFab);

      FArrayBox fFab(bx2, 1);
      FORT_CORIOLISPARAMETER(CHF_CONST_FRA(lonlatFab),
                             CHF_FRA1(fFab, 0),
                             CHF_BOX(bx2));

      const FArrayBox& zetaFab = zeta[dit];
      const FArrayBox& heightFab = heightcen[dit];

      FArrayBox potEnstrophyFab(bx2, 1);
      potEnstrophyFab.copy(zetaFab);
      potEnstrophyFab.plus(fFab);
      potEnstrophyFab.mult(potEnstrophyFab);
      potEnstrophyFab.divide(heightFab, 0, 0);
      potEnstrophyFab *= 0.5;

      const FArrayBox& volumeFab = volume[dit];
      FArrayBox potEnstrophyWeightedFab(bxBase, 1);
      potEnstrophyWeightedFab.copy(potEnstrophyFab);
      potEnstrophyWeightedFab.mult(volumeFab, 0, 0);
      potEnstrophyTotal += potEnstrophyWeightedFab.sum(0);

      a_plot_data[dit].copy(potEnstrophyFab, 0, plot_data_counter, 1);
    }
  plot_data_counter++;

  Real ETotal = 0.;

  // total energy E = 1/2 * h * (v dot v) + 1/2 * G * (H^2 - z^2).
  // Here H = h + z, but we don't have mountains yet, so z = 0.
  // Hence E = 1/2 * h * (v dot v) + 1/2 * G * h^2
  //         = 1/2 * h * ((v dot v) + G*h)
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& bxBase = levelGrids[dit];
      Box bx = grow(bxBase, m_ghostVect);
      Box bx1 = grow(bxBase, ghostVect1);
      Box bx2 = grow(bxBase, ghostVect2);

      FArrayBox xiFab(bx2, SpaceDim);
      xiFab.copy(cellCenters[dit]);

      const FArrayBox& heightFab = heightcen[dit];
      const FArrayBox& velFab = velCScen[dit];

      // E = 1/2 * h * ((v dot v) + G*h)

      FArrayBox vdotvFab(bx2, 1);
      FORT_CUBEDSPHERE2DDOTPROD(CHF_CONST_FRA(xiFab),
                                CHF_CONST_FRA(velFab),
                                CHF_CONST_FRA(velFab),
                                CHF_FRA1(vdotvFab, 0),
                                CHF_BOX(bx2));

      FArrayBox EFab(bx2, 1);
      FORT_SWTOTALENERGY(CHF_CONST_FRA1(vdotvFab, 0),
                         CHF_CONST_FRA1(heightFab, 0),
                         CHF_FRA1(EFab, 0),
                         CHF_BOX(bx2));

      const FArrayBox& volumeFab = volume[dit];
      FArrayBox EWeightedFab(bxBase, 1);
      EWeightedFab.copy(EFab);
      EWeightedFab.mult(volumeFab, 0, 0);
      ETotal += EWeightedFab.sum(0);

      a_plot_data[dit].copy(EFab, 0, plot_data_counter, 1);
    }
  plot_data_counter++;

#ifdef CH_MPI
  Real recvComm;
  int resultComm;

  resultComm = MPI_Allreduce(&potEnstrophyTotal, &recvComm, 1, MPI_CH_REAL,
                             MPI_SUM, Chombo_MPI::comm);
  if (resultComm != MPI_SUCCESS)
    {
      MayDay::Error("Communication error in AMRLevelShallowWaterMapped::getPlotData");
    }
  potEnstrophyTotal = recvComm;

  resultComm = MPI_Allreduce(&ETotal, &recvComm, 1, MPI_CH_REAL,
                             MPI_SUM, Chombo_MPI::comm);
  if (resultComm != MPI_SUCCESS)
    {
      MayDay::Error("Communication error in AMRLevelShallowWaterMapped::getPlotData");
    }
  ETotal = recvComm;
#endif

  pout() << setprecision(16)
         << setiosflags(ios::scientific)
         << "potential enstrophy = " << potEnstrophyTotal
         << " total energy = " << ETotal << endl;

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
          const FArrayBox& UavgFab = Uavg[dit];
          errPhiFab -= UavgFab;
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
AMRLevelShallowWaterMapped::numPlotComps() const
{
  // <UJ>, <J>, original and transformed velocities, zeta, xi, E
  int numcomp = 2*m_numStates + 2*SpaceDim + 3;
  if (s_writeJ) numcomp++;
  if (s_writeVolume) numcomp++;
  if (s_writeError) numcomp += 2*m_numStates;
  return numcomp;
}

#endif

#include "NamespaceFooter.H"
