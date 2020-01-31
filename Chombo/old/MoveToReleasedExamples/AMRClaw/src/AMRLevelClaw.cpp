#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "AMRLevelClaw.H"

#include "CH_HDF5.H"
#include "SPMD.H"
#include "LoadBalance.H"
#include "BoxIterator.H"
#include "computeSum.H"

// Constructor
AMRLevelClaw::AMRLevelClaw()
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelClaw default constructor" << endl;
    }

  m_cfl = 0.8;
  m_domainLength = 1.0;
  m_refineThresh = 0.2;
  m_initial_dt_multiplier = 0.1;
  m_coarse_cfl = 0.0;

  // Default values for Clawpack parameters not set...
}


// Destructor
AMRLevelClaw::~AMRLevelClaw()
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelClaw destructor" << endl;
    }
}


void AMRLevelClaw::clawPatch(const ClawPatch& a_clawPatch)
{
  m_clawPatch.define(a_clawPatch);
}

// Define new AMR level
void AMRLevelClaw::define(AMRLevel*            a_coarserLevelPtr,
                          const ProblemDomain& a_problemDomain,
                          int                  a_level,
                          int                  a_refRatio)

{

  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelClaw::define " << a_level << endl;
    }

  // Call inherited define
  AMRLevel::define(a_coarserLevelPtr,
                   a_problemDomain,
                   a_level,
                   a_refRatio);

  // Get setup information from the next coarser level
  if (a_coarserLevelPtr != NULL)
    {
      AMRLevelClaw* amrClawPtr = dynamic_cast<AMRLevelClaw*>(a_coarserLevelPtr);

      if (amrClawPtr != NULL)
        {
          m_cfl = amrClawPtr->m_cfl;
          m_domainLength = amrClawPtr->m_domainLength;
          m_refineThresh = amrClawPtr->m_refineThresh;
          m_tagBufferSize = amrClawPtr->m_tagBufferSize;
        }
      else
        {
          MayDay::Error("AMRLevelClaw::define: a_coarserLevelPtr is not castable to AMRLevelClaw*");
        }
    }

  CH_assert(m_clawPatch.isDefined());
  // Compute the grid spacing
  m_dx = m_domainLength / a_problemDomain.domainBox().longside();

  m_numGhost   = m_clawPatch.get_mbc();
  m_numStates  = m_clawPatch.get_meqn();
  m_stateNames = m_clawPatch.get_stateNames();
  m_maux = m_clawPatch.get_maux();
  m_mcapa = m_clawPatch.get_mcapa();
  m_hasKappa = m_mcapa > 0;
}

// Advance by one timestep
Real AMRLevelClaw::advance()
{

  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelClaw::advance level " << m_level << " to time " << m_time << endl;
    }

  // Copy the new to the old
  m_phiNew.copyTo(m_phiNew.interval(),
                  m_phiOld,
                  m_phiOld.interval());

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
  // A finer level exists
  if (m_hasFiner)
    {
      // Recall that my flux register goes between my level and the next
      // finer level
      finerFR = &m_fluxRegister;
      finerFR->setToZero();
    }

  // Fill phi's ghost cells using fillInterp
  if (m_hasCoarser)
    {
      // A coarser level exists
      AMRLevelClaw* coarserPtr = getCoarserLevel();

      // Recall that my flux register goes between my level and the next
      // finer level
      coarserFR = &coarserPtr->m_fluxRegister;

      coarserPhiOld = &coarserPtr->m_phiOld;
      coarserPhiNew = &coarserPtr->m_phiNew;

      tCoarserNew = coarserPtr->m_time;
      tCoarserOld = tCoarserNew - coarserPtr->m_dt;

      // Fraction "m_time" falls between the old and the new coarse times
      Real alpha = (m_time - tCoarserOld) / (tCoarserNew - tCoarserOld);

      // Truncate the fraction to the range [0,1] to remove floating-point
      // subtraction roundoff effects
      Real eps = 0.0001 * m_dt / m_ref_ratio;

      if (Abs(alpha) < eps)
        {
          alpha = 0.0;
        }

      if (Abs(1.0-alpha) < eps)
        {
          alpha = 1.0;
        }

      // Current time before old coarse time
      if (alpha < 0.0)
        {
          MayDay::Error( "AMRLevelClaw::advance: alpha < 0.0");
        }

      // Current time after new coarse time
      if (alpha > 1.0)
        {
          MayDay::Error( "AMRLevelClaw::advance: alpha > 1.0");
        }

      // Interpolate ghost cells from next coarser level using both space
      // and time interpolation
      LevelData<FArrayBox>& phiCoarseOld = coarserPtr->m_phiOld;
      LevelData<FArrayBox>& phiCoarseNew = coarserPtr->m_phiNew;
      m_patcher.fillInterp(m_phiNew,
                           phiCoarseOld,
                           phiCoarseNew,
                           alpha,
                           0,0,m_numStates);
    } // done with m_hasCoarser

  // Exchange all the data between grids at this level
  m_phiNew.exchange(m_phiNew.interval());

  // Real maxWaveSpeed = 0.;
  Real maxcfl = 0;

  DataIterator dit = m_phiNew.disjointBoxLayout().dataIterator();

  // Get disjoint box layout from m_phiNew.
  DisjointBoxLayout dblfine = m_phiNew.disjointBoxLayout();
  DisjointBoxLayout dblcoarse;
  DisjointBoxLayout dblfine_aux = m_aux.disjointBoxLayout();
  DisjointBoxLayout dblcoarse_aux;
  LevelData<FArrayBox> phiTmpCoarse, auxTmpCoarse;
  int refRatio;
  if (m_hasCoarser)
    {
      AMRLevelClaw* coarserPtr = getCoarserLevel();
      refRatio = coarserPtr->refRatio();

      // First get coarse q data for qad
      LevelData<FArrayBox>& phiCoarseOld = coarserPtr->m_phiOld;
      coarsen(dblcoarse,dblfine,refRatio);
      phiTmpCoarse.define(dblcoarse,m_numStates,IntVect::Unit);
      Interval interv(0,m_numStates-1);
      // copyTo doesn't copy ghost cell data in non-valid data regions
      phiCoarseOld.copyTo(interv,phiTmpCoarse,interv);

      // now get coarse aux data for qad
      LevelData<FArrayBox>& auxCoarseOld = coarserPtr->m_aux;
      coarsen(dblcoarse_aux,dblfine_aux,refRatio);
      int maux = 1;
      if (m_maux > 0)
        {
          maux = m_maux;
        }
      auxTmpCoarse.define(dblcoarse_aux,maux,IntVect::Unit);
      Interval interv_aux(0,maux-1);
      // copyTo doesn't copy ghost cell data in non-valid data regions
      auxCoarseOld.copyTo(interv_aux,auxTmpCoarse,interv_aux);
    }
  else
    {
      // Don't coarsen if this is the coarsest level - we'll hopefully handle this correctly
      // in ClawPatch;  extra layer of ghost cells won't get used.
      refRatio = 1;

      phiTmpCoarse.define(dblfine,m_numStates,IntVect::Unit);
      Interval interv(0,m_numStates-1);
      m_phiNew.copyTo(interv,phiTmpCoarse,interv);

      // Aux data
      int maux = 1;
      if (m_maux > 0)
        {
          maux = m_maux;
        }
      auxTmpCoarse.define(dblfine_aux,maux,IntVect::Unit);
      Interval interv_aux(0,maux-1);
      m_aux.copyTo(interv_aux,auxTmpCoarse,interv_aux);
    }

  for (dit.begin();dit.ok();++dit)
    {
      FArrayBox& phiPatch = m_phiNew[dit()];
      FArrayBox& auxPatch = m_aux[dit()];
      FArrayBox& phiCoarse = phiTmpCoarse[dit()];
      FArrayBox& auxCoarse = auxTmpCoarse[dit()];

      Box patchBox = m_grids.get(dit());
      FArrayBox fluxp[SpaceDim]; // This will be divided by kapp at this level
      FArrayBox fluxm[SpaceDim];
      FArrayBox fluxpc[SpaceDim]; // This will be divided by coarser kappa
      FArrayBox fluxmc[SpaceDim];
      FArrayBox qadd[SpaceDim];
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          fluxp[idir].resize(surroundingNodes(patchBox,idir),m_numStates);
          fluxm[idir].resize(surroundingNodes(patchBox,idir),m_numStates);
          fluxpc[idir].resize(surroundingNodes(patchBox,idir),m_numStates);
          fluxmc[idir].resize(surroundingNodes(patchBox,idir),m_numStates);
          qadd[idir].resize(surroundingNodes(patchBox,idir),m_numStates);
        }

      //  Call patch integrator.

      // Real maxWaveSpeedLocal =
      Real cflgrid =
        m_clawPatch.ClawPatchIntegrator(phiPatch,auxPatch,fluxp,fluxm,fluxpc, fluxmc,
                                        phiCoarse,auxCoarse,qadd,
                                        patchBox,m_problem_domain,m_time,
                                        m_dt,m_dx,refRatio,m_level);

      // printf("Courant number of grid at level %d is %8.4f\n",m_level, cflgrid);

      // maxWaveSpeed = std::max(maxWaveSpeed,maxWaveSpeedLocal);
      maxcfl = Max(maxcfl, cflgrid);

      // Do flux register updates.

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          // Increment coarse flux register between this level and the next
          // finer level - this level is the next coarser level with respect
          // to the next finer level
          if (m_hasFiner)
            {

              // increment coarse registers between this level and next finer level
              // This one gets divided by kappa at this level, if mcapa > 0
              // fm should get divided by kappa(i-1), (or kappa(j-1) or kappa(k-1))
              // fp should get divided by kappa(i), kappa(j), or kappa(k)
              finerFR->incrementCoarse(fluxm[idir],m_dt,dit(),
                                       phiPatch.interval(),
                                       phiPatch.interval(),idir,Side::Lo);
              finerFR->incrementCoarse(fluxp[idir],m_dt,dit(),
                                       phiPatch.interval(),
                                       phiPatch.interval(),idir,Side::Hi);

            }

          // Increment fine flux registers between this level and the next
          // coarser level - this level is the next finer level with respect
          // to the next coarser level
          if (m_hasCoarser)
            {
              // increment fine registers between this level and next coarser level
              // This one gets divided by kappa at next coarser level.  How do we do this?
              // fm should get divided by m_coarsened_kappa(i-1), etc
              // fp should get divided by m_coarsened_kappa(i), etc
              coarserFR->incrementFine(fluxmc[idir],m_dt,dit(),
                                       phiPatch.interval(),
                                       phiPatch.interval(),idir,Side::Lo);
              coarserFR->incrementFine(fluxpc[idir],m_dt,dit(),
                                       phiPatch.interval(),
                                       phiPatch.interval(),idir,Side::Hi);

              coarserFR->incrementFine(qadd[idir],m_dt,dit(),
                                       phiPatch.interval(),
                                       phiPatch.interval(),idir,Side::Lo);
              coarserFR->incrementFine(qadd[idir],-m_dt,dit(),
                                       phiPatch.interval(),
                                       phiPatch.interval(),idir,Side::Hi);
            }
        } // End of loop over directions for flux register updates.
    } // End of loop over patches in m_phiNew.

  // End application-dependent code - PC.


  // Increment the time.
  if (m_dt < 0.0)
    {
      MayDay::Error("negative time");
    }
  m_time += m_dt;

  // Update the time and store the new timestep
  // Find the minimum of dt's over this level
  // Real dtNew = m_dx / maxWaveSpeed;


  if (s_verbosity >= 1 && m_level == 0)
    {
      pout() << "Courant number = "
             << setw(8)
             << setprecision(3)
             << setiosflags(ios::fixed)
             << setiosflags(ios::showpoint)
             << maxcfl
             << "  dt = "
             << setw(10)
             << setprecision(4)
             << setiosflags(ios::showpoint)
             << setiosflags(ios::scientific)
             << m_dt
             << endl
             << endl;
    }

  // Real dtNew = m_dt*m_cfl/maxcfl;
  Real dtNew = m_dt/maxcfl;
  if (m_level == 0)
    {
      m_coarse_cfl = Max(m_coarse_cfl,maxcfl);
    }

  // Broadcast the minimum to processors and gather dtNew
  Vector<Real> allDt;
  gather(allDt,dtNew,uniqueProc(SerialTask::compute));
  if (procID() == uniqueProc(SerialTask::compute))
    {
      dtNew = allDt[0];
      for (int i = 1; i < allDt.size(); ++i)
        {
          dtNew = Min(dtNew,allDt[i]);
        }
    }
  broadcast(dtNew,uniqueProc(SerialTask::compute));

  m_dtNew = dtNew*m_cfl;
  if (maxcfl > m_clawPatch.get_max_cfl())
    {
      pout() << endl;
      pout() << "Maximum cfl number exceeded; cfl = " << maxcfl << endl;
      pout() << endl;
      m_dtNew = dtNew/2.0;  // take a smaller step next time
    }



  // Return the maximum stable time step

  return m_dtNew;

}
// Things to do after a timestep
void AMRLevelClaw::postTimeStep()
{
  // Used for conservation tests
  static Real orig_integral = 0.0;
  static Real last_integral = 0.0;
  static bool first = true;

  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelClaw::postTimeStep " << m_level << endl;
    }

  if (m_hasFiner)
    {
      // Reflux
      Real scale = -1.0/m_dx;
      m_fluxRegister.reflux(m_phiNew,scale);

      // Average from finer level data
      AMRLevelClaw* amrClawFinerPtr = getFinerLevel();
      if (m_hasKappa)
        {
          amrClawFinerPtr->m_coarseAverage.averageToCoarse(m_phiNew,
                                                           amrClawFinerPtr->m_phiNew,
                                                           amrClawFinerPtr->m_kappa,
                                                           amrClawFinerPtr->m_coarsened_kappa);
        }
      else
        {
          amrClawFinerPtr->m_coarseAverage.averageToCoarse(m_phiNew, amrClawFinerPtr->m_phiNew);
        }
    }

  if (s_verbosity >= 2 && m_level == 0)
    {
      int nRefFine = 1;

      pout() << "AMRLevelClaw::postTimeStep:" << endl;
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
          }

          pout() << endl;

          if (comp == 0)
            {
              if (first)
                {
                  orig_integral = integral;
                  first = false;
                }

              last_integral = integral;
            }
        }
    }
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelClaw::postTimeStep " << m_level << " finished" << endl;
    }
}

// Create tags for regridding
void AMRLevelClaw::tagCells(IntVectSet& a_tags)
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelClaw::tagCells " << m_level << endl;
    }

  // Since tags are calculated using only current time step data, use
  // the same tagging function for initialization and for regridding.
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelClaw::tagCellsInit " << m_level << endl;
    }

  // Begin application-dependent code - PC.

  // Create tags based on undivided gradient of density
  IntVectSet localTags;
  const DisjointBoxLayout& levelDomain = m_phiNew.disjointBoxLayout();
  // If there is a coarser level interpolate undefined ghost cells
  if (m_hasCoarser)
    {
      const AMRLevelClaw* amrClawCoarserPtr = getCoarserLevel();

      PiecewiseLinearFillPatch pwl(levelDomain,
                                   amrClawCoarserPtr->m_phiNew.disjointBoxLayout(),
                                   m_numStates,
                                   amrClawCoarserPtr->m_problem_domain,
                                   amrClawCoarserPtr->m_ref_ratio,
                                   1);

      pwl.fillInterp(m_phiNew,
                     amrClawCoarserPtr->m_phiNew,
                     amrClawCoarserPtr->m_phiNew,
                     1.0,
                     0,
                     0,
                     m_numStates);
    }
  m_phiNew.exchange(Interval(0,m_numStates-1));

  // Compute undivided gradient
  DataIterator dit = levelDomain.dataIterator();   // leveldomain is a disjointBoxLayout
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& b = levelDomain[dit()]; // does this include ghost cells?
      FArrayBox gradFab(b,SpaceDim);
      const FArrayBox& phiFab = m_phiNew[dit()];
      FArrayBox error_estimate(b,1);

      // isBoundary stores a 1 if patch is adjacent to physical boundary and not periodic,
      // and 0 otherwise.
      int isBoundary[2*SpaceDim];
      for (int idir = 0; idir < SpaceDim; ++idir)
        {
          const Box bCenter = b & grow(m_problem_domain,-BASISV(idir));  // shrinks problem domain
          const Box bLo     = b & adjCellLo(bCenter,idir);
          isBoundary[2*idir] = !(bLo.isEmpty() | m_problem_domain.isPeriodic(idir));
          const Box bHi     = b & adjCellHi(bCenter,idir);
          isBoundary[2*idir + 1] = ! (bHi.isEmpty() | m_problem_domain.isPeriodic(idir));
        }

      m_clawPatch.estimateError(phiFab, b, m_problem_domain,m_time,
                                m_dt, m_dx, m_level, isBoundary, m_refineThresh,
                                error_estimate);


      // Left over code from earlier version.  This is now done in call to
      // m_clawPatch.estimateError(...), above.  (DAC, 5/18/04)
      //        FORT_GETGRAD(CHF_FRA1(gradFab,dir),
      //                        CHF_CONST_FRA1(phiFab,0),
      //                        CHF_CONST_INT(dir),
      //                        CHF_BOX(bLo),
      //                        CHF_CONST_INT(hasLo),
      //                        CHF_BOX(bHi),
      //                        CHF_CONST_INT(hasHi),
      //                        CHF_BOX(bCenter));
      //         }
      //
      //       FArrayBox gradMagFab(b,1);
      //       FORT_MAGNITUDE(CHF_FRA1(gradMagFab,0),
      //                      CHF_CONST_FRA(gradFab),
      //                      CHF_BOX(b));


      // Tag where gradient exceeds threshold
      BoxIterator bit(b);
      for (bit.begin(); bit.ok(); ++bit)
        {
          const IntVect& iv = bit();

          if (error_estimate(iv) >= m_refineThresh)
            {
              localTags |= iv;
            }
        }
    }

  localTags.grow(m_tagBufferSize);

  // Need to do this in two steps unless a IntVectSet::operator &=
  // (ProblemDomain) operator is defined
  Box localTagsBox = localTags.minBox();
  localTagsBox &= m_problem_domain;
  localTags &= localTagsBox;

  a_tags = localTags;
}

// Create tags at initialization
void AMRLevelClaw::tagCellsInit(IntVectSet& a_tags)
{
  tagCells(a_tags);
}

// Set up data on this level after regridding
void AMRLevelClaw::regrid(const Vector<Box>& a_newGrids)
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelClaw::regrid " << m_level << endl;
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
  LevelData<FArrayBox> phiOld;
  phiOld.define(m_phiNew);
  m_phiNew.copyTo(m_phiNew.interval(),
                  phiOld,
                  phiOld.interval());

  // Reshape state with new grids
  IntVect ivGhost = m_numGhost*IntVect::Unit;
  m_phiNew.define(m_grids,m_numStates,ivGhost);
  m_phiOld.define(m_grids,m_numStates,ivGhost);

  if (m_maux > 0)
    {
      m_aux.define(m_grids,m_maux,ivGhost);
    }
  else
    {
      // This is just because I am a bit paranoid about sending around empty
      // pointers in a situation where the Clawpack Fortran routines require
      // 'aux' as an input argument.
      m_aux.define(m_grids,1,ivGhost);
    }

  if (m_hasKappa)
    {
      m_kappa.define(m_grids,1,IntVect::Zero);
    }


  // Set initial values for old and new grids,
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      // Probably not necessary to set the values here
      // m_phiNew[dit()].setVal(0.);
      // m_phiOld[dit()].setVal(0.);

      if (m_maux > 0)
        {
          const Box& box = m_grids.get(dit());
          m_clawPatch.setAuxArray(m_aux[dit()],box,m_dx,m_level);
        }
    }

  if (m_hasKappa)
    {
      Interval mcapa_int(m_mcapa-1,m_mcapa-1);
      Interval unit_int(0,0);
      m_aux.copyTo(mcapa_int,m_kappa,unit_int);
    }
  // Set up data structures
  levelSetup();   // this is the regridding call to levelSetup

  // Initialize this series of grids by interpolating from coarser level

  if (m_hasCoarser)
    {
      AMRLevelClaw* amrClawCoarserPtr = getCoarserLevel();

      if (m_hasKappa)
        {
          Interval unit_int(0,0);
          amrClawCoarserPtr->m_kappa.copyTo(unit_int,m_coarsened_kappa,unit_int);
          m_fineInterp.interpToFine(m_phiNew,amrClawCoarserPtr->m_phiNew,
                                    m_kappa, m_coarsened_kappa);
        }
      else
        {
          m_fineInterp.interpToFine(m_phiNew,amrClawCoarserPtr->m_phiNew);
        }

      //       // Copy from old state at this level in case there was already data at a finer level.
      phiOld.copyTo(phiOld.interval(),
                    m_phiNew,
                    m_phiNew.interval());
    }

}

// Initialize grids
void AMRLevelClaw::initialGrid(const Vector<Box>& a_newGrids)
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelClaw::initialGrid " << m_level << endl;
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
  if (m_maux > 0)
    {
      m_aux.define(m_grids,m_maux,ivGhost);
      if (m_hasKappa)
        {
          m_kappa.define(m_grids,1,IntVect::Zero);
        }
    }
  else
    {
      // Have to define it with at least one entry
      m_aux.define(m_grids,1,ivGhost);
    }

  // Set up data structures
  levelSetup();  // This is the initialGrid call to levelSetup
}

// Initialize data
void AMRLevelClaw::initialData()
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelClaw::initialData " << m_level << endl;
    }

  if (m_maux > 0)
    {
      DataIterator dit = m_grids.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          const Box& box = m_aux.getBoxes()[dit()];
          m_clawPatch.setAuxArray(m_aux[dit()],box,m_dx,m_level);
        }
      if (m_hasKappa)
        {
          Interval mcapa_int(m_mcapa-1,m_mcapa-1);
          Interval unit_int(0,0);
          m_aux.copyTo(mcapa_int,m_kappa,unit_int);

          // Have to do this here for first coarseAverage (which happens before regrid is called).
          if (m_hasCoarser)
            {
              Interval unit_int(0,0);
              getCoarserLevel()->m_kappa.copyTo(unit_int,m_coarsened_kappa,unit_int);
            }
        }
    }

  // Iterator of all grids in this level and initialize them
  DataIterator dit = m_grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& box = m_grids.get(dit());
      m_clawPatch.initialize(m_phiNew[dit()],m_aux[dit()],box,m_dx);
    }
}

// Things to do after initialization
void AMRLevelClaw::postInitialize()
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelClaw::postInitialize " << m_level << endl;
    }

}

void AMRLevelClaw::postRegrid(int a_base_level)
{
  if (s_verbosity >= 3)
    {
      pout() << "In postRegrid() at level " << m_level << endl;
    }
}

#ifdef CH_USE_HDF5

// Write checkpoint header
void AMRLevelClaw::writeCheckpointHeader(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelClaw::writeCheckpointHeader" << endl;
    }

  // Setup the number of components
  HDF5HeaderData header;
  header.m_int["num_components"] = m_numStates;

  // Setup the component names
  char compStr[30];
  for (int comp = 0; comp < m_numStates; ++comp)
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
void AMRLevelClaw::writeCheckpointLevel(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelClaw::writeCheckpointLevel" << endl;
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
  header.m_int["auxcomp"]          = m_aux.nComp();
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
  write(a_handle,m_phiNew.boxLayout());
  write(a_handle,m_phiNew,"data");
  write(a_handle,m_aux,"auxdata");
}

// Read checkpoint header
void AMRLevelClaw::readCheckpointHeader(HDF5Handle& a_handle)
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelClaw::readCheckpointHeader" << endl;
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
      MayDay::Error("AMRLevelClaw::readCheckpointHeader: checkpoint file does not have num_components");
    }

  int numStates = header.m_int["num_components"];
  if (numStates != m_numStates)
    {
      MayDay::Error("AMRLevelClaw::readCheckpointHeader: num_components in checkpoint file does not match solver");
    }

  // Get the component names
  std::string stateName;
  char compStr[60];
  for (int comp = 0; comp < m_numStates; ++comp)
    {
      sprintf(compStr,"component_%d",comp);
      if (header.m_string.find(compStr) == header.m_string.end())
        {
          MayDay::Error("AMRLevelClaw::readCheckpointHeader: checkpoint file does not have enough component names");
        }

      stateName = header.m_string[compStr];
      if (stateName != m_stateNames[comp])
        {
          MayDay::Error("AMRLevelClaw::readCheckpointHeader: state_name in checkpoint does not match solver");
        }
    }
}

// Read checkpoint data for this level
void AMRLevelClaw::readCheckpointLevel(HDF5Handle& a_handle)
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelClaw::readCheckpointLevel" << endl;
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
      MayDay::Error("AMRLevelClaw::readCheckpointLevel: file does not contain ref_ratio");
    }
  m_ref_ratio = header.m_int["ref_ratio"];

  if (s_verbosity >= 2)
    {
      pout() << "read ref_ratio = " << m_ref_ratio << endl;
    }

  // Get the tag buffer size
  if (header.m_int.find("tag_buffer_size") == header.m_int.end())
    {
      MayDay::Error("AMRLevelClaw::readCheckpointLevel: file does not contain tag_buffer_size");
    }
  m_tagBufferSize=  header.m_int["tag_buffer_size"];

  if (s_verbosity >= 2)
    {
      pout() << "read tag_buffer_size = " << m_tagBufferSize << endl;
    }

  // Get dx
  if (header.m_real.find("dx") == header.m_real.end())
    {
      MayDay::Error("AMRLevelClaw::readCheckpointLevel: file does not contain dx");
    }
  m_dx = header.m_real["dx"];

  if (s_verbosity >= 2)
    {
      pout() << "read dx = " << m_dx << endl;
    }

  // Get dt
  if (header.m_real.find("dt") == header.m_real.end())
    {
      MayDay::Error("AMRLevelClaw::readCheckpointLevel: file does not contain dt");
    }
  m_dt = header.m_real["dt"];

  if (s_verbosity >= 2)
    {
      pout() << "read dt = " << m_dt << endl;
    }

  // Get time
  if (header.m_real.find("time") == header.m_real.end())
    {
      MayDay::Error("AMRLevelClaw::readCheckpointLevel: file does not contain time");
    }
  m_time = header.m_real["time"];

  if (s_verbosity >= 2)
    {
      pout() << "read time = " << m_time << endl;
    }

  // Get the problem domain
  if (header.m_box.find("prob_domain") == header.m_box.end())
    {
      MayDay::Error("AMRLevelClaw::readCheckpointLevel: file does not contain prob_domain");
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
      MayDay::Error("AMRLevelClaw::readCheckpointLevel: file does not contain a Vector<Box>");
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
      LayoutIterator lit = m_grids.layoutIterator();
      for (lit.begin(); lit.ok(); ++lit)
        {
          const Box& b = m_grids[lit()];
          pout() << lit().intCode() << ": " << b << endl;
        }
      pout() << endl;
    }

  // Reshape state with new grids
  m_phiNew.define(m_grids,m_numStates);
  // need the number of aux components
  int numaux = header.m_int["auxcomp"];

  m_aux.define(m_grids,numaux);
  const int dataStatus = read<FArrayBox>(a_handle,
                                         m_phiNew,
                                         "data",
                                         m_grids);

  const int auxDataStatus = read<FArrayBox>(a_handle,
                                            m_aux,
                                            "auxdata",
                                            m_grids);
  if ((dataStatus != 0) || (auxDataStatus != 0))
    {
      MayDay::Error("AMRLevelClaw::readCheckpointLevel: file does not contain state data");
    }
  m_phiOld.define(m_grids,m_numStates);

  // Set up data structures
  levelSetup();  // This is the readCheckPoint call to levelSetup
}

// Write plotfile header
void AMRLevelClaw::writePlotHeader(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelClaw::writePlotHeader" << endl;
    }


  // Setup the number of components
  HDF5HeaderData header;
  char compStr[30];
  // Begin application-dependent code - PC.

  header.m_int["num_components"] = m_numStates;

  // Setup the component names
  //

  for (int comp = 0; comp < m_numStates; ++comp)
    {
      sprintf(compStr,"component_%d",comp);
      header.m_string[compStr] = m_stateNames[comp];
    }
  // End application-dependent code - PC.

  // Write the header

  header.writeToFile(a_handle);

  if (s_verbosity >= 3)
    {
      pout() << header << endl;
    }
}

// Write plotfile data for this level
void AMRLevelClaw::writePlotLevel(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelClaw::writePlotLevel" << endl;
    }

  /*
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
  */

  writeLevel(a_handle,m_level,m_phiNew,m_dx,m_dtNew,m_time,
             m_problem_domain.domainBox(),m_ref_ratio);


  // write(a_handle,m_phiNew.boxLayout());

  // pout() << "m_phiNew.boxLayout().size() = " << m_phiNew.boxLayout().size() << endl;
  // write(a_handle,m_phiNew,"data");
}

#endif

// Returns the dt computed earlier for this level
Real AMRLevelClaw::computeDt()
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelClaw::computeDt " << m_level << endl;
    }

  Real newDt;
  newDt = m_dtNew;

  return newDt;
}

// Compute dt using initial data
Real AMRLevelClaw::computeInitialDt()
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelClaw::computeInitialDt " << m_level << endl;
    }

  // Real newDT = m_initial_dt_multiplier * m_dx;

  // Added 5/23 :
  Real newDT = m_clawPatch.get_initial_dt();
  return newDT;
}

// Set the CFL number
void AMRLevelClaw::CFL(Real a_cfl)
{
  m_cfl = a_cfl;
}

// Set the physical dimension of the longest side of the domain
void AMRLevelClaw::domainLength(Real a_domainLength)
{
  m_domainLength = a_domainLength;
}

// Set the refinement threshold
void AMRLevelClaw::refinementThreshold(Real a_refineThresh)
{
  m_refineThresh = a_refineThresh;
}

// Set the tag buffer size
void AMRLevelClaw::tagBufferSize(int a_tagBufferSize)
{
  m_tagBufferSize = a_tagBufferSize;
}


// Create a load-balanced DisjointBoxLayout from a collection of Boxes
DisjointBoxLayout AMRLevelClaw::loadBalance(const Vector<Box>& a_grids)
{
  // Load balance and create boxlayout
  Vector<int> procMap;

  // appears to be faster for all procs to do the loadbalance (ndk)
  LoadBalance(procMap,a_grids);

  if (s_verbosity >= 4)
    {
      pout() << "AMRLevelClaw::loadBalance: procesor map: " << endl;
      for (int igrid = 0; igrid < a_grids.size(); ++igrid)
        {
          pout() << igrid << ": " << procMap[igrid] << "  " << endl;
        }
      pout() << endl;
    }

  DisjointBoxLayout dbl(a_grids,procMap,m_problem_domain);
  dbl.close();

  return dbl;
}

// Setup menagerie of data structures
void AMRLevelClaw::levelSetup()
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelClaw::levelSetup " << m_level << endl;
    }

  AMRLevelClaw* amrClawCoarserPtr = getCoarserLevel();
  AMRLevelClaw* amrClawFinerPtr   = getFinerLevel();

  m_hasCoarser = (amrClawCoarserPtr != NULL);
  m_hasFiner   = (amrClawFinerPtr   != NULL);

  if (m_hasCoarser)
    {

      int nRefCrse = amrClawCoarserPtr->refRatio();

      m_coarseAverage.define(m_grids,
                             m_numStates,
                             nRefCrse);
      m_fineInterp.define(m_grids,
                          m_numStates,
                          nRefCrse,
                          m_problem_domain);

      m_patcher.define(m_grids,
                       amrClawCoarserPtr->m_grids,
                       m_numStates,
                       amrClawCoarserPtr->problemDomain(),
                       nRefCrse,
                       m_numGhost);

      // Define m-coarsened_kappa, which will store kappa at coarser level.
      if (m_hasKappa)
        {
          // This is done here, since this routine "LevelSetup" gets called
          // from both InitialGrid, and from regrid.

          // Create coarse grids that outline fine grids.
          DisjointBoxLayout dblcoarse;
          coarsen(dblcoarse,m_grids,nRefCrse);

          // Create coarse kappa
          m_coarsened_kappa.define(dblcoarse,1,IntVect::Zero); // one comp; no ghost cells
        }

      // This may look twisted but you have to do this this way because the
      // coarser levels get setup before the finer levels so, since a flux
      // register lives between this level and the next FINER level, the finer
      // level has to do the setup because it is the only one with the
      // information at the time of construction.

      // Maintain flux registers

      amrClawCoarserPtr->m_fluxRegister.define(m_grids,
                                               amrClawCoarserPtr->m_grids,
                                               m_problem_domain,
                                               amrClawCoarserPtr->m_ref_ratio,
                                               m_numStates);
      amrClawCoarserPtr->m_fluxRegister.setToZero();

    }
}

// Get the next coarser level
AMRLevelClaw* AMRLevelClaw::getCoarserLevel() const
{
  // Begin application-dependent code - PC.

  AMRLevelClaw* coarserPtr = NULL;

  if (m_coarser_level_ptr != NULL)
    {
      coarserPtr = dynamic_cast<AMRLevelClaw*>(m_coarser_level_ptr);

      if (coarserPtr == NULL)
        {
          MayDay::Error("AMRLevelClaw::getCoarserLevel: dynamic cast failed");
        }
    }

  return coarserPtr;

  // End application-dependent code - PC.
}

// Get the next finer level
AMRLevelClaw* AMRLevelClaw::getFinerLevel() const
{
  AMRLevelClaw* amrClawFinerPtr = NULL;

  if (m_finer_level_ptr != NULL)
    {
      amrClawFinerPtr = dynamic_cast<AMRLevelClaw*>(m_finer_level_ptr);

      if (amrClawFinerPtr == NULL)
        {
          MayDay::Error("AMRLevelClaw::getFinerLevel: dynamic cast failed");
        }
    }

  return amrClawFinerPtr;
}
