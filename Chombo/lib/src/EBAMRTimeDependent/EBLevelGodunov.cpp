#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "DebugOut.H"
#include "EBLevelGodunov.H"
#include "BaseIVFactory.H"
#include "BaseIFFactory.H"
#include "BaseIFFAB.H"
#include "EBFluxFAB.H"
#include "FaceIterator.H"
#include "REAL.H"
#include "EBCellFactory.H"
#include "FaceIterator.H"
#include <cstdio>
#include "EBAMRIO.H"
#include "EBArith.H"
#include "memusage.H"
#include "NamespaceHeader.H"
int EBLevelGodunov::s_timestep = 0;
/*****************************/
/*****************************/
EBLevelGodunov::EBLevelGodunov()
{
  m_isDefined = false;
  m_ebPatchGodunovSP = NULL;
}

/*****************************/
/*****************************/
bool
EBLevelGodunov::isDefined() const
{
  return m_isDefined;
}
/*****************************/
/*****************************/
EBLevelGodunov::~EBLevelGodunov()
{
  if(m_isDefined)
    {
      DataIterator dit = m_thisGrids.dataIterator();
      int nbox=dit.size();
      for (int ibox=0; ibox <nbox; ibox++)
        {
          delete(m_ebPatchGodunov[dit[ibox]]);
          m_ebPatchGodunov[dit[ibox]] = NULL;
        }
    }
  if(m_ebPatchGodunovSP != NULL)
    {
      delete(m_ebPatchGodunovSP);
      m_ebPatchGodunovSP = NULL;
    }

}

/*****************************/
/*****************************/
void
EBLevelGodunov::define(const DisjointBoxLayout&      a_thisDBL,
                       const DisjointBoxLayout&      a_coarDBL,
                       const EBISLayout&             a_thisEBISL,
                       const EBISLayout&             a_coarEBISL,
                       const ProblemDomain&          a_domain,
                       const int&                    a_nRefine,
                       const RealVect&               a_dx,
                       const bool&                   a_useMassRedist,
                       const bool&                   a_doSmushing,
                       const bool&                   a_doRZCoords,
                       const bool&                   a_hasSourceTerm,
                       const EBPatchGodunovFactory*  const a_patchGodunov,
                       const bool&                   a_hasCoarser,
                       const bool&                   a_hasFiner,
                       const IntVect &               a_ivGhost,
                       const bool&                   a_forceNoEBCF
                       )
{
  CH_TIME("EBLevelGodunov::define");
  CH_assert(a_dx[0] > 0.0);
  CH_assert(a_nRefine > 0);

  m_ivGhost = a_ivGhost;
  m_thisGrids = a_thisDBL;
  if(m_isDefined)
    {
      
      for (DataIterator dit = m_ebPatchGodunov.dataIterator(); dit.ok(); ++dit)
        {
          delete(m_ebPatchGodunov[dit()]);
          m_ebPatchGodunov[dit()] = NULL;
        }
    }

  m_isDefined = true;
  m_forceNoEBCF = a_forceNoEBCF;
  m_useMassRedist = a_useMassRedist;
  m_thisEBISL = a_thisEBISL;
  m_refRatCrse= a_nRefine;
  m_dx = a_dx;
  m_domain = a_domain;
  m_hasCoarser= a_hasCoarser;
  m_hasFiner= a_hasFiner;
  m_doSmushing = a_doSmushing;
  m_doRZCoords = a_doRZCoords;
  m_hasSourceTerm = a_hasSourceTerm;
  m_ebPatchGodunov.define(m_thisGrids);

  // Determing the number of ghost cells necessary here
  m_nGhost = 4;

  if(m_ebPatchGodunovSP != NULL)
    {
      delete(m_ebPatchGodunovSP);
      m_ebPatchGodunovSP = NULL;
    }
  m_ebPatchGodunovSP = a_patchGodunov->create();

  DataIterator dit = m_thisGrids.dataIterator();
  int nbox=dit.size();
  for (int ibox=0; ibox <nbox; ibox++)
    {
      EBPatchGodunov* ebPatchGodunov = a_patchGodunov->create();
      m_ebPatchGodunov[dit[ibox]] = dynamic_cast<EBPatchGodunov*>(ebPatchGodunov);
      if (m_ebPatchGodunov[dit[ibox]] == NULL)
        {
          MayDay::Error("problem in casting to patch advection class");
        }

      m_ebPatchGodunov[dit[ibox]]->define(m_domain, m_dx);
    }

  if (m_doRZCoords && !m_hasSourceTerm)
    {
      MayDay::Error("LG: RZ implies need of a source term--inconsistent inputs");
    }
  if (m_hasCoarser)
    {
      m_coarGrids = a_coarDBL;
      m_coarEBISL = a_coarEBISL;
    }

  m_nCons = m_ebPatchGodunovSP->numConserved();
  m_nFlux = m_ebPatchGodunovSP->numFluxes();

  //define redistribution object for this level
  //for now set to volume weighting
  //print_memory_line("before redist define");
  m_ebLevelRedist.define(m_thisGrids, m_thisEBISL,
                         m_domain, m_nCons);
  //print_memory_line("after redist define");

  if (m_hasCoarser)
    {
      CH_TIME("coarse_stuff_defs");
      ProblemDomain domainCrse = coarsen(m_domain, m_refRatCrse);

      //patcher is defined with the number of conserved vars.
      //print_memory_line("before patcher define");
      m_patcher.define(m_thisGrids, m_coarGrids,
                       m_thisEBISL, m_coarEBISL,
                       domainCrse, m_refRatCrse, m_nCons,
                       m_nGhost, m_ivGhost, m_forceNoEBCF,
                       Chombo_EBIS::instance());
      //print_memory_line("after patcher define");
    }

  //print_memory_line("before irregsets define");

  m_irregSetsSmall.define(m_thisGrids);

  for (int ibox=0; ibox <nbox; ibox++)
    {
      const EBISBox& ebisBox = m_thisEBISL[dit[ibox]];
      if (!ebisBox.isAllCovered())
        {
          const Box&     thisBox = m_thisGrids.get(dit[ibox]);
          m_irregSetsSmall[dit[ibox]] = ebisBox.getIrregIVS(thisBox);

        }
    }
  //print_memory_line("after irregsets define");

  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      CH_TIME("flux_interpolant_defs");
      //print_memory_line("before flux interpolant");
      EBArith::defineFluxInterpolant(m_fluxInterpolants[faceDir],
                                     m_irregSetsGrown  [faceDir],
                                     m_thisGrids, m_thisEBISL,
                                     m_domain, m_nFlux, faceDir);
      //print_memory_line("after flux_interpolant");
    }
  {
    CH_TIME("EBIrregFlux_defs");
    BaseIVFactory<Real> cellFactorySmall(m_thisEBISL, m_irregSetsSmall);

    //print_memory_line("before dnc");
    m_nonConsDivergence.define(m_thisGrids, m_nCons,
                               IntVect::Zero, cellFactorySmall);
    //print_memory_line("after dnc");
    m_ebIrregFaceFlux.define(m_thisGrids, m_nCons,
                             IntVect::Zero, cellFactorySmall);
    //print_memory_line("after ebirreg flux");
  }

  {
    CH_TIME("coarse-fine_ivs_defs");
    //ebpatch needs coarse-fine IVS to know where to drop order for interpolation
    m_cfIVS.define(m_thisGrids);
    for (DataIterator ditOuter = m_thisGrids.dataIterator();
        ditOuter.ok(); ++ditOuter)
      {
        const EBISBox& ebisBox = m_thisEBISL[ditOuter()];
        if (!ebisBox.isAllCovered())
          {
            // make a grown box of the grid and then subtract off each grid in the
            //domain to make the coarse-fine IVS
            // CFIVS = grow(b, 1) - sum(bi)
            Box grownBox = grow(m_thisGrids.get(ditOuter()), 1);
            grownBox &= m_domain;
            IntVectSet complementIVS(grownBox);

            for (LayoutIterator litInner = m_thisGrids.layoutIterator();
                litInner.ok(); ++litInner)
              {
                Box innerBox = m_thisGrids.get(litInner());
                complementIVS -= innerBox;
              }

            m_cfIVS[ditOuter()] = complementIVS;
          }
      }
  }

  for (int ibox=0; ibox <nbox; ibox++)
    {
      const IntVectSet& cfivs = m_cfIVS[dit[ibox]];
      //IntVectSet cfivs; //not used here.  only used in flux interpolation
      Real time = 0; //needs to get reset later
      Real dt   = 0; //needs to get reset later
      m_ebPatchGodunov[dit[ibox]]->setValidBox(m_thisGrids[dit[ibox]], m_thisEBISL[dit[ibox]], cfivs, time, dt);
    }

  {
    CH_TIME("flattening_defs");
    //create temp data with the correct number of ghost cells
    IntVect ivGhost = m_nGhost*IntVect::Unit;
    EBCellFactory factory(m_thisEBISL);
    IntVect flatGhostIV = 3*IntVect::Unit;
    m_flattening.define(m_thisGrids, 1, flatGhostIV, factory);
  }
}
/*****************************/
void
EBLevelGodunov::
floorConserved(LevelData<EBCellFAB>&         a_consState,
               Real a_time, Real a_dt)
{
  DataIterator dit = m_thisGrids.dataIterator();
  int nbox=dit.size();
  for (int ibox=0; ibox <nbox; ibox++)
    {
      EBCellFAB& consState = a_consState[dit[ibox]];
      //const IntVectSet& cfivs = m_cfIVS[dit[ibox]];
      //const EBISBox& ebisBox = m_thisEBISL[dit[ibox]];
      const Box& cellBox = m_thisGrids.get(dit[ibox]);
      //      m_ebPatchGodunov[dit[ibox]]->setValidBox(cellBox, ebisBox, cfivs, a_time, a_dt);
      m_ebPatchGodunov[dit[ibox]]->floorConserved(consState, cellBox);
    }
}
/*****************************/
void
EBLevelGodunov::
fillConsState(LevelData<EBCellFAB>&         a_consState,
              const LevelData<EBCellFAB>&   a_consStateCoarseOld,
              const LevelData<EBCellFAB>&   a_consStateCoarseNew,
              const Real&                   a_time,
              const Real&                   a_coarTimeOld,
              const Real&                   a_coarTimeNew)
{
  Interval consInterv(0, m_nCons-1);
  Interval fluxInterv(0, m_nFlux-1);

  if (m_hasCoarser)
    {
      //aggsten needs to know ghost up front
      CH_assert(a_consState.ghostVect() == m_ivGhost);
      m_patcher.interpolate(a_consState,
                            a_consStateCoarseOld,
                            a_consStateCoarseNew,
                            a_coarTimeOld,
                            a_coarTimeNew,
                            a_time,
                            consInterv);
    }

  // Exchange all the data between grids
  a_consState.exchange(consInterv);
}
/***************/
void
EBLevelGodunov::
computeFlattening(Real a_time, Real a_dt,
                  LevelData<EBCellFAB>&         a_consState)
{
  CH_TIME("eblevelgodunov::compute_flattening");
  //compute flattening coefficients. this saves a ghost cell. woo hoo.
  int ibox = 0;
  bool verbose = false;
  DataIterator dit = m_thisGrids.dataIterator();
  int nbox=dit.size();
  for (int jbox=0; jbox <nbox; jbox++)
    {
      EBCellFAB& consState = a_consState[dit[jbox]];
      const EBISBox& ebisBox = m_thisEBISL[dit[jbox]];
      if (!ebisBox.isAllCovered())
        {
          const Box& cellBox = m_thisGrids.get(dit[jbox]);
          //const IntVectSet& cfivs = m_cfIVS[dit[jbox]];
          //          m_ebPatchGodunov[dit[jbox]]->setValidBox(cellBox, ebisBox, cfivs, a_time, a_dt);
          m_ebPatchGodunov[dit[jbox]]->setTimeAndDt(a_time, a_dt);
          m_ebPatchGodunov[dit[jbox]]->setCoveredConsVals(consState);
          int nPrim  = m_ebPatchGodunov[dit[jbox]]->numPrimitives();
          EBCellFAB primState(ebisBox, consState.getRegion(), nPrim);
          int logflag = 0;
          //debug
          //verbose = true;
          //end debug
          m_ebPatchGodunov[dit[jbox]]->consToPrim(primState, consState, consState.getRegion(), logflag, verbose);
          EBCellFAB& flatteningFAB = m_flattening[dit[jbox]];
          //this will set the stuff over the coarse-fine interface
          flatteningFAB.setVal(1.);
          if (m_ebPatchGodunov[dit[jbox]]->usesFlattening())
            {
              m_ebPatchGodunov[dit[jbox]]->computeFlattening(flatteningFAB,
                                                  primState,
                                                  cellBox);
            }
          ibox++;
        }
    }
  Interval zerointerv(0,0);
  m_flattening.exchange(zerointerv);
}
/***************/
void
EBLevelGodunov::
doRegularUpdate(EBFluxRegister&               a_fineFluxRegister,
                EBFluxRegister&               a_coarFluxRegister,
                Real a_time, Real a_dt,
                LevelData<EBCellFAB>&         a_consState)
{
  bool verbose = false;
  //int ibox = 0;
  Interval consInterv(0, m_nCons-1);
  Interval fluxInterv(0, m_nFlux-1);
  DataIterator dit = m_thisGrids.dataIterator();
  int nbox=dit.size();
  for (int jbox=0; jbox <nbox; jbox++)
    {
      const Box& cellBox = m_thisGrids.get(dit[jbox]);
      const EBISBox& ebisBox = m_thisEBISL[dit[jbox]];
      if (!ebisBox.isAllCovered())
        {
          //const IntVectSet& cfivs = m_cfIVS[dit[jbox]];

          EBCellFAB& consState = a_consState[dit[jbox]];
          //          m_ebPatchGodunov[dit[jbox]]->setValidBox(cellBox, ebisBox, cfivs, a_time, a_dt);
          m_ebPatchGodunov[dit[jbox]]->setTimeAndDt(a_time, a_dt);

          EBCellFAB source;
          if (m_hasSourceTerm)
            {
              const Box& bigBox = consState.box();
              int nPrim = m_ebPatchGodunov[dit[jbox]]->numPrimitives();
              source.define(ebisBox, bigBox, nPrim);
              //this setval is important
              source.setVal(0.);
              m_ebPatchGodunov[dit[jbox]]->setSource(source, consState, bigBox);
            }

          EBFluxFAB flux(ebisBox, cellBox, m_nFlux);
          BaseIVFAB<Real>& nonConsDiv    = m_nonConsDivergence[dit[jbox]];
          BaseIVFAB<Real>& ebIrregFlux   = m_ebIrregFaceFlux[dit[jbox]];
          flux.setVal(7.89);
          ebIrregFlux.setVal(7.89);
          const IntVectSet& ivsIrreg     = m_irregSetsSmall[dit[jbox]];
          const EBCellFAB& flatteningFAB = m_flattening[dit[jbox]];

          m_ebPatchGodunov[dit[jbox]]->regularUpdate(consState, flux, ebIrregFlux,
                                          nonConsDiv,flatteningFAB,
                                          source, cellBox, ivsIrreg,
                                          dit[jbox],verbose);

          //do fluxregister cha-cha
          /*
            Coarse flux register is flux register with the coarse level.
            Fine flux register is the flux register with the fine level.
            To the finer level FR, this level is the coarse level.
            To the coarser level FR, this level is the fine level.
          */
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              Real scale = a_dt;

              EBFaceFAB fluxRegFlux;
              if (m_doRZCoords)
                {
                  //gather fluxes into flux-register compatible form
                  fluxRegFlux.define(ebisBox, cellBox, idir, m_nCons);
                  m_ebPatchGodunov[dit[jbox]]->assembleFluxReg(fluxRegFlux, flux[idir],
                                                    idir, cellBox);
                }
              if (m_hasFiner)
                {
                  a_fineFluxRegister.incrementCoarseRegular(flux[idir], scale,dit[jbox],
                                                            consInterv, idir);
                }

              if (m_hasCoarser)
                {
                  for (SideIterator sit; sit.ok(); ++sit)
                    {
                      a_coarFluxRegister.incrementFineRegular(flux[idir],scale, dit[jbox],
                                                              consInterv, idir,sit());
                    }
                }
            }

          //copy fluxes into sparse interpolant
          for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
            {
              IntVectSet ivsIrregGrown = m_irregSetsGrown[faceDir][dit[jbox]];
              ivsIrregGrown &= cellBox;
              FaceStop::WhichFaces stopCrit = FaceStop::SurroundingWithBoundary;

              BaseIFFAB<Real>& interpol = m_fluxInterpolants[faceDir][dit[jbox]];
              interpol.setVal(7.7777e7);
              EBFaceFAB& fluxDir = flux[faceDir];
              for (FaceIterator faceit(ivsIrregGrown, ebisBox.getEBGraph(),
                                      faceDir, stopCrit);
                  faceit.ok(); ++faceit)
                {
                  for (int ivar = 0; ivar < m_nFlux; ivar++)
                    {
                      interpol(faceit(), ivar) = fluxDir(faceit(), ivar);
                    }
                }

            }
        }
    }

  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      m_fluxInterpolants[faceDir].exchange(fluxInterv);
    }
}
/***************/
Real
EBLevelGodunov::
doIrregularUpdate(EBFluxRegister&               a_fineFluxRegister,
                  EBFluxRegister&               a_coarFluxRegister,
                  LevelData<BaseIVFAB<Real> >&  a_massDiff,
                  Real a_time, Real a_dt,
                  LevelData<EBCellFAB>&         a_consState)
{
  //now do the irregular update and the max wave speed
  Real maxWaveSpeed = 1.0e-12;
  int ibox = 0;
  Interval consInterv(0, m_nCons-1);
  Interval fluxInterv(0, m_nFlux-1);
  DataIterator dit = m_thisGrids.dataIterator();
  int nbox=dit.size();
  for (int jbox=0; jbox <nbox; jbox++, ibox++)
    //  for (DataIterator dit = m_thisGrids.dataIterator(); dit.ok(); ++dit, ibox++)
    {
      const Box& cellBox = m_thisGrids.get(dit[jbox]);
      const EBISBox& ebisBox = m_thisEBISL[dit[jbox]];
      if (!ebisBox.isAllCovered())
        {
          //const IntVectSet& cfivs = m_cfIVS[dit[jbox]];

          EBCellFAB& consState = a_consState[dit[jbox]];
          BaseIVFAB<Real>& redMass = a_massDiff[dit[jbox]];

          //          m_ebPatchGodunov[dit[jbox]]->setValidBox(cellBox, ebisBox, cfivs, a_time, a_dt);
          m_ebPatchGodunov[dit[jbox]]->setTimeAndDt(a_time, a_dt);

          BaseIFFAB<Real> centroidFlux[SpaceDim];
          const BaseIFFAB<Real>* interpolantGrid[SpaceDim];
          const IntVectSet& ivsIrregSmall = m_irregSetsSmall[dit[jbox]];
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              const BaseIFFAB<Real>& interpol = m_fluxInterpolants[idir][dit[jbox]];
              interpolantGrid[idir] = &interpol;
              BaseIFFAB<Real>& fluxDir= centroidFlux[idir];
              fluxDir.define(ivsIrregSmall, ebisBox.getEBGraph(), idir, m_nFlux);
            }
          m_ebPatchGodunov[dit[jbox]]->interpolateFluxToCentroids(centroidFlux,
                                                       interpolantGrid,
                                                       ivsIrregSmall);

          Real maxWaveSpeedGrid = 0.0;
          //update the state and interpolate the flux
          const BaseIVFAB<Real>& nonConsDiv = m_nonConsDivergence[dit[jbox]];
          const BaseIVFAB<Real>& ebIrregFlux = m_ebIrregFaceFlux[dit[jbox]];
          m_ebPatchGodunov[dit[jbox]]->irregularUpdate(consState,
                                            maxWaveSpeedGrid, redMass,
                                            centroidFlux, ebIrregFlux, nonConsDiv,
                                            cellBox, ivsIrregSmall);

          m_ebLevelRedist.increment(redMass, dit[jbox], consInterv);
          maxWaveSpeed = Max(maxWaveSpeed, maxWaveSpeedGrid);

          //do fluxregister mambo
          /*
            Coarse flux register is flux register with the coarse level.
            Fine flux register is the flux register with the fine level.
            To the finer level FR, this level is the coarse level.
            To the coarser level FR, this level is the fine level.
          */
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              Real scale = a_dt;

              BaseIFFAB<Real> fluxRegFlux;
              if (m_doRZCoords)
                {
                  //gather fluxes into flux-register compatible form
                  fluxRegFlux.define(ivsIrregSmall, ebisBox.getEBGraph(), idir, m_nCons);
                  m_ebPatchGodunov[dit[jbox]]->assembleFluxIrr(fluxRegFlux, centroidFlux[idir],
                                                    idir,  cellBox, ivsIrregSmall);
                }

              if (m_hasFiner)
                {
                  a_fineFluxRegister.incrementCoarseIrregular(centroidFlux[idir],
                                                              scale,dit[jbox],
                                                              consInterv, idir);
                }

              if (m_hasCoarser)
                {
                  for (SideIterator sit; sit.ok(); ++sit)
                    {
                      a_coarFluxRegister.incrementFineIrregular(centroidFlux[idir],
                                                                scale, dit[jbox],
                                                                consInterv, idir,sit());
                    }
                }
            }
        }
    }// end of loop over grids.
  a_consState.exchange(consInterv);
  return maxWaveSpeed;
}
/*****************************/
Real
EBLevelGodunov::
step(LevelData<EBCellFAB>&         a_consState,
     LevelData<BaseIVFAB<Real> >&  a_massDiff,
     EBFluxRegister&               a_fineFluxRegister,
     EBFluxRegister&               a_coarFluxRegister,
     const LevelData<EBCellFAB>&   a_consStateCoarseOld,
     const LevelData<EBCellFAB>&   a_consStateCoarseNew,
     const Real&                   a_time,
     const Real&                   a_coarTimeOld,
     const Real&                   a_coarTimeNew,
     const Real&                   a_dt)
{
  CH_TIME("levelgodunov_step");
  Interval consInterv(0, m_nCons-1);
  Interval fluxInterv(0, m_nFlux-1);
  CH_assert(isDefined());
  CH_assert(a_consState.disjointBoxLayout() == m_thisGrids);
  //pout() << "eblg here1" << endl;
  {
    CH_TIME("fillConsState");
    //cf interpolation and exchange
    fillConsState(a_consState, a_consStateCoarseOld, a_consStateCoarseNew, a_time, a_coarTimeOld, a_coarTimeNew);
  }

  //pout() << "eblg here2" << endl;
  {
    CH_TIME("early_coarse_fine");
    m_ebLevelRedist.setToZero();
    // clear flux registers with fine level
    //remember this level is the coarse level to the fine FR
    if (m_hasFiner)
      {
        a_fineFluxRegister.setToZero();
      }
  }

  //pout() << "eblg here3" << endl;
  {
    CH_TIME("computeFlattening");
    //compute flattening coefficients. this saves a ghost cell. woo hoo.
    computeFlattening(a_time, a_dt, a_consState);
  }

  //pout() << "eblg here3.5" << endl;
  {
    CH_TIME("doRegularUpdate");
    //this includes copying flux into flux interpolant and updating
    //regular grids and incrementing flux registers.
    doRegularUpdate(a_fineFluxRegister, a_coarFluxRegister, a_time, a_dt, a_consState);
  }

  //pout() << "eblg here4" << endl;
  Real maxWaveSpeed;
  {
    CH_TIME("doIrregularUpdate");
    //this does irregular update and deals with flux registers.
    //also computes the mass increment
    maxWaveSpeed = doIrregularUpdate(a_fineFluxRegister, a_coarFluxRegister, a_massDiff, a_time, a_dt, a_consState);
  }

  //pout() << "eblg here5" << endl;
  if (m_doSmushing)
    {
      CH_TIME("smushing");
      if (m_useMassRedist)
        {
          //if use mass weighting, need to
          //fix weights of redistribution object
          DataIterator dit=m_thisGrids.dataIterator();
          int densityIndex = m_ebPatchGodunovSP->densityIndex();
          a_consState.exchange(consInterv);
          m_ebLevelRedist.resetWeights(a_consState, densityIndex);
        }

      //do fluxregister samba
      //redistribute mass at this level.  Redistribution to other
      //levels (including reredist) is handled in EBAMRGodunov::postTimeStep
      m_ebLevelRedist.redistribute(a_consState, consInterv);
    }

  //pout() << "eblg here6" << endl;
  {
    CH_TIME("floors");
    floorConserved(a_consState, a_time, a_dt);
  }

  //max wave speed already gets broadcast
  maxWaveSpeed = getMaxWaveSpeed(a_consState);

  // Find the minimum of dt's over this level
  Real dtNew = m_dx[0] / maxWaveSpeed;

  //pout() << "eblg here7" << endl;
  //return the maximum stable time step
  return dtNew;
}
/*****************************/
void
EBLevelGodunov::getDrhoDtOverRho(LevelData<EBCellFAB>& a_drhoDt,
                                 const LevelData<EBCellFAB>& a_rhoNew,
                                 const LevelData<EBCellFAB>& a_rhoOld,
                                 const Real& a_dt)
{
  DataIterator dit = m_thisGrids.dataIterator();
  int nbox=dit.size();
  for (int jbox=0; jbox <nbox; jbox++)
    {
      EBCellFAB& drhodt = a_drhoDt[dit[jbox]];
      EBCellFAB& rhoold = (EBCellFAB&) a_rhoOld[dit[jbox]];
      EBCellFAB& rhonew = (EBCellFAB&) a_rhoNew[dit[jbox]];

      drhodt.setVal(0.);
      drhodt.setCoveredCellVal(1.0, 0);
      rhonew.setCoveredCellVal(1.0, 0);
      rhoold.setCoveredCellVal(1.0, 0);

      drhodt += rhoold;
      drhodt -= rhonew;
      drhodt /= rhonew;
      drhodt /= a_dt;
    }
}
/*****************************/
Real EBLevelGodunov::getMaxWaveSpeed(const LevelData<EBCellFAB>& a_state)
{
  CH_assert(a_state.disjointBoxLayout() == m_thisGrids);
  Real speed = 0.0;
  DataIterator dit = m_thisGrids.dataIterator();
  int nbox=dit.size();
  for (int jbox=0; jbox <nbox; jbox++)
    {
      const Box& validBox   = m_thisGrids.get(dit[jbox]);
      const EBISBox& ebisBox= m_thisEBISL[dit[jbox]];
      if (!ebisBox.isAllCovered())
        {
          //place holder not used maxWaveSpeed calc
          Real time = 0.0;
          //const IntVectSet& cfivs = m_cfIVS[dit[jbox]];

          //          m_ebPatchGodunov[dit[jbox]]->setValidBox(validBox, ebisBox, cfivs, time, time);
          m_ebPatchGodunov[dit[jbox]]->setTimeAndDt(time,time);
          Real speedOverBox = m_ebPatchGodunov[dit[jbox]]->getMaxWaveSpeed(a_state[dit[jbox]],
                                                                validBox);
          speed = Max(speed,speedOverBox);
        }
    }

  // gather speed
  Vector<Real> all_speeds;
  gather(all_speeds,speed,uniqueProc(SerialTask::compute));
  if (procID() == uniqueProc(SerialTask::compute))
    {
      speed = all_speeds[0];
      for (int i = 1; i < all_speeds.size (); ++i)
        {
          speed = Max(speed,all_speeds[i]);
        }
    }
  broadcast(speed,uniqueProc(SerialTask::compute));

  return speed;
}
/*****************************/
/*****************************/
#include "NamespaceFooter.H"
