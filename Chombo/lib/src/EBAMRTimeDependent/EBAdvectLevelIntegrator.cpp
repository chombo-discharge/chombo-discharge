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
#include "EBAdvectLevelIntegrator.H"
#include "BaseIVFactory.H"
#include "BaseIFFactory.H"
#include "BaseIFFAB.H"
#include "EBFluxFAB.H"
#include "FaceIterator.H"
#include "REAL.H"
#include "EBCellFactory.H"
#include "FaceIterator.H"
#include "EBPatchAdvectF_F.H"
#include <iomanip>
#include <cmath>
#include <cstdio>
#include <string>
#include "EBAMRIO.H"
#include "CH_OpenMP.H"
#include "NamespaceHeader.H"

void
EBAdvectLevelIntegrator::
resetBCs(const RefCountedPtr<EBPhysIBCFactory>&  a_advectBC)
{
  for (DataIterator dit = m_ebPatchAdvect.dataIterator(); dit.ok(); ++dit)
    {
      m_ebPatchAdvect[dit()]->setEBPhysIBC(*a_advectBC);
    }
}
/*****************************/
EBAdvectLevelIntegrator::
EBAdvectLevelIntegrator(const EBLevelGrid &                a_thisEBLG,
                        const EBLevelGrid &                a_coarEBLG,
                        const int&                         a_nRefine,
                        const RealVect&                    a_dx,
                        const bool&                        a_hasCoarser,
                        const bool&                        a_hasFiner,
                        const bool &                       a_forceNoEBCF,
                        const bool &                       a_useSlopeLimiting,
                        const EBIndexSpace*                const a_eb)
{
  CH_TIME("EBAdvectLevelIntegrator::EBAdvectLevelIntegrator");
  CH_assert(a_dx[0] > 0.0);
  CH_assert(a_nRefine > 0);
  m_useSlopeLimiting = a_useSlopeLimiting;
  m_forceNoEBCF  = a_forceNoEBCF;
  m_thisEBLG     = a_thisEBLG;
  m_refRatCrse   = a_nRefine;
  m_dx           = a_dx;
  m_hasCoarser   = a_hasCoarser;
  m_hasFiner     = a_hasFiner;
  if (m_hasCoarser)
    {
      m_coarEBLG = a_coarEBLG;
    }

  m_ebPatchAdvect.define(m_thisEBLG.getDBL());
  m_nVar         =  1;

  for (DataIterator dit = m_ebPatchAdvect.dataIterator(); dit.ok(); ++dit)
    {
      IntVectSet cfivs =  (*(m_thisEBLG.getCFIVS()))[dit()]; 
      m_ebPatchAdvect[dit()] = 
        RefCountedPtr<EBAdvectPatchIntegrator>
        (new          EBAdvectPatchIntegrator(m_thisEBLG.getDBL()[dit()], 
                                              m_thisEBLG.getEBISL()[dit()], 
                                              cfivs, 
                                              m_thisEBLG.getDomain(), m_dx,
                                              m_useSlopeLimiting));
    }
  m_nGhost = 4;

  if (m_hasCoarser)
    {
      CH_TIME("EBAdvectLevelIntegrator::define::fillPatchDefine");
      ProblemDomain domainCrse = coarsen(m_thisEBLG.getDomain(), m_refRatCrse);

      IntVect ivGhost = m_nGhost*IntVect::Unit;
      //patcher is defined with the number of conserved vars.
      m_fillPatch = RefCountedPtr<AggEBPWLFillPatch>
        (new                      AggEBPWLFillPatch(m_thisEBLG, m_coarEBLG,
                                                    m_refRatCrse, m_nVar,
                                                    m_nGhost, ivGhost, m_forceNoEBCF, a_eb));

      m_fillPatchVel= RefCountedPtr<AggEBPWLFillPatch>
        (new                      AggEBPWLFillPatch(m_thisEBLG, m_coarEBLG,
                                                    m_refRatCrse, SpaceDim,
                                                    m_nGhost, ivGhost, m_forceNoEBCF, a_eb));

    }

}
/*****************************/
void
EBAdvectLevelIntegrator::
advectToFacesCol(LevelData< EBFluxFAB >&                         a_extrapState,
                 LayoutData< Vector< BaseIVFAB<Real>* > >&       a_coveredPrimLo,
                 LayoutData< Vector< BaseIVFAB<Real>* > >&       a_coveredPrimHi,
                 const LayoutData< Vector< Vector<VolIndex> > >& a_coveredFaceLo,
                 const LayoutData< Vector< Vector<VolIndex> > >& a_coveredFaceHi,
                 const LayoutData< Vector< IntVectSet> >&        a_coveredSetsLo,
                 const LayoutData< Vector< IntVectSet> >&        a_coveredSetsHi,
                 const LevelData< EBCellFAB >&                   a_consState,
                 const LevelData< EBCellFAB >&                   a_normalVel,
                 const LevelData< EBFluxFAB >&                   a_advectionVel,
                 const LevelData< EBCellFAB >*                   a_consStateCoarseOld,
                 const LevelData< EBCellFAB >*                   a_consStateCoarseNew,
                 const LevelData< EBCellFAB >*                   a_normalVelCoarseOld,
                 const LevelData< EBCellFAB >*                   a_normalVelCoarseNew,
                 const Real&                                     a_timeCoarseOld,
                 const Real&                                     a_timeCoarseNew,
                 const Real&                                     a_timeFine,
                 const Real&                                     a_dt,
                 const LevelData<EBCellFAB>* const               a_source,
                 const LevelData<EBCellFAB>* const               a_sourceCoarOld,
                 const LevelData<EBCellFAB>* const               a_sourceCoarNew)
{
  CH_TIME("EBAdvectLevelIntegrator::advectToFacesCol(level)");

  //create temp data with the correct number of ghost cells
  IntVect ivGhost = m_nGhost*IntVect::Unit;
  Interval consInterv(0, m_nVar-1);
  Interval intervSD(0, SpaceDim-1);

  // LevelData<EBCellFAB>& consTemp = (LevelData<EBCellFAB>&) a_consState;
  // LevelData<EBCellFAB>& veloTemp = (LevelData<EBCellFAB>&) a_normalVel;

  EBCellFactory factory(m_thisEBLG.getEBISL());
  LevelData<EBCellFAB> consTemp(m_thisEBLG.getDBL(), m_nVar, ivGhost, factory);
  LevelData<EBCellFAB> veloTemp(m_thisEBLG.getDBL(), SpaceDim, ivGhost, factory);
  for (DataIterator dit = m_thisEBLG.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      consTemp[dit()].setVal(0.);
    }

  //this should be strictly local
  CH_assert(m_thisEBLG.getDBL() == a_consState.disjointBoxLayout());
  CH_assert(m_thisEBLG.getDBL() == a_normalVel.disjointBoxLayout());
  //if the above asserts are satisfied, copyto is strictly local
  {
    CH_TIME("copies");
    a_consState.localCopyTo(consInterv, consTemp, consInterv);
    a_normalVel.localCopyTo(intervSD  , veloTemp, intervSD);
  }
  // Fill ghost cells using fillInterp, and copyTo.
  if (m_hasCoarser)
    {
      CH_TIME("fillPatch");
      m_fillPatch->interpolate(consTemp,
                               *a_consStateCoarseOld,
                               *a_consStateCoarseNew,
                               a_timeCoarseOld,
                               a_timeCoarseNew,
                               a_timeFine,
                               consInterv);

      m_fillPatchVel->interpolate(veloTemp,
                                  *a_normalVelCoarseOld,
                                  *a_normalVelCoarseNew,
                                  a_timeCoarseOld,
                                  a_timeCoarseNew,
                                  a_timeFine,
                                  intervSD);
    }
  // Exchange all the data between grids
  {
    CH_TIME("initial_exchange");
    consTemp.exchange(consInterv);
    veloTemp.exchange(intervSD  );
  }

  LevelData<EBCellFAB>* srcTmpPtr = NULL;
  if (a_source != NULL)
    {
      // srcTmpPtr = (LevelData<EBCellFAB>*) a_source;
      CH_TIME("source_related_stuff");
      srcTmpPtr = new LevelData<EBCellFAB>(m_thisEBLG.getDBL(), m_nVar, ivGhost, factory);
      for (DataIterator dit = m_thisEBLG.getDBL().dataIterator(); dit.ok(); ++dit)
        {
          (*srcTmpPtr)[dit()].setVal(0.);
        }
      a_source->localCopyTo(consInterv, *srcTmpPtr, consInterv);

      if ( (a_sourceCoarOld != NULL) &&
           (a_sourceCoarNew != NULL) &&
           (m_hasCoarser) )
        {
          CH_TIME("fillPatch");

          m_fillPatch->interpolate(*srcTmpPtr,
                                   *a_sourceCoarOld,
                                   *a_sourceCoarNew,
                                   a_timeCoarseOld,
                                   a_timeCoarseNew,
                                   a_timeFine,
                                   consInterv);
          {
            CH_TIME("initial_exchange");
            srcTmpPtr->exchange(consInterv);
          }
        }
    }

  {
    CH_TIME("advectToFaces");
    DataIterator dit = m_thisEBLG.getDBL().dataIterator();
    int nbox=dit.size();
#pragma omp parallel for
    for (int mybox=0;mybox<nbox; mybox++)
      {
        const Box& cellBox = m_thisEBLG.getDBL().get(dit[mybox]);
        const EBISBox& ebisBox = m_thisEBLG.getEBISL()[dit[mybox]];
        if (!ebisBox.isAllCovered())
          {
            //the viscous term goes into here
            EBCellFAB dummy;
            EBCellFAB* srcPtr = &dummy;
            if (srcTmpPtr != NULL)
              {
                srcPtr = (EBCellFAB*)(&((*srcTmpPtr)[dit[mybox]]));
              }

            const EBCellFAB& source = *srcPtr;

            EBFluxFAB& extrapFAB  = a_extrapState[dit[mybox]];
            advectToFaces(extrapFAB,
                          a_coveredPrimLo[dit[mybox]],
                          a_coveredPrimHi[dit[mybox]],
                          a_coveredFaceLo[dit[mybox]],
                          a_coveredFaceHi[dit[mybox]],
                          a_coveredSetsLo[dit[mybox]],
                          a_coveredSetsHi[dit[mybox]],
                          consTemp[dit[mybox]],
                          veloTemp[dit[mybox]],
                          a_advectionVel[dit[mybox]],
                          cellBox, ebisBox,
                          a_dt,   a_timeFine,
                          source,  dit[mybox]);
          }
      }
  }
  if (srcTmpPtr != NULL)
    {
      delete srcTmpPtr;
    }
}
/*****************************/
void
EBAdvectLevelIntegrator::
advectToFaces(EBFluxFAB&                         a_extrapState,
              Vector< BaseIVFAB<Real>* >&        a_coveredPrimLo,
              Vector< BaseIVFAB<Real>* >&        a_coveredPrimHi,
              const Vector< Vector<VolIndex> >&  a_coveredFaceLo,
              const Vector< Vector<VolIndex> >&  a_coveredFaceHi,
              const Vector< IntVectSet  >&       a_coveredSetsLo,
              const Vector< IntVectSet  >&       a_coveredSetsHi,
              const EBCellFAB &                  a_consState,
              const EBCellFAB &                  a_normalVel,
              const EBFluxFAB &                  a_advectionVel,
              const Box&                         a_cellBox,
              const EBISBox&                     a_ebisBox,
              const Real&                        a_dt,
              const Real&                        a_time,
              const EBCellFAB &                  a_source,
              const DataIndex&                   a_dit)
{
  CH_TIME("EBAdvectLevelIntegrator::advectToFaces");
  IntVectSet cfivs; //not used here.  only used in flux interpolation
  m_ebPatchAdvect[a_dit]->setVelocities(a_normalVel, a_advectionVel);

  m_ebPatchAdvect[a_dit]->extrapolatePrim(a_extrapState,
                                          a_coveredPrimLo,
                                          a_coveredPrimHi,
                                          a_coveredSetsLo,
                                          a_coveredSetsHi,
                                          a_coveredFaceLo,
                                          a_coveredFaceHi,
                                          a_consState,
                                          a_source,
                                          a_dit,
                                          a_time, a_dt);

}
/*****************************/
void
EBAdvectLevelIntegrator::
advectToFacesBCG(LevelData< EBFluxFAB >&                         a_extrapState,
                 const LevelData< EBCellFAB >&                   a_consState,
                 const LevelData< EBCellFAB >&                   a_normalVel,
                 const LevelData< EBFluxFAB >&                   a_advectionVel,
                 const LevelData< EBCellFAB >*                   a_consStateCoarseOld,
                 const LevelData< EBCellFAB >*                   a_consStateCoarseNew,
                 const LevelData< EBCellFAB >*                   a_normalVelCoarseOld,
                 const LevelData< EBCellFAB >*                   a_normalVelCoarseNew,
                 const Real&                                     a_timeCoarseOld,
                 const Real&                                     a_timeCoarseNew,
                 const Real&                                     a_timeFine,
                 const Real&                                     a_dt,
                 const LevelData<EBCellFAB>* const               a_source,
                 const LevelData<EBCellFAB>* const               a_sourceCoarOld,
                 const LevelData<EBCellFAB>* const               a_sourceCoarNew)
{
  CH_TIME("EBAdvectLevelIntegrator::advectToFacesBCG(level)");

  //create temp data with the correct number of ghost cells
  IntVect ivGhost = m_nGhost*IntVect::Unit;
  Interval consInterv(0, m_nVar-1);
  Interval intervSD(0, SpaceDim-1);

  // LevelData<EBCellFAB>& consTemp = (LevelData<EBCellFAB>&) a_consState;
  // LevelData<EBCellFAB>& veloTemp = (LevelData<EBCellFAB>&) a_normalVel;

  EBCellFactory factory(m_thisEBLG.getEBISL());
  LevelData<EBCellFAB> consTemp(m_thisEBLG.getDBL(), m_nVar, ivGhost, factory);
  LevelData<EBCellFAB> veloTemp(m_thisEBLG.getDBL(), SpaceDim, ivGhost, factory);
  for (DataIterator dit = m_thisEBLG.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      consTemp[dit()].setVal(0.);
    }

  a_consState.localCopyTo(consInterv, consTemp, consInterv);
  a_normalVel.localCopyTo(intervSD, veloTemp, intervSD);
  // Fill ghost cells using fillInterp, and copyTo.
  if (m_hasCoarser)
    {
      CH_TIME("fillPatch");
      m_fillPatch->interpolate(consTemp,
                               *a_consStateCoarseOld,
                               *a_consStateCoarseNew,
                               a_timeCoarseOld,
                               a_timeCoarseNew,
                               a_timeFine,
                               consInterv);

      m_fillPatchVel->interpolate(veloTemp,
                                  *a_normalVelCoarseOld,
                                  *a_normalVelCoarseNew,
                                  a_timeCoarseOld,
                                  a_timeCoarseNew,
                                  a_timeFine,
                                  intervSD);
    }
  // Exchange all the data between grids
  {
    CH_TIME("initial_exchange");
    consTemp.exchange(consInterv);
    veloTemp.exchange(intervSD  );
  }

  LevelData<EBCellFAB>* srcTmpPtr = NULL;
  if (a_source != NULL)
    {
      // srcTmpPtr = (LevelData<EBCellFAB>*) a_source;

      srcTmpPtr = new LevelData<EBCellFAB>(m_thisEBLG.getDBL(), m_nVar, ivGhost, factory);
      for (DataIterator dit = m_thisEBLG.getDBL().dataIterator(); dit.ok(); ++dit)
        {
          (*srcTmpPtr)[dit()].setVal(0.);
        }
      a_source->copyTo(consInterv, *srcTmpPtr, consInterv);

      if ( (a_sourceCoarOld != NULL) &&
           (a_sourceCoarNew != NULL) &&
           (m_hasCoarser) )
        {
          CH_TIME("fillPatch");

          m_fillPatch->interpolate(*srcTmpPtr,
                                   *a_sourceCoarOld,
                                   *a_sourceCoarNew,
                                   a_timeCoarseOld,
                                   a_timeCoarseNew,
                                   a_timeFine,
                                   consInterv);
          {
            CH_TIME("initial_exchange");
            srcTmpPtr->exchange(consInterv);
          }
        }
    }

    DataIterator dit = m_thisEBLG.getDBL().dataIterator();
    int nbox=dit.size();
#pragma omp parallel for
    for (int mybox=0;mybox<nbox; mybox++)
      {
        const Box& cellBox = m_thisEBLG.getDBL().get(dit[mybox]);
        const EBISBox& ebisBox = m_thisEBLG.getEBISL()[dit[mybox]];
        if (!ebisBox.isAllCovered())
          {
            //the viscous term goes into here
            EBCellFAB dummy;
            EBCellFAB* srcPtr = &dummy;
            if (srcTmpPtr != NULL)
              {
                srcPtr = (EBCellFAB*)(&((*srcTmpPtr)[dit[mybox]]));
              }

            const EBCellFAB& source = *srcPtr;

            EBFluxFAB& extrapFAB  = a_extrapState[dit[mybox]];
            advectToFacesBCG(extrapFAB,
                             consTemp[dit[mybox]],
                             veloTemp[dit[mybox]],
                             a_advectionVel[dit[mybox]],
                             cellBox,
                             ebisBox,
                             a_dt,
                             a_timeFine,
                             source,
                             dit[mybox]);

          }
      }//dit

  if (srcTmpPtr != NULL)
    {
      delete srcTmpPtr;
    }
}
/*****************************/
void
EBAdvectLevelIntegrator::
advectToFacesBCG(EBFluxFAB&                         a_extrapState,
                 const EBCellFAB &                  a_consState,
                 const EBCellFAB &                  a_normalVel,
                 const EBFluxFAB &                  a_advectionVel,
                 const Box&                         a_cellBox,
                 const EBISBox&                     a_ebisBox,
                 const Real&                        a_dt,
                 const Real&                        a_time,
                 const EBCellFAB &                  a_source,
                 const DataIndex&                   a_dit)
{
  CH_TIME("EBAdvectLevelIntegrator::advectToFacesBCG(fluxfab)");
  //IntVectSet cfivs; //not used here.  only used in flux interpolation

  m_ebPatchAdvect[a_dit]->setVelocities(a_normalVel, a_advectionVel);

  m_ebPatchAdvect[a_dit]->extrapolateBCG(a_extrapState,
                                         a_consState,
                                         a_source,
                                         a_dit, a_time, a_dt);


}
/*****************************/
void
EBAdvectLevelIntegrator::
computeNormalVel(LevelData<EBCellFAB>&                          a_normalVel,
                 const LevelData<EBFluxFAB>&                    a_advectionVel,
                 const LayoutData<Vector<BaseIVFAB<Real> * > >& a_coveredVeloLo,
                 const LayoutData<Vector<BaseIVFAB<Real> * > >& a_coveredVeloHi,
                 const LayoutData<Vector<Vector<VolIndex> > >&  a_coveredFaceLo,
                 const LayoutData<Vector<Vector<VolIndex> > >&  a_coveredFaceHi) const
{
  CH_TIME("EBAdvectLevelIntegrator::computeNormalVel");
  for (DataIterator dit = m_thisEBLG.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      const Box& cellBox = m_thisEBLG.getDBL().get(dit());
      Box  grownBox = grow(cellBox, 1);
      grownBox &= m_thisEBLG.getDomain();
      m_ebPatchAdvect[dit]->averageVelToCC(a_normalVel[dit()],
                                           a_advectionVel[dit()],
                                           a_coveredVeloLo[dit()],
                                           a_coveredVeloHi[dit()],
                                           a_coveredFaceLo[dit()],
                                           a_coveredFaceHi[dit()],
                                           grownBox);
    }

}
/*****************************/
#include "NamespaceFooter.H"
