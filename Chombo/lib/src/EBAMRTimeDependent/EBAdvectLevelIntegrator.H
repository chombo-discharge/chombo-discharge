#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// dtg 12-03-2015

#ifndef _EBADVECTLEVELINTEGRATOR_H_
#define _EBADVECTLEVELINTEGRATOR_H_

#include "EBCellFAB.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "PiecewiseLinearFillPatch.H"
#include "AggEBPWLFillPatch.H"
#include "EBFluxRegister.H"
#include "ProblemDomain.H"
#include "EBLevelGrid.H"
#include "EBAdvectPatchIntegrator.H"
#include "NamespaceHeader.H"

///
/**
   Advection and diffusion on a level using
*/
class EBAdvectLevelIntegrator
{
public:
  /// Destructor
  /**
     destroys all objects created by define(..). Passed in
     data references of define(..) are left alone.
  */
  ~EBAdvectLevelIntegrator()
  {;}

  ///
  /**
     forceNoEBCF = true promises that the embedded boundary 
     will not cross the coarse-fine interface.   This will save you  
     a  lot of compute time if you are certain you are tagging all irregulkar cells.
  */
  EBAdvectLevelIntegrator(const EBLevelGrid      &           a_thisEBLG,
                          const EBLevelGrid      &           a_coarEBLG,
                          const int&                         a_nRefine,
                          const RealVect&                    a_dx,
                          const bool&                        a_hasCoarser,
                          const bool&                        a_hasFiner,
                          const bool&                        a_useSlopeLimiting,
                          const bool&                        a_forceNoEBCF = false,
                          const EBIndexSpace* const a_eb = Chombo_EBIS::instance());


  ///
  /**
     1.   Average advection velocity to cell centers to get unorm.
     2.   Do linear C/F interpolation on a_cellState   and compute
     the diffusion source term if necessary.   Quadratic is unnecessary
     because the laplacian does its own CF interpolation in the
     viscous case.
     3.   Do linear C/F interpolation in space and time
     4.   Extrapolate cellState to to faces using patchgodunov->extrapolatePrim

     Notes:  The advection velocity that is sent in is not used.  It will be
     taken out in later designs
  */
  void
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
                   const LevelData<EBCellFAB>* const               a_source = NULL,
                   const LevelData<EBCellFAB>* const               a_sourceCoarOld = NULL,
                   const LevelData<EBCellFAB>* const               a_sourceCoarNew = NULL);



  void
  advectToFaces(EBFluxFAB&                         a_extrapState,
                Vector< BaseIVFAB<Real>* >&        a_coveredPrimLo,
                Vector< BaseIVFAB<Real>* >&        a_coveredPrimHi,
                const Vector< Vector<VolIndex> >&  a_coveredFaceLo,
                const Vector< Vector<VolIndex> >&  a_coveredFaceHi,
                const Vector< IntVectSet >&        a_coveredSetsLo,
                const Vector< IntVectSet >&        a_coveredSetsHi,
                const EBCellFAB &                  a_consState,
                const EBCellFAB &                  a_normalVel,
                const EBFluxFAB &                  a_advectionVel,
                const Box&                         a_cellBox,
                const EBISBox&                     a_ebisBox,
                const Real&                        a_dt,
                const Real&                        a_time,
                const EBCellFAB &                  a_source,
                const DataIndex&                   a_dit);



  void
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
                   const LevelData<EBCellFAB>* const               a_source = NULL,
                   const LevelData<EBCellFAB>* const               a_sourceCoarOld = NULL,
                   const LevelData<EBCellFAB>* const               a_sourceCoarNew = NULL);



  void
  advectToFacesBCG(EBFluxFAB&                         a_extrapState,
                   const EBCellFAB &                  a_consState,
                   const EBCellFAB &                  a_normalVel,
                   const EBFluxFAB &                  a_advectionVel,
                   const Box&                         a_cellBox,
                   const EBISBox&                     a_ebisBox,
                   const Real&                        a_dt,
                   const Real&                        a_time,
                   const EBCellFAB &                  a_source,
                   const DataIndex&                   a_dit);



  void
  computeNormalVel(LevelData<EBCellFAB>&                          a_normalVel,
                   const LevelData<EBFluxFAB>&                    a_advectionVel,
                   const LayoutData<Vector<BaseIVFAB<Real> * > >& a_coveredVeloLo,
                   const LayoutData<Vector<BaseIVFAB<Real> * > >& a_coveredVeloHi,
                   const LayoutData<Vector<Vector<VolIndex> > >&  a_coveredFaceLo,
                   const LayoutData<Vector<Vector<VolIndex> > >&  a_coveredFaceHi) const;


  void resetBCs(const RefCountedPtr<EBPhysIBCFactory>&  a_advectBC);

  EBAdvectPatchIntegrator& getPatchAdvect(const DataIndex& a_dit)
  {
    return *(m_ebPatchAdvect[a_dit]);
  }

protected:
  
  bool                                                m_hasCoarser;
  bool                                                m_hasFiner;
  bool                                                m_useSlopeLimiting;
  bool                                                m_forceNoEBCF;
  RealVect                                            m_dx;
  int                                                 m_refRatCrse;
  int                                                 m_nVar;
  int                                                 m_nGhost;
                                                      
  EBLevelGrid                                         m_thisEBLG;
  EBLevelGrid                                         m_coarEBLG;
  RefCountedPtr<AggEBPWLFillPatch>                    m_fillPatch;
  RefCountedPtr<AggEBPWLFillPatch>                    m_fillPatchVel;
  ProblemDomain                                       m_domain;
  LayoutData<RefCountedPtr<EBAdvectPatchIntegrator> > m_ebPatchAdvect;

private:

  //disallowed for all the usual reasons
  void operator=(const EBAdvectLevelIntegrator& a_input)
  {
    MayDay::Error("invalid operator");
  }
  EBAdvectLevelIntegrator(const EBAdvectLevelIntegrator& a_input)
  {
    MayDay::Error("invalid operator");
  }
  /// weak construction is bad
  EBAdvectLevelIntegrator()
  {
  }

};

#include "NamespaceFooter.H"
#endif
