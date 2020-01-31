#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBADVECTPATCHINTEGRATOR_H_
#define _EBADVECTPATCHINTEGRATOR_H_

#include "Box.H"
#include "IntVectSet.H"
#include "Vector.H"
#include "CH_HDF5.H"
#include "EBCellFAB.H"
#include "EBFluxFAB.H"
#include "EBISBox.H"
#include "BaseIVFAB.H"
#include "EBPhysIBC.H"
#include "EBPhysIBCFactory.H"
#include "Stencils.H"
#include "AggStencil.H"
#include "NamespaceHeader.H"

#define EBAPI_TOL 1.e-12

///
/**
   This does the same as EBPatchAdvect but presumably with less memory footprint. 
   and with a somewhat more sane interface.   

   I have pared down the EBPatchGodunov interface and tried to optimize for minimal 
   memory without completely destroying performance.
*/
class EBAdvectPatchIntegrator 
{
public:

  ///boundary condions are set via setEBPhysIBC
  EBAdvectPatchIntegrator(const Box&            a_validBox,
                          const EBISBox&        a_ebisBox,
                          const IntVectSet&     a_coarseFineIVS,
                          const ProblemDomain&  a_domain,
                          const RealVect&       a_dx,
                          bool                  a_useSlopeLimiting);
  
  virtual ~EBAdvectPatchIntegrator()
  {;}


  void
  extrapolatePrim(EBFluxFAB&                       a_flux,
                  Vector<BaseIVFAB<Real> * >&      a_coveredFluxMinu,
                  Vector<BaseIVFAB<Real> * >&      a_coveredFluxPlus,
                  const Vector<IntVectSet >&       a_coveredSetsMinu,
                  const Vector<IntVectSet >&       a_coveredSetsPlus,
                  const Vector<Vector<VolIndex> >& a_coveredFaceMinu,
                  const Vector<Vector<VolIndex> >& a_coveredFacePlus,
                  const EBCellFAB&                 a_consState,
                  const EBCellFAB&                 a_source,
                  const DataIndex&                 a_dit,
                  const Real&                      a_time,
                  const Real&                      a_dt);

  ///
  /**
   */
  void
  advectiveDerivative(EBCellFAB&                       a_uDotDelRho,
                      const EBFluxFAB&                 a_faceRho,
                      const EBFluxFAB&                 a_faceVel,
                      const Vector<BaseIVFAB<Real>*> & a_coveredRhoLo,
                      const Vector<BaseIVFAB<Real>*> & a_coveredRhoHi,
                      const Vector<BaseIVFAB<Real>*> & a_coveredVelLo,
                      const Vector<BaseIVFAB<Real>*> & a_coveredVelHi,
                      const Vector<Vector<VolIndex> >& a_coveredFaceLo,
                      const Vector<Vector<VolIndex> >& a_coveredFaceHi,
                      const Box&                       a_box);


  ///
  /**
     Version that leaves out the covered face stuff.  This is wrong near the EB but
     the codes that use it overwrite the EB stuff anyway.
   */
  void
  advectiveDerivative(EBCellFAB&                       a_uDotDelRho,
                      const EBFluxFAB&                 a_faceRho,
                      const EBFluxFAB&                 a_faceVel,
                      const Box&                       a_box);
  ///
  /**
   */
  void
  setEBPhysIBC(const EBPhysIBCFactory& a_bc)
  {
    m_bc = RefCountedPtr<EBPhysIBC>(a_bc.create());
    m_bc->define(m_domain, m_dx);
    m_isBCSet = true;
  }

  void setVelocities(const EBCellFAB& a_normalVel,
                     const EBFluxFAB& a_advectionVel)
  {
    m_normalVelPtr    = &a_normalVel;
    m_advectionVelPtr = &a_advectionVel;
    m_isVelSet = true;
  }
  

  /**
     For when EBFlux is always zero.
   */
  void
  consUndividedDivergence(BaseIVFAB<Real>&       a_divF,
                          const BaseIFFAB<Real>  a_centroidFlux[SpaceDim],
                          const IntVectSet&      a_ivs);


  void
  consUndividedDivergence(BaseIVFAB<Real>&       a_divF,
                          const BaseIFFAB<Real>  a_centroidFlux[SpaceDim],
                          const BaseIVFAB<Real>& a_ebIrregFlux,
                          const IntVectSet&      a_ivs);

  void
  interpolateFluxToCentroids(BaseIFFAB<Real>              a_centroidFlux[SpaceDim],
                             const BaseIFFAB<Real>* const a_fluxInterpolant[SpaceDim],
                             const IntVectSet&            a_irregIVS);

  void
  averageVelToCC(EBCellFAB&                        a_normalVel,
                 const EBFluxFAB&                  a_advectionVel,
                 const Vector<BaseIVFAB<Real> * >& a_coveredVeloLo,
                 const Vector<BaseIVFAB<Real> * >& a_coveredVeloHi,
                 const Vector<Vector<VolIndex> >&  a_coveredFaceLo,
                 const Vector<Vector<VolIndex> >&  a_coveredFaceHi,
                 const Box&                        a_box) const;


  void
  extrapolateBCG(EBFluxFAB&                       a_flux,
                 const EBCellFAB&                 a_consState,
                 const EBCellFAB&                 a_source,
                 const DataIndex&                 a_dit,
                 const Real&                      a_time,
                 const Real&                      a_dt);

  ///here are a couple of awful hooks necessary to get the minutiae of the algorithm correct
  static void setCurComp(int a_curComp)
  {
    s_curComp = a_curComp;
  }
  ///
  static void setDoingVel(int a_yesorno)
  {
    s_doingVel = a_yesorno;
  }

  static int getDoingVel()
  {
    return s_doingVel;
  }
  static int getCurComp()
  {
    return s_curComp;
  }

  ///
  /**
     This is called by EBAMRNoSubCycle.   The insane version that uses cell-centered
     data holders for face data (the plus-minus stuff) is called internally.
  */
  void
  extrapToCoveredFaces(BaseIVFAB<Real>&        a_extendedPrim,
                       const EBFaceFAB&        a_primFace,
                       const EBCellFAB&        a_primState,
                       const Vector<VolIndex>& a_coveredFaces,
                       const int&              a_faceDir,
                       const Side::LoHiSide&   a_sd,
                       const Box&              a_box);

  void setMaxMin(const Real& a_maxVal,
                 const Real& a_minVal)
  {
    m_isMaxMinSet = true;
    m_maxVal = a_maxVal;
    m_minVal = a_minVal;
  }
private:
  const EBFluxFAB* m_advectionVelPtr;
  const EBCellFAB* m_normalVelPtr;
  bool m_isVelSet;
  bool m_isMaxMinSet;
  bool m_useLimiting;
  Real m_maxVal;
  Real m_minVal;

  RefCountedPtr<EBPhysIBC> m_bc;
  bool             m_isBCSet;
  Box              m_validBox;
  ProblemDomain    m_domain;
  EBISBox          m_ebisBox;
  RealVect         m_dx;
  IntVectSet       m_cfivs;
  Box              m_validBoxG4;
  Vector<VolIndex> m_irregVoFs;

  IntVectSet       m_coveredSetsPlusG4[SpaceDim];
  IntVectSet       m_coveredSetsMinuG4[SpaceDim];
  Vector<VolIndex> m_coveredFacePlusG4[SpaceDim];
  Vector<VolIndex> m_coveredFaceMinuG4[SpaceDim];

  ///these exist because special things have to be done for velocity
  static int  s_curComp;
  static int  s_doingVel;

  /// weak construction is bad.
  EBAdvectPatchIntegrator();

  /// internal functions (all the madness below probably needs to get cleaned up)

  ///
  void extrapolatePrim3D(EBCellFAB           a_primMinu[SpaceDim],
                         EBCellFAB           a_primPlus[SpaceDim],
                         const EBCellFAB&    a_consState,
                         const EBCellFAB&    a_source,
                         const DataIndex&    a_dit,
                         const Real&         a_time,
                         const Real&         a_dt);

  ///
  void
  extrapolateBCG(EBCellFAB           a_primMinu[SpaceDim],
                 EBCellFAB           a_primPlus[SpaceDim],
                 const EBCellFAB&    a_consState,
                 const EBCellFAB&    a_source,
                 const DataIndex&    a_dit,
                 const Real&                      a_time,
                 const Real&                      a_dt);

  ///
  void extrapolatePrim2D(EBCellFAB           a_primMinu[SpaceDim],
                         EBCellFAB           a_primPlus[SpaceDim],
                         const EBCellFAB&    a_consState,
                         const EBCellFAB&    a_source,
                         const DataIndex&    a_dit,
                         const Real&         a_time,
                         const Real&         a_dt);

  ///
  void
  riemann(BaseIVFAB<Real>&        a_coveredPrim,
          const BaseIVFAB<Real>&  a_exteState,
          const EBCellFAB&        a_primState,
          const Vector<VolIndex>& a_vofset,
          const int&              a_faceDir,
          const Side::LoHiSide&   a_sd,
          const Box&       a_box);

  ///
  void
  riemann(EBFaceFAB&       a_primGdnv,
          const EBCellFAB& a_primLeft,
          const EBCellFAB& a_primRigh,
          const int&       a_faceDir,
          const Box&       a_box);

  ///
  FaceStencil getInterpStencil(const FaceIndex& a_face) const;


  void
  pointExtrapToCovered2D(Vector<Real>&           a_extrapVal,
                         const EBFaceFAB&        a_primFace,
                         const EBCellFAB&        a_primState,
                         const int&              a_faceDir,
                         const VolIndex&         a_vof,
                         const RealVect&         a_normal,
                         const Side::LoHiSide&   a_sd,
                         const int&              a_numPrim);

  void
  pointExtrapToCovered3D(Vector<Real>&           a_extrapVal,
                         const EBFaceFAB&        a_primFace,
                         const EBCellFAB&        a_primState,
                         const int&              a_faceDir,
                         const VolIndex&         a_vof,
                         const RealVect&         a_normal,
                         const Side::LoHiSide&   a_sd,
                         const int&              a_numPrim);

  /// 
  /**
     This insane version called internally.   The primMinu and primPlus stuff
     are really face centered data but because they *came* from cell-centered data
     they are left there.   
  */
  void
  extrapToCoveredFaces(BaseIVFAB<Real>&        a_extendedPrim,
                       const EBCellFAB&        a_primMinu,
                       const EBCellFAB&        a_primPlus,
                       const EBCellFAB&        a_primState,
                       const Vector<VolIndex>& a_coveredFaces,
                       const int&              a_faceDir,
                       const Side::LoHiSide&   a_sd,
                       const Box&              a_box);

  void
  pointExtrapToCovered2D(Vector<Real>&           a_extrapVal,
                         const EBCellFAB&        a_primMinu,
                         const EBCellFAB&        a_primPlus,
                         const EBCellFAB&        a_primState,
                         const int&              a_faceDir,
                         const VolIndex&         a_vof,
                         const RealVect&         a_normal,
                         const Side::LoHiSide&   a_sd,
                         const int&              a_numPrim);

  void
  pointExtrapToCovered3D(Vector<Real>&           a_extrapVal,
                         const EBCellFAB&        a_primMinu,
                         const EBCellFAB&        a_primPlus,
                         const EBCellFAB&        a_primState,
                         const int&              a_faceDir,
                         const VolIndex&         a_vof,
                         const RealVect&         a_normal,
                         const Side::LoHiSide&   a_sd,
                         const int&              a_numPrim);


  /// floors if m_isMaxMinSet
  void
  floorPrimitives(EBCellFAB&  a_primState,
                  const Box&  a_box);

  /// floors if m_isMaxMinSet
  void
  floorPrimitives(BaseIVFAB<Real>&   a_primState,
                  const IntVectSet&  a_set);

  void 
  coveredExtrapSlopes(Real&            a_dq,
                      const VolIndex&  a_vof,
                      const EBCellFAB& a_primState,
                      const int&       a_dir,
                      const int&       a_ivar);

  void
  pointGetSlopes(Real&            a_dql,
                 Real&            a_dqr,
                 Real&            a_dqc,
                 bool&            a_hasFacesLeft,
                 bool&            a_hasFacesRigh,
                 const VolIndex&  a_vof,
                 const EBCellFAB& a_primState,
                 const int&       a_dir,
                 const int&       a_ivar,
                 const bool&      a_verbose);

  Real
  bilinearFunc(const Real  a_WVal[2][2],
               const Real& a_xd1,
               const Real& a_xd2);

  void
  pointGetSlopesUpwind(Real&            a_dql,
                       Real&            a_dqr,
                       Real&            a_dqc,
                       bool&            a_hasFacesLeft,
                       bool&            a_hasFacesRigh,
                       const VolIndex&  a_vof,
                       const EBCellFAB& a_primState,
                       const int&       a_dir,
                       const int&       a_ivar,
                       const bool&      a_verbose);

  ///and this is the *simplified* version
  void
  doNormalDerivativeExtr2D(EBCellFAB              a_primMinu[SpaceDim],
                           EBCellFAB              a_primPlus[SpaceDim],
                           EBFaceFAB              a_fluxOne[SpaceDim],
                           BaseIVFAB<Real>        a_coveredFluxNormMinu[SpaceDim],
                           BaseIVFAB<Real>        a_coveredFluxNormPlus[SpaceDim],
                           Vector<VolIndex>       a_coveredFaceNormMinu[SpaceDim],
                           Vector<VolIndex>       a_coveredFaceNormPlus[SpaceDim],
                           const EBCellFAB&       a_primState,
                           const EBCellFAB&       a_source,
                           const DataIndex&       a_dit,
                           const Real     &       a_time,
                           const Real     &       a_dt);

  ///options for 4th ordeer slopes and flattening removed
  void
  slope(EBCellFAB&       a_slopePrim,
        const EBCellFAB& a_primState,
        const int&       a_dir,
        const Box&       a_box);


  void
  finalExtrap2D(EBCellFAB              a_primMinu[SpaceDim],
                EBCellFAB              a_primPlus[SpaceDim],
                const BaseIVFAB<Real>  a_coveredFluxNormMinu[SpaceDim],
                const BaseIVFAB<Real>  a_coveredFluxNormPlus[SpaceDim],
                const Vector<VolIndex> a_coveredFaceNormMinu[SpaceDim],
                const Vector<VolIndex> a_coveredFaceNormPlus[SpaceDim],
                const EBFaceFAB        a_fluxOne[SpaceDim],
                const EBCellFAB&       a_primState,
                const Real     &       a_time,
                const Real     &       a_dt);

  //and it goes on like this...
  void
  finalExtrap3D(EBCellFAB              a_primMinu[SpaceDim],
                EBCellFAB              a_primPlus[SpaceDim],
                const BaseIVFAB<Real>  a_coveredFluxMinu3D[SpaceDim][SpaceDim],
                const BaseIVFAB<Real>  a_coveredFluxPlus3D[SpaceDim][SpaceDim],
                const EBFaceFAB        a_fluxTwo[SpaceDim][SpaceDim],
                const EBCellFAB&       a_primState,
                const Real     &       a_time,
                const Real     &       a_dt);

  void
  do111coupling(EBFaceFAB              a_fluxTwo[SpaceDim][SpaceDim],
                BaseIVFAB<Real>        a_coveredFluxMinu3D[SpaceDim][SpaceDim],
                BaseIVFAB<Real>        a_coveredFluxPlus3D[SpaceDim][SpaceDim],
                const EBCellFAB        a_primMinu[SpaceDim],
                const EBCellFAB        a_primPlus[SpaceDim],
                const BaseIVFAB<Real>  a_coveredFluxNormMinu[SpaceDim],
                const BaseIVFAB<Real>  a_coveredFluxNormPlus[SpaceDim],
                const Vector<VolIndex> a_coveredFaceNormMinu[SpaceDim],
                const Vector<VolIndex> a_coveredFaceNormPlus[SpaceDim],
                const EBFaceFAB        a_fluxOne[SpaceDim],
                const EBCellFAB&       a_primState,
                const DataIndex&       a_dit,
                const Real     &       a_time,
                const Real     &       a_dt);

  void
  doNormalDerivativeExtr3D(EBCellFAB              a_primMinu[SpaceDim],
                           EBCellFAB              a_primPlus[SpaceDim],
                           EBFaceFAB              a_fluxOne[SpaceDim],
                           BaseIVFAB<Real>        a_coveredFluxNormMinu[SpaceDim],
                           BaseIVFAB<Real>        a_coveredFluxNormPlus[SpaceDim],
                           Vector<VolIndex>       a_coveredFaceNormMinu[SpaceDim],
                           Vector<VolIndex>       a_coveredFaceNormPlus[SpaceDim],
                           const EBCellFAB&       a_primState,
                           const EBCellFAB&       a_source,
                           const DataIndex&       a_dit,
                           const Real     &       a_time,
                           const Real     &       a_dt);
  

  void
  normalPred(EBCellFAB&       a_primLo,
             EBCellFAB&       a_primHi,
             const EBCellFAB& a_primState,
             const EBCellFAB& a_slopePrim,
             const Real&      a_scale,
             const int&       a_dir,
             const Box&       a_box);

  void
  incrementWithSource(EBCellFAB&       a_primState,
                      const EBCellFAB& a_source,
                      const Real&      a_scale,
                      const Box&       a_box) ;

  void
  updatePrim(EBCellFAB&              a_primMinu,
             EBCellFAB&              a_primPlus,
             const EBFaceFAB&        a_primFace,
             const BaseIVFAB<Real>&  a_coveredPrimMinu,
             const BaseIVFAB<Real>&  a_coveredPrimPlus,
             const Vector<VolIndex>& a_coveredFaceMinu,
             const Vector<VolIndex>& a_coveredFacePlus,
             const int&              a_faceDir,
             const Box&              a_box,
             const Real&             a_scale);

  void
  upwindSlope(EBCellFAB&       a_slopeUpWi,
              const EBCellFAB& a_primState,
              const int&       a_dir,
              const Box&       a_box);

  void
  transversePred(EBCellFAB&       a_rhoLo,
                 EBCellFAB&       a_rhoHi,
                 const EBCellFAB& a_rho,
                 const EBCellFAB& a_dRho,
                 const Real&      a_dtbydx,
                 const int&       a_dir,
                 const Box&       a_box);
};

#include "NamespaceFooter.H"
#endif
