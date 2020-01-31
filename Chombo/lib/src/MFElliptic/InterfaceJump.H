#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _INTERFACEJUMP_H_
#define _INTERFACEJUMP_H_

#include "LevelData.H"
#include "EBCellFAB.H"
#include "MFIndexSpace.H"
#include "Stencils.H"
#include "VoFIterator.H"
#include "BaseBCValue.H"
#include "NamespaceHeader.H"

class InterfaceJump
{
public:

  //constructor
  InterfaceJump();

  //destructor
  ~InterfaceJump();

  void define(const MFIndexSpace     & a_MFIS,
              const ProblemDomain    & a_domain,
              const DisjointBoxLayout& a_dbl,
              const RealVect         & a_vectDx,
              const RealVect         & a_physOrigin,
              const bool               a_isVectorJump = false);

  void setJump(const Real& a_gD,
               const Real& a_gN);

  void setJump(const RealVect& a_gD,
               const RealVect& a_gN);

  void setJump(const Vector< RefCountedPtr<BaseBCValue> >& a_phiValVect,
               const Vector< RefCountedPtr<BaseBCValue> >& a_flxValVect,
               const Vector<Real>&                         a_alpha,
               const Vector<Real>&                         a_beta,
               const Real&                                 a_time = 0.0);

  void resetJump(const Vector<Real>&  a_alpha,
                 const Vector<Real>&  a_beta,
                 const Real&          a_time);

  void getFaceStencil(VoFStencil&          a_istencil,
                      VoFStencil&          a_jstencil,
                      const Vector<Real>&  a_beta,
                      const VolIndex&      a_ivof,
                      const DataIndex&     a_dix,
                      const int&           a_iphase,
                      const int&           a_iface);

  void computeScalarBV(const EBCellFAB       & a_P0,
                       const EBCellFAB       & a_P1,
                       const Vector<Real>    & a_beta,
                       const DataIndex       & a_dataInd,
                       BaseIVFAB<Real>       & a_PBV0,
                       BaseIVFAB<Real>       & a_PBV1,
                       BaseIVFAB<Real>       & a_dpDN0,
                       BaseIVFAB<Real>       & a_dpDN1,
                       bool homogeneous = false);

  void computeVectorBV(const EBCellFAB       & a_P0,
                       const EBCellFAB       & a_P1,
                       const Vector<Real>    & a_beta,
                       const DataIndex       & a_dataInd,
                       BaseIVFAB<Real>       & a_PBV0,
                       BaseIVFAB<Real>       & a_PBV1,
                       BaseIVFAB<Real>       & a_dpDN0,
                       BaseIVFAB<Real>       & a_dpDN1,
                       bool homogeneous = false);

protected:

  //check is defined
  bool m_isDefined;
  // scalar or vector jump conditions
  bool m_isVectorJump;

  RealVect             m_vectDx;
  RealVect             m_d2Vect;
  Real                 m_dxScale;
  RealVect             m_origin;
  ProblemDomain        m_domain;
  int                  m_nComp;

  // jump values. defined on phase 0 VoFs
  LayoutData< BaseIVFAB< Vector<Real> > >       m_scalarGD, m_scalarGN;
  LayoutData< BaseIVFAB< Vector<RealVect> > >   m_vectorGD, m_vectorGN;

  Vector< RefCountedPtr<BaseBCValue> > m_phiValVect, m_flxValVect;

  //used for calculating phase 1 dpDn
  LayoutData< BaseIVFAB< Vector<VoFStencil> > > m_0stencils;
  LayoutData< BaseIVFAB< Vector<Real> > >       m_0inHomDirichWeight;
  LayoutData< BaseIVFAB< Vector<int> > >        m_0dirichletDropOrder;
  LayoutData< BaseIVFAB< Vector<VolIndex> > >   m_0boundaryVofs;
  LayoutData< BaseIVFAB< Vector<int> > >        m_0boundaryFaceIndices;

  //used for calculating phase 2 dpDn
  LayoutData< BaseIVFAB< Vector<VoFStencil> > > m_1stencils;
  LayoutData< BaseIVFAB< Vector<Real> > >       m_1inHomDirichWeight;
  LayoutData< BaseIVFAB< Vector<int> > >        m_1dirichletDropOrder;
  LayoutData< BaseIVFAB< Vector<VolIndex> > >   m_1boundaryVofs;
  LayoutData< BaseIVFAB< Vector<int> > >        m_1boundaryFaceIndices;

  //vofIterators
  LayoutData<VoFIterator > m_vofItBndry;

  //multiPhase ebisl
  EBISLayout m_0ebisl;
  EBISLayout m_1ebisl;

  //one box layout for all phases
  DisjointBoxLayout  m_dbl;

  void calculateJump(const Vector<Real>&  a_alpha,
                     const Vector<Real>&  a_beta,
                     const Real&          a_time);

  //defines the stencil data structures
  void defineLayoutData(EBISLayout&                                     a_ebisl,
                        LayoutData< BaseIVFAB< Vector<VoFStencil> > >&  a_stencils,
                        LayoutData< BaseIVFAB< Vector<Real> > >&        a_inHomDirichWeight,
                        LayoutData< BaseIVFAB< Vector<int> > >&         a_dirichletDropOrder,
                        LayoutData< BaseIVFAB< Vector<VolIndex> > >&    a_boundaryVofs,
                        LayoutData< BaseIVFAB< Vector<int> > >&         a_boundaryFaceIndices);

  // defines the jump for scalar case
  void defineScalarJumpLayoutData(EBISLayout&                               a_ebisl,
                                  LayoutData< BaseIVFAB< Vector<Real> > >&  a_scalarGD,
                                  LayoutData< BaseIVFAB< Vector<Real> > >&  a_scalarGN);

  // defines the jump for vector case
  void defineVectorJumpLayoutData(EBISLayout&                                   a_ebisl,
                                  LayoutData< BaseIVFAB< Vector<RealVect> > >&  a_vectorGD,
                                  LayoutData< BaseIVFAB< Vector<RealVect> > >&  a_vectorGN);

  //build vofIterators
  void cacheVoFIt();

  //builds the stencil per phase
  void definePhase(const MFIndexSpace&                             a_MFIS,
                   const ProblemDomain&                            a_domain,
                   const DisjointBoxLayout&                        a_dbl,
                   EBISLayout&                                     a_ebisl,
                   LayoutData< BaseIVFAB< Vector<VoFStencil> > >&  a_stencils,
                   LayoutData< BaseIVFAB< Vector<Real> > >&        a_inHomDirichWeight,
                   LayoutData< BaseIVFAB< Vector<int> > >&         a_dirichletDropOrder,
                   LayoutData< BaseIVFAB< Vector<VolIndex> > >&    a_boundaryVofs,
                   LayoutData< BaseIVFAB< Vector<int> > >&         a_boundaryFaceIndices,
                   const int&                                      a_numEBISLGhost,
                   const int&                                      a_phase);

  // computes the vofs lying across the boundary from a given vof
  void
  computeBoundaryLayoutData(LayoutData< BaseIVFAB< Vector<VolIndex> > >&  a_0boundaryVofs,
                            LayoutData< BaseIVFAB< Vector<VolIndex> > >&  a_1boundaryVofs,
                            LayoutData< BaseIVFAB< Vector<int> > >&       a_0boundaryFaceIndices,
                            LayoutData< BaseIVFAB< Vector<int> > >&       a_1boundaryFaceIndices);

  //Drop Order stencil if geometry too constricted
  void leastSquaresDirichStencil(VoFStencil&          a_stencil,
                                 Real&                a_weight,
                                 const RealVect&      a_normal,
                                 const RealVect&      a_centroid,
                                 const Real&          a_bndryArea,
                                 const VolIndex&      a_vof,
                                 const EBISBox&       a_ebisBox,
                                 const RealVect&      a_dx,
                                 const ProblemDomain& a_domain,
                                 int                  a_ivar);

  // Johansen stencil
  void johanDirichletStencil(const RealVect& a_normal,
                             const RealVect& a_centroid,
                             const VolIndex& a_vof,
                             const EBISBox&  a_ebisBox,
                             VoFStencil&     a_vofStencil,
                             int&            a_dirichletDropOrder,
                             Real&           a_inHomDirichWeight);

  void boundarySlope(Real*            a_dPhiDn,
                     const RealVect&  a_normal,
                     const VolIndex&  a_vof,
                     const EBCellFAB& a_phiFab,
                     const EBISBox&   a_ebisBox,
                     const int&       a_ivar = 0);

  void setDxConstants();

  void invert2x2(const Real& a_betaZeta1,
                 const Real& a_betaZeta2,
                 const Real& a_r1,
                 const Real& a_r2,
                 Real&       a_pB1,
                 Real&       a_pB2);

  void  normalDeriv(const Real& a_qs,
                    const Real& a_zetas,
                    const Real& a_pBs,
                    Real&       a_dPdns);

  //debug stuff
  void findTangentPlane(RealVect&      a_t1,
                        RealVect&      a_t2,
                        const RealVect a_normal);

  void crossProduct(      Real a_cp[3],
                    const Real a_c1[3],
                    const Real a_c2[3]);

  Real euclidianNorm(const RealVect a_vect);

  Real innerProduct(const RealVect& a_v1,
                    const RealVect& a_v2);

  void AInvTimesVelocity3D(RealVect&            a_answer,
                           const RealVect&      a_velocity,
                           const RealVect&      a_normal,
                           const RealVect&      a_t1,
                           const RealVect&      a_t2,
                           const Vector<Real>&  a_nu,
                           const Vector<Real>&  a_zeta);

  void AInvTimesVelocity2D(RealVect&            a_answer,
                           const RealVect&      a_velocity,
                           const RealVect&      a_normal,
                           const RealVect&      a_t1,
                           const Vector<Real>&  a_nu,
                           const Vector<Real>&  a_zeta);

  void AInvTimesVelocity(RealVect&            a_answer,
                         const RealVect&      a_velocity,
                         const RealVect&      a_normal,
                         const RealVect&      a_t1,
                         const RealVect&      a_t2,
                         const Vector<Real>&  a_nu,
                         const Vector<Real>&  a_zeta);

  void ATimesQ3D(      RealVect& a_answer,
                 const RealVect& a_q,
                 const RealVect& a_normal,
                 const RealVect& a_t1,
                 const RealVect& a_t2,
                 const Real      a_nu);

  void ATimesQ2D(      RealVect& a_answer,
                 const RealVect& a_q,
                 const RealVect& a_normal,
                 const RealVect& a_t1,
                 const Real      a_nu);

  void ATimesQ(      RealVect& a_answer,
               const RealVect& a_q,
               const RealVect& a_normal,
               const RealVect& a_t1,
               const RealVect& a_t2,
               const Real      a_nu);

};
#include "NamespaceFooter.H"
#endif
