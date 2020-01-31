#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MFPOISSONOPFACTORY_H_
#define _MFPOISSONOPFACTORY_H_

#include "REAL.H"
#include "Box.H"
#include "FArrayBox.H"
#include "Vector.H"
#include <map>
#include "RefCountedPtr.H"

#include "AMRMultiGrid.H"

#include "EBIndexSpace.H"
#include "EBCellFAB.H"
#include "EBCellFactory.H"

#include "EBLevelDataOps.H"
#include "BaseEBBC.H"
#include "BaseDomainBC.H"
#include "BaseBCValue.H"
#include "CFIVS.H"
#include "EBFluxRegister.H"
#include "EBMGAverage.H"
#include "EBMGInterp.H"
#include "EBCoarsen.H"
#include "PolyGeom.H"
#include "MFPoissonOp.H"
#include "EBLevelGrid.H"
#include "NamespaceHeader.H"


///
/**
   Factory class to generate MFPoissonOps.  This follows the AMRLevelOpFactory interface.
*/
class MFPoissonOpFactory: public AMRLevelOpFactory<LevelData<MFCellFAB> >
{
public:
  ///
  virtual ~MFPoissonOpFactory();

  ///
  /**
     a_dilbos : layouts at each AMR level \\
     a_domainFactory : domain boundary conditions \\
     a_ebBCFactory:    eb boundary conditions \\
     a_dxCoarse:      grid spacing at coarsest level \\
     a_origin:        offset to lowest corner of the domain \\
     a_refRatio:     refinement ratios. refRatio[i] is between levels i and i+1 \\
     a_preCondIters:  number of iterations to do for pre-conditioning \\
     a_alpha:         coefficent of identity \\
     a_beta:          coefficient of laplacian.\\
     a_time:          time for boundary conditions \\
     a_ghostCellsPhi:  Number of ghost cells in phi, correction (typically one)\\
     a_ghostCellsRhs:  Number of ghost cells in RHS, residual, lphi (typically zero)\\
     Ghost cell arguments are there for caching reasons.  Once you set them, an error is thrown if
     you send in data that does not match.   Use numlevels = -1 if you want to use the
     size of the vectors for numlevels.
  */
  MFPoissonOpFactory(const RefCountedPtr<MFIndexSpace>&            a_mfis,
                     const Vector<DisjointBoxLayout>&              a_dilbos,
                     const Vector<int>&                            a_refRatio,
                     const ProblemDomain&                          a_domainCoar,
                     const RealVect&                               a_dxCoarse,
                     const RealVect&                               a_origin,
                     const Vector<RefCountedPtr<BaseDomainBC> >&   a_bc,
                     const Vector<Real>&                           a_alpha,
                     const Vector<Real>&                           a_beta,
                     const int&                                    a_ncomp,
                     const IntVect&                                a_ghostCellsPhi,
                     const IntVect&                                a_ghostCellsRHS,
                     int a_numLevels = -1);

  ///
  void setJump(const Real& a_gD,
               const Real& a_gN);

  ///
  void setJump(const RealVect& a_gD,
               const RealVect& a_gN);

  ///
  void setJump(const Vector< RefCountedPtr<BaseBCValue> >& a_phiValVect,
               const Vector< RefCountedPtr<BaseBCValue> >& a_flxValVect);

  ///
  virtual MGLevelOp<LevelData<MFCellFAB> >*
  MGnewOp(const ProblemDomain& a_FineindexSpace,
          int                  a_depth,
          bool                 a_homoOnly = true);

  MFPoissonOp* createOperator(const DisjointBoxLayout&             a_dilboMGLevel,
                              const DisjointBoxLayout&             a_dilboCoarMG,
                              const ProblemDomain&           a_domainMGLevel,
                              const bool&                    a_hasMGObjects,
                              const bool&                    a_layoutChanged,
                              const RealVect&                a_dxMGLevel,
                              const RealVect&                a_dxCoar,
                              const int&                     a_whichLevel,
                              const int&                     a_mgLevel);
  ///
  virtual void reclaim(MGLevelOp<LevelData<EBCellFAB> >* a_reclaim);

  ///
  virtual AMRLevelOp<LevelData<MFCellFAB> >*
  AMRnewOp(const ProblemDomain& a_FineindexSpace);

  ///
  virtual void AMRreclaim(MFPoissonOp* a_reclaim);

  ///
  /** Refinement ratio between this level and coarser level.
      Returns 1 when there are no coarser AMRLevelOp objects */
  virtual int refToFiner(const ProblemDomain& a_domain) const;

  ///
  /**
     testRef is the size of the coarsest domain allowed in multigrid. If testRef=2,
     then the coarsest domain in multigrid will be 2x2(x2)
   **/
  static void setTestRef(int a_testRef)
  {
    s_testRef = a_testRef;
  }

  ///
  /**
   **/
  static void setRelaxType(int a_relaxType)
  {
    s_relaxType = a_relaxType;
  }

protected:
  static int s_testRef, s_relaxType;

  bool                                      m_isDefined;
  bool                                      m_isScalarJump, m_analyticJump;
  int                                       m_numLevels;
  int                                       m_ncomp;

  RefCountedPtr<MFIndexSpace>               m_mfis;
  Vector<DisjointBoxLayout>                 m_dilboVec;
  Vector<ProblemDomain>                     m_domainVec;
  Vector< Vector<DisjointBoxLayout> >       m_dilboVecMG;
  Vector< Vector<ProblemDomain> >           m_domainVecMG;
  std::vector<bool>                         m_hasMGObjects;
  std::vector<bool>                         m_layoutChanged;
  std::vector< std::vector<bool> >          m_layoutChangedMG;
  Vector<int>                               m_refRatioVec;
  Vector<RealVect>                          m_dxVec;
  RealVect                                  m_origin;
  Vector<Real>                              m_alpha;
  Vector<Real>                              m_beta;
  const IntVect                             m_ghostCellsPhi;
  const IntVect                             m_ghostCellsRHS;
  Real                                      m_scalarGD, m_scalarGN;
  RealVect                                  m_vectorGD, m_vectorGN;
  Vector< RefCountedPtr<BaseBCValue> >      m_phiValVect, m_flxValVect;
  Vector< RefCountedPtr<BaseDomainBC> >     m_bc;

private:
  ///weak construction bad
  MFPoissonOpFactory()
  {
    MayDay::Error("invalid operator");
  }

  //copy constructor and operator= disallowed for all the usual reasons
  MFPoissonOpFactory(const MFPoissonOpFactory& a_opin)
  {
    MayDay::Error("invalid operator");
  }

  void operator=(const MFPoissonOpFactory& a_opin)
  {
    MayDay::Error("invalid operator");
  }
};

#include "NamespaceFooter.H"
#endif
