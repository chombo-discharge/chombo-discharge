#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// BVS, Nov 3, 2004

#ifndef _MFPOISSONOP_H_
#define _MFPOISSONOP_H_

#include "AMRMultiGrid.H"
#include "REAL.H"
#include "Box.H"
#include "MFCellFAB.H"
#include "EBAMRPoissonOp.H"
#include "MFIndexSpace.H"
#include "MFCellFAB.H"
#include "MFFluxFAB.H"
#include "RefCountedPtr.H"
#include "BaseDomainBC.H"
#include "BaseBCValue.H"
#include "BaseEBBC.H"
#include "BaseIVFAB.H"
#include "InterfaceJump.H"
#include "LevelDataOps.H"
#include "AMRTGA.H"
#include "NamespaceHeader.H"

/// Multifluid poisson operator -- computes (alpha + div(Beta Grad))
/**
 */
class MFPoissonOp : public TGAHelmOp< LevelData<MFCellFAB> >
{
public:
  /// Default constructor
  MFPoissonOp()
  {
  }

  /// destructor
  virtual ~MFPoissonOp();
  /**
     \name LinearOp functions */
  /*@{*/

  /// full define function for AMRLevelOp with both coarser and finer levels
  void define(const MFIndexSpace&      a_mfis,
              int                      a_ncomp,
              const DisjointBoxLayout& a_grids,
              const DisjointBoxLayout& a_gridsCoarMG,
              const bool& a_hasMGObjects,
              const bool& a_layoutChanged,
              const DisjointBoxLayout& a_gridsFiner,
              const DisjointBoxLayout& a_gridsCoarser,
              const RealVect&          a_dxLevel,
              int                      a_refRatio,
              int                      a_refRatioFiner,
              const ProblemDomain&     a_domain,
              const Vector<RefCountedPtr<BaseDomainBC> >& a_bc,
              const IntVect& a_ghostPhi,
              const IntVect& a_ghostRHS,
              bool hasCoarser,
              bool hasFiner  ,
              const Vector<Real>& a_alpha,
              const Vector<Real>& a_beta
              );

  // (TGA) reset constants in the equation
  virtual void setAlphaAndBeta(const Real& a_alpha,
                               const Real& a_beta);

  /// (TGA) set diagonal scale of the operator
  virtual void diagonalScale(LevelData<MFCellFAB>& a_rhs);

  ///TGA --no op in this case
  virtual void divideByIdentityCoef(LevelData<MFCellFAB>& a_rhs)
  {
  }
  /// (TGA) apply operator without any boundary or coarse-fine
  ///  boundary conditions and no finer level
  virtual void applyOpNoBoundary(LevelData<MFCellFAB>&       a_opPhi,
                                 const LevelData<MFCellFAB>& a_phi);

  //conforms to AMRTGA
  virtual void setTime(Real a_time);

  ///
  virtual void residual(  LevelData<MFCellFAB>& a_lhs,
                          const LevelData<MFCellFAB>& a_phi,
                          const LevelData<MFCellFAB>& a_rhs,
                          bool a_homogeneous = false);

  /// Apply the preconditioner
  virtual void preCond(   LevelData<MFCellFAB>& a_correction,
                          const LevelData<MFCellFAB>& a_residual);

  ///
  virtual void applyOp(   LevelData<MFCellFAB>&        a_lhs,
                          const LevelData<MFCellFAB>&  a_phi,
                          DataIterator&                a_dit,
                          bool                         a_homogeneous = false);

  ///
  virtual void applyOp(   LevelData<MFCellFAB>&        a_lhs,
                          const LevelData<MFCellFAB>&  a_phi,
                          bool                         a_homogeneous = false);

  /// create a clone of this MFPoissonOp
  virtual void create(    LevelData<MFCellFAB>& a_lhs,
                          const LevelData<MFCellFAB>& a_rhs);

  ///
  virtual void createCoarsened(    LevelData<MFCellFAB>& a_lhs,
                                   const LevelData<MFCellFAB>& a_rhs,
                                   const int& a_refRat);

  ///
  virtual void assign(    LevelData<MFCellFAB>& a_lhs,
                          const LevelData<MFCellFAB>& a_rhs) ;

  ///
  virtual Real dotProduct(const LevelData<MFCellFAB>& a_1,
                          const LevelData<MFCellFAB>& a_2) ;

  ///
  virtual void incr( LevelData<MFCellFAB>& a_lhs,
                     const LevelData<MFCellFAB>& a_x,
                     Real a_scale) ;

  ///
  virtual void axby( LevelData<MFCellFAB>& a_lhs,
                     const LevelData<MFCellFAB>& a_x,
                     const LevelData<MFCellFAB>& a_y,
                     Real a, Real b) ;

  ///
  virtual void scale(LevelData<MFCellFAB>& a_lhs, const Real& a_scale) ;

  ///
  virtual Real norm(const LevelData<MFCellFAB>& a_x, int a_ord);

  ///
  virtual void setToZero( LevelData<MFCellFAB>& a_x);
  /*@}*/

  /**
     \name MGLevelOp functions */
  /*@{*/

  ///
  virtual void relax(LevelData<MFCellFAB>& a_e,
                     const LevelData<MFCellFAB>& a_residual,
                     int iterations);

  ///
  virtual void createCoarser(LevelData<MFCellFAB>& a_coarse,
                             const LevelData<MFCellFAB>& a_fine,
                             bool ghosted);

  ///
  /**
     calculate restricted residual
     a_resCoarse[2h] = I[h->2h] (rhsFine[h] - L[h](phiFine[h])
  */
  virtual void restrictResidual(LevelData<MFCellFAB>& a_resCoarse,
                                LevelData<MFCellFAB>& a_phiFine,
                                const LevelData<MFCellFAB>& a_rhsFine);


  ///
  /**
     correct the fine solution based on coarse correction
     a_phiThisLevel += I[2h->h](a_correctCoarse)
  */
  virtual void prolongIncrement(LevelData<MFCellFAB>& a_phiThisLevel,
                                const LevelData<MFCellFAB>& a_correctCoarse);

  /*@}*/

  /**
     \name AMRLevelOp functions */
  /*@{*/

  ///
  /** returns 1 when there are no coarser AMRLevelOp objects */
  virtual int refToCoarser()
  {
    return m_refToCoarser;
  }

  ///
  /** a_residual = a_rhs - L(a_phi, a_phiFine, a_phiCoarse) */
  virtual void AMRResidual(LevelData<MFCellFAB>& a_residual,
                           const LevelData<MFCellFAB>& a_phiFine,
                           const LevelData<MFCellFAB>& a_phi,
                           const LevelData<MFCellFAB>& a_phiCoarse,
                           const LevelData<MFCellFAB>& a_rhs,
                           bool a_homogeneousBC,
                           AMRLevelOp<LevelData<MFCellFAB> >* a_finerOp);
  ///
  /** residual assuming no more coarser AMR levels */
  virtual void AMRResidualNC(LevelData<MFCellFAB>& a_residual,
                             const LevelData<MFCellFAB>& a_phiFine,
                             const LevelData<MFCellFAB>& a_phi,
                             const LevelData<MFCellFAB>& a_rhs,
                             bool a_homogeneousBC,
                             AMRLevelOp<LevelData<MFCellFAB> >* a_finerOp);

  /** a_residual = a_rhs - L(a_phi, a_phiCoarse)  */
  virtual void AMRResidualNF(LevelData<MFCellFAB>& a_residual,
                             const LevelData<MFCellFAB>& a_phi,
                             const LevelData<MFCellFAB>& a_phiCoarse,
                             const LevelData<MFCellFAB>& a_rhs, bool a_homogeneousBC);


  ///
  /** apply AMR operator, including coarse-fine matching conditions */
  virtual void AMROperator(LevelData<MFCellFAB>& a_LofPhi,
                           const LevelData<MFCellFAB>& a_phiFine,
                           const LevelData<MFCellFAB>& a_phi,
                           const LevelData<MFCellFAB>& a_phiCoarse,
                           bool a_homogeneousBC,
                           AMRLevelOp<LevelData<MFCellFAB> >* a_finerOp);
  ///
  /** AMR operator assuming no more coarser AMR levels */
  virtual void AMROperatorNC(LevelData<MFCellFAB>& a_LofPhi,
                             const LevelData<MFCellFAB>& a_phiFine,
                             const LevelData<MFCellFAB>& a_phi,
                             bool a_homogeneousBC,
                             AMRLevelOp<LevelData<MFCellFAB> >* a_finerOp);

  /** AMR operator assuming no finer level */
  virtual void AMROperatorNF(LevelData<MFCellFAB>& a_LofPhi,
                             const LevelData<MFCellFAB>& a_phi,
                             const LevelData<MFCellFAB>& a_phiCoarse,
                             bool a_homogeneousBC);




  ///
  /** a_resCoarse = I[h-2h]( a_residual - L(a_correction, a_coarseCorrection))
      it is assumed that a_resCoarse has already been filled in with the coarse
      version of AMRResidualNF and that this operation is free to overwrite
      in the overlap regions.
  */
  virtual void AMRRestrict(LevelData<MFCellFAB>& a_resCoarse,
                           const LevelData<MFCellFAB>& a_residual,
                           const LevelData<MFCellFAB>& a_correction,
                           const LevelData<MFCellFAB>& a_coarseCorrection,
                           bool                        a_skip_res);

  ///
  /** a_correction += I[h->h](a_coarseCorrection) */
  virtual void AMRProlong(LevelData<MFCellFAB>& a_correction,
                          const LevelData<MFCellFAB>& a_coarseCorrection);

  ///
  /** a_residual = a_residual - L(a_correction, a_coarseCorrection) */
  virtual void AMRUpdateResidual(LevelData<MFCellFAB>& a_residual,
                                 const LevelData<MFCellFAB>& a_correction,
                                 const LevelData<MFCellFAB>& a_coarseCorrection);

  ///
  /**
     compute norm over all cells on coarse not covered by finer
   */
  virtual Real AMRNorm(const LevelData<MFCellFAB>& a_coarseResid,
                       const LevelData<MFCellFAB>& a_fineResid,
                       const int&                  a_refRat,
                       const int&                  a_ord);

  /*@}*/

  /// set the jump conditions at the multifluid interface
  void setJump(const Real& gD, const Real& gN);
  void setJump(const RealVect& a_gD, const RealVect& a_gN);

  ///
  /**
      version where we want to specify the interface jumps
      through analytic function
  */
  void setJump(const Vector< RefCountedPtr<BaseBCValue> >& a_phiValVect,
               const Vector< RefCountedPtr<BaseBCValue> >& a_flxValVect);

  void setVal(LevelData<MFCellFAB>& a_phi, const Vector<Real> a_values) const;

  ///
  EBAMRPoissonOp* ebOp(int phase)
  {
    return m_ebops[phase];
  }

  ///
  /**
     return total flux across all irregular faces in a_phase
   */
  Real totalBoundaryFlux(int                          a_phase,
                         const LevelData<MFCellFAB>&  a_phi,
                         Real                         a_factor = 1.0,
                         bool                         a_divideByTotalArea = true);

  void getBoundaryValues(LevelData<MFCellFAB>&  a_phi,
                         LevelData<MFCellFAB>&  a_dPhi_dN,
                         Real                   a_invalidVal = 1.2345678e90);

  Real exactBoundaryFlux(int                         a_phase,
                         RefCountedPtr<BaseBCValue>  a_flxVal,
                         RealVect&                   a_origin,
                         const Real&                 a_time);


  ///
  /**
     dump stencil as matrix to stdout
   */
  class StencilIndex
  {
  public:
    StencilIndex(VolIndex& a_vof,
                 int&      a_phase)
      {
        m_vof   = a_vof;
        m_phase = a_phase;
      }

    StencilIndex()
    {
      MayDay::Error("No weak construction of StencilIndex class.");
    }

    StencilIndex& operator= (const StencilIndex& a_sin)
    {
      return *this;
    }

    VolIndex vof() const
    {
      return m_vof;
    }

    int phase() const
    {
      return m_phase;
    }

  protected:
    VolIndex m_vof;
    int      m_phase;
  };

  class StencilIndexComparator
  {
  public:
    bool operator() (const StencilIndex& a_s1, const StencilIndex& a_s2) const
    {
      /*
      int p1 = a_s1.phase();
      int p2 = a_s2.phase();
      if (p1 == p2)
        {
          return (a_s1.vof()<a_s2.vof());
        }
      else
        {
          const IntVect& iv1 = a_s1.vof().gridIndex();
          const IntVect& iv2 = a_s2.vof().gridIndex();
          if (iv1 == iv2)
            {
              return (p1<p2);
            }
          else
            {
              return (iv1<iv2);
            }
        }
      */
      const IntVect& iv1 = a_s1.vof().gridIndex();
      const IntVect& iv2 = a_s2.vof().gridIndex();
      for (int idir=0; idir<SpaceDim; ++idir)
        {
          if (iv1[idir] != iv2[idir])
            {
              return (iv1[idir]<iv2[idir]);
            }
        }
      int p1 = a_s1.phase();
      int p2 = a_s2.phase();
      if (p1 == p2)
        {
          int ci1 = a_s1.vof().cellIndex();
          int ci2 = a_s2.vof().cellIndex();
          return (ci1<ci2);
        }
      else
        {
          return (p1<p2);
        }
    }
  };

  void dumpStencilMatrix();
  void dumpReferenceStencilMatrix();

  ///
  /**
     return my domain
   */
  ProblemDomain getDomain()
  {
    return m_domain;
  }

  ///
  /**
     return fluxes at cell boundaries. at present only does regular faces.
   */
  void getFlux(MFFluxFAB&                   a_flux,
               const LevelData<MFCellFAB>&  a_data,
               const Box&                   a_grid,
               const DataIndex&             a_dit,
               Real                         a_scale);

  ///
  int m_relax;

  Vector<IntVect> m_colors;

protected:
  Vector<Real> m_beta;
  Vector<Real> m_bCoef;
  Vector<Real> m_alpha;
  Vector<Real> m_aCoef;

  RealVect                m_dx;
  RealVect                m_dxCrse;
  ProblemDomain           m_domain;
  Vector<EBAMRPoissonOp*>          m_ebops;
  Vector<LevelData<EBCellFAB>* >   m_alias;
  LevelDataOps<MFCellFAB> m_ops;
  LevelData<MFCellFAB>    m_tmp; // temp of rhs/residual;
  LevelData<MFCellFAB>    m_weights;

  LevelData<BaseIVFAB<Real> >     m_boundaryD[2];
  LevelData<BaseIVFAB<Real> >     m_boundaryN[2];
  InterfaceJump                   m_jump;// only one for now, just two-phase.

  void computeBoundaryN(const LevelData<MFCellFAB>& a_phi, bool a_homogeneous);

  void levelJacobi( LevelData<MFCellFAB>&       a_phi,
                    const LevelData<MFCellFAB>& a_rhs);

  void levelGSRB( LevelData<MFCellFAB>&       a_phi,
                  const LevelData<MFCellFAB>& a_rhs);

  void levelMulticolorGS( LevelData<MFCellFAB>& a_phi,
                          const LevelData<MFCellFAB>& a_rhs);

  Real kappaNorm(Real&                       a_volume,
                 const LevelData<MFCellFAB>& a_data,
                 int                         a_p) const;

  int m_refToCoarser;
  int m_refToFiner;
  int m_phases;
  int m_ncomp;
  IntVect m_ghostPhi, m_ghostRHS;
};

#include "NamespaceFooter.H"
#endif
