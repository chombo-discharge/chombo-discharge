/*!
  @file mf_helmholtz_op.H
  @brief Declaration of multifluid helmholtz operator for embedment in AMRMultiGrid
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _MF_HELMHOLTZ_OP_
#define _MF_HELMHOLTZ_OP_

#include "type_definitions.H"
#include "mfis.H"
#include "jump_bc.H"
#include "MFLevelGrid.H"
#include "electrode.H"
#include "MFQuadCFInterp.H"

#include <LevelDataOps.H>
#include <RefCountedPtr.H>
#include <MFCellFAB.H>
#include <MFFluxFAB.H>
#include <BaseDomainBC.H>
#include <BaseEBBC.H>
#include <EBConductivityOp.H>
#include <DirichletConductivityEBBC.H>
#include <AMRTGA.H>

/*!
  @brief Operator class for solving the variable-coefficient Helmholtz equation with multigrid in two phases
*/
class mf_helmholtz_op : public TGAHelmOp<LevelData<MFCellFAB> >{
public:

  /*!
    @brief Constructor
  */
  mf_helmholtz_op();
  

  /*!
    @brief Destructor
  */
  virtual ~mf_helmholtz_op();

  /*!
    @brief Define function
  */
  virtual void define(const RefCountedPtr<mfis>&                    a_mfis,
		      const RefCountedPtr<BaseDomainBC>&            a_dombc,
		      const RefCountedPtr<LevelData<MFCellFAB> >&   a_aco,
		      const RefCountedPtr<LevelData<MFFluxFAB> >&   a_bco,
		      const RefCountedPtr<LevelData<MFBaseIVFAB> >& a_bco_irreg,
		      const MFQuadCFInterp&                         a_quadcfi,
		      const MFLevelGrid&                            a_mflg_fine,
		      const MFLevelGrid&                            a_mflg,
		      const MFLevelGrid&                            a_mflg_coar,
		      const MFLevelGrid&                            a_mflg_coar_mg,
		      const ProblemDomain&                          a_domain,
		      const bool&                                   a_layout_changed,
		      const bool&                                   a_has_mg_objects,
		      const bool&                                   a_has_fine,
		      const bool&                                   a_has_coar,
		      const int&                                    a_ref_to_fine,
		      const int&                                    a_ref_to_coar,
		      const int&                                    a_relax_type,
		      const int&                                    a_order_ebbc,
		      const IntVect&                                a_ghost_phi,
		      const IntVect&                                a_ghost_rhs,
		      const Real&                                   a_dx,
		      const Real&                                   a_dx_coar,
		      const Real&                                   a_alpha,
		      const Real&                                   a_beta);



  /*!
    @brief Set jump condition
  */
  virtual void set_jump(const RefCountedPtr<LevelData<BaseIVFAB<Real> > >& a_jump);

  /*!
    @brief Set the non-matched boundary conditions
  */
  virtual void set_electrodes(const Vector<electrode>& a_electrodes);

  /*!
    @brief Update boundary conditions before relaxation
  */
  virtual void update_bc();

  /*!
    @brief Reset for TGA
  */
  virtual void setAlphaAndBeta(const Real& a_alpha, const Real& a_beta);

  /*!
    @brief Set diagonal scale of the operator
  */
  virtual void diagonalScale(LevelData<MFCellFAB>& a_rhs);

  /*!
    @brief TGA --no op in this case
  */
  virtual void divideByIdentityCoef(LevelData<MFCellFAB>& a_rhs);

  /*!
    @brief Apply operator without any boundary or coarse-fine
  */
  virtual void applyOpNoBoundary(LevelData<MFCellFAB>&       a_opPhi,
                                 const LevelData<MFCellFAB>& a_phi);

  /*!
    @brief Set time
  */
  virtual void setTime(Real a_oldTime, Real a_mu, Real a_dt);

  /*!
    @brief Compute residual
  */
  virtual void residual(LevelData<MFCellFAB>&       a_lhs,
                        const LevelData<MFCellFAB>& a_phi,
                        const LevelData<MFCellFAB>& a_rhs,
			bool                        a_homogeneous = false);

  /*!
    @brief Apply preconditioner
  */
  virtual void preCond(LevelData<MFCellFAB>&       a_correction,
                       const LevelData<MFCellFAB>& a_residual);

  /*!
    @brief Apply operator
  */
  virtual void applyOp(LevelData<MFCellFAB>&        a_lhs,
		       const LevelData<MFCellFAB>&  a_phi,
		       bool                         a_homogeneous = false);
  
  /*!
    @brief Create a clone of this operator
  */
  virtual void create(LevelData<MFCellFAB>&       a_lhs,
		      const LevelData<MFCellFAB>& a_rhs);

  /*!
    @brief Create a coarsened version of this operator
  */
  virtual void createCoarsened(LevelData<MFCellFAB>&       a_lhs,
			       const LevelData<MFCellFAB>& a_rhs,
			       const int&                  a_refRat);

  /*!
    @brief Assign operator
  */
  virtual void assign(LevelData<MFCellFAB>&       a_lhs,
                      const LevelData<MFCellFAB>& a_rhs);

  /*!
    @brief Dot product
  */
  virtual Real dotProduct(const LevelData<MFCellFAB>& a_1,
                          const LevelData<MFCellFAB>& a_2);

  /*!
    @brief Incremetn
  */
  virtual void incr(LevelData<MFCellFAB>&       a_lhs,
		    const LevelData<MFCellFAB>& a_x,
		    Real                        a_scale);

  /*!
    @brief Multiply
  */
  virtual void axby(LevelData<MFCellFAB>&       a_lhs,
		    const LevelData<MFCellFAB>& a_x,
		    const LevelData<MFCellFAB>& a_y,
		    Real a,
		    Real b);

  /*!
    @brief Scale
  */
  virtual void scale(LevelData<MFCellFAB>& a_lhs, const Real& a_scale);

  /*!
    @brief Compute norm
  */
  virtual Real norm(const LevelData<MFCellFAB>& a_x, int a_ord);

  /*!
    @brief Set to zero
  */
  virtual void setToZero(LevelData<MFCellFAB>& a_x);

  /*!
    @brief Relax
  */
  virtual void relax(LevelData<MFCellFAB>&       a_e,
                     const LevelData<MFCellFAB>& a_residual,
                     int                         iterations);

  /*!
    @brief Create a coarsened data holder
  */
  virtual void createCoarser(LevelData<MFCellFAB>&       a_coarse,
                             const LevelData<MFCellFAB>& a_fine,
                             bool ghosted);

  /*!
     @brief Calculate restricted residual
     @details a_resCoarse[2h] = I[h->2h] (rhsFine[h] - L[h](phiFine[h])
  */
  virtual void restrictResidual(LevelData<MFCellFAB>&       a_resCoarse,
                                LevelData<MFCellFAB>&       a_phiFine,
                                const LevelData<MFCellFAB>& a_rhsFine);

  /*!
     @brief Correct the fine solution based on coarse correction
     @details a_phiThisLevel += I[2h->h](a_correctCoarse)
  */
  virtual void prolongIncrement(LevelData<MFCellFAB>&       a_phiThisLevel,
                                const LevelData<MFCellFAB>& a_correctCoarse);

  /*!
   @brief Returns 1 when there are no coarser AMRLevelOp objects 
  */
  virtual int refToCoarser();

  /*!
    @brief Compute residual
    @details a_residual = a_rhs - L(a_phi, a_phiFine, a_phiCoarse) 
  */
  virtual void AMRResidual(LevelData<MFCellFAB>&       a_residual,
                           const LevelData<MFCellFAB>& a_phiFine,
                           const LevelData<MFCellFAB>& a_phi,
                           const LevelData<MFCellFAB>& a_phiCoarse,
                           const LevelData<MFCellFAB>& a_rhs,
                           bool                        a_homogeneousBC,
                           AMRLevelOp<LevelData<MFCellFAB> >* a_finerOp);


  /*!
    @brief Compute residual assuming no more coarser AMR levels 
  */
  virtual void AMRResidualNC(LevelData<MFCellFAB>&       a_residual,
                             const LevelData<MFCellFAB>& a_phiFine,
                             const LevelData<MFCellFAB>& a_phi,
                             const LevelData<MFCellFAB>& a_rhs,
                             bool                        a_homogeneousBC,
                             AMRLevelOp<LevelData<MFCellFAB> >* a_finerOp);

  /*!
    @brief Compute residual assuming no finer AMR levels
    @details a_residual = a_rhs - L(a_phi, a_phiCoarse)  
  */
  virtual void AMRResidualNF(LevelData<MFCellFAB>&       a_residual,
                             const LevelData<MFCellFAB>& a_phi,
                             const LevelData<MFCellFAB>& a_phiCoarse,
                             const LevelData<MFCellFAB>& a_rhs,
			     bool                        a_homogeneousBC);

  /*!
    @brief AMR operator assuming no more coarser AMR levels 
  */
  virtual void AMROperatorNC(LevelData<MFCellFAB>&       a_LofPhi,
                             const LevelData<MFCellFAB>& a_phiFine,
                             const LevelData<MFCellFAB>& a_phi,
                             bool                        a_homogeneousBC,
                             AMRLevelOp<LevelData<MFCellFAB> >* a_finerOp);

    /*! 
    @brief AMR operator assuming no finer level 
  */
  virtual void AMROperatorNF(LevelData<MFCellFAB>&       a_LofPhi,
                             const LevelData<MFCellFAB>& a_phi,
                             const LevelData<MFCellFAB>& a_phiCoarse,
                             bool a_homogeneousBC);


  /*!
    @brief Apply AMR operator, including coarse-fine matching conditions
  */
  virtual void AMROperator(LevelData<MFCellFAB>&       a_LofPhi,
                           const LevelData<MFCellFAB>& a_phiFine,
                           const LevelData<MFCellFAB>& a_phi,
                           const LevelData<MFCellFAB>& a_phiCoarse,
                           bool                        a_homogeneousBC,
                           AMRLevelOp<LevelData<MFCellFAB> >* a_finerOp);
  



  /*!
    @brief Restrict residual 
    @details a_resCoarse = I[h-2h]( a_residual - L(a_correction, a_coarseCorrection))
    it is assumed that a_resCoarse has already been filled in with the coarse
    version of AMRResidualNF and that this operation is free to overwrite
    in the overlap regions.
  */
  virtual void AMRRestrict(LevelData<MFCellFAB>&       a_resCoarse,
                           const LevelData<MFCellFAB>& a_residual,
                           const LevelData<MFCellFAB>& a_correction,
                           const LevelData<MFCellFAB>& a_coarseCorrection,
                           bool                        a_skip_res);

  /*!
    @brief Prologn residual 
    @details a_correction += I[h->h](a_coarseCorrection)
  */
  virtual void AMRProlong(LevelData<MFCellFAB>&       a_correction,
                          const LevelData<MFCellFAB>& a_coarseCorrection);


  /*!
    @brief Update the AMR residual
  */
  virtual void AMRUpdateResidual(LevelData<MFCellFAB>&       a_residual,
				 const LevelData<MFCellFAB>& a_correction,
				 const LevelData<MFCellFAB>& a_coarseCorrection);


  /*!
    @brief Get the norm
  */
  virtual Real AMRNorm(const LevelData<MFCellFAB>& a_coar_resid,
		       const LevelData<MFCellFAB>& a_fine_resid,
		       const int&                  a_ref_rat,
		       const int&                  a_ord);

  /*!
    @brief Compute kappa norm
  */
  virtual Real kappaNorm(Real&                       a_volume,
			 const LevelData<MFCellFAB>& a_data,
			 int                         a_p) const;

  /*!
    @brief Jacobi iteration
  */
  virtual void levelJacobi( LevelData<MFCellFAB>&       a_phi,
			    const LevelData<MFCellFAB>& a_rhs);

protected:

  /*!
    @brief EB operators
  */
  Vector<RefCountedPtr<EBConductivityOp> > m_ebops;

  /*!
    @brief a-coefficients for each phase
  */
  Vector<RefCountedPtr<LevelData<EBCellFAB> > > m_acoeffs;

  /*!
    @brief b-coefficients for each phase
  */
  Vector<RefCountedPtr<LevelData<EBFluxFAB> > > m_bcoeffs;

  /*!
    @brief b-coefficients for each phase
  */
  Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > > m_bcoeffs_irr;

  /*!
    @brief Aliasing stuff
  */
  Vector<LevelData<EBCellFAB>* > m_alias;

  /*!
    @brief Surface potential
  */
  LevelData<BaseIVFAB<Real> > m_surfpot;

  /*!
    @brief Arithmetic multifluid operators
  */
  LevelDataOps<MFCellFAB> m_ops;

  /*!
    @brief Some extra storage
  */
  LevelData<MFCellFAB> m_tmp;

  /*!
    @brief Some extra storage
  */
  LevelData<MFCellFAB> m_weights;
  
  /*!
    @brief Class for computing BC matching
  */
  RefCountedPtr<jump_bc> m_jumpbc;

  /*!
    @brief Jump condition. Comes in through factory. 
  */
  RefCountedPtr<LevelData<BaseIVFAB<Real> > > m_jump;

  /*!
    @brief Dirichlet boundary conditions on embedded boundaries
    @note Not sure if this can be used directly for both phases; might need a separate for each
  */
  RefCountedPtr<LevelData<BaseIVFAB<Real> > > m_phibc;

  /*!
    @brief Dirichlet conductivity bc for each phase
  */
  Vector<RefCountedPtr<DirichletConductivityEBBC> > m_ebbc;

  /*!
    @brief Electrodes
  */
  Vector<electrode> m_electrodes;

  /*!
    @brief Problem domain
  */
  ProblemDomain m_domain;

  /*!
    @brief Coarser MFLevelGrid for MG
  */
  MFLevelGrid m_mflg_coar_mg;
  
  /*!
    @brief Number of phases.
  */
  int m_phases;

  /*!
    @brief Number of components
  */
  int m_ncomp;

  /*!
    @brief Relaxation type
  */
  int m_relax;

  /*!
    @brief Refinement factor to coarser level
  */
  int m_ref_to_coarser;

  /*!
    @brief Number of ghost cells
  */
  IntVect m_ghost_rhs;

  /*!
    @brief Number of ghost cells
  */
  IntVect m_ghost_phi;
  
  /*!
    @brief alpha coefficient
  */
  Real m_alpha;

  /*!
    @brief beta-coefficient
  */
  Real m_beta;
  
private:
  
};
#endif
