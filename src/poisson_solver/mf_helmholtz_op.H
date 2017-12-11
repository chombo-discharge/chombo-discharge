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

#include <RefCountedPtr.H>
#include <MFCellFAB.H>
#include <MFFluxFAB.H>
#include <BaseDomainBC.H>
#include <BaseEBBC.H>
#include <EBConductivityOp.H>
#include <DirichletConductivityEBBC.H>

/*!
  @brief Operator class for solving the variable-coefficient Helmholtz equation with multigrid in two phases
*/
class mf_helmholtz_op {
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
		      const RefCountedPtr<MFQuadCFInterp>&          a_quadcfi,
		      const RefCountedPtr<BaseDomainBC>&            a_dombc,
		      const RefCountedPtr<LevelData<MFCellFAB> >&   a_aco,
		      const RefCountedPtr<LevelData<MFFluxFAB> >&   a_bco,
		      const RefCountedPtr<LevelData<MFBaseIVFAB> >& a_bco_irreg,
		      const MFLevelGrid&                            a_mflg_fine,
		      const MFLevelGrid&                            a_mflg,
		      const MFLevelGrid&                            a_mflg_coar,
		      const MFLevelGrid&                            a_mflg_coar_mg,
		      const DisjointBoxLayout&                      a_dbl,
		      const DisjointBoxLayout&                      a_dbl_finer,
		      const DisjointBoxLayout&                      a_dbl_coarser,
		      const DisjointBoxLayout&                      a_dbl_coar_mg,
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
    @brief Reset for TGA
  */
  virtual void setAlphaAndBeta(const Real& a_alpha, const Real& a_beta);

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
  Vector<LevelData<EBCellFAB>* > m_cell_alias;

  /*!
    @brief Aliasing stuff
  */
  Vector<LevelData<EBFluxFAB>* > m_flux_alias;

  /*!
    @brief Aliasing stuff
  */
  Vector<LevelData<BaseIVFAB<Real> >* > m_eb_alias;

  /*!
    @brief Surface potential
  */
  LevelData<BaseIVFAB<Real> > m_surfpot;
  
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
    @brief Number of phases.
  */
  int m_phases;
  
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
