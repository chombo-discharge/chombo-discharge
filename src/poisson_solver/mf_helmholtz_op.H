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

#include <RefCountedPtr.H>
#include <MFCellFAB.H>
#include <MFFluxFAB.H>
#include <BaseDomainBC.H>
#include <BaseEBBC.H>
#include <EBConductivityOp.H>

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
  virtual void define(const RefCountedPtr<mfis>&         a_mfis,
		      const RefCountedPtr<BaseDomainBC>& a_dombc,
		      const RefCountedPtr<jump_bc>&      a_jumpbc,
		      const MFLevelGrid&                 a_mflg_fine,
		      const MFLevelGrid&                 a_mflg,
		      const MFLevelGrid&                 a_mflg_coar,
		      const MFLevelGrid&                 a_mflg_coar_mg,
		      const DisjointBoxLayout&           a_dbl,
		      const DisjointBoxLayout&           a_dbl_finer,
		      const DisjointBoxLayout&           a_dbl_coarser,
		      const DisjointBoxLayout&           a_dbl_coar_mg,
		      const bool&                        a_layout_changed,
		      const bool&                        a_has_mg_objects,
		      const bool&                        a_has_coarser,
		      const bool&                        a_has_finer,
		      const int&                         a_ref_ratio,
		      const int&                         a_ref_ratio_to_finer,
		      const IntVect&                     a_ghost_phi,
		      const IntVect&                     a_ghost_rhs,
		      const Real&                        a_dx,
		      const Real&                        a_alpha,
		      const Real&                        a_beta);

  /*!
    @brief Reset for TGA
  */
  virtual void setAlphaAndBeta(const Real& a_alpha, const Real& a_beta);


  /*!
    @brief Update the AMR residual
  */
  virtual void AMRUpdateResidual(LevelData<MFCellFAB>&       a_residual,
				 const LevelData<MFCellFAB>& a_correction,
				 const LevelData<MFCellFAB>& a_coarseCorrection);


  virtual Real AMRNorm(const LevelData<MFCellFAB>& a_coar_resid,
		       const LevelData<MFCellFAB>& a_fine_resid,
		       const int&                  a_ref_rat,
		       const int&                  a_ord);
protected:

  Vector<RefCountedPtr<EBConductivityOp> > m_ebops;

  /*!
    @brief Aliasing stuff
  */
  Vector<LevelData<EBCellFAB>* > m_alias;

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
