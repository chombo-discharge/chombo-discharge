/*!
  @file nwomfconductivityopfactory.H
  @brief Factory class for creating multifluid helmholtz operators
  @author Robert Marskar
  @todo Move implementations to nwomfconductivityopfactory.H
*/

#ifndef _NWOMFCONDUCTIVITYOPFACTORY_
#define _NWOMFCONDUCTIVITYOPFACTORY_

#include "nwomfconductivityop.H"
#include "mfis.H"
#include "mfalias.H"
#include "MFLevelGrid.H"
#include "NWOMFQuadCFInterp.H"
#include "data_ops.H"
#include "electrode.H" 

#include <EBLevelDataOps.H>
#include <AMRMultiGrid.H>
#include <MFAliasFactory.H>
#include <ebcoarseaverage.H>
#include <BaseDomainBC.H>
#include <BaseEBBC.H>
#include <ProblemDomain.H>
#include <RefCountedPtr.H>
#include <EBIndexSpace.H>
#include <EBQuadCFInterp.H>
#include <MFCellFAB.H>

/*!
  @brief Factory class for generationg nwomfconductivityop's. 
  @details The jump condition is of the Neumann type in the form
  \f[
  a\partial_n\phi_1 - b\partial_n\phi_2 = c\sigma,
  \f]
  where \f$a\f$, \f$b\f$ and \f$c\f$ are variables.  
*/
class nwomfconductivityopfactory : public AMRLevelOpFactory<LevelData<MFCellFAB> > {
public:
  
  /*!
    @brief Factory class for creating multifluid helmholtz operators
    @note This class generates its own EBConductivityOps as well as it's own interpolation utilities. 
    @details This class only supports two-phase solutions with single components
  */
  nwomfconductivityopfactory(const RefCountedPtr<mfis>&             a_mfis,
			     const Vector<MFLevelGrid>&                a_mflg,
			     const Vector<NWOMFQuadCFInterp>&          a_mfquadcfi,
			     const Vector<int>&                        a_ref_rat,
			     const Vector<DisjointBoxLayout>&          a_grids,
			     const MFAMRCellData&                      a_aco,
			     const MFAMRFluxData&                      a_bco,
			     const MFAMRIVData&                        a_bco_irreg,
			     const Real&                               a_alpha,
			     const Real&                               a_beta,
			     const Real&                               a_coarsest_dx,
			     const ProblemDomain&                      a_coarsest_domain,
			     const RefCountedPtr<BaseDomainBCFactory>& a_dombc,
			     const RealVect&                           a_origin,
			     const IntVect&                            a_ghost_phi,
			     const IntVect&                            a_ghost_rhs,
			     const int                                 a_ebbc_order,
			     const int                                 a_drop_bottom = -1,
			     int                                       a_num_levels = -1);

  /*!
    @brief Destructor
  */
  ~nwomfconductivityopfactory();

  /*!
    @brief Set when to drop to the bottom solver
  */
  virtual void set_bottom_drop(const int a_bottom_drop);

  /*!
    @brief Set when to drop to the bottom solver
  */
  virtual void set_relax_type(const int a_relax_type);

  /*!
    @brief MG stuff
  */
  virtual void set_max_box_size(const int a_max_box_size);

  /*!
    @brief Set boundary condition order
  */
  virtual void set_ebbc_order(const int a_order);

  /*!
    @brief Coarsen coefficients for multigrid
  */
  virtual void coarsen_coefficients(LevelData<MFCellFAB>&         a_aco_coar,
				    LevelData<MFFluxFAB>&         a_bco_coar,
				    LevelData<MFBaseIVFAB>&       a_bco_irreg_coar,
				    const MFLevelGrid&            a_mflg_coar,
				    const MFLevelGrid&            a_mflg_fine,
				    const LevelData<MFCellFAB>&   a_aco_fine,
				    const LevelData<MFFluxFAB>&   a_bco_fine,
				    const LevelData<MFBaseIVFAB>& a_bco_irreg_fine,
				    const int&                    a_ref_to_depth);
  
  /*!
    @brief Set jump conditions. 
    @details This call must be performed before each solve; this call updates the surface charge density on the jump interface, 
    including MG levels. 
  */
  virtual void set_jump(const EBAMRIVData& a_sigma, const Real& a_scale);

  /*!
    @brief Set simplified jump condition
  */
  virtual void set_jump(const Real& a_sigma, const Real& a_scale);

  /*!
    @brief Electrodes, this stuff is used for setting appropriate Dirichlet boundary conditions
  */
  virtual void set_electrodes(const Vector<electrode>& a_electrodes, const RefCountedPtr<BaseBCFuncEval> m_potential);

  /*!
    @brief Average down the surface charge on AMR levels
  */
  virtual void average_down_amr();

  /*!
    @brief Average down the surface charge on lower levels of multigrid
  */
  virtual void average_down_mg();

  /*!
    @brief Reclaim
  */
  virtual void reclaim(MGLevelOp<LevelData<EBCellFAB> >* a_reclaim);

  /*!
    @brief Reclaim
  */
  virtual void AMRreclaim(nwomfconductivityop* a_reclaim);

  /*!
    @brief Define things for lower levels of multigrid
  */
  virtual void define_multigrid_stuff();

  /*!
    @brief Allocate internal storage for jump cells on the AMR hierarchy.
  */
  virtual void define_jump_stuff();

  /*!
    @brief Factory method
  */
  virtual MGLevelOp<LevelData<MFCellFAB> >* MGnewOp(const ProblemDomain& a_domain_fine, int a_depth, bool a_homo_only = true);

  /*!
    @brief Factory method
  */
  virtual AMRLevelOp<LevelData<MFCellFAB> >* AMRnewOp(const ProblemDomain& a_fineindexspace);

  /*!
    @brief Refinement ratio between this level and the coarser level. Returns 1 when there are no coarser levels
  */
  virtual int refToFiner(const ProblemDomain& a_domain) const;

protected:

  /*!
    @brief Index spaces
  */
  RefCountedPtr<mfis> m_mfis;

  /*!
    @brief Potential function
  */
  RefCountedPtr<BaseBCFuncEval> m_potential;

  /*!
    @brief a-coefficients
  */
  MFAMRCellData m_aco;

  /*!
    @brief b-coefficients
  */
  MFAMRFluxData m_bco;

  /*!
    @brief b-coefficients
  */
  MFAMRIVData m_bco_irreg;

  /*!
    @brief Jump condition
  */
  EBAMRIVData m_jump;

  /*!
    @brief Electrodes
  */
  Vector<electrode> m_electrodes;

  /*!
    @brief EBLevelGrids (in MFLevelGrid form)
  */
  Vector<MFLevelGrid> m_mflg;

  /*!
    @brief Interpolation utilities
  */
  Vector<NWOMFQuadCFInterp> m_mfquadcfi;

  /*!
    @brief Refinement ratios
  */
  Vector<int> m_ref_rat;

  /*!
    @brief AMR grids
  */
  Vector<DisjointBoxLayout> m_grids;

  /*!
    @brief Problem domains
  */
  Vector<ProblemDomain> m_domains;

  /*!
    @brief Resolutions
  */
  Vector<Real> m_dx;

  /*!
    @brief Interface jump cells
  */
  Vector<RefCountedPtr<LayoutData<IntVectSet> > > m_jumpcells;

  /*!
    @brief Averaging operator. 
  */
  Vector<RefCountedPtr<ebcoarseaverage> > m_aveop;

  /*!
    @brief Averaging operator for lower levels of multigrid. This is needed for coarsening data-based boundary conditions
  */
  Vector<Vector<RefCountedPtr<ebcoarseaverage> > > m_aveop_mg;

  /*!
    @brief Grids for multigrid
  */
  Vector<Vector<DisjointBoxLayout> > m_grids_mg;

  /*!
    @brief Domains for multigrid
  */
  Vector<Vector<ProblemDomain> > m_domains_mg;

  /*!
    @brief Coarsened surface charge for lower levels of multigrid
  */
  Vector<EBAMRIVData> m_jump_mg;

  /*!
    @brief a-coefficient for multigrid
  */
  Vector<MFAMRCellData> m_aco_mg;

  /*!
    @brief b-coefficient for multigrid
  */
  Vector<MFAMRFluxData> m_bco_mg;

  /*!
    @brief Irregular b-coefficient for multigrid
  */
  Vector<MFAMRIVData> m_bco_irreg_mg;

  /*!
    @brief MFLevelGrids for MultiGrid
  */
  Vector<Vector<MFLevelGrid> > m_mflg_mg;

  /*!
    @brief Interface jump cells for lower levels MG
  */
  Vector<Vector<RefCountedPtr<LayoutData<IntVectSet> > > > m_jumpcells_mg;

  /*!
    @brief Has multigrid objects or not
  */
  std::vector<bool> m_has_mg_objects;

  /*!
    @brief Layout changed or not
  */
  std::vector<bool> m_layout_changed;
  
  /*!
    @brief Layout changed or not
  */
  std::vector<std::vector<bool> > m_layout_changed_mg;

  /*!
    @brief Domain boundary conditions
  */
  RefCountedPtr<BaseDomainBCFactory> m_dombc;

  /*!
    @brief Number of AMR levels
  */
  int m_num_levels;

  /*!
    @brief Boundary conditiosn order
  */
  int m_ebbc_order;

  /*!
    @brief Drop to bottom solver
  */
  int m_test_ref;

  /*!
    @brief Internal stuff
  */
  int m_max_box_size;

  /*!
    @brief Relaxation type
  */
  int m_relax_type;

  /*!
    @brief Ghost cells for solution vector
  */
  IntVect m_ghost_phi;

  /*!
    @brief Ghost cells for solution vector
  */
  IntVect m_ghost_rhs;

  /*!
    @brief alpha-coefficient
  */
  Real m_alpha;

  /*!
    @brief beta-coefficient
  */
  Real m_beta;

  /*!
    @brief Origin
  */
  RealVect m_origin;

  
private:
  
  /*!
    @brief Disallowed weak constructor
  */
  nwomfconductivityopfactory(){
    MayDay::Error("nwomfconductivityopfactory:: invalid constructor");
  }
};

#endif
