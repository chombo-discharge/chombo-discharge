/*!
  @file mf_helmholtz_opfactory.H
  @brief Factory class for creating multifluid helmholtz operators
  @author Robert Marskar
  @todo Move implementations to mf_helmholtz_opfactory.H
*/

#ifndef _MF_HELMHOLTZ_OPFACTORY_
#define _MF_HELMHOLTZ_OPFACTORY_

#include "type_definitions.H"
#include "mf_helmholtz_op.H"
#include "mfis.H"
#include "mfalias.H"
#include "MFLevelGrid.H"
#include "MFQuadCFInterp.H"

#include <AMRMultiGrid.H>
#include <MFAliasFactory.H>
#include <EBCoarseAverage.H>
#include <BaseDomainBC.H>
#include <BaseEBBC.H>
#include <ProblemDomain.H>
#include <RefCountedPtr.H>
#include <EBIndexSpace.H>
#include <EBQuadCFInterp.H>
#include <MFCellFAB.H>

/*!
  @brief Factory class for generationg mf_helmholtz_op's. 
  @details The jump condition is of the Neumann type in the form
  \f[
  a\partial_n\phi_1 - b\partial_n\phi_2 = c\sigma,
  \f]
  where \f$a\f$, \f$b\f and \f$c\f$ are variables.  
*/
class mf_helmholtz_opfactory : public AMRLevelOpFactory<LevelData<MFCellFAB> > {
public:
  
  /*!
    @brief Factory class for creating multifluid helmholtz operators
    @note This class generates its own EBConductivityOps as well as it's own interpolation utilities. 
    @details This class only supports two-phase solutions with single components
  */
  mf_helmholtz_opfactory(const RefCountedPtr<mfis>&                a_mfis,
			 const Vector<MFLevelGrid>&                a_mflg,
			 const Vector<MFQuadCFInterp>&             a_mfquadcfi,
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
			 const RefCountedPtr<BaseEBBCFactory>&     a_ebbc,
			 const RealVect&                           a_origin,
			 const IntVect&                            a_ghost_phi,
			 const IntVect&                            a_ghost_rhs,
			 int                                       a_num_levels = -1){

    CH_assert(a_mflg[0].num_phases() == 2); // What, you don't like two-phase?
    
    m_mfis       = a_mfis;
    m_mflg       = a_mflg;
    m_mfquadcfi  = a_mfquadcfi;
    m_aco        = a_aco;
    m_bco        = a_bco;
    m_bco_irreg  = a_bco_irreg;
    m_alpha      = a_alpha;
    m_beta       = a_beta;
    m_ebbc       = a_ebbc;
    m_dombc      = a_dombc;
    m_ref_rat    = a_ref_rat;
    m_grids      = a_grids;
    
    m_num_levels = (a_num_levels > 0) ? a_num_levels : m_grids.size();
    
    m_domains.resize(m_num_levels);
    m_dx.resize(m_num_levels);
    m_dx[0] = a_coarsest_dx;
    m_domains[0] = a_coarsest_domain;
    for (int lvl = 1; lvl < m_num_levels; lvl++){
      m_dx[lvl] = m_dx[lvl-1]/m_ref_rat[lvl-1];

      m_domains[lvl] = m_domains[lvl-1];
      m_domains[lvl].refine(m_ref_rat[lvl-1]);
    }

    m_aco_mg.resize(m_num_levels);
    m_bco_mg.resize(m_num_levels);
    m_bco_irreg_mg.resize(m_num_levels);
    m_mflg_mg.resize(m_num_levels);
    m_grids_mg.resize(m_num_levels);
    m_domains_mg.resize(m_num_levels);
    m_has_mg_objects.resize(m_num_levels);
    m_layout_changed.resize(m_num_levels);
    m_layout_changed_mg.resize(m_num_levels);
    
    for (int lvl = 0; lvl < m_num_levels; lvl++){
      if(lvl == 0 || m_ref_rat[lvl] > 2) { // Must be able to generate MultiGrid objects for bottom level and if ref > 2
	
	m_has_mg_objects[lvl] = true;

	const int mg_refi = 2;             // MultiGrid uses VCycling, refinement of 2. 

	m_aco_mg[lvl].resize(0);
	m_bco_mg[lvl].resize(0);
	m_bco_irreg_mg[lvl].resize(0);
	m_mflg_mg[lvl].resize(0);
	m_grids_mg[lvl].resize(0);
	m_domains_mg[lvl].resize(0);
	m_layout_changed_mg.resize(0);

	m_aco_mg[lvl].push_back(m_aco[lvl]);                         // lvl is AMR level, second index is MG level for current
	m_bco_mg[lvl].push_back(m_bco[lvl]);                         // AMR level *sigh*...
	m_bco_irreg_mg[lvl].push_back(m_bco_irreg[lvl]);
	m_mflg_mg[lvl].push_back(m_mflg[lvl]);
	m_grids_mg[lvl].push_back(m_grids[lvl]);
	m_domains_mg[lvl].push_back(m_domains[lvl]);
	m_layout_changed_mg[lvl].push_back(m_layout_changed[lvl]);

	bool has_coarser = true;
	bool at_amr_lvl  = true;
	ProblemDomain cur_domain = m_domains[lvl];
	while(has_coarser){ 

	  int imgsize = m_grids_mg[lvl].size();
	  const DisjointBoxLayout& fine_grid = m_grids_mg[lvl][imgsize - 1]; // Finer grid for current MG level
	  const MFLevelGrid& mflg_fine       = m_mflg_mg[lvl][imgsize - 1];  // Finer MFLevelGrid for current MG level

	  DisjointBoxLayout grid_coar_mg;
	  ProblemDomain domain_coar_mg;

	  bool layout_changed;

	  // Check if we have coarser stuff
	  has_coarser = EBArith::getCoarserLayouts(grid_coar_mg,   // Grid
						   domain_coar_mg, // Domain  
						   fine_grid,      // Fine level grid
						   cur_domain,     // Current domain
						   mg_refi,        // Refinement factor
						   s_max_box_size, // 
						   layout_changed, //
						   s_test_ref);    //

#if 1 // Debug
	  pout() << "lvl = " << lvl << "\t has_coarser = " << has_coarser << endl;
#endif

	  if(at_amr_lvl && !has_coarser){
	    m_has_mg_objects[lvl] = false;
	  }

	  if(at_amr_lvl){
	    m_layout_changed[lvl] = layout_changed;
	    at_amr_lvl = false;
	  }

	  if(has_coarser){
	    m_grids_mg[lvl].push_back(grid_coar_mg);
	    m_domains_mg[lvl].push_back(domain_coar_mg);
	    m_layout_changed_mg[lvl].push_back(layout_changed);
	    cur_domain.coarsen(mg_refi);

	    const int  ncomps = 1; // Number of components we solve for. Always 1. 
	    const int ebghost = 4; // Ghost cells for MG, using 4 since that allows refinement of 4
	    const int   ghost = 1; // Necessary ghost cells for second order

	    m_mflg_mg[lvl].push_back(MFLevelGrid(grid_coar_mg, domain_coar_mg, ebghost, mflg_fine.get_ebis()));

	    const int img = m_mflg_mg[lvl].size() - 1; // Last one added, i.e. the coarsest that we have so far

	    const MFLevelGrid& mflg_coar = m_mflg_mg[lvl][img  ];    // Coarsened EBLevelGrids for all phases
	    const MFLevelGrid& mflg_fine = m_mflg_mg[lvl][img-1];    // Fine EBLevelGrids for all phases

	    Vector<EBISLayout> ebisl_coar = mflg_coar.get_ebisl();  // Need this stuff for factories. 
	    Vector<int> comps;
	    for (int i = 0; i < ebisl_coar.size(); i++){
	      comps.push_back(ncomps);
	    }

	    MFCellFactory      cellfact(ebisl_coar, comps);
	    MFFluxFactory      fluxfact(ebisl_coar, comps);
	    MFBaseIVFABFactory ivfact(ebisl_coar,   comps);

	    RefCountedPtr<LevelData<MFCellFAB> > aco_coar = RefCountedPtr<LevelData<MFCellFAB> >
	      (new LevelData<MFCellFAB>(grid_coar_mg, ncomps, ghost*IntVect::Unit, cellfact));
	    RefCountedPtr<LevelData<MFFluxFAB> > bco_coar = RefCountedPtr<LevelData<MFFluxFAB> >
	      (new LevelData<MFFluxFAB>(grid_coar_mg, ncomps, ghost*IntVect::Unit, fluxfact));
	    RefCountedPtr<LevelData<MFBaseIVFAB> > bco_irreg_coar = RefCountedPtr<LevelData<MFBaseIVFAB> >
	      (new LevelData<MFBaseIVFAB>(grid_coar_mg, ncomps, ghost*IntVect::Unit, ivfact));
	    

	    this->coarsen_coefficients(*aco_coar, // Coarsen coefficients
	    			       *bco_coar,
	    			       *bco_irreg_coar,
	    			       mflg_coar,
	    			       mflg_fine,
	    			       *m_aco_mg[lvl][img-1],
	    			       *m_bco_mg[lvl][img-1],
	    			       *m_bco_irreg_mg[lvl][img-1],
	    			       mg_refi);

	    m_aco_mg[lvl].push_back(aco_coar);
	    m_bco_mg[lvl].push_back(bco_coar);
	    m_bco_irreg_mg[lvl].push_back(bco_irreg_coar);


	  }
	}
      }
      else {
	m_has_mg_objects[lvl] = false;
      }


    }

    m_defined    = true;
    m_jumpset    = false;

#if 1 // Debugging hook
    for (int lvl = 0; lvl < m_num_levels; lvl++){
      pout() << "lvl = "  << lvl
	     << "\t has MG objects = " << m_has_mg_objects[lvl]
	     << "\t MG levels = " << m_grids_mg[lvl].size()
	     << endl;
    }
#endif
  }

  /*!
    @brief Destructor
  */
  ~mf_helmholtz_opfactory();

  /*!
    @brief Set when to drop to the bottom solver
  */
  static void set_bottom_drop(const int a_bottom_drop);

  /*!
    @brief Set when to drop to the bottom solver
  */
  static void set_relax_type(const int a_relax_type);

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
				    const int&                    a_ref_to_depth){
    CH_assert(a_ref_to_depth > 0);
    
    const int ncomp = 1;
    const Interval interv(0,ncomp - 1);

    if(a_ref_to_depth == 1){
      a_aco_fine.copyTo(interv,       a_aco_coar,       interv);
      a_bco_fine.copyTo(interv,       a_bco_coar,       interv);
      a_bco_irreg_fine.copyTo(interv, a_bco_irreg_coar, interv);
    }
    else {
      for (int i = 0; i < a_mflg_coar.num_phases(); i++){
	const EBLevelGrid& eblg_coar = a_mflg_coar.get_eblg(i);
	const EBLevelGrid& eblg_fine = a_mflg_fine.get_eblg(i);
	EBCoarseAverage aveop(eblg_fine.getDBL(),    eblg_coar.getDBL(),
			      eblg_fine.getEBISL(),  eblg_coar.getEBISL(),
			      eblg_coar.getDomain(), a_ref_to_depth, ncomp,
			      eblg_coar.getEBIS());

	LevelData<EBCellFAB>        aco_coar;
	LevelData<EBCellFAB>        aco_fine;
	LevelData<EBFluxFAB>        bco_coar;
	LevelData<EBFluxFAB>        bco_fine;
	LevelData<BaseIVFAB<Real> > bco_irreg_coar;
	LevelData<BaseIVFAB<Real> > bco_irreg_fine;

	mfalias::aliasMF(aco_coar,       i, a_aco_coar);
	mfalias::aliasMF(aco_fine,       i, a_aco_fine);
	mfalias::aliasMF(bco_coar,       i, a_bco_coar); 
	mfalias::aliasMF(bco_fine,       i, a_bco_fine);       
	mfalias::aliasMF(bco_irreg_coar, i, a_bco_irreg_coar); 
	mfalias::aliasMF(bco_irreg_fine, i, a_bco_irreg_fine);
	
	aveop.average(aco_coar,       aco_fine,       interv);
	aveop.average(bco_coar,       bco_fine,       interv);
	aveop.average(bco_irreg_coar, bco_irreg_fine, interv);

	aco_coar.exchange();
	bco_coar.exchange();
	bco_irreg_coar.exchange();
      }
    }
  }

  /*!
    @brief Coarsen jump conditions; 
    @details Routine is called every time set_jump is called because coarsened BCs are needed for multigrid. 
  */
  virtual void coarsen_stuff(LevelData<BaseIVFAB<Real> >&       a_ajump_coar,
			     LevelData<BaseIVFAB<Real> >&       a_bjump_coar,
			     LevelData<BaseIVFAB<Real> >&       a_cjump_coar,
			     const EBLevelGrid&                 a_eblg_coar,
			     const EBLevelGrid&                 a_eblg_fine,
			     const LevelData<BaseIVFAB<Real> >& a_ajump_fine,
			     const LevelData<BaseIVFAB<Real> >& a_bjump_fine,
			     const LevelData<BaseIVFAB<Real> >& a_cjump_fine,
			     const int a_ref_to_depth){
    CH_assert(a_ref_to_depth > 0);

    const Interval interv(0,0);
    if(a_ref_to_depth == 1){
      a_ajump_fine.copyTo(interv, a_ajump_coar, interv);
      a_bjump_fine.copyTo(interv, a_bjump_coar, interv);
      a_cjump_fine.copyTo(interv, a_cjump_coar, interv);
    }
    else{
      MayDay::Abort("mf_helmholtz_opfactory::coarsen_stuff - EBCoarseAverage is not defined");
      m_aveop[0].average(a_ajump_coar, a_ajump_fine, interv);
      m_aveop[0].average(a_bjump_coar, a_bjump_fine, interv);
      m_aveop[0].average(a_cjump_coar, a_cjump_fine, interv);
    }

    a_ajump_coar.exchange(interv);
    a_bjump_coar.exchange(interv);
    a_cjump_coar.exchange(interv);
  }

  /*!
    @brief Set jump conditions. 
  */
  virtual void set_jump(const EBAMRIVData& a_ajump,
			const EBAMRIVData& a_bjump,
			const EBAMRIVData& a_cjump,
			const EBAMRIVData& a_sigma);

  /*!
    @brief Factory method
  */
  virtual MGLevelOp<LevelData<MFCellFAB> >* MGnewOp(const ProblemDomain& a_domain_fine, int a_depth, bool a_homo_only = true){

    // Find out the starting point. 
    int ref_lvl;
    bool found = false;

    for (int lvl = 0; lvl < m_num_levels && !found; lvl++){
      if(a_domain_fine == m_domains[lvl]){
	found = true;
	ref_lvl = lvl;
      }
    }

    if(!found){
      MayDay::Error("mf_helmholtzopfactory::MGnewOp - no corresponding AMRLevel to starting point of MGnewOp");
    }

    Real           dx_mg_level;
    MFLevelGrid    mflg_level_mg;
    MFLevelGrid    mflg_coar_mg;
    MFQuadCFInterp quadcfi;        // Only defined on amr levels

    RefCountedPtr<LevelData<MFCellFAB> >   aco;
    RefCountedPtr<LevelData<MFFluxFAB> >   bco;
    RefCountedPtr<LevelData<MFBaseIVFAB> > bco_irreg;
    
    if(a_depth == 0){ // Depth = 0 implies an AMR level
      mflg_level_mg = m_mflg[ref_lvl];

      aco       = m_aco[ref_lvl];
      bco       = m_bco[ref_lvl];
      bco_irreg = m_bco_irreg[ref_lvl];
      quadcfi   = m_mfquadcfi[ref_lvl];

      if(m_has_mg_objects[ref_lvl]){ 
	mflg_coar_mg = m_mflg_mg[ref_lvl][1]; 
      }
    }
    else{
      int icoar = 1;
      for (int idep = 0; idep < a_depth; idep++){
	icoar *= 2;
      }

      const ProblemDomain domain_fine = mflg_level_mg.get_eblg(0).getDomain();
      ProblemDomain domain_box_mg_lvl = coarsen(domain_fine, icoar);

      bool found_mg_level = false;
      int num_mg_levels   = m_mflg_mg[ref_lvl].size();

      for (int img = 0; img < num_mg_levels; img++){
	if(m_mflg_mg[ref_lvl][img].get_eblg(0).getDomain() == domain_box_mg_lvl){
	  
	}
      }
    }


    
    return static_cast<MGLevelOp<LevelData<MFCellFAB> >* > (NULL);
  }

  /*!
    @brief Factory method
  */
  virtual AMRLevelOp<LevelData<MFCellFAB> >* AMRnewOp(const ProblemDomain& a_fineindexspace);

  /*!
    @brief Factory method
  */
  virtual mf_helmholtz_op* createOperator(const DisjointBoxLayout&       a_dilboMGLevel,
					  const DisjointBoxLayout&       a_dilboCoarMG,
					  const ProblemDomain&           a_domainMGLevel,
					  const bool&                    a_hasMGObjects,
					  const bool&                    a_layoutChanged,
					  const RealVect&                a_dxMGLevel,
					  const RealVect&                a_dxCoar,
					  const int&                     a_whichLevel,
					  const int&                     a_mgLevel){

    return static_cast<mf_helmholtz_op*> (NULL);
    
  }

  /*!
    @brief Reclaim
  */
  virtual void reclaim(MGLevelOp<LevelData<EBCellFAB> >* a_reclaim);

  /*!
    @brief Reclaim
  */
  virtual void AMRreclaim(mf_helmholtz_op* a_reclaim);

  /*!
    @brief Refinement ratio between this level and the coarser level. Returns 1 when there are no coarser levels
  */
  virtual int refToFiner(const ProblemDomain& a_domain) const;

protected:

  /*!
    @brief Drop to bottom solver
  */
  static int s_test_ref;

  /*!
    @brief Internal stuff
  */
  static int s_max_box_size;

  /*!
    @brief Relaxation type
  */
  static int s_relax_type;

  /*!
    @brief EBIndexSpaces
  */
  RefCountedPtr<mfis> m_mfis;

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
    @brief Jump coefficient
  */
  EBAMRIVData m_ajump;
  
  /*!
    @brief Jump coefficient
  */
  EBAMRIVData m_bjump;

  /*!
    @brief Jump coefficient
  */
  EBAMRIVData m_cjump;

  /*!
    @brief Jump condition
  */
  EBAMRIVData m_sigma;

  /*!
    @brief EBLevelGrids (in MFLevelGrid form)
  */
  Vector<MFLevelGrid> m_mflg;

  /*!
    @brief Interpolation utilities
  */
  Vector<MFQuadCFInterp> m_mfquadcfi;

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
    @brief Averaging operator - not sure how to use this yet :)
  */
  Vector<EBCoarseAverage> m_aveop;

  /*!
    @brief Grids for multigrid
  */
  Vector<Vector<DisjointBoxLayout> > m_grids_mg;

  /*!
    @brief Domains for multigrid
  */
  Vector<Vector<ProblemDomain> > m_domains_mg;

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
    @brief EB boundary conditions
  */
  RefCountedPtr<BaseEBBCFactory> m_ebbc;

  /*!
    @brief Number of AMR levels
  */
  int m_num_levels;

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
    @brief Defined or not
  */
  bool m_defined;

  /*!
    @brief Jump condition set or not
  */
  bool m_jumpset;
  
private:
  
  /*!
    @brief Disallowed weak constructor
  */
  mf_helmholtz_opfactory(){
    MayDay::Error("mf_helmholtz_opfactory:: invalid constructor");
  }
};

#endif
