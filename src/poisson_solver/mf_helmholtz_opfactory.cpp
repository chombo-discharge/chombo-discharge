/*!
  @file mf_helmholtz_opfactory.cpp
  @brief Implementation of mf_helmholtz_opfactory.H
  @author Robert Marskar
*/

#include "mf_helmholtz_opfactory.H"

int mf_helmholtz_opfactory::s_max_box_size = 32;
int mf_helmholtz_opfactory::s_test_ref     = 4;
int mf_helmholtz_opfactory::s_relax_type   = 2;

mf_helmholtz_opfactory::mf_helmholtz_opfactory(const RefCountedPtr<mfis>&                a_mfis,
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
					       int                                       a_num_levels){

  CH_assert(a_mflg[0].num_phases() == 2); // What, you don't like two-phase?

  m_num_levels = (a_num_levels > 0) ? a_num_levels : a_grids.size();
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
    
  m_domains.resize(m_num_levels);
  m_dx.resize(m_num_levels);
    
  m_dx[0]      = a_coarsest_dx;
  m_domains[0] = a_coarsest_domain;
    
  for (int lvl = 1; lvl < m_num_levels; lvl++){
    m_dx[lvl] = m_dx[lvl-1]/m_ref_rat[lvl-1];
    m_domains[lvl] = m_domains[lvl-1];
    m_domains[lvl].refine(m_ref_rat[lvl-1]);
  } 
    
  this->define_jump_stuff();      // Define jump cell stuff
  this->define_multigrid_stuff(); // Define things for lower levels of multigrid. Must happen after define_jump_stuff
  this->set_jump(0.0, 1.0);       // Default, no surface charge. 


#if 0 // Debugging hook
  for (int lvl = 0; lvl < m_num_levels; lvl++){
    pout() << "lvl = "  << lvl
	   << "\t has MG objects = " << m_has_mg_objects[lvl]
	   << "\t MG levels = " << m_grids_mg[lvl].size()
	   << endl;
  }
#endif

}

mf_helmholtz_opfactory::~mf_helmholtz_opfactory(){

}

void mf_helmholtz_opfactory::define_multigrid_stuff(){
  CH_TIME("mf_helmholtz_opfactory::define_multigrid_stuff");
  m_aco_mg.resize(m_num_levels);
  m_bco_mg.resize(m_num_levels);
  m_bco_irreg_mg.resize(m_num_levels);
  m_mflg_mg.resize(m_num_levels);
  m_grids_mg.resize(m_num_levels);
  m_aveop_mg.resize(m_num_levels);
  m_domains_mg.resize(m_num_levels);
  m_has_mg_objects.resize(m_num_levels);
  m_layout_changed.resize(m_num_levels);
  m_layout_changed_mg.resize(m_num_levels);
  m_aveop_mg.resize(m_num_levels);
  m_jump_mg.resize(m_num_levels);
    
  for (int lvl = 0; lvl < m_num_levels; lvl++){
    if(lvl == 0 || m_ref_rat[lvl] > 2) { // Must be able to generate MultiGrid objects for bottom level and if ref > 2
	
      m_has_mg_objects[lvl] = true;

      const int mg_refi = 2;             // MultiGrid uses VCycling, refinement of 2. 

      m_aco_mg[lvl].resize(0);            // aco for all MG levels at this level
      m_bco_mg[lvl].resize(0);            // bco for all MG levels at this level
      m_bco_irreg_mg[lvl].resize(0);      // bco for all MG levels at this level
      m_mflg_mg[lvl].resize(0);           // MFLevelGrids for all MG levels at this level
      m_grids_mg[lvl].resize(0);          // Grids for all MG levels at this level
      m_domains_mg[lvl].resize(0);        // Domains for all MG levels at this level
      m_layout_changed_mg[lvl].resize(0); // Layout changed for all MG levels at this level
      m_aveop_mg[lvl].resize(0);          //
      m_jump_mg[lvl].resize(0);          // sigma for all MG levels at this AMR level

      m_aco_mg[lvl].push_back(m_aco[lvl]);                         // MG depth 0 is an AMR level
      m_bco_mg[lvl].push_back(m_bco[lvl]);                         //  
      m_bco_irreg_mg[lvl].push_back(m_bco_irreg[lvl]);             // 
      m_mflg_mg[lvl].push_back(m_mflg[lvl]);                       // 
      m_grids_mg[lvl].push_back(m_grids[lvl]);                     // 
      m_domains_mg[lvl].push_back(m_domains[lvl]);                 // 
      m_layout_changed_mg[lvl].push_back(m_layout_changed[lvl]);   //
      m_aveop_mg[lvl].push_back(m_aveop[lvl]);                     // This is null for lvl = 0
      m_jump_mg[lvl].push_back(m_jump[lvl]);                       //

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

	  //	    Vector<EBISLayout> ebisl_coar = mflg_coar.get_ebisl();  // Need this stuff for factories. 
	  Vector<int> comps;
	  Vector<EBISLayout> ebisl_coar;
	  for (int i = 0; i < mflg_coar.num_phases(); i++){
	    const EBLevelGrid& eblg = mflg_coar.get_eblg(i);
	    const EBISLayout& ebisl = eblg.getEBISL();
	    ebisl_coar.push_back(ebisl);
	    comps.push_back(ncomps);
	  }

	  // Averaging operator for m_jump_mg
	  const int ncomp      = 1;
	  const int main_phase = 0;
	  const EBLevelGrid& eblg_fine = mflg_fine.get_eblg(main_phase);
	  const EBLevelGrid& eblg_coar = mflg_coar.get_eblg(main_phase);
	  RefCountedPtr<EBCoarseAverage> aveop(new EBCoarseAverage(eblg_fine.getDBL(),    eblg_coar.getDBL(),
								   eblg_fine.getEBISL(),  eblg_coar.getEBISL(),
								   eblg_coar.getDomain(), mg_refi, ncomp,
								   eblg_coar.getEBIS()));

	  // Coarsened m_jump_mg
	  LayoutData<IntVectSet> isect_cells (eblg_coar.getDBL());
	  for (DataIterator dit = isect_cells.dataIterator(); dit.ok(); ++dit){
	    isect_cells[dit()] = m_mfis->interface_region(eblg_coar.getDomain()) & eblg_coar.getDBL().get(dit());
	  }
	  BaseIVFactory<Real> fact(eblg_coar.getEBISL(), isect_cells);

	  RefCountedPtr<LevelData<BaseIVFAB<Real> > > jump_coar = RefCountedPtr<LevelData<BaseIVFAB<Real> > >
	    (new LevelData<BaseIVFAB<Real> > (grid_coar_mg, ncomps, ghost*IntVect::Unit, fact)); 

	  


	  // Coarsened coefficients
	  MFCellFactory      cellfact(ebisl_coar, comps);
	  MFFluxFactory      fluxfact(ebisl_coar, comps);
	  MFBaseIVFABFactory ivfact  (ebisl_coar, comps);

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
	  m_aveop_mg[lvl].push_back(aveop);
	  m_jump_mg[lvl].push_back(jump_coar);
	}
      }
    }
    else {
      m_has_mg_objects[lvl] = false;
    }


  }

}

void mf_helmholtz_opfactory::coarsen_coefficients(LevelData<MFCellFAB>&         a_aco_coar,
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

void mf_helmholtz_opfactory::set_bottom_drop(const int a_bottom_drop){
  s_test_ref = a_bottom_drop;
}

void mf_helmholtz_opfactory::set_relax_type(int a_relax_type){
  s_relax_type = a_relax_type;
}

void mf_helmholtz_opfactory::reclaim(MGLevelOp<LevelData<EBCellFAB> >* a_reclaim){
  delete a_reclaim;
}

void mf_helmholtz_opfactory::AMRreclaim(mf_helmholtz_op* a_reclaim){
  delete a_reclaim;
}

void mf_helmholtz_opfactory::define_jump_stuff(){
  CH_TIME("mf_helmholtz_opfactory::define_jump_cells");

  m_aveop.resize(m_num_levels);
  m_jump.resize(m_num_levels);
  m_jumpcells.resize(m_num_levels);

  const int ghost      = 1; // Using one ghost cell for this
  const int ncomp      = 1; // Only single-component stuff
  const int main_phase = 0; // Interface region is the intersection between gas-side irregular cells and solid phase cells
    
  for (int lvl = 0; lvl < m_num_levels; lvl++){
    const EBLevelGrid& eblg  = m_mflg[lvl].get_eblg(main_phase);
    const EBISLayout& ebisl = eblg.getEBISL();
    const IntVectSet interface_cells = m_mfis->interface_region(m_domains[lvl]);

    m_jumpcells[lvl] = RefCountedPtr<LayoutData<IntVectSet> > (new LayoutData<IntVectSet> (m_grids[lvl]));
    for (DataIterator dit = m_jumpcells[lvl]->dataIterator(); dit.ok(); ++dit){
      (*m_jumpcells[lvl])[dit()] = interface_cells & m_grids[lvl].get(dit());
    }

    BaseIVFactory<Real> fact(ebisl, *m_jumpcells[lvl]);
    m_jump[lvl] = RefCountedPtr<LevelData<BaseIVFAB<Real> > >
      (new LevelData<BaseIVFAB<Real> >(m_grids[lvl], ncomp, ghost*IntVect::Unit, fact));
  }

  for (int lvl = 0; lvl < m_num_levels; lvl++){

    const bool has_coar = lvl > 0;

    if(has_coar){
      const EBLevelGrid& eblg_fine = m_mflg[lvl].get_eblg(main_phase);
      const EBLevelGrid& eblg_coar = m_mflg[lvl-1].get_eblg(main_phase);
      const int ref_ratio          = m_ref_rat[lvl];
	
      m_aveop[lvl] = RefCountedPtr<EBCoarseAverage> (new EBCoarseAverage(eblg_fine.getDBL(),    eblg_coar.getDBL(),
									 eblg_fine.getEBISL(),  eblg_coar.getEBISL(),
									 eblg_coar.getDomain(), ref_ratio, ncomp,
									 eblg_coar.getEBIS()));
    }
  }
}

void mf_helmholtz_opfactory::average_down_amr(){
  CH_TIME("mf_helmholtz_opfactory::average_down_amr");
    
  const int ncomp        = 0;
  const Interval interv  = Interval(0, ncomp -1);
  const int finest_level = m_num_levels - 1;

  for (int lvl = finest_level; lvl > 0; lvl--){ // Average down AMR levels
    m_aveop[lvl]->average(*m_jump[lvl-1], *m_jump[lvl], interv);
#if 1 // Debug
    pout() << "mf_helmholtz_opfactory::average_down_amr from AMR level = " << lvl << " and onto AMR level = " << lvl - 1 << endl;
#endif
  }
}

void mf_helmholtz_opfactory::average_down_mg(){
  CH_TIME("mf_helmholtz_opfactory::average_down_mg");

  const int ncomp        = 0;
  const Interval interv  = Interval(0, ncomp -1);
  
  for (int lvl = 0; lvl < m_num_levels; lvl++){ // Average down the MG stuff
    if(m_has_mg_objects[lvl]){
      EBAMRIVData& jump_mg = m_jump_mg[lvl];

      const int finest_mg_level = jump_mg.size() - 1;

      // img = 0 is the finest MG level, but this has already been averaged down
      for (int img = 1; img <= finest_mg_level; img++){ // This time, 0 is the finest level. m_jump_mg[lvl][0] is the same
#if 1 // DEBUG
	pout() << "mf_helmholtz_opfactory::average_down_mg from AMR level = " << lvl
	       << " from MG level = " << img+1
	       << " to   MG level = " << img << endl;
#endif
	m_aveop_mg[lvl][img]->average(*jump_mg[img], *jump_mg[img+1], interv);
      }
    }
  }
  MayDay::Abort("mf_helmholtz_opfactory::average_down_mg - ensure that everything is averaged down correctly!");
}

void mf_helmholtz_opfactory::set_jump(const Real& a_sigma, const Real& a_scale){
  CH_TIME("mf_helmholtz_opfactory::set_jump(scalar)");
  for (int lvl = 0; lvl < m_num_levels; lvl++){
    EBLevelDataOps::setVal(*m_jump[lvl], a_sigma);
    data_ops::scale(*m_jump[lvl], a_sigma);
  }

  this->average_down_amr();
  this->average_down_mg();
}

void mf_helmholtz_opfactory::set_jump(const EBAMRIVData& a_sigma, const Real& a_scale){
  CH_TIME("mf_helmholtz_opfactory::set_jump(data based)");

  for (int lvl = 0; lvl < m_num_levels; lvl++){
    a_sigma[lvl]->copyTo(*m_jump[lvl]);
    data_ops::scale(*m_jump[lvl], a_scale);
  }

  this->average_down_amr();
  this->average_down_mg();
}

int mf_helmholtz_opfactory::refToFiner(const ProblemDomain& a_domain) const{
  int retval = -1;
  bool found = false;

  for (int lvl = 0; lvl < m_domains.size(); lvl++){
    if(m_domains[lvl] == a_domain){
      retval = m_ref_rat[lvl];
      found  = true;
    }
  }
  
  if(!found){
    MayDay::Error("mf_helmholtz_opfactory::refToFiner - domain not found in AMR hierarchy");
  }
  
  return retval;
}

MGLevelOp<LevelData<MFCellFAB> >* mf_helmholtz_opfactory::MGnewOp(const ProblemDomain& a_domain_fine,
								  int                  a_depth,
								  bool                 a_homo_only){
  CH_TIME("mf_helmholtz_opfactory::MGnewOp");
    
  // Find the AMR level starting point
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

  bool           has_coarser_mg;
  Real           dx_mg_level;
  MFLevelGrid    mflg_level_mg;
  MFLevelGrid    mflg_coar_mg;
  MFQuadCFInterp quadcfi;        // Only defined on amr levels

  RefCountedPtr<LevelData<MFCellFAB> >   aco;
  RefCountedPtr<LevelData<MFFluxFAB> >   bco;
  RefCountedPtr<LevelData<MFBaseIVFAB> > bco_irreg;
    
  if(a_depth == 0){ // This is an AMR level
    mflg_level_mg = m_mflg[ref_lvl];       // Levelgrids    }
    aco           = m_aco[ref_lvl];        // coefficients  }
    bco           = m_bco[ref_lvl];        // coefficient   } Non-coarsened stuff
    bco_irreg     = m_bco_irreg[ref_lvl];  // coefficient   }
    quadcfi       = m_mfquadcfi[ref_lvl];  // QuadCFI       }

    has_coarser_mg = m_has_mg_objects[ref_lvl];
    if(has_coarser_mg){ 
      mflg_coar_mg = m_mflg_mg[ref_lvl][1]; // Coarsened EBLevelGrids
    }
  }
  else{ // Not an AMR level
    int icoar = 1;
    for (int idep = 0; idep < a_depth; idep++){
      icoar *= 2;
    }

    const ProblemDomain domain_fine = mflg_level_mg.get_eblg(0).getDomain();   // 
    ProblemDomain domain_box_mg_lvl = coarsen(domain_fine, icoar);             // Coarsened a_depth times by factor 2

    bool found_mg_level = false;
    int num_mg_levels   = m_mflg_mg[ref_lvl].size();

    for (int img = 0; img < num_mg_levels; img++){
      if(m_mflg_mg[ref_lvl][img].get_eblg(0).getDomain() == domain_box_mg_lvl){ // Check which mg level we're after

	aco            = m_aco_mg[ref_lvl][img];
	bco            = m_bco_mg[ref_lvl][img];
	bco_irreg      = m_bco_irreg_mg[ref_lvl][img];
	found_mg_level = true;

	MayDay::Warning("mf_helmholtz_opfactory::MGnewOp - boundary conditions have not been placed yet");

	has_coarser_mg = img+1 < m_mflg_mg[ref_lvl].size(); // Are there more mg levels below this one?
	if(has_coarser_mg){
	  mflg_coar_mg = m_mflg_mg[ref_lvl][img+1]; // Next coarser MFLevelGrid for MG
	}
	break;
      }
    }
      
    bool coarsenable = found_mg_level;
    dx_mg_level = m_dx[ref_lvl];  // Resolution on AMR level
    dx_mg_level *= Real(icoar);   // Resolution on MG level

    if(!coarsenable) {
      return NULL;
    }
  }

  // Create boundary conditions - should the jump conditions be a part of the BaseEBBC or not?
    
  return static_cast<MGLevelOp<LevelData<MFCellFAB> >* > (NULL);
}

AMRLevelOp<LevelData<MFCellFAB> >* mf_helmholtz_opfactory::AMRnewOp(const ProblemDomain& a_fineindexspace){
  MayDay::Abort("mf_helmholtz_opfactory::AMRnewOp - not implemented");
  return static_cast<AMRLevelOp<LevelData<MFCellFAB> >* > (NULL);
}

mf_helmholtz_op* mf_helmholtz_opfactory::createOperator(const DisjointBoxLayout&       a_dilboMGLevel,
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
