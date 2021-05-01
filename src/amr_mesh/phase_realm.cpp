/*!
  @file   phase_realm.cpp
  @brief  Implementation of phase_realm.H
  @author Robert Marskar
  @date   July 2020
*/

#include "phase_realm.H"
#include "EBFastFineToCoarRedist.H"
#include "EBFastCoarToFineRedist.H"
#include "EBFastCoarToCoarRedist.H"
#include "load_balance.H"
#include "EBFastFluxRegister.H"

#include <EBArith.H>


phase_realm::phase_realm(){
  m_defined   = false;
  m_verbosity = -1;

  // Always do this shit. 
  this->register_operator(s_eb_gradient);
  this->register_operator(s_eb_irreg_interp);
}

phase_realm::~phase_realm(){

}

void phase_realm::define(const Vector<DisjointBoxLayout>& a_grids,
			 const Vector<ProblemDomain>& a_domains,
			 const Vector<int>& a_ref_rat,
			 const Vector<Real>& a_dx,
			 const RealVect a_prob_lo,
			 const int a_finest_level,
			 const int a_ebghost,
			 const int a_num_ghost,
			 const int a_lsf_ghost,
			 const int a_redist_rad,
			 const stencil_type a_centroid_stencil,
			 const stencil_type a_eb_stencil,
			 const bool a_ebcf,
			 const RefCountedPtr<BaseIF>&       a_baseif,
			 const RefCountedPtr<EBIndexSpace>& a_ebis){

  m_ebis = a_ebis;
  m_finest_level = a_finest_level;
  m_grids = a_grids;
  m_domains = a_domains;
  m_ref_ratios = a_ref_rat;
  m_dx = a_dx;
  m_ebcf = a_ebcf;
  m_ebghost = a_ebghost;
  m_num_ghost = a_num_ghost;
  m_lsf_ghost = a_lsf_ghost;
  m_redist_rad = a_redist_rad;
  m_centroid_stencil = a_centroid_stencil;
  m_eb_stencil = a_eb_stencil;
  m_baseif = a_baseif;
  m_prob_lo = a_prob_lo;
  
  if(!m_ebis.isNull()){
    m_defined = true;
  }
}

void phase_realm::set_grids(const Vector<DisjointBoxLayout>& a_grids, const int a_finest_level){
  CH_TIME("phase_realm::set_grids");
  if(m_verbosity > 5){
    pout() << "phase_realm::set_grids" << endl;
  }

  if(m_defined){
    m_grids = a_grids;
    m_finest_level = a_finest_level;
  }
}

void phase_realm::regrid_base(const int a_lmin){
  CH_TIME("phase_realm::regrid_base");
  if(m_verbosity > 5){
    pout() << "phase_realm::regrid_base" << endl;
  }

  if(m_defined){
    this->define_eblevelgrid(a_lmin);
    this->define_neighbors(a_lmin);
    this->define_vofiter(a_lmin);
  }
}

void phase_realm::regrid_operators(const int a_lmin, const int a_lmax, const int a_regsize){
  CH_TIME("phase_realm::regrid_operators_phase");
  if(m_verbosity > 5){
    pout() << "phase_realm::regrid_operators" << endl;
  }

  if(m_defined){
    this->define_eb_coar_ave(a_lmin);             // Define ebcoarseaverage on both phases
    this->define_eb_quad_cfi(a_lmin);             // Define nwoebquadcfinterp on both phases.
    this->define_fillpatch(a_lmin);               // Define operator for piecewise linear interpolation of ghost cells
    this->define_ebpwl_interp(a_lmin);            // Define interpolator for piecewise interpolation of interior points
    this->define_ebmg_interp(a_lmin);             // Define interpolator used for e.g. multigrid (or piecewise constant)
    this->define_flux_reg(a_lmin,a_regsize);      // Define flux register
    this->define_redist_oper(a_lmin, a_regsize);  // Define redistribution (phase::gas only)
    this->define_gradsten(a_lmin);                // Make stencils for computing gradients
    this->define_irreg_sten();                    // Make stencils for doing interpolation to centroids
    this->define_noncons_sten();                  // Make stencils for nonconservative averaging
    this->define_copier(a_lmin);                  // Make stencils for copier
    this->define_ghostcloud(a_lmin);              // Make stencils for ghost clouds with particle depositions
    this->define_levelset(a_lmin, m_lsf_ghost);   // Defining levelset
  }
}

void phase_realm::register_operator(const std::string a_operator){
  CH_TIME("phase_realm::register_operator");
  if(m_verbosity > 5){
    pout() << "phase_realm::register_operator" << endl;
  } 

  // These are the supported operators - issue an error if we ask for something that is not supported. 
  if(!(a_operator.compare(s_eb_coar_ave)     == 0 ||
       a_operator.compare(s_eb_quad_cfi)     == 0 ||
       a_operator.compare(s_eb_fill_patch)   == 0 ||
       a_operator.compare(s_eb_pwl_interp)   == 0 ||
       a_operator.compare(s_eb_flux_reg)     == 0 ||
       a_operator.compare(s_eb_redist)       == 0 ||
       a_operator.compare(s_eb_noncons_div)  == 0 ||
       a_operator.compare(s_eb_copier)       == 0 ||
       a_operator.compare(s_eb_ghostcloud)   == 0 ||
       a_operator.compare(s_eb_gradient)     == 0 ||
       a_operator.compare(s_eb_irreg_interp) == 0 ||
       a_operator.compare(s_eb_mg_interp)    == 0 ||
       a_operator.compare(s_levelset)        == 0 )){

    const std::string str = "phase_realm::register_operator - unknown operator '" + a_operator + "' requested";
    MayDay::Abort(str.c_str());
  }

  if(!this->query_operator(a_operator)){
    m_operator_map.emplace(a_operator, true);
  }
}

bool phase_realm::query_operator(const std::string a_operator) {
  CH_TIME("phase_realm::query_operator");
  if(m_verbosity > 5){
    pout() << "phase_realm::query_operator" << endl;
  }

  //  return m_operator_map[a_operator];

  bool ret = false;
  if(m_defined){
    ret = true;
    
    if(m_operator_map.find(a_operator) == m_operator_map.end()){
      ret = false;
    }
  }

  return ret;
}
  

void phase_realm::define_eblevelgrid(const int a_lmin){
  CH_TIME("phase_realm::define_eblevelgrid");
  if(m_verbosity > 2){
    pout() << "phase_realm::define_eblevelgrid" << endl;
  }

  m_eblg.resize(1 + m_finest_level);
  m_ebisl.resize(1 + m_finest_level);

  for (int lvl = a_lmin; lvl <= m_finest_level; lvl++){
    m_eblg[lvl]  = RefCountedPtr<EBLevelGrid> (new EBLevelGrid(m_grids[lvl], m_domains[lvl], m_ebghost, &(*m_ebis)));
    
    if(lvl > 0) m_eblg[lvl]->setMaxCoarseningRatio(m_ref_ratios[lvl-1], &(*m_ebis));

    m_ebisl[lvl] = m_eblg[lvl]->getEBISL();
  }
}
void phase_realm::define_vofiter(const int a_lmin){
  CH_TIME("phase_realm::define_vofiter");
  if(m_verbosity > 2){
    pout() << "phase_realm::define_vofiter" << endl;
  }

  m_vofiter.resize(1 + m_finest_level);

  for (int lvl = a_lmin; lvl <= m_finest_level; lvl++){
      
    m_vofiter[lvl] = RefCountedPtr<LayoutData<VoFIterator> > (new LayoutData<VoFIterator> (m_grids[lvl]));
    
    for (DataIterator dit = m_grids[lvl].dataIterator(); dit.ok(); ++dit){
      VoFIterator& vofit = (*m_vofiter[lvl])[dit()];

      const Box& box         = m_grids[lvl].get(dit());
      const EBISBox& ebisbox = m_ebisl[lvl][dit()];
      const EBGraph& ebgraph = ebisbox.getEBGraph();
      const IntVectSet& irreg = ebisbox.getIrregIVS(box);

      vofit.define(irreg, ebgraph);
    }
  }
}

void phase_realm::define_neighbors(const int a_lmin){
  CH_TIME("phase_realm::define_neighbors");
  if(m_verbosity > 2){
    pout() << "phase_realm::define_neighbors" << endl;
  }

  m_neighbors.resize(1 + m_finest_level);

  for (int lvl = a_lmin; lvl <= m_finest_level; lvl++){
    m_neighbors[lvl] = RefCountedPtr<LayoutData<Vector<LayoutIndex> > > (new LayoutData<Vector<LayoutIndex> >(m_grids[lvl]));

    const DisjointBoxLayout& dbl = m_grids[lvl];
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const Box& box = dbl.get(dit());

      Vector<LayoutIndex>& curNeighbors = (*m_neighbors[lvl])[dit()];

      for (LayoutIterator lit = dbl.layoutIterator(); lit.ok(); ++lit){
	const Box& otherBox = dbl.get(lit());

	if(otherBox != box){
	  Box grownBox = otherBox;
	  grownBox.grow(1);

	  if(box.intersects(grownBox)){
	    curNeighbors.push_back(lit());
	  }
	}
      }
    }
  }
}

void phase_realm::define_levelset(const int a_lmin, const int a_numGhost){
  CH_TIME("phase_realm::define_levelset");
  if(m_verbosity > 2){
    pout() << "phase_realm::define_leveset" << endl;
  }

  const bool do_this_operator = this->query_operator(s_levelset);

  m_levelset.resize(1 + m_finest_level);

  if(do_this_operator){

    const int comp  = 0;
    const int ncomp = 1;

    for (int lvl = a_lmin; lvl <= m_finest_level; lvl++){
      const Real dx = m_dx[lvl];

      m_levelset[lvl] = RefCountedPtr<LevelData<FArrayBox> > (new LevelData<FArrayBox>(m_grids[lvl], ncomp, a_numGhost*IntVect::Unit));

      for (DataIterator dit(m_grids[lvl]); dit.ok(); ++dit){
	FArrayBox& fab = (*m_levelset[lvl])[dit()];
	const Box bx = fab.box();

	if(!m_baseif.isNull()){
	  for (BoxIterator bit(bx); bit.ok(); ++bit){
	    const IntVect iv = bit();
	    const RealVect pos = m_prob_lo + (0.5*RealVect::Unit + RealVect(iv))*dx;

	    fab(iv, comp) = m_baseif->value(pos); 
	  }
	}
	else{
	  fab.setVal(-1.23456789, comp);
	}
      }
    }
  }
}

void phase_realm::define_eb_coar_ave(const int a_lmin){
  CH_TIME("phase_realm::define_eb_coar_ave");
  if(m_verbosity > 2){
    pout() << "phase_realm::define_eb_coar_ave" << endl;
  }

  const bool do_this_operator = this->query_operator(s_eb_coar_ave);

  m_coarave.resize(1 + m_finest_level);
  
  if(do_this_operator){
    
    const int comps = SpaceDim;

    for (int lvl = a_lmin; lvl <= m_finest_level; lvl++){

      const bool has_coar = lvl > 0;

      if(has_coar){
	m_coarave[lvl] = RefCountedPtr<ebcoarseaverage> (new ebcoarseaverage(m_grids[lvl],
									     m_grids[lvl-1],
									     m_ebisl[lvl],
									     m_ebisl[lvl-1],
									     m_domains[lvl-1],
									     m_ref_ratios[lvl-1],
									     comps,
									     &(*m_ebis)));
      }
    }
  }
}

void phase_realm::define_eb_quad_cfi(const int a_lmin){
  CH_TIME("phase_realm::define_eb_quad_cfi");
  if(m_verbosity > 2){
    pout() << "phase_realm::define_eb_quad_cfi" << endl;
  }

  const bool do_this_operator = this->query_operator(s_eb_quad_cfi);

  m_quadcfi.resize(1 + m_finest_level);
  m_old_quadcfi.resize(1 + m_finest_level);
  
  if(do_this_operator){

    const int ncomps = SpaceDim;

    for (int lvl = a_lmin; lvl <= m_finest_level; lvl++){

      const bool has_coar = lvl > 0;
      
      if(has_coar){
	const LayoutData<IntVectSet>& cfivs = *m_eblg[lvl]->getCFIVS();

	m_quadcfi[lvl] = RefCountedPtr<nwoebquadcfinterp> (new nwoebquadcfinterp(m_grids[lvl],
										 m_grids[lvl-1],
										 m_ebisl[lvl],
										 m_ebisl[lvl-1],
										 m_domains[lvl-1],
										 m_ref_ratios[lvl-1],
										 ncomps,
										 m_dx[lvl],
										 m_num_ghost,
										 cfivs,
										 m_ebis));

	m_old_quadcfi[lvl] = RefCountedPtr<EBQuadCFInterp> (new EBQuadCFInterp(m_grids[lvl],
									       m_grids[lvl-1],
									       m_ebisl[lvl],
									       m_ebisl[lvl-1],
									       m_domains[lvl-1],
									       m_ref_ratios[lvl-1],
									       1,
									       cfivs,
									       &(*m_ebis)));
      }
    }
  }
}


void phase_realm::define_fillpatch(const int a_lmin){
  CH_TIME("phase_realm::define_fillpatch");
  if(m_verbosity > 2){
    pout() << "phase_realm::define_fillpatch" << endl;
  }

  const bool do_this_operator = this->query_operator(s_eb_fill_patch);

  m_pwl_fillpatch.resize(1 + m_finest_level);
  
  if(do_this_operator){
    
    const int comps     = SpaceDim;

    // Should these be input somehow?
    const int radius    = 1;
    const IntVect ghost = m_num_ghost*IntVect::Unit;


    for (int lvl = a_lmin; lvl <= m_finest_level; lvl++){

      const bool has_coar = lvl > 0;

      if(has_coar){
	const LayoutData<IntVectSet>& cfivs = *m_eblg[lvl]->getCFIVS();
	m_pwl_fillpatch[lvl] = RefCountedPtr<AggEBPWLFillPatch> (new AggEBPWLFillPatch(m_grids[lvl],
										       m_grids[lvl-1],
										       m_ebisl[lvl],
										       m_ebisl[lvl-1],
										       m_domains[lvl-1],
										       m_ref_ratios[lvl-1],
										       comps,
										       radius,
										       ghost,
										       !m_ebcf,
										       &(*m_ebis)));
      }
    }
  }
}


void phase_realm::define_ebpwl_interp(const int a_lmin){
  CH_TIME("phase_realm::define_ebpwl_interp");
  if(m_verbosity > 2){
    pout() << "phase_realm::define_ebpwl_interp" << endl;
  }

  const bool do_this_operator = this->query_operator(s_eb_pwl_interp);

  m_pwl_interp.resize(1 + m_finest_level);

  if(do_this_operator){

	
    const int comps     = SpaceDim;

    // Should these be input somehow?
    const int radius    = 1;
    const IntVect ghost = m_num_ghost*IntVect::Unit;

    for (int lvl = a_lmin; lvl <= m_finest_level; lvl++){

      const bool has_coar = lvl > 0;

      if(has_coar){
	m_pwl_interp[lvl] = RefCountedPtr<EBPWLFineInterp> (new EBPWLFineInterp(m_grids[lvl],
										m_grids[lvl-1],
										m_ebisl[lvl],
										m_ebisl[lvl-1],
										m_domains[lvl-1],
										m_ref_ratios[lvl-1],
										comps,
										&(*m_ebis)));
      }
    }
  }
}

void phase_realm::define_ebmg_interp(const int a_lmin){
  CH_TIME("phase_realm::define_ebmg_interp");
  if(m_verbosity > 2){
    pout() << "phase_realm::define_ebmg_interp" << endl;
  }

  const bool do_this_operator = this->query_operator(s_eb_mg_interp);

  m_ebmg_interp.resize(1 + m_finest_level);

  if(do_this_operator){

    
    const int ncomps    = 1;

    // Should these be input somehow?
    const int radius    = 1;
    const IntVect ghost = m_num_ghost*IntVect::Unit;

    for (int lvl = a_lmin; lvl <= m_finest_level; lvl++){

      const bool has_coar = lvl > 0;

      if(has_coar){

	m_ebmg_interp[lvl] = RefCountedPtr<EBMGInterp> (new EBMGInterp(m_grids[lvl],
								       m_grids[lvl-1],
								       m_ebisl[lvl],
								       m_ebisl[lvl-1],
								       m_domains[lvl-1],
								       m_ref_ratios[lvl-1],
								       SpaceDim,
								       &(*m_ebis),
								       m_num_ghost*IntVect::Unit));
      }
    }
  }
}

void phase_realm::define_flux_reg(const int a_lmin, const int a_regsize){
  CH_TIME("phase_realm::define_flux_reg");
  if(m_verbosity > 2){
    pout() << "phase_realm::define_flux_reg" << endl;
  }

  const bool do_this_operator = this->query_operator(s_eb_flux_reg);

  m_flux_reg.resize(1 + m_finest_level);

  if(do_this_operator){
    
    const int comps = a_regsize;
    
    for (int lvl = Max(0,a_lmin-1); lvl <= m_finest_level; lvl++){

      const bool has_fine = lvl < m_finest_level;


      if(has_fine){
	m_flux_reg[lvl] = RefCountedPtr<EBFluxRegister> (new EBFastFluxRegister(m_grids[lvl+1],
										m_grids[lvl],
										m_ebisl[lvl+1],
										m_ebisl[lvl],
										m_domains[lvl].domainBox(),
										m_ref_ratios[lvl],
										comps,
										&(*m_ebis)));
      }
    }
  }
}

void phase_realm::define_redist_oper(const int a_lmin, const int a_regsize){
  CH_TIME("phase_realm::define_redist_oper");
  if(m_verbosity > 2){
    pout() << "phase_realm::define_redist_oper" << endl;
  }

  // TLDR: All these operators either do stuff on an AMR level, or between a coarse and a fine level. The entries
  //       live on these levels:
  //
  //       Oper                        Level
  //       EBFineToCoar [l,  l-1]    l
  //       EBCoarToFine [l-1,l  ] 
  //       when level a_lmin changed we need to update fine-to-coar

  const bool do_this_operator = this->query_operator(s_eb_redist);

  m_level_redist.resize(1 + m_finest_level);
  m_fine_to_coar_redist.resize(1 + m_finest_level);
  m_coar_to_coar_redist.resize(1 + m_finest_level);
  m_coar_to_fine_redist.resize(1 + m_finest_level);

  if(do_this_operator){
    
    const int comps = a_regsize;

    for (int lvl = Max(0, a_lmin-1); lvl <= m_finest_level; lvl++){

      const bool has_coar = lvl > 0;
      const bool has_fine = lvl < m_finest_level;

	      
      if(lvl >= a_lmin){
	m_level_redist[lvl] = RefCountedPtr<EBLevelRedist> (new EBLevelRedist(m_grids[lvl],
									      m_ebisl[lvl],
									      m_domains[lvl],
									      comps,
									      m_redist_rad));
      }

    
      if(m_ebcf){
	if(has_coar){

	  // TLDR: The fine-to-coar redistribution operator that transfers from the fine level to the coar level
	  //       obviously lives on the fine level. But since a_lmin is the coarsest level that changed, we only
	  //       need to update this if lvl >= a_lmin
	  if(lvl >= a_lmin){

	    auto redist = RefCountedPtr<EBFastFineToCoarRedist> (new EBFastFineToCoarRedist());
	    redist->define(*m_eblg[lvl],
			   *m_eblg[lvl-1],
			   *m_neighbors[lvl],
			   *m_neighbors[lvl-1],
			   m_ref_ratios[lvl-1],
			   comps,
			   m_redist_rad);
	    m_fine_to_coar_redist[lvl] = RefCountedPtr<EBFineToCoarRedist> (redist);
	    
	  }
	}

	if(has_fine){
	  // TLDR: The coar-to-fine redistribution operator transfers from the coarse level and to the fine level and
	  //       therefore lives on the coarse level. Since a_lmin is the coarsest level that changed, we need to update
	  //       if lvl >= a_lmin-1
	  if(lvl >= a_lmin-1){

	    auto c2f_redist = RefCountedPtr<EBFastCoarToFineRedist> (new EBFastCoarToFineRedist());
	    c2f_redist->define(*m_eblg[lvl+1],
			       *m_eblg[lvl],
			       *m_neighbors[lvl+1],
			       *m_neighbors[lvl],
			       m_ref_ratios[lvl],
			       comps,
			       m_redist_rad);
	    m_coar_to_fine_redist[lvl] = RefCountedPtr<EBCoarToFineRedist> (c2f_redist);


	    auto c2c_redist = RefCountedPtr<EBFastCoarToCoarRedist> (new EBFastCoarToCoarRedist());
	    c2c_redist->define(*m_eblg[lvl+1],
			       *m_eblg[lvl],
			       *m_neighbors[lvl+1],
			       *m_neighbors[lvl],
			       m_ref_ratios[lvl],
			       comps,
			       m_redist_rad);
	    m_coar_to_coar_redist[lvl] = RefCountedPtr<EBCoarToCoarRedist> (c2c_redist);
	  }
	}
      }
    }
  }
}

void phase_realm::define_gradsten(const int a_lmin){
  CH_TIME("phase_realm::define_gradsten");
  if(m_verbosity > 2){
    pout() << "phase_realm::define_gradsten" << endl;
  }

  const bool do_this_operator = this->query_operator(s_eb_gradient);

  m_gradsten.resize(1 + m_finest_level);

  if(do_this_operator){
    for (int lvl = a_lmin; lvl <= m_finest_level; lvl++){
      const DisjointBoxLayout& dbl = m_grids[lvl];
      const ProblemDomain& domain  = m_domains[lvl];
      const Real dx                = m_dx[lvl];
    
      m_gradsten[lvl] = RefCountedPtr<LayoutData<BaseIVFAB<VoFStencil> > > (new LayoutData<BaseIVFAB<VoFStencil> >(dbl));

      const EBISLayout& ebisl = m_ebisl[lvl];

      LayoutData<IntVectSet>& cfivs = *m_eblg[lvl]->getCFIVS();
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	const Box box          = dbl.get(dit());
	const EBISBox& ebisbox = ebisl[dit()];
	const EBGraph& ebgraph = ebisbox.getEBGraph();

	
	IntVectSet ivs   = ebisbox.getIrregIVS(box);
	for (int dir = 0; dir < SpaceDim; dir++){
	  Box lo_box, hi_box;
	  int has_lo, has_hi;
	
	  EBArith::loHi(lo_box, has_lo, hi_box, has_hi, domain, box, dir);
	
	  if(has_lo){
	    ivs |= IntVectSet(lo_box);
	  }
	  if(has_hi){
	    ivs |= IntVectSet(hi_box);
	  }
	}

	BaseIVFAB<VoFStencil>& vofstencils = (*m_gradsten[lvl])[dit()];
	vofstencils.define(ivs, ebgraph, 1);

	for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
	  const VolIndex& vof = vofit();

	  VoFStencil& sten = vofstencils(vof, 0);
	  sten.clear();
	  for (int dir = 0; dir < SpaceDim; dir++){
	    VoFStencil dirsten;
	    //	    EBArith::getFirstDerivStencil(dirsten, vof, ebisbox, dir, dx, &cfivs[dit()], dir);
	    EBArith::getFirstDerivStencilWidthOne(dirsten, vof, ebisbox, dir, dx, NULL, dir);
	    //	    EBArith::getFirstDerivStencil(dirsten, vof, ebisbox, dir, dx, NULL, dir);
	    sten += dirsten;
	  }
	}
      }
    }
  }
}

void phase_realm::define_copier(const int a_lmin){
  CH_TIME("phase_realm::define_copier");
  if(m_verbosity > 3){
    pout() << "phase_realm::define_copier" << endl;
  }

  const bool do_this_operator = this->query_operator(s_eb_copier);

  m_copier.resize(1 + m_finest_level);
  m_reverse_copier.resize(1 + m_finest_level);

  if(do_this_operator){
    
    for (int lvl = a_lmin; lvl <= m_finest_level; lvl++){
      m_copier[lvl] = RefCountedPtr<Copier> (new Copier(m_grids[lvl],
							m_grids[lvl],
							m_domains[lvl],
							m_num_ghost*IntVect::Unit,
							true));
      
      m_reverse_copier[lvl] = RefCountedPtr<Copier> (new Copier(m_grids[lvl],
								m_grids[lvl],
								m_domains[lvl],
								m_num_ghost*IntVect::Unit,
								true));
      m_reverse_copier[lvl]->reverse();
    }
  }
}


void phase_realm::define_ghostcloud(const int a_lmin){
  CH_TIME("phase_realm::define_ghostcloud");
  if(m_verbosity > 3){
    pout() << "phase_realm::define_ghostcloud" << endl;
  }

  const bool do_this_operator = this->query_operator(s_eb_copier);

  m_ghostclouds.resize(1 + m_finest_level);

  if(do_this_operator){
    for (int lvl = a_lmin; lvl <= m_finest_level; lvl++){
      const bool has_coar = lvl > 0;

      if(has_coar){
	m_ghostclouds[lvl] = RefCountedPtr<EBGhostCloud> (new EBGhostCloud(m_grids[lvl-1],
									   m_grids[lvl],
									   *m_eblg[lvl-1],
									   *m_eblg[lvl],
									   m_domains[lvl-1],
									   m_domains[lvl],
									   m_ref_ratios[lvl-1],
									   1,
									   m_num_ghost));
      }
    }
  }
}


void phase_realm::define_irreg_sten(){
  CH_TIME("phase_realm::define_irreg_sten");
  if(m_verbosity > 2){
    pout() << "phase_realm::define_irreg_sten" << endl;
  }

  const bool do_this_operator = this->query_operator(s_eb_irreg_interp);

  if(do_this_operator){
    
    const int order = 1;
    const int rad   = 1;

    m_centroid_interp = RefCountedPtr<irreg_amr_stencil<centroid_interp> >
      (new irreg_amr_stencil<centroid_interp>(m_grids,
					      m_ebisl,
					      m_domains,
					      m_dx,
					      m_finest_level,
					      order,
					      rad,
					      m_centroid_stencil));
      
    m_eb_centroid_interp = RefCountedPtr<irreg_amr_stencil<eb_centroid_interp> >
      (new irreg_amr_stencil<eb_centroid_interp>(m_grids,
						 m_ebisl,
						 m_domains,
						 m_dx,
						 m_finest_level,
						 order,
						 rad,
						 m_eb_stencil));
  }
}

void phase_realm::define_noncons_sten(){
  CH_TIME("phase_realm::define_noncons_sten");
  if(m_verbosity > 2){
    pout() << "phase_realm::define_noncons_sten" << endl;
  }

  const bool do_this_operator = this->query_operator(s_eb_noncons_div);

  if(do_this_operator){
    const int order = 1; // Dummy argument
    const int rad   = m_redist_rad;

    
    m_noncons_div = RefCountedPtr<irreg_amr_stencil<noncons_div> >
      (new irreg_amr_stencil<noncons_div>(m_grids,
					  m_ebisl,
					  m_domains,
					  m_dx,
					  m_finest_level,
					  order,                 // Dummy argument
					  m_redist_rad,
					  m_centroid_stencil));  // Dummy argumement
  }
}

const RefCountedPtr<EBIndexSpace>& phase_realm::get_ebis() {
  return m_ebis;
}

Vector<int>& phase_realm::get_ref_rat() {
  return m_ref_ratios;
}

Vector<Real>& phase_realm::get_dx() {
  return m_dx;
}

Vector<DisjointBoxLayout>& phase_realm::get_grids() {
  return m_grids;
}

Vector<ProblemDomain>& phase_realm::get_domains() {
  return m_domains;
}

Vector<EBISLayout>& phase_realm::get_ebisl() {
  return m_ebisl;
}

Vector<RefCountedPtr<EBLevelGrid> >& phase_realm::get_eblg() {
  return m_eblg;
}

Vector<RefCountedPtr<LayoutData<Vector<LayoutIndex> > > >& phase_realm::get_neighbors() {
  return m_neighbors;
}
Vector<RefCountedPtr<LayoutData<VoFIterator> > >& phase_realm::get_vofit() {
  return m_vofiter;
}

irreg_amr_stencil<centroid_interp>& phase_realm::get_centroid_interp_stencils() {
  return *m_centroid_interp;
}

irreg_amr_stencil<eb_centroid_interp>& phase_realm::get_eb_centroid_interp_stencils() {
  return *m_eb_centroid_interp;
}

Vector<RefCountedPtr<LayoutData<BaseIVFAB<VoFStencil> > > >& phase_realm::get_gradsten(){
  if(!this->query_operator(s_eb_gradient)) MayDay::Abort("phase_realm::get_gradsten - operator not registered!");
  return m_gradsten;
}

// Throw errors if the operator does not exist

irreg_amr_stencil<noncons_div>& phase_realm::get_noncons_div_stencils() {
  if(!this->query_operator(s_eb_noncons_div)) MayDay::Abort("phase_realm::get_non_cons_div_stencils - operator not registered!");
  
  return *m_noncons_div;
}

EBAMRFAB& phase_realm::get_levelset() {
  if(!this->query_operator(s_levelset)) MayDay::Abort("phase_realm::get_levelset - operator not registered!");

  return m_levelset;
}

Vector<RefCountedPtr<ebcoarseaverage> >& phase_realm::get_coarave() {
  if(!this->query_operator(s_eb_coar_ave)) MayDay::Abort("phase_realm::get_coarave - operator not registered!");
  
  return m_coarave;
}

Vector<RefCountedPtr<EBGhostCloud> >& phase_realm::get_ghostcloud() {
  if(!this->query_operator(s_eb_ghostcloud)) MayDay::Abort("phase_realm::get_ghostcloud - operator not registered!");
  
  return m_ghostclouds;
}

Vector<RefCountedPtr<nwoebquadcfinterp> >& phase_realm::get_quadcfi() {
  if(!this->query_operator(s_eb_quad_cfi)) MayDay::Abort("phase_realm::get_quadcfi - operator not registered!");
  
  return m_quadcfi;
}

Vector<RefCountedPtr<EBQuadCFInterp> >& phase_realm::get_old_quadcfi() {
  if(!this->query_operator(s_eb_quad_cfi)) MayDay::Abort("phase_realm::get_old_quadcfi - operator not registered!");
  
  return m_old_quadcfi;
}

Vector<RefCountedPtr<AggEBPWLFillPatch> >& phase_realm::get_fillpatch() {
  if(!this->query_operator(s_eb_fill_patch)) MayDay::Abort("phase_realm::get_fillpatch - operator not registered!");
  
  return m_pwl_fillpatch;
}

Vector<RefCountedPtr<EBPWLFineInterp> >& phase_realm::get_eb_pwl_interp() {
  if(!this->query_operator(s_eb_pwl_interp)) MayDay::Abort("phase_realm::get_eb_pwl_interp - operator not registered!");
  
  return m_pwl_interp;
}

Vector<RefCountedPtr<EBMGInterp> >& phase_realm::get_eb_mg_interp() {
  if(!this->query_operator(s_eb_mg_interp)) MayDay::Abort("phase_realm::get_eb_mg_interp - operator not registered!");
  
  return m_ebmg_interp;
}

Vector<RefCountedPtr<EBFluxRegister> >&  phase_realm::get_flux_reg() {
  if(!this->query_operator(s_eb_flux_reg)) MayDay::Abort("phase_realm::get_flux_reg - operator not registered!");

  return m_flux_reg;
}

Vector<RefCountedPtr<EBLevelRedist> >&  phase_realm::get_level_redist() {
  if(!this->query_operator(s_eb_redist)) MayDay::Abort("phase_realm::get_level_redist - operator not registered!");

  return m_level_redist;
}

Vector<RefCountedPtr<EBCoarToFineRedist> >&  phase_realm::get_coar_to_fine_redist() {
  if(!this->query_operator(s_eb_redist)) MayDay::Abort("phase_realm::get_coar_to_fine_redist - operator not registered!");

  return m_coar_to_fine_redist;
}

Vector<RefCountedPtr<EBCoarToCoarRedist> >&  phase_realm::get_coar_to_coar_redist() {
  if(!this->query_operator(s_eb_redist)) MayDay::Abort("phase_realm::get_coar_to_coar - operator not registered!");

  return m_coar_to_coar_redist;
}

Vector<RefCountedPtr<EBFineToCoarRedist> >&  phase_realm::get_fine_to_coar_redist() {
  if(!this->query_operator(s_eb_redist)) MayDay::Abort("phase_realm::get_fine_to_coar_redist - operator not registered!");

  return m_fine_to_coar_redist;
}

Vector<RefCountedPtr<Copier> >& phase_realm::get_copier() {
  if(!this->query_operator(s_eb_copier)) MayDay::Abort("phase_realm::get_copier - operator not registered!");

  return m_copier;
}

Vector<RefCountedPtr<Copier> >& phase_realm::get_reverse_copier() {
  if(!this->query_operator(s_eb_copier)) MayDay::Abort("phase_realm::get_reverse_copier - operator not registered!");
  
  return m_reverse_copier;
}
