/*!
  @file   realm.cpp
  @brief  Implementation of realm.H
  @author Robert Marskar
  @date   July 2020
*/

#include "realm.H"
#include "EBFastFineToCoarRedist.H"
#include "EBFastCoarToFineRedist.H"
#include "EBFastCoarToCoarRedist.H"
#include "load_balance.H"

#include <EBArith.H>


realm::realm(){
  m_defined = false;
  m_verbosity = -1;
}

realm::~realm(){

}

void realm::define(const Vector<DisjointBoxLayout>& a_grids,
		   const Vector<ProblemDomain>& a_domains,
		   const Vector<int>& a_ref_rat,
		   const Vector<Real>& a_dx,
		   const int a_finest_level,
		   const int a_ebghost,
		   const int a_num_ghost,
		   const int a_redist_rad,
		   const stencil_type::which_type a_centroid_stencil,
		   const stencil_type::which_type a_eb_stencil,
		   const bool a_ebcf,
		   const RefCountedPtr<EBIndexSpace>& a_ebis){

  m_finest_level = a_finest_level;
  m_grids = a_grids;
  m_domains = a_domains;
  m_ref_ratios = a_ref_rat;
  m_dx = a_dx;
  m_ebis = a_ebis;
  m_ebcf = a_ebcf;
  m_ebghost = a_ebghost;
  m_num_ghost = a_num_ghost;
  m_redist_rad = a_redist_rad;
  m_centroid_stencil = a_centroid_stencil;
  m_eb_stencil = a_eb_stencil;

  m_defined = true;
}

void realm::set_grids(const Vector<DisjointBoxLayout>& a_grids, const int a_finest_level){
  CH_TIME("realm::set_grids");
  if(m_verbosity > 5){
    pout() << "realm::set_grids" << endl;
  }
  
  m_grids = a_grids;
  m_finest_level = a_finest_level;
}

void realm::regrid_base(const int a_lmin, const int a_lmax, const int a_hardcap){
  CH_TIME("realm::regrid_base");
  if(m_verbosity > 5){
    pout() << "realm::regrid_base" << endl;
  }

  this->define_neighbors(a_lmin);
  this->define_eblevelgrid(a_lmin);
  this->define_vofiter(a_lmin);

  // Missing MG grids...
}

void realm::regrid_operators(const int a_lmin, const int a_lmax, const int a_regsize){
  CH_TIME("realm::regrid_operators_phase");
  if(m_verbosity > 5){
    pout() << "realm::regrid_operators" << endl;
  }

  this->define_eb_coar_ave(a_lmin);             // Define ebcoarseaverage on both phases
  this->define_eb_quad_cfi(a_lmin);             // Define nwoebquadcfinterp on both phases.
  this->define_fillpatch(a_lmin);               // Define operator for piecewise linear interpolation of ghost cells
  this->define_ebpwl_interp(a_lmin);            // Define interpolator for piecewise interpolation of interior points
  this->define_ebmg_interp(a_lmin);             // Define interpolator used for e.g. multigrid (or piecewise constant)
  this->define_flux_reg(a_lmin,a_regsize);      // Define flux register (phase::gas only)
  this->define_redist_oper(a_lmin, a_regsize);  // Define redistribution (phase::gas only)
  this->define_gradsten(a_lmin);                // Make stencils for computing gradients
  this->define_irreg_sten();                    // Make stencils for doing interpolation to centroids
  this->define_noncons_sten();                  // Make stencils for nonconservative averaging
  this->define_copier(a_lmin);                  // Make stencils for copier
  this->define_ghostcloud(a_lmin);              // Make stencils for ghost clouds with particle depositions
}

void realm::register_operator(const std::string a_operator){
  CH_TIME("realm::register_operator");
  if(m_verbosity > 5){
    pout() << "realm::register_operator" << endl;
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
       a_operator.compare(s_eb_mg_interp)    == 0)){

    const std::string str = "amr_mesh::register_operator - unknown operator '" + a_operator + "' requested";
    MayDay::Abort(str.c_str());
  }
     

  m_operator_map.emplace(a_operator, true);
}

bool realm::query_operator(const std::string a_operator){
  CH_TIME("realm::query_operator");
  if(m_verbosity > 5){
    pout() << "realm::query_operator" << endl;
  }
  
  return m_operator_map[a_operator];
}
  

void realm::define_eblevelgrid(const int a_lmin){
  CH_TIME("realm::define_eblevelgrid");
  if(m_verbosity > 2){
    pout() << "realm::define_eblevelgrid" << endl;
  }

  m_eblg.resize(1 + m_finest_level);
  m_ebisl.resize(1 + m_finest_level);

  for (int lvl = a_lmin; lvl <= m_finest_level; lvl++){
    m_eblg[lvl]  = RefCountedPtr<EBLevelGrid> (new EBLevelGrid(m_grids[lvl], m_domains[lvl], m_ebghost, &(*m_ebis)));
    
    if(lvl > 0) m_eblg[lvl]->setMaxCoarseningRatio(m_ref_ratios[lvl-1], &(*m_ebis));

    m_ebisl[lvl] = m_eblg[lvl]->getEBISL();
  }
}
void realm::define_vofiter(const int a_lmin){
  CH_TIME("realm::define_vofiter");
  if(m_verbosity > 2){
    pout() << "realm::define_vofiter" << endl;
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

void realm::define_neighbors(const int a_lmin){
  CH_TIME("realm::define_neighbors");
  if(m_verbosity > 2){
    pout() << "realm::define_neighbors" << endl;
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

void realm::define_eb_coar_ave(const int a_lmin){
  CH_TIME("realm::define_eb_coar_ave");
  if(m_verbosity > 2){
    pout() << "realm::define_eb_coar_ave" << endl;
  }

  const bool do_this_operator = this->query_operator(s_eb_coar_ave);
  
  if(do_this_operator){

    m_coarave.resize(1 + m_finest_level);
    
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

void realm::define_eb_quad_cfi(const int a_lmin){
  CH_TIME("realm::define_eb_quad_cfi");
  if(m_verbosity > 2){
    pout() << "realm::define_eb_quad_cfi" << endl;
  }

  const bool do_this_operator = this->query_operator(s_eb_quad_cfi);
  
  if(do_this_operator){

    m_quadcfi.resize(1 + m_finest_level);
    m_old_quadcfi.resize(1 + m_finest_level);

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


void realm::define_fillpatch(const int a_lmin){
  CH_TIME("realm::define_fillpatch");
  if(m_verbosity > 2){
    pout() << "realm::define_fillpatch" << endl;
  }

  const bool do_this_operator = this->query_operator(s_eb_fill_patch);
  
  if(do_this_operator){
    m_pwl_fillpatch.resize(1 + m_finest_level);
    
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


void realm::define_ebpwl_interp(const int a_lmin){
  CH_TIME("amr_mesh::define_ebpwl_interp");
  if(m_verbosity > 2){
    pout() << "amr_mesh::define_ebpwl_interp" << endl;
  }

  const bool do_this_operator = this->query_operator(s_eb_pwl_interp);

  if(do_this_operator){
    m_pwl_interp.resize(1 + m_finest_level);
	
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

void realm::define_ebmg_interp(const int a_lmin){
  CH_TIME("realm::define_ebmg_interp");
  if(m_verbosity > 2){
    pout() << "realm::define_ebmg_interp" << endl;
  }

  const bool do_this_operator = this->query_operator(s_eb_mg_interp);

  if(do_this_operator){
    m_ebmg_interp.resize(1 + m_finest_level);
    
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

void realm::define_flux_reg(const int a_lmin, const int a_regsize){
  CH_TIME("realm::define_flux_reg");
  if(m_verbosity > 2){
    pout() << "realm::define_flux_reg" << endl;
  }

  const bool do_this_operator = this->query_operator(s_eb_flux_reg);

  if(do_this_operator){
    m_flux_reg.resize(1 + m_finest_level);
    
    const int comps = a_regsize;
    
    for (int lvl = Max(0,a_lmin-1); lvl <= m_finest_level; lvl++){

      const bool has_fine = lvl < m_finest_level;


      if(has_fine){
	m_flux_reg[lvl] = RefCountedPtr<EBFluxRegister> (new EBFluxRegister(m_grids[lvl+1],
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


void realm::define_redist_oper(const int a_lmin, const int a_regsize){
  CH_TIME("realm::define_redist_oper");
  if(m_verbosity > 2){
    pout() << "realm::define_redist_oper" << endl;
  }

  // TLDR: All these operators either do stuff on an AMR level, or between a coarse and a fine level. The entries
  //       live on these levels:
  //
  //       Oper                        Level
  //       EBFineToCoar [l,  l-1]    l
  //       EBCoarToFine [l-1,l  ] 
  //       when level a_lmin changed we need to update fine-to-coar

  const bool do_this_operator = this->query_operator(s_eb_redist);

  if(do_this_operator){
    m_level_redist.resize(1 + m_finest_level);
    m_fine_to_coar_redist.resize(1 + m_finest_level);
    m_coar_to_coar_redist.resize(1 + m_finest_level);
    m_coar_to_fine_redist.resize(1 + m_finest_level);
    
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

void realm::define_gradsten(const int a_lmin){
  CH_TIME("realm::define_gradsten");
  if(m_verbosity > 2){
    pout() << "realm::define_gradsten" << endl;
  }

  const bool do_this_operator = this->query_operator(s_eb_gradient);

  if(do_this_operator){
    m_gradsten.resize(1 + m_finest_level);
    
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
	    EBArith::getFirstDerivStencil(dirsten, vof, ebisbox, dir, dx, &cfivs[dit()], dir);
	    sten += dirsten;
	  }
	}
      }
    }
  }
}

void realm::define_copier(const int a_lmin){
  CH_TIME("realm::define_copier");
  if(m_verbosity > 3){
    pout() << "realm::define_copier" << endl;
  }

  const bool do_this_operator = this->query_operator(s_eb_copier);

  if(do_this_operator){
    m_copier.resize(1 + m_finest_level);
    m_reverse_copier.resize(1 + m_finest_level);
    
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


void realm::define_ghostcloud(const int a_lmin){
  CH_TIME("realm::define_ghostcloud");
  if(m_verbosity > 3){
    pout() << "realm::define_ghostcloud" << endl;
  }

  const bool do_this_operator = this->query_operator(s_eb_copier);

  if(do_this_operator){
    m_ghostclouds.resize(1 + m_finest_level);

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


void realm::define_irreg_sten(){
  CH_TIME("realm::define_irreg_sten");
  if(m_verbosity > 2){
    pout() << "realm::define_irreg_sten" << endl;
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

void realm::define_noncons_sten(){
  CH_TIME("realm::define_noncons_sten");
  if(m_verbosity > 2){
    pout() << "realm::define_noncons_sten" << endl;
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
