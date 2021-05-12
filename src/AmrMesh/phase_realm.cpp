/*!
  @file   phase_realm.cpp
  @brief  Implementation of phase_realm.H
  @author Robert Marskar
  @date   July 2020
*/

#include "phase_realm.H"
#include "EBFastFineToCoarRedist.H"
#include "EBFastCoarToFineRedist.H"
#include "CD_EbFastCoarToCoarRedist.H"
#include "load_balance.H"
#include "EBFastFluxRegister.H"

#include <EBArith.H>

#include "CD_NamespaceHeader.H"
phase_realm::phase_realm(){
  m_defined   = false;
  m_verbosity = -1;

  // Always do this shit. 
  this->registerOperator(s_eb_gradient);
  this->registerOperator(s_eb_irreg_interp);
}

phase_realm::~phase_realm(){

}

void phase_realm::define(const Vector<DisjointBoxLayout>& a_grids,
			 const Vector<ProblemDomain>& a_domains,
			 const Vector<int>& a_ref_rat,
			 const Vector<Real>& a_dx,
			 const RealVect a_probLo,
			 const int a_finest_level,
			 const int a_ebghost,
			 const int a_num_ghost,
			 const int a_lsf_ghost,
			 const int a_redist_rad,
			 const IrregStencil::StencilType a_centroid_stencil,
			 const IrregStencil::StencilType a_eb_stencil,
			 const bool a_ebcf,
			 const RefCountedPtr<BaseIF>&       a_baseif,
			 const RefCountedPtr<EBIndexSpace>& a_ebis){

  m_ebis = a_ebis;
  m_finestLevel = a_finest_level;
  m_grids = a_grids;
  m_domains = a_domains;
  m_refinementRatios = a_ref_rat;
  m_dx = a_dx;
  m_hasEbCf = a_ebcf;
  m_numEbGhostsCells = a_ebghost;
  m_numGhostCells = a_num_ghost;
  m_numLsfGhostCells = a_lsf_ghost;
  m_redistributionRadius = a_redist_rad;
  m_centroidStencilType = a_centroid_stencil;
  m_ebCentroidStencilType = a_eb_stencil;
  m_baseif = a_baseif;
  m_probLo = a_probLo;
  
  if(!m_ebis.isNull()){
    m_defined = true;
  }
}

void phase_realm::setGrids(const Vector<DisjointBoxLayout>& a_grids, const int a_finest_level){
  CH_TIME("phase_realm::setGrids");
  if(m_verbosity > 5){
    pout() << "phase_realm::setGrids" << endl;
  }

  if(m_defined){
    m_grids = a_grids;
    m_finestLevel = a_finest_level;
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

void phase_realm::regridOperators(const int a_lmin, const int a_lmax, const int a_regsize){
  CH_TIME("phase_realm::regridOperators_phase");
  if(m_verbosity > 5){
    pout() << "phase_realm::regridOperators" << endl;
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
    this->define_levelset(a_lmin, m_numLsfGhostCells);   // Defining levelset
  }
}

void phase_realm::registerOperator(const std::string a_operator){
  CH_TIME("phase_realm::registerOperator");
  if(m_verbosity > 5){
    pout() << "phase_realm::registerOperator" << endl;
  } 

  // These are the supported operators - issue an error if we ask for something that is not supported. 
  if(!(a_operator.compare(s_eb_coar_ave)     == 0 ||
       a_operator.compare(s_eb_quad_cfi)     == 0 ||
       a_operator.compare(s_eb_fill_patch)   == 0 ||
       a_operator.compare(s_eb_pwl_interp)   == 0 ||
       a_operator.compare(s_eb_flux_reg)     == 0 ||
       a_operator.compare(s_eb_redist)       == 0 ||
       a_operator.compare(s_eb_NonConservativeDivergenceStencil)  == 0 ||
       a_operator.compare(s_eb_copier)       == 0 ||
       a_operator.compare(s_eb_ghostcloud)   == 0 ||
       a_operator.compare(s_eb_gradient)     == 0 ||
       a_operator.compare(s_eb_irreg_interp) == 0 ||
       a_operator.compare(s_eb_mg_interp)    == 0 ||
       a_operator.compare(s_levelset)        == 0 )){

    const std::string str = "phase_realm::registerOperator - unknown operator '" + a_operator + "' requested";
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

  m_eblg.resize(1 + m_finestLevel);
  m_ebisl.resize(1 + m_finestLevel);

  for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++){
    m_eblg[lvl]  = RefCountedPtr<EBLevelGrid> (new EBLevelGrid(m_grids[lvl], m_domains[lvl], m_numEbGhostsCells, &(*m_ebis)));
    
    if(lvl > 0) m_eblg[lvl]->setMaxCoarseningRatio(m_refinementRatios[lvl-1], &(*m_ebis));

    m_ebisl[lvl] = m_eblg[lvl]->getEBISL();
  }
}
void phase_realm::define_vofiter(const int a_lmin){
  CH_TIME("phase_realm::define_vofiter");
  if(m_verbosity > 2){
    pout() << "phase_realm::define_vofiter" << endl;
  }

  m_vofiter.resize(1 + m_finestLevel);

  for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++){
      
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

  m_neighbors.resize(1 + m_finestLevel);

  for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++){
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

  m_levelset.resize(1 + m_finestLevel);

  if(do_this_operator){

    const int comp  = 0;
    const int ncomp = 1;

    for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++){
      const Real dx = m_dx[lvl];

      m_levelset[lvl] = RefCountedPtr<LevelData<FArrayBox> > (new LevelData<FArrayBox>(m_grids[lvl], ncomp, a_numGhost*IntVect::Unit));

      for (DataIterator dit(m_grids[lvl]); dit.ok(); ++dit){
	FArrayBox& fab = (*m_levelset[lvl])[dit()];
	const Box bx = fab.box();

	if(!m_baseif.isNull()){
	  for (BoxIterator bit(bx); bit.ok(); ++bit){
	    const IntVect iv = bit();
	    const RealVect pos = m_probLo + (0.5*RealVect::Unit + RealVect(iv))*dx;

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

  m_coarave.resize(1 + m_finestLevel);
  
  if(do_this_operator){
    
    const int comps = SpaceDim;

    for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++){

      const bool has_coar = lvl > 0;

      if(has_coar){
	m_coarave[lvl] = RefCountedPtr<ebcoarseaverage> (new ebcoarseaverage(m_grids[lvl],
									     m_grids[lvl-1],
									     m_ebisl[lvl],
									     m_ebisl[lvl-1],
									     m_domains[lvl-1],
									     m_refinementRatios[lvl-1],
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

  m_quadcfi.resize(1 + m_finestLevel);
  m_old_quadcfi.resize(1 + m_finestLevel);
  
  if(do_this_operator){

    const int ncomps = SpaceDim;

    for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++){

      const bool has_coar = lvl > 0;
      
      if(has_coar){
	const LayoutData<IntVectSet>& cfivs = *m_eblg[lvl]->getCFIVS();

	m_quadcfi[lvl] = RefCountedPtr<nwoebquadcfinterp> (new nwoebquadcfinterp(m_grids[lvl],
										 m_grids[lvl-1],
										 m_ebisl[lvl],
										 m_ebisl[lvl-1],
										 m_domains[lvl-1],
										 m_refinementRatios[lvl-1],
										 ncomps,
										 m_dx[lvl],
										 m_numGhostCells,
										 cfivs,
										 m_ebis));

	m_old_quadcfi[lvl] = RefCountedPtr<EBQuadCFInterp> (new EBQuadCFInterp(m_grids[lvl],
									       m_grids[lvl-1],
									       m_ebisl[lvl],
									       m_ebisl[lvl-1],
									       m_domains[lvl-1],
									       m_refinementRatios[lvl-1],
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

  m_pwl_fillpatch.resize(1 + m_finestLevel);
  
  if(do_this_operator){
    
    const int comps     = SpaceDim;

    // Should these be input somehow?
    const int radius    = 1;
    const IntVect ghost = m_numGhostCells*IntVect::Unit;


    for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++){

      const bool has_coar = lvl > 0;

      if(has_coar){
	const LayoutData<IntVectSet>& cfivs = *m_eblg[lvl]->getCFIVS();
	m_pwl_fillpatch[lvl] = RefCountedPtr<AggEBPWLFillPatch> (new AggEBPWLFillPatch(m_grids[lvl],
										       m_grids[lvl-1],
										       m_ebisl[lvl],
										       m_ebisl[lvl-1],
										       m_domains[lvl-1],
										       m_refinementRatios[lvl-1],
										       comps,
										       radius,
										       ghost,
										       !m_hasEbCf,
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

  m_pwl_interp.resize(1 + m_finestLevel);

  if(do_this_operator){

	
    const int comps     = SpaceDim;

    // Should these be input somehow?
    const int radius    = 1;
    const IntVect ghost = m_numGhostCells*IntVect::Unit;

    for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++){

      const bool has_coar = lvl > 0;

      if(has_coar){
	m_pwl_interp[lvl] = RefCountedPtr<EBPWLFineInterp> (new EBPWLFineInterp(m_grids[lvl],
										m_grids[lvl-1],
										m_ebisl[lvl],
										m_ebisl[lvl-1],
										m_domains[lvl-1],
										m_refinementRatios[lvl-1],
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

  m_ebmg_interp.resize(1 + m_finestLevel);

  if(do_this_operator){

    
    const int ncomps    = 1;

    // Should these be input somehow?
    const int radius    = 1;
    const IntVect ghost = m_numGhostCells*IntVect::Unit;

    for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++){

      const bool has_coar = lvl > 0;

      if(has_coar){

	m_ebmg_interp[lvl] = RefCountedPtr<EBMGInterp> (new EBMGInterp(m_grids[lvl],
								       m_grids[lvl-1],
								       m_ebisl[lvl],
								       m_ebisl[lvl-1],
								       m_domains[lvl-1],
								       m_refinementRatios[lvl-1],
								       SpaceDim,
								       &(*m_ebis),
								       m_numGhostCells*IntVect::Unit));
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

  m_flux_reg.resize(1 + m_finestLevel);

  if(do_this_operator){
    
    const int comps = a_regsize;
    
    for (int lvl = Max(0,a_lmin-1); lvl <= m_finestLevel; lvl++){

      const bool has_fine = lvl < m_finestLevel;


      if(has_fine){
	m_flux_reg[lvl] = RefCountedPtr<EBFluxRegister> (new EBFastFluxRegister(m_grids[lvl+1],
										m_grids[lvl],
										m_ebisl[lvl+1],
										m_ebisl[lvl],
										m_domains[lvl].domainBox(),
										m_refinementRatios[lvl],
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

  m_level_redist.resize(1 + m_finestLevel);
  m_fine_to_coar_redist.resize(1 + m_finestLevel);
  m_coar_to_coar_redist.resize(1 + m_finestLevel);
  m_coar_to_fine_redist.resize(1 + m_finestLevel);

  if(do_this_operator){
    
    const int comps = a_regsize;

    for (int lvl = Max(0, a_lmin-1); lvl <= m_finestLevel; lvl++){

      const bool has_coar = lvl > 0;
      const bool has_fine = lvl < m_finestLevel;

	      
      if(lvl >= a_lmin){
	m_level_redist[lvl] = RefCountedPtr<EBLevelRedist> (new EBLevelRedist(m_grids[lvl],
									      m_ebisl[lvl],
									      m_domains[lvl],
									      comps,
									      m_redistributionRadius));
      }

    
      if(m_hasEbCf){
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
			   m_refinementRatios[lvl-1],
			   comps,
			   m_redistributionRadius);
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
			       m_refinementRatios[lvl],
			       comps,
			       m_redistributionRadius);
	    m_coar_to_fine_redist[lvl] = RefCountedPtr<EBCoarToFineRedist> (c2f_redist);


	    auto c2c_redist = RefCountedPtr<EbFastCoarToCoarRedist> (new EbFastCoarToCoarRedist());
	    c2c_redist->define(*m_eblg[lvl+1],
			       *m_eblg[lvl],
			       *m_neighbors[lvl+1],
			       *m_neighbors[lvl],
			       m_refinementRatios[lvl],
			       comps,
			       m_redistributionRadius);
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

  m_gradsten.resize(1 + m_finestLevel);

  if(do_this_operator){
    for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++){
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

  m_copier.resize(1 + m_finestLevel);
  m_reverse_copier.resize(1 + m_finestLevel);

  if(do_this_operator){
    
    for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++){
      m_copier[lvl] = RefCountedPtr<Copier> (new Copier(m_grids[lvl],
							m_grids[lvl],
							m_domains[lvl],
							m_numGhostCells*IntVect::Unit,
							true));
      
      m_reverse_copier[lvl] = RefCountedPtr<Copier> (new Copier(m_grids[lvl],
								m_grids[lvl],
								m_domains[lvl],
								m_numGhostCells*IntVect::Unit,
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

  m_ghostclouds.resize(1 + m_finestLevel);

  if(do_this_operator){
    for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++){
      const bool has_coar = lvl > 0;

      if(has_coar){
	m_ghostclouds[lvl] = RefCountedPtr<EBGhostCloud> (new EBGhostCloud(m_grids[lvl-1],
									   m_grids[lvl],
									   *m_eblg[lvl-1],
									   *m_eblg[lvl],
									   m_domains[lvl-1],
									   m_domains[lvl],
									   m_refinementRatios[lvl-1],
									   1,
									   m_numGhostCells));
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

    m_CentroidInterpolationStencil = RefCountedPtr<IrregAmrStencil<CentroidInterpolationStencil> >
      (new IrregAmrStencil<CentroidInterpolationStencil>(m_grids,
					      m_ebisl,
					      m_domains,
					      m_dx,
					      m_finestLevel,
					      order,
					      rad,
					      m_centroidStencilType));
      
    m_EbCentroidInterpolationStencil = RefCountedPtr<IrregAmrStencil<EbCentroidInterpolationStencil> >
      (new IrregAmrStencil<EbCentroidInterpolationStencil>(m_grids,
						 m_ebisl,
						 m_domains,
						 m_dx,
						 m_finestLevel,
						 order,
						 rad,
						 m_ebCentroidStencilType));
  }
}

void phase_realm::define_noncons_sten(){
  CH_TIME("phase_realm::define_noncons_sten");
  if(m_verbosity > 2){
    pout() << "phase_realm::define_noncons_sten" << endl;
  }

  const bool do_this_operator = this->query_operator(s_eb_NonConservativeDivergenceStencil);

  if(do_this_operator){
    const int order = 1; // Dummy argument
    const int rad   = m_redistributionRadius;

    
    m_NonConservativeDivergenceStencil = RefCountedPtr<IrregAmrStencil<NonConservativeDivergenceStencil> >
      (new IrregAmrStencil<NonConservativeDivergenceStencil>(m_grids,
					  m_ebisl,
					  m_domains,
					  m_dx,
					  m_finestLevel,
					  order,                 // Dummy argument
					  m_redistributionRadius,
					  m_centroidStencilType));  // Dummy argumement
  }
}

const RefCountedPtr<EBIndexSpace>& phase_realm::get_ebis() {
  return m_ebis;
}

Vector<int>& phase_realm::getRefinementRatios() {
  return m_refinementRatios;
}

Vector<Real>& phase_realm::getDx() {
  return m_dx;
}

Vector<DisjointBoxLayout>& phase_realm::getGrids() {
  return m_grids;
}

const Vector<DisjointBoxLayout>& phase_realm::getGrids() const {
  return m_grids;
}

Vector<ProblemDomain>& phase_realm::getDomains() {
  return m_domains;
}

Vector<EBISLayout>& phase_realm::getEBISLayout() {
  return m_ebisl;
}

Vector<RefCountedPtr<EBLevelGrid> >& phase_realm::getEBLevelGrid() {
  return m_eblg;
}

Vector<RefCountedPtr<LayoutData<Vector<LayoutIndex> > > >& phase_realm::getNeighbors() {
  return m_neighbors;
}
Vector<RefCountedPtr<LayoutData<VoFIterator> > >& phase_realm::getVofIterator() {
  return m_vofiter;
}

IrregAmrStencil<CentroidInterpolationStencil>& phase_realm::getCentroidInterpolationStencils() {
  return *m_CentroidInterpolationStencil;
}

IrregAmrStencil<EbCentroidInterpolationStencil>& phase_realm::getEbCentroidInterpolationStencilStencils() {
  return *m_EbCentroidInterpolationStencil;
}

Vector<RefCountedPtr<LayoutData<BaseIVFAB<VoFStencil> > > >& phase_realm::get_gradsten(){
  if(!this->query_operator(s_eb_gradient)) MayDay::Abort("phase_realm::get_gradsten - operator not registered!");
  return m_gradsten;
}

// Throw errors if the operator does not exist

IrregAmrStencil<NonConservativeDivergenceStencil>& phase_realm::getNonConservativeDivergenceStencils() {
  if(!this->query_operator(s_eb_NonConservativeDivergenceStencil)) MayDay::Abort("phase_realm::get_non_cons_div_stencils - operator not registered!");
  
  return *m_NonConservativeDivergenceStencil;
}

EBAMRFAB& phase_realm::getLevelset() {
  if(!this->query_operator(s_levelset)) MayDay::Abort("phase_realm::getLevelset - operator not registered!");

  return m_levelset;
}

Vector<RefCountedPtr<ebcoarseaverage> >& phase_realm::getCoarseAverage() {
  if(!this->query_operator(s_eb_coar_ave)) MayDay::Abort("phase_realm::getCoarseAverage - operator not registered!");
  
  return m_coarave;
}

Vector<RefCountedPtr<EBGhostCloud> >& phase_realm::getGhostCloud() {
  if(!this->query_operator(s_eb_ghostcloud)) MayDay::Abort("phase_realm::getGhostCloud - operator not registered!");
  
  return m_ghostclouds;
}

Vector<RefCountedPtr<nwoebquadcfinterp> >& phase_realm::getNWOEBQuadCFInterp() {
  if(!this->query_operator(s_eb_quad_cfi)) MayDay::Abort("phase_realm::getNWOEBQuadCFInterp - operator not registered!");
  
  return m_quadcfi;
}

Vector<RefCountedPtr<EBQuadCFInterp> >& phase_realm::getEBQuadCFInterp() {
  if(!this->query_operator(s_eb_quad_cfi)) MayDay::Abort("phase_realm::getEBQuadCFInterp - operator not registered!");
  
  return m_old_quadcfi;
}

Vector<RefCountedPtr<AggEBPWLFillPatch> >& phase_realm::getFillPatch() {
  if(!this->query_operator(s_eb_fill_patch)) MayDay::Abort("phase_realm::getFillPatch - operator not registered!");
  
  return m_pwl_fillpatch;
}

Vector<RefCountedPtr<EBPWLFineInterp> >& phase_realm::getPwlInterpolator() {
  if(!this->query_operator(s_eb_pwl_interp)) MayDay::Abort("phase_realm::getPwlInterpolator - operator not registered!");
  
  return m_pwl_interp;
}

Vector<RefCountedPtr<EBMGInterp> >& phase_realm::getEBMGInterp() {
  if(!this->query_operator(s_eb_mg_interp)) MayDay::Abort("phase_realm::getEBMGInterp - operator not registered!");
  
  return m_ebmg_interp;
}

Vector<RefCountedPtr<EBFluxRegister> >&  phase_realm::getFluxRegister() {
  if(!this->query_operator(s_eb_flux_reg)) MayDay::Abort("phase_realm::getFluxRegister - operator not registered!");

  return m_flux_reg;
}

Vector<RefCountedPtr<EBLevelRedist> >&  phase_realm::getLevelRedist() {
  if(!this->query_operator(s_eb_redist)) MayDay::Abort("phase_realm::getLevelRedist - operator not registered!");

  return m_level_redist;
}

Vector<RefCountedPtr<EBCoarToFineRedist> >&  phase_realm::getCoarToFineRedist() {
  if(!this->query_operator(s_eb_redist)) MayDay::Abort("phase_realm::getCoarToFineRedist - operator not registered!");

  return m_coar_to_fine_redist;
}

Vector<RefCountedPtr<EBCoarToCoarRedist> >&  phase_realm::getCoarToCoarRedist() {
  if(!this->query_operator(s_eb_redist)) MayDay::Abort("phase_realm::get_coar_to_coar - operator not registered!");

  return m_coar_to_coar_redist;
}

Vector<RefCountedPtr<EBFineToCoarRedist> >&  phase_realm::getFineToCoarRedist() {
  if(!this->query_operator(s_eb_redist)) MayDay::Abort("phase_realm::getFineToCoarRedist - operator not registered!");

  return m_fine_to_coar_redist;
}

Vector<RefCountedPtr<Copier> >& phase_realm::getCopier() {
  if(!this->query_operator(s_eb_copier)) MayDay::Abort("phase_realm::getCopier - operator not registered!");

  return m_copier;
}

Vector<RefCountedPtr<Copier> >& phase_realm::getReverseCopier() {
  if(!this->query_operator(s_eb_copier)) MayDay::Abort("phase_realm::getReverseCopier - operator not registered!");
  
  return m_reverse_copier;
}
#include "CD_NamespaceFooter.H"
