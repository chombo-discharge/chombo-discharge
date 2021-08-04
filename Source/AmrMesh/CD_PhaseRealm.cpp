/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_PhaseRealm.cpp
  @brief  Implementation of CD_PhaseRealm.H
  @author Robert Marskar
*/

// Chombo includes
#include <EBArith.H>
#include <ParmParse.H>

// Our includes
#include <CD_PhaseRealm.H>
#include <CD_LoadBalancing.H>
#include <CD_EbFastFineToCoarRedist.H>
#include <CD_EbFastCoarToFineRedist.H>
#include <CD_EbFastCoarToCoarRedist.H>
#include <CD_EbFastFluxRegister.H>
#include <CD_EBMultigridInterpolator.H>
#include <CD_NamespaceHeader.H>

PhaseRealm::PhaseRealm(){
  // Default settings
  m_isDefined   = false;
  m_verbosity = -1;

  this->registerOperator(s_eb_gradient);
  this->registerOperator(s_eb_irreg_interp);

  // Adding this for debugging purposes. 
  ParmParse pp("PhaseRealm");
  pp.query("verbosity", m_verbosity);
}

PhaseRealm::~PhaseRealm(){

}

void PhaseRealm::define(const Vector<DisjointBoxLayout>&   a_grids,
			const Vector<ProblemDomain>&       a_domains,
			const Vector<int>&                 a_refRat,
			const Vector<Real>&                a_dx,
			const RealVect                     a_probLo,
			const int                          a_finestLevel,
			const int                          a_ebGhost,
			const int                          a_numGhost,
			const int                          a_lsfGhost,
			const int                          a_redistRad,
			const int                          a_mgInterpOrder,
			const int                          a_mgInterpRadius,
			const int                          a_mgInterpWeight,	      
			const IrregStencil::StencilType    a_centroidStencil,
			const IrregStencil::StencilType    a_ebStencil,
			const RefCountedPtr<BaseIF>&       a_baseif,
			const RefCountedPtr<EBIndexSpace>& a_ebis){
  m_grids                        = a_grids;
  m_domains                      = a_domains;
  m_refinementRatios             = a_refRat;
  m_dx                           = a_dx;
  m_probLo                       = a_probLo;
  m_finestLevel                  = a_finestLevel;
  m_numEbGhostsCells             = a_ebGhost;
  m_numGhostCells                = a_numGhost;
  m_numLsfGhostCells             = a_lsfGhost;
  m_redistributionRadius         = a_redistRad;
  m_multigridInterpolationOrder  = a_mgInterpOrder;
  m_multigridInterpolationRadius = a_mgInterpRadius;
  m_multigridInterpolationWeight = a_mgInterpWeight;
  m_centroidStencilType          = a_centroidStencil;
  m_ebCentroidStencilType        = a_ebStencil;
  m_baseif                       = a_baseif;
  m_ebis                         = a_ebis;
  
  if(!m_ebis.isNull()){
    m_isDefined = true;
  }

  m_hasEbCf = true;
}

void PhaseRealm::setGrids(const Vector<DisjointBoxLayout>& a_grids, const int a_finestLevel){
  CH_TIME("PhaseRealm::setGrids");
  if(m_verbosity > 5){
    pout() << "PhaseRealm::setGrids" << endl;
  }

  if(m_isDefined){
    m_grids = a_grids;
    m_finestLevel = a_finestLevel;
  }
}

void PhaseRealm::regridBase(const int a_lmin){
  CH_TIME("PhaseRealm::regridBase");
  if(m_verbosity > 5){
    pout() << "PhaseRealm::regridBase" << endl;
  }

  if(m_isDefined){
    this->defineEBLevelGrid(a_lmin);
    this->defineNeighbors(a_lmin);
    this->defineVofIterator(a_lmin);
  }
}

void PhaseRealm::regridOperators(const int a_lmin){
  CH_TIME("PhaseRealm::regridOperators");
  if(m_verbosity > 5){
    pout() << "PhaseRealm::regridOperators" << endl;
  }

  if(m_isDefined){
    this->defineEbCoarAve(a_lmin);                   
    this->defineEBQuadCFI(a_lmin);                   
    this->defineFillPatch(a_lmin);                   
    this->defineEBPWLInterp(a_lmin);                 
    this->defineMultigridInjection(a_lmin);          
    this->defineFluxReg(a_lmin, 1);                  
    this->defineRedistOper(a_lmin, 1);               
    this->defineGradSten(a_lmin);                    
    this->defineIrregSten();                         
    this->defineNonConsDivSten();                    
    this->defineCopier(a_lmin);                      
    this->defineGhostCloud(a_lmin);                  
    this->defineLevelSet(a_lmin, m_numLsfGhostCells);
    this->defineEBMultigrid(a_lmin);                 
  }
}

void PhaseRealm::registerOperator(const std::string a_operator){
  CH_TIME("PhaseRealm::registerOperator");
  if(m_verbosity > 5){
    pout() << "PhaseRealm::registerOperator" << endl;
  } 

  // These are the supported operators - issue an error if we ask for something that is not supported. 
  if(!(a_operator.compare(s_eb_coar_ave)     == 0 ||
       a_operator.compare(s_eb_quad_cfi)     == 0 ||
       a_operator.compare(s_eb_fill_patch)   == 0 ||
       a_operator.compare(s_eb_pwl_interp)   == 0 ||
       a_operator.compare(s_eb_flux_reg)     == 0 ||
       a_operator.compare(s_eb_redist)       == 0 ||
       a_operator.compare(s_noncons_div)     == 0 ||
       a_operator.compare(s_eb_copier)       == 0 ||
       a_operator.compare(s_eb_ghostcloud)   == 0 ||
       a_operator.compare(s_eb_gradient)     == 0 ||
       a_operator.compare(s_eb_irreg_interp) == 0 ||
       a_operator.compare(s_eb_mg_interp)    == 0 ||
       a_operator.compare(s_eb_multigrid)    == 0 ||
       a_operator.compare(s_levelset)        == 0 )){

    const std::string str = "PhaseRealm::registerOperator - unknown operator '" + a_operator + "' requested";
    MayDay::Abort(str.c_str());
  }

  if(!this->queryOperator(a_operator)){
    m_operatorMap.emplace(a_operator, true);
  }
}

bool PhaseRealm::queryOperator(const std::string a_operator) const {
  CH_TIME("PhaseRealm::queryOperator");
  if(m_verbosity > 5){
    pout() << "PhaseRealm::queryOperator" << endl;
  }

  //  return m_operatorMap[a_operator];

  bool ret = false;
  if(m_isDefined){
    ret = true;
    
    if(m_operatorMap.find(a_operator) == m_operatorMap.end()){
      ret = false;
    }
  }

  return ret;
}
  

void PhaseRealm::defineEBLevelGrid(const int a_lmin){
  CH_TIME("PhaseRealm::defineEBLevelGrid");
  if(m_verbosity > 2){
    pout() << "PhaseRealm::defineEBLevelGrid" << endl;
  }

  m_eblg.resize(1 + m_finestLevel);
  m_ebisl.resize(1 + m_finestLevel);

  for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++){
    m_eblg[lvl]  = RefCountedPtr<EBLevelGrid> (new EBLevelGrid(m_grids[lvl], m_domains[lvl], m_numEbGhostsCells, &(*m_ebis)));

    const bool hasCoar = lvl > 0;
    const bool hasFine = lvl < m_finestLevel;
    
    if(hasCoar) m_eblg[lvl]->setMaxCoarseningRatio(m_refinementRatios[lvl-1], &(*m_ebis));
    if(hasFine) m_eblg[lvl]->setMaxRefinementRatio(m_refinementRatios[lvl]);

    m_ebisl[lvl] = m_eblg[lvl]->getEBISL();
  }
}
void PhaseRealm::defineVofIterator(const int a_lmin){
  CH_TIME("PhaseRealm::defineVofIterator");
  if(m_verbosity > 2){
    pout() << "PhaseRealm::defineVofIterator" << endl;
  }

  m_vofIter.resize(1 + m_finestLevel);

  for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++){
      
    m_vofIter[lvl] = RefCountedPtr<LayoutData<VoFIterator> > (new LayoutData<VoFIterator> (m_grids[lvl]));
    
    for (DataIterator dit = m_grids[lvl].dataIterator(); dit.ok(); ++dit){
      VoFIterator& vofit = (*m_vofIter[lvl])[dit()];

      const Box& box         = m_grids[lvl].get(dit());
      const EBISBox& ebisbox = m_ebisl[lvl][dit()];
      const EBGraph& ebgraph = ebisbox.getEBGraph();
      const IntVectSet& irreg = ebisbox.getIrregIVS(box);

      vofit.define(irreg, ebgraph);
    }
  }
}

void PhaseRealm::defineNeighbors(const int a_lmin){
  CH_TIME("PhaseRealm::defineNeighbors");
  if(m_verbosity > 2){
    pout() << "PhaseRealm::defineNeighbors" << endl;
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

void PhaseRealm::defineLevelSet(const int a_lmin, const int a_numGhost){
  CH_TIME("PhaseRealm::defineLevelSet");
  if(m_verbosity > 2){
    pout() << "PhaseRealm::define_leveset" << endl;
  }

  const bool doThisOperator = this->queryOperator(s_levelset);

  m_levelset.resize(1 + m_finestLevel);

  if(doThisOperator){

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

void PhaseRealm::defineEbCoarAve(const int a_lmin){
  CH_TIME("PhaseRealm::defineEbCoarAve");
  if(m_verbosity > 2){
    pout() << "PhaseRealm::defineEbCoarAve" << endl;
  }

  const bool doThisOperator = this->queryOperator(s_eb_coar_ave);

  m_coarave.resize(1 + m_finestLevel);
  
  if(doThisOperator){
    
    const int comps = SpaceDim;

    for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++){

      const bool has_coar = lvl > 0;

      if(has_coar){
	m_coarave[lvl] = RefCountedPtr<EbCoarAve> (new EbCoarAve(m_grids[lvl],
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

void PhaseRealm::defineEBQuadCFI(const int a_lmin){
  CH_TIME("PhaseRealm::defineEBQuadCFI");
  if(m_verbosity > 2){
    pout() << "PhaseRealm::defineEBQuadCFI" << endl;
  }

  const bool doThisOperator = this->queryOperator(s_eb_quad_cfi);

  m_quadcfi.resize(1 + m_finestLevel);
  
  if(doThisOperator){
    for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++){

      const bool has_coar = lvl > 0;
      
      if(has_coar){
	const LayoutData<IntVectSet>& cfivs = *m_eblg[lvl]->getCFIVS();

	m_quadcfi[lvl] = RefCountedPtr<EBQuadCFInterp> (new EBQuadCFInterp(m_grids[lvl],
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

void PhaseRealm::defineEBMultigrid(const int a_lmin){
  CH_TIME("PhaseRealm::defineEBMultigrid");
  if(m_verbosity > 2){
    pout() << "PhaseRealm::defineEBMultigrid" << endl;
  }

  const bool doThisOperator = this->queryOperator(s_eb_multigrid);

  m_multigridInterpolator.resize(1 + m_finestLevel);
  
  if(doThisOperator){
    for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++){

      const bool hasCoar = lvl > 0;
      
      if(hasCoar){
	const EBLevelGrid& eblgFine         = *m_eblg[lvl  ];
	const EBLevelGrid& eblgCoar         = *m_eblg[lvl-1];
	
	const LayoutData<IntVectSet>& cfivs = *eblgFine.getCFIVS();

	m_multigridInterpolator[lvl] = RefCountedPtr<EBMultigridInterpolator>(new EBMultigridInterpolator(eblgFine,
													  eblgCoar,
													  Location::Cell::Center,
													  m_numGhostCells*IntVect::Unit,
													  m_refinementRatios[lvl-1],
													  1,
													  m_multigridInterpolationRadius,
													  m_multigridInterpolationOrder,
													  m_multigridInterpolationWeight));
      }
    }
  }
}

void PhaseRealm::defineFillPatch(const int a_lmin){
  CH_TIME("PhaseRealm::defineFillPatch");
  if(m_verbosity > 2){
    pout() << "PhaseRealm::defineFillPatch" << endl;
  }

  const bool doThisOperator = this->queryOperator(s_eb_fill_patch);

  m_pwlFillPatch.resize(1 + m_finestLevel);
  
  if(doThisOperator){
    
    const int comps     = SpaceDim;

    // Should these be input somehow?
    const int radius    = 1;
    const IntVect ghost = m_numGhostCells*IntVect::Unit;


    for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++){

      const bool has_coar = lvl > 0;

      if(has_coar){
	const LayoutData<IntVectSet>& cfivs = *m_eblg[lvl]->getCFIVS();
	m_pwlFillPatch[lvl] = RefCountedPtr<AggEBPWLFillPatch> (new AggEBPWLFillPatch(m_grids[lvl],
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


void PhaseRealm::defineEBPWLInterp(const int a_lmin){
  CH_TIME("PhaseRealm::defineEBPWLInterp");
  if(m_verbosity > 2){
    pout() << "PhaseRealm::defineEBPWLInterp" << endl;
  }

  const bool doThisOperator = this->queryOperator(s_eb_pwl_interp);

  m_pwlInterpGhosts.resize(1 + m_finestLevel);

  if(doThisOperator){

	
    const int comps     = SpaceDim;

    // Should these be input somehow?
    const int radius    = 1;
    const IntVect ghost = m_numGhostCells*IntVect::Unit;

    for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++){

      const bool has_coar = lvl > 0;

      if(has_coar){
	m_pwlInterpGhosts[lvl] = RefCountedPtr<EBPWLFineInterp> (new EBPWLFineInterp(m_grids[lvl],
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

void PhaseRealm::defineMultigridInjection(const int a_lmin){
  CH_TIME("PhaseRealm::defineMultigridInjection");
  if(m_verbosity > 2){
    pout() << "PhaseRealm::defineMultigridInjection" << endl;
  }

  const bool doThisOperator = this->queryOperator(s_eb_mg_interp);

  m_mgInjection.resize(1 + m_finestLevel);

  if(doThisOperator){

    
    const int ncomps    = 1;

    // Should these be input somehow?
    const int radius    = 1;
    const IntVect ghost = m_numGhostCells*IntVect::Unit;

    for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++){

      const bool has_coar = lvl > 0;

      if(has_coar){

	m_mgInjection[lvl] = RefCountedPtr<EBMGInterp> (new EBMGInterp(m_grids[lvl],
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

void PhaseRealm::defineFluxReg(const int a_lmin, const int a_regsize){
  CH_TIME("PhaseRealm::defineFluxReg");
  if(m_verbosity > 2){
    pout() << "PhaseRealm::defineFluxReg" << endl;
  }

  const bool doThisOperator = this->queryOperator(s_eb_flux_reg);

  m_fluxReg.resize(1 + m_finestLevel);

  if(doThisOperator){
    
    const int comps = a_regsize;
    
    for (int lvl = Max(0,a_lmin-1); lvl <= m_finestLevel; lvl++){

      const bool has_fine = lvl < m_finestLevel;


      if(has_fine){
#if 1 // Our code
	m_fluxReg[lvl] = RefCountedPtr<EBFluxRegister> (new EbFastFluxRegister(m_grids[lvl+1],
										m_grids[lvl],
										m_ebisl[lvl+1],
										m_ebisl[lvl],
										m_domains[lvl].domainBox(),
										m_refinementRatios[lvl],
										comps,
										&(*m_ebis)));
#else // Chombo code. Leaving this in place in case something breaks. 
	m_fluxReg[lvl] = RefCountedPtr<EBFluxRegister> (new EBFluxRegister(m_grids[lvl+1],
									    m_grids[lvl],
									    m_ebisl[lvl+1],
									    m_ebisl[lvl],
									    m_domains[lvl].domainBox(),
									    m_refinementRatios[lvl],
									    comps,
									    &(*m_ebis)));
#endif
      }
    }
  }
}

void PhaseRealm::defineRedistOper(const int a_lmin, const int a_regsize){
  CH_TIME("PhaseRealm::defineRedistOper");
  if(m_verbosity > 2){
    pout() << "PhaseRealm::defineRedistOper" << endl;
  }

  // TLDR: All these operators either do stuff on an AMR level, or between a coarse and a fine level. The entries
  //       live on these levels:
  //
  //       Oper                        Level
  //       EBFineToCoar [l,  l-1]    l
  //       EBCoarToFine [l-1,l  ] 
  //       when level a_lmin changed we need to update fine-to-coar

  const bool doThisOperator = this->queryOperator(s_eb_redist);

  m_levelRedist.resize(1 + m_finestLevel);
  m_fineToCoarRedist.resize(1 + m_finestLevel);
  m_coarToCoarRedist.resize(1 + m_finestLevel);
  m_coarToFineRedist.resize(1 + m_finestLevel);

  if(doThisOperator){
    
    const int comps = a_regsize;

    for (int lvl = Max(0, a_lmin-1); lvl <= m_finestLevel; lvl++){

      const bool has_coar = lvl > 0;
      const bool has_fine = lvl < m_finestLevel;

	      
      if(lvl >= a_lmin){
	m_levelRedist[lvl] = RefCountedPtr<EBLevelRedist> (new EBLevelRedist(m_grids[lvl],
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

	    auto redist = RefCountedPtr<EbFastFineToCoarRedist> (new EbFastFineToCoarRedist());
	    redist->define(*m_eblg[lvl],
			   *m_eblg[lvl-1],
			   *m_neighbors[lvl],
			   *m_neighbors[lvl-1],
			   m_refinementRatios[lvl-1],
			   comps,
			   m_redistributionRadius);
	    m_fineToCoarRedist[lvl] = RefCountedPtr<EBFineToCoarRedist> (redist);
	    
	  }
	}

	if(has_fine){
	  // TLDR: The coar-to-fine redistribution operator transfers from the coarse level and to the fine level and
	  //       therefore lives on the coarse level. Since a_lmin is the coarsest level that changed, we need to update
	  //       if lvl >= a_lmin-1
	  if(lvl >= a_lmin-1){

	    auto c2f_redist = RefCountedPtr<EbFastCoarToFineRedist> (new EbFastCoarToFineRedist());
	    c2f_redist->define(*m_eblg[lvl+1],
			       *m_eblg[lvl],
			       *m_neighbors[lvl+1],
			       *m_neighbors[lvl],
			       m_refinementRatios[lvl],
			       comps,
			       m_redistributionRadius);
	    m_coarToFineRedist[lvl] = RefCountedPtr<EBCoarToFineRedist> (c2f_redist);


	    auto c2c_redist = RefCountedPtr<EbFastCoarToCoarRedist> (new EbFastCoarToCoarRedist());
	    c2c_redist->define(*m_eblg[lvl+1],
			       *m_eblg[lvl],
			       *m_neighbors[lvl+1],
			       *m_neighbors[lvl],
			       m_refinementRatios[lvl],
			       comps,
			       m_redistributionRadius);
	    m_coarToCoarRedist[lvl] = RefCountedPtr<EBCoarToCoarRedist> (c2c_redist);
	  }
	}
      }
    }
  }
}

void PhaseRealm::defineGradSten(const int a_lmin){
  CH_TIME("PhaseRealm::defineGradSten");
  if(m_verbosity > 2){
    pout() << "PhaseRealm::defineGradSten" << endl;
  }

  const bool doThisOperator = this->queryOperator(s_eb_gradient);

  m_gradsten.resize(1 + m_finestLevel);

  if(doThisOperator){
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

void PhaseRealm::defineCopier(const int a_lmin){
  CH_TIME("PhaseRealm::defineCopier");
  if(m_verbosity > 3){
    pout() << "PhaseRealm::defineCopier" << endl;
  }

  const bool doThisOperator = this->queryOperator(s_eb_copier);

  m_copier.resize(1 + m_finestLevel);
  m_reverseCopier.resize(1 + m_finestLevel);

  if(doThisOperator){
    
    for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++){
      m_copier[lvl] = RefCountedPtr<Copier> (new Copier(m_grids[lvl],
							m_grids[lvl],
							m_domains[lvl],
							m_numGhostCells*IntVect::Unit,
							true));
      
      m_reverseCopier[lvl] = RefCountedPtr<Copier> (new Copier(m_grids[lvl],
								m_grids[lvl],
								m_domains[lvl],
								m_numGhostCells*IntVect::Unit,
								true));
      m_reverseCopier[lvl]->reverse();
    }
  }
}


void PhaseRealm::defineGhostCloud(const int a_lmin){
  CH_TIME("PhaseRealm::defineGhostCloud");
  if(m_verbosity > 3){
    pout() << "PhaseRealm::defineGhostCloud" << endl;
  }

  const bool doThisOperator = this->queryOperator(s_eb_copier);

  m_ghostclouds.resize(1 + m_finestLevel);

  if(doThisOperator){
    for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++){
      const bool has_coar = lvl > 0;

      if(has_coar){
	m_ghostclouds[lvl] = RefCountedPtr<EbGhostCloud> (new EbGhostCloud(m_grids[lvl-1],
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


void PhaseRealm::defineIrregSten(){
  CH_TIME("PhaseRealm::defineIrregSten");
  if(m_verbosity > 2){
    pout() << "PhaseRealm::defineIrregSten" << endl;
  }

  const bool doThisOperator = this->queryOperator(s_eb_irreg_interp);

  if(doThisOperator){
    
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

void PhaseRealm::defineNonConsDivSten(){
  CH_TIME("PhaseRealm::defineNonConsDivSten");
  if(m_verbosity > 2){
    pout() << "PhaseRealm::defineNonConsDivSten" << endl;
  }

  const bool doThisOperator = this->queryOperator(s_noncons_div);

  if(doThisOperator){
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

const RefCountedPtr<EBIndexSpace>& PhaseRealm::getEBIndexSpace() const {
  return m_ebis;
}

const Vector<int>& PhaseRealm::getRefinementRatios() const {
  return m_refinementRatios;
}

const Vector<Real>& PhaseRealm::getDx() const {
  return m_dx;
}

const Vector<DisjointBoxLayout>& PhaseRealm::getGrids() const {
  return m_grids;
}

const Vector<ProblemDomain>& PhaseRealm::getDomains() const {
  return m_domains;
}

const Vector<EBISLayout>& PhaseRealm::getEBISLayout() const {
  return m_ebisl;
}

const Vector<RefCountedPtr<EBLevelGrid> >& PhaseRealm::getEBLevelGrid() const {
  return m_eblg;
}

const Vector<RefCountedPtr<LayoutData<Vector<LayoutIndex> > > >& PhaseRealm::getNeighbors() const {
  return m_neighbors;
}
Vector<RefCountedPtr<LayoutData<VoFIterator> > >& PhaseRealm::getVofIterator() const {
  return m_vofIter;
}

const IrregAmrStencil<CentroidInterpolationStencil>& PhaseRealm::getCentroidInterpolationStencils() const {
  return *m_CentroidInterpolationStencil;
}

const IrregAmrStencil<EbCentroidInterpolationStencil>& PhaseRealm::getEbCentroidInterpolationStencilStencils() const {
  return *m_EbCentroidInterpolationStencil;
}

const Vector<RefCountedPtr<LayoutData<BaseIVFAB<VoFStencil> > > >& PhaseRealm::getGradientStencils() const {
  if(!this->queryOperator(s_eb_gradient)) MayDay::Abort("PhaseRealm::getGradientStencils - operator not registered!");
  
  return m_gradsten;
}

const IrregAmrStencil<NonConservativeDivergenceStencil>& PhaseRealm::getNonConservativeDivergenceStencils() const {
  if(!this->queryOperator(s_noncons_div)) MayDay::Abort("PhaseRealm::get_non_cons_div_stencils - operator not registered!");
  
  return *m_NonConservativeDivergenceStencil;
}

Vector<RefCountedPtr<EbCoarAve> >& PhaseRealm::getCoarseAverage() const {
  if(!this->queryOperator(s_eb_coar_ave)) MayDay::Abort("PhaseRealm::getCoarseAverage - operator not registered!");
  
  return m_coarave;
}

Vector<RefCountedPtr<EbGhostCloud> >& PhaseRealm::getGhostCloud() const {
  if(!this->queryOperator(s_eb_ghostcloud)) MayDay::Abort("PhaseRealm::getGhostCloud - operator not registered!");
  
  return m_ghostclouds;
}

Vector<RefCountedPtr<EBQuadCFInterp> >& PhaseRealm::getEBQuadCFInterp() const {
  if(!this->queryOperator(s_eb_quad_cfi)) MayDay::Abort("PhaseRealm::getEBQuadCFInterp - operator not registered!");
  
  return m_quadcfi;
}

Vector<RefCountedPtr<EBMultigridInterpolator> >& PhaseRealm::getMultigridInterpolator() const {
  if(!this->queryOperator(s_eb_multigrid)) MayDay::Abort("PhaseRealm::getEBMultigridInterpolator - operator not registered!");
  
  return m_multigridInterpolator;
}

Vector<RefCountedPtr<AggEBPWLFillPatch> >& PhaseRealm::getFillPatch() const {
  if(!this->queryOperator(s_eb_fill_patch)) MayDay::Abort("PhaseRealm::getFillPatch - operator not registered!");
  
  return m_pwlFillPatch;
}

Vector<RefCountedPtr<EBPWLFineInterp> >& PhaseRealm::getPwlInterpolator() const {
  if(!this->queryOperator(s_eb_pwl_interp)) MayDay::Abort("PhaseRealm::getPwlInterpolator - operator not registered!");
  
  return m_pwlInterpGhosts;
}

Vector<RefCountedPtr<EBMGInterp> >& PhaseRealm::getEBMGInterp() const {
  if(!this->queryOperator(s_eb_mg_interp)) MayDay::Abort("PhaseRealm::getEBMGInterp - operator not registered!");
  
  return m_mgInjection;
}

Vector<RefCountedPtr<EBFluxRegister> >&  PhaseRealm::getFluxRegister() const {
  if(!this->queryOperator(s_eb_flux_reg)) MayDay::Abort("PhaseRealm::getFluxRegister - operator not registered!");

  return m_fluxReg;
}

Vector<RefCountedPtr<EBLevelRedist> >&  PhaseRealm::getLevelRedist() const {
  if(!this->queryOperator(s_eb_redist)) MayDay::Abort("PhaseRealm::getLevelRedist - operator not registered!");

  return m_levelRedist;
}

Vector<RefCountedPtr<EBCoarToFineRedist> >&  PhaseRealm::getCoarToFineRedist() const {
  if(!this->queryOperator(s_eb_redist)) MayDay::Abort("PhaseRealm::getCoarToFineRedist - operator not registered!");

  return m_coarToFineRedist;
}

Vector<RefCountedPtr<EBCoarToCoarRedist> >&  PhaseRealm::getCoarToCoarRedist() const {
  if(!this->queryOperator(s_eb_redist)) MayDay::Abort("PhaseRealm::get_coar_to_coar - operator not registered!");

  return m_coarToCoarRedist;
}

Vector<RefCountedPtr<EBFineToCoarRedist> >&  PhaseRealm::getFineToCoarRedist() const {
  if(!this->queryOperator(s_eb_redist)) MayDay::Abort("PhaseRealm::getFineToCoarRedist - operator not registered!");

  return m_fineToCoarRedist;
}

Vector<RefCountedPtr<Copier> >& PhaseRealm::getCopier() const {
  if(!this->queryOperator(s_eb_copier)) MayDay::Abort("PhaseRealm::getCopier - operator not registered!");

  return m_copier;
}

Vector<RefCountedPtr<Copier> >& PhaseRealm::getReverseCopier() const {
  if(!this->queryOperator(s_eb_copier)) MayDay::Abort("PhaseRealm::getReverseCopier - operator not registered!");
  
  return m_reverseCopier;
}

const EBAMRFAB& PhaseRealm::getLevelset() const {
  if(!this->queryOperator(s_levelset)) MayDay::Abort("PhaseRealm::getLevelset - operator not registered!");

  return m_levelset;
}

#include <CD_NamespaceFooter.H>
