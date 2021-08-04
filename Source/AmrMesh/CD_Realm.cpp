/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_Realm.cpp
  @brief  Implementation of CD_Realm.H
  @author Robert Marskar
*/

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_Realm.H>
#include <CD_NamespaceHeader.H>

const std::string Realm::Primal = "primal";

template <class T>
void Realm::copy(Vector<RefCountedPtr<LevelData<T> > >& a_dst, const Vector<RefCountedPtr<LevelData<T> > >& a_src){
  for (int lvl = 0; lvl < a_dst.size(); lvl++){
    if(a_src[lvl] != NULL && a_dst[lvl] != NULL){
      a_src[lvl]->copyTo(a_dst[lvl]);
    }
  }
}

Realm::Realm(){
  m_defined = false;
  m_verbosity = -1;

  ParmParse pp("PhaseRealm");

  pp.query("verbosity", m_verbosity);

  // Just empty points until define() is called
  m_realms.emplace(phase::gas,   RefCountedPtr<PhaseRealm> (new PhaseRealm()));
  m_realms.emplace(phase::solid, RefCountedPtr<PhaseRealm> (new PhaseRealm()));
}

Realm::~Realm(){

}

void Realm::define(const Vector<DisjointBoxLayout>& a_grids,
		   const Vector<ProblemDomain>& a_domains,
		   const Vector<int>& a_ref_rat,
		   const Vector<Real>& a_dx,
		   const RealVect a_probLo,
		   const int a_finestLevel,
		   const int a_ebghost,
		   const int a_num_ghost,
		   const int a_lsf_ghost,
		   const int a_redist_rad,
		   const int a_mgInterpOrder,
		   const int a_mgInterpRadius,
		   const int a_mgInterpWeight,
		   const IrregStencil::StencilType a_centroid_stencil,
		   const IrregStencil::StencilType a_eb_stencil,
		   const std::map<phase::which_phase, RefCountedPtr<BaseIF> > a_baseif,
		   const RefCountedPtr<MultiFluidIndexSpace>& a_multiFluidIndexSpace){
  CH_TIME("Realm::define");
  if(m_verbosity > 5){
    pout() << "Realm::define" << endl;
  }



  m_refinementRatios = a_ref_rat;
  m_dx = a_dx;
  m_grids = a_grids;
  m_domains = a_domains;
  m_baseif = a_baseif;

  m_probLo = a_probLo;
  m_multifluidIndexSpace = a_multiFluidIndexSpace;
  m_finestLevel = a_finestLevel;
  
  const RefCountedPtr<EBIndexSpace>& ebis_gas = m_multifluidIndexSpace->getEBIndexSpace(phase::gas);
  const RefCountedPtr<EBIndexSpace>& ebis_sol = m_multifluidIndexSpace->getEBIndexSpace(phase::solid);

  m_realms[phase::gas]->define(a_grids,
			       a_domains,
			       a_ref_rat,
			       a_dx,
			       m_probLo,
			       a_finestLevel,
			       a_ebghost,
			       a_num_ghost,
			       a_lsf_ghost,
			       a_redist_rad,
			       a_mgInterpOrder,
			       a_mgInterpRadius,
			       a_mgInterpWeight,
			       a_centroid_stencil,
			       a_eb_stencil,
			       m_baseif.at(phase::gas),
			       ebis_gas);

  m_realms[phase::solid]->define(a_grids,
				 a_domains,
				 a_ref_rat,
				 a_dx,
				 m_probLo,
				 a_finestLevel,
				 a_ebghost,
				 a_num_ghost,
				 a_lsf_ghost,
				 a_redist_rad,
				 a_mgInterpOrder,
				 a_mgInterpRadius,
				 a_mgInterpWeight,
				 a_centroid_stencil,
				 a_eb_stencil,
				 m_baseif.at(phase::solid),
				 ebis_sol);
}

void Realm::setGrids(const Vector<DisjointBoxLayout>& a_grids, const int a_finestLevel){
  CH_TIME("Realm::setGrids");
  if(m_verbosity > 5){
    pout() << "Realm::setGrids" << endl;
  }

  for (auto& r : m_realms){
    r.second->setGrids(a_grids, a_finestLevel);
  }
}

void Realm::regridBase(const int a_lmin){
  CH_TIME("Realm::regridBase");
  if(m_verbosity > 5){
    pout() << "Realm::regridBase" << endl;
  }
  
  for (auto& r : m_realms){
    r.second->regridBase(a_lmin);
  }
  this->define_mflevelgrid(a_lmin);
}

void Realm::regridOperators(const int a_lmin){
  CH_TIME("Realm::regridOperators");
  if(m_verbosity > 5){
    pout() << "Realm::regridOperators" << endl;
  }

  for (auto& r : m_realms){
    r.second->regridOperators(a_lmin);
  }
  this->define_masks(a_lmin);
}

void Realm::define_masks(const int a_lmin){
  CH_TIME("Realm::define_masks");
  if(m_verbosity > 5){
    pout() << "Realm::define_masks" << endl;
  }

  // Regrid all masks
  this->define_halo_masks(a_lmin);
}

void Realm::define_mflevelgrid(const int a_lmin){
  CH_TIME("Realm::define_mflevelgrid");
  if(m_verbosity > 5){
    pout() << "Realm::define_mflevelgrid" << endl;
  }

  m_mflg.resize(1 + m_finestLevel);

  PhaseRealm& gas = this->getRealm(phase::gas);
  PhaseRealm& sol = this->getRealm(phase::solid);

  const RefCountedPtr<EBIndexSpace>& ebis_gas = gas.getEBIndexSpace();
  const RefCountedPtr<EBIndexSpace>& ebis_sol = sol.getEBIndexSpace();

  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    Vector<EBLevelGrid> eblgs;

    if(!ebis_gas.isNull()) eblgs.push_back(*(gas.getEBLevelGrid()[lvl]));
    if(!ebis_sol.isNull()) eblgs.push_back(*(sol.getEBLevelGrid()[lvl]));

    m_mflg[lvl] = RefCountedPtr<MFLevelGrid> (new MFLevelGrid(m_multifluidIndexSpace, eblgs));
  }
}

void Realm::define_halo_masks(const int a_lmin){
  CH_TIME("Realm::define_halo_masks");
  if(m_verbosity > 5){
    pout() << "Realm::define_halo_masks" << endl;
  }

  // Loop through all masks and do something about the halo masks only. 
  for (auto& m : m_masks){

    // Get mask identifier and buffer. 
    const std::string which_mask = m.first.first;
    const int buffer             = m.first.second;

    if(which_mask == s_particle_halo){
      if(buffer < 0) MayDay::Abort("Realm::define_halo_masks -- cannot have buffer < 0!");

      AMRMask& mask = m.second;

      mask.resize(1 + m_finestLevel);

      for (int lvl = 0; lvl < m_finestLevel; lvl++){
	const DisjointBoxLayout& gridsCoar = m_grids[lvl];
	const DisjointBoxLayout& gridsFine = m_grids[lvl+1];

	const ProblemDomain& domainCoar = m_domains[lvl];
	const ProblemDomain& domainFine = m_domains[lvl+1];

	const int ncomp = 1;

	mask[lvl] = RefCountedPtr<LevelData<BaseFab<bool> > >(new LevelData<BaseFab<bool> >(gridsCoar, ncomp, IntVect::Zero));

	this->define_halo_mask(*mask[lvl],
			       domainCoar,
			       domainFine,
			       gridsCoar,
			       gridsFine,
			       buffer,
			       m_refinementRatios[lvl]);
      }
    }
  }
}

void Realm::define_halo_mask(LevelData<BaseFab<bool> >& a_coarMask,
			     const ProblemDomain&       a_domainCoar,
			     const ProblemDomain&       a_domainFine,
			     const DisjointBoxLayout&   a_gridsCoar,
			     const DisjointBoxLayout&   a_gridsFine,
			     const int                  a_buffer,
			     const int                  a_ref_rat){
  CH_TIME("Realm::define_halo_masks");
  if(m_verbosity > 5){
    pout() << "Realm::define_halo_masks" << endl;
  }

  // TLDR: This function creates a GLOBAL view of the CFIVS, a_buffer cells wide on the coarse side of the refinement boundary. If this is not performant,
  //       there are probably better ways to do this. 
  
  const int comp  = 0;
  const int ncomp = 1;

  // First, reset the mask.
  for (DataIterator dit(a_gridsCoar); dit.ok(); ++dit){
    a_coarMask[dit()].setVal(false);
  }

  // Ok, we need a particle halo. 
  if(a_buffer > 0){

    // Coarsen the fine grid and make a mask on the coarsened fine grid
    DisjointBoxLayout gridsCoFi;
    coarsen(gridsCoFi, a_gridsFine, a_ref_rat);

    IntVectSet halo;                 // Global CFIVS - should find a better way to do this.
    NeighborIterator nit(gridsCoFi); // Neighbor iterator

    // Go through the cofi grid and set the halo to true
    for (DataIterator dit(gridsCoFi); dit.ok(); ++dit){
      const Box coFiBox = gridsCoFi[dit()];

      // Make IntVect set consisting of only ghost cells (a_buffer) for each box. 
      Box grownBox = grow(coFiBox, a_buffer);
      grownBox    &= a_domainCoar;

      // Subtract non-ghosted box. 
      IntVectSet myHalo(grownBox);
      myHalo -= coFiBox;

      // Subtract non-ghosted neighboring boxes. 
      for (nit.begin(dit()); nit.ok(); ++nit){
	const Box neighborBox = gridsCoFi[nit()];
	myHalo -= neighborBox;
      }

      // Add to halo. 
      halo |= myHalo;
    }

    // TLDR: In the above, we found the coarse-grid cells surrounding the fine level, viewed from the fine grids.
    //       Below, we create that view from the coarse grid. We use a BoxLayoutData<FArrayBox> on the coarsened fine grid,
    //       whose "ghost cells" can added to the _actual_ coarse grid. We then loop through those cells and set the mask. 

    // Get some box definitions. 
    Vector<Box> coFiBoxes = gridsCoFi.boxArray();
    Vector<int> coFiProcs = gridsCoFi.procIDs();

    Vector<Box> coarBoxes = a_gridsCoar.boxArray();
    Vector<Box> coarPrioc = a_gridsCoar.boxArray();

    // Grow the coarsened fine boxes
    for (int i = 0; i < coFiBoxes.size(); i++){
      coFiBoxes[i].grow(a_buffer);
      coFiBoxes[i] &= a_domainCoar;
    }

    // Define boxlayouts
    BoxLayout blCoFi(coFiBoxes, coFiProcs);
    BoxLayout blCoar(a_gridsCoar.boxArray(), a_gridsCoar.procIDs());

    BoxLayoutData<FArrayBox> coFiMask(blCoFi, ncomp);
    BoxLayoutData<FArrayBox> coarMask(blCoar, ncomp);
  
    // Reset masks
    for (DataIterator dit(blCoFi); dit.ok(); ++dit){
      coFiMask[dit()].setVal(0.0);
    }
    for (DataIterator dit(a_gridsCoar); dit.ok(); ++dit){
      coarMask[dit()].setVal(0.0);
    }

    // Run through the halo and set the halo cells to 1 on the coarsened fine grids (these are essentially "ghost cells")
    for (DataIterator dit(blCoFi); dit.ok(); ++dit){
      IntVectSet curHalo = halo & blCoFi[dit()];
      for (IVSIterator ivsit(curHalo); ivsit.ok(); ++ivsit){
	coFiMask[dit()](ivsit(), comp) = 1.0;
      }
    }

    // Add the result to the coarse grid. 
    const Interval interv(0,0);
    coFiMask.addTo(interv, coarMask, interv, a_domainCoar.domainBox());

    // Run through the grids and make the boolean mask
    for (DataIterator dit(a_gridsCoar); dit.ok(); ++dit){
      const Box box = a_gridsCoar[dit()];

      BaseFab<bool>& mask     = a_coarMask[dit()];
      const FArrayBox& blMask = coarMask[dit()];

      for (BoxIterator bit(box); bit.ok(); ++bit){
	const IntVect iv = bit();
	if(blMask(iv, comp) > 0.0) mask(iv, comp) = true;
      }
    }
  }
}

void Realm::registerOperator(const std::string a_operator, const phase::which_phase a_phase){
  CH_TIME("Realm::registerOperator(operator, phase)");
  if(m_verbosity > 5){
    pout() << "Realm::registerOperator(operator, phase)" << endl;
  }
  
  m_realms[a_phase]->registerOperator(a_operator);
}

bool Realm::queryOperator(const std::string a_operator, const phase::which_phase a_phase) {
  CH_TIME("Realm::queryOperator");
  if(m_verbosity > 5){
    pout() << "Realm::queryOperator" << endl;
  }
  
  return m_realms[a_phase]->queryOperator(a_operator);
}

void Realm::registerMask(const std::string a_mask, const int a_buffer){
  CH_TIME("Realm::registerMask(mask, buffer)");
  if(m_verbosity > 5){
    pout() << "Realm::registerMask(mask, buffer)" << endl;
  }

  m_masks.emplace(std::pair<string, int>(a_mask, a_buffer), AMRMask());
}

bool Realm::queryMask(const std::string a_mask, const int a_buffer) const {
  CH_TIME("Realm::queryMask(mask, buffer)");
  if(m_verbosity > 5){
    pout() << "Realm::queryMask(mask, buffer)" << endl;
  }

  bool ret;

  if(m_masks.count(std::pair<string, int>(a_mask, a_buffer) ) > 0){
    ret = true;
  }
  else{
    ret = false;
  }

  return ret;
}

PhaseRealm& Realm::getRealm(const phase::which_phase a_phase){
  return *m_realms[a_phase];
}

Vector<int>& Realm::getRefinementRatios() {
  return m_refinementRatios;
}

Vector<Real>& Realm::getDx() {
  return m_dx;
}

const Vector<Real>& Realm::getDx() const {
  return m_dx;
}

Vector<DisjointBoxLayout>& Realm::getGrids() {
  return m_grids;
}

const Vector<DisjointBoxLayout>& Realm::getGrids() const {
  return m_grids;
}

Vector<ProblemDomain>& Realm::getDomains() {
  return m_domains;
}

const Vector<ProblemDomain>& Realm::getDomains() const {
  return m_domains;
}

Vector<RefCountedPtr<MFLevelGrid> >& Realm::getMFLevelGrid(){
  return m_mflg;
}

const RefCountedPtr<EBIndexSpace>& Realm::getEBIndexSpace(const phase::which_phase a_phase) const {
  return m_realms[a_phase]->getEBIndexSpace();
}

const Vector<EBISLayout>& Realm::getEBISLayout(const phase::which_phase a_phase) const {
  return m_realms[a_phase]->getEBISLayout();
}

const Vector<RefCountedPtr<EBLevelGrid> >& Realm::getEBLevelGrid(const phase::which_phase a_phase) const {
  return m_realms[a_phase]->getEBLevelGrid();
}

const Vector<RefCountedPtr<LayoutData<Vector<LayoutIndex> > > >& Realm::getNeighbors(const phase::which_phase a_phase) const {
  return m_realms[a_phase]->getNeighbors();
}

Vector<RefCountedPtr<LayoutData<VoFIterator> > >& Realm::getVofIterator(const phase::which_phase a_phase) const {
  return m_realms[a_phase]->getVofIterator();
}

const IrregAmrStencil<CentroidInterpolationStencil>& Realm::getCentroidInterpolationStencils(const phase::which_phase a_phase) const {
  return m_realms[a_phase]->getCentroidInterpolationStencils();
}

const IrregAmrStencil<EbCentroidInterpolationStencil>& Realm::getEbCentroidInterpolationStencilStencils(const phase::which_phase a_phase) const {
  return m_realms[a_phase]->getEbCentroidInterpolationStencilStencils();
}

const IrregAmrStencil<NonConservativeDivergenceStencil>& Realm::getNonConservativeDivergenceStencils(const phase::which_phase a_phase) const {
  return m_realms[a_phase]->getNonConservativeDivergenceStencils();
}

const Vector<RefCountedPtr<LayoutData<BaseIVFAB<VoFStencil> > > >& Realm::getGradientStencils(const phase::which_phase a_phase) const{
  return m_realms[a_phase]->getGradientStencils();
}

Vector<RefCountedPtr<EbCoarAve> >& Realm::getCoarseAverage(const phase::which_phase a_phase){
  return m_realms[a_phase]->getCoarseAverage();
}

Vector<RefCountedPtr<EbGhostCloud> >& Realm::getGhostCloud(const phase::which_phase a_phase){
  return m_realms[a_phase]->getGhostCloud();
}

Vector<RefCountedPtr<EBQuadCFInterp> >& Realm::getEBQuadCFInterp(const phase::which_phase a_phase){
  return m_realms[a_phase]->getEBQuadCFInterp();
}

Vector<RefCountedPtr<EBMultigridInterpolator> >& Realm::getMultigridInterpolator(const phase::which_phase a_phase){
  return m_realms[a_phase]->getMultigridInterpolator();
}

Vector<RefCountedPtr<AggEBPWLFillPatch> >& Realm::getFillPatch(const phase::which_phase a_phase){
  return m_realms[a_phase]->getFillPatch();
}

Vector<RefCountedPtr<EBPWLFineInterp> >& Realm::getPwlInterpolator(const phase::which_phase a_phase){
  return m_realms[a_phase]->getPwlInterpolator();
}

Vector<RefCountedPtr<EBMGInterp> >& Realm::getEBMGInterp(const phase::which_phase a_phase){
  return m_realms[a_phase]->getEBMGInterp();
}

Vector<RefCountedPtr<EBFluxRegister> >&  Realm::getFluxRegister(const phase::which_phase a_phase){
  return m_realms[a_phase]->getFluxRegister();
}

Vector<RefCountedPtr<EBLevelRedist> >&  Realm::getLevelRedist(const phase::which_phase a_phase) {
  return m_realms[a_phase]->getLevelRedist();
}

Vector<RefCountedPtr<EBCoarToFineRedist> >&  Realm::getCoarToFineRedist(const phase::which_phase a_phase){
  return m_realms[a_phase]->getCoarToFineRedist();
}

Vector<RefCountedPtr<EBCoarToCoarRedist> >&  Realm::getCoarToCoarRedist(const phase::which_phase a_phase){
  return m_realms[a_phase]->getCoarToCoarRedist();
}

Vector<RefCountedPtr<EBFineToCoarRedist> >&  Realm::getFineToCoarRedist(const phase::which_phase a_phase){
  return m_realms[a_phase]->getFineToCoarRedist();
}

Vector<RefCountedPtr<Copier> >& Realm::getCopier(const phase::which_phase a_phase){
  return m_realms[a_phase]->getCopier();
}

Vector<RefCountedPtr<Copier> >& Realm::getReverseCopier(const phase::which_phase a_phase){
  return m_realms[a_phase]->getReverseCopier();
}

const EBAMRFAB& Realm::getLevelset(const phase::which_phase a_phase) const {
  return m_realms[a_phase]->getLevelset();
}

AMRMask& Realm::getMask(const std::string a_mask, const int a_buffer){

  if(!this->queryMask(a_mask, a_buffer)){
    std::string str = "Realm::getMask - could not find mask '" + a_mask + "'";
    MayDay::Abort(str.c_str());
  }

  return m_masks.at(std::pair<std::string, int>(a_mask, a_buffer));
}

#include <CD_NamespaceFooter.H>
