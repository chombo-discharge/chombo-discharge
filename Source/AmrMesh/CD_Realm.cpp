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
#include <NeighborIterator.H>

// Our includes
#include <CD_Realm.H>
#include <CD_BoxLoops.H>
#include <CD_NamespaceHeader.H>

const std::string Realm::Primal = "primal";
const std::string Realm::primal = "primal";

Realm::Realm(){
  m_isDefined = false;
  m_verbosity = -1;

  ParmParse pp("Realm");

  pp.query("verbosity", m_verbosity);

  // Just empty points until define() is called
  m_realms.emplace(phase::gas,   RefCountedPtr<PhaseRealm> (new PhaseRealm()));
  m_realms.emplace(phase::solid, RefCountedPtr<PhaseRealm> (new PhaseRealm()));
}

Realm::~Realm(){

}

void Realm::define(const Vector<DisjointBoxLayout>&                           a_grids,
		   const Vector<ProblemDomain>&                               a_domains,
		   const Vector<int>&                                         a_refRat,
		   const Vector<Real>&                                        a_dx,
		   const RealVect                                             a_probLo,
		   const int                                                  a_finestLevel,
		   const int                                                  a_ebGhost,
		   const int                                                  a_numGhost,
		   const int                                                  a_lsfGhost,
		   const int                                                  a_redistRad,
		   const int                                                  a_mgInterpOrder,
		   const int                                                  a_mgInterpRadius,
		   const int                                                  a_mgInterpWeight,
		   const IrregStencil::StencilType                            a_centroidStencil,
		   const IrregStencil::StencilType                            a_ebStencil,
		   const std::map<phase::which_phase, RefCountedPtr<BaseIF> > a_baseif,
		   const RefCountedPtr<MultiFluidIndexSpace>&                 a_mfis){
  CH_TIME("Realm::define");
  if(m_verbosity > 5){
    pout() << "Realm::define" << endl;
  }

  m_grids                = a_grids;
  m_domains              = a_domains;  
  m_refinementRatios     = a_refRat;
  m_dx                   = a_dx;
  m_probLo               = a_probLo;
  m_finestLevel          = a_finestLevel;
  m_baseif               = a_baseif;
  m_multifluidIndexSpace = a_mfis;
  
  const RefCountedPtr<EBIndexSpace>& ebis_gas = m_multifluidIndexSpace->getEBIndexSpace(phase::gas);
  const RefCountedPtr<EBIndexSpace>& ebis_sol = m_multifluidIndexSpace->getEBIndexSpace(phase::solid);

  m_realms[phase::gas]->define(a_grids,
			       a_domains,
			       a_refRat,
			       a_dx,
			       m_probLo,
			       a_finestLevel,
			       a_ebGhost,
			       a_numGhost,
			       a_lsfGhost,
			       a_redistRad,
			       a_mgInterpOrder,
			       a_mgInterpRadius,
			       a_mgInterpWeight,
			       a_centroidStencil,
			       a_ebStencil,
			       m_baseif.at(phase::gas),
			       ebis_gas);

  m_realms[phase::solid]->define(a_grids,
				 a_domains,
				 a_refRat,
				 a_dx,
				 m_probLo,
				 a_finestLevel,
				 a_ebGhost,
				 a_numGhost,
				 a_lsfGhost,
				 a_redistRad,
				 a_mgInterpOrder,
				 a_mgInterpRadius,
				 a_mgInterpWeight,
				 a_centroidStencil,
				 a_ebStencil,
				 m_baseif.at(phase::solid),
				 ebis_sol);

  m_isDefined = true;
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
  this->defineMFLevelGrid(a_lmin);
}

void Realm::regridOperators(const int a_lmin){
  CH_TIME("Realm::regridOperators");
  if(m_verbosity > 5){
    pout() << "Realm::regridOperators" << endl;
  }

  for (auto& r : m_realms){
    r.second->regridOperators(a_lmin);
  }

  this->defineMasks(a_lmin);
}

void Realm::defineMasks(const int a_lmin){
  CH_TIME("Realm::defineMasks");
  if(m_verbosity > 5){
    pout() << "Realm::defineMasks" << endl;
  }

  // Regrid all masks
  this->defineHaloMasks(a_lmin);
}

void Realm::defineMFLevelGrid(const int a_lmin){
  CH_TIME("Realm::defineMFLevelGrid");
  if(m_verbosity > 5){
    pout() << "Realm::defineMFLevelGrid" << endl;
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

void Realm::defineHaloMasks(const int a_lmin){
  CH_TIME("Realm::defineHaloMasks");
  if(m_verbosity > 5){
    pout() << "Realm::defineHaloMasks" << endl;
  }

  // Loop through all masks and do something about the halo masks only. 
  for (auto& m : m_masks){

    // Get mask identifier and buffer. 
    const std::string which_mask = m.first.first;
    const int buffer             = m.first.second;

    if(which_mask == s_particle_halo){
      if(buffer < 0) MayDay::Abort("Realm::defineHaloMasks -- cannot have buffer < 0!");

      AMRMask& mask = m.second;

      mask.resize(1 + m_finestLevel);

      for (int lvl = 0; lvl < m_finestLevel; lvl++){
	const DisjointBoxLayout& gridsCoar = m_grids[lvl];
	const DisjointBoxLayout& gridsFine = m_grids[lvl+1];

	const ProblemDomain& domainCoar = m_domains[lvl];
	const ProblemDomain& domainFine = m_domains[lvl+1];

	const int ncomp = 1;

	mask[lvl] = RefCountedPtr<LevelData<BaseFab<bool> > >(new LevelData<BaseFab<bool> >(gridsCoar, ncomp, IntVect::Zero));

	this->defineHaloMask(*mask[lvl],
			     domainCoar,
			     domainFine,
			     gridsCoar,
			     gridsFine,
			     buffer,
			     m_refinementRatios[lvl]);
      }

      // Must explicitly set this because we want to the finest level mask to be nullptr, but we could be removing a grid level and in that case
      // the old mask will remain in the vector after resizing. This fixes that. 
      mask[m_finestLevel] = RefCountedPtr<LevelData<BaseFab<bool> > >(nullptr);
    }
  }
}

void Realm::defineHaloMask(LevelData<BaseFab<bool> >& a_coarMask,
			   const ProblemDomain&       a_domainCoar,
			   const ProblemDomain&       a_domainFine,
			   const DisjointBoxLayout&   a_gridsCoar,
			   const DisjointBoxLayout&   a_gridsFine,
			   const int                  a_buffer,
			   const int                  a_refRat){
  CH_TIME("Realm::defineHaloMasks");
  if(m_verbosity > 5){
    pout() << "Realm::defineHaloMasks" << endl;
  }

  // TLDR: This routine defines a "mask" of valid coarse-grid cells (i.e., not covered by a finer level) around a fine grid. The mask is a_buffer wide. It is
  //       created by first fetching the cells around on the fine grid. This is a local operation where we get all the cells around each patch and then subtract
  //       cells that overlap with other boxes. This set of cells is then put on a LevelData<FArrayBox> mask so we can copy the result to the coarse grid and
  //       set the mask. 
  
  constexpr int comp  = 0;
  constexpr int ncomp = 1;

  // First, reset the mask.
  for (DataIterator dit(a_gridsCoar); dit.ok(); ++dit){
    a_coarMask[dit()].setVal(false);
  }

  // Ok, we need a particle halo. 
  if(a_buffer > 0){

    // Coarsen the fine grid and make a mask on the coarsened fine grid
    DisjointBoxLayout dblCoFi;
    coarsen(dblCoFi, a_gridsFine, a_refRat);

    IntVectSet halo;
    
    // Go through the cofi grid and set the halo to true
    for (DataIterator dit(dblCoFi); dit.ok(); ++dit){
      const Box coFiBox = dblCoFi[dit()];

      // Make IntVect set consisting of only ghost cells (a_buffer) for each box. 
      Box grownBox = grow(coFiBox, a_buffer);
      grownBox    &= a_domainCoar;

      // Subtract non-ghosted box. 
      IntVectSet myHalo(grownBox);
      myHalo -= coFiBox;

      // Subtract non-ghosted neighboring boxes.
      NeighborIterator nit(dblCoFi); // Neighbor iterator      
      for (nit.begin(dit()); nit.ok(); ++nit){
	const Box neighborBox = dblCoFi[nit()];
	myHalo -= neighborBox;
      }

      // Add to halo. 
      halo |= myHalo;
    }

    // TLDR: In the above, we found the coarse-grid cells surrounding the fine level, viewed from the fine grids.
    //       Below, we create that view from the coarse grid. We use a BoxLayoutData<FArrayBox> on the coarsened fine grid,
    //       whose "ghost cells" can added to the _actual_ coarse grid. We then loop through those cells and set the mask. 
    LevelData<FArrayBox> coFiMask(dblCoFi,     ncomp, a_buffer*IntVect::Unit);
    LevelData<FArrayBox> coarMask(a_gridsCoar, ncomp,          IntVect::Zero);    
  
    // Reset masks
    for (DataIterator dit(dblCoFi); dit.ok(); ++dit){
      coFiMask[dit()].setVal(0.0);
    }
    for (DataIterator dit(a_gridsCoar); dit.ok(); ++dit){
      coarMask[dit()].setVal(0.0);
    }

    // Run through the halo and set the halo cells to 1 on the coarsened fine grids. Since dblCoFi was a coarsening of the fine
    // grid, the mask value in the valid region (i.e., not including ghosts) is always zero. There used to be a bug here because
    // we only iterated through that region, but obviously the halo masks will always be zero in that case...
    for (DataIterator dit(dblCoFi); dit.ok(); ++dit){
      const Box        region  = coFiMask[dit()].box();
      const IntVectSet curHalo = halo & region;
      for (IVSIterator ivsit(curHalo); ivsit.ok(); ++ivsit){
	coFiMask[dit()](ivsit(), comp) = 1.0;
      }
    }

    // Add the result to the coarse grid.
    Copier copier;
    copier.ghostDefine(dblCoFi, a_gridsCoar, a_domainCoar, a_buffer*IntVect::Unit);
    const Interval interv(0,0);
    coFiMask.copyTo(interv, coarMask, interv, copier, LDaddOp<FArrayBox>());


    // Run through the grids and make the boolean mask
    for (DataIterator dit(a_gridsCoar); dit.ok(); ++dit){
      const Box box = a_gridsCoar[dit()];

      BaseFab<bool>&   boolMask = a_coarMask[dit()];
      const FArrayBox& realMask =   coarMask[dit()];

      auto kernel = [&] (const IntVect& iv) -> void {
	if(realMask(iv, comp) > 0.0){
	  boolMask(iv, comp) = true;
	}
      };

      BoxLoops::loop(box, kernel);
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

bool Realm::queryOperator(const std::string a_operator, const phase::which_phase a_phase) const {
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

const Vector<int>& Realm::getRefinementRatios() const {
  return m_refinementRatios;
}

const Vector<Real>& Realm::getDx() const {
  return m_dx;
}

const Vector<DisjointBoxLayout>& Realm::getGrids() const {
  return m_grids;
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

Vector<RefCountedPtr<EbCoarAve> >& Realm::getCoarseAverage(const phase::which_phase a_phase){
  return m_realms[a_phase]->getCoarseAverage();
}

EBAMRParticleMesh& Realm::getParticleMesh(const phase::which_phase a_phase){
  return m_realms[a_phase]->getParticleMesh();
}

Vector<RefCountedPtr<EBMultigridInterpolator> >& Realm::getMultigridInterpolator(const phase::which_phase a_phase){
  return m_realms[a_phase]->getMultigridInterpolator();
}

const Vector<RefCountedPtr<EBGradient> >& Realm::getGradientOp(const phase::which_phase a_phase) const {
  return m_realms[a_phase]->getGradientOp();
}

Vector<RefCountedPtr<AggEBPWLFillPatch> >& Realm::getFillPatch(const phase::which_phase a_phase){
  return m_realms[a_phase]->getFillPatch();
}

Vector<RefCountedPtr<EBPWLFineInterp> >& Realm::getPwlInterpolator(const phase::which_phase a_phase){
  return m_realms[a_phase]->getPwlInterpolator();
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

const EBAMRFAB& Realm::getLevelset(const phase::which_phase a_phase) const {
  return m_realms[a_phase]->getLevelset();
}

const AMRMask& Realm::getMask(const std::string a_mask, const int a_buffer) const{
  if(!this->queryMask(a_mask, a_buffer)){
    std::string str = "Realm::getMask - could not find mask '" + a_mask + "'";
    MayDay::Abort(str.c_str());
  }

  return m_masks.at(std::pair<std::string, int>(a_mask, a_buffer));
}

#include <CD_NamespaceFooter.H>
