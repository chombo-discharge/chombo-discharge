/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_SigmaSolver.cpp
  @brief  Implementation of CD_SigmaSolver.H
  @author Robert Marskar
*/

// Chombo includes
#include <EBArith.H>

// Our includes
#include <CD_SigmaSolver.H>
#include <CD_DataOps.H>
#include <CD_NamespaceHeader.H>

SigmaSolver::SigmaSolver(){
  this->setVerbosity(-1);
  this->setPhase(phase::gas);
  this->setRealm(Realm::Primal);
}

SigmaSolver::~SigmaSolver(){

}

const std::string SigmaSolver::getRealm() const {
  return m_realm;
}

void SigmaSolver::setRealm(const std::string a_realm){
  m_realm = a_realm;
}

void SigmaSolver::allocateInternals(){
  CH_TIME("SigmaSolver::allocateInternals");
  if(m_verbosity > 5){
    pout() << "SigmaSolver::allocateInternals" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  
  m_amr->allocate(m_phi, m_realm, m_phase, ncomp);
  m_amr->allocate(m_flux,  m_realm, m_phase, ncomp);
}

void SigmaSolver::preRegrid(const int a_lbase, const int a_oldFinestLevel){
  CH_TIME("SigmaSolver::preRegrid");
  if(m_verbosity > 5){
    pout() << "SigmaSolver::preRegrid" << endl;
  }

  const int ncomp = 1;
  const int finest_level = m_amr->getFinestLevel();
  
  m_amr->allocate(m_cache, m_realm, m_phase, ncomp);
  DataOps::setValue(m_cache, 0.0);
  
  for (int lvl = 0; lvl <= a_oldFinestLevel; lvl++){
    m_phi[lvl]->localCopyTo(*m_cache[lvl]);
  }
}

void SigmaSolver::computeRHS(EBAMRIVData& a_rhs){
  CH_TIME("SigmaSolver::computeRHS");
  if(m_verbosity > 5){
    pout() << "SigmaSolver::computeRHS" << endl;
  }
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl ++){
    m_flux[lvl]->localCopyTo(*a_rhs[lvl]);
  }
}

void SigmaSolver::deallocateInternals(){
  CH_TIME("SigmaSolver::deallocateInternals");
  if(m_verbosity > 5){
    pout() << "SigmaSolver::deallocateInternals" << endl;
  }
  
  m_amr->deallocate(m_phi);
  m_amr->deallocate(m_flux);
}

void SigmaSolver::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel){
  CH_TIME("SigmaSolver::regrid");
  if(m_verbosity > 5){
    pout() << "SigmaSolver::regrid" << endl;
  }

  const RefCountedPtr<EBIndexSpace> ebis = m_multifluidIndexSpace->getEBIndexSpace(m_phase);

  const int ebghost = 4; // m_amr->getNumberOfEbGhostCells();
  const int comp  = 0;
  const int ncomp = 1;
  const Interval interv(comp, comp);

  this->allocateInternals();

  DataOps::setValue(m_phi, 0.0);

  // These levels have never changed
  for (int lvl = 0; lvl <= Max(0, a_lmin-1); lvl++){
    m_cache[lvl]->copyTo(*m_phi[lvl]); 
  }

  // These levels have changed
  for (int lvl = Max(1,a_lmin); lvl <= a_newFinestLevel; lvl++){
    const DisjointBoxLayout& fine_grid = m_amr->getGrids(m_realm)[lvl];
    const ProblemDomain& fine_domain   = m_amr->getDomains()[lvl];
    const ProblemDomain& coar_domain   = m_amr->getDomains()[lvl-1];
    const EBISLayout& fine_ebisl       = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const int nref                     = m_amr->getRefinementRatios()[lvl-1];
    
    // Fill a coarsened grid and a coarsened ebisl
    EBISLayout coar_ebisl;
    DisjointBoxLayout coar_grid = DisjointBoxLayout();
    coarsen(coar_grid, fine_grid, nref);
    ebis->fillEBISLayout(coar_ebisl, coar_grid, coar_domain, ebghost);

    // Need extra storage for coarse data
    LayoutData<IntVectSet> sets(coar_grid);
    for (DataIterator dit = coar_grid.dataIterator(); dit.ok(); ++dit){
      Box box          = coar_grid.get(dit());
      box.grow(3);
      box &= coar_domain;
      
      const EBISBox& ebisbox = coar_ebisl[dit()];
      const IntVectSet ivs   = ebisbox.getIrregIVS(box);

      sets[dit()] = ivs;
    }

    // Allocate storage for coarsened data
    BaseIVFactory<Real> ivfact(coar_ebisl, sets);
    LevelData<BaseIVFAB<Real> > coarsened_fine_data(coar_grid, ncomp, m_phi[0]->ghostVect(), ivfact);

    //
    EBLevelDataOps::setVal(coarsened_fine_data, 0.0);
    m_phi[lvl-1]->copyTo(coarsened_fine_data);


    // Loop through coarse grid and interpolate to fine grid
    for (DataIterator dit = coar_grid.dataIterator(); dit.ok(); ++dit){
      BaseIVFAB<Real>& fine_state       = (*m_phi[lvl])[dit()];
      const BaseIVFAB<Real>& coar_state = coarsened_fine_data[dit()];
      const EBISBox& coar_ebisbox       = coar_ebisl[dit()];
      const EBISBox& fine_ebisbox       = fine_ebisl[dit()];

      const IntVectSet ivs   = fine_state.getIVS();
      const EBGraph& ebgraph = fine_state.getEBGraph();

      for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
	const VolIndex& fine_vof = vofit();
	const VolIndex  coar_vof = fine_ebisl.coarsen(fine_vof, nref, dit());
	const Real coarArea      = coar_ebisbox.bndryArea(coar_vof);
	
	// Same sigma for all refined cells such that charge is conserved. 
	const Vector<VolIndex> refined_vofs = coar_ebisl.refine(coar_vof, nref, dit());
	Real fineArea = 0.0;
	for (int i = 0; i < refined_vofs.size(); i++){
	  fineArea += fine_ebisbox.bndryArea(refined_vofs[i]);
	}

	// Initialize conserved charge
	if(fineArea > 0.0 && coarArea > 0.0){
	  fine_state(fine_vof, comp) = coar_state(coar_vof, comp)*nref*coarArea/fineArea;
	}
	else{
	  fine_state(fine_vof, comp) = 0.0;
	}
      }
    }

    // If data already exists, it takes precedence
    if (lvl <= a_oldFinestLevel){
      m_cache[lvl]->copyTo(*m_phi[lvl]);
    }
  }

  this->resetCells(m_phi);
}

void SigmaSolver::registerOperators(){
  CH_TIME("SigmaSolver::registerOperators");
  if(m_verbosity > 5){
    pout() << "SigmaSolver::registerOperators" << endl;
  }

  if(m_amr.isNull()){
    MayDay::Abort("SigmaSolver::registerOperators - need to set AmrMesh!");
  }
  else{
    m_amr->registerOperator(s_eb_coar_ave, m_realm, m_phase);
  }
}

void SigmaSolver::resetCells(EBAMRIVData& a_data){
  CH_TIME("SigmaSolver::resetCells");
  if(m_verbosity > 5){
    pout() << "SigmaSolver::resetCells" << endl;
  }

  const int finest_level = m_amr->getFinestLevel();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const Real dx                = m_amr->getDx()[lvl];
    const MFLevelGrid& mflg      = *m_amr->getMFLevelGrid(m_realm)[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const Box box          = dbl.get(dit());
      BaseIVFAB<Real>& data  = (*a_data[lvl])[dit()];
      IntVectSet ivs         = data.getIVS();
      const EBGraph& ebgraph = data.getEBGraph();
      
      ivs -= mflg.interfaceRegion(box, dit());

      for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();
	for (int comp = 0; comp < data.nComp(); comp++){
	  data(vof, comp) = 0.0;
	}
      }
    }
  }
}

void SigmaSolver::setAmr(const RefCountedPtr<AmrMesh>& a_amr){
  CH_TIME("SigmaSolver::setAmr");
  if(m_verbosity > 5){
    pout() << "SigmaSolver::setAmr" << endl;
  }

  m_amr = a_amr;
}

void SigmaSolver::setComputationalGeometry(const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry){
  CH_TIME("SigmaSolver::setComputationalGeometry");
  if(m_verbosity > 5){
    pout() << "SigmaSolver::setComputationalGeometry" << endl;
  }

  m_computationalGeometry = a_computationalGeometry;
  m_multifluidIndexSpace     = m_computationalGeometry->getMfIndexSpace();
}

void SigmaSolver::setPhase(phase::which_phase a_phase){
  CH_TIME("SigmaSolver::setPhase");
  if(m_verbosity > 5){
    pout() << "SigmaSolver::setPhase" << endl;
  }

  m_phase = a_phase;
}

void SigmaSolver::setSigma(const EBAMRIVData& a_sigma){
  CH_TIME("SigmaSolver::setSigma(ebamrivdata)");
  if(m_verbosity > 5){
    pout() << "SigmaSolver::setSigma(ebamrivdata)" << endl;
  }

  const int finest_level = m_amr->getFinestLevel();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    a_sigma[lvl]->localCopyTo(*m_phi[lvl]);
  }

  this->resetCells(m_phi);
}

void SigmaSolver::setSigma(const Real a_sigma){
  CH_TIME("SigmaSolver::setSigma(constant)");
  if(m_verbosity > 5){
    pout() << "SigmaSolver::setSigma(constant)" << endl;
  }

  const int finest_level = m_amr->getFinestLevel();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    DataOps::setValue(*m_phi[lvl], a_sigma);
  }

  this->resetCells(m_phi);
}

void SigmaSolver::setVerbosity(const int a_verbosity){
  CH_TIME("SigmaSolver::setVerbosity");
  m_verbosity = a_verbosity;
  
  if(m_verbosity > 5){
    pout() << "SigmaSolver::setVerbosity" << endl;
  }
  

}

void SigmaSolver::setTime(const int a_step, const Real a_time, const Real a_dt){
  CH_TIME("SigmaSolver::setTime");
  if(m_verbosity > 5){
    pout() << "SigmaSolver::setTime" << endl;
  }
  
  m_timeStep = a_step;
  m_time = a_time;
  m_dt   = a_dt;
}

#ifdef CH_USE_HDF5
void SigmaSolver::writeCheckpointLevel(HDF5Handle& a_handle, const int a_level) const {
  CH_TIME("SigmaSolver::writeCheckpointLevel");
  if(m_verbosity > 5){
    pout() << "SigmaSolver::writeCheckpointLevel" << endl;
  }

  EBCellFactory fact(m_amr->getEBISLayout(m_realm, phase::gas)[a_level]);
  LevelData<EBCellFAB> scratch(m_amr->getGrids(m_realm)[a_level], 1, 3*IntVect::Unit, fact);
  DataOps::setValue(scratch, 0.0);
  DataOps::incr(scratch, *m_phi[a_level], 1.0);

  // Write state vector
  write(a_handle, scratch, "sigma");
}
#endif

#ifdef CH_USE_HDF5
void SigmaSolver::readCheckpointLevel(HDF5Handle& a_handle, const int a_level){
  CH_TIME("SigmaSolver::readCheckpointLevel");
  if(m_verbosity > 5){
    pout() << "SigmaSolver::readCheckpointLevel" << endl;
  }

  const EBISLayout& ebisl = m_amr->getEBISLayout(m_realm, phase::gas)[a_level];
  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_level];
  
  EBCellFactory fact(ebisl);
  LevelData<EBCellFAB> scratch(dbl, 1, 3*IntVect::Unit, fact);
  DataOps::setValue(scratch, 0.0);
  read<EBCellFAB>(a_handle, scratch, "sigma", dbl, Interval(0,0), false);

		     
  DataOps::setValue(*m_phi[a_level], 0.0);
  DataOps::incr(*m_phi[a_level], scratch, 1.0);
}
#endif

void SigmaSolver::writePlotData(EBAMRCellData& a_output, int& a_comp){
  CH_TIME("SigmaSolver::writePlotData");
  if(m_verbosity > 5){
    pout() << "SigmaSolver::writePlotData" << endl;
  }


  EBAMRCellData scratch;
  m_amr->allocate(scratch, m_realm, m_phase, 1);


  // Write sigma
  DataOps::setValue(scratch, 0.0);
  DataOps::incr(scratch, m_phi, 1.0);
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    if(m_realm == a_output.getRealm()){
      scratch[lvl]->localCopyTo(Interval(0,0), *a_output[lvl], Interval(a_comp, a_comp));
    }
    else{
      scratch[lvl]->copyTo(Interval(0,0), *a_output[lvl], Interval(a_comp, a_comp));
    }
  }
  a_comp++;

  // Write flux
  DataOps::setValue(scratch, 0.0);
  DataOps::incr(scratch, m_flux, 1.0);
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    if(m_realm == a_output.getRealm()){
      scratch[lvl]->localCopyTo(Interval(0,0), *a_output[lvl], Interval(a_comp, a_comp));
    }
    else{
      scratch[lvl]->copyTo(Interval(0,0), *a_output[lvl], Interval(a_comp, a_comp));
    }
  }
  a_comp++;
}

int SigmaSolver::getNumberOfPlotVariables(){
  CH_TIME("SigmaSolver::getNumberOfPlotVariables");
  if(m_verbosity > 5){
    pout() << "SigmaSolver::getNumberOfPlotVariables" << endl;
  }
  
  return 2;
}
  
  
Vector<std::string> SigmaSolver::getPlotVariableNames() const{
  CH_TIME("SigmaSolver::getNumberOfPlotVariables");
  if(m_verbosity > 5){
    pout() << "SigmaSolver::getNumberOfPlotVariables" << endl;
  }
  Vector<std::string> ret(2);

  ret[0] = "surface charge density";
  ret[1] = "surface charge flux";
  
  return ret;
}


Real SigmaSolver::computeCharge(){
  CH_TIME("SigmaSolver::computeCharge");
  if(m_verbosity > 5){
    pout() << "SigmaSolver::computeCharge" << endl;
  }

  m_amr->averageDown(m_phi, m_realm, m_phase);

  Real charge = 0.0;

  const int comp               = 0;
  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[0];
  const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[0];
  const Real dx                = m_amr->getDx()[0];
  
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const EBISBox& ebisbox = ebisl[dit()];
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    const Box box          = dbl.get(dit());
    const IntVectSet irreg = ebisbox.getIrregIVS(box);

    for (VoFIterator vofit(irreg, ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();
      const Real& area    = ebisbox.bndryArea(vof);

      charge += area*(*m_phi[0])[dit()](vof, comp);
    }
  }

  DataOps::sum(charge); // Parallell sum
  
  return charge*dx;
}

EBAMRIVData& SigmaSolver::getPhi(){
  return m_phi;
}

EBAMRIVData& SigmaSolver::getFlux(){
  return m_flux;
}

#include <CD_NamespaceFooter.H>
