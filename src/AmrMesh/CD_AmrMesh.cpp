/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   AmrMesh.cpp
  @brief  Implementation of AmrMesh.H
  @author Robert Marskar
*/

// Chombo includes
#include <BRMeshRefine.H>
#include <EBEllipticLoadBalance.H>
#include <EBLevelDataOps.H>
#include <MFLevelDataOps.H>
#include <EBArith.H>
#include <ParmParse.H>
#include <BaseIFFactory.H>

// Our includes
#include <CD_AmrMesh.H>
#include <mfalias.H>
#include <load_balance.H>
#include <gradientF_F.H>
#include <EBFastFineToCoarRedist.H>
#include <EBFastCoarToFineRedist.H>
#include <EBFastCoarToCoarRedist.H>
#include <DomainFluxIFFABFactory.H>
#include <TiledMeshRefine.H>

#include <CD_NamespaceHeader.H>

AmrMesh::AmrMesh(){
  parseOptions();

  m_finestLevel = 0;
  m_hasGrids    = false;

#if 1
  // I have no idea why this work, and we should check that out. But, for some reason
  // if we use refinement of 4, without a 2 at the end (even if this is NOT used anywhere?!?!?), we get memory
  // leaks through define_eblevelgrid() in the setMaxCoarseningRatio stuff. I don't know if this leak is real,
  // but that's what valgrind tells me....
  //
  // Robert Marskar, Dec 7 (2018)
  //
  m_refinementRatios.resize(m_maxAmrDepth);
  m_refinementRatios.push_back(2);
#endif

}

AmrMesh::~AmrMesh(){
  
}

EBAMRCellData AmrMesh::alias(const phase::which_phase a_phase, const MFAMRCellData& a_mfdata){
  CH_TIME("AmrMesh::alias(phase, mfamrcelldata)");
  if(m_verbosity > 5){
    pout() << "AmrMesh::alias(phase, mfamrcelldata)" << endl;
  }

  EBAMRCellData ret;

  const int finestLevel = a_mfdata.size() - 1;

  this->allocatePointer(ret, finestLevel);

  this->alias(ret, a_phase, a_mfdata, finestLevel);

  return ret;
}

EBAMRFluxData AmrMesh::alias(const phase::which_phase a_phase, const MFAMRFluxData& a_mfdata){
  CH_TIME("AmrMesh::alias(phase, mfamrfluxdata)");
  if(m_verbosity > 5){
    pout() << "AmrMesh::alias(phase, mfamrfluxdata)" << endl;
  }

  EBAMRFluxData ret;

  const int finestLevel = a_mfdata.size() - 1;
  
  this->allocatePointer(ret, finestLevel);

  this->alias(ret, a_phase, a_mfdata, finestLevel);

  return ret;
}

EBAMRIVData AmrMesh::alias(const phase::which_phase a_phase, const MFAMRIVData& a_mfdata){
  CH_TIME("AmrMesh::alias(phase, mfamrivdata)");
  if(m_verbosity > 5){
    pout() << "AmrMesh::alias(phase, mfamrivdata)" << endl;
  }

  EBAMRIVData ret;
  
  this->allocatePointer(ret);

  this->alias(ret, a_phase, a_mfdata);

  return ret;
}

void AmrMesh::alias(EBAMRCellData&           a_data,
		    const phase::which_phase a_phase,
		    const MFAMRCellData&     a_mfdata,
		    const int                a_finestLevel){
  CH_TIME("AmrMesh::alias(hardcap)");
  if(m_verbosity > 5){
    pout() << "AmrMesh::alias(hardcap)" << endl;
  }

  for (int lvl = 0; lvl <= a_finestLevel; lvl++){
    mfalias::aliasMF(*a_data[lvl], a_phase, *a_mfdata[lvl]);
  }
}

void AmrMesh::alias(EBAMRFluxData&           a_data,
		    const phase::which_phase a_phase,
		    const MFAMRFluxData&     a_mfdata,
		    const int                a_finestLevel){
  CH_TIME("AmrMesh::alias(hardcap)");
  if(m_verbosity > 5){
    pout() << "AmrMesh::alias(hardcap)" << endl;
  }

  for (int lvl = 0; lvl <= a_finestLevel; lvl++){
    mfalias::aliasMF(*a_data[lvl], a_phase, *a_mfdata[lvl]);
  }
}

void AmrMesh::alias(EBAMRCellData& a_data, const phase::which_phase a_phase, const MFAMRCellData& a_mfdata){
  CH_TIME("AmrMesh::alias");
  if(m_verbosity > 5){
    pout() << "AmrMesh::alias" << endl;
  }

  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    mfalias::aliasMF(*a_data[lvl], a_phase, *a_mfdata[lvl]);
  }
}

void AmrMesh::alias(EBAMRFluxData& a_data, const phase::which_phase a_phase, const MFAMRFluxData& a_mfdata){
  CH_TIME("AmrMesh::alias");
  if(m_verbosity > 5){
    pout() << "AmrMesh::alias" << endl;
  }

  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    mfalias::aliasMF(*a_data[lvl], a_phase, *a_mfdata[lvl]);
  }
}

void AmrMesh::alias(EBAMRIVData& a_data, const phase::which_phase a_phase, const MFAMRIVData& a_mfdata){
  CH_TIME("AmrMesh::alias");
  if(m_verbosity > 5){
    pout() << "AmrMesh::alias" << endl;
  }

  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    mfalias::aliasMF(*a_data[lvl], a_phase, *a_mfdata[lvl]);
  }
}

void AmrMesh::allocate(EBAMRCellData&           a_data,
		       const std::string        a_realm,
		       const phase::which_phase a_phase,
		       const int                a_ncomp,
		       const int                a_ghost){
  CH_TIME("AmrMesh::allocate(ebamrcell, realm, phase, comp, ghost)");
  if(m_verbosity > 5){
    pout() << "AmrMesh::allocate(ebamrcell, realm, phase, comp, ghost)" << endl;
  }

  if(!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::allocate(ebamcell, realm, phase, comp, ghost) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const int ghost = (a_ghost == -1) ? m_numGhostCells : a_ghost;

  a_data.resize(1 + m_finestLevel);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    const DisjointBoxLayout& dbl = m_realms[a_realm]->getGrids()[lvl];
    const EBISLayout& ebisl      = m_realms[a_realm]->getEBISLayout(a_phase)[lvl];
    
    EBCellFactory fact(ebisl);

    a_data[lvl] = RefCountedPtr<LevelData<EBCellFAB> >
      (new LevelData<EBCellFAB>(dbl, a_ncomp, ghost*IntVect::Unit, fact));

    EBLevelDataOps::setVal(*a_data[lvl], 0.0);
  }

  a_data.set_realm(a_realm);
}

void AmrMesh::allocate(EBAMRFluxData&           a_data,
		       const std::string        a_realm,
		       const phase::which_phase a_phase,
		       const int                a_ncomp,
		       const int                a_ghost){
  CH_TIME("AmrMesh::allocate(ebamrflux, realm, phase, comp, ghost)");
  if(m_verbosity > 5){
    pout() << "AmrMesh::allocate(ebamrflux, realm, phase, comp, ghost)" << endl;
  }

  if(!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::allocate(ebamrflux, realm, phase, comp, ghost) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const int ghost = (a_ghost == -1) ? m_numGhostCells : a_ghost;

  a_data.resize(1 + m_finestLevel);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    const DisjointBoxLayout& dbl = m_realms[a_realm]->getGrids()[lvl];
    const EBISLayout& ebisl      = m_realms[a_realm]->getEBISLayout(a_phase)[lvl];

    EBFluxFactory fact(ebisl);

    a_data[lvl] = RefCountedPtr<LevelData<EBFluxFAB> >
      (new LevelData<EBFluxFAB>(dbl, a_ncomp, ghost*IntVect::Unit, fact));

    EBLevelDataOps::setVal(*a_data[lvl], 0.0);
  }

  a_data.set_realm(a_realm);
}

void AmrMesh::allocate(EBAMRIVData&             a_data,
		       const std::string        a_realm,
		       const phase::which_phase a_phase,
		       const int                a_ncomp,
		       const int                a_ghost){
  CH_TIME("AmrMesh::allocate(ebamriv, realm, phase, comp, ghost)");
  if(m_verbosity > 5){
    pout() << "AmrMesh::allocate(ebamriv, realm, phase, comp, ghost)" << endl;
  }

  if(!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::allocate(ebamriv, realm, phase, comp, ghost) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const int ghost = (a_ghost == -1) ? m_numGhostCells : a_ghost;

  a_data.resize(1 + m_finestLevel);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    const DisjointBoxLayout& dbl = m_realms[a_realm]->getGrids()[lvl];
    const EBISLayout& ebisl      = m_realms[a_realm]->getEBISLayout(a_phase)[lvl];
    const ProblemDomain& domain  = m_realms[a_realm]->getDomains()[lvl];
    
    LayoutData<IntVectSet> irreg_sets(dbl);
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      Box box = dbl.get(dit());
      box.grow(ghost);
      box &= domain;

      irreg_sets[dit()] = ebisl[dit()].getIrregIVS(box);
    }

    BaseIVFactory<Real> fact(ebisl, irreg_sets);

    a_data[lvl] = RefCountedPtr<LevelData<BaseIVFAB<Real> > >
      (new LevelData<BaseIVFAB<Real> >(dbl, a_ncomp, ghost*IntVect::Unit, fact));

    EBLevelDataOps::setVal(*a_data[lvl], 0.0);
  }

  a_data.set_realm(a_realm);
}

void AmrMesh::allocate(EBAMRIFData&             a_data,
		       const std::string        a_realm,
		       const phase::which_phase a_phase,
		       const int                a_ncomp,
		       const int                a_ghost){
  CH_TIME("AmrMesh::allocate(ebamrifdata, realm, phase, comp, ghost)");
  if(m_verbosity > 5){
    pout() << "AmrMesh::allocate(ebamrifdata, realm, phase, comp, ghost)" << endl;
  }

  if(!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::allocate(ebamrifdata, realm, phase, comp, ghost) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const int ghost = (a_ghost == -1) ? m_numGhostCells : a_ghost;

  a_data.resize(1 + m_finestLevel);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    const DisjointBoxLayout& dbl = m_realms[a_realm]->getGrids()[lvl];
    const EBISLayout& ebisl      = m_realms[a_realm]->getEBISLayout(a_phase)[lvl];
    const ProblemDomain& domain  = m_realms[a_realm]->getDomains()[lvl];
    
    DomainFluxIFFABFactory fact(ebisl, domain);

    a_data[lvl] = RefCountedPtr<LevelData<DomainFluxIFFAB> >
      (new LevelData<DomainFluxIFFAB>(dbl, a_ncomp, ghost*IntVect::Unit, fact));
  }

  a_data.set_realm(a_realm);
}

void AmrMesh::allocate(EBAMRBool& a_data, const std::string a_realm, const int a_ncomp, const int a_ghost){
  CH_TIME("AmrMesh::allocate(ebamrbool, realm, comp, ghost)");
  if(m_verbosity > 5){
    pout() << "AmrMesh::allocate(ebamrbool, realm, comp, ghost)" << endl;
  }

  if(!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::allocate(ebamrbool, realm, phase, comp, ghost) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  a_data.resize(1 + m_finestLevel);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    const DisjointBoxLayout& dbl = m_realms[a_realm]->getGrids()[lvl];
    
    a_data[lvl] = RefCountedPtr<LevelData<BaseFab<bool> > >
      (new LevelData<BaseFab<bool> >(dbl, a_ncomp, a_ghost*IntVect::Unit));
  }

  a_data.set_realm(a_realm);
}

void AmrMesh::allocate(MFAMRCellData& a_data, const std::string a_realm, const int a_ncomp, const int a_ghost){
  CH_TIME("AmrMesh::allocate(mfamrcell, realm, comp, ghost)");
  if(m_verbosity > 5){
    pout() << "AmrMesh::allocate(mfamrcell, realm, comp, ghost)" << endl;
  }

  if(!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::allocate(mfamrcell, realm, comp, ghost) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const int ghost   = (a_ghost == -1) ? m_numGhostCells : a_ghost;
  const int ignored = a_ncomp;
  const int nphases = m_multifluidIndexSpace->num_phases();

  a_data.resize(1 + m_finestLevel);

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_realms[a_realm]->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_realms[a_realm]->get_ebis(phase::solid);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    const DisjointBoxLayout& dbl = m_realms[a_realm]->getGrids()[lvl];
    
    Vector<EBISLayout> ebisl(nphases);
    Vector<int>        comps(nphases, a_ncomp);

    if(!ebis_gas.isNull()) ebisl[phase::gas]   = m_realms[a_realm]->getEBISLayout(phase::gas)[lvl];
    if(!ebis_sol.isNull()) ebisl[phase::solid] = m_realms[a_realm]->getEBISLayout(phase::solid)[lvl];
    
    MFCellFactory factory(ebisl, comps);

    a_data[lvl] = RefCountedPtr<LevelData<MFCellFAB> >
      (new LevelData<MFCellFAB>(dbl, ignored, ghost*IntVect::Unit, factory));

    MFLevelDataOps::setVal(*a_data[lvl], 0.0);
  }

  a_data.set_realm(a_realm);
}

void AmrMesh::allocate(MFAMRFluxData& a_data, const std::string a_realm, const int a_ncomp, const int a_ghost){
  CH_TIME("AmrMesh::allocate(mfamrflux, realm, comp, ghost)");
  if(m_verbosity > 5){
    pout() << "AmrMesh::allocate(mfamrflux, realm, comp, ghost)" << endl;
  }

  if(!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::allocate(mfamrflux, realm, comp, ghost) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const int ghost   = (a_ghost == -1) ? m_numGhostCells : a_ghost;
  const int ignored = a_ncomp;
  const int nphases = m_multifluidIndexSpace->num_phases();

  a_data.resize(1 + m_finestLevel);

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_realms[a_realm]->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_realms[a_realm]->get_ebis(phase::solid);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    const DisjointBoxLayout& dbl = m_realms[a_realm]->getGrids()[lvl];
    
    Vector<EBISLayout> ebisl(nphases);
    Vector<int>        comps(nphases, a_ncomp);

    if(!ebis_gas.isNull()) ebisl[phase::gas]   = m_realms[a_realm]->getEBISLayout(phase::gas)[lvl];
    if(!ebis_sol.isNull()) ebisl[phase::solid] = m_realms[a_realm]->getEBISLayout(phase::solid)[lvl];
    
    MFFluxFactory factory(ebisl, comps);

    a_data[lvl] = RefCountedPtr<LevelData<MFFluxFAB> >
      (new LevelData<MFFluxFAB>(dbl, ignored, ghost*IntVect::Unit, factory));
  }

  a_data.set_realm(a_realm);
}

void AmrMesh::allocate(MFAMRIVData& a_data, const std::string a_realm, const int a_ncomp, const int a_ghost){
  CH_TIME("AmrMesh::allocate(mfamrivdata, realm, comp, ghost)");
  if(m_verbosity > 5){
    pout() << "AmrMesh::allocate(mfamrivdata, realm, comp, ghost)" << endl;
  }
  
  if(!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::allocate(mfamrivdata, realm, comp, ghost) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const int ghost   = (a_ghost == -1) ? m_numGhostCells : a_ghost;
  const int ignored = a_ncomp;
  const int nphases = m_multifluidIndexSpace->num_phases();

  a_data.resize(1 + m_finestLevel);

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_realms[a_realm]->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_realms[a_realm]->get_ebis(phase::solid);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    const DisjointBoxLayout& dbl = m_realms[a_realm]->getGrids()[lvl];
    
    Vector<EBISLayout> ebisl(nphases);
    Vector<int>        comps(nphases, a_ncomp);

    if(!ebis_gas.isNull()) ebisl[phase::gas]   = m_realms[a_realm]->getEBISLayout(phase::gas)[lvl];
    if(!ebis_sol.isNull()) ebisl[phase::solid] = m_realms[a_realm]->getEBISLayout(phase::solid)[lvl];
    
    MFBaseIVFABFactory factory(ebisl, comps);

    a_data[lvl] = RefCountedPtr<LevelData<MFBaseIVFAB> >
      (new LevelData<MFBaseIVFAB>(dbl, ignored, ghost*IntVect::Unit, factory));
  }

  a_data.set_realm(a_realm);
}

void AmrMesh::reallocate(EBAMRCellData& a_data, const phase::which_phase a_phase, const int a_lmin){
  CH_TIME("AmrMesh::reallocate(ebamrcell, realm, phase, lmin)");
  if(m_verbosity > 5){
    pout() << "AmrMesh::reallocate(ebamrcell, realm, phase, lmin)" << endl;
  }

  const std::string a_realm = a_data.get_realm();

  if(!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::allocate(ebamrcell, realm, phase, lmin) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const IntVect ghost = a_data[0]->ghostVect();
  const int ncomp     = a_data[0]->nComp();

  a_data.resize(1 + m_finestLevel);

  const int lmin = Max(0, a_lmin+1);
  
  for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++){
    const DisjointBoxLayout& dbl = m_realms[a_realm]->getGrids()[lvl];
    const EBISLayout& ebisl      = m_realms[a_realm]->getEBISLayout(a_phase)[lvl];
    
    EBCellFactory fact(ebisl);

    a_data[lvl] = RefCountedPtr<LevelData<EBCellFAB> >
      (new LevelData<EBCellFAB>(dbl, ncomp, ghost, fact));

    EBLevelDataOps::setVal(*a_data[lvl], 0.0);
  }
}

void AmrMesh::reallocate(EBAMRFluxData& a_data, const phase::which_phase a_phase, const int a_lmin){
  CH_TIME("AmrMesh::reallocate(ebamrflux, realm, phase, lmin)");
  if(m_verbosity > 5){
    pout() << "AmrMesh::reallocate(ebamrflux, realm, phase, lmin)" << endl;
  }

  const std::string a_realm = a_data.get_realm();

  if(!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::reallocate(ebamrflux, realm, phase, lmin) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const IntVect ghost = a_data[0]->ghostVect();
  const int ncomp     = a_data[0]->nComp();

  a_data.resize(1 + m_finestLevel);

  const int lmin = Max(0, a_lmin+1);
  
  for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++){
    const DisjointBoxLayout& dbl = m_realms[a_realm]->getGrids()[lvl];
    const EBISLayout& ebisl      = m_realms[a_realm]->getEBISLayout(a_phase)[lvl];
    
    EBFluxFactory fact(ebisl);

    a_data[lvl] = RefCountedPtr<LevelData<EBFluxFAB> > (new LevelData<EBFluxFAB>(dbl, ncomp, ghost, fact));

    EBLevelDataOps::setVal(*a_data[lvl], 0.0);
  }
}

void AmrMesh::reallocate(EBAMRIVData& a_data, const phase::which_phase a_phase, const int a_lmin){
  CH_TIME("AmrMesh::reallocate(ebamriv, realm, phase, lmin)");
  if(m_verbosity > 5){
    pout() << "AmrMesh::reallocate(ebamriv, realm, phase, lmin)" << endl;
  }

  const std::string a_realm = a_data.get_realm();

  if(!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::reallocate(ebamriv, realm, phase, lmin) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const IntVect ghost = a_data[0]->ghostVect();
  const int ncomp     = a_data[0]->nComp();

  a_data.resize(1 + m_finestLevel);

  const int lmin = Max(0, a_lmin+1);

  for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++){
    const DisjointBoxLayout& dbl = m_realms[a_realm]->getGrids()[lvl];
    const EBISLayout& ebisl      = m_realms[a_realm]->getEBISLayout(a_phase)[lvl];
    const ProblemDomain& domain  = m_realms[a_realm]->getDomains()[lvl];
    
    LayoutData<IntVectSet> irreg_sets(dbl);
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      Box box = dbl.get(dit());
      box.grow(ghost);
      box &= domain;

      irreg_sets[dit()] = ebisl[dit()].getIrregIVS(box);
    }

    BaseIVFactory<Real> fact(ebisl, irreg_sets);

    a_data[lvl] = RefCountedPtr<LevelData<BaseIVFAB<Real> > >
      (new LevelData<BaseIVFAB<Real> >(dbl, ncomp, ghost, fact));

    EBLevelDataOps::setVal(*a_data[lvl], 0.0);
  }
}

void AmrMesh::reallocate(EBAMRIFData& a_data, const phase::which_phase a_phase, const int a_lmin){
  CH_TIME("AmrMesh::reallocate(ebamrifdata, realm, phase, lmin)");
  if(m_verbosity > 5){
    pout() << "AmrMesh::reallocate(ebamrifdata, realm, phase, lmin)" << endl;
  }

  const std::string a_realm = a_data.get_realm(); 

  if(!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::reallocate(ebamrif, realm, phase, lmin) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const IntVect ghost = a_data[0]->ghostVect();
  const int ncomp     = a_data[0]->nComp();

  a_data.resize(1 + m_finestLevel);

  const int lmin = Max(0, a_lmin+1);

  for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++){
    const DisjointBoxLayout& dbl = m_realms[a_realm]->getGrids()[lvl];
    const EBISLayout& ebisl      = m_realms[a_realm]->getEBISLayout(a_phase)[lvl];
    const ProblemDomain& domain  = m_realms[a_realm]->getDomains()[lvl];
    
    DomainFluxIFFABFactory fact(ebisl, domain);

    a_data[lvl] = RefCountedPtr<LevelData<DomainFluxIFFAB> >
      (new LevelData<DomainFluxIFFAB>(dbl, ncomp, ghost, fact));
  }
}

void AmrMesh::reallocate(EBAMRBool& a_data, const int a_lmin){
  CH_TIME("AmrMesh::reallocate(ebamrbool, realm, lmin)");
  if(m_verbosity > 5){
    pout() << "AmrMesh::reallocate(ebamrbool, realm, lmin)" << endl;
  }

  const std::string a_realm = a_data.get_realm();

  if(!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::reallocate(ebamrbool, realm, lmin) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }
  
  const IntVect ghost = a_data[0]->ghostVect();
  const int ncomp     = a_data[0]->nComp();

  a_data.resize(1 + m_finestLevel);

  const int lmin = Max(0, a_lmin+1);

  for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++){
    const DisjointBoxLayout& dbl = m_realms[a_realm]->getGrids()[lvl];

    a_data[lvl] = RefCountedPtr<LevelData<BaseFab<bool> > >
      (new LevelData<BaseFab<bool> >(dbl, ncomp, ghost));
  }
}

void AmrMesh::reallocate(MFAMRCellData& a_data, const int a_lmin){
  CH_TIME("AmrMesh::reallocate(mfamrcell, realm, lmin)");
  if(m_verbosity > 5){
    pout() << "AmrMesh::reallocate(mfamrcell, realm, lmin)" << endl;
  }

  const std::string a_realm = a_data.get_realm();

  if(!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::reallocate(mfamrcell, realm, lmin) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const IntVect ghost = a_data[0]->ghostVect();
  const int ncomp = a_data[0]->nComp();
  const int ignored = ncomp;
  const int nphases = m_multifluidIndexSpace->num_phases();

  a_data.resize(1 + m_finestLevel);

  const int lmin = Max(0, a_lmin+1);

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_realms[a_realm]->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_realms[a_realm]->get_ebis(phase::solid);

  for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++){
    Vector<EBISLayout> ebisl(nphases);
    Vector<int>        comps(nphases, ncomp);

    const DisjointBoxLayout& dbl = m_realms[a_realm]->getGrids()[lvl];

    if(!ebis_gas.isNull()) ebisl[phase::gas]   = m_realms[a_realm]->getEBISLayout(phase::gas)[lvl];
    if(!ebis_sol.isNull()) ebisl[phase::solid] = m_realms[a_realm]->getEBISLayout(phase::solid)[lvl];
    
    MFCellFactory factory(ebisl, comps);

    a_data[lvl] = RefCountedPtr<LevelData<MFCellFAB> > (new LevelData<MFCellFAB>(dbl, ignored, ghost, factory));

    MFLevelDataOps::setVal(*a_data[lvl], 0.0);
  }
}

void AmrMesh::reallocate(MFAMRFluxData& a_data, const int a_lmin){
  CH_TIME("AmrMesh::allocate(mfamrflux, realm, lmin)");
  if(m_verbosity > 5){
    pout() << "AmrMesh::allocate(mfamrflux, realm, lmin)" << endl;
  }

  const std::string a_realm = a_data.get_realm(); 

  if(!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::reallocate(mfamrflux, realm, lmin) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const IntVect ghost = a_data[0]->ghostVect();
  const int ncomp     = a_data[0]->nComp();
  const int ignored   = ncomp;
  const int nphases   = m_multifluidIndexSpace->num_phases();

  a_data.resize(1 + m_finestLevel);

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_realms[a_realm]->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_realms[a_realm]->get_ebis(phase::solid);

  const int lmin = Max(0, a_lmin+1);

  for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++){
    const DisjointBoxLayout& dbl = m_realms[a_realm]->getGrids()[lvl];
    
    Vector<EBISLayout> ebisl(nphases);
    Vector<int>        comps(nphases, ncomp);

    if(!ebis_gas.isNull()) ebisl[phase::gas]   = m_realms[a_realm]->getEBISLayout(phase::gas)[lvl];
    if(!ebis_sol.isNull()) ebisl[phase::solid] = m_realms[a_realm]->getEBISLayout(phase::solid)[lvl];
    
    MFFluxFactory factory(ebisl, comps);

    a_data[lvl] = RefCountedPtr<LevelData<MFFluxFAB> >
      (new LevelData<MFFluxFAB>(dbl, ignored, ghost, factory));
  }
}

void AmrMesh::reallocate(MFAMRIVData& a_data, const int a_lmin){
  CH_TIME("AmrMesh::allocate(mfamrivdata, realm, lmin)");
  if(m_verbosity > 5){
    pout() << "AmrMesh::allocate(mfamrivdata, realm, lmin)" << endl;
  }

  const std::string a_realm = a_data.get_realm();

  if(!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::reallocate(mfamrivdata, realm, lmin) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const IntVect ghost = a_data[0]->ghostVect();
  const int ncomp     = a_data[0]->nComp();
  const int ignored   = ncomp;
  const int nphases   = m_multifluidIndexSpace->num_phases();

  a_data.resize(1 + m_finestLevel);

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_realms[a_realm]->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_realms[a_realm]->get_ebis(phase::solid);

  const int lmin = Max(0, a_lmin+1);

  for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++){
    const DisjointBoxLayout& dbl = m_realms[a_realm]->getGrids()[lvl];
    
    Vector<EBISLayout> ebisl(nphases);
    Vector<int>        comps(nphases, ncomp);

    if(!ebis_sol.isNull()) ebisl[phase::gas]   = m_realms[a_realm]->getEBISLayout(phase::gas)[lvl];
    if(!ebis_sol.isNull()) ebisl[phase::solid] = m_realms[a_realm]->getEBISLayout(phase::solid)[lvl];
    
    MFBaseIVFABFactory factory(ebisl, comps);

    a_data[lvl] = RefCountedPtr<LevelData<MFBaseIVFAB> >
      (new LevelData<MFBaseIVFAB>(dbl, ignored, ghost, factory));
  }
}

void AmrMesh::setMultifluidIndexSpace(const RefCountedPtr<mfis>& a_mfis){
  CH_TIME("AmrMesh::setMultifluidIndexSpace");
  if(m_verbosity > 5){
    pout() << "AmrMesh::setMultifluidIndexSpace" << endl;
  }

  m_multifluidIndexSpace = a_mfis;
}

void AmrMesh::setBaseImplicitFunction(const phase::which_phase a_phase, const RefCountedPtr<BaseIF>& a_baseif){
  CH_TIME("AmrMesh::setBaseImplicitFunction");
  if(m_verbosity > 5){
    pout() << "AmrMesh::setBaseImplicitFunction" << endl;
  }

  m_baseif.emplace(a_phase, a_baseif);
}

void AmrMesh::parseOptions(){
  parseVerbosity();
  parseCoarsestLevelNumCells();
  parseMaxAmrDepth();
  parseMaxSimulationDepth();
  parseEbCf();
  parseRefinementRatios();
  parseBlockingFactor();
  parseMaxBoxSize();
  parseMaxEbisBoxSize();
  parsegridGeneration();
  parseBrBufferSize();
  parseIrregTagGrowth();
  parseBrFillRatio();
  parseRedistributionRadius();;
  parseNumGhostCells();
  parseEbGhostCells();
  parseProbLoHiCorners();
  parseGhostInterpolation();
  parseCentroidStencils();
  parseEbCentroidStencils();
}

void AmrMesh::parseRuntimeOptions(){
  parseVerbosity();
  parseBlockingFactor();
  parseMaxBoxSize();
  parsegridGeneration();
  parseBrBufferSize();
  parseIrregTagGrowth();
  parseBrFillRatio();
}

void AmrMesh::parseProbLoHiCorners(){

  ParmParse pp("AmrMesh");

  Vector<Real> v(SpaceDim);
  pp.getarr("lo_corner", v, 0, SpaceDim); m_probLo  = RealVect(D_DECL(v[0], v[1], v[2]));
  pp.getarr("hi_corner", v, 0, SpaceDim); m_prob_hi = RealVect(D_DECL(v[0], v[1], v[2]));
}

void AmrMesh::parseGhostInterpolation(){

  std::string interp_type;
  ParmParse pp("AmrMesh");
  pp.get("ghost_interp", interp_type);
  if(interp_type == "pwl"){
    m_ghostCellInterpolationMethod = GhostInterpolation::PiecewiseLinear;
  }
  else if(interp_type == "quad"){
    m_ghostCellInterpolationMethod = GhostInterpolation::Quadratic;
  }
  else{
    MayDay::Abort("AmrMesh::parseGhostInterpolation - unknown ghost interpolation requested");
  }
}

void AmrMesh::buildDomains(){
  CH_TIME("AmrMesh::buildDomains");
  if(m_verbosity > 5){
    pout() << "AmrMesh::buildDomains" << endl;
  }

  const int nlevels = 1 + m_maxAmrDepth;
  
  m_domains.resize(nlevels);
  m_dx.resize(nlevels);
  m_grids.resize(nlevels);


  m_dx[0] = (m_prob_hi[0] - m_probLo[0])/m_numCells[0];
  m_domains[0] = ProblemDomain(IntVect::Zero, m_numCells - IntVect::Unit);

  for (int lvl = 1; lvl <= m_maxAmrDepth; lvl++){
    m_dx[lvl]      = m_dx[lvl-1]/m_refinementRatios[lvl-1];
    m_domains[lvl] = m_domains[lvl-1];
    m_domains[lvl].refine(m_refinementRatios[lvl-1]);
  }
}

void AmrMesh::regrid(const Vector<IntVectSet>& a_tags,
		     const int a_lmin,
		     const int a_lmax,
		     const int a_regsize,
		     const int a_hardcap){
  CH_TIME("AmrMesh::regrid");
  if(m_verbosity > 1){
    pout() << "AmrMesh::regrid" << endl;
  }

  this->regridAmr(a_tags, a_lmin, a_lmax, a_hardcap);
  this->regridOperators(a_lmin, a_lmax, a_regsize);
}

void AmrMesh::regridAmr(const Vector<IntVectSet>& a_tags,
			const int a_lmin,
			const int a_lmax,
			const int a_hardcap){
  CH_TIME("AmrMesh::regridAmr(tags, level, level, hardcap)");
  if(m_verbosity > 1){
    pout() << "AmrMesh::regridAmr(tags, level, level, hardcap)" << endl;
  }

  // TLDR: This is the version that reads boxes. AmrMesh makes the grids from the tags and load balances them
  //       by using the patch volume. Those grids are then sent to the various realms. 
  
  Vector<IntVectSet> tags = a_tags; // buildGrids destroys tags, so copy them

  // This is stuff that always gets done
  this->buildGrids(tags, a_lmin, a_lmax, a_hardcap);

  // Define realms with the new grids and redo the realm stuff
  this->defineRealms();

  for (auto& r : m_realms){
    r.second->regrid_base(a_lmin);
  }
}

void AmrMesh::regridOperators(const int a_lmin,
			      const int a_lmax,
			      const int a_regsize){
  CH_TIME("AmrMesh::regridOperators(procs, boxes, level)");
  if(m_verbosity > 1){
    pout() << "AmrMesh::regridOperators(procs, boxes, level)" << endl;
  }
  
  for (auto& r : m_realms){
    r.second->regridOperators(a_lmin, a_lmax, a_regsize);
  }
}

void AmrMesh::buildGrids(Vector<IntVectSet>& a_tags, const int a_lmin, const int a_lmax, const int a_hardcap){
  CH_TIME("AmrMesh::buildGrids");
  if(m_verbosity > 2){
    pout() << "AmrMesh::buildGrids" << endl;
  }

  // TLDR: a_lmin is the coarsest level that changes. A special condition is a_lmin=0 for which we assume
  //       that there are no prior grids.
  //       a_lmax is the finest level that changes. This is typically 

  // base is the coarsest level which does not change. top_level is the finest level where we have tags. We should never
  // have tags on max_amr_depth, and we make that restriction here.
  const int base      = Max(0, a_lmin - 1);
  const int top_level = (m_finestLevel == m_maxAmrDepth) ? m_finestLevel - 1 : a_tags.size() - 1;

  // New and old grid boxes
  Vector<Vector<Box> > new_boxes(1 + top_level);
  Vector<Vector<Box> > old_boxes(1 + top_level);

  // Enforce potential hardcap. 
  const int hardcap = (a_hardcap == -1) ? m_maxAmrDepth : a_hardcap;

  // Inside this loop we make the boxes. 
  if(m_maxAmrDepth > 0 && hardcap > 0){
    domainSplit(m_domains[0], old_boxes[0], m_maxBoxSize, m_blockingFactor);

    if(!m_hasGrids){
      for (int lvl = 1; lvl <= top_level; lvl++){
	old_boxes[lvl].resize(0);
	old_boxes[lvl].push_back(m_domains[lvl].domainBox());
      }
    }
    else{
      for (int lvl = 0; lvl <= top_level; lvl++){ 
	old_boxes[lvl] = m_grids[lvl].boxArray();
      }
    }

    // Berger-Rigoutsos grid generation
    int new_finest_level;
    if(m_gridGenerationMethod == GridGenerationMethod::BergerRigoutsous){
      BRMeshRefine mesh_refine(m_domains[0], m_refinementRatios, m_fillRatioBR, m_blockingFactor, m_bufferSizeBR, m_maxBoxSize);
      new_finest_level = mesh_refine.regrid(new_boxes, a_tags, base, top_level, old_boxes);
    }
    else if (m_gridGenerationMethod == GridGenerationMethod::Tiled){
      TiledMeshRefine mesh_refine(m_domains[0], m_refinementRatios, m_blockingFactor*IntVect::Unit);
      new_finest_level = mesh_refine.regrid(new_boxes, a_tags, base, top_level, old_boxes);
    }
    else{
      MayDay::Abort("AmrMesh::regrid - logic bust, regridding with unknown regrid algorithm");
    }
    
    m_finestLevel = Min(new_finest_level, m_maxAmrDepth); // Don't exceed m_maxAmrDepth
    m_finestLevel = Min(m_finestLevel,   m_maxSimulationDepth); // Don't exceed maximum simulation depth
    m_finestLevel = Min(m_finestLevel,   hardcap);         // Don't exceed hardcap
  }
  else{ // Only end up here if we have a single grid level, i.e. just single-level grid decomposition. 
    new_boxes.resize(1);
    domainSplit(m_domains[0], new_boxes[0], m_maxBoxSize, m_blockingFactor);
    
    m_finestLevel = 0;
  }

  if(a_lmin == 0){ // Coarsest level also changes in this case, but that's not caught by the regridders. 
    domainSplit(m_domains[0], new_boxes[0], m_maxBoxSize, m_blockingFactor);
  }

  // Morton order the boxes.
  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    load_balance::sort(new_boxes[lvl], m_boxSort);
  }

  // Load balance boxes with patch volume as load proxy. 
  Vector<Vector<int> > pid(1 + m_finestLevel);
  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    load_balance::make_balance(pid[lvl], new_boxes[lvl]);
  }

  // Define grids. If a_lmin=0 every grid is new, otherwise keep old grids up to but not including a_lmin
  if(a_lmin == 0){
    m_grids.resize(1 + m_finestLevel);
    for (int lvl = 0; lvl <= m_finestLevel; lvl++){
      m_grids[lvl] = DisjointBoxLayout();
      m_grids[lvl].define(new_boxes[lvl], pid[lvl], m_domains[lvl]);
      m_grids[lvl].close();
    }
  }
  else{
    Vector<DisjointBoxLayout> old_grids = m_grids;
    m_grids.resize(1 + m_finestLevel);
    for (int lvl = 0; lvl < a_lmin; lvl++){ // Copy old grids
      m_grids[lvl] = old_grids[lvl];
    }
    for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++){ // Create new ones from tags
      m_grids[lvl] = DisjointBoxLayout();
      m_grids[lvl].define(new_boxes[lvl], pid[lvl], m_domains[lvl]);
      m_grids[lvl].close();
    }
  }

  m_hasGrids = true;
}

void AmrMesh::computeGradient(LevelData<EBCellFAB>&       a_gradient,
			      const LevelData<EBCellFAB>& a_phi,
			      const std::string           a_realm,
			      const phase::which_phase    a_phase,
			      const int                   a_lvl){
  CH_TIME("AmrMesh::computeGradient(grad, phi, realm, phase, lvl)");
  if(m_verbosity > 5){
    pout() << "AmrMesh::computeGradient(grad, phi, realm, phase,lvl)" << endl;
  }

  CH_assert(a_phi.nComp()      == 1);
  CH_assert(a_gradient.nComp() == SpaceDim);

  if(!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::computeGradient(grad, phi, realm, phase, lvl) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }
    
  const int comp  = 0;
  const int ncomp = 1;
    
  const Real& dx = m_realms[a_realm]->getDx()[a_lvl];
  const DisjointBoxLayout& dbl = m_realms[a_realm]->getGrids()[a_lvl];
  const ProblemDomain& domain  = m_realms[a_realm]->getDomains()[a_lvl];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    EBCellFAB& grad        = a_gradient[dit()];
    const EBCellFAB& phi   = a_phi[dit()];
    const EBISBox& ebisbox = phi.getEBISBox();
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    const Box& region      = dbl.get(dit());

    // For interior cells we do our old friend centered differences. God I hate Chombo Fortran.
    const BaseFab<Real>& phi_fab = phi.getSingleValuedFAB();
    BaseFab<Real>& grad_fab  = grad.getSingleValuedFAB();
    FORT_GRADIENT(CHF_FRA(grad_fab),
		  CHF_CONST_FRA1(phi_fab, comp),
		  CHF_CONST_REAL(dx),
		  CHF_BOX(region));

    BaseIVFAB<VoFStencil>& grad_stencils = (*m_realms[a_realm]->get_gradsten(a_phase)[a_lvl])[dit()];

    for (VoFIterator vofit(grad_stencils.getIVS(), ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof    = vofit();
      const VoFStencil& sten = grad_stencils(vof, 0);

      for (int dir = 0; dir < SpaceDim; dir++){
	grad(vof, dir) = 0.0;
      }

      for (int i = 0; i < sten.size(); i++){
	const VolIndex& ivof = sten.vof(i);
	const Real& iweight  = sten.weight(i);
	const int ivar       = sten.variable(i); // Note: For the gradient stencil the ivar is the direction. 

	grad(vof, ivar) += phi(ivof, comp)*iweight;
      }
    }

    // Set covered to zero
    for (int dir= 0; dir < SpaceDim; dir++){
      grad.setCoveredCellVal(0.0, dir);
    }
  }
}

void AmrMesh::computeGradient(EBAMRCellData&           a_gradient,
			      const EBAMRCellData&     a_phi,
			      const std::string        a_realm,  
			      const phase::which_phase a_phase){
  CH_TIME("AmrMesh::computeGradient(grad, phi, realm, phase)");
  if(m_verbosity > 5){
    pout() << "AmrMesh::computeGradient(grad, phi, realm, phase)" << endl;
  }

  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    this->computeGradient(*a_gradient[lvl], *a_phi[lvl], a_realm, a_phase, lvl);
  }
}

void AmrMesh::computeGradient(MFAMRCellData& a_gradient, const MFAMRCellData& a_phi, const std::string a_realm){
  CH_TIME("AmrMesh::computeGradient(mf grad, mf phi, realm)");
  if(m_verbosity > 5){
    pout() << "AmrMesh::computeGradient(mf grad, mf phi, realm)" << endl;
  }

  for (int iphase = 0; iphase < m_multifluidIndexSpace->num_phases(); iphase++){
    EBAMRCellData alias_grad(1 + m_finestLevel);
    EBAMRCellData alias_phi(1 + m_finestLevel);

    for (int lvl = 0; lvl <= m_finestLevel; lvl++){
      alias_grad[lvl] = RefCountedPtr<LevelData<EBCellFAB> > (new LevelData<EBCellFAB>());
      alias_phi[lvl]  = RefCountedPtr<LevelData<EBCellFAB> > (new LevelData<EBCellFAB>());
      
      mfalias::aliasMF(*alias_grad[lvl], iphase, *a_gradient[lvl]);
      mfalias::aliasMF(*alias_phi[lvl],  iphase, *a_phi[lvl]);
    }

    if(iphase == 0){
      this->computeGradient(alias_grad, alias_phi, a_realm, phase::gas);
    }
    else if(iphase == 1){
      this->computeGradient(alias_grad, alias_phi, a_realm, phase::solid);
    }
  }
}

void AmrMesh::averageDown(EBAMRCellData& a_data, const std::string a_realm, const phase::which_phase a_phase){
  CH_TIME("AmrMesh::averageDown(ebamrcell, realm, phase");
  if(m_verbosity > 3){
    pout() << "AmrMesh::averageDown(ebamrcell, realm, phase)" << endl;
  }

  if(!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::averageDown(ebamrcell, realm, phase) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }
  
  for (int lvl = m_finestLevel; lvl > 0; lvl--){
    ebcoarseaverage& aveOp = *m_realms[a_realm]->getCoarseAverage(a_phase)[lvl];
      
    const int ncomps = a_data[lvl]->nComp();
    const Interval interv (0, ncomps-1);

    aveOp.average(*a_data[lvl-1], *a_data[lvl], interv);
  }
    
  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    a_data[lvl]->exchange();
  }
}

void AmrMesh::averageDown(EBAMRCellData& a_data, const std::string a_realm, const phase::which_phase a_phase, const int a_lvl){
  CH_TIME("AmrMesh::averageDown(ebamrcelldata, realm, phase, level");
  if(m_verbosity > 3){
    pout() << "AmrMesh::averageDown(ebamrcelldata, realm, phase, level)" << endl;
  }

  if(!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::averageDown(ebamrcell, realm, phase, level) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }
  
  const int ncomps = a_data[a_lvl]->nComp();
  const Interval interv (0, ncomps-1);

  ebcoarseaverage& aveOp = *m_realms[a_realm]->getCoarseAverage(a_phase)[a_lvl+1];

  aveOp.average(*a_data[a_lvl], *a_data[a_lvl+1], interv);

  a_data[a_lvl]->exchange();
}

void AmrMesh::averageDown(MFAMRFluxData& a_data, const std::string a_realm){
  CH_TIME("AmrMesh::averageDown(mfamrflux, realm)");
  if(m_verbosity > 3){
    pout() << "AmrMesh::averageDown(mfamrflux, realm)" << endl;
  }

  if(!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::averageDown(mfamrflux, realm) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const RefCountedPtr<EBIndexSpace>& ebis_gas = m_realms[a_realm]->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace>& ebis_sol = m_realms[a_realm]->get_ebis(phase::solid);

  // Alias the data to regular EBFluxFABs
  EBAMRFluxData alias_g(1 + m_finestLevel);
  EBAMRFluxData alias_s(1 + m_finestLevel);
    
  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    alias_g[lvl] = RefCountedPtr<LevelData<EBFluxFAB> > (new LevelData<EBFluxFAB>());
    alias_s[lvl] = RefCountedPtr<LevelData<EBFluxFAB> > (new LevelData<EBFluxFAB>());
      
    if(!ebis_gas.isNull()) mfalias::aliasMF(*alias_g[lvl], phase::gas,   *a_data[lvl]);
    if(!ebis_sol.isNull()) mfalias::aliasMF(*alias_s[lvl], phase::solid, *a_data[lvl]);
  }

  if(!ebis_gas.isNull()) this->averageDown(alias_g, a_realm, phase::gas);
  if(!ebis_sol.isNull()) this->averageDown(alias_s, a_realm, phase::solid);
}

void AmrMesh::averageDown(MFAMRCellData& a_data, const std::string a_realm){
  CH_TIME("AmrMesh::averageDown(mfamrcell, realm)");
  if(m_verbosity > 3){
    pout() << "AmrMesh::averageDown(mfamrcell, realm)" << endl;
  }

  if(!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::averageDown(mfamrcell, realm) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const RefCountedPtr<EBIndexSpace>& ebis_gas = m_realms[a_realm]->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace>& ebis_sol = m_realms[a_realm]->get_ebis(phase::solid);
  
  EBAMRCellData alias_g(1 + m_finestLevel);
  EBAMRCellData alias_s(1 + m_finestLevel);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    alias_g[lvl] = RefCountedPtr<LevelData<EBCellFAB> > (new LevelData<EBCellFAB>());
    alias_s[lvl] = RefCountedPtr<LevelData<EBCellFAB> > (new LevelData<EBCellFAB>());
      
    if(!ebis_gas.isNull()) mfalias::aliasMF(*alias_g[lvl], phase::gas,   *a_data[lvl]);
    if(!ebis_sol.isNull()) mfalias::aliasMF(*alias_s[lvl], phase::solid, *a_data[lvl]);
  }

  if(!ebis_gas.isNull()) this->averageDown(alias_g, a_realm, phase::gas);
  if(!ebis_sol.isNull()) this->averageDown(alias_s, a_realm, phase::solid);
}

void AmrMesh::averageDown(EBAMRFluxData& a_data, const std::string a_realm, const phase::which_phase a_phase){
  CH_TIME("AmrMesh::averageDown(ebamrflux, realm, phase");
  if(m_verbosity > 3){
    pout() << "AmrMesh::averageDown(ebamrflux, realm, phase)" << endl;
  }

  if(!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::averageDown(ebamrflux, realm, phase) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  for (int lvl = m_finestLevel; lvl > 0; lvl--){
    const int ncomps = a_data[lvl]->nComp();
    const Interval interv (0, ncomps-1);

    ebcoarseaverage& aveOp = *m_realms[a_realm]->getCoarseAverage(a_phase)[lvl];
    aveOp.average(*a_data[lvl-1], *a_data[lvl], interv);
  }

  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    a_data[lvl]->exchange();
  }
}

void AmrMesh::averageDown(EBAMRIVData& a_data, const std::string a_realm, const phase::which_phase a_phase){
  CH_TIME("AmrMesh::averageDown(ebamriv, realm, phase)");
  if(m_verbosity > 3){
    pout() << "AmrMesh::averageDown(ebamriv, realm, phase)" << endl;
  }

  if(!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::averageDown(ebamriv, realm, phase) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  for (int lvl = m_finestLevel; lvl > 0; lvl--){
    const int ncomps = a_data[lvl]->nComp();
    const Interval interv (0, ncomps-1);

    ebcoarseaverage& aveOp = *m_realms[a_realm]->getCoarseAverage(a_phase)[lvl];
    aveOp.average(*a_data[lvl-1], *a_data[lvl], interv);
  }

  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    a_data[lvl]->exchange();
  }
}

void AmrMesh::conservativeAverage(EBAMRIVData& a_data, const std::string a_realm, const phase::which_phase a_phase){
  CH_TIME("AmrMesh::conservativeAverage(ebamriv, realm, phase)");
  if(m_verbosity > 3){
    pout() << "AmrMesh::conservativeAverage(ebamriv, realm, phase)" << endl;
  }

  if(!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::conservativeAverage(ebamriv, realm, phase) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  for (int lvl = m_finestLevel; lvl > 0; lvl--){
    const int ncomps = a_data[lvl]->nComp();
    const Interval interv (0, ncomps-1);

    ebcoarseaverage& aveOp = *m_realms[a_realm]->getCoarseAverage(a_phase)[lvl];
    
    aveOp.conservativeAverage(*a_data[lvl-1], *a_data[lvl], interv);
  }

  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    a_data[lvl]->exchange();
  }
}

void AmrMesh::interpGhost(EBAMRCellData& a_data, const std::string a_realm, const phase::which_phase a_phase){
  CH_TIME("AmrMesh::interpGhost(ebamrcell, realm, phase)");
  if(m_verbosity > 3){
    pout() << "AmrMesh::interpGhost(ebamrcell, realm, phase)" << endl;
  }

  if(!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::interpGhost(ebamrcell, realm, phase) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }  

  if(m_ghostCellInterpolationMethod == GhostInterpolation::PiecewiseLinear){
    this->interpGhostPwl(a_data, a_realm, a_phase);
  }
  else if(m_ghostCellInterpolationMethod == GhostInterpolation::Quadratic){
    this->interpGhostQuad(a_data, a_realm, a_phase);
  }
  else{
    MayDay::Abort("AmrMesh::interpGhost - unsupported interpolation type requested");
  }
}

void AmrMesh::interpGhost(LevelData<EBCellFAB>&       a_fineData,
			  const LevelData<EBCellFAB>& a_coarData,
			  const int                   a_fineLevel,
			  const std::string           a_realm,
			  const phase::which_phase    a_phase){
  CH_TIME("AmrMesh::interpGhost(fine, coar, level, realm, phase)");
  if(m_verbosity > 3){
    pout() << "AmrMesh::interpGhost(fine, coar, level, realm, phase)" << endl;
  }

  if(!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::interpGhost(fine, coar, level, realm, phase) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  if(a_fineLevel > 0){

    const int ncomps      = a_fineData.nComp();
    const Interval interv = Interval(0, ncomps-1);
    
    if(m_ghostCellInterpolationMethod == GhostInterpolation::PiecewiseLinear){
      AggEBPWLFillPatch& fillpatch = *m_realms[a_realm]->getFillPatch(a_phase)[a_fineLevel];
    
      fillpatch.interpolate(a_fineData, a_coarData, a_coarData, 0.0, 0.0, 0.0, interv);
    }
    else if(m_ghostCellInterpolationMethod == GhostInterpolation::Quadratic){
      nwoebquadcfinterp& quadcfi = *m_realms[a_realm]->getNWOEBQuadCFInterp(a_phase)[a_fineLevel];
      quadcfi.coarseFineInterp(a_fineData, a_coarData, 0, 0, ncomps);
    }
    else{
      MayDay::Abort("AmrMesh::interpGhost - unsupported interpolation type requested");
    }
  }
}

void AmrMesh::interpGhost(MFAMRCellData& a_data, const std::string a_realm){
  CH_TIME("AmrMesh::interpGhost(mfamrcell, realm)");
  if(m_verbosity > 3){
    pout() << "AmrMesh::interpGhost(mfamrcell, realm)" << endl;
  }

  if(!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::interpGhost(mfamrcell, realm, phase) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  // Do aliasing
  EBAMRCellData alias_g(1 + m_finestLevel);
  EBAMRCellData alias_s(1 + m_finestLevel);

  const RefCountedPtr<EBIndexSpace>& ebis_gas = m_realms[a_realm]->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace>& ebis_sol = m_realms[a_realm]->get_ebis(phase::solid);
  
  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    alias_g[lvl] = RefCountedPtr<LevelData<EBCellFAB> > (new LevelData<EBCellFAB>());
    alias_s[lvl] = RefCountedPtr<LevelData<EBCellFAB> > (new LevelData<EBCellFAB>());
      
    if(!ebis_gas.isNull()) mfalias::aliasMF(*alias_g[lvl], phase::gas,   *a_data[lvl]);
    if(!ebis_sol.isNull()) mfalias::aliasMF(*alias_s[lvl], phase::solid, *a_data[lvl]);
  }

  if(!ebis_gas.isNull()) this->interpGhost(alias_g, a_realm, phase::gas);
  if(!ebis_sol.isNull()) this->interpGhost(alias_s, a_realm, phase::solid);
}

void AmrMesh::interpGhostQuad(EBAMRCellData& a_data, const std::string a_realm, const phase::which_phase a_phase){
  CH_TIME("AmrMesh::interpGhostQuad(ebamrcell, realm, phase)");
  if(m_verbosity > 3){
    pout() << "AmrMesh::interpGhostQuad(ebamrcell, realm, phase)" << endl;
  }

  if(!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::interpGhostQuad(ebamrcell, realm, phase) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  for (int lvl = m_finestLevel; lvl > 0; lvl--){
    const int ncomps = a_data[lvl]->nComp();
    const Interval interv(0, ncomps -1);

    nwoebquadcfinterp& quadcfi = *m_realms[a_realm]->getNWOEBQuadCFInterp(a_phase)[lvl];

    quadcfi.coarseFineInterp(*a_data[lvl], *a_data[lvl-1], 0, 0, ncomps);
  }

  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    a_data[lvl]->exchange();
  }
}

void AmrMesh::interpGhostPwl(EBAMRCellData& a_data, const std::string a_realm, const phase::which_phase a_phase){
  CH_TIME("AmrMesh::interpGhostPwl(ebamrcell, realm, phase)");
  if(m_verbosity > 3){
    pout() << "AmrMesh::interpGhostPwl(ebamrcell, realm, phase)" << endl;
  }

  if(!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::interpGhostPwl(ebamrcell, realm, phase) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }
  
  for (int lvl = m_finestLevel; lvl > 0; lvl--){
    const int ncomps = a_data[lvl]->nComp();
    const Interval interv(0, ncomps -1);

    AggEBPWLFillPatch& fillpatch = *m_realms[a_realm]->getFillPatch(a_phase)[lvl];
    
    fillpatch.interpolate(*a_data[lvl], *a_data[lvl-1], *a_data[lvl-1], 0.0, 0.0, 0.0, interv);
  }

  for (int lvl = 0; lvl <= m_finestLevel; lvl++){
    a_data[lvl]->exchange();
  }
}

void AmrMesh::interpToCentroids(EBAMRCellData& a_data, const std::string a_realm, const phase::which_phase a_phase){
  CH_TIME("AmrMesh::interpToCentroids(ebamrcell, realm, phase)");
  if(m_verbosity > 3){
    pout() << "AmrMesh::interpToCentroids(ebamrcell, realm, phase)" << endl;
  }

  if(!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::interpToCentroids(ebamrcell, realm, phase) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }
  
  irreg_amr_stencil<CentroidInterpolationStencil>& stencil = m_realms[a_realm]->getCentroidInterpolationStencils(a_phase);
  stencil.apply(a_data);
}

void AmrMesh::parseVerbosity(){
  CH_TIME("AmrMesh::parseVerbosity");

  ParmParse pp("AmrMesh");
  pp.get("verbosity", m_verbosity);
}

void AmrMesh::parseCoarsestLevelNumCells(){

  ParmParse pp("AmrMesh");
  Vector<int> cells;
  cells.resize(pp.countval("coarsest_domain"));
  CH_assert(cells.size() >= SpaceDim);
  pp.getarr("coarsest_domain", cells, 0, SpaceDim);

  m_numCells = IntVect(D_DECL(cells[0], cells[1], cells[2]));
}

void AmrMesh::parseMaxAmrDepth(){

  ParmParse pp("AmrMesh");
  int depth;
  pp.get("max_amr_depth", depth);
  if(depth >= 0){
    m_maxAmrDepth = depth;
  }
  else{
    m_maxAmrDepth = 0;
  }
}

void AmrMesh::parseMaxSimulationDepth(){

  ParmParse pp("AmrMesh");
  int depth;
  pp.get("max_sim_depth", depth);
  if(depth >= 0){
    m_maxSimulationDepth = depth;
  }
  else {
    m_maxSimulationDepth = m_maxAmrDepth;
  }
}

void AmrMesh::parseEbCf(){
  ParmParse pp("AmrMesh");
  std::string str;
  pp.get("ebcf", str);
  if(str == "true"){
    m_hasEbCf = true;
  }
  else if(str == "false"){
    m_hasEbCf = false;
  }
}

void AmrMesh::parseRefinementRatios(){

  ParmParse pp("AmrMesh");
  Vector<int> ratios;
  ratios.resize(pp.countval("ref_rat"));
  pp.getarr("ref_rat", ratios, 0, ratios.size());

  // Pad with 2 if user didn't supply enough
  while(ratios.size() < m_maxAmrDepth){
    ratios.push_back(2);
  }
  
  m_refinementRatios = ratios;
}

void AmrMesh::setRefinementRatios(const Vector<int> a_ref_ratios){
  m_refinementRatios = a_ref_ratios;
}

void AmrMesh::parseBrFillRatio(){
  ParmParse pp("AmrMesh");
  Real fill_ratio = 1.0;
  pp.get("fill_ratio", fill_ratio);
  if(fill_ratio > 0.0 && fill_ratio <= 1.0){
    m_fillRatioBR = fill_ratio;
  }
}

void AmrMesh::setFinestLevel(const int a_finest_level){
  m_finestLevel = a_finest_level;
  m_finestLevel = Min(m_finestLevel, m_maxAmrDepth); // Don't exceed m_maxAmrDepth
  m_finestLevel = Min(m_finestLevel, m_maxSimulationDepth); // Don't exceed maximum simulation depth
}

void AmrMesh::setGrids(const Vector<Vector<Box> >& a_boxes, const std::map<std::string, Vector<Vector<long int> > >& a_realms_and_loads){
  CH_TIME("AmrMesh::setGrids(boxes, loads, regsize)");
  if(m_verbosity > 3){
    pout() << "AmrMesh::setGrids(boxes, loads, regsize)" << endl;
  }

  const int lmin = 0;

  for (const auto& r : a_realms_and_loads){
    const std::string&               cur_realm = r.first;
    const Vector<Vector<long int> >& cur_loads = r.second;

    Vector<Vector<int> > pids(1 + m_finestLevel);

    // Do load balancing. 
    for (int lvl = 0; lvl <= m_finestLevel; lvl++){
      load_balance::make_balance(pids[lvl], cur_loads[lvl], a_boxes[lvl]);
    }

    this->regridRealm(cur_realm, pids, a_boxes, lmin);
  }

  // Set the proxy grids, too. These are load balanced using the patch volume. 
  m_grids = m_realms[realm::primal]->getGrids();
  m_hasGrids = true;
}

void AmrMesh::parseMaxBoxSize(){
  ParmParse pp("AmrMesh");
  int box_size;
  pp.get("max_box_size", box_size);
  if(box_size >= 8 && box_size % 2 == 0){
    m_maxBoxSize = box_size;
  }
  else{
    MayDay::Abort("AmrMesh::parseMaxBoxSize - must have box_size >= 8 and divisible by 2");
  }
}

void AmrMesh::parseMaxEbisBoxSize(){

  ParmParse pp("AmrMesh");
  int box_size;
  pp.get("max_ebis_box", box_size);
  if(box_size >= 8 && box_size % 2 == 0){
    m_maxEbisBoxSize = box_size;
  }
  else{
    MayDay::Abort("AmrMesh::parseMaxEbisBoxSize - must have box_size >= 8 and divisible by 2");
  }
}

void AmrMesh::parseBrBufferSize(){

  ParmParse pp("AmrMesh");
  int buffer = 2;
  pp.get("buffer_size", buffer);
  if(buffer > 0){
    m_bufferSizeBR = buffer;
  }
}

void AmrMesh::parsegridGeneration(){

  ParmParse pp("AmrMesh");
  std::string str;
  pp.get("grid_algorithm", str);
  if(str == "br"){
    m_gridGenerationMethod = GridGenerationMethod::BergerRigoutsous;
  }
  else if(str == "tiled"){
    m_gridGenerationMethod = GridGenerationMethod::Tiled;
  }
  else{
    MayDay::Abort("AmrMesh::parsegridGeneration - unknown grid generation method requested");
  }

  pp.get("BoxSorting", str);
  if( str == "none"){
    m_boxSort = BoxSorting::none;
  }
  if( str == "std"){
    m_boxSort = BoxSorting::std;
  }
  else if(str == "shuffle"){
    m_boxSort = BoxSorting::shuffle;
  }
  else if(str == "morton"){
    m_boxSort = BoxSorting::morton;
  }
  else {
    MayDay::Abort("AmrMesh::parsegridGeneration - unknown box sorting method requested");
  }
}

void AmrMesh::parseIrregTagGrowth(){
  ParmParse pp("AmrMesh");
  pp.get("irreg_growth", m_irregTagGrowth);

  m_irregTagGrowth = std::max(m_irregTagGrowth, 0);
}

void AmrMesh::parseBlockingFactor(){
  ParmParse pp("AmrMesh");
  int blocking;
  pp.get("blocking_factor", blocking);
  if(blocking >= 4 && blocking % 2 == 0){
    m_blockingFactor = blocking;
  }
}

void AmrMesh::parseEbGhostCells(){

  ParmParse pp("AmrMesh");
  int ebghost;
  pp.get("eb_ghost", ebghost);
  if(ebghost >= 2){
    m_numEbGhostsCells = ebghost;
  }
}

void AmrMesh::parseNumGhostCells(){
  CH_TIME("AmrMesh::parseNumGhostCells");

  ParmParse pp("AmrMesh");
  int ghost;
  pp.get("num_ghost", m_numGhostCells);
  pp.get("lsf_ghost", m_numLsfGhostCells);
}

void AmrMesh::parseRedistributionRadius(){
  ParmParse pp("AmrMesh");
  int rad = 1;
  pp.get("redist_radius", rad);
  m_redistributionRadius = rad;
}

void AmrMesh::parseCentroidStencils(){
  std::string str = "taylor";
  ParmParse pp("AmrMesh");
  pp.get("centroid_sten", str);


  // Maybe, in the future, we can change these but the user should not care about these (yet)
  m_centroidStencilRadius = 1;
  m_centroidStencilOrder  = 1;

  if(str == "linear"){
    m_centroidStencilType = IrregStencil::StencilType::Linear;
  }
  else if(str == "taylor"){
    m_centroidStencilType = IrregStencil::StencilType::TaylorExtrapolation;
  }
  else if(str == "lsq"){
    m_centroidStencilType = IrregStencil::StencilType::LeastSquares;
  }
  else if(str == "pwl"){
    m_centroidStencilType = IrregStencil::StencilType::PiecewiseLinear;
  }
  else{
    MayDay::Abort("AmrMesh::parseCentroidStencils - unknown stencil requested");
  }
}

void AmrMesh::parseEbCentroidStencils(){
  std::string str = "taylor";
  ParmParse pp("AmrMesh");
  pp.get("eb_sten", str);
  
  // Maybe, in the future, we can change these but the user should not care about these (yet)
  m_ebCentroidStencilRadius = 1;
  m_ebCentroidStencilOrder  = 1;

  if(str == "linear"){
    m_ebCentroidStencilType = IrregStencil::StencilType::Linear;
  }
  else if(str == "taylor"){
    m_ebCentroidStencilType = IrregStencil::StencilType::TaylorExtrapolation;
  }
  else if(str == "lsq"){
    m_ebCentroidStencilType = IrregStencil::StencilType::LeastSquares;
  }
  else if(str == "pwl"){
    m_ebCentroidStencilType = IrregStencil::StencilType::PiecewiseLinear;
  }
  else{
    MayDay::Abort("AmrMesh::parseEbCentroidStencils - unknown stencil requested");
  }
}

void AmrMesh::sanityCheck() const {
  CH_TIME("AmrMesh::sanityCheck");
  if(m_verbosity > 1){
    pout() << "AmrMesh::sanityCheck" << endl;
  }

  CH_assert(m_maxAmrDepth >= 0);
  for (int lvl = 0; lvl < m_refinementRatios.size(); lvl++){
    CH_assert(m_refinementRatios[lvl] == 2 || m_refinementRatios[lvl] == 4);
    CH_assert(m_blockingFactor >= 4 && m_blockingFactor % m_refinementRatios[lvl] == 0);
  }
  CH_assert(m_maxBoxSize >= 8 && m_maxBoxSize % m_blockingFactor == 0);
  CH_assert(m_fillRatioBR > 0. && m_fillRatioBR <= 1.0);
  CH_assert(m_bufferSizeBR > 0);
}

bool AmrMesh::getEbCf() const {
  return m_hasEbCf;
}

RealVect AmrMesh::getProbLo() const {
  return m_probLo;
}

RealVect AmrMesh::getProbHi() const {
  return m_prob_hi;
}

int AmrMesh::getFinestLevel() const {
  return m_finestLevel;
}

int AmrMesh::getIrregTagGrowth() const {
  return m_irregTagGrowth;
}

int AmrMesh::getMaxAmrDepth() const {
  return m_maxAmrDepth;
}

int AmrMesh::getMaxSimulationDepth() const {
  return m_maxSimulationDepth;
}

int AmrMesh::getNumberOfGhostCells() const {
  return m_numGhostCells;
}

int AmrMesh::getNumberOfEbGhostCells() const {
  return m_numEbGhostsCells;
}

int AmrMesh::getRedistributionRadius() const {
  return m_redistributionRadius;
}

int AmrMesh::getBlockingFactor() const {
  return m_blockingFactor;
}

int AmrMesh::getMaxBoxSize() const {
  return m_maxBoxSize;
}

int AmrMesh::getBrBuffer() const {
  return m_bufferSizeBR;
}

int AmrMesh::getMaxEbisBoxSize() const {
  return m_maxEbisBoxSize;
}

int AmrMesh::getRefinementRatio(const int a_level1, const int a_level2) const {
  int coarLevel = Min(a_level1, a_level2);
  int fineLevel = Max(a_level1, a_level2);

  int ref = 1;
  for (int lvl = coarLevel; lvl < fineLevel; lvl++){
    ref = ref*m_refinementRatios[lvl];
  }

  return ref;
}

ProblemDomain AmrMesh::getFinestDomain() const {
  return m_domains[m_maxAmrDepth];
}

Real AmrMesh::getFinestDx() const {
  return m_dx[m_maxAmrDepth];
}

const Vector<Real>& AmrMesh::getDx() const {
  return m_dx;
}

const Vector<int>& AmrMesh::getRefinementRatios() const {
  return m_refinementRatios;
}

const RefCountedPtr<BaseIF>& AmrMesh::getBaseImplicitFunction(const phase::which_phase a_phase) const {
  return m_baseif.at(a_phase);
}

Vector<IntVectSet> AmrMesh::getIrregularTags() const {
  CH_TIME("AmrMesh::getIrregularTags");
  if(m_verbosity > 5){
    pout() << "AmrMesh::getIrregularTags" << endl;
  }

  Vector<IntVectSet> tags(m_maxAmrDepth);

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_multifluidIndexSpace->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_multifluidIndexSpace->get_ebis(phase::solid);

  CH_assert(ebis_gas != NULL);

  for (int lvl = 0; lvl < m_maxAmrDepth; lvl++){ // Don't need tags on maxdepth, we will never generate grids below that.
    const int which_level = ebis_gas->getLevel(m_domains[lvl]);

    tags[lvl] |= ebis_gas->irregCells(which_level);
    if(!ebis_sol.isNull()){
      tags[lvl] |= ebis_sol->irregCells(which_level);
    }
  }

  return tags;
}

const Vector<ProblemDomain>& AmrMesh::getDomains() const {
  return m_domains;
}

const Vector<DisjointBoxLayout>& AmrMesh::getProxyGrids() const {
  return m_grids;
}

const Vector<DisjointBoxLayout>& AmrMesh::getGrids(const std::string a_realm) const {
  if(!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::getGrids - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }
  
  return m_realms[a_realm]->getGrids();
}

const Vector<EBISLayout>& AmrMesh::getEBISLayout(const std::string a_realm, const phase::which_phase a_phase) const{
  return m_realms[a_realm]->getEBISLayout(a_phase);
}

Vector<RefCountedPtr<LayoutData<VoFIterator> > >& AmrMesh::getVofIterator(const std::string a_realm, const phase::which_phase a_phase){
  return m_realms[a_realm]->getVofIterator(a_phase);
}

const Vector<RefCountedPtr<LayoutData<Vector<LayoutIndex> > > >& AmrMesh::getNeighbors(const std::string a_realm,
										       const phase::which_phase a_phase) const{
  return m_realms[a_realm]->getNeighbors(a_phase);
}

const AMRMask& AmrMesh::getMask(const std::string a_mask, const int a_buffer, const std::string a_realm) const {
  return m_realms[a_realm]->getMask(a_mask, a_buffer);
}

const Vector<RefCountedPtr<EBLevelGrid> >& AmrMesh::getEBLevelGrid(const std::string a_realm, const phase::which_phase a_phase) const{
  return m_realms[a_realm]->getEBLevelGrid(a_phase);
}

const Vector<RefCountedPtr<MFLevelGrid> >& AmrMesh::getMFLevelGrid(const std::string a_realm) const {
  return m_realms[a_realm]->getMFLevelGrid();
}

const EBAMRFAB& AmrMesh::getLevelset(const std::string a_realm, const phase::which_phase a_phase) const {
  return m_realms[a_realm]->getLevelset(a_phase);
}

Vector<RefCountedPtr<ebcoarseaverage> >& AmrMesh::getCoarseAverage(const std::string a_realm, const phase::which_phase a_phase){
  return m_realms[a_realm]->getCoarseAverage(a_phase);
}

Vector<RefCountedPtr<EBGhostCloud> >& AmrMesh::getGhostCloud(const std::string a_realm, const phase::which_phase a_phase){
  return m_realms[a_realm]->getGhostCloud(a_phase);
}

Vector<RefCountedPtr<nwoebquadcfinterp> >& AmrMesh::getNWOEBQuadCFInterp(const std::string a_realm, const phase::which_phase a_phase){
  return m_realms[a_realm]->getNWOEBQuadCFInterp(a_phase);
}

Vector<RefCountedPtr<EBQuadCFInterp> >& AmrMesh::getEBQuadCFInterp(const std::string a_realm, const phase::which_phase a_phase){
  return m_realms[a_realm]->getEBQuadCFInterp(a_phase);
}

Vector<RefCountedPtr<AggEBPWLFillPatch> >& AmrMesh::getFillPatch(const std::string a_realm, const phase::which_phase a_phase){
  return m_realms[a_realm]->getFillPatch(a_phase);
}

Vector<RefCountedPtr<EBPWLFineInterp> >& AmrMesh::getPwlInterpolator(const std::string a_realm, const phase::which_phase a_phase){
  return m_realms[a_realm]->getPwlInterpolator(a_phase);
}

Vector<RefCountedPtr<EBMGInterp> >& AmrMesh::getEBMGInterp(const std::string a_realm, const phase::which_phase a_phase){
  return m_realms[a_realm]->getEBMGInterp(a_phase);
}

Vector<RefCountedPtr<EBFluxRegister> >&  AmrMesh::getFluxRegister(const std::string a_realm, const phase::which_phase a_phase){
  return m_realms[a_realm]->getFluxRegister(a_phase);
}

Vector<RefCountedPtr<EBLevelRedist> >& AmrMesh::getLevelRedist(const std::string a_realm, const phase::which_phase a_phase){
  return m_realms[a_realm]->getLevelRedist(a_phase);
}

Vector<RefCountedPtr<EBCoarToFineRedist> >&  AmrMesh::getCoarToFineRedist(const std::string        a_realm,
									  const phase::which_phase a_phase){
  return m_realms[a_realm]->getCoarToFineRedist(a_phase);
}

Vector<RefCountedPtr<EBCoarToCoarRedist> >&  AmrMesh::getCoarToCoarRedist(const std::string        a_realm,
									  const phase::which_phase a_phase){
  return m_realms[a_realm]->getCoarToCoarRedist(a_phase);
}

Vector<RefCountedPtr<EBFineToCoarRedist> >&  AmrMesh::getFineToCoarRedist(const std::string        a_realm,
									  const phase::which_phase a_phase){
  return m_realms[a_realm]->getFineToCoarRedist(a_phase);
}

const irreg_amr_stencil<CentroidInterpolationStencil>& AmrMesh::getCentroidInterpolationStencils(const std::string        a_realm,
										    const phase::which_phase a_phase) const {
  return m_realms[a_realm]->getCentroidInterpolationStencils(a_phase);
}

const irreg_amr_stencil<EbCentroidInterpolationStencil>& AmrMesh::getEbCentroidInterpolationStencilStencils(const std::string        a_realm,
											 const phase::which_phase a_phase) const {
  return m_realms[a_realm]->getEbCentroidInterpolationStencilStencils(a_phase);
}

const irreg_amr_stencil<NonConservativeDivergenceStencil>& AmrMesh::getNonConservativeDivergenceStencils(const std::string a_realm, const phase::which_phase a_phase) const{
  return m_realms[a_realm]->getNonConservativeDivergenceStencils(a_phase);
}

Vector<RefCountedPtr<Copier> >& AmrMesh::getCopier(const std::string a_realm, const phase::which_phase a_phase){
  return m_realms[a_realm]->getCopier(a_phase);
}

Vector<RefCountedPtr<Copier> >& AmrMesh::getReverseCopier(const std::string a_realm, const phase::which_phase a_phase){
  return m_realms[a_realm]->getReverseCopier(a_phase);
}

bool AmrMesh::queryRealm(const std::string a_realm) const {
  CH_TIME("AmrMesh::queryRealm");
  if(m_verbosity > 5){
    pout() << "AmrMesh::queryRealm" << endl;
  }

  bool ret = true;
  
  if(m_realms.find(a_realm) == m_realms.end()){
    ret = false;
  }

  return ret;
}

void AmrMesh::registerRealm(const std::string a_realm){
  CH_TIME("AmrMesh::registerRealm");
  if(m_verbosity > 5){
    pout() << "AmrMesh::registerRealm" << endl;
  }

  if(!this->queryRealm(a_realm)){
    m_realms.emplace(a_realm, RefCountedPtr<realm> (new realm()));
  }
}

void AmrMesh::registerOperator(const std::string a_operator, const std::string a_realm, const phase::which_phase a_phase){
  CH_TIME("AmrMesh::registerOperator(operator, realm, phase)");
  if(m_verbosity > 5){
    pout() << "AmrMesh::registerOperator(operator, realm, phase)" << endl;
  }

  if(!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::registerOperator(operator, realm, phase) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  m_realms[a_realm]->registerOperator(a_operator, a_phase);
}

void AmrMesh::registerMask(const std::string a_mask, const int a_buffer, const std::string a_realm){
  CH_TIME("AmrMesh::registerMask(mask, realm, buffer)");
  if(m_verbosity > 5){
    pout() << "AmrMesh::registerMask(mask, realm, buffer)" << endl;
  }

  if(!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::registerMask(mask, realm, buffer) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  m_realms[a_realm]->registerMask(a_mask, a_buffer);
}

void AmrMesh::defineRealms(){
  CH_TIME("AmrMesh::defineRealms()");
  if(m_verbosity > 5){
    pout() << "AmrMesh::defineRealms()" << endl;
  }

  for (auto& r : m_realms){
    r.second->define(m_grids, m_domains, m_refinementRatios, m_dx, m_probLo, m_finestLevel, m_numEbGhostsCells, m_numGhostCells, m_numLsfGhostCells, m_redistributionRadius,
		     m_hasEbCf, m_centroidStencilType, m_ebCentroidStencilType, m_baseif, m_multifluidIndexSpace);
  }
}

void AmrMesh::regridRealm(const std::string           a_realm,
			  const Vector<Vector<int> >& a_procs,
			  const Vector<Vector<Box> >& a_boxes,
			  const int                   a_lmin){
  CH_TIME("AmrMesh::regridRealm(procs, boxes, level)");
  if(m_verbosity > 1){
    pout() << "AmrMesh::regridRealm(procs, boxes, level)" << endl;
  }

  if(!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::define_realm - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  // Make the dbl
  Vector<DisjointBoxLayout> grids(1 + m_finestLevel);

  // Levels that didn't change. 
  for (int lvl = 0; lvl < a_lmin; lvl++){
    grids[lvl] = this->getGrids(a_realm)[lvl];
  }

  // Levels that did change. 
  for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++){
    grids[lvl] = DisjointBoxLayout();
    grids[lvl].define(a_boxes[lvl], a_procs[lvl], m_domains[lvl]);
    grids[lvl].close();
  }

  m_realms[a_realm]->define(grids, m_domains, m_refinementRatios, m_dx, m_probLo, m_finestLevel, m_numEbGhostsCells, m_numGhostCells, m_numLsfGhostCells, m_redistributionRadius,
			    m_hasEbCf, m_centroidStencilType, m_ebCentroidStencilType, m_baseif, m_multifluidIndexSpace);

  m_realms[a_realm]->regrid_base(a_lmin);
}

std::vector<std::string> AmrMesh::getRealms() const {
  std::vector<std::string> realms;

  for (const auto& r : m_realms){
    realms.push_back(r.first);
  }

  return realms;
}

BoxSorting AmrMesh::getBoxSorting() const{
  return m_boxSort;
}

#include <CD_NamespaceFooter.H>
