/*!
  @file amr_mesh.cpp
  @brief Implementation of amr_mesh.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "amr_mesh.H"
#include "mfalias.H"
#include "load_balance.H"
#include "gradientF_F.H"
#include "EBFastFineToCoarRedist.H"
#include "EBFastCoarToFineRedist.H"
#include "EBFastCoarToCoarRedist.H"
#include "DomainFluxIFFABFactory.H"
#include "TiledMeshRefine.H"

#include <BRMeshRefine.H>
#include <EBEllipticLoadBalance.H>
#include <EBLevelDataOps.H>
#include <MFLevelDataOps.H>
#include <EBArith.H>
#include <ParmParse.H>
#include <BaseIFFactory.H>

#include "CD_NamespaceHeader.H"

amr_mesh::amr_mesh(){
  parse_options();

  m_finest_level = 0;
  m_has_grids    = false;

#if 1
  // This is a fucking HACK! I have no idea why this work, and we should check that out. But, for some reason
  // if we use refinement of 4, without a 2 at the end (even if this is NOT used anywhere?!?!?), we get memory
  // leaks through define_eblevelgrid() in the setMaxCoarseningRatio stuff. I don't know if this leak is real,
  // but that's what valgrind tells me....
  //
  // Robert Marskar, Dec 7 (2018)
  //
  m_ref_ratios.resize(m_max_amr_depth);
  m_ref_ratios.push_back(2);
#endif

}

amr_mesh::~amr_mesh(){
  
}

EBAMRCellData amr_mesh::alias(const phase::which_phase a_phase, const MFAMRCellData& a_mfdata){
  CH_TIME("amr_mesh::alias(phase, mfamrcelldata)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::alias(phase, mfamrcelldata)" << endl;
  }

  EBAMRCellData ret;

  const int finestLevel = a_mfdata.size() - 1;

  allocate_ptr(ret, finestLevel);

  alias(ret, a_phase, a_mfdata, finestLevel);

  return ret;
}

EBAMRFluxData amr_mesh::alias(const phase::which_phase a_phase, const MFAMRFluxData& a_mfdata){
  CH_TIME("amr_mesh::alias(phase, mfamrfluxdata)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::alias(phase, mfamrfluxdata)" << endl;
  }

  EBAMRFluxData ret;

  const int finestLevel = a_mfdata.size() - 1;
  
  allocate_ptr(ret, finestLevel);

  alias(ret, a_phase, a_mfdata, finestLevel);

  return ret;
}

EBAMRIVData amr_mesh::alias(const phase::which_phase a_phase, const MFAMRIVData& a_mfdata){
  CH_TIME("amr_mesh::alias(phase, mfamrivdata)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::alias(phase, mfamrivdata)" << endl;
  }

  EBAMRIVData ret;
  allocate_ptr(ret);

  alias(ret, a_phase, a_mfdata);

  return ret;
}

void amr_mesh::alias(EBAMRCellData&           a_data,
		     const phase::which_phase a_phase,
		     const MFAMRCellData&     a_mfdata,
		     const int                a_finestLevel){
  CH_TIME("amr_mesh::alias(hardcap)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::alias(hardcap)" << endl;
  }

  for (int lvl = 0; lvl <= a_finestLevel; lvl++){
    mfalias::aliasMF(*a_data[lvl], a_phase, *a_mfdata[lvl]);
  }
}

void amr_mesh::alias(EBAMRFluxData&           a_data,
		     const phase::which_phase a_phase,
		     const MFAMRFluxData&     a_mfdata,
		     const int                a_finestLevel){
  CH_TIME("amr_mesh::alias(hardcap)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::alias(hardcap)" << endl;
  }

  for (int lvl = 0; lvl <= a_finestLevel; lvl++){
    mfalias::aliasMF(*a_data[lvl], a_phase, *a_mfdata[lvl]);
  }
}

void amr_mesh::alias(EBAMRCellData& a_data, const phase::which_phase a_phase, const MFAMRCellData& a_mfdata){
  CH_TIME("amr_mesh::alias");
  if(m_verbosity > 5){
    pout() << "amr_mesh::alias" << endl;
  }

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    mfalias::aliasMF(*a_data[lvl], a_phase, *a_mfdata[lvl]);
  }
}

void amr_mesh::alias(EBAMRFluxData& a_data, const phase::which_phase a_phase, const MFAMRFluxData& a_mfdata){
  CH_TIME("amr_mesh::alias");
  if(m_verbosity > 5){
    pout() << "amr_mesh::alias" << endl;
  }

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    mfalias::aliasMF(*a_data[lvl], a_phase, *a_mfdata[lvl]);
  }
}

void amr_mesh::alias(EBAMRIVData& a_data, const phase::which_phase a_phase, const MFAMRIVData& a_mfdata){
  CH_TIME("amr_mesh::alias");
  if(m_verbosity > 5){
    pout() << "amr_mesh::alias" << endl;
  }

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    mfalias::aliasMF(*a_data[lvl], a_phase, *a_mfdata[lvl]);
  }
}

void amr_mesh::allocate(EBAMRCellData&           a_data,
			const std::string        a_realm,
			const phase::which_phase a_phase,
			const int                a_ncomp,
			const int                a_ghost){
  CH_TIME("amr_mesh::allocate(ebamrcell, realm, phase, comp, ghost)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::allocate(ebamrcell, realm, phase, comp, ghost)" << endl;
  }

  if(!this->query_realm(a_realm)) {
    std::string str = "amr_mesh::allocate(ebamcell, realm, phase, comp, ghost) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const int ghost = (a_ghost == -1) ? m_num_ghost : a_ghost;

  a_data.resize(1 + m_finest_level);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_realms[a_realm]->get_grids()[lvl];
    const EBISLayout& ebisl      = m_realms[a_realm]->get_ebisl(a_phase)[lvl];
    
    EBCellFactory fact(ebisl);

    a_data[lvl] = RefCountedPtr<LevelData<EBCellFAB> >
      (new LevelData<EBCellFAB>(dbl, a_ncomp, ghost*IntVect::Unit, fact));

    EBLevelDataOps::setVal(*a_data[lvl], 0.0);
  }

  a_data.set_realm(a_realm);
}

void amr_mesh::allocate(EBAMRFluxData&           a_data,
			const std::string        a_realm,
			const phase::which_phase a_phase,
			const int                a_ncomp,
			const int                a_ghost){
  CH_TIME("amr_mesh::allocate(ebamrflux, realm, phase, comp, ghost)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::allocate(ebamrflux, realm, phase, comp, ghost)" << endl;
  }

  if(!this->query_realm(a_realm)) {
    std::string str = "amr_mesh::allocate(ebamrflux, realm, phase, comp, ghost) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const int ghost = (a_ghost == -1) ? m_num_ghost : a_ghost;

  a_data.resize(1 + m_finest_level);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_realms[a_realm]->get_grids()[lvl];
    const EBISLayout& ebisl      = m_realms[a_realm]->get_ebisl(a_phase)[lvl];

    EBFluxFactory fact(ebisl);

    a_data[lvl] = RefCountedPtr<LevelData<EBFluxFAB> >
      (new LevelData<EBFluxFAB>(dbl, a_ncomp, ghost*IntVect::Unit, fact));

    EBLevelDataOps::setVal(*a_data[lvl], 0.0);
  }

  a_data.set_realm(a_realm);
}

void amr_mesh::allocate(EBAMRIVData&             a_data,
			const std::string        a_realm,
			const phase::which_phase a_phase,
			const int                a_ncomp,
			const int                a_ghost){
  CH_TIME("amr_mesh::allocate(ebamriv, realm, phase, comp, ghost)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::allocate(ebamriv, realm, phase, comp, ghost)" << endl;
  }

  if(!this->query_realm(a_realm)) {
    std::string str = "amr_mesh::allocate(ebamriv, realm, phase, comp, ghost) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const int ghost = (a_ghost == -1) ? m_num_ghost : a_ghost;

  a_data.resize(1 + m_finest_level);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_realms[a_realm]->get_grids()[lvl];
    const EBISLayout& ebisl      = m_realms[a_realm]->get_ebisl(a_phase)[lvl];
    const ProblemDomain& domain  = m_realms[a_realm]->get_domains()[lvl];
    
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

void amr_mesh::allocate(EBAMRIFData&             a_data,
			const std::string        a_realm,
			const phase::which_phase a_phase,
			const int                a_ncomp,
			const int                a_ghost){
  CH_TIME("amr_mesh::allocate(ebamrifdata, realm, phase, comp, ghost)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::allocate(ebamrifdata, realm, phase, comp, ghost)" << endl;
  }

  if(!this->query_realm(a_realm)) {
    std::string str = "amr_mesh::allocate(ebamrifdata, realm, phase, comp, ghost) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const int ghost = (a_ghost == -1) ? m_num_ghost : a_ghost;

  a_data.resize(1 + m_finest_level);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_realms[a_realm]->get_grids()[lvl];
    const EBISLayout& ebisl      = m_realms[a_realm]->get_ebisl(a_phase)[lvl];
    const ProblemDomain& domain  = m_realms[a_realm]->get_domains()[lvl];
    
    DomainFluxIFFABFactory fact(ebisl, domain);

    a_data[lvl] = RefCountedPtr<LevelData<DomainFluxIFFAB> >
      (new LevelData<DomainFluxIFFAB>(dbl, a_ncomp, ghost*IntVect::Unit, fact));
  }

  a_data.set_realm(a_realm);
}

void amr_mesh::allocate(EBAMRBool& a_data, const std::string a_realm, const int a_ncomp, const int a_ghost){
  CH_TIME("amr_mesh::allocate(ebamrbool, realm, comp, ghost)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::allocate(ebamrbool, realm, comp, ghost)" << endl;
  }

  if(!this->query_realm(a_realm)) {
    std::string str = "amr_mesh::allocate(ebamrbool, realm, phase, comp, ghost) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  a_data.resize(1 + m_finest_level);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_realms[a_realm]->get_grids()[lvl];
    
    a_data[lvl] = RefCountedPtr<LevelData<BaseFab<bool> > >
      (new LevelData<BaseFab<bool> >(dbl, a_ncomp, a_ghost*IntVect::Unit));
  }

  a_data.set_realm(a_realm);
}

void amr_mesh::allocate(MFAMRCellData& a_data, const std::string a_realm, const int a_ncomp, const int a_ghost){
  CH_TIME("amr_mesh::allocate(mfamrcell, realm, comp, ghost)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::allocate(mfamrcell, realm, comp, ghost)" << endl;
  }

  if(!this->query_realm(a_realm)) {
    std::string str = "amr_mesh::allocate(mfamrcell, realm, comp, ghost) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const int ghost   = (a_ghost == -1) ? m_num_ghost : a_ghost;
  const int ignored = a_ncomp;
  const int nphases = m_mfis->num_phases();

  a_data.resize(1 + m_finest_level);

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_realms[a_realm]->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_realms[a_realm]->get_ebis(phase::solid);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_realms[a_realm]->get_grids()[lvl];
    
    Vector<EBISLayout> ebisl(nphases);
    Vector<int>        comps(nphases, a_ncomp);

    if(!ebis_gas.isNull()) ebisl[phase::gas]   = m_realms[a_realm]->get_ebisl(phase::gas)[lvl];
    if(!ebis_sol.isNull()) ebisl[phase::solid] = m_realms[a_realm]->get_ebisl(phase::solid)[lvl];
    
    MFCellFactory factory(ebisl, comps);

    a_data[lvl] = RefCountedPtr<LevelData<MFCellFAB> >
      (new LevelData<MFCellFAB>(dbl, ignored, ghost*IntVect::Unit, factory));

    MFLevelDataOps::setVal(*a_data[lvl], 0.0);
  }

  a_data.set_realm(a_realm);
}

void amr_mesh::allocate(MFAMRFluxData& a_data, const std::string a_realm, const int a_ncomp, const int a_ghost){
  CH_TIME("amr_mesh::allocate(mfamrflux, realm, comp, ghost)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::allocate(mfamrflux, realm, comp, ghost)" << endl;
  }

  if(!this->query_realm(a_realm)) {
    std::string str = "amr_mesh::allocate(mfamrflux, realm, comp, ghost) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const int ghost   = (a_ghost == -1) ? m_num_ghost : a_ghost;
  const int ignored = a_ncomp;
  const int nphases = m_mfis->num_phases();

  a_data.resize(1 + m_finest_level);

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_realms[a_realm]->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_realms[a_realm]->get_ebis(phase::solid);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_realms[a_realm]->get_grids()[lvl];
    
    Vector<EBISLayout> ebisl(nphases);
    Vector<int>        comps(nphases, a_ncomp);

    if(!ebis_gas.isNull()) ebisl[phase::gas]   = m_realms[a_realm]->get_ebisl(phase::gas)[lvl];
    if(!ebis_sol.isNull()) ebisl[phase::solid] = m_realms[a_realm]->get_ebisl(phase::solid)[lvl];
    
    MFFluxFactory factory(ebisl, comps);

    a_data[lvl] = RefCountedPtr<LevelData<MFFluxFAB> >
      (new LevelData<MFFluxFAB>(dbl, ignored, ghost*IntVect::Unit, factory));
  }

  a_data.set_realm(a_realm);
}

void amr_mesh::allocate(MFAMRIVData& a_data, const std::string a_realm, const int a_ncomp, const int a_ghost){
  CH_TIME("amr_mesh::allocate(mfamrivdata, realm, comp, ghost)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::allocate(mfamrivdata, realm, comp, ghost)" << endl;
  }
  
  if(!this->query_realm(a_realm)) {
    std::string str = "amr_mesh::allocate(mfamrivdata, realm, comp, ghost) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const int ghost   = (a_ghost == -1) ? m_num_ghost : a_ghost;
  const int ignored = a_ncomp;
  const int nphases = m_mfis->num_phases();

  a_data.resize(1 + m_finest_level);

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_realms[a_realm]->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_realms[a_realm]->get_ebis(phase::solid);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_realms[a_realm]->get_grids()[lvl];
    
    Vector<EBISLayout> ebisl(nphases);
    Vector<int>        comps(nphases, a_ncomp);

    if(!ebis_gas.isNull()) ebisl[phase::gas]   = m_realms[a_realm]->get_ebisl(phase::gas)[lvl];
    if(!ebis_sol.isNull()) ebisl[phase::solid] = m_realms[a_realm]->get_ebisl(phase::solid)[lvl];
    
    MFBaseIVFABFactory factory(ebisl, comps);

    a_data[lvl] = RefCountedPtr<LevelData<MFBaseIVFAB> >
      (new LevelData<MFBaseIVFAB>(dbl, ignored, ghost*IntVect::Unit, factory));
  }

  a_data.set_realm(a_realm);
}

void amr_mesh::reallocate(EBAMRCellData& a_data, const phase::which_phase a_phase, const int a_lmin){
  CH_TIME("amr_mesh::reallocate(ebamrcell, realm, phase, lmin)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::reallocate(ebamrcell, realm, phase, lmin)" << endl;
  }

  const std::string a_realm = a_data.get_realm();

  if(!this->query_realm(a_realm)) {
    std::string str = "amr_mesh::allocate(ebamrcell, realm, phase, lmin) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const IntVect ghost = a_data[0]->ghostVect();
  const int ncomp     = a_data[0]->nComp();

  a_data.resize(1 + m_finest_level);

  const int lmin = Max(0, a_lmin+1);
  
  for (int lvl = a_lmin; lvl <= m_finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_realms[a_realm]->get_grids()[lvl];
    const EBISLayout& ebisl      = m_realms[a_realm]->get_ebisl(a_phase)[lvl];
    
    EBCellFactory fact(ebisl);

    a_data[lvl] = RefCountedPtr<LevelData<EBCellFAB> >
      (new LevelData<EBCellFAB>(dbl, ncomp, ghost, fact));

    EBLevelDataOps::setVal(*a_data[lvl], 0.0);
  }
}

void amr_mesh::reallocate(EBAMRFluxData& a_data, const phase::which_phase a_phase, const int a_lmin){
  CH_TIME("amr_mesh::reallocate(ebamrflux, realm, phase, lmin)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::reallocate(ebamrflux, realm, phase, lmin)" << endl;
  }

  const std::string a_realm = a_data.get_realm();

  if(!this->query_realm(a_realm)) {
    std::string str = "amr_mesh::reallocate(ebamrflux, realm, phase, lmin) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const IntVect ghost = a_data[0]->ghostVect();
  const int ncomp     = a_data[0]->nComp();

  a_data.resize(1 + m_finest_level);

  const int lmin = Max(0, a_lmin+1);
  
  for (int lvl = a_lmin; lvl <= m_finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_realms[a_realm]->get_grids()[lvl];
    const EBISLayout& ebisl      = m_realms[a_realm]->get_ebisl(a_phase)[lvl];
    
    EBFluxFactory fact(ebisl);

    a_data[lvl] = RefCountedPtr<LevelData<EBFluxFAB> > (new LevelData<EBFluxFAB>(dbl, ncomp, ghost, fact));

    EBLevelDataOps::setVal(*a_data[lvl], 0.0);
  }
}

void amr_mesh::reallocate(EBAMRIVData& a_data, const phase::which_phase a_phase, const int a_lmin){
  CH_TIME("amr_mesh::reallocate(ebamriv, realm, phase, lmin)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::reallocate(ebamriv, realm, phase, lmin)" << endl;
  }

  const std::string a_realm = a_data.get_realm();

  if(!this->query_realm(a_realm)) {
    std::string str = "amr_mesh::reallocate(ebamriv, realm, phase, lmin) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const IntVect ghost = a_data[0]->ghostVect();
  const int ncomp     = a_data[0]->nComp();

  a_data.resize(1 + m_finest_level);

  const int lmin = Max(0, a_lmin+1);

  for (int lvl = a_lmin; lvl <= m_finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_realms[a_realm]->get_grids()[lvl];
    const EBISLayout& ebisl      = m_realms[a_realm]->get_ebisl(a_phase)[lvl];
    const ProblemDomain& domain  = m_realms[a_realm]->get_domains()[lvl];
    
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

void amr_mesh::reallocate(EBAMRIFData& a_data, const phase::which_phase a_phase, const int a_lmin){
  CH_TIME("amr_mesh::reallocate(ebamrifdata, realm, phase, lmin)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::reallocate(ebamrifdata, realm, phase, lmin)" << endl;
  }

  const std::string a_realm = a_data.get_realm(); 

  if(!this->query_realm(a_realm)) {
    std::string str = "amr_mesh::reallocate(ebamrif, realm, phase, lmin) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const IntVect ghost = a_data[0]->ghostVect();
  const int ncomp     = a_data[0]->nComp();

  a_data.resize(1 + m_finest_level);

  const int lmin = Max(0, a_lmin+1);

  for (int lvl = a_lmin; lvl <= m_finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_realms[a_realm]->get_grids()[lvl];
    const EBISLayout& ebisl      = m_realms[a_realm]->get_ebisl(a_phase)[lvl];
    const ProblemDomain& domain  = m_realms[a_realm]->get_domains()[lvl];
    
    DomainFluxIFFABFactory fact(ebisl, domain);

    a_data[lvl] = RefCountedPtr<LevelData<DomainFluxIFFAB> >
      (new LevelData<DomainFluxIFFAB>(dbl, ncomp, ghost, fact));
  }
}

void amr_mesh::reallocate(EBAMRBool& a_data, const int a_lmin){
  CH_TIME("amr_mesh::reallocate(ebamrbool, realm, lmin)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::reallocate(ebamrbool, realm, lmin)" << endl;
  }

  const std::string a_realm = a_data.get_realm();

  if(!this->query_realm(a_realm)) {
    std::string str = "amr_mesh::reallocate(ebamrbool, realm, lmin) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }
  
  const IntVect ghost = a_data[0]->ghostVect();
  const int ncomp     = a_data[0]->nComp();

  a_data.resize(1 + m_finest_level);

  const int lmin = Max(0, a_lmin+1);

  for (int lvl = a_lmin; lvl <= m_finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_realms[a_realm]->get_grids()[lvl];

    a_data[lvl] = RefCountedPtr<LevelData<BaseFab<bool> > >
      (new LevelData<BaseFab<bool> >(dbl, ncomp, ghost));
  }
}

void amr_mesh::reallocate(MFAMRCellData& a_data, const int a_lmin){
  CH_TIME("amr_mesh::reallocate(mfamrcell, realm, lmin)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::reallocate(mfamrcell, realm, lmin)" << endl;
  }

  const std::string a_realm = a_data.get_realm();

  if(!this->query_realm(a_realm)) {
    std::string str = "amr_mesh::reallocate(mfamrcell, realm, lmin) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const IntVect ghost = a_data[0]->ghostVect();
  const int ncomp = a_data[0]->nComp();
  const int ignored = ncomp;
  const int nphases = m_mfis->num_phases();

  a_data.resize(1 + m_finest_level);

  const int lmin = Max(0, a_lmin+1);

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_realms[a_realm]->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_realms[a_realm]->get_ebis(phase::solid);

  for (int lvl = a_lmin; lvl <= m_finest_level; lvl++){
    Vector<EBISLayout> ebisl(nphases);
    Vector<int>        comps(nphases, ncomp);

    const DisjointBoxLayout& dbl = m_realms[a_realm]->get_grids()[lvl];

    if(!ebis_gas.isNull()) ebisl[phase::gas]   = m_realms[a_realm]->get_ebisl(phase::gas)[lvl];
    if(!ebis_sol.isNull()) ebisl[phase::solid] = m_realms[a_realm]->get_ebisl(phase::solid)[lvl];
    
    MFCellFactory factory(ebisl, comps);

    a_data[lvl] = RefCountedPtr<LevelData<MFCellFAB> > (new LevelData<MFCellFAB>(dbl, ignored, ghost, factory));

    MFLevelDataOps::setVal(*a_data[lvl], 0.0);
  }
}

void amr_mesh::reallocate(MFAMRFluxData& a_data, const int a_lmin){
  CH_TIME("amr_mesh::allocate(mfamrflux, realm, lmin)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::allocate(mfamrflux, realm, lmin)" << endl;
  }

  const std::string a_realm = a_data.get_realm(); 

  if(!this->query_realm(a_realm)) {
    std::string str = "amr_mesh::reallocate(mfamrflux, realm, lmin) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const IntVect ghost = a_data[0]->ghostVect();
  const int ncomp     = a_data[0]->nComp();
  const int ignored   = ncomp;
  const int nphases   = m_mfis->num_phases();

  a_data.resize(1 + m_finest_level);

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_realms[a_realm]->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_realms[a_realm]->get_ebis(phase::solid);

  const int lmin = Max(0, a_lmin+1);

  for (int lvl = a_lmin; lvl <= m_finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_realms[a_realm]->get_grids()[lvl];
    
    Vector<EBISLayout> ebisl(nphases);
    Vector<int>        comps(nphases, ncomp);

    if(!ebis_gas.isNull()) ebisl[phase::gas]   = m_realms[a_realm]->get_ebisl(phase::gas)[lvl];
    if(!ebis_sol.isNull()) ebisl[phase::solid] = m_realms[a_realm]->get_ebisl(phase::solid)[lvl];
    
    MFFluxFactory factory(ebisl, comps);

    a_data[lvl] = RefCountedPtr<LevelData<MFFluxFAB> >
      (new LevelData<MFFluxFAB>(dbl, ignored, ghost, factory));
  }
}

void amr_mesh::reallocate(MFAMRIVData& a_data, const int a_lmin){
  CH_TIME("amr_mesh::allocate(mfamrivdata, realm, lmin)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::allocate(mfamrivdata, realm, lmin)" << endl;
  }

  const std::string a_realm = a_data.get_realm();

  if(!this->query_realm(a_realm)) {
    std::string str = "amr_mesh::reallocate(mfamrivdata, realm, lmin) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const IntVect ghost = a_data[0]->ghostVect();
  const int ncomp     = a_data[0]->nComp();
  const int ignored   = ncomp;
  const int nphases   = m_mfis->num_phases();

  a_data.resize(1 + m_finest_level);

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_realms[a_realm]->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_realms[a_realm]->get_ebis(phase::solid);

  const int lmin = Max(0, a_lmin+1);

  for (int lvl = a_lmin; lvl <= m_finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_realms[a_realm]->get_grids()[lvl];
    
    Vector<EBISLayout> ebisl(nphases);
    Vector<int>        comps(nphases, ncomp);

    if(!ebis_sol.isNull()) ebisl[phase::gas]   = m_realms[a_realm]->get_ebisl(phase::gas)[lvl];
    if(!ebis_sol.isNull()) ebisl[phase::solid] = m_realms[a_realm]->get_ebisl(phase::solid)[lvl];
    
    MFBaseIVFABFactory factory(ebisl, comps);

    a_data[lvl] = RefCountedPtr<LevelData<MFBaseIVFAB> >
      (new LevelData<MFBaseIVFAB>(dbl, ignored, ghost, factory));
  }
}

void amr_mesh::set_mfis(const RefCountedPtr<mfis>& a_mfis){
  CH_TIME("amr_mesh::set_mfis");
  if(m_verbosity > 5){
    pout() << "amr_mesh::set_mfis" << endl;
  }

  m_mfis = a_mfis;
}

void amr_mesh::set_baseif(const phase::which_phase a_phase, const RefCountedPtr<BaseIF>& a_baseif){
  CH_TIME("amr_mesh::set_baseif");
  if(m_verbosity > 5){
    pout() << "amr_mesh::set_baseif" << endl;
  }

  m_baseif.emplace(a_phase, a_baseif);
}

void amr_mesh::parse_options(){
  parse_verbosity();
  parse_coarsest_num_cells();
  parse_max_amr_depth();
  parse_max_simulation_depth();
  parse_ebcf();
  parse_refinement_ratio();
  parse_blocking_factor();
  parse_max_box_size();
  parse_max_ebis_box_size();
  parse_grid_generation();
  parse_buffer_size();
  parse_irreg_growth();
  parse_fill_ratio();
  parse_redist_rad();;
  parse_num_ghost();
  parse_eb_ghost();
  parse_domain();
  parse_ghost_interpolation();
  parse_centroid_stencils();
  parse_eb_stencils();
}

void amr_mesh::parse_runtime_options(){
  parse_verbosity();
  parse_blocking_factor();
  parse_max_box_size();
  parse_grid_generation();
  parse_buffer_size();
  parse_irreg_growth();
  parse_fill_ratio();
}

void amr_mesh::parse_domain(){

  ParmParse pp("amr_mesh");

  Vector<Real> v(SpaceDim);
  pp.getarr("lo_corner", v, 0, SpaceDim); m_prob_lo = RealVect(D_DECL(v[0], v[1], v[2]));
  pp.getarr("hi_corner", v, 0, SpaceDim); m_prob_hi = RealVect(D_DECL(v[0], v[1], v[2]));
}

void amr_mesh::parse_ghost_interpolation(){

  std::string interp_type;
  ParmParse pp("amr_mesh");
  pp.get("ghost_interp", interp_type);
  if(interp_type == "pwl"){
    m_interp_type = ghost_interpolation::pwl;
  }
  else if(interp_type == "quad"){
    m_interp_type = ghost_interpolation::quad;
  }
  else{
    MayDay::Abort("amr_mesh::parse_ghost_interpolation - unknown ghost interpolation requested");
  }
}

void amr_mesh::build_domains(){
  CH_TIME("amr_mesh::build_domains");
  if(m_verbosity > 5){
    pout() << "amr_mesh::build_domains" << endl;
  }

  const int nlevels = 1 + m_max_amr_depth;
  
  m_domains.resize(nlevels);
  m_dx.resize(nlevels);
  m_grids.resize(nlevels);


  m_dx[0] = (m_prob_hi[0] - m_prob_lo[0])/m_num_cells[0];
  m_domains[0] = ProblemDomain(IntVect::Zero, m_num_cells - IntVect::Unit);

  for (int lvl = 1; lvl <= m_max_amr_depth; lvl++){
    m_dx[lvl]      = m_dx[lvl-1]/m_ref_ratios[lvl-1];
    m_domains[lvl] = m_domains[lvl-1];
    m_domains[lvl].refine(m_ref_ratios[lvl-1]);
  }
}

void amr_mesh::regrid(const Vector<IntVectSet>& a_tags,
		      const int a_lmin,
		      const int a_lmax,
		      const int a_regsize,
		      const int a_hardcap){
  CH_TIME("amr_mesh::regrid");
  if(m_verbosity > 1){
    pout() << "amr_mesh::regrid" << endl;
  }

  this->regrid_amr(a_tags, a_lmin, a_lmax, a_hardcap);
  this->regrid_operators(a_lmin, a_lmax, a_regsize);
}

void amr_mesh::regrid_amr(const Vector<IntVectSet>& a_tags,
			  const int a_lmin,
			  const int a_lmax,
			  const int a_hardcap){
  CH_TIME("amr_mesh::regrid_amr(tags, level, level, hardcap)");
  if(m_verbosity > 1){
    pout() << "amr_mesh::regrid_amr(tags, level, level, hardcap)" << endl;
  }

  // TLDR: This is the version that reads boxes. amr_mesh makes the grids from the tags and load balances them
  //       by using the patch volume. Those grids are then sent to the various realms. 
  
  Vector<IntVectSet> tags = a_tags; // build_grids destroys tags, so copy them

  // This is stuff that always gets done
  this->build_grids(tags, a_lmin, a_lmax, a_hardcap);

  // Define realms with the new grids and redo the realm stuff
  this->define_realms();

  for (auto& r : m_realms){
    r.second->regrid_base(a_lmin);
  }
}


void amr_mesh::regrid_amr(const Vector<Vector<int> >& a_procs, const Vector<Vector<Box> >& a_boxes, const int a_lmin){
  CH_TIME("amr_mesh::regrid_amr(procs, boxes, level)");
  if(m_verbosity > 1){
    pout() << "amr_mesh::regrid_amr(procs, boxes, level)" << endl;
  }

  // TLDR: This is the version that reads boxes. amr_mesh makes the grids by using the patch volume as load,
  //       and those grids are then sent to the various realms. 
  
  for (int lvl = a_lmin; lvl <= m_finest_level; lvl++){
    m_grids[lvl] = DisjointBoxLayout();
    m_grids[lvl].define(a_boxes[lvl], a_procs[lvl], m_domains[lvl]);
    m_grids[lvl].close(); 
  }

  this->define_realms();

  // Regrid the base on every realm. This includes EBLevelGrid,  neighbors, and vof iterators. 
  for (auto& r : m_realms){
    r.second->regrid_base(a_lmin);
  }
}

void amr_mesh::regrid_operators(const int a_lmin,
				const int a_lmax,
				const int a_regsize){
  CH_TIME("amr_mesh::regrid_operators(procs, boxes, level)");
  if(m_verbosity > 1){
    pout() << "amr_mesh::regrid_operators(procs, boxes, level)" << endl;
  }
  
  for (auto& r : m_realms){
    r.second->regrid_operators(a_lmin, a_lmax, a_regsize);
  }
}

void amr_mesh::build_grids(Vector<IntVectSet>& a_tags, const int a_lmin, const int a_lmax, const int a_hardcap){
  CH_TIME("amr_mesh::build_grids");
  if(m_verbosity > 2){
    pout() << "amr_mesh::build_grids" << endl;
  }

  // TLDR: a_lmin is the coarsest level that changes. A special condition is a_lmin=0 for which we assume
  //       that there are no prior grids.
  //       a_lmax is the finest level that changes. This is typically 

  // base is the coarsest level which does not change. top_level is the finest level where we have tags. We should never
  // have tags on max_amr_depth, and we make that restriction here.
  const int base      = Max(0, a_lmin - 1);
  const int top_level = (m_finest_level == m_max_amr_depth) ? m_finest_level - 1 : a_tags.size() - 1;

  // New and old grid boxes
  Vector<Vector<Box> > new_boxes(1 + top_level);
  Vector<Vector<Box> > old_boxes(1 + top_level);

  // Enforce potential hardcap. 
  const int hardcap = (a_hardcap == -1) ? m_max_amr_depth : a_hardcap;

  // Inside this loop we make the boxes. 
  if(m_max_amr_depth > 0 && hardcap > 0){
    domainSplit(m_domains[0], old_boxes[0], m_max_box_size, m_blocking_factor);

    if(!m_has_grids){
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
    if(m_gridgen == grid_generation::berger_rigoustous){
      BRMeshRefine mesh_refine(m_domains[0], m_ref_ratios, m_fill_ratio, m_blocking_factor, m_buffer_size, m_max_box_size);
      new_finest_level = mesh_refine.regrid(new_boxes, a_tags, base, top_level, old_boxes);
    }
    else if (m_gridgen == grid_generation::tiled){
      TiledMeshRefine mesh_refine(m_domains[0], m_ref_ratios, m_blocking_factor*IntVect::Unit);
      new_finest_level = mesh_refine.regrid(new_boxes, a_tags, base, top_level, old_boxes);
    }
    else{
      MayDay::Abort("amr_mesh::regrid - logic bust, regridding with unknown regrid algorithm");
    }
    
    m_finest_level = Min(new_finest_level, m_max_amr_depth); // Don't exceed m_max_amr_depth
    m_finest_level = Min(m_finest_level,   m_max_sim_depth); // Don't exceed maximum simulation depth
    m_finest_level = Min(m_finest_level,   hardcap);         // Don't exceed hardcap
  }
  else{ // Only end up here if we have a single grid level, i.e. just single-level grid decomposition. 
    new_boxes.resize(1);
    domainSplit(m_domains[0], new_boxes[0], m_max_box_size, m_blocking_factor);
    
    m_finest_level = 0;
  }

  if(a_lmin == 0){ // Coarsest level also changes in this case, but that's not caught by the regridders. 
    domainSplit(m_domains[0], new_boxes[0], m_max_box_size, m_blocking_factor);
  }

  // Morton order the boxes.
  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    load_balance::sort(new_boxes[lvl], m_boxsort);
  }

  // Load balance boxes with patch volume as load proxy. 
  Vector<Vector<int> > pid(1 + m_finest_level);
  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    load_balance::make_balance(pid[lvl], new_boxes[lvl]);
  }

  // Define grids. If a_lmin=0 every grid is new, otherwise keep old grids up to but not including a_lmin
  if(a_lmin == 0){
    m_grids.resize(1 + m_finest_level);
    for (int lvl = 0; lvl <= m_finest_level; lvl++){
      m_grids[lvl] = DisjointBoxLayout();
      m_grids[lvl].define(new_boxes[lvl], pid[lvl], m_domains[lvl]);
      m_grids[lvl].close();
    }
  }
  else{
    Vector<DisjointBoxLayout> old_grids = m_grids;
    m_grids.resize(1 + m_finest_level);
    for (int lvl = 0; lvl < a_lmin; lvl++){ // Copy old grids
      m_grids[lvl] = old_grids[lvl];
    }
    for (int lvl = a_lmin; lvl <= m_finest_level; lvl++){ // Create new ones from tags
      m_grids[lvl] = DisjointBoxLayout();
      m_grids[lvl].define(new_boxes[lvl], pid[lvl], m_domains[lvl]);
      m_grids[lvl].close();
    }
  }

  m_has_grids = true;
}

void amr_mesh::compute_gradient(LevelData<EBCellFAB>&       a_gradient,
				const LevelData<EBCellFAB>& a_phi,
				const std::string           a_realm,
				const phase::which_phase    a_phase,
				const int                   a_lvl){
  CH_TIME("amr_mesh::compute_gradient(grad, phi, realm, phase, lvl)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::compute_gradient(grad, phi, realm, phase,lvl)" << endl;
  }

  CH_assert(a_phi.nComp()      == 1);
  CH_assert(a_gradient.nComp() == SpaceDim);

  if(!this->query_realm(a_realm)) {
    std::string str = "amr_mesh::compute_gradient(grad, phi, realm, phase, lvl) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }
    
  const int comp  = 0;
  const int ncomp = 1;
    
  const Real& dx = m_realms[a_realm]->get_dx()[a_lvl];
  const DisjointBoxLayout& dbl = m_realms[a_realm]->get_grids()[a_lvl];
  const ProblemDomain& domain  = m_realms[a_realm]->get_domains()[a_lvl];

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

void amr_mesh::compute_gradient(EBAMRCellData&           a_gradient,
				const EBAMRCellData&     a_phi,
				const std::string        a_realm,  
				const phase::which_phase a_phase){
  CH_TIME("amr_mesh::compute_gradient(grad, phi, realm, phase)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::compute_gradient(grad, phi, realm, phase)" << endl;
  }

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    this->compute_gradient(*a_gradient[lvl], *a_phi[lvl], a_realm, a_phase, lvl);
  }
}

void amr_mesh::compute_gradient(MFAMRCellData& a_gradient, const MFAMRCellData& a_phi, const std::string a_realm){
  CH_TIME("amr_mesh::compute_gradient(mf grad, mf phi, realm)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::compute_gradient(mf grad, mf phi, realm)" << endl;
  }

  for (int iphase = 0; iphase < m_mfis->num_phases(); iphase++){
    EBAMRCellData alias_grad(1 + m_finest_level);
    EBAMRCellData alias_phi(1 + m_finest_level);

    for (int lvl = 0; lvl <= m_finest_level; lvl++){
      alias_grad[lvl] = RefCountedPtr<LevelData<EBCellFAB> > (new LevelData<EBCellFAB>());
      alias_phi[lvl]  = RefCountedPtr<LevelData<EBCellFAB> > (new LevelData<EBCellFAB>());
      
      mfalias::aliasMF(*alias_grad[lvl], iphase, *a_gradient[lvl]);
      mfalias::aliasMF(*alias_phi[lvl],  iphase, *a_phi[lvl]);
    }

    if(iphase == 0){
      this->compute_gradient(alias_grad, alias_phi, a_realm, phase::gas);
    }
    else if(iphase == 1){
      this->compute_gradient(alias_grad, alias_phi, a_realm, phase::solid);
    }
  }
}

void amr_mesh::average_down(EBAMRCellData& a_data, const std::string a_realm, const phase::which_phase a_phase){
  CH_TIME("amr_mesh::average_down(ebamrcell, realm, phase");
  if(m_verbosity > 3){
    pout() << "amr_mesh::average_down(ebamrcell, realm, phase)" << endl;
  }

  if(!this->query_realm(a_realm)) {
    std::string str = "amr_mesh::average_down(ebamrcell, realm, phase) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }
  
  for (int lvl = m_finest_level; lvl > 0; lvl--){
    ebcoarseaverage& aveOp = *m_realms[a_realm]->get_coarave(a_phase)[lvl];
      
    const int ncomps = a_data[lvl]->nComp();
    const Interval interv (0, ncomps-1);

    aveOp.average(*a_data[lvl-1], *a_data[lvl], interv);
  }
    
  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    a_data[lvl]->exchange();
  }
}

void amr_mesh::average_down(EBAMRCellData& a_data, const std::string a_realm, const phase::which_phase a_phase, const int a_lvl){
  CH_TIME("amr_mesh::average_down(ebamrcelldata, realm, phase, level");
  if(m_verbosity > 3){
    pout() << "amr_mesh::average_down(ebamrcelldata, realm, phase, level)" << endl;
  }

  if(!this->query_realm(a_realm)) {
    std::string str = "amr_mesh::average_down(ebamrcell, realm, phase, level) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }
  
  const int ncomps = a_data[a_lvl]->nComp();
  const Interval interv (0, ncomps-1);

  ebcoarseaverage& aveOp = *m_realms[a_realm]->get_coarave(a_phase)[a_lvl+1];

  aveOp.average(*a_data[a_lvl], *a_data[a_lvl+1], interv);

  a_data[a_lvl]->exchange();
}

void amr_mesh::average_down(MFAMRFluxData& a_data, const std::string a_realm){
  CH_TIME("amr_mesh::average_down(mfamrflux, realm)");
  if(m_verbosity > 3){
    pout() << "amr_mesh::average_down(mfamrflux, realm)" << endl;
  }

  if(!this->query_realm(a_realm)) {
    std::string str = "amr_mesh::average_down(mfamrflux, realm) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const RefCountedPtr<EBIndexSpace>& ebis_gas = m_realms[a_realm]->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace>& ebis_sol = m_realms[a_realm]->get_ebis(phase::solid);

  // Alias the data to regular EBFluxFABs
  EBAMRFluxData alias_g(1 + m_finest_level);
  EBAMRFluxData alias_s(1 + m_finest_level);
    
  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    alias_g[lvl] = RefCountedPtr<LevelData<EBFluxFAB> > (new LevelData<EBFluxFAB>());
    alias_s[lvl] = RefCountedPtr<LevelData<EBFluxFAB> > (new LevelData<EBFluxFAB>());
      
    if(!ebis_gas.isNull()) mfalias::aliasMF(*alias_g[lvl], phase::gas,   *a_data[lvl]);
    if(!ebis_sol.isNull()) mfalias::aliasMF(*alias_s[lvl], phase::solid, *a_data[lvl]);
  }

  if(!ebis_gas.isNull()) this->average_down(alias_g, a_realm, phase::gas);
  if(!ebis_sol.isNull()) this->average_down(alias_s, a_realm, phase::solid);
}

void amr_mesh::average_down(MFAMRCellData& a_data, const std::string a_realm){
  CH_TIME("amr_mesh::average_down(mfamrcell, realm)");
  if(m_verbosity > 3){
    pout() << "amr_mesh::average_down(mfamrcell, realm)" << endl;
  }

  if(!this->query_realm(a_realm)) {
    std::string str = "amr_mesh::average_down(mfamrcell, realm) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const RefCountedPtr<EBIndexSpace>& ebis_gas = m_realms[a_realm]->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace>& ebis_sol = m_realms[a_realm]->get_ebis(phase::solid);
  
  EBAMRCellData alias_g(1 + m_finest_level);
  EBAMRCellData alias_s(1 + m_finest_level);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    alias_g[lvl] = RefCountedPtr<LevelData<EBCellFAB> > (new LevelData<EBCellFAB>());
    alias_s[lvl] = RefCountedPtr<LevelData<EBCellFAB> > (new LevelData<EBCellFAB>());
      
    if(!ebis_gas.isNull()) mfalias::aliasMF(*alias_g[lvl], phase::gas,   *a_data[lvl]);
    if(!ebis_sol.isNull()) mfalias::aliasMF(*alias_s[lvl], phase::solid, *a_data[lvl]);
  }

  if(!ebis_gas.isNull()) this->average_down(alias_g, a_realm, phase::gas);
  if(!ebis_sol.isNull()) this->average_down(alias_s, a_realm, phase::solid);
}

void amr_mesh::average_down(EBAMRFluxData& a_data, const std::string a_realm, const phase::which_phase a_phase){
  CH_TIME("amr_mesh::average_down(ebamrflux, realm, phase");
  if(m_verbosity > 3){
    pout() << "amr_mesh::average_down(ebamrflux, realm, phase)" << endl;
  }

  if(!this->query_realm(a_realm)) {
    std::string str = "amr_mesh::average_down(ebamrflux, realm, phase) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  for (int lvl = m_finest_level; lvl > 0; lvl--){
    const int ncomps = a_data[lvl]->nComp();
    const Interval interv (0, ncomps-1);

    ebcoarseaverage& aveOp = *m_realms[a_realm]->get_coarave(a_phase)[lvl];
    aveOp.average(*a_data[lvl-1], *a_data[lvl], interv);
  }

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    a_data[lvl]->exchange();
  }
}

void amr_mesh::average_down(EBAMRIVData& a_data, const std::string a_realm, const phase::which_phase a_phase){
  CH_TIME("amr_mesh::average_down(ebamriv, realm, phase)");
  if(m_verbosity > 3){
    pout() << "amr_mesh::average_down(ebamriv, realm, phase)" << endl;
  }

  if(!this->query_realm(a_realm)) {
    std::string str = "amr_mesh::average_down(ebamriv, realm, phase) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  for (int lvl = m_finest_level; lvl > 0; lvl--){
    const int ncomps = a_data[lvl]->nComp();
    const Interval interv (0, ncomps-1);

    ebcoarseaverage& aveOp = *m_realms[a_realm]->get_coarave(a_phase)[lvl];
    aveOp.average(*a_data[lvl-1], *a_data[lvl], interv);
  }

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    a_data[lvl]->exchange();
  }
}

void amr_mesh::conservative_average(EBAMRIVData& a_data, const std::string a_realm, const phase::which_phase a_phase){
  CH_TIME("amr_mesh::conservative_average(ebamriv, realm, phase)");
  if(m_verbosity > 3){
    pout() << "amr_mesh::conservative_average(ebamriv, realm, phase)" << endl;
  }

  if(!this->query_realm(a_realm)) {
    std::string str = "amr_mesh::conservative_average(ebamriv, realm, phase) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  for (int lvl = m_finest_level; lvl > 0; lvl--){
    const int ncomps = a_data[lvl]->nComp();
    const Interval interv (0, ncomps-1);

    ebcoarseaverage& aveOp = *m_realms[a_realm]->get_coarave(a_phase)[lvl];
    
    aveOp.conservative_average(*a_data[lvl-1], *a_data[lvl], interv);
  }

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    a_data[lvl]->exchange();
  }
}

void amr_mesh::interp_ghost(EBAMRCellData& a_data, const std::string a_realm, const phase::which_phase a_phase){
  CH_TIME("amr_mesh::interp_ghost(ebamrcell, realm, phase)");
  if(m_verbosity > 3){
    pout() << "amr_mesh::interp_ghost(ebamrcell, realm, phase)" << endl;
  }

  if(!this->query_realm(a_realm)) {
    std::string str = "amr_mesh::interp_ghost(ebamrcell, realm, phase) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }  

  if(m_interp_type == ghost_interpolation::pwl){
    this->interp_ghost_pwl(a_data, a_realm, a_phase);
  }
  else if(m_interp_type == ghost_interpolation::quad){
    this->interp_ghost_quad(a_data, a_realm, a_phase);
  }
  else{
    MayDay::Abort("amr_mesh::interp_ghost - unsupported interpolation type requested");
  }
}

void amr_mesh::interp_ghost(LevelData<EBCellFAB>&       a_fineData,
			    const LevelData<EBCellFAB>& a_coarData,
			    const int                   a_fineLevel,
			    const std::string           a_realm,
			    const phase::which_phase    a_phase){
  CH_TIME("amr_mesh::interp_ghost(fine, coar, level, realm, phase)");
  if(m_verbosity > 3){
    pout() << "amr_mesh::interp_ghost(fine, coar, level, realm, phase)" << endl;
  }

  if(!this->query_realm(a_realm)) {
    std::string str = "amr_mesh::interp_ghost(fine, coar, level, realm, phase) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  if(a_fineLevel > 0){

    const int ncomps      = a_fineData.nComp();
    const Interval interv = Interval(0, ncomps-1);
    
    if(m_interp_type == ghost_interpolation::pwl){
      AggEBPWLFillPatch& fillpatch = *m_realms[a_realm]->get_fillpatch(a_phase)[a_fineLevel];
    
      fillpatch.interpolate(a_fineData, a_coarData, a_coarData, 0.0, 0.0, 0.0, interv);
    }
    else if(m_interp_type == ghost_interpolation::quad){
      nwoebquadcfinterp& quadcfi = *m_realms[a_realm]->get_quadcfi(a_phase)[a_fineLevel];
      quadcfi.coarseFineInterp(a_fineData, a_coarData, 0, 0, ncomps);
    }
    else{
      MayDay::Abort("amr_mesh::interp_ghost - unsupported interpolation type requested");
    }
  }
}

void amr_mesh::interp_ghost(MFAMRCellData& a_data, const std::string a_realm){
  CH_TIME("amr_mesh::interp_ghost(mfamrcell, realm)");
  if(m_verbosity > 3){
    pout() << "amr_mesh::interp_ghost(mfamrcell, realm)" << endl;
  }

  if(!this->query_realm(a_realm)) {
    std::string str = "amr_mesh::interp_ghost(mfamrcell, realm, phase) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  // Do aliasing
  EBAMRCellData alias_g(1 + m_finest_level);
  EBAMRCellData alias_s(1 + m_finest_level);

  const RefCountedPtr<EBIndexSpace>& ebis_gas = m_realms[a_realm]->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace>& ebis_sol = m_realms[a_realm]->get_ebis(phase::solid);
  
  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    alias_g[lvl] = RefCountedPtr<LevelData<EBCellFAB> > (new LevelData<EBCellFAB>());
    alias_s[lvl] = RefCountedPtr<LevelData<EBCellFAB> > (new LevelData<EBCellFAB>());
      
    if(!ebis_gas.isNull()) mfalias::aliasMF(*alias_g[lvl], phase::gas,   *a_data[lvl]);
    if(!ebis_sol.isNull()) mfalias::aliasMF(*alias_s[lvl], phase::solid, *a_data[lvl]);
  }

  if(!ebis_gas.isNull()) this->interp_ghost(alias_g, a_realm, phase::gas);
  if(!ebis_sol.isNull()) this->interp_ghost(alias_s, a_realm, phase::solid);
}

void amr_mesh::interp_ghost_quad(EBAMRCellData& a_data, const std::string a_realm, const phase::which_phase a_phase){
  CH_TIME("amr_mesh::interp_ghost_quad(ebamrcell, realm, phase)");
  if(m_verbosity > 3){
    pout() << "amr_mesh::interp_ghost_quad(ebamrcell, realm, phase)" << endl;
  }

  if(!this->query_realm(a_realm)) {
    std::string str = "amr_mesh::interp_ghost_quad(ebamrcell, realm, phase) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  for (int lvl = m_finest_level; lvl > 0; lvl--){
    const int ncomps = a_data[lvl]->nComp();
    const Interval interv(0, ncomps -1);

    nwoebquadcfinterp& quadcfi = *m_realms[a_realm]->get_quadcfi(a_phase)[lvl];

    quadcfi.coarseFineInterp(*a_data[lvl], *a_data[lvl-1], 0, 0, ncomps);
  }

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    a_data[lvl]->exchange();
  }
}

void amr_mesh::interp_ghost_pwl(EBAMRCellData& a_data, const std::string a_realm, const phase::which_phase a_phase){
  CH_TIME("amr_mesh::interp_ghost_pwl(ebamrcell, realm, phase)");
  if(m_verbosity > 3){
    pout() << "amr_mesh::interp_ghost_pwl(ebamrcell, realm, phase)" << endl;
  }

  if(!this->query_realm(a_realm)) {
    std::string str = "amr_mesh::interp_ghost_pwl(ebamrcell, realm, phase) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }
  
  for (int lvl = m_finest_level; lvl > 0; lvl--){
    const int ncomps = a_data[lvl]->nComp();
    const Interval interv(0, ncomps -1);

    AggEBPWLFillPatch& fillpatch = *m_realms[a_realm]->get_fillpatch(a_phase)[lvl];
    
    fillpatch.interpolate(*a_data[lvl], *a_data[lvl-1], *a_data[lvl-1], 0.0, 0.0, 0.0, interv);
  }

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    a_data[lvl]->exchange();
  }
}

void amr_mesh::interpolate_to_centroids(EBAMRCellData& a_data, const std::string a_realm, const phase::which_phase a_phase){
  CH_TIME("amr_mesh::interpolate_to_centroids(ebamrcell, realm, phase)");
  if(m_verbosity > 3){
    pout() << "amr_mesh::interpolate_to_centroids(ebamrcell, realm, phase)" << endl;
  }

  if(!this->query_realm(a_realm)) {
    std::string str = "amr_mesh::interpolate_to_centroids(ebamrcell, realm, phase) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }
  
  irreg_amr_stencil<centroid_interp>& stencil = m_realms[a_realm]->get_centroid_interp_stencils(a_phase);
  stencil.apply(a_data);
}

void amr_mesh::parse_verbosity(){
  CH_TIME("amr_mesh::parse_verbosity");

  ParmParse pp("amr_mesh");
  pp.get("verbosity", m_verbosity);
}

void amr_mesh::parse_coarsest_num_cells(){

  ParmParse pp("amr_mesh");
  Vector<int> cells;
  cells.resize(pp.countval("coarsest_domain"));
  CH_assert(cells.size() >= SpaceDim);
  pp.getarr("coarsest_domain", cells, 0, SpaceDim);

  m_num_cells = IntVect(D_DECL(cells[0], cells[1], cells[2]));
}

void amr_mesh::parse_max_amr_depth(){

  ParmParse pp("amr_mesh");
  int depth;
  pp.get("max_amr_depth", depth);
  if(depth >= 0){
    m_max_amr_depth = depth;
  }
  else{
    m_max_amr_depth = 0;
  }
}

void amr_mesh::parse_max_simulation_depth(){

  ParmParse pp("amr_mesh");
  int depth;
  pp.get("max_sim_depth", depth);
  if(depth >= 0){
    m_max_sim_depth = depth;
  }
  else {
    m_max_sim_depth = m_max_amr_depth;
  }
}

void amr_mesh::parse_ebcf(){
  ParmParse pp("amr_mesh");
  std::string str;
  pp.get("ebcf", str);
  if(str == "true"){
    m_ebcf = true;
  }
  else if(str == "false"){
    m_ebcf = false;
  }
}

void amr_mesh::parse_refinement_ratio(){

  ParmParse pp("amr_mesh");
  Vector<int> ratios;
  ratios.resize(pp.countval("ref_rat"));
  pp.getarr("ref_rat", ratios, 0, ratios.size());

  // Pad with 2 if user didn't supply enough
  while(ratios.size() < m_max_amr_depth){
    ratios.push_back(2);
  }
  
  m_ref_ratios = ratios;
}

void amr_mesh::set_refinement_ratios(const Vector<int> a_ref_ratios){
  m_ref_ratios = a_ref_ratios;
}

void amr_mesh::parse_fill_ratio(){
  ParmParse pp("amr_mesh");
  Real fill_ratio = 1.0;
  pp.get("fill_ratio", fill_ratio);
  if(fill_ratio > 0.0 && fill_ratio <= 1.0){
    m_fill_ratio = fill_ratio;
  }
}

void amr_mesh::set_finest_level(const int a_finest_level){
  m_finest_level = a_finest_level;
  m_finest_level = Min(m_finest_level, m_max_amr_depth); // Don't exceed m_max_amr_depth
  m_finest_level = Min(m_finest_level, m_max_sim_depth); // Don't exceed maximum simulation depth
}

void amr_mesh::set_grids(const Vector<Vector<Box> >& a_boxes, const std::map<std::string, Vector<Vector<long int> > >& a_realms_and_loads){
  CH_TIME("amr_mesh::set_grids(boxes, loads, regsize)");
  if(m_verbosity > 3){
    pout() << "amr_mesh::set_grids(boxes, loads, regsize)" << endl;
  }

  const int lmin = 0;

  for (const auto& r : a_realms_and_loads){
    const std::string&               cur_realm = r.first;
    const Vector<Vector<long int> >& cur_loads = r.second;

    Vector<Vector<int> > pids(1 + m_finest_level);

    // Do load balancing. 
    for (int lvl = 0; lvl <= m_finest_level; lvl++){
      load_balance::make_balance(pids[lvl], cur_loads[lvl], a_boxes[lvl]);
    }

    this->regrid_realm(cur_realm, pids, a_boxes, lmin);
  }

  // Set the proxy grids, too. These are load balanced using the patch volume. 
  m_grids = m_realms[realm::primal]->get_grids();
  m_has_grids = true;
}

void amr_mesh::parse_max_box_size(){
  ParmParse pp("amr_mesh");
  int box_size;
  pp.get("max_box_size", box_size);
  if(box_size >= 8 && box_size % 2 == 0){
    m_max_box_size = box_size;
  }
  else{
    MayDay::Abort("amr_mesh::parse_max_box_size - must have box_size >= 8 and divisible by 2");
  }
}

void amr_mesh::parse_max_ebis_box_size(){

  ParmParse pp("amr_mesh");
  int box_size;
  pp.get("max_ebis_box", box_size);
  if(box_size >= 8 && box_size % 2 == 0){
    m_max_ebis_box_size = box_size;
  }
  else{
    MayDay::Abort("amr_mesh::parse_max_ebis_box_size - must have box_size >= 8 and divisible by 2");
  }
}

void amr_mesh::parse_buffer_size(){

  ParmParse pp("amr_mesh");
  int buffer = 2;
  pp.get("buffer_size", buffer);
  if(buffer > 0){
    m_buffer_size = buffer;
  }
}

void amr_mesh::parse_grid_generation(){

  ParmParse pp("amr_mesh");
  std::string str;
  pp.get("grid_algorithm", str);
  if(str == "br"){
    m_gridgen = grid_generation::berger_rigoustous;
  }
  else if(str == "tiled"){
    m_gridgen = grid_generation::tiled;
  }
  else{
    MayDay::Abort("amr_mesh::parse_grid_generation - unknown grid generation method requested");
  }

  pp.get("box_sorting", str);
  if( str == "none"){
    m_boxsort = box_sorting::none;
  }
  if( str == "std"){
    m_boxsort = box_sorting::std;
  }
  else if(str == "shuffle"){
    m_boxsort = box_sorting::shuffle;
  }
  else if(str == "morton"){
    m_boxsort = box_sorting::morton;
  }
  else {
    MayDay::Abort("amr_mesh::parse_grid_generation - unknown box sorting method requested");
  }
}

void amr_mesh::parse_irreg_growth(){
  ParmParse pp("amr_mesh");
  pp.get("irreg_growth", m_irreg_growth);

  m_irreg_growth = std::max(m_irreg_growth, 0);
}

void amr_mesh::parse_blocking_factor(){
  ParmParse pp("amr_mesh");
  int blocking;
  pp.get("blocking_factor", blocking);
  if(blocking >= 4 && blocking % 2 == 0){
    m_blocking_factor = blocking;
  }
}

void amr_mesh::parse_eb_ghost(){

  ParmParse pp("amr_mesh");
  int ebghost;
  pp.get("eb_ghost", ebghost);
  if(ebghost >= 2){
    m_ebghost = ebghost;
  }
}

void amr_mesh::parse_num_ghost(){
  CH_TIME("amr_mesh::parse_num_ghost");

  ParmParse pp("amr_mesh");
  int ghost;
  pp.get("num_ghost", m_num_ghost);
  pp.get("lsf_ghost", m_lsf_ghost);
}

void amr_mesh::parse_redist_rad(){
  ParmParse pp("amr_mesh");
  int rad = 1;
  pp.get("redist_radius", rad);
  m_redist_rad = rad;
}

void amr_mesh::parse_centroid_stencils(){
  std::string str = "taylor";
  ParmParse pp("amr_mesh");
  pp.get("centroid_sten", str);


  // Maybe, in the future, we can change these but the user should not care about these (yet)
  m_centroid_sten_rad   = 1;
  m_centroid_sten_order = 1;

  if(str == "linear"){
    m_centroid_stencil = stencil_type::linear;
  }
  else if(str == "taylor"){
    m_centroid_stencil = stencil_type::taylor;
  }
  else if(str == "lsq"){
    m_centroid_stencil = stencil_type::lsq;
  }
  else if(str == "pwl"){
    m_centroid_stencil = stencil_type::pwl;
  }
  else{
    MayDay::Abort("amr_mesh::parse_centroid_stencils - unknown stencil requested");
  }
}

void amr_mesh::parse_eb_stencils(){
  std::string str = "taylor";
  ParmParse pp("amr_mesh");
  pp.get("eb_sten", str);
  
  // Maybe, in the future, we can change these but the user should not care about these (yet)
  m_eb_sten_rad   = 1;
  m_eb_sten_order = 1;

  if(str == "linear"){
    m_eb_stencil = stencil_type::linear;
  }
  else if(str == "taylor"){
    m_eb_stencil = stencil_type::taylor;
  }
  else if(str == "lsq"){
    m_eb_stencil = stencil_type::lsq;
  }
  else if(str == "pwl"){
    m_eb_stencil = stencil_type::pwl;
  }
  else{
    MayDay::Abort("amr_mesh::parse_eb_stencils - unknown stencil requested");
  }
}

void amr_mesh::set_irreg_sten_type(const stencil_type a_type){
  CH_TIME("amr_mesh::set_irreg_sten_type");
  m_stencil_type = a_type;

  
  std::string str = "taylor";
  ParmParse pp("amr_mesh");
  pp.query("stencil_type", str);
  if(str == "linear"){
    m_stencil_type = stencil_type::linear;
  }
  else if(str == "taylor"){
    m_stencil_type = stencil_type::taylor;
  }
  else if(str == "lsq"){
    m_stencil_type = stencil_type::lsq;
  }
}

void amr_mesh::set_irreg_sten_order(const int a_irreg_sten_order){
  CH_TIME("amr_mesh::irreg_sten_order");
  m_irreg_sten_order = a_irreg_sten_order;

  ParmParse pp("amr_mesh");
  int order = 1;
  pp.query("stencil_order", order);
  if(order == 1 || order == 2){
    m_irreg_sten_order = order;
  }
}

void amr_mesh::set_irreg_sten_radius(const int a_irreg_sten_radius){
  CH_TIME("amr_mesh::irreg_sten_radius");
  m_irreg_sten_radius = a_irreg_sten_radius;

  int radius = 1;
  ParmParse pp("amr_mesh");
  pp.query("stencil_radius", radius);
  if(radius == 1 || radius == 2){
    m_irreg_sten_radius = radius;
  }
}

void amr_mesh::sanity_check(){
  CH_TIME("amr_mesh::sanity_check");
  if(m_verbosity > 1){
    pout() << "amr_mesh::sanity_check" << endl;
  }

  CH_assert(m_max_amr_depth >= 0);
  for (int lvl = 0; lvl < m_ref_ratios.size(); lvl++){
    CH_assert(m_ref_ratios[lvl] == 2 || m_ref_ratios[lvl] == 4);
    CH_assert(m_blocking_factor >= 4 && m_blocking_factor % m_ref_ratios[lvl] == 0);
  }
  CH_assert(m_max_box_size >= 8 && m_max_box_size % m_blocking_factor == 0);
  CH_assert(m_fill_ratio > 0. && m_fill_ratio <= 1.0);
  CH_assert(m_buffer_size > 0);
}

bool amr_mesh::get_ebcf(){
  return m_ebcf;
}

RealVect amr_mesh::get_prob_lo(){
  return m_prob_lo;
}

RealVect amr_mesh::get_prob_hi(){
  return m_prob_hi;
}

int amr_mesh::get_finest_level(){
  return m_finest_level;
}

int amr_mesh::get_irreg_growth(){
  return m_irreg_growth;
}

int amr_mesh::get_max_amr_depth(){
  return m_max_amr_depth;
}

int amr_mesh::get_max_sim_depth(){
  return m_max_sim_depth;
}

int amr_mesh::get_num_ghost(){
  return m_num_ghost;
}

int amr_mesh::get_eb_ghost(){
  return m_ebghost;
}

int amr_mesh::get_redist_rad(){
  return m_redist_rad;
}

int amr_mesh::get_blocking_factor(){
  return m_blocking_factor;
}

int amr_mesh::get_max_box_size(){
  return m_max_box_size;
}

int amr_mesh::get_buffer(){
  return m_buffer_size;
}

int amr_mesh::get_max_ebis_box_size(){
  return m_max_ebis_box_size;
}

int amr_mesh::get_refinement_ratio(const int a_level1, const int a_level2){
  int coarLevel = Min(a_level1, a_level2);
  int fineLevel = Max(a_level1, a_level2);

  int ref = 1;
  for (int lvl = coarLevel; lvl < fineLevel; lvl++){
    ref = ref*m_ref_ratios[lvl];
  }

  return ref;
}

ProblemDomain amr_mesh::get_finest_domain(){
  return m_domains[m_max_amr_depth];
}

Real amr_mesh::get_finest_dx(){
  return m_dx[m_max_amr_depth];
}

Vector<Real>& amr_mesh::get_dx(){
  return m_dx;
}

Vector<int>& amr_mesh::get_ref_rat(){
  return m_ref_ratios;
}

const RefCountedPtr<BaseIF>& amr_mesh::get_baseif(const phase::which_phase a_phase){
  return m_baseif.at(a_phase);
}

Vector<IntVectSet> amr_mesh::get_irreg_tags() const {
  CH_TIME("amr_mesh::get_irreg_tags");
  if(m_verbosity > 5){
    pout() << "amr_mesh::get_irreg_tags" << endl;
  }

  Vector<IntVectSet> tags(m_max_amr_depth);

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_mfis->get_ebis(phase::solid);

  CH_assert(ebis_gas != NULL);

  for (int lvl = 0; lvl < m_max_amr_depth; lvl++){ // Don't need tags on maxdepth, we will never generate grids below that.
    const int which_level = ebis_gas->getLevel(m_domains[lvl]);

    tags[lvl] |= ebis_gas->irregCells(which_level);
    if(!ebis_sol.isNull()){
      tags[lvl] |= ebis_sol->irregCells(which_level);
    }
  }

  return tags;
}

Vector<ProblemDomain>& amr_mesh::get_domains(){
  return m_domains;
}

Vector<DisjointBoxLayout>& amr_mesh::get_proxy_grids(){
  return m_grids;
}

Vector<DisjointBoxLayout>& amr_mesh::get_grids(const std::string a_realm){
  if(!this->query_realm(a_realm)) {
    std::string str = "amr_mesh::get_grids - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }
  return m_realms[a_realm]->get_grids();
}

Vector<EBISLayout>& amr_mesh::get_ebisl(const std::string a_realm, const phase::which_phase a_phase){
  return m_realms[a_realm]->get_ebisl(a_phase);
}

Vector<RefCountedPtr<LayoutData<VoFIterator> > > amr_mesh::get_vofit(const std::string a_realm, const phase::which_phase a_phase){
  return m_realms[a_realm]->get_vofit(a_phase);
}

Vector<RefCountedPtr<LayoutData<Vector<LayoutIndex> > > >& amr_mesh::get_neighbors(const std::string a_realm,
										   const phase::which_phase a_phase){
  return m_realms[a_realm]->get_neighbors(a_phase);
}

AMRMask& amr_mesh::get_mask(const std::string a_mask, const int a_buffer, const std::string a_realm) {
  return m_realms[a_realm]->get_mask(a_mask, a_buffer);
}

Vector<RefCountedPtr<EBLevelGrid> >& amr_mesh::get_eblg(const std::string a_realm, const phase::which_phase a_phase){
  return m_realms[a_realm]->get_eblg(a_phase);
}

Vector<RefCountedPtr<MFLevelGrid> >& amr_mesh::get_mflg(const std::string a_realm){
  return m_realms[a_realm]->get_mflg();
}

EBAMRFAB& amr_mesh::get_levelset(const std::string a_realm, const phase::which_phase a_phase){
  return m_realms[a_realm]->get_levelset(a_phase);
}

Vector<RefCountedPtr<ebcoarseaverage> >& amr_mesh::get_coarave(const std::string a_realm, const phase::which_phase a_phase){
  return m_realms[a_realm]->get_coarave(a_phase);
}

Vector<RefCountedPtr<EBGhostCloud> >& amr_mesh::get_ghostcloud(const std::string a_realm, const phase::which_phase a_phase){
  return m_realms[a_realm]->get_ghostcloud(a_phase);
}

Vector<RefCountedPtr<nwoebquadcfinterp> >& amr_mesh::get_quadcfi(const std::string a_realm, const phase::which_phase a_phase){
  return m_realms[a_realm]->get_quadcfi(a_phase);
}

Vector<RefCountedPtr<EBQuadCFInterp> >& amr_mesh::get_old_quadcfi(const std::string a_realm, const phase::which_phase a_phase){
  return m_realms[a_realm]->get_old_quadcfi(a_phase);
}

Vector<RefCountedPtr<AggEBPWLFillPatch> >& amr_mesh::get_fillpatch(const std::string a_realm, const phase::which_phase a_phase){
  return m_realms[a_realm]->get_fillpatch(a_phase);
}

Vector<RefCountedPtr<EBPWLFineInterp> >& amr_mesh::get_eb_pwl_interp(const std::string a_realm, const phase::which_phase a_phase){
  return m_realms[a_realm]->get_eb_pwl_interp(a_phase);
}

Vector<RefCountedPtr<EBMGInterp> >& amr_mesh::get_eb_mg_interp(const std::string a_realm, const phase::which_phase a_phase){
  return m_realms[a_realm]->get_eb_mg_interp(a_phase);
}

Vector<RefCountedPtr<EBFluxRegister> >&  amr_mesh::get_flux_reg(const std::string a_realm, const phase::which_phase a_phase){
  return m_realms[a_realm]->get_flux_reg(a_phase);
}

Vector<RefCountedPtr<EBLevelRedist> >& amr_mesh::get_level_redist(const std::string a_realm, const phase::which_phase a_phase){
  return m_realms[a_realm]->get_level_redist(a_phase);
}

Vector<RefCountedPtr<EBCoarToFineRedist> >&  amr_mesh::get_coar_to_fine_redist(const std::string        a_realm,
									       const phase::which_phase a_phase){
  return m_realms[a_realm]->get_coar_to_fine_redist(a_phase);
}

Vector<RefCountedPtr<EBCoarToCoarRedist> >&  amr_mesh::get_coar_to_coar_redist(const std::string        a_realm,
									       const phase::which_phase a_phase){
  return m_realms[a_realm]->get_coar_to_coar_redist(a_phase);
}

Vector<RefCountedPtr<EBFineToCoarRedist> >&  amr_mesh::get_fine_to_coar_redist(const std::string        a_realm,
									       const phase::which_phase a_phase){
  return m_realms[a_realm]->get_fine_to_coar_redist(a_phase);
}

irreg_amr_stencil<centroid_interp>& amr_mesh::get_centroid_interp_stencils(const std::string        a_realm,
									   const phase::which_phase a_phase){
  return m_realms[a_realm]->get_centroid_interp_stencils(a_phase);
}

irreg_amr_stencil<eb_centroid_interp>& amr_mesh::get_eb_centroid_interp_stencils(const std::string        a_realm,
										 const phase::which_phase a_phase){
  return m_realms[a_realm]->get_eb_centroid_interp_stencils(a_phase);
}

irreg_amr_stencil<noncons_div>& amr_mesh::get_noncons_div_stencils(const std::string a_realm, const phase::which_phase a_phase){
  return m_realms[a_realm]->get_noncons_div_stencils(a_phase);
}

Vector<RefCountedPtr<Copier> >& amr_mesh::get_copier(const std::string a_realm, const phase::which_phase a_phase){
  return m_realms[a_realm]->get_copier(a_phase);
}

Vector<RefCountedPtr<Copier> >& amr_mesh::get_reverse_copier(const std::string a_realm, const phase::which_phase a_phase){
  return m_realms[a_realm]->get_reverse_copier(a_phase);
}

Vector<Box> amr_mesh::make_tiles(const Box a_box, const IntVect a_tilesize){
  MayDay::Abort("amr_mesh::make_Tiles - stop, where is this code used...?");
  // Modify the input tilesize to something sensible
  IntVect tilesize = a_tilesize;
  for (int dir = 0; dir < SpaceDim; dir++){
    if(tilesize[dir] > a_box.size(dir)){
      tilesize[dir] = a_box.size(dir);
    }
  }

  // Compute number of tiles along each coordinate
  Vector<int> num_tiles_dir(SpaceDim);
  for (int dir = 0; dir < SpaceDim; dir++){
    if(a_box.size(dir)%tilesize[dir] == 0){
      num_tiles_dir[dir] = a_box.size(dir)/tilesize[dir];
    }
    else{
      MayDay::Abort("amr_mesh::make_tiles - logic bust when building tiles. Tile size was not divisible by blocking factor");
    }
  }

  // Make the tiles
  Vector<Box> tiles;
#if CH_SPACEDIM==3
  for (int k = 0; k < num_tiles_dir[2]; k++){
#endif
    for (int j = 0; j < num_tiles_dir[1]; j++){
      for (int i = 0; i < num_tiles_dir[0]; i++){
	const IntVect ivlo = a_box.smallEnd() + IntVect(D_DECL(i,j,k))*tilesize;
	const IntVect ivhi = a_box.smallEnd() + IntVect(D_DECL(i+1,j+1,k+1))*tilesize - IntVect::Unit;
	tiles.push_back(Box(ivlo, ivhi));
      }
    }
#if CH_SPACEDIM==3
  }
#endif

  return tiles;
}

bool amr_mesh::query_realm(const std::string a_realm){
  CH_TIME("amr_mesh::query_realm");
  if(m_verbosity > 5){
    pout() << "amr_mesh::query_realm" << endl;
  }

  bool ret = true;
  
  if(m_realms.find(a_realm) == m_realms.end()){
    ret = false;
  }

  return ret;
}

void amr_mesh::register_realm(const std::string a_realm){
  CH_TIME("amr_mesh::register_realm");
  if(m_verbosity > 5){
    pout() << "amr_mesh::register_realm" << endl;
  }

  if(!this->query_realm(a_realm)){
    m_realms.emplace(a_realm, RefCountedPtr<realm> (new realm()));
  }
}

void amr_mesh::register_operator(const std::string a_operator, const std::string a_realm, const phase::which_phase a_phase){
  CH_TIME("amr_mesh::register_operator(operator, realm, phase)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::register_operator(operator, realm, phase)" << endl;
  }

  if(!this->query_realm(a_realm)) {
    std::string str = "amr_mesh::register_operator(operator, realm, phase) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  m_realms[a_realm]->register_operator(a_operator, a_phase);
}

void amr_mesh::register_mask(const std::string a_mask, const int a_buffer, const std::string a_realm){
  CH_TIME("amr_mesh::register_mask(mask, realm, buffer)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::register_mask(mask, realm, buffer)" << endl;
  }

  if(!this->query_realm(a_realm)) {
    std::string str = "amr_mesh::register_mask(mask, realm, buffer) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  m_realms[a_realm]->register_mask(a_mask, a_buffer);
}

void amr_mesh::define_realms(){
  CH_TIME("amr_mesh::define_realms()");
  if(m_verbosity > 5){
    pout() << "amr_mesh::define_realms()" << endl;
  }

  for (auto& r : m_realms){
    r.second->define(m_grids, m_domains, m_ref_ratios, m_dx, m_prob_lo, m_finest_level, m_ebghost, m_num_ghost, m_lsf_ghost, m_redist_rad,
		     m_ebcf, m_centroid_stencil, m_eb_stencil, m_baseif, m_mfis);
  }
}

void amr_mesh::regrid_realm(const std::string           a_realm,
			    const Vector<Vector<int> >& a_procs,
			    const Vector<Vector<Box> >& a_boxes,
			    const int                   a_lmin){
  CH_TIME("amr_mesh::regrid_realm(procs, boxes, level)");
  if(m_verbosity > 1){
    pout() << "amr_mesh::regrid_realm(procs, boxes, level)" << endl;
  }

  if(!this->query_realm(a_realm)) {
    std::string str = "amr_mesh::define_realm - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  // Make the dbl
  Vector<DisjointBoxLayout> grids(1 + m_finest_level);

  // Levels that didn't change. 
  for (int lvl = 0; lvl < a_lmin; lvl++){
    grids[lvl] = this->get_grids(a_realm)[lvl];
  }

  // Levels that did change. 
  for (int lvl = a_lmin; lvl <= m_finest_level; lvl++){
    grids[lvl] = DisjointBoxLayout();
    grids[lvl].define(a_boxes[lvl], a_procs[lvl], m_domains[lvl]);
    grids[lvl].close();
  }

  m_realms[a_realm]->define(grids, m_domains, m_ref_ratios, m_dx, m_prob_lo, m_finest_level, m_ebghost, m_num_ghost, m_lsf_ghost, m_redist_rad,
			    m_ebcf, m_centroid_stencil, m_eb_stencil, m_baseif, m_mfis);

  m_realms[a_realm]->regrid_base(a_lmin);
}

std::vector<std::string> amr_mesh::get_realms() const {
  std::vector<std::string> realms;

  for (const auto& r : m_realms){
    realms.push_back(r.first);
  }

  return realms;
}

box_sorting amr_mesh::get_box_sorting() const{
  return m_boxsort;
}
#include "CD_NamespaceFooter.H"
