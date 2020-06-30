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

#define AMR_MESH_DEBUG 0
#define USE_NEW_F2C_REDIST 1
#define USE_NEW_C2F_REDIST 1
#define USE_NEW_C2C_REDIST 1

amr_mesh::amr_mesh(){

  parse_verbosity();
  parse_coarsest_num_cells();
  parse_max_amr_depth();
  parse_max_simulation_depth();
  parse_refine_all_depth();
  parse_ebcf();
  parse_mg_coarsen();
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
  parse_load_balance();
  parse_domain();
  parse_ghost_interpolation();
  parse_centroid_stencils();
  parse_eb_stencils();


  m_finest_level = 0;
  m_has_grids    = false;
  m_has_mg_stuff = false;

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

  m_realm = RefCountedPtr<mf_realm> (new mf_realm());

  // Always do this stuff. 
  this->register_operator(s_eb_gradient,     phase::gas);
  this->register_operator(s_eb_gradient,     phase::solid);
  this->register_operator(s_eb_irreg_interp, phase::gas);
  this->register_operator(s_eb_irreg_interp, phase::solid);

}

amr_mesh::~amr_mesh(){
  
}



void amr_mesh::alias(EBAMRCellData&           a_data,
		     const phase::which_phase a_phase,
		     const MFAMRCellData&     a_mfdata,
		     const int                a_finest_level){
  CH_TIME("amr_mesh::alias(hardcap)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::alias(hardcap)" << endl;
  }

  for (int lvl = 0; lvl <= a_finest_level; lvl++){
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

void amr_mesh::allocate(AMRPVR& a_pvr, const int a_buffer){
  CH_TIME("amr_mesh::allocate(AMR PVR)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::allocate(AMR PVR)" << endl;
  }

  if(m_max_box_size != m_blocking_factor){
    MayDay::Abort("amr_mesh::allocate(particles) - only constant box sizes supported for particle methods");
  }
  
  a_pvr.resize(1 + m_finest_level);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_realm->get_grids()[lvl];
    const ProblemDomain& dom     = m_realm->get_domains()[lvl];

    const bool has_coar = lvl > 0;

    if(!has_coar){
      a_pvr[lvl] = RefCountedPtr<ParticleValidRegion> (new ParticleValidRegion(dbl, NULL, 1, 0));
    }
    else{
      const int ref_ratio = m_ref_ratios[lvl-1];
      a_pvr[lvl] = RefCountedPtr<ParticleValidRegion> (new ParticleValidRegion(dbl, a_pvr[lvl-1]->mask(), ref_ratio, a_buffer));
    }
  }
}

void amr_mesh::allocate(EBAMRCellData& a_data, const phase::which_phase a_phase, const int a_ncomp, const int a_ghost){
  CH_TIME("amr_mesh::allocate(cell)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::allocate(cell)" << endl;
  }

  const int ghost = (a_ghost == -1) ? m_num_ghost : a_ghost;

  a_data.resize(1 + m_finest_level);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_realm->get_grids()[lvl];
    const EBISLayout& ebisl      = m_realm->get_ebisl(a_phase)[lvl];
    
    EBCellFactory fact(ebisl);

    a_data[lvl] = RefCountedPtr<LevelData<EBCellFAB> >
      (new LevelData<EBCellFAB>(dbl, a_ncomp, ghost*IntVect::Unit, fact));

    EBLevelDataOps::setVal(*a_data[lvl], 0.0);
  }
}



void amr_mesh::allocate(EBAMRFluxData& a_data, const phase::which_phase a_phase, const int a_ncomp, const int a_ghost){
  CH_TIME("amr_mesh::allocate(flux)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::allocate(flux)" << endl;
  }

  const int ghost = (a_ghost == -1) ? m_num_ghost : a_ghost;

  a_data.resize(1 + m_finest_level);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_realm->get_grids()[lvl];
    const EBISLayout& ebisl      = m_realm->get_ebisl(a_phase)[lvl];

    EBFluxFactory fact(ebisl);

    a_data[lvl] = RefCountedPtr<LevelData<EBFluxFAB> >
      (new LevelData<EBFluxFAB>(dbl, a_ncomp, ghost*IntVect::Unit, fact));

    EBLevelDataOps::setVal(*a_data[lvl], 0.0);
  }
}

void amr_mesh::allocate(EBAMRIVData& a_data, const phase::which_phase a_phase, const int a_ncomp, const int a_ghost){
  CH_TIME("amr_mesh::allocate(baseiv)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::allocate(baseiv)" << endl;
  }

  const int ghost = (a_ghost == -1) ? m_num_ghost : a_ghost;

  a_data.resize(1 + m_finest_level);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_realm->get_grids()[lvl];
    const EBISLayout& ebisl      = m_realm->get_ebisl(a_phase)[lvl];
    const ProblemDomain& domain  = m_realm->get_domains()[lvl];
    
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
}

void amr_mesh::allocate(EBAMRIFData& a_data, const phase::which_phase a_phase, const int a_ncomp, const int a_ghost){
  CH_TIME("amr_mesh::allocate(EBAMRIFData)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::allocate(EBAMRIFData)" << endl;
  }

  const int ghost = (a_ghost == -1) ? m_num_ghost : a_ghost;

  a_data.resize(1 + m_finest_level);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_realm->get_grids()[lvl];
    const EBISLayout& ebisl      = m_realm->get_ebisl(a_phase)[lvl];
    const ProblemDomain& domain  = m_realm->get_domains()[lvl];
    
    DomainFluxIFFABFactory fact(ebisl, domain);

    a_data[lvl] = RefCountedPtr<LevelData<DomainFluxIFFAB> >
      (new LevelData<DomainFluxIFFAB>(dbl, a_ncomp, ghost*IntVect::Unit, fact));
  }
}

void amr_mesh::allocate(EBAMRBool& a_data, const int a_ncomp, const int a_ghost){
  CH_TIME("amr_mesh::allocate(EBAMRBool)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::allocate(EBAMRBool)" << endl;
  }

  a_data.resize(1 + m_finest_level);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_realm->get_grids()[lvl];
    
    a_data[lvl] = RefCountedPtr<LevelData<BaseFab<bool> > >
      (new LevelData<BaseFab<bool> >(dbl, a_ncomp, a_ghost*IntVect::Unit));
  }
}

void amr_mesh::allocate(MFAMRCellData& a_data, const int a_ncomp, const int a_ghost){
  CH_TIME("amr_mesh::allocate(mf cell)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::allocate(mf cell)" << endl;
  }

  const int ghost   = (a_ghost == -1) ? m_num_ghost : a_ghost;
  const int ignored = a_ncomp;
  const int nphases = m_mfis->num_phases();

  a_data.resize(1 + m_finest_level);

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_realm->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_realm->get_ebis(phase::solid);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_realm->get_grids()[lvl];
    
    Vector<EBISLayout> ebisl(nphases);
    Vector<int>        comps(nphases, a_ncomp);

    ebisl[phase::gas]   = m_realm->get_ebisl(phase::gas)[lvl];
    if(!ebis_sol.isNull()){ // Single phase
      ebisl[phase::solid] = m_realm->get_ebisl(phase::solid)[lvl];
    }
    
    MFCellFactory factory(ebisl, comps);

    a_data[lvl] = RefCountedPtr<LevelData<MFCellFAB> >
      (new LevelData<MFCellFAB>(dbl, ignored, ghost*IntVect::Unit, factory));

    MFLevelDataOps::setVal(*a_data[lvl], 0.0);
  }
  
}

void amr_mesh::allocate(MFAMRFluxData& a_data, const int a_ncomp, const int a_ghost){
  CH_TIME("amr_mesh::allocate(mf flux)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::allocate(mf flux)" << endl;
  }

  const int ghost   = (a_ghost == -1) ? m_num_ghost : a_ghost;
  const int ignored = a_ncomp;
  const int nphases = m_mfis->num_phases();

  a_data.resize(1 + m_finest_level);

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_realm->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_realm->get_ebis(phase::solid);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_realm->get_grids()[lvl];
    
    Vector<EBISLayout> ebisl(nphases);
    Vector<int>        comps(nphases, a_ncomp);

    ebisl[phase::gas]   = m_realm->get_ebisl(phase::gas)[lvl];
    if(!ebis_sol.isNull()){
      ebisl[phase::solid] = m_realm->get_ebisl(phase::solid)[lvl];
    }
    
    MFFluxFactory factory(ebisl, comps);

    a_data[lvl] = RefCountedPtr<LevelData<MFFluxFAB> >
      (new LevelData<MFFluxFAB>(dbl, ignored, ghost*IntVect::Unit, factory));
  }
}
  
void amr_mesh::allocate(MFAMRIVData& a_data, const int a_ncomp, const int a_ghost){
  CH_TIME("amr_mesh::allocate(mf baseivfab)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::allocate(mf baseivfab)" << endl;
  }

  const int ghost   = (a_ghost == -1) ? m_num_ghost : a_ghost;
  const int ignored = a_ncomp;
  const int nphases = m_mfis->num_phases();

  a_data.resize(1 + m_finest_level);

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_mfis->get_ebis(phase::solid);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    Vector<EBISLayout> ebisl(nphases);
    Vector<int>        comps(nphases, a_ncomp);

    ebisl[phase::gas]   = this->get_ebisl(phase::gas)[lvl];
    if(!ebis_sol.isNull()){
      ebisl[phase::solid] = this->get_ebisl(phase::solid)[lvl];
    }
    
    MFBaseIVFABFactory factory(ebisl, comps);

    a_data[lvl] = RefCountedPtr<LevelData<MFBaseIVFAB> >
      (new LevelData<MFBaseIVFAB>(m_grids[lvl], ignored, ghost*IntVect::Unit, factory));
  }
}


void amr_mesh::reallocate(EBAMRCellData& a_data, const phase::which_phase a_phase, const int a_lmin){
  CH_TIME("amr_mesh::reallocate(cell)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::reallocate(cell)" << endl;
  }

  const IntVect ghost = a_data[0]->ghostVect();
  const int ncomp     = a_data[0]->nComp();

  a_data.resize(1 + m_finest_level);

  const int lmin = Max(0, a_lmin+1);
  
  for (int lvl = a_lmin; lvl <= m_finest_level; lvl++){
    EBCellFactory fact(m_ebisl[a_phase][lvl]);

    a_data[lvl] = RefCountedPtr<LevelData<EBCellFAB> >
      (new LevelData<EBCellFAB>(m_grids[lvl], ncomp, ghost, fact));

    EBLevelDataOps::setVal(*a_data[lvl], 0.0);
  }
}

void amr_mesh::reallocate(EBAMRFluxData& a_data, const phase::which_phase a_phase, const int a_lmin){
  CH_TIME("amr_mesh::reallocate(ebamrflux)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::reallocate(ebamrflux)" << endl;
  }

  const IntVect ghost = a_data[0]->ghostVect();
  const int ncomp     = a_data[0]->nComp();

  a_data.resize(1 + m_finest_level);

  const int lmin = Max(0, a_lmin+1);
  
  for (int lvl = a_lmin; lvl <= m_finest_level; lvl++){
    EBFluxFactory fact(m_ebisl[a_phase][lvl]);

    a_data[lvl] = RefCountedPtr<LevelData<EBFluxFAB> > (new LevelData<EBFluxFAB>(m_grids[lvl], ncomp, ghost, fact));

    EBLevelDataOps::setVal(*a_data[lvl], 0.0);
  }
}

void amr_mesh::reallocate(EBAMRIVData& a_data, const phase::which_phase a_phase, const int a_lmin){
  CH_TIME("amr_mesh::allocate(baseiv)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::allocate(baseiv)" << endl;
  }

  const IntVect ghost = a_data[0]->ghostVect();
  const int ncomp     = a_data[0]->nComp();

  a_data.resize(1 + m_finest_level);

  const int lmin = Max(0, a_lmin+1);

  for (int lvl = a_lmin; lvl <= m_finest_level; lvl++){

    LayoutData<IntVectSet> irreg_sets(m_grids[lvl]);
    for (DataIterator dit = m_grids[lvl].dataIterator(); dit.ok(); ++dit){
      Box box = m_grids[lvl].get(dit());
      box.grow(ghost);
      box &= m_domains[lvl];

      irreg_sets[dit()] = m_ebisl[a_phase][lvl][dit()].getIrregIVS(box);
    }

    BaseIVFactory<Real> fact(m_ebisl[a_phase][lvl], irreg_sets);

    a_data[lvl] = RefCountedPtr<LevelData<BaseIVFAB<Real> > >
      (new LevelData<BaseIVFAB<Real> >(m_grids[lvl], ncomp, ghost, fact));

    EBLevelDataOps::setVal(*a_data[lvl], 0.0);
  }
}

void amr_mesh::reallocate(EBAMRIFData& a_data, const phase::which_phase a_phase, const int a_lmin){
  CH_TIME("amr_mesh::reallocate(EBAMRIFData)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::reallocate(EBAMRIFData)" << endl;
  }

  const IntVect ghost = a_data[0]->ghostVect();
  const int ncomp     = a_data[0]->nComp();

  a_data.resize(1 + m_finest_level);

  const int lmin = Max(0, a_lmin+1);

  for (int lvl = a_lmin; lvl <= m_finest_level; lvl++){

    DomainFluxIFFABFactory fact(m_ebisl[a_phase][lvl],m_domains[lvl]);

    a_data[lvl] = RefCountedPtr<LevelData<DomainFluxIFFAB> >
      (new LevelData<DomainFluxIFFAB>(m_grids[lvl], ncomp, ghost, fact));
  }
}

void amr_mesh::reallocate(EBAMRBool& a_data, const int a_ncomp, const int a_lmin){
  CH_TIME("amr_mesh::reallocate(EBAMRBool)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::reallocate(EBAMRBool)" << endl;
  }
  
  const IntVect ghost = a_data[0]->ghostVect();
  const int ncomp     = a_data[0]->nComp();

  a_data.resize(1 + m_finest_level);

  const int lmin = Max(0, a_lmin+1);

  for (int lvl = a_lmin; lvl <= m_finest_level; lvl++){

    a_data[lvl] = RefCountedPtr<LevelData<BaseFab<bool> > >
      (new LevelData<BaseFab<bool> >(m_grids[lvl], ncomp, ghost));
  }
}

void amr_mesh::reallocate(MFAMRCellData& a_data, const int a_lmin){
  CH_TIME("amr_mesh::reallocate(mf cell)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::reallocate(mf cell)" << endl;
  }

  const IntVect ghost = a_data[0]->ghostVect();
  const int ncomp = a_data[0]->nComp();
  const int ignored = ncomp;
  const int nphases = m_mfis->num_phases();

  a_data.resize(1 + m_finest_level);

  const int lmin = Max(0, a_lmin+1);

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_mfis->get_ebis(phase::solid);

  for (int lvl = a_lmin; lvl <= m_finest_level; lvl++){
    Vector<EBISLayout> ebisl(nphases);
    Vector<int>        comps(nphases, ncomp);

    ebisl[phase::gas]   = this->get_ebisl(phase::gas)[lvl];
    if(!ebis_sol.isNull()){ // Single phase
      ebisl[phase::solid] = this->get_ebisl(phase::solid)[lvl];
    }
    
    MFCellFactory factory(ebisl, comps);

    a_data[lvl] = RefCountedPtr<LevelData<MFCellFAB> >
      (new LevelData<MFCellFAB>(m_grids[lvl], ignored, ghost, factory));

    MFLevelDataOps::setVal(*a_data[lvl], 0.0);
  }
}

void amr_mesh::reallocate(MFAMRFluxData& a_data, const int a_lmin){
  CH_TIME("amr_mesh::allocate(mf flux)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::allocate(mf flux)" << endl;
  }

  const IntVect ghost = a_data[0]->ghostVect();
  const int ncomp     = a_data[0]->nComp();
  const int ignored   = ncomp;
  const int nphases   = m_mfis->num_phases();

  a_data.resize(1 + m_finest_level);

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_mfis->get_ebis(phase::solid);

  const int lmin = Max(0, a_lmin+1);

  for (int lvl = a_lmin; lvl <= m_finest_level; lvl++){
    Vector<EBISLayout> ebisl(nphases);
    Vector<int>        comps(nphases, ncomp);

    ebisl[phase::gas]   = this->get_ebisl(phase::gas)[lvl];
    if(!ebis_sol.isNull()){
      ebisl[phase::solid] = this->get_ebisl(phase::solid)[lvl];
    }
    
    MFFluxFactory factory(ebisl, comps);

    a_data[lvl] = RefCountedPtr<LevelData<MFFluxFAB> >
      (new LevelData<MFFluxFAB>(m_grids[lvl], ignored, ghost, factory));
  }
}

void amr_mesh::reallocate(MFAMRIVData& a_data, const int a_lmin){
  CH_TIME("amr_mesh::allocate(mf baseivfab)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::allocate(mf baseivfab)" << endl;
  }

  const IntVect ghost = a_data[0]->ghostVect();
  const int ncomp     = a_data[0]->nComp();
  const int ignored   = ncomp;
  const int nphases   = m_mfis->num_phases();

  a_data.resize(1 + m_finest_level);

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_mfis->get_ebis(phase::solid);

  const int lmin = Max(0, a_lmin+1);

  for (int lvl = a_lmin; lvl <= m_finest_level; lvl++){
    Vector<EBISLayout> ebisl(nphases);
    Vector<int>        comps(nphases, ncomp);

    ebisl[phase::gas]   = this->get_ebisl(phase::gas)[lvl];
    if(!ebis_sol.isNull()){
      ebisl[phase::solid] = this->get_ebisl(phase::solid)[lvl];
    }
    
    MFBaseIVFABFactory factory(ebisl, comps);

    a_data[lvl] = RefCountedPtr<LevelData<MFBaseIVFAB> >
      (new LevelData<MFBaseIVFAB>(m_grids[lvl], ignored, ghost, factory));
  }
}

void amr_mesh::allocate_interface(EBAMRIVData& a_data, phase::which_phase a_phase, const int a_ncomp, const int a_ghost){
  CH_TIME("amr_mesh::allocate_interface(baseiv)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::allocate_interface(baseiv)" << endl;
  }

  MayDay::Abort("amr_mesh::allocate_interface - this is not the way to do it. Please use MFLevelGrid");
  const int ghost = (a_ghost == -1) ? m_num_ghost : a_ghost;

  a_data.resize(1 + m_finest_level);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){

    const IntVectSet interface_ivs = m_mfis->interface_region(m_domains[lvl]);
    LayoutData<IntVectSet> irreg_sets(m_grids[lvl]);
    for (DataIterator dit = m_grids[lvl].dataIterator(); dit.ok(); ++dit){
      Box box = m_grids[lvl].get(dit());
      box.grow(ghost);
      box &= m_domains[lvl];

      irreg_sets[dit()] = interface_ivs & box; 
    }

    BaseIVFactory<Real> fact(m_ebisl[a_phase][lvl], irreg_sets);

    a_data[lvl] = RefCountedPtr<LevelData<BaseIVFAB<Real> > >
      (new LevelData<BaseIVFAB<Real> >(m_grids[lvl], a_ncomp, ghost*IntVect::Unit, fact));
  }
}

void amr_mesh::set_mfis(const RefCountedPtr<mfis>& a_mfis){
  CH_TIME("amr_mesh::set_mfis");
  if(m_verbosity > 5){
    pout() << "amr_mesh::set_mfis" << endl;
  }

  m_mfis = a_mfis;
}

void amr_mesh::parse_domain(){

  std::string balance;
  ParmParse pp("amr_mesh");

  Vector<Real> v(SpaceDim);
  pp.getarr("lo_corner", v, 0, SpaceDim); m_prob_lo = RealVect(D_DECL(v[0], v[1], v[2]));
  pp.getarr("hi_corner", v, 0, SpaceDim); m_prob_hi = RealVect(D_DECL(v[0], v[1], v[2]));
}

void amr_mesh::parse_load_balance(){

  std::string balance;
  ParmParse pp("amr_mesh");
  pp.get("load_balance", balance);

  if(balance == "volume"){
    m_which_balance = load_balance::volume;
  }
  else if(balance == "elliptic"){
    m_which_balance = load_balance::elliptic;
  }
  else{
    MayDay::Abort("amr_mesh::parse_load_balance - unknown load balance requested");
  }
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
  m_mflg.resize(nlevels);
  m_neighbors.resize(nlevels);

  m_eblg.resize(phase::num_phases);
  m_ebisl.resize(phase::num_phases);
  m_vofiter.resize(phase::num_phases);
  m_coarave.resize(phase::num_phases);
  m_quadcfi.resize(phase::num_phases);
  m_flux_reg.resize(phase::num_phases);
  m_old_quadcfi.resize(phase::num_phases);
  m_level_redist.resize(phase::num_phases);
  m_concentration_redist.resize(phase::num_phases);
  m_pwl_fillpatch.resize(phase::num_phases);
  m_pwl_interp.resize(phase::num_phases);
  m_ebmg_interp.resize(phase::num_phases);
  m_centroid_interp.resize(phase::num_phases);
  m_noncons_div.resize(phase::num_phases);
  m_eb_centroid_interp.resize(phase::num_phases);
  m_coar_to_fine_redist.resize(phase::num_phases);
  m_coar_to_coar_redist.resize(phase::num_phases);
  m_fine_to_coar_redist.resize(phase::num_phases);
  m_copier.resize(phase::num_phases);
  m_ghostclouds.resize(phase::num_phases);
  m_reverse_copier.resize(phase::num_phases);

  m_eblg[phase::gas].resize(nlevels);
  m_ebisl[phase::gas].resize(nlevels);
  m_coarave[phase::gas].resize(nlevels);
  m_vofiter[phase::gas].resize(nlevels);
  m_quadcfi[phase::gas].resize(nlevels);
  m_flux_reg[phase::gas].resize(nlevels);
  m_old_quadcfi[phase::gas].resize(nlevels);
  m_pwl_fillpatch[phase::gas].resize(nlevels);
  m_pwl_interp[phase::gas].resize(nlevels);
  m_ebmg_interp[phase::gas].resize(nlevels);
  m_level_redist[phase::gas].resize(nlevels);
  m_coar_to_fine_redist[phase::gas].resize(nlevels);
  m_coar_to_coar_redist[phase::gas].resize(nlevels);
  m_fine_to_coar_redist[phase::gas].resize(nlevels);
  m_concentration_redist[phase::gas].resize(nlevels);
  m_copier[phase::gas].resize(nlevels);
  m_reverse_copier[phase::gas].resize(nlevels);
  m_ghostclouds[phase::gas].resize(nlevels);

  m_eblg[phase::solid].resize(nlevels);
  m_ebisl[phase::solid].resize(nlevels);
  m_coarave[phase::solid].resize(nlevels);
  m_quadcfi[phase::solid].resize(nlevels);
  m_flux_reg[phase::solid].resize(nlevels);
  m_vofiter[phase::solid].resize(nlevels);
  m_old_quadcfi[phase::solid].resize(nlevels);
  m_level_redist[phase::solid].resize(nlevels);
  m_pwl_fillpatch[phase::solid].resize(nlevels);
  m_pwl_interp[phase::solid].resize(nlevels);
  m_ebmg_interp[phase::solid].resize(nlevels);
  m_coar_to_fine_redist[phase::solid].resize(nlevels);
  m_coar_to_coar_redist[phase::solid].resize(nlevels);
  m_fine_to_coar_redist[phase::solid].resize(nlevels);
  m_concentration_redist[phase::solid].resize(nlevels);
  m_copier[phase::solid].resize(nlevels);
  m_reverse_copier[phase::solid].resize(nlevels);
  m_ghostclouds[phase::solid].resize(nlevels);

  m_dx[0] = (m_prob_hi[0] - m_prob_lo[0])/m_num_cells[0];
  m_domains[0] = ProblemDomain(IntVect::Zero, m_num_cells - IntVect::Unit);

  m_gradsten_gas.resize(nlevels);
  m_gradsten_sol.resize(nlevels);

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
  CH_TIME("amr_mesh::regrid_amr");
  if(m_verbosity > 1){
    pout() << "amr_mesh::regrid_amr" << endl;
  }
  
  Vector<IntVectSet> tags = a_tags; // build_grids destroys tags, so copy them

  // This is stuff that always gets done
  this->build_grids(tags, a_lmin, a_lmax, a_hardcap);
  this->define_neighbors(a_lmin);
  this->define_eblevelgrid(a_lmin);  // Define EBLevelGrid objects on both phases
  this->define_vofiter(a_lmin);
  this->define_mflevelgrid(a_lmin);  // Define MFLevelGrid
  if(!m_has_mg_stuff){ // Define MG stuff
    this->define_mg_stuff();
    m_has_mg_stuff = true; // Only needs to be done ONCE per run. 
  }

  this->define_realms();
  m_realm->regrid_base(a_lmin);
}

void amr_mesh::regrid_operators(const int a_lmin,
				const int a_lmax,
				const int a_regsize){
  // Now allocate operators
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

  m_realm->regrid_operators(a_lmin, a_lmax, a_regsize);
}

void amr_mesh::build_grids(Vector<IntVectSet>& a_tags, const int a_lmin, const int a_lmax, const int a_hardcap){
  CH_TIME("amr_mesh::build_grids");
  if(m_verbosity > 2){
    pout() << "amr_mesh::build_grids" << endl;
  }

  // TLDR: a_lmin is the coarsest level that changes. A special condition is a_lmin=0 for which we assume
  //       that there are no prior grids.
  //       a_lmax is the finest level that changes. This is typically 

#if 0 // Original code
  const int base      = 0;                                    // Base level never changes.
#else // Debug code
  const int base = Max(0, a_lmin - 1);
#endif
  const int top_level = (m_finest_level == m_max_amr_depth) ? // top_level is the finest level where we have tags. We should never
    m_finest_level - 1 : a_tags.size() - 1;                   // have tags on max_amr_depth, and we make that restriction here. 
  Vector<Vector<Box> > new_boxes(1 + top_level);  // New boxes to be load balance
  Vector<Vector<Box> > old_boxes(1 + top_level);  // Old grids.

#if 0 // Debug
  pout() << "amr_mesh::build_grids - " << "finest_level = " << m_finest_level
	 << "\t max_depth = " << m_max_amr_depth
    	 << "\t top_level = " << top_level
	 << endl;
#endif

  const int hardcap = (a_hardcap == -1) ? m_max_amr_depth : a_hardcap;

  if(m_max_amr_depth > 0 && hardcap > 0){
    domainSplit(m_domains[0], old_boxes[0], m_max_box_size, m_blocking_factor);

    if(!m_has_grids){

      for (int lvl = 1; lvl <= top_level; lvl++){
#if 1
	domainSplit(m_domains[lvl], old_boxes[lvl], m_max_box_size, m_blocking_factor);
#else
	old_boxes[lvl].resize(0);
	old_boxes[lvl].push_back(m_domains[lvl].domainBox());
#endif
      }
    }
    else{
      for (int lvl = 1; lvl <= top_level; lvl++){ 
	old_boxes[lvl] = m_grids[lvl].boxArray();
      }
    }

    int base_level = m_refine_all_depth;
    if(top_level >= base_level){ // Use tags for regridding
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
    else{ // Tag depth is below uniform refinement depth, create uniform grids everywhere
      m_finest_level = m_refine_all_depth;
      m_finest_level = Min(m_finest_level, m_max_amr_depth); // Don't exceed m_max_amr_depth
      m_finest_level = Min(m_finest_level, m_max_sim_depth); // Don't exceed maximum simulation depth
      m_finest_level = Min(m_finest_level, hardcap);         // Don't exceed hardcap

      new_boxes.resize(1 + m_finest_level);
      for (int lvl = 0; lvl <= m_finest_level; lvl++){
	domainSplit(m_domains[lvl], new_boxes[lvl], m_max_box_size, m_blocking_factor);
      }
    }
  }
  else{
    new_boxes.resize(1);
    domainSplit(m_domains[0], new_boxes[0], m_max_box_size, m_blocking_factor);
    
    m_finest_level = 0;
  }

  // Do Morton ordering
  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    mortonOrdering((Vector<Box>&)new_boxes[lvl]);
  }

  // Load balance boxes
  Vector<Vector<int> > proc_assign(1 + m_finest_level);
  this->loadbalance(proc_assign, new_boxes);

  // Define grids. If a_lmin=0 every grid is new, otherwise keep old grids up to but not including a_lmin
  if(a_lmin == 0){
    for (int lvl = 0; lvl <= m_finest_level; lvl++){
      m_grids[lvl] = DisjointBoxLayout();
      m_grids[lvl].define(new_boxes[lvl], proc_assign[lvl], m_domains[lvl]);
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
      m_grids[lvl].define(new_boxes[lvl], proc_assign[lvl], m_domains[lvl]);
      m_grids[lvl].close();
    }
  }

  m_has_grids = true;
}

void amr_mesh::define_mg_stuff(){
  CH_TIME("amr_mesh::define_mg_stuff");
  if(m_verbosity > 2){
    pout() << "amr_mesh::define_mg_stuff" << endl;
  }

  const RefCountedPtr<EBIndexSpace>& ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace>& ebis_sol = m_mfis->get_ebis(phase::solid);

  const int coar_ref = 2;

  // Redo these
  m_mg_domains.resize(0);
  m_mg_grids.resize(0);
  m_mg_mflg.resize(0);
  
  m_mg_eblg.resize(phase::num_phases);

  m_mg_eblg[phase::gas].resize(0);
  m_mg_eblg[phase::solid].resize(0);

  int num_coar       = 0;
  bool has_coar      = true;
  ProblemDomain fine = m_domains[0];

  // Coarsen problem domains and create grids
  while(num_coar < m_mg_coarsen || !has_coar){

    // Check if we can coarsen
    const ProblemDomain coar = fine.coarsen(coar_ref);
    const Box coar_box       = coar.domainBox();
    for (int dir = 0; dir < SpaceDim; dir++){
      if(coar_box.size()[dir] < m_max_box_size || coar_box.size()[dir]%m_max_box_size != 0){
	has_coar = false;
      }
    }

    if(has_coar){
      // Split the domain into pieces, then order and load balance them
      Vector<Box> boxes;
      Vector<int> proc_assign;
      domainSplit(coar, boxes, m_max_box_size, m_blocking_factor);
      mortonOrdering(boxes);
      load_balance::balance_volume(proc_assign, boxes);

      // Add problem domain and grid
      m_mg_domains.push_back(coar);
      m_mg_grids.push_back(DisjointBoxLayout(boxes, proc_assign, coar));

      // Define the EBLevelGrids
      const int idx = m_mg_grids.size() - 1; // Last element added
      if(!ebis_gas.isNull()){
	m_mg_eblg[phase::gas].push_back(RefCountedPtr<EBLevelGrid> (new EBLevelGrid(m_mg_grids[idx],
										    m_mg_domains[idx],
										    m_ebghost,
										    ebis_gas)));
      }
      if(!ebis_sol.isNull()){
	m_mg_eblg[phase::solid].push_back(RefCountedPtr<EBLevelGrid> (new EBLevelGrid(m_mg_grids[idx],
										      m_mg_domains[idx],
										      m_ebghost,
										      ebis_sol)));
      }

      // Define the MFLevelGrid object
      Vector<EBLevelGrid> eblgs;
      if(!ebis_gas.isNull()){
	eblgs.push_back(*m_mg_eblg[phase::gas][idx]);
      }
      if(!ebis_sol.isNull()){
	eblgs.push_back(*m_mg_eblg[phase::solid][idx]);
      }
      m_mg_mflg.push_back(RefCountedPtr<MFLevelGrid> (new MFLevelGrid(m_mfis, eblgs)));
      

      // Next iterate
      fine = coar;
      num_coar++;
    }
    else{
      break;
    }
  }
}

void amr_mesh::loadbalance(Vector<Vector<int> >& a_procs, Vector<Vector<Box> >& a_boxes){
  CH_TIME("amr_mesh::loadbalance");
  if(m_verbosity > 5){
    pout() << "amr_mesh::loadbalance" << endl;
  }

  // Level-by-level load balancing. This might change in the future. 
  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    if(m_which_balance == load_balance::volume){
      load_balance::balance_volume(a_procs[lvl], a_boxes[lvl]);
    }
    else if(m_which_balance == load_balance::elliptic){
      load_balance::balance_elliptic(a_procs[lvl], a_boxes[lvl], m_mfis->get_ebis(phase::gas), m_domains[lvl], false);
    }
    else if(m_which_balance == load_balance::multifluid){
      load_balance::balance_multifluid(a_procs[lvl], a_boxes[lvl], m_mfis, m_domains[lvl], false);
    }
  }
}

void amr_mesh::compute_gradient(LevelData<EBCellFAB>& a_gradient,
				const LevelData<EBCellFAB>& a_phi,
				const phase::which_phase a_phase,
				const int a_lvl){
  CH_TIME("amr_mesh::compute_gradient(level)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::compute_gradient(level)" << endl;
  }

  CH_assert(a_phi.nComp()      == 1);
  CH_assert(a_gradient.nComp() == SpaceDim);

    
  const int comp  = 0;
  const int ncomp = 1;
    
  const Real& dx = m_dx[a_lvl];
  const DisjointBoxLayout& dbl = m_grids[a_lvl]; // Doing this since I assume everything is defined over m_grids
  const ProblemDomain& domain  = m_domains[a_lvl];

  LayoutData<IntVectSet> cfivs;
  EBArith::defineCFIVS(cfivs, dbl, domain);

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


#if 0 // Old code
    // We can't REALLY trust ghost cells on the boundary. Do the boundary cells using safer stencils.
    IntVectSet bndry_ivs = ebisbox.getIrregIVS(dbl.get(dit()));
    for (int dir = 0; dir < SpaceDim; dir++){
      Box lo_box, hi_box;
      int has_lo, has_hi;

      EBArith::loHi(lo_box, has_lo, hi_box, has_hi, domain, region, dir);

      if(has_lo){
	bndry_ivs |= IntVectSet(lo_box);
      }
      if(has_hi){
	bndry_ivs |= IntVectSet(hi_box);
      }
    }

      
    // Compute stencils for irregular cells
    for (VoFIterator vofit(bndry_ivs, ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();

      for (int dir = 0; dir < SpaceDim; dir++){
	  
	grad(vof, dir) = 0.;

	VoFStencil sten;
	EBArith::getFirstDerivStencil(sten, vof, ebisbox, dir, dx, &cfivs[dit()], 0);
	for (int i = 0; i < sten.size(); i++){
	  const VolIndex& ivof = sten.vof(i);
	  const Real& iweight  = sten.weight(i);
	    
	  grad(vof, dir) += phi(ivof, comp)*iweight;
	}
      }
    }
#else // New code
    BaseIVFAB<VoFStencil>* grad_stencils;
    if(a_phase == phase::gas){
      grad_stencils = &((*m_gradsten_gas[a_lvl])[dit()]);
    }
    else if(a_phase == phase::solid){
      grad_stencils = &((*m_gradsten_sol[a_lvl])[dit()]);
    }


    for (VoFIterator vofit(grad_stencils->getIVS(), ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof    = vofit();
      const VoFStencil& sten = (*grad_stencils)(vof, 0);


      for (int dir = 0; dir < SpaceDim; dir++){
	grad(vof, dir) = 0.0;
      }

      for (int i = 0; i < sten.size(); i++){
	const VolIndex& ivof = sten.vof(i);
	const Real& iweight  = sten.weight(i);
	const int ivar       = sten.variable(i);

	grad(vof, ivar) += phi(ivof, comp)*iweight;
      }
    }
#endif

    // Set covered to zero
    for (int dir= 0; dir < SpaceDim; dir++){
      grad.setCoveredCellVal(0.0, dir);
    }
  }
}

void amr_mesh::compute_gradient(EBAMRCellData& a_gradient, const EBAMRCellData& a_phi, const phase::which_phase a_phase){
  CH_TIME("amr_mesh::compute_gradient(ebamrcell)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::compute_gradient(ebamrcell)" << endl;
  }

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    compute_gradient(*a_gradient[lvl], *a_phi[lvl], a_phase, lvl);
  }
}

void amr_mesh::compute_gradient(MFAMRCellData& a_gradient, const MFAMRCellData& a_phi){
  CH_TIME("amr_mesh::compute_gradient(mfamrcell)");
  if(m_verbosity > 5){
    pout() << "amr_mesh::compute_gradient(mfamrcell)" << endl;
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
      this->compute_gradient(alias_grad, alias_phi, phase::gas);
    }
    else if(iphase == 1){
      this->compute_gradient(alias_grad, alias_phi, phase::solid);
    }
  }
}

void amr_mesh::define_gradsten(const int a_lmin){
  CH_TIME("amr_mesh::define_gradsten");
  if(m_verbosity > 2){
    pout() << "amr_mesh::define_gradsten" << endl;
  }

  const RefCountedPtr<EBIndexSpace>& ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace>& ebis_sol = m_mfis->get_ebis(phase::solid);

  for (int lvl = a_lmin; lvl <= m_finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_grids[lvl];
    const ProblemDomain& domain  = m_domains[lvl];
    const Real dx                = m_dx[lvl];
    
    if(!ebis_gas.isNull()){
      if(this->query_operator(s_eb_gradient, phase::gas)){
	m_gradsten_gas[lvl] = RefCountedPtr<LayoutData<BaseIVFAB<VoFStencil> > > (new LayoutData<BaseIVFAB<VoFStencil> >(dbl));

	const EBISLayout& ebisl = m_ebisl[phase::gas][lvl];

	LayoutData<IntVectSet>& cfivs = *(m_eblg[phase::gas][lvl]->getCFIVS());
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

	  BaseIVFAB<VoFStencil>& vofstencils = (*m_gradsten_gas[lvl])[dit()];
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
    
    if(!ebis_sol.isNull()){
      if(this->query_operator(s_eb_gradient, phase::solid)){
	m_gradsten_sol[lvl] = RefCountedPtr<LayoutData<BaseIVFAB<VoFStencil> > > (new LayoutData<BaseIVFAB<VoFStencil> >(dbl));

	const EBISLayout& ebisl = m_ebisl[phase::solid][lvl];

	LayoutData<IntVectSet>& cfivs = *(m_eblg[phase::solid][lvl]->getCFIVS());
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

	  BaseIVFAB<VoFStencil>& vofstencils = (*m_gradsten_sol[lvl])[dit()];
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
}

void amr_mesh::define_neighbors(const int a_lmin){
  CH_TIME("amr_mesh::define_neighbors");
  if(m_verbosity > 2){
    pout() << "amr_mesh::define_neighbors" << endl;
  }

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

void amr_mesh::define_eblevelgrid(const int a_lmin){
  CH_TIME("amr_mesh::define_eblevelgrid");
  if(m_verbosity > 2){
    pout() << "amr_mesh::define_eblevelgrid" << endl;
  }

  const RefCountedPtr<EBIndexSpace>& ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace>& ebis_sol = m_mfis->get_ebis(phase::solid);

  for (int lvl = a_lmin; lvl <= m_finest_level; lvl++){
    if(!ebis_gas.isNull()){
      m_eblg[phase::gas][lvl]  = RefCountedPtr<EBLevelGrid> (new EBLevelGrid(m_grids[lvl],
									     m_domains[lvl],
									     m_ebghost,
									     ebis_gas));
      if(lvl > 0){
	m_eblg[phase::gas][lvl]->setMaxCoarseningRatio(m_ref_ratios[lvl-1], ebis_gas);
      }
      m_ebisl[phase::gas][lvl] = m_eblg[phase::gas][lvl]->getEBISL();
    }
    if(!ebis_sol.isNull()){
      m_eblg[phase::solid][lvl]  = RefCountedPtr<EBLevelGrid> (new EBLevelGrid(m_grids[lvl],
									       m_domains[lvl],
									       m_ebghost,
									       ebis_sol));
      if(lvl > 0){
	m_eblg[phase::solid][lvl]->setMaxCoarseningRatio(m_ref_ratios[lvl-1], ebis_sol);
      }
      m_ebisl[phase::solid][lvl] = m_eblg[phase::solid][lvl]->getEBISL();
    }
  }
}

void amr_mesh::define_vofiter(const int a_lmin){
  CH_TIME("amr_mesh::define_vofiter");
  if(m_verbosity > 2){
    pout() << "amr_mesh::define_vofiter" << endl;
  }

  const RefCountedPtr<EBIndexSpace>& ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace>& ebis_sol = m_mfis->get_ebis(phase::solid);

  for (int lvl = a_lmin; lvl <= m_finest_level; lvl++){
    if(!ebis_gas.isNull()){
      m_vofiter[phase::gas][lvl] = RefCountedPtr<LayoutData<VoFIterator> > (new LayoutData<VoFIterator> (m_grids[lvl]));

      for (DataIterator dit = m_grids[lvl].dataIterator(); dit.ok(); ++dit){
	VoFIterator& vofit = (*m_vofiter[phase::gas][lvl])[dit()];

	const Box& box         = m_grids[lvl].get(dit());
	const EBISBox& ebisbox = m_ebisl[phase::gas][lvl][dit()];
	const EBGraph& ebgraph = ebisbox.getEBGraph();
	const IntVectSet& irreg = ebisbox.getIrregIVS(box);

	vofit.define(irreg, ebgraph);
      }
    }
    if(!ebis_sol.isNull()){
      m_vofiter[phase::solid][lvl] = RefCountedPtr<LayoutData<VoFIterator> > (new LayoutData<VoFIterator> (m_grids[lvl]));

      for (DataIterator dit = m_grids[lvl].dataIterator(); dit.ok(); ++dit){
	VoFIterator& vofit = (*m_vofiter[phase::solid][lvl])[dit()];

	const Box& box         = m_grids[lvl].get(dit());
	const EBISBox& ebisbox = m_ebisl[phase::solid][lvl][dit()];
	const EBGraph& ebgraph = ebisbox.getEBGraph();
	const IntVectSet& irreg = ebisbox.getIrregIVS(box);

	vofit.define(irreg, ebgraph);
      }
    }
  }
}

void amr_mesh::define_mflevelgrid(const int a_lmin){
  CH_TIME("amr_mesh::define_mflevelgrid");
  if(m_verbosity > 2){
    pout() << "amr_mesh::define_mflevelgrid" << endl;
  }

  const RefCountedPtr<EBIndexSpace>& ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace>& ebis_sol = m_mfis->get_ebis(phase::solid);

  for (int lvl = a_lmin; lvl <= m_finest_level; lvl++){
    Vector<EBLevelGrid> eblgs;
    if(!ebis_gas.isNull()){
      eblgs.push_back(*m_eblg[phase::gas][lvl]);
    }
    if(!ebis_sol.isNull()){
      eblgs.push_back(*m_eblg[phase::solid][lvl]);
    }

    m_mflg[lvl] = RefCountedPtr<MFLevelGrid> (new MFLevelGrid(m_mfis, eblgs));
  }
}

void amr_mesh::define_eb_coar_ave(const int a_lmin){
  CH_TIME("amr_mesh::define_eb_coar_ave");
  if(m_verbosity > 2){
    pout() << "amr_mesh::define_eb_coar_ave" << endl;
  }

  const int comps = SpaceDim;

  const RefCountedPtr<EBIndexSpace>& ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace>& ebis_sol = m_mfis->get_ebis(phase::solid);

  for (int lvl = a_lmin; lvl <= m_finest_level; lvl++){

    const bool has_coar = lvl > 0;

    if(has_coar){
      if(!ebis_gas.isNull()){
	if(this->query_operator(s_eb_coar_ave, phase::gas)){
	  m_coarave[phase::gas][lvl] = RefCountedPtr<ebcoarseaverage> (new ebcoarseaverage(m_grids[lvl],
											   m_grids[lvl-1],
											   m_ebisl[phase::gas][lvl],
											   m_ebisl[phase::gas][lvl-1],
											   m_domains[lvl-1],
											   m_ref_ratios[lvl-1],
											   comps,
											   &(*ebis_gas)));
	}
      }

      if(!ebis_sol.isNull()){
	if(this->query_operator(s_eb_coar_ave, phase::solid)){
	  m_coarave[phase::solid][lvl] = RefCountedPtr<ebcoarseaverage> (new ebcoarseaverage(m_grids[lvl],
											     m_grids[lvl-1],
											     m_ebisl[phase::solid][lvl],
											     m_ebisl[phase::solid][lvl-1],
											     m_domains[lvl-1],
											     m_ref_ratios[lvl-1],
											     comps,
											     ebis_sol));
	}
      }
    }
  }
}

void amr_mesh::define_eb_quad_cfi(const int a_lmin){
  CH_TIME("amr_mesh::define_eb_quad_cfi");
  if(m_verbosity > 2){
    pout() << "amr_mesh::define_eb_quad_cfi" << endl;
  }


    
  const int comps = SpaceDim;

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_mfis->get_ebis(phase::solid);

  for (int lvl = a_lmin; lvl <= m_finest_level; lvl++){

    const bool has_coar = lvl > 0;

    if(has_coar){
      if(!ebis_gas.isNull()){
	if(this->query_operator(s_eb_quad_cfi, phase::gas)){
	  const LayoutData<IntVectSet>& cfivs = *(m_eblg[phase::gas][lvl]->getCFIVS());

	    m_quadcfi[phase::gas][lvl] = RefCountedPtr<nwoebquadcfinterp> (new nwoebquadcfinterp(m_grids[lvl],
												 m_grids[lvl-1],
												 m_ebisl[phase::gas][lvl],
												 m_ebisl[phase::gas][lvl-1],
												 m_domains[lvl-1],
												 m_ref_ratios[lvl-1],
												 comps,
												 m_dx[lvl],
												 m_num_ghost,
												 cfivs,
												 ebis_gas));

	  m_old_quadcfi[phase::gas][lvl] = RefCountedPtr<EBQuadCFInterp> (new EBQuadCFInterp(m_grids[lvl],
											     m_grids[lvl-1],
											     m_ebisl[phase::gas][lvl],
											     m_ebisl[phase::gas][lvl-1],
											     m_domains[lvl-1],
											     m_ref_ratios[lvl-1],
											     1,
											     cfivs,
											     ebis_gas));
	}
      }
      if(!ebis_sol.isNull()){
	if(this->query_operator(s_eb_quad_cfi, phase::solid)){
	  LayoutData<IntVectSet>& cfivs = *(m_eblg[phase::solid][lvl]->getCFIVS());
	    m_quadcfi[phase::solid][lvl] = RefCountedPtr<nwoebquadcfinterp> (new nwoebquadcfinterp(m_grids[lvl],
												   m_grids[lvl-1],
												   m_ebisl[phase::solid][lvl],
												   m_ebisl[phase::solid][lvl-1],
												   m_domains[lvl-1],
												   m_ref_ratios[lvl-1],
												   comps,
												   m_dx[lvl],
												   m_num_ghost,
												   cfivs,
												   ebis_sol));
	
	  m_old_quadcfi[phase::solid][lvl] = RefCountedPtr<EBQuadCFInterp> (new EBQuadCFInterp(m_grids[lvl],
											       m_grids[lvl-1],
											       m_ebisl[phase::solid][lvl],
											       m_ebisl[phase::solid][lvl-1],
											       m_domains[lvl-1],
											       m_ref_ratios[lvl-1],
											       1,
											       cfivs,
											       ebis_sol));
	}
      }
    }
  }
}

void amr_mesh::define_fillpatch(const int a_lmin){
  CH_TIME("amr_mesh::define_fillpatch");
  if(m_verbosity > 2){
    pout() << "amr_mesh::define_fillpatch" << endl;
  }


  const int comps     = SpaceDim;

  // Should these be input somehow?
  const int radius    = 1;
  const IntVect ghost = m_num_ghost*IntVect::Unit;

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_mfis->get_ebis(phase::solid);

  for (int lvl = a_lmin; lvl <= m_finest_level; lvl++){

    const bool has_coar = lvl > 0;

    if(has_coar){
      if(!ebis_gas.isNull()){
	if(this->query_operator(s_eb_fill_patch, phase::gas)){
	  const LayoutData<IntVectSet>& cfivs = *(m_eblg[phase::gas][lvl]->getCFIVS());
	  m_pwl_fillpatch[phase::gas][lvl] = RefCountedPtr<AggEBPWLFillPatch> (new AggEBPWLFillPatch(m_grids[lvl],
												     m_grids[lvl-1],
												     m_ebisl[phase::gas][lvl],
												     m_ebisl[phase::gas][lvl-1],
												     m_domains[lvl-1],
												     m_ref_ratios[lvl-1],
												     comps,
												     radius,
												     ghost,
												     !m_ebcf,
												     ebis_gas));
	}
      }
      if(!ebis_sol.isNull()){
	if(this->query_operator(s_eb_fill_patch, phase::solid)){
	  const LayoutData<IntVectSet>& cfivs = *(m_eblg[phase::solid][lvl]->getCFIVS());
	  m_pwl_fillpatch[phase::solid][lvl] = RefCountedPtr<AggEBPWLFillPatch>
	    (new AggEBPWLFillPatch(m_grids[lvl],
				   m_grids[lvl-1],
				   m_ebisl[phase::solid][lvl],
				   m_ebisl[phase::solid][lvl-1],
				   m_domains[lvl-1],
				   m_ref_ratios[lvl-1],
				   comps,
				   radius,
				   ghost,
				   !m_ebcf,
				   ebis_sol));
	}
      }
    }
  }
}

void amr_mesh::define_ebpwl_interp(const int a_lmin){
  CH_TIME("amr_mesh::define_ebpwl_interp");
  if(m_verbosity > 2){
    pout() << "amr_mesh::define_ebpwl_interp" << endl;
  }

  const int comps     = SpaceDim;

  // Should these be input somehow?
  const int radius    = 1;
  const IntVect ghost = m_num_ghost*IntVect::Unit;

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_mfis->get_ebis(phase::solid);

  for (int lvl = a_lmin; lvl <= m_finest_level; lvl++){

    const bool has_coar = lvl > 0;

    if(has_coar){
      if(!ebis_gas.isNull()){
	if(this->query_operator(s_eb_pwl_interp, phase::gas)){
	  m_pwl_interp[phase::gas][lvl] = RefCountedPtr<EBPWLFineInterp> (new EBPWLFineInterp(m_grids[lvl],
											      m_grids[lvl-1],
											      m_ebisl[phase::gas][lvl],
											      m_ebisl[phase::gas][lvl-1],
											      m_domains[lvl-1],
											      m_ref_ratios[lvl-1],
											      comps,
											      ebis_gas));
	}
      }
      if(!ebis_sol.isNull()){
	if(this->query_operator(s_eb_pwl_interp, phase::solid)){
	  m_pwl_interp[phase::solid][lvl] = RefCountedPtr<EBPWLFineInterp> (new EBPWLFineInterp(m_grids[lvl],
												m_grids[lvl-1],
												m_ebisl[phase::solid][lvl],
												m_ebisl[phase::solid][lvl-1],
												m_domains[lvl-1],
												m_ref_ratios[lvl-1],
												comps,
												ebis_sol));
	}
      }
    }
  }
}

void amr_mesh::define_ebmg_interp(const int a_lmin){
  CH_TIME("amr_mesh::define_ebmg_interp");
  if(m_verbosity > 2){
    pout() << "amr_mesh::define_ebmg_interp" << endl;
  }

  const int ncomps    = 1;

  // Should these be input somehow?
  const int radius    = 1;
  const IntVect ghost = m_num_ghost*IntVect::Unit;

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_mfis->get_ebis(phase::solid);

  for (int lvl = a_lmin; lvl <= m_finest_level; lvl++){

    const bool has_coar = lvl > 0;

    if(has_coar){
      if(!ebis_gas.isNull()){
	if(this->query_operator(s_eb_mg_interp, phase::gas)){
	  m_ebmg_interp[phase::gas][lvl] = RefCountedPtr<EBMGInterp> (new EBMGInterp(m_grids[lvl],
										     m_grids[lvl-1],
										     m_ebisl[phase::gas][lvl],
										     m_ebisl[phase::gas][lvl-1],
										     m_domains[lvl-1],
										     m_ref_ratios[lvl-1],
										     SpaceDim,
										     ebis_gas,
										     m_num_ghost*IntVect::Unit));
	}
      }
      if(!ebis_sol.isNull()){
	if(this->query_operator(s_eb_pwl_interp, phase::solid)){
	  m_ebmg_interp[phase::solid][lvl] = RefCountedPtr<EBMGInterp> (new EBMGInterp(m_grids[lvl],
										       m_grids[lvl-1],
										       m_ebisl[phase::solid][lvl],
										       m_ebisl[phase::solid][lvl-1],
										       m_domains[lvl-1],
										       m_ref_ratios[lvl-1],
										       SpaceDim,
										       ebis_sol,
										       m_num_ghost*IntVect::Unit));
	}
      }
    }
  }
}

void amr_mesh::define_flux_reg(const int a_lmin, const int a_regsize){
  CH_TIME("amr_mesh::define_flux_reg");
  if(m_verbosity > 2){
    pout() << "amr_mesh::define_flux_reg" << endl;
  }

  // TLDR; The flux register between levels [lvl,(lvl+1)] lives on vector entry [lvl]. Since a_lmin is the coarsest
  //       level which has changed, we need to update the flux register in entry a_lmin-1 as well. 

  const int comps = a_regsize;

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_mfis->get_ebis(phase::solid);

  for (int lvl = Max(0,a_lmin-1); lvl <= m_finest_level; lvl++){

    const bool has_fine = lvl < m_finest_level;

    if(!ebis_gas.isNull()){
      if(has_fine){
	if(this->query_operator(s_eb_flux_reg, phase::gas)){
	  m_flux_reg[phase::gas][lvl] = RefCountedPtr<EBFluxRegister> (new EBFluxRegister(m_grids[lvl+1],
											  m_grids[lvl],
											  m_ebisl[phase::gas][lvl+1],
											  m_ebisl[phase::gas][lvl],
											  m_domains[lvl].domainBox(),
											  m_ref_ratios[lvl],
											  comps,
											  ebis_gas));
	}
      }
    }
    if(!ebis_sol.isNull()){
      if(has_fine){
	if(this->query_operator(s_eb_flux_reg, phase::solid)){
	  m_flux_reg[phase::solid][lvl] = RefCountedPtr<EBFluxRegister> (new EBFluxRegister(m_grids[lvl+1],
											    m_grids[lvl],
											    m_ebisl[phase::solid][lvl+1],
											    m_ebisl[phase::solid][lvl],
											    m_domains[lvl].domainBox(),
											    m_ref_ratios[lvl],
											    comps,
											    ebis_sol));
	}
      }
    }
  }
}

void amr_mesh::define_redist_oper(const int a_lmin, const int a_regsize){
  CH_TIME("amr_mesh::define_redist_oper");
  if(m_verbosity > 2){
    pout() << "amr_mesh::define_redist_oper" << endl;
  }

  // TLDR: All these operators either do stuff on an AMR level, or between a coarse and a fine level. The entries
  //       live on these levels:
  //
  //       Oper                        Level
  //       EBFineToCoar [l,  l-1]    l
  //       EBCoarToFine [l-1,l  ] 
  //       when level a_lmin changed we need to update fine-to-coar 

  const int comps = a_regsize;

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_mfis->get_ebis(phase::solid);

  Real t_level = 0.0;
  Real t_coar2fine = 0.0;
  Real t_coar2coar = 0.0;
  Real t_fine2coar = 0.0;
  for (int lvl = Max(0, a_lmin-1); lvl <= m_finest_level; lvl++){

    const bool has_coar = lvl > 0;
    const bool has_fine = lvl < m_finest_level;

    if(!ebis_gas.isNull()){
      if(this->query_operator(s_eb_redist, phase::gas)){
	      
	if(lvl >= a_lmin){
	  t_level -= MPI_Wtime();
	  m_level_redist[phase::gas][lvl] = RefCountedPtr<EBLevelRedist> (new EBLevelRedist(m_grids[lvl],
											    m_ebisl[phase::gas][lvl],
											    m_domains[lvl],
											    comps,
											    m_redist_rad));

	  m_concentration_redist[phase::gas][lvl] = RefCountedPtr<EBLevelConcentrationRedist>
	    (new EBLevelConcentrationRedist(m_grids[lvl],
					    m_ebisl[phase::gas][lvl],
					    m_domains[lvl],
					    comps,
					    m_redist_rad));


	  t_level += MPI_Wtime();
	}

    
	if(m_ebcf){
	  if(has_coar){

	    // TLDR: The fine-to-coar redistribution operator that transfers from the fine level to the coar level
	    //       obviously lives on the fine level. But since a_lmin is the coarsest level that changed, we only
	    //       need to update this if lvl >= a_lmin
	    if(lvl >= a_lmin){
	      t_fine2coar -= MPI_Wtime();
#if USE_NEW_F2C_REDIST
	      auto redist = RefCountedPtr<EBFastFineToCoarRedist> (new EBFastFineToCoarRedist());
	      redist->define(*m_eblg[phase::gas][lvl],
			     *m_eblg[phase::gas][lvl-1],
			     *m_neighbors[lvl],
			     *m_neighbors[lvl-1],
			     m_ref_ratios[lvl-1],
			     comps,
			     m_redist_rad);
	      m_fine_to_coar_redist[phase::gas][lvl] = RefCountedPtr<EBFineToCoarRedist> (redist);
#else
	    
	      m_fine_to_coar_redist[phase::gas][lvl] = RefCountedPtr<EBFineToCoarRedist> (new EBFineToCoarRedist());
	      m_fine_to_coar_redist[phase::gas][lvl]->define(m_grids[lvl],
							     m_grids[lvl-1],
							     m_ebisl[phase::gas][lvl],
							     m_ebisl[phase::gas][lvl-1],
							     m_domains[lvl-1].domainBox(),
							     m_ref_ratios[lvl-1],
							     comps,
							     m_redist_rad,
							     ebis_gas);
#endif
	      t_fine2coar += MPI_Wtime();
	    }
	  }

	  if(has_fine){
	    // TLDR: The coar-to-fine redistribution operator transfers from the coarse level and to the fine level and
	    //       therefore lives on the coarse level. Since a_lmin is the coarsest level that changed, we need to update
	    //       if lvl >= a_lmin-1
	    if(lvl >= a_lmin-1){
	      t_coar2fine -= MPI_Wtime();
#if USE_NEW_C2F_REDIST
	      auto c2f_redist = RefCountedPtr<EBFastCoarToFineRedist> (new EBFastCoarToFineRedist());
	      c2f_redist->define(*m_eblg[phase::gas][lvl+1],
				 *m_eblg[phase::gas][lvl],
				 *m_neighbors[lvl+1],
				 *m_neighbors[lvl],
				 m_ref_ratios[lvl],
				 comps,
				 m_redist_rad);
	      m_coar_to_fine_redist[phase::gas][lvl] = RefCountedPtr<EBCoarToFineRedist> (c2f_redist);
#else
	      m_coar_to_fine_redist[phase::gas][lvl] = RefCountedPtr<EBCoarToFineRedist> (new EBCoarToFineRedist());
	      m_coar_to_fine_redist[phase::gas][lvl]->define(m_grids[lvl+1],
							     m_grids[lvl],
							     m_ebisl[phase::gas][lvl],
							     m_domains[lvl].domainBox(),
							     m_ref_ratios[lvl],
							     comps,
							     m_redist_rad,
							     ebis_gas);
#endif
	      t_coar2fine += MPI_Wtime();

	      // Coarse to coarse redistribution
	      t_coar2coar -= MPI_Wtime();
#if USE_NEW_C2C_REDIST
	      auto c2c_redist = RefCountedPtr<EBFastCoarToCoarRedist> (new EBFastCoarToCoarRedist());
	      c2c_redist->define(*m_eblg[phase::gas][lvl+1],
				 *m_eblg[phase::gas][lvl],
				 *m_neighbors[lvl+1],
				 *m_neighbors[lvl],
				 m_ref_ratios[lvl],
				 comps,
				 m_redist_rad);
	      m_coar_to_coar_redist[phase::gas][lvl] = RefCountedPtr<EBCoarToCoarRedist> (c2c_redist);
#else
	      m_coar_to_coar_redist[phase::gas][lvl] = RefCountedPtr<EBCoarToCoarRedist> (new EBCoarToCoarRedist());
	      m_coar_to_coar_redist[phase::gas][lvl]->define(*m_eblg[phase::gas][lvl+1],
							     *m_eblg[phase::gas][lvl],
							     m_ref_ratios[lvl],
							     comps,
							     m_redist_rad);
#endif
	      t_coar2coar += MPI_Wtime();
	    }
	  }
	}
      }
    }

    if(!ebis_sol.isNull()){
      if(this->query_operator(s_eb_redist, phase::solid)){
	      
	if(lvl >= a_lmin){
	  t_level -= MPI_Wtime();
	  m_level_redist[phase::solid][lvl] = RefCountedPtr<EBLevelRedist> (new EBLevelRedist(m_grids[lvl],
											    m_ebisl[phase::solid][lvl],
											    m_domains[lvl],
											    comps,
											    m_redist_rad));

	  m_concentration_redist[phase::solid][lvl] = RefCountedPtr<EBLevelConcentrationRedist>
	    (new EBLevelConcentrationRedist(m_grids[lvl],
					    m_ebisl[phase::solid][lvl],
					    m_domains[lvl],
					    comps,
					    m_redist_rad));


	  t_level += MPI_Wtime();
	}

    
	if(m_ebcf){
	  if(has_coar){

	    // TLDR: The fine-to-coar redistribution operator that transfers from the fine level to the coar level
	    //       obviously lives on the fine level. But since a_lmin is the coarsest level that changed, we only
	    //       need to update this if lvl >= a_lmin
	    if(lvl >= a_lmin){
	      t_fine2coar -= MPI_Wtime();
#if USE_NEW_F2C_REDIST
	      auto redist = RefCountedPtr<EBFastFineToCoarRedist> (new EBFastFineToCoarRedist());
	      redist->define(*m_eblg[phase::solid][lvl],
			     *m_eblg[phase::solid][lvl-1],
			     *m_neighbors[lvl],
			     *m_neighbors[lvl-1],
			     m_ref_ratios[lvl-1],
			     comps,
			     m_redist_rad);
	      m_fine_to_coar_redist[phase::solid][lvl] = RefCountedPtr<EBFineToCoarRedist> (redist);
#else
	    
	      m_fine_to_coar_redist[phase::solid][lvl] = RefCountedPtr<EBFineToCoarRedist> (new EBFineToCoarRedist());
	      m_fine_to_coar_redist[phase::solid][lvl]->define(m_grids[lvl],
							     m_grids[lvl-1],
							     m_ebisl[phase::solid][lvl],
							     m_ebisl[phase::solid][lvl-1],
							     m_domains[lvl-1].domainBox(),
							     m_ref_ratios[lvl-1],
							     comps,
							     m_redist_rad,
							     ebis_gas);
#endif
	      t_fine2coar += MPI_Wtime();
	    }
	  }

	  if(has_fine){
	    // TLDR: The coar-to-fine redistribution operator transfers from the coarse level and to the fine level and
	    //       therefore lives on the coarse level. Since a_lmin is the coarsest level that changed, we need to update
	    //       if lvl >= a_lmin-1
	    if(lvl >= a_lmin-1){
	      t_coar2fine -= MPI_Wtime();
#if USE_NEW_C2F_REDIST
	      auto c2f_redist = RefCountedPtr<EBFastCoarToFineRedist> (new EBFastCoarToFineRedist());
	      c2f_redist->define(*m_eblg[phase::solid][lvl+1],
				 *m_eblg[phase::solid][lvl],
				 *m_neighbors[lvl+1],
				 *m_neighbors[lvl],
				 m_ref_ratios[lvl],
				 comps,
				 m_redist_rad);
	      m_coar_to_fine_redist[phase::solid][lvl] = RefCountedPtr<EBCoarToFineRedist> (c2f_redist);
#else
	      m_coar_to_fine_redist[phase::solid][lvl] = RefCountedPtr<EBCoarToFineRedist> (new EBCoarToFineRedist());
	      m_coar_to_fine_redist[phase::solid][lvl]->define(m_grids[lvl+1],
							     m_grids[lvl],
							     m_ebisl[phase::solid][lvl],
							     m_domains[lvl].domainBox(),
							     m_ref_ratios[lvl],
							     comps,
							     m_redist_rad,
							     ebis_gas);
#endif
	      t_coar2fine += MPI_Wtime();

	      // Coarse to coarse redistribution
	      t_coar2coar -= MPI_Wtime();
#if USE_NEW_C2C_REDIST
	      auto c2c_redist = RefCountedPtr<EBFastCoarToCoarRedist> (new EBFastCoarToCoarRedist());
	      c2c_redist->define(*m_eblg[phase::solid][lvl+1],
				 *m_eblg[phase::solid][lvl],
				 *m_neighbors[lvl+1],
				 *m_neighbors[lvl],
				 m_ref_ratios[lvl],
				 comps,
				 m_redist_rad);
	      m_coar_to_coar_redist[phase::solid][lvl] = RefCountedPtr<EBCoarToCoarRedist> (c2c_redist);
#else
	      m_coar_to_coar_redist[phase::solid][lvl] = RefCountedPtr<EBCoarToCoarRedist> (new EBCoarToCoarRedist());
	      m_coar_to_coar_redist[phase::solid][lvl]->define(*m_eblg[phase::solid][lvl+1],
							     *m_eblg[phase::solid][lvl],
							     m_ref_ratios[lvl],
							     comps,
							     m_redist_rad);
#endif
	      t_coar2coar += MPI_Wtime();
	    }
	  }
	}
      }
    }
  }

#if AMR_MESH_DEBUG
  pout() << "amr_mesh::define_redist_oper breakdown by operator: " << endl
	 << "t_level = " << t_level << endl
    	 << "t_fine2coar = " << t_fine2coar << endl
    	 << "t_coar2fine = " << t_coar2fine << endl
    	 << "t_coar2coar = " << t_coar2coar << endl;
#endif
}

void amr_mesh::define_irreg_sten(){
  CH_TIME("amr_mesh::define_irreg_sten");
  if(m_verbosity > 2){
    pout() << "amr_mesh::define_irreg_sten" << endl;
  }

  const int order = m_irreg_sten_order;
  const int rad   = m_irreg_sten_radius;

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_mfis->get_ebis(phase::solid);

  if(!ebis_gas.isNull()){
    if(this->query_operator(s_eb_irreg_interp, phase::gas)){
      m_centroid_interp[phase::gas] = RefCountedPtr<irreg_amr_stencil<centroid_interp> >
	(new irreg_amr_stencil<centroid_interp>(m_grids,
						m_ebisl[phase::gas],
						m_domains,
						m_dx,
						m_finest_level,
						m_centroid_sten_order,
						m_centroid_sten_rad,
						m_centroid_stencil));

      m_eb_centroid_interp[phase::gas] = RefCountedPtr<irreg_amr_stencil<eb_centroid_interp> >
	(new irreg_amr_stencil<eb_centroid_interp>(m_grids,
						   m_ebisl[phase::gas],
						   m_domains,
						   m_dx,
						   m_finest_level,
						   m_eb_sten_order,
						   m_eb_sten_rad,
						   m_eb_stencil));
    }
  }

  if(!ebis_sol.isNull()){
    if(this->query_operator(s_eb_irreg_interp, phase::solid)){
      m_centroid_interp[phase::solid] = RefCountedPtr<irreg_amr_stencil<centroid_interp> >
	(new irreg_amr_stencil<centroid_interp>(m_grids,
						m_ebisl[phase::solid],
						m_domains,
						m_dx,
						m_finest_level,
						m_centroid_sten_order,
						m_centroid_sten_rad,
						m_centroid_stencil));
      m_eb_centroid_interp[phase::solid] = RefCountedPtr<irreg_amr_stencil<eb_centroid_interp> >
	(new irreg_amr_stencil<eb_centroid_interp>(m_grids,
						   m_ebisl[phase::solid],
						   m_domains,
						   m_dx,
						   m_finest_level,
						   m_eb_sten_order,
						   m_eb_sten_rad,
						   m_eb_stencil));
    }
  }
}

void amr_mesh::define_noncons_sten(){
  CH_TIME("amr_mesh::define_noncons_sten");
  if(m_verbosity > 2){
    pout() << "amr_mesh::define_noncons_sten" << endl;
  }

  const int order = 1; // Dummy argument
  const int rad   = m_redist_rad;

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_mfis->get_ebis(phase::solid);

  if(!ebis_gas.isNull()){
    if(this->query_operator(s_eb_noncons_div, phase::gas)){
      m_noncons_div[phase::gas] = RefCountedPtr<irreg_amr_stencil<noncons_div> >
	(new irreg_amr_stencil<noncons_div>(m_grids,
					    m_ebisl[phase::gas],
					    m_domains,
					    m_dx,
					    m_finest_level,
					    order,                 // Dummy argument
					    m_redist_rad,
					    m_centroid_stencil));  // Dummy argumement
    }
  }

  if(!ebis_gas.isNull()){
    if(this->query_operator(s_eb_noncons_div, phase::solid)){
      m_noncons_div[phase::solid] = RefCountedPtr<irreg_amr_stencil<noncons_div> >
	(new irreg_amr_stencil<noncons_div>(m_grids,
					    m_ebisl[phase::solid],
					    m_domains,
					    m_dx,
					    m_finest_level,
					    order,                 // Dummy argument
					    m_redist_rad,
					    m_centroid_stencil));  // Dummy argumement
    }
  }


}

void amr_mesh::define_copier(const int a_lmin){
  CH_TIME("amr_mesh::define_copier");
  if(m_verbosity > 3){
    pout() << "amr_mesh::define_copier" << endl;
  }

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_mfis->get_ebis(phase::solid);

  for (int lvl = a_lmin; lvl <= m_finest_level; lvl++){
    if(!ebis_gas.isNull()){
      if(this->query_operator(s_eb_copier, phase::gas)){
	m_copier[phase::gas][lvl] = RefCountedPtr<Copier> (new Copier(m_grids[lvl],
								      m_grids[lvl],
								      m_domains[lvl],
								      m_num_ghost*IntVect::Unit,
								      true));
	m_reverse_copier[phase::gas][lvl] = RefCountedPtr<Copier> (new Copier(m_grids[lvl],
									      m_grids[lvl],
									      m_domains[lvl],
									      m_num_ghost*IntVect::Unit,
									      true));
	m_reverse_copier[phase::gas][lvl]->reverse();
      }
    }

    if(!ebis_sol.isNull()){
      if(this->query_operator(s_eb_copier, phase::solid)){
	m_copier[phase::solid][lvl] = RefCountedPtr<Copier> (new Copier(m_grids[lvl],
									m_grids[lvl],
									m_domains[lvl],
									m_num_ghost*IntVect::Unit,
									true));
	m_reverse_copier[phase::solid][lvl] = RefCountedPtr<Copier> (new Copier(m_grids[lvl],
										m_grids[lvl],
										m_domains[lvl],
										m_num_ghost*IntVect::Unit,
										true));
	m_reverse_copier[phase::solid][lvl]->reverse();
      }
    }
  }
}

void amr_mesh::define_ghostcloud(const int a_lmin){
  CH_TIME("amr_mesh::define_ghostcloud");
  if(m_verbosity > 3){
    pout() << "amr_mesh::define_ghostcloud" << endl;
  }

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_mfis->get_ebis(phase::solid);

  for (int lvl = a_lmin; lvl <= m_finest_level; lvl++){
    const bool has_coar = lvl > 0;

    if(has_coar){
      if(!ebis_gas.isNull()){
	if(this->query_operator(s_eb_copier, phase::gas)){
	  m_ghostclouds[phase::gas][lvl] = RefCountedPtr<EBGhostCloud> (new EBGhostCloud(m_grids[lvl-1],
											 m_grids[lvl],
											 *(m_eblg[phase::gas][lvl-1]),
											 *(m_eblg[phase::gas][lvl]),
											 m_domains[lvl-1],
											 m_domains[lvl],
											 m_ref_ratios[lvl-1],
											 1,
											 m_num_ghost));
	}
      }
      if(!ebis_sol.isNull()){
	if(this->query_operator(s_eb_copier, phase::solid)){
	  m_ghostclouds[phase::solid][lvl] = RefCountedPtr<EBGhostCloud> (new EBGhostCloud(m_grids[lvl-1],
											   m_grids[lvl],
											   *(m_eblg[phase::solid][lvl-1]),
											   *(m_eblg[phase::solid][lvl]),
											   m_domains[lvl-1],
											   m_domains[lvl],
											   m_ref_ratios[lvl-1],
											   1,
											   m_num_ghost));
	}
      }
    }
  }
}

void amr_mesh::average_down(EBAMRCellData& a_data, phase::which_phase a_phase){
  CH_TIME("amr_mesh::average_down");
  if(m_verbosity > 3){
    pout() << "amr_mesh::average_down(ebcell)" << endl;
  }
  const RefCountedPtr<EBIndexSpace>& ebis = m_mfis->get_ebis(a_phase);

  if(!ebis.isNull()){
    for (int lvl = m_finest_level; lvl > 0; lvl--){
      const int ncomps = a_data[lvl]->nComp();
      const Interval interv (0, ncomps-1);

      m_coarave[a_phase][lvl]->average(*a_data[lvl-1], *a_data[lvl], interv);
    }
    for (int lvl = 0; lvl <= m_finest_level; lvl++){
      a_data[lvl]->exchange();
    }
  }
}

void amr_mesh::average_down(EBAMRCellData& a_data, phase::which_phase a_phase, const int a_lvl){
  CH_TIME("amr_mesh::average_down");
  if(m_verbosity > 3){
    pout() << "amr_mesh::average_down(ebcell, level)" << endl;
  }
  
  const RefCountedPtr<EBIndexSpace>& ebis = m_mfis->get_ebis(a_phase);

  if(!ebis.isNull()){
      const int ncomps = a_data[a_lvl]->nComp();
      const Interval interv (0, ncomps-1);

      m_coarave[a_phase][a_lvl+1]->average(*a_data[a_lvl], *a_data[a_lvl+1], interv);
  }

  a_data[a_lvl]->exchange();
}

void amr_mesh::average_down(MFAMRFluxData& a_data){
  CH_TIME("amr_mesh::average_down");
  if(m_verbosity > 3){
    pout() << "amr_mesh::average_down(mfflux)" << endl;
  }

  EBAMRFluxData alias_g(1 + m_finest_level);
  EBAMRFluxData alias_s(1 + m_finest_level);

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_mfis->get_ebis(phase::solid);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    alias_g[lvl] = RefCountedPtr<LevelData<EBFluxFAB> > (new LevelData<EBFluxFAB>());
    alias_s[lvl] = RefCountedPtr<LevelData<EBFluxFAB> > (new LevelData<EBFluxFAB>());
      
    mfalias::aliasMF(*alias_g[lvl], phase::gas,   *a_data[lvl]);
    if(!ebis_sol.isNull()){
      mfalias::aliasMF(*alias_s[lvl], phase::solid, *a_data[lvl]);
    }
  }

  this->average_down(alias_g, phase::gas);
  if(!ebis_sol.isNull()){
    this->average_down(alias_s, phase::solid);
  }
}

void amr_mesh::average_down(MFAMRCellData& a_data){
  CH_TIME("amr_mesh::average_down");
  if(m_verbosity > 3){
    pout() << "amr_mesh::average_down(mfcell)" << endl;
  }
  
  EBAMRCellData alias_g(1 + m_finest_level);
  EBAMRCellData alias_s(1 + m_finest_level);

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_mfis->get_ebis(phase::solid);

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    alias_g[lvl] = RefCountedPtr<LevelData<EBCellFAB> > (new LevelData<EBCellFAB>());
    alias_s[lvl] = RefCountedPtr<LevelData<EBCellFAB> > (new LevelData<EBCellFAB>());
      
    mfalias::aliasMF(*alias_g[lvl], phase::gas,   *a_data[lvl]);
    if(!ebis_sol.isNull()){
      mfalias::aliasMF(*alias_s[lvl], phase::solid, *a_data[lvl]);
    }
  }

  this->average_down(alias_g, phase::gas);
  if(!ebis_sol.isNull()){
    this->average_down(alias_s, phase::solid);
  }
}

void amr_mesh::average_down(EBAMRFluxData& a_data, phase::which_phase a_phase){
  CH_TIME("amr_mesh::average_down");
  if(m_verbosity > 3){
    pout() << "amr_mesh::average_down(face)" << endl;
  }

  for (int lvl = m_finest_level; lvl > 0; lvl--){
    const int ncomps = a_data[lvl]->nComp();
    const Interval interv (0, ncomps-1);

    m_coarave[a_phase][lvl]->average(*a_data[lvl-1], *a_data[lvl], interv);

  }

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    a_data[lvl]->exchange();
  }
}

void amr_mesh::average_down(EBAMRIVData& a_data, phase::which_phase a_phase){
  CH_TIME("amr_mesh::average_down");
  if(m_verbosity > 3){
    pout() << "amr_mesh::average_down(iv)" << endl;
  }

  for (int lvl = m_finest_level; lvl > 0; lvl--){
    const int ncomps = a_data[lvl]->nComp();
    const Interval interv (0, ncomps-1);

    m_coarave[a_phase][lvl]->average(*a_data[lvl-1], *a_data[lvl], interv);

  }

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    a_data[lvl]->exchange();
  }

}

void amr_mesh::conservative_average(EBAMRIVData& a_data, phase::which_phase a_phase){
  CH_TIME("amr_mesh::conservative_average");
  if(m_verbosity > 3){
    pout() << "amr_mesh::conservative_average" << endl;
  }

  for (int lvl = m_finest_level; lvl > 0; lvl--){
    const int ncomps = a_data[lvl]->nComp();
    const Interval interv (0, ncomps-1);

    m_coarave[a_phase][lvl]->conservative_average(*a_data[lvl-1], *a_data[lvl], interv);

  }

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    a_data[lvl]->exchange();
  }

}

void amr_mesh::interp_ghost(EBAMRCellData& a_data, phase::which_phase a_phase){
  CH_TIME("amr_mesh::interp_ghost(eb)");
  if(m_verbosity > 3){
    pout() << "amr_mesh::interp_ghost(eb)" << endl;
  }

  if(m_interp_type == ghost_interpolation::pwl){
    this->interp_ghost_pwl(a_data, a_phase);
  }
  else if(m_interp_type == ghost_interpolation::quad){
    this->interp_ghost_quad(a_data, a_phase);
  }
  else{
    MayDay::Abort("amr_mesh::interp_ghost - unsupported interpolation type requested");
  }
}

void amr_mesh::interp_ghost(MFAMRCellData& a_data){
  CH_TIME("amr_mesh::interp_ghost(mf)");
  if(m_verbosity > 3){
    pout() << "amr_mesh::interp_ghost(mf)" << endl;
  }
  
  EBAMRCellData alias_g(1 + m_finest_level);
  EBAMRCellData alias_s(1 + m_finest_level);

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_mfis->get_ebis(phase::solid);
  
  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    alias_g[lvl] = RefCountedPtr<LevelData<EBCellFAB> > (new LevelData<EBCellFAB>());
    alias_s[lvl] = RefCountedPtr<LevelData<EBCellFAB> > (new LevelData<EBCellFAB>());
      
    mfalias::aliasMF(*alias_g[lvl], phase::gas,   *a_data[lvl]);
    if(!ebis_sol.isNull()){
      mfalias::aliasMF(*alias_s[lvl], phase::solid, *a_data[lvl]);
    }
  }

  this->interp_ghost(alias_g, phase::gas);
  if(!ebis_sol.isNull()){
    this->interp_ghost(alias_s, phase::solid);
  }
}

void amr_mesh::interp_ghost_quad(EBAMRCellData& a_data, phase::which_phase a_phase){
  CH_TIME("amr_mesh::interp_ghost(quad)");
  if(m_verbosity > 3){
    pout() << "amr_mesh::interp_ghost(quad)" << endl;
  }

  for (int lvl = m_finest_level; lvl > 0; lvl--){
    const int ncomps = a_data[lvl]->nComp();
    const Interval interv(0, ncomps -1);

    CH_assert(a_data[lvl]->ghostVect() == m_num_ghost*IntVect::Unit);

    m_quadcfi[a_phase][lvl]->coarseFineInterp(*a_data[lvl], *a_data[lvl-1], 0, 0, ncomps);
  }

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    a_data[lvl]->exchange();
  }
}

void amr_mesh::interp_ghost_pwl(EBAMRCellData& a_data, phase::which_phase a_phase){
  CH_TIME("amr_mesh::interp_ghost(pwl)");
  if(m_verbosity > 3){
    pout() << "amr_mesh::interp_ghost(pwl)" << endl;
  }
  
  for (int lvl = m_finest_level; lvl > 0; lvl--){
    const int ncomps = a_data[lvl]->nComp();
    const Interval interv(0, ncomps -1);

    CH_assert(a_data[lvl]->ghostVect() == m_num_ghost*IntVect::Unit);

    m_pwl_fillpatch[a_phase][lvl]->interpolate(*a_data[lvl], *a_data[lvl-1], *a_data[lvl-1], 0.0, 0.0, 0.0, interv);
  }

  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    a_data[lvl]->exchange();
  }  
}

void amr_mesh::interpolate_to_centroids(EBAMRCellData& a_data, phase::which_phase a_phase){
  m_centroid_interp[a_phase]->apply(a_data);
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

void amr_mesh::parse_mg_coarsen(){

  ParmParse pp("amr_mesh");
  int depth;
  pp.get("mg_coarsen", depth);
  if(depth >= 0){
    m_mg_coarsen = depth;
  }
  else{
    m_mg_coarsen = 0;
  }
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

void amr_mesh::parse_refine_all_depth(){

  ParmParse pp("amr_mesh");
  int depth;
  pp.get("refine_all_lvl", depth);
  if(depth > 0){
    m_refine_all_depth = depth;
  }
  else {
    m_refine_all_depth = 0;
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

void amr_mesh::set_grids(Vector<Vector<Box> >& a_boxes, const int a_regsize){

  // Do Morton ordering
  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    mortonOrdering((Vector<Box>&)a_boxes[lvl]);
  }

  // Load balance boxes
  Vector<Vector<int> > proc_assign(1 + m_finest_level);
  this->loadbalance(proc_assign, a_boxes);

  // Define grids
  m_grids.resize(1 + m_finest_level);
  for (int lvl = 0; lvl <= m_finest_level; lvl++){
    m_grids[lvl] = DisjointBoxLayout();
    m_grids[lvl].define(a_boxes[lvl], proc_assign[lvl], m_domains[lvl]);
    m_grids[lvl].close();
  }

  m_has_grids = true;

  const int a_lmin = 0;
  this->define_neighbors(a_lmin);    
  this->define_eblevelgrid(a_lmin);  // Define EBLevelGrid objects on both phases
  this->define_mflevelgrid(a_lmin);  // Define MFLevelGrid
  this->define_vofiter(a_lmin);      // Define vof iterators
  if(!m_has_mg_stuff){
    this->define_mg_stuff();
    m_has_mg_stuff = true;
  }

#if 0  // Done after the operator shit
  this->define_eb_coar_ave(a_lmin); 
  this->define_eb_quad_cfi(a_lmin); 
  this->define_fillpatch(a_lmin);   
  this->define_ebpwl_interp(a_lmin);
  this->define_flux_reg(a_lmin,a_regsize);    
  this->define_redist_oper(a_lmin, a_regsize);
  this->define_gradsten(a_lmin);
  this->define_irreg_sten();            
  this->define_copier(a_lmin);                
  this->define_ghostcloud(a_lmin);              
#endif

  this->define_realms();
  m_realm->regrid_base(a_lmin);

}

void amr_mesh::parse_max_box_size(){
  ParmParse pp("amr_mesh");
  int box_size;
  pp.get("max_box_size", box_size);
  if(box_size >= 8 && box_size % 2 == 0){
    m_max_box_size = box_size;
  }
  else{
    MayDay::Abort("amr_mesh::parse_max_box_size - must have box_size > 8 and divisible by 2");
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
    MayDay::Abort("amr_mesh::parse_max_ebis_box_size - must have box_size > 8 and divisible by 2");
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
}

void amr_mesh::parse_irreg_growth(){
  ParmParse pp("amr_mesh");
  int buffer;
  pp.get("irreg_growth", buffer);
  if(buffer > 0){
    m_irreg_growth = buffer;
  }
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

void amr_mesh::set_irreg_sten_type(const stencil_type::which_type a_type){
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

int amr_mesh::get_refine_all_depth(){
  return m_refine_all_depth;
}

int amr_mesh::get_num_ghost(){
  return m_num_ghost;
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

int amr_mesh::get_refinement_ratio(const int a_level1, const int a_level2){
  int coarLevel = Min(a_level1, a_level2);
  int fineLevel = Max(a_level1, a_level2);

  int ref = 1;
  for (int lvl = coarLevel; lvl < fineLevel; lvl++){
    ref = ref*m_ref_ratios[lvl];
  }

  return ref;
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

Vector<DisjointBoxLayout>& amr_mesh::get_grids(){
  return m_grids;
}

Vector<DisjointBoxLayout>& amr_mesh::get_mg_grids(){
  return m_mg_grids;
}

Vector<ProblemDomain>& amr_mesh::get_domains(){
  return m_domains;
}

Vector<ProblemDomain>& amr_mesh::get_mg_domains(){
  return m_mg_domains;
}

Vector<EBISLayout>& amr_mesh::get_ebisl(phase::which_phase a_phase){
  return m_ebisl[a_phase];
}

Vector<RefCountedPtr<LayoutData<VoFIterator> > > amr_mesh::get_vofit(phase::which_phase a_phase){
  return m_vofiter[a_phase];
}

Vector<RefCountedPtr<LayoutData<Vector<LayoutIndex> > > >& amr_mesh::get_neighbors(){
  return m_neighbors;
}

Vector<RefCountedPtr<EBLevelGrid> >& amr_mesh::get_eblg(phase::which_phase a_phase){
  return m_eblg[a_phase];
}

Vector<RefCountedPtr<EBLevelGrid> >& amr_mesh::get_mg_eblg(phase::which_phase a_phase){
  return m_mg_eblg[a_phase];
}

Vector<RefCountedPtr<MFLevelGrid> >& amr_mesh::get_mflg(){
  return m_mflg;
}

Vector<RefCountedPtr<MFLevelGrid> >& amr_mesh::get_mg_mflg(){
  return m_mg_mflg;
}

Vector<RefCountedPtr<ebcoarseaverage> >& amr_mesh::get_coarave(phase::which_phase a_phase){
  return m_coarave[a_phase];
}

Vector<RefCountedPtr<EBGhostCloud> >& amr_mesh::get_ghostcloud(phase::which_phase a_phase){
  return m_ghostclouds[a_phase];
}

Vector<RefCountedPtr<nwoebquadcfinterp> >& amr_mesh::get_quadcfi(phase::which_phase a_phase){
  return m_quadcfi[a_phase];
}

Vector<RefCountedPtr<EBQuadCFInterp> >& amr_mesh::get_old_quadcfi(phase::which_phase a_phase){
  return m_old_quadcfi[a_phase];
}

Vector<RefCountedPtr<AggEBPWLFillPatch> >& amr_mesh::get_fillpatch(phase::which_phase a_phase){
  return m_pwl_fillpatch[a_phase];
}

Vector<RefCountedPtr<EBPWLFineInterp> >& amr_mesh::get_eb_pwl_interp(phase::which_phase a_phase){
  return m_pwl_interp[a_phase];
}

Vector<RefCountedPtr<EBMGInterp> >& amr_mesh::get_eb_mg_interp(phase::which_phase a_phase){
  return m_ebmg_interp[a_phase];
}

Vector<RefCountedPtr<EBFluxRegister> >&  amr_mesh::get_flux_reg(phase::which_phase a_phase){
  CH_assert(a_phase == phase::gas); // This is disabled since we only solve cdr in gas phase. 
  return m_flux_reg[a_phase];
}

Vector<RefCountedPtr<EBLevelRedist> >& amr_mesh::get_level_redist(phase::which_phase a_phase){
  CH_assert(a_phase == phase::gas); // This is disabled since we only solve cdr in gas phase. 
  return m_level_redist[a_phase];
}

Vector<RefCountedPtr<EBCoarToFineRedist> >&  amr_mesh::get_coar_to_fine_redist(phase::which_phase a_phase){
  CH_assert(a_phase == phase::gas); // This is disabled since we only solve cdr in gas phase.

  return m_coar_to_fine_redist[a_phase];;
}

Vector<RefCountedPtr<EBCoarToCoarRedist> >&  amr_mesh::get_coar_to_coar_redist(phase::which_phase a_phase){
  CH_assert(a_phase == phase::gas); // This is disabled since we only solve cdr in gas phase.

  return m_coar_to_coar_redist[a_phase];;
}

Vector<RefCountedPtr<EBFineToCoarRedist> >&  amr_mesh::get_fine_to_coar_redist(phase::which_phase a_phase){
  CH_assert(a_phase == phase::gas); // This is disabled since we only solve cdr in gas phase.

  return m_fine_to_coar_redist[a_phase];;
}

Vector<RefCountedPtr<EBLevelConcentrationRedist> >& amr_mesh::get_concentration_redist(phase::which_phase a_phase){
  CH_assert(a_phase == phase::gas); // This is disabled since we only solve cdr in gas phase. 
  return m_concentration_redist[a_phase];
}

irreg_amr_stencil<centroid_interp>& amr_mesh::get_centroid_interp_stencils(phase::which_phase a_phase){
  return *m_centroid_interp[a_phase];
}

irreg_amr_stencil<eb_centroid_interp>& amr_mesh::get_eb_centroid_interp_stencils(phase::which_phase a_phase){
  return *m_eb_centroid_interp[a_phase];
}

irreg_amr_stencil<noncons_div>& amr_mesh::get_noncons_div_stencils(phase::which_phase a_phase){
  return *m_noncons_div[a_phase];
}


Vector<RefCountedPtr<Copier> >& amr_mesh::get_copier(phase::which_phase a_phase){
  return m_copier[a_phase];
}


Vector<RefCountedPtr<Copier> >& amr_mesh::get_reverse_copier(phase::which_phase a_phase){
  return m_reverse_copier[a_phase];
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

void amr_mesh::register_operator(const std::string a_operator, const phase::which_phase a_phase){
  CH_TIME("amr_mesh::register_operator");
  if(m_verbosity > 5){
    pout() << "amr_mesh::register_operator" << endl;
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
     

  m_operator_map.emplace(std::make_pair(a_operator, a_phase), true);

  m_realm->register_operator(a_operator, a_phase);
}

bool amr_mesh::query_operator(const std::string a_operator, const phase::which_phase a_phase) {
  CH_TIME("amr_mesh::query_operator");
  if(m_verbosity > 5){
    pout() << "amr_mesh::query_operator" << endl;
  }

  return m_operator_map[std::make_pair(a_operator, a_phase)];
}

void amr_mesh::define_realms(){
  CH_TIME("amr_mesh::define_realms");
  if(m_verbosity > 5){
    pout() << "amr_mesh::define_realms" << endl;
  }

  m_realm->define(m_grids, m_domains, m_ref_ratios, m_dx, m_finest_level, m_ebghost, m_num_ghost, m_redist_rad,
		  m_ebcf, m_centroid_stencil, m_eb_stencil, m_mfis);
}
