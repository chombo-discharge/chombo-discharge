/*!
  @file poisson_multifluid_gmg.cpp
  @brief Implementation of poisson_multifluid_gmg.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "poisson_multifluid_gmg.H"
#include "data_ops.H"
#include "MFQuadCFInterp.H"

#include <MFCellFAB.H>
#include <LayoutData.H>
#include <MFLevelDataOps.H>

poisson_multifluid_gmg::poisson_multifluid_gmg(){
  m_needs_setup = true;
}

poisson_multifluid_gmg::~poisson_multifluid_gmg(){

}

int poisson_multifluid_gmg::query_ghost() const {
  return 2; // Need this many cells
}

void poisson_multifluid_gmg::solve(){
  if(m_needs_setup){
    this->setup_gmg();
  }

  
}

void poisson_multifluid_gmg::solve(MFAMRCellData& a_state, const MFAMRCellData& a_source){

}

void poisson_multifluid_gmg::setup_gmg(){

  const int finest_level = m_amr->get_finest_level();
  const Vector<DisjointBoxLayout>& grids = m_amr->get_grids();

  m_state.resize(1 + finest_level);
  m_source.resize(1 + finest_level);

  // Create data holders
  for (int lvl = 0; lvl <= finest_level; lvl++){
    Vector<EBISLayout> ebisl_phases(2);
    Vector<int> components(2);
    ebisl_phases[0] = m_amr->get_ebisl(Phase::Gas)[lvl];
    ebisl_phases[1] = m_amr->get_ebisl(Phase::Solid)[lvl];
    components[0] = 1;
    components[1] = 1;
    
    MFCellFactory cellfact(ebisl_phases, components);

    m_state[lvl] = RefCountedPtr<LevelData<MFCellFAB> > (new LevelData<MFCellFAB>(grids[lvl], 1, 2*IntVect::Unit, cellfact));
    m_source[lvl] =RefCountedPtr<LevelData<MFCellFAB> > (new LevelData<MFCellFAB>(grids[lvl], 1, 2*IntVect::Unit, cellfact));

    MFLevelDataOps::setVal(*m_state[lvl], 0.0);
    MFLevelDataOps::setVal(*m_source[lvl], 0.0);
  }

  // Things to pass to the operator factory
  const RefCountedPtr<EBIndexSpace>& ebis_gas = m_mfis->get_ebis(Phase::Gas);
  const RefCountedPtr<EBIndexSpace>& ebis_sol = m_mfis->get_ebis(Phase::Solid);
  const Vector<int>& refinement_ratios  = m_amr->get_ref_rat();
  const Vector<ProblemDomain>& domains  = m_amr->get_domains();
  const Vector<Real>& dx                = m_amr->get_dx();
  const RealVect& origin                = m_physdom->get_prob_lo();

  
  RefCountedPtr<BaseDomainBCFactory> domfact = RefCountedPtr<BaseDomainBCFactory> (NULL);
  RefCountedPtr<BaseEBBCFactory>     ebcfact = RefCountedPtr<BaseEBBCFactory> (NULL);


  // Create two MFBaseIVFABs and copy between them
  for (int lvl = 0; lvl <= finest_level; lvl++){

    Vector<EBISLayout> ebisl(2);
    Vector<int> comps(2);

    ebisl[Phase::Gas] = m_amr->get_ebisl(Phase::Gas)[lvl];
    ebisl[Phase::Solid] = m_amr->get_ebisl(Phase::Solid)[lvl];

    comps[Phase::Gas] = 1;
    comps[Phase::Solid] = 1;

    MFBaseIVFABFactory fact(ebisl, comps);

    const int ignored = 0;
    LevelData<MFBaseIVFAB>* mfiv1 = new LevelData<MFBaseIVFAB> (grids[lvl], ignored, 2*IntVect::Unit, fact);
    LevelData<MFBaseIVFAB>* mfiv2 = new LevelData<MFBaseIVFAB> (grids[lvl], ignored, 2*IntVect::Unit, fact);

    data_ops::set_value(*mfiv1, 1.0);
    data_ops::set_value(*mfiv2, 2.0);

    mfiv1->copyTo(*mfiv2);

  }

  // Allocate coefficients
  EBAMRIVData   ajump;
  EBAMRIVData   bjump;
  EBAMRIVData   cjump;
  EBAMRIVData   sigma;
  Vector<RefCountedPtr<LayoutData<IntVectSet> > > jump_cells(1 + finest_level);
  for (int lvl = 0; lvl <= finest_level; lvl++){
    jump_cells[lvl] = RefCountedPtr<LayoutData<IntVectSet> > (new LayoutData<IntVectSet> (grids[lvl]));
    const IntVectSet interface_ivs = m_mfis->interface_region(domains[lvl]);
    for (DataIterator dit = jump_cells[lvl]->dataIterator(); dit.ok(); ++dit){
      (*jump_cells[lvl])[dit()] = interface_ivs & grids[lvl].get(dit());
    }
  }
  MFAMRCellData aco;
  MFAMRFluxData bco;
  MFAMRIVData   bco_irreg;
  Real          alpha =  0.0;
  Real          beta  = -1.0;

  const int ncomps = 1;
  const int ghosts = 1;
  m_amr->allocate(aco,       ncomps, ghosts);
  m_amr->allocate(bco,       ncomps, ghosts);
  m_amr->allocate(bco_irreg, ncomps, ghosts);
  m_amr->allocate_interface(sigma, Phase::Gas, 1, 0);

  
  Vector<MFLevelGrid> mfeblg(1 + finest_level);
  Vector<MFQuadCFInterp> mfquadcfi(1 + finest_level);
  for (int lvl = 0; lvl <= finest_level; lvl++){
    Vector<EBLevelGrid> eblg_phases(Phase::num_phases);
    Vector<RefCountedPtr<EBQuadCFInterp> > quadcfi_phases(Phase::num_phases);

    
    eblg_phases[Phase::Gas] = *(m_amr->get_eblg(Phase::Gas)[lvl]);
    eblg_phases[Phase::Solid] = *(m_amr->get_eblg(Phase::Solid)[lvl]);

    quadcfi_phases[Phase::Gas] =   (m_amr->get_old_quadcfi(Phase::Gas)[lvl]);
    quadcfi_phases[Phase::Solid] = (m_amr->get_old_quadcfi(Phase::Solid)[lvl]);
    
    mfeblg[lvl].define(eblg_phases);
    mfquadcfi[lvl].define(quadcfi_phases);
  }

  m_opfact = RefCountedPtr<mf_helmholtz_opfactory> (new mf_helmholtz_opfactory(m_mfis,
									       mfeblg,
									       mfquadcfi,
									       refinement_ratios,
									       grids,
									       aco,
									       bco,
									       bco_irreg,
									       alpha,
									       beta,
									       dx[0],
									       domains[0],
  									       domfact,
									       ebcfact,
  									       origin,
  									       2*IntVect::Unit,
  									       2*IntVect::Unit));

  m_opfact->set_jump(1.0, 1.0);
									       
									       
									       
}
