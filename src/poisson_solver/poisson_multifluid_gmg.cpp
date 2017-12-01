/*!
  @file poisson_multifluid_gmg.cpp
  @brief Implementation of poisson_multifluid_gmg.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "poisson_multifluid_gmg.H"

#include <MFCellFAB.H>
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
  
  m_opfact = RefCountedPtr<mf_helmholtz_opfactory> (new mf_helmholtz_opfactory(ebis_gas,
									       ebis_sol,
									       refinement_ratios,
									       domains[0],
									       dx[0],
									       origin,
									       domfact,
									       ebcfact,
									       2*IntVect::Unit,
									       2*IntVect::Unit));
									       
									       
									       
}
