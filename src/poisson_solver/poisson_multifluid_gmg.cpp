/*!
  @file poisson_multifluid_gmgs.cpp
  @brief Implementation of poisson_multifluid_gmg.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "poisson_multifluid_gmg.H"

poisson_multifluid_gmg::poisson_multifluid_gmg(){

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

  const int finest_level          = m_amr->get_finest_level();
  Vector<DisjointBoxLayout> grids = m_amr->get_grids();
  Vector<EBISLayout> ebisl_gas    = m_amr->get_ebisl(Phase::Gas);
  Vector<EBISLayout> ebisl_sol    = m_amr->get_ebisl(Phase::Solid);



  for (int lvl = 0; lvl <= finest_level; lvl++){
    Vector<int> num_comps(Phase::num_phases, 1);
    Vector<EBISLayout> ebisl_phases(Phase::num_phases);
    ebisl_phases[Phase::Gas]   = ebisl_gas[lvl];
    ebisl_phases[Phase::Solid] = ebisl_sol[lvl];
    
    MFCellFactory fact(ebisl_phases, num_comps);

    const int ignored = 0;
    m_state[lvl] = new LevelData<MFCellFAB>(grids[lvl], ignored, this->query_ghost()*IntVect::Unit, fact);
  }
}
