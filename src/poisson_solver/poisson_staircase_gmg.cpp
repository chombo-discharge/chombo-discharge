/*!
  @file poisson_staircase_gmgs.cpp
  @brief Implementation of poisson_staircase_gmg.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "poisson_staircase_gmg.H"
#include "data_ops.H"

#include <DirichletConductivityDomainBC.H>
#include <DirichletConductivityEBBC.H>

poisson_staircase_gmg::poisson_staircase_gmg(){
  CH_TIME("poisson_staircase_gmg::poisson_staircase_gmg");

  
  this->set_verbosity(-1); // Shut up.
  
  m_needs_setup = true;
}

poisson_staircase_gmg::~poisson_staircase_gmg(){

}

int poisson_staircase_gmg::query_ghost() const {
  CH_TIME("poisson_staircase_gmg::query_ghost");
  if(m_verbosity > 5){
    pout() << "poisson_staircase_gmg::query_ghost" << endl;
  }
  
  return 2;
}

void poisson_staircase_gmg::define_coefficients(){
  CH_TIME("poisson_staircase_gmg::define_coefficients");
  if(m_verbosity > 5){
    pout() << "poisson_staircase_gmg::define_coefficients" << endl;
  }

  m_alpha =  0.0;
  m_beta  = -1.0;

  m_amr->allocate(m_aco,       Phase::Gas, 1, m_amr->get_num_ghost());
  m_amr->allocate(m_bco,       Phase::Gas, 1, m_amr->get_num_ghost());
  m_amr->allocate(m_bco_irreg, Phase::Gas, 1, m_amr->get_num_ghost());

  data_ops::set_value(m_aco,       1.0);
  data_ops::set_value(m_bco,       1.0);
  data_ops::set_value(m_bco_irreg, 1.0);
}

void poisson_staircase_gmg::solve(){
  CH_TIME("poisson_staircase_gmg::solve()");
  if(m_verbosity > 5){
    pout() << "poisson_staircase_gmg::solve()" << endl;
  }

  if(m_needs_setup){
    this->define_coefficients();
    this->setup_bc();
    this->setup_gmg();
  }


	      
}

void poisson_staircase_gmg::solve(MFAMRCellData& a_state, const MFAMRCellData& a_source){
  CH_TIME("poisson_staircase_gmg::solve(phi, src)");
  if(m_verbosity > 5){
    pout() << "poisson_staircase_gmg::solve(phi, src)" << endl;
  }

  MayDay::Abort("poisson_staircase_gmg::solve(phi, src) - not implemented");
}

void poisson_staircase_gmg::setup_bc(){
  CH_TIME("poisson_staircase_gmg::setup_bc");
  if(m_verbosity > 5){
    pout() << "poisson_staircase_gmg::setup_bc" << endl;
  }
  
  // Boundary conditions. Dirichlet everywhere for testing purposes
  DirichletConductivityDomainBCFactory* dombc_fact = new DirichletConductivityDomainBCFactory();
  DirichletConductivityEBBCFactory* ebbc_fact      = new DirichletConductivityEBBCFactory();

  dombc_fact->setValue(0.0);
  ebbc_fact->setValue(1.0);
  ebbc_fact->setOrder(2);

  m_domain_bc_factory = RefCountedPtr<BaseDomainBCFactory> (dombc_fact);
  m_eb_bc_factory     = RefCountedPtr<BaseEBBCFactory>     (ebbc_fact);

}

void poisson_staircase_gmg::setup_gmg(){
  CH_TIME("poisson_staircase_gmg::setup_gmg");
  if(m_verbosity > 5){
    pout() << "poisson_staircase_gmg::setup_gmg" << endl;
  }

  const int comps                  = 1;
  const int finest_level           = m_amr->get_finest_level();
  const int ghost                  = m_amr->get_num_ghost();
  Vector<int> ref_ratios           = m_amr->get_ref_rat();
  Vector<Real>& dx                 = m_amr->get_dx();
  Vector<ProblemDomain>& domains   = m_amr->get_domains();
  Vector<DisjointBoxLayout>& grids = m_amr->get_grids();
  Vector<EBISLayout>& ebisl        = m_amr->get_ebisl(Phase::Gas);
  Vector<EBLevelGrid> levelgrids;

  
  for (int lvl = 0; lvl <= finest_level; lvl++){ 
    levelgrids.push_back(*(m_amr->get_eblg(Phase::Gas)[lvl])); // amr_mesh uses RefCounted levelgrids
  }

  // This should be moved to amr_mesh
  Vector<RefCountedPtr<EBQuadCFInterp> > quadcfi(1 + finest_level);
  for (int lvl = 0; lvl <= finest_level; lvl++){
    if(lvl > 0){
      quadcfi[lvl] = RefCountedPtr<EBQuadCFInterp> (new EBQuadCFInterp(grids[lvl],
								       grids[lvl-1],
								       ebisl[lvl],
								       ebisl[lvl-1],
								       domains[lvl],
								       m_amr->get_refinement_ratio(),
								       1,
								       *(levelgrids[lvl].getCFIVS()),
								       m_mfis->get_ebis(Phase::Gas)));
    }
  }
								       

  m_cond_op_fact = RefCountedPtr<EBConductivityOpFactory> (new EBConductivityOpFactory(levelgrids,
										       quadcfi,
										       m_alpha,
										       m_beta,
										       m_aco,
										       m_bco,
										       m_bco_irreg,
										       dx[0],
										       ref_ratios,
										       m_domain_bc_factory,
										       m_eb_bc_factory,
										       ghost*IntVect::Unit,
										       ghost*IntVect::Unit,
										       1));

  
  m_gmg_solver.define(domains[0], *m_cond_op_fact, &m_bicgstab, 1 + finest_level);


  EBAMRCellData phi, src, res;
  m_amr->allocate(phi, Phase::Gas, comps, ghost);
  m_amr->allocate(src, Phase::Gas, comps, ghost);
  m_amr->allocate(res, Phase::Gas, comps, ghost);

  data_ops::set_value(src, 0.0);
  data_ops::set_value(phi, 0.0);
  data_ops::set_value(res, 0.0);

  Vector<LevelData<EBCellFAB>* > phi_ptr, src_ptr, res_ptr;
  for (int lvl = 0; lvl <= finest_level; lvl++){
    phi_ptr.push_back(&(*phi[lvl]));
    src_ptr.push_back(&(*src[lvl]));
    res_ptr.push_back(&(*res[lvl]));
  }


  m_gmg_solver.setSolverParameters(32, 32, 32, 1, 100, 1.E-30, 1.E-30, 1.E-30);
  m_gmg_solver.m_verbosity = 10;
  m_gmg_solver.solve(phi_ptr, src_ptr, finest_level, 0);

  
  // Write data
  Vector<std::string> names(1);
  Vector<Real> covered_values;
  names[0] = "potential";
  m_amr->average_down(phi, Phase::Gas);
  writeEBHDF5("potential.hdf5",
	      grids,
	      phi_ptr,
	      names,
	      domains[0],
	      dx[0],
	      0.0,
	      0.0,
	      m_amr->get_ref_rat(),
	      1 + finest_level,
	      false,
	      covered_values);
  
}
