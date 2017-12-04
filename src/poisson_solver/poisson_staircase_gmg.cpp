/*!
  @file poisson_staircase_gmgs.cpp
  @brief Implementation of poisson_staircase_gmg.H
  @author Robert Marskar
  @date Nov. 2017
  @todo Map boundary conditions onto computational_geometry and plasma_engine. 
*/

#include "poisson_staircase_gmg.H"
#include "data_ops.H"

#include <DirichletConductivityDomainBC.H>
#include <DirichletConductivityEBBC.H>

poisson_staircase_gmg::poisson_staircase_gmg(){
  CH_TIME("poisson_staircase_gmg::poisson_staircase_gmg");

  
  this->set_verbosity(-1); // Shut up.
  this->set_gmg_solver_parameters();
  
  m_needs_setup = true;
}

poisson_staircase_gmg::~poisson_staircase_gmg(){

}

void poisson_staircase_gmg::set_gmg_solver_parameters(relax::which_relax a_relax_type,
						      amrmg::which_mg a_gmg_type,      
						      const int a_verbosity,          
						      const int a_pre_smooth,         
						      const int a_post_smooth,       
						      const int a_bot_smooth,         
						      const int a_max_iter,           
						      const Real a_eps,               
						      const Real a_hang,              
						      const Real a_norm_thresh){
  CH_TIME("poisson_staircase_gmg::set_gmg_solver_parameters");
  if(m_verbosity > 5){
    pout() << "poisson_staircase_gmg::set_gmg_solver_parameters" << endl;
  }

  m_gmg_relax_type  = a_relax_type;
  m_gmg_type        = a_gmg_type;
  m_gmg_verbosity   = a_verbosity;
  m_gmg_pre_smooth  = a_pre_smooth;
  m_gmg_post_smooth = a_post_smooth;
  m_gmg_bot_smooth  = a_bot_smooth;
  m_gmg_max_iter    = a_max_iter;
  m_gmg_eps         = a_eps;
  m_gmg_hang        = a_hang;
  m_gmg_norm_thresh = a_norm_thresh;
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

    m_needs_setup = false;
  }

  
#if 1 // This is a test, GMG is now set up correctly

  const int comps                  = 1;
  const int ghost                  = m_amr->get_num_ghost();
  const int finest_level           = m_amr->get_finest_level();
  Vector<int> ref_ratios           = m_amr->get_ref_rat();
  Vector<Real>& dx                 = m_amr->get_dx();
  Vector<EBISLayout>& ebisl        = m_amr->get_ebisl(Phase::Gas);
  Vector<ProblemDomain>& domains   = m_amr->get_domains();
  Vector<DisjointBoxLayout>& grids = m_amr->get_grids();

  
  EBAMRCellData phi, src, res, E;
  m_amr->allocate(phi, Phase::Gas, comps,    ghost);
  m_amr->allocate(src, Phase::Gas, comps,    ghost);
  m_amr->allocate(res, Phase::Gas, comps,    ghost);
  m_amr->allocate(E,   Phase::Gas, SpaceDim, ghost);

  data_ops::set_value(src, 0.0);
  data_ops::set_value(phi, 0.0);
  data_ops::set_value(res, 0.0);

  Vector<LevelData<EBCellFAB>* > phi_ptr, src_ptr, res_ptr, E_ptr;
  for (int lvl = 0; lvl <= finest_level; lvl++){
    phi_ptr.push_back(&(*phi[lvl]));
    src_ptr.push_back(&(*src[lvl]));
    res_ptr.push_back(&(*res[lvl]));
    E_ptr.push_back(&(*E[lvl]));
  }

  m_gmg_solver.getInfo();

  m_gmg_solver.solve(phi_ptr, src_ptr, finest_level, 0);

  // Compute gradient
  m_amr->average_down(phi, Phase::Gas);
  m_amr->interp_ghost(phi, Phase::Gas);
  m_amr->compute_gradient(E, phi);
  m_amr->average_down(E, Phase::Gas);
  m_amr->interp_ghost(E, Phase::Gas);

  irreg_amr_stencil<centroid_interp> stencils = m_amr->get_centroid_interp_stencils(Phase::Gas);
  stencils.apply(E);
  m_amr->average_down(E, Phase::Gas);
  
  // Write data
  Vector<std::string> names(SpaceDim);
  Vector<Real> covered_values;
  names[0] = "x-E";
  names[1] = "y-E";
  if(SpaceDim == 3){
    names[2] = "z-E";
  }
  m_amr->average_down(phi, Phase::Gas);
  writeEBHDF5("E.hdf5",
	      grids,
	      E_ptr,
	      names,
	      domains[0],
	      dx[0],
	      0.0,
	      0.0,
	      m_amr->get_ref_rat(),
	      1 + finest_level,
	      false,
	      covered_values);
#endif
	      
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
  Vector<EBISLayout>& ebisl        = m_amr->get_ebisl(Phase::Gas);
  Vector<EBLevelGrid> levelgrids;

  
  for (int lvl = 0; lvl <= finest_level; lvl++){ 
    levelgrids.push_back(*(m_amr->get_eblg(Phase::Gas)[lvl])); // amr_mesh uses RefCounted levelgrids. EBConductivityOp does not. 
  }


  m_cond_op_fact = RefCountedPtr<EBConductivityOpFactory> (new EBConductivityOpFactory(levelgrids,
										       m_amr->get_old_quadcfi(Phase::Gas),
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
										       m_gmg_relax_type));

  m_gmg_solver.define(domains[0], *m_cond_op_fact, &m_bicgstab, 1 + finest_level);
  m_gmg_solver.setSolverParameters(m_gmg_pre_smooth, m_gmg_post_smooth, m_gmg_bot_smooth, m_gmg_type, m_gmg_max_iter,
				   m_gmg_eps, m_gmg_hang, m_gmg_norm_thresh);
  m_gmg_solver.m_verbosity = m_gmg_verbosity;
}
