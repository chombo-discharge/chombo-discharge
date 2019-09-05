/*!
  @file   cdr_fhd.cpp
  @brief  Implementation of cdr_fhd.H
  @author Robert Marskar
  @date   Jan. 2018
*/

#include "cdr_fhd.H"
#include "cdr_fhdF_F.H"
#include "data_ops.H"
#include <time.h>
#include <chrono>

#include <ParmParse.H>

cdr_fhd::cdr_fhd() : cdr_gdnv() {
  m_name       = "cdr_fhd";
  m_class_name = "cdr_fhd";
}

cdr_fhd::~cdr_fhd(){
  this->delete_covered();
}

void cdr_fhd::parse_options(){
  CH_TIME("cdr_fhd::parse_options");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_options" << endl;
  }
  
  parse_diffusion();    // Parses stochastic diffusion
  parse_rng_seed();     // Parses RNG seed
  parse_domain_bc();    // Parses domain BC options
  parse_mass_redist();  // Parses mass redistribution
  parse_hybrid_div();   // Parses options for hybrid divergence
  parse_slopelim();     // Parses slope limiter settings
  parse_plot_vars();    // Parses plot variables
  parse_gmg_settings(); // Parses solver parameters for geometric multigrid
  parse_plotmode();    // Parse plot mode
}

void cdr_fhd::advance_euler(EBAMRCellData& a_new_state, const EBAMRCellData& a_old_state, const Real a_dt){
  CH_TIME("cdr_fhd::advance_euler");
  if(m_verbosity > 5){
    pout() << m_name + "::advance_euler (no source)" << endl;
  }
  
  if(!m_stochastic_diffusion){
    cdr_tga::advance_euler(a_new_state, a_old_state, a_dt);
  }
  else{
    bool converged = false;

    const int comp         = 0;
    const int ncomp        = 1;
    const int finest_level = m_amr->get_finest_level();

    EBAMRCellData ransource;
    m_amr->allocate(ransource, m_phase, 1);
    this->GWN_diffusion_source(ransource, a_old_state);


    // Do the regular aliasing stuff for passing into AMRMultiGrid
    Vector<LevelData<EBCellFAB>* > new_state, old_state, source;
    m_amr->alias(new_state, a_new_state);
    m_amr->alias(old_state, a_old_state);
    m_amr->alias(source,    ransource);
      
    const Real alpha = 0.0;
    const Real beta  = 1.0;

    data_ops::set_value(m_diffco_eb, 0.0);
      
    // Multigrid solver
    m_eulersolver->resetAlphaAndBeta(alpha, beta);
    m_eulersolver->oneStep(new_state, old_state, source, a_dt, 0, finest_level, false);
      
    const int status = m_gmg_solver->m_exitStatus;  // 1 => Initial norm sufficiently reduced
    if(status == 1 || status == 8 || status == 9){  // 8 => Norm sufficiently small
      converged = true;
    }
    MayDay::Abort("cdr_fhd::advance_euler - stochastic not implemented");
  }

  data_ops::floor(a_new_state, 0.0);
}

void cdr_fhd::advance_euler(EBAMRCellData&       a_new_state,
			    const EBAMRCellData& a_old_state,
			    const EBAMRCellData& a_source,
			    const Real           a_dt){
  CH_TIME("cdr_fhd::advance_euler");
  if(m_verbosity > 5){
    pout() << m_name + "::advance_euler (with source)" << endl;
  }

  
  if(m_diffusive){
    if(!m_stochastic_diffusion){
      cdr_tga::advance_euler(a_new_state, a_old_state, a_source, a_dt);
    }
    else{
      bool converged = false;

      const int comp         = 0;
      const int ncomp        = 1;
      const int finest_level = m_amr->get_finest_level();

      EBAMRCellData ransource;
      m_amr->allocate(ransource, m_phase, 1);
      this->GWN_diffusion_source(ransource, a_old_state);

      // Make new state = old_state + dt*random_source
#if 0
      data_ops::copy(a_new_state, a_old_state);
      data_ops::incr(a_new_state, ransource, a_dt);
#else
      data_ops::incr(ransource, a_source, 1.0);
#endif
      
      // Do the regular aliasing stuff for passing into AMRMultiGrid
      Vector<LevelData<EBCellFAB>* > new_state, old_state, source;
      m_amr->alias(new_state, a_new_state);
      m_amr->alias(old_state, a_old_state);
      m_amr->alias(source,    ransource);
      
      const Real alpha = 0.0;
      const Real beta  = 1.0;

      data_ops::set_value(m_diffco_eb, 0.0);
      
      // Multigrid solver
      m_eulersolver->resetAlphaAndBeta(alpha, beta);
      m_eulersolver->oneStep(new_state, old_state, source, a_dt, 0, finest_level, false);
      
      const int status = m_gmg_solver->m_exitStatus;  // 1 => Initial norm sufficiently reduced
      if(status == 1 || status == 8 || status == 9){  // 8 => Norm sufficiently small
	converged = true;
      }
    }
  }
  else{
    data_ops::copy(a_new_state, a_old_state);
  }

  data_ops::floor(a_new_state, 0.0);
}

void cdr_fhd::advance_tga(EBAMRCellData& a_new_state, const EBAMRCellData& a_old_state, const Real a_dt){
  CH_TIME("cdr_fhd::advance_tga");
  if(m_verbosity > 5){
    pout() << m_name + "::advance_tga (no source)" << endl;
  }
  
  if(!m_stochastic_diffusion){
    cdr_tga::advance_tga(a_new_state, a_old_state, a_dt);
  }
  else{
    bool converged = false;

    const int comp         = 0;
    const int ncomp        = 1;
    const int finest_level = m_amr->get_finest_level();

    EBAMRCellData ransource;
    m_amr->allocate(ransource, m_phase, 1);
    this->GWN_diffusion_source(ransource, a_old_state);
      
    // Do the regular aliasing stuff for passing into AMRMultiGrid
    Vector<LevelData<EBCellFAB>* > new_state, old_state, source;
    m_amr->alias(new_state, a_new_state);
    m_amr->alias(old_state, a_old_state);
    m_amr->alias(source,    ransource);
      
    const Real alpha = 0.0;
    const Real beta  = 1.0;

    data_ops::set_value(m_diffco_eb, 0.0);
      
    // Multigrid solver
    m_tgasolver->resetAlphaAndBeta(alpha, beta);
    m_tgasolver->oneStep(new_state, old_state, source, a_dt, 0, finest_level, false);
      
    const int status = m_gmg_solver->m_exitStatus;  // 1 => Initial norm sufficiently reduced
    if(status == 1 || status == 8 || status == 9){  // 8 => Norm sufficiently small
      converged = true;
    }
  }
}

void cdr_fhd::advance_tga(EBAMRCellData&       a_new_state,
			  const EBAMRCellData& a_old_state,
			  const EBAMRCellData& a_source,
			  const Real           a_dt){
  CH_TIME("cdr_fhd::advance_tga");
  if(m_verbosity > 5){
    pout() << m_name + "::advance_tga (with source)" << endl;
  }
  
  if(!m_stochastic_diffusion){
    cdr_tga::advance_tga(a_new_state, a_old_state, a_source, a_dt);
  }
  else{
    bool converged = false;

    const int comp         = 0;
    const int ncomp        = 1;
    const int finest_level = m_amr->get_finest_level();

    EBAMRCellData ransource;
    m_amr->allocate(ransource, m_phase, 1);
    this->GWN_diffusion_source(ransource, a_old_state);
    data_ops::incr(ransource, a_source, 1.0);
      
    // Do the regular aliasing stuff for passing into AMRMultiGrid
    Vector<LevelData<EBCellFAB>* > new_state, old_state, source;
    m_amr->alias(new_state, a_new_state);
    m_amr->alias(old_state, a_old_state);
    m_amr->alias(source,    ransource);
      
    const Real alpha = 0.0;
    const Real beta  = 1.0;

    data_ops::set_value(m_diffco_eb, 0.0);
      
    // Multigrid solver
    m_tgasolver->resetAlphaAndBeta(alpha, beta);
    m_tgasolver->oneStep(new_state, old_state, source, a_dt, 0, finest_level, false);
      
    const int status = m_gmg_solver->m_exitStatus;  // 1 => Initial norm sufficiently reduced
    if(status == 1 || status == 8 || status == 9){  // 8 => Norm sufficiently small
      converged = true;
    }
  }
}



void cdr_fhd::compute_divF(EBAMRCellData& a_divF, const EBAMRCellData& a_state, const Real a_extrap_dt, const bool a_redist){
  CH_TIME("cdr_fhd::compute_divF(divF, state)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_divF(divF, state)" << endl;
  }

  if(m_mobile){
    cdr_tga::compute_divF(a_divF, a_state, a_extrap_dt, a_redist); // Deterministic flux
    if(m_stochastic_advection){ // Stochastic flux
      EBAMRCellData ransource;
      m_amr->allocate(ransource, m_phase, 1);
      this->GWN_advection_source(ransource, a_state);
      data_ops::incr(a_divF, ransource, 1.0);
    }
  }
  else{
    data_ops::set_value(a_divF, 0.0);
  }
}
