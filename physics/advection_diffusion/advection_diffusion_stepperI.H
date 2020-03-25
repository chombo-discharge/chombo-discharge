/*!
  @file   advection_diffusion_stepperI.H
  @brief  Implementation of advection_diffusion_stepper
  @author Robert Marskar
  @data   March 2020
*/

#include <ParmParse.H>

#include "advection_diffusion_species.H"

namespace physics::advection_diffusion {

  template <class T>
  advection_diffusion_stepper<T>::advection_diffusion_stepper(){
    ParmParse pp("advection_diffusion_stepper");

    pp.get("diffco",   m_diffco);
    pp.get("omega",    m_omega);
    pp.get("cfl",      m_cfl);
  }

  template <class T>
  advection_diffusion_stepper<T>::~advection_diffusion_stepper(){
  
  }

  template <class T>
  void advection_diffusion_stepper<T>::setup_solvers(){
    m_solver  = RefCountedPtr<cdr_solver>  (new T());
    m_species = RefCountedPtr<cdr_species> (new advection_diffusion_species());

    m_solver->set_species(m_species);

    if(!m_solver->is_mobile() && !m_solver->is_diffusive()){
      MayDay::Abort("advection_diffusion_stepper<T>::setup_solvers - can't turn off both advection AND diffusion");
    }

    m_solver->parse_options();
    m_solver->set_amr(m_amr);
    m_solver->set_computational_geometry(m_compgeom);
    m_solver->set_phase(phase::gas);
    m_solver->sanity_check();
    m_solver->allocate_internals();
  }

  template <class T>
  void advection_diffusion_stepper<T>::initial_data(){
    m_solver->initial_data();       // Fill initial through the cdr species
    m_solver->set_diffco(m_diffco);
    m_solver->set_source(0.0);
    m_solver->set_ebflux(0.0);

    // Set the velocity
    this->set_velocity();
  }

  template <class T>
  void advection_diffusion_stepper<T>::set_velocity(){
    const int finest_level = m_amr->get_finest_level();
    for (int lvl = 0; lvl <= finest_level; lvl++){
      this->set_velocity(lvl);
    }
  }

  template<class T>
  void advection_diffusion_stepper<T>::set_velocity(const int a_level){
    // TLDR: This code goes down to each cell on grid level a_level and sets the velocity to omega*r
    
    const DisjointBoxLayout& dbl = m_amr->get_grids()[a_level];
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const Box& box = dbl.get(dit());

      EBCellFAB& vel = (*(m_solver->get_velo_cell())[a_level])[dit()];
      BaseFab<Real>& vel_reg = vel.getSingleValuedFAB();

      vel.setVal(0.0);
      
      // Regular cells
      for (BoxIterator bit(box); bit.ok(); ++bit){
	const IntVect iv = bit();
	const RealVect pos = m_amr->get_prob_lo() + (RealVect(iv) + 0.5*RealVect::Unit)*m_amr->get_dx()[a_level];

	const Real r     = sqrt(pos[0]*pos[0] + pos[1]*pos[1]);
	const Real theta = atan(pos[1]/pos[0]);

	vel_reg(iv,0) = -r*m_omega*sin(theta);
	vel_reg(iv,1) =  r*m_omega*cos(theta);
      }

      // Irregular and multicells
      const EBISBox& ebisbox = vel.getEBISBox();
      const EBGraph& ebgraph = ebisbox.getEBGraph();
      for (VoFIterator vofit(ebisbox.getIrregIVS(box), ebgraph); vofit.ok(); ++vofit){

	const VolIndex vof = vofit();
	const IntVect iv   = vof.gridIndex();
	const RealVect pos = m_amr->get_prob_lo() + (RealVect(iv) + 0.5*RealVect::Unit)*m_amr->get_dx()[a_level];

	const Real r     = sqrt(pos[0]*pos[0] + pos[1]*pos[1]);
	const Real theta = atan(pos[1]/pos[0]);

	vel(vof,0) = -r*m_omega*sin(theta);
	vel(vof,1) =  r*m_omega*cos(theta);
      }
    }
  }

  template <class T>
  void advection_diffusion_stepper<T>::write_checkpoint_data(HDF5Handle& a_handle, const int a_lvl) const {
    m_solver->write_checkpoint_level(a_handle, a_lvl);
  }

  template <class T>
  void advection_diffusion_stepper<T>::read_checkpoint_data(HDF5Handle& a_handle, const int a_lvl){
    m_solver->read_checkpoint_level(a_handle, a_lvl);
  }

  template <class T>
  void advection_diffusion_stepper<T>::post_checkpoint_setup(){

    MayDay::Abort("advection_diffusion_stepper<T>::post_checkpoint_setup - need to fill velo/diffco");
  }

  template <class T>
  int advection_diffusion_stepper<T>::get_num_plot_vars() const{
    return m_solver->get_num_plotvars();
  }

  template <class T>
  void advection_diffusion_stepper<T>::write_plot_data(EBAMRCellData&       a_output,
						       Vector<std::string>& a_plotvar_names,
						       int&                 a_icomp) const {
    a_plotvar_names.append(m_solver->get_plotvar_names());
    m_solver->write_plot_data(a_output, a_icomp);
  }

  template <class T>
  void advection_diffusion_stepper<T>::compute_dt(Real& a_dt, time_code::which_code& a_timecode){
    MayDay::Abort("advection_diffusion_stepper<T>::post_checkpoint_setup - need to restrict dt");

    // CFL on advection
    if(m_solver->is_mobile()){
      a_dt       = m_solver->compute_cfl_dt();
      a_timecode = time_code::cfl;
    }
    
    // CFL on diffusion, if explicit diffusion
    if(m_solver->is_diffusive()){
      const Real diff_dt = m_solver->compute_diffusive_dt();

      if(diff_dt < a_dt){
	a_dt = diff_dt;
	a_timecode = time_code::diffusion;
      }
    }
  }

  template <class T>
  Real advection_diffusion_stepper<T>::advance(const Real a_dt){
    MayDay::Abort("advection_diffusion_stepper<T>::advance - not implemented");
    return 0.0;
  }

  template <class T>
  void advection_diffusion_stepper<T>::synchronize_solver_times(const int a_step, const Real a_time, const Real a_dt){
    m_solver->set_time(a_step, a_time, a_dt);
  }

  template <class T>
  void advection_diffusion_stepper<T>::print_step_report(){

  }

  // Regrid routines
  template <class T>
  bool advection_diffusion_stepper<T>::need_to_regrid(){
    return false;
  }

  template <class T>
  void advection_diffusion_stepper<T>::cache(){
    m_solver->cache_state();
  }

  template <class T>
  void advection_diffusion_stepper<T>::deallocate(){
    m_solver->deallocate_internals();
  }

  template <class T>
  void advection_diffusion_stepper<T>::regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level){
    m_solver->regrid(a_lmin, a_old_finest_level, a_new_finest_level);
  }
}
