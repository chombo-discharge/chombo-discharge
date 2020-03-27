/*!
  @file   euler_subcycle.cpp
  @brief  Implementation of euler_subcycle
  @author Robert Marskar
  @data   March 2020
*/

#include <ParmParse.H>

#include "euler_subcycle.H"
#include "advection_diffusion_species.H"
#include "data_ops.H"

using namespace physics::advection_diffusion;

euler_subcycle::euler_subcycle(){
  ParmParse pp("advection_diffusion");

  pp.get("diffco",   m_diffco);
  pp.get("omega",    m_omega);
  pp.get("cfl",      m_cfl);
}

euler_subcycle::euler_subcycle(RefCountedPtr<cdr_solver>& a_solver) : euler_subcycle() {
  m_solver = a_solver;
}

euler_subcycle::~euler_subcycle(){
  
}

void euler_subcycle::setup_solvers(){
  advection_diffusion_stepper::setup_solvers();

  // Allocate memory for RK steps
  m_amr->allocate(m_k1,  phase::gas, 1);
  m_amr->allocate(m_k2,  phase::gas, 1);
}

void euler_subcycle::compute_dt(Real& a_dt, time_code::which_code& a_timecode){
  advection_diffusion_stepper::compute_dt(a_dt, a_timecode);

  // I assume that a_dt was restricted on the finest level.
  Real bigstep_dt = a_dt;
  for (int lvl = m_amr->get_finest_level(); lvl > 0; lvl--){
    bigstep_dt *= m_amr->get_ref_rat()[lvl-1];
  }

  //  a_dt = bigstep_dt;
}

Real euler_subcycle::advance(const Real a_dt){

  // a_dt is the coarse step. Each level takes its own steps
  const int finest_level = m_amr->get_finest_level();
  const Vector<int> ref_rat = m_amr->get_ref_rat();
  Vector<int> steps(1+finest_level, 1);
  Vector<Real> dt(1 + finest_level, a_dt);
  for (int lvl = 1; lvl <= m_amr->get_finest_level(); lvl++){
    dt[lvl]    = dt[lvl-1]*m_amr->get_ref_rat()[lvl-1];
    steps[lvl] = steps[lvl-1]*ref_rat[lvl-1];
  }

  // Here are the algorithmic steps for an explicit code
  //  for (int lvl = 0; lvl <= m_amr->get_finest_level()){
  //
  //     1. Compute fluxes
  // 
  //           F^l = -v*phi + D*grad(phi)
  //
  //     2. Compute Dc(F^l) and Dnc(F^l) as usual
  //
  //     3. Do an Euler update
  //          
  //           phi^l = phi^l - dt^l*Dh(F^l)
  //
  //        where Dh is the hybrid divergence
  //
  //     4. Compute redistribution mass
  //  
  //           dM^l = dt^l(1-kappa)*[Dnc(F^l) - Dc(F^l)]
  //
  //     5. Do level distribution on level l
  //
  //           phi^l += dM^l
  //        
  //     6. Increment the flux register between (l,l-1) with the F^l flux, i.e. the fine side flux
  //        from the viewpoint of level (l-1)
  //
  //           dF^(l,l-1) += <F^l>*dt, where <F^l> is the average over fine faces
  //
  //     7. Initialize the flux register between (l+1,l) with the F^l flux. I.e. this is the coarse flux
  //        from the viewpoint of level l,
  //
  //           dF^(l+1,l) = F^l*dt
  //
  //     8. Increment the redistriution register between (l,l-1), i.e. the fine2coar register, i.e.
  //
  //           dMC2F^(l,l-1) += dM^l   (this is mass going from this level and into the coarser level)
  // 
  //     9. Initialize the redistribution register between this level and the next finer lever
  //
  //           dMC2F^(l,l+1) = dM^l  (this is mass moving from l and into l+1 in one substep on level l)
  //           dMF2C^(l+1,l) = 0     (reset register describing mass from the fine level and into this level)
  //           dMC2C^(l,l)   = -dM^l (set the coarse-to-coarse register)
  //
  //     10. Advance level l+1 to the same time as level l
  //
  //     11. Reflux into level l, increment redist register with refluxed portion
  //
  //     12. Redistribute both way into l+1 and l-1, correct level l
  //
  //     13. Average down level l from l+1
  //  }
  //
  // These steps require the following new signatures:
  //
  // * Compute fluxes on level l
  // * Level-wise conservative, nonconservative, and hybrid divergences
  // * Per-level redistribution
  // * Per-level flux registers
  // * Per-level mass redistribution registers
  // * Per-level refluxing, which includes redistribution registers
  // * Per-level redistribution
  // * Per-level averaging
  
  EBAMRCellData& state = m_solver->get_state();
  
  m_solver->compute_divJ(m_k1, state, 0.0);  
  data_ops::incr(state, m_k1, -a_dt);        
  m_solver->compute_divJ(m_k2, state, 0.0);  
  data_ops::incr(state, m_k1,  0.5*a_dt);    
  data_ops::incr(state, m_k2, -0.5*a_dt);

  m_solver->make_non_negative(state);
  
  m_amr->average_down(state, phase::gas);
  m_amr->interp_ghost(state, phase::gas);
  
  return a_dt;
}

void euler_subcycle::regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level){

  // Regrid CDR solver
  m_solver->regrid(a_lmin, a_old_finest_level, a_new_finest_level);
  m_solver->set_source(0.0);
  m_solver->set_ebflux(0.0);
  if(m_solver->is_diffusive()){
    m_solver->set_diffco(m_diffco);
  }
  if(m_solver->is_mobile()){
    this->set_velocity();
  }

  // Allocate memory for RK steps
  m_amr->allocate(m_k1,  phase::gas, 1);
  m_amr->allocate(m_k2,  phase::gas, 1);
}
