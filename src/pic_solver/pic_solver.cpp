/*!
  @file   pic_solver.cpp
  @brief  Implementation of pic_solver.H
  @author Robert Marskar
  @date   June 2019
*/

#include "pic_solver.H"
#include "data_ops.H"
#include "units.H"

#include <ParmParse.H>
#include <EBAlias.H>
#include <EBLevelDataOps.H>
#include <BoxIterator.H>
#include <EBArith.H>
#include <ParmParse.H>
#include <ParticleIO.H>
#include <EBAMRIO.H>

  
pic_solver::pic_solver(){
  m_name = "pic_solver";
  set_verbosity(10);
  set_time(0, 0.0, 0.0);
  set_deposition_type();
  set_phase(phase::gas);
  set_charge(0);
  set_mass(1.0);
}
pic_solver::~pic_solver(){

}

void pic_solver::set_amr(const RefCountedPtr<amr_mesh>& a_amr){
  CH_TIME("pic_solver::set_amr");
  if(m_verbosity > 5){
    pout() << m_name + "::set_amr" << endl;
  }
  
  m_amr = a_amr;
}

void pic_solver::set_computational_geometry(const RefCountedPtr<computational_geometry> a_compgeom){
  CH_TIME("pic_solver::set_computational_geometry");
  if(m_verbosity > 5){
    pout() << m_name + "::set_computational_geometry" << endl;
  }
  m_compgeom = a_compgeom;
}

void pic_solver::set_physical_domain(const RefCountedPtr<physical_domain>& a_physdom){
  CH_TIME("pic_solver::set_physical_domain");
  if(m_verbosity > 5){
    pout() << m_name + "::set_physical_domain" << endl;
  }

  m_physdom = a_physdom;
}

void pic_solver::set_verbosity(const int a_verbosity){
  m_verbosity = a_verbosity;
}

void pic_solver::set_charge(const int a_charge){
  m_charge = a_charge;
}

void pic_solver::set_mass(const Real a_mass){
  m_mass = a_mass;
}

void pic_solver::advance(const Real a_dt){

}

void pic_solver::deposit_charge(EBAMRCellData& a_rho) {
  CH_TIME("pic_solver::deposit_charge");
  if(m_verbosity > 5){
    pout() << m_name + "::deposit_charge" << endl;
  }

  // 1. Deposit particle NUMBERS onto scratch storage
  data_ops::set_value(m_scratch_sca, 0.0);
  
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    const Real dx                = m_amr->get_dx()[lvl];
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const RealVect origin  = m_physdom->get_prob_lo();
    
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const Box box          = dbl.get(dit());
      MeshInterp interp(box, dx*RealVect::Unit, origin);
      interp.deposit((*m_particles[lvl])[dit()].listItems(), (*m_scratch_sca[lvl])[dit()].getFArrayBox(), m_deposition);
    }

    // Add in the contribution from the ghost cells
    const RefCountedPtr<Copier>& reversecopier = m_amr->get_reverse_copier(m_phase)[lvl];
    LDaddOp<FArrayBox> addOp;
    LevelData<FArrayBox> aliasFAB;
    aliasEB(aliasFAB, *m_scratch_sca[lvl]);
    aliasFAB.exchange(Interval(0,0), *reversecopier, addOp);
  }

  // 2. Scale with particle charge
  //  data_ops::scale(m_scratch_sca, m_charge*units::s_Qe);

  // 3. Increment space charge density
  data_ops::incr(a_rho, m_scratch_sca, 1.0);
}

void pic_solver::interpolate_force(const EBAMRCellData& a_E){
  CH_TIME("pic_solver::interpolate_force");
  if(m_verbosity > 5){
    pout() << m_name + "::interpolate_force" << endl;
  }

  // 1. Make a = q*E/m
  data_ops::copy(m_scratch_vec, a_E);
  data_ops::scale(m_scratch_vec, m_charge*units::s_Qe/m_mass);

  // 2. Interpolate a = q*E/m onto particle positions
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    const Real dx                = m_amr->get_dx()[lvl];
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      MeshInterp interp(dbl.get(dit()), m_amr->get_dx()[lvl]*RealVect::Unit, m_physdom->get_prob_lo());
      interp.interpolate((*m_particles[lvl])[dit()].listItems(), (*m_scratch_vec[lvl])[dit()].getFArrayBox(), m_deposition);
    }
  }
}

void pic_solver::initial_data(List<Particle>& a_particles){
  CH_TIME("pic_solver::initial_data");
  if(m_verbosity > 5){
    pout() << m_name + "::initial_data" << endl;
  }

  const RealVect origin  = m_physdom->get_prob_lo();
  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const Real dx                = m_amr->get_dx()[lvl];
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const ProblemDomain& dom     = m_amr->get_domains()[lvl];
    const bool coarsest_level    = (lvl == 0);

    if(coarsest_level){ // On coarsest level, destructively move all particles into their respective boxes
      m_particles[lvl]->clear();
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	ListBox<Particle>& particles = (*m_particles[lvl])[dit()];

	particles.addItemsDestructive(a_particles);
      }
    }
    else { // If there is a coarser level, move particles from the coarser level and onto this one if they overlap with this PVR
      const int ref_ratio = m_amr->get_ref_rat()[lvl-1];
      collectValidParticles(m_particles[lvl]->outcast(),
      			    *m_particles[lvl-1],
      			    m_pvr[lvl]->mask(),
      			    dx*RealVect::Unit,
      			    ref_ratio,
			    false,
			    origin);
      m_particles[lvl]->remapOutcast();
    }

#if 0
    if(!(m_particles[lvl]->outcast().isEmpty())){
      MayDay::Abort("pic_solver::initial_data - shouldn't happen!");
    }
    m_particles[lvl]->outcast().clear();      
#endif
  }


#if 0  // Debug
  std::cout << m_particles[0]->numParticles() << std::endl;
#endif
}

void pic_solver::set_deposition_type(){
  CH_TIME("pic_solver::set_deposition_type");
  if(m_verbosity > 5){
    pout() << m_name + "::set_deposition_type" << endl;
  }

  std::string str;
  ParmParse pp("pic_solver");
  pp.get("deposition", str);

  if(str == "ngp"){
    m_deposition = InterpType::NGP;
  }
  else if(str == "cic"){
    m_deposition = InterpType::CIC;
  }
  else if(str == "tsc"){
    m_deposition = InterpType::TSC;
  }
  else if(str == "w4"){
    m_deposition = InterpType::W4;
  }
  else{
    MayDay::Abort("pic_solver::set_deposition_type - unknown interpolant requested");
  }
}

void pic_solver::insert_particles(EBAMRParticles& a_particles){
  CH_TIME("pic_solver::insert_particles");
  if(m_verbosity > 5){
    pout() << m_name + "::insert_particles" << endl;
  }
}

void pic_solver::allocate_internals(){
  CH_TIME("pic_solver::allocate_internals");
  if(m_verbosity > 5){
    pout() << m_name + "::allocate_internals" << endl;
  }

  const int buffer = 0;
  const int ncomp  = 1;

  m_amr->allocate(m_scratch_sca, m_phase, ncomp);
  m_amr->allocate(m_scratch_vec, m_phase, SpaceDim);

  m_amr->allocate(m_particles);
  m_amr->allocate(m_pvr, buffer);
}

void pic_solver::cache_state(){

}

void pic_solver::set_time(const int a_step, const Real a_time, const Real a_dt){
  m_step = a_step;
  m_time = a_time;
  m_dt   = a_dt;
}

void pic_solver::deallocate_internals(){

}

void pic_solver::regrid(const int a_old_finest_level, const int a_new_finest_level){

}

void pic_solver::write_plot_file() {
  CH_TIME("pic_solver::write_plot_file");
  if(m_verbosity > 5){
    pout() << m_name + "::write_plot_file" << endl;
  }

  // Filename
  char file_char[1000];
  sprintf(file_char, "%s.step%07d.%dd.hdf5", m_name.c_str(), m_step, SpaceDim);

  // Deposit charge to some data holder
  const int ncomps = 1;
  EBAMRCellData output;
  Vector<string> names(ncomps);
  m_amr->allocate(output, m_phase, ncomps);
  data_ops::set_value(output, 0.0);
  deposit_charge(output);
  names[0] = "density";

  Vector<Real> covered_values(ncomps, 0.0);
  string fname(file_char);

  Vector<LevelData<EBCellFAB>* > output_ptr;
  m_amr->alias(output_ptr, output);

  writeEBHDF5(fname,
	      m_amr->get_grids(),
	      output_ptr,
	      names,
	      m_amr->get_domains()[0].domainBox(),
	      m_amr->get_dx()[0],
	      m_dt,
	      m_time,
	      m_amr->get_ref_rat(),
	      m_amr->get_finest_level() + 1,
	      true,
	      covered_values,
	      IntVect::Unit);
}

void pic_solver::set_phase(const phase::which_phase a_phase){
  CH_TIME("pic_solver::set_phase");
  if(m_verbosity > 5){
    pout() << m_name + "::pic_phase" << endl;
  }

  m_phase = a_phase;
}

void pic_solver::leapfrog_rewind(const Real a_dt){
  CH_TIME("pic_solver::leapfrog_rewind");
  if(m_verbosity > 5){
    pout() << m_name + "leapfrog_rewind" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    const Real dx                = m_amr->get_dx()[lvl];
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const RealVect origin  = m_physdom->get_prob_lo();
    
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      List<Particle>& particles = (*m_particles[lvl])[dit()].listItems();

      for (ListIterator<Particle> lit(particles); lit.ok(); ++lit){
	Particle& particle = lit();

	particle.velocity() -= 0.5*particle.acceleration()*a_dt;
      }
    }
  }
}

void pic_solver::advance_leapfrog(const Real a_dt){
  CH_TIME("pic_solver::leapfrog_rewind");
  if(m_verbosity > 5){
    pout() << m_name + "leapfrog_rewind" << endl;
  }

  // 1. Move and accelerate
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    const Real dx                = m_amr->get_dx()[lvl];
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const RealVect origin  = m_physdom->get_prob_lo();
    
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      List<Particle>& particles = (*m_particles[lvl])[dit()].listItems();

      for (ListIterator<Particle> lit(particles); lit.ok(); ++lit){
	Particle& particle = lit();

	
	particle.velocity() += particle.acceleration()*a_dt;
	particle.position() += particle.velocity()*a_dt;
      }
    }
  }

  // 2. Rebin photons
  MayDay::Abort("pic_solver::advance_leapfrog - routine is not finished");
}

Real pic_solver::get_time() const{
  return m_time;
}
  
int pic_solver::get_step() const{
  return m_step;
}

EBAMRParticles& pic_solver::get_particles(){
  return m_particles;
}

EBAMRParticles& pic_solver::get_eb_particles(){
  MayDay::Abort("pic_solver::get_eb_particles - not returning the correct particles. Routine is not done.");
  return m_particles;
}

EBAMRParticles& pic_solver::get_domain_particles(){
  MayDay::Abort("pic_solver::get_domain_particles - not returning the correct particles. Routine is not done.");
  return m_particles;
}
EBAMRParticles& pic_solver::get_coll_particles(){
  MayDay::Abort("pic_solver::get_coll_particles - not returning the correct particles. Routine is not done.");
  return m_particles;
}
