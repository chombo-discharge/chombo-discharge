/*!
  @file   ito_solver.H
  @brief  Declaration of an abstract class for Ito diffusion
  @author Robert Marskar
  @date   April 2020
  @todo   Compute CFL time step
  @todo   Compute diffusion time step, do this based on chance that a particle will move more than 1 cell
  @todo   Boundary conditions on EBs
  @todo   Boundary conditions on domain edges/sides
*/

#include "ito_solver.H"
#include "data_ops.H"
#include "EBAlias.H"

#include <ParmParse.H>

#include <chrono>

ito_solver::ito_solver(){
  m_name       = "ito_solver";
  m_class_name = "ito_solver";
}

ito_solver::~ito_solver(){

}

std::string ito_solver::get_name(){
  return m_name;
}

void ito_solver::parse_options(){
  CH_TIME("ito_solver::parse_options");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_options" << endl;
  }

  this->parse_rng();
  this->parse_plot_vars();
  this->parse_deposition();
  this->parse_bisect_step();
  this->parse_pvr_buffer();
}

void ito_solver::parse_rng(){
  CH_TIME("ito_solver::parse_rng");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_rng" << endl;
  }

    // Seed the RNG
  ParmParse pp(m_class_name.c_str());
  pp.get("seed", m_seed_rng);
  if(m_seed_rng < 0) { // Random seed if input < 0
    m_seed_rng = std::chrono::system_clock::now().time_since_epoch().count();
  }
  
  m_rng = std::mt19937_64(m_seed_rng);

  m_udist01 = std::uniform_real_distribution<Real>( 0.0, 1.0);
  m_udist11 = std::uniform_real_distribution<Real>(-1.0, 1.0);
  m_gauss01 = std::normal_distribution<Real>(0.0, 1.0);
}

void ito_solver::parse_plot_vars(){
  CH_TIME("mc_photo::parse_plot_vars");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_plot_vars" << endl;
  }

  m_plot_phi = false;
  m_plot_vel = false;
  m_plot_dco = false;

  ParmParse pp(m_class_name.c_str());
  const int num = pp.countval("plt_vars");
  Vector<std::string> str(num);
  pp.getarr("plt_vars", str, 0, num);

  for (int i = 0; i < num; i++){
    if(     str[i] == "phi") m_plot_phi = true;
    else if(str[i] == "vel") m_plot_vel = true;
    else if(str[i] == "dco") m_plot_dco = true;
  }
}

void ito_solver::parse_deposition(){
  CH_TIME("ito_solver::parse_rng");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_rng" << endl;
  }

  ParmParse pp(m_class_name.c_str());
  std::string str;

  // Deposition for particle-mesh operations
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
    MayDay::Abort("mc_photo::set_deposition_type - unknown interpolant requested");
  }

  // Deposition for plotting only
  pp.get("plot_deposition", str);

  if(str == "ngp"){
    m_plot_deposition = InterpType::NGP;
  }
  else if(str == "cic"){
    m_plot_deposition = InterpType::CIC;
  }
  else if(str == "tsc"){
    m_plot_deposition = InterpType::TSC;
  }
  else if(str == "w4"){
    m_plot_deposition = InterpType::W4;
  }
  else{
    MayDay::Abort("mc_photo::set_deposition_type - unknown interpolant requested");
  }
}

void ito_solver::parse_bisect_step(){
  CH_TIME("ito_solver::parse_bisect_step");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_bisect_step" << endl;
  }

  ParmParse pp(m_class_name.c_str());
  pp.get("bisect_step", m_bisect_step);
}

void ito_solver::parse_pvr_buffer(){
  CH_TIME("ito_solver::parse_pvr_buffer");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_pvr_buffer" << endl;
  }

  ParmParse pp(m_class_name.c_str());
  pp.get("pvr_buffer", m_pvr_buffer);
}

Vector<std::string> ito_solver::get_plotvar_names() const {
  CH_TIME("ito_solver::get_plotvar_names");
  if(m_verbosity > 5){
    pout() << m_name + "::get_plotvar_names" << endl;
  }

  Vector<std::string> names(0);
  if(m_plot_phi) {
    names.push_back(m_name + " phi");
  }
  if(m_plot_dco && m_diffusive){
    names.push_back(m_name + " diffusion_coefficient");
  }
  if(m_plot_vel && m_mobile){
    names.push_back("x-Velocity " + m_name);
    names.push_back("y-Velocity " + m_name);
    if(SpaceDim == 3){
      names.push_back("z-Velocity " + m_name);
    }
  }

  return names;
}

int ito_solver::get_num_plotvars() const {
  CH_TIME("ito_solver::get_num_plotvars");
  if(m_verbosity > 5){
    pout() << m_name + "::get_num_plotvars" << endl;
  }

  int num_plotvars = 0;
  
  if(m_plot_phi)                num_plotvars += 1;
  if(m_plot_dco && m_diffusive) num_plotvars += 1;
  if(m_plot_vel && m_mobile)    num_plotvars += SpaceDim;

  return num_plotvars;
}

void ito_solver::set_computational_geometry(const RefCountedPtr<computational_geometry> a_compgeom){
  CH_TIME("ito_solver::set_computational_geometry");
  if(m_verbosity > 5){
    pout() << m_name + "::set_computational_geometry" << endl;
  }
  m_compgeom = a_compgeom;
}

void ito_solver::set_amr(const RefCountedPtr<amr_mesh>& a_amr){
  CH_TIME("ito_solver::set_amr");
  if(m_verbosity > 5){
    pout() << m_name + "::set_amr" << endl;
  }

  m_amr = a_amr;
}

void ito_solver::set_phase(phase::which_phase a_phase){
  CH_TIME("ito_solver::set_phase");
  if(m_verbosity > 5){
    pout() << m_name + "::set_phase" << endl;
  }

  m_phase = a_phase;
}

void ito_solver::set_verbosity(const int a_verbosity){
  CH_TIME("ito_solver::set_verbosity");
  m_verbosity = a_verbosity;
  if(m_verbosity > 5){
    pout() << m_name + "::set_verbosity" << endl;
  }
}

void ito_solver::set_time(const int a_step, const Real a_time, const Real a_dt) {
  CH_TIME("ito_solver::set_time");
  if(m_verbosity > 5){
    pout() << m_name + "::set_time" << endl;
  }

  m_step = a_step;
  m_time = a_time;
  m_dt   = a_dt;
}

void ito_solver::initial_data(){
  CH_TIME("ito_solver::initial_data");
  if(m_verbosity > 5){
    pout() << m_name + "::initial_data" << endl;
  }

  // Put the initial particles on the coarsest grid level
  List<ito_particle>& outcastBase = m_particles[0]->outcast();
  outcastBase.catenate(m_species->get_initial_particles()); // This destroys the initial partcies
  m_particles[0]->remapOutcast(); 

  // Move particles to finer levels if they belong there. This piece of code moves particles from lvl-1
  // and onto the outcast list on level lvl. Then, we remap the outcast list
  for (int lvl = 1; lvl <= m_amr->get_finest_level(); lvl++){
    collectValidParticles(m_particles[lvl]->outcast(),
			  *m_particles[lvl-1],
			  m_pvr[lvl]->mask(),
			  m_amr->get_dx()[lvl]*RealVect::Unit,
			  m_amr->get_ref_rat()[lvl-1],
			  false,
			  m_amr->get_prob_lo());
    m_particles[lvl]->remapOutcast();


  }

  // Deposit the particles
  this->deposit_particles(m_state, m_particles, m_deposition);
  
#if 0 // Experimental code. Compute the number of particles in each box
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    Vector<Box> oldBoxes(0);
    Vector<int> oldLoads(0);
    int numParticlesThisProc = 0;
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const Box thisBox = dbl.get(dit());
      const int numPart = (*m_particles[lvl])[dit()].numItems();

      oldBoxes.push_back(thisBox);
      oldLoads.push_back(numPart);

      numParticlesThisProc += numPart;
    }
    pout() << "on level = " << lvl << "\t num particles = " << numParticlesThisProc << endl;


    // Gather boxes and loads
    load_balance::gather_boxes(oldBoxes);
    load_balance::gather_loads(oldLoads);

    Vector<Box> newBoxes = oldBoxes;
    Vector<int> newLoads = oldLoads;

    mortonOrdering(oldBoxes);
    for (int i = 0; i < newBoxes.size(); i++){
      for (int j = 0; j < oldBoxes.size(); j++){
	if(newBoxes[i] == oldBoxes[j]){
	  newLoads[i] = newLoads[j];
	}
      }
    }

    Vector<int> procs;
    load_balance::load_balance_boxes(procs, newLoads, newBoxes);

    int sum = 0;
    for (int i = 0; i < oldLoads.size(); i++){
      sum += oldLoads[i];
    }

    if(procID() == 0){
      std::cout << "num particles on this level = " << sum << std::endl;
      std::cout << "procs" << std::endl;
    }


  }
#endif
}

void ito_solver::regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level){
  CH_TIME("ito_solver::regrid");
  if(m_verbosity > 5){
    pout() << m_name + "::regrid" << endl;
  }

  // This reallocates all internal stuff. After this, the only object with any knowledge of the past
  // is m_particleCache, which has all the photons in the old data holders.
  this->allocate_internals();

  // Here are the steps:
  // 
  // 1. We are regridding a_lmin and above, so we need to move all photons onto level a_lmin-1. There's
  //    probably a better way to do this
  // 
  // 2. Levels 0 through a_lmin-1 did not change. Copy photons back into those levels. 
  // 
  // 3. Levels a_lmin-1 through a_new_finest_level changed. Remap photons into the correct boxes. 
  
  const int base_level  = Max(0, a_lmin-1);
  const RealVect origin = m_amr->get_prob_lo();

  // 1. Move all particles from cache and onto the coarsest level
  List<ito_particle>& base_cache = m_particleCache[0]->outcast();
  base_cache.clear();
  for (int lvl = base_level; lvl <= a_old_finest_level; lvl++){
    for (DataIterator dit = m_particleCache[lvl]->dataIterator(); dit.ok(); ++dit){
      for (ListIterator<ito_particle> li((*m_particleCache[lvl])[dit()].listItems()); li.ok(); ){
	m_particleCache[base_level]->outcast().transfer(li);
      }
    }
  }
    

  // 2. Levels 0 through base_level did not change, so copy photons back onto those levels. 
  for (int lvl = 0; lvl <= base_level; lvl++){
    collectValidParticles(*m_particles[lvl],
			  *m_particleCache[lvl],
			  m_pvr[lvl]->mask(),
			  m_amr->get_dx()[lvl]*RealVect::Unit,
			  1,
			  false, 
			  origin);

    m_particles[lvl]->gatherOutcast();
    m_particles[lvl]->remapOutcast();
  }

  // 3. Levels above base_level may have changed, so we need to move the particles into the correct boxes now
  for (int lvl = base_level; lvl < a_new_finest_level; lvl++){
    collectValidParticles(m_particles[lvl+1]->outcast(),
			  *m_particles[lvl],
			  m_pvr[lvl+1]->mask(),
			  m_amr->get_dx()[lvl+1]*RealVect::Unit,
			  m_amr->get_ref_rat()[lvl],
			  false,
			  origin);
    m_particles[lvl+1]->remapOutcast();
  }
}

void ito_solver::set_species(RefCountedPtr<ito_species>& a_species){
  CH_TIME("ito_solver::set_species");
  if(m_verbosity > 5){
    pout() << m_name + "::set_species" << endl;
  }
  
  m_species   = a_species;
  m_name      = a_species->get_name();
  m_diffusive = m_species->is_diffusive();
  m_mobile    = m_species->is_mobile();
}

void ito_solver::allocate_internals(){
  CH_TIME("ito_solver::allocate_internals");
  if(m_verbosity > 5){
    pout() << m_name + "::allocate_internals" << endl;
  }
  const int ncomp = 1;

  m_amr->allocate(m_state,       m_phase, ncomp);
  m_amr->allocate(m_scratch,     m_phase, ncomp);

  data_ops::set_value(m_state,   0.0);
  data_ops::set_value(m_scratch, 0.0);

  // Only allocate memory if we actually have a mobile solver
  if(m_mobile){
    m_amr->allocate(m_velo_cell, m_phase, SpaceDim);
  }
  else{ 
    m_amr->allocate_ptr(m_velo_cell);
  }

  // Only allocate memory if we actually a diffusion solver
  if(m_diffusive){
    m_amr->allocate(m_diffco_cell, m_phase, 1);
  }
  else{
    m_amr->allocate_ptr(m_diffco_cell);
  }
  
  // This allocates parallel data holders using the load balancing in amr_mesh. This will give very poor
  // load balancing, but we will rectify that by rebalancing later. 
  m_amr->allocate(m_particles);
  m_amr->allocate(m_pvr, m_pvr_buffer);
}


void ito_solver::write_checkpoint_level(HDF5Handle& a_handle, const int a_level) const {
  CH_TIME("ito_solver::write_checkpoint_level");
  if(m_verbosity > 5){
    pout() << m_name + "::write_checkpoint_level" << endl;
  }

  MayDay::Abort("ito_solver::write_checkpoint_level - checkpointing not implemented");
}

void ito_solver::read_checkpoint_level(HDF5Handle& a_handle, const int a_level){
  CH_TIME("ito_solver::read_checkpoint_level");
  if(m_verbosity > 5){
    pout() << m_name + "::read_checkpoint_level" << endl;
  }

  MayDay::Abort("ito_solver::read_checkpoint_level - checkpointing not implemented");
}

void ito_solver::write_plot_data(EBAMRCellData& a_output, int& a_comp){
  CH_TIME("ito_solver::write_plot_data");
  if(m_verbosity > 5){
    pout() << m_name + "::write_plot_data" << endl;
  }

  // Write phi
  if(m_plot_phi){
    this->deposit_particles(m_state, m_particles, m_plot_deposition);
    
    const Interval src(0, 0);
    const Interval dst(a_comp, a_comp);
    
    for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
      m_state[lvl]->localCopyTo(src, *a_output[lvl], dst);
    }
    data_ops::set_covered_value(a_output, 0.0, a_comp);
    a_comp++;
  }

  // Plot diffusion coefficient
  if(m_plot_dco && m_diffusive){
    
    const Interval src(0, 0);
    const Interval dst(a_comp, a_comp);
    
    for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
      m_diffco_cell[lvl]->localCopyTo(src, *a_output[lvl], dst);
    }
    data_ops::set_covered_value(a_output, 0.0, a_comp);
    a_comp++;
  }

  // Write velocities
  if(m_plot_vel && m_mobile){
    const int ncomp = SpaceDim;
    const Interval src(0, ncomp-1);
    const Interval dst(a_comp, a_comp + ncomp-1);

    for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
      m_velo_cell[lvl]->localCopyTo(src, *a_output[lvl], dst);
    }

    for (int c = 0; c < SpaceDim; c++){
      data_ops::set_covered_value(a_output, 0.0, a_comp + c);
    }

    a_comp += ncomp;
  }
}

void ito_solver::deposit_particles(){
  CH_TIME("ito_solver::deposit_particles");
  if(m_verbosity > 5){
    pout() << m_name + "::deposit_particles" << endl;
  }

  this->deposit_particles(m_state, m_particles, m_deposition);
}

void ito_solver::deposit_particles(EBAMRCellData&           a_state,
				   const EBAMRItoParticles& a_particles,
				   const InterpType&        a_deposition){
  CH_TIME("ito_solver::deposit_particles");
  if(m_verbosity > 5){
    pout() << m_name + "::deposit_particles" << endl;
  }

  // TODO: How should we deposit particles near the EB? We should not
  //
  //          1. Deposit into covered cells since this is equivalent to losing mass
  //          2. Deposit into cells that can't be reached through the EBGraph of distance 1.
  //
  //       Should we just take the cut-cell particles and deposit them with an NGP scheme...? 
  MayDay::Warning("ito_solver::deposit_particles - not EB supported yet");
           

  // TLDR: This code deposits on the entire AMR mesh. For all levels l > 0 the data on the coarser grids are interpolated
  //       to the current grid first. Then we deposit.
  //
  //       We always assume that a_state is a scalar, i.e. a_state has only one component. 

  const int comp = 0;
  const Interval interv(comp, comp);

  const RealVect origin  = m_amr->get_prob_lo();
  const int finest_level = m_amr->get_finest_level();

  data_ops::set_value(a_state,    0.0);
  data_ops::set_value(m_scratch,  0.0);

  InterpType deposition = a_deposition;

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const Real dx                = m_amr->get_dx()[lvl];
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const ProblemDomain& dom     = m_amr->get_domains()[lvl];

    const bool has_coar = (lvl > 0);
    const bool has_fine = (lvl < finest_level);

    // 1. If we have a coarser level whose cloud extends beneath this level, interpolate that result here first. 
    if(has_coar){
      RefCountedPtr<EBPWLFineInterp>& interp = m_amr->get_eb_pwl_interp(m_phase)[lvl];
      interp->interpolate(*a_state[lvl], *m_scratch[lvl-1], interv);
    }
    
    // 2. Deposit this levels particles and exchange ghost cells
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const Box box          = dbl.get(dit());
      MeshInterp interp(box, dx*RealVect::Unit, origin);
      interp.deposit((*a_particles[lvl])[dit()].listItems(), (*a_state[lvl])[dit()].getFArrayBox(), deposition);
    }

    // Exchange ghost cells. 
    const RefCountedPtr<Copier>& reversecopier = m_amr->get_reverse_copier(m_phase)[lvl];
    LDaddOp<FArrayBox> addOp;
    LevelData<FArrayBox> aliasFAB;
    aliasEB(aliasFAB, *a_state[lvl]);
    aliasFAB.exchange(Interval(0,0), *reversecopier, addOp);

    // 3. If we have a finer level, copy contributions from this level to the temporary holder that is used for
    //    interpolation of "hanging clouds"
    if(has_fine){
      a_state[lvl]->copyTo(*m_scratch[lvl]);
    }
  }

  // Do a kappa scaling
  data_ops::kappa_scale(a_state);
}

bool ito_solver::is_mobile() const{
  CH_TIME("ito_solver::is_mobile");
  if(m_verbosity > 5){
    pout() << m_name + "::is_mobile" << endl;
  }
  
  return m_mobile;
}
  

bool ito_solver::is_diffusive() const{
  CH_TIME("ito_solver::is_diffusive");
  if(m_verbosity > 5){
    pout() << m_name + "::is_diffusive" << endl;
  }
  return m_diffusive;
}

void ito_solver::cache_state(){
  CH_TIME("ito_solver::cache_state");
  if(m_verbosity > 5){
    pout() << m_name + "::cache_state" << endl;
  }

  m_amr->allocate(m_particleCache);

  // Copy particles onto cache
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    collectValidParticles(m_particleCache[lvl]->outcast(),
			  *m_particles[lvl],
			  m_pvr[lvl]->mask(),
			  m_amr->get_dx()[lvl]*RealVect::Unit,
			  1,
			  false, 
			  m_amr->get_prob_lo());
    m_particleCache[lvl]->remapOutcast();
  }
}

EBAMRItoParticles& ito_solver::get_particles(){
  CH_TIME("ito_solver::get_particles");
  if(m_verbosity > 5){
    pout() << m_name + "::get_particles" << endl;
  }

  return m_particles;
}

EBAMRCellData& ito_solver::get_velo_cell(){
  CH_TIME("ito_solver::get_velo_cell");
  if(m_verbosity > 5){
    pout() << m_name + "::get_velo_cell" << endl;
  }

  return m_velo_cell;
}

EBAMRCellData& ito_solver::get_diffco_cell(){
  CH_TIME("ito_solver::get_diffco_cell");
  if(m_verbosity > 5){
    pout() << m_name + "::get_diffco_cell" << endl;
  }

  return m_diffco_cell;
}

void ito_solver::set_diffco(const Real a_diffco){
  CH_TIME("ito_solver::set_diffco");
  if(m_verbosity > 5){
    pout() << m_name + "::set_diffco" << endl;
  }

  data_ops::set_value(m_diffco_cell, a_diffco);
}

void ito_solver::set_velocity(const RealVect a_vel){
  CH_TIME("ito_solver::set_velocity");
  if(m_verbosity > 5){
    pout() << m_name + "::set_velocity" << endl;
  }

  for (int comp = 0; comp < SpaceDim; comp++){
    data_ops::set_value(m_velo_cell, a_vel[comp], comp);
  }
}

void ito_solver::interpolate_velocities(){
  CH_TIME("ito_solver::interpolate_velocities");
  if(m_verbosity > 5){
    pout() << m_name + "::interpolate_velocities" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      this->interpolate_velocities(lvl, dit());
    }
  }
}

void ito_solver::interpolate_velocities(const int a_lvl, const DataIndex& a_dit){
  CH_TIME("ito_solver::interpolate_velocities");
  if(m_verbosity > 5){
    pout() << m_name + "::interpolate_velocities" << endl;
  }

  const EBCellFAB& velo_cell = (*m_velo_cell[a_lvl])[a_dit];
  const FArrayBox& vel_fab   = velo_cell.getFArrayBox();
  const RealVect dx          = m_amr->get_dx()[a_lvl]*RealVect::Unit;
  const RealVect origin      = m_amr->get_prob_lo();
  const Box box              = m_amr->get_grids()[a_lvl][a_dit];

  List<ito_particle>& particleList = (*m_particles[a_lvl])[a_dit].listItems();

  MeshInterp meshInterp(box, dx, origin);
  meshInterp.interpolate(particleList, vel_fab, m_deposition);
}

void ito_solver::move_particles_eulerf(const Real a_dt){
  CH_TIME("ito_solver::move_particles_eulerf");
  if(m_verbosity > 5){
    pout() << m_name + "::move_particles_eulerf" << endl;
  }

  // Move all particles and rebin them
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      this->move_particles_eulerf(lvl, dit(), a_dt);
    }

    m_particles[lvl]->gatherOutcast();
    m_particles[lvl]->remapOutcast();
  }
}

void ito_solver::move_particles_eulerf(const int a_lvl, const DataIndex& a_dit, const Real a_dt){
  CH_TIME("ito_solver::move_particles_eulerf");
  if(m_verbosity > 5){
    pout() << m_name + "::move_particles_eulerf" << endl;
  }

  List<ito_particle>& particleList = (*m_particles[a_lvl])[a_dit].listItems();
  ListIterator<ito_particle> lit(particleList);
  
  for (lit.rewind(); lit; ++lit){
    ito_particle& p = particleList[lit];

    p.position() += p.velocity()*a_dt;
  }
}
