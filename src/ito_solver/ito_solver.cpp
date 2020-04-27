/*!
  @file   ito_solver.cpp
  @brief  Declaration of an abstract class for Ito diffusion
  @author Robert Marskar
  @date   April 2020
*/

#include "ito_solver.H"
#include "data_ops.H"
#include "EBParticleInterp.H"
#include "units.H"

#include <ParmParse.H>
#include <EBAlias.H>
#include <BaseEBCellFactory.H>

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
  this->parse_diffusion_hop();
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
    m_deposition = DepositionType::NGP;
  }
  else if(str == "cic"){
    m_deposition = DepositionType::CIC;
  }
  else if(str == "tsc"){
    m_deposition = DepositionType::TSC;
  }
  else if(str == "w4"){
    m_deposition = DepositionType::W4;
  }
  else{
    MayDay::Abort("mc_photo::set_deposition_type - unknown interpolant requested");
  }

  // Deposition for plotting only
  pp.get("plot_deposition", str);

  if(str == "ngp"){
    m_plot_deposition = DepositionType::NGP;
  }
  else if(str == "cic"){
    m_plot_deposition = DepositionType::CIC;
  }
  else if(str == "tsc"){
    m_plot_deposition = DepositionType::TSC;
  }
  else if(str == "w4"){
    m_plot_deposition = DepositionType::W4;
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

void ito_solver::parse_diffusion_hop(){
  CH_TIME("ito_solver::parse_diffusion_hop");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_diffusion_hop" << endl;
  }

  ParmParse pp(m_class_name.c_str());

  pp.get("max_diffusion_hop", m_max_diffusion_hop);
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

  // Get initial particles
  m_particles.add_particles(m_species->get_initial_particles());

  // Remove particles that are inside embedded boundaries
#if 1
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[lvl];

    RefCountedPtr<BaseIF> func;
    if(m_phase == phase::gas){
      func = m_compgeom->get_gas_if();
    }
    else{
      func = m_compgeom->get_sol_if();
    }
    
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const EBISBox& ebisbox = ebisl[dit()];
      
      List<ito_particle>& particles = m_particles[lvl][dit()].listItems();

      if(ebisbox.isAllCovered()){ // Box is all covered
	particles.clear();
      }
      else if(ebisbox.isAllRegular()){ // Do nothing
      }
      else{
	List<ito_particle>  particleCopy = List<ito_particle>(particles);
	particles.clear();
	
	ListIterator<ito_particle> lit(particleCopy);
	for (lit.rewind(); lit; ++lit){
	  ito_particle& p = particles[lit];

	  const Real f = func->value(p.position());
	  if(f <= 0.0){
	    particles.add(p);
	  }
	}
      }
    }
  }
#endif

  // Deposit the particles
  this->deposit_particles(m_state, m_particles.get_particles(), m_deposition);
  
#if 0 // Experimental code. Compute the number of particles in each box
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    Vector<Box> oldBoxes(0);
    Vector<int> oldLoads(0);
    int numParticlesThisProc = 0;
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const Box thisBox = dbl.get(dit());
      const int numPart = m_particles[lvl][dit()].numItems();

      oldBoxes.push_back(thisBox);
      oldLoads.push_back(numPart);

      numParticlesThisProc += numPart;
    }
    //    pout() << "on level = " << lvl << "\t num particles = " << numParticlesThisProc << endl;


    // Gather boxes and loads
    load_balance::gather_boxes(oldBoxes);
    load_balance::gather_loads(oldLoads);


    Vector<Box> newBoxes = oldBoxes;
    Vector<int> newLoads = oldLoads;

    mortonOrdering(oldBoxes);
    for (int i = 0; i < newBoxes.size(); i++){
      for (int j = 0; j < oldBoxes.size(); j++){
	if(newBoxes[i] == oldBoxes[j]){
	  newLoads[i] = oldLoads[j];
	}
      }
    }

    Vector<int> procs;
    load_balance::load_balance_boxes(procs, newLoads, newBoxes);



    if(procID() == 0) std::cout << "\nnew load balancing: " << std::endl;
    for (int i = 0; i < newLoads.size(); i++){
      if(procID() == 0){
	std::cout << "box = " << newBoxes[i] << "\t load = " << newLoads[i] << "\t proc = " << procs[i] << std::endl;
      }
    }
  }
#endif
}

void ito_solver::regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level){
  CH_TIME("ito_solver::regrid");
  if(m_verbosity > 5){
    pout() << m_name + "::regrid" << endl;
  }

  // Reallocate mesh data
  m_amr->reallocate(m_state,   m_phase, a_lmin);
  m_amr->reallocate(m_scratch, m_phase, a_lmin);

  // Only allocate memory if we actually have a mobile solver
  if(m_mobile){
    m_amr->reallocate(m_velo_cell, m_phase, a_lmin);
  }
  else{ 
    m_amr->allocate_ptr(m_velo_cell);
  }

  // Only allocate memory if we actually a diffusion solver
  if(m_diffusive){
    m_amr->reallocate(m_diffco_cell, m_phase, a_lmin);
  }
  else{
    m_amr->allocate_ptr(m_diffco_cell);
  }

  // Particle data regrids
  const Vector<DisjointBoxLayout>& grids = m_amr->get_grids();
  const Vector<ProblemDomain>& domains   = m_amr->get_domains();
  const Vector<Real>& dx                 = m_amr->get_dx();
  const Vector<int>& ref_rat             = m_amr->get_ref_rat();

  m_particles.regrid(       grids, domains, dx, ref_rat, a_lmin, a_new_finest_level);
  m_eb_particles.regrid(    grids, domains, dx, ref_rat, a_lmin, a_new_finest_level);
  m_domain_particles.regrid(grids, domains, dx, ref_rat, a_lmin, a_new_finest_level);
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
  
  // This allocates parallel data holders using the load balancing in amr_mesh. This might give poor
  // load balancing, but we will rectify that by rebalancing later.
  m_amr->allocate(m_particles,        m_pvr_buffer);
  m_amr->allocate(m_eb_particles,     m_pvr_buffer);
  m_amr->allocate(m_domain_particles, m_pvr_buffer);

#if 0 // Experimental code for getting the number of particles per cell
  Vector<LevelData<BaseEBCellFAB<int> >* > m_ppc(1 + m_amr->get_finest_level());
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const EBISLayout& ebisl = m_amr->get_ebisl(m_phase)[lvl];

    BaseEBCellFactory<int> fact(ebisl);
    m_ppc[lvl] = new LevelData<BaseEBCellFAB<int> >(dbl, 1, 3*IntVect::Unit, fact);
  }
#endif
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
    const Interval src(0, 0);
    const Interval dst(a_comp, a_comp);
    
    for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
      m_state[lvl]->localCopyTo(src, *a_output[lvl], dst);
    }
    //data_ops::set_covered_value(a_output, a_comp, 0.0);
    a_comp++;
  }

  // Plot diffusion coefficient
  if(m_plot_dco && m_diffusive){
    const Interval src(0, 0);
    const Interval dst(a_comp, a_comp);
    
    for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
      m_diffco_cell[lvl]->localCopyTo(src, *a_output[lvl], dst);
    }
    data_ops::set_covered_value(a_output, a_comp, 0.0);
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
      data_ops::set_covered_value(a_output, a_comp + c, 0.0);
    }

    a_comp += ncomp;
  }
}

void ito_solver::deposit_particles(){
  CH_TIME("ito_solver::deposit_particles");
  if(m_verbosity > 5){
    pout() << m_name + "::deposit_particles" << endl;
  }
  this->deposit_particles(m_state, m_particles.get_particles(), m_deposition);
}

void ito_solver::deposit_particles(EBAMRCellData&           a_state,
				   const AMRParticles<ito_particle>& a_particles,
				   const DepositionType::Which        a_deposition){
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
  //  MayDay::Warning("ito_solver::deposit_particles - not EB supported yet");
           

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

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const Real dx                = m_amr->get_dx()[lvl];
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const ProblemDomain& dom     = m_amr->get_domains()[lvl];
    const RefCountedPtr<EBLevelGrid>& eblg = m_amr->get_eblg(m_phase)[lvl];

    const bool has_coar = (lvl > 0);
    const bool has_fine = (lvl < finest_level);

    // 1. If we have a coarser level whose cloud extends beneath this level, interpolate that result here first. 
    if(has_coar){
      RefCountedPtr<EBPWLFineInterp>& interp = m_amr->get_eb_pwl_interp(m_phase)[lvl];
      interp->interpolate(*a_state[lvl], *m_scratch[lvl-1], interv);
    }
    
    // 2. Deposit this levels particles. Note that this will deposit into ghost cells, which must later
    //    be added to neighboring patches
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const Box box          = dbl.get(dit());
      EBParticleInterp interp(box, dx*RealVect::Unit, origin);
      interp.deposit((*a_particles[lvl])[dit()].listItems(), (*a_state[lvl])[dit()].getFArrayBox(), m_deposition);
    }

    // This code adds contributions from ghost cels
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

#if 1 // Development code
  for (int lvl = 0; lvl <= finest_level; lvl++){
    const Real dx                = m_amr->get_dx()[lvl];
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const ProblemDomain& dom     = m_amr->get_domains()[lvl];
    const EBLevelGrid& eblg      = *m_amr->get_eblg(m_phase)[lvl];

    const bool has_coar = (lvl > 0);
    const bool has_fine = (lvl < finest_level);
    
    if(has_coar){
      const DisjointBoxLayout& gridCoar = m_amr->get_grids()[lvl-1];
      DisjointBoxLayout gridCoFi;
      EBLevelGrid eblgCoFi;

      coarsen(gridCoFi, dbl, 2);
      coarsen(eblgCoFi, eblg, 2);
      
      Vector<Box> coFiBoxes = gridCoFi.boxArray();
      Vector<int> coFiProcs = gridCoFi.procIDs();

      for (int i = 0; i < coFiBoxes.size(); i++){
	coFiBoxes[i].grow(1);
	coFiBoxes[i] &= dom.domainBox();
      }

      BoxLayout coFiBl(coFiBoxes, coFiProcs);
      BoxLayoutData<FArrayBox> dataCoFi(coFiBl, 1);

      // Monkey with the scratch box
      a_state[lvl]->localCopyTo(*m_scratch[lvl]);
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	FArrayBox& scratchBox = (*m_scratch[lvl])[dit()].getFArrayBox();
	scratchBox.setVal(0.0, dbl.get(dit()), 0, 1);
      }
      m_scratch[lvl]->exchange();

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){

	
	Box fineBox = dbl.get(dit());
	Box bx = fineBox;

	bx.coarsen(2);
	
	fineBox.grow(2);
	Box coFiBox = fineBox;
	coFiBox.coarsen(2);

	FArrayBox& fabCoFi        = dataCoFi[dit()];
	const FArrayBox& fabFine = (*m_scratch[lvl])[dit()].getFArrayBox();
	fabCoFi.setVal(0.0);

	// Iterate over the CoFi box
	for (BoxIterator bit(fineBox); bit.ok(); ++bit){
	  IntVect ivFine = bit();
	  IntVect ivCoar = coarsen(ivFine, 2);

	  fabCoFi(ivCoar, 0) += fabFine(ivFine, 0)/4.;
	}

	fabCoFi.setVal(0.0, bx, 0, 1);
      }

      LevelData<FArrayBox> coarAlias;
      aliasEB(coarAlias, *a_state[lvl-1]);
      Interval interv(0,0);
      dataCoFi.addTo(interv, coarAlias, interv, dom.domainBox());
    }
#endif
  }

  // Do a kappa scaling
  data_ops::kappa_scale(a_state);

  m_amr->average_down(a_state, m_phase);
  m_amr->interp_ghost(a_state, m_phase);
}


void ito_solver::deposit_weights(EBAMRCellData& a_state, const AMRParticles<ito_particle>& a_particles){
  CH_TIME("ito_solver::is_mobile");
  if(m_verbosity > 5){
    pout() << m_name + "::is_mobile" << endl;
  }

  this->deposit_particles(a_state, a_particles, DepositionType::NGP);

  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    data_ops::scale(*a_state[lvl], pow(m_amr->get_dx()[lvl], SpaceDim));
  }
}

void ito_solver::add_particles(particle_container<ito_particle>& a_part, const bool a_destructive){
  CH_TIME("ito_solver::add_particles(container");
  if(m_verbosity > 5){
    pout() << m_name + "::add_particles(container)" << endl;
  }

  this->add_particles(a_part.get_particles(), a_destructive);
}

void ito_solver::add_particles(AMRParticles<ito_particle>& a_part, const bool a_destructive){
  CH_TIME("ito_solver::add_particles(amr)");
  if(m_verbosity > 5){
    pout() << m_name + "::add_particles(amr)" << endl;
  }

  const int finest_level = m_amr->get_finest_level();
  
  for (int lvl = 0; lvl <= finest_level; lvl++){
    this->add_particles(*a_part[lvl], lvl, a_destructive);
  }
}

void ito_solver::add_particles(ParticleData<ito_particle>& a_part, const int a_lvl, const bool a_destructive){
  CH_TIME("ito_solver::add_particles(lvl)");
  if(m_verbosity > 5){
    pout() << m_name + "::add_particles(lvl)" << endl;
  }
  const DisjointBoxLayout& dbl = m_amr->get_grids()[a_lvl];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    this->add_particles(a_part[dit()], a_lvl, dit(), a_destructive);
  }
}

void ito_solver::add_particles(ListBox<ito_particle>& a_part, const int a_lvl, const DataIndex a_dit, const bool a_destructive){
  CH_TIME("ito_solver::add_particles(lvl, dit)");
  if(m_verbosity > 5){
    pout() << m_name + "::add_particles(lvl, dit)" << endl;
  }

  ListBox<ito_particle>& my_particles = m_particles[a_lvl][a_dit];
  
  if(a_destructive){
    my_particles.addItemsDestructive(a_part.listItems());
  }
  else{
    my_particles.addItems(a_part.listItems());
  }
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

void ito_solver::pre_regrid(const int a_base, const int a_old_finest_level){
  CH_TIME("ito_solver::pre_regrid");
  if(m_verbosity > 5){
    pout() << m_name + "::pre_regrid" << endl;
  }

  m_particles.pre_regrid(a_base);
}

AMRParticles<ito_particle>& ito_solver::get_particles(){
  CH_TIME("ito_solver::get_particles");
  if(m_verbosity > 5){
    pout() << m_name + "::get_particles" << endl;
  }

  return m_particles.get_particles();
}

AMRParticles<ito_particle>& ito_solver::get_eb_particles(){
  CH_TIME("ito_solver::get_eb_particles");
  if(m_verbosity > 5){
    pout() << m_name + "::get_eb_particles" << endl;
  }

  return m_eb_particles.get_particles();
}

AMRParticles<ito_particle>& ito_solver::get_domain_particles(){
  CH_TIME("ito_solver::get_domain_particles");
  if(m_verbosity > 5){
    pout() << m_name + "::get_domain_particles" << endl;
  }

  return m_domain_particles.get_particles();
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

  List<ito_particle>& particleList = m_particles[a_lvl][a_dit].listItems();

  EBParticleInterp meshInterp(box, dx, origin);
  meshInterp.interpolateVelocity(particleList, vel_fab, m_deposition);
}

void ito_solver::interpolate_diffusion(){
  CH_TIME("ito_solver::interpolate_diffusion");
  if(m_verbosity > 5){
    pout() << m_name + "::interpolate_diffusion" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      this->interpolate_diffusion(lvl, dit());
    }
  }
}

void ito_solver::interpolate_diffusion(const int a_lvl, const DataIndex& a_dit){
  CH_TIME("ito_solver::interpolate_diffusion");
  if(m_verbosity > 5){
    pout() << m_name + "::interpolate_diffusion" << endl;
  }

  const EBCellFAB& dco_cell   = (*m_diffco_cell[a_lvl])[a_dit];
  const FArrayBox& dco_fab   = dco_cell.getFArrayBox();
  const RealVect dx          = m_amr->get_dx()[a_lvl]*RealVect::Unit;
  const RealVect origin      = m_amr->get_prob_lo();
  const Box box              = m_amr->get_grids()[a_lvl][a_dit];

  List<ito_particle>& particleList = m_particles[a_lvl][a_dit].listItems();

  EBParticleInterp meshInterp(box, dx, origin);
  meshInterp.interpolateDiffusion(particleList, dco_fab, m_deposition);
}

RealVect ito_solver::random_gaussian() {
  const double d = m_gauss01(m_rng);
  return d*this->random_direction();
}

RealVect ito_solver::random_direction(){
  const Real EPS = 1.E-8;
#if CH_SPACEDIM==2
  Real x1 = 2.0;
  Real x2 = 2.0;
  Real r  = x1*x1 + x2*x2;
  while(r >= 1.0 || r < EPS){
    x1 = m_udist11(m_rng);
    x2 = m_udist11(m_rng);
    r  = x1*x1 + x2*x2;
  }

  return RealVect(x1,x2)/sqrt(r);
#elif CH_SPACEDIM==3
  Real x1 = 2.0;
  Real x2 = 2.0;
  Real r  = x1*x1 + x2*x2;
  while(r >= 1.0 || r < EPS){
    x1 = m_udist11(m_rng);
    x2 = m_udist11(m_rng);
    r  = x1*x1 + x2*x2;
  }

  const Real x = 2*x1*sqrt(1-r);
  const Real y = 2*x2*sqrt(1-r);
  const Real z = 1 - 2*r;

  return RealVect(x,y,z);
#endif
}

Real ito_solver::compute_min_drift_dt(const Real a_maxCellsToMove) const{
  CH_TIME("ito_solver::compute_drift_dt(allAMRlevels, maxCellsToMove)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_drift_dt(allAMRlevels, maxCellsToMove)" << endl;
  }

  Vector<Real> dt = this->compute_drift_dt(a_maxCellsToMove);

  Real minDt = dt[0];
  for (int lvl = 0; lvl < dt.size(); lvl++){
    minDt = Min(minDt, dt[lvl]);
  }

  return minDt;
}

Vector<Real> ito_solver::compute_drift_dt(const Real a_maxCellsToMove) const {
  CH_TIME("ito_solver::compute_drift_dt(amr, maxCellsToMove)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_drift_dt(amr, maxCellsToMove)" << endl;
  }

  Vector<Real> dt = this->compute_drift_dt();
  for (int lvl = 0; lvl < dt.size(); lvl++){
    dt[lvl] = dt[lvl]*a_maxCellsToMove;
  }

  return dt;

}

Vector<Real> ito_solver::compute_drift_dt() const {
  CH_TIME("ito_solver::compute_drift_dt(amr)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_drift_dt(amr)" << endl;
  }

  const int finest_level = m_amr->get_finest_level();

  Vector<Real> dt(1 + finest_level, 1.2345E6);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    dt[lvl] = this->compute_drift_dt(lvl);
  }

  return dt;
}

Real ito_solver::compute_drift_dt(const int a_lvl) const {
  CH_TIME("ito_solver::compute_drift_dt(level)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_drift_dt(level)" << endl;
  }

  Real dt = 1.E99;
  
  const DisjointBoxLayout& dbl = m_amr->get_grids()[a_lvl];
  const RealVect dx = m_amr->get_dx()[a_lvl]*RealVect::Unit;
  
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Real boxDt = this->compute_drift_dt(a_lvl, dit(), dx);
    dt = Min(dt, boxDt);
  }

#ifdef CH_MPI
    Real tmp = 1.;
    int result = MPI_Allreduce(&dt, &tmp, 1, MPI_CH_REAL, MPI_MIN, Chombo_MPI::comm);
    if(result != MPI_SUCCESS){
      MayDay::Error("ito_solver::compute_drift_level(lvl) - communication error on norm");
    }
    dt = tmp;
#endif  

  return dt;
}

Real ito_solver::compute_drift_dt(const int a_lvl, const DataIndex& a_dit, const RealVect a_dx) const{
  CH_TIME("ito_solver::compute_drift_dt(level, dataindex, dx)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_drift_dt(level, dataindex, dx)" << endl;
  }

  const List<ito_particle>& particleList = m_particles[a_lvl][a_dit].listItems();
  ListIterator<ito_particle> lit(particleList);

  Real dt = 1.E99;
  
  for (lit.rewind(); lit; ++lit){
    const ito_particle& p = particleList[lit];
    const RealVect& v = p.velocity();

    const int maxDir = v.maxDir(true);
    const Real thisDt = a_dx[maxDir]/Abs(v[maxDir]);

    dt = Min(dt, thisDt);
  }

  return dt;
}

Real ito_solver::compute_min_diffusion_dt(const Real a_maxCellsToHop) const{
  CH_TIME("ito_solver::compute_min_diffusion_dt(min, maxCellsToHop)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_min_diffusion_dt(min, maxCellsToHop)" << endl;
  }

  Vector<Real> dt = this->compute_diffusion_dt(a_maxCellsToHop);
  Real minDt = dt[0];
  for (int lvl = 0; lvl < dt.size(); lvl++){
    minDt = Min(minDt, dt[lvl]);
  }

  return minDt;
}

Vector<Real> ito_solver::compute_diffusion_dt(const Real a_maxCellsToHop) const{
  CH_TIME("ito_solver::compute_min_diffusion_dt(maxCellsToHop)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_min_diffusion_dt(maxCellsToHop)" << endl;
  }

  Vector<Real> dt = this->compute_diffusion_dt();
  for (int lvl = 0; lvl < dt.size(); lvl++){
    dt[lvl] = dt[lvl]*a_maxCellsToHop;
  }

  return dt;
}

Vector<Real> ito_solver::compute_diffusion_dt() const{
  CH_TIME("ito_solver::compute_diffusion_dt");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_diffusion_dt" << endl;
  }

  const int finest_level = m_amr->get_finest_level();
    
  Vector<Real> dt(1 + finest_level, 1.2345E6);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    dt[lvl] = this->compute_diffusion_dt(lvl);
  }

  return dt;
}

Real ito_solver::compute_diffusion_dt(const int a_lvl) const{
  CH_TIME("ito_solver::compute_diffusion_dt(level)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_diffusion_dt(level)" << endl;
  }

  
  Real dt = 1.E99;
  
  const DisjointBoxLayout& dbl = m_amr->get_grids()[a_lvl];
  const RealVect dx = m_amr->get_dx()[a_lvl]*RealVect::Unit;
  
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Real boxDt = this->compute_diffusion_dt(a_lvl, dit(), dx);
    dt = Min(dt, boxDt);
  }

#ifdef CH_MPI
  Real tmp = 1.;
  int result = MPI_Allreduce(&dt, &tmp, 1, MPI_CH_REAL, MPI_MIN, Chombo_MPI::comm);
  if(result != MPI_SUCCESS){
    MayDay::Error("ito_solver::compute_diffusion_level(lvl) - communication error on norm");
  }
  dt = tmp;
#endif

  return dt;
}

Real ito_solver::compute_diffusion_dt(const int a_lvl, const DataIndex& a_dit, const RealVect a_dx) const{
  CH_TIME("ito_solver::compute_diffusion_dt(level, dataindex, dx)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_diffusion_dt(level, dataindex, dx)" << endl;
  }
  const List<ito_particle>& particleList = m_particles[a_lvl][a_dit].listItems();
  ListIterator<ito_particle> lit(particleList);

  Real dt = 1.E99;
  
  for (lit.rewind(); lit; ++lit){
    const ito_particle& p = particleList[lit];
    
    const Real thisDt = a_dx[0]/(sqrt(2.0*p.diffusion()));
    
    dt = Min(dt, thisDt);
  }

  return dt;
}

void ito_solver::remap(){
  CH_TIME("ito_solver::remap");
  if(m_verbosity > 5){
    pout() << m_name + "::remap" << endl;
  }

  m_particles.remap();
}

phase::which_phase ito_solver::get_phase() const{
  return m_phase;
}
