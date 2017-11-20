/*!
  @file plasma_engine.cpp
  @brief Implementation of plasma_engine.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "plasma_engine.H"

plasma_engine::plasma_engine(){
  CH_TIME("plasma_engine::plasma_engine(weak)");
  if(m_verbosity > 2){
    pout() << "plasma_engine::plasma_engine(weak)" << endl;
  }
}

plasma_engine::plasma_engine(const RefCountedPtr<physical_domain>&        a_physdom,
			     const RefCountedPtr<computational_geometry>& a_compgeom,
			     const RefCountedPtr<plasma_kinetics>&        a_plaskin,
			     const RefCountedPtr<time_stepper>&           a_timestepper,
			     const RefCountedPtr<amr_mesh>&               a_amr,
			     const RefCountedPtr<cell_tagger>&            a_celltagger){
  CH_TIME("plasma_engine::plasma_engine(full)");
  if(m_verbosity > 2){
    pout() << "plasma_engine::plasma_engine(full)" << endl;
  }

  this->set_physical_domain(a_physdom);         // Set physical domain
  this->set_computational_geometry(a_compgeom); // Set computational geometry
  this->set_plasma_kinetics(a_plaskin);         // Set plasma kinetics
  this->set_time_stepper(a_timestepper);        // Set time stepper
  this->set_amr(a_amr);                         // Set amr
  this->set_cell_tagger(a_celltagger);          // Set cell tagger

  m_amr->sanity_check();  // Sanity check, make sure everything is set up correctly
  m_amr->build_domains(); // Build domains and resolutions, nothing else

  this->define_cell_tagger();
  this->allocate_wall_bc();
}

plasma_engine::~plasma_engine(){
  CH_TIME("plasma_engine::~plasma_engine");
}

void plasma_engine::allocate_wall_bc(){
  m_wallbc.resize(2*SpaceDim);
  for (int i = 0; i < 2*SpaceDim; i++){
    m_wallbc[i] = RefCountedPtr<wall_bc> (NULL);
  }
}

void plasma_engine::set_verbosity(const int a_verbosity){
  CH_TIME("plasma_engine::set_verbosity");
  if(m_verbosity > 5){
    pout() << "plasma_engine::set_verbosity" << endl;
  }
  m_verbosity = a_verbosity;
}

void plasma_engine::set_computational_geometry(const RefCountedPtr<computational_geometry>& a_compgeom){
  CH_TIME("plasma_engine::set_computational_geometry");
  if(m_verbosity > 3){
    pout() << "plasma_engine::set_computational_geometry" << endl;
  }
  m_compgeom = a_compgeom;
  m_mfis     = a_compgeom->get_ebis_mf();
}

void plasma_engine::set_plasma_kinetics(const RefCountedPtr<plasma_kinetics>& a_plaskin){
  m_plaskin = a_plaskin;
}

void plasma_engine::set_time_stepper(const RefCountedPtr<time_stepper>& a_timestepper){
  m_timestepper = a_timestepper;
}

void plasma_engine::set_cell_tagger(const RefCountedPtr<cell_tagger>& a_celltagger){
  m_celltagger = a_celltagger;
}

void plasma_engine::define_cell_tagger(){
  if(!m_celltagger.isNull()){
    m_celltagger->define(m_plaskin, m_timestepper, m_amr, m_compgeom, m_physdom);
  }
}

void plasma_engine::set_amr(const RefCountedPtr<amr_mesh>& a_amr){
  CH_TIME("plasma_engine::set_amr");
  if(m_verbosity > 2){
    pout() << "plasma_engine::set_amr" << endl;
  }

  m_amr = a_amr;
  m_amr->set_mfis(m_compgeom->get_ebis_mf());
}

void plasma_engine::setup_fresh(){
  CH_TIME("plasma_engine::setup_fresh");
  if(m_verbosity > 2){
    pout() << "plasma_engine::setup_fresh" << endl;
  }

  this->sanity_check();  // Sanity check before doing anything expensive
  

  m_compgeom->build_geometries(*m_physdom, m_amr->get_finest_domain(), m_amr->get_finest_dx(), m_amr->get_max_box_size());
  
  this->get_geom_tags();      // Get geometric tags
  m_amr->regrid(m_geom_tags); // Regrid using geometric tags
}

void plasma_engine::setup_for_restart(const std::string a_restart_file){
  
}

void plasma_engine::initial_regrids(const int a_init_regrids){

}

void plasma_engine::set_dirichlet_wall_bc(const int a_dir, Side::LoHiSide a_side, const Potential::GroundLive a_live){
  CH_TIME("plasma_engine::set_dirichlet_wall_bc");
  if(m_verbosity > 3){
    pout() << "plasma_engine::set_dirichlet_wall_bc" << endl;
  }

  const int idx = this->map_bc(a_dir, a_side);
  m_wallbc[idx] = RefCountedPtr<wall_bc> (new wall_bc(a_dir, a_side, WallBC::Dirichlet));
  m_wallbc[idx]->set_live(a_live);
}

void plasma_engine::set_neumann_wall_bc(const int a_dir, Side::LoHiSide a_side, const Real a_value){
  CH_TIME("plasma_engine::set_neumann_wall_bc");
  if(m_verbosity > 3){
    pout() << "plasma_engine::set_neumann_wall_bc" << endl;
  }

  const int idx = this->map_bc(a_dir, a_side);
  m_wallbc[idx] = RefCountedPtr<wall_bc> (new wall_bc(a_dir, a_side, WallBC::Neumann));
  m_wallbc[idx]->set_value(a_value);
}

void plasma_engine::set_physical_domain(const RefCountedPtr<physical_domain>& a_physdom){
  CH_TIME("plasma_engine::set_physical_domain");
  if(m_verbosity > 3){
    pout() << "plasma_engine::set_physical_domain" << endl;
  }
  m_physdom = a_physdom;
}

void plasma_engine::sanity_check(){
  CH_TIME("plasma_engine::sanity_check");
  if(m_verbosity > 2){
    pout() << "plasma_engine::sanity_check" << endl;
  }
  
  for (int dir = 0; dir < SpaceDim; dir++){
    for (SideIterator sideit; sideit.ok(); ++sideit){
      if(m_wallbc[map_bc(dir, sideit())].isNull()){
	MayDay::Abort("computational_geometry::sanity_check() failed. Wall BC has not been set properly");
      }
    }
  }
}

void plasma_engine::regrid(){
  CH_TIME("plasma_engine::regrid");
  if(m_verbosity > 2){
    pout() << "plasma_engine::regrid" << endl;
  }
}

void plasma_engine::write_plot_file(){
  
}

void plasma_engine::write_checkpoint_file(){
  
}

void plasma_engine::read_checkpoint_file(){
  
}

void plasma_engine::get_geom_tags(){
  CH_TIME("plasma_engine::get_geom_tags");
  if(m_verbosity > 2){
    pout() << "plasma_engine::get_geom_tags" << endl;
  }
  

  const int maxdepth = m_amr->get_max_amr_depth();

  m_geom_tags.resize(1 + maxdepth);

  const EBIndexSpace* const ebis_gas = m_mfis->EBIS(Phase::Gas);
  const EBIndexSpace* const ebis_sol = m_mfis->EBIS(Phase::Solid);

  for (int lvl = 0; lvl <= maxdepth; lvl++){
    const ProblemDomain& cur_dom = m_amr->get_domains()[lvl];
    const int which_level = ebis_gas->getLevel(cur_dom);

    m_geom_tags[lvl].makeEmpty();
    m_geom_tags[lvl] |= ebis_gas->irregCells(which_level);
    m_geom_tags[lvl] |= ebis_sol->irregCells(which_level);
  }

  // Grow tags by 2, this is an ad-hoc fix that prevents ugly grid near EBs
  for (int lvl = 0; lvl <= maxdepth; lvl++){
    m_geom_tags[lvl].grow(2);
  }
}

wall_bc& plasma_engine::get_wall_bc(const int a_dir, Side::LoHiSide a_side) const{
  CH_TIME("plasma_engine::get_wall_bc");
  if(m_verbosity > 2){
    pout() << "plasma_engine::get_wall_bc" << endl;
  }
  return *m_wallbc[this->map_bc(a_dir, a_side)];
}

int plasma_engine::map_bc(const int a_dir, const Side::LoHiSide a_side) const {
  const int iside = (a_side == Side::Lo) ? 0 : 1;

  return 2*a_dir + iside;
}
