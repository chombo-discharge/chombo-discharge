/*!
  @file   time_stepper.cpp
  @brief  Implementation of time_stepper.H
  @author Robert Marskar
  @date   March 2020
*/

#include "time_stepper.H"
#include "LoadBalance.H"

void time_stepper::set_amr(const RefCountedPtr<amr_mesh>& a_amr){
  CH_TIME("time_stepper::set_amr");
  if(m_verbosity > 5){
    pout() << "time_stepper::set_amr" << endl;
  }

  m_amr = a_amr;
}

void time_stepper::set_computational_geometry(const RefCountedPtr<computational_geometry>& a_compgeom){
  CH_TIME("time_stepper::set_computational_geometry");
  if(m_verbosity > 5){
    pout() << "time_stepper::set_computational_geometry" << endl;
  }

  m_compgeom = a_compgeom;
}

bool time_stepper::need_to_regrid(){
  CH_TIME("time_stepper::need_to_regrid");
  if(m_verbosity > 5){
    pout() << "time_stepper::need_to_regrid" << endl;
  }

  return false;
}

int time_stepper::get_redistribution_regsize() const {
  CH_TIME("time_stepper::get_redistribution_regsize");
  if(m_verbosity > 5){
    pout() << "time_stepper::get_redistribution_regsize" << endl;
  }

  return 1;
}

Vector<long int> time_stepper::get_checkpoint_loads(const std::string a_realm, const int a_level) const {
  CH_TIME("time_stepper::get_checkpoint_loads");
  if(m_verbosity > 5){
    pout() << "time_stepper::get_checkpoint_loads" << endl;
  }

  const DisjointBoxLayout& dbl = m_amr->get_grids(a_realm)[a_level];
  const Vector<Box>& a_boxes = dbl.boxArray();

  Vector<long int> loads(a_boxes.size(), 0L);
  for (int i = 0; i < a_boxes.size(); i++){
    loads[i] = a_boxes[i].numPts();
  }

  return loads;
}

bool time_stepper::load_balance_realm(const std::string a_realm) const {
  CH_TIME("time_stepper::load_balance_realm");
  if(m_verbosity > 5){
    pout() << "time_stepper::load_balance_realm" << endl;
  }

  return false;
}

void time_stepper::load_balance_boxes(Vector<Vector<int> >&             a_procs,
				      Vector<Vector<Box> >&             a_boxes,
				      const std::string                 a_realm,              
				      const Vector<DisjointBoxLayout>&  a_grids,
				      const int                         a_lmin,
				      const int                         a_finest_level){
  CH_TIME("time_stepper::load_balance_boxes");
  if(m_verbosity > 5){
    pout() << "time_stepper::load_balance_boxes" << endl;
  }

  a_procs.resize(1 + a_finest_level);
  a_boxes.resize(1 + a_finest_level);

  for (int lvl = 0; lvl <= a_finest_level; lvl++){
    a_boxes[lvl] = a_grids[lvl].boxArray();

    LoadBalance(a_procs[lvl], a_boxes[lvl]);
  }
}
