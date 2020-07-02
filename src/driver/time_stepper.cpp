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

void time_stepper::load_balance(Vector<int>& a_procs,
				Vector<Box>& a_boxes,
				const DisjointBoxLayout& a_dbl,
				const int a_level){
  CH_TIME("time_stepper::load_balance");
  if(m_verbosity > 5){
    pout() << "time_stepper::load_balance" << endl;
  }

  Vector<long int> a_loads;
  a_boxes = a_dbl.boxArray();
  a_loads.resize(a_dbl.size(), 0);

  for (DataIterator dit = a_dbl.dataIterator(); dit.ok(); ++dit){
    const int idx = dit().intCode();
    const Box bx  = a_dbl.get(dit());
    a_loads[idx]  = bx.numPts();
  }

  // Gather loads
#ifdef CH_MPI
  int count = a_loads.size();
  Vector<long int> tmp(count);
  MPI_Allreduce(&(a_loads[0]),&(tmp[0]), count, MPI_LONG, MPI_SUM, Chombo_MPI::comm);
  a_loads = tmp;
#endif

  // Load balance
  LoadBalance(a_procs, a_loads, a_boxes);
}
