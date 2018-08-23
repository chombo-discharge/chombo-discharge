/*!
  @file load_balance.cpp
  @details Implementatino of load_balance.H
  @author Robert Marskar
  @date Jan. 2018
*/

#include "load_balance.H"
#include "EBEllipticLoadBalance.H"

void load_balance::balance_volume(Vector<int>& a_procs, const Vector<Box>& a_boxes){
  LoadBalance(a_procs, a_boxes);
}

void load_balance::balance_elliptic(Vector<int>&                       a_procs,
				    const Vector<Box>&                 a_boxes,
				    const RefCountedPtr<EBIndexSpace>& a_ebis,
				    const ProblemDomain&               a_domain,
				    const bool                         a_verbose){

  EBEllipticLoadBalance(a_procs, a_boxes, a_domain, a_verbose, a_ebis);
}

void load_balance::balance_multifluid(Vector<int>&               a_procs,
				      const Vector<Box>&         a_boxes,
				      const RefCountedPtr<mfis>& a_mfis,
				      const ProblemDomain&       a_domain,
				      const bool                 a_verbose){

  MayDay::Abort("load_balance::balance_multifluid - not implemented (yet)");
  
  // Load balance the conventional way
  Vector<Box> boxes = a_boxes;
  LoadBalance(a_procs, a_boxes);

  
  MayDay::Abort("load_balance::balance_multifluid - not implemented");
}

void load_balance::get_mfpoisson_loads(Vector<unsigned long long>& a_loads,
				       Vector<Box>&                a_boxes,
				       RefCountedPtr<mfis>&        a_mfis,
				       const DisjointBoxLayout&    a_dbl,
				       const ProblemDomain&        a_domain){
  
}
