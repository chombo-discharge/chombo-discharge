/*!
  @file mfindexspace.cpp
  @brief Implementation of mfis.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "mfis.H"
#include <AllRegularService.H>

mfis::mfis(){
  m_ebis.resize(phase::num_phases);
  for (int i = 0; i < m_ebis.size(); i++){
    m_ebis[i] = RefCountedPtr<EBIndexSpace> (new EBIndexSpace());
  }
}

mfis::~mfis(){
}

void mfis::define(const Box                     & a_domain,
		  const RealVect                & a_origin,
		  const Real                    & a_dx,
		  const Vector<GeometryService*>& a_geoservers,
		  int                             a_nCellMax,
		  int                             a_max_coar,
		  bool                            a_fix_phase){


  m_ebis[phase::gas]->define(a_domain,   a_origin, a_dx, *a_geoservers[phase::gas],   a_nCellMax, a_max_coar);

  if(a_geoservers[phase::solid] == NULL){
    m_ebis[phase::solid] = RefCountedPtr<EBIndexSpace> (NULL);
  }
  else{
    m_ebis[phase::solid]->define(a_domain, a_origin, a_dx, *a_geoservers[phase::solid], a_nCellMax, a_max_coar);
  }
}
  
const RefCountedPtr<EBIndexSpace>& mfis::get_ebis(phase::which_phase a_whichEBIS) const {
  return m_ebis[a_whichEBIS];
}

const RefCountedPtr<EBIndexSpace>& mfis::get_ebis(const int a_phase) const {
  return m_ebis[a_phase];
}


const int mfis::num_phases() const{
  return m_ebis.size();
}

const IntVectSet mfis::interface_region(const ProblemDomain& a_domain) const {
  CH_TIME("mfis::interface_region");

  const int which_level = m_ebis[0]->getLevel(a_domain);
  
  IntVectSet iface_reg = m_ebis[0]->irregCells(which_level);
  for (int i = 1; i < m_ebis.size(); i++){
    iface_reg &= m_ebis[i]->irregCells(which_level);
  }

  return iface_reg;
}
