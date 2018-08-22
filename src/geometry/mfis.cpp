/*!
  @file   mfis.cpp
  @brief  Implementation of mfis.H
  @author Robert Marskar
  @date   Nov. 2017
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


#if 1
  m_ebis[phase::gas]->define(a_domain,   a_origin, a_dx, *a_geoservers[phase::gas],   a_nCellMax, a_max_coar);

  if(a_geoservers[phase::solid] == NULL){
    m_ebis[phase::solid] = RefCountedPtr<EBIndexSpace> (NULL);
  }
  else{
    m_ebis[phase::solid]->define(a_domain, a_origin, a_dx, *a_geoservers[phase::solid], a_nCellMax, a_max_coar);
  }
#else
  pout() << "defining MFIndexSpace" << endl;
  m_mfis = RefCountedPtr<MFIndexSpace> (new MFIndexSpace());
  m_mfis->define(a_domain, a_origin, a_dx, a_geoservers, a_max_coar, false);

  m_ebis[phase::gas] = RefCountedPtr<EBIndexSpace> (const_cast<EBIndexSpace*> (m_mfis->EBIS(phase::gas)));
  m_ebis[phase::solid] = RefCountedPtr<EBIndexSpace> (const_cast<EBIndexSpace*> (m_mfis->EBIS(phase::solid)));

  pout() << "done defining MFIndexSpace" << endl;
#endif
}
  
const RefCountedPtr<EBIndexSpace>& mfis::get_ebis(const phase::which_phase a_phase) const {
  return m_ebis[a_phase];
}

const RefCountedPtr<EBIndexSpace>& mfis::get_ebis(const int a_phase) const {
  return m_ebis[a_phase];
}

RefCountedPtr<EBIndexSpace>& mfis::get_ebis(const phase::which_phase a_phase){
  return m_ebis[a_phase];
}

RefCountedPtr<EBIndexSpace>& mfis::get_ebis(const int a_phase){
  return m_ebis[a_phase];
}


int mfis::num_phases() const{
  int phases = 0;
  for (int i = 0; i < m_ebis.size(); i++){
    if(!m_ebis[i].isNull()){
      phases++;
    }
  }

  return phases;
}

IntVectSet mfis::interface_region(const ProblemDomain& a_domain) const {
  CH_TIME("mfis::interface_region");

  const int which_level = m_ebis[0]->getLevel(a_domain);
  
  IntVectSet iface_reg = m_ebis[0]->irregCells(which_level);
  for (int i = 1; i < this->num_phases(); i++){
    iface_reg &= m_ebis[i]->irregCells(which_level);
  }

  return iface_reg;
}
