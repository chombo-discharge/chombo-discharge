/*!
  @file mfindexspace.cpp
  @brief Implementation of mfis.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "mfis.H"
#include <AllRegularService.H>

mfis::mfis(){
  m_ebis.resize(Phase::num_phases);
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


  m_ebis[Phase::Gas]->define(a_domain,   a_origin, a_dx, *a_geoservers[Phase::Gas],   a_nCellMax, a_max_coar);

  if(a_geoservers[Phase::Solid] == NULL){
    m_ebis[Phase::Solid] = RefCountedPtr<EBIndexSpace> (NULL);
  }
  else{
    m_ebis[Phase::Solid]->define(a_domain, a_origin, a_dx, *a_geoservers[Phase::Solid], a_nCellMax, a_max_coar);
  }
}
  
const RefCountedPtr<EBIndexSpace>& mfis::get_ebis(Phase::WhichPhase a_whichEBIS) const {
  return m_ebis[a_whichEBIS];
}
