/*!
  @file   mfis.cpp
  @brief  Implementation of mfis.H
  @author Robert Marskar
  @date   Nov. 2017
*/

#include "mfis.H"
#include "computational_geometry.H"
#include "memrep.H"

#include <AllRegularService.H>

#include "CD_NamespaceHeader.H"

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

  // Define the gas geoserver
  if(computational_geometry::s_use_new_gshop){
    m_ebis[phase::gas]->setDistributedData();
  }
  m_ebis[phase::gas]->define(a_domain,   a_origin, a_dx, *a_geoservers[phase::gas],   a_nCellMax, a_max_coar);
#if 1 // Debug
  memrep::get_max_min_memory();
#endif


  // Define the solid state geoserver. This EBIS might not exist. 
  if(a_geoservers[phase::solid] == NULL){
    m_ebis[phase::solid] = RefCountedPtr<EBIndexSpace> (NULL);
  }
  else{
    if(computational_geometry::s_use_new_gshop){
      m_ebis[phase::solid]->setDistributedData();
    }
    m_ebis[phase::solid]->define(a_domain, a_origin, a_dx, *a_geoservers[phase::solid], a_nCellMax, a_max_coar);
#if 1 // Debug
    memrep::get_max_min_memory();
#endif
  }
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
#if 0  
  IntVectSet iface_reg = m_ebis[0]->irregCells(which_level);
  for (int i = 1; i < this->num_phases(); i++){
    iface_reg &= m_ebis[i]->irregCells(which_level);
  }

  return iface_reg;
#else
  IntVectSet irreg0;
  IntVectSet irreg1;
  IntVectSet myIrregCellsPhase0 = m_ebis[0]->irregCells(which_level);
  IntVectSet myIrregCellsPhase1 = m_ebis[1]->irregCells(which_level);

  Vector<IntVectSet> allIrregCellsPhase0(numProc());
  Vector<IntVectSet> allIrregCellsPhase1(numProc());
  
  gather(allIrregCellsPhase0, myIrregCellsPhase0, 0);
  gather(allIrregCellsPhase1, myIrregCellsPhase1, 0);

  IntVectSet ret;
  
  if(procID() == 0){
    for (int i = 0; i < allIrregCellsPhase1.size(); i++){
      irreg0 |= allIrregCellsPhase1[i];
    }

    for (int i = 0; i < allIrregCellsPhase0.size(); i++){
      irreg1 |= allIrregCellsPhase0[i];
    }

    ret = irreg1 & irreg0;

  }

  broadcast(ret, 0);

  return ret;
#endif
}
#include "CD_NamespaceFooter.H"
