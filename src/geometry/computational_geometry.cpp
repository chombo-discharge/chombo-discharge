/*!
  @file computational_geometry.cpp
  @brief Implementation of computational_geometry.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include <MFIndexSpace.H>
#include <IntersectionIF.H>
#include <UnionIF.H>
#include <AllRegularService.H>
#include <GeometryService.H>
#include <WrappedGShop.H>

#include "computational_geometry.H"

int computational_geometry::s_minRef = 0;
int computational_geometry::s_maxRef = 0;

computational_geometry::computational_geometry(){

  //
  m_eps0 = 1.0;
  m_electrodes.resize(0);
  m_dielectrics.resize(0);
  m_wallbc.resize(2*SpaceDim);

  for (int i = 0; i < 2*SpaceDim; i++){
    m_wallbc[i] = RefCountedPtr<wall_bc> (NULL);
  }
  
  m_ebis_mf = RefCountedPtr<MFIndexSpace> (new MFIndexSpace());
}

computational_geometry::~computational_geometry(){
}

const Vector<dielectric>& computational_geometry::get_dielectrics() const  {
  return m_dielectrics;
}

const Vector<electrode>& computational_geometry::get_electrodes() const  {
  return m_electrodes;
}

const Real& computational_geometry::get_eps0() const {
  return m_eps0;
}

const RefCountedPtr<MFIndexSpace> computational_geometry::get_ebis_mf() const {
  return m_ebis_mf;
}

void computational_geometry::set_dielectrics(Vector<dielectric>& a_dielectrics){
  m_dielectrics = a_dielectrics;
}

void computational_geometry::set_electrodes(Vector<electrode>& a_electrodes){
  m_electrodes = a_electrodes;
}

void computational_geometry::set_dirichlet_wall_bc(const int a_dir, Side::LoHiSide a_side, const bool a_live){

  const int idx = this->map_bc(a_dir, a_side);
  m_wallbc[idx] = RefCountedPtr<wall_bc> (new wall_bc(a_dir, a_side, WallBC::Dirichlet));
  m_wallbc[idx]->set_live(a_live);
}

void computational_geometry::set_neumann_wall_bc(const int a_dir, Side::LoHiSide a_side, const Real a_value){

  const int idx = this->map_bc(a_dir, a_side);
  m_wallbc[idx] = RefCountedPtr<wall_bc> (new wall_bc(a_dir, a_side, WallBC::Neumann));
  m_wallbc[idx]->set_value(a_value);
}

const wall_bc& computational_geometry::get_wall_bc(const int a_dir, Side::LoHiSide a_side) const{
  return *m_wallbc[this->map_bc(a_dir, a_side)];
}

void computational_geometry::set_eps0(const Real a_eps0){
  m_eps0 = a_eps0;
}

void computational_geometry::build_geometries(const physical_domain& a_physDom,
					      const ProblemDomain&   a_finestDomain,
					      const Real&            a_finestDx,
					      const int&             a_nCellMax,
					      const int&             a_maxCoarsen){

  this->sanity_check();
  
  // Build geoservers
  Vector<GeometryService*> geoservers(2, NULL);
  this->build_gas_geoserv(geoservers[Phase::Gas],   a_finestDomain, a_physDom.get_prob_lo(), a_finestDx);
  this->build_solid_geoserv(geoservers[Phase::Solid], a_finestDomain, a_physDom.get_prob_lo(), a_finestDx);

  // Define MF
  m_ebis_mf->define(a_finestDomain.domainBox(),
		    a_physDom.get_prob_lo(),
		    a_finestDx,
		    geoservers,
		    a_nCellMax,
		    a_maxCoarsen);
}

void computational_geometry::build_gas_geoserv(GeometryService*&    a_geoserver,
					       const ProblemDomain& a_finestDomain,
					       const RealVect&      a_origin,
					       const Real&          a_dx){
  
  // The gas ebis is the intersection of electrodes and dielectrics
  Vector<BaseIF*> parts;
  for (int i = 0; i < m_dielectrics.size(); i++){
    parts.push_back(&(*m_dielectrics[i].get_function()));
  }
  for (int i = 0; i < m_electrodes.size(); i++){
    parts.push_back(&(*m_electrodes[i].get_function()));
  }

  // Create geoshop; either all regular or an intersection
  if(parts.size() == 0){
    a_geoserver = new AllRegularService();
  }
  else {
    RefCountedPtr<BaseIF> baseif = RefCountedPtr<BaseIF> (new IntersectionIF(parts));

    a_geoserver = static_cast<GeometryService*> (new WrappedGShop(baseif, a_origin, a_dx, a_finestDomain, s_minRef, s_maxRef));
  }
}

void computational_geometry::build_solid_geoserv(GeometryService*&    a_geoserver,
						 const ProblemDomain& a_finestDomain,
						 const RealVect&      a_origin,
						 const Real&          a_dx){

  // The gas ebis is the intersection of dielectrics
  Vector<BaseIF*> parts_dielectrics, parts_electrodes;
  for (int i = 0; i < m_dielectrics.size(); i++){
    parts_dielectrics.push_back(&(*m_dielectrics[i].get_function()));
  }
  for (int i = 0; i < m_electrodes.size(); i++){
    parts_electrodes.push_back(&(*m_electrodes[i].get_function()));
  }

  // Create EBIndexSpace. If there are no solid phases, return ull
  if(parts_dielectrics.size() == 0){
    a_geoserver = new AllRegularService();
  }
  else {
    Vector<BaseIF*> parts;

    //
    RefCountedPtr<BaseIF> diel_baseif = RefCountedPtr<BaseIF> (new IntersectionIF(parts_dielectrics));
    RefCountedPtr<BaseIF> elec_baseif = RefCountedPtr<BaseIF> (new IntersectionIF(parts_electrodes));
    parts.push_back(&(*diel_baseif));
    parts.push_back(&(*elec_baseif));

    // Create union to get rid of electrode part
    RefCountedPtr<BaseIF> baseif = RefCountedPtr<BaseIF> (new UnionIF(parts));

    // Set up the rest
    a_geoserver = static_cast<GeometryService*> (new WrappedGShop(baseif, a_origin, a_dx, a_finestDomain, s_minRef, s_maxRef));

  }
}

void computational_geometry::sanity_check(){
  
  for (int dir = 0; dir < SpaceDim; dir++){
    for (SideIterator sideit; sideit.ok(); ++sideit){
      if(!m_wallbc[map_bc(dir, sideit())].isNull()){
	MayDay::Abort("computational_geometry::sanity_check() failed. Wall BC has not been set properly");
      }
    }
  }
      
}

int computational_geometry::map_bc(const int a_dir, const Side::LoHiSide a_side) const {
  const int iside = (a_side == Side::Lo) ? 0 : 1;

  return 2*a_dir + iside;
}
