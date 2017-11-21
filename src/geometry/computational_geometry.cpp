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
#include <GeometryShop.H>
#include <ComplementIF.H>

#include "computational_geometry.H"

int computational_geometry::s_minRef = 0;
int computational_geometry::s_maxRef = 0;

computational_geometry::computational_geometry(){

  //
  m_eps0 = 1.0;
  m_electrodes.resize(0);
  m_dielectrics.resize(0);
  
  m_mfis = RefCountedPtr<mfis> (new mfis());
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

const RefCountedPtr<mfis> computational_geometry::get_mfis() const {
  return m_mfis;
}

void computational_geometry::set_dielectrics(Vector<dielectric>& a_dielectrics){
  m_dielectrics = a_dielectrics;
}

void computational_geometry::set_electrodes(Vector<electrode>& a_electrodes){
  m_electrodes = a_electrodes;
}

void computational_geometry::set_eps0(const Real a_eps0){
  m_eps0 = a_eps0;
}

void computational_geometry::build_geometries(const physical_domain& a_physDom,
					      const ProblemDomain&   a_finestDomain,
					      const Real&            a_finestDx,
					      const int&             a_nCellMax,
					      const int&             a_maxCoarsen){


  
  // Build geoservers
  Vector<GeometryService*> geoservers(2, NULL);
  this->build_gas_geoserv(geoservers[Phase::Gas],     a_finestDomain, a_physDom.get_prob_lo(), a_finestDx);
  this->build_solid_geoserv(geoservers[Phase::Solid], a_finestDomain, a_physDom.get_prob_lo(), a_finestDx);

  m_mfis->define(a_finestDomain.domainBox(), // Define MF
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
    parts.push_back(&(*(m_dielectrics[i].get_function())));
  }
  for (int i = 0; i < m_electrodes.size(); i++){
    parts.push_back(&(*(m_electrodes[i].get_function())));
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

  // Create EBIndexSpace. If there are no solid phases, return null
  if(parts_dielectrics.size() == 0){
    a_geoserver = NULL;
  }
  else {
    Vector<BaseIF*> union_parts;

    RefCountedPtr<BaseIF> diel_baseif = RefCountedPtr<BaseIF> (new IntersectionIF(parts_dielectrics)); // Gas outside dielectric
    RefCountedPtr<BaseIF> elec_baseif = RefCountedPtr<BaseIF> (new IntersectionIF(parts_electrodes));  // Gas outside electrodes
    RefCountedPtr<BaseIF> diel_compif = RefCountedPtr<BaseIF> (new ComplementIF(*diel_baseif));        // Inside of dielectrics
    
    union_parts.push_back(&(*elec_baseif)); // Parts for union
    union_parts.push_back(&(*diel_compif)); // Parts for union
    
    RefCountedPtr<BaseIF> baseif = RefCountedPtr<BaseIF> (new UnionIF(union_parts)); // Union to get rid of electrode

    a_geoserver = static_cast<GeometryService*> (new WrappedGShop(baseif, a_origin, a_dx, a_finestDomain, s_minRef, s_maxRef));
  }
}
