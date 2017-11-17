/*!
  @file ComputationalGeometry.cpp
  @brief Implementation of ComputationalGeometry.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include <MFIndexSpace.H>
#include <IntersectionIF.H>
#include <UnionIF.H>
#include <AllRegularService.H>
#include <GeometryService.H>
#include <WrappedGShop.H>

#include "ComputationalGeometry.H"

//
ComputationalGeometry::ComputationalGeometry(){

  //
  m_eps0 = 1.0;
  m_electrodes.resize(0);
  m_dielectrics.resize(0);

  m_ebis_mf = RefCountedPtr<MFIndexSpace> (new MFIndexSpace());
}

//
ComputationalGeometry::~ComputationalGeometry(){
}

//
const Vector<Dielectric>& ComputationalGeometry::get_dielectrics() const  {
  return m_dielectrics;
}

//
const Vector<Electrode>& ComputationalGeometry::get_electrodes() const  {
  return m_electrodes;
}

//
const Real& ComputationalGeometry::get_eps0() const {
  return m_eps0;
}

//
const RefCountedPtr<MFIndexSpace> ComputationalGeometry::get_ebis_mf() const {
  return m_ebis_mf;
}

//
void ComputationalGeometry::set_dielectrics(Vector<Dielectric>& a_dielectrics){
  m_dielectrics = a_dielectrics;
}

//
void ComputationalGeometry::set_electrodes(Vector<Electrode>& a_electrodes){
  m_electrodes = a_electrodes;
}

//
void ComputationalGeometry::set_eps0(const Real a_eps0){
  m_eps0 = a_eps0;
}

//
void ComputationalGeometry::build_geometries(const PhysicalDomain& a_physDom,
					     const ProblemDomain&  a_finestDomain,
					     const Real&           a_finestDx,
					     const int&            a_nCellMax,
					     const int&            a_maxCoarsen){

  // Build geoservers
  Vector<GeometryService*> geoservers(2, NULL);
  this->build_gas_geoserv(geoservers[0],   a_finestDomain, a_physDom.get_prob_lo(), a_finestDx);
  this->build_solid_geoserv(geoservers[1], a_finestDomain, a_physDom.get_prob_lo(), a_finestDx);

  // Define MF
  m_ebis_mf->define(a_finestDomain.domainBox(),
		    a_physDom.get_prob_lo(),
		    a_finestDx,
		    geoservers,
		    a_nCellMax,
		    a_maxCoarsen);
}

//
void ComputationalGeometry::build_gas_geoserv(GeometryService*     a_geoserver,
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
    
    a_geoserver = static_cast<GeometryService*> (new WrappedGShop(baseif, a_origin, a_dx, a_finestDomain, 0, 0));
  }
}

//
void ComputationalGeometry::build_solid_geoserv(GeometryService*     a_geoserver,
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
    a_geoserver = static_cast<GeometryService*> (new WrappedGShop(baseif, a_origin, a_dx, a_finestDomain, 0, 0));

  }
}
