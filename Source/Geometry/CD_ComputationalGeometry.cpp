/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_ComputationalGeometry.cpp
  @brief  Implementation of CD_ComputationalGeometry.H
  @author Robert Marskar
*/

// Chombo includes
#include <MFIndexSpace.H>
#include <IntersectionIF.H>
#include <UnionIF.H>
#include <AllRegularService.H>
#include <GeometryService.H>
#include <WrappedGShop.H>
#include <GeometryShop.H>
#include <ComplementIF.H>

// Our includes
#include <CD_ComputationalGeometry.H>
#include <CD_ScanShop.H>
#include <memrep.H>
#include <CD_NamespaceHeader.H>

#if CH_USE_DOUBLE
Real ComputationalGeometry::s_thresh = 1.E-15;
#elif CH_USE_FLOAT
Real ComputationalGeometry::s_thresh = 1.E-6;
#endif

bool          ComputationalGeometry::s_use_new_gshop = false;
ProblemDomain ComputationalGeometry::s_ScanDomain = ProblemDomain(); // Needs to be set

ComputationalGeometry::ComputationalGeometry(){
  m_eps0 = 1.0;
  m_electrodes.resize(0);
  m_dielectrics.resize(0);
  
  m_multifluidIndexSpace = RefCountedPtr<mfis> (new mfis());
}

ComputationalGeometry::~ComputationalGeometry(){}

const Vector<Dielectric>& ComputationalGeometry::getDielectrics() const  {
  return m_dielectrics;
}

const Vector<Electrode>& ComputationalGeometry::getElectrodes() const  {
  return m_electrodes;
}

const RefCountedPtr<BaseIF>& ComputationalGeometry::getGasImplicitFunction() const{
  return m_gas_if;
}

const RefCountedPtr<BaseIF>& ComputationalGeometry::getSolidImplicitFunction() const{
  return m_sol_if;
}

const Real& ComputationalGeometry::getGasPermittivity() const {
  return m_eps0;
}

const RefCountedPtr<mfis>& ComputationalGeometry::getMfIndexSpace() const {
  return m_multifluidIndexSpace;
}

void ComputationalGeometry::setDielectrics(Vector<Dielectric>& a_dielectrics){
  m_dielectrics = a_dielectrics;
}

void ComputationalGeometry::setElectrodes(Vector<Electrode>& a_electrodes){
  m_electrodes = a_electrodes;
}

void ComputationalGeometry::setGasPermittivity(const Real a_eps0){
  m_eps0 = a_eps0;
}

void ComputationalGeometry::buildGeometries(const ProblemDomain   a_finestDomain,
					    const RealVect        a_origin,
					    const Real            a_finestDx,
					    const int             a_nCellMax,
					    const int             a_maxCoarsen){

  // Build geoservers
  Vector<GeometryService*> geoservers(2, NULL);
  this->buildGasGeoServ(geoservers[phase::gas],     a_finestDomain, a_origin, a_finestDx);
  this->buildSolidGeoServ(geoservers[phase::solid], a_finestDomain, a_origin, a_finestDx);

  m_multifluidIndexSpace->define(a_finestDomain.domainBox(), // Define MF
				 a_origin,
				 a_finestDx,
				 geoservers,
				 a_nCellMax,
				 a_maxCoarsen);

  for (int i = 0; i < 2; i++){
    if(geoservers[i] != NULL){
      delete geoservers[i];
    }
  }
}

void ComputationalGeometry::buildGeometriesFromFiles(const std::string&   a_gas_file,
						     const std::string&   a_sol_file){

  RefCountedPtr<EBIndexSpace>& ebis_gas = m_multifluidIndexSpace->getEBIndexSpace(phase::gas);
  RefCountedPtr<EBIndexSpace>& ebis_sol = m_multifluidIndexSpace->getEBIndexSpace(phase::solid);

  // Define gas phase
  HDF5Handle gas_handle(a_gas_file.c_str(), HDF5Handle::OPEN_RDONLY);
  ebis_gas->define(gas_handle);
  gas_handle.close();

  if(m_dielectrics.size() > 0){
    HDF5Handle sol_handle(a_sol_file.c_str(), HDF5Handle::OPEN_RDONLY);
    ebis_sol->define(sol_handle);
    sol_handle.close();
  }
  else {
    ebis_sol = RefCountedPtr<EBIndexSpace> (NULL);
  }
}

void ComputationalGeometry::buildGasGeoServ(GeometryService*&   a_geoserver,
					    const ProblemDomain a_finestDomain,
					    const RealVect      a_origin,
					    const Real          a_dx){
  
  // The gas ebis is the intersection of electrodes and dielectrics
  Vector<BaseIF*> parts;
  for (int i = 0; i < m_dielectrics.size(); i++){
    parts.push_back(&(*(m_dielectrics[i].getImplicitFunction())));
  }
  for (int i = 0; i < m_electrodes.size(); i++){
    parts.push_back(&(*(m_electrodes[i].getImplicitFunction())));
  }

  // Create geoshop; either all regular or an intersection
  // if(parts.size() == 0){ 
  //   a_geoserver = new AllRegularService();
  // }
  // else {
  m_gas_if = RefCountedPtr<BaseIF> (new IntersectionIF(parts));
  if(s_use_new_gshop){
    a_geoserver = static_cast<GeometryService*> (new ScanShop(*m_gas_if,
							      0,
							      a_dx,
							      a_origin,
							      a_finestDomain,
							      s_ScanDomain,
							      s_thresh));
  }
  else{ // Chombo geometry generation
    a_geoserver = static_cast<GeometryService*> (new GeometryShop(*m_gas_if, 0, a_dx*RealVect::Unit, s_thresh));
  }
  //  }
}

void ComputationalGeometry::buildSolidGeoServ(GeometryService*&   a_geoserver,
					      const ProblemDomain a_finestDomain,
					      const RealVect      a_origin,
					      const Real          a_dx){

  // The gas ebis is the intersection of dielectrics
  Vector<BaseIF*> parts_dielectrics, parts_electrodes;
  for (int i = 0; i < m_dielectrics.size(); i++){
    parts_dielectrics.push_back(&(*m_dielectrics[i].getImplicitFunction()));
  }
  for (int i = 0; i < m_electrodes.size(); i++){
    parts_electrodes.push_back(&(*m_electrodes[i].getImplicitFunction()));
  }

  // Create EBIndexSpace. If there are no solid phases, return null
  if(parts_dielectrics.size() == 0){
    a_geoserver = NULL;
  }
  else {
    Vector<BaseIF*> parts;

    RefCountedPtr<BaseIF> diel_baseif = RefCountedPtr<BaseIF> (new IntersectionIF(parts_dielectrics)); // Gas outside dielectric
    RefCountedPtr<BaseIF> elec_baseif = RefCountedPtr<BaseIF> (new IntersectionIF(parts_electrodes));  // Gas outside electrodes
    RefCountedPtr<BaseIF> diel_compif = RefCountedPtr<BaseIF> (new ComplementIF(*diel_baseif));        // Inside of dielectrics

    parts.push_back(&(*diel_compif)); // Parts for intersection
    parts.push_back(&(*elec_baseif)); // Parts for intersection
    
    m_sol_if = RefCountedPtr<BaseIF> (new IntersectionIF(parts)); 

    if(s_use_new_gshop){
      a_geoserver = static_cast<GeometryService*> (new ScanShop(*m_sol_if,
								0,
								a_dx,
								a_origin,
								a_finestDomain,
								s_ScanDomain,
								s_thresh));
    }
    else{ // Chombo geometry generation
      a_geoserver = static_cast<GeometryService*> (new GeometryShop(*m_sol_if, 0, a_dx*RealVect::Unit, s_thresh));
    }
  }
}

#include <CD_NamespaceFooter.H>
