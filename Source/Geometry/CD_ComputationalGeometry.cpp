/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
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
#include <CD_MemoryReport.H>
#include <CD_NamespaceHeader.H>

//Real ComputationalGeometry::s_thresh;

ComputationalGeometry::ComputationalGeometry(){
  CH_TIME("ComputationalGeometry::ComputationalGeometry()");
  
  // Default parameters. 
  m_eps0 = 1.0;
  
  m_electrodes. resize(0);
  m_dielectrics.resize(0);

  m_useScanShop = false;
  m_scanDomain  = ProblemDomain();
  
  m_multifluidIndexSpace = RefCountedPtr<MultiFluidIndexSpace> (new MultiFluidIndexSpace());
}

ComputationalGeometry::~ComputationalGeometry(){
  CH_TIME("ComputationalGeometry::~ComputationalGeometry()");
}

void ComputationalGeometry::useScanShop(const ProblemDomain a_beginDomain){
  CH_TIME("ComputationalGeometry::useScanShop(ProblemDomain)");

  // TLDR: If you called this function you signal that ComputationalGeometry will use ScanShop for geometry generation.
  
  m_useScanShop = true;
  m_scanDomain  = a_beginDomain;
}

void ComputationalGeometry::useChomboShop(){
  CH_TIME("ComputationalGeometry::useChomboShop()");

  // TLDR: If you called this function you signal that ComputationalGeometry will use Chombo's GeometryShop for geometry generation.
  m_useScanShop = false;
  m_scanDomain  = ProblemDomain();
}

const Vector<Dielectric>& ComputationalGeometry::getDielectrics() const  {
  CH_TIME("ComputationalGeometry::getDielectrics()");
  
  return (m_dielectrics);
}

const Vector<Electrode>& ComputationalGeometry::getElectrodes() const  {
  CH_TIME("ComputationalGeometry::getElectrodes()");
  
  return (m_electrodes);
}

const RefCountedPtr<BaseIF>& ComputationalGeometry::getGasImplicitFunction() const{
  CH_TIME("ComputationalGeometry::getGasImplicitFunction()");
  
  return (m_implicitFunctionGas);
}

const RefCountedPtr<BaseIF>& ComputationalGeometry::getSolidImplicitFunction() const{
  CH_TIME("ComputationalGeometry::getSolidImplicitFunction()");
  
  return (m_implicitFunctionSolid);
}

const RefCountedPtr<BaseIF>& ComputationalGeometry::getImplicitFunction(const phase::which_phase a_phase) const {
  CH_TIME("ComputationalGeometry::getImplicitFunction(phase::which_phase)");

  if(a_phase == phase::gas){
    return m_implicitFunctionGas;
  }
  else if(a_phase == phase::solid){
    return m_implicitFunctionSolid;
  }
}

Real ComputationalGeometry::getGasPermittivity() const {
  CH_TIME("ComputationalGeometry::getGasPermittivity()");
  
  return (m_eps0);
}

const RefCountedPtr<MultiFluidIndexSpace>& ComputationalGeometry::getMfIndexSpace() const {
  CH_TIME("ComputationalGeometry::getMfIndexSpace()");
  
  return (m_multifluidIndexSpace);
}

void ComputationalGeometry::setDielectrics(const Vector<Dielectric>& a_dielectrics){
  CH_TIME("ComputationalGeometry::setDielectrics(Vector<Dielectric>)");
  
  m_dielectrics = a_dielectrics;
}

void ComputationalGeometry::setElectrodes(const Vector<Electrode>& a_electrodes){
  CH_TIME("ComputationalGeometry::setElectrodes(Vector<Electrode>)");
  
  m_electrodes = a_electrodes;
}

void ComputationalGeometry::setGasPermittivity(const Real a_eps0){
  CH_TIME("ComputationalGeometry::setGasPermittivity(Real)");
  
  m_eps0 = a_eps0;
}

void ComputationalGeometry::buildGeometries(const ProblemDomain   a_finestDomain,
					    const RealVect        a_probLo,
					    const Real            a_finestDx,
					    const int             a_nCellMax,
					    const int             a_maxGhostEB,
					    const int             a_maxCoarsen){
  CH_TIME("ComputationalGeometry::buildGeometries(ProblemDomain, RealVect, Real, int, int, int)");  

  // Set the default maximum number of EB ghosts that we will ever use. This is needed because ScanShop will look
  // through grown grid patches when it determines if a grid patch is irregular or not. 
  m_maxGhostEB = a_maxGhostEB;
  
  // Build the geoservers. 
  Vector<GeometryService*> geoservers(2, NULL);
  
  this->buildGasGeoServ  (geoservers[phase::gas],   a_finestDomain, a_probLo, a_finestDx);
  this->buildSolidGeoServ(geoservers[phase::solid], a_finestDomain, a_probLo, a_finestDx);

  // Define the multifluid index space.
  const bool useDistributedData = m_useScanShop;
  
  m_multifluidIndexSpace->define(a_finestDomain.domainBox(), // Define MF
				 a_probLo,
				 a_finestDx,
				 geoservers,
				 useDistributedData,
				 a_nCellMax,
				 a_maxCoarsen);

  // Delete temps. 
  for (int i = 0; i < 2; i++){
    if(geoservers[i] != NULL){
      delete geoservers[i];
    }
  }
}

void ComputationalGeometry::buildGasGeoServ(GeometryService*&   a_geoserver,
					    const ProblemDomain a_finestDomain,
					    const RealVect      a_probLo,
					    const Real          a_dx){
  CH_TIME("ComputationalGeometry::buildGasGeoServ(GeometryService, ProblemDomain, RealVect, Real)");
  
  // The gas phase is the intersection of the region outside every object, so IntersectionIF is correct here. We build the
  // various parts and then create the implicit function for the gas-phas using constructive solid geometry.
  Vector<BaseIF*> parts;
  for (int i = 0; i < m_dielectrics.size(); i++){
    parts.push_back(&(*(m_dielectrics[i].getImplicitFunction())));
  }
  for (int i = 0; i < m_electrodes.size(); i++){
    parts.push_back(&(*(m_electrodes[i].getImplicitFunction())));
  }

  m_implicitFunctionGas = RefCountedPtr<BaseIF> (new IntersectionIF(parts));

  // Build the EBIS geometry. Use either ScanShop or Chombo here. 
  if(m_useScanShop){
    a_geoserver = static_cast<GeometryService*> (new ScanShop(*m_implicitFunctionGas,
							      0,
							      a_dx,
							      a_probLo,
							      a_finestDomain,
							      m_scanDomain,
							      m_maxGhostEB,
							      s_thresh));
  }
  else{ // Chombo geometry generation
    a_geoserver = static_cast<GeometryService*> (new GeometryShop(*m_implicitFunctionGas, 0, a_dx*RealVect::Unit, s_thresh));
  }
}

void ComputationalGeometry::buildSolidGeoServ(GeometryService*&   a_geoserver,
					      const ProblemDomain a_finestDomain,
					      const RealVect      a_probLo,
					      const Real          a_dx){
  CH_TIME("ComputationalGeometry::buildSolidGeoServ(GeometryService, ProblemDomain, RealVect, Real)");  
  
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
    
    m_implicitFunctionSolid = RefCountedPtr<BaseIF> (new IntersectionIF(parts)); 

    if(m_useScanShop){
      a_geoserver = static_cast<GeometryService*> (new ScanShop(*m_implicitFunctionSolid,
								0,
								a_dx,
								a_probLo,
								a_finestDomain,
								m_scanDomain,
								m_maxGhostEB,
								s_thresh));
    }
    else{ // Chombo geometry generation
      a_geoserver = static_cast<GeometryService*> (new GeometryShop(*m_implicitFunctionSolid, 0, a_dx*RealVect::Unit, s_thresh));
    }
  }
}

#include <CD_NamespaceFooter.H>
