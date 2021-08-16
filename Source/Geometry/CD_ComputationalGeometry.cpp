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
  
  m_multifluidIndexSpace = RefCountedPtr<MultiFluidIndexSpace> (new MultiFluidIndexSpace());
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

const RefCountedPtr<BaseIF>& ComputationalGeometry::getImplicitFunction(const phase::which_phase a_phase) const {
  if(a_phase == phase::gas){
    return m_gas_if;
  }
  else if(a_phase == phase::solid){
    return m_sol_if;
  }
}

Real ComputationalGeometry::getGasPermittivity() const {
  return m_eps0;
}

const RefCountedPtr<MultiFluidIndexSpace>& ComputationalGeometry::getMfIndexSpace() const {
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

std::pair<Real, Real> ComputationalGeometry::getPrincipalCurvatures(const phase::which_phase a_phase, const RealVect a_pos, const Real a_dx) const {

  const RefCountedPtr<BaseIF>& impFunc = (a_phase == phase::gas) ? m_gas_if : m_sol_if;

  std::pair<Real, Real> ret;
#if CH_SPACEDIM==2
  const Real curv = this->curvature2D(impFunc, a_pos, a_dx);
  ret = std::make_pair(curv, curv);
#elif CH_SPACEDIM==3
  ret = this->curvature3D(impFunc, a_pos, a_dx);
#endif

  return ret;
}

#if CH_SPACEDIM==2
Real ComputationalGeometry::curvature2D(const RefCountedPtr<BaseIF>& a_F, const RealVect a_pos, const Real a_diffDx) const {
  const Real dx      = a_diffDx;

  // Points for centered differences
  const RealVect x0  = a_pos;
  const RealVect xLo = x0 - a_diffDx*BASISREALV(0);
  const RealVect xHi = x0 + a_diffDx*BASISREALV(0);
  const RealVect yLo = x0 - a_diffDx*BASISREALV(1);
  const RealVect yHi = x0 + a_diffDx*BASISREALV(1);

  // Points for mixed derivatives
  const RealVect xll = x0 - a_diffDx*BASISREALV(0) - a_diffDx*BASISREALV(1);
  const RealVect xhh = x0 + a_diffDx*BASISREALV(0) + a_diffDx*BASISREALV(1);
  const RealVect xlh = x0 - a_diffDx*BASISREALV(0) + a_diffDx*BASISREALV(1);
  const RealVect xhl = x0 + a_diffDx*BASISREALV(0) - a_diffDx*BASISREALV(1);

  // Derivs to O(dx^2)
  const Real Fx  = (a_F->value(xHi) - a_F->value(xLo))/(2*dx);
  const Real Fy  = (a_F->value(yHi) - a_F->value(yLo))/(2*dx);

  // Second derivs to O(dx^2)
  const Real Fxx = (a_F->value(xHi) - 2*a_F->value(x0) + a_F->value(xLo))/(dx*dx);
  const Real Fyy = (a_F->value(yHi) - 2*a_F->value(x0) + a_F->value(yLo))/(dx*dx);

  // Mixed deriv to O(dx^2)
  const Real Fxy = (a_F->value(xll) + a_F->value(xhh) - a_F->value(xlh) - a_F->value(xhl))/(4.*dx*dx);

  const Real curv = (-Fy*Fy*Fxx + 2*Fx*Fy*Fxy - Fx*Fx*Fyy)/std::pow(Fx*Fx + Fy*Fy, 3./2.);

  return curv;
}
#endif

#if CH_SPACEDIM==3
std::pair<Real, Real> ComputationalGeometry::curvature3D(const RefCountedPtr<BaseIF>& a_impFunc, const RealVect a_pos, const Real a_diffDx) const {

  // Some shortcuts
  auto phi = [a_impFunc](const RealVect x) -> Real {
    return a_impFunc->value(x);
  };

  const Real d      = a_diffDx;
  const RealVect  x = a_pos;
  const RealVect dx = a_diffDx*BASISREALV(0);
  const RealVect dy = a_diffDx*BASISREALV(1);
  const RealVect dz = a_diffDx*BASISREALV(2);

  // Gradient terms to O(dx^2)
  const Real Fx = (phi(x + dx) - phi(x-dx))/(2*d);
  const Real Fy = (phi(x + dy) - phi(x-dy))/(2*d);
  const Real Fz = (phi(x + dz) - phi(x-dz))/(2*d);

  // Second order derivs to O(dx^2)
  const Real Fxx = (phi(x + dx) - 2*phi(x) + phi(x-dx))/(d*d);
  const Real Fyy = (phi(x + dy) - 2*phi(x) + phi(x-dy))/(d*d);
  const Real Fzz = (phi(x + dz) - 2*phi(x) + phi(x-dz))/(d*d);

  // Mixed derivs
  const Real Fxy = (phi(x+dx+dy) + phi(x-dx-dy) - phi(x+dx-dy) - phi(x-dx+dy))/(4.*d*d);
  const Real Fxz = (phi(x+dx+dz) + phi(x-dx-dz) - phi(x+dx-dz) - phi(x-dx+dz))/(4.*d*d);
  const Real Fyz = (phi(x+dy+dz) + phi(x-dy-dz) - phi(x+dy-dz) - phi(x-dy+dz))/(4.*d*d);

  // Expressions taken from Albin et. al. "Computational assessment of curvatures and principaldirections of implicit surfaces from 3D scalar data"
  const Real denomH = std::pow(Fx*Fx + Fy*Fy + Fz*Fz, 3./2.);
  const Real denomK = std::pow(Fx*Fx + Fy*Fy + Fz*Fz, 2.0);
  const Real kH = (1./(2.*denomH)) * (Fx*Fx*(Fyy + Fzz) + Fy*Fy*(Fxx + Fzz) + Fz*Fz*(Fxx + Fyy))
    - (1./denomH)*(Fx*Fy*Fxy + Fx*Fz*Fxz + Fy*Fz*Fyz);
  const Real kK = (2./denomK)*(Fx*Fy*(Fxz*Fyz - Fxy*Fzz) + Fx*Fz*(Fxy*Fyz - Fxz*Fyy) + Fy*Fz*(Fxy*Fxz - Fyz*Fxx))
    + (1./denomK)*(Fx*Fx*(Fyy*Fzz - Fyz*Fyz) + Fy*Fy*(Fxx*Fzz - Fxz*Fxz) + Fzz*(Fxx*Fyy - Fxy*Fxy));
  

  // Mixed derives
  const Real kmin = kH - sqrt(std::abs(kH*kH - kK));
  const Real kmax = kH + sqrt(std::abs(kH*kH - kK));
  
  return std::make_pair(kmin, kmax);
}
#endif

#include <CD_NamespaceFooter.H>
