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
#include <CD_NewIntersectionIF.H>
#include <CD_ScanShop.H>
#include <CD_MemoryReport.H>
#include <CD_NamespaceHeader.H>

ComputationalGeometry::ComputationalGeometry()
{
  CH_TIME("ComputationalGeometry::ComputationalGeometry()");

  // Default parameters.
  m_eps0 = 1.0;

  m_electrodes.resize(0);
  m_dielectrics.resize(0);

  m_useScanShop = false;
  m_scanDomain  = ProblemDomain();

  m_multifluidIndexSpace = RefCountedPtr<MultiFluidIndexSpace>(new MultiFluidIndexSpace());
}

ComputationalGeometry::~ComputationalGeometry()
{
  CH_TIME("ComputationalGeometry::~ComputationalGeometry()");
}

void
ComputationalGeometry::useScanShop(const ProblemDomain a_beginDomain)
{
  CH_TIME("ComputationalGeometry::useScanShop(ProblemDomain)");

  // TLDR: If you called this function you signal that ComputationalGeometry will use ScanShop for geometry generation.

  m_useScanShop = true;
  m_scanDomain  = a_beginDomain;
}

void
ComputationalGeometry::useChomboShop()
{
  CH_TIME("ComputationalGeometry::useChomboShop()");

  // TLDR: If you called this function you signal that ComputationalGeometry will use Chombo's GeometryShop for geometry generation.
  m_useScanShop = false;
  m_scanDomain  = ProblemDomain();
}

const Vector<Dielectric>&
ComputationalGeometry::getDielectrics() const
{
  CH_TIME("ComputationalGeometry::getDielectrics()");

  return (m_dielectrics);
}

const Vector<Electrode>&
ComputationalGeometry::getElectrodes() const
{
  CH_TIME("ComputationalGeometry::getElectrodes()");

  return (m_electrodes);
}

const RefCountedPtr<BaseIF>&
ComputationalGeometry::getGasImplicitFunction() const
{
  CH_TIME("ComputationalGeometry::getGasImplicitFunction()");

  return (m_implicitFunctionGas);
}

const RefCountedPtr<BaseIF>&
ComputationalGeometry::getSolidImplicitFunction() const
{
  CH_TIME("ComputationalGeometry::getSolidImplicitFunction()");

  return (m_implicitFunctionSolid);
}

const RefCountedPtr<BaseIF>&
ComputationalGeometry::getImplicitFunction(const phase::which_phase a_phase) const
{
  CH_TIME("ComputationalGeometry::getImplicitFunction(phase::which_phase)");

  return (a_phase == phase::gas) ? m_implicitFunctionGas : m_implicitFunctionSolid;
}

Real
ComputationalGeometry::getGasPermittivity() const
{
  CH_TIME("ComputationalGeometry::getGasPermittivity()");

  return (m_eps0);
}

const RefCountedPtr<MultiFluidIndexSpace>&
ComputationalGeometry::getMfIndexSpace() const
{
  CH_TIME("ComputationalGeometry::getMfIndexSpace()");

  return (m_multifluidIndexSpace);
}

void
ComputationalGeometry::setDielectrics(const Vector<Dielectric>& a_dielectrics)
{
  CH_TIME("ComputationalGeometry::setDielectrics(Vector<Dielectric>)");

  m_dielectrics = a_dielectrics;
}

void
ComputationalGeometry::setElectrodes(const Vector<Electrode>& a_electrodes)
{
  CH_TIME("ComputationalGeometry::setElectrodes(Vector<Electrode>)");

  m_electrodes = a_electrodes;
}

void
ComputationalGeometry::setGasPermittivity(const Real a_eps0)
{
  CH_TIME("ComputationalGeometry::setGasPermittivity(Real)");

  m_eps0 = a_eps0;
}

void
ComputationalGeometry::buildGeometries(const ProblemDomain a_finestDomain,
                                       const RealVect      a_probLo,
                                       const Real          a_finestDx,
                                       const int           a_nCellMax,
                                       const int           a_maxGhostEB,
                                       const int           a_maxCoarsen)
{
  CH_TIME("ComputationalGeometry::buildGeometries(ProblemDomain, RealVect, Real, int, int, int)");

  // Set the default maximum number of EB ghosts that we will ever use. This is needed because ScanShop will look
  // through grown grid patches when it determines if a grid patch is irregular or not.
  m_maxGhostEB = a_maxGhostEB;

  // Build the geoservers. This creates the composite implicit functions and the GeometryService* objects which
  // can be passed to Chombo. Note that the
  Vector<GeometryService*> geoServices(2, nullptr);

  this->buildGasGeometry(geoServices[phase::gas], a_finestDomain, a_probLo, a_finestDx);
  this->buildSolidGeometry(geoServices[phase::solid], a_finestDomain, a_probLo, a_finestDx);

  // Define the multifluid index space.
  const bool useDistributedData = m_useScanShop;

  m_multifluidIndexSpace->define(a_finestDomain.domainBox(), // Define MF
                                 a_probLo,
                                 a_finestDx,
                                 geoServices,
                                 useDistributedData,
                                 a_nCellMax,
                                 a_maxCoarsen);

  // Delete temps.
  for (int i = 0; i < 2; i++) {
    if (geoServices[i] != nullptr) {
      delete geoServices[i];
    }
  }
}

void
ComputationalGeometry::buildGasGeometry(GeometryService*&   a_geoserver,
                                        const ProblemDomain a_finestDomain,
                                        const RealVect      a_probLo,
                                        const Real          a_finestDx)
{
  CH_TIME("ComputationalGeometry::buildGasGeometry(GeometryService, ProblemDomain, RealVect, Real)");

  // The gas phase is the intersection of the region outside every object, so IntersectionIF is correct here. We build the
  // various parts and then create the implicit function for the gas-phas using constructive solid geometry.
  Vector<BaseIF*> parts;
  for (int i = 0; i < m_dielectrics.size(); i++) {
    parts.push_back(&(*(m_dielectrics[i].getImplicitFunction())));
  }
  for (int i = 0; i < m_electrodes.size(); i++) {
    parts.push_back(&(*(m_electrodes[i].getImplicitFunction())));
  }

  m_implicitFunctionGas = RefCountedPtr<BaseIF>(new NewIntersectionIF(parts));

  // Build the EBIS geometry. Use either ScanShop or Chombo here.
  if (m_useScanShop) {
    ScanShop* scanShop = new ScanShop(*m_implicitFunctionGas,
                                      0,
                                      a_finestDx,
                                      a_probLo,
                                      a_finestDomain,
                                      m_scanDomain,
                                      m_maxGhostEB,
                                      s_thresh);

    scanShop->setProfileFileName("ScanShopReportGasPhase.dat");

    a_geoserver = static_cast<GeometryService*>(scanShop);
  }
  else { // Chombo geometry generation
    a_geoserver = static_cast<GeometryService*>(
      new GeometryShop(*m_implicitFunctionGas, 0, a_finestDx * RealVect::Unit, s_thresh));
  }
}

void
ComputationalGeometry::buildSolidGeometry(GeometryService*&   a_geoserver,
                                          const ProblemDomain a_finestDomain,
                                          const RealVect      a_probLo,
                                          const Real          a_finestDx)
{
  CH_TIME("ComputationalGeometry::buildSolidGeometry(GeometryService, ProblemDomain, RealVect, Real)");

  // The "solid phase", i.e. the part inside dielectrics is a bit more complicated to compute. We want to get the region
  // outside the electrodes but inside the dielectrics. Fortunately there is a way to do this.

  // Get all the parts (dielectrics/electrodes)
  Vector<BaseIF*> dielectricParts;
  Vector<BaseIF*> electrodeParts;

  for (int i = 0; i < m_dielectrics.size(); i++) {
    dielectricParts.push_back(&(*m_dielectrics[i].getImplicitFunction()));
  }

  for (int i = 0; i < m_electrodes.size(); i++) {
    electrodeParts.push_back(&(*m_electrodes[i].getImplicitFunction()));
  }

  // Create EBIndexSpace. If there are no solid phases, return null
  if (dielectricParts.size() == 0) {
    a_geoserver = nullptr;
  }
  else {
    Vector<BaseIF*> parts;

    RefCountedPtr<BaseIF> dielBaseIF = RefCountedPtr<BaseIF>(
      new NewIntersectionIF(dielectricParts)); // This gives the region outside the dielectrics.
    RefCountedPtr<BaseIF> elecBaseIF = RefCountedPtr<BaseIF>(
      new NewIntersectionIF(electrodeParts)); // This is the region outside the the electrodes.
    RefCountedPtr<BaseIF> dielCompIF = RefCountedPtr<BaseIF>(
      new ComplementIF(*dielBaseIF)); // This is the region inside the dielectrics.

    // We want the function which is the region inside the dielectrics and outside the electrodes, i.e. the intersection of the region "inside"
    // dielectrics and outside the electrods.
    parts.push_back(&(*dielCompIF));
    parts.push_back(&(*elecBaseIF));

    m_implicitFunctionSolid = RefCountedPtr<BaseIF>(new IntersectionIF(parts));

    // Build the EBIS geometry. Use either ScanShop or Chombo here.
    if (m_useScanShop) {
      ScanShop* scanShop = new ScanShop(*m_implicitFunctionSolid,
                                        0,
                                        a_finestDx,
                                        a_probLo,
                                        a_finestDomain,
                                        m_scanDomain,
                                        m_maxGhostEB,
                                        s_thresh);

      scanShop->setProfileFileName("ScanShopReportSolidPhase.dat");

      a_geoserver = static_cast<GeometryService*>(scanShop);
    }
    else { // Chombo geometry generation
      a_geoserver = static_cast<GeometryService*>(
        new GeometryShop(*m_implicitFunctionSolid, 0, a_finestDx * RealVect::Unit, s_thresh));
    }
  }
}

#include <CD_NamespaceFooter.H>
