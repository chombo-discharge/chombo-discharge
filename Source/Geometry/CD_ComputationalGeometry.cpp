/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_ComputationalGeometry.cpp
 * @brief  Implementation of CD_ComputationalGeometry.H
 * @author Robert Marskar
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
#include <CD_TiledMeshRefine.H>
#include <CD_Units.H>

// Std includes
#include <limits>
#include <CD_ComputationalGeometry.H>
#include <CD_NewIntersectionIF.H>
#include <CD_ScanShop.H>
#include <CD_PolyhedralGeometryShop.H>
#include <CD_MemoryReport.H>
#include <CD_NamespaceHeader.H>

ComputationalGeometry::ComputationalGeometry()
  : m_eps0(1.0), m_generator(Generator::GeometryShop), m_geometryRefinement(1)
{
  CH_TIME("ComputationalGeometry::ComputationalGeometry()");

  // Default parameters.

  m_electrodes.resize(0);
  m_dielectrics.resize(0);

  m_scanDomain = ProblemDomain();

  m_multifluidIndexSpace = RefCountedPtr<MultiFluidIndexSpace>(new MultiFluidIndexSpace());
}

ComputationalGeometry::~ComputationalGeometry()
{
  CH_TIME("ComputationalGeometry::~ComputationalGeometry()");
}

void
ComputationalGeometry::useScanShop(const ProblemDomain& a_beginDomain)
{
  CH_TIME("ComputationalGeometry::useScanShop(ProblemDomain)");

  // TLDR: If you called this function you signal that ComputationalGeometry will use ScanShop for geometry generation.

  m_generator  = Generator::ScanShop;
  m_scanDomain = a_beginDomain;
}

void
ComputationalGeometry::usePolyhedralShop(const ProblemDomain& a_beginDomain, const int a_refinement)
{
  CH_TIME("ComputationalGeometry::usePolyhedralShop(ProblemDomain, int)");

  m_generator          = Generator::PolyhedralShop;
  m_scanDomain         = a_beginDomain;
  m_geometryRefinement = a_refinement;
}

void
ComputationalGeometry::useChomboShop()
{
  CH_TIME("ComputationalGeometry::useChomboShop()");

  // TLDR: If you called this function you signal that ComputationalGeometry will use Chombo's GeometryShop for geometry
  // generation.
  m_generator  = Generator::GeometryShop;
  m_scanDomain = ProblemDomain();
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
ComputationalGeometry::buildGeometries(const ProblemDomain& a_finestDomain,
                                       const RealVect&      a_probLo,
                                       const Real           a_finestDx,
                                       const int            a_nCellMax,
                                       const int            a_maxGhostEB,
                                       const int            a_maxCoarsen)
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
  const bool useDistributedData = (m_generator != Generator::GeometryShop);

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

Real
ComputationalGeometry::edgeRoot(const RefCountedPtr<BaseIF>& a_implicitFunction,
                                const RealVect&              a_lowPoint,
                                const int                    a_dir,
                                const Real&                  a_lowValue,
                                const Real&                  a_dx) const
{
  Real lo       = 0.0;
  Real hi       = 1.0;
  Real lowValue = a_lowValue;

  for (int iter = 0; iter < 100; iter++) {
    const Real mid = 0.5 * (lo + hi);

    RealVect x = a_lowPoint;
    x[a_dir] += a_dx * mid;

    const Real value = a_implicitFunction->value(x);

    if (PolyhedralEB::isFluid(value) == PolyhedralEB::isFluid(lowValue)) {
      lo       = mid;
      lowValue = value;
    }
    else {
      hi = mid;
    }

    if (hi - lo < 1.0E-15) {
      break;
    }
  }

  return 0.5 * (lo + hi);
}

void
ComputationalGeometry::tagUnderResolvedCells(IntVectSet&                  a_tags,
                                             const Box&                   a_region,
                                             const ProblemDomain&         a_domain,
                                             const RefCountedPtr<BaseIF>& a_implicitFunction,
                                             const RealVect&              a_probLo,
                                             const Real&                  a_dx,
                                             const Real&                  a_angle) const
{
  CH_TIME("ComputationalGeometry::tagUnderResolvedCells");

  // A cell is compared against every neighbour, so normals are wanted one cell out from the
  // region the tags are for.
  const Box grownRegion = grow(a_region, 1) & a_domain;

  Box nodeBox = grownRegion;
  nodeBox.surroundingNodes();

  BaseFab<Real> nodeValues(nodeBox, 1);

  for (BoxIterator bit(nodeBox); bit.ok(); ++bit) {
    const IntVect node = bit();

    RealVect x = a_probLo;

    for (int d = 0; d < SpaceDim; d++) {
      x[d] += a_dx * static_cast<Real>(node[d]);
    }

    nodeValues(node, 0) = a_implicitFunction->value(x);
  }

  // Edges are shared by the cells meeting along them, so each root is found once and addressed
  // by the edge's own low node.
  BaseFab<Real> intercept[SpaceDim];

  for (int dir = 0; dir < SpaceDim; dir++) {
    Box edgeBox = nodeBox;
    edgeBox.enclosedCells(dir);

    intercept[dir].define(edgeBox, 1);
    intercept[dir].setVal(PolyhedralEB::CutCellSurface::s_noCrossing);
  }

  BaseFab<Real> normal(grownRegion, SpaceDim);
  BaseFab<int>  isCut(grownRegion, 1);

  normal.setVal(0.0);
  isCut.setVal(0);

  for (BoxIterator bit(grownRegion); bit.ok(); ++bit) {
    const IntVect iv = bit();

    PolyhedralEB::CutCellSurface surface;

    bool cut = false;

    for (int c = 0; c < PolyhedralEB::CutCellSurface::s_numCorners; c++) {
      IntVect node = iv;

      for (int d = 0; d < SpaceDim; d++) {
        node[d] += (c >> d) & 1;
      }

      surface.m_corner[c] = nodeValues(node, 0);

      cut = cut || (PolyhedralEB::isFluid(surface.m_corner[c]) != PolyhedralEB::isFluid(surface.m_corner[0]));
    }

    if (!cut) {
      continue;
    }

    for (int e = 0; e < PolyhedralEB::CutCellSurface::s_numEdges; e++) {
      int lo = 0;
      int hi = 0;

      PolyhedralEB::detail::edgeCorners(e, lo, hi);

      if (PolyhedralEB::isFluid(surface.m_corner[lo]) == PolyhedralEB::isFluid(surface.m_corner[hi])) {
        continue;
      }

      const int dir = PolyhedralEB::detail::edgeDirection(e);

      int offset[SpaceDim];
      PolyhedralEB::detail::edgeOrigin(e, offset);

      IntVect edgeIV = iv;

      for (int d = 0; d < SpaceDim; d++) {
        edgeIV[d] += offset[d];
      }

      if (intercept[dir](edgeIV, 0) == PolyhedralEB::CutCellSurface::s_noCrossing) {
        RealVect lowPoint = a_probLo;

        for (int d = 0; d < SpaceDim; d++) {
          lowPoint[d] += a_dx * static_cast<Real>(edgeIV[d]);
        }

        intercept[dir](edgeIV, 0) = this->edgeRoot(a_implicitFunction, lowPoint, dir, surface.m_corner[lo], a_dx);
      }

      surface.m_crossing[e] = intercept[dir](edgeIV, 0);
    }

    const RealVect cellNormal = PolyhedralEB::crossingNormal(surface);

    if (cellNormal.vectorLength() <= 0.0) {
      continue;
    }

    isCut(iv, 0) = 1;

    for (int d = 0; d < SpaceDim; d++) {
      normal(iv, d) = cellNormal[d];
    }
  }

  // A feature thinner than a cell leaves every corner on the same side, so the cell reads as
  // regular and the bend test never sees it. Sampling inside the cell is what catches it, and
  // it is only worth sampling where the interface could plausibly reach: a cell whose smallest
  // corner value exceeds the largest change across any of its edges is further from the surface
  // than the function is seen to move over the cell. That bound comes from the data rather than
  // from assuming a signed distance function, which not all of them are.
  for (BoxIterator bit(a_region); bit.ok(); ++bit) {
    const IntVect iv = bit();

    if (isCut(iv, 0) == 1) {
      continue;
    }

    Real smallest = std::numeric_limits<Real>::max();
    Real spread   = 0.0;

    for (int c = 0; c < PolyhedralEB::CutCellSurface::s_numCorners; c++) {
      IntVect node = iv;

      for (int d = 0; d < SpaceDim; d++) {
        node[d] += (c >> d) & 1;
      }

      smallest = std::min(smallest, std::abs(nodeValues(node, 0)));

      for (int d = 0; d < SpaceDim; d++) {
        if (((c >> d) & 1) == 0) {
          IntVect other = node;
          other[d] += 1;

          spread = std::max(spread, std::abs(nodeValues(other, 0) - nodeValues(node, 0)));
        }
      }
    }

    if (smallest > spread) {
      continue;
    }

    const Real cornerSide = nodeValues(iv, 0);

    bool hidden = false;

    for (int e = 0; e < PolyhedralEB::CutCellSurface::s_numEdges && !hidden; e++) {
      const int dir = PolyhedralEB::detail::edgeDirection(e);

      int offset[SpaceDim];
      PolyhedralEB::detail::edgeOrigin(e, offset);

      RealVect x = a_probLo;

      for (int d = 0; d < SpaceDim; d++) {
        x[d] += a_dx * static_cast<Real>(iv[d] + offset[d]);
      }

      x[dir] += 0.5 * a_dx;

      hidden = PolyhedralEB::isFluid(a_implicitFunction->value(x)) != PolyhedralEB::isFluid(cornerSide);
    }

    if (!hidden) {
      RealVect centre = a_probLo;

      for (int d = 0; d < SpaceDim; d++) {
        centre[d] += a_dx * (static_cast<Real>(iv[d]) + 0.5);
      }

      hidden = PolyhedralEB::isFluid(a_implicitFunction->value(centre)) != PolyhedralEB::isFluid(cornerSide);
    }

    if (hidden) {
      a_tags |= iv;
    }
  }

  const Real cosineThreshold = std::cos(a_angle * Units::pi / 180.0);

  for (BoxIterator bit(a_region); bit.ok(); ++bit) {
    const IntVect iv = bit();

    if (isCut(iv, 0) == 0) {
      continue;
    }

    RealVect here;

    for (int d = 0; d < SpaceDim; d++) {
      here[d] = normal(iv, d);
    }

    bool bends = false;

    for (BoxIterator nit(Box(iv - IntVect::Unit, iv + IntVect::Unit)); nit.ok() && !bends; ++nit) {
      const IntVect other = nit();

      if (other == iv || !grownRegion.contains(other) || isCut(other, 0) == 0) {
        continue;
      }

      RealVect there;

      for (int d = 0; d < SpaceDim; d++) {
        there[d] = normal(other, d);
      }

      bends = here.dotProduct(there) < cosineThreshold;
    }

    if (bends) {
      a_tags |= iv;
    }
  }
}

Vector<IntVectSet>
ComputationalGeometry::getCurvatureTags(const ProblemDomain& a_coarsestDomain,
                                        const Vector<int>&   a_refRatios,
                                        const IntVect&       a_tileSize,
                                        const IntVect&       a_maxBlockSize,
                                        const RealVect&      a_probLo,
                                        const Real&          a_coarsestDx,
                                        const Real&          a_angle,
                                        const int            a_maxDepth) const
{
  CH_TIME("ComputationalGeometry::getCurvatureTags");

  Vector<IntVectSet> tags(std::max(a_maxDepth, 1));

  // The pre-pass descends on its own cap, which has nothing to do with how deep the run is
  // allowed to refine, so the ratios have to reach that far whatever the AMR hierarchy is.
  if (a_refRatios.size() < tags.size() - 1) {
    MayDay::Error("ComputationalGeometry::getCurvatureTags - too few refinement ratios for the requested depth");
  }

  Vector<ProblemDomain> domains(tags.size(), a_coarsestDomain);
  Vector<Real>          dx(tags.size(), a_coarsestDx);

  for (int lvl = 1; lvl < tags.size(); lvl++) {
    domains[lvl] = refine(domains[lvl - 1], a_refRatios[lvl - 1]);
    dx[lvl]      = dx[lvl - 1] / static_cast<Real>(a_refRatios[lvl - 1]);
  }

  const TiledMeshRefine meshRefine(a_coarsestDomain, a_refRatios, a_tileSize, a_maxBlockSize);

  RefCountedPtr<BaseIF> implicitFunctions[2];

  implicitFunctions[0] = m_implicitFunctionGas;
  implicitFunctions[1] = m_implicitFunctionSolid;

  // The coarsest level is swept whole; every finer one only where the level above it tagged.
  Vector<Vector<Box>> regions(tags.size());

  regions[0] = Vector<Box>(1, a_coarsestDomain.domainBox());

  for (int lvl = 0; lvl < tags.size(); lvl++) {
    for (const auto& implicitFunction : implicitFunctions) {
      if (implicitFunction.isNull()) {
        continue;
      }

      for (const auto& region : regions[lvl].stdVector()) {
        this->tagUnderResolvedCells(tags[lvl], region, domains[lvl], implicitFunction, a_probLo, dx[lvl], a_angle);
      }
    }

    if (lvl == tags.size() - 1) {
      break;
    }

    // The cap is what stops a crease, where the turn never falls off under refinement.
    Vector<Vector<Box>> grids;

    meshRefine.regrid(grids, tags);

    regions[lvl + 1] = (lvl + 1 < grids.size()) ? grids[lvl + 1] : Vector<Box>();
  }

  return tags;
}

void
ComputationalGeometry::buildGasGeometry(GeometryService*&    a_geoserver,
                                        const ProblemDomain& a_finestDomain,
                                        const RealVect&      a_probLo,
                                        const Real           a_finestDx)
{
  CH_TIME("ComputationalGeometry::buildGasGeometry(GeometryService, ProblemDomain, RealVect, Real)");

  // The gas phase is the intersection of the region outside every object, so IntersectionIF is correct here. We build
  // the various parts and then create the implicit function for the gas-phas using constructive solid geometry.
  Vector<BaseIF*> parts;
  for (int i = 0; i < m_dielectrics.size(); i++) {
    parts.push_back(&(*(m_dielectrics[i].getImplicitFunction())));
  }
  for (int i = 0; i < m_electrodes.size(); i++) {
    parts.push_back(&(*(m_electrodes[i].getImplicitFunction())));
  }

  m_implicitFunctionGas = RefCountedPtr<BaseIF>(new NewIntersectionIF(parts));

  // Build the EBIS geometry. Use ScanShop, the polyhedral generator, or Chombo here.
  if (m_generator == Generator::PolyhedralShop) {
    auto* shop = new PolyhedralGeometryShop(*m_implicitFunctionGas,
                                            0,
                                            a_finestDx,
                                            a_probLo,
                                            a_finestDomain,
                                            m_scanDomain,
                                            m_maxGhostEB,
                                            s_thresh,
                                            s_strictGeometry,
                                            m_geometryRefinement);

    shop->setProfileFileName("PolyhedralShopReportGasPhase.dat");

    a_geoserver = static_cast<GeometryService*>(shop);
  }
  else if (m_generator == Generator::ScanShop) {
    auto* scanShop = new ScanShop(*m_implicitFunctionGas,
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
ComputationalGeometry::buildSolidGeometry(GeometryService*&    a_geoserver,
                                          const ProblemDomain& a_finestDomain,
                                          const RealVect&      a_probLo,
                                          const Real           a_finestDx)
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

    // We want the function which is the region inside the dielectrics and outside the electrodes, i.e. the intersection
    // of the region "inside" dielectrics and outside the electrods.
    parts.push_back(&(*dielCompIF));
    parts.push_back(&(*elecBaseIF));

    m_implicitFunctionSolid = RefCountedPtr<BaseIF>(new IntersectionIF(parts));

    // Build the EBIS geometry. Use ScanShop, the polyhedral generator, or Chombo here.
    if (m_generator == Generator::PolyhedralShop) {
      auto* shop = new PolyhedralGeometryShop(*m_implicitFunctionSolid,
                                              0,
                                              a_finestDx,
                                              a_probLo,
                                              a_finestDomain,
                                              m_scanDomain,
                                              m_maxGhostEB,
                                              s_thresh,
                                              s_strictGeometry,
                                              m_geometryRefinement);

      shop->setProfileFileName("PolyhedralShopReportSolidPhase.dat");

      a_geoserver = static_cast<GeometryService*>(shop);
    }
    else if (m_generator == Generator::ScanShop) {
      auto* scanShop = new ScanShop(*m_implicitFunctionSolid,
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
