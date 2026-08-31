/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   main.cpp
 * @brief  Unit test for the mirrored-deposition surface data built by PhaseRealm at regrid.
 * @author Robert Marskar
 */

// Std includes
#include <cmath>

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_Driver.H>
#include <CD_GeometryStepper.H>
#include <CD_MirrorDeposition.H>
#include <CD_ParallelOps.H>
#include <CD_SphereSdf.H>
#include <CD_TorusSdf.H>

using namespace ChomboDischarge;
using namespace Physics::Geometry;

namespace {

/**
 * @brief Which analytic surface the test is checking against.
 */
enum class Shape
{
  Sphere,
  Torus
};

/**
 * @brief Test configuration, read once from the input file.
 */
Shape    g_shape         = Shape::Sphere;
Real     g_radius        = 0.4;
Real     g_minorRadius   = 0.15;
bool     g_fluidInside   = false;
RealVect g_centre        = RealVect::Zero;
Real     g_meanTolerance = 0.05;
Real     g_slimTolerance = 0.10;
int      g_numFailures   = 0;

/**
 * @brief Report one check, and remember whether it passed.
 * @param[in] a_what  Human-readable description.
 * @param[in] a_value Measured value.
 * @param[in] a_bound Largest permitted value.
 */
void
expect(const std::string& a_what, const Real a_value, const Real a_bound) noexcept
{
  const bool ok = (a_value <= a_bound);

  if (!ok) {
    g_numFailures++;
  }

  pout() << "  " << (ok ? "ok   " : "FAIL ") << a_what << " = " << a_value << " (bound " << a_bound << ")" << endl;
}

/**
 * @brief Analytic twice-mean-curvature and Gaussian curvature of the test surface at a point.
 * @details Signs follow the convention S = +grad(nhat) with nhat pointing INTO the fluid, so a convex solid has
 * positive mean curvature.
 * @param[in]  a_x    Position on (or near) the surface.
 * @param[out] a_twoH Twice the mean curvature.
 * @param[out] a_K    Gaussian curvature.
 */
void
analyticCurvature(const RealVect& a_x, Real& a_twoH, Real& a_K) noexcept
{
  const Real sign = g_fluidInside ? -1.0 : 1.0;

  switch (g_shape) {
  case Shape::Sphere: {
    // Both principal curvatures are 1/R.
    a_twoH = sign * (SpaceDim - 1) / g_radius;
    a_K    = (SpaceDim == 3) ? 1.0 / (g_radius * g_radius) : 0.0;

    break;
  }
  case Shape::Torus: {
    // Principal curvatures are 1/r (around the tube) and cos(theta)/(R + r*cos(theta)) (around the ring), with
    // theta measured from the equatorial plane and R + r*cos(theta) = rho, the distance from the axis. This is the
    // ONLY surface here that can distinguish a correct anisotropic shape operator from a wrong one: on a sphere
    // every shape operator is isotropic, and therefore frame-invariant, so a sphere cannot see the error.
    const RealVect d   = a_x - g_centre;
    const Real     rho = std::sqrt(d[0] * d[0] + d[1] * d[1]);

    const Real c1 = 1.0 / g_minorRadius;
    const Real c2 = (rho > 0.0) ? (rho - g_radius) / (g_minorRadius * rho) : 0.0;

    a_twoH = sign * (c1 + c2);
    a_K    = c1 * c2;

    break;
  }
  }
}
} // namespace

/**
 * @brief GeometryStepper that additionally registers the mirror-deposition operator and checks its output.
 */
class MirrorSurfaceDataStepper : public GeometryStepper
{
public:
  /**
   * @brief Register the operators this test needs.
   */
  void
  registerOperators() override
  {
    m_amr->registerOperator(s_mirror_deposition, Realm::Primal, phase::gas);
  }

  /**
   * @brief Check the built surface data against the analytic surface.
   */
  void
  postInitialize() override
  {
    this->checkSurfaceData();
  }

  /**
   * @brief Compare the fitted curvature with the analytic value, binned by volume fraction.
   */
  void
  checkSurfaceData() const noexcept;
};

void
MirrorSurfaceDataStepper::checkSurfaceData() const noexcept
{
  const EBAMRCellData& surfaceData = m_amr->getMirrorSurfaceData(Realm::Primal, phase::gas);
  const Vector<int>&   hasCutCells = m_amr->getMirrorHasCutCells(Realm::Primal, phase::gas);

  pout() << "MirrorSurfaceData: checking " << (g_shape == Shape::Sphere ? "sphere" : "torus") << endl;

  // Error accumulators. The 'slim' bin is kappa <= 0.05, reported separately because that is where the mirror
  // matters most and where the discrete normals are worst -- an aggregate hides exactly the cells of interest.
  Real      worstAll   = 0.0;
  Real      worstSlim  = 0.0;
  Real      sumAll     = 0.0;
  Real      worstNormK = 0.0;
  Real      worstSn    = 0.0;
  long long numAll     = 0;
  long long numSlim    = 0;
  long long numFitted  = 0;
  long long numPlanar  = 0;
  long long numNone    = 0;
  long long numAniso   = 0;

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl   = m_amr->getGrids(Realm::Primal)[lvl];
    const EBISLayout&        ebisl = m_amr->getEBISLayout(Realm::Primal, phase::gas)[lvl];
    const DataIterator&      dit   = dbl.dataIterator();
    const Real               dx    = m_amr->getDx()[lvl];

    for (DataIterator it = dit; it.ok(); ++it) {
      const EBCellFAB&     surf    = (*surfaceData[lvl])[it()];
      const BaseFab<Real>& surfReg = surf.getSingleValuedFAB();
      const EBISBox&       ebisbox = surf.getEBISBox();

      for (IVSIterator ivsIt(ebisbox.getIrregIVS(dbl[it()])); ivsIt.ok(); ++ivsIt) {
        const IntVect& iv = ivsIt();

        if (ebisbox.numVoFs(iv) != 1) {
          continue;
        }

        const VolIndex vof    = VolIndex(iv, 0);
        const Real     kappa  = ebisbox.volFrac(vof);
        const int      status = static_cast<int>(std::lround(surfReg(iv, MirrorDeposition::compStatus)));

        if (status == static_cast<int>(MirrorDeposition::Status::None)) {
          numNone++;

          continue;
        }
        if (status == static_cast<int>(MirrorDeposition::Status::Planar)) {
          numPlanar++;

          continue;
        }

        numFitted++;

        RealVect xc;
        RealVect nc;
        for (int dir = 0; dir < SpaceDim; dir++) {
          xc[dir] = surfReg(iv, MirrorDeposition::compCentroid + dir);
          nc[dir] = surfReg(iv, MirrorDeposition::compNormal + dir);
        }

        Real shape[MirrorDeposition::numShapeComp];
        for (int c = 0; c < MirrorDeposition::numShapeComp; c++) {
          shape[c] = surfReg(iv, MirrorDeposition::compShape + c);
        }

        // The world-frame shape operator must annihilate the normal. This is the invariant that a wrong lift or a
        // mismatched tangent basis breaks, and it breaks nothing else visible.
        worstSn = std::max(worstSn, MirrorDeposition::applyShapeOperator(shape, nc).vectorLength());

        Real twoHexact = 0.0;
        Real Kexact    = 0.0;
        analyticCurvature(xc, twoHexact, Kexact);

        const Real twoH = MirrorDeposition::meanCurvatureTimesTwo(shape);

        // det(S_c) is identically zero, because S_c annihilates the normal and so has rank at most SpaceDim-1.
        // That is a trap rather than a curiosity: anyone who "simplifies" gaussianCurvature() to a determinant gets
        // K = 0 everywhere, which silently degrades the exact Jacobian to the linearized one, whose error runs to
        // 48-77%. Assert it rather than describe it.
        //
        // Measured RELATIVE to the curvature scale. det(S_c) is a product of SpaceDim entries each of order 1/R, so
        // its roundoff floor grows as the surface tightens; an absolute bound passes on a gentle surface and fails
        // on a sharp one for no reason but arithmetic.
        Real detS = 0.0;
#if CH_SPACEDIM == 2
        detS = MirrorDeposition::shapeEntry(shape, 0, 0) * MirrorDeposition::shapeEntry(shape, 1, 1) -
               MirrorDeposition::shapeEntry(shape, 0, 1) * MirrorDeposition::shapeEntry(shape, 1, 0);
#else
        {
          const Real a00 = MirrorDeposition::shapeEntry(shape, 0, 0);
          const Real a01 = MirrorDeposition::shapeEntry(shape, 0, 1);
          const Real a02 = MirrorDeposition::shapeEntry(shape, 0, 2);
          const Real a11 = MirrorDeposition::shapeEntry(shape, 1, 1);
          const Real a12 = MirrorDeposition::shapeEntry(shape, 1, 2);
          const Real a22 = MirrorDeposition::shapeEntry(shape, 2, 2);

          detS = a00 * (a11 * a22 - a12 * a12) - a01 * (a01 * a22 - a12 * a02) + a02 * (a01 * a12 - a11 * a02);
        }
#endif

        const Real curvatureScale = std::pow(std::max(std::abs(twoH), 1.E-30), SpaceDim);

        worstNormK = std::max(worstNormK, std::abs(detS) / curvatureScale);

        // Anisotropy check: a torus must actually produce cells with distinct principal curvatures, or a coarsely
        // resolved one has silently degenerated into the sphere case -- which cannot see a frame error at all.
        // Meaningful in 3-D only: a 2-D surface is a curve and has exactly one principal curvature, so the second
        // root below is the trivial one along the normal and every cell would count as "anisotropic".
#if CH_SPACEDIM == 3
        const Real K    = MirrorDeposition::gaussianCurvature(shape);
        const Real disc = 0.25 * twoH * twoH - K;

        if (disc > 0.0) {
          const Real c1 = 0.5 * twoH + std::sqrt(disc);
          const Real c2 = 0.5 * twoH - std::sqrt(disc);

          if (std::abs(c1 - c2) > 0.1 * std::abs(c1 + c2)) {
            numAniso++;
          }
        }
#endif

        if (std::abs(twoHexact) > 0.0) {
          const Real relErr = std::abs(twoH - twoHexact) / std::abs(twoHexact);

          worstAll = std::max(worstAll, relErr);
          sumAll += relErr;
          numAll++;

          if (kappa <= 0.05) {
            worstSlim = std::max(worstSlim, relErr);
            numSlim++;
          }
        }
      }
    }

    if (hasCutCells[lvl] == 0) {
      pout() << "  level " << lvl << " reports no cut cells" << endl;
    }
  }

  numFitted = ParallelOps::sum(numFitted);
  numPlanar = ParallelOps::sum(numPlanar);
  numNone   = ParallelOps::sum(numNone);
  numAll    = ParallelOps::sum(numAll);
  numSlim   = ParallelOps::sum(numSlim);
  numAniso  = ParallelOps::sum(numAniso);

  pout() << "  cut cells: " << numFitted << " fitted, " << numPlanar << " planar fallback, " << numNone << " no data; "
         << numAniso << " anisotropic" << endl;

  if (numAll == 0) {
    pout() << "  FAIL no cut cells were scored at all" << endl;

    g_numFailures++;
  }
  else {
    expect("mean |d(2H)|/|2H| over all cut cells", ParallelOps::sum(sumAll) / numAll, g_meanTolerance);
    expect("worst |S_c n_c|", ParallelOps::max(worstSn), 1.E-12);
    expect("worst |det S_c|/|2H|^D (must be identically zero)", ParallelOps::max(worstNormK), 1.E-12);

    if (numSlim > 0) {
      expect("worst |d(2H)|/|2H| in the kappa <= 0.05 bin", ParallelOps::max(worstSlim), g_slimTolerance);
    }
    else {
      pout() << "  note: no cells in the kappa <= 0.05 bin" << endl;
    }
  }
}

int
main(int argc, char* argv[])
{
  ChomboDischarge::initialize(argc, argv);

  {
    ParmParse pp("MirrorSurfaceData");

    std::string  shape;
    Vector<Real> centre;

    pp.get("shape", shape);
    pp.get("radius", g_radius);
    pp.get("minor_radius", g_minorRadius);
    pp.get("fluid_inside", g_fluidInside);
    pp.get("mean_tolerance", g_meanTolerance);
    pp.get("slim_tolerance", g_slimTolerance);
    pp.getarr("centre", centre, 0, SpaceDim);

    for (int dir = 0; dir < SpaceDim; dir++) {
      g_centre[dir] = centre[dir];
    }

    if (shape == "sphere") {
      g_shape = Shape::Sphere;
    }
    else if (shape == "torus") {
#if CH_SPACEDIM == 2
      MayDay::Error("MirrorSurfaceData - 'shape = torus' is 3-D only; in 2-D TorusSdf is an annulus, which has a "
                    "single principal curvature and so cannot test anisotropy");
#endif
      g_shape = Shape::Torus;
    }
    else {
      MayDay::Error("MirrorSurfaceData - 'shape' must be 'sphere' or 'torus'");
    }
  }

  RefCountedPtr<BaseIF> surface;
  if (g_shape == Shape::Sphere) {
    surface = RefCountedPtr<BaseIF>(new SphereSdf(g_centre, g_radius, g_fluidInside));
  }
  else {
    surface = RefCountedPtr<BaseIF>(new TorusSdf(g_centre, g_radius, g_minorRadius, g_fluidInside));
  }

  auto compgeom = RefCountedPtr<ComputationalGeometry>(new ComputationalGeometry());
  compgeom->setElectrodes(Vector<Electrode>(1, Electrode(surface, true)));

  auto amr         = RefCountedPtr<AmrMesh>(new AmrMesh());
  auto tagger      = RefCountedPtr<CellTagger>(nullptr);
  auto timestepper = RefCountedPtr<MirrorSurfaceDataStepper>(new MirrorSurfaceDataStepper());
  auto engine      = RefCountedPtr<Driver>(new Driver(compgeom, timestepper, amr, tagger));

  engine->setupAndRun();

  const int numFailures = ParallelOps::sum(g_numFailures);

  if (numFailures > 0) {
    MayDay::Error("MirrorSurfaceData - one or more checks failed; see pout.*");
  }

  ChomboDischarge::finalize();
}
