/* chombo-discharge
 * Copyright © 2024 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_RandomInterface.cpp
  @brief  Implementation of CD_RandomInterface.H
  @author Robert Marskar
*/

// Std includes
#include <random>
#include <chrono>

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_EBGeometryIF.H>
#include <CD_RandomInterface.H>
#include <CD_BoundedNoisePlane.H>
#include <CD_NamespaceHeader.H>

RandomInterface::ClampedNoisePlane::ClampedNoisePlane(const Vec3 a_point,
                                                      const Vec3 a_normal,
                                                      const Vec3 a_clampLo,
                                                      const Vec3 a_clampHi,
                                                      const Vec3 a_clampDx,
                                                      const Vec3 a_noiseFrequency,
                                                      const Real a_noiseAmplitude,
                                                      const Real a_noisePersistence,
                                                      const int  a_noiseOctaves) noexcept
{
  m_point   = a_point;
  m_normal  = a_normal;
  m_clampDx = a_clampDx;

  for (int dir = 0; dir < SpaceDim; dir++) {
    m_clampLo[dir] = std::min(a_clampLo[dir], a_clampHi[dir]);
    m_clampHi[dir] = std::max(a_clampLo[dir], a_clampHi[dir]);
  }

  m_plane  = std::make_shared<EBGeometry::PlaneSDF<Real>>(m_point, m_normal);
  m_perlin = std::make_shared<EBGeometry::PerlinSDF<Real>>(a_noiseAmplitude,
                                                           a_noiseFrequency,
                                                           a_noisePersistence,
                                                           a_noiseOctaves);
}

template <class URNG>
void
RandomInterface::ClampedNoisePlane::shuffle(URNG& a_rng) noexcept
{
  m_perlin->shuffle(a_rng);
}

Real
RandomInterface::ClampedNoisePlane::signedDistance(const Vec3& a_point) const noexcept
{
  const Vec3 x0 = m_point;
  const Vec3 x1 = a_point;
  const Vec3 xp = x1 - dot((x1 - x0), m_normal) * m_normal;

  // Clamping function.
  Real clamp = 1.0;
  for (int dir = 0; dir < SpaceDim; dir++) {
    const Real& x  = a_point[dir];
    const Real& xL = m_clampLo[dir];
    const Real& xH = m_clampHi[dir];
    const Real& dx = m_clampDx[dir];

    clamp *= 0.5 * (tanh((x - xL) / dx) - tanh((x - xH) / dx));
  }

  return m_plane->signedDistance(a_point) + clamp * m_perlin->signedDistance(xp);
}

RandomInterface::RandomInterface() noexcept
{
  ParmParse pp("RandomInterface");

  auto getVec = [&](const std::string id) -> Vec3 {
    Vector<Real> v;
    pp.getarr(id.c_str(), v, 0, SpaceDim);
#if CH_SPACEDIM == 2
    return Vec3(v[0], v[1], 0.0);
#else
    return Vec3(v[0], v[1], v[2]);
#endif
  };

  // Read in all input parameters
  Real smoothLength      = 0.0;
  Real gasPermittivity   = 1.0;
  Real solidPermittivity = 1.0;

  bool reseed    = false;
  bool usePlane1 = false;
  bool usePlane2 = false;

  Real noiseAmplitude1 = 0.0;
  Real noiseAmplitude2 = 0.0;

  Real noisePersistence1 = 0.0;
  Real noisePersistence2 = 0.0;

  int noiseOctaves1 = 0;
  int noiseOctaves2 = 0;

  Vec3 noiseFrequency1 = Vec3::zero();
  Vec3 noiseFrequency2 = Vec3::zero();

  Vec3 point1 = Vec3::zero();
  Vec3 point2 = Vec3::zero();

  Vec3 normal1 = Vec3::zero();
  Vec3 normal2 = Vec3::zero();

  Vec3 clampLo1 = Vec3::zero();
  Vec3 clampLo2 = Vec3::zero();

  Vec3 clampHi1 = Vec3::zero();
  Vec3 clampHi2 = Vec3::zero();

  Vec3 clampDx1 = Vec3::zero();
  Vec3 clampDx2 = Vec3::zero();

  pp.get("reseed", reseed);
  pp.get("smooth_len", smoothLength);
  pp.get("gas_permittivity", gasPermittivity);
  pp.get("solid_permittivity", solidPermittivity);

  // Get parameters for the planes
  pp.get("plane1.use", usePlane1);
  pp.get("plane2.use", usePlane2);
  pp.get("plane1.noise_amplitude", noiseAmplitude1);
  pp.get("plane2.noise_amplitude", noiseAmplitude2);
  pp.get("plane1.noise_persistence", noisePersistence1);
  pp.get("plane2.noise_persistence", noisePersistence2);
  pp.get("plane1.noise_octaves", noiseOctaves1);
  pp.get("plane2.noise_octaves", noiseOctaves2);

  noiseFrequency1 = getVec("plane1.noise_frequency");
  noiseFrequency2 = getVec("plane2.noise_frequency");
  point1          = getVec("plane1.point");
  point2          = getVec("plane2.point");
  normal1         = getVec("plane1.normal");
  normal2         = getVec("plane2.normal");
  clampLo1        = getVec("plane1.clamp_lo");
  clampLo2        = getVec("plane2.clamp_lo");
  clampHi1        = getVec("plane1.clamp_hi");
  clampHi2        = getVec("plane2.clamp_hi");
  clampDx1        = getVec("plane1.clamp_dx");
  clampDx2        = getVec("plane2.clamp_dx");

  std::shared_ptr<EBGeometry::ImplicitFunction<Real>> dielectric;
  std::shared_ptr<ClampedNoisePlane>                  plane1;
  std::shared_ptr<ClampedNoisePlane>                  plane2;

  // Planes are constructed along the y-axis
  if (usePlane1) {
    plane1 = std::make_shared<ClampedNoisePlane>(point1,
                                                 normal1,
                                                 clampLo1,
                                                 clampHi1,
                                                 clampDx1,
                                                 noiseFrequency1,
                                                 noiseAmplitude1,
                                                 noisePersistence1,
                                                 noiseOctaves1);

    // Find a seed to use for the RNG
    int seed = reseed ? std::chrono::system_clock::now().time_since_epoch().count() : 0;
#ifdef CH_MPI
    MPI_Bcast(&seed, 1, MPI_INT, 0, Chombo_MPI::comm);
#endif

    std::mt19937 rng(seed);

    plane1->shuffle(rng);
  }

  if (usePlane2) {
    plane2 = std::make_shared<ClampedNoisePlane>(point2,
                                                 normal2,
                                                 clampLo2,
                                                 clampHi2,
                                                 clampDx2,
                                                 noiseFrequency2,
                                                 noiseAmplitude2,
                                                 noisePersistence2,
                                                 noiseOctaves2);

    // Find a seed to use for the RNG
    int seed = reseed ? std::chrono::system_clock::now().time_since_epoch().count() : 1;
#ifdef CH_MPI
    MPI_Bcast(&seed, 1, MPI_INT, 0, Chombo_MPI::comm);
#endif

    std::mt19937 rng(seed);
    plane2->shuffle(rng);
  }

  // Merge the planes if using two planes.
  if (usePlane1 && usePlane2) {
    if (smoothLength > 0.0) {
      dielectric = EBGeometry::SmoothUnion<Real>(plane1, plane2, smoothLength);
    }
    else {
      dielectric = EBGeometry::Union<Real>(plane1, plane2);
    }
  }
  else if (usePlane1) {
    dielectric = plane1;
  }
  else if (usePlane2) {
    dielectric = plane2;
  }

  if (usePlane1 || usePlane2) {
    auto implicitFunction = RefCountedPtr<BaseIF>(new EBGeometryIF<>(dielectric, true));
    m_dielectrics.push_back(Dielectric(implicitFunction, solidPermittivity));
  }

  this->setGasPermittivity(gasPermittivity);
}

RandomInterface::~RandomInterface() noexcept{

}

#include <CD_NamespaceFooter.H>
