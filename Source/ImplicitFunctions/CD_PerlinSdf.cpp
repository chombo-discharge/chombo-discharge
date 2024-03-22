/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_PerlinSdf.cpp
  @brief  Implementation of CD_PerlinSdf.H
  @author Robert Marskar
*/

// Our includes
#include <CD_PerlinSdf.H>
#include <CD_NamespaceHeader.H>

constexpr int PerlinSdf::m_permutationTable[256];

PerlinSdf::PerlinSdf(const Real     a_noiseAmp,
                     const RealVect a_noiseFreq,
                     const Real     a_persistence,
                     const int      a_octaves,
                     const bool     a_reseed)
{

  CH_assert(a_octaves >= 1);

  m_noiseAmp    = a_noiseAmp;
  m_noiseFreq   = a_noiseFreq;
  m_persistence = a_persistence;
  m_octaves     = a_octaves;

  // Use Ken Perlin's original permutation table
  for (int i = 0; i < 256; i++) {
    p[i]       = m_permutationTable[i];
    p[i + 256] = m_permutationTable[i];
  }

  // Reseed the permutation table
  if (a_reseed) {
    this->reseed();
  }
}

PerlinSdf::PerlinSdf(const PerlinSdf& a_inputIF)
{

  m_noiseAmp    = a_inputIF.m_noiseAmp;
  m_noiseFreq   = a_inputIF.m_noiseFreq;
  m_persistence = a_inputIF.m_persistence;
  m_octaves     = a_inputIF.m_octaves;

  for (int i = 0; i < 256; i++) {
    p[i]       = a_inputIF.p[i];
    p[i + 256] = a_inputIF.p[i + 256];
  }
}

PerlinSdf::~PerlinSdf()
{}

Real
PerlinSdf::value(const RealVect& a_pos) const
{
  return this->octaveNoise(a_pos);
}

BaseIF*
PerlinSdf::newImplicitFunction() const
{
  return static_cast<BaseIF*>(new PerlinSdf(*this));
}

int
PerlinSdf::random(const int i)
{
  return std::rand() % i;
}

void
PerlinSdf::reseed()
{

  for (int i = 0; i < 256; i++) {
    p[i] = i;
  }

  // Reseed the RNG
  int seed = time(nullptr);

  // Haven't had trouble with this, but I'm putting this here for safety.
#ifdef CH_MPI
  int result = MPI_Bcast(&seed, 1, MPI_INT, 0, Chombo_MPI::comm);
  if (result != MPI_SUCCESS) {
    MayDay::Error("PerlinSdf::reseed() - broadcast failed");
  }
#endif
  std::srand(seed);

  std::random_shuffle(&p[0], &p[256], random);

  for (int i = 0; i < 256; i++) {
    p[i + 256] = p[i];
  }
}

double
PerlinSdf::noise(const double a_x, const double a_y, const double a_z) const
{

  // Lower cube corner
  const int X = (int)std::floor(a_x) & 255;
  const int Y = (int)std::floor(a_y) & 255;
  const int Z = (int)std::floor(a_z) & 255;

  // Relative distance wrt lower cube corner
  const double x = a_x - std::floor(a_x);
  const double y = a_y - std::floor(a_y);
  const double z = a_z - std::floor(a_z);

  // Fade curves
  const double u = fade(x);
  const double v = fade(y);
  const double w = fade(z);

  // Hash coordinates of 8 cube corners
  const int A  = p[X] + Y;
  const int AA = p[A] + Z;
  const int AB = p[A + 1] + Z;
  const int B  = p[X + 1] + Y;
  const int BA = p[B] + Z;
  const int BB = p[B + 1] + Z;

  // Add blended results from 8 corners of cube
  return lerp(w,
              lerp(v,
                   lerp(u,
                        grad(p[AA], x, y, z),      // AND ADD
                        grad(p[BA], x - 1, y, z)), // BLENDED
                   lerp(u,
                        grad(p[AB], x, y - 1, z),       // RESULTS
                        grad(p[BB], x - 1, y - 1, z))), // FROM  8
              lerp(v,
                   lerp(u,
                        grad(p[AA + 1], x, y, z - 1),      // CORNERS
                        grad(p[BA + 1], x - 1, y, z - 1)), // OF CUBE
                   lerp(u, grad(p[AB + 1], x, y - 1, z - 1), grad(p[BB + 1], x - 1, y - 1, z - 1))));
}

Real
PerlinSdf::noise(const RealVect& a_pos) const
{

  Real x, y, z;
  x = a_pos[0];
  y = a_pos[1];
  z = 0.;

#if CH_SPACEDIM == 3
  z = a_pos[2];
#endif

  // Make noise [0, 1]
  return 0.5 * noise(x, y, z) + 0.5;
}

Real
PerlinSdf::octaveNoise(const RealVect& a_pos) const
{
  Real result = 0.0;

  Real     normFrac = 0.;
  RealVect freq     = m_noiseFreq;
  double   amp      = 1.; //m_noiseAmp;

  // Add noise octaves
  for (int i = 0; i < m_octaves; ++i) {
    result += noise(a_pos * freq) * amp;

    normFrac += amp;
    freq *= 1. / m_persistence;
    amp *= m_persistence;
  }

  // Normalize
  result *= m_noiseAmp; ///normFrac;

  return result;
}

Real
PerlinSdf::lerp(const Real t, const Real a, const Real b) const
{
  return a + t * (b - a);
}

Real
PerlinSdf::fade(const Real t) const
{
  return t * t * t * (t * (t * 6 - 15) + 10);
}

Real
PerlinSdf::grad(const int hash, const double x, const double y, const double z) const
{
  const int    h = hash & 15;
  const double u = h < 8 ? x : y;
  const double v = h < 4 ? y : h == 12 || h == 14 ? x : y;
  return ((h & 1) == 0 ? u : -u) + ((h & 2) == 0 ? v : -v);
}

#include <CD_NamespaceFooter.H>
