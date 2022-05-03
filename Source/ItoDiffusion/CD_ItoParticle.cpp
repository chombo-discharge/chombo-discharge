/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoParticle.cpp
  @brief  Implementation of CD_ItoParticle.H
  @author Robert Marskar
*/

// Our includes
#include <CD_ItoParticle.H>
#include <CD_NamespaceHeader.H>

int ItoParticle::s_numRuntimeScalars = 0;
int ItoParticle::s_numRuntimeVectors = 0;

void
ItoParticle::setNumRuntimeScalars(const int a_numRuntimeScalars)
{
  s_numRuntimeScalars = a_numRuntimeScalars;
}

void
ItoParticle::setNumRuntimeVectors(const int a_numRuntimeVectors)
{
  s_numRuntimeVectors = a_numRuntimeVectors;
}

ItoParticle::ItoParticle() : BinItem() { this->allocateRuntimeBuffers(); }

ItoParticle::ItoParticle(const ItoParticle& a_other)
{
  m_position      = a_other.m_position;
  m_mass          = a_other.m_mass;
  m_velocity      = a_other.m_velocity;
  m_diffusion     = a_other.m_diffusion;
  m_oldPosition   = a_other.m_oldPosition;
  m_mobility      = a_other.m_mobility;
  m_averageEnergy = a_other.m_averageEnergy;

  this->allocateRuntimeBuffers();
}

ItoParticle::ItoParticle(const Real      a_mass,
                         const RealVect& a_position,
                         const RealVect& a_velocity,
                         const Real      a_diffusion,
                         const Real      a_mobility,
                         const Real      a_averageEnergy)
  : BinItem(a_position)
{
  m_mass          = a_mass;
  m_velocity      = a_velocity;
  m_diffusion     = a_diffusion;
  m_oldPosition   = a_position;
  m_mobility      = a_mobility;
  m_averageEnergy = a_averageEnergy;

  this->allocateRuntimeBuffers();
}

ItoParticle::~ItoParticle()
{
  if (m_runtimeScalars != nullptr) {
    delete[] m_runtimeScalars;
    m_runtimeScalars = nullptr;
  }

  if (m_runtimeVectors != nullptr) {
    delete[] m_runtimeVectors;
    m_runtimeVectors = nullptr;
  }
}

void
ItoParticle::allocateRuntimeBuffers()
{
  m_runtimeScalars = nullptr;
  m_runtimeVectors = nullptr;

  if (s_numRuntimeScalars >= 0) {
    m_runtimeScalars = new Real[s_numRuntimeScalars];
  }

  if (s_numRuntimeVectors >= 0) {
    m_runtimeVectors = new RealVect[s_numRuntimeVectors];
  }
}

bool
ItoParticle::operator==(const ItoParticle& a_p) const
{
  return (m_mass == a_p.m_mass && m_position == a_p.m_position && m_velocity == a_p.m_velocity &&
          m_diffusion == a_p.m_diffusion);
}

bool
ItoParticle::operator==(const ItoParticle* a_p) const
{
  return (*this == *a_p);
}

bool
ItoParticle::operator!=(const ItoParticle& a_p) const
{
  return !(*this == a_p);
}

int
ItoParticle::size() const
{
  const int compileSize = BinItem::size() + sizeof(m_mass) + sizeof(m_velocity) + sizeof(m_diffusion) +
                          sizeof(m_oldPosition) + sizeof(m_mobility) + sizeof(m_averageEnergy);
  //  const int runtimeSize = sizeof(*m_runtimeScalars) + sizeof(*m_runtimeVectors);
  const int runtimeSize = s_numRuntimeScalars * sizeof(Real) + s_numRuntimeVectors * sizeof(RealVect);

  //  std::cout << s_numRuntimeVectors*sizeof(RealVect) << std::endl;

  return compileSize + runtimeSize;
}

void
ItoParticle::linearOut(void* buf) const
{
  Real* buffer = (Real*)buf;
  D_TERM6(*buffer++ = m_position[0];, *buffer++ = m_position[1];, *buffer++ = m_position[2];, *buffer++ = m_position[3];
          , *buffer++ = m_position[4];
          , *buffer++ = m_position[5];);

  D_TERM6(*buffer++ = m_velocity[0];, *buffer++ = m_velocity[1];, *buffer++ = m_velocity[2];, *buffer++ = m_velocity[3];
          , *buffer++ = m_velocity[4];
          , *buffer++ = m_velocity[5];);

  D_TERM6(*buffer++ = m_oldPosition[0];, *buffer++ = m_oldPosition[1];, *buffer++ = m_oldPosition[2];
          , *buffer++ = m_oldPosition[3];
          , *buffer++ = m_oldPosition[4];
          , *buffer++ = m_oldPosition[5];);

  *buffer++ = m_diffusion;
  *buffer++ = m_mass;
  *buffer++ = m_mobility;

  // Go through run-time memory
  for (int i = 0; i < s_numRuntimeScalars; i++) {
    *buffer++ = m_runtimeScalars[i];
  }

  for (int i = 0; i < s_numRuntimeVectors; i++) {
    D_TERM6(*buffer++ = m_runtimeVectors[i][0];, *buffer++ = m_runtimeVectors[i][1];
            , *buffer++                                    = m_runtimeVectors[i][2];
            , *buffer++                                    = m_runtimeVectors[i][3];
            , *buffer++                                    = m_runtimeVectors[i][4];
            , *buffer++                                    = m_runtimeVectors[i][5];);
  }

  *buffer = m_averageEnergy;
}

void
ItoParticle::linearIn(void* buf)
{
  Real* buffer = (Real*)buf;
  D_TERM6(m_position[0] = *buffer++;, m_position[1] = *buffer++;, m_position[2] = *buffer++;, m_position[3] = *buffer++;
          , m_position[4] = *buffer++;
          , m_position[5] = *buffer++;);

  D_TERM6(m_velocity[0] = *buffer++;, m_velocity[1] = *buffer++;, m_velocity[2] = *buffer++;, m_velocity[3] = *buffer++;
          , m_velocity[4] = *buffer++;
          , m_velocity[5] = *buffer++;);

  D_TERM6(m_oldPosition[0] = *buffer++;, m_oldPosition[1] = *buffer++;, m_oldPosition[2] = *buffer++;
          , m_oldPosition[3] = *buffer++;
          , m_oldPosition[4] = *buffer++;
          , m_oldPosition[5] = *buffer++;);

  m_diffusion = *buffer++;
  m_mass      = *buffer++;
  m_mobility  = *buffer++;

  // Go through run-time memory
  for (int i = 0; i < s_numRuntimeScalars; i++) {
    m_runtimeScalars[i] = *buffer++;
  }
  for (int i = 0; i < s_numRuntimeVectors; i++) {
    D_TERM6(m_runtimeVectors[i][0] = *buffer++;, m_runtimeVectors[i][1] = *buffer++;
            , m_runtimeVectors[i][2]                                    = *buffer++;
            , m_runtimeVectors[i][3]                                    = *buffer++;
            , m_runtimeVectors[i][4]                                    = *buffer++;
            , m_runtimeVectors[i][5]                                    = *buffer++;);
  }

  m_averageEnergy = *buffer;
}

std::ostream&
operator<<(std::ostream& ostr, const ItoParticle& p)
{
  ostr << " ItoParticle : " << std::endl;
  ostr << " mass " << p.mass() << std::endl;
  ostr << " position ( ";
  for (int i = 0; i < SpaceDim; ++i) {
    ostr << " " << p.position(i);
  }
  ostr << " ) ";
  ostr << std::endl << " velocity ( ";
  for (int i = 0; i < SpaceDim; ++i) {
    ostr << " " << p.velocity(i);
  }
  ostr << " ) ";

  ostr << std::endl << " diffusion " << p.diffusion() << std::endl;
  ostr << std::endl << " mobility " << p.mobility() << std::endl;
  ostr << std::endl << " energy " << p.energy() << std::endl;

  return ostr;
}

#include <CD_NamespaceFooter.H>
