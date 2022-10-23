/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoPlasmaPhysics.cpp
  @brief  Implementation of CD_ItoPlasmaPhysics.H
  @author Robert Marskar
*/

// Std includes
#include <fstream>
#include <sstream>

// Chombo includes
#include <PolyGeom.H>

// Our includes
#include <CD_ItoPlasmaPhysics.H>
#include <CD_Units.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::ItoPlasma;

ItoPlasmaPhysics::ItoPlasmaPhysics()
{

  // Default coupling
  m_coupling = coupling::LFA;

  // Initialize RNGs
  m_udist11 = std::uniform_real_distribution<Real>(-1., 1.);
  m_udist01 = std::uniform_real_distribution<Real>(0., 1.);

  m_reactions.clear();
  m_photoReactions.clear();

  // Default
  m_ppc = 32;

  // Default parameters for hybrid algorithm.
  m_Ncrit     = 10;
  m_eps       = 1.0;
  m_NSSA      = 100;
  m_SSAlim    = 0.1;
  m_algorithm = algorithm::hybrid;

  // Default parameter for lookup tables
  m_table_entries = 1000;
}

ItoPlasmaPhysics::~ItoPlasmaPhysics() {}

const Vector<RefCountedPtr<ItoSpecies>>&
ItoPlasmaPhysics::getItoSpecies() const
{
  return m_ItoSpecies;
}

const Vector<RefCountedPtr<RtSpecies>>&
ItoPlasmaPhysics::getRtSpecies() const
{
  return m_rtSpecies;
}

int
ItoPlasmaPhysics::getNumPlasmaSpecies() const
{
  return m_ItoSpecies.size();
}

int
ItoPlasmaPhysics::getNumPhotonSpecies() const
{
  return m_rtSpecies.size();
}

ItoPlasmaPhysics::coupling
ItoPlasmaPhysics::getCoupling() const
{
  return m_coupling;
}

Real
ItoPlasmaPhysics::initialSigma(const Real a_time, const RealVect a_pos) const
{
  return 0.0;
}

void
ItoPlasmaPhysics::addTable(const std::string a_table_name, const std::string a_file)
{

  LookupTable1D<2> table;

  this->readFile(table, a_file); // Read file
  table.sort();
  table.makeUniform(m_table_entries);    // Make table into a unifom table
  m_tables.emplace(a_table_name, table); // Add table
}

void
ItoPlasmaPhysics::readFile(LookupTable1D<2>& a_table, const std::string a_file)
{

  Real x, y;

  std::ifstream infile(a_file);
  std::string   line;

  while (std::getline(infile, line)) {

    // Trim string
    trim(line);

    std::istringstream iss(line);

    const bool skipline = (line.at(0) == '#') || (line.length() == 0);
    if (!skipline) {
      if (!(iss >> x >> y)) {
        continue;
      }
      a_table.addEntry(x, y);
    }
  }
  infile.close();
}

#include <CD_NamespaceFooter.H>
