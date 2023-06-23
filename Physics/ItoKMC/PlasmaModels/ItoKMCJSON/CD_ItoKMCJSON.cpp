/* chombo-discharge
 * Copyright © 2023 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoKMCJSON.cpp
  @brief  Implementation of CD_ItoKMCJson.cpp
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>
#include <ParmParse.H>

// Our includes
#include <CD_ItoKMCJSON.H>
#include <CD_Units.H>
#include <CD_DataParser.H>
#include <CD_ParticleManagement.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::ItoKMC;

ItoKMCJSON::ItoKMCJSON() noexcept
{
  CH_TIME("ItoKMCJSON::ItoKMCJSON");

  m_className = "ItoKMCJSON";
}

ItoKMCJSON::~ItoKMCJSON() noexcept { CH_TIME("ItoKMCJSON::~ItoKMCJSON"); }

void
ItoKMCJSON::parseRuntimeOptions() noexcept
{
  CH_TIME("ItoKMCJSON::parseRuntimeOptions");

  this->parsePPC();
  this->parseDebug();
  this->parseAlgorithm();
}

bool
ItoKMCJSON::doesFileExist(const std::string a_filename) const noexcept
{
  CH_TIME("ItoKMCJSON::doesFileExist");
  if (m_verbose) {
    pout() << m_className + "::doesFileExist" << endl;
  }

  std::ifstream istream(a_filename);

  return istream.good();
}

void
ItoKMCJSON::parseJSON() noexcept
{
  CH_TIME("ItoKMCJSON::parseJSON");
  if (m_verbose) {
    pout() << m_className + "::parseJSON" << endl;
  }

  if (!(this->doesFileExist(m_jsonFile))) {
    this->throwParserError(m_className + "::parseJSON -- file '" + m_jsonFile + "' does not exist");
  }

  // Parse the JSON file.
  std::ifstream f(m_jsonFile);
  m_json = json::parse(f, nullptr, true, true);
}
}

#include <CD_NamespaceFooter.H>
