/* chombo-discharge
 * Copyright © 2023 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoKMCJSON.cpp
  @brief  Implementation of CD_ItoKMCJSON.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>
#include <ParmParse.H>

// Our includes
#include <CD_ItoKMCJSON.H>
#include <CD_ItoKMCCDRSpecies.H>
#include <CD_Units.H>
#include <CD_DataParser.H>
#include <CD_ParticleManagement.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::ItoKMC;

ItoKMCJSON::ItoKMCJSON() noexcept
{
  CH_TIME("ItoKMCJSON::ItoKMCJSON");

  m_verbose   = false;
  m_className = "ItoKMCJSON";

  this->parseVerbose();
  this->parsePPC();
  this->parseDebug();
  this->parseAlgorithm();
  this->parseJSON();

  // Initialize the gas law and background species
  this->initializeGasLaw();
  this->initializeBackgroundSpecies();

  // Initialize the plasma species
  this->initializePlasmaSpecies();
  this->parseInitialData();

  // Define internals
  this->define();

  // Useful shortcut.
  m_numPlasmaSpecies = this->getNumPlasmaSpecies();
}

ItoKMCJSON::~ItoKMCJSON() noexcept { CH_TIME("ItoKMCJSON::~ItoKMCJSON"); }

void
ItoKMCJSON::parseRuntimeOptions() noexcept
{
  CH_TIME("ItoKMCJSON::parseRuntimeOptions");

  this->parseVerbose();
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

std::string
ItoKMCJSON::trim(const std::string& a_string) const noexcept
{
  CH_TIME("ItoKMCJSON::trim");
  if (m_verbose) {
    pout() << m_className + "::trim" << endl;
  }

  auto ltrim = [](const std::string a_s) -> std::string {
    std::string s = a_s;
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
  };

  auto rtrim = [](std::string a_s) -> std::string {
    std::string s = a_s;
    s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    return s;
  };

  return ltrim(rtrim(a_string));
}

void
ItoKMCJSON::throwParserError(const std::string a_error) const noexcept
{
  CH_TIME("ItoKMCJSON::throwParserError");
  if (m_verbose) {
    pout() << m_className + "::trim" << endl;
  }

  pout() << a_error << endl;

  MayDay::Error(a_error.c_str());
}

void
ItoKMCJSON::throwParserWarning(const std::string a_warning) const noexcept
{
  CH_TIME("ItoKMCJSON::throwParserError");
  if (m_verbose) {
    pout() << m_className + "::trim" << endl;
  }

  pout() << a_warning << endl;

  MayDay::Warning(a_warning.c_str());
}

bool
ItoKMCJSON::containsWildcard(const std::string a_str) const noexcept
{
  CH_TIME("ItoKMCJSON::containsWildcard");
  if (m_verbose) {
    pout() << m_className + "::containsWildcard" << endl;
  }

  return (a_str.find("@") != std::string::npos);
}

void
ItoKMCJSON::checkMolarFraction(const RealVect a_position) const noexcept
{
  CH_TIME("ItoKMCJSON::checkMolarFraction");
  if (m_verbose) {
    pout() << m_className + "::checkMolarFraction" << endl;
  }

  Real sumFractions = 0.0;
  for (const auto& species : m_backgroundSpecies) {
    sumFractions += species.molarFraction(a_position);
  }

  if (std::abs(sumFractions - 1.0) > std::numeric_limits<Real>::epsilon()) {
    MayDay::Warning("ItoKMCJSON::checkMolarFraction -- fractions do not sum to 1");
  }
}

void
ItoKMCJSON::parseVerbose() noexcept
{
  CH_TIME("ItoKMCJSON::parseVerbose");
  if (m_verbose) {
    pout() << m_className + "::parseVerbose" << endl;
  }

  ParmParse pp(m_className.c_str());

  pp.get("verbose", m_verbose);
}

void
ItoKMCJSON::parseJSON()
{
  CH_TIME("ItoKMCJSON::parseJSON");
  if (m_verbose) {
    pout() << m_className + "::parseJSON" << endl;
  }

  const std::string baseError = m_className + "::parseJSON";

  ParmParse pp(m_className.c_str());

  pp.get("chemistry_file", m_jsonFile);

  if (!(this->doesFileExist(m_jsonFile))) {
    const std::string parseError = baseError + " -- file '" + m_jsonFile + "' does not exist";

    this->throwParserError(parseError.c_str());
  }
  else {
    // Parse the JSON file.
    std::ifstream f(m_jsonFile);

    constexpr auto callback        = nullptr;
    constexpr auto allowExceptions = true;
    constexpr auto ignoreComments  = true;

    m_json = nlohmann::json::parse(f, callback, allowExceptions, ignoreComments);

    if (m_verbose) {
      pout() << m_className + "::parseJSON - done reading file!" << endl;
    }
  }
}

void
ItoKMCJSON::initializeGasLaw() noexcept
{
  CH_TIME("ItoKMCJSON::initializeGasLaw");
  if (m_verbose) {
    pout() << m_className + "::initializeGasLaw" << endl;
  }

  const std::string baseError = m_className + "::initializeGasLaw";

  // TLDR: This routine exists for populating m_gasPressure, m_gasTemperature, and m_gasDensity. If you want to add a new gas law, this
  //       routine is where you would do it.

  if (!(m_json["gas"]["law"].contains("id"))) {
    this->throwParserError(baseError + " but field 'gas/law/id' is missing");
  }

  const std::string curLaw = this->trim(m_json["gas"]["law"]["id"].get<std::string>());

  if ((!m_json["gas"]["law"][curLaw].contains("type"))) {
    this->throwParserError(baseError + " but specified gas law '" + curLaw + "' is missing the 'type' specifier");
  }

  const std::string type = this->trim(m_json["gas"]["law"][curLaw]["type"].get<std::string>());

  if (type == "ideal") {
    if (!(m_json["gas"]["law"][curLaw].contains("temperature"))) {
      this->throwParserError(baseError + " and got ideal gas law but field 'temperature' is missing");
    }
    if (!(m_json["gas"]["law"][curLaw].contains("pressure"))) {
      this->throwParserError(baseError + " and got ideal gas law but field 'pressure' is missing");
    }

    const Real T0   = m_json["gas"]["law"][curLaw]["temperature"].get<Real>();
    const Real P0   = m_json["gas"]["law"][curLaw]["pressure"].get<Real>();
    const Real P    = P0 * Units::atm2pascal;
    const Real Rho0 = (P * Units::Na) / (T0 * Units::R);

    m_gasTemperature = [T0](const RealVect a_position) -> Real {
      return T0;
    };
    m_gasPressure = [P](const RealVect a_position) -> Real {
      return P;
    };
    m_gasDensity = [Rho0](const RealVect a_position) -> Real {
      return Rho0;
    };
  }
  else {
    const std::string parseError = baseError + " but gas law '" + type + "' is not supported";

    this->throwParserError(parseError);
  }
}

void
ItoKMCJSON::initializeBackgroundSpecies() noexcept
{
  CH_TIME("ItoKMCJSON::initializeBackgroundSpecies");
  if (m_verbose) {
    pout() << m_className + "::initializeBackgroundSpecies" << endl;
  }

  const std::string baseError = m_className + "::initializeBackgroundSpecies";

  if (!(m_json["gas"].contains("background species"))) {
    this->throwParserError(baseError + " but field 'gas/background species' is missing");
  }
  else {
    const auto backgroundSpecies = m_json["gas"]["background species"];

    for (const auto& species : backgroundSpecies) {
      bool plotSpecies = false;

      if (!(species.contains("id"))) {
        this->throwParserError(baseError + " but species does not contain the 'id' field");
      }
      if (species.contains("plot")) {
        plotSpecies = species["plot"].get<bool>();
      }
      if (!(species["molar fraction"].contains("type"))) {
        this->throwParserError(baseError + " -- but species does not contain the 'molar fraction/type' field");
      }

      const std::string speciesName = this->trim(species["id"].get<std::string>());
      const std::string molFracType = this->trim(species["molar fraction"]["type"].get<std::string>());

      std::function<Real(const RealVect x)> molarFraction;

      if (this->containsWildcard(speciesName)) {
        this->throwParserError(baseError + " but species '" + speciesName + "' should not contain wildcard @");
      }

      if (molFracType == "constant") {
        if (!(species["molar fraction"].contains("value"))) {
          this->throwParserError(baseError + " and got constant molar fraction but field 'value' is missing");
        }

        const Real m = species["molar fraction"]["value"].get<Real>();

        molarFraction = [m](const RealVect x) -> Real {
          return m;
        };
      }
      else if (molFracType == "table vs height") {
        const std::string baseErrorID = baseError + "and got 'table vs height for species '" + speciesName + "'";

        int heightColumn   = 0;
        int fractionColumn = 1;
        int axis           = -1;
        int numPoints      = 1000;

        TableSpacing spacing = TableSpacing::Uniform;

        // Required arguments.
        if (!(species["molar fraction"].contains("file"))) {
          this->throwParserError(baseErrorID + " but 'file' is not specified");
        }
        if (!(species["molar fraction"].contains("axis"))) {
          this->throwParserError(baseErrorID + " but 'axis' is not specified");
        }

        const std::string fileName  = this->trim(species["molar fraction"]["file"].get<std::string>());
        const std::string whichAxis = this->trim(species["molar fraction"]["axis"].get<std::string>());

        if (whichAxis == "x") {
          axis = 0;
        }
        else if (whichAxis == "y") {
          axis = 1;
        }
        else if (whichAxis == "z") {
          axis = 2;
        }
        else {
          this->throwParserError(baseErrorID + " but 'axis' is not 'x', 'y', or 'z'");
        }

        if (!(this->doesFileExist(fileName))) {
          this->throwParserError(baseErrorID + "but file '" + fileName + "' does not exist");
        }

        // Optional arguments
        if (species["molar fraction"].contains("height column")) {
          heightColumn = species["molar fraction"]["height column"].get<int>();

          if (heightColumn < 0) {
            this->throwParserError(baseErrorID + " but can't have 'height column' < 0");
          }
        }
        if (species["molar fraction"].contains("molar fraction column")) {
          fractionColumn = species["molar fraction"]["molar fraction column"].get<int>();
          if (fractionColumn < 0) {
            this->throwParserError(baseErrorID + " but can't have 'molar fraction column' < 0");
          }
        }
        if (species["molar fraction"].contains("num points")) {
          numPoints = species["molar fraction"]["num points"].get<int>();

          if (numPoints < 2) {
            this->throwParserError(baseErrorID + " but can't have 'num points' < 2");
          }
        }

        LookupTable1D<2> table = DataParser::simpleFileReadASCII(fileName, heightColumn, fractionColumn);

        if (species["molar fraction"].contains("height scale")) {
          const Real scaling = species["molar fraction"]["height scale"].get<Real>();
          if (scaling <= 0.0) {
            this->throwParserError(baseErrorID + " but can't have 'height scale' <= 0.0");
          }

          table.scale<0>(scaling);
        }
        if (species["molar fraction"].contains("fraction scale")) {
          const Real scaling = species["molar fraction"]["fraction scale"].get<Real>();
          if (scaling <= 0.0) {
            this->throwParserError(baseErrorID + " but can't have 'fraction scale' <= 0.0");
          }

          table.scale<1>(scaling);
        }
        if (species["molar fraction"].contains("min height")) {
          const Real minHeight = species["molar fraction"]["min height"].get<Real>();
          table.setMinRange(minHeight, 0);
        }
        if (species["molar fraction"].contains("max height")) {
          const Real maxHeight = species["molar fraction"]["max height"].get<Real>();
          table.setMaxRange(maxHeight, 0);
        }
        if (species["molar fraction"].contains("spacing")) {
          const std::string whichSpacing = species["molar fraction"]["spacing"].get<std::string>();

          if (whichSpacing == "linear") {
            spacing = TableSpacing::Uniform;
          }
          else if (whichSpacing == "exponential") {
            spacing = TableSpacing::Exponential;
          }
          else {
            this->throwParserError(baseErrorID + " but spacing '" + whichSpacing + "' is not supported");
          }
        }

        table.setTableSpacing(spacing);
        table.sort(0);
        table.makeUniform(numPoints);

        molarFraction = [table, axis](const RealVect a_position) -> Real {
          return table.getEntry<1>(a_position[axis]);
        };

        if (species["molar fraction"].contains("dump")) {
          const std::string dumpId = species["molar fraction"]["dump"].get<std::string>();

          table.dumpTable(dumpId);
        }
      }
      else {
        const std::string err = baseError + " but got unsupported molar fraction type '" + molFracType + "'";

        MayDay::Error(err.c_str());
      }

      const int idx = m_backgroundSpecies.size();

      m_backgroundSpecies.push_back(ItoKMCBackgroundSpecies(speciesName, molarFraction));
      m_backgroundSpeciesPlot.push_back(plotSpecies);
      m_backgroundSpeciesMap[speciesName] = idx;
      m_backgroundSpeciesMapInverse[idx]  = speciesName;
    }
  }

  // Do a dummy test to see if molar fractions sum to one like they should. This may break if the user inputs function-based values.
  this->checkMolarFraction(RealVect::Zero);
}

void
ItoKMCJSON::initializePlasmaSpecies() noexcept
{
  CH_TIME("ItoKMCJSON::initializePlasmaSpecies");
  if (m_verbose) {
    pout() << m_className + "::initializePlasmaSpecies" << endl;
  }

  const std::string baseError = "ItoKMCJSON::initializePlasmaSpecies";

  if (!(m_json.contains("plasma species"))) {
    this->throwParserError(baseError + " but did not find field 'plasma species'");
  }

  for (const auto& species : m_json["plasma species"]) {
    if (!(species.contains("id"))) {
      this->throwParserError(baseError + " but did not find field 'id' for one of the species");
    }

    const std::string speciesID   = species["id"].get<std::string>();
    const std::string baseErrorID = baseError + " for species '" + speciesID + "'";

    if (this->containsWildcard(speciesID)) {
      this->throwParserError(baseErrorID + " but species name '" + speciesID + "' should not contain wildcard @");
    }
    if (!(species.contains("Z"))) {
      this->throwParserError(baseErrorID + " but did not find 'Z' (must be integer)");
    }
    if (!(species.contains("solver"))) {
      this->throwParserError(baseErrorID + " but did not field find 'solver' (must be 'ito' or 'cdr')");
    }
    if (!(species.contains("mobile"))) {
      this->throwParserError(baseErrorID + " but did not find field 'mobile' (must be true/false)");
    }
    if (!(species.contains("diffusive"))) {
      this->throwParserError(baseErrorID + " but did not find field 'diffusive' (must be true/false)");
    }

    const int         Z         = species["Z"].get<int>();
    const std::string solver    = species["solver"].get<std::string>();
    const bool        mobile    = species["mobile"].get<bool>();
    const bool        diffusive = species["diffusive"].get<bool>();

    if (solver == "ito") {
      m_itoSpecies.push_back(RefCountedPtr<ItoSpecies>(new ItoSpecies(speciesID, Z, mobile, diffusive)));
    }
    else if (solver == "cdr") {
      m_cdrSpecies.push_back(RefCountedPtr<CdrSpecies>(new ItoKMCCDRSpecies(speciesID, Z, mobile, diffusive)));
    }
    else {
      this->throwParserError(baseErrorID + " but 'solver' field must be either 'cdr' or 'ito'");
    }

    if (m_verbose) {
      pout() << "ItoKMCJSON::initializePlasmaSpecies, instantiating species:"
             << "\n"
             << "\tName             = " << speciesID << "\n"
             << "\tZ                = " << Z << "\n"
             << "\tMobile           = " << mobile << "\n"
             << "\tDiffusive        = " << diffusive << "\n"
             << "\tSolver type      = " << solver << "\n"
             << "\n";
    }
  }
}

void
ItoKMCJSON::parseInitialData() noexcept
{
  CH_TIME("ItoKMCJSON::parseInitialData");
  if (m_verbose) {
    pout() << m_className + "::parseInitialData" << endl;
  }
}

Vector<Real>
ItoKMCJSON::computeMobilities(const Real a_time, const RealVect a_pos, const RealVect a_E) const noexcept
{
  CH_TIME("ItoKMCJSON::computeMobilities");
  if (m_verbose) {
    pout() << m_className + "::computeMobilities" << endl;
  }

  Vector<Real> mobilityCoefficients(m_numPlasmaSpecies, 0.0);

  return mobilityCoefficients;
}

Vector<Real>
ItoKMCJSON::computeDiffusionCoefficients(const Real a_time, const RealVect a_pos, const RealVect a_E) const noexcept
{
  CH_TIME("ItoKMCJSON::computeDiffusionCoefficients");
  if (m_verbose) {
    pout() << m_className + "::computeDiffusionCoefficients" << endl;
  }

  Vector<Real> diffusionCoefficients(m_numPlasmaSpecies, 0.0);

  return diffusionCoefficients;
}

#include <CD_NamespaceFooter.H>
