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
#include <CD_ItoKMCPhotonSpecies.H>
#include <CD_Units.H>
#include <CD_DataParser.H>
#include <CD_ParticleManagement.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::ItoKMC;

ItoKMCJSON::ItoKMCJSON()
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

  // Initialize Townsend coefficients
  this->initializeTownsendCoefficient("alpha");
  this->initializeTownsendCoefficient("eta");

  // Initialize the plasma species
  this->initializePlasmaSpecies();
  this->initializeParticles();
  this->initializeMobilities();
  this->initializeDiffusionCoefficients();
  this->initializeTemperatures();

  // Initialize the photon species
  this->initializePhotonSpecies();

  // Initialize reactions
  this->initializePlasmaReactions();
  this->initializePhotoReactions();
  this->initializeSurfaceReactions();

  // Define internals. This includes the KMC solver instantiation.
  this->define();
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
ItoKMCJSON::doesFileExist(const std::string a_filename) const noexcept
{
  CH_TIME("ItoKMCJSON::doesFileExist");
  if (m_verbose) {
    pout() << m_className + "::doesFileExist" << endl;
  }

  std::ifstream istream(a_filename);

  return istream.good();
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

bool
ItoKMCJSON::isBackgroundSpecies(const std::string& a_name) const noexcept
{
  CH_TIME("ItoKMCJSON::isBackgroundSpecies");
  if (m_verbose) {
    pout() << m_className + "::isBackgroundSpecies" << endl;
  }

  return m_backgroundSpeciesMap.count(a_name) != 0;
}

bool
ItoKMCJSON::isPlasmaSpecies(const std::string& a_name) const noexcept
{
  CH_TIME("ItoKMCJSON::isPlasmaSpecies");
  if (m_verbose) {
    pout() << m_className + "::isPlasmaSpecies" << endl;
  }

  return m_plasmaSpeciesTypes.count(a_name) != 0;
}

bool
ItoKMCJSON::isPhotonSpecies(const std::string& a_name) const noexcept
{
  CH_TIME("ItoKMCJSON::isPhotonSpecies");
  if (m_verbose) {
    pout() << m_className + "::isPhotonSpecies" << endl;
  }

  return m_photonIndexMap.count(a_name) != 0;
}

bool
ItoKMCJSON::containsBracket(const std::string a_str) const noexcept
{
  CH_TIME("ItoKMCJSON::containsBracket");
  if (m_verbose) {
    pout() << m_className + "::containsBracket" << endl;
  }

  const std::list<char> bracketList{'(', ')', '[', ']', '{', '}'};

  bool containsBracket = false;

  for (const auto& b : bracketList) {
    if (a_str.find(b) != std::string::npos) {
      containsBracket = true;

      break;
    }
  }

  return containsBracket;
}

bool
ItoKMCJSON::isBracketed(const std::string a_str) const noexcept
{
  return (a_str.front() == '(' && a_str.back() == ')');
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
ItoKMCJSON::initializeGasLaw()
{
  CH_TIME("ItoKMCJSON::initializeGasLaw");
  if (m_verbose) {
    pout() << m_className + "::initializeGasLaw" << endl;
  }

  const std::string baseError = m_className + "::initializeGasLaw";

  // TLDR: This routine exists for populating m_gasPressure, m_gasTemperature, and m_gasNumberDensity. If you want to add a new gas law, this
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
    m_gasNumberDensity = [Rho0](const RealVect a_position) -> Real {
      return Rho0;
    };
  }
  else {
    const std::string parseError = baseError + " but gas law '" + type + "' is not supported";

    this->throwParserError(parseError);
  }
}

void
ItoKMCJSON::initializeBackgroundSpecies()
{
  CH_TIME("ItoKMCJSON::initializeBackgroundSpecies");
  if (m_verbose) {
    pout() << m_className + "::initializeBackgroundSpecies" << endl;
  }

  const std::string baseError = m_className + "::initializeBackgroundSpecies";

  if (!(m_json["gas"].contains("background species"))) {
    this->throwParserWarning(baseError + " but field 'gas/background species' is missing");
  }
  else {
    const auto backgroundSpecies = m_json["gas"]["background species"];

    for (const auto& species : backgroundSpecies) {
      bool plotSpecies = false;

      if (!(species.contains("id"))) {
        this->throwParserError(baseError + "but 'id' field is not specified");
      }
      const std::string speciesName = species["id"].get<std::string>();

      if (this->containsWildcard(speciesName)) {
        this->throwParserError(baseError + " but species '" + speciesName + "' should not contain wildcard @");
      }
      if (this->containsBracket(speciesName)) {
        this->throwParserError(baseError + " but species '" + speciesName + "' should not contain brackets");
      }
      if (m_allSpecies.count(speciesName) != 0) {
        this->throwParserError(baseError + " but species '" + speciesName + "' was already defined elsewhere)");
      }

      if (species.contains("plot")) {
        plotSpecies = species["plot"].get<bool>();
      }
      if (!(species["molar fraction"].contains("type"))) {
        this->throwParserError(baseError + " -- but species does not contain the 'molar fraction/type' field");
      }

      const std::string molFracType = this->trim(species["molar fraction"]["type"].get<std::string>());

      std::function<Real(const RealVect x)> molarFraction;

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
          const std::string whichSpacing = this->trim(species["molar fraction"]["spacing"].get<std::string>());

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
          const std::string dumpId = this->trim(species["molar fraction"]["dump"].get<std::string>());

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
      m_allSpecies.emplace(speciesName);
    }
  }

  // Do a dummy test to see if molar fractions sum to one like they should. This may break if the user inputs function-based values.
  this->checkMolarFraction(RealVect::Zero);
}

void
ItoKMCJSON::initializePlasmaSpecies()
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
      this->throwParserError(baseError + " but one of the species does not contain the 'id' field");
    }
    const std::string speciesID   = species["id"].get<std::string>();
    const std::string baseErrorID = baseError + " for species '" + speciesID + "'";

    if (this->containsWildcard(speciesID)) {
      this->throwParserError(baseErrorID + " but species name '" + speciesID + "' should not contain wildcard @");
    }
    if (this->containsBracket(speciesID)) {
      this->throwParserError(baseError + " but species '" + speciesID + "' should not contain brackets");
    }
    if (m_allSpecies.count(speciesID) != 0) {
      this->throwParserError(baseErrorID + " but species '" + speciesID + "' was already defined elsewhere");
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
      m_plasmaSpeciesTypes[speciesID] = SpeciesType::Ito;
      m_itoSpeciesMap[speciesID]      = m_itoSpecies.size();

      m_itoSpecies.push_back(RefCountedPtr<ItoSpecies>(new ItoSpecies(speciesID, Z, mobile, diffusive)));
    }
    else if (solver == "cdr") {
      m_plasmaSpeciesTypes[speciesID] = SpeciesType::CDR;
      m_cdrSpeciesMap[speciesID]      = m_cdrSpecies.size();

      m_cdrSpecies.push_back(RefCountedPtr<CdrSpecies>(new ItoKMCCDRSpecies(speciesID, Z, mobile, diffusive)));
    }
    else {
      this->throwParserError(baseErrorID + " but 'solver' field must be either 'cdr' or 'ito'");
    }

    m_allSpecies.emplace(speciesID);

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

  m_numPlasmaSpecies = this->getNumPlasmaSpecies();

  // Initialize the solver map, which will assist us when we index into solvers later on.
  int species = 0;
  for (int i = 0; i < m_itoSpecies.size(); i++, species++) {
    m_plasmaIndexMap[m_itoSpecies[i]->getName()] = species;
  }
  for (int i = 0; i < m_cdrSpecies.size(); i++, species++) {
    m_plasmaIndexMap[m_cdrSpecies[i]->getName()] = species;
  }

  if (m_plasmaIndexMap.size() != m_numPlasmaSpecies) {
    this->throwParserError(baseError + " but something went wrong when setting up the solver map");
  }
}

void
ItoKMCJSON::initializeTownsendCoefficient(const std::string a_coeff)
{
  CH_TIME("ItoKMCJSON::initializeTownsendCoefficient");
  if (m_verbose) {
    pout() << m_className + "::initializeTownsendCoefficient" << endl;
  }

  const std::string baseError = "ItoKMCJSON::initializeTownsendCoefficient";

  FunctionEX func;

  if (!(m_json.contains(a_coeff))) {
    this->throwParserError(baseError + " but field '" + a_coeff + "' is missing");
  }
  if (!(m_json[a_coeff].contains("type"))) {
    this->throwParserError(baseError + " but field '" + a_coeff + "' does not contain the required field 'type'");
  }

  const std::string type = m_json[a_coeff]["type"].get<std::string>();

  if (type == "constant") {
    if (!(m_json[a_coeff].contains("value"))) {
      this->throwParserError(baseError + " and got 'constant' but missing the 'value' field");
    }

    const Real value = m_json[a_coeff]["value"].get<Real>();

    if (value < 0.0) {
      this->throwParserError(baseError + " and got 'constant' but can't have negative coefficient");
    }

    func = [value](const Real E, const RealVect x) -> Real {
      return value;
    };
  }
  else if (type == "table vs E/N") {
    const std::string baseErrorTable = baseError + " and got table vs E/N";

    const nlohmann::json& jsonTable = m_json[a_coeff];

    if (!(jsonTable.contains("file"))) {
      this->throwParserError(baseErrorTable + "but 'file' is not specified");
    }

    LookupTable1D<2> tabulatedCoeff = this->parseTableEByN(jsonTable, a_coeff + "/N");

    func = [this, tabulatedCoeff](const Real E, const RealVect x) -> Real {
      const Real N   = m_gasNumberDensity(x);
      const Real Etd = E / (N * Units::Td);

      return tabulatedCoeff.getEntry<1>(Etd) * N;
    };
  }
  else {
    this->throwParserError(baseError + "but type specification '" + type + "' is not supported");
  }

  // Associate with correct internal functions.
  if (a_coeff == "alpha") {
    m_alpha = func;
  }
  else if (a_coeff == "eta") {
    m_eta = func;
  }
}

void
ItoKMCJSON::initializeParticles()
{
  CH_TIME("ItoKMCJSON::initializeParticles");
  if (m_verbose) {
    pout() << m_className + "::initializeParticles" << endl;
  }

  const std::string baseError = "ItoKMCJSON::initializeParticles";

  for (const auto& species : m_json["plasma species"]) {
    const std::string speciesID   = species["id"].get<std::string>();
    const std::string baseErrorID = baseError + " and found 'initial particles' for species '" + speciesID + "'";

    List<PointParticle> initialParticles;

    if (species.contains("initial particles")) {

      for (const auto& initField : species["initial particles"]) {
        const auto obj = initField.get<nlohmann::json::object_t>();

        if (obj.size() != 1) {
          this->throwParserError(baseError + " - logic bust, too many entries!");
        }

        const std::string whichField = (*obj.begin()).first;

        if (whichField == "single particle") {
          unsigned long long weight   = 1ULL;
          RealVect           position = RealVect::Zero;

          if (!(initField["single particle"].contains("position"))) {
            this->throwParserError(baseError + " but 'single particle' field does not contain 'position'");
          }
          if (!(initField["single particle"].contains("weight"))) {
            this->throwParserError(baseError + " but 'single particle' field does not contain 'weight'");
          }

          for (int dir = 0; dir < SpaceDim; dir++) {
            position[dir] = initField["single particle"]["position"][dir].get<Real>();
          }
          weight = initField["single particle"]["weight"].get<unsigned long long>();

#ifdef CH_MPI
          if (procID() == 0) {
            initialParticles.add(PointParticle(position, 1.0 * weight));
          }
#else
          initialParticles.add(PointParticle(position, 1.0 * weight));
#endif
        }
        else if (whichField == "uniform distribution") {
          const nlohmann::json& jsonEntry = initField["uniform distribution"];

          if (!(jsonEntry.contains("low corner"))) {
            this->throwParserError(baseError + "but 'uniform distribution' does not contain 'low corner'");
          }
          if (!(jsonEntry.contains("high corner"))) {
            this->throwParserError(baseError + "but 'uniform distribution' does not contain 'high corner'");
          }
          if (!(jsonEntry.contains("num particles"))) {
            this->throwParserError(baseError + "but 'uniform distribution' does not contain 'num particles'");
          }
          if (!(jsonEntry.contains("weight"))) {
            this->throwParserError(baseError + "but 'uniform distribution' does not contain 'weight'");
          }

          RealVect loCorner = RealVect::Zero;
          RealVect hiCorner = RealVect::Zero;

          for (int dir = 0; dir < SpaceDim; dir++) {
            loCorner[dir] = jsonEntry["low corner"][dir].get<Real>();
            hiCorner[dir] = jsonEntry["high corner"][dir].get<Real>();

            if (hiCorner[dir] < loCorner[dir]) {
              this->throwParserError(baseError + "but 'uniform distribution' can't have 'high corner' < 'low corner'");
            }
          }

          const unsigned long long numParticles   = jsonEntry["num particles"].get<unsigned long long>();
          const unsigned long long particleWeight = jsonEntry["weight"].get<unsigned long long>();

          if (numParticles > 0) {
            List<PointParticle> particles;

            ParticleManagement::drawBoxParticles(particles, numParticles, loCorner, hiCorner);

            for (ListIterator<PointParticle> lit(particles); lit.ok(); ++lit) {
              lit().weight() = 1.0 * particleWeight;
            }

            initialParticles.catenate(particles);
          }
        }
        else if (whichField == "sphere distribution") {
          const nlohmann::json& jsonEntry = initField["sphere distribution"];

          if (!(jsonEntry.contains("center"))) {
            this->throwParserError(baseError + "but 'sphere distribution' does not contain 'center'");
          }
          if (!(jsonEntry.contains("radius"))) {
            this->throwParserError(baseError + "but 'sphere distribution' does not contain 'radius'");
          }
          if (!(jsonEntry.contains("num particles"))) {
            this->throwParserError(baseError + "but 'sphere distribution' does not contain 'num particles'");
          }
          if (!(jsonEntry.contains("weight"))) {
            this->throwParserError(baseError + "but 'sphere distribution' does not contain 'weight'");
          }

          RealVect center;
          for (int dir = 0; dir < SpaceDim; dir++) {
            center[dir] = jsonEntry["center"][dir].get<Real>();
          }

          const Real               radius         = jsonEntry["radius"].get<Real>();
          const unsigned long long numParticles   = jsonEntry["num particles"].get<unsigned long long>();
          const unsigned long long particleWeight = jsonEntry["weight"].get<unsigned long long>();

          if (radius < 0.0) {
            this->throwParserError(baseError + "but 'sphere distribution' can't have 'radius' < 0");
          }

          if (numParticles > 0) {
            List<PointParticle> particles;

            ParticleManagement::drawSphereParticles(particles, numParticles, center, radius);

            for (ListIterator<PointParticle> lit(particles); lit.ok(); ++lit) {
              lit().weight() = 1.0 * particleWeight;
            }

            initialParticles.catenate(particles);
          }
        }
        else if (whichField == "gaussian distribution") {
          const nlohmann::json& jsonEntry = initField["gaussian distribution"];

          if (!(jsonEntry.contains("center"))) {
            this->throwParserError(baseError + "but 'gaussian distribution' does not contain 'center'");
          }
          if (!(jsonEntry.contains("radius"))) {
            this->throwParserError(baseError + "but 'gaussian distribution' does not contain 'radius'");
          }
          if (!(jsonEntry.contains("num particles"))) {
            this->throwParserError(baseError + "but 'gaussian distribution' does not contain 'num particles'");
          }
          if (!(jsonEntry.contains("weight"))) {
            this->throwParserError(baseError + "but 'gaussian distribution' does not contain 'weight'");
          }

          RealVect center;
          for (int dir = 0; dir < SpaceDim; dir++) {
            center[dir] = jsonEntry["center"][dir].get<Real>();
          }

          const Real               radius         = jsonEntry["radius"].get<Real>();
          const unsigned long long numParticles   = jsonEntry["num particles"].get<unsigned long long>();
          const unsigned long long particleWeight = jsonEntry["weight"].get<unsigned long long>();

          if (radius < 0.0) {
            this->throwParserError(baseError + "but 'gaussian distribution' can't have 'radius' < 0");
          }

          if (numParticles > 0) {
            List<PointParticle> particles;

            ParticleManagement::drawGaussianParticles(particles, numParticles, center, radius);

            for (ListIterator<PointParticle> lit(particles); lit.ok(); ++lit) {
              lit().weight() = 1.0 * particleWeight;
            }

            initialParticles.catenate(particles);
          }
        }
        else if (whichField == "list") {
          const nlohmann::json& jsonEntry = initField["list"];

          if (!(jsonEntry.contains("file"))) {
            this->throwParserError(baseError + "but 'list' does not contain 'file'");
          }

          const std::string f = this->trim(jsonEntry["file"].get<std::string>());

          unsigned int xcol = 0;
          unsigned int ycol = 1;
          unsigned int zcol = 2;
          unsigned int wcol = 3;

          if (jsonEntry.contains("x column")) {
            xcol = jsonEntry["x column"].get<unsigned int>();
          }
          if (jsonEntry.contains("y column")) {
            ycol = jsonEntry["y column"].get<unsigned int>();
          }
          if (jsonEntry.contains("z column")) {
            zcol = jsonEntry["z column"].get<unsigned int>();
          }
          if (jsonEntry.contains("w column")) {
            wcol = jsonEntry["w column"].get<unsigned int>();
          }

          List<PointParticle> particles = DataParser::readPointParticlesASCII(f, xcol, ycol, zcol, wcol);

#ifdef CH_MPI
          if (procID() == 0) {
            initialParticles.catenate(particles);
          }
#else
          initialParticles.catenate(particles);
#endif
        }
        else {
          this->throwParserError(baseError + " but specification '" + whichField + "' is not supported");
        }
      }
    }

    // Put the particles in the solvers.
    const SpeciesType& speciesType = m_plasmaSpeciesTypes.at(speciesID);

    switch (speciesType) {
    case SpeciesType::Ito: {
      const int idx = m_itoSpeciesMap.at(speciesID);

      List<ItoParticle>& solverParticles = m_itoSpecies[idx]->getInitialParticles();

      solverParticles.clear();

      for (ListIterator<PointParticle> lit(initialParticles); lit.ok(); ++lit) {
        solverParticles.add(ItoParticle(lit().weight(), lit().position()));
      }

      break;
    }
    case SpeciesType::CDR: {
      const int idx = m_cdrSpeciesMap.at(speciesID);

      List<PointParticle>& solverParticles = m_cdrSpecies[idx]->getInitialParticles();

      solverParticles.clear();

      solverParticles.catenate(initialParticles);

      break;
    }
    default: {
      MayDay::Error("ItoKMCJSON::initializeParticles - logic bust");

      break;
    }
    }
  }
}

void
ItoKMCJSON::initializeMobilities()
{
  CH_TIME("ItoKMCJSON::initializeMobilities");
  if (m_verbose) {
    pout() << m_className + "::initializeMobilities" << endl;
  }

  const std::string baseError = "ItoKMCJSON::initializeMobilities";

  m_mobilityFunctions.resize(m_numPlasmaSpecies);

  // Read in mobilities
  for (const auto& species : m_json["plasma species"]) {
    FunctionEX mobilityFunction = [](const Real E, const RealVect x) -> Real {
      return 0.0;
    };

    const std::string speciesID   = species["id"].get<std::string>();
    const std::string baseErrorID = baseError + " and found mobile species '" + speciesID + "'";

    // Check if the species is mobile.
    const bool isMobile = species["mobile"].get<bool>();
    if (isMobile) {
      if (!(species.contains("mobility"))) {
        this->throwParserError(baseErrorID + " but did not find the required field 'mobility'");
      }

      const nlohmann::json& mobilityJSON = species["mobility"];

      if (!(mobilityJSON).contains("type")) {
        this->throwParserError(baseErrorID + " but 'type' specifier was not found");
      }

      const std::string type = mobilityJSON["type"].get<std::string>();

      if (type == "constant") {
        if (!(mobilityJSON.contains("value"))) {
          this->throwParserError(baseErrorID + " and got constant mobility but did not find the 'value' field");
        }

        const Real mu = mobilityJSON["value"].get<Real>();

        if (mu < 0.0) {
          this->throwParserError(baseErrorID + " and got constant mobility but mobility should not be negative");
        }

        mobilityFunction = [mu](const Real E, const RealVect x) -> Real {
          return mu;
        };
      }
      else if (type == "table vs E/N") {
        const std::string baseErrorTable = baseErrorID + " and also got table vs E/N";

        if (!(mobilityJSON.contains("file"))) {
          this->throwParserError(baseErrorTable + ", but 'file' is not specified");
        }

        LookupTable1D<2> tabulatedCoeff = this->parseTableEByN(mobilityJSON, "mu*N");

        mobilityFunction = [this, tabulatedCoeff](const Real E, const RealVect x) -> Real {
          const Real N   = m_gasNumberDensity(x);
          const Real Etd = E / (N * Units::Td);

          return tabulatedCoeff.getEntry<1>(Etd) / (std::numeric_limits<Real>::epsilon() + N);
        };
      }
      else {
        this->throwParserError(baseErrorID + " but mobility specifier '" + type + "' is not supported");
      }
    }

    // Put the mobility function into the appropriate solver.
    const int idx = m_plasmaIndexMap[speciesID];

    m_mobilityFunctions[idx] = mobilityFunction;
  }
}

void
ItoKMCJSON::initializeDiffusionCoefficients()
{
  CH_TIME("ItoKMCJSON::initializeDiffusionCoefficients");
  if (m_verbose) {
    pout() << m_className + "::initializeDiffusionCoefficients" << endl;
  }

  const std::string baseError = "ItoKMCJSON::initializeDiffusionCoefficients";

  m_diffusionCoefficients.resize(m_numPlasmaSpecies);

  for (const auto& species : m_json["plasma species"]) {
    FunctionEX diffusionCoefficient = [](const Real E, const RealVect x) -> Real {
      return 0.0;
    };

    const std::string speciesID   = species["id"].get<std::string>();
    const std::string baseErrorID = baseError + " and found diffusive species '" + speciesID + "'";

    // Check if the species is diffusive
    const bool isDiffusive = species["diffusive"].get<bool>();
    if (isDiffusive) {
      if (!(species.contains("diffusion"))) {
        this->throwParserError(baseErrorID + " but did not find the required field 'diffusion'");
      }

      const nlohmann::json& diffusionJSON = species["diffusion"];

      if (!(diffusionJSON).contains("type")) {
        this->throwParserError(baseErrorID + " but 'type' specifier was not found");
      }

      const std::string type = diffusionJSON["type"].get<std::string>();

      if (type == "constant") {
        if (!(diffusionJSON.contains("value"))) {
          this->throwParserError(baseErrorID + " and got constant diffusion but did not find the 'value' field");
        }

        const Real D = diffusionJSON["value"].get<Real>();

        if (D < 0.0) {
          this->throwParserError(baseErrorID + " and got constant diffusion but coefficient should not be negative");
        }

        diffusionCoefficient = [D](const Real E, const RealVect x) -> Real {
          return D;
        };
      }
      else if (type == "table vs E/N") {
        const std::string baseErrorTable = baseErrorID + " and also got table vs E/N";

        if (!(diffusionJSON.contains("file"))) {
          this->throwParserError(baseErrorTable + ", but 'file' is not specified");
        }

        LookupTable1D<2> tabulatedCoeff = this->parseTableEByN(diffusionJSON, "D*N");

        diffusionCoefficient = [this, tabulatedCoeff](const Real E, const RealVect x) -> Real {
          const Real N   = m_gasNumberDensity(x);
          const Real Etd = E / (N * Units::Td);

          return tabulatedCoeff.getEntry<1>(Etd) / (std::numeric_limits<Real>::epsilon() + N);
        };
      }
      else {
        this->throwParserError(baseErrorID + " but mobility specifier '" + type + "' is not supported");
      }
    }

    // Put the mobility function into the appropriate solver.
    const int idx = m_plasmaIndexMap[speciesID];

    m_diffusionCoefficients[idx] = diffusionCoefficient;
  }
}

void
ItoKMCJSON::initializeTemperatures()
{
  CH_TIME("ItoKMCJSON::initializeTemperatures");
  if (m_verbose) {
    pout() << m_className + "::initializeTemperatures" << endl;
  }

  const std::string baseError = "ItoKMCJSON::initializeTemperatures";

  m_plasmaTemperatures.resize(m_numPlasmaSpecies);

  for (const auto& species : m_json["plasma species"]) {
    FunctionEX temperature = [](const Real E, const RealVect x) -> Real {
      return 0.0;
    };

    const std::string speciesID   = species["id"].get<std::string>();
    const std::string baseErrorID = baseError + " and found 'temperature' for species '" + speciesID + "'";

    if (species.contains("temperature")) {
      const nlohmann::json& temperatureJSON = species["temperature"];

      if (!(temperatureJSON).contains("type")) {
        this->throwParserError(baseErrorID + " but 'type' specifier was not found");
      }

      const std::string type = temperatureJSON["type"].get<std::string>();

      if (type == "gas") {
        temperature = [T = this->m_gasTemperature](const Real E, const RealVect x) -> Real {
          return T(x);
        };
      }
      else if (type == "constant") {
        if (!(temperatureJSON.contains("T"))) {
          this->throwParserError(baseErrorID + " and got constant temperature but did not find the 'T' field");
        }

        const Real T = temperatureJSON["T"].get<Real>();

        if (T < 0.0) {
          this->throwParserError(baseErrorID + " and got constant temperature but 'T' should not be negative");
        }

        temperature = [T](const Real E, const RealVect x) -> Real {
          return T;
        };
      }
      else if (type == "table vs E/N") {
        const std::string baseErrorTable = baseErrorID + " and also got table vs E/N";

        if (!(temperatureJSON.contains("file"))) {
          this->throwParserError(baseErrorTable + ", but 'file' is not specified");
        }

        LookupTable1D<2> tabulatedCoeff = this->parseTableEByN(temperatureJSON, "eV");

        constexpr Real eVToKelvin = 2.0 * Units::Qe / (3.0 * Units::kb);

        temperature = [this, eVToKelvin, tabulatedCoeff](const Real E, const RealVect x) -> Real {
          const Real N   = m_gasNumberDensity(x);
          const Real Etd = E / (N * Units::Td);

          return eVToKelvin * tabulatedCoeff.getEntry<1>(Etd) / (std::numeric_limits<Real>::epsilon() + N);
        };
      }
      else {
        this->throwParserError(baseErrorID + " but temperature specifier '" + type + "' is not supported");
      }
    }
    else {
      temperature = [T = this->m_gasTemperature](const Real E, const RealVect x) -> Real {
        return T(x);
      };
    }

    // Put the mobility function into the appropriate solver.
    const int idx = m_plasmaIndexMap[speciesID];

    m_plasmaTemperatures[idx] = temperature;
  }
}

void
ItoKMCJSON::initializePhotonSpecies()
{
  CH_TIME("ItoKMCJSON::initializePhotonSpecies");
  if (m_verbose) {
    pout() << m_className + "::initializePhotonSpecies" << endl;
  }

  const std::string baseError = "ItoKMCJSON::initializePhotonSpecies";

  for (const auto& species : m_json["photon species"]) {
    if (!(species.contains("id"))) {
      this->throwParserError(baseError + " but one of the photon species do not contain the 'id' field");
    }

    const std::string speciesID   = species["id"].get<std::string>();
    const std::string baseErrorID = baseError + " for species '" + speciesID + "'";

    FunctionX kappaFunction = [](const RealVect a_pos) -> Real {
      return 0.0;
    };

    if (this->containsWildcard(speciesID)) {
      this->throwParserError(baseErrorID + " but species name '" + speciesID + "' should not contain wildcard @");
    }
    if (this->containsBracket(speciesID)) {
      this->throwParserError(baseError + " but species '" + speciesID + "' should not contain brackets");
    }
    if (m_allSpecies.count(speciesID) != 0) {
      this->throwParserError(baseErrorID + " but species '" + speciesID + "' was already defined elsewhere");
    }

    if (!(species.contains("kappa"))) {
      this->throwParserError(baseErrorID + " but 'kappa' is not specified");
    }

    const nlohmann::json& kappaJSON = species["kappa"];

    const std::string type = this->trim(kappaJSON["type"].get<std::string>());
    if (type == "constant") {
      if (!(kappaJSON.contains("value"))) {
        this->throwParserError(baseErrorID + " and got constant kappa but 'value' field is not specified");
      }

      const Real value = kappaJSON["value"].get<Real>();
      if (value < 0.0) {
        this->throwParserError(baseErrorID + " and got constant kappa but 'value' field can not be negative");
      }

      kappaFunction = [value](const RealVect a_pos) -> Real {
        return value;
      };
    }
    else if (type == "stochastic A") {
      if (!(kappaJSON.contains("f1"))) {
        this->throwParserError(baseErrorID + " and got 'stochastic A' but field 'f1' is not specified");
      }
      if (!(kappaJSON.contains("f2"))) {
        this->throwParserError(baseErrorID + " and got 'stochastic A' but field 'f2' is not specified");
      }
      if (!(kappaJSON.contains("chi min"))) {
        this->throwParserError(baseErrorID + " and got 'stochastic A' but field 'chi min' is not specified");
      }
      if (!(kappaJSON.contains("chi max"))) {
        this->throwParserError(baseErrorID + " and got 'stochastic A' but field 'chi max' is not specified");
      }
      if (!(kappaJSON.contains("neutral"))) {
        this->throwParserError(baseErrorID + " and got 'stochastic A' but field 'neutral' is not specified");
      }

      const Real        f1      = kappaJSON["f1"].get<Real>();
      const Real        f2      = kappaJSON["f2"].get<Real>();
      const Real        chiMin  = kappaJSON["chi min"].get<Real>();
      const Real        chiMax  = kappaJSON["chi max"].get<Real>();
      const std::string neutral = this->trim(kappaJSON["neutral"].get<std::string>());

      if (f1 >= f2) {
        this->throwParserError(baseErrorID + " and got 'stochastic A' but can't have f1 >= f2");
      }
      if (chiMin >= chiMax) {
        this->throwParserError(baseErrorID + " and got 'stochastic A' but can't have 'chi min' >= 'chi max'");
      }
      if (m_backgroundSpeciesMap.count(neutral) != 1) {
        this->throwParserError(baseErrorID + " and got 'stochastic A' but don't now species '" + neutral + "'");
      }

      std::uniform_real_distribution<Real> udist(f1, f2);

      kappaFunction = [f1,
                       f2,
                       udist,
                       x1              = chiMin,
                       x2              = chiMax,
                       &gasPressure    = this->m_gasPressure,
                       &neutralSpecies = this->m_backgroundSpecies[m_backgroundSpeciesMap.at(neutral)]](
                        const RealVect a_position) mutable -> Real {
        const Real m = neutralSpecies.molarFraction(a_position);
        const Real P = gasPressure(a_position);
        const Real p = m * P;

        const Real f  = Random::get(udist);
        const Real a  = (f - f1) / (f2 - f1);
        const Real K1 = x1 * p;
        const Real K2 = x2 * p;

        return K1 * std::pow(K2 / K1, a);
      };
    }
    else {
      this->throwParserError(baseErrorID + " but type specification '" + type + "' is not supported");
    }

    m_photonIndexMap[speciesID] = m_rtSpecies.size();

    m_rtSpecies.push_back(RefCountedPtr<RtSpecies>(new ItoKMCPhotonSpecies(speciesID, kappaFunction)));
  }
}

void
ItoKMCJSON::initializePlasmaReactions()
{
  CH_TIME("ItoKMCJSON::initializePlasmaReactions");
  if (m_verbose) {
    pout() << m_className + "::initializePlasmaReactions" << endl;
  }

  const std::string baseError = "ItoKMCJSON::initializePlasmaReactions";

  for (const auto& reactionJSON : m_json["plasma reactions"]) {
    if (!(reactionJSON.contains("reaction"))) {
      this->throwParserError(baseError + " but one of the reactions is missing the field 'reaction'");
    }
    if (!(reactionJSON.contains("type"))) {
      this->throwParserError(baseError + " but one of the reactions is missing the field 'type'");
    }

    const std::string reaction    = this->trim(reactionJSON["reaction"].get<std::string>());
    const std::string baseErrorID = baseError + " for reaction '" + reaction + "'";

    // Parse the reaction string to figure out the species involved in the reaction. This CAN involve the species wildcard, in which
    // case we also build the reaction superset;
    std::vector<std::string> reactants;
    std::vector<std::string> products;

    this->parseReactionString(reactants, products, reaction);

    // Reactions may have contained a wildcard. Build a superset if it does.
    const auto reactionSets = this->parseReactionWildcards(reactants, products, reactionJSON);

    for (const auto& curReaction : reactionSets) {
      const std::string              wildcard     = std::get<0>(curReaction);
      const std::vector<std::string> curReactants = std::get<1>(curReaction);
      const std::vector<std::string> curProducts  = std::get<2>(curReaction);

      // This is the reaction index for the current index. The reaction we are currently
      // dealing with is put in m_kmcReactions[reactionIdex].
      const int reactionIndex = m_kmcReactions.size();

      // Ignore species which are enclosed by brackets ()
      std::vector<std::string> trimmedProducts;
      for (const auto& p : curProducts) {
        if (!(this->isBracketed(p))) {
          trimmedProducts.emplace_back(p);
        }
      }

      // Make sure the reaction makes sense
      this->sanctifyPlasmaReaction(curReactants, trimmedProducts, reaction);

      // Build the KMC reaction. Note that this does not involve the "background species", which are absorbed into the transition rates.
      std::list<size_t> backgroundReactants;
      std::list<size_t> plasmaReactants;
      std::list<size_t> photonReactants;
      std::list<size_t> backgroundProducts;
      std::list<size_t> plasmaProducts;
      std::list<size_t> photonProducts;

      this->getReactionSpecies(backgroundReactants,
                               plasmaReactants,
                               photonReactants,
                               backgroundProducts,
                               plasmaProducts,
                               photonProducts,
                               curReactants,
                               trimmedProducts);

      // Figure out how to compute the rate for this reaction. The plasma reactants are put in here because
      // we must scale properly against the KMCDualStateReaction method (which operates using the microscopic rates).
      const auto reactionRates      = this->parsePlasmaReactionRate(reactionJSON, backgroundReactants, plasmaReactants);
      const auto reactionPlot       = this->parsePlasmaReactionPlot(reactionJSON);
      const auto gradientCorrection = this->parsePlasmaReactionGradientCorrection(reactionJSON);

      m_kmcReactions.emplace_back(std::make_shared<KMCReaction>(plasmaReactants, plasmaProducts, photonProducts));
      m_kmcReactionRates.emplace_back(reactionRates.first);
      m_fluidRates.emplace_back(reactionRates.second);
      m_kmcReactionRatePlots.emplace_back(reactionPlot);
      m_kmcReactionGradientCorrections.emplace_back(gradientCorrection);
    }
  }
}

void
ItoKMCJSON::initializePhotoReactions()
{
  CH_TIME("ItoKMCJSON::initializePhotoReactions");
  if (m_verbose) {
    pout() << m_className + "::initializePhotoReactions" << endl;
  }

  const std::string baseError = "ItoKMCJSON::initializePhotoReactions";

  Real photoiEfficiency;

  for (const auto& reactionJSON : m_json["photoionization"]) {
    if (!(reactionJSON.contains("reaction"))) {
      this->throwParserError(baseError + " but one of the reactions is missing the field 'reaction'");
    }
    if (!(reactionJSON.contains("efficiency"))) {
      photoiEfficiency = 1.0;
    }
    else {
      photoiEfficiency = reactionJSON["efficiency"].get<Real>();
    }

    const std::string reaction    = this->trim(reactionJSON["reaction"].get<std::string>());
    const std::string baseErrorID = baseError + " for reaction '" + reaction + "'";

    // Parse the reaction string to figure out the species involved in the reaction. This CAN involve the species wildcard, in which
    // case we also build the reaction superset;
    std::vector<std::string> reactants;
    std::vector<std::string> products;

    this->parseReactionString(reactants, products, reaction);

    // Reactions may have contained a wildcard. Build a superset if it does.
    const auto reactionSets = this->parseReactionWildcards(reactants, products, reactionJSON);

    for (const auto& curReaction : reactionSets) {
      const std::string              wildcard     = std::get<0>(curReaction);
      const std::vector<std::string> curReactants = std::get<1>(curReaction);
      const std::vector<std::string> curProducts  = std::get<2>(curReaction);

      // This is the reaction index for the current index. The reaction we are currently
      // dealing with is put in m_kmcReactions[reactionIdex].
      const int reactionIndex = m_kmcReactions.size();

      // Ignore species which are enclosed by brackets ()
      std::vector<std::string> trimmedReactants;
      std::vector<std::string> trimmedProducts;
      for (const auto& r : curReactants) {
        if (!(this->isBracketed(r))) {
          trimmedReactants.emplace_back(r);
        }
      }
      for (const auto& p : curProducts) {
        if (!(this->isBracketed(p))) {
          trimmedProducts.emplace_back(p);
        }
      }

      // Make sure the reaction makes sense
      this->sanctifyPhotoReaction(trimmedReactants, trimmedProducts, reaction);

      // Get all species involved.
      std::list<size_t> backgroundReactants;
      std::list<size_t> plasmaReactants;
      std::list<size_t> photonReactants;
      std::list<size_t> backgroundProducts;
      std::list<size_t> plasmaProducts;
      std::list<size_t> photonProducts;

      this->getReactionSpecies(backgroundReactants,
                               plasmaReactants,
                               photonReactants,
                               backgroundProducts,
                               plasmaProducts,
                               photonProducts,
                               trimmedReactants,
                               trimmedProducts);

      // Add the reaction
      m_photoReactions.emplace_back(
        std::make_shared<ItoKMCPhotoReaction>(photonReactants.front(), plasmaProducts, photoiEfficiency));
    }
  }
}

void
ItoKMCJSON::initializeSurfaceReactions()
{
  CH_TIME("ItoKMCJSON::initializeSurfaceReactions");
  if (m_verbose) {
    pout() << m_className + "::initializeSurfaceReactions" << endl;
  }

  pout() << "ItoKMCJSON::initializeSurfaceReactions - not implemented" << endl;
}

void
ItoKMCJSON::sanctifyPlasmaReaction(const std::vector<std::string>& a_reactants,
                                   const std::vector<std::string>& a_products,
                                   const std::string&              a_reaction) const noexcept
{
  CH_TIME("ItoKMCJSON::sanctifyPlasmaReaction()");
  if (m_verbose) {
    pout() << m_className + "::sanctifyPlasmaReaction()" << endl;
  }

  const std::string baseError = "ItoKMCJSON::sanctifyPlasmaReaction ";

  // All reactants must be in the list of neutral species or in the list of plasma species
  for (const auto& r : a_reactants) {
    const bool isBackground = this->isBackgroundSpecies(r);
    const bool isPlasma     = this->isPlasmaSpecies(r);
    if (!isBackground && !isPlasma) {
      this->throwParserError(baseError + " but reactant '" + r + "' for reaction '" + a_reaction +
                             " should not appear on left hand side");
    }
  }

  // All products should be in the list of plasma or photon species. It's ok if users include a neutral species -- we will ignore it (but tell the user about it).
  for (const auto& p : a_products) {
    const bool isBackground = this->isBackgroundSpecies(p);
    const bool isPlasma     = this->isPlasmaSpecies(p);
    const bool isPhoton     = this->isPhotonSpecies(p);

    if (!isBackground && !isPlasma && !isPhoton) {
      this->throwParserError(baseError + "but I do not know product species '" + p + "' for reaction '" + a_reaction +
                             "'.");
    }
  }

  // Check for charge conservation
  int sumCharge = 0;
  for (const auto& r : a_reactants) {
    if (this->isPlasmaSpecies(r)) {
      const SpeciesType& type = m_plasmaSpeciesTypes.at(r);

      int Z = 0;

      switch (type) {
      case SpeciesType::Ito: {
        const int idx = m_itoSpeciesMap.at(r);

        Z = m_itoSpecies[idx]->getChargeNumber();

        break;
      }
      case SpeciesType::CDR: {
        const int idx = m_cdrSpeciesMap.at(r);

        Z = m_cdrSpecies[idx]->getChargeNumber();

        break;
      }
      default: {
        const std::string err = baseError + " logic bust";

        MayDay::Error(err.c_str());

        break;
      }
      }

      sumCharge -= Z;
    }
  }
  for (const auto& p : a_products) {
    if (this->isPlasmaSpecies(p)) {
      const SpeciesType& type = m_plasmaSpeciesTypes.at(p);

      int Z = 0;

      switch (type) {
      case SpeciesType::Ito: {
        const int idx = m_itoSpeciesMap.at(p);

        Z = m_itoSpecies[idx]->getChargeNumber();

        break;
      }
      case SpeciesType::CDR: {
        const int idx = m_cdrSpeciesMap.at(p);

        Z = m_cdrSpecies[idx]->getChargeNumber();

        break;
      }
      default: {
        const std::string err = baseError + " logic bust";

        MayDay::Error(err.c_str());

        break;
      }
      }

      sumCharge += Z;
    }
  }

  if (sumCharge != 0) {
    this->throwParserWarning(baseError + " but charge is not conserved for reaction '" + a_reaction + "'.");
  }
}

void
ItoKMCJSON::sanctifyPhotoReaction(const std::vector<std::string>& a_reactants,
                                  const std::vector<std::string>& a_products,
                                  const std::string&              a_reaction) const noexcept
{
  CH_TIME("ItoKMCJSON::sanctifyPhotoReaction()");
  if (m_verbose) {
    pout() << m_className + "::sanctifyPhotoReaction()" << endl;
  }

  const std::string baseError = "ItoKMCJSON::sanctifyPhotoReaction ";

  // All reactants must be in the list of neutral species or in the list of plasma species
  for (const auto& r : a_reactants) {
    const bool isBackground = this->isBackgroundSpecies(r);
    const bool isPlasma     = this->isPlasmaSpecies(r);
    const bool isPhoton     = this->isPhotonSpecies(r);

    if (!isBackground && !isPlasma && !isPhoton) {
      this->throwParserError(baseError + " but reactant '" + r + "' for reaction '" + a_reaction +
                             " is not a background/plasma/photon species");
    }

    if (isPlasma) {
      this->throwParserError(baseError + " but reactant '" + r + "' for reaction '" + a_reaction +
                             " should not appear on left hand side");
    }
  }

  // All products should be in the list of plasma. It's ok if users include a neutral species -- we will ignore it (but tell the user about it).
  for (const auto& p : a_products) {
    const bool isBackground = this->isBackgroundSpecies(p);
    const bool isPlasma     = this->isPlasmaSpecies(p);
    const bool isPhoton     = this->isPhotonSpecies(p);

    if (isPhoton) {
      this->throwParserError(baseError + "but photon species '" + p + "' for reaction '" + a_reaction +
                             "' is not allowed on the right-hand side.");
    }

    if (!isBackground && !isPlasma) {
      this->throwParserError(baseError + "but I do not know product species '" + p + "' for reaction '" + a_reaction +
                             "'.");
    }
  }

  // Check for charge conservation
  int sumCharge = 0;

  for (const auto& p : a_products) {
    if (this->isPlasmaSpecies(p)) {
      const SpeciesType& type = m_plasmaSpeciesTypes.at(p);

      int Z = 0;

      switch (type) {
      case SpeciesType::Ito: {
        const int idx = m_itoSpeciesMap.at(p);

        Z = m_itoSpecies[idx]->getChargeNumber();

        break;
      }
      case SpeciesType::CDR: {
        const int idx = m_cdrSpeciesMap.at(p);

        Z = m_cdrSpecies[idx]->getChargeNumber();

        break;
      }
      default: {
        const std::string err = baseError + " logic bust";

        MayDay::Error(err.c_str());

        break;
      }
      }

      sumCharge += Z;
    }
  }

  if (sumCharge != 0) {
    this->throwParserWarning(baseError + " but charge is not conserved for reaction '" + a_reaction + "'.");
  }
}

void
ItoKMCJSON::parseReactionString(std::vector<std::string>& a_reactants,
                                std::vector<std::string>& a_products,
                                const std::string&        a_reaction) const noexcept
{
  CH_TIME("ItoKMCJSON::parseReactionString");
  if (m_verbose) {
    pout() << m_className + "::parseReactionString" << endl;
  }

  const std::string baseError = "ItoKMCJSON::parseReactionString";

  // Parse the string into segments. We split on white-space.
  std::stringstream        ss(a_reaction);
  std::string              segment;
  std::vector<std::string> segments;

  while (std::getline(ss, segment, ' ')) {
    segments.push_back(segment);
  }

  // Discard all whitespace and solitary + from input string
  segments.erase(std::remove(segments.begin(), segments.end(), ""), segments.end());
  segments.erase(std::remove(segments.begin(), segments.end(), "+"), segments.end());

  // Find the element containing "->"
  const auto& it = std::find(segments.begin(), segments.end(), "->");

  // Make sure that -> is in the reaction string.
  if (it == segments.end())
    this->throwParserError(baseError + " -- Reaction '" + a_reaction + "' does not contain '->");

  // Left of "->" are reactants and right of "->" are products
  a_reactants = std::vector<std::string>(segments.begin(), it);
  a_products  = std::vector<std::string>(it + 1, segments.end());
}

void
ItoKMCJSON::getReactionSpecies(std::list<size_t>&              a_backgroundReactants,
                               std::list<size_t>&              a_plasmaReactants,
                               std::list<size_t>&              a_photonReactants,
                               std::list<size_t>&              a_backgroundProducts,
                               std::list<size_t>&              a_plasmaProducts,
                               std::list<size_t>&              a_photonProducts,
                               const std::vector<std::string>& a_reactants,
                               const std::vector<std::string>& a_products) const noexcept
{
  CH_TIME("ItoKMCJSON::getReactionSpecies");
  if (m_verbose) {
    pout() << m_className + "::getReactionSpecies" << endl;
  }

  const std::string baseError = "ItoKMC::getReactionSpecies";

  a_backgroundReactants.clear();
  a_plasmaReactants.clear();
  a_photonReactants.clear();
  a_backgroundProducts.clear();
  a_plasmaProducts.clear();
  a_photonProducts.clear();

  for (const auto& r : a_reactants) {
    const bool isBackground = this->isBackgroundSpecies(r);
    const bool isPlasma     = this->isPlasmaSpecies(r);
    const bool isPhoton     = this->isPhotonSpecies(r);

    if (isBackground && !isPlasma && !isPhoton) {
      a_backgroundReactants.emplace_back(m_backgroundSpeciesMap.at(r));
    }
    else if (!isBackground && isPlasma && !isPhoton) {
      a_plasmaReactants.emplace_back(m_plasmaIndexMap.at(r));
    }
    else if (!isBackground && !isPlasma && isPhoton) {
      a_photonReactants.emplace_back(m_photonIndexMap.at(r));
    }
    else {
      this->throwParserError(baseError + " -- logic bust 1");
    }
  }

  for (const auto& r : a_products) {
    const bool isBackground = this->isBackgroundSpecies(r);
    const bool isPlasma     = this->isPlasmaSpecies(r);
    const bool isPhoton     = this->isPhotonSpecies(r);

    if (isBackground && !isPlasma && !isPhoton) {
      a_backgroundProducts.emplace_back(m_backgroundSpeciesMap.at(r));
    }
    else if (!isBackground && isPlasma && !isPhoton) {
      a_plasmaProducts.emplace_back(m_plasmaIndexMap.at(r));
    }
    else if (!isBackground && !isPlasma && isPhoton) {
      a_photonProducts.emplace_back(m_photonIndexMap.at(r));
    }
    else {
      this->throwParserError(baseError + " -- logic bust 2");
    }
  }
}

std::pair<std::function<Real(const Real E, const Real V, const RealVect x)>,
          std::function<Real(const Real E, const RealVect x)>>
ItoKMCJSON::parsePlasmaReactionRate(const nlohmann::json&    a_reactionJSON,
                                    const std::list<size_t>& a_backgroundReactants,
                                    const std::list<size_t>& a_plasmaReactants) const
{
  CH_TIME("ItoKMCJSON::parsePlasmaReactionRate");
  if (m_verbose) {
    pout() << m_className + "::parsePlasmaReactionRate" << endl;
  }

  // TLDR: ItoKMCPhysics uses KMCDualStateReaction which computes propensities for reactions S + S -> null as a = 0.5 * c * X * (X-1), where c is the
  //       "rate" in the KMC sense. This is correct since there are 0.5 * X * (X-1) distinct pairs of particles. But ItoKMCJSON expects that the input
  //       rates corresponding to the rates in the reaction rate equation, so for S + S -> null we would have dn/dt = -2*k*n*n, or dX/dt = -(2k/dV) * X * X.
  //       On the other hand, the deterministic limit of the tau-leaping equation becomes dX/dt = -2 * a = c*X*(X-1). For consistency we must therefore have
  //       c = 2*k/dV. Likewise, we would have c = k/dV for the reaction S1 + S2 => null (S1 and S2 being different species). In this latter case we do not
  //       need to account for the scaling.
  //
  //       This subtle scaling is important for consistency between the KMC algorithm and the reaction rate equation. Because KMCDualStateReaction operates
  //       using the microscopic rates, we must add this scaling back in. Also note that this scaling does not matter if background species enter on the left
  //       hand side because the rates are simply absorbed into the rate itself.

  FunctionEX fluidRate = [](const Real E, const RealVect x) -> Real {
    return 0.0;
  };

  // Count the number of times each reactant appears on the left hand side.
  std::map<size_t, size_t> reactantNumbers;
  for (const auto& r : a_plasmaReactants) {
    if (reactantNumbers.count(r) != 0) {
      reactantNumbers[r]++;
    }
    else {
      reactantNumbers[r] = 1;
    }
  }

  size_t volumeFactor     = 0;
  Real   propensityFactor = 1.0;
  for (const auto& rn : reactantNumbers) {
    for (size_t i = 2; i <= rn.second; i++) {
      propensityFactor *= i;
    }

    volumeFactor += rn.second;
  }

  const std::string type      = this->trim(a_reactionJSON["type"].get<std::string>());
  const std::string reaction  = this->trim(a_reactionJSON["reaction"].get<std::string>());
  const std::string baseError = "ItoKMC::parsePlasmaReactionRate for reaction '" + reaction + "'";

  if (type == "constant") {
    if (!(a_reactionJSON.contains("value"))) {
      this->throwParserError(baseError + " and got constant rate but 'rate' is not specified");
    }

    const Real value = a_reactionJSON["value"].get<Real>(); // * propensityFactor;
    if (value < 0.0) {
      this->throwParserError(baseError + " and got constant rate but 'value' cannot be negative");
    }

    fluidRate = [value](const Real E, const RealVect x) -> Real {
      return value;
    };
  }
  else if (type == "alpha*v") {
    if (!(a_reactionJSON.contains("species"))) {
      this->throwParserError(baseError + " and got 'alpha*v' but 'species' not not specified!");
    }

    const std::string species = a_reactionJSON["species"].get<std::string>();

    if (!(this->isPlasmaSpecies(species))) {
      this->throwParserError(baseError + "and got 'alpha*v' but species '" + species + "' is not a plasma species");
    }

    const int idx = m_plasmaIndexMap.at(species);

    fluidRate = [&mu = m_mobilityFunctions[idx], &alpha = m_alpha](const Real E, const RealVect x) -> Real {
      return alpha(E, x) * mu(E, x) * E;
    };
  }
  else if (type == "eta*v") {
    if (!(a_reactionJSON.contains("species"))) {
      this->throwParserError(baseError + " and got 'eta*v' but 'species' not not specified!");
    }

    const std::string species = a_reactionJSON["species"].get<std::string>();

    if (!(this->isPlasmaSpecies(species))) {
      this->throwParserError(baseError + "and got 'eta*v' but species '" + species + "' is not a plasma species");
    }

    const int idx = m_plasmaIndexMap.at(species);

    fluidRate = [&mu = m_mobilityFunctions[idx], &eta = m_eta](const Real E, const RealVect x) -> Real {
      return eta(E, x) * mu(E, x) * E;
    };
  }
  else if (type == "table vs E/N") {
    const std::string baseErrorTable = baseError + " and got table vs E/N";

    // Read in the table.
    if (!(a_reactionJSON.contains("file"))) {
      this->throwParserError(baseErrorTable + "but 'file' is not specified");
    }

    LookupTable1D<2> tabulatedCoeff = this->parseTableEByN(a_reactionJSON, "rate/N");

    fluidRate = [&N = this->m_gasNumberDensity, tabulatedCoeff](const Real E, const RealVect x) -> Real {
      const Real Etd = E / (N(x) * Units::Td);

      return tabulatedCoeff.getEntry<1>(Etd);
    };
  }
  else if (type == "function T A") {
    if (!(a_reactionJSON.contains("T"))) {
      this->throwParserError(baseError + "and got 'functionT A' but field 'T' was not found");
    }
    if (!(a_reactionJSON.contains("c1"))) {
      this->throwParserError(baseError + "and got 'functionT A' but field 'c2' was not found");
    }
    if (!(a_reactionJSON.contains("c2"))) {
      this->throwParserError(baseError + "and got 'functionT A' but field 'c2' was not found");
    }

    const std::string species = this->trim(a_reactionJSON["T"].get<std::string>());
    const Real        c1      = a_reactionJSON["c1"].get<Real>();
    const Real        c2      = a_reactionJSON["c2"].get<Real>();

    const bool isBackground = this->isBackgroundSpecies(species);
    const bool isPlasma     = this->isPlasmaSpecies(species);

    if (!isBackground || !isPlasma) {
      this->throwParserError(baseError + " but species '" + species + "' is not a background or plasma species");
    }

    FunctionEX speciesTemperature;
    if (isBackground) {
      speciesTemperature = [&backgroundTemperature = this->m_gasTemperature](const Real E, const RealVect x) -> Real {
        return backgroundTemperature(x);
      };
    }
    else if (isPlasma) {
      speciesTemperature = m_plasmaTemperatures[m_plasmaIndexMap.at(species)];
    }

    fluidRate = [c1, c2, T = speciesTemperature](const Real E, const RealVect x) -> Real {
      return c1 * std::pow(T(E, x), c2);
    };
  }
  else if (type == "function TT A") {
    if (!(a_reactionJSON.contains("T1"))) {
      this->throwParserError(baseError + " and got 'functionT1T2 A' but field 'T1' was not found");
    }
    if (!(a_reactionJSON.contains("T2"))) {
      this->throwParserError(baseError + " and got 'functionT1T2 A' but field 'T2' was not found");
    }
    if (!(a_reactionJSON.contains("c1"))) {
      this->throwParserError(baseError + " and got 'functionT1T2 A' but field 'c1' was not found");
    }
    if (!(a_reactionJSON.contains("c2"))) {
      this->throwParserError(baseError + " and got 'functionT1T2 A' but field 'c2' was not found");
    }

    const std::string speciesT1 = this->trim(a_reactionJSON["T1"].get<std::string>());
    const std::string speciesT2 = this->trim(a_reactionJSON["T2"].get<std::string>());
    const std::string err       = " and got function 'functionT1T2 A' but unrecognized species '";

    const bool isPlasmaT1 = this->isPlasmaSpecies(speciesT1);
    const bool isPlasmaT2 = this->isPlasmaSpecies(speciesT2);

    const bool isBackgroundT1 = this->isBackgroundSpecies(speciesT1);
    const bool isBackgroundT2 = this->isBackgroundSpecies(speciesT2);

    // Make sure that the specified species exist.
    if (!isPlasmaT1 && !isBackgroundT1) {
      this->throwParserError(baseError + err + speciesT1 + "'");
    }
    if (!isPlasmaT2 && !isBackgroundT2) {
      this->throwParserError(baseError + err + speciesT2 + "'");
    }

    FunctionEX speciesTemperature1;
    FunctionEX speciesTemperature2;

    if (isBackgroundT1) {
      speciesTemperature1 = [&backgroundTemperature = this->m_gasTemperature](const Real E, const RealVect x) -> Real {
        return backgroundTemperature(x);
      };
    }
    else {
      speciesTemperature1 = m_plasmaTemperatures[m_plasmaIndexMap.at(speciesT1)];
    }
    if (isBackgroundT2) {
      speciesTemperature2 = [&backgroundTemperature = this->m_gasTemperature](const Real E, const RealVect x) -> Real {
        return backgroundTemperature(x);
      };
    }
    else {
      speciesTemperature2 = m_plasmaTemperatures[m_plasmaIndexMap.at(speciesT2)];
    }

    const Real c1 = a_reactionJSON["c1"].get<Real>();
    const Real c2 = a_reactionJSON["c1"].get<Real>();

    fluidRate = [c1, c2, T1 = speciesTemperature1, T2 = speciesTemperature2](const Real E, const RealVect x) -> Real {
      return c1 * std::pow(T1(E, x) / T2(E, x), c2);
    };
  }
  else if (type == "function E/N exp A") {
    const std::string err1 = baseError + " and got 'function E/N exp A' but field '";
    const std::string err2 = "' is not specified";

    if (!(a_reactionJSON.contains("c1"))) {
      this->throwParserError(err1 + "c1" + err2);
    }
    if (!(a_reactionJSON.contains("c2"))) {
      this->throwParserError(err1 + "c2" + err2);
    }
    if (!(a_reactionJSON.contains("c3"))) {
      this->throwParserError(err1 + "c3" + err2);
    }
    if (!(a_reactionJSON.contains("c4"))) {
      this->throwParserError(err1 + "c4" + err2);
    }
    if (!(a_reactionJSON.contains("c5"))) {
      this->throwParserError(err1 + "c5" + err2);
    }

    // Get the constants
    const Real c1 = a_reactionJSON["c1"].get<Real>();
    const Real c2 = a_reactionJSON["c2"].get<Real>();
    const Real c3 = a_reactionJSON["c3"].get<Real>();
    const Real c4 = a_reactionJSON["c4"].get<Real>();
    const Real c5 = a_reactionJSON["c5"].get<Real>();

    fluidRate = [c1, c2, c3, c4, c5, &N = this->m_gasNumberDensity](const Real E, const RealVect x) -> Real {
      const Real ETd = E / (N(x) * Units::Td);

      return c1 * exp(-std::pow(c2 / (c3 + c4 * ETd), c5));
    };
  }
  else {
    this->throwParserError(baseError + " but 'type' specifier '" + type + "' is not supported");
  }

  // Scale the reaction according to various input variables.

  if (a_reactionJSON.contains("scale")) {
    const Real scale = a_reactionJSON["scale"].get<Real>();

    fluidRate = [fluidRate, scale](const Real E, const RealVect x) {
      return fluidRate(E, x) * scale;
    };
  }
  if (a_reactionJSON.contains("efficiency")) {
    const Real efficiency = a_reactionJSON["efficiency"].get<Real>();

    fluidRate = [fluidRate, efficiency](const Real E, const RealVect x) {
      return fluidRate(E, x) * efficiency;
    };
  }
  if (a_reactionJSON.contains("quenching pressure")) {
    const Real pq = a_reactionJSON["quenching pressure"].get<Real>();

    fluidRate = [fluidRate, pq, &p = this->m_gasPressure](const Real E, const RealVect x) {
      return fluidRate(E, x) * pq / (pq + p(x));
    };
  }
  if (a_reactionJSON.contains("quenching rates")) {
    const std::string derivedError = baseError + " and got 'quenching rates' but";

    if (!(a_reactionJSON["quenching rates"].contains("kr"))) {
      this->throwParserError(derivedError + " radiative rate 'kr' is missing");
    }
    if (!(a_reactionJSON["quenching rates"].contains("kp"))) {
      this->throwParserError(derivedError + " predissociation rate 'kp' is missing");
    }
    if (!(a_reactionJSON["quenching rates"].contains("kq/N"))) {
      this->throwParserError(derivedError + " quenching factor 'kq/N' is missing");
    }

    const Real kr    = a_reactionJSON["quenching rates"]["kr"].get<Real>();
    const Real kp    = a_reactionJSON["quenching rates"]["kp"].get<Real>();
    const Real kqByN = a_reactionJSON["quenching rates"]["kq/N"].get<Real>();

    fluidRate = [fluidRate, kr, kp, kqByN, &N = this->m_gasNumberDensity](const Real E, const RealVect x) -> Real {
      const Real kq = kqByN * N(x);

      return fluidRate(E, x) * (kr / (kr + kp + kq));
    };
  }

  // This is the KMC rate -- note that it absorbs the background species.
  FunctionEVX kmcRate = [fluidRate,
                         volumeFactor,
                         propensityFactor,
                         a_backgroundReactants,
                         &S = this->m_backgroundSpecies,
                         &N = this->m_gasNumberDensity](const Real E, const Real V, const RealVect x) -> Real {
    Real k = fluidRate(E, x);

    // Multiply by neutral densities
    for (const auto& idx : a_backgroundReactants) {
      const Real n = S[idx].molarFraction(x) * N(x);

      k *= n;
    }

    // Multiply by propensity factor (because of ItoKMCDualStateReaction)
    k *= propensityFactor;

    // Multiply by volume factor (for higher-order reactions)
    if (volumeFactor > 0) {
      k *= 1. / std::pow(V, volumeFactor - 1);
    }

    return k;
  };

  return std::make_pair(kmcRate, fluidRate);
}

std::pair<bool, std::string>
ItoKMCJSON::parsePlasmaReactionPlot(const nlohmann::json& a_reactionJSON) const
{
  CH_TIME("ItoKMCJSON::parsePlasmaReactionPlot");
  if (m_verbose) {
    pout() << m_className + "::parsePlasmaReactionPlot" << endl;
  }

  bool        plot = false;
  std::string id   = this->trim(a_reactionJSON["reaction"].get<std::string>());

  if (a_reactionJSON.contains("plot")) {
    plot = a_reactionJSON["plot"].get<bool>();
  }

  if (a_reactionJSON.contains("description")) {
    id = this->trim(a_reactionJSON["description"].get<std::string>());
  }

  return std::make_pair(plot, id);
}

std::pair<bool, std::string>
ItoKMCJSON::parsePlasmaReactionGradientCorrection(const nlohmann::json& a_reactionJSON) const
{
  CH_TIME("ItoKMCJSON::parsePlasmaReactionGradientCorrection");
  if (m_verbose) {
    pout() << m_className + "::parsePlasmaReactionGradientCorrection" << endl;
  }

  const std::string baseError = "ItoKMCJSON::parsePlasmaReactionGradientCorrection";

  std::pair<bool, std::string> ret = std::make_pair(false, "invalid");

  if (a_reactionJSON.contains("gradient correction")) {
    const std::string species = this->trim(a_reactionJSON["gradient correction"].get<std::string>());

    if (m_plasmaSpeciesTypes.count(species) == 0) {
      this->throwParserError(baseError + " but species '" + species + " is not a plasma species");
    }

    const SpeciesType type = m_plasmaSpeciesTypes.at(species);

    bool isMobile;
    bool isDiffusive;

    switch (type) {
    case SpeciesType::Ito: {
      const int idx = m_itoSpeciesMap.at(species);

      isMobile    = m_itoSpecies[idx]->isMobile();
      isDiffusive = m_itoSpecies[idx]->isDiffusive();

      break;
    }
    case SpeciesType::CDR: {
      const int idx = m_cdrSpeciesMap.at(species);

      isMobile    = m_cdrSpecies[idx]->isMobile();
      isDiffusive = m_cdrSpecies[idx]->isDiffusive();

      break;
    }
    default: {
      const std::string err = baseError + " - logic bust";

      this->throwParserError(err.c_str());

      break;
    }
    };

    if (isMobile && isDiffusive) {
      ret = std::make_pair(true, species);
    }
    else {
      this->throwParserError(baseError + " but species '" + species + "' is not mobile and diffusive!");
    }
  }

  return ret;
}

LookupTable1D<2>
ItoKMCJSON::parseTableEByN(const nlohmann::json& a_tableEntry, const std::string& a_dataID) const
{
  CH_TIME("ItoKMCJSON::parseTableEByN");
  if (m_verbose) {
    pout() << m_className + "::parseTableEByN" << endl;
  }

  const std::string preError  = "ItoKMCJSON::parseTableEByN";
  const std::string postError = "for dataID = " + a_dataID;

  if (!(a_tableEntry.contains("file"))) {
    this->throwParserError(preError + " but could not find the 'file' specifier " + postError);
  }

  const std::string fileName = this->trim(a_tableEntry["file"].get<std::string>());
  if (!(this->doesFileExist(fileName))) {
    this->throwParserError(preError + " but file '" + fileName + "' " + postError + " was not found");
  }

  int columnEbyN  = 0;
  int columnCoeff = 1;
  int numPoints   = 1000;

  TableSpacing spacing = TableSpacing::Exponential;

  if (a_tableEntry.contains("EbyN column")) {
    columnEbyN = a_tableEntry["EbyN column"].get<int>();
  }
  if (a_tableEntry.contains(a_dataID + " column")) {
    columnCoeff = a_tableEntry[a_dataID + " column"].get<int>();
  }
  if (a_tableEntry.contains("num points")) {
    numPoints = a_tableEntry["num points"];
  }
  if (a_tableEntry.contains("spacing")) {
    const std::string whichSpacing = this->trim(a_tableEntry["spacing"].get<std::string>());

    if (whichSpacing == "linear") {
      spacing = TableSpacing::Uniform;
    }
    else if (whichSpacing == "exponential") {
      spacing = TableSpacing::Exponential;
    }
    else {
      this->throwParserError(preError + " but spacing '" + whichSpacing + "' is not supported");
    }
  }

  LookupTable1D<2> tabulatedCoefficient;

  if (a_tableEntry.contains("header")) {
    const std::string header = this->trim(a_tableEntry["header"].get<std::string>());

    tabulatedCoefficient = DataParser::fractionalFileReadASCII(fileName, header, "", columnEbyN, columnCoeff);
  }
  else {
    tabulatedCoefficient = DataParser::simpleFileReadASCII(fileName, columnEbyN, columnCoeff);
  }

  // Scale table if asked
  if (a_tableEntry.contains("scale E/N")) {
    const Real scaling = a_tableEntry["scale E/N"].get<Real>();
    if (scaling <= 0.0) {
      this->throwParserWarning(preError + " but shouldn't have 'scale E/N' <= 0.0 for " + postError);
    }

    tabulatedCoefficient.scale<0>(scaling);
  }
  if (a_tableEntry.contains("scale " + a_dataID)) {
    const Real scaling = a_tableEntry["scale " + a_dataID].get<Real>();
    if (scaling <= 0.0) {
      this->throwParserWarning(preError + " but shouldn't have 'scale E/N' <= 0.0 for " + postError);
    }

    tabulatedCoefficient.scale<1>(scaling);
  }

  // Set min/max range for internal table
  if (a_tableEntry.contains("min E/N")) {
    const Real minEbyN = a_tableEntry["min E/N"].get<Real>();

    tabulatedCoefficient.setMinRange(minEbyN, 0);
  }
  if (a_tableEntry.contains("max E/N")) {
    const Real maxEbyN = a_tableEntry["max E/N"].get<Real>();

    tabulatedCoefficient.setMaxRange(maxEbyN, 0);
  }

  // Make the table uniform and meaningful
  tabulatedCoefficient.setTableSpacing(spacing);
  tabulatedCoefficient.sort(0);
  tabulatedCoefficient.makeUniform(numPoints);

  if (a_tableEntry.contains("dump")) {
    const std::string dumpFile = this->trim(a_tableEntry["dump"].get<std::string>());

    tabulatedCoefficient.dumpTable(dumpFile);
  }

  return tabulatedCoefficient;
}

std::list<std::tuple<std::string, std::vector<std::string>, std::vector<std::string>>>
ItoKMCJSON::parseReactionWildcards(const std::vector<std::string>& a_reactants,
                                   const std::vector<std::string>& a_products,
                                   const nlohmann::json&           a_reactionJSON) const noexcept
{
  CH_TIME("ItoKMCJSON::parseReactionWildcards()");
  if (m_verbose) {
    pout() << m_className + "::parseReactionWildcards()" << endl;
  }

  // This is what we return. A horrific creature.
  std::list<std::tuple<std::string, std::vector<std::string>, std::vector<std::string>>> reactionSets;

  // This is the reaction name.
  const std::string reaction  = a_reactionJSON["reaction"].get<std::string>();
  const std::string baseError = "ItoKMCJSON::parseReactionWildcards for reaction '" + reaction + "'";

  // Check if reaction string had a wildcard '@'. If it did we replace the wildcard with the corresponding species. This means that we need to
  // build additional reactions.
  const bool containsWildcard = this->containsWildcard(reaction);

  if (containsWildcard) {
    if (!(a_reactionJSON.contains("@"))) {
      this->throwParserError(baseError + "got reaction wildcard '@' but array '@:' was not specified");
    }

    // Get the wildcards array.
    const std::vector<std::string> wildcards = a_reactionJSON["@"].get<std::vector<std::string>>();

    // Go through the wildcards and replace appropriately.
    for (const auto& w : wildcards) {
      std::vector<std::string> curReactants;
      std::vector<std::string> curProducts;

      // Replace by wildcard in reactants.
      for (const auto& r : a_reactants) {
        if (this->containsWildcard(r)) {
          curReactants.emplace_back(w);
        }
        else {
          curReactants.emplace_back(r);
        }
      }

      // Replace by wildcard in reactants.
      for (const auto& p : a_products) {
        if (this->containsWildcard(p)) {
          curProducts.emplace_back(w);
        }
        else {
          curProducts.emplace_back(p);
        }
      }

      reactionSets.emplace_back(w, curReactants, curProducts);
    }
  }
  else {
    reactionSets.emplace_back("", a_reactants, a_products);
  }

  return reactionSets;
}

Real
ItoKMCJSON::computeDt(const RealVect a_E, const RealVect a_pos, const Vector<Real> a_densities) const noexcept
{
  CH_TIME("ItoKMCJSON::computeDt");
  if (m_verbose) {
    pout() << m_className + "::computeDt" << endl;
  }

  return std::numeric_limits<Real>::infinity();
}

Real
ItoKMCJSON::computeAlpha(const Real a_E, const RealVect a_pos) const noexcept
{
  CH_TIME("ItoKMCJSON::computeAlpha");
  if (m_verbose) {
    pout() << m_className + "::computeAlpha" << endl;
  }

  return m_alpha(a_E, a_pos);
}

Real
ItoKMCJSON::computeEta(const Real a_E, const RealVect a_pos) const noexcept
{
  CH_TIME("ItoKMCJSON::computeEta");
  if (m_verbose) {
    pout() << m_className + "::computeEta" << endl;
  }

  return m_eta(a_E, a_pos);
}

Vector<Real>
ItoKMCJSON::computeMobilities(const Real a_time, const RealVect a_pos, const RealVect a_E) const noexcept
{
  CH_TIME("ItoKMCJSON::computeMobilities");
  if (m_verbose) {
    pout() << m_className + "::computeMobilities" << endl;
  }

  const Real E = a_E.vectorLength();

  Vector<Real> mobilityCoefficients(m_numPlasmaSpecies, 0.0);
  for (int i = 0; i < m_numPlasmaSpecies; i++) {
    mobilityCoefficients[i] = m_mobilityFunctions[i](E, a_pos);
  }

  return mobilityCoefficients;
}

Vector<Real>
ItoKMCJSON::computeDiffusionCoefficients(const Real a_time, const RealVect a_pos, const RealVect a_E) const noexcept
{
  CH_TIME("ItoKMCJSON::computeDiffusionCoefficients");
  if (m_verbose) {
    pout() << m_className + "::computeDiffusionCoefficients" << endl;
  }

  const Real E = a_E.vectorLength();

  Vector<Real> diffusionCoefficients(m_numPlasmaSpecies, 0.0);
  for (int i = 0; i < m_numPlasmaSpecies; i++) {
    diffusionCoefficients[i] = m_diffusionCoefficients[i](E, a_pos);
  }

  return diffusionCoefficients;
}

void
ItoKMCJSON::updateReactionRates(const RealVect          a_E,
                                const RealVect          a_pos,
                                const Vector<Real>&     a_phi,
                                const Vector<RealVect>& a_gradPhi,
                                const Real              a_dx,
                                const Real              a_kappa) const noexcept
{
  CH_TIME("ItoKMCJSON::updateReactionRates");
  if (m_verbose) {
    pout() << m_className + "::updateReactionRates" << endl;
  }

#ifndef NDEBUG
  this->checkMolarFraction(a_pos);
#endif

  // Update basic reaction rates.
  const Real E = a_E.vectorLength();
  const Real V = std::pow(a_dx, SpaceDim);

  for (int i = 0; i < m_kmcReactions.size(); i++) {
    m_kmcReactions[i]->rate() = m_kmcReactionRates[i](E, V, a_pos);

    // Add gradient correction if the user has asked for it.
    const std::pair<bool, std::string> gradientCorrection = m_kmcReactionGradientCorrections[i];

    if (std::get<0>(gradientCorrection)) {
      const int idx = m_plasmaIndexMap.at(std::get<1>(gradientCorrection));

      const Real     n  = a_phi[idx];
      const Real     mu = m_mobilityFunctions[idx](E, a_pos);
      const Real     D  = m_diffusionCoefficients[idx](E, a_pos);
      const RealVect g  = a_gradPhi[idx];

      constexpr Real safety = std::numeric_limits<Real>::epsilon();

      Real fcorr = 1.0 + a_E.dotProduct(D * g) / (safety + n * mu * E * E);

      fcorr = std::max(fcorr, 0.0);
      fcorr = std::min(fcorr, 1.0);

      m_kmcReactions[i]->rate() *= fcorr;
    }
  }
}

void
ItoKMCJSON::secondaryEmissionEB(Vector<List<ItoParticle>>&       a_secondaryParticles,
                                Vector<Real>&                    a_secondaryCDRFluxes,
                                Vector<List<Photon>>&            a_secondaryPhotons,
                                const Vector<List<ItoParticle>>& a_primaryParticles,
                                const Vector<Real>&              a_primaryCDRFluxes,
                                const Vector<List<Photon>>&      a_primaryPhotons,
                                const RealVect&                  a_E,
                                const RealVect&                  a_cellCenter,
                                const RealVect&                  a_cellCentroid,
                                const RealVect&                  a_bndryCentroid,
                                const RealVect&                  a_bndryNormal,
                                const Real                       a_bndryArea,
                                const Real                       a_dx,
                                const Real                       a_dt,
                                const bool                       a_isDielectric,
                                const int                        a_matIndex) const noexcept
{
  CH_TIME("ItoKMCJSON::secondaryEmissionEB");
  if (m_verbose) {
    pout() << m_className + "::secondaryEmissionEB" << endl;
  }
}

int
ItoKMCJSON::getNumberOfPlotVariables() const noexcept
{
  CH_TIME("ItoKMCJSON::getNumberOfPlotVariables");
  if (m_verbose) {
    pout() << m_className + "::getNumberOfPlotVariables" << endl;
  }

  int numPlots = 0;

  if (m_plotGas) {
    numPlots += 3;
  }

  for (int i = 0; i < m_backgroundSpeciesPlot.size(); i++) {
    if (m_backgroundSpeciesPlot[i]) {
      numPlots++;
    }
  }

  for (int i = 0; i < m_kmcReactionRatePlots.size(); i++) {
    if (std::get<0>(m_kmcReactionRatePlots[i])) {
      numPlots++;
    }
  }

  return numPlots;
}

Vector<std::string>
ItoKMCJSON::getPlotVariableNames() const noexcept
{
  CH_TIME("ItoKMCJSON::getPlotVariableNames");
  if (m_verbose) {
    pout() << m_className + "::getPlotVariableNames" << endl;
  }

  Vector<std::string> plotVariableNames;

  if (m_plotGas) {
    plotVariableNames.push_back("gas pressure");
    plotVariableNames.push_back("gas temperature");
    plotVariableNames.push_back("gas number density");
  }

  for (int i = 0; i < m_backgroundSpeciesPlot.size(); i++) {
    if (m_backgroundSpeciesPlot[i]) {
      plotVariableNames.push_back(m_backgroundSpecies[i].getName() + " molar fraction");
    }
  }

  for (int i = 0; i < m_kmcReactionRatePlots.size(); i++) {
    if (std::get<0>(m_kmcReactionRatePlots[i])) {
      plotVariableNames.push_back(std::get<1>(m_kmcReactionRatePlots[i]));
    }
  }

  return plotVariableNames;
}

Vector<Real>
ItoKMCJSON::getPlotVariables(const RealVect          a_E,
                             const RealVect          a_pos,
                             const Vector<Real>&     a_phi,
                             const Vector<RealVect>& a_gradPhi,
                             const Real              a_dx,
                             const Real              a_kappa) const noexcept
{

  CH_TIME("ItoKMCJSON::getPlotVariables");
  if (m_verbose) {
    pout() << m_className + "::getPlotVariables" << endl;
  }

  Vector<Real> plotVars;

  // Update basic reaction rates.
  const Real E = a_E.vectorLength();

  if (m_plotGas) {
    plotVars.push_back(m_gasPressure(a_pos));
    plotVars.push_back(m_gasTemperature(a_pos));
    plotVars.push_back(m_gasNumberDensity(a_pos));
  }

  for (int i = 0; i < m_backgroundSpeciesPlot.size(); i++) {
    if (m_backgroundSpeciesPlot[i]) {
      plotVars.push_back(m_backgroundSpecies[i].molarFraction(a_pos));
    }
  }

  for (int i = 0; i < m_kmcReactions.size(); i++) {
    if (std::get<0>(m_kmcReactionRatePlots[i])) {
      plotVars.push_back(m_fluidRates[i](E, a_pos));
    }
  }

  return plotVars;
}

#include <CD_NamespaceFooter.H>
