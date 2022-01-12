/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaJSON.cpp
  @brief  Implementation of CD_CdrPlasmaJSON.H
  @author Robert Marskar
*/

// Std includees
#include <iostream>
#include <fstream>

// Chombo includes
#include <ParmParse.H>
#include <CH_Timer.H>

// Our includes
#include <CD_CdrPlasmaJSON.H>
#include <CD_DataParser.H>
#include <CD_Units.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::CdrPlasma;

CdrPlasmaJSON::CdrPlasmaJSON(){
  CH_TIME("CdrPlasmaJSON::CdrPlasmaJSON()");

  // Parse options and read JSON file.
  this->parseOptions();
  this->parseJSON();

  // Parse initial data.
  this->initializeNeutralSpecies();
  this->initializePlasmaSpecies();
  this->initializeSigma();

  // Parse CDR mobilities and diffusion coefficients. 
  this->parseMobilities();

  // Populate the stuff that is needed by CdrPlasmaPhysics
  m_RtSpecies.resize(0);

  m_numCdrSpecies = m_CdrSpecies.size();
  m_numRtSpecies  = m_RtSpecies. size();
}

CdrPlasmaJSON::~CdrPlasmaJSON(){
  CH_TIME("CdrPlasmaJSON::~CdrPlasmaJSON()");
}

void CdrPlasmaJSON::parseOptions() {
  CH_TIME("CdrPlasmaJSON::parseOptions");
  
  ParmParse pp("CdrPlasmaJSON");

  pp.get("verbose",        m_verbose );
  pp.get("chemistry_file", m_jsonFile);
}

void CdrPlasmaJSON::parseJSON() {
  CH_TIME("CdrPlasmaJSON::parseJSON");
  if(m_verbose){
    pout() << "CdrPlasmaJSON::parseJSON -- file is = " << m_jsonFile << endl;
  }

  // Parse the JSON file
  std::ifstream istream(m_jsonFile);
  istream >> m_json;  
}

void CdrPlasmaJSON::throwParserError(const std::string a_error) const {
  pout() << a_error << endl;

  MayDay::Error(a_error.c_str());
}

void CdrPlasmaJSON::throwParserWarning(const std::string a_warning) const {
  pout() << a_warning << endl;

  MayDay::Warning(a_warning.c_str());
}

void CdrPlasmaJSON::initializeSigma() {
  CH_TIME("CdrPlasmaJSON::initializeSigma");
  if(m_verbose){
    pout() << "CdrPlasmaJSON::initializeSigma()" << endl;
  }

  m_initialSigma = [](const Real a_time, const RealVect a_pos) -> Real {
    return 0.0;
  };

  if(m_json.contains("sigma")){
    if(m_json["sigma"].contains("initial_density")){
      const Real sigma = m_json["sigma"]["initial_density"].get<Real>();
      
      m_initialSigma = [sigma] (const Real a_time, const RealVect a_pos) -> Real {
	return sigma;
      };
    }
  }
}

void CdrPlasmaJSON::initializeNeutralSpecies() {
  CH_TIME("CdrNeutralJSON::initializeNeutralSpecies");
  if(m_verbose){
    pout() << "CdrNeutralJSON::initializeNeutralSpecies()" << endl;
  }

  // These fields are required
  if(!(m_json["gas"].contains("temperature"    ))) this->throwParserError("In JSON field 'gas' - field 'temperature' is missing"    );
  if(!(m_json["gas"].contains("pressure"       ))) this->throwParserError("In JSON field 'gas' - field 'pressure' is missing"       );
  if(!(m_json["gas"].contains("law"            ))) this->throwParserError("In JSON field 'gas' - field 'law' is missing"            );
  if(!(m_json["gas"].contains("neutral_species"))) this->throwParserError("In JSON field 'gas' - field 'neutral_species' is missing");

  const auto referenceTemperature = m_json["gas"]["temperature"].get<Real       >();
  const auto referencePressure    = m_json["gas"]["pressure"   ].get<Real       >();
  const auto gasLaw               = m_json["gas"]["law"        ].get<std::string>();

  // Instantiate the pressure, density, and temperature of the gas. Note: The density  is the NUMBER density. 
  if(gasLaw == "ideal"){
    // Set the gas temperature, density, and pressure from the ideal gas law. No extra parameters needed and no variation in space either. 
      
    const Real referenceDensity = (referencePressure * Units::atm2pascal * Units::Na)/ (referenceTemperature * Units::R);	

    m_gasTemperature = [T   = referenceTemperature] (const RealVect a_position) -> Real { return T;   };
    m_gasPressure    = [P   = referencePressure   ] (const RealVect a_position) -> Real { return P;   };
    m_gasDensity     = [Rho = referenceDensity    ] (const RealVect a_position) -> Real { return Rho; };
  }
  else if(gasLaw == "troposphere"){

    // These fields are required.  
    if(!(m_json["gas"].contains("molar_mass"      ))) this->throwParserError("gas_law is 'troposphere' but I did not find field 'molar_mass'"      );
    if(!(m_json["gas"].contains("gravity"         ))) this->throwParserError("gas_law is 'troposphere' but I did not find field 'gravity'"         );
    if(!(m_json["gas"].contains("lapse_rate"      ))) this->throwParserError("gas_law is 'troposphere' but I did not find field 'lapse_rate'"      );

    const Real g    = m_json["gas"]["gravity"   ].get<Real>();
    const Real L    = m_json["gas"]["lapse_rate"].get<Real>();
    const Real M    = m_json["gas"]["molar_mass"].get<Real>();
    const Real gMRL = (g * M) / (Units::R * L);

    // Temperature is T = T0 - L*(z-h0)
    m_gasTemperature = [T  = referenceTemperature, L] (const RealVect a_position) -> Real {
      return T - L * a_position[SpaceDim-1];
    };

    // Pressure is p = p0 * (1 - L*h/T0)^(g*M/(R*L)
    m_gasPressure = [T  = referenceTemperature, P = referencePressure, L, gMRL ] (const RealVect a_position) -> Real {
      return P * std::pow(( 1 - L*a_position[SpaceDim-1]/T), gMRL);
    };

    // Density is rho = P*Na/(T*R)
    m_gasDensity = [P = this->m_gasPressure, T = this->m_gasTemperature] (const RealVect a_position) -> Real {
      return (P(a_position) * Units::atm2pascal * Units::Na)/ (T(a_position) * Units::R);          
    };
  }
  else{
    this->throwParserError("CdrPlasmaJSON::initializeNeutralSpecies gas law '" + gasLaw + "' not recognized.");
  }

  // Instantiate the species densities. Note that we need to go through this twice because we need to normalize the molar fractions in case users
  // were a bit inconsiderate when setting them. 
  Real molarSum = 0.0;
  for (const auto& species : m_json["gas"]["neutral_species"]){
    if(!(species.contains("name")))           this->throwParserError("In JSON field 'neutral_species' - field 'name' is required"          );
    if(!(species.contains("molar_fraction"))) this->throwParserError("In JSON field 'neutral_species' - field 'molar_fraction' is required");

    molarSum += species["molar_fraction"].get<Real>();
  }

  // Initialize the species
  for (const auto& species : m_json["gas"]["neutral_species"]){
    const std::string speciesName     = species["name"          ].get<std::string>();
    const Real        speciesFraction = species["molar_fraction"].get<Real>() / molarSum;

    // Set the species density function. 
    const std::function<Real(const RealVect)> speciesDensity  = [f = speciesFraction, N = this->m_gasDensity] (const RealVect a_position) {
      return f * N(a_position);
    };

    // Add the species. Make sure the maps are consist.
    const int idx = m_neutralSpecies.size();
    
    // Create the neutral species.
    m_neutralSpecies.push_back(std::shared_ptr<NeutralSpeciesJSON>((new NeutralSpeciesJSON(speciesName, speciesDensity))));

    // Create the string-int maps
    m_neutralSpeciesMap.       insert(std::make_pair(speciesName, idx        ));
    m_neutralSpeciesInverseMap.insert(std::make_pair(idx ,        speciesName));
  }

  // Figure out if we should plot the gas quantities
  m_plotGas = false;
  if(m_json["gas"].contains("plot")){
    m_plotGas = m_json["gas"]["plot"].get<bool>();
  }
}

void CdrPlasmaJSON::initializePlasmaSpecies() {
  CH_TIME("CdrPlasmaJSON::initializePlasmaSpecies");
  if(m_verbose){
    pout() << "CdrPlasmaJSON::initializePlasmaSpecies()" << endl;
  }

  if(!(m_json.contains("plasma_species"))) this->throwParserWarning("CdrPlasmaJSON::initializeSpecies -- did not find any plasma species");

  // Iterate through all species defined in the JSON file. 
  for (const auto& species : m_json["plasma_species"]){
    if(!(species.contains("name"     ))) this->throwParserError("CdrPlasmaJSON::initializeSpecies -- 'plasma_species' must have field 'name'"     );
    if(!(species.contains("Z"        ))) this->throwParserError("CdrPlasmaJSON::initializeSpecies -- 'plasma_species' must have field 'Z'"        );
    if(!(species.contains("mobile"   ))) this->throwParserError("CdrPlasmaJSON::initializeSpecies -- 'plasma_species' must have field 'mobile'"   );
    if(!(species.contains("diffusive"))) this->throwParserError("CdrPlasmaJSON::initializeSpecies -- 'plasma_species' must have field 'diffusive'");

    const auto name        = species["name"].     get<std::string>();
    const auto Z           = species["Z"].        get<int        >();
    const auto mobile      = species["mobile"].   get<bool       >();
    const auto diffusive   = species["diffusive"].get<bool       >();
    
    const bool hasInitData = species.contains("initial_data");

    // Get the initial data. 
    std::function<Real(const RealVect, const Real)> initFunc;

    if(hasInitData){
      initFunc = [data = species["initial_data"]] (const RealVect a_point, const Real a_time) -> Real {
	Real ret = 0.0;

	// Add uniform density. 
	if(data.contains("uniform")){
	  ret += data["uniform"].get<Real>();
	}

	// Add gaussian seed.
	if(data.contains("gauss2")){
	  const Real     amplitude = data["gauss2"]["amplitude"].get<Real>();	  
	  const Real     radius    = data["gauss2"]["radius"   ].get<Real>();
	  const RealVect center    = RealVect(D_DECL(data["gauss2"]["position"][0].get<Real>(),
						     data["gauss2"]["position"][1].get<Real>(),
						     data["gauss2"]["position"][2].get<Real>()));
	  const RealVect delta      = center - a_point;

	  ret += amplitude * exp(-delta.dotProduct(delta)/(2*std::pow(radius,2)));
	}

	// Add super-Gaussian seed.
	if(data.contains("gauss4")){
	  const Real     amplitude = data["gauss4"]["amplitude"].get<Real>();	  
	  const Real     radius    = data["gauss4"]["radius"   ].get<Real>();
	  const RealVect center    = RealVect(D_DECL(data["gauss4"]["position"][0].get<Real>(),
						     data["gauss4"]["position"][1].get<Real>(),
						     data["gauss4"]["position"][2].get<Real>()));
	  const RealVect delta      = center - a_point;

	  ret += amplitude * exp(-std::pow(delta.dotProduct(delta),2)/(2*std::pow(radius, 4)));
	}

	return ret;
      };
    }
    else{
      initFunc = [](const RealVect a_point, const Real a_time) -> Real {
	return 0.0;
      };
    }

    // Print out a message if we're verbose.
    if(m_verbose){
      pout() << "CdrPlasmaJSON::initializeSpecies: instantiating species" << "\n"
	     << "\tName        = " << name        << "\n"
	     << "\tZ           = " << Z           << "\n"
	     << "\tMobile      = " << mobile      << "\n"
	     << "\tDiffusive   = " << diffusive   << "\n"
	     << "\tInitialData = " << hasInitData << "\n";    	
    }

    // Initialize the species.
    const int num = m_CdrSpecies.size();

    // Make the string-int map encodings. 
    m_cdrSpeciesMap.       emplace(std::make_pair(name, num ));
    m_cdrSpeciesInverseMap.emplace(std::make_pair(num , name));

    // Push the JSON entry and the new CdrSpecies to corresponding vectors. 
    m_CdrSpecies.    push_back(RefCountedPtr<CdrSpecies> (new CdrSpeciesJSON(name, Z, diffusive, mobile, initFunc)));
    m_cdrSpeciesJSON.push_back(species);
  }
}

void CdrPlasmaJSON::parseMobilities() {
  CH_TIME("CdrPlasmaJSON::parseMobilities");
  if(m_verbose){
    pout() << "CdrPlasmaJSON::parseMobilities - file is = " << m_jsonFile << endl;
  }

  // Iterate through all tracked species.
  for (const auto& species : m_cdrSpeciesJSON){
    const std::string name = species["name"].get<std::string>();
    const int         idx  = m_cdrSpeciesMap.at(name);

    if(m_CdrSpecies[idx]->isMobile()){

      // This is a required field. We use it for specifying the mobility. 
      if(!(species.contains("mobility"))) this->throwParserError("species " + name + " is mobile but file does not contain field 'mobility'");
      const json& mobilityJSON = species["mobility"];

      // Get the mobility lookup method. This must either be a constant, a function, or a table. We parse these cases differently.
      const std::string lookup = mobilityJSON["lookup"].get<std::string>();
	
      if(lookup == "constant"){
	// User specified a constant mobility. We look for a field 'value' in the JSON file and set the mobility from that. If the
	// field does not exist then it's an error. 
	if(!(mobilityJSON.contains("value"))) this->throwParserError("CdrPlasmaJSON::parseMobilities -- constant lookup was specified but I can't find the field 'value'");

	const Real value = mobilityJSON["value"].get<Real>();

	m_mobilityLookup.  emplace(std::make_pair(idx, LookupMethod::Constant));
	m_mobilityConstant.emplace(std::make_pair(idx, value                 ));
      }
      else if(lookup == "table E/N"){
	if(!(mobilityJSON.contains("file"      ))) this->throwParserError("CdrPlasmaJSON::parseMobilities -- tabulated mobility was specified but field 'file' was not"   );
	if(!(mobilityJSON.contains("header"    ))) this->throwParserError("CdrPlasmaJSON::parseMobilities -- tabulated mobility was specified but field 'header' was not" );
	if(!(mobilityJSON.contains("E/N"       ))) this->throwParserError("CdrPlasmaJSON::parseMobilities -- tabulated mobility was specified but field 'E/N' was not"    );
	if(!(mobilityJSON.contains("mu*N"      ))) this->throwParserError("CdrPlasmaJSON::parseMobilities -- tabulated mobility was specified but field 'mu*N' was not"   );
	if(!(mobilityJSON.contains("min E/N"   ))) this->throwParserError("CdrPlasmaJSON::parseMobilities -- tabulated mobility was specified but field 'min E/N' was not");
	if(!(mobilityJSON.contains("max E/N"   ))) this->throwParserError("CdrPlasmaJSON::parseMobilities -- tabulated mobility was specified but field 'max E/N' was not");
	if(!(mobilityJSON.contains("num_points"))) this->throwParserError("CdrPlasmaJSON::parseMobilities -- tabulated mobility was specified but field 'num_points' was not");
	
	const std::string filename  = mobilityJSON["file"  ].get<std::string>();
	const std::string startRead = mobilityJSON["header"].get<std::string>();
	const std::string stopRead  = "";

	const int xColumn   = mobilityJSON["E/N"       ].get<int>();
	const int yColumn   = mobilityJSON["mu*N"      ].get<int>();
	const int numPoints = mobilityJSON["num_points"].get<int>();

	const Real minEN = mobilityJSON["min E/N"].get<Real>();
	const Real maxEN = mobilityJSON["max E/N"].get<Real>();	

	// Read the table and format it. We happen to know that this function reads data into the approprate columns. So if
	// the user specified the correct E/N column then that data will be put in the first column. The data for mu*N will be in the
	// second column. 
	LookupTable<2> mobilityTable = DataParser::fractionalFileReadASCII(filename, startRead, stopRead, xColumn, yColumn);

	// Format the table
	mobilityTable.setRange(minEN, maxEN, 0);
	mobilityTable.sort(0);
	mobilityTable.makeUniform(numPoints);

	m_mobilityLookup.emplace(std::make_pair(idx, LookupMethod::TableLFA));
	m_mobilityTables.emplace(std::make_pair(idx, mobilityTable         ));
      }		
      else if (lookup == "function E/N"){
	if(!(mobilityJSON.contains("function"))) this->throwParserError("CdrPlasmaJSON::parseMobilities -- using function E/N but did not find function specification");

	FunctionEN func;
	
	const std::string whichFunction = mobilityJSON["function"].get<std::string>();
	if(whichFunction == "ABC"){
	  if(!(mobilityJSON.contains("A"))) this->throwParserError("CdrPlasmaJSON::parseMobilities -- using function 'ABC'  did not find 'A'");
	  if(!(mobilityJSON.contains("B"))) this->throwParserError("CdrPlasmaJSON::parseMobilities -- using function 'ABC'  did not find 'B'");
	  if(!(mobilityJSON.contains("C"))) this->throwParserError("CdrPlasmaJSON::parseMobilities -- using function 'ABC'  did not find 'C'");

	  const Real A = mobilityJSON["A"].get<Real>();
	  const Real B = mobilityJSON["B"].get<Real>();
	  const Real C = mobilityJSON["C"].get<Real>();

	  func = [A, B, C](const Real a_E, const Real a_N) -> Real {return A * std::pow(a_E, B) / std::pow(a_N, C);};
	}
	else{
	  this->throwParserError("CdrPlasmaJSON::parseMobilities -- using function E/N but do not recognize specification '" + whichFunction + "'");
	}

	m_mobilityLookup.     emplace(std::make_pair(idx, LookupMethod::FunctionLFA));
	m_mobilityFunctionsEN.emplace(std::make_pair(idx, func                     ));	
      }
      else{
	this->throwParserError("CdrPlasmaJSON::parseMobilities -- lookup = '" + lookup + "' was specified but this is not 'constant', 'function', or 'table'");
      }
    }
  }
}

int CdrPlasmaJSON::getNumberOfPlotVariables() const {
  int ret = 0;

  if(m_plotGas){
    ret = 3;
  }

  return ret;
}

Vector<std::string> CdrPlasmaJSON::getPlotVariableNames() const {
  Vector<std::string> ret(0);

  if(m_plotGas){
    ret.push_back("gas pressure"   );
    ret.push_back("gas temperature");
    ret.push_back("gas density"    );
  }

  return ret;
}

Vector<Real> CdrPlasmaJSON::getPlotVariables(const Vector<Real> a_cdrDensities,
					     const Vector<Real> a_rteDensities,
					     const RealVect     a_E,
					     const RealVect     a_position,
					     const Real         a_time) const {
  Vector<Real> ret(0);

  if(m_plotGas){
    ret.push_back(m_gasPressure   (a_position));
    ret.push_back(m_gasTemperature(a_position));
    ret.push_back(m_gasDensity    (a_position));
  }

  return ret;
}

Real CdrPlasmaJSON::getNeutralSpeciesDensity(const std::string a_name, const RealVect& a_position) const {
  return (*m_neutralSpecies[m_neutralSpeciesMap.at(a_name)])(a_position);
}

Real CdrPlasmaJSON::computeAlpha(const RealVect a_E) const {
  MayDay::Warning("CdrPlasmaJSON::computeAlpha -- don't know how to do this yet");

  return 0.0;
}

void CdrPlasmaJSON::advanceReactionNetwork(Vector<Real>&          a_cdrSources,
					   Vector<Real>&          a_rteSources,
					   const Vector<Real>     a_cdrDensities,
					   const Vector<RealVect> a_cdrGradients,
					   const Vector<Real>     a_rteDensities,
					   const RealVect         a_E,
					   const RealVect         a_pos,
					   const Real             a_dx,
					   const Real             a_dt,
					   const Real             a_time,
					   const Real             a_kappa) const {
  //  MayDay::Warning("CdrPlasmaJSON::advanceReactionnetwork -- don't know how to do this yet");

  for (int i = 0; i < a_cdrSources.size(); i++){
    a_cdrSources[i] = 0.0;
  }

  for (int i = 0; i < a_rteSources.size(); i++){
    a_rteSources[i] = 0.0;
  }

  return;
}

Vector<RealVect> CdrPlasmaJSON::computeCdrDriftVelocities(const Real         a_time,
							  const RealVect     a_position,
							  const RealVect     a_E,
							  const Vector<Real> a_cdrDensities) const {
  Vector<RealVect> velocities(m_numCdrSpecies, RealVect::Zero);

  for (int i = 0; i < a_cdrDensities.size(); i++){
    const bool isMobile = m_CdrSpecies[i]->isMobile();
    const int  Z        = m_CdrSpecies[i]->getChargeNumber();
    
    if(isMobile && Z != 0){
      // Figure out the mobility.
      const LookupMethod& method = m_mobilityLookup.at(i);

      Real mobility = 0.0;
      
      switch(method) {
      case LookupMethod::Constant:
	{
	  mobility = m_mobilityConstant.at(i);
	  break;
	}
      case LookupMethod::FunctionLFA:
	{
	  
	  const Real E = a_E.vectorLength();	  
	  const Real N = m_gasDensity(a_position);

	  mobility = m_mobilityFunctionsEN.at(i)(E, N);

	  break;
	}
      case LookupMethod::TableLFA:
	{
	  // Recall; the mobility tables are stored as (E/N, mu*N) so we need to extract mu from that. 
	  const LookupTable<2>& mobilityTable = m_mobilityTables.at(i);

	  const Real E = a_E.vectorLength();
	  const Real N = m_gasDensity(a_position);

	  mobility  = mobilityTable.getEntry<1>(E/N); // Get mu*N
	  mobility /= N;                              // Get mu

	  break;
	}
      default:
	{
	  MayDay::Error("CdrPlasmaJSON::computeCdrDriftVelocities -- logic bust");
	}
      }

      int sgn = 0;
      if (Z > 0){
	sgn = 1;
      }
      else if (Z < 0){
	sgn = -1;
      }

      mobility *= sgn;

      velocities[i] = mobility * a_E;
    }
  }

  return velocities;
}

Vector<Real> CdrPlasmaJSON::computeCdrDiffusionCoefficients(const Real         a_time,
							    const RealVect     a_pos,
							    const RealVect     a_E,
							    const Vector<Real> a_cdrDensities) const {
  //  MayDay::Warning("CdrPlasmaJSON::computeDiffusionCoefficients -- don't know how to do this yet");

  return Vector<Real>(m_numCdrSpecies, 0.0);
}

Vector<Real> CdrPlasmaJSON::computeCdrElectrodeFluxes(const Real         a_time,
						      const RealVect     a_pos,
						      const RealVect     a_normal,
						      const RealVect     a_E,
						      const Vector<Real> a_cdrDensities,
						      const Vector<Real> a_cdrVelocities,
						      const Vector<Real> a_cdrGradients,
						      const Vector<Real> a_rteFluxeses,
						      const Vector<Real> a_extrapCdrFluxes) const {
  //  MayDay::Warning("CdrPlasmaJSON::computeCdrElectrodeFluxes -- don't know how to do this yet");

  return Vector<Real>(m_numCdrSpecies, 0.0);  
}

Vector<Real> CdrPlasmaJSON::computeCdrDielectricFluxes(const Real         a_time,
						       const RealVect     a_pos,
						       const RealVect     a_normal,
						       const RealVect     a_E,
						       const Vector<Real> a_cdrDensities,
						       const Vector<Real> a_cdrVelocities,
						       const Vector<Real> a_cdrGradients,
						       const Vector<Real> a_rteFluxeses,
						       const Vector<Real> a_extrapCdrFluxes) const {
  //  MayDay::Warning("CdrPlasmaJSON::computeCdrDielectricFluxes -- don't know how to do this yet");

  return Vector<Real>(m_numCdrSpecies, 0.0);    
}

Vector<Real> CdrPlasmaJSON::computeCdrDomainFluxes(const Real           a_time,
						   const RealVect       a_pos,
						   const int            a_dir,
						   const Side::LoHiSide a_side,
						   const RealVect       a_E,
						   const Vector<Real>   a_cdrDensities,
						   const Vector<Real>   a_cdrVelocities,
						   const Vector<Real>   a_cdrGradients,
						   const Vector<Real>   a_rteFluxeses,
						   const Vector<Real>   a_extrapCdrFluxes) const {
  //  MayDay::Warning("CdrPlasmaJSON::computeCdrDomainFluxes -- don't know how to do this yet");
  
  return Vector<Real>(m_numCdrSpecies, 0.0);      
}

Real CdrPlasmaJSON::initialSigma(const Real a_time, const RealVect a_pos) const {
  return m_initialSigma(a_time, a_pos);
}


#include <CD_NamespaceFooter.H>
