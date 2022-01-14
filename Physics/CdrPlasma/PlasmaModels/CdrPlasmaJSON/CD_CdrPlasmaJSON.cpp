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
#include <sstream>
#include <algorithm>

// Chombo includes
#include <ParmParse.H>
#include <CH_Timer.H>

// Our includes
#include <CD_CdrPlasmaJSON.H>
#include <CD_DataParser.H>
#include <CD_Random.H>
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
  this->initializePhotonSpecies();  
  this->initializeSigma();

  // Do a sanity check. 
  this->sanityCheckSpecies();

  // Parse CDR mobilities and diffusion coefficients. 
  this->parseMobilities();
  this->parseDiffusion();

  // Parse plasma-reactions and photo-reactions
  this->parsePlasmaReactions();

  // Parse boundary conditions and surface reactions


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

void CdrPlasmaJSON::sanityCheckSpecies() const {
  CH_TIME("CdrPlasmaJSON::sanityCheckSpecies");
  if(m_verbose){
    pout() << "CdrPlasmaJSON::sanityCheckSpecies" << endl;
  }

  std::vector<std::string> allSpecies;

  for (const auto& s : m_neutralSpeciesMap) {
    allSpecies.emplace_back(s.first);
  }

  for (const auto& s : m_cdrSpeciesMap) {
    allSpecies.emplace_back(s.first);
  }

  for (const auto& s : m_rteSpeciesMap) {
    allSpecies.emplace_back(s.first);
  }


  const int numSpecies = allSpecies.size();

  // Sort the container. 
  std::sort(allSpecies.begin(), allSpecies.end());
  for (int i = 0; i < allSpecies.size() - 1; i++){
    if(allSpecies[i] == allSpecies[i+1]){
      this->throwParserError("CdrPlasmaJSON::sanityCheckSpecies -- species '" + allSpecies[i] + "' was defined more than once. Check your neutral, plasma, and photon species");
    }
  }
}

std::string CdrPlasmaJSON::trim(const std::string& a_string) const {

  auto ltrim = [](const std::string a_s) -> std::string {
    std::string s = a_s;
    s.erase(s.begin(), std::find_if(s.begin(), s.end(),
				    std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
  };


  auto rtrim = [](std::string a_s) -> std::string {
    std::string s = a_s;    
    s.erase(std::find_if(s.rbegin(), s.rend(),
			 std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    return s;
  };

  return ltrim(rtrim(a_string));
}

void CdrPlasmaJSON::parseReactionString(std::vector<std::string>& a_reactants,
					std::vector<std::string>& a_products,
					const std::string&        a_reaction) const {
  CH_TIME("CdrPlasmaJSON::parseReactionString");
  if(m_verbose){
    pout() << "CdrPlasmaJSON::parseReactionString()" << endl;
  }

  // Parse the string into segments. We split on white-space. 
  std::stringstream ss(a_reaction);
  std::string segment;
  std::vector<std::string> segments;
    
  while(std::getline(ss, segment, ' ')){
    segments.push_back(segment);
  }

  // Discard all whitespace and solitary + from input string
  segments.erase(std::remove(segments.begin(), segments.end(), "" ), segments.end());
  segments.erase(std::remove(segments.begin(), segments.end(), "+"), segments.end());

  // Find the element containing "->"
  const auto& it = std::find(segments.begin(), segments.end(), "->");  

  // Make sure that -> is in the reaction string.
  if(it == segments.end()) this->throwParserError("CdrPlasmaJSON::parseReactionString -- -" + a_reaction + "' does not contain '->");

  // Left of "->" are reactants and right of "->" are products
  a_reactants = std::vector<std::string>(segments.begin(), it);
  a_products  = std::vector<std::string>(it + 1, segments.end());
}

void CdrPlasmaJSON::initializeNeutralSpecies() {
  CH_TIME("CdrPlasmaJSON::initializeNeutralSpecies");
  if(m_verbose){
    pout() << "CdrPlasmaJSON::initializeNeutralSpecies()" << endl;
  }

  // These fields are required
  if(!(m_json["gas"].contains("temperature"    ))) this->throwParserError("In JSON field 'gas' - field 'temperature' is missing"    );
  if(!(m_json["gas"].contains("pressure"       ))) this->throwParserError("In JSON field 'gas' - field 'pressure' is missing"       );
  if(!(m_json["gas"].contains("law"            ))) this->throwParserError("In JSON field 'gas' - field 'law' is missing"            );
  if(!(m_json["gas"].contains("neutral_species"))) this->throwParserError("In JSON field 'gas' - field 'neutral_species' is missing");

  const auto referenceTemperature =      m_json["gas"]["temperature"].get<Real       >() ;
  const auto referencePressure    =      m_json["gas"]["pressure"   ].get<Real       >() ;
  const auto gasLaw               = trim(m_json["gas"]["law"        ].get<std::string>());

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

  // Initialize the species. This will iterate through the neutral species in the JSON input file. 
  for (const auto& species : m_json["gas"]["neutral_species"]){
    const std::string speciesName     = trim(species["name"          ].get<std::string>());
    const Real        speciesFraction =      species["molar_fraction"].get<Real>() / molarSum;

    // It's an error if a species was defined twice. 
    if(m_neutralSpeciesMap.find(speciesName) != m_neutralSpeciesMap.end()){
      this->throwParserError("Neutral species '" + speciesName + "' was defined more than once");
    }

    // Set the species density function. 
    const std::function<Real(const RealVect)> speciesDensity  = [f = speciesFraction, N = this->m_gasDensity] (const RealVect a_position) {
      return f * N(a_position);
    };

    // Add the species. Make sure the maps are consist.
    const int idx = m_neutralSpecies.size();
    
    // Create the neutral species (and the mapped background density).
    m_neutralSpecies.         push_back(std::shared_ptr<NeutralSpeciesJSON>((new NeutralSpeciesJSON(speciesName, speciesFraction, speciesDensity))));
    m_neutralSpeciesDensities.push_back(speciesDensity);
    
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

  if(!(m_json.contains("plasma_species"))) this->throwParserWarning("CdrPlasmaJSON::initializePlasmaSpecies -- did not find any plasma species");

  // Iterate through all species defined in the JSON file. 
  for (const auto& species : m_json["plasma_species"]){
    if(!(species.contains("name"     ))) this->throwParserError("CdrPlasmaJSON::initializePlasmaSpecies -- 'plasma_species' must have field 'name'"     );
    if(!(species.contains("Z"        ))) this->throwParserError("CdrPlasmaJSON::initializePlasmaSpecies -- 'plasma_species' must have field 'Z'"        );
    if(!(species.contains("mobile"   ))) this->throwParserError("CdrPlasmaJSON::initializePlasmaSpecies -- 'plasma_species' must have field 'mobile'"   );
    if(!(species.contains("diffusive"))) this->throwParserError("CdrPlasmaJSON::initializePlasmaSpecies -- 'plasma_species' must have field 'diffusive'");

    const auto name        = trim(species["name"].     get<std::string>());
    const auto Z           =      species["Z"].        get<int        >() ;
    const auto mobile      =      species["mobile"].   get<bool       >() ;
    const auto diffusive   =      species["diffusive"].get<bool       >() ;


    // It's an error if a species was defined twice. 
    if(m_cdrSpeciesMap.find(name) != m_cdrSpeciesMap.end()){
      this->throwParserError("Plasma species '" + name + "' was defined more than once");
    }    
    
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
      pout() << "CdrPlasmaJSON::initializePlasmaSpecies: instantiating species" << "\n"
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

void CdrPlasmaJSON::initializePhotonSpecies() {
  CH_TIME("CdrPlasmaJSON::initializePhotonSpecies");
  if(m_verbose){
    pout() << "CdrPlasmaJSON::initializePhotonSpecies - file is = " << m_jsonFile << endl;
  }

  if(!(m_json.contains("photon_species"))) this->throwParserWarning("CdrPlasmaJSON::initializePhotonSpecies -- did not find any photon species");

  for (const auto& species : m_json["photon_species"]){
    if(!(species.contains("name" ))) this->throwParserError("CdrPlasmaJSON::initializePhotonSpecies -- 'photon_species' must have field 'name'"     );
    if(!(species.contains("kappa"))) this->throwParserError("CdrPlasmaJSON::initializePhotonSpecies -- 'photon_species' must have field 'mobile'"   );

    const auto name  = trim(species["name" ].get<std::string>());
    const auto kappa = trim(species["kappa"].get<std::string>());

    // Set the kappa-function needed by RteSpeciesJSON.
    std::function<Real(const RealVect a_position)> kappaFunction = [](const RealVect a_position) -> Real {return 1.0;};
    
    if(kappa == "constant"){
      if(!(species.contains("value"))) this->throwParserError("CdrPlasmaJSON::initializePhotonSpecies -- got 'constant' but field 'value' is missing");

      const Real value = species["value"].get<Real>();

      kappaFunction = [value](const RealVect a_position) {return value;};
    }
    else if (kappa == "helmholtz"){
      // For the Bourdon model the translation between the Eddington and Bourdon model is
      //
      //    kappa = pX*lambda/sqrt(3).
      //
      // where pX is the partial pressure for a species. This requires that we are able to compute the partial pressure. Fortunately, we have m_gasPressure
      // which computes the pressure, and m_neutralSpecies also stores the molar fraction for each species so this is comparatively easy to reconstruct.     

      // Make sure that 'lambda' is found in the photon species and that O2 is found in the neutral species. 
      if(!(species.contains("lambda" ))) this->throwParserError("CdrPlasmaJSON::initializePhotonSpecies -- got 'helmholtz' but field 'lambda' is missing" );
      if(!(species.contains("neutral"))) this->throwParserError("CdrPlasmaJSON::initializePhotonSpecies -- got 'helmholtz' but field 'neutral' is missing");

      // Get lambda and the neutral species that we base the partial pressure upon. 
      const auto neutral = trim(species["neutral"].get<std::string>());
      const auto lambda  =      species["lambda" ].get<Real       >() ;      

      // Make sure that the neutral is in the list of species. 
      if(m_neutralSpeciesMap.find(neutral) == m_neutralSpeciesMap.end()) {
	this->throwParserError("CdrPlasmaJSON::initializePhotonSpecies -- got 'helmholtz' but '" + neutral + "' is not in the list of neutral species");
      }

      // Get the molar fraction for this specific neutral
      const Real molarFraction = m_neutralSpecies[m_neutralSpeciesMap.at(neutral)]->getMolarFraction();

      kappaFunction = [lambda, P = this->m_gasPressure, m = molarFraction](const RealVect a_position) -> Real {
	return m*P(a_position)*Units::atm2pascal*lambda/sqrt(3.0);
      };
    }
    else if (kappa == "stochastic_chanrion"){
      // For the Chanrion stochastic model we compute the absorption length as
      //
      //    kappa = k1 * (k2/k1)^((f-f1)/(f2-f1))
      //
      // where k1 = chi_min * p(neutral) and p(neutral) is the partial pressure (in Pa) for some species.
      //
      // Note that the frequency f is sampled stochastically. 
      //
      // This requires that we are able to compute the partial pressure.
      // Fortunately, we have m_gasPressure which computes the pressure, and m_neutralSpecies also stores the molar fraction for each species so this is comparatively
      // easy to reconstruct.
      
      if(!(species.contains("f1"     ))) this->throwParserError("CdrPlasmaJSON::initializePhotonSpecies -- got 'stochastic_chanrion' but field 'f1' is missing"     );
      if(!(species.contains("f2"     ))) this->throwParserError("CdrPlasmaJSON::initializePhotonSpecies -- got 'stochastic_chanrion' but field 'f2' is missing"     );
      if(!(species.contains("chi_min"))) this->throwParserError("CdrPlasmaJSON::initializePhotonSpecies -- got 'stochastic_chanrion' but field 'chi_max' is missing");
      if(!(species.contains("chi_min"))) this->throwParserError("CdrPlasmaJSON::initializePhotonSpecies -- got 'stochastic_chanrion' but field 'chi_min' is missing");
      if(!(species.contains("neutral"))) this->throwParserError("CdrPlasmaJSON::initializePhotonSpecies -- got 'stochastic_chanrion' but field 'neutral' is missing");

      const auto f1      =      species["f1"     ].get<Real       >() ;
      const auto f2      =      species["f1"     ].get<Real       >() ;
      const auto chi_min =      species["chi_min"].get<Real       >() ;
      const auto chi_max =      species["chi_max"].get<Real       >() ;
      const auto neutral = trim(species["neutral"].get<std::string>());

      // Make sure that the neutral is in the list of species. 
      if(m_neutralSpeciesMap.find(neutral) == m_neutralSpeciesMap.end()) {
	this->throwParserError("CdrPlasmaJSON::initializePhotonSpecies -- got 'bourdon' but '" + neutral + "' is not in the list of neutral species");
      }

      // Get the molar fraction for the specified species. 
      const Real molarFraction = m_neutralSpecies[m_neutralSpeciesMap.at(neutral)]->getMolarFraction();

      // Create the absorption function. 
      kappaFunction = [f1, f2, x1=chi_min, x2=chi_max, m=molarFraction, P=this->m_gasPressure](const RealVect a_position) -> Real {

	// Create a uniform distribution on the range [f1,f2]	
	std::uniform_real_distribution<Real> udist = std::uniform_real_distribution<Real>(f1, f2);	

	const Real f = Random::get(udist);
	const Real a = (f-f1)/(f2-f1);

	const Real p  = m * P(a_position);
	const Real K1 = x1*p;
	const Real K2 = x2*p;

	return K1*std::pow(K2/K1, a);
      };
    }
    else{
      this->throwParserError("CdrPlasmaJSON::initializePhotonSpecies -- 'kappa ' is neither 'constant', 'helmholtz', 'stochastic_chanrion'");
    }

    // Initialize the species.
    const int num = m_RtSpecies.size();

    // Make the string-int map encodings. 
    m_rteSpeciesMap.       emplace(std::make_pair(name, num ));
    m_rteSpeciesInverseMap.emplace(std::make_pair(num , name));

    // Push the JSON entry and the new CdrSpecies to corresponding vectors. 
    m_RtSpecies.     push_back(RefCountedPtr<RtSpecies> (new RteSpeciesJSON(name, kappaFunction)));
    m_rteSpeciesJSON.push_back(species);    
  }
}

void CdrPlasmaJSON::initializeSigma() {
  CH_TIME("CdrPlasmaJSON::initializeSigma");
  if(m_verbose){
    pout() << "CdrPlasmaJSON::initializeSigma()" << endl;
  }

  m_initialSigma = [](const RealVect a_pos, const Real a_time) -> Real {
    return 0.0;
  };

  if(m_json.contains("sigma")){
    if(m_json["sigma"].contains("initial_density")){
      const Real sigma = m_json["sigma"]["initial_density"].get<Real>();
      
      m_initialSigma = [sigma] (const RealVect a_position, const Real a_time) -> Real {
	return sigma;
      };
    }
  }
}

void CdrPlasmaJSON::parseMobilities() {
  CH_TIME("CdrPlasmaJSON::parseMobilities");
  if(m_verbose){
    pout() << "CdrPlasmaJSON::parseMobilities - file is = " << m_jsonFile << endl;
  }

  // Iterate through all tracked species.
  for (const auto& species : m_cdrSpeciesJSON){
    const std::string name = trim(species["name"].get<std::string>());
    const int         idx  = m_cdrSpeciesMap.at(name);

    if(m_CdrSpecies[idx]->isMobile()){

      // This is a required field. We use it for specifying the mobility. 
      if(!(species.contains("mobility"))) this->throwParserError("species " + name + " is mobile but file does not contain field 'mobility'");
      const json& mobilityJSON = species["mobility"];

      // Get the mobility lookup method. This must either be a constant, a function, or a table. We parse these cases differently.
      const std::string lookup = trim(mobilityJSON["lookup"].get<std::string>());
	
      if(lookup == "constant"){
	// User specified a constant mobility. We look for a field 'value' in the JSON file and set the mobility from that. If the
	// field does not exist then it's an error. 
	if(!(mobilityJSON.contains("value"))) this->throwParserError("CdrPlasmaJSON::parseMobilities -- constant lookup was specified but I can't find the field 'value'");

	const Real value = mobilityJSON["value"].get<Real>();

	m_mobilityLookup.   emplace(std::make_pair(idx, LookupMethod::Constant));
	m_mobilityConstants.emplace(std::make_pair(idx, value                 ));
      }
      else if(lookup == "table E/N"){
	if(!(mobilityJSON.contains("file"      ))) this->throwParserError("CdrPlasmaJSON::parseMobilities -- tabulated mobility was specified but field 'file' was not"      );
	if(!(mobilityJSON.contains("header"    ))) this->throwParserError("CdrPlasmaJSON::parseMobilities -- tabulated mobility was specified but field 'header' was not"    );
	if(!(mobilityJSON.contains("E/N"       ))) this->throwParserError("CdrPlasmaJSON::parseMobilities -- tabulated mobility was specified but field 'E/N' was not"       );
	if(!(mobilityJSON.contains("mu*N"      ))) this->throwParserError("CdrPlasmaJSON::parseMobilities -- tabulated mobility was specified but field 'mu*N' was not"      );
	if(!(mobilityJSON.contains("min E/N"   ))) this->throwParserError("CdrPlasmaJSON::parseMobilities -- tabulated mobility was specified but field 'min E/N' was not"   );
	if(!(mobilityJSON.contains("max E/N"   ))) this->throwParserError("CdrPlasmaJSON::parseMobilities -- tabulated mobility was specified but field 'max E/N' was not"   );
	if(!(mobilityJSON.contains("num_points"))) this->throwParserError("CdrPlasmaJSON::parseMobilities -- tabulated mobility was specified but field 'num_points' was not");
	
	const std::string filename  = trim(mobilityJSON["file"  ].get<std::string>());
	const std::string startRead = trim(mobilityJSON["header"].get<std::string>());
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

	m_mobilityLookup.  emplace(std::make_pair(idx, LookupMethod::TableLFA));
	m_mobilityTablesEN.emplace(std::make_pair(idx, mobilityTable         ));
      }		
      else if (lookup == "function E/N"){
	if(!(mobilityJSON.contains("function"))) this->throwParserError("CdrPlasmaJSON::parseMobilities -- using function E/N but did not find function specification");

	FunctionEN func;
	
	const std::string whichFunction = trim(mobilityJSON["function"].get<std::string>());
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

void CdrPlasmaJSON::parseDiffusion() {
  CH_TIME("CdrPlasmaJSON::parseDiffusion");
  if(m_verbose){
    pout() << "CdrPlasmaJSON::parseDiffusion - file is = " << m_jsonFile << endl;
  }

  // Iterate through all tracked species.
  for (const auto& species : m_cdrSpeciesJSON){
    const std::string name = trim(species["name"].get<std::string>());
    const int         idx  = m_cdrSpeciesMap.at(name);

    if(m_CdrSpecies[idx]->isDiffusive()){

      // This is a required field. We use it for specifying the mobility. 
      if(!(species.contains("diffusion"))) this->throwParserError("species " + name + " is diffusive but file does not contain field 'diffusion'");
      const json& diffusionJSON = species["diffusion"];

      // Get the mobility lookup method. This must either be a constant, a function, or a table. We parse these cases differently.
      const std::string lookup = trim(diffusionJSON["lookup"].get<std::string>());
	
      if(lookup == "constant"){
	// User specified a constant mobility. We look for a field 'value' in the JSON file and set the mobility from that. If the
	// field does not exist then it's an error. 
	if(!(diffusionJSON.contains("value"))) this->throwParserError("CdrPlasmaJSON::parseDiffusion -- constant lookup was specified but I can't find the field 'value'");

	const Real value = diffusionJSON["value"].get<Real>();

	m_diffusionLookup.  emplace(std::make_pair(idx, LookupMethod::Constant));
	m_diffusionConstants.emplace(std::make_pair(idx, value                 ));
      }
      else if(lookup == "table E/N"){
	if(!(diffusionJSON.contains("file"      ))) this->throwParserError("CdrPlasmaJSON::parseDiffusion -- tabulated diffusion was specified but field 'file' was not."      );
	if(!(diffusionJSON.contains("header"    ))) this->throwParserError("CdrPlasmaJSON::parseDiffusion -- tabulated diffusion was specified but field 'header' was not."    );
	if(!(diffusionJSON.contains("E/N"       ))) this->throwParserError("CdrPlasmaJSON::parseDiffusion -- tabulated diffusion was specified but field 'E/N' was not."       );
	if(!(diffusionJSON.contains("D*N"       ))) this->throwParserError("CdrPlasmaJSON::parseDiffusion -- tabulated diffusion was specified but field 'D*N' was not"        );
	if(!(diffusionJSON.contains("min E/N"   ))) this->throwParserError("CdrPlasmaJSON::parseDiffusion -- tabulated diffusion was specified but field 'min E/N' was not."   );
	if(!(diffusionJSON.contains("max E/N"   ))) this->throwParserError("CdrPlasmaJSON::parseDiffusion -- tabulated diffusion was specified but field 'max E/N' was not."   );
	if(!(diffusionJSON.contains("num_points"))) this->throwParserError("CdrPlasmaJSON::parseDiffusion -- tabulated diffusion was specified but field 'num_points' was not.");
	
	const std::string filename  = trim(diffusionJSON["file"  ].get<std::string>());
	const std::string startRead = trim(diffusionJSON["header"].get<std::string>());
	const std::string stopRead  = "";

	const int xColumn   = diffusionJSON["E/N"       ].get<int>();
	const int yColumn   = diffusionJSON["D*N"       ].get<int>();
	const int numPoints = diffusionJSON["num_points"].get<int>();

	const Real minEN = diffusionJSON["min E/N"].get<Real>();
	const Real maxEN = diffusionJSON["max E/N"].get<Real>();	

	// Read the table and format it. We happen to know that this function reads data into the approprate columns. So if
	// the user specified the correct E/N column then that data will be put in the first column. The data for D*N will be in the
	// second column. 
	LookupTable<2> diffusionTable = DataParser::fractionalFileReadASCII(filename, startRead, stopRead, xColumn, yColumn);

	// Format the table
	diffusionTable.setRange(minEN, maxEN, 0);
	diffusionTable.sort(0);
	diffusionTable.makeUniform(numPoints);

	m_diffusionLookup.  emplace(std::make_pair(idx, LookupMethod::TableLFA));
	m_diffusionTablesEN.emplace(std::make_pair(idx, diffusionTable         ));
      }		
      else if (lookup == "function E/N"){
	if(!(diffusionJSON.contains("function"))) this->throwParserError("CdrPlasmaJSON::parseDiffusion -- using function E/N but did not find function specification");

	FunctionEN func;
	
	const std::string whichFunction = trim(diffusionJSON["function"].get<std::string>());
	if(whichFunction == "ABC"){
	  if(!(diffusionJSON.contains("A"))) this->throwParserError("CdrPlasmaJSON::parseDiffusion -- using function 'ABC'  did not find 'A'");
	  if(!(diffusionJSON.contains("B"))) this->throwParserError("CdrPlasmaJSON::parseDiffusion -- using function 'ABC'  did not find 'B'");
	  if(!(diffusionJSON.contains("C"))) this->throwParserError("CdrPlasmaJSON::parseDiffusion -- using function 'ABC'  did not find 'C'");

	  const Real A = diffusionJSON["A"].get<Real>();
	  const Real B = diffusionJSON["B"].get<Real>();
	  const Real C = diffusionJSON["C"].get<Real>();

	  func = [A, B, C](const Real a_E, const Real a_N) -> Real {return A * std::pow(a_E, B) / std::pow(a_N, C);};
	}
	else{
	  this->throwParserError("CdrPlasmaJSON::parseDiffusion -- using function E/N but do not recognize specification '" + whichFunction + "'");
	}

	m_diffusionLookup.     emplace(std::make_pair(idx, LookupMethod::FunctionLFA));
	m_diffusionFunctionsEN.emplace(std::make_pair(idx, func                     ));	
      }
      else{
	this->throwParserError("CdrPlasmaJSON::parseDiffusion -- lookup = '" + lookup + "' was specified but this is not 'constant', 'function', or 'table'");
      }
    }
  }
}

void CdrPlasmaJSON::parsePlasmaReactions() {
  CH_TIME("CdrPlasmaJSON::parsePlasmaReactions");
  if(m_verbose){
    pout() << "CdrPlasmaJSON::parsePlasmaReactions - file is = " << m_jsonFile << endl;
  }

  if(!(m_json.contains("plasma_reactions"))) this->throwParserWarning("CdrPlasmaJSON::parsePlasmaReactions -- did not find any reactions");  

  for (const auto& R : m_json["plasma_reactions"]){
    if(!(R.contains("reaction"))) this->throwParserError("CdrPlasmaJSON::parsePlasmaReactions -- field 'reaction' is missing");
    if(!(R.contains("rate"    ))) this->throwParserError("CdrPlasmaJSON::parsePlasmaReactions -- field 'rate' is missing");

    const std::string reaction = trim(R["reaction"].get<std::string>());
    const std::string lookup   = trim(R["lookup"  ].get<std::string>());

    // Parse the reaction string to figure out the species involved in the reaction.
    std::vector<std::string> reactants;
    std::vector<std::string> products ;

    this->parseReactionString   (reactants, products, reaction);
    this->sanctifyPlasmaReaction(reactants, products, reaction);

    // Make the string-int encoding so we can encode the reaction properly.
    std::list<int> plasmaReactants ;
    std::list<int> neutralReactants;
    std::list<int> plasmaProducts  ;
    std::list<int> photonProducts  ;

    this->getPlasmaReactionProducts(plasmaReactants,
				    neutralReactants,
				    plasmaProducts,
				    photonProducts,
				    reactants,
				    products);


    //
    const int reactionIndex = m_plasmaReactions.size();


    // Start monkeying with the rate.
    //    LookupMethod 
    if(lookup == "constant"){

      // Constant reaction rates are easy, just fetch it and put it where it belongs. 
      if(!(R.contains("rate"))) this->throwParserError("CdrPlasmaJSON::parsePlasmaReactions -- got 'constant' but did not get 'rate'");

      const Real k = R["rate"].get<Real>();

      /// Add the rate and lookup method. 
      m_plasmaReactionLookup.   emplace(std::make_pair(reactionIndex, LookupMethod::Constant));
      m_plasmaReactionConstants.emplace(std::make_pair(reactionIndex, k                     ));
    }
    else if(lookup == "function E/N"){
      MayDay::Error("CdrPlasmaJSON::parsePlasmaReaction -- 'function E/N' not supported yet");
    }
    else if (lookup == "table E/N"){
      MayDay::Error("CdrPlasmaJSON::parsePlasmaReaction -- 'table E/N' not supported yet");      
    }
    else{
      this->throwParserError("CdrPlasmaJSON::parsePlasmaReactions -- lookup = '" + lookup + "' was specified but this is not 'constant', 'function', or 'table'");      
    }


    // Add the reaction to the list of reactions.
    m_plasmaReactions.emplace_back(plasmaReactants, neutralReactants, plasmaProducts, photonProducts);    
  }
}

void CdrPlasmaJSON::sanctifyPlasmaReaction(const std::vector<std::string>& a_reactants,
					   const std::vector<std::string>& a_products,
					   const std::string               a_reaction) const {
  CH_TIME("CdrPlasmaJSON::sanctifyPlasmaReaction");
  if(m_verbose){
    pout() << "CdrPlasmaJSON::sanctifyPlasmaReaction" << m_jsonFile << endl;
  }

    // All reactants must be in the list of neutral species or in the list of plasma species
    for (const auto& r : a_reactants){
      if(m_cdrSpeciesMap.    find(r) == m_cdrSpeciesMap.    end() &&
	 m_neutralSpeciesMap.find(r) == m_neutralSpeciesMap.end()) {
	this->throwParserError("CdrPlasmaJSON::sanctifyPlasmaReaction -- I do not know reacting species '" + r + "' for reaction '" + a_reaction + "'.");
      }
    }

    // All products should be in the list of plasma or photon species. It's ok if users include a neutral species -- we will ignore it (but tell the user about it).
    for (const auto& p : a_products){
      if(m_cdrSpeciesMap.find(p) == m_cdrSpeciesMap.end() && m_rteSpeciesMap.find(p) == m_rteSpeciesMap.end()) {
	if(m_neutralSpeciesMap.find(p) == m_neutralSpeciesMap.end()) {
	  this->throwParserError("CdrPlasmaJSON::parsePlasmaReactions -- I do not know product species '" + p + "' for reaction '" + a_reaction + "'.");
	}
	else{
	  this->throwParserWarning("CdrPlasmaJSON::parsePlasmaReactions -- neutral species '" + p + "' was specified in reaction but will be ignored");
	}
      }
    }

    // Check for charge conservation
    int sumCharge;
    for (const auto& r : a_reactants){
      if (m_cdrSpeciesMap.find(r) != m_cdrSpeciesMap.end()){
	sumCharge -= m_CdrSpecies[m_cdrSpeciesMap.at(r)]->getChargeNumber();
      }
    }
    for (const auto& p : a_products){
      if (m_cdrSpeciesMap.find(p) != m_cdrSpeciesMap.end()){
	sumCharge += m_CdrSpecies[m_cdrSpeciesMap.at(p)]->getChargeNumber();
      }
    }

    if(sumCharge != 0) this->throwParserError("CdrPlasmaJSON::parsePlasmaReaction -- charge not conserved for reaction '" + a_reaction + "'.");  
}

void CdrPlasmaJSON::getPlasmaReactionProducts(std::list<int>&                 a_plasmaReactants,
					      std::list<int>&                 a_neutralReactants,
					      std::list<int>&                 a_plasmaProducts,
					      std::list<int>&                 a_photonProducts,
					      const std::vector<std::string>& a_reactants,
					      const std::vector<std::string>& a_products) const {

  CH_TIME("CdrPlasmaJSON::getPlasmaReactionProducts");
  if(m_verbose){
    pout() << "CdrPlasmaJSON::getPlasmaReactionProducts" << endl;
  }

  a_plasmaReactants. clear();
  a_neutralReactants.clear();
  a_plasmaProducts.  clear();
  a_photonProducts.  clear();  

  // Figure out the reactants. 
  for (const auto& r : a_reactants){

    // Figure out if the reactants is a plasma species or a neutral species.
    const bool isPlasmaSpecies  = m_cdrSpeciesMap.    find(r) != m_cdrSpeciesMap.    end();
    const bool isNeutralSpecies = m_neutralSpeciesMap.find(r) != m_neutralSpeciesMap.end();

    if(isPlasmaSpecies && !isNeutralSpecies){
      a_plasmaReactants.emplace_back(m_cdrSpeciesMap.at(r));
    }
    else if (!isPlasmaSpecies && isNeutralSpecies){
      a_neutralReactants.emplace_back(m_neutralSpeciesMap.at(r));      
    }
    else {
      this->throwParserError("CdrPlasmaJSON::getPlasmaReactionProducts - logic bust 1");
    }
  }

  // Figure out the products.
  for (const auto& p : a_products){
    const bool isPlasmaSpecies  = m_cdrSpeciesMap.    find(p) != m_cdrSpeciesMap.    end();
    const bool isNeutralSpecies = m_neutralSpeciesMap.find(p) != m_neutralSpeciesMap.end();
    const bool isPhotonSpecies  = m_rteSpeciesMap.    find(p) != m_rteSpeciesMap.    end();

    if(isPlasmaSpecies && !isPhotonSpecies && !isNeutralSpecies){
      a_plasmaProducts.emplace_back(m_cdrSpeciesMap.at(p));
    }
    else if(!isPlasmaSpecies && isPhotonSpecies && !isNeutralSpecies){
      a_photonProducts.emplace_back(m_rteSpeciesMap.at(p));
    }
    else if(!isPlasmaSpecies && !isPhotonSpecies && isNeutralSpecies){
      // do nothing
    }
    else{
      this->throwParserError("CdrPlasmaJSON::getPlasmaReactionProducts - logic bust 2");
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

  // I really hate Chombo sometimes. 
  std::vector<Real>& cdrSources = a_cdrSources.stdVector();
  std::vector<Real>& rteSources = a_rteSources.stdVector();

  const std::vector<Real    >& cdrDensities = ((Vector<Real    >&) a_cdrDensities).stdVector();
  const std::vector<RealVect>& cdrGradients = ((Vector<RealVect>&) a_cdrGradients).stdVector();
  const std::vector<Real    >& rteDensities = ((Vector<Real    >&) a_rteDensities).stdVector();


  // Set all sources to zero. 
  for (auto& S : cdrSources){
    S = 0.0;
  }

  for (auto& S : rteSources){
    S = 0.0;
  }

  // Here is the lambda for firing a plasma reaction that has a rate constant k. 
  auto FirePlasmaReaction = [&position = a_pos,
			     &CdrSrc   = cdrSources,
			     &RteSrc   = rteSources,
			     &CdrPhi   = cdrDensities,
			     &Neutrals = this->m_neutralSpeciesDensities] (const Real a_rate,
									   const CdrPlasmaReactionJSON& a_reaction) -> void {
    Real volumetricRate = a_rate;

    // These are the various species involved. 
    const std::list<int>& plasmaReactants  = a_reaction.getPlasmaReactants ();    
    const std::list<int>& neutralReactants = a_reaction.getNeutralReactants();
    const std::list<int>& plasmaProducts   = a_reaction.getPlasmaProducts  ();    
    const std::list<int>& photonProducts   = a_reaction.getPhotonProducts  ();

    // Compute the total consumption. 
    for (const auto& n : neutralReactants){
      volumetricRate *= (Neutrals[n])(position);
    }

    for (const auto& r : plasmaReactants){
      volumetricRate *= CdrPhi[r];
    }

    // Remove consumption on the left-hand side.
    for (const auto& r : plasmaReactants){
      CdrSrc[r] -= volumetricRate;
    }

    // Add mass on the right-hand side.
    for (const auto& p : plasmaProducts){
      CdrSrc[p] += volumetricRate;
    }

    for (const auto& p : plasmaProducts){
      RteSrc[p] += volumetricRate;
    }
  };

  // Here is the lambda that fires a photon reaction that has a rate constant k. 


  // Iterate through plasma reactions
  for (int i = 0; i < m_plasmaReactions.size(); i++){
    
    // Figure out the reaction rate for this reaction. 
    Real k = 0.0;

    const LookupMethod& method = m_plasmaReactionLookup.at(i);

    switch(method){
      case LookupMethod::Constant:
	{
	  k = m_plasmaReactionConstants.at(i);

	  break;
	}
      case LookupMethod::FunctionLFA:
	{
	  MayDay::Error("CdrPlasmaJSON::advanceReactionNetwork -- FunctionLFA not supported yet");
	  
	  break;
	}
      case LookupMethod::TableLFA:
	{
	  MayDay::Error("CdrPlasmaJSON::advanceReactionNetwork -- Table not supported yet");
	  
	  break;
	}
    default:
      {
      MayDay::Error("CdrPlasmaJSON::advanceReactionNetwork -- logic bust");
      break;
      }
    }

    // Fire the reaction.
    FirePlasmaReaction(k, m_plasmaReactions[i]);
  }


  return;
}

Vector<RealVect> CdrPlasmaJSON::computeCdrDriftVelocities(const Real         a_time,
							  const RealVect     a_position,
							  const RealVect     a_E,
							  const Vector<Real> a_cdrDensities) const {
  Vector<RealVect> velocities(m_numCdrSpecies, RealVect::Zero);

  const Real E   = a_E.vectorLength();
  const Real N   = m_gasDensity(a_position);
  const Real Etd = (E/(N * Units::Td));  

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
	  mobility = m_mobilityConstants.at(i);
	  break;
	}
      case LookupMethod::FunctionLFA:
	{
	  mobility = m_mobilityFunctionsEN.at(i)(E, N);
	  break;
	}
      case LookupMethod::TableLFA:
	{
	  // Recall; the mobility tables are stored as (E/N, mu*N) so we need to extract mu from that. 
	  const LookupTable<2>& mobilityTable = m_mobilityTablesEN.at(i);

	  mobility  = mobilityTable.getEntry<1>(Etd); // Get mu*N
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
							    const RealVect     a_position,
							    const RealVect     a_E,
							    const Vector<Real> a_cdrDensities) const {
  Vector<Real> diffusionCoefficients(m_numCdrSpecies, 0.0);

  const Real E   = a_E.vectorLength();
  const Real N   = m_gasDensity(a_position);
  const Real Etd = (E/(N * Units::Td));    

  for (int i = 0; i < a_cdrDensities.size(); i++){
    if(m_CdrSpecies[i]->isDiffusive()){
      
      // Figure out how we compute the diffusion coefficient for this species. 
      const LookupMethod& method = m_diffusionLookup.at(i);

      Real Dco = 0.0;
      
      switch(method) {
      case LookupMethod::Constant:
	{
	  Dco = m_diffusionConstants.at(i);
	  break;
	}
      case LookupMethod::FunctionLFA:
	{
	  const Real E = a_E.vectorLength();	  
	  const Real N = m_gasDensity(a_position);

	  Dco = m_diffusionFunctionsEN.at(i)(E, N);

	  break;
	}
      case LookupMethod::TableLFA:
	{
	  // Recall; the diffusion tables are stored as (E/N, D*N) so we need to extract D from that. 
	  const LookupTable<2>& diffusionTable = m_diffusionTablesEN.at(i);

	  Dco  = diffusionTable.getEntry<1>(Etd); // Get D*N
	  Dco /= N;                               // Get D

	  break;
	}
      default:
	{
	  MayDay::Error("CdrPlasmaJSON::computeCdrDiffusionCoefficients -- logic bust");
	}
      }

      diffusionCoefficients[i] = Dco;
    }
  }

  return diffusionCoefficients;  
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
  return m_initialSigma(a_pos, a_time);
}


#include <CD_NamespaceFooter.H>
