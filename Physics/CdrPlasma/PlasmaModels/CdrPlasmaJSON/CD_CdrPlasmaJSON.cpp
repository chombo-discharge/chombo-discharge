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
  this->parseTemperatures();

  // Parse plasma-reactions and photo-reactions
  this->parsePlasmaReactions();
  this->parsePhotoReactions();

  m_numCdrSpecies = m_CdrSpecies.size();
  m_numRtSpecies  = m_RtSpecies. size();
}

CdrPlasmaJSON::~CdrPlasmaJSON(){
  CH_TIME("CdrPlasmaJSON::~CdrPlasmaJSON()");
}

void CdrPlasmaJSON::parseOptions() {
  CH_TIME("CdrPlasmaJSON::parseOptions");
  
  ParmParse pp("CdrPlasmaJSON");

  pp.get("verbose",          m_verbose );
  pp.get("chemistry_file",   m_jsonFile);
  pp.get("discrete_photons", m_discretePhotons);
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

bool CdrPlasmaJSON::containsWildcard(const std::string a_str) const {
  return (a_str.find("@") != std::string::npos);
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
  if(!(m_json["gas"].contains("neutral species"))) this->throwParserError("In JSON field 'gas' - field 'neutral species' is missing");

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
    if(!(m_json["gas"].contains("molar mass"))) this->throwParserError("gas_law is 'troposphere' but I did not find field 'molar mass'");
    if(!(m_json["gas"].contains("gravity"   ))) this->throwParserError("gas_law is 'troposphere' but I did not find field 'gravity'"   );
    if(!(m_json["gas"].contains("lapse rate"))) this->throwParserError("gas_law is 'troposphere' but I did not find field 'lapse rate'");

    const Real g    = m_json["gas"]["gravity"   ].get<Real>();
    const Real L    = m_json["gas"]["lapse rate"].get<Real>();
    const Real M    = m_json["gas"]["molar mass"].get<Real>();
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
    m_gasDensity = [&P = this->m_gasPressure, &T = this->m_gasTemperature] (const RealVect a_position) -> Real {
      return (P(a_position) * Units::atm2pascal * Units::Na)/ (T(a_position) * Units::R);          
    };
  }
  else{
    this->throwParserError("CdrPlasmaJSON::initializeNeutralSpecies gas law '" + gasLaw + "' not recognized.");
  }

  // Instantiate the species densities. Note that we need to go through this twice because we need to normalize the molar fractions in case users
  // were a bit inconsiderate when setting them. 
  Real molarSum = 0.0;
  for (const auto& species : m_json["gas"]["neutral species"]){
    if(!(species.contains("name")))           this->throwParserError("In JSON field 'neutral species' - field 'name' is required"          );
    if(!(species.contains("molar fraction"))) this->throwParserError("In JSON field 'neutral species' - field 'molar fraction' is required");

    molarSum += species["molar fraction"].get<Real>();
  }

  // Initialize the species. This will iterate through the neutral species in the JSON input file. 
  for (const auto& species : m_json["gas"]["neutral species"]){
    const std::string speciesName     = trim(species["name"          ].get<std::string>());
    const Real        speciesFraction =      species["molar fraction"].get<Real>() / molarSum;

    // Names can not contain the at letter.
    if(containsWildcard(speciesName)) this->throwParserError("CdrPlasmaJSON::initializeNeutralSpecies -- species name must not contain '@' letter");

    // It's an error if a species was defined twice.
    if(isNeutralSpecies(speciesName)) this->throwParserError("Neutral species '" + speciesName + "' was defined more than once");        

    // Set the species density function. 
    const std::function<Real(const RealVect)> speciesDensity  = [f = speciesFraction, &N = this->m_gasDensity] (const RealVect a_position) {
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

  if(!(m_json.contains("plasma species"))) this->throwParserWarning("CdrPlasmaJSON::initializePlasmaSpecies -- did not find any plasma species");

  // Iterate through all species defined in the JSON file. 
  for (const auto& species : m_json["plasma species"]){
    if(!(species.contains("name"     ))) this->throwParserError("CdrPlasmaJSON::initializePlasmaSpecies -- 'plasma species' must have field 'name'"     );
    if(!(species.contains("Z"        ))) this->throwParserError("CdrPlasmaJSON::initializePlasmaSpecies -- 'plasma species' must have field 'Z'"        );
    if(!(species.contains("mobile"   ))) this->throwParserError("CdrPlasmaJSON::initializePlasmaSpecies -- 'plasma species' must have field 'mobile'"   );
    if(!(species.contains("diffusive"))) this->throwParserError("CdrPlasmaJSON::initializePlasmaSpecies -- 'plasma species' must have field 'diffusive'");

    const auto name        = trim(species["name"].     get<std::string>());
    const auto Z           =      species["Z"].        get<int        >() ;
    const auto mobile      =      species["mobile"].   get<bool       >() ;
    const auto diffusive   =      species["diffusive"].get<bool       >() ;

    // Does not get to contain at letter
    if(containsWildcard(name)) this->throwParserError("CdrPlasmaJSON::initializePlasmaSpecies -- species name must not contain '@' letter");    

    // It's an error if the species was already defined. 
    if(isPlasmaSpecies(name)) this->throwParserError("Plasma species '" + name + "' was defined more than once");
    
    const bool hasInitData = species.contains("initial data");

    // Get the initial data. 
    std::function<Real(const RealVect, const Real)> initFunc;

    if(hasInitData){
      initFunc = [data = species["initial data"]] (const RealVect a_point, const Real a_time) -> Real {
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

  if(!(m_json.contains("photon species"))) this->throwParserWarning("CdrPlasmaJSON::initializePhotonSpecies -- did not find any photon species");

  for (const auto& species : m_json["photon species"]){
    if(!(species.contains("name" ))) this->throwParserError("CdrPlasmaJSON::initializePhotonSpecies -- 'photon species' must have field 'name'"     );
    if(!(species.contains("kappa"))) this->throwParserError("CdrPlasmaJSON::initializePhotonSpecies -- 'photon species' must have field 'mobile'"   );

    const auto name  = trim(species["name" ].get<std::string>());
    const auto kappa = trim(species["kappa"].get<std::string>());

    // Does not get to contain at letter.
    if(containsWildcard(name)) this->throwParserError("CdrPlasmaJSON::initializePhotonSpecies -- species name cannot contain '@'");

    // It's an error if the species is already defined.
    if(isPhotonSpecies(name)) this->throwParserError("Photon species '" + name + "' was defined more than once");    

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
    if(m_json["sigma"].contains("initial density")){
      const Real sigma = m_json["sigma"]["initial density"].get<Real>();
      
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
	if(!(mobilityJSON.contains("num points"))) this->throwParserError("CdrPlasmaJSON::parseMobilities -- tabulated mobility was specified but field 'num points' was not");
	
	const std::string filename  = trim(mobilityJSON["file"  ].get<std::string>());
	const std::string startRead = trim(mobilityJSON["header"].get<std::string>());
	const std::string stopRead  = "";

	const int xColumn   = mobilityJSON["E/N"       ].get<int>();
	const int yColumn   = mobilityJSON["mu*N"      ].get<int>();
	const int numPoints = mobilityJSON["num points"].get<int>();

	const Real minEN = mobilityJSON["min E/N"].get<Real>();
	const Real maxEN = mobilityJSON["max E/N"].get<Real>();	

	// Read the table and format it. We happen to know that this function reads data into the approprate columns. So if
	// the user specified the correct E/N column then that data will be put in the first column. The data for mu*N will be in the
	// second column. 
	LookupTable<2> mobilityTable = DataParser::fractionalFileReadASCII(filename, startRead, stopRead, xColumn, yColumn);

      // If the table is empty then it's an error.
	if(mobilityTable.getNumEntries() == 0){
	  this->throwParserError("Mobility table '" + startRead + "' in file '" + filename + "'is empty. This is probably an error");
	}	

	// Format the table
	mobilityTable.setRange(minEN, maxEN, 0);
	mobilityTable.sort(0);
	mobilityTable.makeUniform(numPoints);

	m_mobilityLookup.  emplace(std::make_pair(idx, LookupMethod::TableEN));
	m_mobilityTablesEN.emplace(std::make_pair(idx, mobilityTable         ));
      }		
      else if (lookup == "functionEN_A"){
	FunctionEN func;
	
	if(!(mobilityJSON.contains("c1"))) this->throwParserError("CdrPlasmaJSON::parseMobilities -- using 'functionEN_A' did not find 'c1'");
	if(!(mobilityJSON.contains("c2"))) this->throwParserError("CdrPlasmaJSON::parseMobilities -- using 'functionEN_A' did not find 'c2'");
	if(!(mobilityJSON.contains("c3"))) this->throwParserError("CdrPlasmaJSON::parseMobilities -- using 'functionEN_A' did not find 'c3'");

	const Real A = mobilityJSON["c1"].get<Real>();
	const Real B = mobilityJSON["c2"].get<Real>();
	const Real C = mobilityJSON["c3"].get<Real>();

	func = [A, B, C](const Real a_E, const Real a_N) -> Real {return A * std::pow(a_E, B) / std::pow(a_N, C);};

	m_mobilityLookup.     emplace(std::make_pair(idx, LookupMethod::FunctionEN));
	m_mobilityFunctionsEN.emplace(std::make_pair(idx, func                     ));	
      }
      else{
	this->throwParserError("CdrPlasmaJSON::parseMobilities - logic bust");
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
	if(!(diffusionJSON.contains("file"      ))) this->throwParserError("CdrPlasmaJSON::parseDiffusion - tabulated diffusion was specified but field 'file' was not."      );
	if(!(diffusionJSON.contains("header"    ))) this->throwParserError("CdrPlasmaJSON::parseDiffusion - tabulated diffusion was specified but field 'header' was not."    );
	if(!(diffusionJSON.contains("E/N"       ))) this->throwParserError("CdrPlasmaJSON::parseDiffusion - tabulated diffusion was specified but field 'E/N' was not."       );
	if(!(diffusionJSON.contains("D*N"       ))) this->throwParserError("CdrPlasmaJSON::parseDiffusion - tabulated diffusion was specified but field 'D*N' was not"        );
	if(!(diffusionJSON.contains("min E/N"   ))) this->throwParserError("CdrPlasmaJSON::parseDiffusion - tabulated diffusion was specified but field 'min E/N' was not."   );
	if(!(diffusionJSON.contains("max E/N"   ))) this->throwParserError("CdrPlasmaJSON::parseDiffusion - tabulated diffusion was specified but field 'max E/N' was not."   );
	if(!(diffusionJSON.contains("num points"))) this->throwParserError("CdrPlasmaJSON::parseDiffusion - tabulated diffusion was specified but field 'num points' was not.");
	
	const std::string filename  = trim(diffusionJSON["file"  ].get<std::string>());
	const std::string startRead = trim(diffusionJSON["header"].get<std::string>());
	const std::string stopRead  = "";

	const int xColumn   = diffusionJSON["E/N"       ].get<int>();
	const int yColumn   = diffusionJSON["D*N"       ].get<int>();
	const int numPoints = diffusionJSON["num points"].get<int>();

	const Real minEN = diffusionJSON["min E/N"].get<Real>();
	const Real maxEN = diffusionJSON["max E/N"].get<Real>();	

	// Read the table and format it. We happen to know that this function reads data into the approprate columns. So if
	// the user specified the correct E/N column then that data will be put in the first column. The data for D*N will be in the
	// second column. 
	LookupTable<2> diffusionTable = DataParser::fractionalFileReadASCII(filename, startRead, stopRead, xColumn, yColumn);

	// If the table is empty then it's an error.
	if(diffusionTable.getNumEntries() == 0){
	  this->throwParserError("Diffusion table '" + startRead + "' in file '" + filename + "'is empty. This is probably an error");	  
	}		

	// Format the table
	diffusionTable.setRange(minEN, maxEN, 0);
	diffusionTable.sort(0);
	diffusionTable.makeUniform(numPoints);

	m_diffusionLookup.  emplace(std::make_pair(idx, LookupMethod::TableEN));
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

	m_diffusionLookup.     emplace(std::make_pair(idx, LookupMethod::FunctionEN));
	m_diffusionFunctionsEN.emplace(std::make_pair(idx, func                     ));	
      }
      else{
	this->throwParserError("CdrPlasmaJSON::parseDiffusion -- lookup = '" + lookup + "' was specified but this is not 'constant', 'function', or 'table'");
      }
    }
  }
}

void CdrPlasmaJSON::parseTemperatures() {
  CH_TIME("CdrPlasmaJSON::parseTemperatures");
  if(m_verbose){
    pout() << "CdrPlasmaJSON::parseTemperatures - file is = " << m_jsonFile << endl;
  }

  for (const auto& species : m_cdrSpeciesJSON){
    const std::string name = trim(species["name"].get<std::string>());
    const int         idx  = m_cdrSpeciesMap.at(name);

    // If temperature is not specified -- we use the background gas temperature. 
    if(!(species.contains("temperature"))){
      m_temperatureLookup.   emplace(idx, LookupMethod::FunctionX);
      m_temperatureConstants.emplace(idx, m_gasTemperature       );
    }
    else{
      const json& S = species["temperature"];

      // We MUST have a lookup field in order to determine how we compute the temperature for a species. 
      if(!(S.contains("lookup"))) this->throwParserError("CdrPlasmaJSON::parseTemperatures -- temperature specified but field 'lookup' is missing");


      // Figure out the lookup method. 
      const std::string lookup = trim(S["lookup"].get<std::string>());
      if(lookup == "constant"){
	if(!(S.contains("value"))) this->throwParserError("CdrPlasmaJSON::parseTemperatures -- 'constant' specified but field 'value' is missing");

	const Real value = S["value"].get<Real>();

	// Create a function which returns a constant value everywhere. 
	m_temperatureLookup.   emplace(idx, LookupMethod::FunctionX);
	m_temperatureConstants.emplace(idx, [value] (const RealVect a_psition) -> Real {return value;});
      }
      else if (lookup == "table E/N"){
	if(!(S.contains("file"      ))) this->throwParserError("CdrPlasmaJSON::parseTemperatures - tabulated diffusion was specified but field 'file' was not."      );
	if(!(S.contains("header"    ))) this->throwParserError("CdrPlasmaJSON::parseTemperatures - tabulated diffusion was specified but field 'header' was not."    );
	if(!(S.contains("E/N"       ))) this->throwParserError("CdrPlasmaJSON::parseTemperatures - tabulated diffusion was specified but field 'E/N' was not."       );
	if(!(S.contains("eV"        ))) this->throwParserError("CdrPlasmaJSON::parseTemperatures - tabulated diffusion was specified but field 'eV' was not"         );
	if(!(S.contains("min E/N"   ))) this->throwParserError("CdrPlasmaJSON::parseTemperatures - tabulated diffusion was specified but field 'min E/N' was not."   );
	if(!(S.contains("max E/N"   ))) this->throwParserError("CdrPlasmaJSON::parseTemperatures - tabulated diffusion was specified but field 'max E/N' was not."   );
	if(!(S.contains("num points"))) this->throwParserError("CdrPlasmaJSON::parseTemperatures - tabulated diffusion was specified but field 'num points' was not.");

	const std::string filename  = trim(S["file"  ].get<std::string>());
	const std::string startRead = trim(S["header"].get<std::string>());
	const std::string stopRead  = "";

	const int xColumn   = S["E/N"       ].get<int>();
	const int yColumn   = S["eV"        ].get<int>();
	const int numPoints = S["num points"].get<int>();

	const Real minEN = S["min E/N"].get<Real>();
	const Real maxEN = S["max E/N"].get<Real>();	

	// Read the table and format it. We happen to know that this function reads data into the approprate columns. So if
	// the user specified the correct E/N column then that data will be put in the first column. The data for D*N will be in the
	// second column. 
	LookupTable<2> temperatureTable = DataParser::fractionalFileReadASCII(filename, startRead, stopRead, xColumn, yColumn);

	// If the table is empty then it's an error.
	if(temperatureTable.getNumEntries() == 0){
	  this->throwParserError("Temperature table '" + startRead + "' in file '" + filename + "'is empty. This is probably an error");	  
	}		

	// Format the table
	temperatureTable.setRange(minEN, maxEN, 0);
	temperatureTable.sort(0);
	temperatureTable.makeUniform(numPoints);

	// Conversion factor is eV to Kelvin.
	temperatureTable.scale<1>( (2.0*Units::Qe) / (3.0*Units::kb) );

	m_temperatureLookup.  emplace(std::make_pair(idx, LookupMethod::TableEN));
	m_temperatureTablesEN.emplace(std::make_pair(idx, temperatureTable      ));	
      }
      else{
	this->throwParserError("CdrPlasmaJSON::parseTemperatures -- I do not know the field '" + lookup + "'. Supported are 'constant' and 'table E/N'");
      }
    }
  }
}

void CdrPlasmaJSON::parsePlasmaReactions() {
  CH_TIME("CdrPlasmaJSON::parsePlasmaReactions");
  if(m_verbose){
    pout() << "CdrPlasmaJSON::parsePlasmaReactions - file is = " << m_jsonFile << endl;
  }

  for (const auto& R : m_json["plasma reactions"]){
    if(!(R.contains("reaction"))) this->throwParserError("CdrPlasmaJSON::parsePlasmaReactions -- field 'reaction' is missing from one of the reactions");
    if(!(R.contains("lookup"  ))) this->throwParserError("CdrPlasmaJSON::parsePlasmaReactions -- field 'lookup' is missing from one of the reactions");

    const std::string reaction  = trim(R["reaction"].get<std::string>());
    const std::string baseError = "CdrPlasmaJSON::parsePlasmaReactions for reaction '" + reaction + "' ";    

    // Parse the reaction string to figure out the species involved in the reaction. Note that this CAN involve the species
    // wildcard @, in which case we need need to build a superset of reaction strings that we parse. We do that below. 
    std::vector<std::string> reactants;
    std::vector<std::string> products ;
    
    this->parseReactionString(reactants, products, reaction);

    // Parse wildcards and make the reaction sets. std::list<std::pair<std::vector<std::string>, std::vector<std::string> > >
    const auto reactionSets = this->parsePlasmaReactionWildcards(reactants, products, R);

    // Now run through the superset of reactions encoded by the input string. Because of the @ character some of the rate functions
    // need special handling. 
    for (const auto& curReaction : reactionSets){
      const std::vector<std::string> curReactants = curReaction.first ;
      const std::vector<std::string> curProducts  = curReaction.second;

      // This is the reaction index for the current index. The reaction we are currently
      // dealing with is put in m_plasmaReactions[reactionIdex]. 
      const int reactionIndex = m_plasmaReactions.size();      

      // Make sure reaction string makes sense. 
      this->sanctifyPlasmaReaction(curReactants, curProducts, reaction);      

      // Parse the reaction parameters. 
      this->parsePlasmaReactionRate       (reactionIndex, R);
      this->parsePlasmaReactionScaling    (reactionIndex, R);
      this->parsePlasmaReactionPlot       (reactionIndex, R);
      this->parsePlasmaReactionDescription(reactionIndex, R);
      this->parsePlasmaReactionSoloviev   (reactionIndex, R);

      // Make the string-int encoding so we can encode the reaction properly. Then add the reaction to the pile. 
      std::list<int> plasmaReactants ;
      std::list<int> neutralReactants;
      std::list<int> plasmaProducts  ;
      std::list<int> photonProducts  ;

      this->getPlasmaReactionProducts(plasmaReactants,
				      neutralReactants,
				      plasmaProducts,
				      photonProducts,
				      curReactants,
				      curProducts);      

      m_plasmaReactions.emplace_back(plasmaReactants, neutralReactants, plasmaProducts, photonProducts);
    }
  }
}

std::list<std::pair<std::vector<std::string>, std::vector<std::string> > > CdrPlasmaJSON::parsePlasmaReactionWildcards(const std::vector<std::string>& a_reactants,
														       const std::vector<std::string>& a_products,
														       const json& a_R){
  CH_TIME("CdrPlasmaJSON::parsePlasmaReactionWildcards");
  if(m_verbose){
    pout() << "CdrPlasmaJSON::parsePlasmaReactionWildcards" << endl;
  }

  // This is what we return. A horrific creature.
  std::list<std::pair<std::vector<std::string>, std::vector<std::string> > > reactionSets;  

  // This is the reaction name. 
  const std::string reaction  = a_R["reaction"].get<std::string>();
  const std::string baseError = "CdrPlasmaJSON::parsePlasmaReactionWildcards for reaction '" + reaction + "'";

  // Check if reaction string had a wildcard '@'. If it did we replace the wildcard with the corresponding species. This means that we need to
  // build additional reactions. 
  const bool containsWildcard = this->containsWildcard(reaction);
  
  if(containsWildcard){
    if(!(a_R.contains("@"))) {
      this->throwParserError(baseError + "got reaction wildcard '@' but array '@:' was specified");
    }

    // Get the wildcards array. 
    const std::vector<std::string> wildcards = a_R["@"].get<std::vector<std::string> >();

    // Go through the wildcards and replace appropriately. 
    for (const auto& w : wildcards){
      std::vector<std::string> curReactants;
      std::vector<std::string> curProducts;

      // Replace by wildcard in reactants.  	
      for (const auto& r : a_reactants){
	if(this->containsWildcard(r)){
	  curReactants.emplace_back(w);
	}
	else{
	  curReactants.emplace_back(r);
	}
      }

      // Replace by wildcard in reactants.  	
      for (const auto& p : a_products){
	if(this->containsWildcard(p)){
	  curProducts.emplace_back(w);
	}
	else{
	  curProducts.emplace_back(p);
	}
      }

      reactionSets.emplace_back(curReactants, curProducts);
    }
  }
  else{
    reactionSets.emplace_back(a_reactants, a_products);
  }

  return reactionSets;
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
    if(!isPlasmaSpecies(r) && !isNeutralSpecies(r))
      this->throwParserError("CdrPlasmaJSON::sanctifyPlasmaReaction -- I do not know reacting species '" + r + "' for reaction '" + a_reaction + "'");
  }

  // All products should be in the list of plasma or photon species. It's ok if users include a neutral species -- we will ignore it (but tell the user about it).
  for (const auto& p : a_products){
    if(!isPlasmaSpecies(p) && !isPhotonSpecies(p) && !isNeutralSpecies(p)){
      this->throwParserError("CdrPlasmaJSON::parsePlasmaReactions -- I do not know product species '" + p + "' for reaction '" + a_reaction + "'.");
    }
  }

  // Check for charge conservation
  int sumCharge = 0;
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

  if(sumCharge != 0) {
    this->throwParserError("CdrPlasmaJSON::sanctifyPlasmaReaction -- charge not conserved for reaction '" + a_reaction + "'.");
  }
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

void CdrPlasmaJSON::parsePlasmaReactionRate(const int a_reactionIndex, const json& a_R){
  CH_TIME("CdrPlasmaJSON::parsePlasmaReactionRate");
  if(m_verbose){
    pout() << "CdrPlasmaJSON::parsePlasmaReactionRate" << endl;
  }

  // These MUST exist because they are checked for in parsePlasmaReactions().
  const std::string reaction  = trim(a_R["reaction"].get<std::string>());
  const std::string lookup    = trim(a_R["lookup"  ].get<std::string>());
  const std::string baseError = "CdrPlasmaJSON::parsePlasmaReactionRate for reaction '" + reaction + "' ";
      
  // Figure out how we compute the reaction rate for this reaction. 
  if(lookup == "constant"){

    // Constant reaction rates are easy, just fetch it and put it where it belongs. 
    if(!(a_R.contains("rate"))) this->throwParserError(baseError + " and got 'constant' but did not get 'rate'");

    const Real k = a_R["rate"].get<Real>();

    /// Add the rate and lookup method. 
    m_plasmaReactionLookup.   emplace(std::make_pair(a_reactionIndex, LookupMethod::Constant));
    m_plasmaReactionConstants.emplace(std::make_pair(a_reactionIndex, k                     ));
  }
  else if(lookup == "function T1T2_A"){
    if(!(a_R.contains("T1"))) this->throwParserError(baseError + "and got 'function T1T2_A' but field 'T1' was not found");
    if(!(a_R.contains("T2"))) this->throwParserError(baseError + "and got 'function T1T2_A' but field 'T2' was not found");
    if(!(a_R.contains("c1"))) this->throwParserError(baseError + "and got 'function T1T2_A' but field 'c1' was not found");
    if(!(a_R.contains("c2"))) this->throwParserError(baseError + "and got 'function T1T2_A' but field 'c2' was not found");    

    const std::string speciesT1 = trim(a_R["T1"].get<std::string>());
    const std::string speciesT2 = trim(a_R["T2"].get<std::string>());

    const bool isPlasmaT1  = this->isPlasmaSpecies(speciesT1);
    const bool isPlasmaT2  = this->isPlasmaSpecies(speciesT2);

    const bool isNeutralT1 = this->isNeutralSpecies(speciesT1);
    const bool isNeutralT2 = this->isNeutralSpecies(speciesT2);      

    // Make sure that the specified species exist. 
    if(!isPlasmaT1 && !isNeutralT1) this->throwParserError(baseError + "and got function 'function T1T2_A' but do not know species '" + speciesT1 + "'");
    if(!isPlasmaT2 && !isNeutralT2) this->throwParserError(baseError + "and got function 'function T1T2_A' but do not know species '" + speciesT2 + "'");    

    // This syntax may look weird, but we need to know precisely which temperatures are involved in the reaction. In general the reaction rate
    // is just a function k = f(T1, T2) but T1 and T2 could be the temperatures for either a plasma or a neutral species. So, the plasmaReactionFunctionsT1T2 map
    // is a very weird one because only the temperatures are passed into the functions and not the species themselves. This also happens to be the correct
    // design because the temperatures themselves can be computed in many formats, and we just don't have a way of reconstructing them here. 
    int firstIndex  = -1;
    int secondIndex = -1;

    if(isPlasmaT1) firstIndex  = m_cdrSpeciesMap.at(speciesT1);
    if(isPlasmaT2) secondIndex = m_cdrSpeciesMap.at(speciesT2);

    const Real c1 = a_R["c1"].get<Real>();
    const Real c2 = a_R["c2"].get<Real>();

    FunctionTT functionT1T2 = [=](const Real a_T1, const Real a_T2) -> Real {
      return c1*std::pow(a_T1/a_T2, c2);
    };

    m_plasmaReactionFunctionsTT.emplace(a_reactionIndex, std::make_tuple(firstIndex, secondIndex, functionT1T2));
    m_plasmaReactionLookup.     emplace(a_reactionIndex, LookupMethod::FunctionTT);
  }
  else if (lookup == "table E/N"){
    if(!(a_R.contains("file"      ))) this->throwParserError(baseError + "and got 'table E/N' but field 'file' was not found"      );
    if(!(a_R.contains("header"    ))) this->throwParserError(baseError + "and got 'table E/N' but field 'header' was not found"    );
    if(!(a_R.contains("E/N"       ))) this->throwParserError(baseError + "and got 'table E/N' but field 'E/N' was not found"       );
    if(!(a_R.contains("rate/N"    ))) this->throwParserError(baseError + "and got 'table E/N' but field 'rate/N' was not found"    );
    if(!(a_R.contains("min E/N"   ))) this->throwParserError(baseError + "and got 'table E/N' but field 'min E/N' was not found"   );
    if(!(a_R.contains("max E/N"   ))) this->throwParserError(baseError + "and got 'table E/N' but field 'max E/N' was not found"   );
    if(!(a_R.contains("num points"))) this->throwParserError(baseError + "and got 'table E/N' but field 'num points' was not found");

    const std::string filename  = trim(a_R["file"  ].get<std::string>());
    const std::string startRead = trim(a_R["header"].get<std::string>());
    const std::string stopRead  = "";

    const int xColumn   = a_R["E/N"       ].get<int>();
    const int yColumn   = a_R["rate/N"    ].get<int>();
    const int numPoints = a_R["num points"].get<int>();

    const Real minEN = a_R["min E/N"].get<Real>();
    const Real maxEN = a_R["max E/N"].get<Real>();	

    // Read the table and format it. We happen to know that this function reads data into the approprate columns. So if
    // the user specified the correct E/N column then that data will be put in the first column. The data for D*N will be in the
    // second column. 
    LookupTable<2> reactionTable = DataParser::fractionalFileReadASCII(filename, startRead, stopRead, xColumn, yColumn);

    // If the table is empty then it's an error.
    if(reactionTable.getNumEntries() == 0){
      this->throwParserError(baseError + "and got 'table E/N' but table is empty. This is probably an error");
    }

    // Format the table. 
    reactionTable.setRange(minEN, maxEN, 0);
    reactionTable.sort(0);
    reactionTable.makeUniform(numPoints);

    // Add the tabulated rate and identifier. 
    m_plasmaReactionLookup.  emplace(std::make_pair(a_reactionIndex, LookupMethod::TableEN));
    m_plasmaReactionTablesEN.emplace(std::make_pair(a_reactionIndex, reactionTable         ));      
  }
  else{
    this->throwParserError(baseError + "but lookup = '" + lookup + "' is not recognized");
  }  
}

void CdrPlasmaJSON::parsePlasmaReactionScaling(const int a_index, const json& a_R) {
  CH_TIME("CdrPlasmaJSON::parsePlasmaReactionScaling");
  if(m_verbose){
    pout() << "CdrPlasmaJSON::parsePlasmaReactionScaling" << endl;
  }

  // TLDR: This routine tries to figure out how we should scale reactions. In general, reactions are scaled parametrically by a function f = f(E,x)
  //       where E is the electric field magnitude and x is the physical position. This is completely general, but there are many, many ways that
  //       such functions could be formed. Currently, we put these functions as one of the following:
  //
  //       f = s
  //         s -> scaling
  // 
  //       f = s * pq/(p(x) + pq)
  //         s  -> scaling
  //         pq -> quenching pressure (in atm)
  //      
  //
  //       f = s * kr/(kr + kp + kq(N)) * nu * mu
  //          s  -> scaling
  //          kr -> radiation rate
  //          kp -> predissociation rate
  //          kq -> quenching rate
  //          nu -> photoionization efficiency
  //          mu -> excitation efficiency efficiency  

  // This is the reaction name. We need it for debugging.
  const std::string reaction  = a_R["reaction"].get<std::string>();
  const std::string baseError = "CdrPlasmaJSON::parsePlasmaReactionScaling for reaction '" + reaction + "' ";

  // Scaling factor for reaction. 
  Real scale = 1.0;

  // Quenching pressure in case that is specified. 
  Real pq    = 0.0;

  // Stuff for when "photoionization" is specified. 
  Real kr        = 0.0;
  Real kp        = 0.0;
  Real kqN       = 0.0;
  Real photoiEff = 0.0;
  Real exciteEff = 0.0;      

  // Flags for determining how to scale the reaction. 
  bool doPressureQuenching = false;
  bool doPhotoIonization   = false;

  // Get the scaling factor.
  if(a_R.contains("scale")) {
    scale = a_R["scale"].get<Real>();
  }      

  if(a_R.contains("quenching pressure")){
    doPressureQuenching = true;
    
    pq = a_R["quenching pressure"].get<Real>();
  }

  // Parse the photoionization field. 
  if(a_R.contains("photoionization")){
    const json& photoi = a_R["photoionization"];
	
    // These fields are required.
    if(!(photoi.contains("kr"  )))       this->throwParserError(baseError + "got 'photoionization but field 'kr' is missing"        );
    if(!(photoi.contains("kp"  )))       this->throwParserError(baseError + "got 'photoionization but field 'kp' is missing"        );
    if(!(photoi.contains("kq/N")))       this->throwParserError(baseError + "got 'photoionization but field 'kq/N' is missing"      );
    if(!(photoi.contains("photoi eff"))) this->throwParserError(baseError + "got 'photoionization but field 'photoi eff' is missing");
    if(!(photoi.contains("excite eff"))) this->throwParserError(baseError + "got 'photoionization but field 'excite eff' is missing");
	
    kr        = photoi["kr"        ].get<Real>();
    kp        = photoi["kp"        ].get<Real>();
    kqN       = photoi["kq/N"      ].get<Real>();
    photoiEff = photoi["photoi eff"].get<Real>();
    exciteEff = photoi["excite eff"].get<Real>();	

    doPhotoIonization = true;
  }

  // This is the scaling function. It will be stored in m_plasmaReactionEfficiencies. 
  FunctionEX func;

  // Now make ourselves a lambda that we can use for scaling the reactions. 
  if(doPhotoIonization && doPressureQuenching) {
    this->throwParserError(baseError + "- cannot specify both 'photoionization' and 'quenching pressure'");
  }
  else if(doPressureQuenching && !doPhotoIonization){
    func = [scale, pq, p = this->m_gasPressure](const Real E, const RealVect x){
      return scale * pq/(pq + p(x));
    };
  }
  else if(!doPressureQuenching && doPhotoIonization){
    func = [scale, kr, kp, kqN, photoiEff, exciteEff, &N = this->m_gasDensity](const Real E, const RealVect x) -> Real {
      const Real kq = kqN * N(x);

      return scale * kr/(kr + kp + kq) * photoiEff * exciteEff;
    };
  }
  else {
    func = [scale](const Real E, const RealVect x) -> Real {
      return scale;
    };
  }

  // Add it to the pile. 
  m_plasmaReactionEfficiencies.emplace(a_index, func);
}

void CdrPlasmaJSON::parsePlasmaReactionPlot(const int a_reactionIndex, const json& a_R) {
  CH_TIME("CdrPlasmaJSON::parsePlasmaReactionPlot");
  if(m_verbose){
    pout() << "CdrPlasmaJSON::parsePlasmaReactionPlot" << endl;
  }

  bool plot = false;

  if(a_R.contains("plot")){
    plot = a_R.get<bool>();
  }

  m_plasmaReactionPlot.emplace(a_reactionIndex, plot);  
}

void CdrPlasmaJSON::parsePlasmaReactionDescription(const int a_reactionIndex, const json& a_R) {
  CH_TIME("CdrPlasmaJSON::parsePlasmaReactionDescription");
  if(m_verbose){
    pout() << "CdrPlasmaJSON::parsePlasmaReactionDescription" << endl;
  }

  // Reaction name. 
  std::string str = a_R["reaction"].get<std::string>();  

  // Determine if the reaction had a field "description". If it did, we will use that description in I/O files
  if(a_R.contains("description")){
    str = a_R["description"].get<std::string>();
  }

  m_plasmaReactionDescriptions.emplace(a_reactionIndex, str);        
}

void CdrPlasmaJSON::parsePlasmaReactionSoloviev(const int a_reactionIndex, const json& a_R) {
  CH_TIME("CdrPlasmaJSON::parsePlasmaReactionScaling");
  if(m_verbose){
    pout() << "CdrPlasmaJSON::parsePlasmaReactionScaling" << endl;
  }

  // TLDR: This will look through reactions and check if we should use the Soloviev energy correction for LFA-based models. 
  if(a_R.contains("soloviev")){
    
    const std::string reaction  = a_R["reaction"].get<std::string>();
    const std::string baseError = "CdrPlasmaJSON::parsePlasmaReactionSoloviev for '" + reaction + "' ";

    const json& solo = a_R["soloviev"];

    if(!solo.contains("correction")) this->throwParserError(baseError + "but did not find field 'correction'");
    if(!solo.contains("species"   )) this->throwParserError(baseError + "but did not find field 'species'"   );

    const auto correct = solo["correction"].get<bool       >();
    const auto species = solo["species"   ].get<std::string>();

    if(correct){
      if(!isPlasmaSpecies(species)) this->throwParserError(baseError + "but '" + species + "' is not a plasma species");

      const int plasmaSpecies = m_cdrSpeciesMap.at(species);

      const bool isMobile    = m_CdrSpecies[plasmaSpecies]->isMobile   ();
      const bool isDiffusive = m_CdrSpecies[plasmaSpecies]->isDiffusive();

      if(!isMobile   ) this->throwParserError(baseError + "but species  + '" + species + "' isn't mobile."   );
      if(!isDiffusive) this->throwParserError(baseError + "but species  + '" + species + "' isn't diffusive.");

      m_plasmaReactionSolovievCorrection.emplace(a_reactionIndex, std::make_pair(true, plasmaSpecies));
    }
  }
  else{
    m_plasmaReactionSolovievCorrection.emplace(a_reactionIndex, std::make_pair(false, -1));
  }
}


void CdrPlasmaJSON::sanctifyPhotoReaction(const std::vector<std::string>& a_reactants,
					  const std::vector<std::string>& a_products,
					  const std::string               a_reaction) const {
  CH_TIME("CdrPlasmaJSON::sanctifyPhotoReaction");
  if(m_verbose){
    pout() << "CdrPlasmaJSON::sanctifyPhotoReaction" << m_jsonFile << endl;
  }

  const std::string baseError = "CdrPlasmaJSON::sanctifyPhotoReaction for reaction '" + a_reaction + "' ";

  // All reactants must be in the list of neutral species or in the list of photon species
  int numPhotonSpecies = 0;  
  for (const auto& r : a_reactants){
    const bool isPlasma  = this->isPlasmaSpecies (r);
    const bool isNeutral = this->isNeutralSpecies(r);
    const bool isPhoton  = this->isPhotonSpecies (r);

    if( isNeutral) this->throwParserError(baseError + "neutral species (" + r + ") not allowed on left-hand side");
    if( isPlasma ) this->throwParserError(baseError + "plasma species (" + r + ") not allowed on left-hand side");
    if(!isPhoton ) this->throwParserError(baseError + "I do not know species species '" + r + "' on left hand side");
    if( isPhoton ) numPhotonSpecies++;
  }

  // There can only be one photon species on the left-hand side of the reaction. 
  if(numPhotonSpecies != 1){
    this->throwParserError("CdrPlasmaJSON::sanctifyPhotoReaction -- only one photon species allowed on left-hand side of photo-reaction '" + a_reaction + "'");
  }

  // All products should be in the list of plasma, neutral, or photon species. 
  for (const auto& p : a_products){
    const bool isPlasma  = this->isPlasmaSpecies (p);
    const bool isNeutral = this->isNeutralSpecies(p);
    const bool isPhoton  = this->isPhotonSpecies (p);

    if( isPhoton              ) this->throwParserError(baseError + "photon species '" + p + "' not allowed on right hand side"); 
    if(!isPlasma && !isNeutral) this->throwParserError(baseError + "I do not know species '" + p + "' on right hand side"     );
  }

  // Check for charge conservation
  int sumCharge = 0;
  for (const auto& r : a_reactants){
    if(this->isPlasmaSpecies(r)) sumCharge -= m_CdrSpecies[m_cdrSpeciesMap.at(r)]->getChargeNumber();
  }
  for (const auto& p : a_products){
    if(this->isPlasmaSpecies(p)) sumCharge += m_CdrSpecies[m_cdrSpeciesMap.at(p)]->getChargeNumber();
  }

  if(sumCharge != 0) {
    this->throwParserError(baseError + "charge is not conserved!");
  }
}

void CdrPlasmaJSON::parsePhotoReactions(){
  CH_TIME("CdrPlasmaJSON::parsePhotoReactions");
  if(m_verbose){
    pout() << "CdrPlasmaJSON::parsePhotoReactions file is = " << m_jsonFile << endl;
  }

  for (const auto& R : m_json["photo reactions"]){
    if(!(R.contains("reaction"))) this->throwParserError("CdrPlasmaJSON::parsePhotoReactions -- field 'reaction' is missing from one of the reactions");

    const std::string reaction = trim(R["reaction"].get<std::string>());

    std::vector<std::string> reactants;
    std::vector<std::string> products;

    // Index.
    const int reactionIndex = m_photoReactions.size();

    this->parseReactionString  (reactants, products, reaction);
    this->sanctifyPhotoReaction(reactants, products, reaction);


    // Parse if we use Helmholtz reconstruction for this photoionization method.
    this->parsePhotoReactionScaling(reactionIndex, R);

    // Make the string-int encoding so we can encode the reaction properly. Then add the reaction to the pile. 
    std::list<int> plasmaReactants ;
    std::list<int> neutralReactants;
    std::list<int> photonReactants ;
    std::list<int> plasmaProducts  ;
    std::list<int> neutralProducts ;

    this->getPhotoReactionProducts(plasmaReactants,
				   neutralReactants,
				   photonReactants,
				   plasmaProducts,
				   neutralProducts,
				   reactants,
				   products);

    // Add it to the pile.
    m_photoReactions.emplace_back(plasmaReactants, neutralReactants, photonReactants, plasmaProducts, neutralProducts);
  }
}

void CdrPlasmaJSON::parsePhotoReactionScaling(const int a_reactionIndex, const json& a_R) {
  CH_TIME("CdrPlasmaJSON::parsePhotoReactionScale");
  if(m_verbose){
    pout() << "CdrPlasmaJSON::parsePhotoReactionScale" << endl;
  }

  // TLDR: This routine tries to figure out how we should scale photo-reactions. In general, reactions are scaled parametrically by a function f = f(E,x)
  //       where E is the electric field magnitude and x is the physical position. This is completely general, but there are many, many ways that
  //       such functions could be formed. Currently, we put these functions as one of the following:
  //
  //       f = s
  //         s -> scaling
  // 
  //       f = s * c * A * pX/(sqrt(3) * l)
  //         s  -> scaling
  //         c  -> Speed of light
  //         A  -> Helmholtz constants
  //         pX -> Partial pressure
  //         pq -> quenching pressure (in atm)

  // This is the reaction name. We need it for debugging.
  const std::string reaction  = a_R["reaction"].get<std::string>();
  const std::string baseError = "CdrPlasmaJSON::parsePhotoReactionScaling for reaction '" + reaction + "' ";

  // Regular scale
  Real scale = 1.0;

  // Parameters for the Helmholtz reconstruction. 
  Real helmholtzFactor = 0.0;
  bool doHelmholtz = false;

  if(a_R.contains("scale")){
    scale = a_R["scale"].get<Real>();
  }

  if(a_R.contains("helmholtz")){
    const json& helm = a_R["helmholtz"];

    if(!(helm.contains("A")))       this->throwParserError(baseError + " and got Helmholtz but field 'A' is missing"      );
    if(!(helm.contains("lambda")))  this->throwParserError(baseError + " and got Helmholtz but field 'lambda' is missing" );
    if(!(helm.contains("species"))) this->throwParserError(baseError + " and got Helmholtz but field 'neutral' is missing");

    const auto A       = helm["A"      ].get<Real       >();
    const auto lambda  = helm["lambda" ].get<Real       >();
    const auto neutral = helm["species"].get<std::string>();
    
    // Get the specified species and it's molar fraction. 
    // Make sure that it's a neutral species.
    if(!(this->isNeutralSpecies(neutral))) {
      this->throwParserError(baseError + " and got Helmholtz but species '" + neutral + "' is not a neutral species");
    }

    // Make the Helmholtz factor
    const Real molarFraction = m_neutralSpecies[m_neutralSpeciesMap.at(neutral)]->getMolarFraction();

    helmholtzFactor = Units::c * A * molarFraction * Units::atm2pascal/(sqrt(3.0) * lambda);

    doHelmholtz = true;
  }


  // Make the function object that we scale with. 
  FunctionEX func;

  if(doHelmholtz){
    func = [scale, helmholtzFactor, &p = this->m_gasPressure](const Real E, const RealVect x){
      return scale * helmholtzFactor * p(x);
    };
  }
  else{
    func = [scale] (const Real E, const RealVect x){
      return scale;
    };
  }

  m_photoReactionEfficiencies.emplace(a_reactionIndex, func       );
  m_photoReactionUseHelmholtz.emplace(a_reactionIndex, doHelmholtz);
}

void CdrPlasmaJSON::getPhotoReactionProducts(std::list<int>&                 a_plasmaReactants,
					     std::list<int>&                 a_neutralReactants,
					     std::list<int>&                 a_photonReactants,					    
					     std::list<int>&                 a_plasmaProducts,
					     std::list<int>&                 a_neutralProducts,
					     const std::vector<std::string>& a_reactants,
					     const std::vector<std::string>& a_products) const {
  CH_TIME("CdrPlasmaJSON::getPhotoReactionProducts");
  if(m_verbose){
    pout() << "CdrPlasmaJSON::getPhotoReactionProducts" << endl;
  }

  a_plasmaReactants. clear();
  a_neutralReactants.clear();
  a_photonReactants. clear();    
  a_plasmaProducts.  clear();
  a_neutralProducts. clear();  

  // Figure out the reactants. 
  for (const auto& r : a_reactants){

    // Figure out if the reactants is a plasma species or a neutral species.
    const bool isPlasma  = this->isPlasmaSpecies (r);
    const bool isNeutral = this->isNeutralSpecies(r);
    const bool isPhoton  = this->isPhotonSpecies (r);

    if(isPlasma && !isNeutral && !isPhoton){
      a_plasmaReactants.emplace_back(m_cdrSpeciesMap.at(r));
    }
    else if(!isPlasma && isNeutral && !isPhoton){
      a_neutralReactants.emplace_back(m_neutralSpeciesMap.at(r));      
    }
    else if(!isPlasma && !isNeutral && isPhoton){
      a_photonReactants.emplace_back(m_rteSpeciesMap.at(r));      
    }    
    else {
      this->throwParserError("CdrPlasmaJSON::getPhotoReactionProducts - logic bust 1");
    }
  }

  // Figure out the reactants. 
  for (const auto& p : a_products){

    // Figure out if the reactants is a plasma species or a neutral species.
    const bool isPlasma  = this->isPlasmaSpecies (p);
    const bool isNeutral = this->isNeutralSpecies(p);
    const bool isPhoton  = this->isPhotonSpecies (p);

    if(isPlasma && !isNeutral && !isPhoton){
      a_plasmaProducts.emplace_back(m_cdrSpeciesMap.at(p));
    }
    else if(!isPlasma && isNeutral && !isPhoton){
      a_neutralProducts.emplace_back(m_neutralSpeciesMap.at(p));      
    }
    else {
      this->throwParserError("CdrPlasmaJSON::getPhotoReactionProducts - logic bust 1");
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

bool CdrPlasmaJSON::isNeutralSpecies(const std::string& a_name) const {
  bool found = false;

  if(m_neutralSpeciesMap.find(a_name) != m_neutralSpeciesMap.end()){
    found = true;
  }
  else {
    found = false;
  }

  return found;
}

bool CdrPlasmaJSON::isPlasmaSpecies(const std::string& a_name) const {
  bool found = false;

  if(m_cdrSpeciesMap.find(a_name) != m_cdrSpeciesMap.end()){
    found = true;
  }
  else {
    found = false;
  }

  return found;
}

bool CdrPlasmaJSON::isPhotonSpecies(const std::string& a_name) const {
  bool found = false;

  if(m_rteSpeciesMap.find(a_name) != m_rteSpeciesMap.end()){
    found = true;
  }
  else {
    found = false;
  }

  return found;
}

std::vector<Real> CdrPlasmaJSON::computePlasmaSpeciesMobilities(const RealVect&          a_position,
								const RealVect&          a_E,
								const std::vector<Real>& a_cdrDensities) const {

  // Get E/N .
  const Real E   = a_E.vectorLength();
  const Real N   = m_gasDensity(a_position);
  const Real Etd = (E/(N * Units::Td));

  // vector of mobilities
  std::vector<Real> mu(m_numCdrSpecies, 0.0);  

  // Go through each species. 
  for (int i = 0; i < a_cdrDensities.size(); i++){
    const bool isMobile = m_CdrSpecies[i]->isMobile();
    const int  Z        = m_CdrSpecies[i]->getChargeNumber();

    // Figure out how to compute the moiblity. 
    if(isMobile && Z != 0){

      const LookupMethod& method = m_mobilityLookup.at(i);
      
      switch(method) {
      case LookupMethod::Constant:
	{
	  mu[i] = m_mobilityConstants.at(i);
	  break;
	}
      case LookupMethod::FunctionEN:
	{
	  mu[i] = m_mobilityFunctionsEN.at(i)(E, N);
	  break;
	}
      case LookupMethod::TableEN:
	{
	  // Recall; the mobility tables are stored as (E/N, mu*N) so we need to extract mu from that. 
	  const LookupTable<2>& mobilityTable = m_mobilityTablesEN.at(i);

	  mu[i]  = mobilityTable.getEntry<1>(Etd); // Get mu*N
	  mu[i] /= N;                              // Get mu

	  break;
	}
      default:
	{
	  MayDay::Error("CdrPlasmaJSON::computePlasmaSpeciesMobilities -- logic bust when computing the mobility. ");
	}
      }
    }
  }

  return mu;
}

std::vector<Real> CdrPlasmaJSON::computePlasmaSpeciesTemperatures(const RealVect&          a_position,
								  const RealVect&          a_E,
								  const std::vector<Real>& a_cdrDensities) const {

  // Electric field and neutral density. 
  const Real N   = m_gasDensity(a_position);
  const Real E   = a_E.vectorLength();
  const Real Etd = (E/(N * Units::Td));

  // Return vector of temperatures. 
  std::vector<Real> T(m_numCdrSpecies, 0.0);
  
  for (int i = 0; i < m_numCdrSpecies; i++){
    const LookupMethod lookup = m_temperatureLookup.at(i);

    // Switch between various lookup methods. 
    switch(lookup) {
    case LookupMethod::FunctionX:
      {
	T[i] = (m_temperatureConstants.at(i))(a_position);
	
	break;
      }
    case LookupMethod::TableEN:
      {
	// Recall; the temperature tables are stored as (E/N, K) so we can fetch the temperature immediately. 
	const LookupTable<2>& temperatureTable = m_temperatureTablesEN.at(i);

	T[i] = temperatureTable.getEntry<1>(Etd);

	break;
      }
    default:
      {
	MayDay::Error("CdrPlasmaJSON::computePlasmaSpeciesTemperatures -- logic bust when computing species temperature");
	
	break;
      }
    }
  }

  return T;
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
  // I really hate Chombo sometimes. 
  std::vector<Real>& cdrSources = a_cdrSources.stdVector();
  std::vector<Real>& rteSources = a_rteSources.stdVector();

  const std::vector<Real    >& cdrDensities = ((Vector<Real    >&) a_cdrDensities).stdVector();
  const std::vector<RealVect>& cdrGradients = ((Vector<RealVect>&) a_cdrGradients).stdVector();
  const std::vector<Real    >& rteDensities = ((Vector<Real    >&) a_rteDensities).stdVector();

  // These may or may not be needed.
  const std::vector<Real    >  cdrMobilities            = this->computePlasmaSpeciesMobilities  (        a_pos, a_E,   cdrDensities)            ;  
  const std::vector<Real    >& cdrDiffusionCoefficients = this->computeCdrDiffusionCoefficients (a_time, a_pos, a_E, a_cdrDensities).stdVector();
  const std::vector<Real    >  cdrTemperatures          = this->computePlasmaSpeciesTemperatures(        a_pos, a_E,   cdrDensities)            ;

  // Electric field and reduce electric field. 
  const Real E   = a_E.vectorLength();
  const Real N   = m_gasDensity(a_pos);
  const Real Etd = (E/(N * Units::Td));      

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

    for (const auto& p : photonProducts){
      RteSrc[p] += volumetricRate;
    }
  };

  // Here is the lambda that fires a photon reaction that has a rate constant k. Note that this is different
  // if we use discrete photons!


  // Iterate through plasma reactions
  for (int i = 0; i < m_plasmaReactions.size(); i++) {
    
    // Figure out the reaction rate for this reaction. 
    Real k = 0.0;

    const LookupMethod& method = m_plasmaReactionLookup.at(i);

    switch(method){
    case LookupMethod::Constant:
      {
	k = m_plasmaReactionConstants.at(i);

	break;
      }
    case LookupMethod::FunctionEN:
      {
	k = m_plasmaReactionFunctionsEN.at(i)(E, N);
	  
	break;
      }
    case LookupMethod::FunctionTT:
      {
	const std::tuple<int, int, FunctionTT>& tup = m_plasmaReactionFunctionsTT.at(i);

	const int         idx1 = std::get<0>(tup);
	const int         idx2 = std::get<1>(tup);
	const FunctionTT& func = std::get<2>(tup);

	const Real T1 = (idx1 < 0) ? m_gasTemperature(a_pos) : cdrTemperatures[idx1];
	const Real T2 = (idx2 < 0) ? m_gasTemperature(a_pos) : cdrTemperatures[idx2];
	
	k = func(T1, T2);

	break;
      }
    case LookupMethod::TableEN:
      {
	// Recall; the reaction tables are stored as (E/N, rate/N) so we need to extract mu from that. 
	const LookupTable<2>& reactionTable = m_plasmaReactionTablesEN.at(i);

	k  = reactionTable.getEntry<1>(Etd); // Get rate/N
	k *= N;                              // Get rate
	  
	break;
      }
    default:
      {
	MayDay::Error("CdrPlasmaJSON::advanceReactionNetwork -- logic bust");
	break;
      }
    }

    // Modify by other parameters.
    k *= m_plasmaReactionEfficiencies.at(i)(E, a_pos);

    // This is a hook that uses the Soloviev correction. It modifies the reaction rate according to k = k * (1 + (E.D*grad(n))/(K * n * E^2) where
    // K is the electron mobility. 
    if( (m_plasmaReactionSolovievCorrection.at(i)).first){
      const int species = (m_plasmaReactionSolovievCorrection.at(i)).second;

      const Real&     N  = cdrDensities            [species];
      const Real&     mu = cdrMobilities           [species];
      const Real&     D  = cdrDiffusionCoefficients[species];
      const RealVect& g  = cdrGradients            [species];

      // Compute correction factor 1 + E.(D*grad(n))/(K * n * E^2). But since v = mu*E we put it as 1+ E.(D*grad(n)) / (n * |v| * |E|).
      constexpr Real safety = 1.0;
      const     Real fcorr  = 1.0 + (a_E.dotProduct(D*g)) / ( safety + N * mu * E * E);

      // No negative rates please. 
      k *= std::max(fcorr, 0.0);
    }

    // Fire the reaction.
    FirePlasmaReaction(k, m_plasmaReactions[i]);
  }

  // Iterate through the photon reactions. 


  // If using stochastic photons -- then we need to run Poisson sampling of the photons.
  if(m_discretePhotons){
    MayDay::Warning("CdrPlasmaJSON::advanceReactionNetwork -- discrete photons not yet done");
  }


  return;
}



Vector<RealVect> CdrPlasmaJSON::computeCdrDriftVelocities(const Real         a_time,
							  const RealVect     a_position,
							  const RealVect     a_E,
							  const Vector<Real> a_cdrDensities) const {

  // I really hate Chombo sometimes.
  const std::vector<Real>& cdrDensities = ((Vector<Real>&) a_cdrDensities).stdVector();  

  // Compute the mobilities for each species. 
  const std::vector<Real> mu = this->computePlasmaSpeciesMobilities(a_position, a_E, cdrDensities);

  // Return vector
  Vector<RealVect> velocities(m_numCdrSpecies, RealVect::Zero);

  // Make sure v = +/- mu*E depending on the sign charge. 
  for (int i = 0; i < a_cdrDensities.size(); i++){
    const int Z = m_CdrSpecies[i]->getChargeNumber();

    if(Z > 0){
      velocities[i] = + mu[i] * a_E;
    }
    else if(Z < 0){
      velocities[i] = - mu[i] * a_E;
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
      case LookupMethod::FunctionEN:
	{
	  const Real E = a_E.vectorLength();	  
	  const Real N = m_gasDensity(a_position);

	  Dco = m_diffusionFunctionsEN.at(i)(E, N);

	  break;
	}
      case LookupMethod::TableEN:
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
						      const Vector<Real> a_rteFluxes,
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
						       const Vector<Real> a_rteFluxes,
						       const Vector<Real> a_extrapCdrFluxes) const {
  //  MayDay::Warning("CdrPlasmaJSON::computeCdrDielectricFluxes -- don't know how to do this yet");

  Vector<Real> fluxes(m_numCdrSpecies, 0.0);
}

Vector<Real> CdrPlasmaJSON::computeCdrDomainFluxes(const Real           a_time,
						   const RealVect       a_pos,
						   const int            a_dir,
						   const Side::LoHiSide a_side,
						   const RealVect       a_E,
						   const Vector<Real>   a_cdrDensities,
						   const Vector<Real>   a_cdrVelocities,
						   const Vector<Real>   a_cdrGradients,
						   const Vector<Real>   a_rteFluxes,
						   const Vector<Real>   a_extrapCdrFluxes) const {
  //  MayDay::Warning("CdrPlasmaJSON::computeCdrDomainFluxes -- don't know how to do this yet");
  
  return Vector<Real>(m_numCdrSpecies, 0.0);      
}

Real CdrPlasmaJSON::initialSigma(const Real a_time, const RealVect a_pos) const {
  return m_initialSigma(a_pos, a_time);
}


#include <CD_NamespaceFooter.H>
