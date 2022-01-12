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
#include <CD_Units.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::CdrPlasma;

CdrPlasmaJSON::CdrPlasmaJSON(){
  CH_TIME("CdrPlasmaJSON::CdrPlasmaJSON()");
  
  this->parseOptions();
  this->parseJSON();
  this->initializeSigma();  
  this->instantiateCdrSpecies();
  this->parseMobilities();

  // Populate the stuff that is needed by CdrPlasmaPhysics
  m_RtSpecies.resize(0);
  m_CdrSpecies = m_trackedCdrSpecies;

  m_numCdrSpecies = m_CdrSpecies.size();
  m_numRtSpecies  = 0;
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

void CdrPlasmaJSON::instantiateCdrSpecies() {
  CH_TIME("CdrPlasmaJSON::instantiateCdrSpecies");
  if(m_verbose){
    pout() << "CdrPlasmaJSON::instantiateCdrSpecies()" << endl;
  }

  // Iterate through all species defined in the JSON file. 
  for (const auto& species : m_json["species"]){

    std::string name        = "default_name";
    int         Z           = 0;    
    bool        tracked     = false;
    bool        mobile      = false;
    bool        diffusive   = false;
    bool        hasInitData = false;
    bool        plot        = false; // Irrelevant for tracked species

    std::function<Real(const RealVect, const Real)> initFunc = [](const RealVect, const Real) -> Real {
      return 0.0;
    };

    // Find the species name.
    if(species.contains("name")){
      name = species["name"].get<std::string>();
    }
    else{
      const std::string err = "CdrPlasmaJSON::instantiateSpecies -- name is missing for field " + species.dump();
      pout() << err << endl;
      MayDay::Error(err.c_str());
    }

    // Figure out if the species is tracked or not. 
    if(species.contains("tracked")){
      tracked = species["tracked"].get<bool>();
    }
    else{
      const std::string err = "CdrPlasmaJSON::instantiateSpecies -- field 'tracked' is missing for species " + name;
      pout() << err << endl;     
      MayDay::Error(err.c_str());
    }

    // Find the charge number.
    if(species.contains("Z")){
      Z = species["Z"].get<int>();
    }
    else{
      const std::string err = "CdrPlasmaJSON::instantiateSpecies -- field 'Z' is missing for species " + name;
      pout() << err << endl;
      MayDay::Error(err.c_str());
    }

    // Figure out if we should plot the species
    if(species.contains("plot")){
      plot = species["plot"].get<bool>();
    }
    else{
      plot = false;
    }

    // Find whether or not the species is mobile. If it contains a mobility field, it must also contain the 'mobile' field.
    if(species.contains("mobile")){
      mobile = species["mobile"].get<bool>();
    }
    else{
      mobile = false;
    }

    // Find whether or not the species is diffusive. If it contains a 'diffusion' field, it must also contain the 'diffusive' field. 
    if(species.contains("diffusive")){
      diffusive = species["diffusive"].get<bool>();
    }
    else{
      diffusive = false;
    }
    
    // Initial data -- we build an incrementing function. 
    if(species.contains("initial_data")){
      hasInitData = true;

      initFunc = [data = species["initial_data"]] (const RealVect a_point, const Real a_time) -> Real {
	Real ret = 0.0;

	// Add uniform density. 
	if(data.contains("uniform")){
	  ret += data["uniform"].get<Real>();
	}

	// Add gaussian seed.
	if(data.contains("gauss2")){
	  const Real     amplitude = data["gauss2"]["amplitude"].get<Real>();	  
	  const Real     radius    = data["gauss2"]["radius"].   get<Real>();
	  const RealVect center    = RealVect(D_DECL(data["gauss2"]["position"][0].get<Real>(),
						     data["gauss2"]["position"][1].get<Real>(),
						     data["gauss2"]["position"][2].get<Real>()));
	  const RealVect delta      = center - a_point;

	  ret += amplitude * exp(-delta.dotProduct(delta)/(2*std::pow(radius,2)));
	}

	// Add super-Gaussian seed.
	if(data.contains("gauss4")){
	  const Real     amplitude = data["gauss4"]["amplitude"].get<Real>();	  
	  const Real     radius    = data["gauss4"]["radius"].   get<Real>();
	  const RealVect center    = RealVect(D_DECL(data["gauss4"]["position"][0].get<Real>(),
						     data["gauss4"]["position"][1].get<Real>(),
						     data["gauss4"]["position"][2].get<Real>()));
	  const RealVect delta      = center - a_point;


	  ret += amplitude * exp(-std::pow(delta.dotProduct(delta),2)/(2*std::pow(radius, 4)));
	}

	// Ideal gas. 
	if(data.contains("ideal_gas")){
	  const Real T = data["ideal_gas"]["T"].             get<Real>();
	  const Real P = data["ideal_gas"]["P"].             get<Real>();
	  const Real f = data["ideal_gas"]["molar_fraction"].get<Real>();

	  ret += f * (P * Units::atm2pascal * Units::Na)/ (T * Units::R);
	}

	return ret;
      };

      initFunc(RealVect::Zero, 0.0);
    }
    else{
      initFunc = [](const RealVect a_point, const Real a_time) -> Real {
	return 0.0;
      };
    }

    // Print out a message if we're verbose.
    if(m_verbose){
      pout() << "CdrPlasmaJSON::instantiateSpecies: instantiating species" << "\n"
	     << "\tName        = " << name        << "\n"
	     << "\tTracked     = " << tracked     << "\n"
	     << "\tZ           = " << Z           << "\n"
	     << "\tMobile      = " << mobile      << "\n"
	     << "\tDiffusive   = " << diffusive   << "\n"
	     << "\tInitialData = " << hasInitData << "\n";    	
    }

    // Instantiate the species.
    if(tracked){
      const int num = m_trackedCdrSpecies.size();

      // Make the string-int map encodings. 
      m_trackedCdrSpeciesMap.       emplace(std::make_pair(name, num ));
      m_trackedCdrSpeciesInverseMap.emplace(std::make_pair(num , name));

      // Tracked/plot maps.
      m_trackedMap.emplace(std::make_pair(name, true));      
      m_plotCdr.   emplace(std::make_pair(name, false));

      // Push the JSON entry and the new CdrSpecies to corresponding vectors. 
      m_trackedCdrSpecies.emplace_back(RefCountedPtr<CdrSpecies> (new CdrSpeciesJSON(name, Z, diffusive, mobile, initFunc)));
      m_trackedCdrJSON.   emplace_back(species);

    }
    else {
      const int num = m_untrackedCdrSpecies.size();
      
      // Make the name-int encoding. 
      m_untrackedCdrSpeciesMap.       emplace(std::make_pair(name, num));
      m_untrackedCdrSpeciesInverseMap.emplace(std::make_pair(num , name));

      // Tracked/plot maps.
      m_trackedMap.emplace(std::make_pair(name, false));      
      m_plotCdr.   emplace(std::make_pair(name, plot ));

      // Push the JSON entry and the new CdrSpecies to corresponding vectors. 
      m_untrackedCdrSpecies.emplace_back(RefCountedPtr<CdrSpecies> (new CdrSpeciesJSON(name, Z, false, false, initFunc)));      
      m_untrackedCdrJSON.   emplace_back(species);      
    }
  }
}

void CdrPlasmaJSON::parseMobilities() {
  CH_TIME("CdrPlasmaJSON::parseMobilities");
  if(m_verbose){
    pout() << "CdrPlasmaJSON::parseMobilities - file is = " << m_jsonFile << endl;
  }

  // Iterate through all tracked species.
  for (const auto& species : m_trackedCdrJSON){
    const std::string name = species["name"].get<std::string>();
    const int         idx  = m_trackedCdrSpeciesMap.at(name);

    // Only do tracked species.
    bool reallyMobile = false;

    if(species.contains("mobile")){
      if(species["mobile"].get<bool>()) {
	reallyMobile = true;
      }
    }

    // Create the map over mobile species
    m_isMobile.emplace(std::make_pair(idx, reallyMobile));

    if(reallyMobile){
      if(!(species.contains("mobility"))){ // This is an error -- if a species is tracked AND is mobile it must contain the field 'mobility'
	const std::string err = "species " + name + " is mobile but file does not contain field 'mobility'";
	pout() << err << endl;
	MayDay::Error(err.c_str());	
      }
      else{
	const json& mobilityJSON = species["mobility"];

	// Get the mobility lookup method. This must either be a constant, a function, or a table. We parse these cases differently.
	
	const std::string lookup = mobilityJSON["lookup"].get<std::string>();
	if(lookup == "constant"){

	  // User specified a constant mobility. We look for a field 'value' in the JSON file and set the mobility from that. If the
	  // field does not exist then it's an error. 
	  if(!(mobilityJSON.contains("value"))) {
	    const std::string err = "CdrPlasmaJSON::parseMobilities -- constant lookup was specified but I can't find the field 'value'";
	    pout() << err << endl;
	    MayDay::Error(err.c_str());
	  }
	  else{
	    const Real value = mobilityJSON["value"].get<Real>();

	    m_mobilityLookup.  emplace(std::make_pair(idx, LookupMethod::Constant));
	    m_mobilityConstant.emplace(std::make_pair(idx, value                 ));
	  }
	}
	else if (lookup == "function"){
	  pout() << "using function" << endl;
	}
	else if(lookup == "table"){
	  pout() << "using table" << endl;
	}	
	else{
	  const std::string err = "CdrPlasmaJSON::parseMobilities -- lookup = '" + lookup + "' was specified but this is not 'constant', 'function', or 'table'";
	  pout() << err << endl;
	  MayDay::Error(err.c_str());
	}
      }
    }
  }
}

int CdrPlasmaJSON::getNumberOfPlotVariables() const {
  int ret = 0;

  for (const auto& plotCdr : m_plotCdr){
    if(plotCdr.second){
      ret++;
    }
  }

  return ret;
}

Vector<std::string> CdrPlasmaJSON::getPlotVariableNames() const {
  Vector<std::string> ret;

  for (const auto& plotCdr : m_plotCdr){
    if(plotCdr.second){
      ret.push_back(plotCdr.first + " untracked density");
    }
  }

  return ret;
}

Vector<Real> CdrPlasmaJSON::getPlotVariables(const Vector<Real> a_cdrDensities,
					     const Vector<Real> a_rteDensities,
					     const RealVect     a_E,
					     const RealVect     a_position,
					     const Real         a_time) const {
  Vector<Real> ret;

  for (const auto & plotCdr : m_plotCdr){
    if(plotCdr.second){
      ret.push_back(this->getUntrackedCdrDensity(plotCdr.first, a_position, a_time));
    }
  }
  
  return ret;
}

const RefCountedPtr<CdrSpecies>& CdrPlasmaJSON::getCdrSpecies(const std::string a_name) const {
  return m_trackedMap.at(a_name) ? this->getTrackedCdrSpecies(a_name) : this->getUntrackedCdrSpecies(a_name);
}

const RefCountedPtr<CdrSpecies>& CdrPlasmaJSON::getUntrackedCdrSpecies(const std::string a_name) const {
  return m_untrackedCdrSpecies[m_untrackedCdrSpeciesMap.at(a_name)];
}

const RefCountedPtr<CdrSpecies>& CdrPlasmaJSON::getTrackedCdrSpecies(const std::string a_name) const {
  return m_trackedCdrSpecies[m_trackedCdrSpeciesMap.at(a_name)];
}

Real CdrPlasmaJSON::getCdrDensity(const std::string& a_name, const Vector<Real>& a_cdrDensities, const RealVect& a_pos, const Real& a_time) const {
  return m_trackedMap.at(a_name) ? this->getTrackedCdrDensity(a_name, a_cdrDensities) : this->getUntrackedCdrDensity(a_name, a_pos, a_time);
}

Real CdrPlasmaJSON::getUntrackedCdrDensity(const std::string& a_name, const RealVect& a_pos, const Real& a_time) const {
  const RefCountedPtr<CdrSpecies>& species = this->getUntrackedCdrSpecies(a_name);
  
  return species->initialData(a_pos, a_time);
}

Real CdrPlasmaJSON::getTrackedCdrDensity(const std::string& a_name, const Vector<Real>& a_cdrDensities) const {
  return a_cdrDensities[m_trackedCdrSpeciesMap.at(a_name)];
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
							  const RealVect     a_pos,
							  const RealVect     a_E,
							  const Vector<Real> a_cdrDensities) const {
  Vector<RealVect> velocities(m_numCdrSpecies, RealVect::Zero);

  for (const auto& mobileSpecies : m_isMobile){
    if(mobileSpecies.second){ // Ok, species is mobile. 
      const int idx = mobileSpecies.first;
      const int Z   = m_trackedCdrSpecies[idx]->getChargeNumber();

      // Figure out the mobility.
      const LookupMethod& method = m_mobilityLookup.at(idx);

      Real mobility = 0.0;
      
      switch(method) {
      case LookupMethod::Constant:
	{
	  mobility = m_mobilityConstant.at(idx);
	  break;
	}
      case LookupMethod::Function:
	{
	  MayDay::Error("CdrPlasmaJSON::computeCdrDriftVelocities -- function not implemented yet");
	  break;
	}
      case LookupMethod::Table:
	{
	  MayDay::Error("CdrPlasmaJSON::computeCdrDriftVelocities -- table not implemented yet");
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

      velocities[idx] = mobility * a_E;
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
