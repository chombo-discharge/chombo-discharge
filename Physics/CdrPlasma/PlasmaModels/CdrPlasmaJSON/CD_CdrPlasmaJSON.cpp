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
#include <CD_NamespaceHeader.H>

using namespace Physics::CdrPlasma;

CdrPlasmaJSON::CdrPlasmaJSON(){
  CH_TIME("CdrPlasmaJSON::CdrPlasmaJSON()");
  
  this->parseOptions();
  this->parseJSON();
  this->instantiateCdrSpecies();
  this->initializeSigma();

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
      name = species["name"];
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
      Z = species["Z"];
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
    if(species.contains("mobility")){
      if(species["mobility"].contains("mobile")) {
	mobile = species["mobility"]["mobile"].get<bool>();
      }
      else{
	const std::string err = "CdrPlasmaJSON::instantiateSpecies -- field 'mobile' is missing for species " + name;
	pout() << err << endl;	
	MayDay::Error(err.c_str());
      }
    }

    // Find whether or not the species is diffusive. If it contains a 'diffusion' field, it must also contain the 'diffusive' field. 
    if(species.contains("mobility")){
      if(species["diffusion"].contains("diffusive")) {
	mobile = species["diffusion"]["diffusive"].get<bool>();
      }
      else{
	const std::string err = "CdrPlasmaJSON::instantiateSpecies -- field 'diffusive' is missing for species " + name;
	pout() << err << endl;
	MayDay::Error(err.c_str());
      }
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

	// Add exponential fall-off
	if(data.contains("exp_falloff")){
	  const Real     amplitude = data["exp_falloff"]["amplitude"].get<Real>();	  
	  const Real     l         = data["exp_falloff"]["l"].        get<Real>();
	  const Real     h0        = data["exp_falloff"]["h0"].       get<Real>();	  
	  const RealVect axis      = RealVect(D_DECL(data["exp_falloff"]["axis"][0].get<Real>(),
						     data["exp_falloff"]["axis"][1].get<Real>(),
						     data["exp_falloff"]["axis"][2].get<Real>()));

	  const RealVect n          = axis/axis.vectorLength();
	  const RealVect delta      = a_point - h0*n;


	  ret += amplitude * exp(-delta.dotProduct(n)/l);
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
      m_trackedCdrSpeciesMap.emplace(std::make_pair(name, m_trackedCdrSpecies.size()));
      m_trackedCdrSpecies.emplace_back(RefCountedPtr<CdrSpecies> (new CdrSpeciesJSON(name, Z, diffusive, mobile, initFunc)));
      m_trackedMap.emplace(std::make_pair(name, true));

      m_plotCdr.emplace(std::make_pair(name, false));
    }
    else {
      m_untrackedCdrSpeciesMap.emplace(std::make_pair(name, m_untrackedCdrSpecies.size()));
      m_untrackedCdrSpecies.emplace_back(RefCountedPtr<CdrSpecies> (new CdrSpeciesJSON(name, Z, false, false, initFunc)));
      m_trackedMap.emplace(std::make_pair(name, false));      

      m_plotCdr.emplace(std::make_pair(name, plot));
    }
  }
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
      ret.push_back(plotCdr.first + " density (untracked)");
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

const RefCountedPtr<CdrSpecies>& CdrPlasmaJSON::getUntrackedCdrSpecies(const std::string a_name) const {
  return m_untrackedCdrSpecies[m_untrackedCdrSpeciesMap.at(a_name)];
}

Real CdrPlasmaJSON::getUntrackedCdrDensity(const std::string& a_name, const RealVect& a_pos, const Real& a_time) const {
  const RefCountedPtr<CdrSpecies>& species = this->getUntrackedCdrSpecies(a_name);
  
  return species->initialData(a_pos, a_time);
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
  //  MayDay::Warning("CdrPlasmaJSON::computeCdrDriftVelocities -- don't know how to do this yet");

  return Vector<RealVect>(m_numCdrSpecies, RealVect::Zero);
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
