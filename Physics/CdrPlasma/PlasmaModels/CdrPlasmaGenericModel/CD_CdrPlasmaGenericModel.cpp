/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaGenericModel.cpp
  @brief  Implementation of CD_CdrPlasmaGenericModel.H
  @author Robert Marskar
*/

// Std includees
#include <iostream>
#include <fstream>

// Chombo includes
#include <ParmParse.H>
#include <CH_Timer.H>

// Our includes
#include <CD_CdrPlasmaGenericModel.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::CdrPlasma;

CdrPlasmaGenericModel::CdrPlasmaGenericModel(){
  CH_TIME("CdrPlasmaGenericModel::CdrPlasmaGenericModel()");
  
  this->parseOptions();
  this->parseJSON();
  this->instantiateCdrSpecies();
}

CdrPlasmaGenericModel::~CdrPlasmaGenericModel(){
  CH_TIME("CdrPlasmaGenericModel::~CdrPlasmaGenericModel()");
}

void CdrPlasmaGenericModel::parseOptions() {
  CH_TIME("CdrPlasmaGenericModel::parseOptions");
  
  ParmParse pp("CdrPlasmaGenericModel");

  pp.get("verbose",        m_verbose);
  pp.get("chemistry_file", m_jsonFile);
}

void CdrPlasmaGenericModel::parseJSON() {
  CH_TIME("CdrPlasmaGenericModel::parseJSON");
  if(m_verbose){
    pout() << "CdrPlasmaGenericModel::parseJSON -- file is = " << m_jsonFile << endl;
  }

  // Parse the JSON file
  std::ifstream istream(m_jsonFile);
  istream >> m_json;  
}

void CdrPlasmaGenericModel::instantiateCdrSpecies() {
  CH_TIME("CdrPlasmaGenericModel::instantiateCdrSpecies");
  if(m_verbose){
    pout() << "CdrPlasmaGenericModel::instantiateCdrSpecies()" << endl;
  }

  // Iterate through all species defined in the JSON file. 
  for (const auto& species : m_json["species"]){

    std::string name        = "default_name";
    int         Z           = 0;    
    bool        tracked     = false;
    bool        mobile      = false;
    bool        diffusive   = false;
    bool        hasInitData = false;

    std::function<Real(const RealVect, const Real)> initFunc = [](const RealVect, const Real) -> Real {
      return 0.0;
    };

    // Find the species name.
    if(species.contains("name")){
      name = species["name"];
    }
    else{
      const std::string err = "CdrPlasmaGenericModel::instantiateSpecies -- name is missing for field " + species.dump();
      pout() << err << endl;
      MayDay::Error(err.c_str());
    }

    // Figure out if the species is tracked or not. 
    if(species.contains("tracked")){
      tracked = species["tracked"].get<bool>();
    }
    else{
      const std::string err = "CdrPlasmaGenericModel::instantiateSpecies -- field 'tracked' is missing for species " + name;
      pout() << err << endl;     
      MayDay::Error(err.c_str());
    }

    // Find the charge number.
    if(species.contains("Z")){
      Z = species["Z"];
    }
    else{
      const std::string err = "CdrPlasmaGenericModel::instantiateSpecies -- field 'Z' is missing for species " + name;
      pout() << err << endl;
      MayDay::Error(err.c_str());
    }

    // Find whether or not the species is mobile. If it contains a mobility field, it must also contain the 'mobile' field. 
    if(species.contains("mobility")){
      if(species["mobility"].contains("mobile")) {
	mobile = species["mobility"]["mobile"].get<bool>();
      }
      else{
	const std::string err = "CdrPlasmaGenericModel::instantiateSpecies -- field 'mobile' is missing for species " + name;
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
	const std::string err = "CdrPlasmaGenericModel::instantiateSpecies -- field 'diffusive' is missing for species " + name;
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
	  const Real     h0        = data["exp_falloff"]["h0"].        get<Real>();	  
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
      pout() << "CdrPlasmaGenericModel::instantiateSpecies: instantiating species" << "\n"
	      << "\tName        = " << name        << "\n"
	      << "\tTracked     = " << tracked     << "\n"
	      << "\tZ           = " << Z           << "\n"
	      << "\tMobile      = " << mobile      << "\n"
	      << "\tDiffusive   = " << diffusive   << "\n"
	      << "\tInitialData = " << hasInitData << "\n";    	
    }

    // Instantiate the species.
    if(tracked){
      m_trackedCdrSpeciesMap.emplace(std::make_pair(m_trackedCdrSpecies.size(), name));
      m_trackedCdrSpecies.emplace_back(RefCountedPtr<CdrSpecies> (new CdrPlasmaSpecies(name, Z, diffusive, mobile, initFunc)));
    }
    else {
      m_untrackedCdrSpeciesMap.emplace(std::make_pair(m_untrackedCdrSpecies.size(), name));
      m_untrackedCdrSpecies.emplace_back(RefCountedPtr<CdrSpecies> (new CdrPlasmaSpecies(name, Z, false, false, initFunc)));
    }
  }
}

#include <CD_NamespaceFooter.H>
