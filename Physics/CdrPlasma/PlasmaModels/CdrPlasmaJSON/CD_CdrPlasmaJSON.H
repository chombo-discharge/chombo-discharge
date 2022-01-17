/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaJSON.H
  @brief  Declaration of a generic CdrPlasma model.
  @author Robert Marskar
*/

#ifndef CD_CdrPlasmaJSON_H
#define CD_CdrPlasmaJSON_H

// Std includes
#include <map>
#include <memory>
#include <string>

// Third-party includes
#include <nlohmann/json.hpp>

// Our includes
#include <CD_CdrPlasmaPhysics.H>
#include <CD_NeutralSpeciesJSON.H>
#include <CD_CdrSpeciesJSON.H>
#include <CD_RteSpeciesJSON.H>
#include <CD_CdrPlasmaReactionJSON.H>
#include <CD_LookupTable.H>
#include <CD_NamespaceHeader.H>

using json = nlohmann::json;

namespace Physics {
  namespace CdrPlasma {

    /*!
      @brief Class which implements CdrPlasmaPhysics and parses plasma chemistry from a JSON input file
    */
    class CdrPlasmaJSON : public CdrPlasmaPhysics {
    public:

      /*!
	@brief Function for encapsulating operations f = f(E,N). 
	field in Townsend units
	@param[in] a_E Electric field in SI units.
	@param[in] a_N Neutral density
	@return Returns f = f(E,N). 
	@note a_E is in physical units. 
      */
      using FunctionEN = std::function<Real(const Real a_E, const Real a_N)>;

      /*!
	@brief Function for encapsulating a function f = f(x) where x is the physical coordinates
	@param[in] a_position Physical coordinates
	@return Returns f(x)
      */
      using FunctionX = std::function<Real(const RealVect a_position)>;

      /*!
	@brief Function for encapsulating a function f = f(T) where T is the temperature of some species
	@param[in] a_T Some temperature. 
      */
      using FunctionT = std::function<Real(const Real a_T)>;

      /*!
	@brief Function for encapsulating a function f = f(T1, T2) where T1/T2 are temperatures of two species. 
	@param[in] a_T1 Some species temperature. 
	@param[in] a_T2 Some other species temperature. 
      */
      using FunctionTT = std::function<Real(const Real a_T1, const Real a_T2)>;      

      /*!
	@brief Default constructor.
      */
      CdrPlasmaJSON();

      /*!
	@brief Destructor
      */
      virtual ~CdrPlasmaJSON();

      /*!
	@brief Get number of plot variables for this physics class. 
	@details This is used by CdrPlasmaStepper for pre-allocating data that will be put in a plot file. The current implementation 
	plots the gas pressure, temperature, and density
      */
      virtual int getNumberOfPlotVariables() const override;

      /*!
	@brief Get plot variable names. The names and positions between this routine and getPlotVariables must be consistent!
      */
      virtual Vector<std::string> getPlotVariableNames() const override;

      /*!
	@brief Provide plot variables. This is used by CdrPlasmaStepper when writing plot files. 
	@details The plot variables should have length this->getNumberOfPlotVariables (as should getPlotVariableNames)
	@param[in] a_cdrDensities CDR densities
	@param[in] a_rteDensities RTE densities
	@param[in] a_position     Physical coordinates
	@param[in] a_time         Time. 
      */
      virtual Vector<Real> getPlotVariables(const Vector<Real> a_cdrDensities,
					    const Vector<Real> a_rteDensities,
					    const RealVect     a_E,
					    const RealVect     a_position,
					    const Real         a_time) const override;

      /*!
	@brief Compute alpha. Should return Townsend ionization coefficient. 
	@details This function is mostly used for the cell tagging classes
	@param[in] a_E Electric field.
      */
      virtual Real computeAlpha(const RealVect a_E) const override;
      
      /*!
	@brief Routine intended for advancing a reaction network over a time a_dt. 
	@details This routine assumes that the subsequent advance is in the form phi^(k+1) = phi^k + S*a_dt. Thus, this routine exists such that users can EITHER
	fill a_cdrSources and a_rteSources directly with an explicit rule, OR they can perform a fully implicit advance within this routine and set S from that. 
	@param[out] a_cdrSources    Source terms for CDR equations.
	@param[out] a_rteSources    Source terms for RTE equations.
	@param[in]  a_cdrDensities  Grid-based density for particle species.
	@param[in]  a_cdrGradients  Grid-based gradients for particle species.
	@param[in]  a_rteDensities  Grid-based densities for photons.
	@param[in]  a_E             Electric field.
	@param[in]  a_pos           Position in space.
	@param[in]  a_dx            Grid resolution. 
	@param[in]  a_dt            Advanced time.
	@param[in]  a_time          Current time.
	@param[in]  a_kappa         Grid cell unit volume. 
      */
      virtual void advanceReactionNetwork(Vector<Real>&          a_cdrSources,
					  Vector<Real>&          a_rteSources,
					  const Vector<Real>     a_cdrDensities,
					  const Vector<RealVect> a_cdrGradients,
					  const Vector<Real>     a_rteDensities,
					  const RealVect         a_E,
					  const RealVect         a_pos,
					  const Real             a_dx,
					  const Real             a_dt,
					  const Real             a_time,
					  const Real             a_kappa) const override;

      /*!
	@brief Compute velocities for the CDR equations
	@param[in] a_time         Time
	@param[in] a_pos          Position
	@param[in] a_E            Electric field
	@param[in] a_cdrDensities CDR densities
	@return Returns the drift velocities for each CDR species. The vector ordering is the same as m_CdrSpecies. 
      */
      virtual Vector<RealVect> computeCdrDriftVelocities(const Real         a_time,
							 const RealVect     a_pos,
							 const RealVect     a_E,
							 const Vector<Real> a_cdrDensities) const override;

      /*!
	@brief Compute diffusion coefficients for the CDR equations. 
	@param[in] a_time         Time
	@param[in] a_pos          Position
	@param[in] a_E            Electric field
	@param[in] a_cdrDensities CDR densities
	@return Returns the diffusion coefficients for each CDR species. The vector ordering is the same as m_CdrSpecies. 
      */
      virtual Vector<Real> computeCdrDiffusionCoefficients(const Real         a_time,
							   const RealVect     a_pos,
							   const RealVect     a_E,
							   const Vector<Real> a_cdrDensities) const override;


      /*!
	@brief Compute CDR fluxes on electrode-gas interfaces. This is used as a boundary condition in the CDR equations. 
	@param[in] a_time            Time
	@param[in] a_pos             Position
	@param[in] a_normal          Boundary normal vector. This points into the gas phase. 
	@param[in] a_E               Electric field
	@param[in] a_cdrVelocities   CDR velocities. Normal component only. 
	@param[in] a_cdrDensities    CDR densities. 
	@param[in] a_cdrGradients    Normal gradients of cdr densities
	@param[in] a_rteFluxeses     RTE fluxes (normal component only)
	@param[in] a_extrapCdrFluxes Extrapolated fluxes from the gas side. 
	@return Returns the flux on an electrode interface cell. The vector ordering must be the same as m_CdrSpecies. 
      */
      virtual Vector<Real> computeCdrElectrodeFluxes(const Real         a_time,
						     const RealVect     a_pos,
						     const RealVect     a_normal,
						     const RealVect     a_E,
						     const Vector<Real> a_cdrDensities,
						     const Vector<Real> a_cdrVelocities,
						     const Vector<Real> a_cdrGradients,
						     const Vector<Real> a_rteFluxeses,
						     const Vector<Real> a_extrapCdrFluxes) const override;

      /*!
	@brief Compute CDR fluxes on dielectric-gas interfaces. This is used as a boundary condition in the CDR equations. 
	@param[in] a_time            Time
	@param[in] a_pos             Position
	@param[in] a_normal          Normal vector. This points into the gas phase. 
	@param[in] a_E Electric      field
	@param[in] a_cdrDensities    CDR densities (on the EB)
	@param[in] a_cdrVelocities   Normal component of CDR velocities (on the EB).
	@param[in] a_cdrGradients    Normal gradients of cdr densities
	@param[in] a_rteFluxeses     RTE fluxes (normal component only)
	@param[in] a_extrapCdrFluxes Extrapolated fluxes from the gas side. 
	@return Returns the flux on a dielectric interface cell. The vector ordering must be the same as m_CdrSpecies. 
      */
      virtual Vector<Real> computeCdrDielectricFluxes(const Real         a_time,
						      const RealVect     a_pos,
						      const RealVect     a_normal,
						      const RealVect     a_E,
						      const Vector<Real> a_cdrDensities,
						      const Vector<Real> a_cdrVelocities,
						      const Vector<Real> a_cdrGradients,
						      const Vector<Real> a_rteFluxeses,
						      const Vector<Real> a_extrapCdrFluxes) const override;

      /*!
	@brief Compute CDR fluxes through domain sides. This is used as a boundary condition in the CDR equations.
	@param[in] a_time            Time
	@param[in] a_pos             Position
	@param[in] a_dir             Direction (0 = x, 1=y etc)
	@param[in] a_side            Side (low or high side)
	@param[in] a_E               Electric field
	@param[in] a_cdrDensities    CDR densities. 
	@param[in] a_cdrVelocities   CDR velocities (normal component only). 
	@param[in] a_cdrGradients    CDR gradients (normal component only)
	@param[in] a_rteFluxeses     RTE fluxes (normal component only)
	@param[in] a_extrapCdrFluxes Extrapolated fluxes from the gas side. 
      */
      virtual Vector<Real> computeCdrDomainFluxes(const Real           a_time,
						  const RealVect       a_pos,
						  const int            a_dir,
						  const Side::LoHiSide a_side,
						  const RealVect       a_E,
						  const Vector<Real>   a_cdrDensities,
						  const Vector<Real>   a_cdrVelocities,
						  const Vector<Real>   a_cdrGradients,
						  const Vector<Real>   a_rteFluxeses,
						  const Vector<Real>   a_extrapCdrFluxes) const override;

      /*!
	@brief Set the initial surface charge
	@param[in] a_time Time
	@param[in] a_pos Position
      */
      virtual Real initialSigma(const Real a_time, const RealVect a_pos) const override;
      
    protected:

      /*!
	@brief Enum class for distinguishing types of computation methods when computing transport data stuff. 
      */
      enum class LookupMethod {
	Constant,
	FunctionX,
	FunctionT,
	FunctionTT,
	FunctionEN,
	TableEN
      };

      /*!
	@brief Verbose or not
      */
      bool m_verbose;

      /*!
	@brief Plot gas pressure, density, and temperature
      */
      bool m_plotGas;

      /*!
	@brief JSON definition. This is populated when calling parseJSON.
      */
      json m_json;

      /*!
	@brief Input JSON file name
      */
      std::string m_jsonFile;

      /*!
	@brief JSON entries for species in plasma_species
      */
      std::vector<json> m_cdrSpeciesJSON;

      /*!
	@brief JSON entries for species in photon_species
      */
      std::vector<json> m_rteSpeciesJSON;

      /*!
	@brief Initial surface charge
      */
      std::function<Real(const RealVect a_position, const Real a_time)> m_initialSigma;      

      /*!
	@brief Gas pressure (in atmospheres)
      */
      FunctionX m_gasPressure;

      /*!
	@brief Gas temperature (in Kelvin)
      */
      FunctionX m_gasTemperature;

      /*!
	@brief Gas number density (in m^(-3))
      */
      FunctionX m_gasDensity;

      /*!
	@brief Neutral species densities
      */
      std::vector<FunctionX> m_neutralSpeciesDensities;

      /*!
	@brief These are the neutral species.
      */
      std::vector<std::shared_ptr<NeutralSpeciesJSON> > m_neutralSpecies;

      /*!
	@brief Map for figuring out which where in m_neutralSpecies a neutral species is found. 
      */
      std::map<std::string, int> m_neutralSpeciesMap;

      /*!
	@brief Inverse of m_neutralSpeciesMap
      */
      std::map<int, std::string> m_neutralSpeciesInverseMap;

      /*!
	@brief string-int encoding of the CDr species. 
	@details This is needed because we sometimes need to use the species name for indexing in the vector. 
      */
      std::map<std::string, int> m_cdrSpeciesMap;

      /*!
	@brief int-string encoding of the CDR species. 
	@details This is needed because we sometimes need to use the species name for indexing in the vector. 
	@note This is the inverse of m_cdrSpeciesMap
      */
      std::map<int, std::string> m_cdrSpeciesInverseMap;

      /*!
	@brief string-int encoding of the RTE species. 
	@details This is needed because we sometimes need to use the species name for indexing in the vector. 
      */
      std::map<std::string, int> m_rteSpeciesMap;

      /*!
	@brief int-string encoding of the RTE species. 
	@details This is needed because we sometimes need to use the species name for indexing in the vector. 
	@note This is the inverse of m_rteSpeciesMap
      */
      std::map<int, std::string> m_rteSpeciesInverseMap;

      // MOBILITY QUANTITIES BEGIN HERE      

      /*!
	@brief Mobility lookup method for each species. 
      */
      std::map<int, LookupMethod> m_mobilityLookup;

      /*!
	@brief Map for constant mobilities. 
      */
      std::map<int, Real> m_mobilityConstants;

      /*!
	@brief Map for table-based mobilities.
      */
      std::map<int, FunctionEN> m_mobilityFunctionsEN;

      /*!
	@brief Map for function-based mobilities. 
      */
      std::map<int, LookupTable<2> > m_mobilityTablesEN;

      // DIFFUSION  QUANTITIES BEGIN HERE

      /*!
	@brief Diffusion lookup method. 
      */
      std::map<int, LookupMethod> m_diffusionLookup;

      /*!
	@brief Map for constant diffusion coefficients. 
      */
      std::map<int, Real> m_diffusionConstants;

      /*!
	@brief Map for function-based diffusion coefficients. .
      */
      std::map<int, FunctionEN> m_diffusionFunctionsEN;

      /*!
	@brief Map for table-based diffusion coefficients. 
      */
      std::map<int, LookupTable<2> > m_diffusionTablesEN;

      // TEMPERATURE DATA BEGINS HERE

      /*!
	@brief Temperature lookup method
      */
      std::map<int, LookupMethod> m_temperatureLookup;

      /*!
	@brief Constant temperatures.
      */
      std::map<int, FunctionX> m_temperatureConstants;

      /*!
	@brief Temperatures as functions of E/N. 
      */
      std::map<int, LookupTable<2> > m_temperatureTablesEN;

      // REACTION DATA BEGINS HERE

      /*!
	@brief Description of plasma reactions. Only used for I/O
      */
      std::map<int, std::string> m_plasmaReactionDescriptions;

      /*!
	@brief Map for figuring out how to look up the rate for a certain plasma reaction. 
      */
      std::map<int, LookupMethod> m_plasmaReactionLookup;

      /*!
	@brief Constant plasma reaction rates.
      */
      std::map<int, Real> m_plasmaReactionConstants;

      /*!
	@brief Maps for functions of the type k = f(T1,T2) where T1 and T2 are the temperatures of some species. 
	@details This signature is absolute horrific -- what it means is that we have a reaction which should be evaluated
	as f(T1, T2), but we need to know which species are involved. By design, this should be general so that T1 and T2 
	can be the temperatures of any species, including neutral species. So, we make a tuple for indicating which species
	we are talking about. The first index in the tuple is the first species, the second is the second species and the third
	entry in the tuple is the actual function. A special rule occurs if one of the first two indices is < 0 in which case
	the temperature is replaced by the background gas temperature. 
      */
      std::map<int, std::tuple<int, int, FunctionTT> > m_plasmaReactionFunctionsTT;;

      /*!
	@brief Function-based plasma reaction rates. 
      */
      std::map<int, FunctionEN> m_plasmaReactionFunctionsEN;

      /*!
	@brief Map for table-based reaction coefficients (using reduced electric field)
      */
      std::map<int, LookupTable<2> > m_plasmaReactionTablesEN;

      /*!
	@brief Plasma reactions. 
      */
      std::vector<CdrPlasmaReactionJSON> m_plasmaReactions;

      /*!
	@brief Flag for whether or not reaction includes Soloviev energy correction. 
	@details If this is true, the rate for a reaction (in the local field approximation) will be modified as k * (1 + E.(D * grad(phi))/n * 
      */
      std::map<int, std::pair<bool, int> > m_plasmaReactionSolovievCorrection;      

      /*!
	@brief Parse class options
      */
      virtual void parseOptions();

      /*!
	@brief Parse the JSON file
      */
      virtual void parseJSON();

      /*!
	@brief Initialize surface charge
	@details This will initialize the surface charge based on the "sigma" field in the JSON file. If you want more complex initial surface charges, 
	you will have to either extend this routine or overwrite it. 
      */
      virtual void initializeSigma();

      /*!
	@brief Initialize neutral species
      */
      virtual void initializeNeutralSpecies();
      
      /*!
	@brief Initialize species. 
	@details This will initialize the species based on the "plasma_species" field in the JSON file. 
      */
      virtual void initializePlasmaSpecies();

      /*!
	@brief Initialize photon species.
	@details This will initialize the radiative transfer species absed on the "photon_species" field in the JSON file. 
      */
      virtual void initializePhotonSpecies();

      /*!
	@brief Initialize species mobilities.
      */
      virtual void parseMobilities();

      /*!
	@brief Initialize species diffusion coefficients.
      */
      virtual void parseDiffusion();

      /*!
	@brief Initialize species temperatures.
      */
      virtual void parseTemperatures();      

      /*!
	@brief Parse plasma reactions
      */
      virtual void parsePlasmaReactions();

      /*!
	@brief Parses a reaction string into reactangs and products
	@param[out] a_reactants Left-hand side of reaction
	@param[out] a_products  Right-hand side of reaction
	@param[in]  a_reaction  Reaction string. Must be in format "a + b + c -> e + f + g". 
      */
      virtual void parseReactionString(std::vector<std::string>& a_reactants,
				       std::vector<std::string>& a_products,
				       const std::string&        a_reaction) const;

      /*!
	@brief Do a species sanity check
	@details This will make sure that species names are not duplicates. 
      */
      virtual void sanityCheckSpecies() const;

      /*!
	@brief Check if a plasma-reaction makes sense in terms of the species that have been defined.
	@details This will throw errors if the species do not exist or charge is not conserved.
	@param[in] a_reactants Reactants
	@param[in] a_products  Reaction products
	@param[in] a_reaction  Reaction string
      */
      virtual void sanctifyPlasmaReaction(const std::vector<std::string>& a_reactants,
					  const std::vector<std::string>& a_products,
					  const std::string               a_reaction) const;

      /*!
	@brief Get the int-encoding corresponding to species involved in a plasma reaction.
	@param[out] a_plasmaReactants   Plasma reactants
	@param[out] a_neutralReactants  Neutral reactants
	@param[out] a_plasmaProducts    Plasma products
	@param[out] a_photonProducts    Photon products
	@param[in]  a_reactants         Reactants
	@param[in]  a_products          Reaction products
      */
      virtual void getPlasmaReactionProducts(std::list<int>&                 a_plasmaReactants,
					     std::list<int>&                 a_neutralReactants,
					     std::list<int>&                 a_plasmaProducts,
					     std::list<int>&                 a_photonProducts,
					     const std::vector<std::string>& a_reactants,
					     const std::vector<std::string>& a_products) const;


      /*!
	@brief Compute the various plasma species temperatures.
	@param[in] a_position     Physical coordinates
	@param[in] a_E            Electric field (SI units)
	@param[in] a_cdrDensities List of plasma species densities.
      */
      virtual std::vector<Real> computePlasmaSpeciesTemperatures(const RealVect&          a_position,
								 const RealVect&          a_E,
								 const std::vector<Real>& a_cdrDensities) const;

      /*!
	@brief Compute the various plasma species mobilities. 
	@param[in] a_position     Physical coordinates
	@param[in] a_E            Electric field (SI units)
	@param[in] a_cdrDensities List of plasma species densities.
      */
      virtual std::vector<Real> computePlasmaSpeciesMobilities(const RealVect&          a_position,
							       const RealVect&          a_E,
							       const std::vector<Real>& a_cdrDensities) const;      

      /*!
	@brief Throw a parser error
	@param[in] a_error Error code.
      */
      void throwParserError(const std::string a_error) const;

      /*!
	@brief Throw a parser wearning
	@param[in] a_error Warning
      */
      void throwParserWarning(const std::string a_warning) const;

      /*!
	@brief Protect the @ character in a string
	@param[in] a_str String If it contains the at character we throw an error. 
      */
      bool containsAt(const std::string a_str) const;

      /*!
	@brief Remove whitespace from string
      */
      std::string trim(const std::string& a_string) const;

      /*!
	@brief Return true if species exists in map and false otherwise
	@param[in] a_name Neutral species name
      */
      bool isNeutralSpecies(const std::string& a_name) const;            

      /*!
	@brief Return true if species exists in map and false otherwise
	@param[in] a_name Cdr species name
      */
      bool isPlasmaSpecies(const std::string& a_name) const;

      /*!
	@brief Return true if species exists in map and false otherwise
	@param[in] a_name Rte species name
      */
      bool isPhotonSpecies(const std::string& a_name) const;
    };
  }
}

#include <CD_NamespaceFooter.H>

#endif
