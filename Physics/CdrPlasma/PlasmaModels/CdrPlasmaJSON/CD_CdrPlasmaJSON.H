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
#include <CD_CdrPlasmaReaction.H>
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
	@brief Enum class for distinguishing types of mobility lookups. 
      */
      enum class LookupMethod {
	Constant,
	FunctionLFA,
	TableLFA
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
	@brief Gas pressure (in atmospheres)
      */
      std::function<Real(const RealVect a_position)> m_gasPressure;

      /*!
	@brief Gas temperature (in Kelvin)
      */
      std::function<Real(const RealVect a_position)> m_gasTemperature;

      /*!
	@brief Gas number density (in m^(-3))
      */
      std::function<Real(const RealVect a_position)> m_gasDensity;

      /*!
	@brief Initial surface charge
      */
      std::function<Real(const Real a_time, const RealVect a_pos) > m_initialSigma;

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

      /*!
	@brief Map for finding a plasma reaction with a particular description
      */
      std::map<std::string, int> m_plasmaReactionsMap;

      /*!
	@brief Plasma reactions. 
      */
      std::vector<CdrPlasmaReaction> m_plasmaReactions;
      
      /*!
	@brief Get a species. This will return the species whether it is tracked or not. This yields an error if the species does not exist. 
	@param[in] a_name Species name. 
      */
      const RefCountedPtr<CdrSpecies>& getCdrSpecies(const std::string a_name) const;

      /*!
	@brief Get a CdrSpecies density.
	@details This will return the density for the species a_name. If the species is actively tracked, then it will return the entry
	in a_cdrDensities that corresponds to the species name. If the species is NOT tracked, it will return the initial density of the species. 
	@param[in] a_name         Species name. 
	@param[in] a_cdrDensities List of CDR densities for tracked species. 
	@param[in] a_pos          Physical coordinates
	@param[in] a_time         Time
      */
      Real getCdrDensity(const std::string& a_name, const Vector<Real>& a_cdrDensities, const RealVect& a_pos, const Real& a_time) const;

      /*!
	@brief Get neutral species density
	@param[in] a_name Species name
	@param[in] a_pos  Physical coordinates
      */      
      virtual Real getNeutralSpeciesDensity(const std::string a_name, const RealVect& a_position) const;

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
	@brief Initialize photon species
	@details This will initialize the radiative transfer species absed on the "photon_species" field in the JSON file. 
      */
      virtual void initializePhotonSpecies();

      /*!
	@brief Initialize species mobilities
      */
      virtual void parseMobilities();

      /*!
	@brief Initialize species diffusion coefficients
      */
      virtual void parseDiffusion();

      /*!
	@brief Parse plasma reactions
      */
      virtual void parsePlasmaReactions();

      /*!
	@brief Parses a reaction string into reactangs and products
	@param[out] a_reactants Left-hand side of reaction
	@param[out] a_products  Right-hand side of reaction
	@param[in]  _reaction   Reaction string. Must be in format "a + b + c -> e + f + g". 
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
	@brief Remove whitespace from string
      */
      std::string trim(const std::string& a_string) const;
    };
  }
}

#include <CD_NamespaceFooter.H>

#endif
