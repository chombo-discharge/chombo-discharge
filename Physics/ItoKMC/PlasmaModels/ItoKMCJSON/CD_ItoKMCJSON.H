/* chombo-discharge
 * Copyright © 2023 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoKMCJSON.H
  @brief  Declaration of a ItoKMC plasma model which reads reactions from a JSON file.
  @author Robert Marskar
*/

#ifndef CD_ItoKMCJSON_H
#define CD_ItoKMCJSON_H

// Std includes
#include <map>
#include <memory>
#include <string>

// Third-party includes
#include <nlohmann/json.hpp>

// Our includes
#include <CD_ItoKMCBackgroundSpecies.H>
#include <CD_ItoKMCPhysics.H>
#include <CD_LookupTable.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace ItoKMC {

    /*!
      @brief Implementation of ItoKMCPhysics which parses input data from a JSON file
    */
    class ItoKMCJSON : public ItoKMCPhysics
    {
    public:
      /*!
	@brief Function alias for e.g. initial data.
	@param[in] a_position Physical coordinates
	@param[in] a_time     Time
      */
      using FunctionXt = std::function<Real(const RealVect a_position, const Real a_time)>;

      /*!
	@brief Function for encapsulating operations f = f(E,N). 
	field in Townsend units
	@param[in] a_E Electric field in SI units.
	@param[in] a_N Neutral density.
	@return Returns f = f(E,N). 
	@note a_E is in SI units and the neutral density is in m^-3.
      */
      using FunctionEN = std::function<Real(const Real a_E, const Real a_N)>;

      /*!
	@brief Function for encapsulating a function f = f(x) where x is the physical coordinates
	@param[in] a_position Physical coordinates
	@return Returns f(x)
      */
      using FunctionX = std::function<Real(const RealVect a_position)>;

      /*!
	@brief Function for encapsulating a function f = f(E, x) where E is the electric field at physical coordinates x.
	@param[in] E Electric field magnitude in SI units. 
	@param[in] x Physical coordinates
      */
      using FunctionEX = std::function<Real(const Real E, const RealVect x)>;

      /*!
	@brief Function for encapsulating a function f = f(E, V, x) where E is the electric field at physical coordinates x, and
	V is a grid cell volume.
	@param[in] E Electric field magnitude in SI units. 
	@param[in] V Volume
	@param[in] dx Grid resolution
	@param[in] dt Time step	
	@param[in] x Physical coordinate
	@param[in] phi List of plasma species densities
      */
      using FunctionEVXTP = std::function<
        Real(const Real E, const Real V, const Real dx, const Real dt, const RealVect x, const Vector<Real>& phi)>;

      /*!
	@brief Function which encapsulates 
	@param[in] dx Grid resolution
	@param[in] dx List of plasma species densities.
      */
      using FunctionDXP = std::function<Real(const Real dx, const Vector<Real>& phi)>;

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
	@brief Default constructor
      */
      ItoKMCJSON();

      /*!
	@brief Destructor
      */
      virtual ~ItoKMCJSON() noexcept;

      /*!
	@brief Parse run-time options
      */
      virtual void
      parseRuntimeOptions() noexcept override;

      /*!
	@brief Get the neutral density at a position in space
	@param[in] a_pos Physical position
	@return Neutral density
      */
      virtual Real
      getNeutralDensity(const RealVect a_pos) const noexcept override;

      /*!
	@brief Compute Townsend ionization coefficient
	@param[in] a_E   Electric field. 
	@param[in] a_pos Physical coordinates
      */
      virtual Real
      computeAlpha(const Real a_E, const RealVect a_pos) const noexcept override;

      /*!
	@brief Compute Townsend attachment coefficient
	@param[in] a_E   Electric field. 
	@param[in] a_pos Physical coordinates
      */
      virtual Real
      computeEta(const Real a_E, const RealVect a_pos) const noexcept override;

      /*!
	@brief Get number of plot variables
      */
      virtual Vector<std::string>
      getPlotVariableNames() const noexcept override;

      /*!
	@brief Get plot variables
	@param[in] a_E       Electric field
	@param[in] a_pos     Physical position
	@param[in] a_phi     Plasma species densities
	@param[in] a_gradPhi Density gradients for plasma species. 
	@param[in] a_dx      Grid resolution
	@param[in] a_kappa   Cut-cell volume fraction
      */
      virtual Vector<Real>
      getPlotVariables(const RealVect          a_E,
                       const RealVect          a_pos,
                       const Vector<Real>&     a_phi,
                       const Vector<RealVect>& a_gradPhi,
                       const Real              a_dx,
                       const Real              a_kappa) const noexcept override;

      /*!
	@brief Get number of plot variables
      */
      virtual int
      getNumberOfPlotVariables() const noexcept override;

      /*!
	@brief Return true/false if physics model needs species gradients.
      */
      virtual bool
      needGradients() const noexcept override;

      /*!
	@brief Compute the Ito solver mobilities.
	@param[in] a_time Time
	@param[in] a_pos  Position
	@param[in] a_E    Electric field
      */
      virtual Vector<Real>
      computeMobilities(const Real a_time, const RealVect a_pos, const RealVect a_E) const noexcept override;

      /*!
	@brief Compute the Ito solver diffusion coefficients
	@param[in] a_time      Time
	@param[in] a_pos       Position
	@param[in] a_E         Electric field
      */
      virtual Vector<Real>
      computeDiffusionCoefficients(const Real a_time, const RealVect a_pos, const RealVect a_E) const noexcept override;

      /*!
	@brief Resolve secondary emission at the EB.
	@details Routine is here to handle charge injection, secondary emission etc.
	@param[out] a_secondaryParticles Outgoing plasma species particles.
	@param[out] a_secondaryCDRFluxes Outgoing plasma species CDR fluxes
	@param[out] a_secondaryPhotons   Photons injected through the EB
	@param[in]  a_primaryParticles   Particles that left the computational domain through the EB
	@param[in]  a_primaryCDRFluxes   CDR fluxes leaving the computational domain through the EB
	@param[in]  a_primaryPhotons     Photons that left the computational domain through the EB
	@param[in]  a_newNumParticles    Total number of particles in the cut-cell AFTER the transport step
	@param[in]  a_oldNumParticles    Total number of particles in the cut-cell BEFORE the transport step
	@param[in]  a_electricField      Electric field                  
	@param[in]  a_cellCenter         Physical position of the cell center. 
	@param[in]  a_cellCentroid       Cell centroid relative to the cell center (not multiplied by dx)
	@param[in]  a_bndryCentroid      EB face centroid relative to the cell center (not multiplied by dx)
	@param[in]  a_bndryNormal        Cut-cell normal vector. 
	@param[in]  a_bndryArea          Cut-cell boundary area - not multiplied by dx (2D) or dx^2 (3D)
	@param[in]  a_dx                 Grid resolution on this level. 
	@param[in]  a_dt                 Time step
	@param[in]  a_isDielectric       Dielectric or electrode. 
	@param[in]  a_matIndex           Material index (taken from computationalGeometry)
      */
      virtual void
      secondaryEmissionEB(Vector<List<ItoParticle>>&       a_secondaryParticles,
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
                          const int                        a_matIndex) const noexcept override;

    protected:
      /*!
	@brief Verbose or not
      */
      bool m_verbose;

      /*!
	@brief Plot gas pressure, density, and temperature
      */
      bool m_plotGas;

      /*!
	@brief Plot Townsend ionization coefficient
      */
      bool m_plotAlpha;

      /*!
	@brief Plot Townsend attachment coefficient
      */
      bool m_plotEta;

      /*!
	@brief A flag for skipping reactions completely. 
      */
      bool m_skipReactions;

      /*!
	@brief Preview rates or not
      */
      bool m_previewRates;

      /*!
	@brief Total number of plasma species
      */
      int m_numPlasmaSpecies;

      /*!
	@brief Total number of photon species
      */
      int m_numPhotonSpecies;

      /*!
	@brief Nohmanns implementation of JSON files
      */
      nlohmann::json m_json;

      /*!
	@brief Input JSON file name
      */
      std::string m_jsonFile;

      /*!
	@brief All species that have been defined.
	@details Only used to make sure we're not doubly defining some species
      */
      std::set<std::string> m_allSpecies;

      /*!
	@brief Map of string-int identifiers for background species
      */
      std::map<std::string, int> m_backgroundSpeciesMap;

      /*!
	@brief Map of int-string identifiers for background species
      */
      std::map<int, std::string> m_backgroundSpeciesMapInverse;

      /*!
	@brief Map of background species. 
	@details Built in initialization order but can look up arbitrary species name but using the maps. 
      */
      std::vector<ItoKMCBackgroundSpecies> m_backgroundSpecies;

      /*!
	@brief Plot background species or not
	@details Built in initialization order but can look up arbitrary species name but using the maps. 
      */
      std::vector<bool> m_backgroundSpeciesPlot;

      /*!
	@brief Map for identifying a species name with a solver type
      */
      std::map<std::string, SpeciesType> m_plasmaSpeciesTypes;

      /*!
	@brief Map for indexing from an ID to an Ito solver index
      */
      std::map<std::string, int> m_itoSpeciesMap;

      /*!
	@brief Map for indexing from an ID to an Ito solver index
      */
      std::map<std::string, int> m_cdrSpeciesMap;

      /*!
	@brief Map for indexing a plasma species name to a global index in the solver vectors
      */
      std::map<std::string, int> m_plasmaIndexMap;

      /*!
	@brief Map for indexing a photon species name to a global index in the solver vectors
      */
      std::map<std::string, int> m_photonIndexMap;

      /*!
	@brief Mobilities for the various species
      */
      std::vector<FunctionEX> m_mobilityFunctions;

      /*!
	@brief Diffusion coefficients for the various species
      */
      std::vector<FunctionEX> m_diffusionCoefficients;

      /*!
	@brief Plasma species temperatures
      */
      std::vector<FunctionEX> m_plasmaTemperatures;

      /*!
	@brief Function-based plasma reaction rates. Same index as the actual reactions.
      */
      std::vector<FunctionEVXTP> m_kmcReactionRates;

      /*!
	@brief Function-based plasma reaction rates as if they appeared in the reaction rate equation. Same index as actual reactions. 
	@details If S = k * n_1 * n_2 * n_3 where n_1, n_2, n_3 is any species, then this contains 'k'.
      */
      std::vector<FunctionEX> m_fluidRates;

      /*!
	@brief List of background species involved in a reaction
      */
      std::vector<std::list<size_t>> m_plasmaReactionBackgroundReactants;

      /*!
	@brief List of plasma species on the left-hand side of a reaction.
      */
      std::vector<std::list<size_t>> m_plasmaReactionPlasmaReactants;

      /*!
	@brief List of plasma species on the right-hand side of a reaction.
      */
      std::vector<std::list<size_t>> m_plasmaReactionPlasmaProducts;

      /*!
	@brief List of photon species on the right-hand side of a reaction.
      */
      std::vector<std::list<size_t>> m_plasmaReactionPhotonProducts;

      /*!
	@brief Whether or not a reaction uses a gradient correction for the rate (Soloviev et. al. paper)
      */
      std::vector<std::pair<bool, std::string>> m_kmcReactionGradientCorrections;

      /*!
	@brief Map of reactions that will be plotted.
      */
      std::vector<std::pair<bool, std::string>> m_kmcReactionRatePlots;

      /*!
	@brief Map of reactions that will be plotted.
      */
      std::map<int, std::pair<bool, int>> m_plasmaReactionSolovievCorrection;

      /*!
	@brief Flag if the user wants to compute the alpha coefficient from the reactions 
      */
      std::pair<bool, std::string> m_autoAlpha;

      /*!
	@brief Flag if the user wants to compute the alpha coefficient from the reactions 
      */
      std::pair<bool, std::string> m_autoEta;

      /*!
	@brief List of dielectric field emission reactions.
      */
      std::vector<std::pair<int, FunctionEN>> m_dielectricFieldEmission;

      /*!
	@brief List of electrode field emission reactions.
      */
      std::vector<std::pair<int, FunctionEN>> m_electrodeFieldEmission;

      /*!
	@brief Townsend ionization coefficient. Defined during parseAlpha
      */
      FunctionEX m_alpha;

      /*!
	@brief Townsend attachment coefficient. Defined during parseEta
      */
      FunctionEX m_eta;

      /*!
	@brief Background gas pressure
      */
      FunctionX m_gasPressure;

      /*!
	@brief Background gas temperature
      */
      FunctionX m_gasTemperature;

      /*!
	@brief Background gas number density
      */
      FunctionX m_gasNumberDensity;

      /*!
	@brief Trim a string. This removes whitespace before/after
	@param[in] a_string String to be trimmed
      */
      virtual std::string
      trim(const std::string& a_string) const noexcept;

      /*!
	@brief Parse the JSON file
      */
      virtual void
      parseJSON();

      /*!
	@brief Parse chattiness
      */
      virtual void
      parseVerbose() noexcept;

      /*!
	@brief Throw a parser error
	@param[in] a_error Error code.
      */
      virtual void
      throwParserError(const std::string a_error) const noexcept;

      /*!
	@brief Throw a parser wearning
	@param[in] a_error Warning
      */
      virtual void
      throwParserWarning(const std::string a_warning) const noexcept;

      /*!
	@brief Check if a file exists
	@param[in] a_filename File name
      */
      virtual bool
      doesFileExist(const std::string a_filename) const noexcept;

      /*!
	@brief Return true if species is a background species
	@param[in] a_name Species name
      */
      virtual bool
      isBackgroundSpecies(const std::string& a_name) const noexcept;

      /*!
	@brief Return true if species is a plasma species
	@param[in] a_name Species name
      */
      virtual bool
      isPlasmaSpecies(const std::string& a_name) const noexcept;

      /*!
	@brief Return true if species is a photon species
	@param[in] a_name Species name
      */
      virtual bool
      isPhotonSpecies(const std::string& a_name) const noexcept;

      /*!
	@brief Check if a string contains the wildcard @ and return true if it does
	@param[in] a_str Input string
      */
      virtual bool
      containsWildcard(const std::string a_str) const noexcept;

      /*!
	@brief Check if a string contains standard brackets
	@param[in] a_str Input string
      */
      virtual bool
      containsBracket(const std::string a_str) const noexcept;

      /*!
	@brief Check if a string contains starts with '(' and ends with ')'
	@param[in] a_str Input string
      */
      virtual bool
      isBracketed(const std::string a_str) const noexcept;

      /*!
	@brief Check that molar fraction is one. Throws a warning if it isn't.
	@param[in] a_position Physical coordinates
      */
      virtual void
      checkMolarFraction(const RealVect a_position) const noexcept;

      /*!
	@brief Check that a plasma reaction makes sense. I.e. that all species are properly defined and that it conserves charge.
	@param[in] a_reactants Reactants
	@param[in] a_products  Reaction products
	@param[in] a_reaction  Reaction string
      */
      virtual void
      sanctifyPlasmaReaction(const std::vector<std::string>& a_reactants,
                             const std::vector<std::string>& a_products,
                             const std::string&              a_reaction) const noexcept;

      /*!
	@brief Check that a photo-reaction makes sense. I.e. that all species are properly defined and that it conserves charge.
	@param[in] a_reactants Reactants
	@param[in] a_products  Reaction products
	@param[in] a_reaction  Reaction string
      */
      virtual void
      sanctifyPhotoReaction(const std::vector<std::string>& a_reactants,
                            const std::vector<std::string>& a_products,
                            const std::string&              a_reaction) const noexcept;

      /*!
	@brief Initialize gas law.
      */
      virtual void
      initializeGasLaw();

      /*!
	@brief Initialize the background species.
      */
      virtual void
      initializeBackgroundSpecies();

      /*!
	@brief Parse one of the Townsend coefficients. 
	@param[in] a_coeff Which coefficient. Must be "alpha" or "eta"
      */
      virtual void
      initializeTownsendCoefficient(const std::string a_coeff);

      /*!
	@brief Initialize automatic Townsend coefficient (ionizatino or attachment) from the reaction rates.
	@param[in] a_coeff Which coefficient. Must be "alpha" or "eta"
      */
      virtual void
      initializeAutomaticTownsend(const std::string a_coeff);

      /*!
	@brief Simple function for previewing a rate/coefficient defined as a function f = f(E,x)
	@param[in] a_json JSON entry containing the previewing information.
	@param[in] a_function Function to be previewed.
      */
      virtual void
      previewFunctionEX(const nlohmann::json& a_json, const FunctionEX& a_function) const;

      /*!
	@brief Print the fluid-representation of the reaction rates.
      */
      virtual void
      printFluidRates() const noexcept;

      /*!
	@brief Initialize the plasma species
      */
      virtual void
      initializePlasmaSpecies();

      /*!
	@brief Parse initial particles
      */
      virtual void
      initializeParticles();

      /*!
	@brief Parse initial densities for CDR and Ito species. 
      */
      virtual void
      initializeDensities();

      /*!
	@brief Initialize mobility functions
      */
      virtual void
      initializeMobilities();

      /*!
	@brief Initialize diffusion coefficients
      */
      virtual void
      initializeDiffusionCoefficients();

      /*!
	@brief Initialize the plasma species temperatures
      */
      virtual void
      initializeTemperatures();

      /*!
	@brief Initialize the photon species
      */
      virtual void
      initializePhotonSpecies();

      /*!
	@brief Initialize plasma reactions
      */
      virtual void
      initializePlasmaReactions();

      /*!
	@brief Initialize photo-reactions
      */
      virtual void
      initializePhotoReactions();

      /*!
	@brief Initialize surface reactions
	@param[in] a_surface Surface type
      */
      virtual void
      initializeSurfaceEmission(const std::string a_surface);

      /*!
	@brief Initialize field emission reactions
      */
      virtual void
      initializeFieldEmission();

      /*!
	@brief Initialize the particle placement algorithm
      */
      virtual void
      initializeParticlePlacement();

      /*!
	@brief Parses a reaction string into reactangs and products
	@param[out] a_reactants Left-hand side of reaction
	@param[out] a_products  Right-hand side of reaction
	@param[in]  a_reaction  Reaction string. Must be in format "a + b + c -> e + f + g". 
      */
      virtual void
      parseReactionString(std::vector<std::string>& a_reactants,
                          std::vector<std::string>& a_products,
                          const std::string&        a_reaction) const noexcept;

      /*!
	@brief Turn reactants and products strings into indices
	@param[out] a_backgroundReactants Background reactant indices
	@param[out] a_plasmaReactants Plasma reactant indices
	@param[out] a_photonReactants Photon reactant indices
	@param[out] a_backgroundProducts Background product indices
	@param[out] a_plasmaProducts Plasma product indices
	@param[out] a_photonProducts Photon product indices
	@param[in]  a_reactants Reactants 
	@param[in]  a_products Products
      */
      virtual void
      getReactionSpecies(std::list<size_t>&              a_backgroundReactants,
                         std::list<size_t>&              a_plasmaReactants,
                         std::list<size_t>&              a_photonReactants,
                         std::list<size_t>&              a_backgroundProducts,
                         std::list<size_t>&              a_plasmaProducts,
                         std::list<size_t>&              a_photonProducts,
                         const std::vector<std::string>& a_reactants,
                         const std::vector<std::string>& a_products) const noexcept;

      /*!
	@brief Parse a plasma reaction rate and turn it into a function
	@param[in] a_reactionJSON Reaction entry in JSON file
	@param[in] a_backgroundReactants List of background species appearing on the left hand side of the reaction
	@param[in] a_plasmaReactants List of plasma species appearing on the left-hand side of the reaction
	@return This should return a tuple of KMC and fluid rates. The KMC rate will also be multiplied by background species on the 
	left hand side of the reaction.
      */
      virtual std::pair<FunctionEVXTP, FunctionEX>
      parsePlasmaReactionRate(const nlohmann::json&    a_reactionJSON,
                              const std::list<size_t>& a_backgroundReactants,
                              const std::list<size_t>& a_plasmaReactants) const;

      /*!
	@brief Parse whether or not a reaction rate should be plotted
	@param[in] a_reactionJSON Reaction entry in JSON file
      */
      virtual std::pair<bool, std::string>
      parsePlasmaReactionPlot(const nlohmann::json& a_reactionJSON) const;

      /*!
	@brief Parse whether or not a reaction uses a gradient correction for the rate. 
	@param[in] a_reactionJSON Reaction entry in JSON file
      */
      virtual std::pair<bool, std::string>
      parsePlasmaReactionGradientCorrection(const nlohmann::json& a_reactionJSON) const;

      /*!
	@brief Parse a table which is stored in E/N format
	@param[in] a_tableEntry Table entry in JSON format. 
	@param[in] a_dataID     Data identifier
      */
      virtual LookupTable1D<Real, 1>
      parseTableEByN(const nlohmann::json& a_tableEntry, const std::string& a_dataID) const;

      /*!
	@brief Make a reaction set into a superset. This parses wildcards '@' in reaction string.
	@param[in] a_reactants List of reactants. Can contain wildcard.
	@param[in] a_products List of products. Can contain wildcard.
	@param[in] a_reactionJSON  JSON reaction entry. 
	@return Returns a list of reactants and products with wilcards replaced. The tuple entries are: <wildcard, reactants, products>
      */
      virtual std::vector<std::tuple<std::string, std::vector<std::string>, std::vector<std::string>>>
      parseReactionWildcards(const std::vector<std::string>& a_reactants,
                             const std::vector<std::string>& a_products,
                             const nlohmann::json&           a_reactionJSON) const noexcept;

      /*!
	@brief Update reaction rates
	@param[out] a_kmcReactions Reaction rates to be set. 
	@param[in] a_E       Electric field
	@param[in] a_pos     Physical position
	@param[in] a_phi     Plasma species densities
	@param[in] a_gradPhi Density gradients for plasma species.
	@param[in] a_dt      Time step
	@param[in] a_dx      Grid resolution
	@param[in] a_kappa   Cut-cell volume fraction
      */
      virtual void
      updateReactionRates(std::vector<std::shared_ptr<const KMCReaction>>& a_kmcReactions,
                          const RealVect                                   a_E,
                          const RealVect                                   a_pos,
                          const Vector<Real>&                              a_phi,
                          const Vector<RealVect>&                          a_gradPhi,
                          const Real                                       a_dt,
                          const Real                                       a_dx,
                          const Real                                       a_kappa) const noexcept override;

      /*!
	@brief Sample a multinomial distribution with N samples 
	@param[in] a_N Number of samples
	@param[in] a_distribution Discrete distribution
      */
      std::vector<size_t>
      multinomial(const size_t N, const std::discrete_distribution<size_t>& a_distribution) const noexcept;
    };
  } // namespace ItoKMC
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
