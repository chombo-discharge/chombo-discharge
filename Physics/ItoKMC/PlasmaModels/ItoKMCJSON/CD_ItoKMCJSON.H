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
#include <CD_LookupTable1D.H>
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
	@brief Compute a physics-based time step. 
	@param[in] a_E            Electric field. 
	@param[in] a_pos          Position
	@param[in] a_numParticles Number of particles per cell
      */
      virtual Real
      computeDt(const RealVect a_E, const RealVect a_pos, const Vector<Real> a_densities) const noexcept override;

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
	@brief Update reaction rates
	@param[in] a_E       Electric field
	@param[in] a_pos     Physical position
	@param[in] a_phi     Plasma species densities
	@param[in] a_gradPhi Density gradients for plasma species. 
	@param[in] a_dx      Grid resolution
	@param[in] a_kappa   Cut-cell volume fraction
      */
      virtual void
      updateReactionRates(const RealVect          a_E,
                          const RealVect          a_pos,
                          const Vector<Real>&     a_phi,
                          const Vector<RealVect>& a_gradPhi,
                          const Real              a_dx,
                          const Real              a_kappa) const noexcept override
      {}

      /*!
	@brief Compute the Ito solver mobilities.
	@param[in] a_time Time
	@param[in] a_pos  Position
	@param[in] a_E    Electric field
      */
      Vector<Real>
      computeMobilities(const Real a_time, const RealVect a_pos, const RealVect a_E) const noexcept override;

      /*!
	@brief Compute the Ito solver diffusion coefficients
	@param[in] a_time      Time
	@param[in] a_pos       Position
	@param[in] a_E         Electric field
      */
      Vector<Real>
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
                          const int                        a_matIndex) const noexcept override
      {}

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
	@brief Total number of plasma species
      */
      int m_numPlasmaSpecies;

      /*!
	@brief Nohmanns implementation of JSON files
      */
      nlohmann::json m_json;

      /*!
	@brief Input JSON file name
      */
      std::string m_jsonFile;

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
	@brief Map for indexing a species name to a global index in the solver vectors
      */
      std::map<std::string, int> m_solverIndexMap;

      /*!
	@brief Mobilities for the various species
      */
      std::vector<FunctionEX> m_mobilityFunctions;

      /*!
	@brief Diffusion coefficients for the various species
      */
      std::vector<FunctionEX> m_diffusionCoefficients;

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
	@brief Check if a file exists
	@param[in] a_filename File name
      */
      virtual bool
      doesFileExist(const std::string a_filename) const noexcept;

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
	@brief Check if a string contains the wildcard @ and return true if it does
      */
      virtual bool
      containsWildcard(const std::string a_str) const noexcept;

      /*!
	@brief Check that molar fraction is one. Throws a warning if it isn't.
	@param[in] a_position Physical coordinates
      */
      virtual void
      checkMolarFraction(const RealVect a_position) const noexcept;

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
      initializeTownsendCoefficient(const std::string a_str);

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
	@brief Parse a table which is stored in E/N format
	@param[in] a_tableEntry Table entry in JSON format. 
	@param[in] a_dataID     Data identifier
      */
      virtual LookupTable1D<2>
      parseTableEByN(const nlohmann::json& a_tableEntry, const std::string& a_dataID) const;
    };
  } // namespace ItoKMC
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
