/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Copyright © 2022 NTNU.
 * Copyright © 2022 Fanny Skirbekk. 
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaJSON.H
  @brief  Declaration of a generic CdrPlasma model.
  @author Robert Marskar, Fanny Skirbekk
  @todo   Energy stuff on domain sides is missing. 
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
#include <CD_CdrPlasmaPhotoReactionJSON.H>
#include <CD_CdrPlasmaSurfaceReactionJSON.H>
#include <CD_LookupTable.H>
#include <CD_NamespaceHeader.H>

using json = nlohmann::json;

namespace Physics {
  namespace CdrPlasma {

    /*!
      @brief Class which implements CdrPlasmaPhysics and parses plasma chemistry from a JSON input file
    */
    class CdrPlasmaJSON : public CdrPlasmaPhysics
    {
    public:
      /*!
	@brief Function alias for initial data function. 
	@param[in] a_position Physical coordinates
	@param[in] a_time     Time
      */
      using InitialDataFunction = std::function<Real(const RealVect a_position, const Real a_time)>;

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
	@brief Default constructor.
      */
      CdrPlasmaJSON();

      /*!
	@brief Destructor
      */
      virtual ~CdrPlasmaJSON();

      /*!
	@brief Parse run-time class options.
      */
      virtual void
      parseRuntimeOptions();

      /*!
	@brief Get number of plot variables for this physics class. 
	@details This is used by CdrPlasmaStepper for pre-allocating data that will be put in a plot file. The current implementation 
	plots the gas pressure, temperature, and density
      */
      virtual int
      getNumberOfPlotVariables() const override;

      /*!
	@brief Get plot variable names. The names and positions between this routine and getPlotVariables must be consistent!
      */
      virtual Vector<std::string>
      getPlotVariableNames() const override;

      /*!
	@brief Provide plot variables. This is used by CdrPlasmaStepper when writing plot files. 
	@details The plot variables should have length this->getNumberOfPlotVariables (as should getPlotVariableNames)
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
      virtual Vector<Real>
      getPlotVariables(const Vector<Real>     a_cdrDensities,
                       const Vector<RealVect> a_cdrGradients,
                       const Vector<Real>     a_rteDensities,
                       const RealVect         a_E,
                       const RealVect         a_position,
                       const Real             a_dx,
                       const Real             a_dt,
                       const Real             a_time,
                       const Real             a_kappa) const override;

      /*!
	@brief Compute alpha. Should return Townsend ionization coefficient. 
	@details This function is mostly used for the cell tagging classes, but can be used for reactions as well. 
	@param[in] a_E Electric field.
      */
      virtual Real
      computeAlpha(const Real E, const RealVect a_position) const override;

      /*!
	@brief Compute eta. Should return Townsend attachment coefficient. 
	@details This function is mostly used for the cell tagging classes, but can be used for reactions as well. 
	@param[in] a_E        Electric field.
	@param[in] a_position Position
      */
      virtual Real
      computeEta(const Real a_E, const RealVect a_position) const override;

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
      virtual void
      advanceReactionNetwork(Vector<Real>&          a_cdrSources,
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
	@return Returns the drift velocities for each CDR species. The vector ordering is the same as m_cdrSpecies. 
      */
      virtual Vector<RealVect>
      computeCdrDriftVelocities(const Real         a_time,
                                const RealVect     a_pos,
                                const RealVect     a_E,
                                const Vector<Real> a_cdrDensities) const override;

      /*!
	@brief Compute diffusion coefficients for the CDR equations. 
	@param[in] a_time         Time
	@param[in] a_pos          Position
	@param[in] a_E            Electric field
	@param[in] a_cdrDensities CDR densities
	@return Returns the diffusion coefficients for each CDR species. The vector ordering is the same as m_cdrSpecies. 
      */
      virtual Vector<Real>
      computeCdrDiffusionCoefficients(const Real         a_time,
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
	@param[in] a_rteFluxes     RTE fluxes (normal component only)
	@param[in] a_extrapCdrFluxes Extrapolated fluxes from the gas side. 
	@return Returns the flux on an electrode interface cell. The vector ordering must be the same as m_cdrSpecies. 
      */
      virtual Vector<Real>
      computeCdrElectrodeFluxes(const Real         a_time,
                                const RealVect     a_pos,
                                const RealVect     a_normal,
                                const RealVect     a_E,
                                const Vector<Real> a_cdrDensities,
                                const Vector<Real> a_cdrVelocities,
                                const Vector<Real> a_cdrGradients,
                                const Vector<Real> a_rteFluxes,
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
	@param[in] a_rteFluxes     RTE fluxes (normal component only)
	@param[in] a_extrapCdrFluxes Extrapolated fluxes from the gas side. 
	@return Returns the flux on a dielectric interface cell. The vector ordering must be the same as m_cdrSpecies. 
      */
      virtual Vector<Real>
      computeCdrDielectricFluxes(const Real         a_time,
                                 const RealVect     a_pos,
                                 const RealVect     a_normal,
                                 const RealVect     a_E,
                                 const Vector<Real> a_cdrDensities,
                                 const Vector<Real> a_cdrVelocities,
                                 const Vector<Real> a_cdrGradients,
                                 const Vector<Real> a_rteFluxes,
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
	@param[in] a_rteFluxes     RTE fluxes (normal component only)
	@param[in] a_extrapCdrFluxes Extrapolated fluxes from the gas side. 
      */
      virtual Vector<Real>
      computeCdrDomainFluxes(const Real           a_time,
                             const RealVect       a_pos,
                             const int            a_dir,
                             const Side::LoHiSide a_side,
                             const RealVect       a_E,
                             const Vector<Real>   a_cdrDensities,
                             const Vector<Real>   a_cdrVelocities,
                             const Vector<Real>   a_cdrGradients,
                             const Vector<Real>   a_rteFluxes,
                             const Vector<Real>   a_extrapCdrFluxes) const override;

      /*!
	@brief Set the initial surface charge
	@param[in] a_time Time
	@param[in] a_pos Position
      */
      virtual Real
      initialSigma(const Real a_time, const RealVect a_pos) const override;

    protected:
      /*!
	@brief Enum class for distinguishing types of computation methods when computing transport data stuff. 
      */
      enum class LookupMethod
      {
        Constant,
        FunctionX,
        FunctionT,
        FunctionTT,
        FunctionEN,
        FunctionEX,
        TableEN,
        TableEnergy,
        AlphaV,
        EtaV
      };

      /*!
	@brief Enum for distinguishing integration methods. Note that 'None' just fills directly with source terms. 
      */
      enum class ReactionIntegrator
      {
        None,
        ExplicitEuler,
        ExplicitTrapezoidal,
        ExplicitMidpoint,
        ExplicitRK4
      };

      /*!
	@brief Enum class for distinguishing how we add/lose energy when running LEA-based models. The user will
	specify how to do this in the chemistry file, where 'AddMean' means that we add the mean energy to the equation.
	Likewise, 'SubtractMean' implies removing the mean energy, and 'External' is a user-specified energy to be added/removed.
      */
      enum class ReactiveEnergyLoss
      {
        AddMean,
        SubtractMean,
        AddDirect,
        SubtractDirect,
        External
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
	@brief Plot Townsend ionization coefficient
      */
      bool m_plotAlpha;

      /*!
	@brief Plot Townsend attachment coefficient
      */
      bool m_plotEta;

      /*!
	@brief Using discrete photons or not.
      */
      bool m_discretePhotons;

      /*!
	@brief A flag for skipping reactions completely. 
      */
      bool m_skipReactions;

      /*!
	@brief JSON definition. This is populated when calling parseJSON.
      */
      json m_json;

      /*!
	@brief Reaction integrator
      */
      ReactionIntegrator m_reactionIntegrator;

      /*!
	@brief Input JSON file name
      */
      std::string m_jsonFile;

      /*!
	@brief JSON entries for species in the defined field 'plasma species'. 
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
	@brief Gas pressure (in Pascal).
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
	@brief Chemistry time step. 
      */
      Real m_chemistryDt;

      /*!
	@brief Neutral species densities
      */
      std::vector<FunctionX> m_neutralSpeciesDensities;

      /*!
	@brief These are the neutral species.
      */
      std::vector<std::shared_ptr<NeutralSpeciesJSON>> m_neutralSpecies;

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
	@brief int-bool encoding for determining if a solver is an energy solver. 
      */
      std::map<int, bool> m_cdrIsEnergySolver;

      /*!
	@brief int-bool encoding for determining if a CDR solver HAS an associated energy solver. 
      */
      std::map<int, bool> m_cdrHasEnergySolver;

      /*!
	@brief Parameters for computing the mean energy from energy density and density. 
	@details The energy is computed as e = max(Emin, min(Emax, E)) where E = n_energy/(min(n_density, safety));
      */
      std::map<int, std::tuple<Real, Real, Real>> m_cdrEnergyComputation;

      /*!
	@brief int-int encoding for associating a transport solver with an energy solver. 
	@details This is only defined for the species that has/is an energy solver. The first index is the transport solver index
	and the second index is the energy solver index. 
      */
      std::map<int, int> m_cdrTransportEnergyMap;

      /*!
	@brief Map of the species masses. This is needed for imposing BCs on the energy equations.
	@details This is populated in initializePlasmaSpecies.
      */
      std::map<int, Real> m_cdrMasses;

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

      // =====================================================
      // TOWNSEND IONIZATION COEFFICIENT QUANTITIES BEGIN HERE
      // =====================================================

      /*!
	@brief Lookup method for Townsend ionization coefficient
      */
      LookupMethod m_alphaLookup;

      /*!
	@brief Lookup method for Townsend attachment coefficient
      */
      LookupMethod m_etaLookup;

      /*!
	@brief For when we can use alpha = constant
      */
      Real m_alphaConstant;

      /*!
	@brief For when we can use eta = constant
      */
      Real m_etaConstant;

      /*!
	@brief For when we can put alpha = alpha(E,N) as an analytic function. 
      */
      FunctionEN m_alphaFunctionEN;

      /*!
	@brief For when we can put eta = eta(E,N) as an analytic function. 
      */
      FunctionEN m_etaFunctionEN;

      /*!
	@brief For when we can put alpha = table(E,N)
      */
      LookupTable1D<Real, 1> m_alphaTableEN;

      /*!
	@brief For when we can put eta = table(E,N)
      */
      LookupTable1D<Real, 1> m_etaTableEN;

      // ================================
      // MOBILITY QUANTITIES BEGIN HERE
      // ================================

      /*!
	@brief Mobility lookup method for each species. 
      */
      std::map<int, LookupMethod> m_mobilityLookup;

      /*!
	@brief Map for constant mobilities. 
      */
      std::map<int, Real> m_mobilityConstants;

      /*!
	@brief Map for function-based mobilities mu = mu(E,N)
      */
      std::map<int, FunctionEN> m_mobilityFunctionsEN;

      /*!
	@brief Map for function-based mobilities mu = mu(E,x)
      */
      std::map<int, FunctionEX> m_mobilityFunctionsEX;

      /*!
	@brief Map for table-based mobilities. Stored as tables (E/N, mu*N)
      */
      std::map<int, LookupTable1D<Real, 1>> m_mobilityTablesEN;

      /*!
	@brief Map for table-based mobilities as function of energy. 
	@details Stored as tables (eV, mu*N)
      */
      std::map<int, LookupTable1D<Real, 1>> m_mobilityTablesEnergy;

      // ================================
      // DIFFUSION  QUANTITIES BEGIN HERE
      // ================================

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
	@brief Map for table-based diffusion coefficients D = D(E,N).
	@details Tables stored as (E/N, D*N)
      */
      std::map<int, LookupTable1D<Real, 1>> m_diffusionTablesEN;

      /*!
	@brief Map for table-based diffusion coefficients as function of energy. 
	@details Stored as tables (eV, D*N)
      */
      std::map<int, LookupTable1D<Real, 1>> m_diffusionTablesEnergy;

      // ===========================================
      // PLASMA SPECIES TEMPERATURE DATA BEGINS HERE
      // ===========================================

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
      std::map<int, LookupTable1D<Real, 1>> m_temperatureTablesEN;

      // ================================
      // REACTION DATA BEGINS HERE
      // ================================

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
	@brief Plasma reaction rates that are alpha*|v|
      */
      std::map<int, int> m_plasmaReactionAlphaV;

      /*!
	@brief Plasma reaction rates that are eta*|v|
      */
      std::map<int, int> m_plasmaReactionEtaV;
      /*!
	@brief Maps for functions of the type k = f(T) where T is the temperature of some species.
	@details In the tuple the first index in the pair indicates the species we are talking about here. A special rule is enforced
	when this index is < 0 in which case the temperature is replaced by the background gas temperature. 
      */
      std::map<int, std::pair<int, FunctionT>> m_plasmaReactionFunctionsT;

      /*!
	@brief Maps for functions of the type k = f(T1,T2) where T1 and T2 are the temperatures of some species. 
	@details This signature is absolute horrific -- what it means is that we have a reaction which should be evaluated
	as f(T1, T2), but we need to know which species are involved. By design, this should be general so that T1 and T2 
	can be the temperatures of any species, including neutral species. So, we make a tuple for indicating which species
	we are talking about. The first index in the tuple is the first species, the second is the second species and the third
	entry in the tuple is the actual function. A special rule occurs if one of the first two indices is < 0 in which case
	the temperature is replaced by the background gas temperature. 
      */
      std::map<int, std::tuple<int, int, FunctionTT>> m_plasmaReactionFunctionsTT;

      /*!
	@brief Function-based plasma reaction rates. 
      */
      std::map<int, FunctionEN> m_plasmaReactionFunctionsEN;

      /*!
	@brief Map for table-based reaction coefficients, where k = k(E,N).
      */
      std::map<int, LookupTable1D<Real, 1>> m_plasmaReactionTablesEN;

      /*!
	@brief Map for table-based reaction coefficients where k = k(energy).
	@details The first index is the reaction index, while the pair describes which species energy and the tabulated data. 
      */
      std::map<int, std::pair<int, LookupTable1D<Real, 1>>> m_plasmaReactionTablesEnergy;

      /*!
	@brief Scaled plasma reactions. These account for e.g. reaction efficiencies, collisional quenching, etc. 
      */
      std::map<int, FunctionEX> m_plasmaReactionEfficiencies;

      /*!
	@brief Plasma reactions. 
      */
      std::vector<CdrPlasmaReactionJSON> m_plasmaReactions;

      /*!
	@brief Plot plasma reaction or not
      */
      std::map<int, bool> m_plasmaReactionPlot;

      /*!
	@brief Flag for whether or not reaction includes Soloviev energy correction. 
	@details If this is true, the rate for a reaction (in the local field approximation) will be modified as k * (1 + E.(D * grad(phi))/n * 
      */
      std::map<int, std::pair<bool, int>> m_plasmaReactionSolovievCorrection;

      /*!
	@brief For mapping reactive energy losses for all reactions. 
	@details The index in the first map is the reaction (i.e., index in m_plasmaReactions).
	The second map determines how reactive energy losses/gains are computed for each species. The first index in the second map
	indicates the species which will add/lose energy, and the std::pair indicates how this loss is computed. 
      */
      std::map<int, std::map<int, std::pair<ReactiveEnergyLoss, Real>>> m_plasmaReactionEnergyLosses;

      /*!
	@brief Associative container for determining if a reaction is associated with an energy loss/gain.
      */
      std::map<int, bool> m_plasmaReactionHasEnergyLoss;

      // ================================
      // PHOTO-REACTIONS BEGIN HERE
      // ================================

      /*!
	@brief Flag for photo-reaction efficiencies. Includes Helmholtz corrections, if present. 
      */
      std::map<int, FunctionEX> m_photoReactionEfficiencies;

      /*!
	@brief Map over the Helmholtz reconstructions. 
      */
      std::map<int, bool> m_photoReactionUseHelmholtz;

      /*!
	@brief Photo-reactions
      */
      std::vector<CdrPlasmaPhotoReactionJSON> m_photoReactions;

      /*!
	@brief Associated energy losses for a photo-reaction.
	@details The list is a list of species and the corresponding energy losses for a reaction. Note that although 
	we write 'loss', this can also be an energy 'gain'. 
      */
      std::map<int, std::list<std::pair<int, Real>>> m_photoReactionEnergyLosses;

      /*!
	@brief Associative container for determining if a reaction is associated with an energy loss/gain.
      */
      std::map<int, bool> m_photoReactionHasEnergyLoss;

      // ===========================================
      // ELECTRODE SURFACE REACTIONS BEGIN HERE
      // ===========================================

      /*!
	@brief Lookup method for the electrode surface reaction rates
      */
      std::map<int, LookupMethod> m_electrodeReactionLookup;

      /*!
	@brief Constant electrode reaction rate. 
      */
      std::map<int, Real> m_electrodeReactionConstants;

      /*!
	@brief Electrode reaction effiencies. Used for scaling reactions on electrodes in a "generic" way. 
      */
      std::map<int, FunctionEX> m_electrodeReactionEfficiencies;

      /*!
	@brief List of electrode reactions
      */
      std::vector<CdrPlasmaSurfaceReactionJSON> m_electrodeReactions;

      /*!
	@brief Associated energy losses for a surface reaction on electrodes. 
	@details The list is a list of species and the corresponding energy losses for a reaction. Note that although 
	we write 'loss', this can also be an energy 'gain'. 
      */
      std::map<int, std::list<std::pair<int, Real>>> m_electrodeReactionEnergyLosses;

      /*!
	@brief Associative container for determining if a reaction is associated with an energy loss/gain.
      */
      std::map<int, bool> m_electrodeReactionHasEnergyLoss;

      /*!
	@brief A container which determines if we should add the extrapolated flux as an inflow condition. 
      */
      std::map<int, bool> m_electrodeExtrapBC;

      // ===========================================
      // DIELECTRIC SURFACE REACTIONS BEGIN HERE
      // ===========================================

      /*!
	@brief Lookup method for the dielectric surface reaction rates
      */
      std::map<int, LookupMethod> m_dielectricReactionLookup;

      /*!
	@brief Constant dielectric reaction rate. 
      */
      std::map<int, Real> m_dielectricReactionConstants;

      /*!
	@brief Dielectric reaction effiencies. Used for scaling reactions on dielectrics in a "generic" way. 
      */
      std::map<int, FunctionEX> m_dielectricReactionEfficiencies;

      /*!
	@brief List of dielectric reactions
      */
      std::vector<CdrPlasmaSurfaceReactionJSON> m_dielectricReactions;

      /*!
	@brief Associated energy losses for a surface reaction on dielectrics. 
	@details The list is a list of species and the corresponding energy losses for a reaction. Note that although 
	we write 'loss', this can also be an energy 'gain'. 
      */
      std::map<int, std::list<std::pair<int, Real>>> m_dielectricReactionEnergyLosses;

      /*!
	@brief Associative container for determining if a reaction is associated with an energy loss/gain.
      */
      std::map<int, bool> m_dielectricReactionHasEnergyLoss;

      /*!
	@brief A container which determines if we should add the extrapolated flux as an inflow condition. 
      */
      std::map<int, bool> m_dielectricExtrapBC;

      // ===========================================
      // DOMAIN SURFACE REACTIONS BEGIN HERE
      // ===========================================

      /*!
	@brief Lookup method for the domain reaction rates. The pair is made up of an int representing direction (0=x, 1=y, 2=z) and a Side::LoHiSide
	representing side (Side::Lo, Side::Hi)
      */
      std::map<std::pair<int, Side::LoHiSide>, std::map<int, LookupMethod>> m_domainReactionLookup;

      /*!
	@brief Constant domain reaction rate. The pair is made up of an int representing direction (0=x, 1=y, 2=z) and a Side::LoHiSide
	representing side (Side::Lo, Side::Hi)
      */
      std::map<std::pair<int, Side::LoHiSide>, std::map<int, Real>> m_domainReactionConstants;

      /*!
	@brief Domain reaction effiencies. Used for scaling reactions on domains in a "generic" way. 
	The pair is made up of an int representing direction (0=x, 1=y, 2=z) and a Side::LoHiSide
	representing side (Side::Lo, Side::Hi)
      */
      std::map<std::pair<int, Side::LoHiSide>, std::map<int, FunctionEX>> m_domainReactionEfficiencies;

      /*! 
	@brief List of domain reactions. The pair is made up of an int representing direction (0=x, 1=y, 2=z) and a Side::LoHiSide
	representing side (Side::Lo, Side::Hi)
      */
      std::map<std::pair<int, Side::LoHiSide>, std::vector<CdrPlasmaSurfaceReactionJSON>> m_domainReactions;

      /*! 
	@brief map to translate dir from char to int
      */
      const std::map<char, int> m_dirCharToInt{{'x', 0}, {'y', 1}, {'z', 2}};

      /*! 
	@brief map to translate side from std::string to Side::LoHiSide
      */
      const std::map<std::string, Side::LoHiSide> m_sideStringToSide{{"lo", Side::Lo}, {"hi", Side::Hi}};

      /*!
	@brief A container which determines if we should add the extrapolated flux as an inflow condition. 
	@details The first entry in the map is the species integer index. The tuple indicates (coordinate_direction, side, doExtrap).
      */
      std::map<std::tuple<int, Side::LoHiSide, int>, bool> m_domainExtrapBC;

      // ===========================================
      // VARIOUS MEMBER FUNCTIONS BEGIN HERE
      // ===========================================

      /*!
	@brief Parse class options
      */
      virtual void
      parseOptions();

      /*!
	@brief Parse the reactive integrator
      */
      virtual void
      parseIntegrator();

      /*!
	@brief Parse the JSON file
      */
      virtual void
      parseJSON();

      /*!
	@brief Initialize surface charge
	@details This will initialize the surface charge based on the "sigma" field in the JSON file. If you want more complex initial surface charges, 
	you will have to either extend this routine or overwrite it. 
      */
      virtual void
      initializeSigma();

      /*!
	@brief Initialize neutral species
      */
      virtual void
      initializeNeutralSpecies();

      /*!
	@brief Initialize species. 
	@details This will initialize the species based on the "plasma_species" field in the JSON file. 
      */
      virtual void
      initializePlasmaSpecies();

      /*!
	@brief Initialize photon species.
	@details This will initialize the radiative transfer species absed on the "photon_species" field in the JSON file. 
      */
      virtual void
      initializePhotonSpecies();

      /*!
	@brief Initialize species mobilities.
      */
      virtual void
      parseMobilities();

      /*!
	@brief Initialize species diffusion coefficients.
      */
      virtual void
      parseDiffusion();

      /*!
	@brief Initialize species temperatures.
      */
      virtual void
      parseTemperatures();

      /*!
	@brief Parse the Townsend ionization coefficient
      */
      virtual void
      parseAlpha();

      /*!
	@brief Parse the Townsend attachment coefficient
      */
      virtual void
      parseEta();

      /*!
	@brief Parse plasma reactions
      */
      virtual void
      parsePlasmaReactions();

      /*!
	@brief Generate an initial data function for a given plasma species
	@param[in] a_json JSON field, usually (always?) describing one of the objects in the 'plasma species' field.
      */
      virtual InitialDataFunction
      parsePlasmaSpeciesInitialData(const json& a_json) const;

      /*!
	@brief Generate initial particles for a given plasma species
	@param[in] a_json JSON field, usually (always?) describing one of the objects in the 'plasma species' field.
      */
      virtual List<PointParticle>
      parsePlasmaSpeciesInitialParticles(const json& a_json) const;

      /*!
	@brief Make a reaction set into a superset. This parses wildcards '@' in reaction string.
	@param[in] a_reactants List of reactants. Can contain wildcard.
	@param[in] a_products  List of products. Can contain wildcard.
	@param[in] a_reaction  JSON reaction entry. 
	@return Returns a list of reactants and products with wilcards replaced. The tuple entries are: <wildcard, reactants, products>
      */
      virtual std::list<std::tuple<std::string, std::vector<std::string>, std::vector<std::string>>>
      parseReactionWildcards(const std::vector<std::string>& a_reactants,
                             const std::vector<std::string>& a_products,
                             const json&                     a_reaction);
      /*!
	@brief Parses a reaction string into reactangs and products
	@param[out] a_reactants Left-hand side of reaction
	@param[out] a_products  Right-hand side of reaction
	@param[in]  a_reaction  Reaction string. Must be in format "a + b + c -> e + f + g". 
      */
      virtual void
      parseReactionString(std::vector<std::string>& a_reactants,
                          std::vector<std::string>& a_products,
                          const std::string&        a_reaction) const;

      /*!
	@brief Check if a plasma-reaction makes sense in terms of the species that have been defined.
	@details This will throw errors if the species do not exist or charge is not conserved.
	@param[in] a_reactants Reactants
	@param[in] a_products  Reaction products
	@param[in] a_reaction  Reaction string
      */
      virtual void
      sanctifyPlasmaReaction(const std::vector<std::string>& a_reactants,
                             const std::vector<std::string>& a_products,
                             const std::string               a_reaction) const;

      /*!
	@brief Get the int-encoding corresponding to species involved in some reaction. 
	@param[out] a_plasmaReactants   Plasma reactants
	@param[out] a_neutralReactants  Neutral reactants
	@param[out] a_photonReactants   Photon reactants
	@param[out] a_plasmaProducts    Plasma products
	@param[out] a_neutralProducts   Neutral products
	@param[out] a_photonProducts    Photon products
	@param[in]  a_reactants         Reactant names
	@param[in]  a_products          Reaction names
      */
      virtual void
      getReactionSpecies(std::list<int>&                 a_plasmaReactants,
                         std::list<int>&                 a_neutralReactants,
                         std::list<int>&                 a_photonReactants,
                         std::list<int>&                 a_plasmaProducts,
                         std::list<int>&                 a_neutralProducts,
                         std::list<int>&                 a_photonProducts,
                         const std::vector<std::string>& a_reactants,
                         const std::vector<std::string>& a_products) const;

      /*!
	@brief Parse reaction rate for plasma reaction.
	@param[in] a_reactionIndex  Reaction index. 
	@param[in] a_reactionJSON   Input reaction. Must be one of the entries in the 'plasma reactions' array. 
      */
      virtual void
      parsePlasmaReactionRate(const int a_reactionIndex, const json& a_reactionJSON);

      /*!
	@brief Parse scaling factors for reactions. 
	@param[in] a_reactionIndex  Reaction index. 
	@param[in] a_reactionJSON   Input reaction. Must be one of the entries in the 'plasma reactions' array. 
      */
      virtual void
      parsePlasmaReactionScaling(const int a_reactionIndex, const json& a_reactionJSON);

      /*!
	@brief Parse reaction plotting. 
	@param[in] a_reactionIndex  Reaction index. 
	@param[in] a_reactionJSON   Input reaction. Must be one of the entries in the 'plasma reactions' array. 
      */
      virtual void
      parsePlasmaReactionPlot(const int a_reactionIndex, const json& a_reactionJSON);

      /*!
	@brief Parse plasma reaction descriptions. 
	@param[in] a_reactionIndex  Reaction index. 
	@param[in] a_reactionJSON   Input reaction. Must be one of the entries in the 'plasma reactions' array. 
      */
      virtual void
      parsePlasmaReactionDescription(const int         a_reactionIndex,
                                     const json&       a_reactionJSON,
                                     const std::string a_wildcard);

      /*!
	@brief Parse plasma reaction energy correction. 
	@param[in] a_reactionIndex  Reaction index. 
	@param[in] a_reactionJSON   Input reaction. Must be one of the entries in the 'plasma reactions' array. 
      */
      virtual void
      parsePlasmaReactionSoloviev(const int a_reactionIndex, const json& a_reactionJSON);

      /*!
	@brief Parse plasma reaction energy losses.
	@param[in] a_reactionIndex  Reaction index. 
	@param[in] a_reactionJSON   Input reaction. Must be one of the entries in the 'plasma reactions' array. 
      */
      virtual void
      parsePlasmaReactionEnergyLosses(const int a_reactionIndex, const json& a_reactionJSON);

      /*!
	@brief Parse photo-reactions
      */
      virtual void
      parsePhotoReactions();

      /*!
	@brief Parse scaling for photo-reactions. Includes Helmholtz corrections if doing Helmholtz reconstruction of photoionization profiles. 
	@param[in] a_reactionIndex  Reaction index. 
	@param[in] a_reactionJSON   Input reaction. Must be one of the entries in the 'photo reactions' array. 
      */
      virtual void
      parsePhotoReactionScaling(const int a_reactionIndex, const json& a_reactionJSON);

      /*!
	@brief Parse photo-reaction energy losses.
	@param[in] a_reactionIndex  Reaction index. 
	@param[in] a_reactionJSON   Input reaction. Must be one of the entries in the 'photo reactions' array. 
      */
      virtual void
      parsePhotoReactionEnergyLosses(const int a_reactionIndex, const json& a_reactionJSON);

      /*!
	@brief Do a species sanity check
	@details This will make sure that species names are not duplicates. 
      */
      virtual void
      sanityCheckSpecies() const;

      /*!
	@brief Check if a photo-reaction makes sense in terms of the species that have been defined.
	@details This will throw errors if the species do not exist or charge is not conserved.
	@param[in] a_reactants Reactants
	@param[in] a_products  Reaction products
	@param[in] a_reaction  Reaction string
      */
      virtual void
      sanctifyPhotoReaction(const std::vector<std::string>& a_reactants,
                            const std::vector<std::string>& a_products,
                            const std::string               a_reaction) const;

      /*!
	@brief Parse secondary emission on electrodes
      */
      virtual void
      parseElectrodeReactions();

      /*!
	@brief Parse reaction rate for electrode surface reactions. 
	@param[in] a_reactionIndex  Reaction index. 
	@param[in] a_reactionJSON   Input reaction. Must be one of the entries in the 'electrode reactions' array. 
      */
      virtual void
      parseElectrodeReactionRate(const int a_reactionIndex, const json& a_reactionJSON);

      /*!
	@brief Parse electrode reaction scaling for a specific reaction.
	@details This can be used to adjust for e.g. positional or field dependence.
	@param[in] a_reactionIndex Reaction index (this is the array index in the reaction array supplied in the JSON file)
	@param[in] a_reactionJSON  JSON entry for the reaction in the JSON input file.
      */
      virtual void
      parseElectrodeReactionScaling(const int a_reactionIndex, const json& a_reactionJSON);

      /*!
	@brief Parse electrode-reaction energy losses.
	@param[in] a_reactionIndex  Reaction index. 
	@param[in] a_reactionJSON   Input reaction. Must be one of the entries in the 'electrode reactions' array. 
      */
      virtual void
      parseElectrodeReactionEnergyLosses(const int a_reactionIndex, const json& a_reactionJSON);

      /*!
	@brief Parse secondary emission on dielectrics
      */
      virtual void
      parseDielectricReactions();

      /*!
	@brief Parse reaction rate for dielectric surface reactions. 
	@param[in] a_reactionIndex  Reaction index. 
	@param[in] a_reactionJSON   Input reaction. Must be one of the entries in the 'electrode reactions' array. 
      */
      virtual void
      parseDielectricReactionRate(const int a_reactionIndex, const json& a_reactionJSON);

      /*!
	@brief Parse dielectric electrode reaction scaling for a specific reaction.
	@details This can be used to adjust for e.g. positional or field dependence.
	@param[in] a_reactionIndex Reaction index (this is the array index in the dielectric reaction array supplied in the JSON file)
	@param[in] a_reactionJSON  JSON entry for the reaction in the JSON input file.
      */
      virtual void
      parseDielectricReactionScaling(const int a_reactionIndex, const json& a_reactionJSON);

      /*!
	@brief Parse dielectric-reaction energy losses.
	@param[in] a_reactionIndex  Reaction index. 
	@param[in] a_reactionJSON   Input reaction. Must be one of the entries in the 'dielectric reactions' array. 
      */
      virtual void
      parseDielectricReactionEnergyLosses(const int a_reactionIndex, const json& a_reactionJSON);

      /*!
	@brief Parse secondary emission on domain
      */
      virtual void
      parseDomainReactions();

      /*!
	@brief Parse reaction rate for domain reactions. 
	@param[in] a_reactionIndex  Reaction index. 
	@param[in] a_reactionJSON   Input reaction. Must be one of the entries in the 'domain reactions' array.
	@param[in] a_sides          List of the domain sides that this occurs at
      */
      virtual void
      parseDomainReactionRate(const int                       a_reactionIndex,
                              const json&                     a_reactionJSON,
                              const std::vector<std::string>& a_sides);

      /*!
	@brief Parse domain reaction scaling for a specific reaction.
	@details This can be used to adjust for e.g. positional or field dependence.
	@param[in] a_reactionIndex Reaction index (this is the array index in the domain reaction array supplied in the JSON file)
	@param[in] a_reactionJSON  JSON entry for the reaction in the JSON input file.
	@param[in] a_sides          List of the domain sides that this occurs at
      */
      virtual void
      parseDomainReactionScaling(const int                       a_reactionIndex,
                                 const json&                     a_reactionJSON,
                                 const std::vector<std::string>& a_sides);

      /*!
	@brief Check if a surface-reaction makes sense in terms of the species that have been defined.
	@param[in] a_reactants Reactants
	@param[in] a_products  Reaction products
	@param[in] a_reaction  Reaction string
      */
      virtual void
      sanctifySurfaceReaction(const std::vector<std::string>& a_reactants,
                              const std::vector<std::string>& a_products,
                              const std::string               a_reaction) const;

      /*!
	@brief Compute the various plasma species temperatures.
	@param[in] a_position     Physical coordinates
	@param[in] a_E            Electric field (SI units)
	@param[in] a_cdrDensities List of plasma species densities.
      */
      virtual std::vector<Real>
      computePlasmaSpeciesTemperatures(const RealVect&          a_position,
                                       const RealVect&          a_E,
                                       const std::vector<Real>& a_cdrDensities) const;

      /*!
	@brief Compute the various plasma species energies. Returns a list of energies in electron-volts. We assume that 
	temperature-energy relations are e = 3/2 * kB * T. 
	@details If a transported species is associated with an energy solver the energy is computed from e = (n_eps)/ne. Otherwise, we compute
	the energy and returns 3/2 * kB * T.
	@param[in] a_position     Physical coordinates
	@param[in] a_E            Electric field (SI units)
	@param[in] a_cdrDensities List of plasma species densities.
      */
      virtual std::vector<Real>
      computePlasmaSpeciesEnergies(const RealVect&          a_position,
                                   const RealVect&          a_E,
                                   const std::vector<Real>& a_cdrDensities) const;

      /*!
	@brief Compute the various plasma species mobilities. 
	@param[in] a_position     Physical coordinates
	@param[in] a_E            Electric field (SI units)
	@param[in] a_cdrDensities List of plasma species densities.
      */
      virtual std::vector<Real>
      computePlasmaSpeciesMobilities(const RealVect&          a_position,
                                     const RealVect&          a_E,
                                     const std::vector<Real>& a_cdrDensities) const;

      /*!
	@brief Compute the various plasma species diffusion coefficients. 
	@param[in] a_position     Physical coordinates
	@param[in] a_E            Electric field (SI units)
	@param[in] a_cdrDensities List of plasma species densities
	@return Returns the diffusion coefficients for each CDR species. The vector ordering is the same as m_cdrSpecies. 
      */
      virtual std::vector<Real>
      computePlasmaSpeciesDiffusion(const RealVect          a_position,
                                    const RealVect          a_E,
                                    const std::vector<Real> a_cdrDensities) const;

      /*!
	@brief Compute the reaction rate for a plasma reaction.
	@details This routine exists because we need to compute the rates both in advanceReactionNetwork and getPlotVariables. This function
	reduces code duplication. 
	@param[in] a_reactionIndex            Reaction index
	@param[in] a_cdrDensities             Plasma species densities. 
	@param[in] a_cdrMobilities            Plasma species mobilities. 
	@param[in] a_cdrDiffusionCoefficients Plasma species diffusion coefficients. 
	@param[in] a_cdrTemperatures          Plasma species temperatures. 
	@param[in] a_cdrEnergies              Plasma species energies.
	@param[in] a_cdrGradients             Plasma species gradients. 
	@param[in] a_pos                      Position (physical coordinates)
	@param[in] a_vectorE                  Electric field (vector)
	@param[in] a_E                        Electric field magnitude (SI units)
	@param[in] a_Etd                      Electric field magnitude (Townsend units)
	@param[in] a_N                        Neutral density
	@param[in] a_alpha                    Townsend ionization coefficient
	@param[in] a_eta                      Townsend attachment coefficient
	@param[in] a_time                     Time
      */
      virtual Real
      computePlasmaReactionRate(const int&                   a_reactionIndex,
                                const std::vector<Real>&     a_cdrDensities,
                                const std::vector<Real>&     a_cdrMobilities,
                                const std::vector<Real>&     a_cdrDiffusionCoefficients,
                                const std::vector<Real>&     a_cdrTemperatures,
                                const std::vector<Real>&     a_cdrEnergies,
                                const std::vector<RealVect>& a_cdrGradients,
                                const RealVect&              a_pos,
                                const RealVect&              a_vectorE,
                                const Real&                  a_E,
                                const Real&                  a_Etd,
                                const Real&                  a_N,
                                const Real&                  a_alpha,
                                const Real&                  a_eta,
                                const Real&                  a_time) const;

      /*!
	@brief Throw a parser error
	@param[in] a_error Error code.
      */
      void
      throwParserError(const std::string a_error) const;

      /*!
	@brief Throw a parser wearning
	@param[in] a_error Warning
      */
      void
      throwParserWarning(const std::string a_warning) const;

      /*!
	@brief Protect the @ character in a string
	@param[in] a_str Input string. If it contains the at character we throw an error. 
      */
      bool
      containsWildcard(const std::string a_str) const;

      /*!
	@brief Protect all kinds of brackets in a string.
	@param[in] a_str Input string. If it contains any bracket we throw an error.
      */
      bool
      containsBracket(const std::string a_str) const;

      /*!
	@brief Return true if string starts and ends with a paranthesis.
	@param[in] a_str Input string. 
	@return True if the first character in the string is ( and the last character is )
      */
      bool
      isBracketed(const std::string a_str) const;

      /*!
	@brief Remove whitespace from string
      */
      std::string
      trim(const std::string& a_string) const;

      /*!
	@brief Return true if species exists in map and false otherwise
	@param[in] a_name Neutral species name
      */
      bool
      isNeutralSpecies(const std::string& a_name) const;

      /*!
	@brief Return true if species exists in map and false otherwise
	@param[in] a_name Cdr species name
      */
      bool
      isPlasmaSpecies(const std::string& a_name) const;

      /*!
	@brief Return true if species exists in map and false otherwise
	@param[in] a_name Rte species name
      */
      bool
      isPhotonSpecies(const std::string& a_name) const;

      /*!
	@brief Check if file exists
	@param[in] a_filename File name
      */
      bool
      doesFileExist(const std::string a_filename) const;

      /*!
	@brief Add photoionization products to transport equations source terms.
	@param[inout] a_cdrSources   Source terms for CDR densities.
	@param[in]    a_rteDensities RTE mesh densities.
	@param[in]    a_position     Physical coordinates
	@param[in]    a_E            Electric field (SI units)
	@param[in]    a_dt           Time step
	@param[in]    a_dx           Grid resolution. 
      */
      virtual void
      addPhotoIonization(std::vector<Real>&       a_cdrSources,
                         const std::vector<Real>& a_rteDensities,
                         const RealVect           a_position,
                         const Real               a_E,
                         const Real               a_dt,
                         const Real               a_dx) const;

      /*!
	@brief Routine for integrating the reactive-only problem using various algorithms. 
	@param[inout] a_cdrDensities     On input, contains n(t). On output it contains n(t+dt).
	@param[out]   a_photonProduction On input, should be equal to zero. On output it will contain the number of photons produced during the time step. 
	@param[in]    a_cdrGradients     CDR gradients at time a_time
	@param[in]    a_E                Electric field
	@param[in]    a_pos              Physical coordinates
	@param[in]    a_dx               Grid resolution
	@param[in]    a_dt               Time step
	@param[in]    a_kappa            Volume fraction 
      */
      virtual void
      integrateReactions(std::vector<Real>&          a_cdrDensities,
                         std::vector<Real>&          a_photonProduction,
                         const std::vector<RealVect> a_cdrGradients,
                         const RealVect              a_E,
                         const RealVect              a_pos,
                         const Real                  a_dx,
                         const Real                  a_dt,
                         const Real                  a_time,
                         const Real                  a_kappa) const;

      /*!
	@brief Routine for filling the source terms in the reactive problem.
	@param[out]   a_cdrSources       Contains source term for CDR equations. 
	@param[out]   a_rteSources       Contains source terms for RTE equations. 
	@param[in]    a_cdrDensities     CDR densities
	@param[in]    a_cdrGradients     CDR gradients at time a_time
	@param[in]    a_E                Electric field
	@param[in]    a_pos              Physical coordinates
	@param[in]    a_dx               Grid resolution
	@param[in]    a_time             Time
	@param[in]    a_kappa            Volume fraction 
      */
      void
      fillSourceTerms(std::vector<Real>&          a_cdrSources,
                      std::vector<Real>&          a_rteSources,
                      const std::vector<Real>     a_cdrDensities,
                      const std::vector<RealVect> a_cdrGradients,
                      const RealVect              a_E,
                      const RealVect              a_pos,
                      const Real                  a_dx,
                      const Real                  a_time,
                      const Real                  a_kappa) const;

      /*!
	@brief Routine for integrating the reactive-only problem using the explicit Euler rule. 
	@param[inout] a_cdrDensities     On input, contains n(t). On output it contains n(t+dt).
	@param[out]   a_photonProduction On input, should be equal to zero. On output it will contain the number of photons produced during the time step. 
	@param[in]    a_cdrGradients     CDR gradients at time a_time
	@param[in]    a_E                Electric field
	@param[in]    a_pos              Physical coordinates
	@param[in]    a_dx               Grid resolution
	@param[in]    a_dt               Time step
	@param[in]    a_time             Time
	@param[in]    a_kappa            Volume fraction 
      */
      void
      integrateReactionsExplicitEuler(std::vector<Real>&          a_cdrDensities,
                                      std::vector<Real>&          a_photonProduction,
                                      const std::vector<RealVect> a_cdrGradients,
                                      const RealVect              a_E,
                                      const RealVect              a_pos,
                                      const Real                  a_dx,
                                      const Real                  a_dt,
                                      const Real                  a_time,
                                      const Real                  a_kappa) const;

      /*!
	@brief Routine for integrating the reactive-only problem using the implicit Euler rule. 
	@param[inout] a_cdrDensities     On input, contains n(t). On output it contains n(t+dt).
	@param[out]   a_photonProduction On input, should be equal to zero. On output it will contain the number of photons produced during the time step. 
	@param[in]    a_cdrGradients     CDR gradients at time a_time
	@param[in]    a_E                Electric field
	@param[in]    a_pos              Physical coordinates
	@param[in]    a_dx               Grid resolution
	@param[in]    a_dt               Time step
	@param[in]    a_time             Time
	@param[in]    a_kappa            Volume fraction 
      */
      void
      integrateReactionsImplicitEuler(std::vector<Real>&          a_cdrDensities,
                                      std::vector<Real>&          a_photonProduction,
                                      const std::vector<RealVect> a_cdrGradients,
                                      const RealVect              a_E,
                                      const RealVect              a_pos,
                                      const Real                  a_dx,
                                      const Real                  a_dt,
                                      const Real                  a_time,
                                      const Real                  a_kappa) const;

      /*!
	@brief Routine for integrating the reactive-only problem using a second order Runge-Kutta method. 
	@details This includes the tableu through the a_tableuAlpha parameter. 
	@param[inout] a_cdrDensities     On input, contains n(t). On output it contains n(t+dt).
	@param[out]   a_photonProduction On input, should be equal to zero. On output it will contain the number of photons produced during the time step. 
	@param[in]    a_cdrGradients     CDR gradients at time a_time
	@param[in]    a_E                Electric field
	@param[in]    a_pos              Physical coordinates
	@param[in]    a_dx               Grid resolution
	@param[in]    a_dt               Time step
	@param[in]    a_time             Time
	@param[in]    a_kappa            Volume fraction 
	@param[in]    a_tableuAlpha      RK2 tableu alpha. Use 0.5 for midpoint and 1.0 for trapezoidal (Heun's method)
      */
      void
      integrateReactionsExplicitRK2(std::vector<Real>&          a_cdrDensities,
                                    std::vector<Real>&          a_photonProduction,
                                    const std::vector<RealVect> a_cdrGradients,
                                    const RealVect              a_E,
                                    const RealVect              a_pos,
                                    const Real                  a_dx,
                                    const Real                  a_dt,
                                    const Real                  a_time,
                                    const Real                  a_kappa,
                                    const Real                  a_tableuAlpha) const;

      /*!
	@brief Routine for integrating the reactive-only problem using the foruth order Runge-Kutta method. 
	@param[inout] a_cdrDensities     On input, contains n(t). On output it contains n(t+dt).
	@param[out]   a_photonProduction On input, should be equal to zero. On output it will contain the number of photons produced during the time step. 
	@param[in]    a_cdrGradients     CDR gradients at time a_time
	@param[in]    a_E                Electric field
	@param[in]    a_pos              Physical coordinates
	@param[in]    a_dx               Grid resolution
	@param[in]    a_dt               Time step
	@param[in]    a_time             Time
	@param[in]    a_kappa            Volume fraction 
      */
      void
      integrateReactionsExplicitRK4(std::vector<Real>&          a_cdrDensities,
                                    std::vector<Real>&          a_photonProduction,
                                    const std::vector<RealVect> a_cdrGradients,
                                    const RealVect              a_E,
                                    const RealVect              a_pos,
                                    const Real                  a_dx,
                                    const Real                  a_dt,
                                    const Real                  a_time,
                                    const Real                  a_kappa) const;
    };
  } // namespace CdrPlasma
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
