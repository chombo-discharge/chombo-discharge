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
#include <string>

// Third-party includes
#include <nlohmann/json.hpp>

// Our includes
#include <CD_CdrPlasmaPhysics.H>
#include <CD_CdrSpeciesJSON.H>
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
	@brief Default constructor.
      */
      CdrPlasmaJSON();

      /*!
	@brief Destructor
      */
      virtual ~CdrPlasmaJSON();

      /*!
	@brief Get number of plot variables for this physics class. 
	@details This is used by CdrPlasmaStepper for pre-allocating data that will be put in a plot file. The CdrPlasmaJSON version will plot the untracked
	species if the user adds "plot": true to the species
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
	@brief Verbose or not
      */
      bool m_verbose;

      /*!
	@brief JSON definition. This is populated when calling parseJSON.
      */
      json m_json;

      /*!
	@brief Input JSON file name
      */
      std::string m_jsonFile;      

      /*!
	@brief Initial surface charge
      */
      std::function<Real(const Real a_time, const RealVect a_pos) > m_initialSigma;

      /*!
	@brief int-string encoding of the tracked species. 
	@details This is needed because we sometimes need to use the species name for indexing in the vector. 
      */
      std::map<std::string, int> m_trackedCdrSpeciesMap;

      /*!
	@brief int-string encoding of the untracked species. 
	@details This is needed because we sometimes need to use the species name for indexing in the vector. 
      */
      std::map<std::string, int> m_untrackedCdrSpeciesMap;

      /*!
	@brief Map for checking if a species is tracked or not. It returns true if a species is tracked by a CDR solver. 
	@details For given a species name we can check whether or not it is tracked through a CDR solver or just included parametrically. In reactions, this gives us
	the opportunity to check if we should get the density through a CDR solver OR through an initial data function. 
      */
      std::map<std::string, bool> m_trackedMap;

      /*!
	@brief Map over which species that will be plotted AND are not tracked. I.e. all the densities that are only parametrically defined. 
      */
      std::map<std::string, bool> m_plotCdr;

      /*!
	@brief Tracked CDR species. Note that this is the one that is needed by CdrPlasmaPhysics, i.e. they become CdrSolvers. 
      */
      std::vector<RefCountedPtr<CdrSpecies> > m_trackedCdrSpecies;

      /*!
	@brief Untracked species. These do not become CdrSolvers. 
      */
      std::vector<RefCountedPtr<CdrSpecies> > m_untrackedCdrSpecies;

      /*!
	@brief Get a species. This will return the species whether it is tracked or not. This yields an error if the species does not exist. 
	@param[in] a_name Species name. 
      */
      const RefCountedPtr<CdrSpecies>& getCdrSpecies(const std::string a_name) const;

      /*!
	@brief Get an untracked species. This will yield an error if you try to obtain a species which does not exist, or it is tracked. 
	@param[in] a_name Species name. 
      */
      const RefCountedPtr<CdrSpecies>& getUntrackedCdrSpecies(const std::string a_name) const;

      /*!
	@brief Get tracked species. This will yield an error if you try to obtain a species which does not exist, or if it is untracked. 
	@param[in] a_name Species name. 
      */
      const RefCountedPtr<CdrSpecies>& getTrackedCdrSpecies(const std::string a_name) const;

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
	@brief Get untracked species density
	@param[in] a_name Untracked species name
	@param[in] a_pos  Physical coordinates
	@param[in] a_time Time
      */
      Real getUntrackedCdrDensity(const std::string& a_name, const RealVect& a_pos, const Real& a_time) const;

      /*!
	@brief Get a tracked species density
	@param[in] a_name         Tracked species name
	@param[in] a_cdrDensities List of CDR densities for tracked species. 
      */
      Real getTrackedCdrDensity(const std::string& a_name, const Vector<Real>& a_cdrDensities) const;

      /*!
	@brief Parse class options
      */
      virtual void parseOptions();

      /*!
	@brief Parse the JSON file
      */
      virtual void parseJSON();
      
      /*!
	@brief Instantiate species. 
	@details This will instantiate the species based on the "species" field in the JSON file. If you want e.g. new boundary conditions, 
	you will have either extend this routine or overwrite it. 
      */
      virtual void instantiateCdrSpecies();

      /*!
	@brief Initialize surface charge
	@details This will instantiate the surface charge based on the "sigma" field in the JSON file. If you want more complex initial surface charges, 
	you will have to either extend this routine or overwrite it. 
      */
      virtual void initializeSigma();
    };
  }
}

#include <CD_NamespaceFooter.H>

#endif
