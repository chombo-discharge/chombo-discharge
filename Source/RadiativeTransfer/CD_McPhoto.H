/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_McPhoto.H
  @brief  Declaration of a monte carlo class for doing radiative transfer
  @author Robert Marskar
*/

#ifndef CD_McPhoto_H
#define CD_McPhoto_H

// Std includes
#include <random>

// Chombo includes
#include <Particle.H>
#include <ParticleData.H>
#include <BinItem.H>

// Our includes
#include <CD_RtSolver.H>
#include <CD_EbParticleInterp.H>
#include <CD_ParticleContainer.H>
#include <CD_WallBc.H>
#include <CD_Photon.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Radiative tranfer equation solver using Monte-Carlo simulation.
  @details This class is, by default, a non-stationary radiative transfer solver which simulates computational photons that can also be deposited on the mesh. Various
  options are available for specifying how the photons are generated, and how they are propagated (e.g. with speed of light or instantaneously). 
*/
class McPhoto : public RtSolver {
public:

  /*!
    @brief Enum class for identifying various containers. Only used for interface reasons. 
  */
  enum class WhichContainer {
    Photons,
    Bulk,
    EB,
    Domain,
    Source
  };

  /*!
    @brief Constructor
  */
  McPhoto();

  /*!
    @brief Destructor
  */
  virtual ~McPhoto();

  /*!
    @brief Advance RTE and deposit photon particles on the mesh
    @param[in]  a_dt Time step
    @param[out] a_phi     Photon density on the mesh
    @param[in]  a_source  Source term, i.e. number of photons produced per unit time. The photons will be generated
    based on what the user specifies in the input parameters for this class. I.e., whether or not a_source will generate stochastic or deterministic photons, and
    scaled based on what a_source actually contains (number of photons, number/volume, number/second, number/(volume*second)). 
    @param[in]  a_zeroPhi Dead parameter not used in this implementation. 
    @note This is the "fluid interface" to the radiative transfer solver in which the solver uses discrete photons for sampling the isotropic density. Particle codes
    will probably use advancePhotonsInstantaneous or advancePhotonsTransient .
    @details This routine switches between stationary and transient solvers depending on whether or not the solver is stationary. 
  */
  virtual bool advance(const Real a_dt, EBAMRCellData& a_phi, const EBAMRCellData& a_source, const bool a_zeroPhi = false) override;

  /*!
    @brief Instantaneous solver or not. 
  */
  virtual bool isInstantaneous();

  /*!
    @brief Generate computational photon sfrom a source term 
    @param[out] a_photons Computational photons. 
    @param[in]  a_source  Source term, i.e. number of photons produced per unit time. The photons will be generated
    based on what the user specifies in the input parameters for this class. I.e., whether or not a_source will generate stochastic or deterministic photons, and
    scaled based on what a_source actually contains (number of photons, number/volume, number/second, number/(volume*second)). 
    @param[in] a_dt Time step
    @details This will generate photons and put them in a_particles
  */
  virtual void generatePhotons(ParticleContainer<Photon>& a_photons, const EBAMRCellData& a_source, const Real a_dt);

  /*!
    @brief Move photons and absorb them on various objects
    @param[out]   a_bulkPhotons   Photons absorbed on the mesh
    @param[out]   a_ebPhotons     Photons absorbed on the EB
    @param[out]   a_domainPhotons Photons absorbed on the domain edges (faces)
    @param[inout] a_photons       Original photons
    @note This routine moves photons with an instantaneous kernel, i.e. all photons are always absorbed and none are left behind as free-flight photons (i.e. a_photons will
    be empty on output)
  */
  virtual void advancePhotonsInstantaneous(ParticleContainer<Photon>& a_bulkPhotons,
					   ParticleContainer<Photon>& a_ebPhotons,
					   ParticleContainer<Photon>& a_domainPhotons,
					   ParticleContainer<Photon>& a_photons);

  /*!
    @brief Move photons and absorb them on various objects
    @param[out]   a_bulkPhotons   Photons absorbed on the mesh
    @param[out]   a_ebPhotons     Photons absorbed on the EB
    @param[out]   a_domainPhotons Photons absorbed on the domain edges (faces)
    @param[inout] a_photons       Original photons
    @param[in]    a_dt            Time step
    @note This routine moves photons with a transient kernel, i.e. all photons propagate at most c*dt and checks are made to determine if they are absorbed on the mesh
    or various objects (EB, domain). a_photons may be non-empty on output, in which case there are still free-flight photons in the solver. 
  */
  virtual void advancePhotonsTransient(ParticleContainer<Photon>& a_bulkPhotons,
				       ParticleContainer<Photon>& a_ebPhotons,
				       ParticleContainer<Photon>& a_domainPhotons,
				       ParticleContainer<Photon>& a_photons,
				       const Real                 a_dt);

  /*!
    @brief Remap computational particles. This remaps m_photons
  */
  virtual void remap();

  /*!
    @brief Remap computational particles
    @param[in] a_photons Computational particles to be remapped.
  */
  virtual void remap(ParticleContainer<Photon>& a_photons);

  /*!
    @brief Deposit photons on the mesh. 
    @details This deposits m_photons onto m_phi
  */
  virtual void depositPhotons();

  /*!
    @brief Deposit photons
    @param[out] a_phi        Deposited mesh density
    @param[in]  a_photons    Computational photons to be deposited
    @param[in]  a_deposition Deposition method
    @note Calls the other version
  */
  virtual void depositPhotons(EBAMRCellData&                   a_phi,
			      const ParticleContainer<Photon>& a_particles,
			      const DepositionType&            a_deposition);

  /*!
    @brief Deposit on the mesh. 
    @param[out] a_phi        Deposited mesh density
    @param[in]  a_photons    Computational photons to be deposited
    @param[in]  a_deposition Deposition method
  */
  virtual void depositPhotons(EBAMRCellData&               a_phi,
			      const AMRParticles<Photon>&  a_photons,
			      const DepositionType&        a_deposition);

  /*!
    @brief Sort container by cell
    @details WhichContainer::Photon = m_photons, WhichContainer::EB = m_ebPhotons and so on. 
    @param[in] a_which Which container to sort. 
  */
  virtual void sortPhotonsByCell(const WhichContainer& a_which);

  /*!
    @brief Sort container by patch
    @details WhichContainer::Photon = m_photons, WhichContainer::EB = m_ebPhotons and so on. 
    @param[in] a_which Which container to sort. 
  */
  virtual void sortPhotonsByPatch(const WhichContainer& a_which);

  /*!
    @brief Sort the m_photons container by cell
  */
  virtual void sortPhotonsByCell();

  /*!
    @brief Sort the m_photons container by patch
  */
  virtual void sortPhotonsByPatch();

  /*!
    @brief Sort the m_bulkPhotons container by cell
  */
  virtual void sortBulkPhotonsByCell();

  /*!
    @brief Sort the m_bulkPhotons container by patch
  */
  virtual void sortBulkPhotonsByPatch();

  /*!
    @brief Sort the m_sourcePhotons container by cell
  */
  virtual void sortSourcePhotonsByCell();

  /*!
    @brief Sort the m_sourcePhotons container by patch
  */
  virtual void sortSourcePhotonsByPatch();

  /*!
    @brief Sort the m_sourcePhotons container by cell
  */
  virtual void sortEbPhotonsByCell();

  /*!
    @brief Sort the m_sourcePhotons container by patch
  */
  virtual void sortEbPhotonsByPatch();

  /*!
    @brief Sort the m_sourcePhotons container by cell
  */
  virtual void sortDomainPhotonsByCell();

  /*!
    @brief Sort the m_sourcePhotons container by patch
  */
  virtual void sortDomainPhotonsByPatch();

  /*!
    @brief Parse class options
  */
  virtual void parseOptions() override;

  /*!
    @brief Parse runtime options
  */
  virtual void parseRuntimeOptions() override;
  
  /*!
    @brief Allocate internal storage
  */
  virtual void allocateInternals() override;

  /*!
    @brief preRegrid operations
    @param[in] a_base Coarsest level which changes during regrid
    @param[in] a_oldFinestLevel Finest grid level before regrid
  */
  virtual void preRegrid(const int a_base, const int a_oldFinestLevel) override;

  /*!
    @brief Deallocate internal storage
  */
  virtual void deallocateInternals() override;

  /*!
    @brief Regrid function for this class
    @param[in] a_lmin Coarsest level that changed during regrid
    @param[i]  a_oldFinestLevel Finest grid level before the regrid
    @param[i]  a_newFinestLevel Finest grid level after the regrid
  */
  virtual void regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) override;

  /*!
    @brief Register operators that this solver needs.
  */
  virtual void registerOperators() override;

  /*!
    @brief Compute the boundary flux. 
    @note This sets the boundary flux to zero. If you want to use a boundary flux, use it through m_ebPhotons
  */
  virtual void computeBoundaryFlux(EBAMRIVData& a_ebFlux, const EBAMRCellData& a_phi) override;

  /*!
    @brief Compute the domain flux. 
    @note This sets the boundary flux to zero. If you want to use a domain flux, use it through m_domainPhotons
  */
  virtual void computeDomainFlux(EBAMRIFData& a_domainflux, const EBAMRCellData& a_phi) override;

  /*!
    @brief Compute the flux
    @note Issues an error -- currently don't know how to compute the mesh flux from computational photons....
  */
  virtual void computeFlux(EBAMRCellData& a_flux, const EBAMRCellData& a_phi) override;

  /*!
    @brief Compute isotropic radiative density from mesh solution
    @note Issues an error - this routine doesn't make a lot of sense, really. 
  */
  virtual void computeDensity(EBAMRCellData& a_isotropic, const EBAMRCellData& a_phi) override;

  /*!
    @brief Clear data holder
  */
  virtual void clear();

  /*!
    @brief Clear data holder
  */
  virtual void clear(ParticleContainer<Photon>& a_Photon);

  /*!
    @brief Clear data holder
  */
  virtual void clear(AMRParticles<Photon>& a_photons);

  /*!
    @brief Write plot file
  */
  virtual void writePlotFile();

  /*!
    @brief Write plot data
  */
  virtual void writePlotData(EBAMRCellData& a_output, int& a_comp);

  /*!
    @brief Write checkpoint data into handle
  */
  virtual void writeCheckpointLevel(HDF5Handle& a_handle, const int a_level) const;

  /*!
    @brief Read checkpoint data from handle
  */
  virtual void readCheckpointLevel(HDF5Handle& a_handle, const int a_level);

  /*!
    @brief Get output plot names
  */
  virtual Vector<std::string> getPlotVariableNames() const;

  /*!
    @brief Get PVR buffer
  */
  virtual int getPVRBuffer() const;

  /*!
    @brief Get the halo buffer
  */
  virtual int getHaloBuffer() const;

  /*!
    @brief Set the PVR buffer
  */
  virtual void setPVRBuffer(const int a_buffer);

  /*!
    @brief Set the halo buffer
  */
  virtual void setHalobuffer(const int a_buffer);

  /*!
    @brief Get number of output fields
  */
  virtual int getNumberOfPlotVariables() const;

  /*!
    @brief Return required number of ghost cells
  */
  virtual int queryGhost() const;

  /*!
    @brief Count number of Photons in particle list
  */
  virtual int countPhotons(const AMRParticles<Photon>& a_photons) const;

  /*!
    @brief Count number of outcast Photons in particle list
  */
  virtual int countOutcast(const AMRParticles<Photon>& a_photons) const;

  /*!
    @brief Get Photons
  */
  virtual ParticleContainer<Photon>& getPhotons();

  /*!
    @brief Get bulk Photons
  */
  virtual ParticleContainer<Photon>& getBulkPhotons();

  /*!
    @brief Get EB Photons
  */
  virtual ParticleContainer<Photon>& getEbPhotons();

  /*!
    @brief Get domain Photons
  */
  virtual ParticleContainer<Photon>& getDomainPhotons();

  /*!
    @brief Get source Photons
  */
  virtual ParticleContainer<Photon>& getSourcePhotons();
  
protected:

  /*!
    @brief Enum for interpreting how photons are generated when using fluid codes
  */
  enum class PhotonGeneration {
    Deterministic,
    Stochastic
  };

  /*!
    @brief Enum for adding flexibility in what the fluid source term contains. 
  */
  enum class SourceType {
    Number,
    PerVol,
    PerVolSecond,
    PerSecond
  };

  /*!
    @brief Instantaneous transport or not
  */
  bool m_instantaneous;

  /*!
    @brief Flag for blending the deposition clouds with the nonconservative divergence
  */
  bool m_blendConservation;

  /*!
    @brief If true, the NUMBER of of Photons will be deposited in each cell
  */
  bool m_deposit_numbers;

  /*!
    @brief Switch for plotting numbers or densities
  */
  bool m_plotNumbers;

  /*!
    @brief Check if m_photons should be plotted
  */
  bool m_plot_phot;

  /*!
    @brief Check if m_bulkPhotons should be plotted
  */
  bool m_plot_bulk_phot;

  /*!
    @brief Check if source_bulkPhotons should be plotted
  */
  bool m_plot_source_phot;

  /*!
    @brief Check if m_ebPhotons should be plotted
  */
  bool m_plot_eb_phot;

  /*!
    @brief Check if m_ebPhotons should be plotted
  */
  bool m_plot_dom_phot;

  /*!
    @brief Number of Photons per pseudoPhoton
  */
  int m_max_photons;

  /*!
    @brief RNG seed
  */
  int m_seed;

  /*!
    @brief Threshold for swapping Poisson and exponential distributions
  */
  int m_poiss_exp_swap;

  /*!
    @brief PVR buffer
  */
  int m_pvr_buffer;

  /*!
    @brief Halo size for particles
  */
  int m_halo_buffer;

  /*!
    @brief Number of Photons
  */
  long long m_num_photons;

  /*!
    @brief Photon generation type
  */
  PhotonGeneration m_photogen;

  /*!
    @brief Source type
  */
  SourceType m_src_type;

  /*!
    @brief Bisection step length for boundary intersections
  */
  Real m_bisect_step;

  /*!
    @brief RNG engine
  */
  std::mt19937_64* m_rng; // Engine

  /*!
    @brief Distribution on [0,1]
  */
  std::uniform_real_distribution<Real>* m_udist01; 

  /*!
    @brief Distribution on [-1,1]
  */
  std::uniform_real_distribution<Real>* m_udist11; 

  /*!
    @brief Deposition type
  */
  DepositionType m_deposition;

  /*!
    @brief Plot deposition type
  */
  DepositionType m_plot_deposition;

  /*!
    @brief Coarse data for interpolation of deposition clouds
  */
  EBAMRCellData m_scratch;

  /*!
    @brief Scratch storage for holding the non-conservative deposition
  */
  EBAMRIVData m_depositionNC;

  /*!
    @brief Scratch storage for holding the mass difference when using hybrid deposition
  */
  EBAMRIVData m_massDiff;

  /*!
    @brief All particles
  */
  ParticleContainer<Photon> m_photons;

  /*!
    @brief Photons absorbed in the volume
  */
  ParticleContainer<Photon> m_bulkPhotons;

  /*!
    @brief This is a particle container for Photons that crossed EBs. It is filled during the advance step. 
  */
  ParticleContainer<Photon> m_ebPhotons;

  /*!
    @brief This is a particle container for Photons that crossed boundaries. It is filled during the advance step. 
  */
  ParticleContainer<Photon> m_domainPhotons;

  /*!
    @brief This is a particle container that can be used to add Photons directly. 
  */
  ParticleContainer<Photon> m_sourcePhotons;

  /*!
    @brief Domain boundary conditions
  */
  Vector<wallbc::which_bc> m_domainbc;

  /*!
    @brief Draw Photons
  */
  int drawPhotons(const Real a_source, const Real a_volume, const Real a_dt);

  /*!
    @brief Mapping function for domain boundary conditions
  */
  int domainBcMap(const int a_dir, const Side::LoHiSide a_side);

  /*!
    @brief Random Poisson trial
  */
  int randomPoisson(const Real a_mean);

  /*!
    @brief Random exponential trial
  */
  Real randomExponential(const Real a_mean);

  /*!
    @brief Compute and return a random direction in 2D/3D space. 
  */
  RealVect randomDirection();

#if CH_SPACEDIM==2
  /*!
    @brief Random direction in 2D
    @return Returns a random direction (in 2D)
  */
  RealVect randomDirection2D();
#elif CH_SPACEDIM==3
  /*!
    @brief Random direction in 3D
    @return Returns a random direction (in 3D)
  */
  RealVect randomDirection3D();
#endif

  /*!
    @brief This computes the "conservative" deposition, multiplied by kappa
    @param[out] a_phi       Mesh density. 
    @param[in] a_particles  Particles to be deposited
    @param[in] a_deposition Deposition method
  */
  void depositKappaConservative(EBAMRCellData&              a_phi,
				const AMRParticles<Photon>& a_particles,
				const DepositionType        a_deposition);

  /*!
    @brief Make the "non-conservative" kappa deposition
  */
  void depositNonConservative(EBAMRIVData& a_depositionNC, const EBAMRCellData& a_depositionKappaC);

  /*!
    @brief Make the hybrid deposition. Also compute the mass difference
  */
  void depositHybrid(EBAMRCellData& a_depositionH, EBAMRIVData& a_massDifference, const EBAMRIVData& a_depositionNC);

  /*!
    @brief Increment redistribution objects with mass to be redistributed
  */
  void incrementRedist(const EBAMRIVData& a_massDifference);

  /*!
    @brief Do redistribution on all levels, no coarse-fine stuff
  */
  void levelRedist(EBAMRCellData& a_phi);

  /*!
    @brief Increment the coarse-fine redist with mass difference (only for EB x CF)
  */
  void coarseFineIncrement(const EBAMRIVData& m_massDifference);

  /*!
    @brief Do coarse-fine redistribution
  */
  void coarseFineRedistribution(EBAMRCellData& a_phi);

  /*!
    @brief Parse RNG options
  */
  void parseRng();

  /*!
    @brief Parse the conservation
  */
  void parseDivergenceComputation();

  /*!
    @brief Parse pseudoPhotons
  */
  void parsePseudoPhotons();

  /*!
    @brief Parse photogeneration type
  */
  void parsePhotoGeneration();

  /*!
    @brief Parse whether stationary or not
  */
  void parseInstantaneous();

  /*!
    @brief Parse source type
  */
  void parseSourceType();

  /*!
    @brief Parse deposition
  */
  void parseDeposition();

  /*!
    @brief Parse domain BCs
  */
  void parseDomainBc();

  /*!
    @brief Parse PVR buffer
  */
  void parsePvrBuffer();

  /*!
    @brief Parse plot variables
  */
  void parsePlotVariables();

};

#include <CD_NamespaceFooter.H>

#endif
