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
  @brief Radiative tranfer equation solver using stationary Monte-Carlo
*/
class McPhoto : public RtSolver {
public:

  /*!
    @brief Constructor
  */
  McPhoto();

  /*!
    @brief Destructor
  */
  virtual ~McPhoto();

  /*!
    @brief Advance RTE and deposit Photons on a_phi
    @details This routine switches between stationary and transient solvers depending on what you ask for. 
  */
  virtual bool advance(const Real a_dt, EBAMRCellData& a_phi, const EBAMRCellData& a_source, const bool a_zerophi = false);

  /*!
    @brief Instantaneous solver or not
  */
  virtual bool isInstantaneous();

  /*!
    @brief Generate Photons
  */
  virtual void generatePhotons(ParticleContainer<Photon>& a_particles, const EBAMRCellData& a_source, const Real a_dt);

  /*!
    @brief Move Photons and absorb them. 
    @details This is the routine that is used to move Photons when the solver is stationary. 
  */
  virtual void advancePhotonsStationary(ParticleContainer<Photon>& a_bulkPhotons,
					ParticleContainer<Photon>& a_ebPhotons,
					ParticleContainer<Photon>& a_domainPhotons,
					ParticleContainer<Photon>& a_photons);

  /*!
    @brief This moves the Photons over a time step dt. If a Photon is absorbed in the gas, it is put in
    a_bulkPhotons. If it is absorbed on the EB, it is put in a_ebPhotons etc. 
    @details This is the routine that is used to move Photons when the solver is transient. 
  */
  virtual void advancePhotonsTransient(ParticleContainer<Photon>& a_bulkPhotons,
				       ParticleContainer<Photon>& a_ebPhotons,
				       ParticleContainer<Photon>& a_domainPhotons,
				       ParticleContainer<Photon>& a_oPhotons,
				       const Real                 a_dt);

  /*!
    @brief Remap Photons
  */
  virtual void remap();

  /*!
    @brief Remap Photons
  */
  virtual void remap(ParticleContainer<Photon>& a_photons);

  /*!
    @brief Deposit Photons
  */
  virtual void depositPhotons();

  /*!
    @brief Deposit Photons
  */
  virtual void depositPhotons(EBAMRCellData&                    a_phi,
			      const ParticleContainer<Photon>&  a_particles,
			      const DepositionType::Which&      a_deposition);

  /*!
    @brief Deposit Photons
  */
  virtual void depositPhotons(EBAMRCellData&               a_phi,
			      const AMRParticles<Photon>&  a_photons,
			      const DepositionType::Which& a_deposition);

  /*!
    @brief Sort the m_photons container by cell
  */
  virtual void sortPhotonsByCell();

  /*!
    @brief Sort the m_photons container by cell
  */
  virtual void sortPhotonsByPatch();

  /*!
    @brief Sort the m_bulkPhotons container by cell
  */
  virtual void sortBulkPhotonsByCell();

  /*!
    @brief Sort the m_bulkPhotons container by cell
  */
  virtual void sortBulkPhotonsByPatch();

  /*!
    @brief Sort the m_sourcePhotons container by cell
  */
  virtual void sortSourcePhotonsByCell();

  /*!
    @brief Sort the m_sourcePhotons container by cell
  */
  virtual void sortSourcePhotonsByPatch();

  /*!
    @brief Sort the m_sourcePhotons container by cell
  */
  virtual void sortEbPhotonsByCell();

  /*!
    @brief Sort the m_sourcePhotons container by cell
  */
  virtual void sortEbPhotonsByPatch();

  /*!
    @brief Sort the m_sourcePhotons container by cell
  */
  virtual void sortDomainPhotonsByCell();

  /*!
    @brief Sort the m_sourcePhotons container by cell
  */
  virtual void sortDomainPhotonsByPatch();

  /*!
    Parse options
  */
  virtual void parseOptions();

  /*!
    Parse options
  */
  virtual void parseRuntimeOptions();
  
  /*!
    @brief Allocate internal storage
  */
  virtual void allocateInternals();

  /*!
    @brief Stuff to do before a regrid
  */
  virtual void preRegrid(const int a_base, const int a_oldFinestLevel);

  /*!
    @brief Deallocate internal storage
  */
  virtual void deallocateInternals();

  /*!
    @brief Regrid function for this class
  */
  virtual void regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel);

  /*!
    @brief Register operators
  */
  virtual void registerOperators();

  /*!
    @brief Compute the boundary flux. For Eddington, the boundary flux is = c*phi/2
  */
  virtual void computeBoundaryFlux(EBAMRIVData& a_ebFlux, const EBAMRCellData& a_phi);

  /*!
    @brief Compute the domain flux. For Eddington, the domain flux is = c*phi/2
  */
  virtual void computeDomainFlux(EBAMRIFData& a_domainflux, const EBAMRCellData& a_phi);

  /*!
    @brief Compute the flux
    @details For diffusive models, the flux will be something like grad(a_phi). For higher-order models, the flux 
    will be contained in a_phi (somehow);
  */
  virtual void computeFlux(EBAMRCellData& a_flux, const EBAMRCellData& a_phi);

  /*!
    @brief Get isotropic part. 
  */
  virtual void computeDensity(EBAMRCellData& a_isotropic, const EBAMRCellData& a_phi);

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
  DepositionType::Which m_deposition;

  /*!
    @brief Plot deposition type
  */
  DepositionType::Which m_plot_deposition;

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
    @brief Random direction
  */
  RealVect randomDirection();

#if CH_SPACEDIM==2
  /*!
    @brief Random direction in 2D
  */
  RealVect randomDirection2D();
#elif CH_SPACEDIM==3
  /*!
    @brief Random direction in 2D
  */
  RealVect randomDirection3D();
#endif

  /*!
    @brief This computes the "conservative" deposition, multiplied by kappa
  */
  void depositKappaConservative(EBAMRCellData&              a_phi,
				const AMRParticles<Photon>& a_particles,
				const DepositionType::Which a_deposition);

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
    @brief Parse bisection step
  */
  void parseBisectStep();

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
