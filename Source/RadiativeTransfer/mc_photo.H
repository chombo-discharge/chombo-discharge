/*!
  @file   mc_photo.H
  @author Robert Marskar
  @date   Jan. 2018
  @todo   Particle deposition should use the new EBMeshInterp and not MeshInterp
*/

#ifndef _MC_PHOTO_
#define _MC_PHOTO_

#include <CD_RtSolver.H>
#include <CD_EbParticleInterp.H>
#include <CD_ParticleContainer.H>
#include <CD_WallBc.H>
#include <CD_Photon.H>

#include <Particle.H>
#include <ParticleData.H>
#include <BinItem.H>

#include <random>

#include "CD_NamespaceHeader.H"

/*!
  @brief Radiative tranfer equation solver using stationary Monte-Carlo
*/
class mc_photo : public RtSolver {
public:

  /*!
    @brief Constructor
  */
  mc_photo();

  /*!
    @brief Destructor
  */
  virtual ~mc_photo();

  /*!
    @brief Advance RTE and deposit Photons on a_phi
    @details This routine switches between stationary and transient solvers depending on what you ask for. 
  */
  virtual bool advance(const Real a_dt, EBAMRCellData& a_phi, const EBAMRCellData& a_source, const bool a_zerophi = false);

  /*!
    @brief Instantaneous solver or not
  */
  virtual bool is_instantaneous();

  /*!
    @brief Generate Photons
  */
  virtual void generate_Photons(ParticleContainer<Photon>& a_particles, const EBAMRCellData& a_source, const Real a_dt);

  /*!
    @brief Move Photons and absorb them. 
    @details This is the routine that is used to move Photons when the solver is stationary. 
  */
  virtual void advance_Photons_stationary(ParticleContainer<Photon>& a_bulk_Photons,
					  ParticleContainer<Photon>& a_eb_Photons,
					  ParticleContainer<Photon>& a_domain_Photons,
					  ParticleContainer<Photon>& a_Photons);

  /*!
    @brief This moves the Photons over a time step dt. If a Photon is absorbed in the gas, it is put in
    a_bulk_Photons. If it is absorbed on the EB, it is put in a_eb_Photons etc. 
    @details This is the routine that is used to move Photons when the solver is transient. 
  */
  virtual void advance_Photons_transient(ParticleContainer<Photon>& a_bulk_Photons,
					 ParticleContainer<Photon>& a_eb_Photons,
					 ParticleContainer<Photon>& a_domain_Photons,
					 ParticleContainer<Photon>& a_oPhotons,
					 const Real                  a_dt);

  /*!
    @brief Remap Photons
  */
  virtual void remap();

  /*!
    @brief Remap Photons
  */
  virtual void remap(ParticleContainer<Photon>& a_Photons);

  /*!
    @brief Deposit Photons
  */
  virtual void deposit_Photons();

  /*!
    @brief Deposit Photons
  */
  virtual void deposit_Photons(EBAMRCellData&                    a_phi,
			       const ParticleContainer<Photon>& a_particles,
			       const DepositionType::Which&      a_deposition);

  /*!
    @brief Deposit Photons
  */
  virtual void deposit_Photons(EBAMRCellData&               a_phi,
			       const AMRParticles<Photon>&  a_Photons,
			       const DepositionType::Which& a_deposition);

  /*!
    @brief Sort the m_Photons container by cell
  */
  virtual void sort_Photons_by_cell();

  /*!
    @brief Sort the m_Photons container by cell
  */
  virtual void sort_Photons_by_patch();

  /*!
    @brief Sort the m_bulk_Photons container by cell
  */
  virtual void sort_bulk_Photons_by_cell();

  /*!
    @brief Sort the m_bulk_Photons container by cell
  */
  virtual void sort_bulk_Photons_by_patch();

  /*!
    @brief Sort the m_source_Photons container by cell
  */
  virtual void sort_source_Photons_by_cell();

  /*!
    @brief Sort the m_source_Photons container by cell
  */
  virtual void sort_source_Photons_by_patch();

  /*!
    @brief Sort the m_source_Photons container by cell
  */
  virtual void sort_eb_Photons_by_cell();

  /*!
    @brief Sort the m_source_Photons container by cell
  */
  virtual void sort_eb_Photons_by_patch();

  /*!
    @brief Sort the m_source_Photons container by cell
  */
  virtual void sort_domain_Photons_by_cell();

  /*!
    @brief Sort the m_source_Photons container by cell
  */
  virtual void sort_domain_Photons_by_patch();

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
  virtual void clear(AMRParticles<Photon>& a_Photons);

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
  virtual int getPVR_buffer() const;

  /*!
    @brief Get the halo buffer
  */
  virtual int get_halo_buffer() const;

  /*!
    @brief Set the PVR buffer
  */
  virtual void set_pvr_buffer(const int a_buffer);

  /*!
    @brief Set the halo buffer
  */
  virtual void set_halo_buffer(const int a_buffer);

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
  virtual int count_Photons(const AMRParticles<Photon>& a_Photons) const;

  /*!
    @brief Count number of outcast Photons in particle list
  */
  virtual int count_outcast(const AMRParticles<Photon>& a_Photons) const;

  /*!
    @brief Get Photons
  */
  virtual ParticleContainer<Photon>& get_Photons();

  /*!
    @brief Get bulk Photons
  */
  virtual ParticleContainer<Photon>& get_bulk_Photons();

  /*!
    @brief Get EB Photons
  */
  virtual ParticleContainer<Photon>& get_eb_Photons();

  /*!
    @brief Get domain Photons
  */
  virtual ParticleContainer<Photon>& get_domain_Photons();

  /*!
    @brief Get source Photons
  */
  virtual ParticleContainer<Photon>& getSource_Photons();
  
protected:

  enum Photon_generation {
    deterministic = 0,
    stochastic = 1,
  };

  
  enum source_type {
    number   = 0,
    per_vols = 1, 
    per_vol  = 2,
    per_s    = 3,
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
    @brief Check if m_Photons should be plotted
  */
  bool m_plot_phot;

  /*!
    @brief Check if m_bulk_Photons should be plotted
  */
  bool m_plot_bulk_phot;

  /*!
    @brief Check if source_bulk_Photons should be plotted
  */
  bool m_plotSource_phot;

  /*!
    @brief Check if m_eb_Photons should be plotted
  */
  bool m_plot_eb_phot;

  /*!
    @brief Check if m_eb_Photons should be plotted
  */
  bool m_plot_dom_phot;

  /*!
    @brief Number of Photons per pseudoPhoton
  */
  int m_max_Photons;

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
  long long m_num_Photons;

  /*!
    @brief Photon generation type
  */
  Photon_generation m_photogen;

  /*!
    @brief Source type
  */
  source_type m_src_type;

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
  ParticleContainer<Photon> m_Photons;

  /*!
    @brief Photons absorbed in the volume
  */
  ParticleContainer<Photon> m_bulk_Photons;

  /*!
    @brief This is a particle container for Photons that crossed EBs. It is filled during the advance step. 
  */
  ParticleContainer<Photon> m_eb_Photons;

  /*!
    @brief This is a particle container for Photons that crossed boundaries. It is filled during the advance step. 
  */
  ParticleContainer<Photon> m_domain_Photons;

  /*!
    @brief This is a particle container that can be used to add Photons directly. 
  */
  ParticleContainer<Photon> m_source_Photons;

  /*!
    @brief Domain boundary conditions
  */
  Vector<wallbc::which_bc> m_domainbc;

  /*!
    @brief Draw Photons
  */
  int draw_Photons(const Real a_source, const Real a_volume, const Real a_dt);

  /*!
    @brief Mapping function for domain boundary conditions
  */
  int domainbc_map(const int a_dir, const Side::LoHiSide a_side);

  /*!
    @brief Random Poisson trial
  */
  int random_fieldSolver(const Real a_mean);

  /*!
    @brief Random exponential trial
  */
  Real random_exponential(const Real a_mean);

  /*!
    @brief Random direction
  */
  RealVect random_direction();

#if CH_SPACEDIM==2
  /*!
    @brief Random direction in 2D
  */
  RealVect random_direction2D();
#elif CH_SPACEDIM==3
  /*!
    @brief Random direction in 2D
  */
  RealVect random_direction3D();
#endif

  /*!
    @brief This computes the "conservative" deposition, multiplied by kappa
  */
  void deposit_kappaConservative(EBAMRCellData&              a_phi,
				 const AMRParticles<Photon>& a_particles,
				 const DepositionType::Which a_deposition);

  /*!
    @brief Make the "non-conservative" kappa deposition
  */
  void deposit_nonConservative(EBAMRIVData& a_depositionNC, const EBAMRCellData& a_depositionKappaC);

  /*!
    @brief Make the hybrid deposition. Also compute the mass difference
  */
  void deposit_hybrid(EBAMRCellData& a_depositionH, EBAMRIVData& a_massDifference, const EBAMRIVData& a_depositionNC);

  /*!
    @brief Increment redistribution objects with mass to be redistributed
  */
  void incrementRedist(const EBAMRIVData& a_massDifference);

  /*!
    @brief Do redistribution on all levels, no coarse-fine stuff
  */
  void level_redistribution(EBAMRCellData& a_phi);

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
  void parse_rng();

  /*!
    @brief Parse the conservation
  */
  void parseDivergenceComputation();

  /*!
    @brief Parse pseudoPhotons
  */
  void parse_pseudoPhotons();

  /*!
    @brief Parse photogeneration type
  */
  void parse_photogen();

  /*!
    @brief Parse whether stationary or not
  */
  void parse_instantaneous();

  /*!
    @brief Parse source type
  */
  void parse_source_type();

  /*!
    @brief Parse deposition
  */
  void parse_deposition();

  /*!
    @brief Parse bisection step
  */
  void parse_bisect_step();

  /*!
    @brief Parse domain BCs
  */
  void parseDomainBc();

  /*!
    @brief Parse PVR buffer
  */
  void parse_pvr_buffer();

  /*!
    @brief Parse plot variables
  */
  void parsePlotVariables();

};
#include "CD_NamespaceFooter.H"

#endif
