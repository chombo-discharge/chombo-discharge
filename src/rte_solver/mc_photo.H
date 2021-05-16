/*!
  @file   mc_photo.H
  @author Robert Marskar
  @date   Jan. 2018
  @todo   Particle deposition should use the new EBMeshInterp and not MeshInterp
*/

#ifndef _MC_PHOTO_
#define _MC_PHOTO_

#include "rte_solver.H"
#include "EBParticleInterp.H"
#include "particle_container.H"
#include "wall_bc.H"
#include "photon.H"

#include <Particle.H>
#include <ParticleData.H>
#include <BinItem.H>

#include <random>

#include "CD_NamespaceHeader.H"

/*!
  @brief Radiative tranfer equation solver using stationary Monte-Carlo
*/
class mc_photo : public rte_solver {
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
    @brief Advance RTE and deposit photons on a_phi
    @details This routine switches between stationary and transient solvers depending on what you ask for. 
  */
  virtual bool advance(const Real a_dt, EBAMRCellData& a_phi, const EBAMRCellData& a_source, const bool a_zerophi = false);

  /*!
    @brief Instantaneous solver or not
  */
  virtual bool is_instantaneous();

  /*!
    @brief Generate photons
  */
  virtual void generate_photons(particle_container<photon>& a_particles, const EBAMRCellData& a_source, const Real a_dt);

  /*!
    @brief Move photons and absorb them. 
    @details This is the routine that is used to move photons when the solver is stationary. 
  */
  virtual void advance_photons_stationary(particle_container<photon>& a_bulk_photons,
					  particle_container<photon>& a_eb_photons,
					  particle_container<photon>& a_domain_photons,
					  particle_container<photon>& a_photons);

  /*!
    @brief This moves the photons over a time step dt. If a photon is absorbed in the gas, it is put in
    a_bulk_photons. If it is absorbed on the EB, it is put in a_eb_photons etc. 
    @details This is the routine that is used to move photons when the solver is transient. 
  */
  virtual void advance_photons_transient(particle_container<photon>& a_bulk_photons,
					 particle_container<photon>& a_eb_photons,
					 particle_container<photon>& a_domain_photons,
					 particle_container<photon>& a_ophotons,
					 const Real                  a_dt);

  /*!
    @brief Remap photons
  */
  virtual void remap();

  /*!
    @brief Remap photons
  */
  virtual void remap(particle_container<photon>& a_photons);

  /*!
    @brief Deposit photons
  */
  virtual void deposit_photons();

  /*!
    @brief Deposit photons
  */
  virtual void deposit_photons(EBAMRCellData&                    a_phi,
			       const particle_container<photon>& a_particles,
			       const DepositionType::Which&      a_deposition);

  /*!
    @brief Deposit photons
  */
  virtual void deposit_photons(EBAMRCellData&               a_phi,
			       const AMRParticles<photon>&  a_photons,
			       const DepositionType::Which& a_deposition);

  /*!
    @brief Sort the m_photons container by cell
  */
  virtual void sort_photons_by_cell();

  /*!
    @brief Sort the m_photons container by cell
  */
  virtual void sort_photons_by_patch();

  /*!
    @brief Sort the m_bulk_photons container by cell
  */
  virtual void sort_bulk_photons_by_cell();

  /*!
    @brief Sort the m_bulk_photons container by cell
  */
  virtual void sort_bulk_photons_by_patch();

  /*!
    @brief Sort the m_source_photons container by cell
  */
  virtual void sort_source_photons_by_cell();

  /*!
    @brief Sort the m_source_photons container by cell
  */
  virtual void sort_source_photons_by_patch();

  /*!
    @brief Sort the m_source_photons container by cell
  */
  virtual void sort_eb_photons_by_cell();

  /*!
    @brief Sort the m_source_photons container by cell
  */
  virtual void sort_eb_photons_by_patch();

  /*!
    @brief Sort the m_source_photons container by cell
  */
  virtual void sort_domain_photons_by_cell();

  /*!
    @brief Sort the m_source_photons container by cell
  */
  virtual void sort_domain_photons_by_patch();

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
  virtual void compute_boundary_flux(EBAMRIVData& a_ebFlux, const EBAMRCellData& a_phi);

  /*!
    @brief Compute the domain flux. For Eddington, the domain flux is = c*phi/2
  */
  virtual void compute_domain_flux(EBAMRIFData& a_domainflux, const EBAMRCellData& a_phi);

  /*!
    @brief Compute the flux
    @details For diffusive models, the flux will be something like grad(a_phi). For higher-order models, the flux 
    will be contained in a_phi (somehow);
  */
  virtual void computeFlux(EBAMRCellData& a_flux, const EBAMRCellData& a_phi);

  /*!
    @brief Get isotropic part. 
  */
  virtual void compute_density(EBAMRCellData& a_isotropic, const EBAMRCellData& a_phi);

  /*!
    @brief Clear data holder
  */
  virtual void clear();

  /*!
    @brief Clear data holder
  */
  virtual void clear(particle_container<photon>& a_photon);

  /*!
    @brief Clear data holder
  */
  virtual void clear(AMRParticles<photon>& a_photons);

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
  virtual Vector<std::string> get_plotVariableNames() const;

  /*!
    @brief Get PVR buffer
  */
  virtual int get_pvr_buffer() const;

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
    @brief Count number of photons in particle list
  */
  virtual int count_photons(const AMRParticles<photon>& a_photons) const;

  /*!
    @brief Count number of outcast photons in particle list
  */
  virtual int count_outcast(const AMRParticles<photon>& a_photons) const;

  /*!
    @brief Get photons
  */
  virtual particle_container<photon>& get_photons();

  /*!
    @brief Get bulk photons
  */
  virtual particle_container<photon>& get_bulk_photons();

  /*!
    @brief Get EB photons
  */
  virtual particle_container<photon>& get_eb_photons();

  /*!
    @brief Get domain photons
  */
  virtual particle_container<photon>& get_domain_photons();

  /*!
    @brief Get source photons
  */
  virtual particle_container<photon>& getSource_photons();
  
protected:

  enum photon_generation {
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
    @brief If true, the NUMBER of of photons will be deposited in each cell
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
    @brief Check if m_bulk_photons should be plotted
  */
  bool m_plot_bulk_phot;

  /*!
    @brief Check if source_bulk_photons should be plotted
  */
  bool m_plotSource_phot;

  /*!
    @brief Check if m_eb_photons should be plotted
  */
  bool m_plot_eb_phot;

  /*!
    @brief Check if m_eb_photons should be plotted
  */
  bool m_plot_dom_phot;

  /*!
    @brief Number of photons per pseudophoton
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
    @brief Number of photons
  */
  long long m_num_photons;

  /*!
    @brief Photon generation type
  */
  photon_generation m_photogen;

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
  particle_container<photon> m_photons;

  /*!
    @brief Photons absorbed in the volume
  */
  particle_container<photon> m_bulk_photons;

  /*!
    @brief This is a particle container for photons that crossed EBs. It is filled during the advance step. 
  */
  particle_container<photon> m_eb_photons;

  /*!
    @brief This is a particle container for photons that crossed boundaries. It is filled during the advance step. 
  */
  particle_container<photon> m_domain_photons;

  /*!
    @brief This is a particle container that can be used to add photons directly. 
  */
  particle_container<photon> m_source_photons;

  /*!
    @brief Domain boundary conditions
  */
  Vector<wallbc::which_bc> m_domainbc;

  /*!
    @brief Draw photons
  */
  int draw_photons(const Real a_source, const Real a_volume, const Real a_dt);

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
				 const AMRParticles<photon>& a_particles,
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
    @brief Parse pseudophotons
  */
  void parse_pseudophotons();

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
