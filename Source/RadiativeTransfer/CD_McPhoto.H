/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_McPhoto.H
  @brief  Declaration of a radiative transfer solver which uses Monte Carlo sampling of computational or real photons.
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
#include <CD_EBParticleMesh.H>
#include <CD_ParticleContainer.H>
#include <CD_Photon.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Radiative tranfer equation solver using Monte-Carlo simulation.
  @details This class is, by default, a non-stationary radiative transfer solver which simulates computational photons that can also be deposited on the mesh. Various
  options are available for specifying how the photons are generated, and how they are propagated (e.g. with speed of light or instantaneously). Note that this solver
  runs its own intersection tests with the EB due to dependencies on how the photons are moved (we store free-flight photons, absorbed photons, intersected photons, etc).
  @note The solver is missing appropriate boundary conditions on domain faces. 
*/
class McPhoto : public RtSolver
{
public:
  /*!
    @brief Enum class for identifying various containers. Only used for interface reasons. 
  */
  enum class WhichContainer
  {
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
    @brief Disallowed copy constructor
  */
  McPhoto(const McPhoto& a_other);

  /*!
    @brief Disallowed move constructor
  */
  McPhoto(const McPhoto&& a_other);

  /*!
    @brief Disallowed assignment operator
  */
  McPhoto&
  operator=(const McPhoto& a_other);

  /*!
    @brief Disallowed move assignement operator
  */
  McPhoto&
  operator=(const McPhoto&& a_other);

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
  virtual bool
  advance(const Real a_dt, EBAMRCellData& a_phi, const EBAMRCellData& a_source, const bool a_zeroPhi = false) override;

  /*!
    @brief Instantaneous solver or not. 
  */
  virtual bool
  isInstantaneous();

  /*!
    @brief Generate computational photon sfrom a source term 
    @param[out] a_photons Computational photons. 
    @param[in]  a_source  Source term, i.e. number of photons produced per unit time. The photons will be generated
    based on what the user specifies in the input parameters for this class. I.e., whether or not a_source will generate stochastic or deterministic photons, and
    scaled based on what a_source actually contains (number of photons, number/volume, number/second, number/(volume*second)). 
    @param[in] a_dt Time step
    @details This will generate photons and put them in a_particles
    @todo Implementation could use a lambda to avoid the double definition. 
  */
  virtual void
  generatePhotons(ParticleContainer<Photon>& a_photons, const EBAMRCellData& a_source, const Real a_dt);

  /*!
    @brief Move photons and absorb them on various objects
    @param[out]   a_bulkPhotons   Photons absorbed on the mesh
    @param[out]   a_ebPhotons     Photons absorbed on the EB
    @param[out]   a_domainPhotons Photons absorbed on the domain edges (faces)
    @param[inout] a_photons       Original photons
    @note This routine moves photons with an instantaneous kernel, i.e. all photons are always absorbed and none are left behind as free-flight photons (i.e. a_photons will
    be empty on output)
  */
  virtual void
  advancePhotonsInstantaneous(ParticleContainer<Photon>& a_bulkPhotons,
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
  virtual void
  advancePhotonsTransient(ParticleContainer<Photon>& a_bulkPhotons,
                          ParticleContainer<Photon>& a_ebPhotons,
                          ParticleContainer<Photon>& a_domainPhotons,
                          ParticleContainer<Photon>& a_photons,
                          const Real                 a_dt);

  /*!
    @brief Remap computational particles. This remaps m_photons
  */
  virtual void
  remap();

  /*!
    @brief Remap computational particles
    @param[in] a_photons Computational particles to be remapped.
  */
  virtual void
  remap(ParticleContainer<Photon>& a_photons);

  /*!
    @brief Deposit photons on the mesh. 
    @details This deposits m_photons onto m_phi
  */
  virtual void
  depositPhotons();

  /*!
    @brief Deposit photons
    @param[out]   a_phi        Deposited mesh density
    @param[inout] a_photons    Computational photons to be deposited
    @param[in]    a_deposition Deposition method
  */
  virtual void
  depositPhotons(EBAMRCellData& a_phi, ParticleContainer<Photon>& a_particles, const DepositionType& a_deposition);

  /*!
    @brief Sort container by cell
    @details WhichContainer::Photon = m_photons, WhichContainer::EB = m_ebPhotons and so on. 
    @param[in] a_which Which container to sort. 
  */
  virtual void
  sortPhotonsByCell(const WhichContainer& a_which);

  /*!
    @brief Sort container by patch
    @details WhichContainer::Photon = m_photons, WhichContainer::EB = m_ebPhotons and so on. 
    @param[in] a_which Which container to sort. 
  */
  virtual void
  sortPhotonsByPatch(const WhichContainer& a_which);

  /*!
    @brief Parse class options
  */
  virtual void
  parseOptions() override;

  /*!
    @brief Parse runtime options
  */
  virtual void
  parseRuntimeOptions() override;

  /*!
    @brief Allocate internal storage
  */
  virtual void
  allocate() override;

  /*!
    @brief preRegrid operations
    @param[in] a_base Coarsest level which changes during regrid
    @param[in] a_oldFinestLevel Finest grid level before regrid
  */
  virtual void
  preRegrid(const int a_base, const int a_oldFinestLevel) override;

  /*!
    @brief Deallocate internal storage
  */
  virtual void
  deallocate() override;

  /*!
    @brief Regrid function for this class
    @param[in] a_lmin Coarsest level that changed during regrid
    @param[i]  a_oldFinestLevel Finest grid level before the regrid
    @param[i]  a_newFinestLevel Finest grid level after the regrid
  */
  virtual void
  regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) override;

  /*!
    @brief Register operators that this solver needs.
  */
  virtual void
  registerOperators() override;

  /*!
    @brief Compute the boundary flux. 
    @param[out] a_ebFlux Flux on the EB
    @param[in]  a_phi    RTE solution
    @note This sets the boundary flux to zero. If you want to use a boundary flux, use it through m_ebPhotons
  */
  virtual void
  computeBoundaryFlux(EBAMRIVData& a_ebFlux, const EBAMRCellData& a_phi) override;

  /*!
    @brief Compute the domain flux. 
    @param[out] a_domainFlux Domain flux
    @param[in]  a_phi        RTE solution
    @note This sets the boundary flux to zero. If you want to use a domain flux, use it through m_domainPhotons
  */
  virtual void
  computeDomainFlux(EBAMRIFData& a_domainFlux, const EBAMRCellData& a_phi) override;

  /*!
    @brief Compute the flux
    @param[out] a_flux Flux
    @param[in]  a_phi  RTE solution
    @note Issues an error -- currently don't know how to compute the mesh flux from computational photons....
  */
  virtual void
  computeFlux(EBAMRCellData& a_flux, const EBAMRCellData& a_phi) override;

  /*!
    @brief Compute isotropic radiative density from mesh solution
    @param[out] a_isotropic Isotropic part of RTE solution
    @param[in]  a_phi       RTE solution
    @note Issues an error - this routine doesn't make a lot of sense, really. 
  */
  virtual void
  computeDensity(EBAMRCellData& a_isotropic, const EBAMRCellData& a_phi) override;

  /*!
    @brief Clear data holder
    @param[in] a_which Which container to clear
  */
  virtual void
  clear(const WhichContainer& a_which);

  /*!
    @brief Clear data holder
  */
  virtual void
  clear();

  /*!
    @brief Clear data holder
    @param[in] a_photons Which container to clear
  */
  virtual void
  clear(ParticleContainer<Photon>& a_photons);

  /*!
    @brief Clear data holder
    @param[in] a_photons Which container to clear
  */
  virtual void
  clear(AMRParticles<Photon>& a_photons);

  /*!
    @brief Write plot file
    @note Currently calls an error because it is not implemented. 
  */
  virtual void
  writePlotFile() override;

  /*!
    @brief Write plot data
    @param[inout] a_ouput Output data holder
    @param[inout] a_icomp Starting component where we begin to write data to a_output
  */
  virtual void
  writePlotData(EBAMRCellData& a_output, int& a_icomp) override;

  /*!
    @brief Write checkpoint data into handle
    @param[inout] a_handle HDF5 file
    @param[in]    a_level Grid level
  */
#ifdef CH_USE_HDF5
  virtual void
  writeCheckpointLevel(HDF5Handle& a_handle, const int a_level) const override;
#endif

  /*!
    @brief Read checkpoint data from handle
    @param[inout] a_handle HDF5 file
    @param[in]    a_level Grid level
  */
#ifdef CH_USE_HDF5
  virtual void
  readCheckpointLevel(HDF5Handle& a_handle, const int a_level) override;
#endif

  /*!
    @brief Get output plot names
    @return Returns a list of plot variable names. 
  */
  virtual Vector<std::string>
  getPlotVariableNames() const override;

  /*!
    @brief Get number of output variables
    @return Returns number of plot variables that will be written to file in writePlotData
  */
  virtual int
  getNumberOfPlotVariables() const;

  /*!
    @brief Count number of photons in particle list
    @param[in] a_photons Particle list
  */
  virtual int
  countPhotons(const AMRParticles<Photon>& a_photons) const;

  /*!
    @brief Count number of outcast photons in particle list
    @param[in] a_photons Particle list
  */
  virtual int
  countOutcast(const AMRParticles<Photon>& a_photons) const;

  /*!
    @brief Get m_photons
    @return m_photons
  */
  virtual ParticleContainer<Photon>&
  getPhotons();

  /*!
    @brief Get bulk photons, i.e. photons absorbed on the mesh
    @return m_bulkPhotons
  */
  virtual ParticleContainer<Photon>&
  getBulkPhotons();

  /*!
    @brief Get eb Photons, i.e. photons absorbed on the EB
    @return m_ebPhotons
  */
  virtual ParticleContainer<Photon>&
  getEbPhotons();

  /*!
    @brief Get domain photons, i.e. photons absorbed on domain edges/faces
    @return m_domainPhotonsl
  */
  virtual ParticleContainer<Photon>&
  getDomainPhotons();

  /*!
    @brief Get source photons
    @return m_sourcePhotons
  */
  virtual ParticleContainer<Photon>&
  getSourcePhotons();

protected:
  /*!
    @brief Enum for interpreting how photons are generated when using fluid codes
  */
  enum class PhotonGeneration
  {
    Deterministic,
    Stochastic
  };

  /*!
    @brief Enum for adding flexibility in what the fluid source term contains. 
  */
  enum class SourceType
  {
    Number,
    PerVol,
    PerVolSecond,
    PerSecond
  };

  /*!
    @brief An enum for switching between various types of EB intersection algorithms when intersecting photons with the EB
    @details Raycast means ray-casting algorithm. Bisection means that the traveled path is divided into intervals and we apply a bisection algorithm
    for computing the intersection point. LSF is just ray-casting, but it uses the implicit function on the mesh rather than calling it directly. 
  */
  enum class IntersectionEB
  {
    Raycast,
    Bisection,
  };

  /*!
    @brief Turn on/off transparent boundaries
  */
  bool m_transparentEB;

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
  bool m_depositNumber;

  /*!
    @brief Switch for plotting numbers or densities
  */
  bool m_plotNumbers;

  /*!
    @brief Check if m_photons should be plotted
  */
  bool m_plotPhotons;

  /*!
    @brief Check if m_bulkPhotons should be plotted
  */
  bool m_plotBulkPhotons;

  /*!
    @brief Check if source_bulkPhotons should be plotted
  */
  bool m_plotSourcePhotons;

  /*!
    @brief Check if m_ebPhotons should be plotted
  */
  bool m_plotEBPhotons;

  /*!
    @brief Check if m_ebPhotons should be plotted
  */
  bool m_plotDomainPhotons;

  /*!
    @brief Number of Photons per pseudoPhoton
  */
  size_t m_maxPhotonsGeneratedPerCell;

  /*!
    @brief RNG seed
  */
  int m_seed;

  /*!
    @brief Halo size for particles
  */
  int m_haloBuffer;

  /*!
    @brief Bisection step for when we use the bisection intersection algorithm
  */
  Real m_bisectStep;

  /*!
    @brief Photon generation type
  */
  PhotonGeneration m_photoGenerationMethod;

  /*!
    @brief Source type
  */
  SourceType m_sourceType;

  /*!
    @brief
  */
  IntersectionEB m_intersectionEB;

  /*!
    @brief Deposition type
  */
  DepositionType m_deposition;

  /*!
    @brief Coarse-fine deposition strategy
  */
  CoarseFineDeposition m_coarseFineDeposition;

  /*!
    @brief Plot deposition type
  */
  DepositionType m_plotDeposition;

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
    @brief Draw photons in a cell and volume
    @param[in] a_source Source term, i.e. number of photons generated/x
    @param[in] a_volume Grid cell volume
    @param[in] a_dt     Time step
    @details The return result of this function depends on the input parameters to this class. 
  */
  size_t
  drawPhotons(const Real a_source, const Real a_volume, const Real a_dt);

  /*!
    @brief Mapping function for domain boundary conditions
    @param[in] a_dir  Coordinate direction
    @param[in] a_side Coordinate side
  */
  int
  domainBcMap(const int a_dir, const Side::LoHiSide a_side);

  /*!
    @brief Random exponential trial
    @param[in] a_mean Mean value in exponential distribution. 
  */
  Real
  randomExponential(const Real a_mean);

  /*!
    @brief This computes the "conservative" deposition, multiplied by kappa
    @param[out] a_phi       Mesh density. 
    @param[in] a_particles  Particles to be deposited
    @param[in] a_deposition Deposition method
    @param[in] a_coarseFineDeposition Coarse-fine deposition strategy.
  */
  void
  depositKappaConservative(EBAMRCellData&             a_phi,
                           ParticleContainer<Photon>& a_particles,
                           const DepositionType       a_deposition,
                           const CoarseFineDeposition a_coarseFineDeposition);

  /*!
    @brief Make the "non-conservative" kappa deposition
    @param[out] a_depositionNC     Non-conservative deposition in cut-cells
    @param[in]  a_depositionKappaC Conservative deposition
  */
  void
  depositNonConservative(EBAMRIVData& a_depositionNC, const EBAMRCellData& a_depositionKappaC);

  /*!
    @brief Make the hybrid deposition. Also compute the mass difference
    @param[inout] a_depositionH    On input this should be the conservative deposition
    @param[out]   a_massDifference Mass difference between conservative and non-conservative
    @param[in]    a_depositioNC    Non-conservative deposition
  */
  void
  depositHybrid(EBAMRCellData& a_depositionH, EBAMRIVData& a_massDifference, const EBAMRIVData& a_depositionNC);

  /*!
    @brief Increment redistribution objects with mass to be redistributed
    @param[in] a_massDifference Mass difference
  */
  void
  incrementRedist(const EBAMRIVData& a_massDifference);

  /*!
    @brief Do redistribution on all levels, no coarse-fine EB magic. 
    @param[in] a_phi Data to be redistributed (with whatever has been put in the redistribution registers)
  */
  void
  levelRedist(EBAMRCellData& a_phi);

  /*!
    @brief Increment the coarse-fine redist with mass difference (only for EB x CF)
    @param[in] m_massDifference Mass difference
  */
  void
  coarseFineIncrement(const EBAMRIVData& m_massDifference);

  /*!
    @brief Do coarse-fine redistribution
    @param[inout] Data to be redistributed (with whatever has been put in the redistribution registers)
  */
  void
  coarseFineRedistribution(EBAMRCellData& a_phi);

  /*!
    @brief Turn on/off transparent boundaries
  */
  void
  parseTransparentBoundaries();

  /*!
    @brief Parse the divergence computation, i.e. if we blend with non-conservative divergence or not
  */
  void
  parseDivergenceComputation();

  /*!
    @brief Parse pseudophotons, i.e. parse how many photons can be generated per cell. 
  */
  void
  parsePseudoPhotons();

  /*!
    @brief Parse photogeneration type
  */
  void
  parsePhotoGeneration();

  /*!
    @brief Parse whether instantaneous propagation or not
  */
  void
  parseInstantaneous();

  /*!
    @brief Parse source term type
  */
  void
  parseSourceType();

  /*!
    @brief Parse EB intersection algorithm
  */
  void
  parseIntersectionEB();

  /*!
    @brief Parse deposition method
  */
  void
  parseDeposition();

  /*!
    @brief Parse plot variables
  */
  void
  parsePlotVariables();
};

#include <CD_NamespaceFooter.H>

#endif
