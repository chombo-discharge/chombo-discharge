/*!
  @file CdrSolver.H
  @brief Declaration of an abstract class for evolving cdr equations
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _CDR_SOLVER_
#define _CDR_SOLVER_

#include "computational_geometry.H"
#include "CD_AmrMesh.H"
#include "cdr_species.H"

#include <random>
#include <time.h>
#include <chrono>

#include "CD_NamespaceHeader.H"

/*!
  @brief Boundary condition types for the CDR solvers
*/
namespace cdr_bc {
  enum which_bc {
    external = 0,
    wall     = 1,
    outflow  = 2,
    extrap   = 3,
  };
};

/*!
  @brief Abstract class for cdr-equations. This is (almost) a whole solver. 
*/
class CdrSolver {
public:

  /*!
    @brief Default constructor. Stay away.
  */
  CdrSolver();

  /*!
    @brief Constructor
  */
  virtual ~CdrSolver();

  /*!
    @brief Parse class options from the command line
  */
  virtual void parseOptions() = 0;

  /*!
    @brief Parse runtime options
  */
  virtual void parseRuntimeOptions() = 0;

  /*!
    @brief Implicit diffusion TGA advance without source terms
  */
  virtual void advanceTGA(EBAMRCellData& a_newPhi, const EBAMRCellData& a_oldPhi, const Real a_dt) = 0;

  /*!
    @brief Implicit diffusion TGA advance with source terms
  */
  virtual void advanceTGA(EBAMRCellData&       a_newPhi,
			   const EBAMRCellData& a_oldPhi,
			   const EBAMRCellData& a_source,
			   const Real           a_dt) = 0;

  /*!
    @brief Implicit diffusion Euler advance without source terms
  */
  virtual void advanceEuler(EBAMRCellData& a_newPhi, const EBAMRCellData& a_oldPhi, const Real a_dt) = 0;

  /*!
    @brief Implicit diffusion Euler advance with source terms
  */
  virtual void advanceEuler(EBAMRCellData&       a_newPhi,
			     const EBAMRCellData& a_oldPhi,
			     const EBAMRCellData& a_source,
			     const Real           a_dt) = 0;

  /*!
    @brief Compute div(nv - D*grad(n)) explicitly
    @note a_phi is non-const because ghost cells will be re-filled
  */
  virtual void computeDivJ(EBAMRCellData& a_divJ,
			    EBAMRCellData& a_phi,
			    const Real     a_extrapDt,
			    const bool     a_ebFlux = true) = 0;

  /*!
    @brief Compute div(v*phi) explicitly
    @note a_phi is non-const because ghost cells will be re-filled
  */
  virtual void computeDivF(EBAMRCellData&  a_divF,
			    EBAMRCellData&  a_phi,
			    const Real      a_extrapDt,
			    const bool      a_ebFlux = true) = 0;

  /*!
    @brief Compute div(D*grad(phi)) explicitly
    @note a_phi is non-const because ghost cells will be re-filled
  */
  virtual void computeDivD(EBAMRCellData& a_divD,
			    EBAMRCellData& a_phi,
			    const bool     a_ebFlux = false) = 0;

  /*!
    @brief Compute div(G) where G is a general face-centered flux on face centers and EB centers
    @param[in]    a_divG div(G). This includes redistribution
    @param[inout] a_G Must contain face-centered fluxes on input. Contains face-centroid fluxes on output. 
    @param[in]    a_ebFlux Flux on the EB centroids
    @note If you also want charge injection through the domain sides, you must put those in a_G before entering this routine
  */
  virtual void computeDivG(EBAMRCellData& a_divG, EBAMRFluxData& a_G, const EBAMRIVData& a_ebFlux);

  /*!
    @brief Random source for diffusion
  */
  virtual void gwnDiffusionSource(EBAMRCellData& a_noiseSource, const EBAMRCellData& a_cellPhi);

  /*!
    @brief Smoothed face average
  */
  virtual void smoothHeavisideFaces(EBAMRFluxData& a_facePhi, const EBAMRCellData& a_cellPhi);

  /*!
    @brief Gaussian noise field
  */
  virtual void fillGwn(EBAMRFluxData& a_noise, const Real a_mean);

  /*!
    @brief Average velocities to faces
  */
  virtual void averageVelocityToFaces();

  /*!
    @brief Pre-regrid
  */
  virtual void preRegrid(const int a_lbase, const int a_oldFinestLevel);

  /*!
    @brief Deallocate internal storage
  */
  virtual void deallocateInternals();

  /*!
    @brief Register operators for AMR operations
  */
  virtual void registerOperators();

  /*!
    @brief Sanity check
  */
  virtual void sanityCheck();

  /*!
    @brief Set the Realm
  */
  virtual void setRealm(const std::string a_realm);

  /*!
    @brief Set species
  */
  virtual void setSpecies(const RefCountedPtr<cdr_species> a_species);
  
  /*!
    @brief Set computational geometry
  */
  virtual void setComputationalGeometry(const RefCountedPtr<computational_geometry> a_computationalGeometry);

  /*!
    @brief Set the amr object
  */
  virtual void setAmr(const RefCountedPtr<AmrMesh>& a_amr);

  /*!
    @brief Set domain bc type
  */
  virtual void setDomainBc(const cdr_bc::which_bc a_bctype);

  /*!
    @brief Set phase
    @details This must be done BEFORE callilng setComputationalGeometry
  */
  virtual void setPhase(phase::which_phase a_phase = phase::gas);

  /*!
    @brief Set verbosity.
  */
  virtual void setVerbosity(const int a_verbosity);

  /*!
    @brief Get current time
  */
  virtual void setTime(const int a_step, const Real a_time, const Real a_dt);

  /*!
    @brief Set velocity
  */
  virtual void setVelocity(const EBAMRCellData& a_velo);

  /*!
    @brief Set constant velocity
  */
  virtual void setVelocity(const RealVect a_velo);

  /*!
    @brief Set diffusion coefficients
  */
  virtual void setDiffusionCoefficient(const EBAMRFluxData& a_diffusionCoefficient, const EBAMRIVData& a_ebDiffusionCoefficient);

  /*!
    @brief Set constant
  */
  virtual void setDiffusionCoefficient(const Real a_diffusionCoefficient);

  /*!
    @brief Set source terms
  */
  virtual void setSource(const EBAMRCellData& a_source);

  /*!
    @brief Set constant source terms
  */
  virtual void setSource(const Real a_source);

  /*!
    @brief Set the EB flux
  */
  virtual void setEbFlux(const EBAMRIVData& a_ebFlux);

  /*!
    @brief Set the EB flux
  */
  virtual void setEbFlux(const Real a_ebFlux);

  /*!
    @brief Set the domain flux
  */
  virtual void setDomainFlux(const Real a_domainFlux);

  /*!
    @brief Fill internal state with initial data
  */
  virtual void initialData();

  /*!
    @brief Fill initial data from a distribution function
  */
  virtual void initialDataDistribution();

  /*!
    @brief Fill initial data from particles
  */
  virtual void initialDataParticles();

  /*!
    @brief Inject the EB flux. 
    @details This routine is equivalent to an update phi^(k+1) = phi^k - dt*div(G)
    where G is non-zero only on the EB
  */
  virtual void injectEbFlux(EBAMRCellData& a_phi, const EBAMRIVData& a_ebFlux, const Real a_dt);

  /*!
    @brief Allocate internal storage
  */
  virtual void allocateInternals();

  /*!
    @brief Write plot file
  */
  virtual void writePlotFile();

  /*!
    @brief Write output data to a_output
  */
  virtual void writePlotData(EBAMRCellData& a_output, int& a_comp);

  /*!
    @brief Write data to output. Convenience function. 
  */
  virtual void writeData(EBAMRCellData& a_output, int& a_comp, const EBAMRCellData& a_data, const bool a_interp);

  /*!
    @brief Write checkpoint data into handle
  */
  virtual void writeCheckpointLevel(HDF5Handle& a_handle, const int a_level) const;

  /*!
    @brief Read checkpoint data from handle
  */
  virtual void readCheckpointLevel(HDF5Handle& a_handle, const int a_level);

  /*!
    @brief Regrid internals
  */
  virtual void regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel);

  /*!
    @brief Set which variables to output
  */
  virtual void set_plot_variables();

  /*!
    @brief Reset redistribution weights
  */
  virtual void setRedistWeights(const EBAMRCellData& a_phi);

  /*!
    @brief Get the Realm
  */
  virtual const std::string getRealm() const;

  /*!
    @brief Get solver name
  */
  virtual std::string getName();

  /*!
    @brief Get output plot names
  */
  virtual Vector<std::string> get_plotVariableNames() const;

  /*!
    @brief Get the ghost cell requirement
  */
  virtual int queryGhost() const = 0;

  /*!
    @brief Get number of output fields
  */
  virtual int getNumberOfPlotVariables() const;

  /*!
    @brief Compute the largest possible diffusive time step (for explicit methods)
  */
  virtual Real computeAdvectionDt();

  /*!
    @brief Compute the largest possible diffusive time step (for explicit methods)
  */
  virtual Real computeDiffusionDt();

  /*!
    @brief Compute the largest possible diffusive time step (for explicit methods)
  */
  virtual Real computeAdvectionDiffusionDt();

  /*!
    @brief Compute the largest possible source time step (for explicit methods
  */
  virtual Real computeSourceDt(const Real a_max, const Real a_tolerance);
  
  /*!
    @brief Compute the total mass in m_phi
  */
  virtual Real computeMass();

  /*!
    @brief Compute the total charge in m_phi
  */
  virtual Real computeCharge();

  /*!
    @brief Diffusive or not
  */
  virtual bool isDiffusive();

  /*!
    @brief Mobile or not
  */
  virtual bool isMobile();

  /*!
    @brief Check if we can extrapolate
  */
  virtual bool extrapolateSourceTerm() const;
  
  /*!
    @brief Get the internal state
  */
  virtual EBAMRCellData& getPhi();

  /*!
    @brief Get the source term
  */
  virtual EBAMRCellData& getSource();

  /*!
    @brief Get the cell-centered velocity
  */
  virtual EBAMRCellData& getCellCenteredVelocity();

  /*!
    @brief Get the face-centered velocities
  */
  virtual EBAMRFluxData& getFaceCenteredVelocity();

  /*!
    @brief Get the eb-centered velocities
  */
  virtual EBAMRIVData& getEbCenteredVelocity();

  /*!
    @brief Get the face-centered diffusion coefficient
  */
  virtual EBAMRFluxData& getFaceCenteredDiffusionCoefficient();

  /*!
    @brief Get the EB-centered diffusion coefficient
  */
  virtual EBAMRIVData& get_diffco_eb();

  /*!
    @brief Get the eb flux data holder
  */
  virtual EBAMRIVData& getEbFlux();

  /*!
    @brief Get the domain flux data holder
  */
  virtual EBAMRIFData& getDomainFlux();


protected:
  
  /*!
    @brief Phase
  */
  RefCountedPtr<EBIndexSpace> m_ebis;

  /*!
    @brief Species
  */
  RefCountedPtr<cdr_species> m_species;
  
  /*!
    @brief Computational geometry
  */
  RefCountedPtr<computational_geometry> m_computationalGeometry;

  /*!
    @brief AMR; needed for grid stuff
  */
  RefCountedPtr<AmrMesh> m_amr;

  /*!
    @brief Stencils for interpolating flux to centroids
  */
  Vector<RefCountedPtr<LayoutData<BaseIFFAB<FaceStencil> > > > m_interpStencils[SpaceDim];

  /*!
    @brief Stencils for computing the non-conservative divergence
  */
  Vector<RefCountedPtr<LayoutData<BaseIVFAB<VoFStencil> > > > m_stencils_nc;

  /*!
    @brief Storage for flux interpolant
  */
  Vector<RefCountedPtr<LevelData<BaseIFFAB<Real> > > > m_interpolant[SpaceDim];

  /*!
    @brief Cells on which to compute the flux interpolant
  */
  Vector<RefCountedPtr<LayoutData<IntVectSet> > > m_interpSets[SpaceDim];

  /*!
    @brief Phase
  */
  phase::which_phase m_phase;

  /*!
    @brief Domain bc condition
  */
  cdr_bc::which_bc m_domainBc;

  /*!
    @brief Name
  */
  std::string m_name;

  /*!
    @brief Class name. 
  */
  std::string m_className;

  /*!
    @brief Realm
  */
  std::string m_Realm;

  /*!
    @brief State data (i.e. the ion density)
  */
  EBAMRCellData m_phi;

  /*!
    @brief Source term 
  */
  EBAMRCellData m_source;

  /*!
    @brief Cell-centered velocities
  */
  EBAMRCellData m_cellVelocity;

  /*!
    @brief Holder for face centered states
  */
  EBAMRFluxData m_faceStates;

  /*!
    @brief Scratch storage. 
  */
  EBAMRCellData m_scratch;

  /*!
    @brief Scratch storage
  */
  EBAMRFluxData m_scratchFluxOne;

  /*!
    @brief Scratch storage
  */
  EBAMRFluxData m_scratchFluxTwo;

  /*!
    @brief Scratch storage
  */
  EBAMRIVData m_nonConservativeDivG;

  /*!
    @brief Scratch storage
  */
  EBAMRIVData m_massDifference;

  /*!
    @brief Zero flux
  */
  EBAMRIVData m_EbZero;

  /*!
    @brief Cached state vector for regrid. DO NOT TOUCH!
  */
  EBAMRCellData m_cachePhi;

  /*!
    @brief Cached source term for regrids. DO NOT TOUCH!
  */
  EBAMRCellData m_cacheSource;

  /*!
    @brief Face-centered velocities (only normal components)
  */
  EBAMRFluxData m_faceVelocity;

  /*!
    @brief EB-centered velocities
  */
  EBAMRIVData m_ebVelocity;

  /*!
    @brief Flux through the embedded boundary
  */
  EBAMRIVData m_ebFlux;

  /*!
    @brief Domain flux
  */
  EBAMRIFData m_domainFlux;

  /*!
    @brief Dummy data
  */
  EBAMRCellData m_aCoefficient;

  /*!
    @brief Diffusion coefficients (can
  */
  EBAMRFluxData m_faceCenteredDiffusionCoefficient;

  /*!
    @brief Diffusion coefficients (can
  */
  EBAMRIVData m_faceCenteredDiffusionCoefficient_eb;

  /*!
    @brief Verbosity
  */
  int m_verbosity;

  /*!
    @brief Time step
  */
  int m_timeStep;

  /*!
    @brief Current time
  */
  Real m_time;

  /*!
    @brief Time step
  */
  Real m_dt;

  /*!
    @brief Flag for mass-weighted redistribution
  */
  bool m_useMassWeightedRedistribution;

  /*!
    @brief Flag for blending the hybrid divergence
  */
  bool m_blendConservation;

  /*!
    @brief Also solve for diffusion
  */
  bool m_isDiffusive;

  /*!
    @brief Also solve for advection/convection
  */
  bool m_isMobile;
  
  /*!
    @brief Turn on/off source terms when time-extrapolating
  */
  bool m_extrapolateSourceTerm;

  /*!
    @brief Output state
  */
  bool m_plotPhi;

  /*!
    @brief Output velocities
  */
  bool m_plotVelocity;

  /*!
    @brief Output diffusion coefficients
  */
  bool m_plotDiffusionCoefficient;

  /*!
    @brief Output EB fluxes
  */
  bool m_plotEbFlux;

  /*!
    @brief Output source term
  */
  bool m_plotSource;

  /*!
    @brief Plot numbers or densities
  */
  bool m_plotNumbers;

  /*!
    @brief RNG engine
  */
  std::mt19937_64 m_rng;

  /*!
    @brief RNG seed
  */
  int m_seed;

  /*!
    @brief Average cell-centered velocities to face centers
  */
  virtual void averageVelocityToFaces(EBAMRFluxData& a_faceVelocity, const EBAMRCellData& a_cellVelocity);

  /*!
    @brief Set up face-centered advection flux
  */
  virtual void computeFlux(EBAMRFluxData&       a_flux,
			    const EBAMRFluxData& a_facePhi,
			    const EBAMRFluxData& a_faceVelocity,
			    const EBAMRIFData&   a_domainFlux);

  /*!
    @brief Set up face-centered advection flux
  */
  virtual void computeFlux(LevelData<EBFluxFAB>&              a_flux,
			    const LevelData<EBFluxFAB>&        a_facePhi,
			    const LevelData<EBFluxFAB>&        a_faceVelocity,
			    const LevelData<DomainFluxIFFAB>&  a_domainFlux,
			    const int                          a_lvl);

  /*!
    @brief Compute the face-centered diffusion flux
  */
  virtual void computeDiffusionFlux(EBAMRFluxData& a_flux, const EBAMRCellData& a_phi);

  /*!
    @brief Compute the face-centered difufsion flux
  */
  virtual void computeDiffusionFlux(LevelData<EBFluxFAB>& a_flux, const LevelData<EBCellFAB>& a_phi, const int a_lvl);

  /*!
    @brief Compute conservative divergence from fluxes. 
  */
  virtual void conservativeDivergenceNoKappaDivision(EBAMRCellData& a_conservativeDivergence, EBAMRFluxData& a_flux, const EBAMRIVData& a_ebFlux);

  /*!
    @brief Compute the conservative divergence with the assumption that only the eb-flux is non-zero
  */
  virtual void conservativeDivergenceNoKappaDivision_eb(EBAMRCellData& a_conservativeDivergence, const EBAMRIVData& a_ebFlux);


  /*!
    @brief Compute the non-conservative divergence
  */
  virtual void nonconservativeDivergenceNoKappaDivision(EBAMRIVData& a_nonConservativeDivergence, const EBAMRCellData& a_divG);

  /*!
    @brief Use the non-conservative divergence to make the conservative divergence hold the hybrid divergence
  */
  virtual void hybridDivergence(EBAMRCellData&       a_divF,
				 EBAMRIVData&         a_massDifference,
				 const EBAMRIVData&   a_NonConservativeDivergenceStencil);

  
  /*!
    @brief Conservative divergence over regular cells
  */
  virtual void consdiv_regular(LevelData<EBCellFAB>& a_divJ, const LevelData<EBFluxFAB>& a_flux, const int a_lvl);

  /*!
    @brief Make the hybrid divergence. On the way in, a_hybridDivergence must hold the conservative divergence
  */
  virtual void hybridDivergence(LevelData<EBCellFAB>&              a_hybridDivergence,
				 LevelData<BaseIVFAB<Real> >&       a_massDifference,
				 const LevelData<BaseIVFAB<Real> >& a_nonConservativeDivergence,
				 const int                          a_lvl);

  /*!
    @brief Compute flux interpolant
  */
  virtual void setupFluxInterpolant(const LevelData<EBFluxFAB>& a_flux, const int a_lvl);

  /*!
    @brief Interpolate flux to centroids
  */
  virtual void interpolateFluxToFaceCentroids(LevelData<EBFluxFAB>& a_flux, const int a_lvl);

  /*!
    @brief Compute div(F) on irregular cells
  */
  virtual void computeDivergenceIrregular(LevelData<EBCellFAB>&              a_divG,
				  const LevelData<BaseIVFAB<Real> >& a_ebFlux,
				  const int                          a_lvl);

  /*!
    @brief Increment flux registers. 
  */
  virtual void incrementFluxRegister(const EBAMRFluxData& a_facePhi, const EBAMRFluxData& a_faceVelocity);

  /*!
    @brief Increment flux registers
  */
  virtual void incrementFluxRegister(const EBAMRFluxData& a_flux);

  /*!
    @brief Reset the flux register
  */
  virtual void reset_flux_register();

  /*!
    @brief Increment the coarse-fine redist with mass difference (only for EB x CF)
  */
  virtual void coarseFineIncrement(const EBAMRIVData& m_massDifference);

  /*!
    @brief Do mass redistribution
  */
  virtual void hyperbolicRedistribution(EBAMRCellData& a_stableDivergence, const EBAMRIVData&   a_massDifference);

  /*!
    @brief Increment redistribution flux registers
  */
  virtual void incrementRedistFlux();

  /*!
    @brief Increment redistribution objects with mass to be redistributed
  */
  virtual void incrementRedist(const EBAMRIVData& a_massDifference);

  /*!
    @brief Do coarse-fine redistribution
  */
  virtual void coarseFineRedistribution(EBAMRCellData& a_phi);

  /*!
    @brief Reflux state
  */
  virtual void reflux(EBAMRCellData& a_phi);

  /*!
    @brief Set ebis
  */
  virtual void setEbIndexSpace(const RefCountedPtr<EBIndexSpace>& a_ebis);

  /*!
    @brief Define stencils
  */
  virtual void defineInterpolationStencils();

  /*!
    @brief Define interpolant
  */
  virtual void defineInterpolant();

  /*!
    @brief Parses domain BC options
  */
  virtual void parseDomainBc();    // Parses domain BC options

  /*!
    @brief Parses plot variables
  */
  virtual void parsePlotVariables();    // Parses plot variables

  /*!
    @brief Parse random number generator
  */
  virtual void parseRngSeed();

  /*!
    @brief Parse plot mode
  */
  virtual void parsePlotMode();

  /*!
    @brief Parse the conservation
  */
  virtual void parseDivergenceComputation();

  /*!
    ®brief Parse source term stuff
  */
  virtual void parseExtrapolateSourceTerm();
};
#include "CD_NamespaceFooter.H"

#endif
