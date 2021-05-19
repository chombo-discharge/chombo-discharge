/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_CdrSolver.H
  @brief  Declaration of an abstract class for solving scalar convection-diffusion-reaction problems. 
  @author Robert Marskar
*/

#ifndef CD_CdrSolver_H
#define CD_CdrSolver_H

// Std includes
#include <random>
#include <time.h>
#include <chrono>

// Our includes
#include <CD_AmrMesh.H>
#include <CD_ComputationalGeometry.H>
#include <CD_CdrSpecies.H>
#include <CD_CdrBc.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Abstract class for cdr-equations. 
*/
class CdrSolver {
public:

  /*!
    @brief Default constructor. 
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
    @brief Implicit diffusion TGA advance without source term.
    @param[inout] a_newPhi Solution at time t + dt
    @param[in]    a_oldPhi Solution at time t
    @param[in]    a_dt     Time step
    @note This version does not use a source term. 
  */
  virtual void advanceTGA(EBAMRCellData& a_newPhi, const EBAMRCellData& a_oldPhi, const Real a_dt) = 0;

  /*!
    @brief Implicit diffusion TGA advance with source term.
    @param[inout] a_newPhi Solution at time t + dt
    @param[in]    a_oldPhi Solution at time t
    @param[in]    a_source Source (at t+dt/2)
    @param[in]    a_dt     Time step
  */
  virtual void advanceTGA(EBAMRCellData&       a_newPhi,
			  const EBAMRCellData& a_oldPhi,
			  const EBAMRCellData& a_source,
			  const Real           a_dt) = 0;

  /*!
    @brief Implicit diffusion Euler advance without source term.
    @param[inout] a_newPhi Solution at time t + dt
    @param[in]    a_oldPhi Solution at time t
    @param[in]    a_dt     Time step
  */
  virtual void advanceEuler(EBAMRCellData& a_newPhi, const EBAMRCellData& a_oldPhi, const Real a_dt) = 0;

  /*!
    @brief Implicit diffusion Euler advance with source term.
    @param[inout] a_newPhi Solution at time t + dt
    @param[in]    a_oldPhi Solution at time t
    @param[in]    a_source Source term. 
    @param[in]    a_dt     Time step
    @details This solves the implicit diffusino equation equation a_newPhi - a_oldPhi = dt*Laplacian(a_newPhi) + a_source. 
  */
  virtual void advanceEuler(EBAMRCellData&       a_newPhi,
			    const EBAMRCellData& a_oldPhi,
			    const EBAMRCellData& a_source,
			    const Real           a_dt) = 0;

  /*!
    @brief Compute div(J) where J = nV - D*grad(n) explicitly
    @param[out] a_divJ      Divergence term, including redistribution magic. 
    @param[in]  a_phi       Cell-centered state
    @param[in]  a_extrapDt  Extrapolation in time, i.e. shifting of div(J) towards e.g. half time step. Only affects the advective term. 
    @param[in]  a_useEbFlux If true, the embedded boundary flux will be included in div(J)
    @note a_phi is non-const because ghost cells will be re-filled
  */
  virtual void computeDivJ(EBAMRCellData& a_divJ,
			   EBAMRCellData& a_phi,
			   const Real     a_extrapDt,
			   const bool     a_useEbFlux = true) = 0;

  /*!
    @brief Compute div(v*phi) explicitly
    @param[out] a_divF      Divergence term, including redistribution magic. 
    @param[in]  a_phi       Cell-centered state
    @param[in]  a_extrapDt  Extrapolation in time, i.e. shifting of div(F) towards e.g. half time step. Only affects the advective term. 
    @param[in]  a_useEbFlux If true, the embedded boundary flux will be included in div(F)
    @note a_phi is non-const because ghost cells will be re-filled
  */
  virtual void computeDivF(EBAMRCellData&  a_divF,
			   EBAMRCellData&  a_phi,
			   const Real      a_extrapDt,
			   const bool      a_ebFlux = true) = 0;

  /*!
    @brief Compute div(D*grad(phi)) explicitly
    @param[out] a_divF      Divergence term, including redistribution magic. 
    @param[in]  a_phi       Cell-centered state
    @param[in]  a_useEbFlux If true, the embedded boundary flux will be included in div(D)
    @note a_phi is non-const because ghost cells will be re-filled
  */
  virtual void computeDivD(EBAMRCellData& a_divD,
			   EBAMRCellData& a_phi,
			   const bool     a_ebFlux = false) = 0;

  /*!
    @brief Compute div(G) where G is a general face-centered flux on face centers and EB centers
    @param[in]    a_divG div(G). This includes redistribution
    @param[inout] a_G Vector field which contains face-centered fluxes on input. Contains face-centroid fluxes on output. 
    @param[in]    a_ebFlux Flux on the EB centroids
    @note If you also want charge injection through the domain sides, you must put those in a_G before entering this routine
  */
  virtual void computeDivG(EBAMRCellData& a_divG, EBAMRFluxData& a_G, const EBAMRIVData& a_ebFlux);

  /*!
    @brief Compute a random gaussian white noise source term. 
    @param[out] a_noiseSource Source term
    @param[in]  a_cellPhi     Cell-centered states
  */
  virtual void gwnDiffusionSource(EBAMRCellData& a_noiseSource, const EBAMRCellData& a_cellPhi);

  /*!
    @brief Use Heaviside smoothing for computing face-centered states. 
    @param[out] a_facePhi Face centered state
    @param[in]  a_cellPhi Cell-centered states
    @note This is only intended to be used with FHD functionality. 
  */
  virtual void smoothHeavisideFaces(EBAMRFluxData& a_facePhi, const EBAMRCellData& a_cellPhi);

  /*!
    @brief Gaussian noise field
    @param[out] a_noise Gaussian white nosie
    @param[out] a_sigma Standard deviation
  */
  virtual void fillGwn(EBAMRFluxData& a_noise, const Real a_sigma);

  /*!
    @brief Average velocities to faces
  */
  virtual void averageVelocityToFaces();

  /*!
    @brief Perform pre-regrid operations. 
    @param[in] a_lbase          Coarsest level that changed.
    @param[in] a_oldFinestLevel Previous finest level.
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
    @brief Do a sanity check
  */
  virtual void sanityCheck();

  /*!
    @brief Set the realm for this solver.
    @param[in] a_realm Realm identifier
  */
  virtual void setRealm(const std::string a_realm);

  /*!
    @brief Set species
    @param[in] a_species Species. 
  */
  virtual void setSpecies(const RefCountedPtr<CdrSpecies> a_species);
  
  /*!
    @brief Set computational geometry
    @param[in] a_computationalGeometry The computational geometry.
  */
  virtual void setComputationalGeometry(const RefCountedPtr<ComputationalGeometry> a_computationalGeometry);

  /*!
    @brief Set the amr object
    @param[in] a_amr AmrMesh object. 
  */
  virtual void setAmr(const RefCountedPtr<AmrMesh>& a_amr);

  /*!
    @brief Set domain bc type.
    @param[in] a_bcType Boundary condition type. 
  */
  virtual void setDomainBc(const CdrBc a_bcType);

  /*!
    @brief Set phase.
    @param[in] a_phase Phse.
  */
  virtual void setPhase(phase::which_phase a_phase = phase::gas);

  /*!
    @brief Set verbosity.
    @param[in] a_verbosity Verbosity level. 
  */
  virtual void setVerbosity(const int a_verbosity);

  /*!
    @brief Set the time for this solver. 
    @param[in] a_step Time step number
    @param[in] a_time Time (in seconds)
    @param[in] a_dt   Time step increment
  */
  virtual void setTime(const int a_step, const Real a_time, const Real a_dt);

  /*!
    @brief Set velocity
    @param[in] a_velo Cell centered velocities. 
    @details This copies a_velo onto m_cellCenteredVelocities. 
  */
  virtual void setVelocity(const EBAMRCellData& a_velo);

  /*!
    @brief Set constant velocity
    @param[in] a_velo Velocity
    @details This sets the cell centered velocity to be a_velo in every grid cell. 
  */
  virtual void setVelocity(const RealVect a_velo);

  /*!
    @brief Set diffusion coefficients
    @param[in] a_diffusionCoefficient   Face-centered diffusion coefficient. 
    @param[in] a_ebDiffusionCoefficient EB-centered diffusion coefficient. 
  */
  virtual void setDiffusionCoefficient(const EBAMRFluxData& a_diffusionCoefficient, const EBAMRIVData& a_ebDiffusionCoefficient);

  /*!
    @brief Set diffusion coefficient to be constant
    @param[in] a_diffusionCoefficient Global diffusion coefficient. 
    @details This sets the diffusion coefficient to be constant on both faces and EB faces. 
  */
  virtual void setDiffusionCoefficient(const Real a_diffusionCoefficient);

  /*!
    @brief Set source terms
    @param[in] a_source Source term.
    @details This sets m_source to be equal to a_source
  */
  virtual void setSource(const EBAMRCellData& a_source);

  /*!
    @brief Set constant source terms
    @param[in] a_source Source term.
  */
  virtual void setSource(const Real a_source);

  /*!
    @brief Set the EB flux
    @param[in] a_ebFlux Flux on EB centroids. 
  */
  virtual void setEbFlux(const EBAMRIVData& a_ebFlux);

  /*!
    @brief Set the EB flux.
    @param[in] a_ebFlux Flux on EB centroids. 
    @details This sets a constant EB-flux on all EBs. 
  */
  virtual void setEbFlux(const Real a_ebFlux);

  /*!
    @brief Set the domain flux
    @param[in] a_domainFlux Domain flux
    @details This sets a constant domain flux on all domain faces. 
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
    @brief Inject the EB flux. This includes redistribution magic. 
    @param[inout] a_phi    Cell-centered states
    @param[in]    a_ebFlux EB flux. 
    @param[in]    a_dt     Time step. 
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
    @param[inout] a_output Output data holder.
    @param[inout] a_icomp  Starting component where this solver begins writing the output. 
  */
  virtual void writePlotData(EBAMRCellData& a_output, int& a_comp);

  /*!
    @brief Write data to output. Convenience function. 
    @param[inout] a_output Output data holder.
    @param[inout] a_icomp  Starting component where this solver begins writing the output. 
    @param[in]    a_data   Data to write. 
    @param[in]    a_interp If true, a_data will be interpolated to cell centroids before writing to a_output. 
  */
  virtual void writeData(EBAMRCellData& a_output, int& a_comp, const EBAMRCellData& a_data, const bool a_interp);

  /*!
    @brief Write checkpoint data into handle
    @paramo[out] a_handle HDF5 handle.
    @param[in]   a_level Grid level
  */
  virtual void writeCheckpointLevel(HDF5Handle& a_handle, const int a_level) const;

  /*!
    @brief Read checkpoint data from handle
    @param[in] a_handle HDF5 handle.
    @param[in] const int a_level Grid level
  */
  virtual void readCheckpointLevel(HDF5Handle& a_handle, const int a_level);

  /*!
    @brief Regrid this solver. 
    @param[in] a_lmin           Coarsest level where grids did not change. 
    @param[in] a_oldFinestLevel Previous finest AMR level
    @param[in] a_newFinestLevel New finest AMR level. 
    @details This linearly interpolates (with limiters) m_potential to the new grids and recomputes the electric field (from the interpolated potential). 
  */
  virtual void regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel);

  /*!
    @brief Set redistribution weights
    @param[in] a_weights Weights to use for redistribution (e.g. mass)
  */
  virtual void setRedistWeights(const EBAMRCellData& a_weights);

  /*!
    @brief Get the realm where this solver is registered. 
  */
  virtual const std::string getRealm() const;

  /*!
    @brief Get solver name
    @note Not necessarily equal to class name (we can have many CdrSolvers instantiated). 
  */
  virtual std::string getName();

  /*!
    @brief Get output plot names
  */
  virtual Vector<std::string> getPlotVariableNames() const;

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
    @param[in] a_max       Maximum value of m_phi
    @param[in] a_tolerance Tolerance
    @details This computes dt = phi/source, but only for cells where phi lies within a_tolerance*a_max (and source > 0)
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
    @brief Get the cell-centered phi
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
  virtual EBAMRIVData& getEbCenteredDiffusionCoefficient();

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
  RefCountedPtr<CdrSpecies> m_species;
  
  /*!
    @brief Computational geometry
  */
  RefCountedPtr<ComputationalGeometry> m_computationalGeometry;

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
  CdrBc m_domainBc;

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
  std::string m_realm;

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
  EBAMRIVData m_ebZero;

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
  EBAMRCellData m_aCoef;

  /*!
    @brief Diffusion coefficients (can
  */
  EBAMRFluxData m_faceCenteredDiffusionCoefficient;

  /*!
    @brief Diffusion coefficients (can
  */
  EBAMRIVData m_ebCenteredDiffusionCoefficient;

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
    @param[out] a_faceVelocity Face-centered velocities
    @param[in]  a_cellVelocity Cell-centered velocities
  */
  virtual void averageVelocityToFaces(EBAMRFluxData& a_faceVelocity, const EBAMRCellData& a_cellVelocity);

  /*!
    @brief Set up face-centered advection flux.
    @param[out] a_flux         Face-centered fluxes
    @param[in]  a_facePhi      Face-centered state
    @param[in]  a_faceVelocity Face-centered velocities
    @param[in]  a_domainFlux   Domain flux
  */
  virtual void computeFlux(EBAMRFluxData&       a_flux,
			   const EBAMRFluxData& a_facePhi,
			   const EBAMRFluxData& a_faceVelocity,
			   const EBAMRIFData&   a_domainFlux);

  /*!
    @brief Set up face-centered advection flux on a grid level
    @param[out] a_flux         Face-centered fluxes
    @param[in]  a_facePhi      Face-centered state
    @param[in]  a_faceVelocity Face-centered velocities
    @param[in]  a_domainFlux   Domain flux
    @param[in]  a_level        Grid level
  */
  virtual void computeFlux(LevelData<EBFluxFAB>&              a_flux,
			   const LevelData<EBFluxFAB>&        a_facePhi,
			   const LevelData<EBFluxFAB>&        a_faceVelocity,
			   const LevelData<DomainFluxIFFAB>&  a_domainFlux,
			   const int                          a_lvl);

  /*!
    @brief Compute the face-centered diffusion flux
    @param[out] a_flux Face-centered flux
    @param[in]  a_phi  Cell centered phi
  */
  virtual void computeDiffusionFlux(EBAMRFluxData& a_flux, const EBAMRCellData& a_phi);

  /*!
    @brief Compute the face-centered diffusion flux
    @param[out] a_flux Face-centered flux
    @param[in]  a_phi  Cell centered phi
    @param[in]  a_lvl  Grid level
  */
  virtual void computeDiffusionFlux(LevelData<EBFluxFAB>& a_flux, const LevelData<EBCellFAB>& a_phi, const int a_lvl);

  /*!
    @brief Compute conservative divergence from fluxes. 
    @param[out] a_conservativeDivergence Conservative divergence computed as div(F) using finite volumes. Not divided by kappa. 
    @param[in]  a_flux                   Face-centered fluxes. Includes domain fluxes.
    @param[in]  a_ebFlux                 EB flux. 
  */
  virtual void conservativeDivergenceNoKappaDivision(EBAMRCellData& a_conservativeDivergence, EBAMRFluxData& a_flux, const EBAMRIVData& a_ebFlux);

  /*!
    @brief Compute the conservative divergence with the assumption that only the eb-flux is non-zero
    @param[out] a_conservativeDivergence Conservative divergence.
    @param[in]  a_ebFlux                 EB flux
  */
  virtual void conservativeDivergenceNoKappaDivisionOnlyEbFlux(EBAMRCellData& a_conservativeDivergence, const EBAMRIVData& a_ebFlux);

  /*!
    @brief Compute the non-conservative divergence
    @param[out] a_nonConservativeDivergence Non-conservative divergence.
    @param[in]  a_divG                      Conservative divergence.
  */
  virtual void nonConservativeDivergence(EBAMRIVData& a_nonConservativeDivergence, const EBAMRCellData& a_divG);

  /*!
    @brief Use the non-conservative divergence to make the conservative divergence hold the hybrid divergence
    @param[inout] a_hybridDivergence          On input, contains the conservative divergence (without kappa division). Contains hybrid divergence on output.
    @param[out]   a_massDifference            Mass difference between updating with hybrid divergence and true divergence
    @param[in]    a_nonConservativeDivergence Non-conservative divergence.
  */
  virtual void hybridDivergence(EBAMRCellData&       a_hybridDivergence,
				EBAMRIVData&         a_massDifference,
				const EBAMRIVData&   a_nonConservativeDivergence);

  /*!
    @brief Make the hybrid divergence. On the way in, a_hybridDivergence must hold the conservative divergence
    @param[inout] a_hybridDivergence          On input, contains the conservative divergence (without kappa division). Contains hybrid divergence on output.
    @param[out]   a_massDifference            Mass difference between updating with hybrid divergence and true divergence
    @param[in]    a_nonConservativeDivergence Non-conservative divergence.
    @param[in]    a_lvl                       Grid level. 
  */
  virtual void hybridDivergence(LevelData<EBCellFAB>&              a_hybridDivergence,
				LevelData<BaseIVFAB<Real> >&       a_massDifference,
				const LevelData<BaseIVFAB<Real> >& a_nonConservativeDivergence,
				const int                          a_lvl);
  
  /*!
    @brief Compute the conservative divergence over regular cells
    @param[out] a_divJ Conservative divergence. Not kappa-divided. 
    @param[in]  a_flux Face-centered fluxes
    @param[in]  a_lvl  Grid level
  */
  virtual void conservativeDivergenceRegular(LevelData<EBCellFAB>& a_divJ, const LevelData<EBFluxFAB>& a_flux, const int a_lvl);

  /*!
    @brief Compute flux interpolant
    @param[in] a_flux The flux to be interpolated to face centroids
    @param[in] a_lvl  Grid level
    @details This fills m_interpolant
  */
  virtual void setupFluxInterpolant(const LevelData<EBFluxFAB>& a_flux, const int a_lvl);

  /*!
    @brief Interpolate flux to centroids
    @param[out] a_flux Centroid fluxes
    @param[in]  a_lvl  Grid level
    @details This interpolates face fluxes to centroids and puts the result in a_flux and m_interpolant. 
  */
  virtual void interpolateFluxToFaceCentroids(LevelData<EBFluxFAB>& a_flux, const int a_lvl);

  /*!
    @brief Compute conservative divergence on irregular cells (not kappa divided)
    @param[out] a_divG   Conservative divergence (not kappa divided)
    @param[in]  a_ebFlux EB flux
    @param[in]  a_lvl    Grid level
  */
  virtual void computeDivergenceIrregular(LevelData<EBCellFAB>&              a_divG,
					  const LevelData<BaseIVFAB<Real> >& a_ebFlux,
					  const int                          a_lvl);

  /*!
    @brief Increment flux registers
    @param[in] a_flux Face centered fluxes
  */
  virtual void incrementFluxRegister(const EBAMRFluxData& a_flux);

  /*!
    @brief Reset the flux register
  */
  virtual void resetFluxRegister();

  /*!
    @brief Increment the coarse-fine redist with mass difference (only for EB x CF)
    @param[in] a_massDifference Mass difference
  */
  virtual void coarseFineIncrement(const EBAMRIVData& m_massDifference);

  /*!
    @brief Do mass redistribution on each level. 
    @param[inout] a_stableDivergence Hybrid divergence on input, stable and conservative divergence on output. 
    @param[in]    a_massDifference Mass difference
    @details This does not account for distribution across refinement boundaries (need the other routines for that). 
  */
  virtual void hyperbolicRedistribution(EBAMRCellData& a_stableDivergence, const EBAMRIVData& a_massDifference);

  /*!
    @brief Increment redistribution flux registers
  */
  virtual void incrementRedistFlux();

  /*!
    @brief Increment redistribution objects with mass to be redistributed.
    @param[in] a_massDifference Mass to be redistributed. 
  */
  virtual void incrementRedist(const EBAMRIVData& a_massDifference);

  /*!
    @brief Do coarse-fine redistribution
    @param[inout] a_phi On output, this will be incremented with redistributed mass across refinement boundaries. 
  */
  virtual void coarseFineRedistribution(EBAMRCellData& a_phi);

  /*!
    @brief Reflux across refinement boundaries, using whatever is inside the flux registers. 
    @param[inout] a_phi Cell-centered state
  */
  virtual void reflux(EBAMRCellData& a_phi);

  /*!
    @brief Set ebis
    @param[in] a_ebis EBIndexspace
  */
  virtual void setEbIndexSpace(const RefCountedPtr<EBIndexSpace>& a_ebis);

  /*!
    @brief Define stencils for doing face-centered to face-centroid-centered states. 
  */
  virtual void defineInterpolationStencils();

  /*!
    @brief Define interpolant for holding the face-interpolated results. 
  */
  virtual void defineInterpolant();

  /*!
    @brief Parses domain BC options
  */
  virtual void parseDomainBc();    

  /*!
    @brief Parses plot variables
  */
  virtual void parsePlotVariables();    

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

#include <CD_NamespaceFooter.H>

#endif
