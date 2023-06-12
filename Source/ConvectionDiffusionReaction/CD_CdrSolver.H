/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
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
#include <functional>

// Our includes
#include <CD_AmrMesh.H>
#include <CD_ComputationalGeometry.H>
#include <CD_CdrSpecies.H>
#include <CD_CdrDomainBC.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Base class for solving convection-diffusion-reaction equations. 
  @details On embedded boundaries, users must always set the flux directly (but possibly via a function, if desired). 
  @details On domain edges/faces, the user can select between various ways of setting the domain boundary conditions (e.g. wall, data-based, function-based, outflow). 
*/
class CdrSolver
{
public:
  /*!
    @brief Default constructor. 
    @note This sets the realm to primal and populates the domain flux functions (with zero-fluxes, i.e. wall BCs)
  */
  CdrSolver();

  /*!
    @brief Disallowed copy constructor
    @param[in] a_other Other solver
  */
  CdrSolver(const CdrSolver& a_other) = delete;

  /*!
    @brief Disallowed move constructor
    @param[in] a_other Other solver
  */
  CdrSolver(const CdrSolver&& a_other) = delete;

  /*!
    @brief Disallowed assignment operator
    @param[in] a_other Other solver
  */
  CdrSolver&
  operator=(const CdrSolver& a_other) = delete;

  /*!
    @brief Disallowed move assignement operator
    @param[in] a_other Other solver
  */
  CdrSolver&
  operator=(const CdrSolver&& a_other) = delete;

  /*!
    @brief Constructor
  */
  virtual ~CdrSolver();

  /*!
    @brief Parse class options
  */
  virtual void
  parseOptions() = 0;

  /*!
    @brief Parse runtime options
  */
  virtual void
  parseRuntimeOptions() = 0;

  /*!
    @brief This sets default boundary conditions (wall type). 
    @details Called in constructor for setting domain BCs.
    @note This sets the domain boundary condition to be a wall BC (no incoming/outgoing mass). This is the default
    behavior -- if the user wants different BCs he can call setDomainBcType and setDomainBcFunction. 
  */
  void
  setDefaultDomainBC();

  /*!
    @brief Set domain bc type on domain side. 
    @param[in] a_domainSide Domain side
    @param[in] a_bcType     Boundary condition type
  */
  void
  setDomainBcType(const CdrDomainBC::DomainSide a_domainSide, const CdrDomainBC::BcType a_bcType);

  /*!
    @brief Set domain bc function on particular domain side
    @param[in] a_domainSide   Domain side
    @param[in] a_fluxFunction Flux function
    @note This only associates a function with a domain side -- to use it one must also call setDomainBcType(a_domainSide, BcType::Function)
  */
  void
  setDomainBcFunction(const CdrDomainBC::DomainSide a_domainSide, const CdrDomainBC::FluxFunction a_fluxFunction);

  /*!
    @brief Implicit diffusion Euler advance without source term.
    @param[inout] a_newPhi Solution at time t + dt
    @param[in]    a_oldPhi Solution at time t
    @param[in]    a_dt     Time step
    @note This calls the other version with a zero source term. 
  */
  virtual void
  advanceEuler(EBAMRCellData& a_newPhi, const EBAMRCellData& a_oldPhi, const Real a_dt);

  /*!
    @brief Implicit diffusion Euler advance with source term.
    @param[inout] a_newPhi Solution at time t + dt
    @param[in]    a_oldPhi Solution at time t
    @param[in]    a_source Source term. 
    @param[in]    a_dt     Time step
    @note For purely implicit Euler the source term should be centered at t+dt (otherwise it's an implicit-explicit method)
    @details This solves the implicit diffusion equation equation a_newPhi - a_oldPhi = dt*Laplacian(a_newPhi) + dt*a_source. 
  */
  virtual void
  advanceEuler(EBAMRCellData&       a_newPhi,
               const EBAMRCellData& a_oldPhi,
               const EBAMRCellData& a_source,
               const Real           a_dt) = 0;

  /*!
    @brief Implicit diffusion Crank-Nicholson advance without source term.
    @param[inout] a_newPhi Solution at time t + dt
    @param[in]    a_oldPhi Solution at time t
    @param[in]    a_dt     Time step
    @note This calls the other version with a zero source term. 
  */
  virtual void
  advanceCrankNicholson(EBAMRCellData& a_newPhi, const EBAMRCellData& a_oldPhi, const Real a_dt);

  /*!
    @brief Implicit diffusion Crank-Nicholson advance with source term.
    @param[inout] a_newPhi Solution at time t + dt
    @param[in]    a_oldPhi Solution at time t
    @param[in]    a_source Source term. 
    @param[in]    a_dt     Time step
  */
  virtual void
  advanceCrankNicholson(EBAMRCellData&       a_newPhi,
                        const EBAMRCellData& a_oldPhi,
                        const EBAMRCellData& a_source,
                        const Real           a_dt) = 0;

  /*!
    @brief Compute div(J) explicitly, where J = nV - D*grad(n)
    @param[out] a_divJ             Divergence term, i.e. finite volume approximation to 
    @param[in]  a_phi              Cell-centered state
    @param[in]  a_extrapDt         Extrapolation in time, i.e. shifting of div(J) towards e.g. half time step. Only affects the advective term. 
    @param[in]  a_conservativeOnly If true, we compute div(J) = 1/dx*sum(fluxes), which does not involve redistribution.
    @param[in]  a_ebFlux           If true, the embedded boundary flux will injected and included in div(J)
    @param[in]  a_domainFlux       If true, the domain flux will injected and included in div(J)
    @note a_phi is non-const because ghost cells will be re-filled
  */
  virtual void
  computeDivJ(EBAMRCellData& a_divJ,
              EBAMRCellData& a_phi,
              const Real     a_extrapDt,
              const bool     a_conservativeOnly,
              const bool     a_ebFlux,
              const bool     a_domainFlux) = 0;

  /*!
    @brief Compute div(v*phi) explicitly
    @param[out] a_divF             Divergence term, i.e. finite volume approximation to Div(v*phi), including redistribution magic. 
    @param[in]  a_phi              Cell-centered state
    @param[in]  a_extrapDt         Extrapolation in time, i.e. shifting of div(F) towards e.g. half time step. Only affects the advective term. 
    @param[in]  a_conservativeOnly If true, we compute div(F)=1/dx*sum(fluxes), which does not involve redistribution.
    @param[in]  a_domainBc         How to set domain fluxes
    @param[in]  a_ebFlux           If true, the embedded boundary flux will injected be included in div(F)
    @param[in]  a_domainFlux       If true, the domain flux will injected and included in div(F)
    @note a_phi is non-const because ghost cells will be interpolated in this routine. Valid data in a_phi is not touched. 
  */
  virtual void
  computeDivF(EBAMRCellData& a_divF,
              EBAMRCellData& a_phi,
              const Real     a_extrapDt,
              const bool     a_conservativeOnly,
              const bool     a_ebFlux,
              const bool     a_domainFlux) = 0;

  /*!
    @brief Compute div(D*grad(phi)) explicitly
    @param[out] a_divF             Divergence term, i.e. finite volume approximation to Div(D*Grad(phi)). 
    @param[in]  a_phi              Cell-centered state
    @param[in]  a_domainBc         Flag for setting domain boundary conditions
    @param[in]  a_conservativeOnly If true, we compute div(D) = 1/dx*sum(fluxes), which does not involve redistribution.
    @param[in]  a_useEbFlux        If true, the embedded boundary flux will be injected and included in div(D)
    @param[in]  a_domainFlux       If true, the domain flux will injected and included in div(D)
    @note a_phi is non-const because ghost cells will be interpolated in this routine. Valid data in a_phi is not touched. 
  */
  virtual void
  computeDivD(EBAMRCellData& a_divD,
              EBAMRCellData& a_phi,
              const bool     a_conservativeOnly,
              const bool     a_ebFlux,
              const bool     a_domainFlux) = 0;

  /*!
    @brief Compute div(G) where G is a general face-centered flux on face centers and EB centers. This can involve mass redistribution.
    @param[in]    a_divG             div(G) or kappa*div(G).
    @param[inout] a_G                Vector field which contains face-centered fluxes on input. Contains face-centroid fluxes on output. 
    @param[in]    a_ebFlux           Flux on the EB centroids
    @param[in]    a_conservativeOnly If true, we compute div(G)=1/dx*sum(fluxes), which does not involve redistribution.
  */
  virtual void
  computeDivG(EBAMRCellData& a_divG, EBAMRFluxData& a_G, const EBAMRIVData& a_ebFlux, const bool a_conservativeOnly);

  /*!
    @brief Compute a random gaussian white noise source term. 
    @param[out] a_noiseSource Source term
    @param[in]  a_cellPhi     Cell-centered states
    @note Used e.g. for fluctuating hydrodynamics. 
  */
  virtual void
  gwnDiffusionSource(EBAMRCellData& a_noiseSource, const EBAMRCellData& a_cellPhi);

  /*!
    @brief Average velocities to faces
    @note This takes the cell-centered velocities in m_cellVelocity and averages them onto m_faceVelocity. Only the velocity component normal
    to the face is stored in m_faceVelocity. 
  */
  virtual void
  averageVelocityToFaces();

  /*!
    @brief Perform pre-regrid operations. 
    @param[in] a_lbase          Coarsest level that changed during regrid. 
    @param[in] a_oldFinestLevel Finest grid level before the regrid operation. 
    @note This stores m_phi and m_source. 
  */
  virtual void
  preRegrid(const int a_lbase, const int a_oldFinestLevel);

  /*!
    @brief Deallocate internal storage
    @note This deallocates a bunch of storage which is not needed (during regrids). It can be used to trim memory during Berger-Rigoutsous regrids (which EAT memory). 
  */
  virtual void
  deallocate();

  /*!
    @brief Register operators for AMR operations
    @note This includes operators for redistribution, flux registers, regridding, ghost cell interpolation, and conservative coarsening. 
  */
  virtual void
  registerOperators();

  /*!
    @brief Set the realm for this solver.
    @param[in] a_realm Realm identifier
  */
  virtual void
  setRealm(const std::string a_realm);

  /*!
    @brief Get the CDR species
  */
  virtual RefCountedPtr<CdrSpecies>&
  getSpecies() noexcept;

  /*!
    @brief Get the CDR species
  */
  virtual const RefCountedPtr<CdrSpecies>&
  getSpecies() const noexcept;

  /*!
    @brief Set species
    @param[in] a_species Species (where this solver gets its name, charge, initial conditions etc). 
  */
  virtual void
  setSpecies(const RefCountedPtr<CdrSpecies>& a_species);

  /*!
    @brief Set computational geometry
    @param[in] a_computationalGeometry The computational geometry.
  */
  virtual void
  setComputationalGeometry(const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry);

  /*!
    @brief Set the amr object
    @param[in] a_amr AmrMesh object. 
  */
  virtual void
  setAmr(const RefCountedPtr<AmrMesh>& a_amr);

  /*!
    @brief Set phase.
    @param[in] a_phase Phase.
  */
  virtual void
  setPhase(phase::which_phase a_phase);

  /*!
    @brief Set verbosity.
    @param[in] a_verbosity Verbosity level. 
  */
  virtual void
  setVerbosity(const int a_verbosity);

  /*!
    @brief Set the time for this solver. 
    @param[in] a_step Time step number
    @param[in] a_time Time (in seconds)
    @param[in] a_dt   Time step increment
    @note This sets m_step=a_step, m_time=a_time, m_dt=a_dt
  */
  virtual void
  setTime(const int a_step, const Real a_time, const Real a_dt);

  /*!
    @brief Set velocity from data holder. 
    @param[in] a_velocity Cell centered velocities. 
    @details This copies a_velo onto m_cellVelocity. 
  */
  virtual void
  setVelocity(const EBAMRCellData& a_velocity);

  /*!
    @brief Set constant velocity
    @param[in] a_velocity Velocity
    @details This sets the cell centered velocity to be a_velo in every grid cell. 
  */
  virtual void
  setVelocity(const RealVect a_velocity);

  /*!
    @brief Set velocity using a polymorphic function. 
    @param[in] a_velocity Velocity to be set. You can use e.g. a lambda to set it. 
  */
  virtual void
  setVelocity(const std::function<RealVect(const RealVect a_pos)>& a_velocity);

  /*!
    @brief Data-based version of setting diffusion coefficients (which are stored on faces)
    @param[in] a_diffusionCoefficient   Face-centered diffusion coefficient. 
    @param[in] a_ebDiffusionCoefficient EB-centered diffusion coefficient. 
    @note Realms do not have to be the same. 
  */
  virtual void
  setDiffusionCoefficient(const EBAMRFluxData& a_diffusionCoefficient, const EBAMRIVData& a_ebDiffusionCoefficient);

  /*!
    @brief Set diffusion coefficient to be constant everywhere (they are stored on faces)
    @param[in] a_diffusionCoefficient Global diffusion coefficient. 
    @details This sets the diffusion coefficient to be constant on both faces and EB faces. 
  */
  virtual void
  setDiffusionCoefficient(const Real a_diffusionCoefficient);

  /*!
    @brief Polymorphic way of setting diffusion coefficients every. 
    @param[in] a_diffusionCoefficient Diffusion coefficient. 
    @note This will evaluate the function in all grid cells and set the diffusion coefficient 
  */
  virtual void
  setDiffusionCoefficient(const std::function<Real(const RealVect a_position)>& a_diffusionCoefficient);

  /*!
    @brief Data based version of setting source terms. 
    @param[in] a_source Source term.
    @details This sets m_source to be equal to a_source
  */
  virtual void
  setSource(const EBAMRCellData& a_source);

  /*!
    @brief Set constant source terms everywhere
    @param[in] a_source Source term.
    @note This sets m_source to a_source everywhere
  */
  virtual void
  setSource(const Real a_source);

  /*!
    @brief Polymorphic way of setting source terms
    @param[in] a_source Source term.
    @note This fills m_source with the evaluation of a_source
  */
  virtual void
  setSource(const std::function<Real(const RealVect a_position)> a_source);

  /*!
    @brief Data-based version of setting the EB flux
    @param[in] a_ebFlux Flux on EB centroids. 
    @note a_ebFlux does not have to be the same realm as in the solver. 
  */
  virtual void
  setEbFlux(const EBAMRIVData& a_ebFlux);

  /*!
    @brief Set the eb flux to a constant. 
    @param[in] a_ebFlux Flux on EB centroids. 
    @details This sets a constant EB-flux on all EBs. 
  */
  virtual void
  setEbFlux(const Real a_ebFlux);

  /*!
    @brief Fill m_phi state with initial data from m_species
    @note This will call initialDataParticles() and initialDataFunction().
  */
  virtual void
  initialData();

  /*!
    @brief Allocate internal storage
    @note This allocates a bunch of storage -- if the solver is not diffusive or mobile then we do our best to trim memory. 
  */
  virtual void
  allocate();

  /*!
    @brief Write plot file
    @details The name of the plot file is m_name.stepXXXXX.DIM.hdf5
    @note This calls writePlotData(...)
  */
  virtual void
  writePlotFile();

  /*!
    @brief Write output data to a_output
    @param[inout] a_output Output data holder.
    @param[inout] a_icomp  Starting component where this solver begins writing the output. 
    @param[in]    a_outputRealm Realm where a_output belongs
    @param[in]    a_level  Grid level
    @note This will write the plot data in this solver to a_output, starting on a_comp
    @note This routine writes m_phi on centroids (and not cell-centers). 
  */
  virtual void
  writePlotData(LevelData<EBCellFAB>& a_output, int& a_icomp, const std::string a_outputRealm, const int a_level) const noexcept;

  /*!
    @brief Write checkpoint data into HDF5 file. 
    @paramo[out] a_handle HDF5 file. 
    @param[in]   a_level Grid level
  */
#ifdef CH_USE_HDF5
  virtual void
  writeCheckpointLevel(HDF5Handle& a_handle, const int a_level) const;
#endif

  /*!
    @brief Read checkpoint data from HDF5 file. 
    @param[in] a_handle HDF5 handle.
    @param[in] const int a_level Grid level
  */
#ifdef CH_USE_HDF5
  virtual void
  readCheckpointLevel(HDF5Handle& a_handle, const int a_level);
#endif

  /*!
    @brief Regrid this solver. 
    @param[in] a_lmin           Coarsest level where grids did not change. 
    @param[in] a_oldFinestLevel Finest AMR level before the regrid. 
    @param[in] a_newFinestLevel Finest AMR level after the regrid. 
    @details This linearly interpolates (with limiters) m_state to the new grids. Velocities are left undefined. 
  */
  virtual void
  regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel);

  /*!
    @brief Get the realm where this solver is registered. 
    @return Returns realm name. 
  */
  virtual std::string
  getRealm() const;

  /*!
    @brief Get solver name
    @note Not necessarily equal to class name (we can have many CdrSolvers instantiated). 
  */
  virtual std::string
  getName() const;

  /*!
    @brief Get output plot names
    @return Return a list of plot variable names. 
  */
  virtual Vector<std::string>
  getPlotVariableNames() const;

  /*!
    @brief Get number of output fields
    @return Returns number of plot variables include in writePlotFile() and writePlotData()
  */
  virtual int
  getNumberOfPlotVariables() const;

  /*!
    @brief Compute the largest possible diffusive time step (for explicit methods)
    @details This computes dt = dx/sum(|vx| + |vy| + |vz|), minimized over all grid levels and patches.  
    @note This is the appropriate time step routine for explicit advection solvers. 
  */
  virtual Real
  computeAdvectionDt();

  /*!
    @brief Compute the largest possible diffusive time step (for explicit methods)
    @details This computes dt = (dx*dx)/(2*D*d) where D is the diffusion coefficient. The result is minimized over all grid levels and patches. 
    @note This is the appropriate time step routine for explicit diffusion solvers. 
  */
  virtual Real
  computeDiffusionDt();

  /*!
    @brief Compute the largest possible diffusive time step (for explicit methods)
    @details This computes dt = 1/[ dx/(|vx|+|vy|+|vz|) + (dx*dx)/(2*D*d) ] and minimizes the result over all grid levels and patches. 
    @note This is the appropriate time step routine for explicit advection-diffusion solvers. 
  */
  virtual Real
  computeAdvectionDiffusionDt();

  /*!
    @brief Compute the largest possible source time step (for explicit methods
    @param[in] a_max       Maximum value of m_phi
    @param[in] a_tolerance Tolerance
    @details This computes dt = phi/source, but only for cells where phi lies within a_tolerance*a_max (and source > 0)
    @todo An old routine which we could (probably) remove. 
  */
  virtual Real
  computeSourceDt(const Real a_max, const Real a_tolerance);

  /*!
    @brief Compute an upwind-weighted version of phi
    @details This is equal to sum(alpha * v^p * phi)/sum(alpha * v^p) where alpha is the face aperture and v is the face velocity. phi is
    taken as the upwind value. 
    @param[out] a_weightedUpwindPhi Weighted upwind phi
    @param[in]  a_pow               Weighting factor
  */
  virtual void
  weightedUpwind(EBAMRCellData& a_weightedUpwindPhi, const int a_pow);

  /*!
    @brief Compute the "physical mass" in m_phi
    @details This conservatively coarsens the solution and computes on the coarsest level only
  */
  virtual Real
  computeMass();

  /*!
    @brief Compute the "physical mass" in the input argument. 
    @details This conservatively coarsens the solution and computes on the coarsest level only
    @param[in] a_phi        Cell-centered data. 
    @param[in] a_kappaScale Multiply by kappa in irregular cells.
  */
  virtual Real
  computeMass(const EBAMRCellData& a_phi, const bool a_kappaScale = true);

  /*!
    @brief Compute the total charge in m_phi
    @details Calls computeMass and multiplies by charge
  */
  virtual Real
  computeCharge();

  /*!
    @brief Return true if the solver is diffusive and false otherwise
  */
  virtual bool
  isDiffusive();

  /*!
    @brief Return true if the solver is mobile and false otherwise
  */
  virtual bool
  isMobile();

  /*!
    @brief Get the cell-centered phi
    @return m_phi
  */
  virtual EBAMRCellData&
  getPhi();

  /*!
    @brief Get the source term
    @return m_source
  */
  virtual EBAMRCellData&
  getSource();

  /*!
    @brief Get the cell-centered velocity
    @return m_cellVelocity
  */
  virtual EBAMRCellData&
  getCellCenteredVelocity();

  /*!
    @brief Get the face-centered velocities
    @return m_faceVelocity
  */
  virtual EBAMRFluxData&
  getFaceCenteredVelocity();

  /*!
    @brief Get the eb-centered velocities
    @return m_ebVelocity
  */
  virtual EBAMRIVData&
  getEbCenteredVelocity();

  /*!
    @brief Get the cell-centered diffusion coefficient
  */
  virtual EBAMRCellData&
  getCellCenteredDiffusionCoefficient();

  /*!
    @brief Get the face-centered diffusion coefficient
    @return m_faceCentereDiffusionCoefficient
  */
  virtual EBAMRFluxData&
  getFaceCenteredDiffusionCoefficient();

  /*!
    @brief Get the EB-centered diffusion coefficient
    @return m_ebCentereDiffusionCoefficient
  */
  virtual EBAMRIVData&
  getEbCenteredDiffusionCoefficient();

  /*!
    @brief Get the eb flux data holder
    @return m_ebFlux
  */
  virtual EBAMRIVData&
  getEbFlux();

  /*!
    @brief Get the domain flux data holder
    @return m_domainFlux
  */
  virtual EBAMRIFData&
  getDomainFlux();

  /*!
    @brief Extrapolate advective flux to EB. 
    @details This does not use any slope limiting (derived classes might want to redo this routine). The returned flux
    is equal to -n.(phi * v) where n is the normal vector pointing into the solution domain.
    @note Calls the other version with m_ebFlux as argument
  */
  virtual void
  extrapolateAdvectiveFluxToEB() noexcept;

  /*!
    @brief Extrapolate advective flux to EB. 
    @details This does not use any slope limiting (derived classes might want to redo this routine). The returned flux
    is equal to -n.(phi * v) where n is the normal vector pointing into the solution domain.
    @param[out] a_ebFlux Projected flux on the EB. 
  */
  virtual void
  extrapolateAdvectiveFluxToEB(EBAMRIVData& a_ebFlux) const noexcept;

  /*!
    @brief Add data through redistribution into cell-centered holders
    @param[inout] a_phi   Target data
    @param[in]    a_delta Data to be redistributed
  */
  virtual void
  redistribute(EBAMRCellData& a_phi, const EBAMRIVData& a_delta) const noexcept;

protected:
  /*!
    @brief Redistribution method
  */
  enum class Redistribution
  {
    VolumeWeighted,
    None
  };

  /*!
    @brief Component number in data holder
  */
  static constexpr int m_comp = 0;

  /*!
    @brief Number of components that this solver solves for
  */
  static constexpr int m_nComp = 1;

  /*!
    @brief Species through which e.g. mobility/diffusion and initial conditions is passed. 
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
    @brief Stencils for interpolating face-centered fluxes to face centroids
  */
  Vector<RefCountedPtr<LayoutData<BaseIFFAB<FaceStencil>>>> m_interpStencils[SpaceDim];

  /*!
    @brief Phase
  */
  phase::which_phase m_phase;

  /*!
    @brief Solver name
  */
  std::string m_name;

  /*!
    @brief Class name. 
    @details This will be different for different implementations of CdrSolver
  */
  std::string m_className;

  /*!
    @brief Realm where this solver is registered
  */
  std::string m_realm;

  /*!
    @brief Domain BCs
  */
  CdrDomainBC m_domainBC;

  /*!
    @brief Cell-centered data (i.e. the advected-diffused quantity)
  */
  EBAMRCellData m_phi;

  /*!
    @brief Source term 
  */
  EBAMRCellData m_source;

  /*!
    @brief Cell-centered velocities
    @note Memory is only allocated if the solver is mobile. 
  */
  EBAMRCellData m_cellVelocity;

  /*!
    @brief Holder for face centered states
    @note Memory is only allocated if the solver is mobile. 
  */
  EBAMRFluxData m_faceStates;

  /*!
    @brief Face-centered velocities (only normal components)
    @note Memory is only allocated if the solver is mobile. 
  */
  EBAMRFluxData m_faceVelocity;

  /*!
    @brief Scratch storage for the non-conservative divergence
  */
  EBAMRIVData m_nonConservativeDivG;

  /*!
    @brief Scratch storage for the mass difference
  */
  EBAMRIVData m_massDifference;

  /*!
    @brief Scratch storage for the EB flux
  */
  EBAMRIVData m_ebZero;

  /*!
    @brief Cached state vector for regrid. 
  */
  EBAMRCellData m_cachePhi;

  /*!
    @brief Cached source term for regrids. 
  */
  EBAMRCellData m_cacheSource;

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
    @brief Diffusion coefficients on cell centers
    @note Memory is only allocated if the solver is diffusive
  */
  EBAMRCellData m_cellCenteredDiffusionCoefficient;

  /*!
    @brief Diffusion coefficients on face centers
    @note Memory is only allocated if the solver is diffusive
  */
  EBAMRFluxData m_faceCenteredDiffusionCoefficient;

  /*!
    @brief Diffusion coefficients on EB faces
    @note Memory is only allocated if the solver is diffusive
  */
  EBAMRIVData m_ebCenteredDiffusionCoefficient;

  /*!
    @brief Domain flux functions 
    @details Used when the user specifies that the domain flux is function-based
  */
  std::map<CdrDomainBC::DomainSide, CdrDomainBC::FluxFunction> m_domainFluxFunctions;

  /*!
    @brief Solver verbosity. 
    @details A high number gives lots of output 
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
    @brief Last time step increment. 
  */
  Real m_dt;

  /*!
    @brief Which type of redistribution to use
  */
  Redistribution m_whichRedistribution;

  /*!
    @brief Flag for blending the hybrid divergence
  */
  bool m_blendConservation;

  /*!
    @brief Is the solver diffusive or not
  */
  bool m_isDiffusive;

  /*!
    @brief Solve for advection/convection or not
  */
  bool m_isMobile;

  /*!
    @brief If true, m_phi is added to plot files. 
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
    @brief Use slopes when regridding
  */
  bool m_regridSlopes;

  /*!
    @brief RNG seed
  */
  int m_seed;

  /*!
    @brief Average cell-centered velocities to face centers
    @param[out] a_faceVelocity Face-centered velocities
    @param[in]  a_cellVelocity Cell-centered velocities
    @note This sets the velocity component normal to the face to be the arithmetic average of the cell center velocities on the low/high side of the face. 
    @note Virtual because some solvers might want to do this differently. 
  */
  virtual void
  averageVelocityToFaces(EBAMRFluxData& a_faceVelocity, const EBAMRCellData& a_cellVelocity);

  /*!
    @brief Advection-only extrapolation to faces.
    @param[out] a_facePhi  Phi on faces
    @param[in]  a_phi      Phi on cell center
    @param[in]  a_extrapDt Time centering/extrapolation (if the advective integrator can do it)
  */
  virtual void
  advectToFaces(EBAMRFluxData& a_facePhi, const EBAMRCellData& a_phi, const Real a_extrapDt) = 0;

  /*!
    @brief Set up face-centered advection flux.
    @param[out] a_flux          Face-centered fluxes
    @param[in]  a_facePhi       Face-centered state
    @param[in]  a_faceVelocity  Face-centered velocities
    @param[in]  a_addDomainFlux Add the domain flux 
    @note This computes flux = phi*vel on all face centers. If a_addDomainFlux is true, we enforce BC fluxes on the domain sides. Otherwise, the flux on the domain
    sides is set to zero. 
  */
  virtual void
  computeAdvectionFlux(EBAMRFluxData&       a_flux,
                       const EBAMRFluxData& a_facePhi,
                       const EBAMRFluxData& a_faceVelocity,
                       const bool           a_addDomainFlux = true);

  /*!
    @brief Set up face-centered advection flux on a grid level
    @param[out] a_flux          Face-centered fluxes
    @param[in]  a_facePhi       Face-centered state
    @param[in]  a_faceVelocity  Face-centered velocities
    @param[in]  a_domainFlux    Domain flux
    @param[in]  a_level         Grid level
    @note This computes flux = v*phi on all face centers. If a_addDomainFlux is true, we enforce BC fluxes on the domain sides. Otherwise, the flux on the domain
    sides is set to zero. 
  */
  virtual void
  computeAdvectionFlux(LevelData<EBFluxFAB>&       a_flux,
                       const LevelData<EBFluxFAB>& a_facePhi,
                       const LevelData<EBFluxFAB>& a_faceVelocity,
                       const int                   a_lvl);

  /*!
    @brief Compute the face-centered diffusion flux
    @param[out] a_flux          Face-centered flux
    @param[in]  a_phi           Cell centered phi
    @param[in]  a_addDomainFlux Add domain flux or not
    @note This computes flux = D*grad(phi) on all face centers. If a_addDomainFlux is true, we enforce BC fluxes on the domain sides. Otherwise, the flux on the domain
    sides is set to zero. 
  */
  virtual void
  computeDiffusionFlux(EBAMRFluxData& a_flux, const EBAMRCellData& a_phi, const bool a_addDomainFlux);

  /*!
    @brief Compute the face-centered diffusion flux
    @param[out] a_flux Face-centered flux
    @param[in]  a_phi  Cell centered phi
    @param[in]  a_lvl  Grid level
    @note This computes flux = D*grad(phi) on all face centers. The flux is set to zero on domain edges/faces
  */
  virtual void
  computeDiffusionFlux(LevelData<EBFluxFAB>& a_flux, const LevelData<EBCellFAB>& a_phi, const int a_lvl);

  /*!
    @brief Compute the full advection-diffusion flux. This assumes that the solver is mobile and diffusive.
    @param[out] a_flux           Face-centered flux. 
    @param[in]  a_cellStates     Cell-centered states. 
    @param[in]  a_faceVelocities Face-centered velocities.
    @param[in]  a_faceStates     Face-centered (advected) states. 
    @param[in]  a_faceDiffCo     Face-centered diffusion coefficients. 
    @param[in]  a_addDomainFlux  Add domain flux or not. 
  */
  virtual void
  computeAdvectionDiffusionFlux(EBAMRFluxData&       a_flux,
                                const EBAMRCellData& a_cellStates,
                                const EBAMRFluxData& a_faceStates,
                                const EBAMRFluxData& a_faceVelocities,
                                const EBAMRFluxData& a_faceDiffCo,
                                const bool           a_addDomainFlux);

  /*!
    @brief Set flux to zero on domain boundaries
    @param[in] a_flux Flux data holder -- on domain edges this is modified so the flux is zero. 
  */
  virtual void
  resetDomainFlux(EBAMRFluxData& a_flux);

  /*!
    @brief Set domain in data holder. This sets the flux on the boundary to either zero or to m_domainFlux
    @param[inout] a_flux Flux to be modified. 
  */
  virtual void
  fillDomainFlux(EBAMRFluxData& a_flux);

  /*!
    @brief Set domain in data holder. This sets the flux on the boundary to either zero or to m_domainFlux
    @param[inout] a_flux  Flux to be modified. 
    @param[in]    a_level Grid level
  */
  virtual void
  fillDomainFlux(LevelData<EBFluxFAB>& a_flux, const int a_level);

  /*!
    @brief Compute conservative divergence from fluxes. 
    @param[out] a_conservativeDivergence Conservative divergence computed as div(F) using finite volumes. Not divided by kappa. 
    @param[in]  a_flux                   Face-centered fluxes. Includes domain fluxes.
    @param[in]  a_ebFlux                 EB flux. 
    @note This computes the conservative divergence kappa*div(F) = sum(fluxes). 
  */
  virtual void
  conservativeDivergenceNoKappaDivision(EBAMRCellData&     a_conservativeDivergence,
                                        EBAMRFluxData&     a_flux,
                                        const EBAMRIVData& a_ebFlux);

  /*!
    @brief Compute the non-conservative divergence
    @param[out] a_nonConservativeDivergence Non-conservative divergence.
    @param[in]  a_divG                      Conservative divergence.
    @note This will compute the non-conservative divergence as divNC(G) = sum(kappa*div(G))/sum(kappa) in the cut-cell neighborhood. 
  */
  virtual void
  nonConservativeDivergence(EBAMRIVData& a_nonConservativeDivergence, const EBAMRCellData& a_divG);

  /*!
    @brief Use the non-conservative divergence to make the conservative divergence hold the hybrid divergence. 
    @param[inout] a_hybridDivergence          On input, contains the conservative divergence (without kappa division). Contains hybrid divergence on output.
    @param[out]   a_massDifference            Mass difference between updating with hybrid divergence and true divergence
    @param[in]    a_nonConservativeDivergence Non-conservative divergence.
  */
  virtual void
  hybridDivergence(EBAMRCellData&     a_hybridDivergence,
                   EBAMRIVData&       a_massDifference,
                   const EBAMRIVData& a_nonConservativeDivergence);

  /*!
    @brief Make the hybrid divergence. On the way in, a_hybridDivergence must hold the conservative divergence
    @param[inout] a_hybridDivergence          On input, contains the conservative divergence (without kappa division). Contains hybrid divergence on output.
    @param[out]   a_massDifference            Mass difference between updating with hybrid divergence and true divergence
    @param[in]    a_nonConservativeDivergence Non-conservative divergence.
    @param[in]    a_lvl                       Grid level. 
  */
  virtual void
  hybridDivergence(LevelData<EBCellFAB>&             a_hybridDivergence,
                   LevelData<BaseIVFAB<Real>>&       a_massDifference,
                   const LevelData<BaseIVFAB<Real>>& a_nonConservativeDivergence,
                   const int                         a_lvl);

  /*!
    @brief Compute the conservative divergence over regular cells
    @param[out] a_divJ Conservative divergence. Not kappa-divided. 
    @param[in]  a_flux Face-centered fluxes
    @param[in]  a_lvl  Grid level
    @note This computes the finite volume approximation to div(J) in regular cells -- cut-cells (in a_divJ are set to zero)
  */
  virtual void
  conservativeDivergenceRegular(LevelData<EBCellFAB>& a_divJ, const LevelData<EBFluxFAB>& a_flux, const int a_lvl);

  /*!
    @brief Interpolate flux to centroids
    @param[inout] a_flux On input, contains centered fluxes. Output contains centroid fluxes
    @details This interpolates face fluxes to centroids. 
  */
  virtual void
  interpolateFluxToFaceCentroids(EBAMRFluxData& a_flux);

  /*!
    @brief Interpolate flux to centroids
    @param[inout] a_flux On input, contains centered fluxes. Output contains centroid fluxes
    @param[in]    a_lvl  Grid level
    @details This interpolates face fluxes to centroids. 
  */
  virtual void
  interpolateFluxToFaceCentroids(LevelData<EBFluxFAB>& a_flux, const int a_lvl);

  /*!
    @brief Compute conservative divergence on irregular cells (not kappa divided)
    @param[out] a_divG         Conservative divergence (not kappa divided)
    @param[out] a_centroidFlux Face centroid centered fluxes
    @param[in]  a_ebFlux       EB flux
    @param[in]  a_lvl          Grid level
    @note This overwrites the result of divG in cut-cells, using the eb flux and face centroid fluxes to compute kappa*div(G) = sum(fluxes)/dx
  */
  virtual void
  computeDivergenceIrregular(LevelData<EBCellFAB>&             a_divG,
                             const LevelData<EBFluxFAB>&       a_centroidFluxes,
                             const LevelData<BaseIVFAB<Real>>& a_ebFlux,
                             const int                         a_lvl);

  /*!
    @brief Fill initial data from a distribution function
    @note This increments m_phi on the mesh with data from m_species->initialData().
  */
  virtual void
  initialDataDistribution();

  /*!
    @brief Fill initial data from particles
    @note This will set m_phi to zero before depositing the particles (which are deposited with an NGP scheme, ignoring modifactions near cut-cells)
  */
  virtual void
  initialDataParticles();

  /*!
    @brief Define stencils for doing face-centered to face-centroid-centered states. 
    @note This computes standard finite-difference stencils for interpolating from face centers to face centroids. 
  */
  virtual void
  defineInterpolationStencils();

  /*!
    @brief Parses domain BC options
  */
  virtual void
  parseDomainBc();

  /*!
    @brief Parses plot variables
  */
  virtual void
  parsePlotVariables();

  /*!
    @brief Parse plot mode
  */
  virtual void
  parsePlotMode();

  /*!
    @brief Parse the conservation
  */
  virtual void
  parseDivergenceComputation();

  /*!
    @brief Parse slope regrid
  */
  virtual void
  parseRegridSlopes();

  /*!
    @brief Shortcut for making a boundary condition string. 
    @param[in] a_dir  Direction.
    @param[in] a_side Coordinate side. 
    @return Returns string of type m_className.bc.direction.side. 
  */
  virtual std::string
  makeBcString(const int a_dir, const Side::LoHiSide a_side) const;

  /*!
    @brief Gaussian noise field
    @param[out] a_noise Gaussian white nosie
    @param[out] a_sigma Standard deviation
  */
  virtual void
  fillGwn(EBAMRFluxData& a_noise, const Real a_sigma);

  /*!
    @brief Use Heaviside smoothing for computing face-centered states. 
    @param[out] a_facePhi Face centered state
    @param[in]  a_cellPhi Cell-centered states
    @note This is only intended to be used with FHD functionality. 
  */
  virtual void
  smoothHeavisideFaces(EBAMRFluxData& a_facePhi, const EBAMRCellData& a_cellPhi);

  /*!
    @brief Write data to output. Convenience function. 
    @param[inout] a_output Output data holder.
    @param[inout] a_icomp  Starting component where this solver begins writing the output. 
    @param[in]    a_data   Data to write. 
    @param[in]    a_level  Grid level
    @param[in]    a_interpToCentroids If true, a_data will be interpolated to cell centroids before writing to a_output. 
    @param[in]    a_interpGhost If true, interpolate ghost cells
  */
  virtual void
  writeData(LevelData<EBCellFAB>& a_output,
            int&                  a_comp,
            const EBAMRCellData&  a_data,
            const int             a_level,
            const bool            a_interpToCentroids,
            const bool            a_interpGhost) const noexcept;
};

#include <CD_NamespaceFooter.H>

#endif
