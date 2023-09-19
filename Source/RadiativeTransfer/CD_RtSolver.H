/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_RtSolver.H
  @brief  Abstract parent class for various radiative transfer solvers
  @author Robert Marskar
*/

#ifndef CD_RtSolver_H
#define CD_RtSolver_H

// Std includes
#include <functional>

// Our includes
#include <CD_ComputationalGeometry.H>
#include <CD_AmrMesh.H>
#include <CD_RtSpecies.H>
#include <CD_Location.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Abstract RTE solver class for doing various kinds of radiative transfer equations. This class
  is a pure class and it is basically just an interface.
*/
class RtSolver
{
public:
  /*!
    @brief Constructor
  */
  RtSolver();

  /*!
    @brief Constructor (does nothing)
  */
  virtual ~RtSolver();

  /*!
    @brief Get solver name
  */
  virtual std::string
  getName();

  /*!
    @brief Get the realm where the solver lives
  */
  virtual const std::string
  getRealm() const;

  /*!
    @brief Parse options
  */
  virtual void
  parseOptions() = 0;

  /*!
    @brief Parse runtime options
  */
  virtual void
  parseRuntimeOptions() = 0;

  /*!
    @brief Advance equation one time step
    @param[in] a_dt      Time step
    @param[in] a_zeroPhi An optional argument for e.g. elliptic solves which can start from previous solution or zero.
  */
  virtual bool
  advance(const Real a_dt, const bool a_zeroPhi = false);

  /*!
    @brief Advance method. Advances one time step.
    @param[in] a_dt      Time step
    @param[in] a_phi     Mesh solution
    @param[in] a_zeroPhi An optional argument for e.g. elliptic solves which can start from previous solution or zero.
  */
  virtual bool
  advance(const Real a_dt, EBAMRCellData& a_phi, const bool a_zeroPhi = false);

  /*!
    @brief Advance method. Advances one time step.
    @param[in] a_dt      Time step
    @param[in] a_phi     Mesh solution
    @param[in] a_source  Source term
    @param[in] a_zeroPhi An optional argument for e.g. elliptic solves which can start from previous solution or zero.
  */
  virtual bool
  advance(const Real a_dt, EBAMRCellData& a_phi, const EBAMRCellData& a_source, const bool a_zeroPhi = false) = 0;

  /*!
    @brief Set realm where this solver lives.
    @param[in] a_realm Realm
  */
  virtual void
  setRealm(const std::string a_realm);

  /*!
    @brief Set the radiative transfer species (RtSpecies)
    @param[in] a_species Species (see RtSpecies.H)
  */
  virtual void
  setRtSpecies(const RefCountedPtr<RtSpecies>& a_species);

  /*!
    @brief Set computational geometry
    @param[in] a_computationalGeometry Computational geometry
  */
  virtual void
  setComputationalGeometry(const RefCountedPtr<ComputationalGeometry> a_computationalGeometry);

  /*!
    @brief Get computational loads for a specific grid level
    @param[out] a_loads Grid loads for this level. 
    @param[in] a_dbl Grids on input level
    @param[in] a_level Input level
    @return Loads for each box on a grid level. 
    @note The return vector should have the same order as the boxes in a_dbl. E.g. ret[0] must be the load for a_dbl.boxArray()[0];
    @details The default implementation returns the number of cells in the grid patch as a proxy for the load. 
  */
  virtual void
  computeLoads(Vector<long long>& a_loads, const DisjointBoxLayout& a_dbl, const int a_level) const noexcept;

  /*!
    @brief Set the amr object
    @param[in] a_amr AmrMesh object.
  */
  virtual void
  setAmr(const RefCountedPtr<AmrMesh>& a_amr);

  /*!
    @brief Set phase
    @param[in] a_phase Phase (gas/solid) where the solver lives.
    @details This must be done BEFORE callilng setComputationalGeometry
  */
  virtual void
  setPhase(phase::which_phase a_phase = phase::gas);

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
    @brief Set stationary solver or not
    @param[in] a_stationary If true, the solver is set to stationary mode.
  */
  virtual void
  setStationary(const bool a_stationary);

  /*!
    @brief Sanity check
  */
  virtual void
  sanityCheck();

  /*!
    @brief Check if solver is stationary
  */
  virtual bool
  isStationary();

  /*!
    @brief Compute the boundary flux given a state (a_phi will be different for different RTE approximations)
    @param[inout] a_ebFlux The flux on the EB
    @param[in]    a_phi    Cell-centered solution for the RTE state (e.g. for Eddington this is just the isotropic part).
  */
  virtual void
  computeBoundaryFlux(EBAMRIVData& a_ebFlux, const EBAMRCellData& a_phi) = 0;

  /*!
    @brief Compute the domain flux given a state (a_phi will be different for different RTE approximations)
    @param[inout] a_domainFlux The flux on the domain wall
    @param[in]    a_phi        Cell-centered solution for the RTE state (e.g. for Eddington this is just the isotropic part).
  */
  virtual void
  computeDomainFlux(EBAMRIFData& a_domainflux, const EBAMRCellData& a_phi) = 0;

  /*!
    @brief Compute the flux
    @param[inout] a_flux The RTE flux on the domain wall
    @param[in]    a_phi  Cell-centered solution for the RTE state (e.g. for Eddington this is just the isotropic part).
    @details For diffusive models, the flux will be something like grad(a_phi). For higher-order models, the flux 
    will be contained in a_phi (somehow).
  */
  virtual void
  computeFlux(EBAMRCellData& a_flux, const EBAMRCellData& a_phi) = 0;

  /*!
    @brief Get isotropic part. 
    @param[inout] a_isotropic Isotropic part of the RTE solution
    @param[in]    a_phi       Full RTE solution.
    @details This is normally just a copy, but I need the function signature in order to have a clear plasma-coupling interface
  */
  virtual void
  computeDensity(EBAMRCellData& a_isotropic, const EBAMRCellData& a_phi) = 0;

  /*!
    @brief Allocate internal storage
  */
  virtual void
  allocate() = 0;

  /*!
    @brief Deallocate internal storage
  */
  virtual void
  deallocate() = 0;

  /*!
    @brief Perform pre-regrid operations. 
    @param[in] a_lbase          Coarsest level that changed during regrid. 
    @param[in] a_oldFinestLevel Finest grid level before the regrid operation. 
  */
  virtual void
  preRegrid(const int a_lbase, const int a_oldFinestLevel) = 0;

  /*!
    @brief Regrid this solver. 
    @param[in] a_lmin           Coarsest level where grids did not change. 
    @param[in] a_oldFinestLevel Finest AMR level before the regrid. 
    @param[in] a_newFinestLevel Finest AMR level after the regrid. 
  */
  virtual void
  regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) = 0;

  /*!
    @brief Register operators
  */
  virtual void
  registerOperators() = 0;

  /*!
    @brief Fill solver with initial data. By default, this sets internal data to zero. 
    @details Default implementation sets to zero -- override if you want to have different initial conditions. 
  */
  virtual void
  initialData();

  /*!
    @brief Set source term
    @param[in] a_source Source term
  */
  virtual void
  setSource(const EBAMRCellData& a_source);

  /*!
    @brief Set source
    @param[in] a_source Source term
  */
  virtual void
  setSource(const Real a_source);

  /*!
    @brief Set source
    @param[in] a_source Source term (varies in space)
  */
  virtual void
  setSource(const std::function<Real(const RealVect a_pos)> a_source);

  /*!
    @brief Write plot file
  */
  virtual void
  writePlotFile() = 0;

  /*!
    @brief Get number of output fields
    @return Returns number of variables that will be plotted to file. 
  */
  virtual int
  getNumberOfPlotVariables() const;

  /*!
    @brief Get output plot names
  */
  virtual Vector<std::string>
  getPlotVariableNames() const;

  /*!
    @brief Write output data to a_output
    @param[inout] a_output      Output data holder.
    @param[inout] a_icomp       Starting component where this solver begins writing the output. 
    @param[in]    a_outputRealm Realm where a_output belongs
    @param[in]    a_level       Grid level
    @note This routine writes m_phi on centroids (and not cell-centers). 
  */
  virtual void
  writePlotData(LevelData<EBCellFAB>& a_output,
                int&                  a_comp,
                const std::string     a_outputRealm,
                const int             a_level) const noexcept;

#ifdef CH_USE_HDF5
  /*!
    @brief Write checkpoint data into HDF5 file. 
    @paramo[out] a_handle HDF5 file. 
    @param[in]   a_level Grid level
  */
  virtual void
  writeCheckpointLevel(HDF5Handle& a_handle, const int a_level) const = 0;
#endif

#ifdef CH_USE_HDF5
  /*!
    @brief Read checkpoint data from HDF5 file. 
    @param[in] a_handle HDF5 handle.
    @param[in] const int a_level Grid level
  */
  virtual void
  readCheckpointLevel(HDF5Handle& a_handle, const int a_level) = 0;
#endif

  /*!
    @brief Get current time
    @return Returns m_time
  */
  virtual Real
  getTime() const;

  /*!
    @brief Get the RTE phase
    @return Returns m_phase
  */
  virtual phase::which_phase
  getPhase();

  /*!
    @brief Get solver state.
    @return Returns m_phi. For diffusive solves this will have a single component. 
  */
  virtual EBAMRCellData&
  getPhi();

  /*!
    @brief Get multifluid source
    @return Returns m_source. For diffusive solves this will have a single component. 
  */
  virtual EBAMRCellData&
  getSource();

  /*!
    @brief Get the absorption length.
    @return m_kappa
  */
  virtual EBAMRFluxData&
  getKappa();

  /*!
    @brief Get the absorption coefficient on irregular EB faces
    @return m_kappaEB
  */
  virtual EBAMRIVData&
  getKappaEb();

  /*!
    @brief Get species
  */
  virtual RefCountedPtr<RtSpecies>&
  getSpecies();

protected:
  /*!
    @brief Default component that we solve for
  */
  static constexpr int m_comp = 0;

  /*!
    @brief Default number of components
  */
  static constexpr int m_nComp = 1;

  /*!
    @brief Data location
  */
  Location::Cell m_dataLocation;

  /*!
    @brief Realm where this solver lives.
  */
  std::string m_realm;

  /*!
    @brief EBIndexSpace for this solver.
  */
  RefCountedPtr<EBIndexSpace> m_ebis;

  /*!
    @brief Radiative transfer species (contains meta-information like initial conditions)
  */
  RefCountedPtr<RtSpecies> m_rtSpecies;

  /*!
    @brief Computational geometry. 
  */
  RefCountedPtr<ComputationalGeometry> m_computationalGeometry;

  /*!
    @brief AMR; needed for grid stuff
  */
  RefCountedPtr<AmrMesh> m_amr;

  /*!
    @brief Phase
  */
  phase::which_phase m_phase;

  /*!
    @brief Name for this solver. 
  */
  std::string m_name;

  /*!
    @brief Class name -- needed because inherited classes will be named different. 
  */
  std::string m_className;

  /*!
    @brief Cached state used for regridding
  */
  EBAMRCellData m_cachePhi;

  /*!
    @brief Internal state. 
    @details For diffusive models, this will only contain the isotropic density. For higher order models, this should also
    include the flux, Eddigton tensor etc.
  */
  EBAMRCellData m_phi;

  /*!
    @brief Source term 
    @details For diffusive models, this will only contain the isotropic source. For higher order models, this should also 
    include beam sources. 
  */
  EBAMRCellData m_source;

  /*!
    @brief Absorption coefficient. 
  */
  EBAMRFluxData m_kappa;

  /*!
    @brief Absorption coefficient on EB faces. 
  */
  EBAMRIVData m_kappaEB;

  /*!
    @brief Time
  */
  Real m_time;

  /*!
    @brief Time increment
  */
  Real m_dt;

  /*!
    @brief Stationary solver or not
  */
  bool m_stationary;

  /*!
    @brief Output state
  */
  bool m_plotPhi;

  /*!
    @brief Output source term
  */
  bool m_plotSource;

  /*!
    @brief
  */
  int m_verbosity;

  /*!
    @brief Time step
  */
  int m_timeStep;

  /*!
    @brief Set ebis
  */
  void
  setEbIndexSpace(const RefCountedPtr<EBIndexSpace>& a_ebis);

  /*!
    @brief Parse verbosity
  */
  void
  parseVerbosity() noexcept;

  /*!
    @brief Write data to output. Convenience function. 
    @param[inout] a_output Output data holder.
    @param[inout] a_icomp  Starting component where this solver begins writing the output. 
    @param[in]    a_data   Data to write. 
    @param[in]    a_level  Grid level
    @param[in]    a_outputRealm Realm to which a_output belongs
    @param[in]    a_interpToCentroids If true, a_data will be interpolated to cell centroids before writing to a_output. 
    @param[in]    a_interpGhost If true, interpolate ghost cells
  */
  virtual void
  writeData(LevelData<EBCellFAB>& a_output,
            int&                  a_comp,
            const EBAMRCellData&  a_data,
            const std::string     a_outputRealm,
            const int             a_level,
            const bool            a_interpToCentroids,
            const bool            a_interpGhost) const noexcept;
};

#include <CD_NamespaceFooter.H>

#endif
