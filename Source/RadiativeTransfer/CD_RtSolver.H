/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_RtSolver.H
  @brief  Abstract parent class for various radiative transfer solvers
  @author Robert Marskar
*/

#ifndef CD_RtSolver_H
#define CD_RtSolver_H

// Our includes
#include <computational_geometry.H>
#include <CD_AmrMesh.H>
#include <rte_species.H>
#include <CD_NamespaceHeader.H>
  
/*!
  @brief Abstract RTE solver class. This is essentially just an interface. 
*/
class RtSolver {
public:
  
  /*!
    @brief Constructor
  */
  RtSolver();

  /*!
    @brief Constructor
  */
  virtual ~RtSolver();

  /*!
    @brief Get solver name
  */
  virtual std::string getName();

  /*!
    @brief Get Realm
  */
  virtual const std::string getRealm() const;

  /*!
    @brief Parse options
  */
  virtual void parseOptions() = 0;

  /*!
    @brief Parse options
  */
  virtual void parseRuntimeOptions() = 0;

  /*!
    @brief Get output plot names
  */
  virtual Vector<std::string> getPlotVariableNames() const;
  
  /*!
    @brief Advance equation one time step
  */
  virtual bool advance(const Real a_dt, const bool a_zerophi = false);

  /*!
    @brief Advance method
  */
  virtual bool advance(const Real a_dt, EBAMRCellData& a_phi, const bool a_zerophi = false);

  /*!
    @brief Advance RTE onto isotropic state a_phi
  */
  virtual bool advance(const Real a_dt, EBAMRCellData& a_phi, const EBAMRCellData& a_source, const bool a_zerophi = false) = 0;

  /*!
    @brief Set Realm
  */
  virtual void setRealm(const std::string a_realm);

  /*!
    @brief Set the photon group
  */
  virtual void setRtSpecies(const RefCountedPtr<rte_species> a_rte_species);

  /*!
    @brief Set computational geometry
  */
  virtual void setComputationalGeometry(const RefCountedPtr<computational_geometry> a_computationalGeometry);

  /*!
    @brief Set the amr object
  */
  virtual void setAmr(const RefCountedPtr<AmrMesh>& a_amr);

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
    @brief Set stationary solver or not
  */
  virtual void setStationary(const bool a_stationary);

  /*!
    @brief Sanity check
  */
  virtual void sanityCheck();

  /*!
    @brief Check if solver is stationary
  */
  virtual bool isStationary();

  /*!
    @brief Compute the boundary flux given a state (a_phi will be different for different RTE approximations)
  */
  virtual void computeBoundaryFlux(EBAMRIVData& a_ebFlux, const EBAMRCellData& a_phi) = 0;

  /*!
    @brief Compute the domain flux given a state (a_phi will be different for different RTE approximations)
  */
  virtual void computeDomainFlux(EBAMRIFData& a_domainflux, const EBAMRCellData& a_phi) = 0;

  /*!
    @brief Compute the flux
    @details For diffusive models, the flux will be something like grad(a_phi). For higher-order models, the flux 
    will be contained in a_phi (somehow);
  */
  virtual void computeFlux(EBAMRCellData& a_flux, const EBAMRCellData& a_phi) = 0;

  /*!
    @brief Get isotropic part. 
    @details This is normally just a copy, but I need the function signature in order to have a clear plasma-coupling interface
  */
  virtual void computeDensity(EBAMRCellData& a_isotropic, const EBAMRCellData& a_phi) = 0;

  /*!
    @brief Allocate internal storage
  */
  virtual void allocateInternals() = 0;

  /*!
    @brief Deallocate internal storage
  */
  virtual void deallocateInternals() = 0;

  /*!
    @brief Pre regrid
  */
  virtual void preRegrid(const int a_lbase, const int a_oldFinestLevel) = 0;

  /*!
    @brief Regrid function
  */
  virtual void regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) = 0;

  /*!
    @brief Register operators
  */
  virtual void registerOperators() = 0;

  /*!
    @brief Fill solver with initial data. By default, this sets internal data to zero. 
    @details Only relevant for non-stationary solvers
  */
  virtual void initialData();

  /*!
    @brief Set source term
  */
  virtual void setSource(const EBAMRCellData& a_source);

  /*!
    @brief Set source
  */
  virtual void setSource(const Real a_source);

  /*!
    @brief Write plot file
  */
  virtual void writePlotFile() = 0;

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
  virtual void writeCheckpointLevel(HDF5Handle& a_handle, const int a_level) const = 0;

  /*!
    @brief Read checkpoint data from handle
  */
  virtual void readCheckpointLevel(HDF5Handle& a_handle, const int a_level) = 0;

  /*!
    @brief Return required number of ghost cells
  */
  virtual int queryGhost() const = 0;

  /*!
    @brief Get number of output fields
  */
  virtual int getNumberOfPlotVariables() const;

  /*!
    @brief Get time
  */
  virtual Real getTime() const;

  /*!
    @brief Get the RTE phase
  */
  virtual phase::which_phase get_phase();

  /*!
    @brief Get multifluid state
  */
  virtual EBAMRCellData& getPhi();

  /*!
    @brief Get multifluid source
  */
  virtual EBAMRCellData& getSource();

  /*!
    @brief Set the absorption coefficient
  */
  virtual EBAMRFluxData& getKappa();

  /*!
    @brief Set the absorption coefficient on irregular EB faces
  */
  virtual EBAMRIVData& getKappaEb();

  /*!
    @brief Get species
  */
  virtual RefCountedPtr<rte_species>& getSpecies();
  
protected:

  /*!
    @brief Realm
  */
  std::string m_realm;

  /*!
    @brief Phase 
  */
  RefCountedPtr<EBIndexSpace> m_ebis;

  /*!
    @brief Photon group
  */
  RefCountedPtr<rte_species> m_rte_species;

  /*!
    @brief Computational geometry. 
  */
  RefCountedPtr<computational_geometry> m_computationalGeometry;

  /*!
    @brief AMR; needed for grid stuff
  */
  RefCountedPtr<AmrMesh> m_amr;

  /*!
    @brief Phase
  */
  phase::which_phase m_phase;

  /*!
    @brief Name
  */
  std::string m_name;

  /*!
    @brief Class name
  */
  std::string m_className;

  /*!
    @brief Cached state used for regridding
  */
  EBAMRCellData m_cache;

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
  EBAMRIVData m_kappa_eb;

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
  void setEbIndexSpace(const RefCountedPtr<EBIndexSpace>& a_ebis);
};

#include <CD_NamespaceFooter.H>

#endif
