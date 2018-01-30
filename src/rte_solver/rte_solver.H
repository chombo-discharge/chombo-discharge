/*!
  @file rte_solver.H
  @brief Abstract multi-phase Poisson solver
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _RTE_SOLVER_
#define _RTE_SOLVER_

#include "type_definitions.H"
#include "computational_geometry.H"
#include "amr_mesh.H"
#include "physical_domain.H"
#include "photon_group.H"

/*!
  @brief Abstract RTE solver class. This is essentially just an interface. 
*/
class rte_solver {
public:
  
  /*!
    @brief Constructor
  */
  rte_solver();

  /*!
    @brief Constructor
  */
  virtual ~rte_solver();

  /*!
    @brief Set the photon group
  */
  virtual void set_photon_group(const RefCountedPtr<photon_group> a_photon_group);

  /*!
    @brief Set computational geometry
  */
  virtual void set_computational_geometry(const RefCountedPtr<computational_geometry> a_compgeom);

  /*!
    @brief Set the physical domain
  */
  virtual void set_physical_domain(const RefCountedPtr<physical_domain>& a_physdom);

  /*!
    @brief Set the amr object
  */
  virtual void set_amr(const RefCountedPtr<amr_mesh>& a_amr);

  /*!
    @brief Set phase
    @details This must be done BEFORE callilng set_computational_geometry
  */
  virtual void set_phase(phase::which_phase a_phase = phase::gas);

    /*!
    @brief Set verbosity.
  */
  virtual void set_verbosity(const int a_verbosity);

  /*!
    @brief Get current time
  */
  virtual void set_time(const Real a_time);

  /*!
    @brief Set stationary solver or not
  */
  virtual void set_stationary(const bool a_stationary);

  /*!
    @brief Sanity check
  */
  virtual void sanity_check();

  /*!
    @brief Advance equation one time step
  */
  virtual bool advance(const Real a_dt, const bool a_zerophi = false);

  /*!
    @brief Advance RTE onto isotropic state a_state
  */
  virtual bool advance(const Real a_dt, EBAMRCellData& a_state, const bool a_zerophi = false) = 0;

  /*!
    @brief Compute the boundary flux given a state (a_state will be different for different RTE approximations)
  */
  virtual void compute_boundary_flux(EBAMRIVData& a_ebflux, const EBAMRCellData& a_state) = 0;

  /*!
    @brief Compute the flux
    @details For diffusive models, the flux will be something like grad(a_state). For higher-order models, the flux 
    will be contained in a_state (somehow);
  */
  virtual void compute_flux(EBAMRCellData& a_flux, const EBAMRCellData& a_state) = 0;

  /*!
    @brief Get isotropic part. 
    @details This is normally just a copy, but I need the function signature in order to have a clear plasma-coupling interface
  */
  virtual void compute_density(EBAMRCellData& a_isotropic, const EBAMRCellData& a_state) = 0;

  /*!
    @brief Allocate internal storage
  */
  virtual void allocate_internals() = 0;

  /*!
    @brief Set source term
  */
  virtual void set_source(EBAMRCellData& a_source);

  /*!
    @brief Set source
  */
  virtual void set_source(const Real a_source);

#ifdef CH_USE_HDF5
  /*!
    @brief Write plot file
  */
  virtual void write_plot_file(const int a_step = 0) = 0;
#endif

  /*!
    @brief Return required number of ghost cells
  */
  virtual int query_ghost() const = 0;

  /*!
    @brief Get time
  */
  virtual Real get_time() const;

  /*!
    @brief Get multifluid state
  */
  EBAMRCellData& get_state();

  /*!
    @brief Get multifluid source
  */
  EBAMRCellData& get_source();

  /*!
    @brief Set the absorption coefficient
  */
  EBAMRFluxData& get_kappa();

  /*!
    @brief Set the absorption coefficient on irregular EB faces
  */
  EBAMRIVData& get_kappa_eb();
  
protected:

  /*!
    @brief Phase 
  */
  RefCountedPtr<EBIndexSpace> m_ebis;

  /*!
    @brief Photon group
  */
  RefCountedPtr<photon_group> m_photon_group;

  /*!
    @brief Computational geometry. 
  */
  RefCountedPtr<computational_geometry> m_compgeom;

  /*!
    @brief Physical domain
  */
  RefCountedPtr<physical_domain> m_physdom;

  /*!
    @brief AMR; needed for grid stuff
  */
  RefCountedPtr<amr_mesh> m_amr;

  /*!
    @brief Phase
  */
  phase::which_phase m_phase;

  /*!
    @brief Name
  */
  std::string m_name;

  /*!
    @brief Internal state. 
    @details For diffusive models, this will only contain the isotropic density. For higher order models, this should also
    include the flux, Eddigton tensor etc.
  */
  EBAMRCellData m_state;

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
    @brief Stationary solver or not
  */
  bool m_stationary;

  /*!
    @brief
  */
  int m_verbosity;

  /*!
    @brief Set ebis
  */
  virtual void set_ebis(const RefCountedPtr<EBIndexSpace>& a_ebis);

};
#endif
