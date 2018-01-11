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

/*!
  @brief Abstract RTE solver class. This is essentially just an interface. 
*/
class rte_solver {
public:

  /*!
    @brief Full constructor
  */
  rte_solver(const RefCountedPtr<computational_geometry> a_compgeom);

  /*!
    @brief Constructor
  */
  virtual ~rte_solver();

  /*!
    @brief Advance equation one time step
  */
  virtual void advance(const Real& a_dt) = 0;

  /*!
    @brief Advance RTE onto isotropic state a_state
  */
  virtual void advance(const Real& a_dt, EBAMRCellData& a_state) = 0;

  /*!
    @brief Advance RTE onto state a_state
  */
  virtual void advance(const Real& a_dt, EBAMRCellData& a_state, const EBAMRCellData& a_source) = 0;

  /*!
    @brief Compute the boundary flux given a state (a_state will be different for different RTE approximations)
  */
  virtual void compute_boundary_flux(EBAMRIVData& a_bc_flux, const EBAMRCellData& a_state) = 0;

  /*!
    @brief Return required number of ghost cells
  */
  virtual int query_ghost() const = 0;

  /*!
    @brief Get multifluid state
  */
  EBAMRCellData& get_isotropic_state();

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

  /*!
    @brief Get current time
  */
  Real get_time() const;

  /*!
    @brief Check if solver is stationary
  */
  bool is_stationary() const;

  /*!
    @brief Get current time
  */
  void set_time(const Real a_time);
  
protected:
  
  /*!
    @brief Constructor
  */
  rte_solver();

  /*!
    @brief Set computational geometry
  */
  void set_computational_geometry(const RefCountedPtr<computational_geometry> a_compgeom);

  /*!
    @brief Computational geometry. 
  */
  RefCountedPtr<computational_geometry> m_compgeom;

  /*!
    @brief Isotropic photon density
  */
  EBAMRCellData m_isotropic_state;

  /*!
    @brief Source term 
  */
  EBAMRCellData m_source;

  /*!
    @brief Absorption coefficient. I can't image anything else than storing this on faces. 
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
};
#endif
