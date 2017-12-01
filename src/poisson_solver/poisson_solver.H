/*!
  @file poisson_solver.H
  @brief Abstract multi-phase Poisson solver
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _POISSON_SOLVER_
#define _POISSON_SOLVER_

#include "type_definitions.H"
#include "computational_geometry.H"
#include "amr_mesh.H"
#include "wall_bc.H"

namespace Potential {
  enum GroundLive {
    Ground = false,
    Live   = true
  };
};

/*!
  @brief Abstract Poisson solver class. Just an interface, so pretty lightweight stuff.
*/
class poisson_solver {
public:

  /*!
    @brief Constructor
  */
  poisson_solver();

  /*!
    @brief Constructor
  */
  virtual ~poisson_solver();

  /*!
    @brief Set the computational geometry
  */
  virtual void set_computational_geometry(const RefCountedPtr<computational_geometry>& a_compgeom);

    /*!
    @brief Set the physical domain
  */
  virtual void set_physical_domain(const RefCountedPtr<physical_domain>& a_physdom);

  /*!
    @brief Set the amr object
  */
  virtual void set_amr(const RefCountedPtr<amr_mesh>& a_amr);

  /*!
    @brief Solve Poisson equation
  */
  virtual void solve();

  /*!
    @brief Solve Poisson onto state 
  */
  virtual void solve(MFAMRCellData& a_state);

  /*!
    @brief Solve Poisson onto state 
  */
  virtual void solve(MFAMRCellData& a_state, const MFAMRCellData& a_source) = 0;

  /*!
    @brief Query the required number of ghost cells
  */
  virtual int query_ghost() const = 0;

  /*!
    @brief Set Dirichlet wall BC
  */
  virtual void set_dirichlet_wall_bc(const int a_dir, Side::LoHiSide a_side, const Potential::GroundLive a_live);

  /*!
    @brief Set Neumann wall BC
  */
  virtual void set_neumann_wall_bc(const int a_dir, Side::LoHiSide a_side, const Real a_value);

  /*!
    @brief Get Poisson BC at specified domain wall
  */
  virtual wall_bc& get_wall_bc(const int a_dir, Side::LoHiSide a_side) const;

  /*!
    @brief Sanity check
  */
  virtual void sanity_check();
  
  /*!
    @brief Get multifluid state
  */
  MFAMRCellData& get_state();

  /*!
    @brief Get multifluid source
  */
  MFAMRCellData& get_source();

  /*!
    @brief Get state on specified phase. This uses aliasing
  */
  EBAMRCellData& get_state_phase(Phase::WhichPhase a_phase);

  /*!
    @brief Get source on specified phase. This uses aliasing
  */
  EBAMRCellData& get_source_phase(Phase::WhichPhase a_phase);

  /*!
    @brief Get the jump state
  */
  EBAMRIVData& get_jump();

  /*!
    @brief Get current time
  */
  Real get_time() const;

  /*!
    @brief Get current time
  */
  void set_time(const Real a_time);

  /*!
    @brief Set verbosity.
  */
  void set_verbosity(const int a_verbosity);

  /*!
    @brief Set potential at current time. 
  */
  void set_potential(const Real a_potential);
  
protected:

  /*!
    @brief MFindex space
  */
  RefCountedPtr<mfis> m_mfis;

  /*!
    @brief Computational geometry. Needed for passing in BCs on EBs. 
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
    @brief Wall boundary conditions for Poisson equation. Must be set by the plasma engine. 
  */
  Vector<RefCountedPtr<wall_bc> > m_wallbc;

  /*!
    @brief State data (i.e. the potential)
  */
  MFAMRCellData m_state;

  /*!
    @brief State data (i.e. the 
  */
  MFAMRCellData m_source;

  /*!
    @brief State at Phase::Gas
  */
  EBAMRCellData m_state_gas;

  /*!
    @brief State at
  */
  EBAMRCellData m_state_solid;

    /*!
    @brief State at Phase::Gas
  */
  EBAMRCellData m_source_gas;

  /*!
    @brief State at
  */
  EBAMRCellData m_source_solid;
  
  /*!
    @brief Jump state on interface - this gives the BC on the gas-solid interfaces
  */
  EBAMRIVData m_jump;

  /*!
    @brief Function pointer to potential function - this gives the BC on EBs (where we have Dirichlet)
  */
  Real m_potential;

  /*!
    @brief Time
  */
  Real m_time;

  /*!
    @brief Current potential. For passing into BC objects. 
  */
  Real m_current_potential;

  /*!
    @brief Verbosity
  */
  int m_verbosity;

  /*!
    @brief Alias internal states so ensure pointers point where they are supposed to
  */
  virtual void alias_internals();

  /*!
    @brief Set the mfis
  */
  virtual void set_mfis(const RefCountedPtr<mfis>& a_mfis);
  
  /*!
    @brief Allocate data holders for wall bcs so they can be set
  */
  virtual void allocate_wall_bc();
};
#endif
