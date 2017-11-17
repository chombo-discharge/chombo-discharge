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
    @brief Set up Poisson solver
  */
  virtual void setup_base(const computational_geometry* a_compgeom);

  /*!
    @brief Solve Poisson equation
  */
  virtual void solve() = 0;

  /*!
    @brief Solve Poisson onto state 
  */
  virtual void solve(MFAMRCellData& a_state) = 0;

  /*!
    @brief Solve Poisson onto state 
  */
  virtual void solve(MFAMRCellData& a_state, const MFAMRCellData& a_source) = 0;

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
    @brief Set potential at current time. 
  */
  void set_potential(const Real a_potential);
  
protected:

  /*!
    @brief Computational geometry. Needed for passing in BCs on EBs. 
  */
  const computational_geometry* m_compgeom;

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
    @brief Alias internal states so ensure pointers point where they are supposed to
  */
  void alias_internals();
};
#endif
