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

namespace potential {
  enum ground_live {
    ground = 0,
    live   = 1,
    mapped = 2, // Not currently supported. But it will. 
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
    @brief Solve Poisson equation
  */
  virtual bool solve(const bool a_zerophi = false);

  /*!
    @brief Solve Poisson onto state 
  */
  virtual bool solve(MFAMRCellData& a_state, const bool a_zerophi = false);

  /*!
    @brief Solve Poisson onto state 
  */
  virtual bool solve(MFAMRCellData&       a_state,
		     const MFAMRCellData& a_rho,
		     const EBAMRIVData&   a_sigma,
		     const bool           a_zerophi = false) = 0;

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
    @brief Allocate internal storage
  */
  virtual void allocate_internals();

  /*!
    @brief Set Dirichlet wall BC
  */
  virtual void set_dirichlet_wall_bc(const int a_dir, Side::LoHiSide a_side, const potential::ground_live a_live);

  /*!
    @brief Set Neumann wall BC
  */
  virtual void set_neumann_wall_bc(const int a_dir, Side::LoHiSide a_side, const Real a_value);

#ifdef CH_USE_HDF5
  /*!
    @brief Write plot file
  */
  virtual void write_plot_file();
#endif

  /*!
    @brief Sanity check
  */
  virtual void sanity_check();

  /*!
    @brief Get current time
  */
  virtual void set_time(const int a_step, const Real a_time, const Real a_dt);

  /*!
    @brief Set verbosity.
  */
  virtual void set_verbosity(const int a_verbosity);

  /*!
    @brief Set potential at current time. 
  */
  virtual void set_potential(Real (*a_potential)(const Real a_time, const RealVect a_pos));

  /*!
    @brief Query the required number of ghost cells
  */
  virtual int query_ghost() const = 0;

  /*!
    @brief Get current time
  */
  virtual Real get_time() const;

  /*!
    @brief Get Poisson BC at specified domain wall
  */
  virtual wall_bc& get_wall_bc(const int a_dir, Side::LoHiSide a_side) const;
  
  /*!
    @brief Get multifluid state
  */
  MFAMRCellData& get_state();

  /*!
    @brief Get multifluid source
  */
  MFAMRCellData& get_source();

  /*!
    @brief Get the residue
  */
  MFAMRCellData& get_resid();

  
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
    @brief State data 
  */
  MFAMRCellData m_source;

  /*!
    @brief Surface source
  */
  EBAMRIVData m_sigma;

  /*!
    @brief Final residue
  */
  MFAMRCellData m_resid;

    /*!
    @brief Verbosity
  */
  int m_verbosity;

  /*!
    @brief Time step
  */
  int m_step;

  /*!
    @brief Time step increment
  */
  Real m_dt;

  /*!
    @brief Time
  */
  Real m_time;

  /*!
    @brief Set potential function pointer
  */
  Real (*m_potential)(const Real a_time, const RealVect a_pos);

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
