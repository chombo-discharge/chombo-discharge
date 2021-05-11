/*!
  @file FieldSolver.H
  @brief Abstract multi-phase Poisson solver
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _FIELD_SOLVER_
#define _FIELD_SOLVER_

#include "computational_geometry.H"
#include "amr_mesh.H"
#include "CD_ElectrostaticDomainBc.H"
#include "CD_ElectrostaticEbBc.H"
#include "CD_NamespaceHeader.H"

/*!
  @brief Abstract Poisson solver class. Boundary conditions are passed through this class. 
*/
class FieldSolver {
public:

  
  /*!
    @brief Silly enum for describing the potential at domain faces
  */
  enum class potential {
    ground,
    live,
    mapped // Not supported yet...?
  };

  /*!
    @brief Default function for space-time dependence of domain boundary conditions. 
  */
  static Real s_defaultDomainBcFunction(const RealVect a_position, const Real a_time);

  /*!
    @brief Potential of one.
  */
  static Real s_potential_one(const Real a_time);

  /*!
    @brief Constructor
  */
  FieldSolver();

  /*!
    @brief Constructor
  */
  virtual ~FieldSolver();

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
    @brief Cell-centered E-field computation
  */
  virtual void compute_E();

  /*!
    @brief Cell-centered E-field computation
  */
  virtual void compute_E(MFAMRCellData& a_E, const MFAMRCellData& a_potential);
  
  /*!
    @brief Parse options
  */
  virtual void parse_options() = 0;

  /*!
    @brief Parse runtime options
  */
  virtual void parse_runtime_options() = 0;
  
  /*!
    @brief Allocate internal storage
  */
  virtual void allocate_internals();

  /*!
    @brief Automatic tuning of solver
  */
  virtual void auto_tune() = 0;

  /*!
    @brief Cache state before regridding
  */
  virtual void pre_regrid(const int a_lbase, const int a_old_finest_level);

  /*!
    @brief Compute displacement field from E
  */
  virtual void compute_D(MFAMRCellData& a_D, const MFAMRCellData& a_E);

  /*!
    @brief Deallocate internal storage
  */
  virtual void deallocate_internals();

  /*!
    @brief Regrid method. 
  */
  virtual void regrid(const int a_lmin, const int a_old_finest, const int a_new_finest);

  /*!
    @brief Register operators for AMR
  */
  virtual void register_operators() = 0;

  /*!
    @brief Set the computational geometry
  */
  virtual void set_computational_geometry(const RefCountedPtr<computational_geometry>& a_compgeom);
  
  /*!
    @brief Set the amr object
  */
  virtual void set_amr(const RefCountedPtr<amr_mesh>& a_amr);

  /*!
    @brief Set which variables to output
  */
  virtual void set_output_variables();

  /*!
    @brief Write plot file
  */
  virtual void write_plot_file();

  /*!
    @brief Write checkpoint data into handle
  */
  virtual void write_checkpoint_level(HDF5Handle& a_handle, const int a_level) const;

  /*!
    @brief Read checkpoint data from handle
  */
  virtual void read_checkpoint_level(HDF5Handle& a_handle, const int a_level);

  /*!
    @brief Post checkpoint things
  */
  virtual void post_checkpoint();

  /*!
    @brief Write output data to a_output
  */
  virtual void write_plot_data(EBAMRCellData& a_output, int& a_comp);

  /*!
    @brief Write mf-data to output
  */
  virtual void write_mfdata(EBAMRCellData& a_output, int& a_comp, const MFAMRCellData& a_data, const bool a_interp);

  /*!
    @brief Set realm
  */
  virtual void set_realm(const std::string a_realm);

  /*!
    @brief Get current time
  */
  virtual void set_time(const int a_step, const Real a_time, const Real a_dt);

  /*!
    @brief Set verbosity.
  */
  virtual void set_verbosity(const int a_verbosity);

  /*!
    @brief Set potential dependence in time. 
  */
  virtual void set_potential(std::function<Real(const Real a_time)> a_potential);
  
  /*!
    @brief Boundary condition function on a wall. 
    @param[in] a_dir Coordinate direction.
    @param[in] a_side Side in the coordinate direction.
    @param[in] a_function Boundary condition function.
  */
  virtual void set_domain_bc_wall_function(const int a_dir, const Side::LoHiSide a_side, const ElectrostaticDomainBc::BcFunction& a_function);

  /*!
    @brief Set embedded boundary Dirichlet function on a specific electrode 
  */
  virtual void set_electrode_dirichlet_function(const int a_electrode, const ElectrostaticEbBc::BcFunction& a_function);

  /*!
    @brief Set potential in covered electrode cells
  */
  virtual void set_covered_potential(EBAMRCellData& a_phi, const int a_comp, const Real a_time);

  /*!
    @brief Query the required number of ghost cells
  */
  virtual int query_ghost() const = 0;

  /*!
    @brief Get number of output fields
  */
  virtual int get_num_plotvars() const;

  /*!
    @brief Get current time
  */
  virtual Real get_time() const;

  /*!
    @brief Compute the capacitance. This assumes that m_state is a current solution to your Poisson problem!
  */
  virtual Real compute_capacitance();

  /*!
    @brief Compute energy density U = 0.5*int(E.dot.D dV)
  */
  virtual Real compute_U(const MFAMRCellData& a_E);

  /*!
    @brief Get output plot names
  */
  virtual Vector<std::string> get_plotvar_names() const;

  /*!
    @brief Get the realm
  */
  const std::string get_realm() const;
  
  /*!
    @brief Get multifluid state
  */
  MFAMRCellData& get_state();

  /*!
    @brief Get electric field
  */
  MFAMRCellData& get_E();

  /*!
    @brief Get E-field on specific phase
  */
  EBAMRCellData& get_E(const phase::which_phase a_phase);

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
    @brief Realm
  */
  std::string m_realm;
  
  /*!
    @brief Class name
  */
  std::string m_class_name;

  /*!
    @brief MFindex space
  */
  RefCountedPtr<mfis> m_mfis;

  /*!
    @brief Computational geometry. Needed for passing in BCs on EBs. 
  */
  RefCountedPtr<computational_geometry> m_compgeom;

  /*!
    @brief AMR; needed for grid stuff
  */
  RefCountedPtr<amr_mesh> m_amr;

  /*!
    @brief Cached state used for regridding. 
  */
  MFAMRCellData m_cache;

  /*!
    @brief State data (i.e. the potential)
  */
  MFAMRCellData m_state;

  /*!
    @brief E-field
  */
  MFAMRCellData m_E;

  /*!
    @brief State data 
  */
  MFAMRCellData m_source;

  /*!
    @brief Gas-side E-field
  */
  EBAMRCellData m_E_gas;

  /*!
    @brief Solid-side E-field
  */
  EBAMRCellData m_E_sol;

  /*!
    @brief Surface source
  */
  EBAMRIVData m_sigma;

  /*!
    @brief Final residue
  */
  MFAMRCellData m_resid;

  /*!
    @brief Automatic tuning
  */
  bool m_autotune;

  /*!
    @brief Output state
  */
  bool m_plot_phi;

  /*!
    @brief Output rhs
  */
  bool m_plot_rho;

  /*!
    @brief Output electric field
  */
  bool m_plot_E;

  /*!
    @brief Output GMG residue
  */
  bool m_plot_res;

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
    @brief Domain boundary conditions for FieldSolver
  */
  ElectrostaticDomainBc m_domainBc;

  /*!
    @brief Dirichlet boundary conditions on electrodes
  */
  ElectrostaticEbBc m_ebBc;

  /*!
    @brief Domain BC functions
  */
  std::map<ElectrostaticDomainBc::Wall, ElectrostaticDomainBc::BcFunction> m_domainBcFunctions;

  /*!
    @brief BC functions (Dirichlet) on the electrodes
  */
  std::vector<std::pair<electrode, ElectrostaticEbBc::BcFunction> > m_electrodeBcFunctions;

  /*!
    @brief Potential function. 
  */
  std::function<Real(const Real a_time)> m_potential;

  /*!
    @brief Set the mfis
  */
  virtual void set_mfis(const RefCountedPtr<mfis>& a_mfis);

  /*!
    @brief Parse domain BCs
  */
  virtual void parse_domain_bc();

  /*!
    @brief Set default BC functions. This sets all the m_domainBcFunction objects to s_defaultDomainBcFunction (returns 1 everywhere). 
  */
  virtual void set_default_domain_bc_functions();

  /*!
    @brief Set default Dirichlet boundary conditions on the embedded boundaries
  */
  virtual void set_default_eb_bc_functions();

  /*!
    @brief Returns BC type based on string
  */
  virtual ElectrostaticDomainBc::BcType parse_bc_string(const std::string a_str) const;

  /*!
    @brief Shortcut for making a boundary condition string. This just returns m_class_name.bc.direction.side
  */
  virtual std::string make_bc_string(const int a_dir, const Side::LoHiSide a_side) const;
};
#include "CD_NamespaceFooter.H"

#endif
