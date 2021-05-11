/*!
  @file FieldSolver.H
  @brief Abstract multi-phase Poisson solver
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef FieldSolver_H
#define FieldSolver_H

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
    @brief Default function for space-time dependence of domain boundary conditions. 
  */
  static Real s_defaultDomainBcFunction(const RealVect a_position, const Real a_time);

  /*!
    @brief Potential of one.
  */
  static Real s_voltageOne(const Real a_time);

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
  virtual void computeElectricField();

  /*!
    @brief Cell-centered E-field computation
  */
  virtual void computeElectricField(MFAMRCellData& a_E, const MFAMRCellData& a_potential);
  
  /*!
    @brief Parse options
  */
  virtual void parseOptions() = 0;

  /*!
    @brief Parse runtime options
  */
  virtual void parseRuntimeOptions() = 0;
  
  /*!
    @brief Allocate internal storage
  */
  virtual void allocateInternals();

  /*!
    @brief Cache state before regridding
  */
  virtual void preRegrid(const int a_lbase, const int a_old_finest_level);

  /*!
    @brief Compute displacement field from E
  */
  virtual void computeDisplacementField(MFAMRCellData& a_D, const MFAMRCellData& a_E);

  /*!
    @brief Deallocate internal storage
  */
  virtual void deallocateInternals();

  /*!
    @brief Regrid method. 
  */
  virtual void regrid(const int a_lmin, const int a_old_finest, const int a_new_finest);

  /*!
    @brief Register operators for AMR
  */
  virtual void registerOperators() = 0;

  /*!
    @brief Set the computational geometry
  */
  virtual void setComputationalGeometry(const RefCountedPtr<computational_geometry>& a_compgeom);
  
  /*!
    @brief Set the amr object
  */
  virtual void setAmr(const RefCountedPtr<amr_mesh>& a_amr);

  /*!
    @brief Write plot file
  */
  virtual void writePlotFile();

  /*!
    @brief Write checkpoint data into handle
  */
  virtual void writeCheckpointLevel(HDF5Handle& a_handle, const int a_level) const;

  /*!
    @brief Read checkpoint data from handle
  */
  virtual void readCheckpointLevel(HDF5Handle& a_handle, const int a_level);

  /*!
    @brief Post checkpoint things
  */
  virtual void postCheckpoint();

  /*!
    @brief Write output data to a_output
  */
  virtual void writePlotData(EBAMRCellData& a_output, int& a_comp);

  /*!
    @brief Write mf-data to output
  */
  virtual void writeMultifluidData(EBAMRCellData& a_output, int& a_comp, const MFAMRCellData& a_data, const bool a_interp);

  /*!
    @brief Set realm
  */
  virtual void setRealm(const std::string a_realm);

  /*!
    @brief Get current time
  */
  virtual void setTime(const int a_step, const Real a_time, const Real a_dt);

  /*!
    @brief Set verbosity.
  */
  virtual void setVerbosity(const int a_verbosity);

  /*!
    @brief Set potential dependence in time. 
  */
  virtual void setVoltage(std::function<Real(const Real a_time)> a_voltage);
  
  /*!
    @brief Boundary condition function on a wall. 
    @param[in] a_dir Coordinate direction.
    @param[in] a_side Side in the coordinate direction.
    @param[in] a_function Boundary condition function.
  */
  virtual void setDomainBcWallFunction(const int a_dir, const Side::LoHiSide a_side, const ElectrostaticDomainBc::BcFunction& a_function);

  /*!
    @brief Set embedded boundary Dirichlet function on a specific electrode 
  */
  virtual void setElectrodeDirichletFunction(const int a_electrode, const ElectrostaticEbBc::BcFunction& a_function);

  /*!
    @brief Set potential in covered electrode cells
  */
  virtual void setCoveredPotential(EBAMRCellData& a_phi, const int a_comp, const Real a_time);

  /*!
    @brief Get number of output fields
  */
  virtual int getNumberOfPlotVariables() const;

  /*!
    @brief Get current time
  */
  virtual Real getTime() const;

  /*!
    @brief Compute the capacitance. This assumes that m_state is a current solution to your Poisson problem!
  */
  virtual Real computeCapacitance();

  /*!
    @brief Compute energy density U = 0.5*int(E.dot.D dV)
  */
  virtual Real computeEnergyDensity(const MFAMRCellData& a_electricField);

  /*!
    @brief Get output plot names
  */
  virtual Vector<std::string> getPlotVariableNames() const;

  /*!
    @brief Get the realm
  */
  const std::string getRealm() const;
  
  /*!
    @brief Get multifluid state
  */
  MFAMRCellData& getPotential();

  /*!
    @brief Get electric field
  */
  MFAMRCellData& getElectricField();

  /*!
    @brief Get E-field on specific phase
  */
  EBAMRCellData& getElectricField(const phase::which_phase a_phase);

  /*!
    @brief Get space charge density. 
  */
  MFAMRCellData& getRho();

  /*!
    @brief Get the residue
  */
  MFAMRCellData& getResidue();
  
protected:

  /*!
    @brief Realm
  */
  std::string m_realm;
  
  /*!
    @brief Class name
  */
  std::string m_className;

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
  MFAMRCellData m_potential;

  /*!
    @brief E-field
  */
  MFAMRCellData m_electricField;

  /*!
    @brief Space charge density. 
  */
  MFAMRCellData m_rho;

  /*!
    @brief Gas-side E-field
  */
  EBAMRCellData m_electricFieldGasSide;

  /*!
    @brief Solid-side E-field
  */
  EBAMRCellData m_electricFieldSolidSide;

  /*!
    @brief Surface charge density. 
  */
  EBAMRIVData m_sigma;

  /*!
    @brief Residue after solving Poisson equation. 
  */
  MFAMRCellData m_residue;

  /*!
    @brief Output state
  */
  bool m_plotPotential;

  /*!
    @brief Output rhs
  */
  bool m_plotRho;

  /*!
    @brief Output electric field
  */
  bool m_plotElectricField;

  /*!
    @brief Output GMG residue
  */
  bool m_plotResidue;

  /*!
    @brief Verbosity
  */
  int m_verbosity;

  /*!
    @brief Time step
  */
  int m_timeStep;

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
  std::function<Real(const Real a_time)> m_voltage;

  /*!
    @brief Set which variables to output
  */
  virtual void parsePlotVariables();

  /*!
    @brief Parse domain BCs
  */
  virtual void parseDomainBc();

  /*!
    @brief Set default BC functions. This sets all the m_domainBcFunction objects to s_defaultDomainBcFunction (returns 1 everywhere). 
  */
  virtual void setDefaultDomainBcFunctions();

  /*!
    @brief Set default Dirichlet boundary conditions on the embedded boundaries
  */
  virtual void setDefaultEbBcFunctions();

  /*!
    @brief Returns BC type based on string
  */
  virtual ElectrostaticDomainBc::BcType parseBcString(const std::string a_str) const;

  /*!
    @brief Shortcut for making a boundary condition string. This just returns m_class_name.bc.direction.side
  */
  virtual std::string makeBcString(const int a_dir, const Side::LoHiSide a_side) const;
};

#include "CD_NamespaceFooter.H"

#endif
