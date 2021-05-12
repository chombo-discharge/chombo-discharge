/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   FieldSolver.H
  @brief  Contains declaration of an abstract electrostatic field solver. 
  @author Robert Marskar
*/

#ifndef CD_FieldSolver_H
#define CD_FieldSolver_H

// Our includes
#include <computational_geometry.H>
#include <CD_AmrMesh.H>
#include <CD_ElectrostaticDomainBc.H>
#include <CD_ElectrostaticEbBc.H>

#include <CD_NamespaceHeader.H>

/*!
  @brief Abstract electrostatic field solver class. This class contains an interface to a solve routine, data holders, and interface to boundary conditions.
*/
class FieldSolver {
public:
  
  /*!
    @brief Default function for space-time dependence of domain boundary conditions. 
    @details This is a utility function 
    @param[in] a_position Position in space
    @param[in] a_time     Time
    @return This returns one (1.0)
  */
  static Real s_defaultDomainBcFunction(const RealVect a_position, const Real a_time);

  /*!
    @brief Potential of one.
    @todo Move to FieldStepper
  */
  static Real s_voltageOne(const Real a_time);

  /*!
    @brief Constructor
  */
  FieldSolver();

  /*!
    @brief Constructor
    @details This sets the class name to "FieldSolver" and the Realm to Realm::Primal. 
  */
  virtual ~FieldSolver();

  /*!
    @brief Solve Poisson equation using m_potential, m_rho, and m_sigma. 
    @param[in] a_zeroPhi Set m_potential to zero before calling other function. 
    @return True if we found a solution and false otherwise. 
    @note This is just a wrapper which calls the other version. 
  */
  virtual bool solve(const bool a_zeroPhi = false);

  /*!
    @brief Solve Poisson equation onto a_state using m_rho, and m_sigma as right-hand sides. 
    @param[inout] a_potential Potential
    @param[in]    a_zeroPhi Set a_potential to zero before calling other function. 
    @return True if we found a solution and false otherwise. 
    @note This is just a wrapper which calls the other version. 
  */
  virtual bool solve(MFAMRCellData& a_potential, const bool a_zerophi = false);

  /*!
    @brief Solves Poisson equation onto a_state using a_rho and a_sigma as right-hand sides. 
    @param[inout] a_potential Potential
    @param[in]    a_rho Space charge density
    @param[in]    a_sigma Surface charge density. 
    @param[in]    a_zeroPhi Set a_potential to zero first. 
    @return True if we found a solution and false otherwise. 
    @note a_sigma must be defined on the gas phase. 
  */
  virtual bool solve(MFAMRCellData&       a_state,
		     const MFAMRCellData& a_rho,
		     const EBAMRIVData&   a_sigma,
		     const bool           a_zerophi = false) = 0;

  /*!
    @brief Compute the cell-centered electric field. 
    @details This uses m_potential for computing the electric field (put into m_electricField). 
  */
  virtual void computeElectricField();

  /*!
    @brief Compute the cell-centered electric field. 
    @param[out] a_electricField Cell-centered electric field
    @param[in]  a_potential     Cell-centered potential
  */
  virtual void computeElectricField(MFAMRCellData& a_E, const MFAMRCellData& a_potential);
  
  /*!
    @brief Parse options (for derived class)
  */
  virtual void parseOptions() = 0;

  /*!
    @brief Parse runtime options (for derived class)
  */
  virtual void parseRuntimeOptions() = 0;
  
  /*!
    @brief Allocates internal storage for FieldSolver. Derived classes may want to overwrite. 
  */
  virtual void allocateInternals();

  /*!
    @brief Cache state before regridding
    @param[in] a_lbase          Coarsest level which changes during regrid. 
    @param[in] a_oldFinestLevel Finest AMR level before regrid. 
    @details This allocates m_cache which gets a copy of m_potential
  */
  virtual void preRegrid(const int a_lbase, const int a_oldFinestLevel);

  /*!
    @brief Compute displacement field from the electric field
    @param[out] a_displacementField Displacement field
    @param[in]  a_electricField     Electric field
  */
  virtual void computeDisplacementField(MFAMRCellData& a_displacementField, const MFAMRCellData& a_electricField);

  /*!
    @brief Deallocate internal storage
  */
  virtual void deallocateInternals();

  /*!
    @brief Regrid method. 
    @param[in] a_lmin           Coarsest level where grids did not change. 
    @param[in] a_oldFinestLevel Previous finest AMR level
    @param[in] a_newFinestLevel New finest AMR level. 
    @details This linearly interpolates (with limiters) m_potential to the new grids and recomputes the electric field (from the interpolated potential). 
  */
  virtual void regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel);

  /*!
    @brief Register operators for AMR. Derived classes have to implement these. 
  */
  virtual void registerOperators() = 0;

  /*!
    @brief Set the computational geometry
    @param[in] a_compgeom Computational geometry. 
  */
  virtual void setComputationalGeometry(const RefCountedPtr<computational_geometry>& a_compgeom);
  
  /*!
    @brief Set the amr object. 
    @param[in] a_amr AmrMesh object. 
  */
  virtual void setAmr(const RefCountedPtr<AmrMesh>& a_amr);

  /*!
    @brief Write plot file
    @details This writes a plot file named FieldSolver.stepXXXXXXX.2d.hdf5 using the specified plot variables. 
  */
  virtual void writePlotFile();

  /*!
    @brief Write checkpoint data for a level
    @paramo[out] a_handle HDF5 handle.
    @param[in]   a_level Grid level
    @details This writes m_potential[a_level] to the checkpoint file. 
  */
  virtual void writeCheckpointLevel(HDF5Handle& a_handle, const int a_level) const;

  /*!
    @brief Read checkpoint data onto a level
    @param[in] a_handle HDF5 handle.
    @param[in] const int a_level Grid level
    @details This fills m_potential[a_level] with data from a_handle. 
  */
  virtual void readCheckpointLevel(HDF5Handle& a_handle, const int a_level);

  /*!
    @brief Post checkpoint. This does not do anything but derived classes will usually need it. 
  */
  virtual void postCheckpoint();

  /*!
    @brief Write output data to a_output
    @param[inout] a_output Output data holder
    @param[inout] a_comp   Current component in a_output. 
    @return On output, the solver will have written its plot variables to a_output and increment a_comp by the number of plotted variables. 
  */
  virtual void writePlotData(EBAMRCellData& a_output, int& a_comp);

  /*!
    @brief Write multifluid data to single-fluid data holders. 
    @details This takes the valid data from each phase and writes it to a_output. 
    @param[inout] a_output Output data holder
    @param[inout] a_comp   Current component in a_output. 
    @param[in]    a_data   Multifluid data holder
    @param[in]    a_interp Special flag which, if true, tells AmrMesh to interpolate the data to centroid (in irregular cells). 
  */
  virtual void writeMultifluidData(EBAMRCellData& a_output, int& a_comp, const MFAMRCellData& a_data, const bool a_interp);

  /*!
    @brief Set Realm
    @param[in] a_Realm Realm identifier. 
  */
  virtual void setRealm(const std::string a_Realm);

  /*!
    @brief Set time for this solver.
    @param[in] a_timeStep Time step
    @param[in] a_time     Time (in seconds). 
    @param[in] a_dt       Time step size (in seconds). 
    @details This sets m_timeStep, m_time, and m_dt.
  */
  virtual void setTime(const int a_timeStep, const Real a_time, const Real a_dt);

  /*!
    @brief Set verbosity.
    @param[in] a_verbosity Verbosity factor (lower yields less printed output). 
  */
  virtual void setVerbosity(const int a_verbosity);

  /*!
    @brief Set potential dependence in time. 
    @param[in] a_voltage Function pointer which sets the voltage travel curve. 
    @details If you want something more complex, the voltage can be set individually for each electrode. 
  */
  virtual void setVoltage(std::function<Real(const Real a_time)> a_voltage);
  
  /*!
    @brief Boundary condition function on a wall. 
    @param[in] a_dir Coordinate direction.
    @param[in] a_side Side in the coordinate direction.
    @param[in] a_function Boundary condition function.
    @details This sets a boundary condition for a particular side on a wall. The user must also specify how to use this BC in the input script. 
  */
  virtual void setDomainBcWallFunction(const int a_dir, const Side::LoHiSide a_side, const ElectrostaticDomainBc::BcFunction& a_function);

  /*!
    @brief Set embedded boundary Dirichlet function on a specific electrode.
    @param[in] a_electrode Electrode index. Follows the same order as ComputationalGeometry. 
    @param[in] a_function  Voltage on the electrode. 
  */
  virtual void setElectrodeDirichletFunction(const int a_electrode, const ElectrostaticEbBc::BcFunction& a_function);

  /*!
    @brief Get number of output fields.
  */
  virtual int getNumberOfPlotVariables() const;

  /*!
    @brief Get time. 
  */
  virtual Real getTime() const;

  /*!
    @brief Compute the capacitance. 
    @details This will first obtain a solution without any sources, and then compute the energy density. The capacitance is then C = 2*EnergyDensity/(V*V)
    @return Capacitance for system. 
  */
  virtual Real computeCapacitance();

  /*!
    @brief Compute energy density U = 0.5*int(E.dot.D dV)
    @param[in] a_electricField The electric field. 
    @return Energy density
  */
  virtual Real computeEnergyDensity(const MFAMRCellData& a_electricField);

  /*!
    @brief Get output plot names.
    @return All the plot variable names. 
  */
  virtual Vector<std::string> getPlotVariableNames() const;

  /*!
    @brief Get the Realm where this solver is registered. 
  */
  const std::string getRealm() const;
  
  /*!
    @brief Get potential on both phases. 
  */
  MFAMRCellData& getPotential();

  /*!
    @brief Get electric field on both phases. 
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
  std::string m_Realm;
  
  /*!
    @brief Class name
  */
  std::string m_className;

  /*!
    @brief MFindex space
  */
  RefCountedPtr<mfis> m_multifluidIndexSpace;

  /*!
    @brief Computational geometry. Needed for passing in BCs on EBs. 
  */
  RefCountedPtr<computational_geometry> m_compgeom;

  /*!
    @brief AMR; needed for grid stuff
  */
  RefCountedPtr<AmrMesh> m_amr;

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
    @brief Parses which plot variables to write to plot files
  */
  virtual void parsePlotVariables();

  /*!
    @brief Parses domain boundary conditions. 
  */
  virtual void parseDomainBc();

  /*!
    @brief Set default BC functions. This sets all the m_domainBcFunction objects to s_defaultDomainBcFunction, which return 1 everywhere. 
  */
  virtual void setDefaultDomainBcFunctions();

  /*!
    @brief Set default Dirichlet boundary conditions on the embedded boundaries. 
    @details For each electrode the default boundary condition is set to m_voltage*electroc.get_fraction(). 
  */
  virtual void setDefaultEbBcFunctions();

  /*!
    @brief Returns BC type based on string.
    @param[in] a_str Boundary condition string. Must be "dirichlet number", "neumann number", "dirichlet_custom" or "neumann_custom"
  */
  virtual ElectrostaticDomainBc::BcType parseBcString(const std::string a_str) const;

  /*!
    @brief Shortcut for making a boundary condition string. 
    @param[in] a_dir  Direction.
    @param[in] a_side Coordinate side. 
    @return Returns string of type m_class_name.bc.direction.side. 
  */
  virtual std::string makeBcString(const int a_dir, const Side::LoHiSide a_side) const;
};

#include <CD_NamespaceFooter.H>

#endif
