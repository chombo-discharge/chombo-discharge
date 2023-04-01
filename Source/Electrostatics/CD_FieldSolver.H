/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_FieldSolver.H
  @brief  Contains declaration of a base electrostatics solver class. 
  @author Robert Marskar
*/

#ifndef CD_FieldSolver_H
#define CD_FieldSolver_H

// Std includes
#include <functional>

// Our includes
#include <CD_Location.H>
#include <CD_ComputationalGeometry.H>
#include <CD_AmrMesh.H>
#include <CD_ElectrostaticDomainBc.H>
#include <CD_ElectrostaticEbBc.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Base class for electrostatic solvers. 
  @details This class contains an interface to a solve routine, data holders, and interface to boundary conditions. Solver routines are
  supposed to be implemented by derived classes. 
*/
class FieldSolver
{
public:
  /*!
    @brief Constructor
  */
  FieldSolver();

  /*!
    @brief Disallowed copy constructor
    @param[in] a_other Other solver
  */
  FieldSolver(const FieldSolver& a_other) = delete;

  /*!
    @brief Disallowed move constructor
    @param[in] a_other Other solver
  */
  FieldSolver(const FieldSolver&& a_other) = delete;

  /*!
    @brief Disallowed copy assignment operator
    @param[in] a_other Other solver
  */
  FieldSolver&
  operator=(const FieldSolver& a_other) = delete;

  /*!
    @brief Disallowed move assignment operator
    @param[in] a_other Other solver
  */
  FieldSolver&
  operator=(const FieldSolver&& a_other) = delete;

  /*!
    @brief Constructor
  */
  virtual ~FieldSolver();

  /*!
    @brief Set up solver routines
  */
  virtual void
  setupSolver() = 0;

  /*!
    @brief A special routine for when solver permittivities need to change but solver does not.
    @param[in] a_permittivityCell Permittivity on cell center
    @param[in] a_permittivityFace Permittivity on faces
    @param[in] a_permittivityEB   Permittivity on EB
    @note Used by FieldSolverMultigrid to change permittivities/conductivities under the hood of the multigrid solver.
  */
  virtual void
  setSolverPermittivities(const MFAMRCellData& a_permittivityCell,
                          const MFAMRFluxData& a_permittivityFace,
                          const MFAMRIVData&   a_permittivityEB);

  /*!
    @brief Solve Poisson equation using m_potential, m_rho, and m_sigma. 
    @param[in] a_zeroPhi Set m_potential to zero before calling other function. 
    @return    True if we found a solution and false otherwise. 
    @details   This calls the other version. 
  */
  virtual bool
  solve(const bool a_zeroPhi = false);

  /*!
    @brief Solve Poisson equation onto a_phi using m_rho, and m_sigma as right-hand sides. 
    @param[inout] a_potential Potential
    @param[in]    a_zeroPhi Set a_potential to zero before calling other function. 
    @return       True if we found a solution and false otherwise. 
    @note         This calls the other (pure) version.
  */
  virtual bool
  solve(MFAMRCellData& a_potential, const bool a_zerophi = false);

  /*!
    @brief Solves Poisson equation onto a_phi using a_rho and a_sigma as right-hand sides. 
    @param[inout] a_potential Potential
    @param[in]    a_rho Space charge density
    @param[in]    a_sigma Surface charge density. 
    @param[in]    a_zeroPhi Set a_potential to zero first. 
    @return True if we found a solution and false otherwise. 
    @note a_sigma must be defined on the gas phase. 
  */
  virtual bool
  solve(MFAMRCellData& a_phi, const MFAMRCellData& a_rho, const EBAMRIVData& a_sigma, const bool a_zerophi = false) = 0;

  /*!
    @brief   Compute the cell-centered electric field. 
    @details This uses m_potential for computing the electric field and puts the result into m_electricField. 
  */
  virtual void
  computeElectricField();

  /*!
    @brief Compute the cell-centered electric field. 
    @param[out] a_electricField Cell-centered electric field
    @param[in]  a_potential     Cell-centered potential
  */
  virtual void
  computeElectricField(MFAMRCellData& a_E, const MFAMRCellData& a_potential) const = 0;

  /*!
    @brief Compute the face-centered electric field. 
    @param[out] a_electricField Face-centered electric field
    @param[in]  a_potential     Cell-centered potential
  */
  virtual void
  computeElectricField(MFAMRFluxData& a_E, const MFAMRCellData& a_potential) const = 0;

  /*!
    @brief Compute the cell-centered electric field on a specific phase. 
    @param[out] a_electricField Face-centered electric field
    @param[in]  a_phase         Phase
    @param[in]  a_potential     Cell-centered potential
  */
  virtual void
  computeElectricField(EBAMRCellData&           a_E,
                       const phase::which_phase a_phase,
                       const MFAMRCellData&     a_potential) const = 0;

  /*!
    @brief Compute the face-centered electric field on a specific phase. 
    @param[out] a_electricField Face-centered electric field
    @param[in]  a_phase         Phase
    @param[in]  a_potential     Cell-centered potential
  */
  virtual void
  computeElectricField(EBAMRFluxData&           a_E,
                       const phase::which_phase a_phase,
                       const MFAMRCellData&     a_potential) const = 0;

  /*!
    @brief Parse options (for derived class)
  */
  virtual void
  parseOptions() = 0;

  /*!
    @brief Parse runtime options (for derived class)
  */
  virtual void
  parseRuntimeOptions() = 0;

  /*!
    @brief Allocates internal storage for FieldSolver. Derived classes may want to overwrite. 
  */
  virtual void
  allocate();

  /*!
    @brief Cache state before regridding
    @param[in] a_lbase          Coarsest level which changes during regrid. 
    @param[in] a_oldFinestLevel Finest AMR level before regrid. 
    @details This allocates m_cache which gets a copy of m_potential
  */
  virtual void
  preRegrid(const int a_lbase, const int a_oldFinestLevel);

  /*!
    @brief Compute displacement field from the electric field
    @param[out] a_displacementField Displacement field
    @param[in]  a_electricField     Electric field
  */
  virtual void
  computeDisplacementField(MFAMRCellData& a_displacementField, const MFAMRCellData& a_electricField);

  /*!
    @brief Deallocate internal storage
  */
  virtual void
  deallocate();

  /*!
    @brief Regrid method. 
    @param[in] a_lmin           Coarsest level allowed to change. 
    @param[in] a_oldFinestLevel Finest level before the regrid. 
    @param[in] a_newFinestLevel Finest level after the regrid. 
    @details This linearly interpolates (with limiters) m_potential to the new grids and recomputes the electric field (from the interpolated potential). 
  */
  virtual void
  regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel);

  /*!
    @brief Register operators for AMR. Derived classes have to implement these. 
  */
  virtual void
  registerOperators() = 0;

  /*!
    @brief Set space charge to constant value everywhere
    @param[in] a_rho Value of space charge
  */
  void
  setRho(const Real a_rho);

  /*!
    @brief Set space charge to spatially dependent function
    @param[in] a_rho Spatially dependent value of space charge
  */
  void
  setRho(const std::function<Real(const RealVect)>& a_rho);

  /*!
    @brief Set surface charge to specified value
  */
  void
  setSigma(const Real a_sigma);

  /*!
    @brief Set the computational geometry
    @param[in] a_computationalGeometry Computational geometry. 
  */
  void
  setComputationalGeometry(const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry);

  /*!
    @brief Set the amr object. 
    @param[in] a_amr AmrMesh object. 
  */
  void
  setAmr(const RefCountedPtr<AmrMesh>& a_amr);

  /*!
    @brief Set the permittivities
    @details This sets m_permittivityCell, m_permittivityFace, and m_permittivityEB.
  */
  virtual void
  setPermittivities();

  /*!
    @brief Write plot file
    @details This writes a plot file named FieldSolver.stepXXXXXXX.2d.hdf5 using the specified plot variables. 
  */
  virtual void
  writePlotFile();

  /*!
    @brief Write checkpoint data for a level
    @paramo[out] a_handle HDF5 handle.
    @param[in]   a_level Grid level
    @details This writes m_potential[a_level] to the checkpoint file. 
  */
#ifdef CH_USE_HDF5
  virtual void
  writeCheckpointLevel(HDF5Handle& a_handle, const int a_level) const;
#endif

  /*!
    @brief Read checkpoint data onto a level
    @param[in] a_handle HDF5 handle.
    @param[in] const int a_level Grid level
    @details This fills m_potential[a_level] with data from a_handle. 
  */
#ifdef CH_USE_HDF5
  virtual void
  readCheckpointLevel(HDF5Handle& a_handle, const int a_level);
#endif

  /*!
    @brief Post checkpoint method. 
    @details Default implementation does not do anything, but derived classes will usually need it. 
  */
  virtual void
  postCheckpoint();

  /*!
    @brief Write output data to a_output
    @param[inout] a_output Output data holder
    @param[inout] a_comp   Current component in a_output. 
    @param[in]    a_forceNoInterp Forces data to be written using native centering
    @return On output, the solver will have written its plot variables to a_output and increment a_comp by the number of plotted variables. 
  */
  virtual void
  writePlotData(EBAMRCellData& a_output, int& a_comp, const bool a_forceNoInterp = false);

  /*!
    @brief Set Realm
    @param[in] a_realm Realm identifier. 
  */
  void
  setRealm(const std::string a_realm);

  /*!
    @brief Set time for this solver.
    @param[in] a_timeStep Time step
    @param[in] a_time     Time (in seconds). 
    @param[in] a_dt       Time step size (in seconds). 
    @details This sets m_timeStep, m_time, and m_dt.
  */
  void
  setTime(const int a_timeStep, const Real a_time, const Real a_dt);

  /*!
    @brief Set verbosity.
    @param[in] a_verbosity Verbosity factor (lower yields less printed output). 
  */
  void
  setVerbosity(const int a_verbosity);

  /*!
    @brief Set potential dependence in time. 
    @param[in] a_voltage Function pointer which sets the voltage travel curve. 
    @details If you want something more complex, the voltage can be set individually for each electrode using setElectrodeDirichletFunction.
  */
  virtual void
  setVoltage(std::function<Real(const Real a_time)> a_voltage);

  /*!
    @brief Boundary condition function on a domain side
    @param[in] a_dir Coordinate direction.
    @param[in] a_side Side in the coordinate direction.
    @param[in] a_function Boundary condition function.
    @details This sets a boundary condition for a particular side on a wall. The user must also specify how to use this BC in the input script. 
  */
  virtual void
  setDomainSideBcFunction(const int                                a_dir,
                          const Side::LoHiSide                     a_side,
                          const ElectrostaticDomainBc::BcFunction& a_function);

  /*!
    @brief Set embedded boundary Dirichlet function on a specific electrode.
    @param[in] a_electrode electrode index. Follows the same order as ComputationalGeometry. 
    @param[in] a_function  Voltage on the electrode. 
  */
  virtual void
  setElectrodeDirichletFunction(const int a_electrode, const ElectrostaticEbBc::BcFunction& a_function);

  /*!
    @brief Get number of output fields.
    @return Number of plot variables
  */
  virtual int
  getNumberOfPlotVariables() const;

  /*!
    @brief Get voltage function
    @return Returns m_voltage
  */
  const std::function<Real(const Real a_time)>&
  getVoltage() const;

  /*!
    @brief Get current voltage
    @return Evaluates m_voltage(m_time) and returns the result. 
  */
  Real
  getCurrentVoltage() const;

  /*!
    @brief Get time. 
    @return m_time
  */
  Real
  getTime() const;

  /*!
    @brief Compute the capacitance. 
    @details This will first obtain a solution without any sources, and then compute the energy density. The capacitance is then C = 2*EnergyDensity/(V*V)
    @return Capacitance for system. 
  */
  Real
  computeCapacitance();

  /*!
    @brief Compute energy density U = 0.5*int(E.dot.D dV)
    @param[in] a_electricField The electric field. 
    @return Energy density
  */
  Real
  computeEnergy(const MFAMRCellData& a_electricField);

  /*!
    @brief Get output plot names.
    @return All the plot variable names. 
  */
  virtual Vector<std::string>
  getPlotVariableNames() const;

  /*!
    @brief Get computational loads for a level
    @return Loads for each box on a grid level. 
    @note The return vector should have the same order as the boxes in a_dbl. E.g. ret[0] must be the load for a_dbl.boxArray()[0];
    @details The default implementation returns the number of cells in the grid patch as a proxy for the load. 
  */
  virtual Vector<long long>
  computeLoads(const DisjointBoxLayout& a_dbl, const int a_level);

  /*!
    @brief Set the data location for the solver
    @param[in] a_dataLocation Desired data location. 
    @details Only Location::Cell::Center and Location::Cell::Centroid are accepted arguments. 
  */
  void
  setDataLocation(const Location::Cell a_dataLocation);

  /*!
    @brief Get the Realm where this solver is registered. 
    @return Returns realm (as string)
  */
  std::string
  getRealm() const;

  /*!
    @brief Get potential on both phases. 
  */
  MFAMRCellData&
  getPotential();

  /*!
    @brief Get electric field on both phases. 
  */
  MFAMRCellData&
  getElectricField();

  /*!
    @brief Get storage for the space charge density. 
    @details This lives in the FieldSolver so users don't need to allocate their own storage for the space charge. 
  */
  MFAMRCellData&
  getRho();

  /*!
    @brief Get the residue
  */
  MFAMRCellData&
  getResidue();

  /*!
    @brief Get cell-centered permittivity
  */
  MFAMRCellData&
  getPermittivityCell();

  /*!
    @brief Get face-centered permittivity
  */
  MFAMRFluxData&
  getPermittivityFace();

  /*!
    @brief Get irregular b coefficient
  */
  MFAMRIVData&
  getPermittivityEB();

  /*!
    @brief Get m_sigma
  */
  EBAMRIVData&
  getSigma();

  /*!
    @brief Get data location
  */
  Location::Cell
  getDataLocation() const;

  /*!
    @brief Write multifluid data to single-fluid data holders. 
    @details This takes the valid data from each phase and writes it to a_output. On cut-cells, the a_phase flag determines which phase will be plotted.
    @param[inout] a_output Output data holder
    @param[inout] a_comp   Current component in a_output. 
    @param[in]    a_data   Multifluid data holder
    @param[in]    a_phase  Main phase in cut-cells. 
    @param[in]    a_interp Special flag which, if true, tells AmrMesh to interpolate the data to centroids (in irregular cells). 
  */
  virtual void
  writeMultifluidData(EBAMRCellData&           a_output,
                      int&                     a_comp,
                      const MFAMRCellData&     a_data,
                      const phase::which_phase a_phase,
                      const bool               a_interp) const;

protected:
  /*!
    @brief Component number where data is stored
  */
  constexpr static int m_comp = 0;

  /*!
    @brief Number of components in data holders
  */
  constexpr static int m_nComp = 1;

  /*!
    @brief Flag which specifies that data location. 
  */
  Location::Cell m_dataLocation;

  /*!
    @brief Flag which specifies where the permittivities are stored
  */
  Location::Face m_faceLocation;

  /*!
    @brief Realm where this solver is registered. 
  */
  std::string m_realm;

  /*!
    @brief Class name (i.e., "FieldSolver" for the base class)
  */
  std::string m_className;

  /*!
    @brief Multifluid index space. 
  */
  RefCountedPtr<MultiFluidIndexSpace> m_multifluidIndexSpace;

  /*!
    @brief Computational geometry. 
  */
  RefCountedPtr<ComputationalGeometry> m_computationalGeometry;

  /*!
    @brief AMR - needed for pretty much everything. 
  */
  RefCountedPtr<AmrMesh> m_amr;

  /*!
    @brief Cached state used for regridding. 
  */
  MFAMRCellData m_cache;

  /*!
    @brief State data, i.e. the potential. The centering of this is the same as m_dataLocation. 
  */
  MFAMRCellData m_potential;

  /*!
    @brief Electric field. The centering of this is the same as m_dataLocation. 
  */
  MFAMRCellData m_electricField;

  /*!
    @brief Storage for space charge density.
  */
  MFAMRCellData m_rho;

  /*!
    @brief Storage for surface charge density. 
  */
  EBAMRIVData m_sigma;

  /*!
    @brief Residue, e.g. used after solving Poisson equation. 
  */
  MFAMRCellData m_residue;

  /*!
    @brief Cell permittivity
  */
  MFAMRCellData m_permittivityCell;

  /*!
    @brief Face permittivity. 
  */
  MFAMRFluxData m_permittivityFace;

  /*!
    @brief EB permittivity.
  */
  MFAMRIVData m_permittivityEB;

  /*!
    @brief Flag for checking if voltage has been set. 
  */
  bool m_isVoltageSet;

  /*!
    @brief If true, potential will be added to plot files. 
  */
  bool m_plotPotential;

  /*!
    @brief If true, space charge will be added to plot files. 
  */
  bool m_plotRho;

  /*!
    @brief If true, the electric field will be added to plot files.
  */
  bool m_plotElectricField;

  /*!
    @brief If true, the electric field on the inside of dielectrics will be added to plot files. 
  */
  bool m_plotElectricFieldSolid;

  /*!
    @brief If true, the residue will be added to plot files. 
  */
  bool m_plotResidue;

  /*!
    @brief If true, the permittivity will be added to plot files. 
  */
  bool m_plotPermittivity;

  /*!
    @brief If true, m_sigma will be added to plot files. 
  */
  bool m_plotSigma;

  /*!
    @brief Use slopes when regridding or ont
  */
  bool m_regridSlopes;

  /*!
    @brief Verbosity for this calss. 
  */
  int m_verbosity;

  /*!
    @brief Current time step
  */
  int m_timeStep;

  /*!
    @brief Last time step increment
  */
  Real m_dt;

  /*!
    @brief Current time
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
    @brief Domain BC functions. This is used to map space/time to a voltage/field at the domain faces. 
  */
  std::map<ElectrostaticDomainBc::DomainSide, ElectrostaticDomainBc::BcFunction> m_domainBcFunctions;

  /*!
    @brief BC functions (Dirichlet) on the electrodes. Used to map space/time to a voltage on the electrodes. 
  */
  std::vector<std::pair<Electrode, ElectrostaticEbBc::BcFunction>> m_electrodeBcFunctions;

  /*!
    @brief Voltage function. 
  */
  std::function<Real(const Real a_time)> m_voltage;

  /*!
    @brief Parse solver class verbosity
  */
  virtual void
  parseVerbosity();

  /*!
    @brief Function which parses which plot variables to write to plot files
  */
  virtual void
  parsePlotVariables();

  /*!
    @brief Parse domain boundary conditions. 
  */
  virtual void
  parseDomainBc();

  /*!
    @brief Parse slope regrid
  */
  virtual void
  parseRegridSlopes();

  /*!
    @brief Set default BC functions. This sets all the m_domainBcFunction objects to s_defaultDomainBcFunction, which return 1 everywhere. 
  */
  virtual void
  setDefaultDomainBcFunctions();

  /*!
    @brief Set default Dirichlet boundary conditions on the embedded boundaries. 
    @details For each electrode the default boundary condition is set to m_voltage*electrode.getFraction(). 
  */
  virtual void
  setDefaultEbBcFunctions();

  /*!
    @brief Get relative permittivity at some point in space
    @param[in] a_position Physical position
    @note This routine is used for computing the permittivity inside dielectrics. 
  */
  inline Real
  getDielectricPermittivity(const RealVect& a_pos) const;

  /*!
    @brief Set cell-centered permittivities
    @param[out] a_perm    Permittivity (on either cell center or centroid)
    @param[in]  a_cellBox Computation box
    @param[in]  a_ebisbox EBIS box
    @param[in]  a_probLo  Lower-left corner off comptuational domain
    @param[in]  a_dx      Resolution
  */
  virtual void
  setCellPermittivities(EBCellFAB&      a_perm,
                        const Box&      a_cellBox,
                        const EBISBox&  a_ebisbox,
                        const RealVect& a_probLo,
                        const Real&     a_dx);

  /*!
    @brief Set face-centered permittivities
    @param[out] a_perm Permittivity (on face centroid or center)
    @param[in]  a_cellBox Computation box
    @param[in]  a_ebisbox EBIS box
    @param[in]  a_probLo  Lower-left corner off comptuational domain
    @param[in]  a_dx       Resolution
  */
  virtual void
  setFacePermittivities(EBFluxFAB&      a_perm,
                        const Box&      a_cellBox,
                        const EBISBox&  a_ebisbox,
                        const RealVect& a_probLo,
                        const Real&     a_dx);

  /*!
    @brief Set EB-centered permittivities
    @param[out] a_perm Permittivity (on face centroid or center)
    @param[in]  a_cellBox Computation box
    @param[in]  a_ebisbox EBIS box
    @param[in]  a_probLo  Lower-left corner off comptuational domain
    @param[in]  a_dx       Resolution
  */
  virtual void
  setEbPermittivities(BaseIVFAB<Real>& a_perm,
                      const Box&       a_cellBox,
                      const EBISBox&   a_ebisbox,
                      const RealVect&  a_origin,
                      const Real&      a_dx);

  /*!
    @brief Write surface data to volume data holder. 
    @details This takes the valid data from a_data and puts in onto a_output. 
    @param[inout] a_output Output data holder
    @param[inout] a_comp   Current component in a_output. 
    @param[in]    a_data   Surface data. 
    @param[in]    a_interp Special flag which, if true, tells AmrMesh to interpolate the data to centroids (in irregular cells). 
  */
  virtual void
  writeSurfaceData(EBAMRCellData& a_output, int& a_comp, const EBAMRIVData& a_data);

  /*!
    @brief Returns BC type based on string.
    @param[in] a_str Boundary condition string. Must be "dirichlet number", "neumann number", "dirichlet_custom" or "neumann_custom"
    @return This returns the boundary condition type, either ElectrostaticDomainBc::Dirichlet or ElectrostaticDomainBc::Neumann
  */
  virtual ElectrostaticDomainBc::BcType
  parseBcString(const std::string a_str) const;

  /*!
    @brief Shortcut for making a boundary condition string. 
    @param[in] a_dir  Direction.
    @param[in] a_side Coordinate side. 
    @return Returns string of type m_className.bc.direction.side. 
  */
  virtual std::string
  makeBcString(const int a_dir, const Side::LoHiSide a_side) const;
};

#include <CD_NamespaceFooter.H>

#include <CD_FieldSolverImplem.H>

#endif
