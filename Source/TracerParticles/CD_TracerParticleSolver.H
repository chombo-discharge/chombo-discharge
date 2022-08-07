/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_TracerParticleSolver.H
  @brief  Declaration of a solver class that advances tracer particles.
  @author Robert Marskar
*/

#ifndef CD_TracerParticleSolver_H
#define CD_TracerParticleSolver_H

// Our includes
#include <CD_AmrMesh.H>
#include <CD_ComputationalGeometry.H>
#include <CD_TracerParticle.H>
#include <CD_ParticleContainer.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Base class for a tracer particle solver. This solver can advance particles in a pre-defined velocity field.
  @details The user can set the velocity field through public member functions. This class is templated so the user can 
  switch tracer particle implementations.

  This is a single-phase solver -- i.e. the particles only live on one of the phases. Extensions to multiphase
  is certainly possible. 

  The template requirements on the particle type P are: 
  1) It must contain a function RealVect& position() (derived from BinItem)
  2) It must contain a function const Real& weight() const.
  3) It must contain a function RealVect& velocity().
*/
template <typename P>
class TracerParticleSolver
{
public:
  /*!
    @brief Default constructor
  */
  TracerParticleSolver();

  /*!
    @brief Disallowed constructor. 
    @param[in] a_other Other solver. 
  */
  TracerParticleSolver(const TracerParticleSolver& a_other) = delete;

  /*!
    @brief Disallowed move constructor. 
    @param[in] a_other Other solver. 
  */
  TracerParticleSolver(const TracerParticleSolver&& a_other) = delete;

  /*!
    @brief Full contructor
    @param[in] a_amr      Handle to AmrMesh.
    @param[in] a_compGeom Computational geometry. 
  */
  TracerParticleSolver(const RefCountedPtr<AmrMesh>& a_amr, const RefCountedPtr<ComputationalGeometry> a_compGeom);

  /*!
    @brief Disallowed assignment operator. 
    @param[in] a_other Other solver. 
  */
  TracerParticleSolver&
  operator=(const TracerParticleSolver& a_other) = delete;

  /*!
    @brief Disallowed move assignment operator. 
    @param[in] a_other Other solver. 
  */
  TracerParticleSolver&
  operator=(const TracerParticleSolver&& a_other) = delete;

  /*!
    @brief Destructor
  */
  virtual ~TracerParticleSolver();

  /*!
    @brief Allocate storage for this solver. 
  */
  virtual void
  allocate();

  /*!
    @brief Register operators needed for AMR core functionality. 
  */
  virtual void
  registerOperators() const;

  /*!
    @brief Set AmrMesh.
    @param[in] a_amrMesh AmrMesh
  */
  virtual void
  setAmr(const RefCountedPtr<AmrMesh>& a_amrMesh);

  /*!
    @brief Set the computational geometry.
    @param[in] a_compGeom Computational geometry
  */
  virtual void
  setComputationalGeometry(const RefCountedPtr<ComputationalGeometry>& a_compGeom);

  /*!
    @brief Set the solver realm. 
    @param[in] a_realm Realm name. 
  */
  virtual void
  setRealm(const std::string& a_realm);

  /*!
    @brief Set the solver phase.
    @param[in] a_phase Phase
  */
  virtual void
  setPhase(const phase::which_phase& a_phase);

  /*!
    @brief Set the time for this solver. 
    @param[in] a_step Time step number
    @param[in] a_time Time (in seconds)
    @param[in] a_dt   Time step increment
    @note This sets m_step=a_step, m_time=a_time, m_dt=a_dt
  */
  virtual void
  setTime(const int a_step, const Real a_time, const Real a_dt);

  /*!
    @brief Set the tracer particle velocity field.
    @param[in] a_velocityField Velocity field.
  */
  virtual void
  setVelocity(const EBAMRCellData& a_velocityField);

  /*!
    @brief Perform pre-regrid operations. 
    @param[in] a_lbase          Coarsest level that changed during regrid. 
    @param[in] a_oldFinestLevel Finest grid level before the regrid operation. 
    @note This stores m_phi and m_source. 
  */
  virtual void
  preRegrid(const int a_lbase, const int a_oldFinestLevel);

  /*!
    @brief Regrid this solver. 
    @param[in] a_lmin           Coarsest level where grids did not change. 
    @param[in] a_oldFinestLevel Finest AMR level before the regrid. 
    @param[in] a_newFinestLevel Finest AMR level after the regrid. 
    @note The velocity is NOT remeshed -- the user is responsible for setting this after the regrid. 
  */
  virtual void
  regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel);

  /*!
    @brief Remap particles.
  */
  virtual void
  remap();

  /*!
    @brief Deposit particle weight on mesh. 
    @param[out] a_phi Deposited weight. 
  */
  virtual void
  deposit(EBAMRCellData& a_phi) const;

  /*!
    @brief Interpolate a scalar field onto the particle weight
  */
  virtual void
  interpolateWeight(const EBAMRCellData& a_scalar) noexcept;

  /*!
    @brief Interpolate particles velocities.
  */
  virtual void
  interpolateVelocities();

  /*!
    @brief Parse solver options.
  */
  virtual void
  parseOptions();

  /*!
    @brief Parse solver run-time options.
  */
  virtual void
  parseRuntimeOptions();

  /*!
    @brief Get the number of plot variables
    @return Returns number of plot variables. 
  */
  virtual int
  getNumberOfPlotVariables() const;

  /*!
    @brief Get plot variable names
    @return List of plot variable names. 
  */
  virtual Vector<std::string>
  getPlotVariableNames() const;

  /*!
    @brief Write plot file
    @details The name of the plot file is m_name.stepXXXXX.DIM.hdf5
    @note This calls writePlotData(...)
  */
  virtual void
  writePlotFile();

  /*!
    @brief Write output data to a_output
    @param[inout] a_output Output data holder.
    @param[inout] a_icomp  Starting component where this solver begins writing the output. 
    @note This will write the plot data in this solver to a_output, starting on a_comp
    @note This routine writes m_phi on centroids (and not cell-centers). 
  */
  virtual void
  writePlotData(EBAMRCellData& a_output, int& a_icomp);

  /*!
    @brief Write checkpoint data into HDF5 file. 
    @paramo[out] a_handle HDF5 file. 
    @param[in]   a_level Grid level
  */
#ifdef CH_USE_HDF5
  virtual void
  writeCheckpointLevel(HDF5Handle& a_handle, const int a_level) const;
#endif

  /*!
    @brief Read checkpoint data from HDF5 file. 
    @param[in] a_handle HDF5 handle.
    @param[in] const int a_level Grid level
  */
#ifdef CH_USE_HDF5
  virtual void
  readCheckpointLevel(HDF5Handle& a_handle, const int a_level);
#endif

  /*!
    @brief Compute dt = dx/max(v_x,v_y, v_z) minimized over all particles
  */
  virtual Real
  computeDt() const;

  /*!
    @brief Get all particles.
    @return m_particles
  */
  virtual ParticleContainer<P>&
  getParticles();

  /*!
    @brief Get all particles. Const version. 
    @return m_particles
  */
  virtual const ParticleContainer<P>&
  getParticles() const;

  /*!
    @brief Return the velocity field
    @return *m_velocityField
  */
  const EBAMRCellData&
  getVelocityField() const;

protected:
  /*!
    @brief Handle to AMR mesh. 
  */
  RefCountedPtr<AmrMesh> m_amr;

  /*!
    @brief Handle to computational geometry.
  */
  RefCountedPtr<ComputationalGeometry> m_computationalGeometry;

  /*!
    @brief Deposition type when depositing particles to grid
  */
  DepositionType m_deposition;

  /*!
    @brief Coarse-fine deposition strategy. 
  */
  CoarseFineDeposition m_coarseFineDeposition;

  /*!
    @brief Realm where this solver lives
  */
  std::string m_realm;

  /*!
    @brief Solver name
  */
  std::string m_name;

  /*!
    @brief Class name
  */
  std::string m_className;

  /*!
    @brief Phase where this solver lives
  */
  phase::which_phase m_phase;

  /*!
    @brief Time step
  */
  Real m_dt;

  /*!
    @brief Time
  */
  Real m_time;

  /*!
    @brief Time step
  */
  int m_timeStep;

  /*!
    @brief Verbosity level
  */
  int m_verbosity;

  /*!
    @brief PVR buffer for particle deposition. 
  */
  int m_pvrBuffer;

  /*!
    @brief Halo buffer for particle deposition. 
  */
  int m_haloBuffer;

  /*!
    @brief Hook for adding weight to plot file. 
  */
  bool m_plotWeight;

  /*!
    @brief Hook for adding velocity field to plot file.
  */
  bool m_plotVelocity;

  /*!
    @brief Do volumetric scaling when depositing (or not).
  */
  bool m_volumeScale;

  /*!
    @brief Mesh-based velocity field. Must be set through setVelocity
  */
  const EBAMRCellData* m_velocityField;

  /*!
    @brief Computational particles. 
  */
  ParticleContainer<P> m_particles;

  /*!
    @brief Parse deposition method
  */
  void
  parseDeposition();

  /*!
    @brief Parse plot variables.
  */
  void
  parsePlotVariables();

  /*!
    @brief Parse solver verbosity.
  */
  void
  parseVerbosity();

  /*!
    @brief Generic particle deposition method for putting a scalar field onto the mesh.
    @param[out] a_phi                  Deposited quantity.
    @param[out] a_particles            Particles to be deposited. 
    @param[in]  a_baseDeposition       Base deposition method.
    @param[in]  a_coarseFineDeposition Coarse-fine deposition strategy.
  */
  template <typename T, const Real& (T::*particleScalarFunction)() const>
  void
  depositParticles(EBAMRCellData&             a_phi,
                   const ParticleContainer<T> a_particles,
                   const DepositionType       a_baseDeposition,
                   const CoarseFineDeposition a_coarseFineDeposition) const;
};

#include <CD_NamespaceFooter.H>

#include <CD_TracerParticleSolverImplem.H>

#endif
