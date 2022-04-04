/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_TracerParticleSolver.H
  @brief  Declaration of a solver class that advances tracer particles
  @author Robert Marskar
*/

#ifndef CD_TracerParticleSolver_H
#define CD_TracerParticleSolver_H

// Our includes
#include <CD_AmrMesh.H>
#include <CD_ComputationalGeometry.H>
#include <CD_TracerParticle.H>

/*!
  @brief Base class for a tracer particle solver. This solver can advance particles in a pre-defined velocity field.
  @details The user can set the velocity field through public member functions. This class is templated because the
  tracer particles can contain an arbitrary number of traced fields. 

  This is also a single-phase solver -- i.e. the particles only live on one of the phases. Extensions to multiphase
  is certainly possible. 
*/
template <size_t M, size_t N> 
class TracerParticleSolver {
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
  TracerParticleSolver(const RefCountedPtr<AmrMesh>& a_amr,
		       const RefCountedPtr<ComputationalGeometry> m_compGeom);

  /*!
    @brief Disallowed assignment operator. 
    @param[in] a_other Other solver. 
  */
  TracerParticleSolver& operator=(const TracerParticleSolver& a_other) = delete;

  /*!
    @brief Disallowed move assignment operator. 
    @param[in] a_other Other solver. 
  */
  TracerParticleSolver& operator=(const TracerParticleSolver&& a_other) = delete;

  /*!
    @brief Destructor
  */
  virtual ~TracerParticleSolver();

  /*!
    @brief Allocate internal storage
    @note This allocates a bunch of storage -- if the solver is not diffusive or mobile then we do our best to trim memory. 
  */
  virtual void allocateInternals();



  /*!
    @brief Set the tracer particle velocity field
    @param[in] a_velocityField Velocity field.
  */
  virtual void setVelocityField(const EBAMRCellData& a_velocityField);

  /*!
    @brief Set the tracer particle velocity field
    @param[in] a_velocityField Velocity field.
  */  
  virtual void setVelocityField(const std::function<RealVect(const RealVect& a_position)>& a_velocityField);

  /*!
    @brief Perform pre-regrid operations. 
    @param[in] a_lbase          Coarsest level that changed during regrid. 
    @param[in] a_oldFinestLevel Finest grid level before the regrid operation. 
    @note This stores m_phi and m_source. 
  */
  virtual void preRegrid(const int a_lbase, const int a_oldFinestLevel);

  /*!
    @brief Regrid this solver. 
    @param[in] a_lmin           Coarsest level where grids did not change. 
    @param[in] a_oldFinestLevel Finest AMR level before the regrid. 
    @param[in] a_newFinestLevel Finest AMR level after the regrid. 
    @note The velocity is NOT remeshed -- the user is responsible for setting this after the regrid. 
  */
  virtual void regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel);

  /*!
    @brief Write plot file
    @details The name of the plot file is m_name.stepXXXXX.DIM.hdf5
    @note This calls writePlotData(...)
  */
  virtual void writePlotFile();

  /*!
    @brief Write output data to a_output
    @param[inout] a_output Output data holder.
    @param[inout] a_icomp  Starting component where this solver begins writing the output. 
    @note This will write the plot data in this solver to a_output, starting on a_comp
    @note This routine writes m_phi on centroids (and not cell-centers). 
  */
  virtual void writePlotData(EBAMRCellData& a_output, int& a_icomp);

  /*!
    @brief Write checkpoint data into HDF5 file. 
    @paramo[out] a_handle HDF5 file. 
    @param[in]   a_level Grid level
  */
#ifdef CH_USE_HDF5
  virtual void writeCheckpointLevel(HDF5Handle& a_handle, const int a_level) const;
#endif

  /*!
    @brief Read checkpoint data from HDF5 file. 
    @param[in] a_handle HDF5 handle.
    @param[in] const int a_level Grid level
  */
#ifdef CH_USE_HDF5  
  virtual void readCheckpointLevel(HDF5Handle& a_handle, const int a_level);
#endif

  /*!
    @brief Get particles
    @return m_particles. 
  */
  virtual ParticleContainer<TracerParticle<M, N> >& getParticles();

  /*!
    @brief Get particles. Const version. 
    @return m_particles. 
  */
  virtual const ParticleContainer<TracerParticle<M, N> >& getParticles() const;  
  
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
    @brief Verbosity level
  */
  int m_verbosity;

  /*!
    @brief Mesh-based velocity field
  */
  EBAMRCellData m_velocityField;

  /*!
    @brief Particle data.
  */
  ParticleContainer<TracerParticle<M, N> > m_particles;
};

#endif
