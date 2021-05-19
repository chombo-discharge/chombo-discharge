/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_SigmaSolver.H
  @brief  Declaration of a surface charge solver
  @author Robert Marskar
*/

#ifndef CD_SigmaSolver_H
#define CD_SigmaSolver_H

// Our includes
#include <CD_ComputationalGeometry.H>
#include <CD_AmrMesh.H>
#include <CD_NamespaceHeader.H>
  
/*!
  @brief Surface charge solver
*/
class SigmaSolver {
public:

  /*!
    @brief Constructor
  */
  SigmaSolver();

  /*!
    @brief Destructor
  */
  ~SigmaSolver();

  /*!
    @brief Get the Realm
  */
  const std::string getRealm() const;

  /*!
    @brief Set the Realm
  */
  virtual void setRealm(const std::string a_realm); 

  /*!
    @brief Allocate internal storage
  */
  virtual void allocateInternals();

  /*!
    @brief Cache state for regridding
  */
  virtual void preRegrid(const int a_lbase, const int a_oldFinestLevel);

  /*!
    @brief Compute the right-hand side
  */
  virtual void computeRHS(EBAMRIVData& a_rhs);

  /*!
    @brief Deallocate internal storage
  */
  virtual void deallocateInternals();

  /*!
    @brief Regrid method
  */
  virtual void regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel);

  /*!
    @brief Register operators
  */
  virtual void registerOperators();
  
  /*!
    @brief Reset data holder on cells
  */
  virtual void resetCells(EBAMRIVData& a_data);

  /*!
    @brief Set the amr object
  */
  virtual void setAmr(const RefCountedPtr<AmrMesh>& a_amr);

  /*!
    @brief Set computational geometry
  */
  virtual void setComputationalGeometry(const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry);

  /*!
    @brief Set phase
    @details This must be done BEFORE callilng setComputationalGeometry
  */
  virtual void setPhase(phase::which_phase a_phase = phase::gas);

  /*!
    @brief Convenience function. Set the surface charge
  */
  virtual void setSigma(const EBAMRIVData& a_sigma);

  /*!
    @brief Convenience function. Set the surface charge. 
  */
  virtual void setSigma(const Real a_sigma);

  /*!
    @brief Set verbosity.
  */
  virtual void setVerbosity(const int a_verbosity);

  /*!
    @brief Get current time
  */
  virtual void setTime(const int a_step, const Real a_time, const Real a_dt);

  /*!
    @brief Write checkpoint data into handle
  */
  virtual void writeCheckpointLevel(HDF5Handle& a_handle, const int a_level) const;

  /*!
    @brief Read checkpoint data from handle
  */
  virtual void readCheckpointLevel(HDF5Handle& a_handle, const int a_level);

  /*!
    @brief Write output data to a_output
  */
  virtual void writePlotData(EBAMRCellData& a_output, int& a_comp);

  /*!
    @brief Number of variables to be plotted
  */
  virtual int getNumberOfPlotVariables();

  /*!
    @brief Get output plot names
  */
  virtual Vector<std::string> getPlotVariableNames() const;

  /*!
    @brief Compute total charge
  */
  virtual Real computeCharge();

  /*!
    @brief Get internal state
  */
  virtual EBAMRIVData& getPhi();

  /*!
    @brief Get internal state
  */
  virtual EBAMRIVData& getFlux();
  
protected:

  /*!
    @brief Realm
  */
  std::string m_realm;

  /*!
    @brief Index spaces
  */
  RefCountedPtr<mfis> m_multifluidIndexSpace;

  /*!
    @brief Computational geometry
  */
  RefCountedPtr<ComputationalGeometry> m_computationalGeometry;

  /*!
    @brief AMR; needed for grid stuff
  */
  RefCountedPtr<AmrMesh> m_amr;

  /*!
    @brief State
  */
  EBAMRIVData m_phi;

  /*!
    @brief Cached state for regridding. DO NOT TOUCH!
  */
  EBAMRIVData m_cache;

  /*!
    @brief Flux
  */
  EBAMRIVData m_flux;

  /*!
    @brief Phase
  */
  phase::which_phase m_phase;

  /*!
    @brief Solver verbosity
  */
  int m_verbosity;

  /*!
    @brief Time step
  */
  int m_timeStep;

  /*!
    @brief Current time
  */
  Real m_time;

  /*!
    @brief Time increment
  */
  Real m_dt;

  /*!
    @brief Output sigma
  */
  bool m_plotPhi;
  
};

#include <CD_NamespaceFooter.H>

#endif

