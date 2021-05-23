/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrGodunov.H
  @brief  Godunov advection interface to CdrSolver
  @author Robert Marskar
*/

#ifndef CD_CdrGodunov_H
#define CD_CdrGodunov_H

// Chombo includes
#include <EBAdvectLevelIntegrator.H>

// Our includes
#include <CD_CdrTGA.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Godunov implementation for advection. 
*/
class CdrGodunov : public CdrTGA {
public:

  /*!
    @brief Constructor
  */
  CdrGodunov();

  /*!
    @brief Destructor
  */
  virtual ~CdrGodunov();

  /*!
    @brief Parse options
  */
  virtual void parseOptions() override;

  /*!
    @brief Parse options
  */
  virtual void parseRuntimeOptions() override;

  /*!
    @brief Query number of ghost cells
  */
  virtual int queryGhost() const override;

  /*!
    @brief Allocate internal stuff
  */
  virtual void allocateInternals() override;

  /*!
    @brief Average velocities to faces
  */
  virtual void averageVelocityToFaces() override;

  /*!
    @brief Call parent function: Then fix no influx
  */
  virtual void averageVelocityToFaces(EBAMRFluxData& a_faceVelocity, const EBAMRCellData& a_cellVelocity) override;

protected:

  /*!
    @brief Godunov face extrapolation method for advection
  */
  virtual void advectToFaces(EBAMRFluxData& a_facePhi, const EBAMRCellData& a_phi, const Real a_extrapDt) override;

  /*!
    @brief Level advect integrator
  */
  Vector<RefCountedPtr<EBAdvectLevelIntegrator> > m_level_advect;

  /*!
    @brief Turn on/off slope limiting
  */
  bool m_slopelim;

  /*!
    @brief Parses slope limiter options
  */
  virtual void parseSlopeLimiter();     // Parses slope limiter settings
};

#include <CD_NamespaceFooter.H>

#endif
