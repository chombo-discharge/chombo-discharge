/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrGodunov.H
  @brief  Declaration of a class which implements CdrMultigrid using Godunov-type advection. 
  @author Robert Marskar
*/

#ifndef CD_CdrGodunov_H
#define CD_CdrGodunov_H

// Chombo includes
#include <EBAdvectLevelIntegrator.H>

// Our includes
#include <CD_CdrMultigrid.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Godunov implementation for advection. 
  @details The CdrGodunov class perform Godunov-type advection with time-extrapolated slope-limiting to faces. Note that this class implements CdrMultigrid which
  is the base class that implements second order diffusion (either explicit or implicit). The rationale for doing time-extrapolation ala Godunov is that one can use
  centered differencing phi^(k+1) = phi^k - dt*div[v^(k+1/2) * phi^(k+1/2)] for a second order accurate discretization in time. This class uses Chombo functionality
  from EBAdvectLevelIntegrator which uses the Graves/Trebotich discretization in "An adaptive finite volume method for the incompressible Navier–Stokes equations 
  in complex geometries" (DOI: 10.2140/camcos.2015.10.43). Note that if source terms are involved the extrapolation is more involved. 
*/
class CdrGodunov : public CdrMultigrid
{
public:
  /*!
    @brief Constructor
  */
  CdrGodunov();

  /*!
    @brief Destructor (does nothing)
  */
  virtual ~CdrGodunov();

  /*!
    @brief Parse class options to put object in usable state. 
  */
  virtual void
  parseOptions() override;

  /*!
    @brief Parse runtime options
  */
  virtual void
  parseRuntimeOptions() override;

  /*!
    @brief Allocate internal data holders in class object. 
  */
  virtual void
  allocate() override;

  /*!
    @brief Compute the largest possible advective time step (for explicit methods)
    @details This computes dt = dx/max(|vx|,|vy|,|vz|), minimized over all grid levels and patches. 
    @note This is the appropriate time step routine for the BCG reconstruction. 
  */
  virtual Real
  computeAdvectionDt() override;

protected:
  /*!
    @brief Advection object -- this implements the Graves/Trebotich discretization. 
  */
  Vector<RefCountedPtr<EBAdvectLevelIntegrator>> m_levelAdvect;

  /*!
    @brief If true, slopes are limited (always use limiting). 
  */
  bool m_limitSlopes;

  /*!
    @brief Turn on/off source terms when time-extrapolating
  */
  bool m_extrapolateSourceTerm;

  /*!
    @brief Parses slope limiter options
  */
  virtual void
  parseSlopeLimiter();

  /*!
    ®brief Parse whether or not to include a source-term in the Godunov extrapolation.
  */
  virtual void
  parseExtrapolateSourceTerm();

  /*!
    @brief Godunov face extrapolation method for advection
    @param[out] a_facePhi  Phi on face centers
    @param[in]  a_cellPhi  Phi on cell centers
    @param[in]  a_extrapDt Time centering (i.e. extrapolation) of the face-centered states. 
    @details This functions calls Chombo code which implements the Graves/Trebotich discretization. 
  */
  virtual void
  advectToFaces(EBAMRFluxData& a_facePhi, const EBAMRCellData& a_phi, const Real a_extrapDt) override;
};

#include <CD_NamespaceFooter.H>

#endif
