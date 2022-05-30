/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrCTU.H
  @brief  Declaration of a class which implements CdrMultigrid using MUSCL for advection. 
  @author Robert Marskar
*/

#ifndef CD_CdrCTU_H
#define CD_CdrCTU_H

// Our includes
#include <CD_CdrMultigrid.H>
#include <CD_CdrGodunov.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class that uses a slope limited method for advection, in combination with corner transport upwind (CTU). 
*/
class CdrCTU : public CdrMultigrid
{
public:
  /*!
    @brief Default constructor. 
  */
  CdrCTU();

  /*!
    @brief Destructor (does nothing)
  */
  virtual ~CdrCTU();

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
    @brief MUSCL advection to faces
    @param[out] a_facePhi  Phi on face centers
    @param[in]  a_cellPhi  Phi on cell centers
    @param[in]  a_dt       Time step (i.e. extrapolation) of the face-centered states. 
  */
  virtual void
  advectToFaces(EBAMRFluxData& a_facePhi, const EBAMRCellData& a_cellPhi, const Real a_dt) override;

  /*!
    @brief Compute the largest possible advective time step (for explicit methods)
    @details This computes dt = dx/max(|vx|,|vy|,|vz|), minimized over all grid levels and patches. 
  */
  virtual Real
  computeAdvectionDt() override;

protected:
  /*!
    @brief Supported limiters
  */
  enum class Limiter
  {
    None,
    MinMod,
    Superbee,
    MonotonizedCentral,
  };

  /*!
    @brief If true, slopes are limited in the normal dircetion. 
  */
  Limiter m_limiter;

  /*!
    @brief If true, the discretization uses CTU. Otherwise it's DCU (donor cell upwind)
  */
  bool m_useCTU;

  /*!
    @brief Compute slopes. 
    @param[out] a_normalSlopes Cell slopes. Must have SpaceDim components (for slopes in each direction). 
    @param[in]  a_cellPhi      Cell-centered states
    @param[in]  a_cellBox      Cell-centered box identifying the cells where we compute the slopes. 
    @param[in]  a_domain       Domain on grid level
    @param[in]  a_level        AMR level
    @param[in]  a_dit          Grid index
  */
  virtual void
  computeNormalSlopes(EBCellFAB&           a_normalSlopes,
                      const EBCellFAB&     a_cellPhi,
                      const Box&           a_cellBox,
                      const ProblemDomain& a_domain,
                      const int            a_level,
                      const DataIndex&     a_dit);

  /*!
    @brief Upwind/Riemann solve. This extrapolates the cell-centered data to the inside face centers and selects the upwind side.
    @details This uses the Colella CTU method where we extrapolate with the corner cells as well. a_normalSlopes holds
    the slopes normal to the face direction, and the transverse slopes are added inside this method. This routine also
    solves the corresponding Riemann problem (which is just the upwind face). 
    @param[out] a_facePhi      Face-centered states.
    @param[in]  a_normalSlopes Normal cell slopes
    @param[in]  a_cellPhi      Cell-centered states
    @param[in]  a_cellVel      Cell-centered velocities.
    @param[in]  a_faceVel      Face-centered velocities.
    @param[in]  a_domain       Domain on grid level.
    @param[in]  a_cellBox      Grid box/patch
    @param[in]  a_level        AMR level
    @param[in]  a_dit          Grid index. 
    @param[in]  a_dt           Time step. 
  */
  virtual void
  upwind(EBFluxFAB&           a_facePhi,
         const EBCellFAB&     a_normalSlopes,
         const EBCellFAB&     a_cellPhi,
         const EBCellFAB&     a_cellVel,
         const EBFluxFAB&     a_faceVel,
         const ProblemDomain& a_domain,
         const Box&           a_cellBox,
         const int&           a_level,
         const DataIndex&     a_dit,
         const Real&          a_dt);

  /*!
    @brief Parse slope limiting on/off
  */
  virtual void
  parseSlopeLimiter();

  /*!
    @brief minmod slope function. 
    @param[in] dwl Left slope
    @param[in] dwr Right slope
  */
  Real
  minmod(const Real& dwl, const Real& dwr) const noexcept;

  /*!
    @brief Superbee slope limiter. 
    @param[in] dwl Left slope
    @param[in] dwr Right slope
  */
  Real
  superbee(const Real& dwl, const Real& dwr) const noexcept;

  /*!
    @brief Monotonized central difference slope limiter.
    @param[in] dwl Left slope
    @param[in] dwr Right slope
  */
  Real
  monotonizedCentral(const Real& dwl, const Real& dwr) const noexcept;
};

#include <CD_NamespaceFooter.H>

#endif
