/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrMuscl.H
  @brief  Declaration of a class which implements CdrMultigrid using MUSCL for advection. 
  @author Robert Marskar
*/

#ifndef CD_CdrMuscl_H
#define CD_CdrMuscl_H

// Our includes
#include <CD_CdrMultigrid.H>
#include <CD_CdrGodunov.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class that uses MUSCL discretization advection.
  @details This class performs MUSCL-type advection using a van Leer limiter for obtaining the face-centered states. 
*/
class CdrMuscl : public CdrMultigrid
{
public:
  /*!
    @brief Default constructor. 
  */
  CdrMuscl();

  /*!
    @brief Destructor (does nothing)
  */
  virtual ~CdrMuscl();

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
    @param[in]  a_extrapDt Time centering (i.e. extrapolation) of the face-centered states. 
    @details This compute the face-centered states using van Leer limiting. It calls computeSlopes(...) and upwind(...). 
  */
  virtual void
  advectToFaces(EBAMRFluxData& a_facePhi, const EBAMRCellData& a_cellPhi, const Real a_extrapDt) override;

protected:
  /*!
    @brief If true, slopes are limited (always use limiting). 
  */
  bool m_limitSlopes;

  /*!
    @brief Compute slopes. 
    @param[out] a_deltaC  Cell slopes. Must have SpaceDim components (for slopes in each direction). 
    @param[in]  a_cellPhi Cell-centered states
    @param[in]  a_cellBox Cell-centered box identifying the cells where we compute the slopes. 
    @param[in]  a_domain  Domain on grid level
    @param[in]  a_level   AMR level
    @param[in]  a_dit     Grid index
  */
  virtual void
  computeSlopes(EBCellFAB&           a_deltaC,
                const EBCellFAB&     a_cellPhi,
                const Box&           a_cellBox,
                const ProblemDomain& a_domain,
                const int            a_level,
                const DataIndex&     a_dit);

  /*!
    @brief Upwind/Riemann solve. This extrapolates the cell-centered data to the (inside) face centers and chooses the upwind side. 
    @param[out] a_facePhi Face-centered states.
    @param[in]  a_slopes  Cell slopes
    @param[in]  a_cellPhi Cell-centered states
    @param[in]  a_faceVel Face-centered velocities.
    @param[in]  a_domain  Domain on grid level.
    @param[in]  a_cellBox Grid box/patch
    @param[in]  a_level   AMR level
    @param[in]  a_dit     Grid index. 
  */
  virtual void
  upwind(EBFluxFAB&           a_facePhi,
         const EBCellFAB&     a_slopes,
         const EBCellFAB&     a_cellPhi,
         const EBFluxFAB&     a_faceVel,
         const ProblemDomain& a_domain,
         const Box&           a_cellBox,
         const int&           a_level,
         const DataIndex&     a_dit);

  /*!
    @brief Parse slope limiting on/off
  */
  virtual void
  parseSlopeLimiter();
};

#include <CD_NamespaceFooter.H>

#endif
