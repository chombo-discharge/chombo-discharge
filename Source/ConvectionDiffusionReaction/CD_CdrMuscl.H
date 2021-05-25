/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrMuscl.H
  @brief  Declaration of a muscl class for the convection diffusion reaction equations
  @author Robert Marskar
*/

#ifndef CD_CdrMuscl_H
#define CD_CdrMuscl_H

// Our includes
#include <CD_CdrTGA.H>
#include <CD_CdrGodunov.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class that uses MUSCL discretization for the advective term in the cdr equation
*/
class CdrMuscl : public CdrTGA {
public:

  /*!
    @brief Constructor
  */
  CdrMuscl();

  /*!
    @brief Destructor
  */
  virtual ~CdrMuscl();

  /*!
    @brief Query number of ghost cells
  */
  virtual int queryGhost() const override;

  /*!
    @brief Parse options
  */
  virtual void parseOptions() override;

  /*!
    @brief Parse options
  */
  virtual void parseRuntimeOptions() override;

  /*!
    @brief Pure advective advance. Calling this is an error - SG currently only works with MOL methods.
  */
  virtual void advance_advect(EBAMRCellData& a_phi, const Real a_dt);

  /*!
    @brief MUSCL advection to faces
  */
  virtual void advectToFaces(EBAMRFluxData& a_facePhi, const EBAMRCellData& a_phi, const Real a_extrapDt) override;

  /*!
    @brief Allocate internal stuff
  */
  virtual void allocateInternals() override;

protected:

  /*!
    @brief
  */
  bool m_slopelim;

  /*!
    @brief Compute slopes
  */
  virtual void computeSlopes(EBCellFAB&           a_deltaC,
			      const EBCellFAB&     a_phi,
			      const Box&           a_box,
			      const ProblemDomain& a_domain,
			      const int            a_level,
			      const DataIndex&     a_dit);

  /*!
    @brief Compute outflow fluxes at the boundary
  */
  virtual void computeBoundaryOutflow(LevelData<EBFluxFAB>& a_flux, const int a_lvl);

  /*!
    @brief Riemann solve (just choose the upwind side)
  */
  virtual void upwind(EBFluxFAB&           a_facePhi,
		      const EBCellFAB&     a_slopes,
		      const EBCellFAB&     a_phi,
		      const EBFluxFAB&     a_velo,
		      const ProblemDomain& a_domain,
		      const Box&           a_box,
		      const int&           a_level,
		      const DataIndex&     a_dit);

  virtual void parseSlopeLimiter();     // Parses slope limiter settings
};

#include <CD_NamespaceFooter.H>

#endif
