/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrMuscl.H
  @brief  Declaration of a class which implements CdrTGA using MUSCL for advection. 
  @author Robert Marskar
*/

#ifndef CD_CdrMuscl_H
#define CD_CdrMuscl_H

// Our includes
#include <CD_CdrTGA.H>
#include <CD_CdrGodunov.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class that uses MUSCL discretization advection.
  @details This class performs MUSCL-type 
*/
class CdrMuscl : public CdrTGA {
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
  virtual void parseOptions() override;

  /*!
    @brief Parse runtime options
  */
  virtual void parseRuntimeOptions() override;

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
    @brief If true, slopes are limited (always use limiting). 
  */
  bool m_limitSlopes;

  /*!
    @brief Compute slopes. 
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
