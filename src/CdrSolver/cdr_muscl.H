/*!
  @file   cdr_muscl.H
  @brief  Declaration of a muscl class for the convection diffusion reaction equations
  @author Robert Marskar
  @date   Feb. 2018
*/

#ifndef _CDR_MUSCL_
#define _CDR_MUSCL_

#include "cdr_tga.H"
#include "cdr_gdnv.H"

#include "CD_NamespaceHeader.H"

/*!
  @brief Class that uses MUSCL discretization for the advective term in the cdr equation
*/
class cdr_muscl : public cdr_tga {
public:

  /*!
    @brief Constructor
  */
  cdr_muscl();

  /*!
    @brief Destructor
  */
  virtual ~cdr_muscl();

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
  virtual void advect_to_faces(EBAMRFluxData& a_facePhi, const EBAMRCellData& a_phi, const Real a_extrapDt) override;

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
  virtual void compute_slopes(EBCellFAB&           a_deltaC,
			      const EBCellFAB&     a_phi,
			      const Box&           a_box,
			      const ProblemDomain& a_domain,
			      const int            a_level,
			      const DataIndex&     a_dit);

  /*!
    @brief Compute outflow fluxes at the boundary
  */
  virtual void compute_bndry_outflow(LevelData<EBFluxFAB>& a_flux, const int a_lvl);

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

  virtual void parse_slopelim();     // Parses slope limiter settings
};
#include "CD_NamespaceFooter.H"
#endif
