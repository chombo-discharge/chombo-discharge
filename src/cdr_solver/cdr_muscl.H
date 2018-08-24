/*!
  @file   cdr_muscl.H
  @brief  Declaration of a muscl class for the convection diffusion reaction equations
  @author Robert Marskar
  @date   Feb. 2018
*/

#ifndef _CDR_MUSCL_
#define _CDR_MUSCL_

#include "cdr_tga.H"

/*!
  @brief Class that uses MUSCL discretization for the advective term in the cdr equation

  @details 

  Class options
  -------------

  The user may modify the cell slopes through

      cdr_muscl.slope_limiter = superbee     # Valid options are 'koren', 'minmod', 'superbee', 'van_leer'
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
  virtual int query_ghost() const;

  /*!
    @brief Pure advective advance. Calling this is an error - SG currently only works with MOL methods.
  */
  virtual void advance_advect(EBAMRCellData& a_state, const Real a_dt);

  /*!
    @brief MUSCL advection to faces
  */
  virtual void advect_to_faces(EBAMRFluxData& a_face_state, const EBAMRCellData& a_state, const Real a_extrap_dt);

  /*!
    @brief Allocate internal stuff
  */
  virtual void allocate_internals();

protected:

  /*!
    @brief
  */
  bool m_use_slope_limiter;

  /*!
    @brief Compute slopes
  */
  virtual void compute_slopes(EBCellFAB& a_deltaC, const EBCellFAB& a_state, const Box& a_box, const ProblemDomain& a_domain);

  /*!
    @brief Compute outflow fluxes at the boundary
  */
  virtual void compute_bndry_outflow(LevelData<EBFluxFAB>& a_flux, const int a_lvl);

  /*!
    @brief Riemann solve (just choose the upwind side)
  */
  virtual void upwind(EBFluxFAB&           a_face_states,
		      const EBCellFAB&     a_slopes,
		      const EBCellFAB&     a_state,
		      const EBFluxFAB&     a_velo,
		      const ProblemDomain& a_domain,
		      const Box&           a_box);

};
#endif
