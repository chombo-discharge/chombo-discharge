/*!
  @file   cdr_sg.H
  @brief  Declaration of a Scharfetter-Gummel implementation of cdr_solver
  @author Robert Marskar
  @date   Jan. 2018
*/

#include "cdr_solver.H"

#ifndef _CDR_SG_
#define _CDR_SG_

class cdr_sg : public cdr_solver {
public:

  /*!
    @brief Constructor
  */
  cdr_sg();

  /*!
    @brief Destructor
  */
  virtual ~cdr_sg();

  /*!
    @brief Required number of ghost cells
  */
  virtual int query_ghost() const;

  /*!
    @brief Pure advective advance. Calling this is an error - SG currently only works with MOL methods.
  */
  virtual void advance_advect(EBAMRCellData& a_state, const Real a_dt);

  /*!
    @brief Pure diffusion advance. Calling this is an error - SG currently only works with MOL methods. 
  */
  virtual void advance_diffusion(EBAMRCellData& a_state, EBAMRCellData& a_error, const Real a_dt);

  /*!
    @brief Compute advective derivative. Calling this is an error - SG currently only works with MOL methods. 
  */
  virtual void compute_divF(EBAMRCellData& a_divF, const EBAMRCellData& a_state, const Real a_extrap_dt, const bool a_redist);

  /*!
    @brief Forward Euler advance with subcycling
  */
  virtual void eulerF_subcycle(EBAMRCellData& a_state, const Real a_dt, const bool a_redist);
    
  /*!
    @brief Compute diffusion term. Calling this is an error - SG currently only works with MOL methods. 
  */
  virtual void compute_divD(EBAMRCellData& a_diffusive_term, const EBAMRCellData& a_state);

protected:
  
  /*!
    @brief Compute div(nv - D*grad(n))
  */
  virtual void compute_divJ(EBAMRCellData& a_divJ, const EBAMRCellData& a_state, const Real a_extrap_dt);

  /*!
    @brief Compute the Scharfetter-Gummel flux
  */
  virtual void compute_sg_flux(EBAMRFluxData&       a_flux,
			       const EBAMRCellData& a_state,
			       const EBAMRFluxData& a_velo,
			       const EBAMRFluxData& a_diffco);

  /*!
    @brief Compute the Scharfetter-Gummel flux
  */
  virtual void compute_sg_flux(LevelData<EBFluxFAB>&       a_flux,
			       const LevelData<EBCellFAB>& a_state,
			       const LevelData<EBFluxFAB>& a_velo,
			       const LevelData<EBFluxFAB>& a_diffco,
			       const int                   a_lvl);

  /*!
    @brief Set boundary outflow flux
  */
  virtual void compute_bndry_outflow(LevelData<EBFluxFAB>&       a_flux,
				     const int                   a_lvl);
};

#endif
