/*!
  @file robinconductivityebbc.H
  @brief Declaration of a factory class for Robin type EB boundary condition
  @author Robert Marskar
  @date Jan. 2018
*/

#ifndef _ROBINCONDUCTIVITYEBBCFACTORY_
#define _ROBINCONDUCTIVITYEBBCFACTORY_

#include "EBAMRData.H"
#include "irreg_stencil.H"
#include "robinconductivityebbc.H"

#include <BaseEBBC.H>

namespace ChomboDischarge {

  /*!
    @brief Factory class for generating Robin boundary conditions \f$a\phi + b\partial_n\phi = g\f$, where \f$n\f$ is the EB unit normal. 
  */
  class robinconductivityebbcfactory : public BaseEBBCFactory {
  public:

    /*!
      @brief Constructor
    */
    robinconductivityebbcfactory(const RealVect a_origin);

    /*!
      @brief Destructor
    */
    ~robinconductivityebbcfactory();

    /*!
      @brief Set bc order
    */
    virtual void set_type(const stencil_type a_type);

    /*!
      @brief Function based coefficients and right hand sides
    */
    virtual void set_coefs(RefCountedPtr<robin_coef> a_robinco);

    /*!
      @brief Set constant coefficients
    */
    virtual void set_coefs(const Real a_aco, const Real b_bco, const Real a_rhs);

    /*!
      @brief Set data-based coefficients
    */
    virtual void set_coefs(const EBAMRIVData& a_aco,
			   const EBAMRIVData& a_bco,
			   const EBAMRIVData& a_rhs);

    /*!
      @brief Create Robin type boundary condition
    */
    robinconductivityebbc* create(const ProblemDomain& a_domain,
				  const EBISLayout&    a_ebisl,
				  const RealVect&      a_dx,
				  const IntVect*       a_ghost_phi = 0,
				  const IntVect*       a_ghost_rhs = 0);

  protected:

    stencil_type m_type;

    RealVect m_origin;

    RefCountedPtr<robin_coef> m_robinco;

    Real m_aco;

    Real m_bco;
  
    Real m_rhs;

    EBAMRIVData m_acodata;

    EBAMRIVData m_bcodata;
  
    EBAMRIVData m_rhsdata;

    bool m_const_coeff;

    bool m_data_coeff;

    bool m_func_coeff;
  };
}

#endif
