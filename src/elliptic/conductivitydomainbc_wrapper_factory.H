/*!
  @file   conductivitydomainbc_wrapper_factory.H
  @brief  Declaration of a factory class for conductivitydomainbc_wrapper
  @author Robert Marskar
  @date   Jan. 2018
*/

#ifndef _CONDUCTIVITYDOMAINBC_WRAPPER_FACTORY_
#define _CONDUCTIVITYDOMAINBC_WRAPPER_FACTORY_

#include "conductivitydomainbc_wrapper.H"
#include "robin_coef.H"

#include <ConductivityBaseDomainBC.H>

namespace ChomboDischarge {

  /*!
    @brief Class for generating conductivitydomainbc_wrapper objects
    @details This class does not (currently) do Robin boundary conditions, although it can probably be extended to take care of 
    those cases as well. 
  */
  class conductivitydomainbc_wrapper_factory : public BaseDomainBCFactory {
  public:
  
    /*!
      @brief Constructor
    */
    conductivitydomainbc_wrapper_factory();

    /*!
      @brief Destructor
    */
    virtual ~conductivitydomainbc_wrapper_factory();

    virtual void set_wallbc(const Vector<RefCountedPtr<wall_bc> >& a_wallbc);
    virtual void set_potentials(const Vector<RefCountedPtr<BaseBCFuncEval> >& a_potentials);
    virtual void set_robin_coefs(const Vector<RefCountedPtr<robin_coef> >& a_robinco);

    /*!
      @brief Factory method
    */
    virtual conductivitydomainbc_wrapper* create(const ProblemDomain& a_domain,
						 const EBISLayout&    a_ebisl,
						 const RealVect&      a_dx);
  
  protected:

    bool m_hasbc;

    Vector<RefCountedPtr<wall_bc> >        m_wallbc;      // BC identifier
    Vector<RefCountedPtr<BaseBCFuncEval> > m_potentials;  // Wall potential
    Vector<RefCountedPtr<robin_coef> >     m_robinco;     // Robin coefficients for Robin

  };
}
#endif
