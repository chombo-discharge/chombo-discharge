/*!
  @file conductivitydomainbc_wrapper_factory.H
  @brief Declaration of a factory class for conductivitydomainbc_wrapper
  @author Robert Marskar
  @date Jan. 2018
  @todo Segregate implementation file
*/

#ifndef _CONDUCTIVITYDOMAINBC_WRAPPER_FACTORY_
#define _CONDUCTIVITYDOMAINBC_WRAPPER_FACTORY_

#include "conductivitydomainbc_wrapper.H"

#include <ConductivityBaseDomainBC.H>

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
  conductivitydomainbc_wrapper_factory(){
    m_hasbc = false;
  }

  /*!
    @brief Destructor
  */
  virtual ~conductivitydomainbc_wrapper_factory(){}


  /*!
    @brief Set wall bc
  */
  virtual void set_wallbc(const Vector<RefCountedPtr<wall_bc> >& a_wallbc){
    m_wallbc = a_wallbc;

    m_hasbc = true;
  }

  virtual void set_potential(const RefCountedPtr<BaseBCFuncEval>& a_potential){
    m_potential = a_potential;
  }

  /*!
    @brief Factory method
  */
  virtual conductivitydomainbc_wrapper* create(const ProblemDomain& a_domain,
					       const EBISLayout&    a_ebisl,
					       const RealVect&      a_dx){
    CH_assert(m_hasbc);
    
    conductivitydomainbc_wrapper* fresh = new conductivitydomainbc_wrapper();

    fresh->set_potential(m_potential); 
    fresh->set_wallbc(m_wallbc);
    
    return fresh;
  }
  
protected:

  /*!
    @brief Make sure BCs are set
  */
  bool m_hasbc;

  /*!
    @brief Wall boundary conditions
  */
  Vector<RefCountedPtr<wall_bc> > m_wallbc;

  /*!
    @brief Potential function
  */
  RefCountedPtr<BaseBCFuncEval> m_potential;

};
#endif
