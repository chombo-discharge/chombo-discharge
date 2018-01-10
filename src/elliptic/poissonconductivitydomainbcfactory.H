/*!
  @file poissonconductivitydomainbcfactory.H
  @brief Declaration of a factory class for poissonconductivitydomainbc
  @author Robert Marskar
  @date Jan. 2018
  @todo Segregate implementation file
*/

#ifndef _POISSONCONDUCTIVITYDOMAINBCFACTORY_
#define _POISSONCONDUCTIVITYDOMAINBCFACTORY_

#include "poissonconductivitydomainbc.H"

#include <ConductivityBaseDomainBC.H>

/*!
  @brief Class for generating poissonconductivitydomainbc objects
  @details This class does not (currently) do Robin boundary conditions, although it can probably be extended to take care of 
  those cases as well. 
*/
class poissonconductivitydomainbcfactory : public BaseDomainBCFactory {
public:
  
  /*!
    @brief Constructor
  */
  poissonconductivitydomainbcfactory(){
  }

  /*!
    @brief Destructor
  */
  ~poissonconductivitydomainbcfactory(){}


  /*!
    @brief Set function
  */
  virtual void set_function(RefCountedPtr<BaseBCFuncEval> a_func){
    m_flux = a_func;
  }

  /*!
    @brief Set value
  */
  virtual void set_value(Real a_value){
    m_value = a_value;
  }

  /*!
    @brief Factory method
  */
  virtual poissonconductivitydomainbc* create(const ProblemDomain& a_domain,
					      const EBISLayout&    a_ebisl,
					      const RealVect&      a_dx){

    poissonconductivitydomainbc* fresh = new poissonconductivitydomainbc();

    if(m_isfunction){
      fresh->set_function(m_flux);
    }
    else{
      fresh->set_value(m_value);
    }
    
    return fresh;
  }
  
protected:

  /*!
    @brief Check if function
  */
  bool m_isfunction;

  /*!
    @brief Constant value
  */
  Real m_value;

  /*!
    @brief Function
  */
  RefCountedPtr<BaseBCFuncEval> m_flux;

};
#endif
