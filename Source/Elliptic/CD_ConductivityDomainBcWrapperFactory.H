/*!
  @file   ConductivityDomainBcWrapperFactory.H
  @brief  Declaration of a factory class for ConductivityDomainBcWrapper
  @author Robert Marskar
  @date   Jan. 2018
*/

#ifndef _CONDUCTIVITYDOMAINBC_WRAPPER_FACTORY_
#define _CONDUCTIVITYDOMAINBC_WRAPPER_FACTORY_

#include <CD_ConductivityDomainBcWrapper.H>
#include <CD_RobinCoefficients.H>

#include <ConductivityBaseDomainBC.H>

#include "CD_NamespaceHeader.H"

/*!
  @brief Class for generating ConductivityDomainBcWrapper objects
  @details This class does not (currently) do Robin boundary conditions, although it can probably be extended to take care of 
  those cases as well. 
*/
class ConductivityDomainBcWrapperFactory : public BaseDomainBCFactory {
public:
  
  /*!
    @brief Constructor
  */
  ConductivityDomainBcWrapperFactory();

  /*!
    @brief Destructor
  */
  virtual ~ConductivityDomainBcWrapperFactory();

  virtual void setWallBc(const Vector<RefCountedPtr<WallBc> >& a_wallbc);
  virtual void setPotentials(const Vector<RefCountedPtr<BaseBCFuncEval> >& a_potentials);
  virtual void setRobinCoefficients(const Vector<RefCountedPtr<RobinCoefficients> >& a_robinco);

  /*!
    @brief Factory method
  */
  virtual ConductivityDomainBcWrapper* create(const ProblemDomain& a_domain,
					       const EBISLayout&    a_ebisl,
					       const RealVect&      a_dx);
  
protected:

  bool m_hasbc;

  Vector<RefCountedPtr<WallBc> >        m_wallbc;      // BC identifier
  Vector<RefCountedPtr<BaseBCFuncEval> > m_potentials;  // Wall potential
  Vector<RefCountedPtr<RobinCoefficients> >     m_robinco;     // Robin coefficients for Robin

};
#include "CD_NamespaceFooter.H"
#endif
