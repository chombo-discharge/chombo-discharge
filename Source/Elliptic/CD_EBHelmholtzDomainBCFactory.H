/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzDomainBCFactory.H
  @brief  Declaration of a factory class for making domain boundary condition objects for EBHelmholtzOp
  @author Robert Marskar
*/

#ifndef CD_EBHelmholtzDomainBCFactory_H
#define CD_EBHelmholtzDomainBCFactory_H

// Chombo includes
#include <EBLevelGrid.H>
#include <RefCountedPtr.H>

// Our includes
#include <CD_EBHelmholtzDomainBC.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Base class for passing domain boundary conditions into EBHelmholtzOp. 
*/
class EBHelmholtzDomainBCFactory {
public:

  /*!
    @brief Constructor. 
  */
  EBHelmholtzDomainBCFactory();

  /*!
    @brief Destructor
  */
  ~EBHelmholtzDomainBCFactory();

  /*!
    @brief Disallowed -- don't see why you would need it.
  */
  EBHelmholtzDomainBCFactory(const EBHelmholtzDomainBCFactory& a_other) = delete;

  /*!
    @brief Disallowed - don't see why you would need it.
  */
  EBHelmholtzDomainBCFactory& operator=(const EBHelmholtzDomainBCFactory& a_other) = delete;

  /*!
    @brief Factory method. Implementations must override this. 
  */
  virtual RefCountedPtr<EBHelmholtzDomainBC> create() const = 0;

};

#include <CD_NamespaceFooter.H>

#endif
