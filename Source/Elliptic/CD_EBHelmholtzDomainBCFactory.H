/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzDomainBCFactory.H
  @brief  Declaration of a factory class for EBHelmholtzDomainBCFactory
  @author Robert Marskar
*/

#ifndef CD_EBHelmholtzDomainBCFactory_H
#define CD_EBHelmholtzDomainBCFactory_H

// Chombo includes
#include <EBLevelGrid.H>

// Our includes
#include <CD_EBHelmholtzDomainBC.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Base class for passing EB boundary conditions into EBHelmholtzOp. 
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
    @brief Factory method
  */
  virtual EBHelmholtzDomainBC* create(const EBLevelGrid& a_eblg, const Real a_dx) = 0;

};

#include <CD_NamespaceFooter.H>

#endif
