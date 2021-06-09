/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzEBBCFactoryFactory.H
  @brief  Declaration of a factory class for EBHelmholtzEBBCFactory
  @author Robert Marskar
*/

#ifndef CD_EBHelmholtzEBBCFactoryFactory_H
#define CD_EBHelmholtzEBBCFactoryFactory_H

// Chombo includes
#include <EBLevelGrid.H>

// Our includes
#include <CD_EBHelmholtzEBBC.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Base class for passing EB boundary conditions into EBHelmholtzOp. 
*/
class EBHelmholtzEBBCFactory {
public:

  /*!
    @brief Constructor. 
  */
  EBHelmholtzEBBCFactory();

  /*!
    @brief Destructor
  */
  ~EBHelmholtzEBBCFactory();

  /*!
    @brief Disallowed -- don't see why you would need it.
  */
  EBHelmholtzEBBCFactory(const EBHelmholtzEBBCFactory& a_other) = delete;

  /*!
    @brief Disallowed - don't see why you would need it.
  */
  EBHelmholtzEBBCFactory& operator=(const EBHelmholtzEBBCFactory& a_other) = delete;

  /*!
    @brief Factory method
  */
  virtual EBHelmholtzEBBC* create(const EBLevelGrid& a_eblg, const Real a_dx) = 0;

};

#include <CD_NamespaceFooter.H>

#endif
