/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzEBBCFactory.H
  @brief  Declaration of a factory class for EBHelmholtzEBBCFactory
  @author Robert Marskar
*/

#ifndef CD_EBHelmholtzEBBCFactory_H
#define CD_EBHelmholtzEBBCFactory_H

// Chombo includes
#include <EBLevelGrid.H>

// Our includes
#include <CD_EBHelmholtzEBBC.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Factory class for making EBHelmholtzEBBCs. 
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
    @brief Factory method which returns new BC object. 
    @note  Defines are done in EBHelmholtzOp so no need to worry about that here. 
  */
  virtual RefCountedPtr<EBHelmholtzEBBC> create() = 0;
};

#include <CD_NamespaceFooter.H>

#endif
