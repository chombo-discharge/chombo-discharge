/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/*
  @file   CD_EBHelmholtzDomainBCFactory.H
  @brief  Declaration of a factory class for making domain boundary condition objects for EBHelmholtzOp
  @author Robert Marskar
*/

#ifndef CD_EBHELMHOLTZDOMAINBCFACTORY_H
#define CD_EBHELMHOLTZDOMAINBCFACTORY_H

// Chombo includes
#include <EBLevelGrid.H>
#include <RefCountedPtr.H>

// Our includes
#include <CD_EBHelmholtzDomainBC.H>
#include <CD_NamespaceHeader.H>

/**
  @brief Base class for passing domain boundary conditions into EBHelmholtzOp. 
*/
class EBHelmholtzDomainBCFactory
{
public:
  /**
    @brief Constructor. 
  */
  EBHelmholtzDomainBCFactory();

  /**
    @brief Destructor
  */
  virtual ~EBHelmholtzDomainBCFactory();

  /**
    @brief Disallowed -- don't see why you would need it.
  */
  EBHelmholtzDomainBCFactory(const EBHelmholtzDomainBCFactory& a_other) = delete;

  /**
    @brief Disallowed - don't see why you would need it.
  */
  EBHelmholtzDomainBCFactory&
  operator=(const EBHelmholtzDomainBCFactory& a_other) = delete;

  /**
    @brief Factory method. Implementations must override this. 
    @return Return value
  */
  virtual RefCountedPtr<EBHelmholtzDomainBC>
  create() const = 0;
};

#include <CD_NamespaceFooter.H>

#endif
