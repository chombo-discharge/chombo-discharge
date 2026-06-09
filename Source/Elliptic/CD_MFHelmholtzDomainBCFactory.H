/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/*
  @file   CD_MFHelmholtzDomainBCFactory.H
  @brief  Declaration of a factory class for making boundary condition objects for MFHelmholtzOp
  @author Robert Marskar
*/

#ifndef CD_MFHELMHOLTZDOMAINBCFACTORY_H
#define CD_MFHELMHOLTZDOMAINBCFACTORY_H

// Chombo includes
#include <RefCountedPtr.H>

// Our includes
#include <CD_EBHelmholtzDomainBC.H>
#include <CD_NamespaceHeader.H>

/**
  @brief Base class for creating domain boundary condition objects in MFHelmholtzOp. 
*/
class MFHelmholtzDomainBCFactory
{
public:
  /**
    @brief Constructor. 
  */
  MFHelmholtzDomainBCFactory();

  /**
    @brief Destructor
  */
  virtual ~MFHelmholtzDomainBCFactory();

  /**
    @brief Disallowed -- don't see why you would need it.
  */
  MFHelmholtzDomainBCFactory(const MFHelmholtzDomainBCFactory& a_other) = delete;

  /**
    @brief Disallowed - don't see why you would need it.
  */
  MFHelmholtzDomainBCFactory&
  operator=(const MFHelmholtzDomainBCFactory& a_other) = delete;

  /**
    @brief Factory method
  */
  virtual RefCountedPtr<EBHelmholtzDomainBC>
  create(const int a_iphase) const = 0;
};

#include <CD_NamespaceFooter.H>

#endif
