/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_MFHelmholtzEBBCFactory.H
  @brief  Declaration of an EB boundary condition factory class for MFHelmholtzOp.
  @author Robert Marskar
*/

#ifndef CD_MFHELMHOLTZEBBCFACTORY_H
#define CD_MFHELMHOLTZEBBCFACTORY_H

// Our includes
#include <CD_EBHelmholtzEBBC.H>
#include <CD_MFHelmholtzJumpBC.H>
#include <CD_NamespaceHeader.H>

/**
  @brief Base class for making boundary conditions on the EB in a multifluid context. 
*/
class MFHelmholtzEBBCFactory
{
public:
  /**
    @brief Weak constructor. Must subsequently set order, weight, and value. 
  */
  MFHelmholtzEBBCFactory();

  /**
    @brief Default constructor
  */
  virtual ~MFHelmholtzEBBCFactory();

  /**
    @brief Factory method for 
  */
  virtual RefCountedPtr<EBHelmholtzEBBC>
  create(const int a_iphase, const RefCountedPtr<MFHelmholtzJumpBC>& a_jumpBC) const = 0;
};

#include <CD_NamespaceFooter.H>

#endif
