/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFHelmholtzEBBCFactory.H
  @brief  Declaration of an EB boundary condition factory class for MFHelmholtzOp.
  @author Robert Marskar
*/

#ifndef CD_MFHelmholtzEBBCFactory_H
#define CD_MFHelmholtzEBBCFactory_H

// Our includes
#include <CD_EBHelmholtzEBBC.H>
#include <CD_JumpBC.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for making  boundary conditions on the EB in a multifluid context. 
  @note This class reuses EBHelmholtzEBBCFactory, but adjusts these on a subset of cells (the "matching" cells). 
*/
class MFHelmholtzEBBCFactory {
public:

  /*!
    @brief Weak constructor. Must subsequently set order, weight, and value. 
  */
  MFHelmholtzEBBCFactory();

  /*!
    @brief Default constructor
  */
  ~MFHelmholtzEBBCFactory();

  /*!
    @brief Factory method for 
  */
  virtual RefCountedPtr<EBHelmholtzEBBC> create(const int a_iphase, const RefCountedPtr<JumpBC>& a_jumpBC) const = 0;

};

#include <CD_NamespaceFooter.H>

#endif
