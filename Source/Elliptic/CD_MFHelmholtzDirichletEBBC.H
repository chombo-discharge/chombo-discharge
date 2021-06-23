/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFHelmholtzDirichletEBBC.H
  @brief  Declaration of a Dirichlet boundary condition class for MFHelmholtzOp.
  @author Robert Marskar
*/

#ifndef CD_MFHelmholtzDirichletEBBC_H
#define CD_MFHelmholtzDirichletEBBC_H

// Our includes
#include <CD_EBHelmholtzDirichletEBBC.H>
#include <CD_JumpBC.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for making Dirichlet boundary conditions on the EB in a multifluid context. 
  @note This class reuses the stencils from EBHelmholtzDirichletEBBC, but adjusts these on a subset of cells (the "matching" cells). 
*/
class MFHelmholtzDirichletEBBC : public EBHelmholtzDirichletEBBC {
public:

  /*!
    @brief Disallowed weak constructor. Factory must call the other version. 
  */
  MFHelmholtzDirichletEBBC() = delete;

  /*!
    @brief Full constructor
  */
  MFHelmholtzDirichletEBBC(const RefCountedPtr<JumpBC>& a_jumpBC);

  /*!
    @brief Default constructor
  */
  ~MFHelmholtzDirichletEBBC();

protected:

  /*!
    @brief Jump bc object
  */
  RefCountedPtr<JumpBC> m_jumpBC;

  /*!
    @brief Phase
  */
  int m_phase;
};

#include <CD_NamespaceFooter.H>

#endif
