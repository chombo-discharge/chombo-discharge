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
    @brief Disallowed weak constructor. 
  */
  MFHelmholtzDirichletEBBC();

  /*!
    @brief Default constructor
  */
  ~MFHelmholtzDirichletEBBC();

  /*!
    @brief Set jump bc. Must be done prior to calling define
    @param[in] a_jumpBC Jump bc object. 
  */
  void setJumpBC(const RefCountedPtr<JumpBC>& a_jumpBC);

protected:

  /*!
    @brief Jump bc object
  */
  RefCountedPtr<JumpBC> m_jumpBC;

  /*!
    @brief Define function. 
    @note Overriding this because this BC class must adjust the stencils for the "matching condition". 
  */
  virtual void define() override;
};

#include <CD_NamespaceFooter.H>

#endif
