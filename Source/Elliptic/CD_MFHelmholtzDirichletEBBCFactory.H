/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFHelmholtzDirichletEBBCFactory.H
  @brief  Declaration of a Dirichlet boundary condition factory class for MFHelmholtzOp.
  @author Robert Marskar
*/

#ifndef CD_MFHelmholtzDirichletEBBCFactory_H
#define CD_MFHelmholtzDirichletEBBCFactory_H

// Our includes
#include <CD_EBHelmholtzDirichletEBBCFactory.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for making Dirichlet boundary conditions on the EB in a multifluid context. 
  @note This class reuses EBHelmholtzDirichletEBBCFactory, but adjusts these on a subset of cells (the "matching" cells). 
*/
class MFHelmholtzDirichletEBBCFactory : public EBHelmholtzDirichletEBBCFactory {
public:

  /*!
    @brief Weak constructor. Must subsequently set order, weight, and value. 
  */
  MFHelmholtzDirichletEBBCFactory();

  /*!
    @brief Full constructor which sets order and value on EB. 
    @param[in] a_order  Approximation order
    @param[in] a_weight Weighting factor for stencils
    @param[in] a_value  Constant value on EB
    @note The order is also the attempted radius of the stencil. 
  */
  MFHelmholtzDirichletEBBCFactory(const int a_order, const int a_weight, const Real a_value);

  /*!
    @brief Full constructor which sets order and spatially varying value on EBs
    @param[in] a_order Approximation order
    @param[in] a_value Spatially varying value on EBs
  */
  MFHelmholtzDirichletEBBCFactory(const int a_order, const int a_weight, const std::function<Real(const RealVect& a_pos)>& a_value);

  /*!
    @brief Default constructor
  */
  ~MFHelmholtzDirichletEBBCFactory();

  /*!
    @brief Factory method
  */
  RefCountedPtr<EBHelmholtzEBBC> create() override;

};

#include <CD_NamespaceFooter.H>

#endif
