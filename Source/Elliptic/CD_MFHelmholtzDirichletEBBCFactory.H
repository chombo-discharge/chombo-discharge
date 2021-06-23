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
#include <CD_MFHelmholtzEBBCFactory.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for making Dirichlet boundary conditions on the EB in a multifluid context. 
  @note This class reuses EBHelmholtzDirichletEBBCFactory, but adjusts these on a subset of cells (the "matching" cells). 
*/
class MFHelmholtzDirichletEBBCFactory : public MFHelmholtzEBBCFactory {
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
    @brief Set BC order
    @param[in] a_order Approximation order (1 and 2 are accepted)
    @note For least squares, the order is also the attempted radius of the stencil. 
  */
  virtual void setOrder(const int a_order);

  /*!
    @brief Set BC order
    @param[in] a_weight Weight for least squares stencils. 
    @note a_weight = 0 can lead to ill-conditioned least-squares systems. 
  */
  virtual void setWeight(const int a_weight);

  /*!
    @brief Set the value on the EB. 
    @param[in] a_value Constant value on EB
  */
  virtual void setValue(const Real a_value);

  /*!
    @brief Set the spatially varying value on the EB. 
    @param[in] a_value Spatially dependent value on EB
  */
  virtual void setValue(const std::function<Real(const RealVect& a_pos)>& a_value);

  /*!
    @brief Factory method
    @param[in] a_jumpBC Jump bc object. 
  */
  RefCountedPtr<EBHelmholtzEBBC> create(const int a_iphase, const RefCountedPtr<JumpBC>& a_jumpBC) const override;

protected:

  /*!
    @brief Use constant for BC
  */
  bool m_useConstant;

  /*!
    @brief Use function for BC value
  */
  bool m_useFunction;

  /*!
    @brief Stencil approximation order,
  */
  int m_order;

  /*!
    @brief Stencil weighting
  */
  int m_weight;

  /*!
    @brief Value on EB
  */
  Real m_constantValue;

  /*!
    @brief Value on the EB. 
  */
  std::function<Real(const RealVect& a_pos) > m_functionValue;
};

#include <CD_NamespaceFooter.H>

#endif
