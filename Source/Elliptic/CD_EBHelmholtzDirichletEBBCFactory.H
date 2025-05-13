/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzDirichletEBBCFactory.H
  @brief  Declaration of a factory class for making Dirichlet boundary conditions on EBs for EBHelmholtzOp
  @author Robert Marskar
*/

#ifndef CD_EBHelmholtzDirichletEBBCFactory_H
#define CD_EBHelmholtzDirichletEBBCFactory_H

// Our includes
#include <CD_EBHelmholtzEBBCFactory.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Factory class for making Dirichlet BCs on EBs in EBHelmholtzOp
*/
class EBHelmholtzDirichletEBBCFactory : public EBHelmholtzEBBCFactory
{
public:
  /*!
    @brief Constructor. 
  */
  EBHelmholtzDirichletEBBCFactory();

  /*!
    @brief Full constructor which sets order and value on EB. 
    @param[in] a_order  Approximation order
    @param[in] a_weight Weighting factor for stencils
    @param[in] a_value  Constant value on EB
    @note The order is also the attempted radius of the stencil. 
  */
  EBHelmholtzDirichletEBBCFactory(const int a_order, const int a_weight, const Real a_value);

  /*!
    @brief Full constructor which sets order and spatially varying value on EBs
    @param[in] a_order Approximation order
    @param[in] a_value Spatially varying value on EBs
  */
  EBHelmholtzDirichletEBBCFactory(const int                                         a_order,
                                  const int                                         a_weight,
                                  const std::function<Real(const RealVect& a_pos)>& a_value);

  /*!
    @brief Destructor
  */
  virtual ~EBHelmholtzDirichletEBBCFactory();

  /*!
    @brief Set BC order
    @param[in] a_order Approximation order (1 and 2 are accepted)
    @note For least squares, the order is also the attempted radius of the stencil. 
  */
  virtual void
  setOrder(const int a_order);

  /*!
    @brief Set BC order
    @param[in] a_weight Weight for least squares stencils. 
    @note a_weight = 0 can lead to ill-conditioned least-squares systems. 
  */
  virtual void
  setWeight(const int a_weight);

  /*!
    @brief Set the value on the EB. 
    @param[in] a_value Constant value on EB
  */
  virtual void
  setValue(const Real a_value);

  /*!
    @brief Drop BC order if domain size is equal or below this.
  */
  virtual void
  setDomainDropOrder(const int a_domainSize);

  /*!
    @brief Signal that it should drop the stencil order if the cell is a coarsened cell
    @param[in] a_dropOrder
  */
  virtual void
  setCoarseGridDropOrder(const bool a_dropOrder);

  /*!
    @brief Set the spatially varying value on the EB. 
    @param[in] a_value Spatially dependent value on EB
  */
  virtual void
  setValue(const std::function<Real(const RealVect& a_pos)>& a_value);

  /*!
    @brief Factory method
  */
  virtual RefCountedPtr<EBHelmholtzEBBC>
  create() override;

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
    @brief Special flag for dropping stencil order when domains become coarser than this.
  */
  int m_domainDropOrder;

  /*!
    @brief Special flag for dropping stencil order if the cell is a coarsened cell
  */
  bool m_dropOrder;

  /*!
    @brief Value on EB
  */
  Real m_constantValue;

  /*!
    @brief Value on the EB. 
  */
  std::function<Real(const RealVect& a_pos)> m_functionValue;
};

#include <CD_NamespaceFooter.H>

#endif
