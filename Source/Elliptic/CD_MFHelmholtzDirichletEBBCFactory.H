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
*/
class MFHelmholtzDirichletEBBCFactory : public MFHelmholtzEBBCFactory
{
public:
  /*!
    @brief Weak constructor. Must subsequently set order, weight, and value. 
  */
  MFHelmholtzDirichletEBBCFactory() = delete;

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
  MFHelmholtzDirichletEBBCFactory(const int                                         a_order,
                                  const int                                         a_weight,
                                  const std::function<Real(const RealVect& a_pos)>& a_value);

  /*!
    @brief Default constructor
  */
  ~MFHelmholtzDirichletEBBCFactory();

  /*!
    @brief Set the value on the EB. 
    @param[in] a_value Constant value on EB
  */
  virtual void
  setValue(const Real a_value);

  /*!
    @brief Set the spatially varying value on the EB. 
    @param[in] a_value Spatially dependent value on EB
  */
  virtual void
  setValue(const std::function<Real(const RealVect& a_pos)>& a_value);

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
    @brief Factory method
    @param[in] a_jumpBC Jump bc object. 
  */
  RefCountedPtr<EBHelmholtzEBBC>
  create(const int a_iphase, const RefCountedPtr<MFHelmholtzJumpBC>& a_jumpBC) const override;

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
    @brief Value on EB
  */
  Real m_constantValue;

  /*!
    @brief Stencil order
  */
  int m_order;

  /*!
    @brief Stencil weight (for weighted least squares)
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
    @brief Value on the EB. 
  */
  std::function<Real(const RealVect& a_pos)> m_functionValue;
};

#include <CD_NamespaceFooter.H>

#endif
