/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzRobinEBBCFactory.H
  @brief  Declaration of a factory class for making Robin EB boundary conditions with EBHelmholtzOp
  @author Robert Marskar
*/

#ifndef CD_EBHelmholtzRobinEBBCFactory_H
#define CD_EBHelmholtzRobinEBBCFactory_H

// Std includes
#include <functional>

// Our includes
#include <CD_EBHelmholtzEBBCFactory.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Base class for making Robin boundary conditions with EBHelmholtzOp. 
  @details Robin bcs are in the form A*phi + B*dphi/dn = C. We impose dphi/dn = (C - A*phi)/B where phi is extrapolated to the EB to first order. 
  @note The B-coefficient is not the same as in the Helmholtz equation...
  @note In our notation, n points into the computational domain. 
*/
class EBHelmholtzRobinEBBCFactory : public EBHelmholtzEBBCFactory
{
public:
  /*!
    @brief Constructor. 
  */
  EBHelmholtzRobinEBBCFactory();

  /*!
    @brief Full constructor which sets constant coefficients
    @param[in] a_order  Stencil order. Must be > 0
    @param[in] a_weight Equation weights. Must be >=0
    @param[in] a_A      Robin A-coefficient
    @param[in] a_B      Robin B-coefficient
    @param[in] a_C      Robin C-coefficient
  */
  EBHelmholtzRobinEBBCFactory(const int a_order, const int a_weight, const Real a_A, const Real a_B, const Real a_C);

  /*!
    @brief Full constructor which sets variable coefficients
    @param[in] a_order  Stencil order. Must be > 0
    @param[in] a_weight Equation weights. Must be >=0
    @param[in] a_A      Robin A-coefficient
    @param[in] a_B      Robin B-coefficient
    @param[in] a_C      Robin C-coefficient
  */
  EBHelmholtzRobinEBBCFactory(const int                                         a_order,
                              const int                                         a_weight,
                              const std::function<Real(const RealVect& a_pos)>& a_A,
                              const std::function<Real(const RealVect& a_pos)>& a_B,
                              const std::function<Real(const RealVect& a_pos)>& a_C);

  /*!
    @brief Destructor (does nothing). 
  */
  virtual ~EBHelmholtzRobinEBBCFactory();

  /*!
    @brief Set BC order
    @param[in] a_order Order
  */
  void
  setOrder(const int a_order);

  /*!
    @brief Set equation weights for least squares reconstruction. 
    @param[in] a_weight Weights
  */
  void
  setWeight(const int a_weight);

  /*!
    @brief Set constant coefficients
    @param[in] a_A A-coefficient
    @param[in] a_B B-coefficient
    @param[in] a_C C-coefficient
  */
  void
  setCoefficients(const Real a_A, const Real a_B, const Real a_C);

  /*!
    @brief Set variable coefficients
    @param[in] a_A A-coefficient
    @param[in] a_B B-coefficient
    @param[in] a_C C-coefficient
  */
  void
  setCoefficients(const std::function<Real(const RealVect& a_pos)>& a_A,
                  const std::function<Real(const RealVect& a_pos)>& a_B,
                  const std::function<Real(const RealVect& a_pos)>& a_C);

  /*!
    @brief Drop BC order if domain size is equal or below this.
  */
  virtual void
  setDomainDropOrder(const int a_domainSize);

  /*!
    @brief Factory method. 
    @return An EBHelmholtzRobinEBBC object with specified order and coefficients. 
  */
  RefCountedPtr<EBHelmholtzEBBC>
  create() override;

protected:
  /*
    @brief Use constant for BC
  */
  bool m_useConstant;

  /*!
    @brief Use function for BC value
  */
  bool m_useFunction;

  /*!
    @brief (Attempted) stencil order
  */
  int m_order;

  /*!
    @brief Weight for least squares gradient reconstruction. 
  */
  int m_weight;

  /*!
    @brief Special flag for dropping stencil order when domains become coarser than this.
  */
  int m_domainDropOrder;

  /*!
    @brief Constant A-coefficient
  */
  Real m_constantA;

  /*!
    @brief Constant B-coefficient
  */
  Real m_constantB;

  /*!
    @brief Constant C-coefficient
  */
  Real m_constantC;

  /*!
    @brief Function-based A-coefficient
  */
  std::function<Real(const RealVect& a_pos)> m_functionA;

  /*!
    @brief Function-based B-coefficient
  */
  std::function<Real(const RealVect& a_pos)> m_functionB;

  /*!
    @brief Function-based C-coefficient
  */
  std::function<Real(const RealVect& a_pos)> m_functionC;
};

#include <CD_NamespaceFooter.H>

#endif
