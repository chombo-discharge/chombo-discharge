/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFHelmholtzRobinEBBCFactory.H
  @brief  Declaration of a Robin boundary condition factory class for MFHelmholtzOp.
  @author Robert Marskar
*/

#ifndef CD_MFHelmholtzRobinEBBCFactory_H
#define CD_MFHelmholtzRobinEBBCFactory_H

// Std includes
#include <functional>

// Our includes
#include <CD_MFHelmholtzEBBCFactory.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for making Robin boundary conditions on the EB in a multifluid context. 
*/
class MFHelmholtzRobinEBBCFactory : public MFHelmholtzEBBCFactory
{
public:
  /*!
    @brief Disallowed weak construction. Use the full constructors. 
  */
  MFHelmholtzRobinEBBCFactory() = delete;

  /*!
    @brief Full constructor which sets constant coefficients
    @param[in] a_order  Stencil order. Must be > 0
    @param[in] a_weight Equation weights. Must be >=0
    @param[in] a_A      Robin A-coefficient
    @param[in] a_B      Robin B-coefficient
    @param[in] a_C      Robin C-coefficient
  */
  MFHelmholtzRobinEBBCFactory(const int a_order, const int a_weight, const Real a_A, const Real a_B, const Real a_C);

  /*!
    @brief Full constructor which sets variable coefficients
    @param[in] a_order  Stencil order. Must be > 0
    @param[in] a_weight Equation weights. Must be >=0
    @param[in] a_A      Robin A-coefficient
    @param[in] a_B      Robin B-coefficient
    @param[in] a_C      Robin C-coefficient
  */
  MFHelmholtzRobinEBBCFactory(const int                                         a_order,
                              const int                                         a_weight,
                              const std::function<Real(const RealVect& a_pos)>& a_A,
                              const std::function<Real(const RealVect& a_pos)>& a_B,
                              const std::function<Real(const RealVect& a_pos)>& a_C);

  /*!
    @brief Default constructor (does nothing)
  */
  ~MFHelmholtzRobinEBBCFactory();

  /*!
    @brief Set BC order
    @param[in] a_order Order. Must be > 0
  */
  void
  setOrder(const int a_order);

  /*!
    @brief Set equation weights for least squares reconstruction. 
    @param[in] a_weight Weights. Must be >= 0
  */
  void
  setWeight(const int a_weight);

  /*!
    @brief Drop BC order if domain size is equal or below this.
  */
  void
  setDomainDropOrder(const int a_domainSize);

  /*!
    @brief Set constant Robin BC coefficients
    @param[in] a_A A-coefficient
    @param[in] a_B B-coefficient
    @param[in] a_C C-coefficient
  */
  void
  setCoefficients(const Real a_A, const Real a_B, const Real a_C);

  /*!
    @brief Set variable Robin BC coefficients
    @param[in] a_A A-coefficient
    @param[in] a_B B-coefficient
    @param[in] a_C C-coefficient
  */
  void
  setCoefficients(const std::function<Real(const RealVect& a_pos)>& a_A,
                  const std::function<Real(const RealVect& a_pos)>& a_B,
                  const std::function<Real(const RealVect& a_pos)>& a_C);

  /*!
    @brief Factory method
    @param[in] a_jumpBC Jump bc object. 
  */
  RefCountedPtr<EBHelmholtzEBBC>
  create(const int a_iphase, const RefCountedPtr<MFHelmholtzJumpBC>& a_jumpBC) const override;

protected:
  /*
    @brief Hook when using constant coefficients on the EB
  */
  bool m_useConstant;

  /*!
    @brief Hook when using variable coefficients on the EB
  */
  bool m_useFunction;

  /*!
    @brief Stencil order
  */
  int m_order;

  /*!
    @brief Equation weighting (in least squares reconstruction)
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
