/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_MFHelmholtzRobinDomainBCFactory.H
  @brief  Declaration of a factory class for making Robin domain boundary conditions with MFHelmholtzOp
  @author Robert Marskar
*/

#ifndef CD_MFHelmholtzRobinDomainBCFactory_H
#define CD_MFHelmholtzRobinDomainBCFactory_H

// Std includes
#include <functional>

// Our includes
#include <CD_MFHelmholtzDomainBCFactory.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Base class for making Robin domain boundary conditions with MFHelmholtzOp. 
  @details Robin bcs are in the form A*phi + B*dphi/dn = C. We impose dphi/dn = (C - A*phi)/B where phi is extrapolated to the EB to first order. 
  @note The B-coefficient is not the same as in the Helmholtz equation...
  @note In our notation, n points into the computational domain. 
*/
class MFHelmholtzRobinDomainBCFactory : public MFHelmholtzDomainBCFactory
{
public:
  /*!
    @brief Constructor. 
  */
  MFHelmholtzRobinDomainBCFactory();

  /*!
    @brief Full constructor which sets constant coefficients
    @param[in] a_A A-coefficient
    @param[in] a_B B-coefficient
    @param[in] a_C C-coefficient
  */
  MFHelmholtzRobinDomainBCFactory(const Real a_A, const Real a_B, const Real a_C);

  /*!
    @brief Full constructor which sets variable coefficients
    @param[in] a_A A-coefficient
    @param[in] a_B B-coefficient
    @param[in] a_C C-coefficient
  */
  MFHelmholtzRobinDomainBCFactory(const std::function<Real(const RealVect& a_pos)>& a_A,
                                  const std::function<Real(const RealVect& a_pos)>& a_B,
                                  const std::function<Real(const RealVect& a_pos)>& a_C);

  /*!
    @brief Destructor
  */
  virtual ~MFHelmholtzRobinDomainBCFactory();

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
    @brief Factory method
  */
  RefCountedPtr<EBHelmholtzDomainBC>
  create(const int a_iphase) const override;

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
