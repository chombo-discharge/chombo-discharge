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
  @note This class reuses EBHelmholtzRobinEBBCFactory, but adjusts these on a subset of cells (the "matching" cells). 
*/
class MFHelmholtzRobinEBBCFactory : public MFHelmholtzEBBCFactory {
public:

  /*!
    @brief Weak constructor. Must subsequently set order, weight, and value. 
  */
  MFHelmholtzRobinEBBCFactory() = delete;

  /*!
    @brief Full constructor which sets order and value on EB. 
    @param[in] a_A      Constant A-coefficient
    @param[in] a_B      Constant B-coefficient
    @param[in] a_C      Constant C-coefficient
    @note The order is also the attempted radius of the stencil. 
  */
  MFHelmholtzRobinEBBCFactory(const Real a_A,
			      const Real a_B,
			      const Real a_C);

  /*!
    @brief Full constructor which sets order and spatially varying value on EBs
    @param[in] a_order  Approximation order
    @param[in] a_weight Weighting factor for stencils
    @param[in] a_A      Spatially varying A-coefficient
    @param[in] a_B      Spatially varying B-coefficient
    @param[in] a_C      Spatially avrying C-coefficient

  */
  MFHelmholtzRobinEBBCFactory(const std::function<Real(const RealVect& a_pos)>& a_A,
			      const std::function<Real(const RealVect& a_pos)>& a_B,
			      const std::function<Real(const RealVect& a_pos)>& a_C);

  /*!
    @brief Default constructor
  */
  ~MFHelmholtzRobinEBBCFactory();

  /*!
    @brief Set constant coefficients
    @param[in] a_A A-coefficient
    @param[in] a_B B-coefficient
    @param[in] a_C C-coefficient
  */
  void setCoefficients(const Real a_A, const Real a_B, const Real a_C);

  /*!
    @brief Set variable coefficients
    @param[in] a_A A-coefficient
    @param[in] a_B B-coefficient
    @param[in] a_C C-coefficient
  */
  void setCoefficients(const std::function<Real(const RealVect& a_pos) >& a_A,
		       const std::function<Real(const RealVect& a_pos) >& a_B,
		       const std::function<Real(const RealVect& a_pos) >& a_C);

  /*!
    @brief Factory method
    @param[in] a_jumpBC Jump bc object. 
  */
  RefCountedPtr<EBHelmholtzEBBC> create(const int a_iphase, const RefCountedPtr<JumpBC>& a_jumpBC) const override;

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
