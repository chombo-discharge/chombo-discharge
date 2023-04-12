/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzRobinEBBC.H
  @brief  Declaration of a boundary condition class for EBHelmholtzOp
  @author Robert Marskar
*/

#ifndef CD_EBHelmholtzRobinEBBC_H
#define CD_EBHelmholtzRobinEBBC_H

// Std includes
#include <functional>

// Our includes
#include <CD_VofUtils.H>
#include <CD_EBHelmholtzEBBC.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Base class for passing Robin EB boundary conditions into EBHelmholtzOp. 
  @details Robin bcs are in the form A*phi + B*dphi/dn = C. We impose dphi/dn = (C - A*phi)/B where phi is extrapolated to the EB to first order. 
  @note The B-coefficient is not the same as in the Helmholtz equation...
  @note In our notation, n points into the computational domain. 
*/
class EBHelmholtzRobinEBBC : public EBHelmholtzEBBC
{
public:
  /*!
    @brief Default constructor. Must subsequently set coefficients.
  */
  EBHelmholtzRobinEBBC();

  /*!
    @brief Full constructor which sets constant coefficients
    @param[in] a_order  Stencil order. Must be > 0
    @param[in] a_weight Equation weights. Must be >=0
    @param[in] a_A      Robin A-coefficient
    @param[in] a_B      Robin B-coefficient
    @param[in] a_C      Robin C-coefficient
  */
  EBHelmholtzRobinEBBC(const int a_order, const int a_weight, const Real a_A, const Real a_B, const Real a_C);

  /*!
    @brief Full constructor which sets variable coefficients
    @param[in] a_order  Stencil order. Must be > 0
    @param[in] a_weight Equation weights. Must be >=0
    @param[in] a_A      Robin A-coefficient
    @param[in] a_B      Robin B-coefficient
    @param[in] a_C      Robin C-coefficient
  */
  EBHelmholtzRobinEBBC(const int                                         a_order,
                       const int                                         a_weight,
                       const std::function<Real(const RealVect& a_pos)>& a_A,
                       const std::function<Real(const RealVect& a_pos)>& a_B,
                       const std::function<Real(const RealVect& a_pos)>& a_C);

  /*!
    @brief Destructor
  */
  virtual ~EBHelmholtzRobinEBBC();

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
    @brief Drop BC order if domain size is equal or below this.
  */
  virtual void
  setDomainDropOrder(const int a_domainSize);

  /*!
    @brief Set constant coefficients
    @param[in] a_A Constant A-coefficient
    @param[in] a_B Constant B-coefficient
    @param[in] a_C Constant C-coefficient
  */
  void
  setCoefficients(const Real a_A, const Real a_B, const Real a_C);

  /*!
    @brief Set variable coefficients
    @param[in] a_A Spatially dependent A-coefficient
    @param[in] a_B Spatially dependent B-coefficient
    @param[in] a_C Spatially dependent C-coefficient
  */
  void
  setCoefficients(const std::function<Real(const RealVect& a_pos)>& a_A,
                  const std::function<Real(const RealVect& a_pos)>& a_B,
                  const std::function<Real(const RealVect& a_pos)>& a_C);

  /*!
    @brief Apply the EB flux. This is the version that is called by EBHelmholtzOp.
    @param[inout] a_vofit             Iterator for cut-cells
    @param[in]    a_Lphi              Operator kappa*L(phi)
    @param[in]    a_phi               Phi, cell-centered. 
    @param[in]    a_Bcoef             EB-centered B-coefficient
    @param[in]    a_dit               Data index    
    @param[in]    a_beta              Beta. Passed in from operator. 
    @param[in]    a_homogeneousPhysBC Homogeneous BC or not. 
  */
  virtual void
  applyEBFluxRelax(VoFIterator&           a_vofit,
                   EBCellFAB&             a_Lphi,
                   const EBCellFAB&       a_phi,
                   const BaseIVFAB<Real>& a_Bcoef,
                   const DataIndex&       a_dit,
                   const Real&            a_beta,
                   const bool&            a_homogeneousPhysBC) const override;

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
    @brief Stencil order
  */
  int m_order;

  /*!
    @brief Stencil weight (for weighting equations in least squares reconstruction)
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

  /*!
    @brief User define function. 
  */
  virtual void
  define() override;

  /*!
    @brief Get a least-squares based interpolation stencil 
    @param[in] a_vof          Grid cell
    @param[in] a_dit          Data index. 
    @param[in] a_neighborhood Vof neighborhood (either quadrant-based or radius)
    @param[in] a_order        Stencil order
    @note This extrapolates from cell-centers to the EB using a first order expression. 
  */
  VoFStencil
  getInterpolationStencil(const VolIndex&              a_vof,
                          const DataIndex&             a_dit,
                          const VofUtils::Neighborhood a_neighborhood,
                          const int                    a_order) const;
};

#include <CD_NamespaceFooter.H>

#endif
