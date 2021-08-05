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
#include <CD_EBHelmholtzEBBC.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Base class for passing Robin EB boundary conditions into EBHelmholtzOp. 
  @details Robin bcs are in the form A*phi + B*dphi/dn = C. We impose dphi/dn = (C - A*phi)/B where phi is extrapolated to the EB to first order. 
  @note The B-coefficient is not the same as in the Helmholtz equation...
  @note In our notation, n points into the computational domain. 
*/
class EBHelmholtzRobinEBBC : public EBHelmholtzEBBC {
public:

  /*!
    @brief Default constructor. 
  */
  EBHelmholtzRobinEBBC();

  /*!
    @brief Destructor
  */
  ~EBHelmholtzRobinEBBC();
  
  /*!
    @brief Set constant coefficients
    @param[in] a_A Constant A-coefficient
    @param[in] a_B Constant B-coefficient
    @param[in] a_C Constant C-coefficient
  */
  void setCoefficients(const Real a_A, const Real a_B, const Real a_C);

  /*!
    @brief Set variable coefficients
    @param[in] a_A Spatially dependent A-coefficient
    @param[in] a_B Spatially dependent B-coefficient
    @param[in] a_C Spatially dependent C-coefficient
  */
  void setCoefficients(const std::function<Real(const RealVect& a_pos) >& a_A,
		       const std::function<Real(const RealVect& a_pos) >& a_B,
		       const std::function<Real(const RealVect& a_pos) >& a_C);


  /*!
    @brief Apply the EB flux. This is the version that is called by EBHelmholtzOp.
    @param[inout] a_vofit             Iterator for cut-cells
    @param[in]    a_Lphi              Operator kappa*L(phi)
    @param[in]    a_phi               Phi, cell-centered. 
    @param[in]    a_dit               Data index    
    @param[in]    a_beta              Beta. Passed in from operator. 
    @param[in]    a_homogeneousPhysBC Homogeneous BC or not. 
  */
  virtual void applyEBFlux(VoFIterator&       a_vofit,
			   EBCellFAB&         a_Lphi,
			   const EBCellFAB&   a_phi,
			   const DataIndex&   a_dit,
			   const Real&        a_beta,
			   const bool&        a_homogeneousPhysBC) const override;
  
protected:

  /*!
    @brief Stencil order
  */
  int m_order;

  /*!
    @brief Stencil radius
  */
  int m_radius;

  /*!
    @brief Stencil weight
  */
  int m_weight;  

  /*!
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
    
  /*!
    @brief User define function. 
  */
  virtual void define() override;

  /*!
    @brief Get an interpolation stencil 
    @param[in] a_vof Grid cell
    @param[in] a_dit Data index. 
    @note This extrapolates from cell-centers to the EB using a first order expression. 
  */
  VoFStencil getInterpolationStencil(const VolIndex& a_vof, const DataIndex& a_dit) const;

  /*!
    @brief Get an interpolation stencil using all vofs in a monotone path (avoids reaching over EBs)
    @param[in] a_vof Grid cell
    @param[in] a_dit Data index. 
  */
  VoFStencil getMonoPathStencil(const VolIndex& a_vof, const DataIndex& a_dit) const;

  /*!
    @brief Get an interpolation stencil using Taylor extrapolation. 
    @param[in] a_vof Grid cell
    @param[in] a_dit Data index. 
  */
  VoFStencil getTaylorStencil(const VolIndex& a_vof, const DataIndex& a_dit) const;
};

#include <CD_NamespaceFooter.H>

#endif

