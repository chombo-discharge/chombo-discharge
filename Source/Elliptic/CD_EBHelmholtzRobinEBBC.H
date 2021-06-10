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
    @brief Apply the EB flux in a VoF. 
  */
  virtual void applyEBFlux(EBCellFAB&         a_Lphi,
			   const EBCellFAB&   a_phi,
			   const VolIndex&    a_vof,
			   const DataIndex&   a_dit,
			   const Real&        a_beta) const override;

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

