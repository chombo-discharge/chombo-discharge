/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzNeumannEBBC.H
  @brief  Declaration of a boundary condition class for EBHelmholtzOp
  @author Robert Marskar
*/

#ifndef CD_EBHelmholtzNeumannEBBC_H
#define CD_EBHelmholtzNeumannEBBC_H

// Std includes
#include <functional>

// Our includes
#include <CD_EBHelmholtzEBBC.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Base class for passing Neumann EB boundary conditions into EBHelmholtzOp. 
  @details This allows the user to specify dphi/dn on the boundary or B*dphi/dn, either using a constant value
  or with a function. 
*/
class EBHelmholtzNeumannEBBC : public EBHelmholtzEBBC
{
public:
  /*!
    @brief Default constructor. 
  */
  EBHelmholtzNeumannEBBC();

  /*!
    @brief Destructor
  */
  virtual ~EBHelmholtzNeumannEBBC();

  /*!
    @brief Set value of dphi/dn on the EB. 
    @param[in] a_DphiDn Value of dphi/dn. 
    @note A positive values gives a flux INTO the domain. 
  */
  void
  setDphiDn(const Real a_DphiDn);

  /*!
    @brief Set spatially dependent value of dphi/dn on the EB. 
    @param[in] a_DphiDn Value of dphi/dn. 
    @note A positive values gives a flux INTO the domain. 
  */
  void
  setDphiDn(const std::function<Real(const RealVect& a_pos)>& a_DphiDn);

  /*!
    @brief Set the value of B*dphi/dn on the EB.
    @param[in] a_BxDphiDn Value of B*dphi/dn. 
    @note A positive values gives a flux INTO the domain. 
  */
  void
  setBxDphiDn(const Real a_BxDphiDn);

  /*!
    @brief Set the spatially dependent value of B*dphi/dn on the EB.
    @param[in] a_BxDphiDn Value of B*dphi/dn. 
    @note A positive values gives a flux INTO the domain. 
  */
  void
  setBxDphiDn(const std::function<Real(const RealVect& a_pos)>& a_BxDphiDn);

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
    @brief Multiply by b-coefficient or not
  */
  bool m_multByBco;

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
  Real m_constantDphiDn;

  /*!
    @brief Value on the EB. 
  */
  std::function<Real(const RealVect& a_pos)> m_functionDphiDn;

  /*!
    @brief User define function. 
  */
  virtual void
  define() override;
};

#include <CD_NamespaceFooter.H>

#endif
