/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzNeumannDomainBC.H
  @brief  Declaration of a Neumann boundary condition class for EBHelmholtzOp
  @author Robert Marskar
*/

#ifndef CD_EBHelmholtzNeumannDomainBC_H
#define CD_EBHelmholtzNeumannDomainBC_H

// Std includes
#include <functional>

// Our includes
#include <CD_EBHelmholtzDomainBC.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Boundary condition class for passing Neumann domain boundary conditions into EBHelmholtzOp
*/
class EBHelmholtzNeumannDomainBC : public EBHelmholtzDomainBC {
public:

  /*!
    @brief Constructor. Must subsequently call setValue. 
  */
  EBHelmholtzNeumannDomainBC();

  /*!
    @brief Constructor. 
  */
  EBHelmholtzNeumannDomainBC(const Real a_value);

    /*!
    @brief Constructor. 
  */
  EBHelmholtzNeumannDomainBC(const std::function<Real(const RealVect& a_pos)>& a_value);

  /*!
    @brief Destructor
  */
  ~EBHelmholtzNeumannDomainBC();

  /*!
    @brief Set the value of dphi/dn on the EB. 
  */
  void setDphiDn(const int a_DphiDn);

  /*!
    @brief Set the value of dphi/dn on the EB. 
  */
  void setDphiDn(const std::function<Real(const RealVect& a_pos)>& a_DphiDn);

  /*!
    @brief Set the value of B*dphi/dn on the EB.
  */
  void setBxDphiDn(const int a_BxDphiDn);

  /*!
    @brief Set the value of B*dphi/dn on the EB.
  */
  void setBxDphiDn(const std::function<Real(const RealVect& a_pos)>& a_BxDphiDn);

  /*!
    @brief Get face flux. This is for regular cells
  */
  virtual void getFaceFlux(BaseFab<Real>&        a_faceFlux,
			   const BaseFab<Real>&  a_phi,
			   const int&            a_dir,
			   const Side::LoHiSide& a_side,
			   const DataIndex&      a_dit,
			   const bool            a_useHomogeneous) const override;

  /*!
    @brief Get face flux. This is for irregular cells. 
  */
  virtual Real getFaceFlux(const VolIndex&       a_vof,
			   const EBCellFAB&      a_phi,
			   const int&            a_dir,
			   const Side::LoHiSide& a_side,
			   const DataIndex&      a_dit,
			   const bool            a_useHomogeneous) const override;

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
  std::function<Real(const RealVect& a_pos) > m_functionDphiDn;
};

#include <CD_NamespaceFooter.H>

#endif
