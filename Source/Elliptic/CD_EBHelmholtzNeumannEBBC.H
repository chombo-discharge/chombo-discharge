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
  @details This allows the user to specify dphi/dn on the boundary, either using a constant value
  or with a function. 
*/
class EBHelmholtzNeumannEBBC : public EBHelmholtzEBBC {
public:

  /*!
    @brief Default constructor. 
  */
  EBHelmholtzNeumannEBBC();

  /*!
    @brief Destructor
  */
  ~EBHelmholtzNeumannEBBC();

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
  Real m_constantDphiDn;;

  /*!
    @brief Value on the EB. 
  */
  std::function<Real(const RealVect& a_pos) > m_functionDphiDn;

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
  
};

#include <CD_NamespaceFooter.H>

#endif
