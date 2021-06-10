/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzDirichletEBBC.H
  @brief  Declaration of a boundary condition class for EBHelmholtzOp
  @author Robert Marskar
*/

#ifndef CD_EBHelmholtzDirichletEBBC_H
#define CD_EBHelmholtzDirichletEBBC_H

// Std includes
#include <functional>

// Our includes
#include <CD_EBHelmholtzEBBC.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Base class for passing Dirichlet EB boundary conditions into EBHelmholtzOp. 
*/
class EBHelmholtzDirichletEBBC : public EBHelmholtzEBBC {
public:

  /*!
    @brief Default constructor. 
  */
  EBHelmholtzDirichletEBBC();

  /*!
    @brief Destructor
  */
  ~EBHelmholtzDirichletEBBC();

  /*!
    @brief Set BC order
  */
  void setOrder(const int a_order);

  /*!
    @brief Set the value on the EB. 
  */
  void setValue(const int a_value);

  /*!
    @brief Set the value on the EB. 
  */
  void setValue(const std::function<Real(const RealVect& a_pos)>& a_value);

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
    @brief Stencil approximation order,
  */
  int m_order;

  /*!
    @brief Value on EB
  */
  Real m_constantValue;

  /*!
    @brief Value on the EB. 
  */
  std::function<Real(const RealVect& a_pos) > m_functionValue;
  
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
			   const Real&        a_factor,
			   const bool&        a_useHomogeneous) const override;
};

#include <CD_NamespaceFooter.H>

#endif
