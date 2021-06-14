/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzDirichletDomainBC.H
  @brief  Declaration of a Dirichlet boundary condition class for EBHelmholtzOp
  @author Robert Marskar
*/

#ifndef CD_EBHelmholtzDirichletDomainBC_H
#define CD_EBHelmholtzDirichletDomainBC_H

// Std includes
#include <functional>

// Our includes
#include <CD_EBHelmholtzDomainBC.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Boundary condition class for passing Dirichlet domain boundary conditions into EBHelmholtzOp
*/
class EBHelmholtzDirichletDomainBC : public EBHelmholtzDomainBC {
public:

  /*!
    @brief Constructor. Must subsequently call setValue. 
  */
  EBHelmholtzDirichletDomainBC();

  /*!
    @brief Constructor. 
  */
  EBHelmholtzDirichletDomainBC(const Real a_value);

    /*!
    @brief Constructor. 
  */
  EBHelmholtzDirichletDomainBC(const std::function<Real(const RealVect& a_pos)>& a_value);

  /*!
    @brief Destructor
  */
  ~EBHelmholtzDirichletDomainBC();

  /*!
    @brief Set the value on the domain side. 
  */
  void setValue(const Real a_value);

  /*!
    @brief Set the value on the domain side. 
  */
  void setValue(const std::function<Real(const RealVect& a_pos)>& a_value);

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
  Real m_constantValue;

  /*!
    @brief Value on the EB. 
  */
  std::function<Real(const RealVect& a_pos) > m_functionValue;
};

#include <CD_NamespaceFooter.H>

#endif
