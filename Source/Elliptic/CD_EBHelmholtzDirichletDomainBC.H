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
class EBHelmholtzDirichletDomainBC : public EBHelmholtzDomainBC
{
public:
  /*!
    @brief Constructor. Must subsequently call setValue. 
  */
  EBHelmholtzDirichletDomainBC();

  /*!
    @brief Constructor. 
    @param[in] a_value Constant value on all domain sides. 
  */
  EBHelmholtzDirichletDomainBC(const Real a_value);

  /*!
    @brief Full constructor. 
    @param[in] a_value Spatially dependent value on all domain sides. 
  */
  EBHelmholtzDirichletDomainBC(const std::function<Real(const RealVect& a_pos)>& a_value);

  /*!
    @brief Destructor
  */
  virtual ~EBHelmholtzDirichletDomainBC();

  /*!
    @brief Set the constant value on the domain side. 
    @param[in] a_value Constant value on all domain sides. 
  */
  void
  setValue(const Real a_value);

  /*!
    @brief Set the spatially dependent value on the domain side. 
    @param[in] a_value Spatially dependent value on all domain sides. 
  */
  void
  setValue(const std::function<Real(const RealVect& a_pos)>& a_value);

  /*!
    @brief Get face flux. This is for regular cells, and it uses first order finite differences for filling the flux. 
    @param[inout] a_faceFlux       Face-centered flux. 
    @param[in]    a_phi            Cell-centered data
    @param[in]    a_Bcoef          Face-centered B-coefficient
    @param[in]    a_dir            Coordinate direction
    @param[in]    a_side           High or low side
    @param[in]    a_dit            Data index
    @param[in]    a_useHomogeneous Homogeneous physical BCs or not. 
    @note a_faceFlux is defined in the first strip of interior cells on the (dir,side) edge of the domain. 
  */
  virtual void
  getFaceFlux(BaseFab<Real>&        a_faceFlux,
              const BaseFab<Real>&  a_phi,
              const BaseFab<Real>&  a_Bcoef,
              const int&            a_dir,
              const Side::LoHiSide& a_side,
              const DataIndex&      a_dit,
              const bool            a_useHomogeneous) const override;

  /*!
    @brief Get face flux. This is for irregular cells. 
    @param[in] a_vof            Grid cell
    @param[in] a_phi            Cell-centered data
    @param[in]    a_Bcoef          Face-centered B-coefficient
    @param[in] a_dir            Coordinate direction
    @param[in] a_side           High or low side
    @param[in] a_dit            Data index
    @param[in] a_useHomogeneous Homogeneous physical BCs or not. 
  */
  virtual Real
  getFaceFlux(const VolIndex&       a_vof,
              const EBCellFAB&      a_phi,
              const EBFaceFAB&      a_Bcoef,
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
    @brief Constant value on domain sides. 
  */
  Real m_constantValue;

  /*!
    @brief Spatially dependent value on domain sides
  */
  std::function<Real(const RealVect& a_pos)> m_functionValue;
};

#include <CD_NamespaceFooter.H>

#endif
