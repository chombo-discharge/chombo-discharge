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
  @details This class can set either dphi/dn or B*dphi/dn on the boundary. 
*/
class EBHelmholtzNeumannDomainBC : public EBHelmholtzDomainBC
{
public:
  /*!
    @brief Constructor. Must subsequently call setDphiDn or setBxDphiDn
  */
  EBHelmholtzNeumannDomainBC();

  /*!
    @brief Constructor. 
    @param[in] a_DphiDn Constant value of dphi/dn on all domain sides. 
    @note A positive values gives a flux INTO the domain. 
  */
  EBHelmholtzNeumannDomainBC(const Real a_DphiDn);

  /*!
    @brief Full constructor. 
    @param[in] a_DphiDn Spatially dependent value of dphi/dn on all domain sides. 
    @note A positive values gives a flux INTO the domain. 
  */
  EBHelmholtzNeumannDomainBC(const std::function<Real(const RealVect& a_pos)>& a_value);

  /*!
    @brief Destructor
  */
  virtual ~EBHelmholtzNeumannDomainBC();

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
    @brief Get face flux. This is for regular cells
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
    @param[in] a_Bcoef          Face-centered B-coefficient
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
    @brief Hook for multiplying by b-coefficient or not
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
    @brief Value on the domain side. 
  */
  Real m_constantDphiDn;

  /*!
    @brief Value on the domain side. 
  */
  std::function<Real(const RealVect& a_pos)> m_functionDphiDn;
};

#include <CD_NamespaceFooter.H>

#endif
