/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_MFHelmholtzNeumannEBBC.H
  @brief  Declaration of a boundary condition class for MFHelmholtzOp
  @author Robert Marskar
*/

#ifndef CD_MFHelmholtzNeumannEBBC_H
#define CD_MFHelmholtzNeumannEBBC_H

// Std includes
#include <functional>

// Our includes
#include <CD_MFHelmholtzEBBC.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Base class for passing Neumann MF boundary conditions into MFHelmholtzOp. 
  @details This allows the user to specify dphi/dn on the boundary or B*dphi/dn, either using a constant value
  or with a function. 
*/
class MFHelmholtzNeumannEBBC : public MFHelmholtzEBBC
{
public:
  /*!
    @brief Full constructor
  */
  MFHelmholtzNeumannEBBC(const int a_phase, const RefCountedPtr<MFHelmholtzJumpBC>& a_jumpBC);

  /*!
    @brief Destructor
  */
  virtual ~MFHelmholtzNeumannEBBC();

  /*!
    @brief Set value of dphi/dn on the MF. 
    @param[in] a_DphiDn Value of dphi/dn. 
    @note A positive values gives a flux INTO the domain. 
  */
  void
  setDphiDn(const Real a_DphiDn);

  /*!
    @brief Set spatially dependent value of dphi/dn on the MF. 
    @param[in] a_DphiDn Value of dphi/dn. 
    @note A positive values gives a flux INTO the domain. 
  */
  void
  setDphiDn(const std::function<Real(const RealVect& a_pos)>& a_DphiDn);

  /*!
    @brief Set the value of B*dphi/dn on the MF.
    @param[in] a_BxDphiDn Value of B*dphi/dn. 
    @note A positive values gives a flux INTO the domain. 
  */
  void
  setBxDphiDn(const Real a_BxDphiDn);

  /*!
    @brief Set the spatially dependent value of B*dphi/dn on the MF.
    @param[in] a_BxDphiDn Value of B*dphi/dn. 
    @note A positive values gives a flux INTO the domain. 
  */
  void
  setBxDphiDn(const std::function<Real(const RealVect& a_pos)>& a_BxDphiDn);

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
    @brief Value on MF
  */
  Real m_constantDphiDn;

  /*!
    @brief Value on the MF. 
  */
  std::function<Real(const RealVect& a_pos)> m_functionDphiDn;

  /*!
    @brief User define function. 
  */
  virtual void
  defineSinglePhase() override;

  /*!
    @brief Apply flux single phase
    @param[inout] a_singlePhaseVofs   Iterator for single-phase cut-cells
    @param[inout] a_Lphi              Operator kappa*L(phi)
    @param[in]    a_phi               Phi, cell-centered. 
    @param[in]    a_Bcoef             EB-centered B-coefficient
    @param[in]    a_dit               Data index    
    @param[in]    a_beta              Beta. Passed in from operator. 
    @param[in]    a_homogeneousPhysBC Homogeneous BC or not. 
  */
  virtual void
  applyEBFluxSinglePhase(VoFIterator&           a_singlePhaseVofs,
                         EBCellFAB&             a_Lphi,
                         const EBCellFAB&       a_phi,
                         const BaseIVFAB<Real>& a_Bcoef,
                         const DataIndex&       a_dit,
                         const Real&            a_beta,
                         const bool&            a_homogeneousPhysBC) const override;
};

#include <CD_NamespaceFooter.H>

#endif
