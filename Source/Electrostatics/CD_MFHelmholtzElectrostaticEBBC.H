/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFHelmholtzElectrostaticEBBC.H
  @brief  Declaration of a Electrostatic boundary condition class for MFHelmholtzOp.
  @author Robert Marskar
*/

#ifndef CD_MFHelmholtzElectrostaticEBBC_H
#define CD_MFHelmholtzElectrostaticEBBC_H

// Our includes
#include <CD_MFHelmholtzEBBC.H>
#include <CD_ElectrostaticEbBc.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for making Electrostatic boundary conditions on the EB in a multifluid context. 
  @details This implements MFHelmholtzEBBC, which is essentially an EBHelmholtzEBBC-derived class that performs automatic "matching" on cells
  with jumps in the b-coefficients. Users will never have to interact directly with this object. 
*/
class MFHelmholtzElectrostaticEBBC : public MFHelmholtzEBBC
{
public:
  /*!
    @brief Full constructor
    @param[in] a_phase            Phase
    @param[in] a_electrostaticBCs Electrostatic boundary conditions on EBs
    @param[in] a_jumpBC           Jumpc BC object. 
  */
  MFHelmholtzElectrostaticEBBC(const int                               a_phase,
                               const ElectrostaticEbBc&                a_electrostaticBCs,
                               const RefCountedPtr<MFHelmholtzJumpBC>& a_jumpBC);

  /*!
    @brief No weak construction.
  */
  MFHelmholtzElectrostaticEBBC() = delete;

  /*!
    @brief No copy construction
    @param[in] a_other Other object
  */
  MFHelmholtzElectrostaticEBBC(const MFHelmholtzElectrostaticEBBC& a_other) = delete;

  /*!
    @brief No move construction
    @param[in] a_other Other object
  */
  MFHelmholtzElectrostaticEBBC(const MFHelmholtzElectrostaticEBBC&& a_other) = delete;

  /*!
    @brief Default constructor (does nothing). 
  */
  virtual ~MFHelmholtzElectrostaticEBBC();

  /*!
    @brief No copy assignment
    @param[in] a_other Other object
  */
  MFHelmholtzElectrostaticEBBC&
  operator=(const MFHelmholtzElectrostaticEBBC& a_other) = delete;

  /*!
    @brief No move assignment
    @param[in] a_other Other object
  */
  MFHelmholtzElectrostaticEBBC&
  operator=(const MFHelmholtzElectrostaticEBBC&& a_other) = delete;

  /*!
    @brief Set stencil order
    @details EBHelmholtzEBBC uses least-squares based stencils for gradient reconstruction. This routine sets the
    order of that reconstruction. 
    @param[in] a_order. Must be > 0. 
    @note The radius of the stencil is the same as the order -- do not use an order larger than the number of ghost cells that you will use. 
  */
  void
  setOrder(const int a_order);

  /*!
    @brief Set stencil weight
    @param[in] a_weight Equation weights. 
    @details EBHelmholtzEBBC uses least-squares based stencils for gradient reconstruction. This routine sets the equation weights in that reconstruction. 
  */
  void
  setWeight(const int a_weight);

  /*!
    @brief Drop BC order if domain size is equal or below this.
  */
  void
  setDomainDropOrder(const int a_domainSize);

protected:
  /*!
    @brief Stencil order
  */
  int m_order;

  /*!
    @brief Stencil weight
  */
  int m_weight;

  /*!
    @brief Special flag for dropping stencil order when domains become coarser than this.
  */
  int m_domainDropOrder;

  /*!
    @brief Electrostatic boundary conditions
  */
  ElectrostaticEbBc m_electrostaticBCs;

  /*!
    @brief User define function. 
  */
  virtual void
  defineSinglePhase() override;

  /*!
    @brief Apply flux single phase
    @details This iterates over single-phase cells only, multi-phase cells are automatically enforced (in MFHelmholtzEBBC) by the jump boundary conditions. 
    @param[inout] a_singlePhaseVofs   Iterator for single-phase cut-cells
    @param[inout] a_Lphi              Operator kappa*L(phi)
    @param[in]    a_phi               Phi, cell-centered.
    @param[in]    a_Bcoef             EB-centered B-coefficient
    @param[in]    a_dit               Data index    
    @param[in]    a_beta              Beta. Passed in from operator. 
    @param[in]    a_homogeneousPhysBC Homogeneous BC or not. 
  */
  virtual void
  applyEBFluxRelaxSinglePhase(VoFIterator&           a_singlePhaseVofs,
                              EBCellFAB&             a_Lphi,
                              const EBCellFAB&       a_phi,
                              const BaseIVFAB<Real>& a_Bcoef,
                              const DataIndex&       a_dit,
                              const Real&            a_beta,
                              const bool&            a_homogeneousPhysBC) const override;

  /*!
    @brief Get electrode potential
    @param[in] a_pos Position in space
  */
  inline Real
  getElectrodePotential(const RealVect& a_pos) const;
};

#include <CD_NamespaceFooter.H>

#include <CD_MFHelmholtzElectrostaticEBBCImplem.H>

#endif
