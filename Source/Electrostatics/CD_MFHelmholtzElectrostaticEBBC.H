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
*/
class MFHelmholtzElectrostaticEBBC : public MFHelmholtzEBBC {
public:

  /*!
    @brief Full constructor
  */
  MFHelmholtzElectrostaticEBBC(const int a_phase, const ElectrostaticEbBc& a_electrostaticBCs, const RefCountedPtr<JumpBC>& a_jumpBC);

  /*!
    @brief No weak construction
  */
  MFHelmholtzElectrostaticEBBC() = delete;

  /*!
    @brief No copy construction
  */
  MFHelmholtzElectrostaticEBBC(const MFHelmholtzElectrostaticEBBC& a_other) = delete;

  /*!
    @brief No move construction
  */
  MFHelmholtzElectrostaticEBBC(const MFHelmholtzElectrostaticEBBC&& a_other) = delete;

  /*!
    @brief Default constructor
  */
  ~MFHelmholtzElectrostaticEBBC();

  /*!
    @brief No copy assignment
  */
  MFHelmholtzElectrostaticEBBC& operator=(const MFHelmholtzElectrostaticEBBC& a_other) = delete;

  /*!
    @brief No move assignment
  */
  MFHelmholtzElectrostaticEBBC& operator=(const MFHelmholtzElectrostaticEBBC&& a_other) = delete;

  /*!
    @brief Set stencil order
  */
  void setOrder(const int a_order);

  /*!
    @brief Set stencil weight
  */
  void setWeight(const int a_weight);

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
    @brief Electrostatic boundary conditions
  */
  ElectrostaticEbBc m_electrostaticBCs;

  /*!
    @brief User define function. 
  */
  virtual void defineSinglePhase() override;

  /*!
    @brief Apply flux single phase
    @param[inout] a_singlePhaseVofs   Iterator for single-phase cut-cells
    @param[inout] a_Lphi              Operator kappa*L(phi)
    @param[in]    a_phi               Phi, cell-centered. 
    @param[in]    a_dit               Data index    
    @param[in]    a_beta              Beta. Passed in from operator. 
    @param[in]    a_homogeneousPhysBC Homogeneous BC or not. 
  */
  virtual void applyEBFluxSinglePhase(VoFIterator&       a_singlePhaseVofs,
				      EBCellFAB&         a_Lphi,
				      const EBCellFAB&   a_phi,
				      const DataIndex&   a_dit,
				      const Real&        a_beta,
				      const bool&        a_homogeneousPhysBC) const override;

  /*!
    @brief Get electrode potential
    @param[in] a_pos Position in space
  */
  inline
  Real getElectrodePotential(const RealVect& a_pos) const;

};

#include <CD_NamespaceFooter.H>

#include <CD_MFHelmholtzElectrostaticEBBCImplem.H>

#endif
