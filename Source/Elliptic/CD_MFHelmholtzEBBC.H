/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFHelmholtzEBBC.H
  @brief  Declaration of an EB boundary condition class for MFHelmholtzOp.
  @author Robert Marskar
  @todo   Consider using an interface here where we do matching in multiphase cells here in the base class (rather than giving the responsibility to derived classes)
*/

#ifndef CD_MFHelmholtzEBBC_H
#define CD_MFHelmholtzEBBC_H

// Our includes
#include <CD_EBHelmholtzEBBC.H>
#include <CD_JumpBC.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for making boundary conditions on the EB in a multifluid context. 
*/
class MFHelmholtzEBBC : public EBHelmholtzEBBC {
public:

  /*!
    @brief Disallowed weak constructor. Factory must call the other version. 
  */
  MFHelmholtzEBBC() = delete;

  /*!
    @brief Disallowed copy constructor. Factory must call the other version. 
  */
  MFHelmholtzEBBC(const MFHelmholtzEBBC& a_other) = delete;

  /*!
    @brief Full constructor
  */
  MFHelmholtzEBBC(const int a_phase, const RefCountedPtr<JumpBC>& a_jumpBC);

  /*!
    @brief Default constructor
  */
  ~MFHelmholtzEBBC();
  
  /*!
    @brief Apply the EB flux. This is the version that is called by EBHelmholtzOp.
    @param[inout] a_vofit             Iterator for cut-cells
    @param[inout] a_Lphi              Operator kappa*L(phi)
    @param[in]    a_phi               Phi, cell-centered. 
    @param[in]    a_dit               Data index    
    @param[in]    a_beta              Beta. Passed in from operator. 
    @param[in]    a_homogeneousPhysBC Homogeneous BC or not. 
  */
  virtual void applyEBFlux(VoFIterator&       a_vofit,
			   EBCellFAB&         a_Lphi,
			   const EBCellFAB&   a_phi,
			   const DataIndex&   a_dit,
			   const Real&        a_beta,
			   const bool&        a_homogeneousPhysBC) const override final;

protected:

  /*!
    @brief Jump bc object
  */
  RefCountedPtr<JumpBC> m_jumpBC;

  /*!
    @brief Phase
  */
  int m_phase;

  /*!
    @brief User define function. 
  */
  virtual void define() override = 0;

  /*!
    @brief Flux application function.
    @param[inout] a_singlePhaseVofs   Iterator for single-phase cut-cells
    @param[inout] a_multiPhaseVofs    Iterator for multi-phase cut-cells
    @param[inout] a_Lphi              Operator kappa*L(phi)
    @param[in]    a_phi               Phi, cell-centered. 
    @param[in]    a_dit               Data index    
    @param[in]    a_beta              Beta. Passed in from operator. 
    @param[in]    a_homogeneousPhysBC Homogeneous BC or not. 
  */
  virtual void applyEBFlux(VoFIterator&       a_singlePhaseVofs,
			   VoFIterator&       a_multiPhaseVofs,
			   EBCellFAB&         a_Lphi,
			   const EBCellFAB&   a_phi,
			   const DataIndex&   a_dit,
			   const Real&        a_beta,
			   const bool&        a_homogeneousPhysBC) const = 0;
};

#include <CD_NamespaceFooter.H>

#endif
