/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFHelmholtzElectrostaticEBBCFactory.H
  @brief  Declaration of an electrostatic boundary condition factory class for MFHelmholtzOp.
  @author Robert Marskar
*/

#ifndef CD_MFHelmholtzElectrostaticEBBCFactory_H
#define CD_MFHelmholtzElectrostaticEBBCFactory_H

// Our includes
#include <CD_MFHelmholtzEBBCFactory.H>
#include <CD_ElectrostaticEbBc.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for making electrostatic boundary conditions on EBs in an MFHelmholtzOp context. 
  @details Users will pass an object of this class to MFHelmholtzOpFactory. The input arguments in the constructor will set the 
  stencil order and weights for the least squares reconstruction that is used in EBHelmholtzEBBC. The actual BC values are passed
  in through the a_electrostaticBCs argument in the constructor. 
*/
class MFHelmholtzElectrostaticEBBCFactory : public MFHelmholtzEBBCFactory {
public:

  /*!
    @brief Weak constructor is bad, use the full constructor. 
  */
  MFHelmholtzElectrostaticEBBCFactory() = delete;

  /*!
    @brief Full constructor
    @param[in] a_order            Stencil order for least squares reconstruction
    @param[in] a_weight           Equation weights for least squares reconstruction
    @param[in] a_electrostaticBCs Electrostatic boundary conditions. 
    @note EBHelmholtzEBBC will use a stencil radius equal to the order -- do not use a stencil order higher than your number of ghost cells. 
  */
  MFHelmholtzElectrostaticEBBCFactory(const int a_order, const int a_weight, const ElectrostaticEbBc& a_electrostaticBCs);

  /*!
    @brief No copy construction
    @param[in] a_other Other factory
  */
  MFHelmholtzElectrostaticEBBCFactory(MFHelmholtzElectrostaticEBBCFactory& a_other) = delete;

  /*!
    @brief No move construction
    @param[in] a_other Other factory
  */
  MFHelmholtzElectrostaticEBBCFactory(MFHelmholtzElectrostaticEBBCFactory&& a_other) = delete;

  /*!
    @brief No copy assignment
    @param[in] a_other Other factory
  */
  MFHelmholtzElectrostaticEBBCFactory& operator=(const MFHelmholtzElectrostaticEBBCFactory&) = delete;

  /*!
    @brief No move assignment
    @param[in] a_other Other factory
  */
  MFHelmholtzElectrostaticEBBCFactory& operator=(const MFHelmholtzElectrostaticEBBCFactory&&) = delete;

  /*!
    @brief Default constructor (does nothing)
  */
  virtual ~MFHelmholtzElectrostaticEBBCFactory();
  
  /*!
    @brief Factory method
    @param[in] a_iphase Phase (0=gas and 1=solid). 
    @param[in] a_jumpBC Jump bc object. 
  */
  RefCountedPtr<EBHelmholtzEBBC> create(const int a_iphase, const RefCountedPtr<MFHelmholtzJumpBC>& a_jumpBC) const override;

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
};

#include <CD_NamespaceFooter.H>

#endif
