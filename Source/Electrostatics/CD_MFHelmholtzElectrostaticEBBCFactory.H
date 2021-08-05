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
  @brief Class for making electrostatic boundary conditions in an MFHelmholtzOp context
*/
class MFHelmholtzElectrostaticEBBCFactory : public MFHelmholtzEBBCFactory {
public:

  /*!
    @brief Weak constructor is bad, use the full constructor. 
  */
  MFHelmholtzElectrostaticEBBCFactory() = delete;

  /*!
    @brief Full constructor
    @param[in] a_electrostaticBCs Electrostatic boundary conditions
  */
  MFHelmholtzElectrostaticEBBCFactory(const int a_order, const int a_weight, const ElectrostaticEbBc& a_electrostaticBCs);

  /*!
    @brief No copy construction
  */
  MFHelmholtzElectrostaticEBBCFactory(MFHelmholtzElectrostaticEBBCFactory& a_other) = delete;

  /*!
    @brief No move construction
  */
  MFHelmholtzElectrostaticEBBCFactory(MFHelmholtzElectrostaticEBBCFactory&& a_other) = delete;

  /*!
    @brief No copy assignemnt
  */
  MFHelmholtzElectrostaticEBBCFactory& operator=(const MFHelmholtzElectrostaticEBBCFactory&) = delete;

  /*!
    @brief No move assignement
  */
  MFHelmholtzElectrostaticEBBCFactory& operator=(const MFHelmholtzElectrostaticEBBCFactory&&) = delete;

  /*!
    @brief Default constructor
  */
  ~MFHelmholtzElectrostaticEBBCFactory();
  
  /*!
    @brief Factory method
    @param[in] a_jumpBC Jump bc object. 
  */
  RefCountedPtr<EBHelmholtzEBBC> create(const int a_iphase, const RefCountedPtr<JumpBC>& a_jumpBC) const override;

protected:

  /*!
    @brief Electrostatic boundary conditions
  */
  ElectrostaticEbBc m_electrostaticBCs;
};

#include <CD_NamespaceFooter.H>

#endif
