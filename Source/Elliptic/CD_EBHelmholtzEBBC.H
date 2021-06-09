/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzEBBC.H
  @brief  Declaration of a boundary condition class for EBHelmholtzOp
  @author Robert Marskar
*/

#ifndef CD_EBHelmholtzEBBC_H
#define CD_EBHelmholtzEBBC_H

// Chombo includes
#include <EBLevelGrid.H>
#include <RefCountedPtr.H>
#include <LevelData.H>
#include <BaseIVFAB.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Base class for passing EB boundary conditions into EBHelmholtzOp. 
*/
class EBHelmholtzEBBC {
public:

  /*!
    @brief Constructor. 
  */
  EBHelmholtzEBBC() = delete;

  /*!
    @brief Destructor
  */
  ~EBHelmholtzEBBC();

  /*!
    @brief Disallowed -- don't see why you would need it.
  */
  EBHelmholtzEBBC(const EBHelmholtzEBBC& a_other) = delete;

  /*!
    @brief Disallowed - don't see why you would need it.
  */
  EBHelmholtzEBBC& operator=(const EBHelmholtzEBBC& a_other) = delete;

  /*!
    @brief Define function
  */
  virtual void define(const LayoutData<IntVectSet>& a_cfivs, const Real a_factor) = 0;

  /*!
    @brief Set the B-coefficient. This needs to be a separate function because the factories will not know about the coefficients. 
  */
  virtual void setBcoef(EBLevelGrid& a_eblg, const RefCountedPtr<LevelData<BaseIVFAB<Real> > >& a_Bcoef);

protected:

  /*!
    @brief Level grid
  */
  EBLevelGrid m_eblg;

  /*!
    @brief B-coefficent. 
  */
  RefCountedPtr<LevelData<BaseIVFAB<Real> > > m_Bcoef;

  /*!
    @brief Stencils for computing the flux
  */
  LayoutData<BaseIVFAB<VoFStencil> > m_stencils;

};

#include <CD_NamespaceFooter.H>

#endif
