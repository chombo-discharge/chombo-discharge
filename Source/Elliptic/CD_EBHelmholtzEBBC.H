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
    @brief Default constructor. 
  */
  EBHelmholtzEBBC();

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
    @brief Define function that is called by EBHelmholtzOp
  */
  void define(const EBLevelGrid& a_eblg, const RefCountedPtr<LevelData<BaseIVFAB<Real> > >& a_Bcoef, const Real a_factor);

  /*!
    @brief Apply the EB flux. This is the version that is called by EBHelmholtzOp.
  */
  void applyEBFlux(VoFIterator&       a_vofit,
			   EBCellFAB&         a_Lphi,
			   const EBCellFAB&   a_phi,
			   const RealVect&    a_probLo,
			   const RealVect&    a_dx,
			   const Real&        a_factor,
			   const bool&        a_useHomogeneous) const;

  /*!
    @brief Get the flux stencil. 
    @note Maybe this should be virtual....
  */
  const LayoutData<BaseIVFAB<VoFStencil> >& getFluxStencil() const;

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
  LayoutData<BaseIVFAB<VoFStencil> > m_fluxStencil;

  /*!
    @brief User define function. 
  */
  virtual void define(const EBLevelGrid& a_eblg, const Real a_factor) = 0;

  /*!
    @brief Apply the EB flux in a VoF. 
  */
  virtual void applyEBFlux(EBCellFAB&         a_Lphi,
			   const EBCellFAB&   a_phi,
			   const VolIndex&    a_vof,
			   const RealVect&    a_probLo,
			   const RealVect&    a_dx,
			   const Real&        a_factor,
			   const bool&        a_useHomogeneous) const = 0;  

};

#include <CD_NamespaceFooter.H>

#endif
