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
    @param[in] a_eblg  Grid where this object is defined. 
    @param[in] a_Bcoef B-coefficient. You will need this for the flux stencil. 
    @param[in] a_dx    Grid resolution. 
  */
  void define(const EBLevelGrid& a_eblg, const RefCountedPtr<LevelData<BaseIVFAB<Real> > >& a_Bcoef, const RealVect& m_probLo, const Real& a_dx);

  /*!
    @brief Apply the EB flux. This is the version that is called by EBHelmholtzOp.
    @param[in] a_vofit          Iterator for cut-cells
    @param[in] a_Lphi           Operator kappa*L(phi)
    @param[in] a_phi            Phi, cell-centered. 
    @param[in] a_factor         Multiplication factor. 
    @param[in] a_useHomogeneous Use homogeneous BCs or not. 
  */
  void applyEBFlux(VoFIterator&       a_vofit,
		   EBCellFAB&         a_Lphi,
		   const EBCellFAB&   a_phi,
		   const Real&        a_factor,
		   const bool&        a_useHomogeneous) const;

  /*!
    @brief Get the stencil for computing the finite-volume approximation to kappa*Div(F). 
    @details The finite volume approximation to kappa*Div(F) is sum(Fluxes)/dx. This routine should return the flux on the EB, i.e. (B*dphi/dn)/dx. The derivative
             can be written dphi/dn = w0*phi0 + sum(wi*phi(i)) where phi0 is a boundary weight and phi(i) are grid points. This routine should return the stencil
	     B*wi*phi(i).
    @note    Multiplication by beta is done in the operator.
  */
  const LayoutData<BaseIVFAB<VoFStencil> >& getKappaDivFStencil() const;

protected:

  /*!
    @brief Grid resolution
  */
  Real m_dx;

  /*!
    @brief Lower-left corner of computational domain
  */
  RealVect m_probLo;

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
  LayoutData<BaseIVFAB<VoFStencil> > m_kappaDivFStencil;

  /*!
    @brief User define function. 
  */
  virtual void define() = 0;

  /*!
    @brief Apply the EB flux in a VoF. 
  */
  virtual void applyEBFlux(EBCellFAB&         a_Lphi,
			   const EBCellFAB&   a_phi,
			   const VolIndex&    a_vof,
			   const Real&        a_factor,
			   const bool&        a_useHomogeneous) const = 0;

  /*!
    @brief Returns physical position at the boundary
  */
  RealVect position(VolIndex& a_vof, const EBCellFAB& a_phi) const;
};

#include <CD_NamespaceFooter.H>

#endif
