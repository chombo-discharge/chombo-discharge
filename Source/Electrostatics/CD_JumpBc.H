/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   JumpBc.H
  @brief  Declaration of a multifluid boundary condition class for the Poisson equation
  @author Robert Marskar
*/

#ifndef CD_JumpBc_H
#define CD_JumpBc_H

// Chombo includes
#include <MFCellFAB.H>
#include <Stencils.H>
#include <EBISBox.H>

// Our includes
#include <CD_MultiFluidIndexSpace.H>
#include <CD_MFLevelGrid.H>
#include <CD_MFInterfaceFAB.H>
#include <CD_MFBaseIVFAB.H>
#include <CD_NamespaceHeader.H>
  
/*!
  @brief Class matching the surface potential with multifluid poisson
*/
class JumpBc {
public:

  /*!
    @brief Empty constructor, must call define
  */
  JumpBc();

  /*!
    @brief Full constructor
  */
  JumpBc(const MFLevelGrid&            a_mflg,
	 const LevelData<MFBaseIVFAB>& a_bco,
	 const Real&                   a_dx,
	 const int                     a_order,
	 const LayoutData<IntVectSet>* a_cfivs);

  /*!
    @brief Destructor
  */
  virtual ~JumpBc();

  /*!
    @brief Switch between quadrant-based and non-quadrant based least squares stencils
  */
  static bool s_quadrant_based;

  /*!
    @brief Radius for least squares stencils
  */
  static int s_lsq_radius;

  /*!
    @brief Define function
  */
  virtual void define(const MFLevelGrid&            a_mflg,
		      const LevelData<MFBaseIVFAB>& a_bco,
		      const Real&                   a_dx,
		      const int                     a_order,
		      const LayoutData<IntVectSet>* a_cfivs);

  /*!
    @brief Match the boundary condition
  */
  virtual void matchBc(LevelData<BaseIVFAB<Real> >&       a_phibc,
		       const LevelData<BaseIVFAB<Real> >& a_jump,
		       const LevelData<MFCellFAB>&        a_phi,
		       const bool                         a_homogeneous);

  /*!
    @brief Match the boundary condition
  */
  virtual void matchBc(LevelData<BaseIVFAB<Real> >&       a_phibc,
		       const LevelData<BaseIVFAB<Real> >& a_jump,
		       const LevelData<MFCellFAB>&        a_phi,
		       DataIterator&                      a_dit,
		       const bool                         a_homogeneous);

  /*!
    @brief Compute dphidn for all phases
  */
  virtual void computeDnPhi(Vector<LevelData<BaseIVFAB<Real> > >&        a_dphidn,
			    const Vector<LevelData<BaseIVFAB<Real> > >&  a_phibc,
			    const LevelData<MFCellFAB>&                  a_phi);
  /*!
    @brief Compute dphidn for specified phase usign precomputed phibc
  */
  virtual void computeDnPhi(LevelData<BaseIVFAB<Real> >&       a_dphidn,
			    const LevelData<BaseIVFAB<Real> >& a_phibc,
			    const LevelData<MFCellFAB>&        a_phi,
			    const int                          a_phase);
  
  /*!
    @brief Compute dphidn for specified phase usign precomputed phibc
  */
  virtual void computeDnPhi(LevelData<BaseIVFAB<Real> >&       a_dphidn,
			    const LevelData<BaseIVFAB<Real> >& a_phibc,
			    const LevelData<EBCellFAB>&        a_phi,
			    const int                          a_phase);

  /*!
    @brief Compute dphidn for specified phase usign precomputed phibc
  */
  virtual void computeDnPhi(BaseIVFAB<Real>&       a_dphidn,
			    const BaseIVFAB<Real>& a_phibc,
			    const EBCellFAB&       a_phi,
			    const DataIndex&       a_dit,
			    const int              a_phase);

  /*!
    @brief Get stencils
  */
  virtual LayoutData<MFInterfaceFAB<VoFStencil> >& getStencils();

  /*!
    @brief Get stencils
  */
  virtual LayoutData<MFInterfaceFAB<VoFStencil> >& getAvgStencils();

  /*!
    @brief Get weights
  */
  virtual LayoutData<MFInterfaceFAB<Real> >& getWeights();

  /*!
    @brief Get weights
  */
  virtual LayoutData<MFInterfaceFAB<Real> >& getAvgWeights();

  /*!
    @brief Get b-coefficient on each side
  */
  virtual LayoutData<MFInterfaceFAB<Real> >& getAvgBco();

  /*!
    @brief Get b-coefficient on each side
  */
  virtual LayoutData<MFInterfaceFAB<Real> >& get_avgBco();

  /*!
    @brief Get inhomogeneous contribution on each side
  */
  virtual LayoutData<MFInterfaceFAB<Real> >& getInhomogeneous();

  /*!
    @brief Homogeneous contribution
  */
  virtual LayoutData<MFInterfaceFAB<Real> >& getHomogeneous();


protected:

  RefCountedPtr<MultiFluidIndexSpace> m_multifluidIndexSpace;

  /*!
    @brief Level grid
  */
  MFLevelGrid m_mflg;

  /*!
    @brief Grid
  */
  DisjointBoxLayout m_grids;

  /*!
    @brief Problem domain
  */
  ProblemDomain m_domain;

  /*!
    @brief Stencil order
  */
  int m_order;

  /*!
    @brief Resolution at level
  */
  Real m_dx;

  /*!
    @brief Is defined or not
  */
  bool m_defined;

  /*!
    @brief Coarse-fine interface
  */
  const LayoutData<IntVectSet>* m_cfivs;

  /*!
    @brief Cells on which a jump condition is defined
  */
  LayoutData<IntVectSet> m_ivs;

  /*!
    @brief b-coefficient on each side
  */
  LayoutData<MFInterfaceFAB<Real> > m_bco;

  /*!
    @brief Weights on each side
  */
  LayoutData<MFInterfaceFAB<Real> > m_weights;

  /*!
    @brief Inhomogeneous contribution to bc
  */
  LayoutData<MFInterfaceFAB<Real> > m_inhomo;

  /*!
    @brief Homogeneous contribution to bc
  */
  LayoutData<MFInterfaceFAB<Real> > m_homog;

  /*!
    @brief Stencils on each side
  */
  LayoutData<MFInterfaceFAB<VoFStencil> > m_stencils;

  /*!
    @brief Average derivative stencil
  */
  LayoutData<MFInterfaceFAB<VoFStencil> > m_avgStencils;
  LayoutData<MFInterfaceFAB<Real> >       m_avgWeights;
  LayoutData<MFInterfaceFAB<Real> >       m_avgBco;
  LayoutData<MFInterfaceFAB<Real> >       m_avgJump;
  LayoutData<MFInterfaceFAB<Real> >       m_avgFactor;

  /*!
    @brief VoFIterator for gas-side cells
  */
  LayoutData<VoFIterator> m_vofit_gas;

  /*!
    @brief VoFIterator for solid-side cells
  */
  LayoutData<VoFIterator> m_vofit_sol;

  /*!
    @brief VoFIterator for bcs
  */
  LayoutData<VoFIterator> m_vofit_bc;

  /*!
    @brief Define VoFIterators
  */
  virtual void defineVofIterator();

  /*!
    @brief Set b-coefficient
  */
  virtual void setBcoefficient(const LevelData<MFBaseIVFAB>& a_bco);

  /*!
    @brief Build stencils for matching
  */
  virtual void buildStencils();

  /*!
    @brief Get a second order stencil 
  */
  virtual bool getSecondOrderStencil(Real&             a_weight,
				     VoFStencil&       a_stencil,
				     const VolIndex&   a_vof,
				     const EBISBox&    a_ebisbox,
				     const IntVectSet& a_cfivs);

  /*!
    @brief Get a first order stencil 
  */
  virtual void getFirstOrderStencil(Real&             a_weight,
				    VoFStencil&       a_stencil,
				    const VolIndex&   a_vof,
				    const EBISBox&    a_ebisbox,
				    const IntVectSet& a_cfivs);

  /*!
    @brief Match the boundary condition
  */
  inline
  void matchBc(BaseIVFAB<Real>&                  a_phibc,
	       const BaseIVFAB<Real>&            a_jump,
	       const MFCellFAB&                  a_phi,
	       const DataIndex&                  a_dit);

  /*!
    @brief Compute the average jump condition
  */
  inline
  void computeAvgJump(const BaseIVFAB<Real>& a_jump, const MFCellFAB& a_phi, const DataIndex& a_dit);
};

#include <CD_NamespaceFooter.H>

#include <CD_JumpBcImplem.H>

#endif
