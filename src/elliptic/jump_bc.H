/*!
  @file jump_bc.H
  @brief Declaration of a multifluid boundary condition class for the poisson equation
  @author Robert Marskar
  @date Dec. 2017
*/

#ifndef _JUMP_BC_
#define _JUMP_BC_

#include "mfis.H"
#include "MFLevelGrid.H"
#include "MFInterfaceFAB.H"
#include "MFBaseIVFAB.H"

#include <Stencils.H>
#include <EBISBox.H>

/*!
  @brief Class matching the surface potential with multifluid poisson
*/
class jump_bc {
public:

  /*!
    @brief Empty constructor, must call define
  */
  jump_bc();

  /*!
    @brief Full constructor
  */
  jump_bc(const MFLevelGrid&            a_mflg,
	  const LevelData<MFBaseIVFAB>& a_bco,
	  const Real&                   a_dx,
	  const int                     a_order,
	  const LayoutData<IntVectSet>* a_cfivs);

  /*!
    @brief Destructor
  */
  virtual ~jump_bc();

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
  virtual void match_bc(LevelData<BaseIVFAB<Real> >&       a_phibc,
			const LevelData<BaseIVFAB<Real> >& a_jump,
			const LevelData<MFCellFAB>&        a_phi,
			const bool                         a_homogeneous);

  /*!
    @brief Compute dphidn for all phases
  */
  virtual void compute_dphidn(Vector<LevelData<BaseIVFAB<Real> > >&        a_dphidn,
			      const Vector<LevelData<BaseIVFAB<Real> > >&  a_phibc,
			      const LevelData<MFCellFAB>&                  a_phi);
  /*!
    @brief Compute dphidn for specified phase usign precomputed phibc
  */
  virtual void compute_dphidn(LevelData<BaseIVFAB<Real> >&       a_dphidn,
			      const LevelData<BaseIVFAB<Real> >& a_phibc,
			      const LevelData<MFCellFAB>&        a_phi,
			      const int                          a_phase);
  
  /*!
    @brief Compute dphidn for specified phase usign precomputed phibc
  */
  virtual void compute_dphidn(LevelData<BaseIVFAB<Real> >&       a_dphidn,
			      const LevelData<BaseIVFAB<Real> >& a_phibc,
			      const LevelData<EBCellFAB>&        a_phi,
			      const int                          a_phase);

  /*!
    @brief Compute dphidn for specified phase usign precomputed phibc
  */
  virtual void compute_dphidn(BaseIVFAB<Real>&       a_dphidn,
			      const BaseIVFAB<Real>& a_phibc,
			      const EBCellFAB&       a_phi,
			      const DataIndex&       a_dit,
			      const int              a_phase);

  /*!
    @brief Get stencils
  */
  virtual LayoutData<MFInterfaceFAB<VoFStencil> >& get_stencils();

  /*!
    @brief Get stencils
  */
  virtual LayoutData<MFInterfaceFAB<VoFStencil> >& get_avgStencils();

  /*!
    @brief Get weights
  */
  virtual LayoutData<MFInterfaceFAB<Real> >& get_weights();

    /*!
    @brief Get weights
  */
  virtual LayoutData<MFInterfaceFAB<Real> >& get_avgWeights();

  /*!
    @brief Get b-coefficient on each side
  */
  virtual LayoutData<MFInterfaceFAB<Real> >& get_bco();

    /*!
    @brief Get b-coefficient on each side
  */
  virtual LayoutData<MFInterfaceFAB<Real> >& get_avgBco();

  /*!
    @brief Get inhomogeneous contribution on each side
  */
  virtual LayoutData<MFInterfaceFAB<Real> >& get_inhomo();

  /*!
    @brief Homogeneous contribution
  */
  virtual LayoutData<MFInterfaceFAB<Real> >& get_homog();


protected:

  RefCountedPtr<mfis> m_mfis;

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
  virtual void define_vofiter();

  /*!
    @brief Set b-coefficient
  */
  virtual void set_bco(const LevelData<MFBaseIVFAB>& a_bco);

  /*!
    @brief Build stencils for matching
  */
  virtual void build_stencils();

  /*!
    @brief Get a second order stencil 
  */
  virtual bool get_second_order_sten(Real&             a_weight,
				     VoFStencil&       a_stencil,
				     const VolIndex&   a_vof,
				     const EBISBox&    a_ebisbox,
				     const IntVectSet& a_cfivs);

  /*!
    @brief Get a first order stencil 
  */
  virtual void get_first_order_sten(Real&             a_weight,
				    VoFStencil&       a_stencil,
				    const VolIndex&   a_vof,
				    const EBISBox&    a_ebisbox,
				    const IntVectSet& a_cfivs);

  /*!
    @brief Match the boundary condition
  */
  inline void match_bc(BaseIVFAB<Real>&                  a_phibc,
		       const BaseIVFAB<Real>&            a_jump,
		       const MFCellFAB&                  a_phi,
		       const DataIndex&                  a_dit);

  /*!
    @brief Compute the average jump condition
  */
  inline void compute_avg_jump(const BaseIVFAB<Real>& a_jump, const MFCellFAB& a_phi, const DataIndex& a_dit);
};

#endif
