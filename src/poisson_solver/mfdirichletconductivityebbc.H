/*!
  @file mfdirichletconductivityebbc.H
  @brief Declaration of a boundary condition class for doing convential Dirichlet BCs as well as multifluid matching bcs 
  on selected cells
  @author Robert Marskar
  @date Dec. 2017
*/

#ifndef _MFDIRICHLETCONDUCTIVITYEBBC_
#define _MFDIRICHLETCONDUCTIVITYEBBC_

#include "MFLevelGrid.H"
#include "mfis.H"

#include <DirichletConductivityEBBC.H>

/*!
  @brief Class for doing regular Dirichlet boundary conditions on EBConductivityOps with selected multifluid matching 
  boundary conditions.
  @details This is a simple class that essentially just ignores calls to homogeneous BCs on interface cells. 
*/
class mfdirichletconductivityebbc : public DirichletConductivityEBBC {
public:
  /*!
    @brief Constructor. Implements DirichletConductivityEBBC
  */
  mfdirichletconductivityebbc(const ProblemDomain& a_domain,
			      const EBISLayout&    a_ebisl,
			      const RealVect&      a_dx,
			      const IntVect*       a_ghost_phi,
			      const IntVect*       a_ghost_rhs);
  
  /*!
    @brief Destructor
  */
  virtual ~mfdirichletconductivityebbc();

  /*!
    @brief Switch between quadrant-based and non-quadrant based least squares stencils
  */
  static bool s_quadrant_based;

  /*!
    @brief Radius for least squares stencils
  */
  static int s_lsq_radius;

  /*!
    @brief Define interface cells
  */
  virtual void define_ivs(const MFLevelGrid& a_mflg);

  /*!
    @brief Define stencils and stuff
  */
  virtual void define(const LayoutData<IntVectSet>& a_cfivs, const Real& a_factor);

#if 1
  /*!
    @brief Apply EB flux
  */
  virtual void applyEBFlux(EBCellFAB&                    a_lphi,
			   const EBCellFAB&              a_phi,
			   VoFIterator&                  a_vofit,
			   const LayoutData<IntVectSet>& a_cfivs,
			   const DataIndex&              a_dit,
			   const RealVect&               a_probLo,
			   const RealVect&               a_dx,
			   const Real&                   a_factor,
			   const bool&                   a_useHomogeneous,
			   const Real&                   a_time);
#endif

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
    @brief Set bc order
  */
  virtual void setOrder(int a_order);


#if 1 // As long as we overwrite this, we can not expect convergence unless EVERYTHING works
  /*!
    @brief Get stencils for computing the flux
  */
  virtual LayoutData<BaseIVFAB<VoFStencil> >* getFluxStencil(int ivar);
#endif


protected:

  /*!
    @brief Problem domain
  */
  ProblemDomain m_domain;
  
  /*!
    @brief EBIS layout
  */
  EBISLayout m_ebisl;

  /*!
    @brief Resolution
  */
  RealVect m_dx;

  /*!
    @brief Boundary condition order
  */
  int m_order;

  /*!
    @brief Interface cells
  */
  LayoutData<IntVectSet> m_ivs;

  /*!
    @brief Weights for irregular cells
  */
  LayoutData<BaseIVFAB<Real> > m_irreg_weights;

  /*!
    @brief Weights for matching cells
  */
  LayoutData<BaseIVFAB<Real> > m_matching_weights;

  /*!
    @brief Stencils on irregular cells
  */
  LayoutData<BaseIVFAB<VoFStencil> > m_irreg_stencils;

  /*!
    @brief Stencils on matching cells
  */
  LayoutData<BaseIVFAB<VoFStencil> > m_matching_stencils;

  /*!
    @brief Flux to be smooshed into applyEBFlux for multi-fluid cells
  */
  LevelData<BaseIVFAB<Real> >* m_dphidn;

  /*!
    @brief VoFIterator for boundary condition matching cells
  */
  VoFIterator m_vofiter_matching;

  /*!
    @brief Make sure multifluid cells have been set
  */
  bool m_definedmf;
};
#endif
