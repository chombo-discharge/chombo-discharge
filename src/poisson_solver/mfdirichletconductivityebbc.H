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
			      const IntVect*       a_ghost_rhs) : DirichletConductivityEBBC(a_domain,
											    a_ebisl,
											    a_dx,
											    a_ghost_phi,
											    a_ghost_rhs){

  }
  
  /*!
    @brief Destructor
  */
  virtual ~mfdirichletconductivityebbc();

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


protected:

  /*!
    @brief Interface cells
  */
  LayoutData<IntVectSet> m_ivs;

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
