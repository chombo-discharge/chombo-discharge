
/*!
  @file irreg_stencil.H
  @brief Abstract class for holding stencils on irregular cells over a single AMR level
  @author Robert Marskar
  @date Nov. . 2017
*/

#ifndef _IRREG_STENCIL_
#define _IRREG_STENCIL_

#include <DisjointBoxLayout.H>
#include <EBISLayout.H>
#include <BaseIVFAB.H>
#include <RefCountedPtr.H>
#include <Stencils.H>
#include <LayoutData.H>

#include "CD_NamespaceHeader.H"

/*!
  @brief Class for holding stencils on irregular cells over a single AMR level
*/
class irreg_stencil {
public:

  enum class StencilType {
    Linear,
    TaylorExtrapolation,
    LeastSquarse,
    PiecewiseLinear
  };

  /*!
    @brief Constructor. You must subseqently call define. No strong construction for this class since buildStencils is pure virtual
  */
  irreg_stencil();

  /*!
    @brief Full constructor
  */
  irreg_stencil(const DisjointBoxLayout&       a_dbl,
		const EBISLayout&              a_ebisl,
		const ProblemDomain&           a_domain,
		const Real&                    a_dx,
		const int                      a_order,
		const int                      a_radius,
		const stencil_type             a_type);

  /*!
    @brief Destructor
  */
  ~irreg_stencil();

  /*!
    @brief Define function
  */
  virtual void define(const DisjointBoxLayout&       a_dbl,
		      const EBISLayout&              a_ebisl,
		      const ProblemDomain&           a_domain,
		      const Real&                    a_dx,
		      const int                      a_order,
		      const int                      a_radius,
		      const stencil_type             a_type);

  /*!
    @brief Get the stencils over a single box
  */
  const BaseIVFAB<VoFStencil>& operator[](const DataIndex& a_dit) const {
    return *m_stencils[a_dit];
  }
  
  /*!
    @brief Get the stencils over a single box
  */
  BaseIVFAB<VoFStencil>& operator[](const DataIndex& a_dit) {
    return *m_stencils[a_dit];
  }

protected:

  LayoutData<RefCountedPtr<BaseIVFAB<VoFStencil> > > m_stencils;

  /*!
    @brief Stencils defined or not
  */
  bool m_defined;;

  /*!
    @brief Grids
  */
  DisjointBoxLayout m_dbl;

  /*!
    @brief EBIS layout
  */
  EBISLayout m_ebisl;

  /*!
    @brief Level resolution
  */
  Real m_dx;

  /*!
    @brief Stencil radius
  */
  int m_radius;

  /*!
    @brief Stencil order
  */
  int m_order;

  /*!
    @brief Stencil type
  */
  stencil_type m_stencilType;

  /*!
    @brief Domain
  */
  ProblemDomain m_domain;
  
  /*!
    @brief Build the desired stencil
  */
  virtual void build_stencil(VoFStencil&              a_sten,
			     const VolIndex&          a_vof,
			     const DisjointBoxLayout& a_dbl,
			     const ProblemDomain&     a_domain,
			     const EBISBox&           a_ebisbox,
			     const Box&               a_box,
			     const Real&              a_dx,
			     const IntVectSet&        a_cfivs) = 0;
};
#include "CD_NamespaceFooter.H"
#endif
