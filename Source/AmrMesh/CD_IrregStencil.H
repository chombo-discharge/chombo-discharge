/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   IrregStencil.H
  @brief  Abstract class for building stencils on cut-cells. 
  @author Robert Marskar
*/

#ifndef CD_IrregStencil_H
#define CD_IrregStencil_H

// Chombo includes
#include <DisjointBoxLayout.H>
#include <EBISLayout.H>
#include <BaseIVFAB.H>
#include <RefCountedPtr.H>
#include <Stencils.H>
#include <LayoutData.H>
#include <VoFIterator.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for holding stencils on irregular cells over a single AMR level
*/
class IrregStencil {
public:

  /*!
    @brief Enum for identifying stencil -- only meant for enhancing code visibility. 
  */
  enum class StencilType {
    NoType,
    Linear,
    TaylorExtrapolation,
    LeastSquares,
    PiecewiseLinear
  };

  /*!
    @brief Constructor, but weak construction is not allowed. 
  */
  IrregStencil();

  /*!
    @brief Full constructor
    @param[in] a_dbl    Grids
    @param[in] a_ebisl  EBIS layout
    @param[in] a_domain Problem domain
    @param[in] a_dx     Resolutions
    @param[in] a_order  Interpolation order
    @param[in] a_radius Radius for least squares
    @param[in] a_type   Stencil type
  */
  IrregStencil(const DisjointBoxLayout&        a_dbl,
	       const EBISLayout&               a_ebisl,
	       const ProblemDomain&            a_domain,
	       const Real&                     a_dx,
	       const int                       a_order,
	       const int                       a_radius,
	       const IrregStencil::StencilType a_type);

  /*!
    @brief Destructor
  */
  ~IrregStencil();



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

  /*!
    @brief Get the VoFIterator
  */
  VoFIterator& getVofIterator(const DataIndex& a_dit) const {
    return m_vofIter[a_dit];
  }

protected:

  /*!
    @brief Stencil
  */
  LayoutData<RefCountedPtr<BaseIVFAB<VoFStencil> > > m_stencils;

  /*!
    @brief VoFIterators
  */
  mutable LayoutData<VoFIterator> m_vofIter;

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
  IrregStencil::StencilType m_stencilType;

  /*!
    @brief Domain
  */
  ProblemDomain m_domain;

  /*!
    @brief Define function
  */
  virtual void define(const DisjointBoxLayout&        a_dbl,
		      const EBISLayout&               a_ebisl,
		      const ProblemDomain&            a_domain,
		      const Real&                     a_dx,
		      const int                       a_order,
		      const int                       a_radius,
		      const IrregStencil::StencilType a_type);
  
  /*!
    @brief Build the desired stencil
  */
  virtual void buildStencil(VoFStencil&              a_sten,
			    const VolIndex&          a_vof,
			    const DisjointBoxLayout& a_dbl,
			    const ProblemDomain&     a_domain,
			    const EBISBox&           a_ebisbox,
			    const Box&               a_box,
			    const Real&              a_dx,
			    const IntVectSet&        a_cfivs) = 0;
};

#include <CD_NamespaceFooter.H>

#endif
