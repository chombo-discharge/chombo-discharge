/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_IrregStencil.H
  @brief  Abstract class for holding cut-cell stencils on an AMR level. 
  @author Robert Marskar
*/

#ifndef CD_IrregStencil_H
#define CD_IrregStencil_H

// Chombo includes
#include <DisjointBoxLayout.H>
#include <EBISLayout.H>
#include <BaseIVFAB.H>
#include <EBCellFAB.H>
#include <RefCountedPtr.H>
#include <Stencils.H>
#include <LayoutData.H>
#include <VoFIterator.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for holding stencils on irregular cells over a single AMR level
  @details Users can use this class for holding VoFStencils on an AMR level, and to embed it in an AMR context
  through IrregAmrStencil. Applications functions do not use stencil aggregration, so do not use this class
  in performance-critical stencil code (like multigrid). To use this class, override the buildStencil function
  which will build the stencil in each cut-cell. 
  @note By default, a single stencil is allocated in the cut-cells. If you need more, override the define and
  apply functions. 
*/
class IrregStencil
{
public:
  /*!
    @brief Enum for identifying stencil -- only meant for enhancing code visibility. 
  */
  enum class StencilType
  {
    NoType,
    Linear,
    TaylorExtrapolation,
    LeastSquares,
    PiecewiseLinear
  };

  /*!
    @brief Weak constructor (does nothing)
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
  virtual ~IrregStencil();

  /*!
    @brief Get the stencils over a single box
  */
  const BaseIVFAB<VoFStencil>&
  operator[](const DataIndex& a_dit) const;

  /*!
    @brief Get the stencils over a single box
  */
  BaseIVFAB<VoFStencil>&
  operator[](const DataIndex& a_dit);

  /*!
    @brief Apply the stencil
  */
  virtual void
  apply(EBCellFAB& a_dst, const EBCellFAB& a_src, const DataIndex& a_dit) const;

  /*!
    @brief Apply the stencil
  */
  virtual void
  apply(BaseIVFAB<Real>& a_dst, const EBCellFAB& a_src, const DataIndex& a_dit) const;

protected:
  /*!
    @brief Default stencil component
  */
  static constexpr int m_defaultStenComp = 0;

  /*!
    @brief Default stencil component
  */
  static constexpr int m_defaultNumSten = 1;

  /*!
    @brief Stencil
  */
  LayoutData<RefCountedPtr<BaseIVFAB<VoFStencil>>> m_stencils;

  /*!
    @brief VoFIterators
  */
  mutable LayoutData<VoFIterator> m_vofIter;

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
  virtual void
  define(const DisjointBoxLayout&        a_dbl,
         const EBISLayout&               a_ebisl,
         const ProblemDomain&            a_domain,
         const Real&                     a_dx,
         const int                       a_order,
         const int                       a_radius,
         const IrregStencil::StencilType a_type);

  /*!
    @brief Build the desired stencil
  */
  virtual void
  buildStencil(VoFStencil&              a_sten,
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
