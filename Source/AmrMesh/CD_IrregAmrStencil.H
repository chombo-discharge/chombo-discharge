/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_IrregAmrStencil.H
  @brief  Class for holding stencils on irregular over an entire AMR hierarchy
  @author Robert Marskar
*/

#ifndef CD_IrregAmrStencil_H
#define CD_IrregAmrStencil_H

// Our includes
#include <CD_EBAMRData.H>
#include <CD_IrregStencil.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for holding VoFStencils on irregular cells over an entire AMR hierarchy. The template
  parameter must inherit from IrregStencil. The stencil can only be applied to EBCellFABs. 
  @details This class holds a single stencil on each irregular cells. When you apply 
  the stencil to a data holder, is will be applie to all the components of that data
  holder. 
  @note This does not use aggregate stencil functionality. Do not use this class inside multigrid!
*/
template <class IrregSten>
class IrregAmrStencil {
public:

  /*!
    @brief Empty constructor. You must subsequently call define. 
  */
  IrregAmrStencil();

  /*!
    @brief Full constructor. Subsequently calls define
    @param[in] a_grids       AMR grids
    @param[in] a_ebisl       EBIS layouts on each level
    @param[in] a_domains     Domains on each level
    @param[in] a_dx          Resolutions on each level
    @param[in] a_finestLevel Finest AMR level
    @param[in] a_order       Stencil order
    @param[in] a_radius      Stencil radius
    @param[in] a_type        Stencil type
  */
  IrregAmrStencil(const Vector<DisjointBoxLayout>& a_grids,
		  const Vector<EBISLayout>&        a_ebisl,
		  const Vector<ProblemDomain>&     a_domains,
		  const Vector<Real>&              a_dx,
		  const int                        a_finestLevel,
		  const int                        a_order,
		  const int                        a_radius,
		  const IrregStencil::StencilType  a_type);

  /*!
    @brief Destructor
  */
  virtual ~IrregAmrStencil();

  /*!
    @brief Get the stencils over an AMR level
    @param[in] a_lvl AMR level
  */
  const IrregStencil& operator[](const int& a_lvl) const;

  /*!
    @brief Get the stencils over an AMR level
    @param[in] a_lvl AMR level
  */
  IrregStencil& operator[](const int& a_lvl);

  /*!
    @brief Define function
    @param[in] a_grids       AMR grids
    @param[in] a_ebisl       EBIS layouts on each level
    @param[in] a_domains     Domains on each level
    @param[in] a_dx          Resolutions on each level
    @param[in] a_finestLevel Finest AMR level
    @param[in] a_order       Stencil order
    @param[in] a_radius      Stencil radius
    @param[in] a_type        Stencil type
  */
  virtual void define(const Vector<DisjointBoxLayout>& a_grids,
		      const Vector<EBISLayout>&        a_ebisl,
		      const Vector<ProblemDomain>&     a_domains,
		      const Vector<Real>&              a_dx,
		      const int                        a_finestLevel,
		      const int                        a_order,
		      const int                        a_radius,
		      const IrregStencil::StencilType  a_type);

  /*!
    @brief Apply the stencils to an existing data holder. 
    @details This routine will perform a copy of all regular cells, and apply stencils to irregular cells only. 
    @param[out] a_dst Destination data. 
    @param[in]  a_src Source data
    @param[in]  a_lvl Level
    @param[in]  a_limit Limiter. This disallows sign changes in the interpolation (floored to zero)
  */
  virtual void apply(LevelData<EBCellFAB>&       a_dst,
		     const LevelData<EBCellFAB>& a_src,
		     const int                   a_lvl,
		     const bool                  a_limit = false) const;

  /*!
    @brief Apply stencils to data holder. This replaces original data
    @param[inout] a_data Transformed data
    @param[in]    a_lvl Level
    @param[in]    a_limit Limiter. This disallows sign changes in the interpolation (floored to zero)
  */
  virtual void apply(LevelData<EBCellFAB>& a_data,
		     const int             a_lvl,
		     const bool            a_limit = false) const;

  /*!
    @brief Apply the stencils to an existing data holder. 
    @details This routine will perform a copy of all regular cells, and apply stencils to irregular cells only. 
    @param[out] a_dst Destination data. 
    @param[in]  a_src Source data
    @param[in]  a_lvl Level
    @param[in]  a_limit Limiter. This disallows sign changes in the interpolation (floored to zero)
  */
  virtual void apply(LevelData<BaseIVFAB<Real> >& a_dst,
		     const LevelData<EBCellFAB>&  a_src,
		     const int                    a_lvl,
		     const bool                   a_limit = false) const;
  
  /*!
    @brief Apply the stencils to an existing data holder. 
    @details This routine will perform a copy of all regular cells, and apply stencils to irregular cells only. 
    @param[out] a_dst Destination data. 
    @param[in]  a_src Source data
    @param[in]  a_limit Limiter. This disallows sign changes in the interpolation (floored to zero)
  */
  virtual void apply(EBAMRCellData&       a_dst,
		     const EBAMRCellData& a_src,
		     const bool         a_limit = false) const;

  /*!
    @brief Apply the stencils to an existing data holder. This replaces original data.
    @param[inout] a_dst Data
    @param[in]  a_limit Limiter. This disallows sign changes in the interpolation (floored to zero)
  */
  virtual void apply(EBAMRCellData& a_dst,
		     const bool   a_limit = false) const;

  /*!
    @brief Apply the stencils to an existing data holder
    @param[out] a_dst Destination data. 
    @param[in]  a_src Source data
    @param[in]  a_limit Limiter. This disallows sign changes in the interpolation (floored to zero)
  */
  virtual void apply(EBAMRIVData&         a_dst,
		     const EBAMRCellData& a_src,
		     const bool         a_limit = false) const;

  /*!
    @brief Apply stencils to an existing data holders. This replaces original data
    @param[inout] a_data
    @param[in]  a_limit Limiter. This disallows sign changes in the interpolation (floored to zero)
  */
  virtual void apply(Vector<EBAMRCellData*>& a_data,
		     const bool            a_limit = false) const;

  /*!
    @brief Apply stencils to existing data holders
    @param[out] a_dst Destination data holders. 
    @param[in]  a_src Source data holders
    @param[in]  a_limit Limiter. This disallows sign changes in the interpolation (floored to zero)
  */
  virtual void apply(Vector<EBAMRCellData*>&       a_dst,
		     const Vector<EBAMRCellData*>& a_src,
		     const bool                  a_limit = false) const;

  /*!
    @brief Apply stencils to existing data holders
    @param[out] a_dst Destination data holders. 
    @param[in]  a_src Source data holders
    @param[in]  a_limit Limiter. This disallows sign changes in the interpolation (floored to zero)
  */
  virtual void apply(Vector<EBAMRIVData*>&         a_dst,
		     const Vector<EBAMRCellData*>& a_src,
		     const bool                  a_limit = false) const;



protected:

  /*!
    @brief AMR grids
  */
  Vector<DisjointBoxLayout> m_grids;

  /*!
    @brief EBIS layouts
  */
  Vector<EBISLayout> m_ebisl;

  /*!
    @brief Problem domains
  */
  Vector<ProblemDomain> m_domains;

  /*!
    @brief Resolutions
  */
  Vector<Real> m_dx;
  
  /*!
    @brief Stencils. 
  */
  Vector<RefCountedPtr<IrregStencil> > m_stencils;

  /*!
    @brief Stencils defined or not
  */
  bool m_defined;

  /*!
    @brief Order
  */
  int m_order;

  /*!
    @brief Radius
  */
  int m_radius;
  
  /*!
    @brief Finest level
  */
  int m_finestLevel;

  /*!
    @brief Stencil type
  */
  IrregStencil::StencilType m_stencilType;
};

#include <CD_NamespaceFooter.H>

#include <CD_IrregAmrStencilImplem.H>

#endif
