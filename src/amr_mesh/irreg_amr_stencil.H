/*!
  @file irreg_amr_stencil.H
  @brief Class for holding stencils on irregular over an entire AMR hierarchy
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _IRREG_AMR_STENCIL_
#define _IRREG_AMR_STENCIL_

#include "EBAMRData.H"
#include "irreg_stencil.H"
#include "type_definitions.H"

/*!
  @brief Class for holding stencils on irregular cells over an entire AMR hierarchy
  @note This class holds a single stencil on each irregular cells. When you apply the stencil to a data holder, the same stencil will be applied to all the components of that holder. 
*/
template <class T> class irreg_amr_stencil {
public:

  /*!
    @brief Empty constructor. You must subsequently call define. 
  */
  irreg_amr_stencil();

  /*!
    @brief Full constructor. Subsequently calls define
  */
  irreg_amr_stencil(const Vector<DisjointBoxLayout>& a_grids,
		    const Vector<EBISLayout>&        a_ebisl,
		    const Vector<ProblemDomain>&     a_domains,
		    const Vector<Real>&              a_dx,
		    const int                        a_finest_level,
		    const int                        a_order,
		    const int                        a_radius,
		    const stencil_type               a_type);

  /*!
    @brief Destructor
  */
  virtual ~irreg_amr_stencil();
  

  /*!
    @brief Define function
  */
  virtual void define(const Vector<DisjointBoxLayout>& a_grids,
		      const Vector<EBISLayout>&        a_ebisl,
		      const Vector<ProblemDomain>&     a_domains,
		      const Vector<Real>&              a_dx,
		      const int                        a_finest_level,
		      const int                        a_order,
		      const int                        a_radius,
		      const stencil_type               a_type);

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

  /*!
    @brief Get the stencils over an AMR level
    @param[in] a_lvl AMR level
  */
  const irreg_stencil& operator[](const int& a_lvl) const {
    return *m_stencils[a_lvl];
  }

  /*!
    @brief Get the stencils over an AMR level
    @param[in] a_lvl AMR level
  */
  irreg_stencil& operator[](const int& a_lvl) {
    return *m_stencils[a_lvl];
  }

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
    @brief Stencils
  */
  Vector<RefCountedPtr<irreg_stencil> > m_stencils;

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
  int m_finest_level;

  /*!
    @brief Stencil type
  */
  stencil_type m_stencil_type;
};

#include "irreg_amr_stencilI.H"

#endif
