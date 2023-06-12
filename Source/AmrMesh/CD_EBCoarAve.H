/* chombo-discharge
 * Copyright © 2023 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBCoarAve.H
  @brief  Declaration of conservative coarsening utility 
  @author Robert Marskar
*/

#ifndef CD_EBCoarAve_H
#define CD_EBCoarAve_H

// Chombo includes
#include <AggStencil.H>
#include <EBLevelGrid.H>
#include <VoFIterator.H>
#include <FaceIterator.H>

// Our includes
#include <CD_Average.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class which replaces data at coarse level of refinement with average
  at fine level of refinement. 
  @details When defining this class, use the EBLevelGrid-based constructors if the coarsened grids are available. 
*/
class EBCoarAve
{
public:
  /*!
    @brief Type of coarsening supported
  */
  enum class Type
  {
    Arithmetic,
    Conservative,
    Harmonic
  };

  /*!
    @brief Default constructor. Must call define afterwards. 
  */
  EBCoarAve() noexcept;

  /*!
    @brief Copy constructor not allowed
  */
  EBCoarAve(const EBCoarAve& a_other) = delete;

  /*!
    @brief Defining constructor.
    @details This version creates it's own coarsened grids
    @param[in] a_dblFine    Fine grids
    @param[in] a_dblCoar    Coarse grids
    @param[in] a_ebislFine  Fine EBISL
    @param[in] a_ebislCoar  Coarse EBISL
    @param[in] a_domainCoar Coarse domain
    @param[in] a_refRat     Refinement ratio
    @param[in] a_ebisPtr    EBIS pointer.
  */
  EBCoarAve(const DisjointBoxLayout& a_dblFine,
            const DisjointBoxLayout& a_dblCoar,
            const EBISLayout&        a_ebislFine,
            const EBISLayout&        a_ebislCoar,
            const ProblemDomain&     a_domainCoar,
            const int&               a_refRat,
            const EBIndexSpace*      a_ebisPtr) noexcept;

  /*!
    @brief Defining constructor. Calls the define function. 
    @param[in] a_eblgFine  Fine grids
    @param[in] a_eblgCoar  Coarse grids
    @param[in] a_eblgCoFi  Coarsened fine grids
    @param[in] a_refRat    Refinement ratio between coarse and level
  */
  EBCoarAve(const EBLevelGrid& a_eblgFine,
            const EBLevelGrid& a_eblgCoar,
            const EBLevelGrid& a_eblgCoFi,
            const int&         a_refRat) noexcept;

  /*!
    @brief Destructor (does nothing)
  */
  virtual ~EBCoarAve() noexcept;

  /*!
    @brief Define function -- puts operator in usable state.
    @param[in] a_eblgFine Fine grids
    @param[in] a_eblgCoar Coarse grids
    @param[in] a_eblgCoFi Coarsened fine grids
    @param[in] a_refRat   Refinement ratio between coarse and level
  */
  virtual void
  define(const EBLevelGrid& a_eblgFine,
         const EBLevelGrid& a_eblgCoar,
         const EBLevelGrid& a_eblgCoFi,
         const int&         a_refRat) noexcept;

  /*!
    @brief Assignement not allowed.
  */
  virtual void
  operator=(const EBCoarAve& fabin) = delete;

  /*!
    @brief Do an average of cell data. 
    @param[inout] a_coarData  Coarse data
    @param[in]    a_fineData  Fine data
    @param[in]    a_variables Variables to coarsen
    @param[in]    a_average   Averaging method
  */
  virtual void
  averageData(LevelData<EBCellFAB>&       a_coarData,
              const LevelData<EBCellFAB>& a_fineData,
              const Interval&             a_variables,
              const Average&              a_average) const noexcept;

  /*!
    @brief Do an average of face data. 
    @param[inout] a_coarData  Coarse data
    @param[in]    a_fineData  Fine data
    @param[in]    a_variables Variables to coarsen
    @param[in]    a_average   Averaging method
  */
  virtual void
  averageData(LevelData<EBFluxFAB>&       a_coarData,
              const LevelData<EBFluxFAB>& a_fineData,
              const Interval&             a_variables,
              const Average&              a_average) const noexcept;

  /*!
    @brief Do an average of EB data. 
    @param[inout] a_coarData  Coarse data
    @param[in]    a_fineData  Fine data
    @param[in]    a_variables Variables to coarsen
    @param[in]    a_average   Averaging method
  */
  virtual void
  averageData(LevelData<BaseIVFAB<Real>>&       a_coarData,
              const LevelData<BaseIVFAB<Real>>& a_fineData,
              const Interval&                   a_variables,
              const Average&                    a_average) const noexcept;

protected:
  /*!
    @brief Defined or not
  */
  bool m_isDefined;

  /*!
    @brief Refinement ratio
  */
  int m_refRat;

  /*!
    @brief Fine grids
  */
  EBLevelGrid m_eblgFine;

  /*!
    @brief Coarse grids
  */
  EBLevelGrid m_eblgCoar;

  /*!
    @brief Coarsened fine grids
  */
  EBLevelGrid m_eblgCoFi;

  /*!
    @brief Pre-built copier for copying from m_eblgCoFi to m_eblgCoar
  */
  Copier m_cellCopier;

  /*!
    @brief Pre-built copier for copying from m_eblgCoFi to m_eblgCoar
  */
  Copier m_faceCopier;

  /*!
    @brief Pre-built copier for copying from m_eblgCoFi to m_eblgCoar
  */
  Copier m_ebCopier;

  /*!
    @brief Irregular cells on the coarsened fine layout
  */
  LayoutData<IntVectSet> m_irregSetsCoFi;

  /*!
    @brief Irregular cells on the coarsened fine layout
  */
  mutable LayoutData<VoFIterator> m_irregCellsCoFi;

  /*!
    @brief Irregular faces on the coarsened fine layout
  */
  mutable LayoutData<std::map<int, FaceIterator>> m_irregFacesCoFi;

  /*!
    @brief Stencils for arithmetic coarsening of volume data
  */
  LayoutData<BaseIVFAB<VoFStencil>> m_cellArithmeticStencils;

  /*!
    @brief Stencils for harmonic coarsening of volume data
    @note This is different than the cellConservative and cellArithmetic stencils because its weights are all one and
    the harmonic average is applied differently. 
  */
  LayoutData<BaseIVFAB<VoFStencil>> m_cellHarmonicStencils;

  /*!
    @brief Stencils for conservative coarsening of volume data
  */
  LayoutData<BaseIVFAB<VoFStencil>> m_cellConservativeStencils;

  /*!
    @brief Stencils for arithmetic coarsening of face data
  */
  LayoutData<std::map<int, BaseIFFAB<FaceStencil>>> m_faceArithmeticStencils;

  /*!
    @brief Stencils for harmonic coarsening of face data
  */
  LayoutData<std::map<int, BaseIFFAB<FaceStencil>>> m_faceHarmonicStencils;

  /*!
    @brief Stencils for conservative coarsening of face data
  */
  LayoutData<std::map<int, BaseIFFAB<FaceStencil>>> m_faceConservativeStencils;

  /*!
    @brief Stencils for arithmetic coarsening of EB data
  */
  LayoutData<BaseIVFAB<VoFStencil>> m_ebArithmeticStencils;

  /*!
    @brief Stencils for harmonic coarsening of EB data
  */
  LayoutData<BaseIVFAB<VoFStencil>> m_ebHarmonicStencils;

  /*!
    @brief Stencils for conservative coarsening of EB data
  */
  LayoutData<BaseIVFAB<VoFStencil>> m_ebConservativeStencils;

  /*!
    @brief Define stencils for the various coarsening functions of cell data
  */
  virtual void
  defineCellStencils() noexcept;

  /*!
    @brief Define stencils for the various coarsening functions of face data
  */
  virtual void
  defineFaceStencils() noexcept;

  /*!
    @brief Define stencils for the various coarsening functions of EB data
  */
  virtual void
  defineEBStencils() noexcept;

  /*!
    @brief Define buffer storage
  */
  virtual void
  defineBuffers() noexcept;

  /*!
    @brief Do an arithmetic average of cell-centered data when coarsening. 
    @param[inout] a_coarData Coarse data
    @param[in]    a_fineData Fine data
    @param[in]    a_datInd Grid index
    @param[in]    a_coarVar Coarse variables 
    @param[in]    a_fineVar Fine variable
  */
  virtual void
  arithmeticAverage(EBCellFAB&       a_coarData,
                    const EBCellFAB& a_fineData,
                    const DataIndex& a_datInd,
                    const int&       a_coarVar,
                    const int&       a_fineVar) const noexcept;
  /*!
    @brief Do a harmonic average of cell-centered data when coarsening. 
    @param[inout] a_coarData Coarse data
    @param[in]    a_fineData Fine data
    @param[in]    a_datInd Grid index
    @param[in]    a_coarVar Coarse variables 
    @param[in]    a_fineVar Fine variable
  */
  virtual void
  harmonicAverage(EBCellFAB&       a_coarData,
                  const EBCellFAB& a_fineData,
                  const DataIndex& a_datInd,
                  const int&       a_coarVar,
                  const int&       a_fineVar) const noexcept;

  /*!
    @brief Do a conservative average of cell-centered data when coarsening. 
    @param[inout] a_coarData Coarse data
    @param[in]    a_fineData Fine data
    @param[in]    a_datInd Grid index
    @param[in]    a_coarVar Coarse variables 
    @param[in]    a_fineVar Fine variable
  */
  virtual void
  conservativeAverage(EBCellFAB&       a_coarData,
                      const EBCellFAB& a_fineData,
                      const DataIndex& a_datInd,
                      const int&       a_coarVar,
                      const int&       a_fineVar) const noexcept;

  /*!
    @brief Do an arithmetic average of face-centered data when coarsening. 
    @param[inout] a_coarData Coarse data
    @param[in]    a_fineData Fine data
    @param[in]    a_datInd Grid index
    @param[in]    a_coarVar Coarse variable
    @param[in]    a_fineVar Fine variable
    @param[in]    a_dir Face direction
  */
  virtual void
  arithmeticAverage(EBFaceFAB&       a_coarData,
                    const EBFaceFAB& a_fineData,
                    const DataIndex& a_datInd,
                    const int&       a_fineVar,
                    const int&       a_coarVar,
                    const int&       a_dir) const noexcept;

  /*!
    @brief Do a harmonic average of face-centered data when coarsening. 
    @param[inout] a_coarData Coarse data
    @param[in]    a_fineData Fine data
    @param[in]    a_datInd Grid index
    @param[in]    a_coarVar Coarse variable
    @param[in]    a_fineVar Fine variable
    @param[in]    a_dir Face direction
  */
  virtual void
  harmonicAverage(EBFaceFAB&       a_coarData,
                  const EBFaceFAB& a_fineData,
                  const DataIndex& a_datInd,
                  const int&       a_coarVar,
                  const int&       a_fineVar,
                  const int&       a_dir) const noexcept;

  /*!
    @brief Do a conservative average of face-centered data when coarsening. 
    @param[inout] a_coarData Coarse data
    @param[in]    a_fineData Fine data
    @param[in]    a_datInd Grid index
    @param[in]    a_coarVar Coarse variable
    @param[in]    a_fineVar Fine variable
    @param[in]    a_dir Face direction
  */
  virtual void
  conservativeAverage(EBFaceFAB&       a_coarData,
                      const EBFaceFAB& a_fineData,
                      const DataIndex& a_datInd,
                      const int&       a_coarVar,
                      const int&       a_fineVar,
                      const int&       a_dir) const noexcept;

  /*!
    @brief Arithmetic average data defined on irregular cells
    @param[inout] a_coarData Coarse data
    @param[in]    a_fineData Fine data
    @param[in]    a_datInd Grid index
    @param[in]    a_coarVar Coarse variable
    @param[in]    a_finevar Fine Variable
  */
  virtual void
  arithmeticAverage(BaseIVFAB<Real>&       a_coarData,
                    const BaseIVFAB<Real>& a_fineData,
                    const DataIndex&       a_datInd,
                    const int&             a_coarVar,
                    const int&             a_fineVar) const noexcept;

  /*!
    @brief Harmonic average data defined on irregular cells
    @param[inout] a_coarData Coarse data
    @param[in]    a_fineData Fine data
    @param[in]    a_datInd Grid index
    @param[in]    a_coarVar Coarse variable
    @param[in]    a_finevar Fine Variable
  */
  virtual void
  harmonicAverage(BaseIVFAB<Real>&       a_coarData,
                  const BaseIVFAB<Real>& a_fineData,
                  const DataIndex&       a_datInd,
                  const int&             a_coarVar,
                  const int&             a_fineVar) const noexcept;

  /*!
    @brief Conservatively average data defined on irregular cells
    @param[inout] a_coarData Coarse data
    @param[in]    a_fineData Fine data
    @param[in]    a_datInd Grid index
    @param[in]    a_coarVar Coarse variable
    @param[in]    a_finevar Fine Variable
  */
  virtual void
  conservativeAverage(BaseIVFAB<Real>&       a_coarData,
                      const BaseIVFAB<Real>& a_fineData,
                      const DataIndex&       a_datInd,
                      const int&             a_coarVar,
                      const int&             a_fineVar) const noexcept;
};

#include <CD_NamespaceFooter.H>

#endif
