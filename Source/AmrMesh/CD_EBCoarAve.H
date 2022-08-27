/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
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
#include <EBLevelGrid.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Various averaging methods. 
  @details 'Arithmetic'/'Harmonic' mean what you think, while 'Conservative' also incorporate volume, face, 
  and EB face fractions such that the total quantity is conserved. 
*/
enum class Average
{
  Arithmetic,
  Harmonic,
  Conservative
};

/*!
  @brief Class which replaces data at coarse level of refinement with average
  at fine level of refinement. 
*/
class EBCoarAve
{
public:
  /*!
    @brief Default constructor. Must call define afterwards. 
  */
  EBCoarAve();

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
            const EBIndexSpace*      a_ebisPtr);

  /*!
    @brief Defining constructor. Calls the define function. 
    @param[in] a_eblgFine Fine grids
    @param[in] a_eblgCoar Coarse grids
    @param[in] a_eblgCoFi Coarsened fine grids
    @param[in] a_refRat   Refinement ratio between coarse and level
  */
  EBCoarAve(const EBLevelGrid& a_eblgFine,
            const EBLevelGrid& a_eblgCoar,
            const EBLevelGrid& a_eblgCoFi,
            const int&         a_refRat);

  /*!
    @brief Destructor (does nothing)
  */
  ~EBCoarAve();

  /*!
    @brief Define function -- puts operator in usable state.
    @param[in] a_eblgFine Fine grids
    @param[in] a_eblgCoar Coarse grids
    @param[in] a_eblgCoFi Coarsened fine grids
    @param[in] a_refRat   Refinement ratio between coarse and level
  */
  void
  define(const EBLevelGrid& a_eblgFine,
         const EBLevelGrid& a_eblgCoar,
         const EBLevelGrid& a_eblgCoFi,
         const int&         a_refRat);

  /*!
    @brief Assignement not allowed.
  */
  void
  operator=(const EBCoarAve& fabin) = delete;

  /*!
    @brief Do an average of cell data. 
    @param[inout] a_coarData  Coarse data
    @param[in]    a_fineData  Fine data
    @param[in]    a_variables Variables to coarsen
    @param[in]    a_average   Averaging method
    @note Input data holders should be defined over the input grids. 
  */
  void
  averageData(LevelData<EBCellFAB>&       a_coarData,
              const LevelData<EBCellFAB>& a_fineData,
              const Interval&             a_variables,
              const Average&              a_average);

  /*!
    @brief Do an average of face data. 
    @param[inout] a_coarData  Coarse data
    @param[in]    a_fineData  Fine data
    @param[in]    a_variables Variables to coarsen
    @param[in]    a_average   Averaging method
    @note Input data holders should be defined over the input grids. 
  */
  void
  averageData(LevelData<EBFluxFAB>&       a_coarData,
              const LevelData<EBFluxFAB>& a_fineData,
              const Interval&             a_variables,
              const Average&              a_average);

  /*!
    @brief Do an average of EB data. 
    @param[inout] a_coarData  Coarse data
    @param[in]    a_fineData  Fine data
    @param[in]    a_variables Variables to coarsen
    @param[in]    a_average   Averaging method
    @note Input data holders should be defined over the input grids. 
  */
  void
  averageData(LevelData<BaseIVFAB<Real>>&       a_coarData,
              const LevelData<BaseIVFAB<Real>>& a_fineData,
              const Interval&                   a_variables,
              const Average&                    a_average);

protected:
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
    @brief Irregular cells on the coarsened fine layout
  */
  LayoutData<IntVectSet> m_irregSetsCoFi;

  /*!
    @brief Refinement ratio
  */
  int m_refRat;

  /*!
    @brief Defined or not
  */
  bool m_isDefined;

  /*!
    @brief Do an arithmetic average of cell-centered data when coarsening. 
    @param[inout] a_coarData Coarse data
    @param[in]    a_fineData    Fine data
    @param[in]    a_datInd      Grid index
    @param[in]    a_fineInterv  Fine interval
    @param[in]    a_coarInterv  Coarse interval
  */
  void
  arithmeticAverage(EBCellFAB&       a_coarData,
                    const EBCellFAB& a_fineData,
                    const DataIndex& a_datInd,
                    const Interval&  a_coarInterv,
                    const Interval&  a_fineInterv);
  /*!
    @brief Do a harmonic average of cell-centered data when coarsening. 
    @param[inout] a_coarData Coarse data
    @param[in]    a_fineData    Fine data
    @param[in]    a_datInd      Grid index
    @param[in]    a_fineInterv  Fine interval
    @param[in]    a_coarInterv  Coarse interval
  */
  void
  harmonicAverage(EBCellFAB&       a_coarData,
                  const EBCellFAB& a_fineData,
                  const DataIndex& a_datInd,
                  const Interval&  a_coarInterv,
                  const Interval&  a_fineInterv);

  /*!
    @brief Do a conservative average of cell-centered data when coarsening. 
    @param[inout] a_coarData Coarse data
    @param[in]    a_fineData    Fine data
    @param[in]    a_datInd      Grid index
    @param[in]    a_fineInterv  Fine interval
    @param[in]    a_coarInterv  Coarse interval
  */
  void
  conservativeAverage(EBCellFAB&       a_coarData,
                      const EBCellFAB& a_fineData,
                      const DataIndex& a_datInd,
                      const Interval&  a_coarInterv,
                      const Interval&  a_fineInterv);

  /*!
    @brief Do an arithmetic average of face-centered data when coarsening. 
    @param[inout] a_coarData Coarse data
    @param[in]    a_fineData    Fine data
    @param[in]    a_datInd      Grid index
    @param[in]    a_fineInterv  Fine interval
    @param[in]    a_coarInterv  Coarse interval
    @param[in]    a_dir         Face direction
  */
  void
  arithmeticAverage(EBFaceFAB&       a_coarData,
                    const EBFaceFAB& a_fineData,
                    const DataIndex& a_datInd,
                    const Interval&  a_fineInterv,
                    const Interval&  a_coarInterv,
                    const int&       a_dir);

  /*!
    @brief Do a harmonic average of face-centered data when coarsening. 
    @param[inout] a_coarData Coarse data
    @param[in]    a_fineData    Fine data
    @param[in]    a_datInd      Grid index
    @param[in]    a_coarInterv  Coarse interval
    @param[in]    a_fineInterv  Fine interval
    @param[in]    a_dir         Face direction
  */
  void
  harmonicAverage(EBFaceFAB&       a_coarData,
                  const EBFaceFAB& a_fineData,
                  const DataIndex& a_datInd,
                  const Interval&  a_coarInterv,
                  const Interval&  a_fineInterv,
                  const int&       a_dir);

  /*!
    @brief Do a conservative average of face-centered data when coarsening. 
    @param[inout] a_coarData Coarse data
    @param[in]    a_fineData    Fine data
    @param[in]    a_datInd      Grid index
    @param[in]    a_coarInterv  Coarse interval
    @param[in]    a_fineInterv  Fine interval
    @param[in]    a_dir         Face direction
  */
  void
  conservativeAverage(EBFaceFAB&       a_coarData,
                      const EBFaceFAB& a_fineData,
                      const DataIndex& a_datInd,
                      const Interval&  a_coarInterv,
                      const Interval&  a_fineInterv,
                      const int&       a_dir);

  /*!
    @brief Arithmetic average data defined on irregular cells
    @param[inout] a_coarData   Coarse data
    @param[in]    a_fineData   Fine data
    @param[in]    a_datInd     Grid index
    @param[in]    a_coarInterv Coarse variables
    @param[in]    a_fineInterv Fine variables
  */
  void
  arithmeticAverage(BaseIVFAB<Real>&       a_coarData,
                    const BaseIVFAB<Real>& a_fineData,
                    const DataIndex&       a_datInd,
                    const Interval&        a_coarInterv,
                    const Interval&        a_fineINterv) const;

  /*!
    @brief Harmonic average data defined on irregular cells
    @param[inout] a_coarData   Coarse data
    @param[in]    a_fineData   Fine data
    @param[in]    a_datInd     Grid index
    @param[in]    a_coarInterv Coarse variables
    @param[in]    a_fineInterv Fine variables
  */
  void
  harmonicAverage(BaseIVFAB<Real>&       a_coarData,
                  const BaseIVFAB<Real>& a_fineData,
                  const DataIndex&       a_datInd,
                  const Interval&        a_coarInterv,
                  const Interval&        a_fineInterv) const;

  /*!
    @brief Conservatively average data defined on irregular cells
    @param[inout] a_coarData   Coarse data
    @param[in]    a_fineData   Fine data
    @param[in]    a_datInd     Grid index
    @param[in]    a_coarInterv Coarse variables
    @param[in]    a_fineInterv Fine variables
  */
  void
  conservativeAverage(BaseIVFAB<Real>&       a_coarData,
                      const BaseIVFAB<Real>& a_fineData,
                      const DataIndex&       a_datInd,
                      const Interval&        a_coarInterv,
                      const Interval&        a_fineInterv) const;
};

#include <CD_NamespaceFooter.H>

#endif
