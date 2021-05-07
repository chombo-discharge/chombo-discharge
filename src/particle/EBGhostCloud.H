/*!
  @file   EBGhostCloud
  @brief  A class for adding coarsened ghost cells to a coarse level
  @author Robert Marskar
  @date   April 2020
*/

#ifndef _EBGHOSTCLOUD_
#define _EBGHOSTCLOUD_

#include "DisjointBoxLayout.H"
#include "EBLevelGrid.H"
#include "EBCellFAB.H"
#include "ProblemDomain.H"

namespace ChomboDischarge {

  class EBGhostCloud {
  public:

    /*!
      @brief Default constructor. Needs to call define
    */
    EBGhostCloud();

    /*!
      @brief Full constructor. Calls define. 
    */
    EBGhostCloud(const DisjointBoxLayout& a_gridsCoar,
		 const DisjointBoxLayout& a_gridsFine,
		 const EBLevelGrid&       a_eblgCoar,
		 const EBLevelGrid&       a_eblgFine,
		 const ProblemDomain&     a_domainCoar,
		 const ProblemDomain&     a_domainFine,
		 const int                a_ref_rat,
		 const int                a_nComp,
		 const int                a_ghost);

    /*!
      @brief Dtor
    */
    ~EBGhostCloud();

    /*!
      @brief Define function
    */
    void define(const DisjointBoxLayout& a_gridsCoar,
		const DisjointBoxLayout& a_gridsFine,
		const EBLevelGrid&       a_eblgCoar,
		const EBLevelGrid&       a_eblgFine,
		const ProblemDomain&     a_domainCoar,
		const ProblemDomain&     a_domainFine,
		const int                a_ref_rat,
		const int                a_nComp,
		const int                a_ghost);

    /*!
      @brief Coarsen the fine-level ghost cells and add them to the coarse level
    */
    void addFineGhostsToCoarse(LevelData<EBCellFAB>& a_coarData, const LevelData<EBCellFAB>& a_fineData);

    /*!
      @brief Add the ghost cells on the refined coarse level to the fine level. 
    */
    void addFiCoDataToFine(LevelData<EBCellFAB>& a_fineData, const BoxLayoutData<FArrayBox>& a_fiCoData);

    /*!
      @brief Get the refined coarse buffer
    */
    BoxLayoutData<FArrayBox>& getFiCoBuffer();

    /*!
      @brief Get the refined coarse buffer
    */
    const BoxLayoutData<FArrayBox>& getFiCoBuffer() const;

    /*!
      @brief Get the refined version of the coarse eblevelgrid
    */
    const EBLevelGrid& getEblgFiCo() const;

  protected:

    /*!
      @brief Is defined or not
    */
    bool m_isDefined;

    /*!
      @brief Coarse grids
    */
    DisjointBoxLayout m_gridsCoar;

    /*!
      @brief Fine grids
    */
    DisjointBoxLayout m_gridsFine;

    /*!
      @brief Coarse EBLevelGrid
    */
    EBLevelGrid m_eblgCoar;

    /*!
      @brief Fine EBLevelGrid
    */
    EBLevelGrid m_eblgFine;

    /*!
      @brief Refined coarse EBLevelGrid
    */
    EBLevelGrid m_eblgFiCo;

    /*!
      @brief Coarse domain
    */
    ProblemDomain m_domainCoar;

    /*!
      @brief Fine domain
    */
    ProblemDomain m_domainFine;

    /*!
      @brief Coarsened grids, grown by the required amount of ghost cells. 
    */
    BoxLayout m_gridsCoFi;

    /*!
      @brief Refined coarse grids.
    */
    BoxLayout m_gridsFiCo;

    /*!
      @brief Data defined on the coarsened fine grid
    */
    BoxLayoutData<FArrayBox> m_dataCoFi;

    /*!
      @brief Data defined on the refined coarse grid
    */
    BoxLayoutData<FArrayBox> m_dataFiCo;

    /*!
      @brief Scratch storage for fine grid data
    */
    LevelData<EBCellFAB> m_scratchFine;

    /*!
      @brief Refinement ratio between the two levels
    */
    int m_refRat;

    /*!
      @brief Maximum number of components
    */
    int m_nComp;

    /*!
      @brief Required number of ghost cells
    */
    int m_ghost;

    /*!
      @brief Make the coarsened fine grids and data
    */
    void makeCoFiStuff();

    /*!
      @brief Make the refined coarse grids and data
    */
    void makeFiCoStuff();
  };
}
#endif
