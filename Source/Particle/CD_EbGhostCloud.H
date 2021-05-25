/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EbGhostCloud.H
  @brief  A class for adding ghost cells from the fine level to the coarse level
  @author Robert Marskar
*/

#ifndef CD_EbGhostCloud_H
#define CD_EbGhostCloud_H

// Chombo includes
#include <DisjointBoxLayout.H>
#include <EBLevelGrid.H>
#include <EBCellFAB.H>
#include <ProblemDomain.H>

// Our includes
#include <CD_NamespaceHeader.H>

class EbGhostCloud {
public:

  /*!
    @brief Default constructor. Needs to call define
  */
  EbGhostCloud();

  /*!
    @brief Full constructor. Calls define. 
  */
  EbGhostCloud(const DisjointBoxLayout& a_gridsCoar,
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
  ~EbGhostCloud();

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

#include <CD_NamespaceFooter.H>

#endif
