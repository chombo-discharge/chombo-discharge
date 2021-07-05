/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBMultigridProlong.H
  @brief  Declaration of a class for doing multigrid restriction
  @author Robert Marskar
*/

#ifndef CD_EBMultigridProlong_H
#define CD_EBMultigridProlong_H

// Chombo includes
#include <EBLevelGrid.H>

// Our includes
#include <CD_NamespaceHeader.H>

class EBMultigridProlong {
public:

  /*!
    @brief Disallowed default constructor
  */
  EBMultigridProlong();

  /*!
    @brief Full constructor
    @param[in] a_eblgFine       Fine levelgrid
    @param[in] a_eblgCoar       Coarse grid
    @param[in] a_refRat         Refinement factor
    @param[in] a_volumeWeighted Use volume-weighted restriction or not
  */
  EBMultigridProlong(const EBLevelGrid& a_eblgFine,
		     const EBLevelGrid& a_eblgCoar,
		     const int&         a_refRat,
		     const bool&        a_volumeWeighted);

  /*!
    @brief Disallowed copy constructor
  */
  EBMultigridProlong(const EBMultigridProlong& a_other) = delete;

  /*!
    @brief Destructor.
  */  
  ~EBMultigridProlong();

  /*!
    @brief Disallowed assignment constructor
  */
  EBMultigridProlong& operator=(const EBMultigridProlong& a_other) = delete;

  /*!
    @brief Define function
    @param[in] a_eblgFine       Fine grids
    @param[in] a_eblgCoar       Coarse grids
    @param[in] a_refRat         Refinement factor
    @param[in] a_volumeWeighted Use volume-weighted prolongation or not
  */
  void define(const EBLevelGrid& a_eblgFine,
	      const EBLevelGrid& a_eblgCoar,
	      const int&         a_refRat,
	      const bool&        a_volumeWeighted);

  /*!
    @brief Prolong onto coarse
    @param[inout] a_fineData  Fine data to be interpolated.
    @param[in]    a_coarData  Coarse data 
    @param[in]    a_variables Variables to restrict. 
  */
  void prolong(LevelData<EBCellFAB>& a_fineData, const LevelData<EBCellFAB>& a_coarData, const Interval& a_variables) const;
  
protected:

  /*!
    @brief Internal component number (for storage). 
    @details Does not matter outside of class.
  */
  static constexpr int m_comp = 0;

  /*!
    @brief Internal number of components (for storage). 
    @details Does not matter outside of class.
  */
  static constexpr int m_nComp = 1;

  /*!
    @brief Fine grids
  */
  EBLevelGrid m_eblgFine;

  /*!
    @brief Coarse grids
  */
  EBLevelGrid m_eblgCoar;

  /*!
    @brief Refine coarse grids
  */
  EBLevelGrid m_eblgFiCo;

  /*!
    @brief Refinement ratio between fine/coarse
  */
  int m_refRat;

  /*!
    @brief Use volume weighted restriction or not
  */
  bool m_volumeWeighted;

  /*!
    @brief Is defined or not
  */
  bool m_isDefined;

  /*!
    @brief Refined coarse data. 
  */
  mutable LevelData<EBCellFAB> m_fineData;
  
  /*!
    @brief VoF iterator for irregular data
  */
  mutable LayoutData<VoFIterator> m_vofitIrregFine;

  /*!
    @brief Averaging stencils
  */
  LayoutData<BaseIVFAB<VoFStencil> > m_prolongStencils;

  /*!
    @brief Define stencils
  */
  void defineStencils();
};

#include <CD_NamespaceFooter.H>

#endif
