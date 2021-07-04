/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBMultigridRestrict.H
  @brief  Declaration of a class for doing multigrid restriction
  @author Robert Marskar
*/

#ifndef CD_EBMultigridRestrict_H
#define CD_EBMultigridRestrict_H

// Chombo includes
#include <EBLevelGrid.H>

// Our includes
#include <CD_NamespaceHeader.H>

class EBMultigridRestrict {
public:

  /*!
    @brief Disallowed default constructor
  */
  EBMultigridRestrict();

  /*!
    @brief Full constructor
    @param[in] a_eblgFine       Fine levelgrid
    @param[in] a_eblgCoar       Coarse grid
    @param[in] a_refRat         Refinement factor
    @param[in] a_volumeWeighted Use volume-weighted restriction or not
  */
  EBMultigridRestrict(const EBLevelGrid& a_eblgFine,
		      const EBLevelGrid& a_eblgCoar,
		      const int&         a_refRat,
		      const bool&        a_volumeWeighted);

  /*!
    @brief Disallowed copy constructor
  */
  EBMultigridRestrict(const EBMultigridRestrict& a_other) = delete;

  /*!
    @brief Destructor.
  */  
  ~EBMultigridRestrict();

  /*!
    @brief Disallowed assignment constructor
  */
  EBMultigridRestrict& operator=(const EBMultigridRestrict& a_other) = delete;

  /*!
    @brief Define function
  */
  void define(const EBLevelGrid& a_eblgFine,
	      const EBLevelGrid& a_eblgCoar,
	      const int&         a_refRat,
	      const bool&        a_volumeWeighted);

  /*!
    @brief Restrict onto coarse
    @param[inout] a_coarData  Coarse data to be restricted. 
    @param[in]    a_fineData  Fine data
    @param[in]    a_variables Variables to restrict. 
  */
  void restrict(LevelData<EBCellFAB>& a_coarData, const LevelData<EBCellFAB>& a_fineData, const Interval& a_variables) const;
  
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
    @brief Coarsened fine grids
  */
  EBLevelGrid m_eblgCoFi;

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
    @brief Buffer for coarsened fine data. 
  */
  mutable LevelData<EBCellFAB> m_coFiData;
  
  /*!
    @brief VoF iterator for irregular data
  */
  mutable LayoutData<VoFIterator> m_vofitIrreg;

  /*!
    @brief Averaging stencils
  */
  LayoutData<BaseIVFAB<VoFStencil> > m_avgStencils;

  /*!
    @brief Define stencils
  */
  void defineStencils();
};

#include <CD_NamespaceFooter.H>

#endif
