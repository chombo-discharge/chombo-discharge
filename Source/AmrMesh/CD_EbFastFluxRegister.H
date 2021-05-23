/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EbFastFluxRegister.H
  @brief  Declaration of a new EBFluxRegister that uses some optimization
  @author Robert Marskar
*/

#ifndef CD_EbFastFluxRegister_H
#define CD_EbFastFluxRegister_H

// Chombo includes
#include <EBFluxRegister.H>

// Our includes
#include <CD_NamespaceHeader.H>
  
/*!
  @brief Version of EBFluxRegister with faster define function. 
*/
class EbFastFluxRegister : public EBFluxRegister { 
public:

  
  /*!
    @brief Default constructor. Leaves object in undefined state.
  */
  EbFastFluxRegister();
  
  /*!
    @brief Full constructor.  Calls the define function which creates
    a levels worth of flux registers.
    @param[in] a_eblgFine    The fine level EBLevelGrid.
    @param[in] a_eblgCoar    The coarse level EBLevelGrid.
    @param[in] a_nRef        Refinement ratio between levels.
    @param[in] a_nVar        Number of variables. 
    @param[in] a_forceNoEBCF Forces EBCF functionality to be turned off if true.
  */
  EbFastFluxRegister(const EBLevelGrid& a_eblgFine,
		     const EBLevelGrid& a_eblgCoar,
		     const int&         a_nRef,
		     const int&         a_nVar,
		     const bool         a_forceNoEBCF = false);

  /*!
    @brief Full constructor.  Calls the define function which creates
    a levels worth of flux registers.
    @param[in] a_dblFine     Fine grids
    @param[in] a_dblCoar     Coarse grids
    @param[in] a_ebislFine   EBISLayout on the fine level
    @param[in] a_ebislCoar   EBISLayout on the coarse level
    @param[in] a_domainCoar  ProblemDomain on the coarse level.
    @param[in] a_nRef        Refinement ratio between levels.
    @param[in] a_nVar        Number of variables. 
    @param[in] a_ebisPtr     EBIndexSpace 
    @param[in] a_forceNoEBCF Forces EBCF functionality to be turned off if true.
  */
  EbFastFluxRegister(const DisjointBoxLayout& a_dblFine,
		     const DisjointBoxLayout& a_dblCoar,
		     const EBISLayout&        a_ebislFine,
		     const EBISLayout&        a_ebislCoar,
		     const Box&               a_domainCoar,
		     const int&               a_nRef,
		     const int&               a_nVar,
		     const EBIndexSpace*      ebisPtr,
		     const bool               a_forceNoEBCF = false);

  /*!
    @brief Dtor
  */
  virtual ~EbFastFluxRegister();

  /*!
    @brief Full define function.  Creates a levels worth of flux registers.
    @param[in] a_dblFine     Fine grids
    @param[in] a_dblCoar     Coarse grids
    @param[in] a_ebislFine   EBISLayout on the fine level
    @param[in] a_ebislCoar   EBISLayout on the coarse level
    @param[in] a_domainCoar  ProblemDomain on the coarse level.
    @param[in] a_nRef        Refinement ratio between levels.
    @param[in] a_nVar        Number of variables. 
    @param[in] a_ebisPtr     EBIndexSpace 
    @param[in] a_forceNoEBCF Forces EBCF functionality to be turned off if true.
  */
  void define(const DisjointBoxLayout& a_dblFine,
	      const DisjointBoxLayout& a_dblCoar,
	      const EBISLayout&        a_ebislFine,
	      const EBISLayout&        a_ebislCoar,
	      const ProblemDomain&     a_domainCoar,
	      const int&               a_nRef,
	      const int&               a_nVar,
	      const EBIndexSpace*      ebisPtr,
	      const bool               a_forceNoEBCF = false);

  /*!
    @brief Define function.
    @param[in] a_eblgFine    The fine level EBLevelGrid.
    @param[in] a_eblgCoar    The coarse level EBLevelGrid.
    @param[in] a_nRef        Refinement ratio between levels.
    @param[in] a_nVar        Number of variables. 
    @param[in] a_forceNoEBCF Forces EBCF functionality to be turned off if true.
  */
  void define(const EBLevelGrid&       a_eblgFine,
	      const EBLevelGrid&       a_eblgCoar,
	      const int&               a_refRat,
	      const int&               a_nVar,
	      const bool               a_forceNoEBCF = false);

protected:
  
  /*!
    @brief Mucho faster version than the one in EBFastFR
  */
  void fastDefineSetsAndIterators();

  /*!
    @brief Define masks
  */
  void defineMasks();

  /*!
    @brief This is a mask for describing which coar cells boarder the CFIVS
  */
  LayoutData<IntVectSet> m_coarCFIVS[2*SpaceDim];

private:

  /*!
    @brief Assignment operator - disallowed for all the usual reasons
  */
  void operator=(const EbFastFluxRegister& out) = delete;

  /*!
    @brief Copy constructor is not allowed.
  */
  EbFastFluxRegister(const EbFastFluxRegister& out) = delete;
};
  
#include <CD_NamespaceFooter.H>

#endif

