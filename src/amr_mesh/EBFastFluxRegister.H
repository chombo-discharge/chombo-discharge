/*!
  @file   EBFastFluxRegister.H
  @brief  Declaration of a new EBFluxRegister that uses some optimization
  @author Robert Marskar
  @date   Aug. 2020
*/

#ifndef _EBFastFluxRegister_H_
#define _EBFastFluxRegister_H_

#include "EBFastFluxRegister.H"
#include "EBFluxRegister.H"

class EBCoarToFineRedist;
class EBCoarToCoarRedist;

/*!
  @brief Class that does refluxing and interaction with redistribution classes. 
*/
class EBFastFluxRegister : public EBFluxRegister { // Forgive me Father, for I have sinned. 
public:

  
  /*!
    @brief Default constructor.  Leaves object undefined.
  */
  EBFastFluxRegister();
  
  /*!
    @brief Full constructor.  Calls the define function which creates
    a levels worth of flux registers.
  */
  EBFastFluxRegister(const EBLevelGrid& a_eblgFine,
		     const EBLevelGrid& a_eblgCoar,
		     const int&         a_nref,
		     const int&         a_nvar,
		     const bool         a_forceNoEBCF = false);

  /*!
    @brief Full constructor.  Calls the define function which creates
    a levels worth of flux registers.
  */
  EBFastFluxRegister(const DisjointBoxLayout& a_dblFine,
		     const DisjointBoxLayout& a_dblCoar,
		     const EBISLayout&        a_ebislFine,
		     const EBISLayout&        a_ebislCoar,
		     const Box&               a_domainCoar,
		     const int&               a_nref,
		     const int&               a_nvar,
		     const EBIndexSpace*      ebisPtr,
		     const bool               a_forceNoEBCF = false);

  /*!
    @brief Dtor
  */
  virtual ~EBFastFluxRegister();

  /*!
    @brief Full define function.  Creates a levels worth of flux registers.
    @details Calls other version
  */
  void define(const DisjointBoxLayout& a_dblFine,
	      const DisjointBoxLayout& a_dblCoar,
	      const EBISLayout&        a_ebislFine,
	      const EBISLayout&        a_ebislCoar,
	      const ProblemDomain&     a_domainCoar,
	      const int&               a_nref,
	      const int&               a_nvar,
	      const EBIndexSpace*      ebisPtr,
	      const bool               a_forceNoEBCF = false);

  /*!
    @brief Define function.
  */
  void define(const EBLevelGrid&       a_eblgFine,
	      const EBLevelGrid&       a_eblgCoar,
	      const int&               a_refRat,
	      const int&               a_nvar,
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
  
  void operator= (const EBFastFluxRegister& out) = delete;
  
  EBFastFluxRegister(const EBFastFluxRegister& out) = delete;
};

#endif

