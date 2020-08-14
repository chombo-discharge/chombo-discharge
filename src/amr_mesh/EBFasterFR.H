/*!
  @brief   EBFasterFR.H
  @details Declaration of a new flux register class with optimizations. 
  @author  Robert Maskar
  @date    Aug. 2020
*/

#ifndef _EBFasterFR_H_
#define _EBFasterFR_H_

#include "EBFastFR.H"

/*!
  @brief Class that optimizes some EB-related operations of EBFastFR
*/

class EBFasterFR : public EBFastFR {
public:
  
  /*!
    @brief Default constructor. Leaves object undefined.
  */
  EBFasterFR();

  /*!
    @brief Full constructor. Class the define function. 
  */
  EBFasterFR(const EBLevelGrid&       a_eblgFine,
	     const EBLevelGrid&       a_eblgCoar,
	     const int&               a_nref,
	     const int&               a_nvar,
	     const bool               a_forceNoEBCF = false);

  /*!
    @brief Dtor
  */
  virtual ~EBFasterFR();

  /*
    @brief Full define function.  Creates a levels worth of flux registers.
  */
  void define(const EBLevelGrid&       a_eblgFine,
	      const EBLevelGrid&       a_eblgCoar,
	      const int&               a_nref,
	      const int&               a_nvar,
	      const bool               a_forceNoEBCF = false);

private:
  
  void operator= (const EBFasterFR& out){
    MayDay::Error("EBFasterFR::invalid operator");
  }
  
  EBFasterFR(const EBFasterFR& out){
    MayDay::Error("EBFasterFR::invalid operator");
  }
};

#endif
