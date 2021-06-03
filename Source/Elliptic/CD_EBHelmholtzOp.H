/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */
/*!
  @file   CD_EBHelmholtzOp.H
  @brief  Declaration of Helmholtz multigrid operators. 
  @author Robert Marskar
*/

#ifndef CD_EBHelmholtzOp_H
#define CD_EBHelmholtzOp_H

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Factory class for making variable-coefficient Helmholtz operators. This is meant to be used with LevelTGA so the operator can be used for TGA time stepping. 
*/
class EBHelmholtzOp : public LevelTGAHelmOp<LevelData<EBCellFAB>, EBFluxFAB> {
public:

  /*!
    @brief Disallowed default constructor
  */
  EBHelmholtzOp() = delete;

  /*!
    @brief Disallowed copy constructor
  */
  EBHelmholtzOp(const EBHelmholtzOp& a_other) = delete;

  /*!
    @brief Full constructor
  */
  EBHelmholtzOp(const int a){

  }

  
protected:

};

#include <CD_NamespaceFooter.H>

#endif
