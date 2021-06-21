/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFHelmholtzOp.H
  @brief  Declaration of a factory class for making multifluid Helmholtz operators for multigrid. 
  @author Robert Marskar
*/

#ifndef CD_MFHelmholtzOp_H
#define CD_MFHelmholtzOp_H

// Chombo includes
#include <MFCellFAB.H>
#include <MFFluxFAB.H>

// Our includes
#include <CD_MFBaseIVFAB.H>
#include <CD_EBHelmholtzOp.H>
#include <CD_EBHelmholtzDomainBCFactory.H>
#include <CD_EBHelmholtzEBBCFactory.H>
#include <CD_NamespaceHeader.H>

class MFHelmholtzOp : public AMRLevelOp<LevelData<MFCellFAB> > {
public:

  /*!
    @brief Constructor. Must subsequently call define(...)
  */
  MFHelmholtzOp();

  /*!
    @brief Destructor. 
  */
  ~MFHelmholtzOp();

  /*!
    @brief Define function. Not implemented (yet). 
  */
  void define();

  /*!
    @brief Return coarsening factor to coarser level (1 if there is no coarser level);
  */
  int refToCoarser() override final;

  /*!
    @brief Weird but true, coarsening between VCycles give better convergence
  */
  unsigned int orderOfAccuracy(void) const override final;

  /*!
    @brief Enforce CF consistency
  */
  void enforceCFConsistency(LevelData<MFCellFAB>& a_coarCorr, const LevelData<MFCellFAB>& a_fineCorr) override final;

protected:
  
};

#include <CD_NamespaceFooter.H>

#endif
