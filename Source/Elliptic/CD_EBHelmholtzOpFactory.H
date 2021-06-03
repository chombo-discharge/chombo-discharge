/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */
/*!
  @file   CD_EBHelmholtzOpFactory.H
  @brief  Declaration of a factory class for making Poisson operators for multigrid. 
  @author Robert Marskar
*/

#ifndef CD_EBHelmholtzOpFactory_H
#define CD_EBHelmholtzOpFactory_H

// Chombo includes
#include <RealVect.H>
#include <EBLevelGrid.H>
#include <EBFluxRegister.H>
#include <EBQuadCFInterp.H>

// Our includes
#include <CD_EBHelmholtzOp.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Factory class for making variable-coefficient Helmholtz operators. 
*/
class EBHelmholtzOpFactory : public AMRLevelOpFactory<LevelData<EBCellFAB> > {
public:

  /*!
    @brief Disallowed constructor. Use the full constructor
  */
  EBHelmholtzOpFactory() = delete;

  /*!
    @brief Disallowed constructor. Use the full constructor
  */
  EBHelmholtzOpFactory(const EBHelmholtzOpFactory& a_otherFactory) = delete;

  /*!
    @brief Destructor. Does nothing.
  */
  ~EBHelmholtzOpFactory();

  /*!
    @brief Disallowed assignment operator
  */
  void operator=(const EBHelmholtzOpFactory& a_opin) = delete;

  /*!
    @brief Create multigrid operator
    @param[in] a_fineDomain      Domain
    @param[in] a_depth           Depth. This specifies that the operator will be created at depth coarsen(a_fineDomain, 2^a_depth);
    @param[in] a_homogeneousOnly If true, only homogeneous boundary conditions will be needed. 
  */
  EBHelmholtzOp* MGnewOp(const ProblemDomain& a_fineDomain, int a_depth, bool a_homogeneousOnly = true) override final;

  /*!
    @brief Create AMR operator for specified domain
    @param[in] a_domain Domain
  */
  EBHelmholtzOp* AMRnewOp(const ProblemDomain& a_domain) override final;

  /*!
    @brief Get refinement ratio to next finest level.
    @details Return 1 when there are not finer levels. 
  */
  int refToFiner(const ProblemDomain& a_indexspace) const override final;

protected:

};

#include <CD_NamespaceFooter.H>

#endif
