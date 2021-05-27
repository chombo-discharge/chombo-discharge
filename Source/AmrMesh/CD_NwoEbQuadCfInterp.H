/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_NwoEbQuadCfInterp.H
  @brief  This is just a copy of Chombos NWOEBQuadCFInterp which also takes an ebis pointer argument. 
  @author Robert Marskar
*/

#ifndef CD_NwoEbQuadCfInterp_H
#define CD_NwoEbQuadCfInterp_H

// Chombo includes
#include <NWOEBQuadCFInterp.H>

// Our includes
#include <CD_NamespaceHeader.H>
  
/*!
  @brief Reimplementation of NWOEBQuadCFInterp which takes an additional argument; the EBIndexSpace
*/
class NwoEbQuadCfInterp : public NWOEBQuadCFInterp {
public:

  /*!
    @brief Don't see why anyone would need this constructor
  */
  NwoEbQuadCfInterp() = delete;

  /*!
    @brief Don't see why anyone would need this constructor
  */
  //  NwoEbQuadCfInterp(NwoEbQuadCfInterp& a_input) = delete;

  /*!
    @brief Constructor
    @param[in] a_gridsFine  Fine grids
    @param[in] a_gridsCoar  Coarse grids
    @param[in] a_ebislFine  Fine EBISLayout
    @param[in] a_ebislCoar  Coarse EBISLayout
    @param[in] a_domainCoar Coarse level domainw
    @param[in] a_nref       Refinement ratio
    @param[in] a_nvar       Number of variables
    @param[in] a_dxFine     Grid resolution on fine level
    @param[in] a_ghost      Ghost cells. 
    @param[in] a_cfivs      Coarse-fine interface
    @param[in] a_ebisPtr    EBIndexSpace pointer. 
   */
  NwoEbQuadCfInterp(const DisjointBoxLayout&           a_gridsFine,
		    const DisjointBoxLayout&           a_gridsCoar,
		    const EBISLayout&                  a_ebislFine,
		    const EBISLayout&                  a_ebislCoar,
		    const ProblemDomain&               a_domainCoar,
		    const int&                         a_nref,
		    const int&                         a_nvar,
		    const Real&                        a_dxFine,
		    const int&                         a_ghost,
		    const LayoutData<IntVectSet>&      a_cfivs,
		    const RefCountedPtr<EBIndexSpace>& a_ebisPtr);

  /*!
    @brief Destructor
  */
  ~NwoEbQuadCfInterp();

  /*!
    @brief Nope, no assignement.
  */
  void operator=(NwoEbQuadCfInterp& a_input) = delete;

protected:
  
  void defineInternals(const LayoutData<IntVectSet>& a_cfivs, const RefCountedPtr<EBIndexSpace>& a_ebisPtr);
  
};

#include <CD_NamespaceFooter.H>

#endif
