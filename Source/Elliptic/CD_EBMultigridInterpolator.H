/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBMultigridInterpolator.H
  @brief  Declaration of a factory class for interpolating ghost cells in a multigrid context. 
  @author Robert Marskar
*/

#ifndef CD_EBMultigridInterpolator_H
#define CD_EBMultigridInterpolator_H

// Chombo includes
#include <EBLevelGrid.H>
#include <EBQuadCFInterp.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Multigrid interpolator class. 
  @note Currently just a shell of EBQuadCFInterp (so we have something to work with). 
*/
class EBMultigridInterpolator : public EBQuadCFInterp {
public:

  /*!
    @brief Default constructor. Must subsequently call define
  */
  EBMultigridInterpolator();

  /*!
    @brief Disallowed copy constructor.
  */
  EBMultigridInterpolator(const EBMultigridInterpolator& a_other) = delete;

  /*!
    @brief Full constructor
  */
  EBMultigridInterpolator(const EBLevelGrid&            a_eblgFine,
			  const EBLevelGrid&            a_eblgCoar,
			  const int                     a_nRef,
			  const int                     a_nVar,
			  const LayoutData<IntVectSet>& a_ghostCells);


  /*!
    @brief Dtor. Releases stencils. 
  */
  ~EBMultigridInterpolator();

protected:
};

#include <CD_NamespaceFooter.H>

#endif
