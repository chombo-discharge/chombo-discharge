/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MultigridInterpolator.H
  @brief  Declaration of a factory class for interpolating ghost cells in a multigrid context. 
  @author Robert Marskar
*/

#ifndef CD_MultigridInterpolator_H
#define CD_MultigridInterpolator_H

// Chombo includes
#include <EBLevelGrid.H>
#include <EBQuadCFInterp.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Multigrid interpolator class. 
  @note Currently just a shell of EBQuadCFInterp (so we have something to work with). 
*/
class MultigridInterpolator : public EBQuadCFInterp {
public:

  /*!
    @brief Default constructor. Must subsequently call define
  */
  MultigridInterpolator();

  /*!
    @brief Disallowed copy constructor.
  */
  MultigridInterpolator(const MultigridInterpolator& a_other) = delete;

  /*!
    @brief Full constructor
  */
  MultigridInterpolator(const EBLevelGrid&            a_eblgFine,
			const EBLevelGrid&            a_eblgCoar,
			const int                     a_nRef,
			const int                     a_nVar,
			const LayoutData<IntVectSet>& a_ghostCells);


  /*!
    @brief Dtor. Releases stencils. 
  */
  ~MultigridInterpolator();

protected:
};

#include <CD_NamespaceFooter.H>

#endif
