/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MultigridInterpolator.cpp
  @brief  Implementation of CD_MultigridInterpolator.H
  @author Robert Marskar
*/

// Our includes
#include <CD_MultigridInterpolator.H>
#include <CD_NamespaceHeader.H>

MultigridInterpolator::MultigridInterpolator(){
  m_isDefined = false;
}

MultigridInterpolator::MultigridInterpolator(const EBLevelGrid&            a_eblgFine,
					     const EBLevelGrid&            a_eblgCoar,
					     const int                     a_nRef,
					     const int                     a_nVar,
					     const LayoutData<IntVectSet>& a_ghostCells)
  : EBQuadCFInterp(a_eblgFine.getDBL(),
		   a_eblgCoar.getDBL(),
		   a_eblgFine.getEBISL(),
		   a_eblgCoar.getEBISL(),
		   a_eblgCoar.getDomain(),
		   a_nRef,
		   a_nVar,
		   a_ghostCells,
		   a_eblgFine.getEBIS(),
		   true){

  
}

MultigridInterpolator::~MultigridInterpolator(){

}

#include <CD_NamespaceFooter.H>
