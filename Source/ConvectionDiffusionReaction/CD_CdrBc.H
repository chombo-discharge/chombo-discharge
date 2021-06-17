/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrBc.H
  @brief  Declaration of boundary condition enum for cdr solver. 
  @author Robert Marskar
*/

#ifndef CD_CdrBc_H
#define CD_CdrBc_H

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Boundary condition types for the CDR solvers
*/
enum class CdrBc {
  Data,
  Wall,
  Function,
  Outflow,
  Solver
};

#include <CD_NamespaceFooter.H>

#endif
