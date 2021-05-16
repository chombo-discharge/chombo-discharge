/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   cdr_bc.H
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
namespace cdr_bc {
  enum which_bc {
    external = 0,
    wall     = 1,
    outflow  = 2,
    extrap   = 3,
  };
};

#include <CD_NamespaceFooter.H>

#endif
