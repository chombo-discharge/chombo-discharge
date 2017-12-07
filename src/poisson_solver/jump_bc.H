/*!
  @file jump_bc.H
  @brief Declaration of a multifluid boundary condition class for the poisson equation
  @author Robert Marskar
  @date Dec. 2017
*/

#ifndef _JUMP_BC_
#define _JUMP_BC_

#include <BaseEBBC.H>


/*!
  @brief Class for switching between multifluid-based Dirichlet and pure Dirichelt boundary conditions
*/
class jump_bc : public ConductivityBaseEBBC {
public:

  jump_bc(){

  }

  ~jump_bc(){

  }
  
protected:
  
};

#endif
