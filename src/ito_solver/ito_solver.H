/*!
  @file   ito_solver.H
  @brief  Declaration of an abstract class for Ito diffusion
  @author Robert Marskar
  @date   April 2020
*/

#ifndef _ITO_SOLVER_
#define _ITO_SOLVER_

/*!
  @brief Base class for Ito diffusion particle models
*/
class ito_solver {
public:

  ito_solver();

  virtual ~ito_solver();

protected:
};
#endif
