/*!
  @file mf_helmholtz_op.H
  @brief Declaration of multifluid helmholtz operator for embedment in AMRMultiGrid
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _MF_HELMHOLTZ_OP_
#define _MF_HELMHOLTZ_OP_

/*!
  @brief Operator class for solving the variable-coefficient Helmholtz equation with multigrid in two phases
*/
class mf_helmholtz_op {
public:

  mf_helmholtz_op();

  virtual ~mf_helmholtz_op();

protected:

private:
};
#endif
