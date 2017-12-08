/*!
  @file jump_bc.H
  @brief Declaration of a multifluid boundary condition class for the poisson equation
  @author Robert Marskar
  @date Dec. 2017
*/

#ifndef _JUMP_BC_
#define _JUMP_BC_

#include "MFLevelGrid.H"


/*!
  @brief Class for deriving the surface potential for multifluid poisson equations
*/
class jump_bc {
public:

  /*!
    @brief Empty constructor, must call define
  */
  jump_bc(){}

  /*!
    @brief Destructor
  */
  virtual ~jump_bc(){}

  /*!
    @brief Define function
  */
  virtual void define(const ProblemDomain& a_domain, const MFLevelGrid& a_mflg, const Real& a_dx){}

  /*!
    @brief Build stencils for matching
  */
  virtual void build_stencils(){}

protected:

  /*!
    @brief B-coefficient on each side
  */
  MFBaseIVFAB m_bco;

  /*!
    @brief Weights at the boundary
  */
  MFBaseIVFAB m_wb;
    
};

#endif
