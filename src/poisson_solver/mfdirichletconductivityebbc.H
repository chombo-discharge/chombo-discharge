/*!
  @file mfdirichletconductivityebbc.H
  @brief Declaration of a boundary condition class for doing convential Dirichlet BCs as well as multifluid matching bcs 
  on selected cells
  @author Robert Marskar
  @date Dec. 2017
*/

#ifndef _MFDIRICHLETCONDUCTIVITYEBBC_
#define _MFDIRICHLETCONDUCTIVITYEBBC_

#include <DirichletConductivityEBBC.H>

/*!
  @brief Class for doing regular Dirichlet boundary conditions on EBConductivityOps with selected multifluid matching 
  boundary conditions.
  @details This is a simple class that essentially just ignores calls to homogeneous BCs on interface cells. 
*/
class mfdirichletconductivityebbc : public DirichletConductivityEBBC {
public:
  /*!
    @brief Constructor. Implements DirichletConductivityEBBC
  */
  mfdirichletconductivityebbc(const ProblemDomain& a_domain,
			      const EBISLayout&    a_ebisl,
			      const RealVect&      a_dx,
			      const IntVect*       a_ghost_phi,
			      const IntVect*       a_ghost_rhs) : DirichletConductivityEBBC(a_domain,
											    a_ebisl,
											    a_dx,
											    a_ghost_phi,
											    a_ghost_rhs){

  }
  
  /*!
    @brief Destructor
  */
  virtual ~mfdirichletconductivityebbc();

};
#endif
