/*!
  @file domain_bc.H
  @brief Declaration of class for setting boundary conditions on individual faces
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _DOMAIN_BC_
#define _DOMAIN_BC_


#include <BaseBCValue.H>

#include "wall_bc.H"

/*!
  @brief Class for setting Dirichlet or Neumann conditions on individual faces
*/
class domain_bc : public BaseBCValue {
public:

  /*!
    @brief Full constructor. Subsequently calls define
  */
  domain_bc(Vector<RefCountedPtr<wall_bc> >& a_bc);

  /*!
    @brief Domain BC function
  */
  virtual ~domain_bc();

  /*!
    @brief Set potential. 
  */
  virtual void set_potential_ptr(Real (*a_ptr)(const Real a_time));

  /*!
    @brief Define function
  */
  virtual void define(Vector<RefCountedPtr<wall_bc> >& a_bc);

  /*!
    @brief Value function
  */
  virtual Real value(const RealVect& a_point, const RealVect& a_normal, const Real& a_time, const int& a_comp) const;

  /*!
    @brief Value function
  */
  virtual Real value(const FaceIndex&      a_face, 
		     const Side::LoHiSide& a_side, 
		     const DataIndex&      a_dit, 
		     const RealVect&       a_point,
		     const RealVect&       a_normal, 
		     const Real&           a_time, 
		     const int&            a_comp) const;

protected:

  /*!
    @brief Domain bc function
  */
  domain_bc();

  /*!
    @brief Potential 
  */
  Real (*m_ptr)(const Real a_time);

  /*!
    @brief Boundary conditions
  */
  Vector<RefCountedPtr<wall_bc> > m_bc;
  
};

#endif
