/*!
  @file box_if.H
  @brief Declare an interface for defining a two- or three-dimensional box with arbitrary centroid and orientation.
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _BOX_IF_
#define _BOX_IF_

#include <BaseIF.H>

#include "CD_NamespaceHeader.H"
  
/*!
  @brief Class for defining a two- or three-dimensional box with arbitrary centroid and orientation.
*/
class box_if: public BaseIF {
public:

  /*!
    @brief Base constructor
    @param[in] a_loCorner Lower left corner
    @param[in] a_hiCorner Upper right corner
    @param[in] a_inside Object on inside or outside
  */
  box_if(const RealVect& a_loCorner, 
	 const RealVect& a_hiCorner,
	 const bool&     a_fluidInside);
    
  /*!
    @brief Copy constructor
    @param[in] a_inputIF Another box_if
  */
  box_if(const box_if& a_inputIF);

  /*!
    @brief Destructor
  */
  virtual ~box_if();

  /*!
    @brief Test if object is on inside or outside
  */
  virtual Real value(const RealVect& a_pos) const;

  /*!
    @brief IF Factory method
  */
  virtual BaseIF* newImplicitFunction() const;

protected:
  
  /*!
    @brief Lower left corner
  */
  RealVect m_loCorner;

  /*!
    @brief Upper right corner
  */
  RealVect m_hiCorner;

  /*!
    @brief Fluid inside box or not. 
  */
  bool m_fluidInside;
};
#include "CD_NamespaceFooter.H"

#endif
