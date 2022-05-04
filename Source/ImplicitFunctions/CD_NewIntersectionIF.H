/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_NewIntersectionIF.H
  @brief  Declaration of a signed distance function which is infinitely far away from everything. 
  @author Robert Marskar
*/

#ifndef CD_NewIntersectionIF_H
#define CD_NewIntersectionIF_H

// Chombo includes
#include <BaseIF.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief New intersection IF which does not mess up the return value function when there are no implicit functions.
*/
class NewIntersectionIF : public BaseIF
{
public:
  /*!
    @brief Weak constructor. Does nothing. 
  */
  NewIntersectionIF();

  /*!
    @brief Constructor taking several functions. 
    @param[in] a_impFuncs Implicit functions
  */
  NewIntersectionIF(const Vector<BaseIF*>& a_impFuncs);

  /*!
    @brief Destructor
  */
  virtual ~NewIntersectionIF();

  /*!
    @brief Get distance to objects. 
    @param[in] a_point Physical position. 
  */
  virtual Real
  value(const RealVect& a_point) const override;

  /*!
    @brief Factory method
  */
  virtual BaseIF*
  newImplicitFunction() const override;

protected:
  /*!
    @brief Number of implicit functions
  */
  int m_numFuncs; // number of implicit functions

  /*!
    @brief Implicit functions
  */
  Vector<BaseIF*> m_impFuncs;
};

#include <CD_NamespaceFooter.H>

#endif
