/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_GECReferenceCell.H
  @brief  Simplified GEC reference cell geometry.
  @author Robert Marskar
*/

#ifndef CD_GECReferenceCell_H
#define CD_GECReferenceCell_H

// our includes
#include <CD_ComputationalGeometry.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Simplified GEC reference cell geometry.
*/
class GECReferenceCell : public ComputationalGeometry
{
public:
  /*!
    @brief Constructor
  */
  GECReferenceCell();

  /*!
    @brief Destructor
  */
  ~GECReferenceCell();
};

#include <CD_NamespaceFooter.H>

#endif
