/* chombo-discharge
 * Copyright © 2025 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_Rod.H
  @brief  File containing headers for a rod-like geometry.
  @author Robert Marskar
*/

#ifndef CD_Rod_H
#define CD_Rod_H

// Our includes
#include <CD_ComputationalGeometry.H>
#include <CD_EBGeometryIF.H>
#include <CD_NamespaceHeader.H>

class Rod : public ComputationalGeometry
{
public:
  /*!
    @brief Full constructor.
    @details This function parses all input parameters and creates the geometry
  */
  Rod() noexcept;

  /*!
    @brief Destructor (does nothing)
  */
  virtual ~Rod() noexcept = default;
};

#include <CD_NamespaceFooter.H>

#endif
