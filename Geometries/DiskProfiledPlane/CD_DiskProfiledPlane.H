/* chombo-discharge
 * Copyright © SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @brief CD_DiskProfiledPlane.H
  @brief Declaration of a geometry describing a disk electrode and a profiled dielectric surface. 
  @author Robert Marskar
*/

#ifndef CD_DiskProfiledPlane_H
#define CD_DiskProfiledPlane_H

// Our includes
#include <CD_ComputationalGeometry.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Geometry class for disk electrode and profiled surface experiments. 
  @details This geometry consists of a disk electrode oriented along +y (user can shift it), and a profiled
  plane with specified profiles: squares or square channels, spheres, or cylinders. Profiled are organized
  along a user-specified lattice (coinciding with the coordinate directions), and can be shifted. 
*/
class DiskProfiledPlane : public ComputationalGeometry
{
public:
  /*!
    @brief Constructor - reads in all parameters
  */
  DiskProfiledPlane() noexcept;

  /*!
    @brief Destructor - does nothing
  */
  virtual ~DiskProfiledPlane() = default;

protected:
  /*!
    @brief Define the electrode
  */
  virtual void
  defineElectrode() noexcept;

  /*!
    @brief Define the dielectric plate.
  */
  virtual void
  defineDielectric() noexcept;
};

#include <CD_NamespaceFooter.H>

#endif
