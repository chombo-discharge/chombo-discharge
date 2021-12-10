/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_Vessel.H
  @brief  High-voltage vessel geometry
  @author Robert Marskar
*/

#ifndef CD_Vessel_H
#define CD_Vessel_H

// our includes
#include <CD_ComputationalGeometry.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief High-voltage vessel geometry
  @details This consists of a rod-plane gap. The "plane" is an elevated plate. In 2D the geometry is oriented along +y and in 3D it is oriented
  along +z. 

  To use this class, simply call the constructor (it will fill all internal data from the input script).
*/
class Vessel : public ComputationalGeometry {
public:

  /*!
    @brief Constructor
  */
  Vessel();

  /*!
    @brief Destructor
  */  
  ~Vessel();
  
protected:

  /*!
    @brief Center point for the rod
  */
  RealVect m_rod_center;

  /*!
    @brief Center point for the elevated plate
  */
  RealVect m_shroom_center;

  /*!
    @brief Rod radius
  */
  Real m_rod_R;

  /*!
    @brief Elevated plate radius
  */
  Real m_shroom_R;

  /*!
    @brief Elevated plate support radius
  */
  Real m_shroom_r;

  /*!
    @brief Elevated plate thickness
  */
  Real m_shroom_d;

  /*!
    @brief Elevated plate corner curvature
  */
  Real m_shroom_c;

  /*!
    @brief If true, rod is at live voltage
  */
  bool m_rod_live;

  /*!
    @brief If true, plate is at live voltage
  */  
  bool m_shroom_live;

  /*!
    @brief If false, elevated plate is removed from the geometry
  */
  bool m_shroom;

  /*!
    @brief If false, rod is removed from the geometry
  */  
  bool m_rod;
};

#include <CD_NamespaceFooter.H>

#endif
