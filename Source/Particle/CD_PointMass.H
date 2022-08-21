/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_PointMass.H
  @brief  Declaration of a simple class for storing a point with mass.
  @author Robert Marskar
*/

#ifndef CD_PointMass_H
#define CD_PointMass_H

// Chombo includes
#include <RealVect.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Declaration of a simple class for storing a point with mass. 
  @details This is not a computational particle.
*/
class PointMass
{
public:
  /*!
    @brief Constructor. Should subsequently set parameters or call define
  */
  inline PointMass();

  /*!
    @brief Full constructor. 
    @param[in] a_pos    Physical position
    @param[in] a_mass   Computational weight
    @param[in] a_energy Particle energy
  */
  inline PointMass(const RealVect a_pos, const Real a_mass, const Real a_energy);

  /*!
    @brief Alternative constructor which computes an "average" PointMass.
    @details This computes the average position as x = sum_i(m_i*x_i)/sum(m_i) etc. 
    @param[in] a_pointMasses list of point masses.
  */
  inline PointMass(const std::vector<PointMass>& a_pointMasses);

  /*!
    @brief Destructor (does nothing)
  */
  inline ~PointMass();

  /*!
    @brief Define function. 
    @param[in] a_pos    Physical position
    @param[in] a_mass   Computational weight
    @param[in] a_energy Particle energy
  */
  inline void
  define(const RealVect& a_pos, const Real& a_mass, const Real& a_energy);

  /*!
    @brief Position operator. Get the physical position of the input coordinate
    @param[in] a_dir Coordinate
    @return Returns m_pos[a_dir]
  */
  inline Real&
  operator[](const int a_dir);

  /*!
    @brief Position operator. Get the physical position of the input coordinate
    @param[in] a_dir Coordinate
    @return Returns m_pos[a_dir]
  */
  inline const Real&
  operator[](int a_dir) const;

  /*!
    @brief Get physical position
    @return Returns m_pos
  */
  inline RealVect&
  pos();

  /*!
    @brief Get physical position
    @return Returns m_pos
  */
  inline const RealVect&
  pos() const;

  /*!
    @brief Get mass (i.e., computational weight)
    @return Returns m_mass
  */
  inline Real&
  mass();

  /*!
    @brief Get mass (i.e., computational weight)
    @return Returns m_mass
  */
  inline const Real&
  mass() const;

  /*!
    @brief Get energy
    @return Returns m_energy
  */
  inline Real&
  energy();

  /*!
    @brief Get energy
    @return Returns m_energy
  */
  inline const Real&
  energy() const;

  /*!
    @brief Return true if this PointMass can be split into finer PointMass-es.
  */
  inline bool
  canSplit() const;

  /*!
    @brief Physical position of this PointMass
  */
  RealVect m_pos;

  /*!
    @brief Computational weight of this PointMass
  */
  Real m_mass;

  /*!
    @brief Energy of this PointMass
  */
  Real m_energy;
};

#include <CD_NamespaceFooter.H>

#include <CD_PointMassImplem.H>

#endif
