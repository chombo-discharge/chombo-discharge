/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_SimpleItoParticleParticle.H
  @brief  Declaration of a particle which only holds position, weight, and energy. 
  @author Robert Marskar
*/

#ifndef CD_SimpleItoParticleParticle_H
#define CD_SimpleItoParticleParticle_H

// Std includes
#include <string>
#include <typeinfo>
#include <cstdio>
#include <cmath>

// Our includes
#include <CD_GenericParticle.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief A particle class for reducing HDF5 I/O file sizes when using ItoSolver. 
  @details This is a computational particle class which only stores the weight, position and energy of a particle. It is used
  with ItoSolver in order to reduce the file size of HDF5 files. 
*/
class SimpleItoParticle : public GenericParticle<2, 0>
{
public:
  /*!
    @brief Default constructor. Should subsequently set weight/position/energy.
  */
  inline SimpleItoParticle();

  /*!
    @brief Constructor. Sets weight, position, and energy.
    @param[in] a_weight   Computational weight
    @param[in] a_position Position in physical space. 
    @param[in] a_energy   Average particle energy
  */
  inline SimpleItoParticle(const Real a_weight, const RealVect a_position, const Real a_energy);

  /*!
    @brief Destructor (does nothing)
  */
  inline virtual ~SimpleItoParticle();

  /*!
    @brief Define function. Sets weight, position, and energy.
    @param[in] a_weight   Computational weigth
    @param[in] a_position Position in physical space. 
    @param[in] a_energy   Average particle energy
  */
  inline void
  define(const Real a_weight, const RealVect a_position, const Real a_energy);

  /*!
    @brief Get weight
    @return Returns m_weight
  */
  inline Real&
  weight();

  /*! 
    @brief Const retrieve weight
    @return Returns m_weight
  */
  inline const Real&
  weight() const;

  /*!
    @brief Get energy
    @return Returns m_energy
  */
  inline Real&
  energy();

  /*! 
    @brief Const retrieve energy
    @return Returns m_energy
  */
  inline const Real&
  energy() const;
};

#include <CD_NamespaceFooter.H>

#include <CD_SimpleItoParticleImplem.H>

#endif
