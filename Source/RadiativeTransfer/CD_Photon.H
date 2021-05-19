/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_Photon.H
  @brief  Declaration of a photon class for particle methods
  @author Robert Marskar
*/

#ifndef CD_Photon_H
#define CD_Photon_H

// Std includes
#include <string>
#include <typeinfo>
#include <cstdio>
#include <cmath>

// Chombo includes
#include <BinItem.H>

// Our includes
#include <CD_NamespaceHeader.H>
  
/*!
  @brief Photon class
*/
class Photon : public BinItem {
public:

  Photon();
  Photon(const RealVect& a_position, const RealVect& a_velocity, const Real& a_kappa, const Real a_mass = 1.0);

  virtual ~Photon();

  void define(const RealVect& a_position, const RealVect& a_velocity, const Real& a_kappa, const Real a_mass = 1.0);

  void setKappa(const Real a_kappa);
  void setMass(const Real a_mass);
  void setVelocity(const RealVect& a_velocity);

  Real&     mass();
  Real&     kappa();
  RealVect& velocity();

  const Real&     kappa()    const;
  const Real&     mass()   const;
  const RealVect& velocity() const;

  virtual bool operator==(const Photon& a_p) const;
  virtual bool operator==(const Photon* a_p) const;
  virtual bool operator!=(const Photon& a_p) const;

  virtual int  size() const;
  virtual void linearOut(void* buf) const;
  virtual void linearIn(void* buf);

protected:

  Real     m_kappa;
  Real     m_mass;
  RealVect m_velocity;
};

#include <CD_NamespaceFooter.H>

#endif 
