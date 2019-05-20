/*!
  @file   photon.H
  @brief  Declaration of a photon class for particle methods
  @author Robert Marskar
  @date   May 2019
*/

#ifndef _PHOTON_H_
#define _PHOTON_H_

#include <string>
#include <typeinfo>
#include <cstdio>
#include <cmath>

#include "BinItem.H"

/*!
  @brief photon class
*/
class photon : public BinItem {
public:

  photon();
  photon(const RealVect& a_position, const RealVect& a_velocity, const Real& a_kappa, const Real a_mass = 1.0);

  virtual ~photon();

  void define(const RealVect& a_position, const RealVect& a_velocity, const Real& a_kappa, const Real a_mass = 1.0);

  void set_kappa(const Real a_kappa);
  void set_mass(const Real a_mass);
  void set_position(const RealVect& a_position);
  void set_velocity(const RealVect& a_velocity);


  Real&     mass();
  Real&     kappa();
  RealVect& position();
  RealVect& velocity();

  const Real&     kappa()    const;
  const Real&     mass()   const;
  const RealVect& position() const ;
  const RealVect& velocity() const;

  virtual bool operator==(const photon& a_p) const;
  virtual bool operator==(const photon* a_p) const;
  virtual bool operator!=(const photon& a_p) const;

  virtual int  size() const;
  virtual void linearOut(void* buf) const;
  virtual void linearIn(void* buf);

protected:

  Real     m_kappa;
  Real     m_mass;
  RealVect m_velocity;
};


#endif 
