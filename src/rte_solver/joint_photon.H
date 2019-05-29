/*!
  @file   joint_photon.H
  @brief  Declaration of a joint_photon class for particle methods
  @author Robert Marskar
  @date   May 2019
*/

#ifndef _JOINT_PHOTON_H_
#define _JOINT_PHOTON_H_

#include <string>
#include <typeinfo>
#include <cstdio>
#include <cmath>

#include "BinItem.H"

/*!
  @brief   joint_photon class
  @details This class is used for deposition of photons on a finer AMR level in order to correctly capture fine-level clouds that
           hand into a coarse level
*/
class joint_photon : public BinItem {
public:

  joint_photon();
  joint_photon(const Real& a_mass, const RealVect& a_position, const size_t& a_num_photons = 0);

  virtual ~joint_photon();

  void define(const Real& a_mass, const RealVect& a_position, const size_t& a_num_photons = 0);
  void add_photon(const photon* const a_photon);
  void clear();

  size_t&   num_photons();
  Real&     mass();
  RealVect& position();

  const size_t&   num_photons() const;
  const Real&     mass()          const;
  const RealVect& position()      const;

  virtual bool operator==(const joint_photon& a_p) const;
  virtual bool operator==(const joint_photon* a_p) const;
  virtual bool operator!=(const joint_photon& a_p) const;

  virtual int  size() const;
  virtual void linearOut(void* buf) const;
  virtual void linearIn(void* buf);

protected:

  Real     m_mass;
  size_t   m_num_photons;
};


#endif 
