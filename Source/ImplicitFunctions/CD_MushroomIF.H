/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file    CD_MushroomIF.H
  @brief   Declaration of a "mushroom", used in HV Vessel
  @author  Robert Marskar
*/

#ifndef CD_MushroomIF_H
#define CD_MushroomIF_H

// Chombo includes
#include <BaseIF.H>
#include <RealVect.H>

// Our includes
#include <CD_NamespaceHeader.H>

class MushroomIF : public BaseIF {
public:  

  MushroomIF(const RealVect a_center,
	      const Real     a_R,
	      const Real     a_r,
	      const Real     a_L,
	      const Real     a_d,
	      const Real     a_curv,
	      const bool     a_fluidInside);

  MushroomIF(const MushroomIF& a_inputIF);
  virtual ~MushroomIF();
  virtual Real value(const RealVect& a_pos) const override;
  virtual BaseIF* newImplicitFunction() const override;
  
protected:

  bool m_fluidInside;
  RefCountedPtr<BaseIF> m_baseif;

    
};

#include <CD_NamespaceFooter.H>

#endif
