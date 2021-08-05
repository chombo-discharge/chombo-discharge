/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFHelmholtzElectrostaticEBBCImplem.H
  @brief  Implementation of CD_MFHelmholtzElectrostaticEBBC.H
  @author Robert Marskar
*/

#ifndef CD_MFHelmholtzElectrostaticEBBCImplem_H
#define CD_MFHelmholtzElectrostaticEBBCImplem_H

//
#include <limits>

// Our includes
#include <CD_MFHelmholtzElectrostaticEBBC.H>
#include <CD_NamespaceHeader.H>

Real MFHelmholtzElectrostaticEBBC::getElectrodePotential(const RealVect& a_pos) const {


  // Find closest electrode
  int closestElectrode;
  Real minDist = std::numeric_limits<Real>::infinity();

  const std::vector<std::pair<Electrode, ElectrostaticEbBc::BcFunction> >& electrodeBCs = m_electrostaticBCs.getBcs();
  
  for (int i = 0; i < electrodeBCs.size(); i++){
    const RefCountedPtr<BaseIF>& impFunc = electrodeBCs[i].first.getImplicitFunction();

    const Real curDist = impFunc->value(a_pos);

    if(std::abs(curDist) < minDist){
      closestElectrode = i;
      minDist          = curDist;
    }
  }

  // Return potential of closest electrode
  const Real dummyDt = 0.0;
  
  return electrodeBCs[closestElectrode].second(a_pos, dummyDt);
}

#include <CD_NamespaceFooter.H>

#endif
