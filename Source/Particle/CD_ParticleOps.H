/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_ParticleOps.H
  @brief  Declaration of some common useful particle routines that would otherwise be replicated across classes. 
  @author Robert Marskar
*/

#ifndef CD_ParticleOps_H
#define CD_ParticleOps_H

// Our includes
#include <CD_NamespaceHeader.H>
  
class ParticleOps {
public:

  static
  inline bool domainIntersection(const RealVect& a_oldPos,
				 const RealVect& a_newPos,
				 const RealVect& a_path,
				 const RealVect& a_probLo,
				 const RealVect& a_prob_hi,
				 Real&           a_s);

  static
  inline bool ebIntersectionBisect(const RefCountedPtr<BaseIF>& a_impfunc,
				   const RealVect&              a_oldPos,
				   const RealVect&              a_newPos,
				   const Real&                  a_pathLen,
				   const Real&                  a_bisectStep,
				   Real&                        a_s);

  static
  inline bool ebIntersectionRaycast(const RefCountedPtr<BaseIF>& a_impfunc,
				    const RealVect&              a_x0,
				    const RealVect&              a_x1,
				    const Real&                  a_tolerance,
				    Real&                        a_s);
  
};

#include <CD_NamespaceFooter.H>

#include <CD_ParticleOpsImplem.H>

#endif
