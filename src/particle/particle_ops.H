/*!
  @file   particle_ops.H
  @brief  Declaration of some common useful particle routines that would otherwise be replicated across classes. 
  @author Robert Marskar
  @date   June 2020
*/

#ifndef _PARTICLE_OPS_H
#define _PARTICLE_OPS_H

class particle_ops {
public:

  static inline bool domain_bc_intersection(const RealVect& a_oldPos,
					    const RealVect& a_newPos,
					    const RealVect& a_path,
					    const RealVect& a_prob_lo,
					    const RealVect& a_prob_hi,
					    Real&           a_s);

  static inline bool eb_bc_intersection(const RefCountedPtr<BaseIF>& a_impfunc,
					const RealVect&              a_oldPos,
					const RealVect&              a_newPos,
					const Real&                  a_pathLen,
					const Real&                  a_bisectStep,
					Real&                        a_s);
  
};

#include "particle_opsI.H"
#endif
