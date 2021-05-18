/*!
  @file   particle_ops.H
  @brief  Declaration of some common useful particle routines that would otherwise be replicated across classes. 
  @author Robert Marskar
  @date   June 2020
*/

#ifndef _PARTICLE_OPSI_H
#define _PARTICLE_OPSI_H

#include "particle_ops.H"
#include "poly.H"

#include <PolyGeom.H>

#include "CD_NamespaceHeader.H"

inline
bool particle_ops::domain_intersection(const RealVect& a_oldPos,
				       const RealVect& a_newPos,
				       const RealVect& a_path,
				       const RealVect& a_probLo,
				       const RealVect& a_prob_hi,
				       Real&           a_s){

  // TLDR: This code does a boundary intersection test and returns where on the interval [oldPos, newPos] the intersection
  //       happened.

  a_s = 1.E99;
  bool retval = false;

  for (int dir = 0; dir < SpaceDim; dir++){
    for (SideIterator sit; sit.ok(); ++sit){
      const Side::LoHiSide side = sit();
      const RealVect wallPoint  = (side == Side::Lo) ? a_probLo : a_prob_hi; // A point on the domain side
      const RealVect n0         = sign(side)*RealVect(BASISV(dir));           // Normal vector pointing OUT of the domain
      const Real norm_path      = PolyGeom::dot(n0, a_path);                  // Component relative to wall

      if(norm_path > 0.0){ 
	const Real s = PolyGeom::dot(wallPoint-a_oldPos, n0)/norm_path;  
	if(s >= 0.0 && s <= 1.0){
	  retval        = true;
	  if(s < a_s){
	    a_s = s;
	  }
	}
      }
    }
  }

  return retval;
}


inline
bool particle_ops::eb_intersection_bisect(const RefCountedPtr<BaseIF>& a_impfunc,
					  const RealVect&              a_oldPos,
					  const RealVect&              a_newPos,
					  const Real&                  a_pathLen,
					  const Real&                  a_bisect_step,
					  Real&                        a_s){

  bool retval = false;
  
  const int nsteps      = ceil(a_pathLen/a_bisect_step);
  const RealVect dxStep = (a_newPos - a_oldPos)/nsteps;
	    
  // Check each interval
  RealVect cur_pos  = a_oldPos;
  for (int istep = 0; istep < nsteps; istep++){
    const Real fa = a_impfunc->value(cur_pos);
    const Real fb = a_impfunc->value(cur_pos + dxStep);
    
    if(fa*fb <= 0.0){ 
      // We happen to know that f(pos+dxStep) > 0.0 and f(pos) < 0.0 so we must now compute the precise location
      // where the photon crossed the EB. For that we use a Brent root finder on the interval [pos, pos+dxStep].
      const RealVect xcol = poly::brent_root_finder(a_impfunc, cur_pos, cur_pos + dxStep);
      a_s = (xcol - a_oldPos).vectorLength()/a_pathLen;
      retval = true;

      break;
    }
    else{ // Move to next interval
      cur_pos += dxStep;
    }
  }

  return retval;
}

inline
bool particle_ops::eb_intersection_raycast(const RefCountedPtr<BaseIF>& impfunc,
					   const RealVect&              x0,
					   const RealVect&              x1,
					   const Real&                  tol,
					   Real&                        s){
  
  bool ret = false;

  auto dist = [&](const RealVect& x) -> Real {
    return std::abs(impfunc->value(x));
  };
  
  const Real  D = (x1-x0).vectorLength();  
  const Real D0 = dist(x0);

  if(D > D0){

    const RealVect t = (x1-x0)/D;

    // Move x0 along +t towards the boundary. If we end up too close to the boundary the particle has intersected the BC. 
    RealVect xa = x0;
    Real r  = D;
    Real d  = dist(xa);
    
    while(d < r){
      if(d < tol) {
	s   = (xa-x0).vectorLength()/D;
	ret = true;
	break;
      }
      else{
	xa += t*d;
	r  -= d;
	d   = dist(xa);
      }
    }
  }

  return ret;
}
#include "CD_NamespaceFooter.H"

#endif
