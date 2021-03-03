/*!
  @file   perlin_slab.cpp
  @brief  Implementation of perlin_slab_if.H
  @author Robert Marskar
  @date   March 2021
*/

#include "perlin_slab_if.H"
#include "perlin_plane_if.H"

#include <SmoothUnion.H>
#include <PlaneIF.H>
#include <TransformIF.H>


perlin_slab_if::perlin_slab_if(const RealVect a_ccPoint,
			       const RealVect a_normal,
			       const RealVect a_xyz,
			       const RealVect a_noiseFreq,
			       const int      a_octaves,
			       const Real     a_noiseAmp,
			       const Real     a_persistence,
			       const Real     a_cornerCurv,
			       const bool     a_reseed,
			       const bool     a_fluidInside){

  // By default, we make the slab such that the noisy side normal vector is along +z (+y in 2D). Then we
  // rotate to a_normal


  // Make a slab whose center is at zero and whose widths are given by a_xyz
  Vector<BaseIF*> parts;
  for (int dir = 0; dir < SpaceDim; dir++){
    for (SideIterator sit; sit.ok(); ++sit){
      const int s      = sign(sit());
      const RealVect n = s*BASISREALV(dir);
      const RealVect p = s*0.5*a_xyz[dir]*n;

      if(dir == SpaceDim && sit() == Side::Hi){
	BaseIF* baseif = (BaseIF*) new perlin_plane_if(n, p, false, a_noiseAmp, a_noiseFreq, a_persistence, a_octaves, a_reseed);
	parts.push_back(baseif);
      }
      else{
	BaseIF* baseif = (BaseIF*) new PlaneIF(n, p, false);
	parts.push_back(baseif);
      }
    }
  }

  // Make the smooth union
  BaseIF* bif = (BaseIF*) new SmoothUnion(parts, a_cornerCurv);

  // Rotate and translate into place
  TransformIF* tif = new TransformIF(*bif);
  tif->translate(-BASISREALV(CH_SPACEDIM)*a_xyz[CH_SPACEDIM]); // Move so that "top" point is at RealVect::Zero. This is the rotation point
  tif->rotate(BASISREALV(SpaceDim), a_normal);                 // Rotate so that +z/+y points along a_normal
  tif->translate(a_ccPoint);                                   // Translate to point

  // Done. Delete the rest. 
  m_baseif = RefCountedPtr<BaseIF> (tif);

  for (int i = 0; i < parts.size(); i++){
    delete parts[i];
  }
  delete bif;
}

perlin_slab_if::perlin_slab_if(const perlin_slab_if& a_inputIF){
  m_baseif      = a_inputIF.m_baseif;
  m_fluidInside = a_inputIF.m_fluidInside;
}

perlin_slab_if::~perlin_slab_if(){

}

Real perlin_slab_if::value(const RealVect& a_pos) const {
  Real retval = m_baseif->value(a_pos);

  if(m_fluidInside){
    retval = -retval;
  }

  return retval;
}

BaseIF* perlin_slab_if::newImplicitFunction() const {
  return (BaseIF*) new perlin_slab_if(*this);
}
