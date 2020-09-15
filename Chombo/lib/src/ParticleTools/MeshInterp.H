#ifndef _MESHINTERP_H_
#define _MESHINTERP_H_

#include "Box.H"
#include "FArrayBox.H"
#include "REAL.H"
#include "Vector.H"
#include "RealVect.H"
#include "MeshInterpF_F.H"

#include "NamespaceHeader.H"

typedef enum {NGP, CIC, TSC, W4} InterpType;

/**
   This class is used for handling the interaction between the particles and the mesh
   in PIC calculations. It handles depositing particle masses or charges onto the 
   mesh, and interpolating the result of the force solve back to the particle positions.
 **/
class MeshInterp
{
public:

  MeshInterp();

  /// Full constructor 
  MeshInterp(const Box&  a_domain,
	     const RealVect& a_dx,
	     const RealVect& m_domainLeftEdge);

  /// Define function 
  void define(const Box&  a_domain,
	      const RealVect& a_dx,
	      const RealVect& m_domainLeftEdge);

  /// Deposit the particles in a_particleList onto a_rho.
  /// The type of interpolation to perform is passed in as an argument.
  template <class P>
  void deposit(const List<P>& a_particleList,
	       FArrayBox& a_rho,
	       InterpType& a_interpType);
  
  /// Interpolate the field given by a_field onto the particles in a_particleList.
  /// a_field must have SpaceDim components. The type of interpolation to perform
  /// is passed in as an argument.
  template <class P>
  void interpolate(List<P>& a_particleList,
		   const FArrayBox& a_field,
		   InterpType& a_interpType);
    
private:

  /// a wrapper function that deposits a single particle using the 
  /// interpolating functions defined in MeshInterpF.ChF
  void depositParticle(FArrayBox& a_rho,
		       const RealVect& a_domainLeftEdge,
		       const RealVect&     a_dx,
		       const RealVect& a_position,
		       const Real&     a_strength,
		       const InterpType a_interpType);
  
  /// a wrapper function that interpolates the field to a single
  /// particle position using the interpolating functions defined in
  /// MeshInterpF.ChF
  void interpolateParticle(RealVect& a_particleField,
			   const FArrayBox& a_field,
			   const RealVect& a_domainLeftEdge,
			   const RealVect& a_dx,
			   const RealVect& a_position,
			   const InterpType& a_interpType);
  
  Box m_domain;
  RealVect m_dx;
  RealVect m_domainLeftEdge;
};

#include "NamespaceFooter.H"

// implementation
#include "MeshInterpI.H"

#endif // include guard
