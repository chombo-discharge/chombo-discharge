#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PATCHPARTICLE_H_
#define _PATCHPARTICLE_H_

#include "Box.H"
#include "Vector.H"
#include "ListBox.H"
#include "ListBoxFactory.H"
#include "MeshInterp.H"
#include "Cosmology.H"

#include <string>
using std::string;

#include "NamespaceHeader.H"

/**
   The base class PatchParticle provides an implementation of a particle
   method for advancing particle solution on a grid patch. It consists of
   three major parts: a method to project the particles mass on the grid,
   a method to interpolate the grid forces to the particle positions and
   a method to push the particles.

*/

class PatchParticle
{
public:
  /// null Constructor
  PatchParticle();

  /// full constructor
  PatchParticle(const ProblemDomain&     a_domain,
                const Real&              a_dx,
                const InterpType&        a_poissonInterpType,
                const Cosmology&         a_cosmology)
  {
    define(a_domain, a_dx, a_poissonInterpType, a_cosmology);
  }

  /// Destructor
  virtual ~PatchParticle();

  /// same as full constructor
  virtual void define(const ProblemDomain&     a_domain,
                      const Real&              a_dx,
                      const InterpType&        a_poissonInterpType,
                      const Cosmology&         a_cosmology);

  // set parameters
  virtual void setParameters(const InterpType& a_poissonInterpType);

  /// Set the MeshInterp object
  virtual void setMeshInterp(const InterpType& a_poissonInterpType);

  /// get MeshInterp object
  MeshInterp* getMeshInterp() const;

  /// Factory method - this object is its own factory
  virtual PatchParticle* new_patchParticle() const;

  /// set coefficients for particle pusher
  void setPusherCoeff(const Real& a_time, const Real& a_dt);

  // set m_isPusherCoeffSet=false
  virtual void resetPusherCoeff();

  // return interpolation scheme used by the patcher
  virtual const InterpType& interpOrder() const;

  /// Return true if the object is completely defined.
  virtual bool isDefined() const;

  /// Compute the mass distribution given a particle ListBox
  template <class P>
  void deposit(FArrayBox&        a_rho,
               const ListBox<P>& a_PH,
               const Box&        a_box);

  /// Update/correct the particle variables and return the max
  /// particle speed.
  template <class P, class FTypeBox>
  void updateState(ListBox<P>&     a_PH,
                   Real&           a_maxPartSpeed,
                   const FTypeBox& a_field,
                   const Real&     a_dt,
                   const Box&      a_box);

  /// Implement 2nd order corrections
  template <class P, class FTypeBox>
  void correctState(ListBox<P>&     a_PH,
                    const FTypeBox& a_field,
                    const Real&     a_dt,
                    const Box&      a_box);

  /// Compute the maximum wave speed
  template <class P>
  Real getMaxPartSpeed(const ListBox<P>& a_PH,
                       const Box&        a_box);

protected:
  /// KDK particle pusher: advances particles x,v by dt
  template <class P, class FTypeBox>
  Real KDK_Pusher(ListBox<P>&      a_PH,
                  const FTypeBox&  a_F,
                  const Box&       a_box);

  /// KDK particle pusher: applies 2nd order corrections to particles x,v
  template <class P, class FTypeBox>
  void KDK_Corrector(ListBox<P>&      a_PH,
                     const FTypeBox&  a_F,
                     const Box&       a_box);

  // set coefficients for KDK pusher
  void KDK_PusherCoeff(const Real& a_time, const Real& a_dt);

protected:

  // Has define() been called?
  bool m_isDefined;

  // Problem domain and grid spacing
  ProblemDomain m_domain;
  Real m_dx;
  RealVect m_origin;

  // Handles cosmic expansion
  Cosmology m_cosmology;

  // Interpolation order
  InterpType m_poissonInterpType;

  // was setParameters called?
  bool m_isParametersSet;

  // timestep, pusher time coeff. and have they been set
  Real m_dt;
  Real m_coef_exp_pred, m_coef_vel_pred, m_coef_acc_pred;
  Real m_coef_exp, m_coef_vel, m_coef_acc;
  bool m_isPusherCoeffSet;

  // MeshInterp object: manages particle-grid operations
  MeshInterp* m_meshInterp;
  bool m_isMeshInterpSet;

private:
  // Disallowed for all the usual reasons
  void operator=(const PatchParticle& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  PatchParticle(const PatchParticle& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#include "NamespaceFooter.H"

#include "PatchParticleI.H"

#endif
