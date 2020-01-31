#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _ADVECTCOSINEBELLIBC_H_
#define _ADVECTCOSINEBELLIBC_H_

#include <iostream>

#include "LevelData.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "RealVect.H"
//#include "AMRIO.H"
#include "PhysAdvectMappedIBC.H"

#include "NamespaceHeader.H"

/// Example of a cosine bell distribution
/**
   Parameters:
     a_ambientHeight  - Ambient height
     a_deltaHeight    - Height of the cosine bell
     a_center         - Center of the cosine bell (vector)
     a_size           - Standard deviation of the gaussian
     a_advectVelocity - Advection velocity
     a_advectAngle    - Advection angle
         a_evalTime       - Time to evaluate initial profile
 */
class AdvectCosineBellIBC: public PhysAdvectMappedIBC
{
public:

 // type of velocity field
  enum velType
    {
      UNIFORM,
      SOLIDBODY
    };

  /// Null constructor
  /**
   */
  AdvectCosineBellIBC();

  /// Constructor which defines parameters
  /**
   */
  AdvectCosineBellIBC(
                    const Real&     a_ambientDensity,
                    const Real&     a_deltaDensity,
                    const RealVect& a_center,
                    const Real&     a_size,
                    const Real&     a_advectVelocity,
                    const Real&     a_advectAngle,
                    const Real&     a_evalTime);

  /// Destructor
  virtual ~AdvectCosineBellIBC();

  /// set uniform velocity field
  void setUniformVel();

  /// set solid-body rotation velocity field (this is default)
  void setSolidBodyRotation();

  PhysMappedIBC *new_physIBC();

  /// set uniform velocity field m_velUniform to a_velUniform
  void setVelUniform(const RealVect& a_velUniform);

  /// Set up initial conditions
  /**
   */
  void initialize(LevelData<FArrayBox>& a_U);

 ///
  /* Return advection velocity.
   */
  virtual void advVel(FArrayBox& a_advVelFab,
                      Real a_time = 0.);

  virtual void primBC(FArrayBox&            a_WGdnv,
                      const FArrayBox&      a_Wextrap,
                      const FArrayBox&      a_W,
                      const int&            a_dir,
                      const Side::LoHiSide& a_side,
                      const Real&           a_time);

protected:

  virtual Real cosineBell(Real a_dist);

  //  bool m_params_are_set;

  // either UNIFORM or SOLIDBODY
  int m_velType;

  /// Ambient density of the advected field
  Real m_ambientDensity;

  /// Height of the cosine bell
  Real m_deltaDensity;

  /// Center of the cosine bell as (latitude, longitude)
  RealVect m_center;

  /// Initial time to evaluate profile
  Real m_evalTime;

  /// Size of the cosine bell
  Real m_size;

  /// Advection velocity
  Real m_advectVelocity;

  /// Advection angle
  Real m_advectAngle;

  /// Uniform advection velocity, if any; it will be MULTIPLIED by m_advectVelocity
  RealVect m_velUniform;

  /// whether to use m_velUniform
  bool m_useVelUniform;

  //  LevelData<FArrayBox> m_advVel;
  //
  //  /// Is m_advVel filled in?
  //  bool m_haveAdvVel;
};

#include "NamespaceFooter.H"

#endif
