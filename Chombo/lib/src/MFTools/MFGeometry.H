#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MFGEOMETRY_H_
#define _MFGEOMETRY_H_

#include "MFIndexSpace.H"
#include "PolyGeom.H"
#include "GeometryShop.H"
#include "PlaneIF.H"
#include "SphereIF.H"
#include "EllipsoidIF.H"
#include "TiltedCylinderIF.H"

#include "NamespaceHeader.H"
///
/**
   class to encapsulate the common geometric operations
   for MF
 */
class MFGeometry
{
public:
  ///
  /**
   */
  static
  void makeLevelSet(      Vector<LevelData<FArrayBox>* >& a_levelSet,
                    const BaseIF&                         a_if,
                    const RealVect&                       a_origin,
                    const RealVect&                       a_coarsestDx,
                    const Vector<int>&                    a_refRatio);

  ///
  /**
   */
  static
  void defineGeometry(      RefCountedPtr<MFIndexSpace>&  a_mfis,
                      const BaseIF&                       a_if,
                      const ProblemDomain&                a_domain,
                      const RealVect&                     a_dx,
                      const RealVect&                     a_origin,
                            int                           a_maxCells = 1024);

  ///
  /**
   */
  static
  RefCountedPtr<PlaneIF>
  definePlaneGeometry(      RefCountedPtr<MFIndexSpace>&  a_mfis,
                      const ProblemDomain&                a_domain,
                      const RealVect&                     a_dx,
                      const RealVect&                     a_origin,
                      const RealVect&                     a_normal,
                      const RealVect&                     a_point,
                            bool                          a_inside = true,
                            int                           a_maxCells = 1024);

  ///
  /**
   */
  static
  RefCountedPtr<SphereIF>
  defineSphereGeometry(      RefCountedPtr<MFIndexSpace>&  a_mfis,
                       const ProblemDomain&                a_domain,
                       const RealVect&                     a_dx,
                       const RealVect&                     a_origin,
                       const RealVect&                     a_center,
                       const Real                          a_radius,
                             bool                          a_inside = true,
                             int                           a_maxCells = 1024);

  ///
  /**
   */
  static
  RefCountedPtr<EllipsoidIF>
  defineEllipsoidGeometry(      RefCountedPtr<MFIndexSpace>&  a_mfis,
                          const ProblemDomain&                a_domain,
                          const RealVect&                     a_dx,
                          const RealVect&                     a_origin,
                          const RealVect&                     a_center,
                          const RealVect&                     a_radii,
                                bool                          a_inside = true,
                          const int                           a_maxCells = 1024);

  ///
  /**
   */
  static
  RefCountedPtr<TiltedCylinderIF>
  defineTiltedCylinderGeometry(      RefCountedPtr<MFIndexSpace>&  a_mfis,
                               const ProblemDomain&                a_domain,
                               const RealVect&                     a_dx,
                               const RealVect&                     a_origin,
                               const Real                          a_radius,
                               const RealVect&                     a_axis,
                               const RealVect&                     a_center,
                                     bool                          a_inside = true,
                                     int                           a_maxCells = 1024);
};
#include "NamespaceFooter.H"
#endif
