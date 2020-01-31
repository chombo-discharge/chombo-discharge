#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "MFGeometry.H"
#include "ComplementIF.H"
#include "NamespaceHeader.H"

///
/**
 */
void MFGeometry::makeLevelSet(      Vector<LevelData<FArrayBox>* >& a_levelSet,
                              const BaseIF&                         a_if,
                              const RealVect&                       a_origin,
                              const RealVect&                       a_coarsestDx,
                              const Vector<int>&                    a_refRatio)
{
  int nlevels = a_levelSet.size();
  Real dx = a_coarsestDx[0];
  for (int ilev=0; ilev<nlevels; ilev++)
    {
      LevelData<FArrayBox>& data = *(a_levelSet[ilev]);
      for (DataIterator dit = data.dataIterator(); dit.ok(); ++dit)
        {
          FArrayBox& fab = data[dit];
          Box fbox(fab.box());
          for (BoxIterator bit(fbox); bit.ok(); ++bit)
            {
              RealVect loc(bit());
              loc += 0.5;
              loc *= dx;
              loc += a_origin;
              Real val = a_if.value(loc);
              fab.set(bit(), 0, val);
            }
        }
      dx /= a_refRatio[ilev];
    }
}

///
/**
 */
void MFGeometry::defineGeometry(      RefCountedPtr<MFIndexSpace>&  a_mfis,
                                const BaseIF&                       a_if,
                                const ProblemDomain&                a_domain,
                                const RealVect&                     a_dx,
                                const RealVect&                     a_origin,
                                      int                           a_maxCells)
{
  int verbosity = 0;
  Vector<GeometryService*> geometry(2, NULL);

  ComplementIF inside(a_if, false);
  GeometryShop workshop0(inside, verbosity, a_dx);
  workshop0.m_phase = 0;
  geometry[0] = &workshop0;

  ComplementIF outside(a_if, true);
  GeometryShop workshop1(outside, verbosity, a_dx);
  workshop1.m_phase=1;
  geometry[1] = &workshop1;

  a_mfis->define(a_domain.domainBox(), a_origin, a_dx[0], geometry, a_maxCells);
}

///
/**
 */
RefCountedPtr<PlaneIF>
MFGeometry::definePlaneGeometry(      RefCountedPtr<MFIndexSpace>&  a_mfis,
                                const ProblemDomain&                a_domain,
                                const RealVect&                     a_dx,
                                const RealVect&                     a_origin,
                                const RealVect&                     a_normal,
                                const RealVect&                     a_point,
                                      bool                          a_inside,
                                      int                           a_maxCells)
{
  RefCountedPtr<PlaneIF> ifunc = RefCountedPtr<PlaneIF>(new PlaneIF(a_normal,
                                                                    a_point,
                                                                    a_inside));
  defineGeometry(a_mfis,
                 *ifunc,
                 a_domain,
                 a_dx,
                 a_origin,
                 a_maxCells);
  return ifunc;
}

///
/**
 */
RefCountedPtr<SphereIF>
MFGeometry::defineSphereGeometry(      RefCountedPtr<MFIndexSpace>&  a_mfis,
                                 const ProblemDomain&                a_domain,
                                 const RealVect&                     a_dx,
                                 const RealVect&                     a_origin,
                                 const RealVect&                     a_center,
                                 const Real                          a_radius,
                                       bool                          a_inside,
                                       int                           a_maxCells)
{
  RefCountedPtr<SphereIF> ifunc = RefCountedPtr<SphereIF>(new SphereIF(a_radius,
                                                                       a_center,
                                                                       a_inside));
  defineGeometry(a_mfis,
                 *ifunc,
                 a_domain,
                 a_dx,
                 a_origin,
                 a_maxCells);
  return ifunc;
}

///
/**
 */
RefCountedPtr<EllipsoidIF>
MFGeometry::defineEllipsoidGeometry(      RefCountedPtr<MFIndexSpace>&  a_mfis,
                                    const ProblemDomain&                a_domain,
                                    const RealVect&                     a_dx,
                                    const RealVect&                     a_origin,
                                    const RealVect&                     a_center,
                                    const RealVect&                     a_radii,
                                          bool                          a_inside,
                                          int                           a_maxCells)
{
  RefCountedPtr<EllipsoidIF> ifunc = RefCountedPtr<EllipsoidIF>(new EllipsoidIF(a_radii,
                                                                                a_center,
                                                                                a_inside));
  defineGeometry(a_mfis,
                 *ifunc,
                 a_domain,
                 a_dx,
                 a_origin,
                 a_maxCells);
  return ifunc;
}

///
/**
 */
RefCountedPtr<TiltedCylinderIF>
MFGeometry::defineTiltedCylinderGeometry(      RefCountedPtr<MFIndexSpace>&  a_mfis,
                                         const ProblemDomain&                a_domain,
                                         const RealVect&                     a_dx,
                                         const RealVect&                     a_origin,
                                         const Real                          a_radius,
                                         const RealVect&                     a_axis,
                                         const RealVect&                     a_center,
                                               bool                          a_inside,
                                               int                           a_maxCells)
{
  RefCountedPtr<TiltedCylinderIF> ifunc = RefCountedPtr<TiltedCylinderIF>
    (new TiltedCylinderIF(a_radius,
                          a_axis,
                          a_center,
                          a_inside));
  defineGeometry(a_mfis,
                 *ifunc,
                 a_domain,
                 a_dx,
                 a_origin,
                 a_maxCells);
  return ifunc;
}

#include "NamespaceFooter.H"
