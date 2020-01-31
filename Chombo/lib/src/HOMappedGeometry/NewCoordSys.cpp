#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "NewCoordSys.H"
#include "BoxIterator.H"
#include "IntVect.H"
#include "MBMiscUtil.H"
#include "SimpleDivergence.H"

#include "NamespaceHeader.H"

// given coordinates in mapped space, return locations in real space
void
NewCoordSys::realCoord(FArrayBox& a_x, const FArrayBox& a_Xi,
                       const Box& a_box) const
{
  // default implementation is probably inefficient, but it should work
  BoxIterator bit(a_box);
  for (bit.begin(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      RealVect xi(D_DECL6(a_Xi(iv,0),
                          a_Xi(iv,1),
                          a_Xi(iv,2),
                          a_Xi(iv,3),
                          a_Xi(iv,4),
                          a_Xi(iv,5)));
      RealVect realLoc = realCoord(xi);
      D_TERM6(a_x(iv,0) = realLoc[0];,
              a_x(iv,1) = realLoc[1];,
              a_x(iv,2) = realLoc[2];,
              a_x(iv,3) = realLoc[3];,
              a_x(iv,4) = realLoc[4];,
              a_x(iv,5) = realLoc[5];)
        }
}

/// given coordinate in real space, return its location in the mapped space
void
NewCoordSys::mappedCoord(FArrayBox& a_Xi, const FArrayBox& a_x,
                         const Box& a_box) const
{
  // default implementation is probably inefficient, but it should work
  BoxIterator bit(a_box);
  for (bit.begin(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      RealVect x(D_DECL6(a_x(iv,0),
                         a_x(iv,1),
                         a_x(iv,2),
                         a_x(iv,3),
                         a_x(iv,4),
                         a_x(iv,5)));

      RealVect mappedLoc = mappedCoord(x);
      D_TERM6(a_Xi(iv,0) = mappedLoc[0];,
              a_Xi(iv,1) = mappedLoc[1];,
              a_Xi(iv,2) = mappedLoc[2];,
              a_Xi(iv,3) = mappedLoc[3];,
              a_Xi(iv,4) = mappedLoc[4];,
              a_Xi(iv,5) = mappedLoc[5];)
        }

}

/// given indices, return centers of cells in mapped coordinates
void
NewCoordSys::getCenterMappedCoordinates(FArrayBox& a_Xi,
                                        const Box& a_box) const
{
  CH_TIME("NewCoordSys::getCenterMappedCoordinates");

  CH_assert(a_Xi.box().type() == a_box.type());
  RealVect thisdx = dx();
  RealVect offset = centeringOffset(a_box, thisdx);

  // default implementation is probably inefficient, but it should work
  BoxIterator bit(a_box);
  for (bit.begin(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      RealVect mappedLoc(iv);
      mappedLoc *= thisdx;
      mappedLoc += offset;
      D_TERM6(a_Xi(iv,0) = mappedLoc[0];,
              a_Xi(iv,1) = mappedLoc[1];,
              a_Xi(iv,2) = mappedLoc[2];,
              a_Xi(iv,3) = mappedLoc[3];,
              a_Xi(iv,4) = mappedLoc[4];,
              a_Xi(iv,5) = mappedLoc[5];)
        }
}

/// given indices, return centers of cells in mapped coordinates
RealVect
NewCoordSys::centerMappedCoordinates(const IntVect& a_iv) const
{
  RealVect thisdx = dx();
  RealVect mappedLoc(a_iv);
  mappedLoc += 0.5;
  mappedLoc *= thisdx;
  return mappedLoc;
}

/// note that a_X is in mapped space.
/** fills the destComp component of a_dxdXi
    with the derivative of x w/ respect to Xi
    in the dirX direction
    (default implementation is there in case derived class doesn't
    provide it)
*/
void
NewCoordSys::dXdXi(FArrayBox& a_dxdXi,
                   const FArrayBox& a_xi,
                   int a_destComp,
                   int a_dirX,
                   int a_dirXi,
                   const Box& a_box) const
{
  // default implementation is probably inefficient, but it should work
  BoxIterator bit(a_box);
  for (bit.begin(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      RealVect xi(D_DECL6(a_xi(iv,0),
                          a_xi(iv,1),
                          a_xi(iv,2),
                          a_xi(iv,3),
                          a_xi(iv,4),
                          a_xi(iv,5)));

      Real dxdXi = dXdXi(xi,a_dirX,a_dirXi);
      a_dxdXi(iv,a_destComp) = dxdXi;

    }
}

/// Jacobian evaluated at locations Xi in mapped space
void
NewCoordSys::pointwiseJ(FArrayBox& a_J,
                        const FArrayBox& a_Xi,
                        const Box& a_box) const
{
  BoxIterator bit(a_box);
  for (bit.begin(); bit.ok(); ++bit)
    {
      IntVect iv = bit();

      RealVect xi(D_DECL6(a_Xi(iv,0),
                          a_Xi(iv,1),
                          a_Xi(iv,2),
                          a_Xi(iv,3),
                          a_Xi(iv,4),
                          a_Xi(iv,5)));

      Real J = pointwiseJ(xi);
      a_J(iv,0) = J;
    }
}


/// returns integral of divergence over mapped-grid cells
void
NewCoordSys::computeDivergence(FArrayBox& a_divF,
                               const FluxBox& a_F,
                               const FluxBox& a_N,
                               const Box& a_box,
                               Interval& divInterval) const
{
  // compute N^T*F
  //RealVect thisdx = dx();
  RealVect fakeDx = RealVect::Unit;
  FluxBox NTF(a_box, divInterval.size());
  computeMetricTermProductAverage(NTF, a_F, a_N, a_box, true);
  simpleDivergence(a_divF, NTF, a_box, fakeDx);


}


/// computes 4th-order average of product = N^T*F
/** if a_fourthOrder is false, then only do a second-order dot product
 */
void
NewCoordSys::computeMetricTermProductAverage(FluxBox& a_product,
                                             const FluxBox& a_F,
                                             const FluxBox& a_N,
                                             const Box& a_box,
                                             bool a_fourthOrder) const
{
  computeMetricTermProductAverage(a_product, a_F, a_N, a_F, a_box,
                                  a_fourthOrder);
}


VectorTransformation
NewCoordSys::vectorTransformationMappedToReal(const RealVect& a_Xi) const
{
  VectorTransformation tfm;
  for (int idirMapped = 0; idirMapped < SpaceDim; idirMapped++)
    {
      // The components of vec will be set to
      // the components in each real coordinate direction
      // of the unit vector in the idirMapped mapped coordinate direction
      // at the point a_Xi in mapped space.
      RealVect vec;
      for (int idirReal = 0; idirReal < SpaceDim; idirReal++)
        {
          vec[idirReal] = dXdXi(a_Xi, idirReal, idirMapped);
        }
      // normalize vec
      Real vecLen = vec.vectorLength();
      vec /= vecLen;
      for (int idirReal = 0; idirReal < SpaceDim; idirReal++)
        {
          tfm.setComponent(idirReal, idirMapped, vec[idirReal]);
        }
    }
  return tfm;
}


VectorTransformation
NewCoordSys::vectorTransformationRealToMapped(const RealVect& a_Xi) const
{
  VectorTransformation tfmMapped2Real = vectorTransformationMappedToReal(a_Xi);
  VectorTransformation tfmReal2Mapped = tfmMapped2Real.inverse();
  return tfmReal2Mapped;
}

VectorTransformation
NewCoordSys::vectorTransformationMappedToRealCenter(const IntVect& a_iv) const
{
  RealVect xi = centerMappedCoordinates(a_iv);
  VectorTransformation tfm = vectorTransformationMappedToReal(xi);
  return tfm;
}

VectorTransformation
NewCoordSys::vectorTransformationRealToMappedCenter(const IntVect& a_iv) const
{
  RealVect xi = centerMappedCoordinates(a_iv);
  VectorTransformation tfm = vectorTransformationRealToMapped(xi);
  return tfm;
}

void
NewCoordSys::vectorTransformMappedToRealCenterFab(FArrayBox& a_vectorFab) const
{
  CH_assert(a_vectorFab.nComp() == SpaceDim);
  BoxIterator bit(a_vectorFab.box());
  for (bit.begin(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      VectorTransformation tfm = vectorTransformationMappedToRealCenter(iv);
      RealVect vecMapped(D_DECL6(a_vectorFab(iv, 0),
                                 a_vectorFab(iv, 1),
                                 a_vectorFab(iv, 2),
                                 a_vectorFab(iv, 3),
                                 a_vectorFab(iv, 4),
                                 a_vectorFab(iv, 5)));
      RealVect vecReal = tfm.transform(vecMapped);
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          a_vectorFab(iv, idir) = vecReal[idir];
        }
    }
}

void
NewCoordSys::vectorTransformRealToMappedCenterFab(FArrayBox& a_vectorFab) const
{
  CH_assert(a_vectorFab.nComp() == SpaceDim);
  BoxIterator bit(a_vectorFab.box());
  for (bit.begin(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      VectorTransformation tfm = vectorTransformationRealToMappedCenter(iv);
      RealVect vecReal(D_DECL6(a_vectorFab(iv, 0),
                               a_vectorFab(iv, 1),
                               a_vectorFab(iv, 2),
                               a_vectorFab(iv, 3),
                               a_vectorFab(iv, 4),
                               a_vectorFab(iv, 5)));
      RealVect vecMapped = tfm.transform(vecReal);
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          a_vectorFab(iv, idir) = vecMapped[idir];
        }
    }
}

void
NewCoordSys::vectorTransformRealToMappedAvgFab(FArrayBox& a_vectorFab,
                                               const FArrayBox& a_gradVectorFab) const
{
  // Transform a_vectorFab, containing average of vector over each cell,
  // from basis in real coordinates to basis in mapped coordinates.
  // If T is the transformation and vec is the vector, then:
  // <T(vec)> = <T>(<vec>) + h^2/12 * (grad T) dot (grad vec)
  // So we need:  <T>, grad T, grad <vec>.
  CH_assert(a_vectorFab.nComp() == SpaceDim);
  const Box& bx = a_vectorFab.box();
  Box bx1 = grow(bx, 1);
  BaseFab<VectorTransformation> tfmCenterFab(bx1, 1);
  BoxIterator bit1(bx1);
  for (bit1.begin(); bit1.ok(); ++bit1)
    {
      IntVect iv = bit1();
      tfmCenterFab(iv, 0) = vectorTransformationRealToMappedCenter(iv);
    }
  vectorTransformInternalAvgFab(a_vectorFab,
                                a_gradVectorFab,
                                tfmCenterFab);
}

void
NewCoordSys::vectorTransformMappedToRealAvgFab(FArrayBox& a_vectorFab,
                                               const FArrayBox& a_gradVectorFab) const
{
  // Transform a_vectorFab, containing average of vector over each cell,
  // from basis in mapped coordinates to basis in real coordinates.
  // If T is the transformation and vec is the vector, then:
  // <T(vec)> = <T>(<vec>) + h^2/12 * (grad T) dot (grad vec)
  // So we need:  <T>, grad T, grad <vec>.
  CH_assert(a_vectorFab.nComp() == SpaceDim);
  const Box& bx = a_vectorFab.box();
  Box bx1 = grow(bx, 1);
  BaseFab<VectorTransformation> tfmCenterFab(bx1, 1);
  BoxIterator bit1(bx1);
  for (bit1.begin(); bit1.ok(); ++bit1)
    {
      IntVect iv = bit1();
      tfmCenterFab(iv, 0) = vectorTransformationMappedToRealCenter(iv);
    }
  vectorTransformInternalAvgFab(a_vectorFab,
                                a_gradVectorFab,
                                tfmCenterFab);
}

void
NewCoordSys::vectorTransformInternalAvgFab(FArrayBox& a_vectorFab,
                                           const FArrayBox& a_gradVectorFab,
                                           const BaseFab<VectorTransformation>& a_tfmCenterFab) const
{
  Real frac = 1./24.;
  Real centerScaling = 1. - 2 * SpaceDim * frac;
  Real scaling = 1./12.;
  const Box& bx = a_vectorFab.box();
  BoxIterator bit(bx);
  for (bit.begin(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      // Initialize tfmAvg to cell-centered T, and then we'll apply
      // convolution so that tfm ends up being <T>.
      VectorTransformation tfmAvg = a_tfmCenterFab(iv, 0);
      tfmAvg *= centerScaling;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          tfmAvg += (a_tfmCenterFab(iv - BASISV(idir)) +
                     a_tfmCenterFab(iv + BASISV(idir))) * frac;
        }

      /*
        Now we have the vector transformations we need.
        Set vecNew to transformation of vecOld.
      */

      RealVect vecOld(D_DECL6(a_vectorFab(iv, 0),
                              a_vectorFab(iv, 1),
                              a_vectorFab(iv, 2),
                              a_vectorFab(iv, 3),
                              a_vectorFab(iv, 4),
                              a_vectorFab(iv, 5)));
      RealVect vecNew = tfmAvg.transform(vecOld);
      int indGrad = 0;
      for (int idir = 0; idir < SpaceDim; idir++)
        { // derivative with respect to dx[idir]
          const VectorTransformation& tfmLo = a_tfmCenterFab(iv - BASISV(idir));
          const VectorTransformation& tfmHi = a_tfmCenterFab(iv + BASISV(idir));
          VectorTransformation tfmGradDir = (tfmHi - tfmLo) / 2.;
          // Order of components within gradData, in 2D:
          // dData[0]/dx[0], dData[1]/dx[0], dData[0]/dx[1], dData[1]/dx[1].
          RealVect gradVecDir(D_DECL6(a_gradVectorFab(iv, indGrad),
                                      a_gradVectorFab(iv, indGrad + 1),
                                      a_gradVectorFab(iv, indGrad + 2),
                                      a_gradVectorFab(iv, indGrad + 3),
                                      a_gradVectorFab(iv, indGrad + 4),
                                      a_gradVectorFab(iv, indGrad + 5)));
          // gradVecDirTransformed is component idir of (grad T) dot (grad vec).
          RealVect gradVecDirTransformed = tfmGradDir.transform(gradVecDir);
          vecNew += gradVecDirTransformed * scaling;
          indGrad += SpaceDim;
        }
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          a_vectorFab(iv, idir) = vecNew[idir];
        }
    }
}

#include "NamespaceFooter.H"
