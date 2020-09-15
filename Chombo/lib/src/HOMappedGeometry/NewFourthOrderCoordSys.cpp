#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "NewFourthOrderCoordSys.H"
#include "BoxIterator.H"
#include "LoHiCenter.H"
#include "FourthOrderUtil.H"
#include "MBMiscUtil.H"

#include "FourthOrderCoordSysF_F.H"
#include "PointwiseDotProdF_F.H"
#include "FourthOrderUtilF_F.H"
#include "FCDivergenceF_F.H"
#include "FaceGradientF_F.H"
#include "SimpleDivergence.H"

#include "GaussianQuadrature.H"
#include "NewtonCotesQuadrature.H"

#include "DebugOut.H"
#include "NamespaceHeader.H"

inline Real computePointwiseJ( const Vector<Real>& Jmtrx );

  /// default constructor
NewFourthOrderCoordSys::NewFourthOrderCoordSys() : m_dx(-1*IntVect::Unit)
{
  m_quadraturePtr = NULL;
  // set default quadrature to be 4-pt Newton-Cotes
  int numPts = 4;
  //m_quadraturePtr = new GaussianQuadrature(numPts);
  m_quadraturePtr = new NewtonCotesQuadrature(numPts);

  // default volume interval is all components
  m_volInterval.define(0, SpaceDim-1);
}

NewFourthOrderCoordSys::~NewFourthOrderCoordSys()
{
  if (m_quadraturePtr != NULL)
    {
      delete m_quadraturePtr;
    }
}

/// given coordinates in mapped space, return locations in real space
void
NewFourthOrderCoordSys::realCoord(FArrayBox& a_x, const FArrayBox& a_Xi,
                                  const Box& a_box) const
{
  NewCoordSys::realCoord(a_x, a_Xi, a_box);
}

/// given coordinate in real space, return its location in the mapped space
void
NewFourthOrderCoordSys::mappedCoord(FArrayBox& a_Xi, const FArrayBox& a_x,
                                    const Box& a_box) const
{
  NewCoordSys::mappedCoord(a_Xi, a_x, a_box);
}

/// return Cartesian XYZ locations of nodes
/** nodeCoords should have dimension returned by dimension()
 */
void
NewFourthOrderCoordSys::getNodeRealCoordinates(FArrayBox& a_nodeCoords,
                                               const Box& a_box) const
{
  BoxIterator bit(a_box);
  for (bit.begin(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      // node centering
      RealVect X = m_dx * iv;

      // default implementation just calls RealCoord
      RealVect nodeCoord = realCoord(X);
      D_TERM6(a_nodeCoords(iv,0) = nodeCoord[0];,
              a_nodeCoords(iv,1) = nodeCoord[1];,
              a_nodeCoords(iv,2) = nodeCoord[2];,
              a_nodeCoords(iv,3) = nodeCoord[3];,
              a_nodeCoords(iv,4) = nodeCoord[4];,
              a_nodeCoords(iv,5) = nodeCoord[5];)
        } // end loop over nodes
}


/// note that a_Xi is in mapped space.
/** fills the destComp component of a_dxdXi
    with the derivative of x w/ respect to Xi
    in the dirX direction
    (default implementation is there in case derived class doesn't
    provide it)
*/
void
NewFourthOrderCoordSys::dXdXi(FArrayBox& a_dxdXi,
                              const FArrayBox& a_Xi,
                              int a_destComp,
                              int a_dirX,
                              int a_dirXi,
                              const Box& a_box) const
{
  NewCoordSys::dXdXi(a_dxdXi, a_Xi, a_destComp, a_dirX,
                     a_dirXi, a_box);

}

/*--------------------------------------------------------------------*/
///  Compute the average volume flux on a face
/**  \param[out] a_volFlux
 *                      The volume flux on each face
 *   \param[in]  a_Nt   Metric terms on the faces (only need terms
 *                      normal to each face.  Must have 1 extra ghost
 *                      cells over 'a_box'.
 *   \param[in]  a_box  Cell-centered box.  The volume fluxes are
 *                      required on each face.
 *   \note
 *   <ul>
 *     <li> After testing, cellVol should use this routine.
 *   </ul>
 *//*-----------------------------------------------------------------*/

void
NewFourthOrderCoordSys::volFlux(FluxBox& a_volFlux, const FluxBox& a_Nt,
                                const Box& a_box) const
{
  // Here the physical flux is x/D given by the mapping
  // \tilde(X(\xi)) = X(\xi)/D.  Note that X needs to have two more ghost cells
  // cells than a_volFlux and a_Nt needs to have one more ghost cell:
  // +1 one for tangential gradients (of X and a_Nt) required for the product
  // and then +1 again to get the average of X.
  CH_assert(a_Nt.box().contains(grow(a_box, 1)));
  FluxBox X(grow(a_box, 2), SpaceDim);

  // Compute the point values of X
  X.setVal(0.0);
  for (int dir=0; dir<SpaceDim; dir++)
    {
      RealVect offset = 0.5*m_dx;
      offset[dir] = 0.0;

      FArrayBox& thisXdir = X[dir];
      BoxIterator bit(thisXdir.box());
      // this is going to be slow, but we can
      // eventually move this into fortran
      for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect iv = bit();
          RealVect mappedLoc = m_dx*iv;
          mappedLoc += offset;
          RealVect realLoc = realCoord(mappedLoc);

          for (int comp=m_volInterval.begin(); comp<=m_volInterval.end();
               comp++)
            {
              // only use components in volInterval
              thisXdir(iv,comp) = realLoc[comp];
            }
        }
      thisXdir /= m_volInterval.size();
    } // end loop over directions

  // Convert point values to 4th-order face averages
  fourthOrderAverageFace(X);

  // Compute the volume flux on the faces.
  computeMetricTermProductAverage(a_volFlux, X, a_Nt, SpaceDim, X, a_box);
}

/// computes cell volumes
void
NewFourthOrderCoordSys::cellVol(FArrayBox& a_vol, const FluxBox& a_N,
                                const Box& a_box) const
{

  // note that F needs to have two more ghost cells than cellVolumes
  // because we need to take derivatives of <F>, and because we
  // need to take derivatives in order to compute 4th-order avg of F
  Box fbox(a_box);
  fbox.grow(2);
  FluxBox F(fbox, SpaceDim);


  F.setVal(0.0);
  for (int dir=0; dir<SpaceDim; dir++)
    {
      RealVect offset = 0.5*m_dx;
      offset[dir] = 0.0;

      FArrayBox& thisFdir = F[dir];
      BoxIterator bit(thisFdir.box());
      // this is going to be slow, but we can
      // eventually move this into fortran
      for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect iv = bit();
          RealVect mappedLoc = m_dx*iv;
          mappedLoc += offset;
          RealVect realLoc = realCoord(mappedLoc);

          for (int comp=m_volInterval.begin(); comp<=m_volInterval.end(); comp++)
            {
              // only use components in volInterval
              thisFdir(iv,comp) = realLoc[comp];
            }
#if 0
          D_TERM6(thisFdir(iv,0) = realLoc[0];,
                  thisFdir(iv,1) = realLoc[1];,
                  thisFdir(iv,2) = realLoc[2];,
                  thisFdir(iv,3) = realLoc[3];,
                  thisFdir(iv,4) = realLoc[4];,
                  thisFdir(iv,5) = realLoc[5]);
#endif
        }
    } // end loop over directions


  // convert point values to 4th-order face averages
  fourthOrderAverageFace(F);

  Interval divInterval = a_vol.interval();
  computeDivergence(a_vol, F, a_N, a_box, divInterval);

  // now divide by number of components
  //       a_cellVolumes[dit] /= SpaceDim;
  a_vol /= m_volInterval.size();

}

/// computes integral of N over each face of a_box
void
NewFourthOrderCoordSys::getN(FluxBox& a_N, const Box& a_box) const
{
  a_N.setVal(0.0);

  for (int faceDir=0; faceDir<SpaceDim; faceDir++)
    {
      FArrayBox& thisDirData = a_N[faceDir];

      // loop over components on this face
      for (int d=0; d<SpaceDim; d++)
        {
          // normal component -- use Stoke's theorem, integrate
          // on edges
          if (d == faceDir)
            {
              // loop over high-low for edges
              SideIterator sit;
              for (sit.begin(); sit.ok(); ++sit)
                {
                  // loop over alternate directions
                  for (int edgeDir = 0; edgeDir<SpaceDim; edgeDir++)
                    {
                      if (edgeDir != faceDir)
                        {
                          incrementFaceMetricWithEdgeTerm(thisDirData,
                                                          faceDir, edgeDir,
                                                          thisDirData.box(),
                                                          sit());
                        }
                    } // end loop over edge directions
                } // end loop over high-low
            }
          else
            {
              // transverse direction...
              computeTransverseFaceMetric(thisDirData,
                                          thisDirData.box(),
                                          faceDir,
                                          d);

            }  // end if transverse direction
        } // end loop over components on this face
    } // end loop over face directions
}

/// Compute cell-average J from the volume flux
void
NewFourthOrderCoordSys::getAvgJ(FArrayBox&     a_avgJ,
                                const FluxBox& a_volFlux,
                                const Box&     a_box) const
{
  simpleDivergence(a_avgJ, a_volFlux, a_box, m_dx);
}

/// computes cell-averaged J
void
NewFourthOrderCoordSys::getAvgJ(FArrayBox& a_avgJ, const Box& a_box) const
{

  // avgJ is simply cell volume divided by volume in mapped space
  // this may be slow. Will probably want to have a version which
  // passes in a precomputed N matrix

  // N needs an extra ghost in order to do fourth-order multiplication
  int nGhostN = 1;
  Box Nbox(a_box);
  Nbox.grow(nGhostN);
  FluxBox N(Nbox, getNumN());
  getN(N, Nbox);

  // now compute cell volumes
  cellVol(a_avgJ, N, a_box);

  Real mappedVol = D_TERM6(m_dx[0],
                           *m_dx[1],
                           *m_dx[2],
                           *m_dx[3],
                           *m_dx[4],
                           *m_dx[5]);
  a_avgJ /= mappedVol;

}

/// computes face-averaged 1/J
void
NewFourthOrderCoordSys::getAvgJinverse(FluxBox& a_avgJinverse,
                                       const Box& a_box) const
{
  // first compute cell-centered J
  // need 2 extra ghost cells around a_box
  Box Jbox(a_box);
  Jbox.grow(2);
  FArrayBox cellAvgJ(Jbox, 1);
  getAvgJ(cellAvgJ, Jbox);
  getAvgJinverse(a_avgJinverse, cellAvgJ, a_box);
}

/// computes face-averaged 1/J
void
NewFourthOrderCoordSys::getAvgJinverse(FluxBox& a_avgJinverse,
                                       const FArrayBox& a_avgJ,
                                       const Box& a_box) const
{
  // need a data holder with an extra ghost cell in order to
  // be able to compute tangential gradients to compute
  Box tempBox(a_box);
  tempBox.grow(1);
  FluxBox tempJinverse(tempBox, 1);
  fourthOrderCellToFace(tempJinverse, a_avgJ);

  // now compute tanGrad(J)
  int numTanGradComps = SpaceDim-1;
  // probably actually want to reduce this by one
  Box tanGradBox(tempBox);
  tanGradBox.grow(-1);
  FluxBox tanGradJ(tanGradBox, numTanGradComps);

  computeTangentialGrad(tanGradJ, tempJinverse, tanGradBox);

  // now compute <J^(-1)> = 1/<J>[<1> + h^2/24(grad <J>^2)/<J>^2]
  FluxBox& thisJinverse = a_avgJinverse;
  FluxBox& thisTempJinverse = tempJinverse;
  FluxBox& thisTanGradJ = tanGradJ;

  //thisTanGradJ.setVal(0.0);
  for (int dir=0; dir<SpaceDim; dir++)
    {
      // compute box over which we can compute this
      // note that we don't need dx because the h^2 cancels with the
      // denominator in the (undivided) gradient
      thisJinverse[dir].copy(thisTempJinverse[dir], thisJinverse[dir].box());

      Box thisFaceBox(thisJinverse[dir].box());
      thisFaceBox &= thisTanGradJ[dir].box();
      FORT_COMPUTEJINVERSE(CHF_FRA1(thisJinverse[dir],0),
                           CHF_FRA(thisTanGradJ[dir]),
                           CHF_BOX(thisFaceBox));
    }

}


/// given N and Jinverse, computes N/J
/** provided for convenience, since the multicomponent N/J product
    requires a bit of processing to work with the existing
    fourth-order multiplication function
*/
void
NewFourthOrderCoordSys::computeNJinverse(FluxBox& a_NJinverse,
                                         const FluxBox& a_Jinverse,
                                         const FluxBox& a_N,
                                         const Box& a_cellBox) const
{
  // essentially a multicomponent wrapper around metricTermProductAverage

  // the absolute simplest way to do this (minimizing code duplication) is
  // to copy 1/J to a multi-component temp data holder, and then call
  // computeProductAverage
  FluxBox multiCompJinverse(a_Jinverse.box(),
                            a_NJinverse.nComp());

  // single component holder for NJinverse
  FluxBox tempStorage(a_cellBox, 1);

  // need to do this one component at a time
  for (int sComp=0; sComp<SpaceDim; sComp++)
    {
      for (int dComp=0; dComp<SpaceDim; dComp++)
        {
          // now do a copy of a_JInverse into each
          // component of the temporary

          // want to set only the dir-th component to be Jinverse, and
          // the rest to be zero, in order to isolate this component
          multiCompJinverse.setVal(0.0);

          // now copy Jinverse into the comp-th component
          // loop over faceDirs
          for (int dir=0; dir<SpaceDim; dir++)
            {
              const FArrayBox& thisJinverseDir = a_Jinverse[dir];
              FArrayBox& thisMultiCompDir = multiCompJinverse[dir];

              // copy from 0th comp of Jinverse->comp-th
              // component of multiCompJinverse
              int copyComp = getNcomponent(sComp,dComp);
              thisMultiCompDir.copy(thisJinverseDir, 0, copyComp, 1);
            } // end loop over face directions

          // can now do this by using the computeProductAverage function
          computeMetricTermProductAverage(tempStorage, multiCompJinverse,
                                          a_N, a_cellBox, true);

          // now copy tempStorage into the relevent component in a_NJinverse
          //FluxBox& thisTemp = tempStorage[dit];
          //FluxBox& thisNJinverse = a_NJinverse[dit];
          int comp = getNcomponent(sComp,dComp);
          for (int dir=0; dir<SpaceDim; dir++)
            {
              a_NJinverse[dir].copy(tempStorage[dir], 0, comp, 1);
            }
        } // end loop over d
    } // end loop over components

}



/// Jacobian evaluated at locations Xi in mapped space
void
NewFourthOrderCoordSys::pointwiseJ(FArrayBox& a_J,
                                   const FArrayBox& a_Xi,
                                   const Box& a_box) const
{
  BoxIterator bit(a_box);
  for (bit.begin(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      RealVect mappedLoc(D_DECL6(a_Xi(iv,0),
                                 a_Xi(iv,1),
                                 a_Xi(iv,2),
                                 a_Xi(iv,3),
                                 a_Xi(iv,4),
                                 a_Xi(iv,5)));
      a_J(iv,0) = pointwiseJ(mappedLoc);
    }

}

//**FIXME This is independent of coordSys and should probably be in
//**FourthOrderUtil instead -- except for getNcomponent which
//**seemingly should  be constant anyways.

//  Computes 4th-order average of product = N^T*F.  F can be a flux dyad
/** This general implementation supports computing an average of a product where
    F is a flux dyad.

    Typically N only contains the metrics components normal to the face and F is
    a flux; for each variable in F, a scalar product is returned on the face.

    In an alternate configuration, all components in N are provided and F
    has a_NNumComp terms for each variable.  But again, for each variable in
    F, a scalar product is returned on the face.

    If F is a flux dyad, it can be stored as (var, space) or (space, var) where
    the number of space dimensions is taken to be a_NNumComp (to support the
    alternate configuration described above).  If a_fluxSpaceContiguous == 1,
    the former storage is assumed (var, space) where indexing in space is
    contiguous.  Otherwise the stride in space is determined from the number
    of variables in F.  If F is just a vector (1 variable) then
    a_fluxSpaceContiguous should be set to 1

    Intervals can be used to specify where the results are placed in a_product
    and which components of the flux to use.  The interval for the flux is
    only in terms of variables indexes, not components in the flux dyad.
    An assertion confirms that
    (a_varIntervalProduct.size() == a_varIntervalF.size())

    If (a_FforGrad < grow(a_box, 1)), one-sided differences may be used to
    compute the gradients.  However, this is only allowed at the domain boundary
    and a_problemDomainPtr must be specified

    \param[out] a_product
                        The average of the product
    \param[in]  a_F     Flux component
    \param[in]  a_N     Metrics component.  Stored on getNcomponent(0, dir)
                        if a_NNumComp == SpaceDim and everywhere otherwise.
    \param[in]  a_NNumComp
                        Specifies the configuration.  If == SpaceDim, then only
                        the normal components of a_N are used and F is a normal
                        flux dyad.  If == SpaceDim^2, the alternate
                        configuration is specified
    \param[in]  a_FforGrad
                        Used for computing gradients of F
    \param[in]  a_box   The box on whose faces the product is computed.  Note
                        the following requirements:
                          a_product  >= a_box
                        2nd-order:
                          a_F        >= a_box
                          a_N        >= a_box
                        4th-order:
                          a_F        >= a_box
                          a_N        >= grow(a_box, 1)
                          a_FforGrad >= grow(a_box, 1) unless one-sided used
    \param[in]  a_fourthOrder (default T)
                        T - 4th order average
                        F - 2nd order average
    \param[in]  a_varIntervalProduct
                        The components of a_product in which to store the
                        results for each variable.
                        (default has size == 0 and the interval is specified as
                        0:a_product.nComp() - 1
    \param[in]  a_varIntervalF
                        The interval of variables to be used in F (not the
                        interval of components if F is a dyad).  Think of it
                        as begin() giving a start location for the first
                        variable and end() providing the number of variables.
                        **FIXME This shouldn't be obtained by an interval.
                        (default has size == 0 and the interval is specified as
                        0:a_F.nComp()/a_NNumComp - 1
    \param[in]  a_fluxSpaceContiguous
                        1 - Indexing in space is contiguous in the flux dyad
                        (default 1)
    \param[in]  a_problemDomainPtr
                        Pointer to the problem domain, only required if one-
                        sided difference are to be used to compute grad F.
*/
void
NewFourthOrderCoordSys::computeMetricTermProductAverage(
  FluxBox&                   a_product,
  const FluxBox&             a_F,
  const FluxBox&             a_N,
  const int                  a_NNumComp,
  const FluxBox&             a_FforGrad,
  const Box&                 a_box,
  bool                       a_fourthOrder,
  Interval                   a_varIntervalProduct,
  Interval                   a_varIntervalF,
  int                        a_fluxSpaceContiguous,
  const ProblemDomain *const a_problemDomainPtr) const
{
  CH_TIME("NewFourthOrderCoordSys::computeMetricTermProductAverage");

  // Confirm that we have sufficient data
  CH_assert(a_product.box().contains(a_box));
  CH_assert(a_F.box().contains(a_box));
  bool useCenterGradF = true;  // If can used a centered difference everywhere
  if (a_fourthOrder)
    {
      // We need an extra layer for computing gradients
      Box grownBox = grow(a_box, 1);
      CH_assert(a_N.box().contains(grownBox));
      CH_assert(a_FforGrad.box().contains(a_box));
      if (a_FforGrad.box().contains(grownBox))
        {
          useCenterGradF = true;
        }
      else
        {
          useCenterGradF = false;
          CH_assert(a_problemDomainPtr != NULL);
          if (a_problemDomainPtr->domainBox().contains(grownBox))
            {
              MayDay::Error("From NewFourthOrderCoordSys::computeMetricTerm"
                            "ProductAverage:  You are only supposed to use "
                            "one-sided differences of a_FforGrad at the "
                            "boundary but a_box is away from the boundary");
              // This means that the size of a_FforGrad relative to a_box is
              // used to determine when one-sided differences are required.  But
              // we only expect this to be required at domain boundaries.  What
              // this means is that you should have made a_FforGrad bigger here.
            }
        }
    }
  else
    {
      CH_assert(a_N.box().contains(a_box));
    }

  // If the size of the intervals are < 0, (probably by default), set according
  // to the fabs
  if (a_varIntervalProduct.size() == 0)
    {
      a_varIntervalProduct.define(0, a_product.nComp()-1);
    }
  if (a_varIntervalF.size() == 0)
    {
      a_varIntervalF.define(0, a_F.nComp()/a_NNumComp - 1);
    }
  // Number of solution variables
  const int numVar = a_varIntervalProduct.size();
  // Start index in product fab
  const int varProductStart =  a_varIntervalProduct.begin();
  // Start index in F fab (in terms of solution variable index)
  const int varFStart = a_varIntervalF.begin();
  // Sanity checks
  CH_assert(numVar == a_varIntervalF.size());
  CH_assert(a_F.nComp() >= a_NNumComp*numVar);

  // This is the stride through the space components of the flux dyad
  int FSpaceStride = 1;
  if (a_fluxSpaceContiguous != 1)
    {
      // The variables are contiguous.  What is the step to advance to the
      // next space component of the same variable?  This tries to determine
      // the total number of variables in a_F which may be different from
      // a_product and different from numVar.
      FSpaceStride = a_F.nComp()/a_NNumComp;
    }

  // Another check that we have enough components in a_F,
  // in which the max index is a_F.nComp()-1.
  // Within FORT_POINTFDOTN:
  // iVar - varProdStart has range varFStart + 0:(varNum-1).
  int iFluxInitMax = varFStart + numVar-1;
  if (FSpaceStride == 1)
    {
      iFluxInitMax *= SpaceDim;
    }
  int iFluxMax = iFluxInitMax + (a_NNumComp - 1) * FSpaceStride;
  CH_assert(a_F.nComp()-1 >= iFluxMax);

  // Check that we have enough components in a_product,
  // in which the max index is a_product.nComp()-1.
  int iVarMax = varProductStart + numVar - 1;
  CH_assert(a_product.nComp()-1 >= iVarMax);

//--Loop through the face directions

  for (int dir=0; dir<SpaceDim; dir++)
    {
      FArrayBox& thisProductDir = a_product[dir];
      const FArrayBox& thisFDir = a_F[dir];
      const FArrayBox& thisNDir = a_N[dir];

      Box faceBox(a_box);
      faceBox.surroundingNodes(dir);

//--First part is easy -- dot product of F and N

      int NStartComp = 0;
      // N is always size SpaceDim*SpaceDim.  If we are only using the normal
      // components (size SpaceDim), find the start index.
      if (a_NNumComp == SpaceDim)
        {
          NStartComp = getNcomponent(0, dir);
        }

      // Check that we have enough components in thisNDir,
      // in which the max index is thisNDir.nComp()-1.
      int iNMax = NStartComp + a_NNumComp - 1;
      CH_assert(thisNDir.nComp()-1 >= iNMax);

      FORT_POINTFDOTN(CHF_FRA(thisProductDir),
                      CHF_CONST_INT(varProductStart),
                      CHF_CONST_INT(numVar),
                      CHF_CONST_FRA(thisFDir),
                      CHF_CONST_INT(varFStart),
                      CHF_CONST_INT(FSpaceStride),
                      CHF_CONST_FRA(thisNDir),
                      CHF_CONST_INT(NStartComp),
                      CHF_CONST_INT(a_NNumComp),
                      CHF_BOX(faceBox));

//--If fourth order, add in the products of the gradients

      if (a_fourthOrder)
        {
          const FArrayBox& thisFforGradDir = a_FforGrad[dir];

          // Temp storage for tangential gradients.  Stored in C as
          // (comp, tanDir) where comp is over all components.

          // For the dyad F, the storage is specifically
          // (0:numVar-1, 0:a_NNumComp-1, 0:nTanDirs-1)
          FArrayBox tanGradF(faceBox, numVar*a_NNumComp*(SpaceDim-1));
          // For N, this is (0:a_NNumComp-1, 0:nTanDirs-1)
          FArrayBox tanGradN(faceBox, a_NNumComp*(SpaceDim-1));

          int tanDirs[SpaceDim] =
          {
            D_DECL6(0, 1, 2, 3, 4, 5)
          };
          const int nTanDirs = SpaceDim - 1;
          for (int i = dir; i != nTanDirs; ++i)
            {
              ++tanDirs[i];
            }

          // Want to compute undivided differences here
          Real fakeDx = 1.0;

          if (useCenterGradF)
            {
              FORT_GRADFACETANDIRS(CHF_FRA(tanGradF),
                                   CHF_CONST_FRA(thisFforGradDir),
                                   CHF_CONST_INT(varFStart),
                                   CHF_CONST_INT(numVar),
                                   CHF_CONST_INT(FSpaceStride),
                                   CHF_CONST_INT(a_NNumComp),
                                   CHF_BOX(faceBox),
                                   CHF_I1D(tanDirs, nTanDirs),
                                   CHF_CONST_REAL(fakeDx));
            }
          else
            {
              Box loBox, hiBox, centerBox, entireBox;
              int hasLo, hasHi;
              for (int i = 0; i != nTanDirs; ++i)
                {
                  const int tanDir = tanDirs[i];
                  Box stencilBox(a_box);
                  stencilBox.grow(tanDir, 1);
                  loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
                             stencilBox, *a_problemDomainPtr, tanDir);
                  // All boxes index face in direction 'dir'
                  centerBox.surroundingNodes(dir);
                  if (hasLo) loBox.surroundingNodes(dir); // lowest layer
                  if (hasHi) hiBox.surroundingNodes(dir); // highest layer
                  FORT_GRADFACETANDIR(CHF_FRA(tanGradF),
                                      CHF_CONST_FRA(thisFforGradDir),
                                      CHF_CONST_INT(varFStart),
                                      CHF_CONST_INT(numVar),
                                      CHF_CONST_INT(FSpaceStride),
                                      CHF_CONST_INT(a_NNumComp),
                                      CHF_BOX(loBox),
                                      CHF_CONST_INT(hasLo),
                                      CHF_BOX(hiBox),
                                      CHF_CONST_INT(hasHi),
                                      CHF_BOX(centerBox),
                                      CHF_CONST_INT(i),
                                      CHF_CONST_INT(tanDir),
                                      CHF_CONST_REAL(fakeDx));
                }
            }

          // For gradients of N, we pretend we are working on one variable of
          // a flux dyad with the dyad in C storage as (var, space).
          const int NStartVar = (a_NNumComp == SpaceDim) ? dir : 0;
          const int NNumVar = 1;
          const int NSpaceStride = 1;
          FORT_GRADFACETANDIRS(CHF_FRA(tanGradN),
                               CHF_CONST_FRA(thisNDir),
                               CHF_CONST_INT(NStartVar),
                               CHF_CONST_INT(NNumVar),
                               CHF_CONST_INT(NSpaceStride),
                               CHF_CONST_INT(a_NNumComp),
                               CHF_BOX(faceBox),
                               CHF_I1D(tanDirs, nTanDirs),
                               CHF_CONST_REAL(fakeDx));

          // We now have nicely ordered gradients.  Take the dot product
          // and add to the product

          // Number of terms in the dot product for each solution variables
          const int numDotTerms = nTanDirs*a_NNumComp;
          // Temporary storage of the dot product
          FArrayBox dotTanGrads(faceBox, 1);
          FORT_REDUCEDFDOTN(CHF_FRA(thisProductDir),
                            CHF_FRA1(dotTanGrads, 0),
                            CHF_CONST_INT(varProductStart),
                            CHF_CONST_INT(numVar),
                            CHF_CONST_FRA(tanGradF),
                            CHF_CONST_FRA(tanGradN),
                            CHF_CONST_INT(numDotTerms),
                            CHF_BOX(faceBox));
          // dummy statement in order to get around gdb bug
          int dummy_unused = 0; dummy_unused = 0;
        } // end if we're doing the 4th-order corrections
    } // end loop over face directions
}

/// computes 4th-order average of product = N^T*F
/** if a_fourthOrder is false, then only do a second-order dot product
 */
void
NewFourthOrderCoordSys::computeMetricTermProductAverage(FluxBox& a_product,
                                                        const FluxBox& a_F,
                                                        const FluxBox& a_N,
                                                        const FluxBox& a_FforGrad,
                                                        const Box& a_box,
                                                        bool a_fourthOrder) const
{
  CH_TIME("NewFourthOrderCoordSys::computeMetricTermProductAverage");

  int NNumComp = 0;
  // two possibilities here
  if (a_F.nComp() == SpaceDim)
    {
      // one vector component of F
      // in this case, we only use the face-normal terms here if
      // there are only SpaceDim components to F.
      NNumComp = SpaceDim;
    }
  else if (a_F.nComp() == a_N.nComp())
    {
      // however, if there are the same number of components in
      // F as in the faceMetricTerms, then do all of the components
      NNumComp = a_N.nComp();
    }
  // In the old routine, we could only process on variables
  Interval varInterval(0, 0);
  // An therefore, the stride through the space components of the flux dyad
  // (really a vector), is 1.
  int fluxSpaceContiguous = 1;
  computeMetricTermProductAverage(a_product,
                                  a_F,
                                  a_N,
                                  NNumComp,
                                  a_FforGrad,
                                  a_box,
                                  a_fourthOrder,
                                  varInterval,
                                  varInterval,
                                  fluxSpaceContiguous);
}

/// computes 4th-order average of product = N^T*F
/** if a_fourthOrder is false, then only do a second-order dot product
 */
void
NewFourthOrderCoordSys::computeMetricTermProductAverage(FluxBox& a_product,
                                                        const FluxBox& a_F,
                                                        const FluxBox& a_N,
                                                        const Box& a_box,
                                                        bool a_fourthOrder) const
    {
  computeMetricTermProductAverage(a_product, a_F, a_N, a_F, a_box,
                                  a_fourthOrder);
}

//**FIXME This is independent of coordSys and should probably be in
//**FourthOrderUtil instead -- except for getNcomponent which
//**seemingly should  be constant anyways.
/*--------------------------------------------------------------------*/
//   Computes magnitude of <N^T> on each face
/**  The normal components (row) of \f$<N^T>\f$ are included on each
 *   face resulting in the ratio of the area of the face in physical
 *   space to the area in computational space.
 *   \param[out] a_NMag Magnitude of N on each face
 *   \param[in]  a_N    Metrics component.  Stored on
 *                      getNcomponent(0, dir). for each direction.
 *   \param[in]  a_box  The box on whose faces the magnitude is
 *                      computed.
 *//*-----------------------------------------------------------------*/

void
NewFourthOrderCoordSys::magnitudeN(FluxBox&       a_NMag,
                                   const FluxBox& a_N,
                                   const Box&     a_box) const
{
  CH_assert(a_NMag.box().contains(a_box));
  CH_assert(a_N.box().contains(a_box));

//--Loop through the face directions

  for (int dir=0; dir != SpaceDim; ++dir)
    {
      Box faceBox(a_box);
      faceBox.surroundingNodes(dir);
      const int beginNDirComp = getNcomponent(0, dir);
      FORT_MAGNITUDEN(CHF_FRA1(a_NMag[dir], 0),
                      CHF_CONST_FRA(a_N[dir]),
                      CHF_CONST_INT(beginNDirComp),
                      CHF_BOX(faceBox));
    }
}

/// Evaluate \f$(D-1)\mathcal{N}_{d,d'}^s = \det\left[\left(\nabla_{\vec{\xi}}\vec{x}\right)^T\left(\vec{x}|d'\right)\left(s|\hat{e}^d\right)\right]\f$ at a point.
// this evaluates the script N values from equation 12 in Phil's
// notes; for the most general case, we do this by actually evaluating
// the determinant, which we do by doing Gaussian elimination to make the
// matrix trangular.  This is fairly inefficient, and should be revisited
// when we actually care about performance
Real
NewFourthOrderCoordSys::getN(const RealVect& a_Xi,
                             int a_s, int a_d, int a_d1) const
{
  // allocate temporaries
   Real polarity = 1.0;
  // this keeps track of rows before/after swapping
  int rowIndex[SpaceDim];
  D_TERM6(rowIndex[0] = 0;,
          rowIndex[1] = 1;,
          rowIndex[2] = 2;,
          rowIndex[3] = 3;,
          rowIndex[4] = 4;,
          rowIndex[5] = 5);

  // A[i][j] = A[row][column]
  Real A[SpaceDim][SpaceDim];

  // first, since we know that the s_th row is the dth unit vector, swap
  // the sth row with the dth row
  if (a_s != a_d)
    {
      rowIndex[a_s] = a_d;
      rowIndex[a_d] = a_s;
      polarity *= -1;
    }

  // now fill in diagonal values, since we know that we'll need them
  D_TERM6(
          A[0][0] = getNMatrixEntry(a_Xi, a_s, a_d, a_d1,
                                    rowIndex[0], 0);,
          A[1][1] = getNMatrixEntry(a_Xi, a_s, a_d, a_d1,
                                    rowIndex[1], 1);,
          A[2][2] = getNMatrixEntry(a_Xi, a_s, a_d, a_d1,
                                    rowIndex[2], 2);,
          A[3][3] = getNMatrixEntry(a_Xi, a_s, a_d, a_d1,
                                    rowIndex[3], 3);,
          A[4][4] = getNMatrixEntry(a_Xi, a_s, a_d, a_d1,
                                    rowIndex[4], 4);,
          A[5][5] = getNMatrixEntry(a_Xi, a_s, a_d, a_d1,
                                    rowIndex[5], 5) );

  // if we're in 2D, then no need to do any more elimination
#if CH_SPACEDIM > 2
  //  if (SpaceDim > 2)
    {
      // need to come back and redo this if SpaceDim > 3
      CH_assert(SpaceDim == 3);

      Real zeroVal = 1.0e-12;
      if (a_d == 0)
        {
          // eliminate A12
          A[1][2] = getNMatrixEntry(a_Xi, a_s, a_d, a_d1,
                                    rowIndex[1], 2);

          // if A22 is already zero, then can just swap rows 1 and 2
          if (abs(A[2][2]) < zeroVal)
            {
              A[2][2] = A[1][2];
              A[1][1] = getNMatrixEntry(a_Xi, a_s, a_d, a_d1,
                                        rowIndex[2], 1);

              polarity *= -1;
            }

          // don't need to do elimination if already zero
          else if (abs(A[1][2]) > zeroVal)
            {
              A[2][1] = getNMatrixEntry(a_Xi, a_s, a_d, a_d1,
                                        rowIndex[2], 1);


              A[1][1] = A[1][1] - A[2][1]*A[1][2]/A[2][2];
            }
        }
      else if (a_d == 1)
        {
          // eliminate A20
          // (note that since row 1 is (0 1 0), eliminating A21 is trivial
          A[2][0] = getNMatrixEntry(a_Xi, a_s, a_d, a_d1,
                                    rowIndex[2], 0);

          // if A00 is already zero, then can swap rows 0 and 2
          if (abs(A[0][0]) < zeroVal)
            {
              // note that A02 = -A20
              A[0][0] = A[2][0];
              A[2][2] = getNMatrixEntry(a_Xi, a_s, a_d, a_d1,
                                        rowIndex[0], 2);
              polarity *= -1;
            }
          // don't need to do elimination if already zero
          else if (abs(A[2][0]) > zeroVal)
            {
              A[0][2] = getNMatrixEntry(a_Xi, a_s, a_d, a_d1,
                                        rowIndex[0], 2);

              A[2][2] = A[2][2] - A[0][2]*A[2][0]/A[0][0];
            }
        }
      else if (a_d == 2)
        {
          // eliminate A10
          A[1][0] = getNMatrixEntry(a_Xi, a_s, a_d, a_d1,
                                    rowIndex[1], 0);

          // if A00 is already zero, then can just swap rows 0 and 1
          if (abs(A[0][0]) < zeroVal)
            {
              A[0][0] = A[1][0];
              A[1][1] = getNMatrixEntry(a_Xi, a_s, a_d, a_d1,
                                        rowIndex[0], 1);
              polarity *= -1;
            }
          // don't need to do elimination if already zero
          else if (abs(A[1][0]) > zeroVal)
            {
              A[0][1] = getNMatrixEntry(a_Xi, a_s, a_d, a_d1,
                                        rowIndex[0], 1);

              A[1][1] = A[1][1] - A[0][1]*A[1][0]/A[0][0];
            }
        }
      else
        {
          MayDay::Error("NewFourthOrderCoordSys::getN -- Bad value for a_d");
        }

    } // end if SpaceDim > 2
#endif

  // note that N blows up for Dim = 1
  // (come back and fix this later, if required)
#if CH_SPACEDIM > 1
  Real N = polarity*D_TERM6(A[0][0],*A[1][1],*A[2][2],*A[3][3],*A[4][4],*A[5][5])/(SpaceDim - 1);
#else
  Real N = 1.0;
  MayDay::Error("NewFourthOrderCoordSys::getN not defined for 1D");
#endif

  return N;
}

Real
NewFourthOrderCoordSys::getNMatrixEntry(const RealVect& a_Xi,
                                        int a_s, int a_d, int a_d1,
                                        int a_row, int a_column) const
{
  Real entry = 0.0;
  // simplest thing -- if row = s, then A_ij = delta_dj
  if (a_row == a_s)
    {
      if (a_column == a_d) entry = 1.0;
    }
  else if (a_column == a_d1)
    {
      entry = realCoord(a_Xi)[a_row];
    }
  else
    {
      entry = dXdXi(a_Xi, a_row, a_column);
    }

  return entry;
}

/*--------------------------------------------------------------------*/
///  Integrates \f$\mathcal(N)^s\f$ on a specific set of codimension 2
///  hyperedges that are orthogonal to the given directions.
/**  \param[out] a_scrN FArrayBox orientated on the hyperedge.  The
 *                      components are the quadratures of
 *                      \f$\mathcal(N)^s\f$ for each \f$s\f$.
 *   \param[in]  a_dir0 The first orthogonal direction specifying the
 *                      orientation of the hyperedge
 *   \param[in]  a_dir0 The second orthogonal direction specifying the
 *                      orientation of the hyperedge
 *   \param[in]  a_box  A box for iterating over a_scrN.  It must have
 *                      the same orientation as a_scrN to iterate over
 *                      the hyperedges orthogonal to a_dir0 and a_dir1
 *   \note
 *   <ul>
 *     <li> This is similar to routine
 *          'incrementFaceMetricWithEdgeTerm' except that the
 *          quadratures of \f$\mathcal(N)^s\f$ are saved on each
 *          hyperedge and there is a more general treatment of
 *          directions.
 *     <li> HE means hyperedge.
 *     <li> Ready for SpaceDim > 3 if the quadratures exist.
 *   </ul>
 *//*-----------------------------------------------------------------*/

void
NewFourthOrderCoordSys::integrateScriptN(FArrayBox& a_scrN,
                                         const int a_dir0,
                                         const int a_dir1,
                                         const Box& a_box) const
{
  BoxIterator bit(a_box);
  if (SpaceDim == 2)
    {
      for (bit.begin(); bit.ok(); ++bit)
        {
          const IntVect iv = bit();
          const RealVect HECenter = m_dx*iv;
          for (int s = 0; s != SpaceDim; ++s)
            a_scrN(iv, s) = getN(HECenter, s, a_dir0, a_dir1);
        }
    }
  else if (SpaceDim > 3)
    {
      MayDay::Error("NewFourthOrderCoordSys::integrateScriptN - Confirm "
                    "quadratures available for SpaceDim > 3");
    }
  else
    {
      // The offset from the iterator location to the HE center
      RealVect HEOffset = 0.5*m_dx;
      HEOffset[a_dir0] = 0.;
      HEOffset[a_dir1] = 0.;
      const Vector<QuadratureElement> quadPts =
        m_quadraturePtr->coefficients(a_dir0, a_dir1);
      Real weightMult = m_quadraturePtr->weightMult(m_dx, a_dir0, a_dir1);
      for (bit.begin(); bit.ok(); ++bit)
        {
          const IntVect iv = bit();
          const RealVect HECenter = m_dx*iv + HEOffset;
          for (int s = 0; s != SpaceDim; ++s)
            {
              a_scrN(iv, s) = 0.;
              for (int i = 0; i < quadPts.size(); ++i)
                {
                  RealVect mappedLoc(HECenter);
                  mappedLoc += quadPts[i].location*(m_dx/2.0);
                  Real scrNPoint = getN(mappedLoc, s, a_dir0, a_dir1);
                  a_scrN(iv, s) += scrNPoint*quadPts[i].weight;
                }
              a_scrN(iv, s) *= weightMult;
            }
        }
    }
}

void
NewFourthOrderCoordSys::incrementFaceMetricWithEdgeTerm(FArrayBox& a_N,
                                                        int a_faceDir,
                                                        int a_edgeDir,
                                                        const Box& a_box,
                                                        const Side::LoHiSide& a_side) const
{
  // make sure we have enough components
  CH_assert(a_N.nComp() >= SpaceDim);

  // this is where we do the integration along the edges
  Real mult = 1.0;
  if (a_side == Side::Lo) mult = -1.0;

  // this is the offset dx*(i,j,k) for the faces in a_faceMetrics
  RealVect faceOffset = m_dx;
  faceOffset *= 0.5;
  faceOffset[a_faceDir] = 0.0;

  RealVect edgeOffset = faceOffset;
  edgeOffset[a_edgeDir] += mult*0.5*m_dx[a_edgeDir];

  const Vector<QuadratureElement > quadPts = m_quadraturePtr->coefficients(a_faceDir,
                                                                           a_edgeDir);
  Real weightMult = m_quadraturePtr->weightMult(m_dx, a_faceDir, a_edgeDir);

  // this is gonna be slow, but we can hopefully make it faster at
  // some point down the road
  BoxIterator bit(a_box);

  for (bit.begin(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      RealVect faceCenter = m_dx*iv + faceOffset;

      RealVect edgeCenter = m_dx*iv + edgeOffset;

      for (int sDir = 0; sDir<SpaceDim; sDir++)
        {
          Real edgeVal = 0.0;

          // since 2D is just nodes, do this separately
          if (SpaceDim == 2)
            {
              //RealVect realLoc = realCoord(edgeCenter);
              edgeVal = getN(edgeCenter, sDir, a_faceDir, a_edgeDir);
            }
          else if (SpaceDim == 3)
            {
              // integrate along tangentDir using quadrature
              for (int i=0; i<quadPts.size(); i++)
                {
                  RealVect mappedLoc(edgeCenter);
                  mappedLoc += (quadPts[i].location)*m_dx/2.0;
                  // convert from mapped->real space
                  //RealVect realLoc = realCoord(mappedLoc);
                  Real Nvalue = getN(mappedLoc, sDir, a_faceDir, a_edgeDir);
                  edgeVal += Nvalue*quadPts[i].weight*weightMult;
                }

            }
          else
            {
              MayDay::Error("NewFourthOrderCoordSys::faceMetrics not defined for SpaceDim > 3");
            }

          // since this is for normal components, faceDir = d
          int component = getNcomponent(sDir,a_faceDir);
          a_N(iv, component) += mult*edgeVal;

        } // end loop over s directions
    } // end loop over faces in faceMetrics

}


void
NewFourthOrderCoordSys::computeTransverseFaceMetric(FArrayBox& a_faceMetrics,
                                                    const Box& a_box,
                                                    int a_faceDir,
                                                    int a_dDir) const
{
  // first, compute face area
  Real faceArea = m_dx[0];
  for (int dir=1; dir<SpaceDim; dir++)
    {
      faceArea *= m_dx[dir];
    }
  faceArea /= m_dx[a_faceDir];


  // temp storage for nodal values
  Box nodeBox(a_box);
  Box LapBox(nodeBox);
  // switch box to node centering and grow by one in the
  // transverse directions (Laplacian doesn't need the extra cells)
  for (int dir=0; dir<SpaceDim; dir++)
    {
      if (dir != a_faceDir)
        {
          nodeBox.surroundingNodes(dir);
          nodeBox.grow(dir, 1);
          LapBox.surroundingNodes(dir);
        }
    }

  FArrayBox nodeN(nodeBox, SpaceDim);


  // compute N at nodes
  BoxIterator bit(nodeBox);
  for (bit.begin(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      RealVect mappedLoc = m_dx*iv;
      RealVect realLoc = realCoord(mappedLoc);

      // now, compute derivatives for Jacobian matrix.  Note that
      // we will not need anything in the dDir column. so don't
      // bother to compute them. All other entries wind up being used 2x
      // (which is why we're precomputing them)
      // do this with an eye toward using the pointwiseJ function
      // NOTE -- since eqn(6) in Phil's notes is actually for N^T instead
      // of N, use transpose of Jacobian matrix (dXdXi(row,col) instead
      // of dxdXi(col,row) as one might expect.  This results in computing
      // (N^T)^T, which is, of course, N.
      Vector<Real> Jmatrix(SpaceDim*SpaceDim, 0.0);

      for (int col=0; col<SpaceDim; col++)
        {
          if (col != a_dDir)
            {
              for (int row=0; row<SpaceDim; row++)
                {
                  Jmatrix[row + SpaceDim*col] = dXdXi(mappedLoc, row, col);
                }
            }
        } // end loop over columns for jacobian matrix computation


      for (int sDir=0; sDir<SpaceDim; sDir++)
        {
          // make use of the computePointwiseJ function to take the
          // determinant
          // first, copy the derivative values
          Vector<Real> jmtrx = Jmatrix;
          // now, replace the sDir'th row with the dth unit vector
          for (int j=0; j<SpaceDim; j++)
            {
              if (j == a_dDir)
                {
                  jmtrx[sDir + SpaceDim*j] = 1.0;
                }
              else
                {
                  jmtrx[sDir + SpaceDim*j] = 0.0;
                }
            }
          Real determinant = computePointwiseJ(jmtrx);
          nodeN(iv,sDir) = determinant;
        }
    } // end loop over cells


  // compute node-centered tangential Laplacian in fortran
  FArrayBox lapN(LapBox,SpaceDim);
  lapN.setVal(0.0);

  for (int dir=0; dir<SpaceDim; dir++)
    {
      if (dir != a_faceDir)
        {
          // simple undivided difference approximation to Laplacian
          Real factor = 1.0;
          FORT_INCREMENTLAPLACIAN(CHF_FRA(lapN),
                                  CHF_CONST_FRA(nodeN),
                                  CHF_BOX(LapBox),
                                  CHF_CONST_INT(dir),
                                  CHF_CONST_REAL(factor));
        }
    }

  // now average to cell centers from nodes.

  // this box represents which nodes will be averaged to
  // the face center
  Box averageBox(IntVect::Zero, IntVect::Unit);
  averageBox.setBig(a_faceDir, 0);
  BoxIterator averageIterator(averageBox);

  Real averageFactor = D_TERM6(1.0, *0.5, *0.5, *0.5, *0.5, *0.5);

  // this is another candidate for fortran...
  BoxIterator faceBit(a_faceMetrics.box());
  for (faceBit.begin(); faceBit.ok(); ++faceBit)
    {
      IntVect faceIV = faceBit();
      // start by just computing simple average to face
      // do this by looping through the averageBox
      for (averageIterator.begin(); averageIterator.ok(); ++averageIterator)
        {
          IntVect nodeLoc = faceIV + averageIterator();
          for (int sDir=0; sDir<SpaceDim; sDir++)
            {
              int metricComp = getNcomponent(sDir, a_dDir);
              a_faceMetrics(faceIV,metricComp) += averageFactor*nodeN(nodeLoc,sDir);
              // increment with Laplacian term as well
              // this is a fourth-order point value
              //a_faceMetrics(faceIV,metricComp) -= 0.125*averageFactor*lapN(nodeLoc,sDir);
              // this is a fourth-order face average
              a_faceMetrics(faceIV,metricComp) -= averageFactor*lapN(nodeLoc,sDir)/12.0;

            }
        }

    }
  // multiply by face area to get the integral over the face
  for (int sDir=0; sDir<SpaceDim; sDir++)
    {
      int metricComp = getNcomponent(sDir, a_dDir);
      //
      a_faceMetrics.mult(faceArea, metricComp, 1);
    }
}

inline Real computePointwiseJ( const Vector<Real>& Jmtrx )
{
   CH_assert(Jmtrx.size()==SpaceDim*SpaceDim);
   Real detJ = 0.0;
   if ( SpaceDim == 1 )
   {
      detJ = Jmtrx[0];
   }
   else if ( SpaceDim == 2 )
   {
      detJ = Jmtrx[0] * Jmtrx[3] - Jmtrx[1] * Jmtrx[2];
   }
   else if ( SpaceDim == 3 )
   {
      detJ = Jmtrx[0] * ( Jmtrx[4] * Jmtrx[8] - Jmtrx[5] * Jmtrx[7] )
           + Jmtrx[1] * ( Jmtrx[5] * Jmtrx[6] - Jmtrx[3] * Jmtrx[8] )
           + Jmtrx[2] * ( Jmtrx[3] * Jmtrx[7] - Jmtrx[4] * Jmtrx[6] );
   }
   else
   {
      MayDay::Error("NewFourthOrderCoordSys::pointwiseJ not defined for SpaceDim > 3");
   }
   return detJ;
}

void
NewFourthOrderCoordSys::computeTangentialGrad(FluxBox& a_gradPhi,
                                              const FluxBox& a_phiFace,
                                              const Box& a_box) const
{
  // want to compute undivided differences here
  Real fakeDx = 1.0;
  int nPhiComp = a_phiFace.nComp();
  CH_assert(a_gradPhi.nComp() >= (SpaceDim-1)*nPhiComp);

  //FluxBox& thisGradPhi = a_gradPhi[dit];
  //const FluxBox& thisPhi = a_phiFace[dit];
  for (int faceDir=0; faceDir<SpaceDim; faceDir++)
    {
      FArrayBox& thisGradPhiDir = a_gradPhi[faceDir];
      const FArrayBox& thisPhiDir = a_phiFace[faceDir];

      for (int n=0; n<nPhiComp; n++)
        {
          // now loop over tangential directions
          for (int gradDir=0; gradDir<SpaceDim; gradDir++)
            {
              if (gradDir != faceDir)
                {
                  int gradComp = tanGradComp(faceDir, gradDir, n);

                  Box gradBox(thisPhiDir.box());
                  gradBox.grow(gradDir, -1);
                  gradBox &= thisGradPhiDir.box();
                  FORT_TANGRADFACE(CHF_FRA1(thisGradPhiDir, gradComp),
                                   CHF_CONST_FRA1(thisPhiDir, n),
                                   CHF_BOX(gradBox),
                                   CHF_INT(gradDir),
                                   CHF_CONST_REAL(fakeDx));

                } // end if gradDir is a tangential dir
            } // end loop over grad directions
        } // end loop over components
    } // end loop over face directions

}


Real
NewFourthOrderCoordSys::pointwiseJ(const RealVect& a_Xi) const
{
   // this needs to be improved and generalized
   Vector<Real> Jmtrx(SpaceDim*SpaceDim);
   for (int j=0; j<SpaceDim; j++)
   {
      for (int i=0; i<SpaceDim; i++)
      {
         Jmtrx[i+SpaceDim*j] = dXdXi( a_Xi, i, j );
      }
   }
   return computePointwiseJ( Jmtrx );
}


int
NewFourthOrderCoordSys::tanGradComp(const int a_faceDir,
                                    const int a_tanDir,
                                    const int a_comp) const
{
  CH_assert(a_faceDir != a_tanDir);

  int indexVal;
  if (a_tanDir < a_faceDir)
    {
      indexVal = (SpaceDim-1)*a_comp + a_tanDir;
    }
  else
    {
      indexVal = (SpaceDim-1)*a_comp + a_tanDir -1;
    }
  return indexVal;
}


void
NewFourthOrderCoordSys::contravariantMetric(FArrayBox& a_metric,
                                            int a_dir) const
{
  MayDay::Error("NewFourthOrderCoordSys::contravariantMetric not implemented yet");
}


void
NewFourthOrderCoordSys::orthonormalize(FluxBox& a_flux,
                                       const Interval& a_vectorIntv) const
{
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      orthonormalize(a_flux[idir], a_vectorIntv);
    }
}


void
NewFourthOrderCoordSys::orthonormalize(FArrayBox& a_fluxFab,
                                       const Interval& a_vectorIntv) const
{
  if (a_vectorIntv.size() > 0)
    {
      CH_assert(a_vectorIntv.size() == SpaceDim);
      IntVect csComps, orthoComps;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          csComps[idir] = a_vectorIntv.begin() + idir;
          orthoComps[idir] = idir;
        }
      const Box& bx = a_fluxFab.box();
      int idir = faceDimension(bx);
      FArrayBox orthoFab(bx, SpaceDim);
      orthonormalize(a_fluxFab, orthoFab, bx, idir, csComps, orthoComps);
      // Overwrite the vector components with deorthonormalized FAB.
      a_fluxFab.copy(orthoFab, 0, a_vectorIntv.begin(), SpaceDim);
    }
}


void
NewFourthOrderCoordSys::orthonormalize(
                                       const FArrayBox& a_csFab,
                                       FArrayBox& a_orthoFab,
                                       const Box& a_box,
                                       int a_idir,
                                       const IntVect& a_csComps,
                                       const IntVect& a_orthoComps) const
{
  // Convert over a_box, on lines where direction a_idir is constant,
  // from a_csFab[a_csComps[d]] in coordinate-system basis
  // to a_orthoFab[a_orthoComps[d]] in orthonormal basis.

  // FIXME: Now this is just a copy; correct it later.
  for (int d = 0; d < SpaceDim; d++)
    {
      a_orthoFab.copy(a_csFab, a_box, a_csComps[d], a_box, a_orthoComps[d], 1);
    }
}


void
NewFourthOrderCoordSys::deorthonormalize(FluxBox& a_flux,
                                         const Interval& a_vectorIntv) const
{
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      deorthonormalize(a_flux[idir], a_vectorIntv);
    }
}


void
NewFourthOrderCoordSys::deorthonormalize(FArrayBox& a_fluxFab,
                                         const Interval& a_vectorIntv) const
{
  if (a_vectorIntv.size() > 0)
    {
      CH_assert(a_vectorIntv.size() == SpaceDim);
      IntVect orthoComps, csComps;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          orthoComps[idir] = a_vectorIntv.begin() + idir;
          csComps[idir] = idir;
        }
      const Box& bx = a_fluxFab.box();
      int idir = faceDimension(bx);
      FArrayBox csFab(bx, SpaceDim);
      deorthonormalize(a_fluxFab, csFab, bx, idir, orthoComps, csComps);
      // Overwrite the vector components with deorthonormalized FAB.
      a_fluxFab.copy(csFab, 0, a_vectorIntv.begin(), SpaceDim);
    }
}


void
NewFourthOrderCoordSys::deorthonormalize(
                                         const FArrayBox& a_orthoFab,
                                         FArrayBox& a_csFab,
                                         const Box& a_box,
                                         int a_idir,
                                         const IntVect& a_orthoComps,
                                         const IntVect& a_csComps) const
{
  // Convert over a_box, on lines where direction a_idir is constant,
  // from a_orthoFab[a_orthoComps[d]] in orthonormal basis
  // to a_csFab[a_csComps[d]] in coordinate-system basis.

  // FIXME: Now this is just a copy; correct it later.
  for (int d = 0; d < SpaceDim; d++)
    {
      a_csFab.copy(a_orthoFab, a_box, a_orthoComps[d], a_box, a_csComps[d], 1);
    }
}

#include "NamespaceFooter.H"

