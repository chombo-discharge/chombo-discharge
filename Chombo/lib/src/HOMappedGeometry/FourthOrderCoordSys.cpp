#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "FourthOrderCoordSys.H"
#include "FourthOrderCoordSysF_F.H"
#include "FourthOrderUtil.H"
#include "FourthOrderUtilF_F.H"
#include "FCDivergenceF_F.H"
#include "FaceGradientF_F.H"
#include "BoxIterator.H"
#include "GaussianQuadrature.H"
#include "NewtonCotesQuadrature.H"
#include "PointwiseDotProdF_F.H"
#include "FourthOrderUtil.H"
#include "NamespaceHeader.H"

inline Real computePointwiseJ( const Vector<Real>& Jmtrx );

FourthOrderCoordSys::FourthOrderCoordSys()
{
  m_quadraturePtr = NULL;
  // for now, set default quadrature to be 3-pt Gaussian
  //int numPts = 1;
  int numPts = 4;
  //m_quadraturePtr = new GaussianQuadrature(numPts);
  m_quadraturePtr = new NewtonCotesQuadrature(numPts);

  // default volume interval is all components
  m_volInterval.define(0, SpaceDim-1);

}

FourthOrderCoordSys::~FourthOrderCoordSys()
{
  if (m_quadraturePtr != NULL)
    {
      delete m_quadraturePtr;
      m_quadraturePtr = NULL;
    }
}


void
FourthOrderCoordSys::define(const DisjointBoxLayout& a_grids,
                            const ProblemDomain& a_domain,
                            const RealVect& a_cellSpacing,
                            const IntVect& a_ghostVect)
{
  m_grids = a_grids;
  m_dx = a_cellSpacing;
  m_domain = a_domain;
  m_grids = a_grids;

  m_ghostVect = a_ghostVect;

  // in order to compute fourth-order face-averaged Jinverse*N,
  // need to have Jinverse with one extra ghost cell, which means
  // that we need to have cell volumes, etc with 2 extra ghost cells.
  // (or could use a different averaging formula)
  IntVect grownGhost1(a_ghostVect);
  grownGhost1 += IntVect::Unit;
  IntVect grownGhost2(a_ghostVect);
  grownGhost2 += 2*IntVect::Unit;
  IntVect grownGhost3(a_ghostVect);
  grownGhost3 += 3*IntVect::Unit;

  int nFaceMetricTerms = SpaceDim*SpaceDim;
  m_faceMetricTerms.define(a_grids, nFaceMetricTerms, grownGhost3);
  //int nGradNTerms = (SpaceDim)*(SpaceDim -1);
  int nGradNTerms = (nFaceMetricTerms)*(SpaceDim -1);
  // kludge to prevent crashing in 1D, not that m_tanGradNTerms will be used.
  if (SpaceDim == 1) nGradNTerms = 1;
  m_tanGradNTerms.define(a_grids, nGradNTerms, grownGhost2);

  // need cell volumes defined with at least two ghost cells in order to compute
  // fourth-order face-centered Jinverse
  m_cellVolumes.define(a_grids, 1, grownGhost2);
  m_J.define(a_grids, 1, grownGhost2);

  // m_Jinverse needs an extra ghost cell in order to compute gradients for NinverseJ
  m_JInverse.define(a_grids, 1, grownGhost1);

  m_NinverseJ.define(a_grids, nFaceMetricTerms, a_ghostVect);

  // this is where we compute the face-centered metric terms
  // from eqn 11 in phil's document

  defineFaceMetricTerms(m_faceMetricTerms);

  // now compute face-centered tangential gradients of the metric terms,
  // G_perp<N_d>

  computeTangentialGrad(m_tanGradNTerms, m_faceMetricTerms);

  computeCellVolumes(m_cellVolumes);

  computeJ(m_J);

  computeJinverse(m_JInverse);

  computeNJinverse(m_NinverseJ);
}

// compute realCoord on an FArrayBox. Default implementation
// the intent here is that derived classes can build more efficient
// implementations which call fortran
void
FourthOrderCoordSys::realCoord(FArrayBox& a_x,
                               const FArrayBox& a_Xi) const
{
  // start with assumption of a cell-centred fab
  RealVect offset = 0.5*RealVect::Unit;
  offset *= m_dx;
  for (int dir=0; dir<SpaceDim; dir++)
    {
      if (a_Xi.box().type(dir) == IndexType::NODE)
        {
          offset[dir] = 0.0;
        }
    }

  BoxIterator bit(a_x.box());
  for (bit.begin(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      RealVect mappedLoc = m_dx*iv + offset;
      RealVect realLoc = realCoord(mappedLoc);
      D_TERM6(a_x(iv,0) = realLoc[0];,
              a_x(iv,1) = realLoc[1];,
              a_x(iv,2) = realLoc[2];,
              a_x(iv,3) = realLoc[3];,
              a_x(iv,4) = realLoc[4];,
              a_x(iv,5) = realLoc[5];)
        }


}

void
FourthOrderCoordSys::mappedCoord(FArrayBox& a_Xi, const FArrayBox& a_x) const
{
  MayDay::Error("FArraybox version of FourthOrderCoorSys::mappedCoord not implemented");
}

void
FourthOrderCoordSys::pointwiseJ(FArrayBox& a_J) const
{
  const Box& Jbox = a_J.box();
  // start with assumption of a cell-centred fab
  RealVect offset = 0.5*RealVect::Unit;
  offset *= m_dx;
  for (int dir=0; dir<SpaceDim; dir++)
    {
      if (Jbox.type(dir) == IndexType::NODE)
        {
          offset[dir] = 0.0;
        }
    }

  BoxIterator bit(Jbox);
  for (bit.begin(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      RealVect mappedLoc = m_dx*iv + offset;
      RealVect realLoc;
      realLoc = realCoord(mappedLoc);
      a_J(iv,0) = pointwiseJ(realLoc);
    }

}

void
FourthOrderCoordSys::pointwiseJ(FArrayBox& a_J, const FArrayBox& a_Xi) const
{
  CH_assert(a_Xi.nComp() == SpaceDim);

  // start with assumption of a cell-centred fab
  RealVect offset = 0.5*RealVect::Unit;
  offset *= m_dx;
  for (int dir=0; dir<SpaceDim; dir++)
    {
      if (a_Xi.box().type(dir) == IndexType::NODE)
        {
          offset[dir] = 0.0;
        }
    }

  BoxIterator bit(a_Xi.box());
  for (bit.begin(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      RealVect mappedLoc, realLoc;
      D_TERM6(mappedLoc[0] = a_Xi(iv,0);,
              mappedLoc[1] = a_Xi(iv,1);,
              mappedLoc[2] = a_Xi(iv,2);,
              mappedLoc[3] = a_Xi(iv,3);,
              mappedLoc[4] = a_Xi(iv,4);,
              mappedLoc[5] = a_Xi(iv,5);)

      realLoc = realCoord(mappedLoc);
      a_J(iv,0) = pointwiseJ(realLoc);
    }

}

void
FourthOrderCoordSys::mappedGridDivergence(LevelData<FArrayBox>& a_divF,
                                          const LevelData<FluxBox>& a_F) const
{
  int nFluxComp = a_F.nComp()/SpaceDim;
  LevelData<FluxBox> tempFlux(a_F.getBoxes(), nFluxComp,
                              a_F.ghostVect());

  computeMetricTermProductAverage(tempFlux, a_F);

  simpleDivergence(a_divF, tempFlux);
}

void
FourthOrderCoordSys::simpleDivergence(LevelData<FArrayBox>& a_divF,
                                      const LevelData<FluxBox>& a_F) const
{
  // now compute divergence in the usual way
  //const DisjointBoxLayout& grids = a_divF.getBoxes();

  DataIterator dit = a_divF.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& thisDiv = a_divF[dit];
      const FluxBox& thisFlux = a_F[dit];

      // box over which we can compute div
      Box thisBox = thisDiv.box();
      thisBox &= thisFlux.box();

      // this is designed to ensure that any cells in divF which
      // don't contain the divergence are set to a bogus value
      thisDiv.setVal(1.23456e10);

      // first, set divF to 0
      thisDiv.setVal(0.0, thisBox, 0, thisDiv.nComp());
      // since we're computing volume*div, and the fluxes
      // are multiplied by the appropriate areas, don't need
      // a dx here.
      Real fakeDx = 1.0;

      //const Box& thisBox = thisDiv.box();
#if 0
      const Box& thisBox = grids[dit];
      if (thisBox != thisDiv.box())
        {
          MayDay::Warning("computing Div over grid box instead of divBox");
        }
#endif

      // now loop over directions and increment with directional
      // derivative
      for (int dir=0; dir<SpaceDim; dir++)
        {
          // use fortran from divergence operator here
          FORT_FCDIVERGENCE(CHF_CONST_FRA(thisFlux[dir]),
                            CHF_FRA(thisDiv),
                            CHF_BOX(thisBox),
                            CHF_CONST_REAL(fakeDx),
                            CHF_INT(dir));

        }
    }

}

void
FourthOrderCoordSys::setQuadrature(const EdgeQuadrature* a_quadrature)
{
  if (m_quadraturePtr != NULL)
    {
      delete m_quadraturePtr;
      m_quadraturePtr = NULL;
    }
  m_quadraturePtr = a_quadrature->new_quadrature();
}

const LevelData<FluxBox>&
FourthOrderCoordSys::getFaceMetricTerms() const
{
  return m_faceMetricTerms;
}

const
LevelData<FArrayBox>&
FourthOrderCoordSys::getCellVolumes() const
{
  return m_cellVolumes;
}

const
LevelData<FArrayBox>&
FourthOrderCoordSys::getJ() const
{
  return m_J;
}

const
LevelData<FluxBox>&
FourthOrderCoordSys::getJInverse() const
{
  return m_JInverse;
}

const
LevelData<FluxBox>&
FourthOrderCoordSys::getNJinverse() const
{
  return m_NinverseJ;
}

// this evaluates the script N values from equation 12 in Phil's
// notes; for the most general case, we do this by actually evaluating
// the determinant, which we do by doing Gaussian elimination to make the
// matrix trangular.  This is fairly inefficient, and should be revisited
// when we actually care about performance
Real
FourthOrderCoordSys::getN(const RealVect& a_X,
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
          A[0][0] = getNMatrixEntry(a_X, a_s, a_d, a_d1,
                                    rowIndex[0], 0);,
          A[1][1] = getNMatrixEntry(a_X, a_s, a_d, a_d1,
                                    rowIndex[1], 1);,
          A[2][2] = getNMatrixEntry(a_X, a_s, a_d, a_d1,
                                    rowIndex[2], 2);,
          A[3][3] = getNMatrixEntry(a_X, a_s, a_d, a_d1,
                                    rowIndex[3], 3);,
          A[4][4] = getNMatrixEntry(a_X, a_s, a_d, a_d1,
                                    rowIndex[4], 4);,
          A[5][5] = getNMatrixEntry(a_X, a_s, a_d, a_d1,
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
          A[1][2] = getNMatrixEntry(a_X, a_s, a_d, a_d1,
                                    rowIndex[1], 2);

          // if A22 is already zero, then can just swap rows 1 and 2
          if (abs(A[2][2]) < zeroVal)
            {
              A[2][2] = A[1][2];
              A[1][1] = getNMatrixEntry(a_X, a_s, a_d, a_d1,
                                        rowIndex[2], 1);

              polarity *= -1;
            }

          // don't need to do elimination if already zero
          else if (abs(A[1][2]) > zeroVal)
            {
              A[2][1] = getNMatrixEntry(a_X, a_s, a_d, a_d1,
                                        rowIndex[2], 1);


              A[1][1] = A[1][1] - A[2][1]*A[1][2]/A[2][2];
            }
        }
      else if (a_d == 1)
        {
          // eliminate A20
          // (note that since row 1 is (0 1 0), eliminating A21 is trivial
          A[2][0] = getNMatrixEntry(a_X, a_s, a_d, a_d1,
                                    rowIndex[2], 0);

          // if A00 is already zero, then can swap rows 0 and 2
          if (abs(A[0][0]) < zeroVal)
            {
              // note that A02 = -A20
              A[0][0] = A[2][0];
              A[2][2] = getNMatrixEntry(a_X, a_s, a_d, a_d1,
                                        rowIndex[0], 2);
              polarity *= -1;
            }
          // don't need to do elimination if already zero
          else if (abs(A[2][0]) > zeroVal)
            {
              A[0][2] = getNMatrixEntry(a_X, a_s, a_d, a_d1,
                                        rowIndex[0], 2);

              A[2][2] = A[2][2] - A[0][2]*A[2][0]/A[0][0];
            }
        }
      else if (a_d == 2)
        {
          // eliminate A10
          A[1][0] = getNMatrixEntry(a_X, a_s, a_d, a_d1,
                                    rowIndex[1], 0);

          // if A00 is already zero, then can just swap rows 0 and 1
          if (abs(A[0][0]) < zeroVal)
            {
              A[0][0] = A[1][0];
              A[1][1] = getNMatrixEntry(a_X, a_s, a_d, a_d1,
                                        rowIndex[0], 1);
              polarity *= -1;
            }
          // don't need to do elimination if already zero
          else if (abs(A[1][0]) > zeroVal)
            {
              A[0][1] = getNMatrixEntry(a_X, a_s, a_d, a_d1,
                                        rowIndex[0], 1);

              A[1][1] = A[1][1] - A[0][1]*A[1][0]/A[0][0];
            }
        }
      else
        {
          MayDay::Error("FourthOrderCoordSys::getN -- Bad value for a_d");
        }

    } // end if SpaceDim > 2
#endif

  // note that N blows up for Dim = 1
  // (come back and fix this later, if required)
#if CH_SPACEDIM > 1
  Real N = polarity*D_TERM6(A[0][0],*A[1][1],*A[2][2],*A[3][3],*A[4][4],*A[5][5])/(SpaceDim - 1);
#else
  Real N = 1.0;
  MayDay::Error("FourthOrderCoordSys::getN not defined for 1D");
#endif

  return N;
}

Real
FourthOrderCoordSys::getNMatrixEntry(const RealVect& a_X,
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
      entry = a_X[a_row];
    }
  else
    {
      entry = dXdXi(a_X, a_row, a_column);
    }

  return entry;
}


void
FourthOrderCoordSys::computeProductAverage(LevelData<FluxBox>& a_Product,
                                           const LevelData<FluxBox>& a_F,
                                           const LevelData<FluxBox>& a_G,
                                           bool a_fourthOrder) const
{

  int nFlux = a_Product.nComp();
  const DisjointBoxLayout& grids = a_F.getBoxes();

  CH_assert(a_G.getBoxes() == grids);
  CH_assert(a_F.nComp() == nFlux);
  CH_assert(a_G.nComp() == nFlux);

  int numGradComps = nFlux*(SpaceDim-1);
  IntVect grownGhost2(m_ghostVect);
  grownGhost2 += 2*IntVect::Unit;

  // temp storage -- same as metric term gradients
  // probably eventually want to only do this if a_fourthOrder is true


  LevelData<FluxBox> tanGradF(grids,
                              numGradComps,
                              grownGhost2);

  LevelData<FluxBox> tanGradG(grids,
                              numGradComps,
                              grownGhost2);

  LevelData<FluxBox> dotTanGrads(grids,
                                 nFlux,
                                 grownGhost2);

  if (a_fourthOrder)
    {
      // compute tangential gradients of F and G
      computeTangentialGrad(tanGradF, a_F);

      computeTangentialGrad(tanGradG, a_G);
    }

  DataIterator dit = grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const FluxBox& thisF = a_F[dit];
      const FluxBox& thisG = a_G[dit];
      FluxBox& thisFlux = a_Product[dit];

      FluxBox& thisGradF = tanGradF[dit];
      FluxBox& thisGradG = tanGradG[dit];
      FluxBox& thisDotGrad = dotTanGrads[dit];

      for (int dir=0; dir<SpaceDim; dir++)
        {
          FArrayBox& thisFluxDir = thisFlux[dir];
          const FArrayBox& thisFDir = thisF[dir];
          const FArrayBox& thisGDir = thisG[dir];

          Box intersectBox(thisFluxDir.box());
          intersectBox &= thisGDir.box();
          FArrayBox& thisGradFDir = thisGradF[dir];
          FArrayBox& thisGradGDir = thisGradG[dir];
          FArrayBox& thisDotGradDir = thisDotGrad[dir];

          // first part is easy -- dot product of F and G
          for (int comp=0; comp<nFlux; comp++)
            {
              int nDotProd = 1;
              FORT_POINTDOTPROD(CHF_FRA1(thisFluxDir, comp),
                                CHF_CONST_FRA(thisFDir),
                                CHF_CONST_INT(comp),
                                CHF_CONST_FRA(thisGDir),
                                CHF_CONST_INT(comp),
                                CHF_CONST_INT(nDotProd),
                                CHF_BOX(intersectBox) );

              if (a_fourthOrder)
                {
                  // now compute correction term
                  // first do pointwise dot products of gradF and gradG

                  // this little bit of wierdness is due to the fact that
                  // the tanGrad direction is not allowed to be the same as the
                  // face direction (i.e. dir can't equal nGradDir)
                  int nGradDir = 0;
                  if (dir == 0) nGradDir = 1;
                  int fstartcomp = tanGradComp(dir,nGradDir,comp);
                  int numcomp = (SpaceDim-1);

                  FORT_POINTWISEREDUCEDDOTPROD(CHF_FRA1(thisDotGradDir,comp),
                                               CHF_CONST_FRA(thisGradFDir),
                                               CHF_CONST_INT(fstartcomp),
                                               CHF_CONST_FRA(thisGradGDir),
                                               CHF_CONST_INT(fstartcomp),
                                               CHF_CONST_INT(numcomp),
                                               CHF_BOX(thisDotGradDir.box()) );

                } // end if incorporating 4th-order corrections
            } // end loop over components

          if (a_fourthOrder)
             {
                // since we computed undivided differences for the tanGrad's,
                // dx should also be one
                Real mult = 1.0/12.0;

                thisFluxDir.plus(thisDotGradDir, mult);
             } // end if incorporating 4th-order corrections

        } // end loop over face directions
    } // end loop over grids
}


// this function computes face-centered average of N^T*F
void
FourthOrderCoordSys::computeMetricTermProductAverage(LevelData<FluxBox>& a_Product,
                                                     const LevelData<FluxBox>& a_F,
                                                     bool a_fourthOrder) const
{

  int nFlux = a_Product.nComp();
  const DisjointBoxLayout& grids = a_F.getBoxes();

  const LevelData<FluxBox>& faceMetricTerms = getFaceMetricTerms();


  // temp storage
  LevelData<FluxBox> tanGradF(m_tanGradNTerms.getBoxes(),
                              m_tanGradNTerms.nComp(),
                              m_tanGradNTerms.ghostVect());

  LevelData<FluxBox> dotTanGrads(m_tanGradNTerms.getBoxes(),
                                 nFlux,
                                 m_tanGradNTerms.ghostVect());

#if 0
  // for debugging purposes
  for (DataIterator dit=grids.dataIterator(); dit.ok(); ++dit)
    {
      tanGradF[dit].setVal(0.0);
      dotTanGrads[dit].setVal(0.0);
    }
#endif

  if (a_fourthOrder)
    {
      // compute tangential gradients of F
      computeTangentialGrad(tanGradF, a_F);
    }

  DataIterator dit = grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const FluxBox& thisF = a_F[dit];
      FluxBox& thisFlux = a_Product[dit];

      // precomputed face-centered metric terms
      const FluxBox& thisFaceMetric = faceMetricTerms[dit];

      const FluxBox& thisGradN = m_tanGradNTerms[dit];
      FluxBox& thisGradF = tanGradF[dit];
      FluxBox& thisDotGrad = dotTanGrads[dit];


      for (int dir=0; dir<SpaceDim; dir++)
        {
          FArrayBox& thisFluxDir = thisFlux[dir];
          const FArrayBox& thisFDir = thisF[dir];
          const FArrayBox& thisFaceMetricDir = thisFaceMetric[dir];

          Box intersectBox(thisFluxDir.box());
          intersectBox &= thisFaceMetricDir.box();

          // first part is easy -- dot product of F and metrics
          int metricStartComp = -1;
          int fStartComp = 0;
          int ncomp = -1;

          // two possibilities here
          if (thisFDir.nComp() == SpaceDim)
            {
              // one vector component of F
              // in this case, we only use the face-normal terms here if
              // there are only SpaceDim components to F.
              // note that since we're multiplying by N^T instead of N.
              // we're doing (0,dir) instead of (dir,0).
              // because of the way that the metric terms are indexed,
              // metricStartComp will multiply by F1,
              // metricStartComp+1 will multiply by F2, etc...
              metricStartComp = getMetricTermComponent(0,dir);
              ncomp = SpaceDim;
            }
          // however, if there are the same number of components in
          // F as in the faceMetricTerms, then do all of the components
          else if (thisFDir.nComp() == thisFaceMetricDir.nComp())
            {
              metricStartComp = 0;
              ncomp = thisFDir.nComp();
            }

          FORT_POINTDOTPROD(CHF_FRA1(thisFluxDir, 0),
                            CHF_CONST_FRA(thisFDir),
                            CHF_CONST_INT(fStartComp),
                            CHF_CONST_FRA(thisFaceMetricDir),
                            CHF_CONST_INT(metricStartComp),
                            CHF_CONST_INT(ncomp),
                            CHF_BOX(intersectBox) );

          if (a_fourthOrder)
            {
              // now compute correction term
              // first do pointwise dot products of gradF and gradN
              const FArrayBox& thisGradNDir = thisGradN[dir];
              FArrayBox& thisGradFDir = thisGradF[dir];
              FArrayBox& thisDotGradDir = thisDotGrad[dir];
              int fGradStartComp = 0;
              // this little bit of wierdness is due to the fact that
              // the tanGrad direction is not allowed to be the same as the
              // face direction (i.e. dir can't equal nGradDir)
              int nGradDir = 0;
              if (dir == 0) nGradDir = 1;
              int nGradStartComp = tanGradComp(dir, nGradDir, metricStartComp);

              int nGradComp = ncomp*(SpaceDim-1);

              FORT_POINTWISEREDUCEDDOTPROD(CHF_FRA1(thisDotGradDir, 0),
                                           CHF_CONST_FRA(thisGradNDir),
                                           CHF_CONST_INT(nGradStartComp),
                                           CHF_CONST_FRA(thisGradFDir),
                                           CHF_CONST_INT(fGradStartComp),
                                           CHF_CONST_INT(nGradComp),
                                           CHF_BOX(thisDotGradDir.box()) );

              // since we computed undivided differences for the tanGrad's,
              //dx should also be one
              Real mult = 1.0/12.0;

              thisFluxDir.plus(thisDotGradDir, mult);
            } // end if we're doing the 4th-order corrections
        } // end loop over face directions
    } // end loop over grids
}


/// note that a_X is in real (not mapped) space.
/** fills the destComp component of a_dxdXi
    with the derivative of x w/ respect to Xi
    in the dirX direction
    Note that the default implementation of this in
    FourthOrderCoordSys class throws an error, since there is no
    way to get real=space coordinate x from index location except in
    the derived class..
    (default implementation is there in case derived class doesn't
    provide it)
*/
void
FourthOrderCoordSys::dXdXi(FArrayBox& a_dxdXi,
                           const FArrayBox& a_X,
                           int a_destComp,
                           int a_dirX,
                           int a_dirXi) const
{
  // if you're here, somebody did something wrong...
  MayDay::Error("default FourthOrderCoordsys::dXdXi not implemented");
}



void
FourthOrderCoordSys::defineFaceMetricTerms(LevelData<FluxBox>& a_faceMetricTerms)
{
  // this is where we compute the face-centered metric terms
  // from eqn 11 in phil's document

  CH_assert(a_faceMetricTerms.isDefined());

  DataIterator dit = a_faceMetricTerms.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FluxBox& thisFluxBox = a_faceMetricTerms[dit];
      thisFluxBox.setVal(0.0);

      for (int faceDir=0; faceDir<SpaceDim; faceDir++)
        {
          FArrayBox& thisDirData = thisFluxBox[faceDir];

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
                                                              sit());
                            }
                        } // end loop over edge directions
                    } // end loop over high-low
                }
              else
                {
                  // transverse direction...
                  computeTransverseFaceMetric(thisDirData,
                                              faceDir,
                                              d);

                }  // end if transverse direction
            } // end loop over components on this face
        } // end loop over face directions
    } // end loop over boxes
}



void
FourthOrderCoordSys::incrementFaceMetricWithEdgeTerm(FArrayBox& a_faceMetrics,
                                                     int a_faceDir,
                                                     int a_edgeDir,
                                                     const Side::LoHiSide& a_side)
{
  // make sure we have enough components
  CH_assert(a_faceMetrics.nComp() >= SpaceDim);

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
  BoxIterator bit(a_faceMetrics.box());
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
              RealVect realLoc = realCoord(edgeCenter);
              edgeVal = getN(realLoc, sDir, a_faceDir, a_edgeDir);
            }
          else if (SpaceDim == 3)
            {
              // integrate along tangentDir using quadrature
              for (int i=0; i<quadPts.size(); i++)
                {
                  RealVect mappedLoc(edgeCenter);
                  mappedLoc += (quadPts[i].location)*m_dx/2.0;
                  // convert from mapped->real space
                  RealVect realLoc = realCoord(mappedLoc);
                  Real Nvalue = getN(realLoc, sDir, a_faceDir, a_edgeDir);
                  edgeVal += Nvalue*quadPts[i].weight*weightMult;
                }

            }
          else
            {
              MayDay::Error("FourthOrderCoordSys::faceMetrics not defined for SpaceDim > 3");
            }

          // since this is for normal components, faceDir = d
          int component = getMetricTermComponent(sDir,a_faceDir);
          a_faceMetrics(iv, component) += mult*edgeVal;

        } // end loop over s directions
    } // end loop over faces in faceMetrics

}

void
FourthOrderCoordSys::computeTransverseFaceMetric(FArrayBox& a_faceMetrics,
                                                 int a_faceDir,
                                                 int a_dDir)
{
  // first, compute face area
  Real faceArea = m_dx[0];
  for (int dir=1; dir<SpaceDim; dir++)
    {
      faceArea *= m_dx[dir];
    }
  faceArea /= m_dx[a_faceDir];


  // temp storage for nodal values
  Box nodeBox(a_faceMetrics.box());
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
                  Jmatrix[row + SpaceDim*col] = dXdXi(realLoc, row, col);
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
                                  CHF_FRA(nodeN),
                                  CHF_BOX(LapBox),
                                  CHF_INT(dir),
                                  CHF_REAL(factor));
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
              int metricComp = getMetricTermComponent(sDir, a_dDir);
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
      int metricComp = getMetricTermComponent(sDir, a_dDir);
      //
      a_faceMetrics.mult(faceArea, metricComp, 1);
    }
}


void
FourthOrderCoordSys::computeTangentialGrad(LevelData<FluxBox>& a_gradPhi,
                                           const LevelData<FluxBox>& a_phiFace) const
{
  // want to compute undivided differences here
  Real fakeDx = 1.0;
  int nPhiComp = a_phiFace.nComp();
  CH_assert(a_gradPhi.nComp() >= (SpaceDim-1)*nPhiComp);
  DataIterator dit = a_gradPhi.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FluxBox& thisGradPhi = a_gradPhi[dit];
      const FluxBox& thisPhi = a_phiFace[dit];
      for (int faceDir=0; faceDir<SpaceDim; faceDir++)
        {
          FArrayBox& thisGradPhiDir = thisGradPhi[faceDir];
          const FArrayBox& thisPhiDir = thisPhi[faceDir];

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
    } // end loop over boxes

}

#if 0
void
FourthOrderCoordSys::fourthOrderAverage(LevelData<FluxBox>& a_phi) const
{
  DataIterator dit = a_phi.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FluxBox& thisPhi = a_phi[dit];

      for (int dir=0; dir<SpaceDim; dir++)
        {
          FArrayBox& thisPhiDir = thisPhi[dir];
          // derivBox shrunk by 1 in tangential directions
          Box derivBox(thisPhiDir.box());
          derivBox.grow(-1);
          derivBox.grow(dir, 1);
          FArrayBox tempLap(derivBox, thisPhiDir.nComp());
          tempLap.setVal(0.0);

          for (int tanDir=0; tanDir <SpaceDim; tanDir++)
            {
              if (tanDir != dir)
                {
                  Real factor = 1.0/24.0;
                  FORT_INCREMENTLAPLACIAN(CHF_FRA(tempLap),
                                          CHF_FRA(thisPhiDir),
                                          CHF_BOX(derivBox),
                                          CHF_INT(tanDir),
                                          CHF_REAL(factor));

                } // end if tangential direction
            } // end loop over tangential directions
          thisPhiDir.plus(tempLap, derivBox, 0, 0, tempLap.nComp());
        } // end loop over FluxBox face directions
    } // end loop over boxes
}
#endif

int
FourthOrderCoordSys::tanGradComp(const int a_faceDir,
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
FourthOrderCoordSys::computeCellVolumes(LevelData<FArrayBox>& a_cellVolumes)
{
  const DisjointBoxLayout& grids = a_cellVolumes.getBoxes();
  // note that F needs to have two more ghost cells than cellVolumes
  // because we need to take derivatives of F, and because we
  // need to take derivatives in order to compute 4th-order avg of F
  IntVect Fghost = a_cellVolumes.ghostVect() + 2*IntVect::Unit;
  LevelData<FluxBox> F(grids, SpaceDim, Fghost);


  // initialize F to be spatial position
  DataIterator dit = grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FluxBox& thisF = F[dit];
      thisF.setVal(0.0);
      for (int dir=0; dir<SpaceDim; dir++)
        {
          RealVect offset = 0.5*m_dx;
          offset[dir] = 0.0;

          FArrayBox& thisFdir = thisF[dir];
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
    }

  // convert point values to 4th-order face averages
  fourthOrderAverage(F);

#if CH_SPACEDIM > 1
  mappedGridDivergence(a_cellVolumes, F);
#else
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& cellVolFab = a_cellVolumes[dit];
      const FluxBox& thisF = F[dit];
      const FArrayBox& thisFdir = thisF[0];
      BoxIterator bit(cellVolFab.box());
      // this is going to be slow, but we can
      // eventually move this into fortran
      for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect iv = bit();
          IntVect ivNext = iv + BASISV(0);
          cellVolFab(iv, 0) = abs(thisFdir(ivNext, 0) - thisFdir(iv));
        }
    }
#endif

  // now divide by number of components
  for (dit.begin(); dit.ok(); ++dit)
    {
       //       a_cellVolumes[dit] /= SpaceDim;
       a_cellVolumes[dit] /= m_volInterval.size();
    }

}

void
FourthOrderCoordSys::computeJ(LevelData<FArrayBox>& a_J)
{
  // compute cell-centered J, which is just cell-volume divided
  // by volume in mapped space
  const DisjointBoxLayout& grids = a_J.getBoxes();

#if 0
  Real mappedGridVol = D_TERM6(m_dx[0],*m_dx[1],*m_dx[2],*m_dx[3],*m_dx[4],*m_dx[5]);

  DataIterator dit = grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& thisJ = a_J[dit];
      thisJ.copy(m_cellVolumes[dit]);
      // compute J = cellVolume/mappedGridVol
      thisJ.divide(mappedGridVol);
    }
#else

  // note that Jpw needs to have one more ghost cell than J because
  // we need to take derivatives in order to compute 4th-order avg
  IntVect ghost = a_J.ghostVect() + IntVect::Unit;
  LevelData<FArrayBox> Jpw(grids, 1, ghost);

  // initialize Jpw to be cell-centered Jacobian
  DataIterator dit = grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
  {
     FArrayBox& thisJpw = Jpw[dit];
     thisJpw.setVal(0.0);

     BoxIterator bit(thisJpw.box());
     for (bit.begin();bit.ok();++bit)
     {
        IntVect iv = bit();
        RealVect xi;
        for (int idir=0; idir<SpaceDim; idir++)
        {
           xi[idir] = ( 0.5 + iv[idir] ) * m_dx[idir];
        }
        RealVect x( realCoord( xi ) );
        thisJpw(iv) = pointwiseJ( x );
     }
  }

  // convert point values to 4th-order cell averages
  fourthOrderAverage( Jpw );

  // copy Jpw into a_J
  for (dit.begin(); dit.ok(); ++dit)
  {
     Box intersectBox = Jpw[dit].box();
     intersectBox &= a_J[dit].box();
     a_J[dit].copy( Jpw[dit], intersectBox );
  }
#endif
}

void
FourthOrderCoordSys::computeJinverse(LevelData<FluxBox>& a_Jinverse)
{
  // first, compute cell-centered J, which is just cell-volume divided
  // by volume in mapped space
  const DisjointBoxLayout& grids = a_Jinverse.getBoxes();

  // need a data holder with an extra ghost cell in order to
  // be able to compute tangential gradients to compute
  IntVect grownGhost = a_Jinverse.ghostVect() + IntVect::Unit;
  LevelData<FluxBox> tempJinverse(grids, 1, grownGhost);
  fourthOrderCellToFace(tempJinverse, m_J);

  // now compute tanGrad(J)
  int numTanGradComps = SpaceDim-1;
  // kludge to prevent crashing in 1D, not that m_tanGradNTerms will be used.
  if (SpaceDim == 1) numTanGradComps = 1;
  // probably actually want to reduce this by one
  IntVect numTanGradGhost = m_cellVolumes.ghostVect() - IntVect::Unit;
  LevelData<FluxBox> tanGradJ(grids, numTanGradComps, numTanGradGhost);

  computeTangentialGrad(tanGradJ, tempJinverse);

  // now compute <J^(-1)> = 1/<J>[<1> + h^2/24(grad <J>^2)/<J>^2]
  DataIterator dit = grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FluxBox& thisJinverse = a_Jinverse[dit];
      FluxBox& thisTempJinverse = tempJinverse[dit];
      FluxBox& thisTanGradJ = tanGradJ[dit];
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
}

void
FourthOrderCoordSys::computeNJinverse(LevelData<FluxBox>& a_NJinverse)
{
  // do this in pretty much the same way we compute fluxes for divergence

  const DisjointBoxLayout& grids = a_NJinverse.getBoxes();
  CH_assert (grids == m_JInverse.getBoxes());

  // the absolute simplest way to do this (minimizing code duplication) is
  // to copy 1/J to a multi-component temp data holder, and then call
  // computeProductAverage
  LevelData<FluxBox> multiCompJinverse(grids, a_NJinverse.nComp(),
                                       m_JInverse.ghostVect());

  // single component holder for NJinverse
  LevelData<FluxBox> tempStorage(grids, 1, a_NJinverse.ghostVect());

  // need to do this one component at a time
  for (int sComp=0; sComp<SpaceDim; sComp++)
    {
      for (int dComp=0; dComp<SpaceDim; dComp++)
      {
        // now do a fab-by-fab copy of m_JInverse into each
        // component of the temporary
        DataIterator dit = grids.dataIterator();
        for (dit.begin(); dit.ok(); ++dit)
          {
            FluxBox& thisMultiComp = multiCompJinverse[dit];
            const FluxBox& thisJinverse = m_JInverse[dit];

            // want to set only the dir-th component to be Jinverse, and
            // the rest to be zero, in order to isolate this component
            thisMultiComp.setVal(0.0);

            // now copy Jinverse into the comp-th component
            // loop over faceDirs
            for (int dir=0; dir<SpaceDim; dir++)
              {
                const FArrayBox& thisJinverseDir = thisJinverse[dir];
                FArrayBox& thisMultiCompDir = thisMultiComp[dir];

                // copy from 0th comp of Jinverse->comp-th
                // component of multiCompJinverse
                int copyComp = getMetricTermComponent(sComp,dComp);
                thisMultiCompDir.copy(thisJinverseDir, 0, copyComp, 1);
              } // end loop over face directions
          } // end loop over boxes

        // can now do this by using the computeProductAverage function
        computeMetricTermProductAverage(tempStorage, multiCompJinverse);

        // now copy tempStorage into the relevent component in a_NJinverse
        for (dit.begin(); dit.ok(); ++dit)
          {
            FluxBox& thisTemp = tempStorage[dit];
            FluxBox& thisNJinverse = a_NJinverse[dit];
            int comp = getMetricTermComponent(sComp,dComp);
            for (int dir=0; dir<SpaceDim; dir++)
              {
                thisNJinverse[dir].copy(thisTemp[dir], 0, comp, 1);
              }
          }
      } // end loop over d
    } // end loop over components

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
      MayDay::Error("FourthOrderCoordSys::pointwiseJ not defined for SpaceDim > 3");
   }
   return detJ;
}


Real FourthOrderCoordSys::pointwiseJ(const RealVect& a_X) const
{
   // this needs to be improved and generalized
   Vector<Real> Jmtrx(SpaceDim*SpaceDim);
   for (int j=0; j<SpaceDim; j++)
   {
      for (int i=0; i<SpaceDim; i++)
      {
         Jmtrx[i+SpaceDim*j] = dXdXi( a_X, i, j );
      }
   }
   return computePointwiseJ( Jmtrx );
}

#if 0
void FourthOrderCoordSys::pointwiseJ(FArrayBox& a_pointwiseJ,
                                     const FArrayBox& a_X) const
{
   Box intersectBox(a_poitwiseJ.box());
   intersectBox &= a_X.box();

   // get gradX on this box
   FArrayBox gradX( intersectBox, SpaceDim * SpaceDim );
   getJacobianOnBox( gradX );

   // this is going to be slow, but we can
   // eventually move this into fortran
   BoxIterator bit(intersectBox);
   for (bit.begin(); bit.ok(); ++bit)
   {
      IntVect iv = bit();
      Vector<Real> gX(SpaceDim * SpaceDim);
      gradX( &(gX[0]), iv, 0, SpaceDim * SpaceDim );
      computePointwiseJ( a_pointwiseJ(iv), gX );
   }
}
#endif



#include "NamespaceFooter.H"
