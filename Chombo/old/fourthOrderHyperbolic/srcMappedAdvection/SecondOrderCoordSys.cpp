#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "SecondOrderCoordSys.H"
#include "SecondOrderCoordSysF_F.H"
#include "DivergenceF_F.H"
#include "BoxIterator.H"
#include "GaussianQuadrature.H"
#include "NewtonCotesQuadrature.H"
#include "PointwiseDotProdF_F.H"

#include "NamespaceHeader.H"

void secondOrderCellToFace(LevelData<FluxBox>& a_faceData,
                           LevelData<FArrayBox>& a_cellData)
{
  const DisjointBoxLayout& grids = a_faceData.getBoxes();

  DataIterator dit = grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& thisCellData = a_cellData[dit];
      FluxBox& thisFaceData = a_faceData[dit];
      for (int dir=0; dir<SpaceDim; dir++)
        {
          FArrayBox& thisFaceDataDir = thisFaceData[dir];
          // compute box over which we can do the averaging
          Box testBox(thisCellData.box());
          // need two cells in each direction for averaging
          testBox.grow(dir, -2);
          testBox.surroundingNodes(dir);
          testBox &= thisFaceDataDir.box();
          FORT_CELLTOFACE2NDORDER(CHF_FRA(thisFaceDataDir),
                                  CHF_FRA(thisCellData),
                                  CHF_BOX(testBox),
                                  CHF_INT(dir));
        }
    } // end loop over grid boxes
}

SecondOrderCoordSys::SecondOrderCoordSys()
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

SecondOrderCoordSys::~SecondOrderCoordSys()
{
  if (m_quadraturePtr != NULL)
    {
      delete m_quadraturePtr;
      m_quadraturePtr = NULL;
    }
}

void
SecondOrderCoordSys::define(const DisjointBoxLayout& a_grids,
                            const ProblemDomain& a_domain,
                            const RealVect& a_cellSpacing,
                            const IntVect& a_ghostVect)
{
  m_dx = a_cellSpacing;
  m_domain = a_domain;

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

  int nFaceMetricTerms = SpaceDim;
  m_faceMetricTerms.define(a_grids, nFaceMetricTerms, grownGhost3);

  // need cell volumes defined with at least two ghost cells in order to compute
  // second-order face-centered Jinverse
  m_cellVolumes.define(a_grids, 1, grownGhost2);

  // m_Jinverse needs an extra ghost cell in order to compute gradients for NinverseJ
  m_JInverse.define(a_grids, 1, grownGhost1);

  m_NinverseJ.define(a_grids, SpaceDim, a_ghostVect);

  // this was originally in the derived class's define function.
  defineMetricTerms();

}

void
SecondOrderCoordSys::mappedGridDivergence(LevelData<FArrayBox>& a_divF,
                                          const LevelData<FluxBox>& a_F)
{
  int nFluxComp = a_F.nComp()/SpaceDim;
  LevelData<FluxBox> tempFlux(a_F.getBoxes(), nFluxComp,
                              a_F.ghostVect());

  computeFaceFlux(tempFlux, a_F);

  // now compute divergence in the usual way
  DataIterator dit = a_divF.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& thisDiv = a_divF[dit];
      FluxBox& thisFlux = tempFlux[dit];
      // first, set divF to 0
      thisDiv.setVal(0.0);
      // since we're computing volume*div, and the fluxes
      // are multiplied by the appropriate areas, don't need
      // a dx here.
      Real fakeDx = 1.0;

      // now loop over directions and increment with directional
      // derivative
      for (int dir=0; dir<SpaceDim; dir++)
        {
          const Box& thisBox = thisDiv.box();
          // use fortran from divergence operator here
          FORT_DIVERGENCE(CHF_CONST_FRA(thisFlux[dir]),
                          CHF_FRA(thisDiv),
                          CHF_BOX(thisBox),
                          CHF_CONST_REAL(fakeDx),
                          CHF_INT(dir));

        }
    }

}

void
SecondOrderCoordSys::setQuadrature(const EdgeQuadrature* a_quadrature)
{
  if (m_quadraturePtr != NULL)
    {
      delete m_quadraturePtr;
      m_quadraturePtr = NULL;
    }
  m_quadraturePtr = a_quadrature->new_quadrature();
}

const LevelData<FluxBox>&
SecondOrderCoordSys::getFaceMetricTerms() const
{
  return m_faceMetricTerms;
}

const
LevelData<FArrayBox>&
SecondOrderCoordSys::getCellVolumes() const
{
  return m_cellVolumes;
}

const
LevelData<FluxBox>&
SecondOrderCoordSys::getJInverse() const
{
  return m_JInverse;
}

const
LevelData<FluxBox>&
SecondOrderCoordSys::getNJinverse() const
{
  return m_NinverseJ;
}

Real
SecondOrderCoordSys::getN(const RealVect& a_X,
                          int a_s, int a_d, int a_d1) const
{
  // allocate temporaries
  Real polarity = 1.0;
  // this keeps track of rows before/after swapping
  int rowIndex[SpaceDim];
  D_TERM(rowIndex[0] = 0;,
         rowIndex[1] = 1;,
         rowIndex[2] = 2;);

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
  D_TERM(
         A[0][0] = getNMatrixEntry(a_X, a_s, a_d, a_d1,
                                   rowIndex[0], 0);,
         A[1][1] = getNMatrixEntry(a_X, a_s, a_d, a_d1,
                                   rowIndex[1], 1);,
         A[2][2] = getNMatrixEntry(a_X, a_s, a_d, a_d1,
                                   rowIndex[2], 2) );

  // if we're in 2D, then no need to do any more elimination
  if (SpaceDim > 2)
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
          MayDay::Error("SecondOrderCoordSys::getN -- Bad value for a_d");
        }

    } // end if SpaceDim > 2

  // note that N blows up for Dim = 1
  // (come back and fix this later, if required)
  CH_assert(SpaceDim > 1);
  Real N = polarity*D_TERM(A[0][0], * A[1][1], *A[2][2])/(SpaceDim - 1);

  return N;
}

Real
SecondOrderCoordSys::getNMatrixEntry(const RealVect& a_X,
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
SecondOrderCoordSys::computeFaceFlux(LevelData<FluxBox>& a_Flux,
                                     const LevelData<FluxBox>& a_F)
{

  //int nFlux = a_Flux.nComp();
  const DisjointBoxLayout& grids = a_F.getBoxes();

  const LevelData<FluxBox>& faceMetricTerms = getFaceMetricTerms();

  DataIterator dit = grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const FluxBox& thisF = a_F[dit];
      FluxBox& thisFlux = a_Flux[dit];

      // precomputed face-centered metric terms
      const FluxBox& thisFaceMetric = faceMetricTerms[dit];

      for (int dir=0; dir<SpaceDim; dir++)
        {
          FArrayBox& thisFluxDir = thisFlux[dir];
          const FArrayBox& thisFDir = thisF[dir];
          const FArrayBox& thisFaceMetricDir = thisFaceMetric[dir];

          Box intersectBox(thisFluxDir.box());
          intersectBox &= thisFaceMetricDir.box();

          int fstart = 0;
          int metricStart = 0;
          int ncomp = thisFDir.nComp();

          // first part is easy -- dot product of F and metrics
          FORT_POINTDOTPROD(CHF_FRA1(thisFluxDir, 0),
                            CHF_CONST_FRA(thisFDir),
                            CHF_CONST_INT(fstart),
                            CHF_CONST_FRA(thisFaceMetricDir),
                            CHF_CONST_INT(metricStart),
                            CHF_CONST_INT(ncomp),
                            CHF_BOX(intersectBox) );

        } // end loop over face directions
    } // end loop over grids
}



void
SecondOrderCoordSys::defineFaceMetricTerms(LevelData<FluxBox>& a_faceMetricTerms)
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
        } // end loop over face directions
    } // end loop over boxes
}

void
SecondOrderCoordSys::defineMetricTerms()
{
  // this is where we compute the face-centered metric terms
  // from eqn 11 in phil's document

  CH_assert(m_faceMetricTerms.isDefined());

  defineFaceMetricTerms(m_faceMetricTerms);

  computeCellVolumes(m_cellVolumes);

  computeJinverse(m_JInverse);

  computeNJinverse(m_NinverseJ);
}


void
SecondOrderCoordSys::incrementFaceMetricWithEdgeTerm(FArrayBox& a_faceMetrics,
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
              MayDay::Error("SecondOrderCoordSys::faceMetrics not defined for SpaceDim > 3");
            }

          a_faceMetrics(iv, sDir) += mult*edgeVal;

        } // end loop over s directions
    } // end loop over faces in faceMetrics

}


void
SecondOrderCoordSys::computeCellVolumes(LevelData<FArrayBox>& a_cellVolumes)
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
              D_TERM(thisFdir(iv,0) = realLoc[0];,
                     thisFdir(iv,1) = realLoc[1];,
                     thisFdir(iv,2) = realLoc[2];)
#endif
            }
        } // end loop over directions
    }

  mappedGridDivergence(a_cellVolumes, F);

  // now divide by number of components
  for (dit.begin(); dit.ok(); ++dit)
    {
      //      a_cellVolumes[dit] /= SpaceDim;
      a_cellVolumes[dit] /= m_volInterval.size();
    }

}

void
SecondOrderCoordSys::computeJinverse(LevelData<FluxBox>& a_Jinverse)
{
  // first, compute cell-centered J, which is just cell-volume divided
  // by volume in mapped space
  Real mappedGridVol = D_TERM(m_dx[0],*m_dx[1],*m_dx[2]);

  const DisjointBoxLayout& grids = a_Jinverse.getBoxes();

  LevelData<FArrayBox> tempJ(grids, 1, m_cellVolumes.ghostVect());
  // need a data holder with an extra ghost cell in order to
  // be able to compute tangential gradients to compute
  IntVect grownGhost = a_Jinverse.ghostVect() + IntVect::Unit;
  LevelData<FluxBox> tempJinverse(grids, 1, grownGhost);

  DataIterator dit = grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& thisTempJ = tempJ[dit];
      thisTempJ.copy(m_cellVolumes[dit]);
      // compute J = cellVolume/mappedGridVol
      thisTempJ.divide(mappedGridVol);
    }

  secondOrderCellToFace(tempJinverse, tempJ);

  // now compute <J^(-1)> = 1/<J>[<1> + h^2/24(grad <J>^2)/<J>^2]
  for (dit.begin(); dit.ok(); ++dit)
    {
      FluxBox& thisJinverse = a_Jinverse[dit];
      FluxBox& thisTempJinverse = tempJinverse[dit];
      for (int dir=0; dir<SpaceDim; dir++)
        {
          // compute box over which we can compute this
          // note that we don't need dx because the h^2 cancels with the
          // denominator in the (undivided) gradient
          thisJinverse[dir].copy(thisTempJinverse[dir], thisJinverse[dir].box());
        }
    }
}

void
SecondOrderCoordSys::computeNJinverse(LevelData<FluxBox>& a_NJinverse)
{
  // do this in pretty much the same way we compute fluxes for divergence

  const DisjointBoxLayout& grids = a_NJinverse.getBoxes();
  CH_assert (grids == m_JInverse.getBoxes());

  // the absolute simplest way to do this (minimizing code duplication) is
  // to copy 1/J to a multi-component temp data holder, and then call
  // computeFaceFlux
  LevelData<FluxBox> multiCompJinverse(grids, a_NJinverse.nComp(),
                                       m_JInverse.ghostVect());

  // single component holder for NJinverse
  LevelData<FluxBox> tempStorage(grids, 1, a_NJinverse.ghostVect());

  // need to do this one component at a time
  for (int comp=0; comp<a_NJinverse.nComp(); comp++)
    {
      // now do a fab-by-fab copy of m_JInverse into each component of the temporary
      DataIterator dit = grids.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          FluxBox& thisMultiComp = multiCompJinverse[dit];
          const FluxBox& thisJinverse = m_JInverse[dit];

          // want to set only the dir-th component to be Jinverse, and
          // the rest to be zero, in order to isolate this component
          thisMultiComp.setVal(0.0);

          // now copy Jinverse into the comp-th component
          for (int dir=0; dir<SpaceDim; dir++)
            {
              const FArrayBox& thisJinverseDir = thisJinverse[dir];
              FArrayBox& thisMultiCompDir = thisMultiComp[dir];

              // copy from 0th comp of Jinverse->comp-th
              // component of multiCompJinverse
              thisMultiCompDir.copy(thisJinverseDir, 0, comp, 1);
            } // end loop over face directions
        } // end loop over boxes

      // can now do this by using the computeFaceFlux function
      computeFaceFlux(tempStorage, multiCompJinverse);

      // now copy tempStorage into the relevent component in a_NJinverse
      for (dit.begin(); dit.ok(); ++dit)
        {
          FluxBox& thisTemp = tempStorage[dit];
          FluxBox& thisNJinverse = a_NJinverse[dit];
          for (int dir=0; dir<SpaceDim; dir++)
            {
              thisNJinverse[dir].copy(thisTemp[dir], 0, comp, 1);
            }
        }
    } // end loop over components

}

#include "NamespaceFooter.H"
