#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


#include "mappedLimiter.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "LevelData.H"
#include "FourthOrderCoordSys.H"
#include "mappedLimiterF_F.H"

//#include <sstream>
//#include <fstream>

//#include <complex>
#include <algorithm>


#include "NamespaceHeader.H"

using std::max;

void
applyMappedLimiter(LevelData<FluxBox>& a_facePhi,
                   const LevelData<FArrayBox>& a_cellPhi,
                   const LevelData<FluxBox>& a_faceVel,
                   const FourthOrderCoordSys* a_FOCS,
                   const RealVect& a_dx,
                   const bool velocities_are_physical)
{
  int nComp = a_facePhi.nComp();

  // this specifies the number of ghost faces we will need in
  // each transverse and normal direction
  int transverseGrow = 2;
  int normalGrow = 0;

  // use C value from Colella and Sekora
  Real limiterConst =  1.25;

  // may need to do exchange on facePhi
  a_facePhi.exchange();

  const DisjointBoxLayout& grids = a_facePhi.getBoxes();
  DataIterator dit = grids.dataIterator();

  // in order to do upwinding, need normal velocities in
  // computational space. computing Fourth-order face
  // averages just in order to upwinding is over-working the
  // issue, but it's the most convenient thing to do here.
  LevelData<FluxBox> normalVel(grids, 1, a_faceVel.ghostVect());

  if (!velocities_are_physical)
     {
        a_FOCS->computeMetricTermProductAverage(normalVel, a_faceVel);
     }
  else
     {
        const RealVect& dx = a_FOCS->dx();

        RealVect face_area;
        for (int dir=0; dir<SpaceDim; ++dir)
        {
           face_area[dir] = 1.;
           for (int tdir=0; tdir<SpaceDim; ++tdir)
           {
              if (tdir != dir) face_area[dir] *= dx[tdir];
           }
        }

        for (dit.begin(); dit.ok(); ++dit)
        {
           const FluxBox& thisFaceVel = a_faceVel[dit];
           FluxBox& thisNormalVel = normalVel[dit];
           for (int dir=0; dir<SpaceDim; ++dir)
           {
              FArrayBox& thisNormalVel_dir = thisNormalVel[dir];
              thisNormalVel_dir.copy(thisFaceVel[dir],dir,0,1);
              thisNormalVel_dir *= face_area[dir];
           }
        }
     }

  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& gridBox = grids[dit];
      Box lapBox(gridBox);
      // need two laplacian ghost cells for this limiter
      int lapBoxGrow = max(normalGrow+2, transverseGrow);
      lapBox.grow(lapBoxGrow);
      FArrayBox ccLaplacians(lapBox, nComp);
      FluxBox& thisFacePhi = a_facePhi[dit];
      const FArrayBox& thisCellPhi = a_cellPhi[dit];
      const FluxBox& thisNormalVel = normalVel[dit];
      // check to be sure that we have enough phi to compute
      // all of these Laplacians
      {
        Box LapPhiBox(lapBox);
        LapPhiBox.grow(1);
        CH_assert(thisCellPhi.box().contains(LapPhiBox));
      }

      for (int dir=0; dir<SpaceDim; dir++)
        {
          // box of valid edges for this grid
          Box faceBox(gridBox);
          // for 4th order, we need extra faces in the transverse
          // directions. (for mapped-grid case, we actually need
          // _2_ transverse faces

          // need extra faces in the tangential directions in order
          // to handle 4th-order multiplies
          faceBox.grow(transverseGrow);
          // however, need different number of extra faces in normal dir
          faceBox.grow(dir,normalGrow - transverseGrow);

          faceBox.surroundingNodes(dir);
          FArrayBox& thisFacePhiDir = thisFacePhi[dir];
          const FArrayBox& thisNormalVelDir = thisNormalVel[dir];
          {
            // context for computing limited face values
            Box grownFaceBox(faceBox);
            // need an extra face's worth of the FC laplacians
            grownFaceBox.grow(dir,1);

            FArrayBox centeredLaplacian(grownFaceBox, nComp);

            // compute centered Laplacian
            centeredLaplacian.setVal(0.0);
            FORT_CENTEREDLAPLACIAN(CHF_FRA(centeredLaplacian),
                                   CHF_CONST_FRA(thisFacePhiDir),
                                   CHF_CONST_FRA(thisCellPhi),
                                   CHF_BOX(grownFaceBox),
                                   CHF_CONST_INT(dir),
                                   CHF_CONST_REAL(a_dx[dir]));

            // compute cell-centered Laplacians
            ccLaplacians.setVal(0.0);
            FORT_CCLAPLACIAN(CHF_FRA(ccLaplacians),
                             CHF_CONST_FRA(thisCellPhi),
                             CHF_BOX(lapBox),
                             CHF_CONST_INT(dir),
                             CHF_CONST_REAL(a_dx[dir]));


            // now compute limited face value
            FORT_LIMITFACEVALUES(CHF_FRA(thisFacePhiDir),
                                 CHF_CONST_FRA(thisCellPhi),
                                 CHF_CONST_FRA(centeredLaplacian),
                                 CHF_CONST_FRA(ccLaplacians),
                                 CHF_BOX(grownFaceBox),
                                 CHF_CONST_INT(dir),
                                 CHF_CONST_REAL(a_dx[dir]),
                                 CHF_CONST_REAL(limiterConst));
            // end context for computing limited face values
          } // (centeredLaplacian goes out of scope)

          // now compute parabolic interpolants
          // storage for cell-centered D^2a
          Box lapBoxDir(gridBox);
          // need transverse faces for 4th order
          lapBoxDir.grow(transverseGrow);
          // need this to be grown by one in normal dir
          lapBoxDir.grow(dir,normalGrow-transverseGrow+1);

          FArrayBox D2a(lapBoxDir,nComp);
          FArrayBox D2aLim(lapBoxDir,nComp);

          // initialize D2a to be -a_6/(h^2)
          // first compute a_6...
          FORT_COMPUTEA6(CHF_FRA(D2a),
                         CHF_CONST_FRA(thisCellPhi),
                         CHF_CONST_FRA(thisFacePhiDir),
                         CHF_BOX(lapBoxDir),
                         CHF_CONST_INT(dir));

          // then multiply by -2/(h^2)
          // (dfm 5/5/09 -- was missing the factor of 2)
          Real mult = -2.0/(a_dx[dir]*a_dx[dir]);
          D2a *= mult;

          // now limit D2a w/r/t other 2nd-derivatives
          FORT_LIMITD2A(CHF_FRA(D2aLim),
                        CHF_CONST_FRA(D2a),
                        CHF_CONST_FRA(ccLaplacians),
                        CHF_BOX(lapBoxDir),
                        CHF_CONST_INT(dir),
                        CHF_CONST_REAL(limiterConst));


          // storage for left and right interpolants
          //   note that while leftPhi and rightPhi are associated
          //   with face indicies, it is more natural to compute the
          //   PPM limited face values on a per-cell basis.  Thus,
          //   each of these FArrayBoxes is one cell longer than it
          //   need to be
          Box growBox( faceBox );
          growBox.grow(dir,1);
          FArrayBox leftPhi(growBox, nComp);
          FArrayBox rightPhi(growBox, nComp);

          // We operate on the cells of the domain plus one ghost on each
          // end in the current direction
          Box cellBox( growBox );
          cellBox.growHi(dir,-1);
          FORT_LEFTRIGHTSTATES(CHF_FRA(leftPhi),
                               CHF_FRA(rightPhi),
                               CHF_CONST_FRA(thisFacePhiDir),
                               CHF_CONST_FRA(thisCellPhi),
                               CHF_CONST_FRA(D2aLim),
                               CHF_CONST_FRA(D2a),
                               CHF_BOX(cellBox),
                               CHF_CONST_INT(dir));

          // need to do something about boundary values here

          // now pick upwind state
          // note that normalVel only has a single component, so we
          // use 0 as the component argument for the CHF_FRA1 macro
          // as a debugging check, setThisFacePhiDir to a bogus value
          // first
          thisFacePhiDir.setVal(1000000000);
          FORT_SELECTUPWIND(CHF_FRA(thisFacePhiDir),
                            CHF_CONST_FRA(leftPhi),
                            CHF_CONST_FRA(rightPhi),
                            CHF_CONST_FRA1(thisNormalVelDir,0),
                            CHF_BOX(faceBox));

        } // end loop over directions

    } // end loop over grids
}

#include "NamespaceFooter.H"
