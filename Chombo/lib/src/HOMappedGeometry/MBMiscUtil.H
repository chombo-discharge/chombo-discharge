#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MBMISCUTIL_H_
#define _MBMISCUTIL_H_

#include "FArrayBox.H"
#include "LevelData.H"
#include "RealVect.H"

#include "NamespaceHeader.H"

/// return a_inputBox intersected with a version of a_domainBox extended indefinitely below if offset is -1 and above if offset is +1
Box intersectOffset(const Box& a_inputBox,
                    const Box& a_domainBox,
                    const IntVect& a_offset);

///
/**
   Returns dx * gradient of a_dataFab, using a second-order formula.
   Takes data in a_dataFab available on a_dataBox
   and returns dx * gradient in a_gradFab on a_gradBox.

   Order of components within a_gradFab, with 2D vector data:
   difference in data[0] in direction 0;
   difference in data[1] in direction 0;
   difference in data[0] in direction 1;
   difference in data[1] in direction 1.

   That is:  dData[0]/dx[0], dData[1]/dx[0], dData[0]/dx[1], dData[1]/dx[1].
 */
void order2gradient(FArrayBox&          a_gradFab,
                    const Box&          a_gradBox,
                    const FArrayBox&    a_dataFab,
                    const Box&          a_dataBox,
                    const Vector<int>&  a_interpDimsVect);

/// returns a_gradData = gradient(a_data) on wherever a_data is defined minus 1 layer
void order2gradient(LevelData<FArrayBox>&        a_gradData,
                    const LevelData<FArrayBox>&  a_data,
                    Vector<int>                  a_interpDimsVect);

/// returns the dimension of the faces of a_box; crashes if not on faces of codimension 1
int faceDimension(const Box& a_box);

/// returns the faces of a_gridBox that are on the (a_idir, a_side) face of a_blockBox
Box faceSlice(const Box& a_gridBox,
              int a_idir,
              Side::LoHiSide a_side,
              const Box& a_blockBox);

/// component d of return value is a_dxVect[d]/2 if a_box is cell-centered in dimension d, and 0 if node-centered in dimension d.  For a d-face-centered box, returns 0 in dimension d, and a_dxVect/2 in other dimensions.
RealVect centeringOffset(const Box& a_box,
                         const RealVect& a_dxVect);

/// convert Vector<RealVect> to Vector<Real>
Vector<Real> VectorRealVectToReal(const Vector<RealVect>& a_rv);

/// convert Vector<Real> to Vector<RealVect>
Vector<RealVect> VectorRealToRealVect(const Vector<Real>& a_r);

/// change a_bx to its minBox across all processors
void minBoxGlobal(Box& a_bx);

int minBufferSizeAMR(const Vector<int>& a_refRatio,
                     int                a_maxLevel,
                     int                a_numGhost,
                     int                a_order = 4,
                     bool               a_fromJU = true);

/// returns ceil(a_num / a_denom)
int ceilRatio(int a_num,
              int a_denom);

/// returns ceil(a_num / a_denom) component by component
IntVect ceilRatio(const IntVect& a_num,
                  const IntVect& a_denom);

Real dotSubvectors(const Vector<Real>& a_vec1,
                   int a_start1,
                   const Vector<Real>& a_vec2,
                   int a_start2,
                   int a_length);

void addPowersPoint(Vector<Real>& a_powers,
                    const RealVect& a_x,
                    const Vector<IntVect>& a_exponents,
                    Real a_weight);

void addPowersVector(Vector<Real>& a_powers,
                     const RealVect& a_x,
                     const Vector<IntVect>& a_exponents,
                     Real a_weight,
                     const Real* const a_vecTransform);

#include "NamespaceFooter.H"

#endif
