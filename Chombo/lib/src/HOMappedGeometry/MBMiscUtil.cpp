#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "MBMiscUtil.H"
#include "SecondOrderGradientF_F.H"
#include "SPMD.H"

#include "NamespaceHeader.H"

Box intersectOffset(const Box& a_inputBox,
                    const Box& a_domainBox,
                    const IntVect& a_offset)
{
  // added by petermc, 29 Sep 2011
  if ( a_domainBox.isEmpty() ) return Box();
  Box returnBox(a_inputBox);
  // Intersect returnBox with a_domainBox.
  // But a_domainBox has no outer limit, only an inner limit;
  // that is why we can't just do returnBox &= a_domainBox.
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      int biggest = a_domainBox.bigEnd(idir); // useful iff offset -1 or 0
      int smallest = a_domainBox.smallEnd(idir); // useful iff offset 0 or +1
      switch (a_offset[idir])
        {
        case -1:
          if (returnBox.smallEnd(idir) <= biggest)
            { // returnBox within a_domainBox; remove the part outside
              returnBox.setBig(idir, biggest);
            }
          else
            { // returnBox does NOT intersect a_domainBox, so discard
              return Box();
            }
          break;
        case +1:
          if (returnBox.bigEnd(idir) >= smallest)
            { // returnBox within a_domainBox; remove the part outside
              returnBox.setSmall(idir, smallest);
            }
          else
            { // returnBox does NOT intersect a_domainBox, so discard
              return Box();
            }
          break;
        case 0:
          if (returnBox.smallEnd(idir) < smallest)
            returnBox.setSmall(idir, smallest);
          if (returnBox.bigEnd(idir) > biggest)
            returnBox.setBig(idir, biggest);
          break;
        default:
          MayDay::Error("intersectOffset() must have offset in [-1:1]^D");
        }
    }
  return returnBox;
}


void order2gradient(FArrayBox&          a_gradFab,
                    const Box&          a_gradBox,
                    const FArrayBox&    a_dataFab,
                    const Box&          a_dataBox,
                    const Vector<int>&  a_interpDimsVect)
{
  int ncomp = a_dataFab.nComp();
  int indGrad = 0;
  for (int ind = 0; ind < a_interpDimsVect.size(); ind++)
    {
      int idir = a_interpDimsVect[ind];
      // idir specifies which derivative of data[comp]
      Box loBox, hiBox;
      int hasLo, hasHi;
      Box centerBox = a_gradBox;

      int smallData = a_dataBox.smallEnd(idir);
      int smallGrad = a_gradBox.smallEnd(idir);
      if (smallData > smallGrad)
        {
          MayDay::Error("order2gradient: need gradient on bigger box than data");
        }
      else if (smallData < smallGrad)
        {
          // Got enough data to use centered difference on low side.
          hasLo = 0;
        }
      else
        {
          // smallData == smallGrad;
          hasLo = 1;
          // Set loBox to sliver along low side of a_gradBox.
          loBox = a_gradBox;
          loBox.setBig(idir, smallData);
          centerBox.setSmall(idir, smallData+1);
        }

      int bigData = a_dataBox.bigEnd(idir);
      int bigGrad = a_gradBox.bigEnd(idir);
      if (bigData < bigGrad)
        {
          MayDay::Error("order2gradient: need gradient on bigger box than data");
        }
      else if (bigData > bigGrad)
        {
          // Got enough data to use centered difference on high side.
          hasHi = 0;
        }
      else
        {
          // bigData == bigGrad;
          hasHi = 1;
          // Set hiBox to sliver along high side of a_gradBox.
          hiBox = a_gradBox;
          hiBox.setSmall(idir, bigData);
          centerBox.setBig(idir, bigData-1);
        }

      for (int comp = 0; comp < ncomp; comp++)
        { // comp specifies which component of Data
          FORT_SECONDORDERGRAD(CHF_FRA1(a_gradFab, indGrad),
                               CHF_CONST_FRA1(a_dataFab, comp),
                               CHF_CONST_INT(idir),
                               CHF_BOX(loBox),
                               CHF_CONST_INT(hasLo),
                               CHF_BOX(hiBox),
                               CHF_CONST_INT(hasHi),
                               CHF_BOX(centerBox));
          // calcVectorFunGradFab.divide(dxDir, indGrad); // component indGrad only
          indGrad++;
        }
    }
}


void order2gradient(LevelData<FArrayBox>&        a_gradData,
                    const LevelData<FArrayBox>&  a_data,
                    Vector<int>                  a_interpDimsVect)
{
  DataIterator dit = a_data.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const FArrayBox& dataFab = a_data[dit];
      const Box& dataBox = dataFab.box();
      FArrayBox& gradFab = a_gradData[dit];
      const Box& gradBox = gradFab.box();
      order2gradient(gradFab, gradBox,
                     dataFab, dataBox,
                     a_interpDimsVect);
    }
  // We're finding gradient of a_data on all cells (not just ghosts),
  // so we don't lose anything in exchange().
  a_gradData.exchange();
}


int faceDimension(const Box& a_box)
{
  const IntVect& tp = a_box.type();
  int idirFace;
  bool found = false;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (tp[idir] == 1)
        {
          if (found)
            { // already found
              MayDay::Error("faceDimension takes a Box with faces of codimension 1 only");
            }
          else
            {
              found = true;
              idirFace = idir;
            }
        }
    }
  CH_assert(found);
  return idirFace;
}


Box faceSlice(const Box& a_gridBox,
              int a_idir,
              Side::LoHiSide a_side,
              const Box& a_blockBox)
{
  Box blockBoxFaces = surroundingNodes(a_blockBox, a_idir);
  int blockExtremity = blockBoxFaces.sideEnd(a_side)[a_idir];

  Box gridBoxFaces = surroundingNodes(a_gridBox, a_idir);
  int gridBoxLo = gridBoxFaces.smallEnd(a_idir);
  int gridBoxHi = gridBoxFaces.bigEnd(a_idir);
  Box returnBoxFaces;
  if ((gridBoxLo <= blockExtremity) && (blockExtremity <= gridBoxHi))
    { // Otherwise, return empty box.
      // Note that returnBoxFaces can be bigger than face of blockBoxFaces, so
      // we don't just return intersection of gridBoxFaces with blockBoxFaces.
      returnBoxFaces = gridBoxFaces;
      returnBoxFaces.setRange(a_idir, blockExtremity);
    }
  return returnBoxFaces;
}


RealVect centeringOffset(const Box& a_box,
                         const RealVect& a_dxVect)
{
  IntVect tp = a_box.type(); // 0 for cell centers, 1 for nodes
  IntVect offHalf = IntVect::Unit - tp; // 1 for cell centers, 0 for nodes
  RealVect offset = 0.5 * a_dxVect * RealVect(offHalf);
  // offset == dx/2 for cell centers, 0 for nodes.
  return offset;
}


Vector<Real> VectorRealVectToReal(const Vector<RealVect>& a_rv)
{
  int len = a_rv.size();
  Vector<Real> returnVec(len * SpaceDim);
  int arri = 0;
  for (int i = 0; i < len; i++)
    {
      const RealVect& rvThis = a_rv[i];
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          returnVec[arri + idir] = rvThis[idir];
        }
      arri += SpaceDim;
    }
  return returnVec;
}


Vector<RealVect> VectorRealToRealVect(const Vector<Real>& a_r)
{
  int len = a_r.size() / SpaceDim;
  Vector<RealVect> returnVec(len);
  int arri = 0;
  for (int i = 0; i < len; i++)
    {
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          returnVec[i][idir] = a_r[arri + idir];
        }
      arri += SpaceDim;
    }
  return returnVec;
}


void minBoxGlobal(Box& a_bx)
{
#ifdef CH_MPI
  int iprocdest = 0;
  Vector<Box> allBoxes(numProc());
  gather(allBoxes, a_bx, iprocdest);

  // Processor iprocdest now knows Vector<Box> allBoxes.
  // Now broadcast it.
  for (int iproc = 0; iproc < numProc(); iproc++)
    {
      // void broadcast(T& a_inAndOut, int a_src) requires for T:
      // linearSize<T>, linearIn<T>, linearOut<T>.
      broadcast(allBoxes[iproc], iprocdest);
    }

  // Now every processor has Vector<Box> allBoxes.
  Box returnBox;
  for (int ibox = 0; ibox < allBoxes.size(); ibox++)
    {
      returnBox.minBox(allBoxes[ibox]);
    }

  a_bx = returnBox;
#endif
}


int minBufferSizeAMR(const Vector<int>& a_refRatio,
                     int                a_maxLevel,
                     int                a_numGhost,
                     int                a_order,
                     bool               a_fromJU)
{
  if (a_maxLevel == 0)
    {
      if (a_fromJU)
        return 1;
      else
        return 0;
    }

  // Find the minimum in a_refRatio.
  int arraySize = a_refRatio.size();
  CH_assert(a_maxLevel <= arraySize);
  CH_assert(arraySize > 0);
  int minRefRatio = a_refRatio[0];
  for (int i = 1; i < a_maxLevel; i++)
    {
      if (a_refRatio[i] < minRefRatio) minRefRatio = a_refRatio[i];
    }
  CH_assert(minRefRatio > 0);

  // Compute the buffer size.
  int buffer =
    ceilRatio(a_numGhost, minRefRatio) // coarse cells under fine cells
    + (a_order + 1)/2; // for stencil
  if (a_fromJU) buffer++; // for computing <U> from <JU>
  return buffer;
}


int ceilRatio(int a_num,
              int a_denom)
{
  int retval = a_num / a_denom;
  if (retval * a_denom < a_num) retval++;
  return retval;
}


/// returns ceil(a_num / a_denom) component by component
IntVect ceilRatio(const IntVect& a_num,
                  const IntVect& a_denom)
{
  IntVect retVect;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      retVect[idir] = ceilRatio(a_num[idir], a_denom[idir]);
    }
  return retVect;
}


Real dotSubvectors(const Vector<Real>& a_vec1,
                   int a_start1,
                   const Vector<Real>& a_vec2,
                   int a_start2,
                   int a_length)
{
  Real sum = 0.;
  for (int i1 = a_start1, i2 = a_start2, ilen = 0;
       ilen < a_length;
       i1++, i2++, ilen++)
    {
      sum += a_vec1[i1] * a_vec2[i2];
    }
  return sum;
}


void addPowersPoint(Vector<Real>& a_powers,
                    const RealVect& a_x,
                    const Vector<IntVect>& a_exponents,
                    Real a_weight)
{
  int numPowers = a_exponents.size();
  CH_assert(a_powers.size() == numPowers);

  IntVect exponentsPrev = IntVect::Zero;
  RealVect curPower = RealVect::Unit;
  for (int ipow = 0; ipow < numPowers; ipow++)
    {
      IntVect exponentsThis = a_exponents[ipow];
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          if (exponentsThis[idir] < exponentsPrev[idir])
            {
              // exponent in dimension idir decreases,
              // so reset current power in dimension idir to 1.
              curPower[idir] = 1.;
            }
          else if (exponentsThis[idir] == exponentsPrev[idir] + 1)
            {
              // exponent in dimension idir increases by 1,
              // so multiply current power in dimension idir by a_x[idir].
              curPower[idir] *= a_x[idir];
              // and make no changes to current powers in other dimensions
              break;
            }
        }
      // save exponents for next iteration
      exponentsPrev = exponentsThis;
      a_powers[ipow] += a_weight * curPower.product();
    }
}


void addPowersVector(Vector<Real>& a_powers,
                     const RealVect& a_x,
                     const Vector<IntVect>& a_exponents,
                     Real a_weight,
                     const Real* const a_vecTransform)
{
  int numPowers = a_exponents.size();
  CH_assert(a_powers.size() == (SpaceDim*SpaceDim) * numPowers);

  IntVect exponentsPrev = IntVect::Zero;
  RealVect curPower = RealVect::Unit;
  for (int ipow = 0; ipow < numPowers; ipow++)
    {
      IntVect exponentsThis = a_exponents[ipow];
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          if (exponentsThis[idir] < exponentsPrev[idir])
            {
              // exponent in dimension idir decreases,
              // so reset current power in dimension idir to 1.
              curPower[idir] = 1.;
            }
          else if (exponentsThis[idir] == exponentsPrev[idir] + 1)
            {
              // exponent in dimension idir increases by 1,
              // so multiply current power in dimension idir by a_x[idir].
              curPower[idir] *= a_x[idir];
              // and make no changes to current powers in other dimensions
              break;
            }
        }
      // save exponents for next iteration
      exponentsPrev = exponentsThis;
      Real weightedCurrentProd = a_weight * curPower.product();
      for (int ivec = 0, ind = ipow;
           ivec < SpaceDim*SpaceDim;
           ivec++, ind += numPowers)
        {
          a_powers[ind] += weightedCurrentProd * a_vecTransform[ivec];
        }
    }
}

#include "NamespaceFooter.H"
