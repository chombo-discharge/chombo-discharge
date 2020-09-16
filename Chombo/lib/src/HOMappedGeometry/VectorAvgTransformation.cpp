#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "VectorAvgTransformation.H"
#include "MatrixVectorTransformF_F.H"
// #include "UGIO.H"
// #include "CH_HDF5.H"

#include "NamespaceHeader.H"

// ---------------------------------------------------------
VectorAvgTransformation::VectorAvgTransformation(const VectorTransformation& a_tfmAvg,
                                                 const Vector<VectorTransformation>& a_tfmGrad)
{
  m_interpDims = a_tfmGrad.size();
  int ncoeffs1 = SpaceDim*SpaceDim;
  int ncoeffsAll = (1 + m_interpDims) * ncoeffs1;
  m_coeffsAll.resize(ncoeffsAll);

  Real scaling = 1./12.;
  const Real* coeffsAvg = a_tfmAvg.m_coeffs;
  for (int i = 0; i < ncoeffs1; i++)
    {
      m_coeffsAll[i] = coeffsAvg[i];
    }
  int icoeff = ncoeffs1;
  for (int idirGrad = 0; idirGrad < m_interpDims; idirGrad++)
    {
      const Real* gradCoeffsDir = a_tfmGrad[idirGrad].m_coeffs;
      for (int i = 0; i < ncoeffs1; i++)
        {
          m_coeffsAll[icoeff] = scaling * gradCoeffsDir[i];
          icoeff++;
        }
    }
}

// ---------------------------------------------------------
RealVect
VectorAvgTransformation::transform(const RealVect& a_vec,
                                   const Vector<RealVect>& a_gradVec) const
{
  RealVect rv;
  int icoeff = 0;
  for (int idirDst = 0; idirDst < SpaceDim; idirDst++)
    {
      rv[idirDst] = D_TERM6(m_coeffsAll[icoeff] * a_vec[0],
                            + m_coeffsAll[icoeff+1] * a_vec[1],
                            + m_coeffsAll[icoeff+2] * a_vec[2],
                            + m_coeffsAll[icoeff+3] * a_vec[3],
                            + m_coeffsAll[icoeff+4] * a_vec[4],
                            + m_coeffsAll[icoeff+5] * a_vec[5]);
      icoeff += SpaceDim;
    }
  for (int idirGrad = 0; idirGrad < m_interpDims; idirGrad++)
    {
      const RealVect& grad = a_gradVec[idirGrad];
      for (int idirDst = 0; idirDst < SpaceDim; idirDst++)
        {
          rv[idirDst] += D_TERM6(m_coeffsAll[icoeff] * grad[0],
                                 + m_coeffsAll[icoeff+1] * grad[1],
                                 + m_coeffsAll[icoeff+2] * grad[2],
                                 + m_coeffsAll[icoeff+3] * grad[3],
                                 + m_coeffsAll[icoeff+4] * grad[4],
                                 + m_coeffsAll[icoeff+5] * grad[5]);
          icoeff += SpaceDim;
        }
    }
  return rv;
}


// ---------------------------------------------------------
void
VectorAvgTransformation::addTransformFab(FArrayBox& a_dstFab,
                                         const FArrayBox& a_srcFab,
                                         const FArrayBox& a_gradFab,
                                         const Box& a_bx,
                                         Real a_weight) const
{
  FORT_VECTORTRANSFORMFABAVG(CHF_FRA(a_dstFab),
                             CHF_CONST_FRA(a_srcFab),
                             CHF_CONST_FRA(a_gradFab),
                             CHF_BOX(a_bx),
                             CHF_CONST_INT(m_interpDims),
                             CHF_CONST_VR(m_coeffsAll),
                             CHF_CONST_REAL(a_weight));
}


#include "NamespaceFooter.H"
