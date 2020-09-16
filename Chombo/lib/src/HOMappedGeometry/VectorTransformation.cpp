#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "VectorTransformation.H"
#include "MatrixVectorTransformF_F.H"
// #include "UGIO.H"
// #include "CH_HDF5.H"

#include "NamespaceHeader.H"

// ---------------------------------------------------------
RealVect
VectorTransformation::transform(const RealVect& a_vec) const
{
  RealVect rv;
  int icoeff = 0;
  for (int idirDst = 0; idirDst < SpaceDim; idirDst++)
    {
      rv[idirDst] = D_TERM6(m_coeffs[icoeff] * a_vec[0],
                            + m_coeffs[icoeff+1] * a_vec[1],
                            + m_coeffs[icoeff+2] * a_vec[2],
                            + m_coeffs[icoeff+3] * a_vec[3],
                            + m_coeffs[icoeff+4] * a_vec[4],
                            + m_coeffs[icoeff+5] * a_vec[5]);
      icoeff += SpaceDim;
    }
  return rv;
}


// ---------------------------------------------------------
void
VectorTransformation::addTransformFab(FArrayBox& a_dstFab,
                                      const FArrayBox& a_srcFab,
                                      const Box& a_bx,
                                      Real a_weight) const
{
  FORT_VECTORTRANSFORMFABCEN(CHF_FRA(a_dstFab),
                             CHF_CONST_FRA(a_srcFab),
                             CHF_BOX(a_bx),
                             CHF_CONST_R1D(m_coeffs, s_ncoeffs),
                             CHF_CONST_REAL(a_weight));
}


// ---------------------------------------------------------
void
VectorTransformation::addTransformFabAvg(FArrayBox& a_dstFab,
                                         const FArrayBox& a_srcFab,
                                         const FArrayBox& a_gradFab,
                                         const Vector<VectorTransformation>& a_vecTransformGrad,
                                         const Box& a_bx,
                                         Real a_weight) const
{
  FORT_VECTORTRANSFORMFABCEN(CHF_FRA(a_dstFab),
                             CHF_CONST_FRA(a_srcFab),
                             CHF_BOX(a_bx),
                             CHF_CONST_R1D(m_coeffs, s_ncoeffs),
                             CHF_CONST_REAL(a_weight));

  Real scaling = 1./12.;
  Real scaledWeight = a_weight * scaling;
  int nInterpDims = a_vecTransformGrad.size();
  FArrayBox& gradFab = const_cast<FArrayBox&>(a_gradFab);
  for (int idirGrad = 0; idirGrad < nInterpDims; idirGrad++)
    {
      int intvlDirBegin = idirGrad * SpaceDim;
      int intvlDirEnd = intvlDirBegin + SpaceDim-1;
      Interval intvlDir(intvlDirBegin, intvlDirEnd);
      const Real* gradCoeffsDir = a_vecTransformGrad[idirGrad].m_coeffs;
      FArrayBox gradDirFab(intvlDir, gradFab);
      FORT_VECTORTRANSFORMFABCEN(CHF_FRA(a_dstFab),
                                 CHF_CONST_FRA(gradDirFab),
                                 CHF_BOX(a_bx),
                                 CHF_CONST_R1D(gradCoeffsDir, s_ncoeffs),
                                 CHF_CONST_REAL(scaledWeight));
    }
}                          


// ---------------------------------------------------------
VectorTransformation
VectorTransformation::inverse() const
{
  Vector< Vector<Real> > matrixVV = vecvec(m_coeffs);
  Real det = determinant(matrixVV);
  // inv(A) = Adj(A) / det(A)
  // where Adj(A) = C^T, and C[i, j] = (-1)^{i+j} * M[i, j]
  // where M[i, j] = det(A without row i and column j).
  VectorTransformation tfm;
  int isignRow = 1;
  for (int i = 0; i < SpaceDim; i++)
    {
      int isign = isignRow;
      for (int j = 0; j < SpaceDim; j++)
        {
          Real val = isign * matrixMinor(matrixVV, j, i) / det;
          tfm.setComponent(i, j, val);
          isign *= -1;
        }
      isignRow *= -1;
    }
  return tfm;
}


// ---------------------------------------------------------
Vector< Vector<Real> >
VectorTransformation::vecvec(const Real* a_coeffs) const
{
  // change form of matrix to Vector< Vector<Real> >
  Vector< Vector<Real> > matrixVV(SpaceDim);
  int icoeff = 0;
  for (int row = 0; row < SpaceDim; row++)
    {
      for (int col = 0; col < SpaceDim; col++)
        {
          matrixVV[row].push_back(a_coeffs[icoeff]);
          icoeff++;
        }
    }
  return matrixVV;
}


// ---------------------------------------------------------
Real
VectorTransformation::determinant(const Real* a_matrix) const
{
  Vector< Vector<Real> > matrixVV = vecvec(a_matrix);
  Real det = determinant(matrixVV);
  return det;
}


// ---------------------------------------------------------
Real
VectorTransformation::determinant(const Vector< Vector<Real> >& a_matrix) const
{
  int matrixSize = a_matrix.size();
  Real det;
  if (matrixSize == 1)
    {
      det = a_matrix[0][0];
    }
  else
    {
      // Calculate determinant recursively using Laplace expansion:
      // det(A) = sum_{j} { (-1)^j * A[0, j] * det(A without row 0 and column j) }
      det = 0.;
      int isign = 1;
      for (int j = 0; j < matrixSize; j++)
        {
          det += isign * a_matrix[0][j] * matrixMinor(a_matrix, 0, j);
          isign *= -1;
        }
    }
  return det;
}


// ---------------------------------------------------------
Real
VectorTransformation::matrixMinor(const Vector< Vector<Real> >& a_matrix,
                                  int a_row,
                                  int a_col) const
{
  // Return determinant of a_matrix without specified row and column.
  int matrixSize = a_matrix.size();

  Real det;
  if (matrixSize == 1)
    {
      det = a_matrix[0][0];
    }
  else
    {
      // Set reducedMatrix to a_matrix without row a_row and column a_col.
      Vector< Vector<Real> > reducedMatrix(matrixSize - 1);
      int rowNew = 0;
      for (int i = 0; i < matrixSize; i++)
        {
          if (i != a_row)
            {
              const Vector<Real>& matrixRow = a_matrix[i];
              for (int j = 0; j < matrixSize; j++)
                {
                  if (j != a_col)
                    {
                      // rowNew is either i or i-1
                      reducedMatrix[rowNew].push_back(matrixRow[j]);
                    }
                }
              rowNew++;
            }
        }
      det = determinant(reducedMatrix);
    }
  return det;
}


// ---------------------------------------------------------
VectorTransformation
VectorTransformation::compose(const VectorTransformation& a_next) const
{
  VectorTransformation tfmComposite;
  int icoeffComposite = 0;
  for (int idirDst = 0; idirDst < SpaceDim; idirDst++)
    {
      for (int idirSrc = 0; idirSrc < SpaceDim; idirSrc++)
        {
          tfmComposite.m_coeffs[icoeffComposite] = 0.;
          for (int idirMid = 0; idirMid < SpaceDim; idirMid++)
            {
              tfmComposite.m_coeffs[icoeffComposite] +=
                a_next.m_coeffs[matrixIndex(idirDst, idirMid)] *
                m_coeffs[matrixIndex(idirMid, idirSrc)];
            }
          icoeffComposite++;
        }
    }
  return tfmComposite;
}

// ---------------------------------------------------------
int
VectorTransformation::InitStatics()
{
  // taken from IntVect.cpp
  VectorTransformation* ptrID =
    const_cast<VectorTransformation*>( &VectorTransformation::Identity );
  // *ptrID = VectorTransformation();
  // Set identity transformation.
  Real* ptrCoeffs = ptrID->m_coeffs;
  int icoeff = 0;
  for (int idirDst = 0; idirDst < SpaceDim; idirDst++)
    {
      for (int idirSrc = 0; idirSrc < SpaceDim; idirSrc++)
        {
          ptrCoeffs[icoeff] = (idirSrc == idirDst) ? 1. : 0.;
          icoeff++;
        }
    }
  return 0; // arbitrary
}


// ---------------------------------------------------------

VectorTransformation&
VectorTransformation::operator+= (const VectorTransformation& p)
{
  for (int i = 0; i < s_ncoeffs; i++)
    {
      m_coeffs[i] += p.m_coeffs[i];
    }
  return *this;
}

VectorTransformation // friend
operator+ (const VectorTransformation& s,
           const VectorTransformation& p)
{
  VectorTransformation tfm;
  for (int i = 0; i < p.s_ncoeffs; i++)
    {
      tfm.m_coeffs[i] = s.m_coeffs[i] + p.m_coeffs[i];
    }
  return tfm;
}

VectorTransformation&
VectorTransformation::operator-= (const VectorTransformation& p)
{
  for (int i = 0; i < s_ncoeffs; i++)
    {
      m_coeffs[i] -= p.m_coeffs[i];
    }
  return *this;
}

VectorTransformation // friend
operator- (const VectorTransformation& s,
           const VectorTransformation& p)
{
  VectorTransformation tfm;
  for (int i = 0; i < p.s_ncoeffs; i++)
    {
      tfm.m_coeffs[i] = s.m_coeffs[i] - p.m_coeffs[i];
    }
  return tfm;
}

VectorTransformation&
VectorTransformation::operator*= (Real s)
{
  for (int i = 0; i < s_ncoeffs; i++)
    {
      m_coeffs[i] *= s;
    }
  return *this;
}

VectorTransformation
VectorTransformation::operator* (Real s) const
{
  VectorTransformation tfm;
  for (int i = 0; i < s_ncoeffs; i++)
    {
      tfm.m_coeffs[i] = s * m_coeffs[i];
    }
  return tfm;
}

VectorTransformation // friend
operator* (Real            s,
           const VectorTransformation& p)
{
  VectorTransformation tfm;
  for (int i = 0; i < p.s_ncoeffs; i++)
    {
      tfm.m_coeffs[i] = s * p.m_coeffs[i];
    }
  return tfm;
}

VectorTransformation&
VectorTransformation::operator/= (Real s)
{
  for (int i = 0; i < s_ncoeffs; i++)
    {
      m_coeffs[i] /= s;
    }
  return *this;
}

VectorTransformation
VectorTransformation::operator/ (Real s) const
{
  VectorTransformation tfm;
  for (int i = 0; i < s_ncoeffs; i++)
    {
      tfm.m_coeffs[i] = m_coeffs[i] / s;
    }
  return tfm;
}

VectorTransformation // friend
operator/ (Real            s,
           const VectorTransformation& p)
{
  VectorTransformation tfm;
  for (int i = 0; i < p.s_ncoeffs; i++)
    {
      tfm.m_coeffs[i] = p.m_coeffs[i] / s;
    }
  return tfm;
}

// ---------------------------------------------------------

const VectorTransformation VectorTransformation::Identity;

static int s_dummyForVectorTransformationCpp( VectorTransformation::InitStatics() );

#include "NamespaceFooter.H"
