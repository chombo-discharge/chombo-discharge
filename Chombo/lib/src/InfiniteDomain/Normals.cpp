#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// Normals.cpp
// petermc, 19 Aug 2004

#include "Normals.H"
#include "InfiniteDomainF.H"
#include "GetNormalF_F.H"

#include "NamespaceHeader.H"

// ---------------------------------------------------------
// default constructor
Normals::Normals()
{
  m_isDefined = false;
}


// ---------------------------------------------------------
Normals::Normals(const Box&        a_bx,
                 const RealVect&   a_dx,
                 int               a_degree)
{
  define(a_bx, a_dx, a_degree);
}


// ---------------------------------------------------------
Normals::~Normals()
{
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      delete m_coeffs[idir];
      delete m_coeffsExtra[idir];
    }
}


// ---------------------------------------------------------
void
Normals::define(const Box&        a_bx,
                const RealVect&   a_dx,
                int               a_degree)
{
  m_bx = a_bx;
  m_dx = a_dx;
  m_degree = a_degree;
  m_verbose = 0;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      m_coeffs[idir] = new Real[a_degree + 1];
      m_coeffsExtra[idir] = new Real[a_degree + 1];
      Real dxDir = m_dx[idir];
      // FORT_COEFFS1SIDEDDERIV(m_coeffs[idir], &m_degree, &dxDir);
      Real zeroPt = 0.;
      FORT_COEFFSDERIV(m_coeffs[idir], &m_degree, &dxDir, &zeroPt, &zeroPt);
      Real extraPt = -dxDir;
      FORT_COEFFSDERIV(m_coeffsExtra[idir], &m_degree, &dxDir, &extraPt, &zeroPt);
    }
  m_isDefined = true;
}


// ---------------------------------------------------------
bool
Normals::isDefined() const
{
  return m_isDefined;
}


// ---------------------------------------------------------
void
Normals::setVerbose(int a_verbose)
{
  m_verbose = a_verbose;
}


// ---------------------------------------------------------
void
Normals::evalNormal(Vector<FArrayBox*>&   a_deriv,
                    const FArrayBox&      a_phi)
{
  int faceID = 0;
  for (SideIterator sit; sit.ok(); sit.next())
    {
      int iside = sign(sit());
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          FArrayBox& deriv = *(a_deriv[faceID]);
          Box faceBox = deriv.box();
          // usedBox is the box from which data are taken to get deriv.
          Box usedBox = faceBox;
          usedBox.growDir(idir, flip(sit()), m_degree);
          CH_assert(a_phi.box().contains(usedBox));
          FORT_GETNORMAL(CHF_FRA1(deriv, 0),
                         CHF_CONST_FRA1(a_phi, 0),
                         CHF_BOX(faceBox),
                         CHF_CONST_R1D(m_coeffs[idir], m_degree+1),
                         CHF_CONST_INT(idir),
                         CHF_CONST_INT(iside));
          faceID++;
        }
    }
}


// ---------------------------------------------------------
void
Normals::evalNormal(FArrayBox&              a_deriv,
                    const FArrayBox&        a_phi,
                    int                     a_idir,
                    const Side::LoHiSide&   a_side)
{
  int iside = sign(a_side);
  Box faceBox = a_deriv.box();
  // usedBox is the box from which data are taken to get a_deriv.
  Box usedBox = faceBox;
  usedBox.growDir(a_idir, flip(a_side), m_degree);
  CH_assert(a_phi.box().contains(usedBox));
  FORT_GETNORMAL(CHF_FRA1(a_deriv, 0),
                 CHF_CONST_FRA1(a_phi, 0),
                 CHF_BOX(a_deriv.box()),
                 CHF_CONST_R1D(m_coeffs[a_idir], m_degree+1),
                 CHF_CONST_INT(a_idir),
                 CHF_CONST_INT(iside));
}


// ---------------------------------------------------------
void
Normals::evalNormalExtra(Vector<FArrayBox*>&         a_deriv,
                         const FArrayBox&            a_phi,
                         const Vector<FArrayBox*>&   a_phiExtra)
{
  int faceID = 0;
  for (SideIterator sit; sit.ok(); sit.next())
    {
      int iside = sign(sit());
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          FArrayBox& deriv = *(a_deriv[faceID]);
          Box faceBox = deriv.box();
          // usedBox is the box from which data are taken to get deriv.
          Box usedBox = faceBox;
          usedBox.growDir(idir, flip(sit()), m_degree);
          CH_assert(a_phi.box().contains(usedBox));
          FArrayBox& faceExtra = *(a_phiExtra[faceID]);
          FORT_GETNORMALEXTRA(CHF_FRA1(deriv, 0),
                              CHF_CONST_FRA1(a_phi, 0),
                              CHF_CONST_FRA1(faceExtra, 0),
                              CHF_BOX(faceBox),
                              CHF_CONST_R1D(m_coeffsExtra[idir], m_degree+1),
                              CHF_CONST_INT(idir),
                              CHF_CONST_INT(iside));
          faceID++;
        }
    }
}


// ---------------------------------------------------------
void
Normals::evalNormalExtra(FArrayBox&              a_deriv,
                         const FArrayBox&        a_phi,
                         const FArrayBox&        a_phiExtra,
                         int                     a_idir,
                         const Side::LoHiSide&   a_side)
{
  int iside = sign(a_side);
  Box faceBox = a_deriv.box();
  // usedBox is the box from which data are taken to get a_deriv.
  Box usedBox = faceBox;
  usedBox.growDir(a_idir, flip(a_side), m_degree);
  CH_assert(a_phi.box().contains(usedBox));
  FORT_GETNORMALEXTRA(CHF_FRA1(a_deriv, 0),
                      CHF_CONST_FRA1(a_phi, 0),
                      CHF_CONST_FRA1(a_phiExtra, 0),
                      CHF_BOX(a_deriv.box()),
                      CHF_CONST_R1D(m_coeffsExtra[a_idir], m_degree+1),
                      CHF_CONST_INT(a_idir),
                      CHF_CONST_INT(iside));
}

#include "NamespaceFooter.H"
