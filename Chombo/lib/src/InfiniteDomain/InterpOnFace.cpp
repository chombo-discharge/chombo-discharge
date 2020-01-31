#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// InterpOnFace.cpp
// petermc, 20 Aug 2004

#include "InterpOnFace.H"
// #include "BoxIterator.H"
// #include "InfiniteDomainF.H"
#include "Projections.H"
#include "InterpOnFaceF_F.H"

#include "NamespaceHeader.H"

// ---------------------------------------------------------
// default constructor
InterpOnFace::InterpOnFace()
{
  m_isDefined = false;
}


// ---------------------------------------------------------
InterpOnFace::InterpOnFace(/// NODE-centered face
                           const Box&   a_fineFaceBox,
                           /// normal direction of the face
                           int          a_idir,
                           /// coarsening ratio
                           int          a_coarsening,
                           /// width of border around face
                           int          a_interpBorder)
{
  define(a_fineFaceBox, a_idir, a_coarsening, a_interpBorder);
}


// ---------------------------------------------------------
InterpOnFace::~InterpOnFace()
{
}


// ---------------------------------------------------------
void
InterpOnFace::define(/// NODE-centered face in fine coordinates
                     const Box&   a_fineFaceBox,
                     /// normal direction of the face
                     int          a_idir,
                     /// coarsening ratio
                     int          a_coarsening,
                     /// width of border around face
                     int          a_interpBorder)
{
  m_coarsening = a_coarsening;
  m_fineFaceBox = a_fineFaceBox;
  m_fineCoarseShift = m_fineFaceBox.smallEnd();
  // NODE-centered
  // Box crseFaceBoxBasic = coarsen(m_fineFaceBox - m_fineCoarseShift,
  // m_coarsening);

  // Because of the shift by m_fineCoarseShift,
  // crseFaceBoxBasic will have lower corner at (0, 0, 0).
  Box crseFaceBoxBasic = coarsenOuter(m_fineFaceBox - m_fineCoarseShift,
                                      m_coarsening);

  // The lower corner of crseFaceBoxBasicRefined will be at (0, 0, 0),
  // but the upper corner may be beyond
  // that of m_fineFaceBox - m_fineCoarseShift.
  m_crseFaceBoxBasicRefined = refine(crseFaceBoxBasic, a_coarsening);

  // Check that m_fineFaceBox - m_fineCoarseShift is
  // coarsenable by a_coarsening.
  // petermc, 3 Sep 2004
  // Previously, required m_fineFaceBox to be coarsenable by a_coarsening,
  // but that is not necessary, and disallows the case of the corners
  // not being on coarse points.
  m_coarsenedExact = (m_crseFaceBoxBasicRefined + m_fineCoarseShift ==
                      m_fineFaceBox);

  if (! m_coarsenedExact)
    {
      Box fineShifted = m_fineFaceBox - m_fineCoarseShift;
      // Example:  (6:19)/4
      // m_fineCoarseShift = 6
      // fineShifted = (0:13)
      // crseFaceBoxBasic = (0:4)
      // m_crseFaceBoxBasicRefined = (0:16)
      // remainder = 16 - 13 = 3
      // m_offLoExt = 3 / 2 = 1
      IntVect remainder =
        m_crseFaceBoxBasicRefined.bigEnd() - fineShifted.bigEnd();
      m_offLoExt = remainder / 2;
      // Then when you want to interpolate from (0:4) to (6:19):
      // Fill in (0:4) to (0:16).
      // Then shift by -m_offLoExt = -1 so that it's (-1:15).
      // Then shift by m_fineCoarseShift = 6 so that it's (5:21).
      // Then copy the same indices from (5:21) to (6:19).
    }

  m_idir = a_idir;
  m_ipar1 = (m_idir == 0) ? 1 : 0;
  m_ipar2 = (m_idir == 2) ? 1 : 2;
  m_interpBorder = a_interpBorder;

  // We'll interpolate from coarse data on m_coarseFaceBox.
  IntVect ivGrowParallel = m_interpBorder *
    (BASISV(m_ipar1) + BASISV(m_ipar2));
  m_coarseFaceBox = grow(crseFaceBoxBasic, ivGrowParallel);

  if (m_coarsening != 1)
    {
      m_icbx1 = Box(-m_interpBorder*BASISV(m_ipar1),
                    (m_interpBorder+1)*BASISV(m_ipar1));
      m_weights1.define(m_icbx1, m_coarsening-1);
      FORT_GETWEIGHTSINTERP(CHF_FRA(m_weights1),
                            CHF_BOX(m_icbx1),
                            CHF_CONST_INT(m_coarsening),
                            CHF_CONST_INT(m_interpBorder),
                            CHF_CONST_INT(m_ipar1));
      m_ifbx1 = Box(BASISV(m_ipar1), (m_coarsening-1)*BASISV(m_ipar1));
      m_cbx1n = crseFaceBoxBasic;
      m_cbx1n.grow(m_ipar2, m_interpBorder);
      m_cbx1e = enclosedCells(m_cbx1n, m_ipar1);

      m_icbx2 = Box(-m_interpBorder*BASISV(m_ipar2),
                    (m_interpBorder+1)*BASISV(m_ipar2));
      m_weights2.define(m_icbx2, m_coarsening-1);
      FORT_GETWEIGHTSINTERP(CHF_FRA(m_weights2),
                            CHF_BOX(m_icbx2),
                            CHF_CONST_INT(m_coarsening),
                            CHF_CONST_INT(m_interpBorder),
                            CHF_CONST_INT(m_ipar2));
      m_ifbx2 = Box(BASISV(m_ipar2), (m_coarsening-1)*BASISV(m_ipar2));
      // m_cfbx2e is crseFaceBoxBasic, with
      // NODES refined by m_coarsening in direction m_ipar1,
      // CELLS in direction m_ipar2.
      IntVect ivRefine2 = BASISV(m_idir) + BASISV(m_ipar2) +
        m_coarsening * BASISV(m_ipar1);
      m_cfbx2n = refine(crseFaceBoxBasic, ivRefine2);
      m_cfbx2n.shift(m_ipar1, m_fineCoarseShift[m_ipar1]);
      m_cfbx2e = enclosedCells(m_cfbx2n, m_ipar2);

      IntVect ivGrow2 = m_interpBorder * BASISV(m_ipar2);
      Box m_workbox = grow(m_cfbx2n, ivGrow2);
      m_work.define(m_workbox, m_coarsening-1);
    }
  m_verbose = 0;
  m_isDefined = true;
}


// ---------------------------------------------------------
bool
InterpOnFace::isDefined() const
{
  return m_isDefined;
}


// ---------------------------------------------------------
void
InterpOnFace::setVerbose(int a_verbose)
{
  m_verbose = a_verbose;
}


// ---------------------------------------------------------
void
InterpOnFace::interpolate(FArrayBox&         a_fine,
                          const FArrayBox&   a_coarse)
{
  CH_assert(a_fine.box() == m_fineFaceBox);
  CH_assert(a_coarse.box() == m_coarseFaceBox);
  if (m_coarsening == 1)
    {
      a_fine.shift(-m_fineCoarseShift);
      a_fine.copy(a_coarse);
      a_fine.shift(m_fineCoarseShift);
    }
  else if (m_coarsenedExact)
    {
      FORT_INTERP2DFACE(CHF_FRA1(a_fine, 0),
                        CHF_CONST_FRA1(a_coarse, 0),
                        CHF_CONST_INT(m_coarsening),
                        CHF_CONST_INT(m_idir),
                        CHF_BOX(m_cbx1e),
                        CHF_BOX(m_cbx1n),
                        CHF_BOX(m_icbx1),
                        CHF_BOX(m_ifbx1),
                        CHF_CONST_FRA(m_weights1),
                        CHF_FRA1(m_work, 0),
                        CHF_BOX(m_cfbx2e),
                        CHF_BOX(m_cfbx2n),
                        CHF_BOX(m_icbx2),
                        CHF_BOX(m_ifbx2),
                        CHF_CONST_FRA(m_weights2),
                        CHF_CONST_INTVECT(m_fineCoarseShift));
    }
  else
    {
      FArrayBox fineExtended(m_crseFaceBoxBasicRefined + m_fineCoarseShift, 1);
      FORT_INTERP2DFACE(CHF_FRA1(fineExtended, 0),
                        CHF_CONST_FRA1(a_coarse, 0),
                        CHF_CONST_INT(m_coarsening),
                        CHF_CONST_INT(m_idir),
                        CHF_BOX(m_cbx1e),
                        CHF_BOX(m_cbx1n),
                        CHF_BOX(m_icbx1),
                        CHF_BOX(m_ifbx1),
                        CHF_CONST_FRA(m_weights1),
                        CHF_FRA1(m_work, 0),
                        CHF_BOX(m_cfbx2e),
                        CHF_BOX(m_cfbx2n),
                        CHF_BOX(m_icbx2),
                        CHF_BOX(m_ifbx2),
                        CHF_CONST_FRA(m_weights2),
                        CHF_CONST_INTVECT(m_fineCoarseShift));
      // CHF_CONST_INTVECT(IntVect::Zero));
      fineExtended.shift(-m_offLoExt);
      a_fine.copy(fineExtended);
    }
}

#include "NamespaceFooter.H"
