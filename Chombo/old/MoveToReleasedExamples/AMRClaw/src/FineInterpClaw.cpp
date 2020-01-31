#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "BaseFab.H"
#include "REAL.H"
#include "DataIterator.H"
#include "Tuple.H"
#include "InterpF_F.H"

#include "FineInterpClaw.H"

#include "UsingNamespace.H"

FineInterpClaw::FineInterpClaw()
  :
  is_defined(false)
{
}


FineInterpClaw::~FineInterpClaw()
{
}


FineInterpClaw::FineInterpClaw(const DisjointBoxLayout& a_fine_domain,
                       const int&  a_numcomps,
                       const int& a_ref_ratio,
                       const Box& a_fine_problem_domain)
  :
  is_defined(false)
{
  ProblemDomain fineProbDomain(a_fine_problem_domain);
  define(a_fine_domain, a_numcomps, a_ref_ratio, fineProbDomain);
}

FineInterpClaw::FineInterpClaw(const DisjointBoxLayout& a_fine_domain,
                       const int&  a_numcomps,
                       const int& a_ref_ratio,
                       const ProblemDomain& a_fine_problem_domain)
  :
  is_defined(false)
{
  define(a_fine_domain, a_numcomps, a_ref_ratio, a_fine_problem_domain);
}


void
FineInterpClaw::define(const DisjointBoxLayout& a_fine_domain,
                   const int& a_numcomps,
                   const int& a_ref_ratio,
                   const Box& a_fine_problem_domain)
{
  ProblemDomain fineProbDomain(a_fine_problem_domain);
  define(a_fine_domain, a_numcomps, a_ref_ratio, fineProbDomain);
}

void
FineInterpClaw::define(const DisjointBoxLayout& a_fine_domain,
                   const int& a_numcomps,
                   const int& a_ref_ratio,
                   const ProblemDomain& a_fine_problem_domain)
{
  // check for consistency
 CH_assert (a_fine_domain.checkPeriodic(a_fine_problem_domain));
  m_ref_ratio = a_ref_ratio;
  m_coarse_problem_domain = coarsen(a_fine_problem_domain, m_ref_ratio);
  //
  // create the work array
  DisjointBoxLayout coarsened_fine_domain;
  coarsen ( coarsened_fine_domain,
            a_fine_domain,
            m_ref_ratio );
  m_coarsened_fine_data.define ( coarsened_fine_domain,
                                 a_numcomps,
                                 IntVect::Unit );
  is_defined = true;
}



bool
FineInterpClaw::isDefined() const
{
  return ( is_defined );
}


// interpolate from coarse level to fine level
void FineInterpClaw::interpToFine(LevelData<FArrayBox>& a_fine_data,
                         const LevelData<FArrayBox>& a_coarse_data)
{
 CH_assert(is_defined);
#ifndef NDEBUG
  // debugging check
  {
    DataIterator crseDit = m_coarsened_fine_data.dataIterator();
    for (crseDit.reset(); crseDit.ok(); ++crseDit)
      {
        m_coarsened_fine_data[crseDit()].setVal(1.0e9);
      }
  }
#endif

  // this should handle all the periodic BCs as well,
  // by filling in the ghost cells in an appropriate way
  a_coarse_data.copyTo(a_coarse_data.interval(),
                       m_coarsened_fine_data,
                       m_coarsened_fine_data.interval() );

  const BoxLayout fine_domain = a_fine_data.boxLayout();
  DataIterator dit = fine_domain.dataIterator();

  for (dit.begin(); dit.ok(); ++dit)
    {
      const FArrayBox& coarsened_fine = m_coarsened_fine_data[dit()];
      const Box& coarsened_fine_box = m_coarsened_fine_data.getBoxes()[dit()];
      FArrayBox& fine = a_fine_data[dit()];
      // interpGridData interpolates from an entire coarse grid onto an
      // entire fine grid.
      interpGridData(fine,
                     coarsened_fine,
                     coarsened_fine_box,
                     m_ref_ratio);
    }
}


void FineInterpClaw::interpToFine(LevelData<FArrayBox>& a_fine_data,
                                  const LevelData<FArrayBox>& a_coarse_data,
                                  const LevelData<FArrayBox>& a_fine_kappa,
                                  const LevelData<FArrayBox>& a_coarsened_fine_kappa)
{
 CH_assert(is_defined);
#ifndef NDEBUG
  // debugging check
  {
    DataIterator crseDit = m_coarsened_fine_data.dataIterator();
    for (crseDit.reset(); crseDit.ok(); ++crseDit)
      {
        m_coarsened_fine_data[crseDit()].setVal(1.0e9);
      }
  }
#endif

  // this should handle all the periodic BCs as well,
  // by filling in the ghost cells in an appropriate way
  a_coarse_data.copyTo(a_coarse_data.interval(),
                       m_coarsened_fine_data,
                       m_coarsened_fine_data.interval() );

  const BoxLayout fine_domain = a_fine_data.boxLayout();
  DataIterator dit = fine_domain.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& fine = a_fine_data[dit()];
      const FArrayBox& coarsened_fine = m_coarsened_fine_data[dit()];
      // const Box& coarsened_fine_box = m_coarsened_fine_data.getBoxes()[dit()];
      // interpGridData interpolates from an entire coarse grid onto an
      // entire fine grid.
      const Box& coarsened_fine_box = m_coarsened_fine_data.getBoxes()[dit()];
      const Box& b = a_fine_data.disjointBoxLayout().get(dit());
      interpGridData(fine,
                     coarsened_fine,
                     coarsened_fine_box,
                     m_ref_ratio);

      // Do the fixup needed when capacity function is used.
      const FArrayBox& fineKappa = a_fine_kappa[dit()];
      const FArrayBox& coarseKappa = a_coarsened_fine_kappa[dit()];
      // const Box& b = a_fine_data.disjointBoxLayout().get(dit());

      interpFixup(fine,fineKappa,coarsened_fine,coarseKappa,b,m_ref_ratio);
    }
}


// interpolate from fine grid to coarse grid.  prerequisite:
// coarsened.box contains coarsen(fine.box).
//
// uses piecewise bilinear interpolation with multidimensional-limited
// slopes.  see design document for details.
void FineInterpClaw::interpGridData(FArrayBox& a_fine,
                           const FArrayBox& a_coarse,
                           const Box& a_coarsened_fine_box,
                           int a_ref_ratio)
  const
{
  // fill fine data with piecewise constant coarse data
  const Box& b = a_coarsened_fine_box;
  const int num_comp = a_fine.nComp ();
  Box refbox(IntVect::Zero,
             (a_ref_ratio-1)*IntVect::Unit);

  FORT_INTERPCONSTANT ( CHF_FRA(a_fine),
                        CHF_CONST_FRA(a_coarse),
                        CHF_BOX(b),
                        CHF_CONST_INT(a_ref_ratio),
                        CHF_BOX(refbox)
                        );
  //  Tuple<BaseFab<Real>, SpaceDim> slopes;
  //  for (int dir = 0; dir < SpaceDim; ++dir)
  // hardwired to 3 due to lack of variable number of arguments in chfpp
  BaseFab<Real> slopes[3];
  for (int dir = 0; dir < 3; ++dir)
    {
      BaseFab<Real>& dir_slope = slopes[dir];
      dir_slope.resize(b, num_comp);
    }
  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      BaseFab<Real>& dir_slope = slopes[dir];

      const Box bcenter = grow(m_coarse_problem_domain,-BASISV(dir)) & b;
      if (!bcenter.isEmpty())
        {
          FORT_INTERPCENTRALSLOPE ( CHF_FRA ( dir_slope ),
                                    CHF_CONST_FRA ( a_coarse ),
                                    CHF_BOX ( bcenter ),
                                    CHF_CONST_INT ( dir )
                                    );
        }
      const Box blo = b & adjCellLo(grow(m_coarse_problem_domain,-BASISV(dir)),dir);
      if (!blo.isEmpty())
        {
          FORT_INTERPHISIDESLOPE ( CHF_FRA ( dir_slope ),
                                   CHF_CONST_FRA ( a_coarse ),
                                   CHF_BOX ( blo ),
                                   CHF_CONST_INT ( dir )

                                   );
        }
      const Box bhi = b & adjCellHi(grow(m_coarse_problem_domain,-BASISV(dir)),dir);
      if (!bhi.isEmpty())
        {
          FORT_INTERPLOSIDESLOPE ( CHF_FRA ( dir_slope ),
                                   CHF_CONST_FRA ( a_coarse ),
                                   CHF_BOX ( bhi ),
                                   CHF_CONST_INT ( dir )
                                   );
        }
    }

  // to do limits, we need to have a box which includes
  // the neighbors of a given point (to check for the
  // local maximum...
  Box neighborBox(-1*IntVect::Unit,
                  IntVect::Unit);

  // GHM 7/12/01
  // interplimit iterates over box b_mod (was b), but cells within
  // 1 of the physical boundary never enter result (and this
  // wasted calculation may call upon uninitialized memory).
  // DFM 10/8/01
  // note that this turns off slope limiting for cells adjacent to the
  // boundary -- may want to revisit this in the future
  Box b_mod(b);
  b_mod.grow(1);
  b_mod = m_coarse_problem_domain & b_mod;
  b_mod.grow(-1);

  // create a box grown big enough to remove periodic BCs from domain
  Box domBox = grow(b, 2);
  domBox = m_coarse_problem_domain & domBox;

  FORT_INTERPLIMIT ( CHF_FRA ( slopes[0] ),
                     CHF_FRA ( slopes[1] ),
                     CHF_FRA ( slopes[2] ),
                     CHF_CONST_FRA ( a_coarse ),
                     CHF_BOX ( b_mod ),
                     CHF_BOX ( neighborBox ),
                     CHF_BOX (domBox)
                     );

  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      BaseFab<Real>& dir_slope = slopes[dir];

      FORT_INTERPLINEAR ( CHF_FRA ( a_fine ),
                          CHF_CONST_FRA ( dir_slope ),
                          CHF_BOX ( b ),
                          CHF_CONST_INT ( dir ),
                          CHF_CONST_INT ( a_ref_ratio ),
                          CHF_BOX ( refbox )
                          );
    }
}

void FineInterpClaw::interpFixup(FArrayBox& a_fine_data,
                                 const FArrayBox& a_fine_kappa,
                                 const FArrayBox& a_coarse_data,
                                 const FArrayBox& a_coarse_kappa,
                                 const Box& a_fine_box,
                                 const int& a_refRatio)
{

    int mxf = a_fine_box.bigEnd(0) - a_fine_box.smallEnd(0) + 1;
    int myf = a_fine_box.bigEnd(1) - a_fine_box.smallEnd(1) + 1;
    int mxc = mxf/a_refRatio;
    int myc = myf/a_refRatio;
#if CH_SPACEDIM == 3
    int mzf = a_fine_box.bigEnd(2) - a_fine_box.smallEnd(2) + 1;
    int mzc = mzf/a_refRatio;
#endif

    int meqn = a_fine_data.nComp();
    const IntVect& vf = a_fine_data.size(); // Includes ghost cells
    int mbc = (vf[0] - mxf)/2; // For fine grid only

          Real* q_fine       = a_fine_data.dataPtr();
    const Real* q_coarse     = a_coarse_data.dataPtr();
    const Real* kappa_fine   = a_fine_kappa.dataPtr();
    const Real* kappa_coarse = a_coarse_kappa.dataPtr();

#if CH_SPACEDIM == 2
    fixcapaq2_(mxc,myc,mxf,myf, mbc,meqn, q_coarse,
               kappa_coarse,q_fine,kappa_fine,a_refRatio);
#elif CH_SPACEDIM == 3
        fixcapaq3_(mxc,myc,mzc,mxf,myf,mzf,mbc,meqn,q_coarse,
               kappa_coarse,q_fine,kappa_fine,a_refRatio);
#endif

}
