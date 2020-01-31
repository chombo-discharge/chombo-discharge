#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>
using std::endl;

#include "DisjointBoxLayout.H"
#include "LevelData.H"
// #include "BaseFab.H"
// #include "REAL.H"
#include "FArrayBox.H"
#include "DataIterator.H"
#include "LayoutIterator.H"
#include "parstream.H"
#include "AverageF_F.H"

#include "CoarseAverageClaw.H"

#include "UsingNamespace.H"

CoarseAverageClaw::CoarseAverageClaw()
  :
  is_defined(false)
{
}


CoarseAverageClaw::~CoarseAverageClaw()
{
}


CoarseAverageClaw::CoarseAverageClaw(const DisjointBoxLayout& a_fine_domain,
                             int a_numcomps,
                             int a_ref_ratio)
  :
  is_defined(false), m_is_copier_defined(false)
{
  define(a_fine_domain, a_numcomps, a_ref_ratio);
}


CoarseAverageClaw::CoarseAverageClaw(const DisjointBoxLayout& a_fine_domain,
                             const DisjointBoxLayout& a_crse_domain,
                             int a_numcomps,
                             int a_ref_ratio)
  :
  is_defined(false), m_is_copier_defined(false)
{
  define(a_fine_domain, a_crse_domain, a_numcomps, a_ref_ratio);
}


void
CoarseAverageClaw::define(const DisjointBoxLayout& a_fine_domain,
                          int a_numcomps,
                          int a_ref_ratio)
{
    m_ref_ratio = a_ref_ratio;
    DisjointBoxLayout coarsened_fine_domain;
    coarsen (coarsened_fine_domain, a_fine_domain, m_ref_ratio);
    m_coarsened_fine_data.define (coarsened_fine_domain, a_numcomps);
    is_defined = true;
}


void
CoarseAverageClaw::define(const DisjointBoxLayout& a_fine_domain,
                      const DisjointBoxLayout& a_crse_domain,
                      int a_numcomps,
                      int a_ref_ratio)
{
  m_ref_ratio = a_ref_ratio;
  DisjointBoxLayout coarsened_fine_domain;
  coarsen (coarsened_fine_domain, a_fine_domain, m_ref_ratio);
  m_coarsened_fine_data.define (coarsened_fine_domain, a_numcomps);

  // also can pre-define copier here
  // note that since CoarseAverage only operates on interior
  // data, there is no need to worry about whether the domain
  // is periodic.
  m_copier.define(coarsened_fine_domain, a_crse_domain);
  m_is_copier_defined = true;

  is_defined = true;
}



bool
CoarseAverageClaw::isDefined() const
{
  return ( is_defined );
}


void CoarseAverageClaw::averageToCoarse(LevelData<FArrayBox>& a_coarse_data,
                                        const LevelData<FArrayBox>& a_fine_data,
                                        const LevelData<FArrayBox>& a_fine_kappa,
                                        const LevelData<FArrayBox>& a_coarsened_fine_kappa)
{
    CH_assert(is_defined);
    // it would be nice if this could check for validity of a_fine_data.
    // this could be done with a redundant DisjointBoxLayout

    DataIterator dit = a_fine_data.boxLayout().dataIterator();
    for (dit.begin(); dit.ok(); ++dit)

    {
        // coarsenGridData coarsens from the entire fine grid onto the entire
        // coarse grid.
        FArrayBox &coarsened_fine_data = m_coarsened_fine_data[dit()];
        const FArrayBox& fine_data = a_fine_data[dit()];
        const FArrayBox&  coarsened_fine_kappa = a_coarsened_fine_kappa[dit()];
        const FArrayBox& fine_kappa = a_fine_kappa[dit()];
        const Box& b = a_fine_data.disjointBoxLayout().get(dit());

        averageGridData(coarsened_fine_data, fine_data, coarsened_fine_kappa, fine_kappa, b, m_ref_ratio);
    }

    if (m_is_copier_defined)
    {
        // we can use the pre-defined copier to make things faster
        m_coarsened_fine_data.copyTo(m_coarsened_fine_data.interval(),
                                     a_coarse_data,
                                     a_coarse_data.interval(),
                                     m_copier);
    }
    else
    {
        m_coarsened_fine_data.copyTo(m_coarsened_fine_data.interval(),
                                     a_coarse_data,
                                     a_coarse_data.interval() );
    }
}


void CoarseAverageClaw::averageGridData(FArrayBox& a_coarsened_fine_data,
                                        const FArrayBox& a_fine_data,
                                        const FArrayBox& a_coarsened_fine_kappa,
                                        const FArrayBox& a_fine_kappa,
                                        const Box& a_fine_box,
                                        int a_ref_ratio) const
{
    const Box &b = a_fine_box;  // doesn't include ghost cells
    int mxf = b.bigEnd(0) - b.smallEnd(0) + 1;
    int myf = b.bigEnd(1) - b.smallEnd(1) + 1;
    int mxc = mxf/a_ref_ratio;
    int myc = myf/a_ref_ratio;
#if CH_SPACEDIM == 3
    int mzf = b.bigEnd(2) - b.smallEnd(2) + 1;
    int mzc = mzf/a_ref_ratio;
#endif
    int meqn = a_coarsened_fine_data.nComp();
    const IntVect& vf = a_fine_data.size(); // Includes ghost cells
    int mbc = (vf[0] - mxf)/2; // For fine grid only

    Real* qcoarse = a_coarsened_fine_data.dataPtr();  // doesn't have ghost cells!!!
    const Real* qfine = a_fine_data.dataPtr();
    const Real* kappacoarse = a_coarsened_fine_kappa.dataPtr();
    const Real* kappafine = a_fine_kappa.dataPtr();

#if CH_SPACEDIM == 2
    average2_(mxc,myc,mxf,myf,mbc,meqn,qcoarse,qfine,kappacoarse,
              kappafine,a_ref_ratio);
#elif CH_SPACEDIM == 3
    average3_(mxc,myc,mzc,mxf,myf,mzf,mbc,meqn,qcoarse,qfine,kappacoarse,
              kappafine,a_ref_ratio);
#endif
}


void CoarseAverageClaw::averageToCoarse(LevelData<FArrayBox>& a_coarse_data,
                                        const LevelData<FArrayBox>& a_fine_data)
{
    CH_assert(is_defined);
    // it would be nice if this could check for validity of a_fine_data.
    // this could be done with a redundant DisjointBoxLayout
    DataIterator dit = a_fine_data.boxLayout().dataIterator();
    for (dit.begin(); dit.ok(); ++dit)

    {
        // coarsenGridData coarsens from the entire fine grid onto the entire
        // coarse grid.
        FArrayBox &coarsened_fine = m_coarsened_fine_data[dit()];
        const FArrayBox& fine = a_fine_data[dit()];
        averageGridData(coarsened_fine, fine, m_ref_ratio);
    }

    if (m_is_copier_defined)
    {
        // we can use the pre-defined copier to make things faster
        m_coarsened_fine_data.copyTo(m_coarsened_fine_data.interval(),
                                     a_coarse_data,
                                     a_coarse_data.interval(),
                                     m_copier);
    }
    else
    {
        m_coarsened_fine_data.copyTo(m_coarsened_fine_data.interval(),
                                     a_coarse_data,
                                     a_coarse_data.interval() );
    }
}


void CoarseAverageClaw::averageGridData(FArrayBox& a_coarse,
                                        const FArrayBox& a_fine,
                                        int a_ref_ratio) const
{
    const Box& b = a_coarse.box();
    Box refbox(IntVect::Zero,
               (a_ref_ratio-1)*IntVect::Unit);
    FORT_AVERAGE( CHF_FRA(a_coarse),
                  CHF_CONST_FRA(a_fine),
                  CHF_BOX(b),
                  CHF_CONST_INT(a_ref_ratio),
                  CHF_BOX(refbox));
}
