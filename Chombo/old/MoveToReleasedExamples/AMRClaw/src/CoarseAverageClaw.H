#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _COARSEAVERAGECLAW_H_
#define _COARSEAVERAGECLAW_H_

#include "REAL.H"
// #include "BaseFab.H"
#include "FArrayBox.H"
#include "LevelData.H"

#include "NamespaceHeader.H"

extern "C"
{
#if CH_SPACEDIM == 2
void    average2_(const int& mxc,
                  const int& myc,
                  const int& mxf,
                  const int& myf,
                  const int& mbc,
                  const int& meqn,
                  Real* qcoarse,
                  const Real* qfine,
                  const Real* kappacoarse,
                  const Real* kappafine,
                  const int& refratio);

#elif CH_SPACEDIM == 3
void    average3_(const int& mxc,
                  const int& myc,
                  const int& mzc,
                  const int& mxf,
                  const int& myf,
                  const int& mzf,
                  const int& mbc,
                  const int& meqn,
                  Real* qcoarse,
                  const Real *qfine,
                  const Real *kappacoarse,
                  const Real *kappafine,
                  const int& refratio);

#endif
}  // end of extern "C"

class DisjointBoxLayout;

/// replaces coarse level data with an average of fine level data.
/**
   This class replaces data at a coarse level of refinement with an
   average of data at a finer level of refinement, in areas where fine
   data is present.  Coarse level data is not modified where fine
   level data is not present.

  */
class CoarseAverageClaw
{
public:
  ///
  /**
     Default constructor.  User must subsequently call define().
  */
  CoarseAverageClaw();

  ///
  /**
     Destructor.
  */
  ~CoarseAverageClaw();

  ///
  /**
     Defining constructor.  Constructs a valid object.
     Equivalent to default construction followed by define().

     {\bf Arguments:}\\
     a_fine_domain (not modified): the fine level domain.\\
     a_numcomps (not modified): the number of components.\\
     a_ref_ratio (not modified): the refinement ratio.\\

  */
  CoarseAverageClaw(const DisjointBoxLayout& a_fine_domain,
                int a_numcomps,
                int a_ref_ratio);

  ///
  /**
     Defining constructor.  Constructs a valid object.
     Equivalent to default construction followed by define().
     This version takes the coarser level as well, results in
     a faster averaging operation, since Copier can be pre-constructed.

     {\bf Arguments:}\\
     a_fine_domain (not modified): the fine level domain.\\
     a_crse_domain (not modified): the crse level domain.\\
     a_numcomps (not modified): the number of components.\\
     a_ref_ratio (not modified): the refinement ratio.\\

  */
  CoarseAverageClaw(const DisjointBoxLayout& a_fine_domain,
                const DisjointBoxLayout& a_crse_domain,
                int a_numcomps,
                int a_ref_ratio);

  ///
  /**
     Defines this object.  Existing information is overriden.

     {\bf Arguments:}\\
     a_fine_domain (not modified): the fine level domain.\\
     a_numcomps (not modified): the number of components.\\
     a_ref_ratio (not modified): the refinement ratio.\\

     {\bf This:}\\
     ---This object is modified.---

  */
  void
  define(const DisjointBoxLayout& a_fine_domain, // the fine level domain
         int a_numcomps,                         // the number of components
         int a_ref_ratio);                       // the refinement ratio


  ///
  /**
     Defines this object.  Existing information is overriden.
     This version takes the coarser-level grids as well, which
     allows for a faster averaging operation, since the Copier
     can be pre-defined

     {\bf Arguments:}\\
     a_fine_domain (not modified): the fine level domain.\\
     a_crse_domain (not modified): the coarse level domain.\\
     a_numcomps (not modified): the number of components.\\
     a_ref_ratio (not modified): the refinement ratio.\\

     {\bf This:}\\
     ---This object is modified.---

  */
  void
  define(const DisjointBoxLayout& a_fine_domain, // the fine level domain
         const DisjointBoxLayout& a_crse_domain, // the crse level domain
         int a_numcomps,                         // the number of components
         int a_ref_ratio);                       // the refinement ratio


  ///
  /**
     Returns true if this object was created with the defining
     constructor or if define() has been called.

     {\bf This:}\\
     This object is not modified.
  */
  bool
  isDefined() const;

  ///
  /**
     Replaces a_coarse_data with the average of a_fine_data within the
     coarsening of the domain of the fine level.  Elsewhere, the
     coarse data is unchanged.  It is an error to call if not
     this->isDefined(). The domain of a_fine_data should be
     the same as the fine domain specified in the most recent call to
     define().  It is expected that the coarse and fine level's
     domains are properly nested.  Both a_coarse_data and a_fine_data
     should have the same number of components specified in the most
     recent call to define().

     {\bf Arguments:}\\
     a_coarse_data (modified): coarse data. \\
     a_fine_data (not modified): fine data. \\

     {\bf This:}\\
     Well, it's complicated.  As far as the user is concerned, this object
     is not modified.  See the design document if you care for details.

  */
  //     this method would like to be const, but the work array is changed.
  //     this suggests that the work array should not be persistent.
  void
  averageToCoarse(LevelData<FArrayBox>& a_coarse_data,
                  const LevelData<FArrayBox>& a_fine_data,
                  const LevelData<FArrayBox>& a_coarse_kappa,
                  const LevelData<FArrayBox>& a_coarsened_fine_kappa);

  void
  averageToCoarse(LevelData<FArrayBox>& a_coarse_data,
                  const LevelData<FArrayBox>& a_fine_data);

protected:
    void averageGridData(FArrayBox& a_coarse,
                         const FArrayBox& a_fine,
                         int a_ref_ratio) const;

    void averageGridData(FArrayBox& a_coarse,
                         const FArrayBox& a_fine,
                         const FArrayBox& a_coarse_kappa,
                         const FArrayBox& a_fine_kappa,
                         const Box& b,
                         int a_ref_ratio) const;


protected:
    bool is_defined;

    // the refinement ratio
    int m_ref_ratio;
    // work array for the coarsening of the fine data, of the same "shape"
    // as the fine data.
    LevelData<FArrayBox> m_coarsened_fine_data;

    // has a copier been defined to transfer data to coarse-grid layout?
    bool m_is_copier_defined;

    // cached copier to handle transfer to coarse-grid layout.
    Copier m_copier;

};

#include "NamespaceFooter.H"

#endif
