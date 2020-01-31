#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


#ifndef _POSITIVITYPOSTPROCESSOR_H_
#define _POSITIVITYPOSTPROCESSOR_H_

#include "LevelData.H"
#include "FArrayBox.H"

#include "NamespaceHeader.H"

  /// Applies Hilditch and Colella redistribution scheme to enforce positivity
  /**
  */

class PositivityPostProcessor
{
   public:

      /// Default Constructor.
      /**
       */
      PositivityPostProcessor() :
         m_maximum_number_of_passes(0),
         m_verbose(false)
      {;}

      /// Constructor.
      /**
       * @param[in] a_halo IntVect of number of cells in each direction over
       *   which to redistribute value deficit.  Must be greater than one and
       *   less than or equal to the number of ghost cells for the LevelData
       *   to which the enforcer is applied.
       * @param[in] a_maximum_number_of_passes Maximum number of passes to
       *   attempt in order to redistribute deficit.
       * @param[in] a_verbose Provide verbose details of progress.
       */
      PositivityPostProcessor(
         const IntVect& a_halo,
         const int&     a_maximum_number_of_passes,
         const bool&    a_verbose = false);

      /// Destructor.
      /**
       */
      ~PositivityPostProcessor()
      {
      }

      /// Define the enforcer.
      /**
       * @param[in] a_halo IntVect of number of cells in each direction over
       *   which to redistribute value deficit.  Must be greater than one and
       *   less than or equal to the number of ghost cells for the LevelData
       *   to which the enforcer is applied.
       * @param[in] a_maximum_number_of_passes Maximum number of passes to
       *   attempt in order to redistribute deficit.
       * @param[in] a_verbose Provide verbose details of progress.
       */
      void define(
         const IntVect& a_halo,
         const int&     a_maximum_number_of_passes,
         const bool&    a_verbose = false);

      /// Enforces the ref value through local redistribution of deficit.
      /**
       * @param[in,out] a_phi Scalar or vector cell-averaged field on which
       *   to enforce the ref value.
       * @param[in] a_ref_val Reference value used to estimate zero, i.e.
       *   the value that is to be considered O(1).
       */
      void enforce(
         LevelData<FArrayBox>& a_phi,
         const Real&           a_ref_val = 1.0 ) const;

   private:

      void computeRedistribution( LevelData<FArrayBox>& a_deltaPhi,
                                  const LevelData<FArrayBox>& a_phi,
                                  const Real&    a_ref_val ) const;

      void accumulateRedistribution( LevelData<FArrayBox>& a_deltaPhi,
                                     const Copier& a_reverseCopier ) const;

      void applyRedistribution( LevelData<FArrayBox>& a_phi,
                                const LevelData<FArrayBox>& a_deltaPhi ) const;

      int numberOfNegativeValues( const LevelData<FArrayBox>& a_phi ) const;

      void checkNeighborhood( const IntVect& a_nghosts ) const;

      void checkIfDefined() const;

      PositivityPostProcessor( const PositivityPostProcessor& foo );

      PositivityPostProcessor& operator=( const PositivityPostProcessor& foo );

      Box m_neighborhood;
      int m_maximum_number_of_passes;
      bool m_verbose;
};

#include "NamespaceFooter.H"

#endif
