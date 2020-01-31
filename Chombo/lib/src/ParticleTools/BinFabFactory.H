#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//  ANAG, LBNL

#ifndef _BINFABFACTORY_H_
#define _BINFABFACTORY_H_

#include "BoxLayoutData.H"
#include "BinFab.H"
#include "DataIndex.H"
#include "NamespaceHeader.H"

/// Factory class to produce one or more instances of \p BinFab on the same level.
/**
   The mesh and \p ProblemDomain will be the same for all instances created
   by the factory.  This is needed for LevelData construction.
*/
template <class T>
class BinFabFactory
  : public DataFactory<BinFab<T> >
{
public:
  /// Factory function.
  /**
      Creates a new \p BinFab object and returns a pointer to it.
      Responsiblitly for calling operator 'delete' on this pointer
      is passed to the user.

      Assumes \a a_ncomps is 1 and ignores \a a_dit.
  */
  virtual BinFab<T>* create(const Box&       a_box,
                            int              a_numComps,
                            const DataIndex& a_dit) const;

  /// Use these parameters for all manufactured \p BinFab instances.
  /**
     Fills internal data with inputs
  */
  BinFabFactory(const RealVect&      a_meshSpacing,
                const RealVect&      a_origin)
  {
    define(a_meshSpacing,a_origin);
  }

  /// Use these parameters for all manufactured \p BinFab instances, default the rest
  /**
     Fills internal data with inputs, defaulting origin to 0 and promoting meshspacing
     from scalar to RealVect
  */
  BinFabFactory(const Real&          a_meshSpacing)
  {
    define(RealVect::Unit*a_meshSpacing,RealVect::Zero);
  };

  /// Use these parameters for all manufactured \p BinFab instances.
  /**
     Fills internal data with inputs
  */
  void define(const RealVect&      a_meshSpacing,
              const RealVect&      a_origin);

  /// Destructor
  virtual ~BinFabFactory();

private:
  /// cell size in physical coordinates
  RealVect m_mesh_spacing;

  /// location in physical space of the lower corner of the 0 cell
  RealVect m_origin;

};

#include "NamespaceFooter.H"

// Implementation
#include "BinFabFactoryI.H"

#endif
