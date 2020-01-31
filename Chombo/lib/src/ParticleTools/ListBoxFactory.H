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

#ifndef _LISTBOXFACTORY_H_
#define _LISTBOXFACTORY_H_

#include "BoxLayoutData.H"
#include "ListBox.H"
#include "DataIndex.H"
#include "NamespaceHeader.H"

/// Factory class to produce one or more instances of \p ListBox on the
/// same level.
/**
   The mesh and \p ProblemDomain will be the same for all instances created
   by the factory.  This is needed for LevelData construction.
*/
template <class T>
class ListBoxFactory : public DataFactory<ListBox<T> >
{
public:

  /// Weak constructor
  ListBoxFactory();

  /// Constructor. Use these parameters for all manufactured \p ListBox instances.
  ListBoxFactory(const RealVect& a_meshSpacing,
                 const RealVect& a_origin);

  /// Alternative constructor that assumes the mesh spacing is the same in all 
  /// directions and that the domain offset is zero.
  ListBoxFactory(const Real& a_meshSpacing);

  /// Use these parameters for all manufactured \p ListBox instances.
  void define(const RealVect& a_meshSpacing,
              const RealVect& a_origin);

  /// Factory function. Note that both a_numComps and a_dit are ignored.
  virtual ListBox<T>* create(const Box& a_box, int a_numComps, 
			     const DataIndex& a_datInd) const;
  /// Destructor
  virtual ~ListBoxFactory();

  /// Retrieve the cell size
  RealVect meshSpacing() const;

  /// Retrieve the origin
  RealVect origin() const;

private:
  /// cell size in physical coordinates
  RealVect m_mesh_spacing;

  /// location in physical space of the lower corner of the 0 cell
  RealVect m_origin;
};

#include "NamespaceFooter.H"

// Implementation
#include "ListBoxFactoryI.H"

#endif
