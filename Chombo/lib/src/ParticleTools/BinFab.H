#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _BINFAB_H_
#define _BINFAB_H_

#include "BaseFab.H"
#include "RealVect.H"
#include "FluxBox.H"
#include "List.H"
#include "Box.H"
#include "NamespaceHeader.H"


/// Base class for particle data on a \p Box
/** This is a class for holding and sorting particle items.  Class <T>
    must have a \p RealVect \b <T>::position() const function which is used
    to place the item in an appropriate bin.  See class \p BinItem for a base
    class for particles.  All BinFabs are assumed to have 1 component.
*/

template<class T>
class BinFab : public BaseFab<List<T> >
{
public:
  /// Null constructor, copy constructor and operator= can be compiler defined.
  BinFab();

  /// Constructs an empty \p BinFab on \a a_domain.
  BinFab(const Box&           a_domain,
         const RealVect&      a_mesh_spacing,
         const RealVect&      a_origin);

  /// Copy constructor
  BinFab(const BinFab& a_binfab);

#ifdef __IBMCPP__
  //[NOTE: this is here only because the IBM xlC compiler tries to compile
  //       DefaultDataFactory<BinFab<BinItem>> in particleTest.cpp, even
  //       though it isn't used anywhere.  Stupid compiler. <dbs> July2005]
  BinFab(const Box& a_domain, int a_ncomps);
#endif

  /// Destructor
  virtual ~BinFab();

  /// Same as the constructor
  virtual void define(const Box&           a_domain,
                      const RealVect&      a_mesh_spacing,
                      const RealVect&      a_origin);

  /// Sort all the items into the appropriate bins, losing items not in the FAB's domain.
  /** This is used after "pushing" particles (i.e. updating physical positions).
   */
  virtual void reBin();

  /// Sort all the items into the appropriate bins and return items in (or not in) the \p Box
  /** This is used after "pushing" particles, to find which have left the domain of the \b BinFab.
   */
  virtual void reBin(List<T>& a_lost, const Box & a_valid = Box(), bool a_in = true);

  /// Copy one item into the appropriate bin.
  virtual void addItem(const T& a_item);

  /// Copy one item into the specified bin.
  virtual void addItem(const T& a_item, const IntVect& a_bin);

  /// Copy items into the appropriate bins.
  virtual void addItems(const List<T>& a_list);

  /// Copy items into the appropriate bins and remove them from the \p List.
  virtual void addItemsDestructive(List<T>& a_list);

  /// Copy items into the appropriate bins.  Remove items that are in (or not in) the \p Box from the \p List.
  virtual void addItemsDestructive(List<T>& a_list, const Box& a_valid, bool a_in = true);

  /// Move particles from a \p BinFab to the overlapping cells in this one.
  /** Performs a point-wise transfer of the whole \p List<T> in src cells
      to the corresponding cells in *this.  Notice that as a result both
      *this and a_src are modified.
  */
  virtual void transfer(BinFab<T>&  a_src,
                        const Box&  a_srcBox,
                        const Box&  a_destBox);

//XXX  //XXX testing if copy() can be overloaded to handle transfer() <dbs>
//XXX  BinFab<T> & copy(BinFab<T>& a_src,
//XXX                   const Box& a_srcBox,
//XXX                   const Box& a_destBox);
//XXX

  /// Delete all the items in this \p BinFab and reset it to the null \p Box.
  virtual void clear();

  /// Count the number of items in this \p BinFab in the specified \p Box.
  /** The default is to do the whole FAB.
      This is slow, so use it sparingly.
   */
  virtual int numItems(const Box& a_box = Box()) const;

  /// Retrieve the mesh cell size, in physical coordinates.
  virtual RealVect meshSpacing() const;

  /// Set the mesh cell size, in physical coordinates.
  virtual void meshSpacing(const RealVect& a_mesh_spacing) ;

  /// Retrieve the physical coordinates of the lower corner of cell zero.
  RealVect origin() const;

  /// Set the physical coordinates of the lower corner of cell zero.
  void origin(const RealVect& a_origin);

  // Linearization functions

  /// The number of bytes used by \b linearIn()/linearOut().
  /** Returns the size, in number of bytes, of a flat linear
      representation of the items in the specified \p Box.
  */
  virtual size_t size(const Box& a_box, const Interval& a_comps = Interval(0,0)) const;

  /// Write a serialized (packed) representation into \a a_buf.
  /** Write a linear representation of the items in the specified \p Box.
      Assumes that sufficient memory for the buffer has already been
      allocated by the caller.
  */
  virtual void linearOut(void* a_buf, const Box& a_box,
                         const Interval& a_comps) const;

  /// Write a serialized (packed) representation into \a a_buf and remove from *this.
  /** Write a linear representation of the items in the specified \p Box and
      removes them from this \p BinFab.  Assumes that sufficient memory for
      the buffer has already been allocated by the caller.
  */
  virtual void linearOutDestructive(void* buf, const Box& a_box,
                                    const Interval& a_comps = Interval(0,0));

  /// Extract a serialized (packed) representation from \a a_buf.
  /** Given a linear representation of the class data previously made
      using linearIn() in \a a_buf, set the data for this class on the
      specified \p Box.
  */
  virtual void linearIn(void* a_buf, const Box& a_box,
                        const Interval& a_comps = Interval(0,0));

  ///
  static int preAllocatable()
  {
    return 2;
  }

  /// Should be inherited from BaseFab but it isn't, sigh; code is the same.
  /** needed for loop-unroll macro */
  IntVect size () const
  {
    return this->m_domain.size();
  }

protected:
  /// location in physical space of the 0 index
  RealVect m_origin;

  /// cell size in physical space
  RealVect m_mesh_spacing;

  /// compute the cell index containing the physical position of the item
  inline virtual IntVect locateBin(const T& a_item) const;
};

// Implementation

#include "NamespaceFooter.H"
#include "BinFabImplem.H"

#endif
