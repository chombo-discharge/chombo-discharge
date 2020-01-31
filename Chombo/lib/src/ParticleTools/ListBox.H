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

#ifndef _LISTBOX_H_
#define _LISTBOX_H_

#include "Arena.H"
#include "Interval.H"
#include "RealVect.H"
#include "List.H"
#include "Box.H"
#include "NamespaceHeader.H"

/** A class for storing and sorting particle items. Particles are held in a List 
    and belong to this ListBox if their position is within the domain represented 
    by the object Box..  Class <T> must have a \p RealVect \b <T>::position() const 
    function which is used to assign the item to the ListBox.  See class \p BinItem
    for a base class for particles.
*/

template<class T>
class ListBox
{
public:

  /// Null constructor.
  ListBox();

  /// Constructs an empty \p ListBox on \a a_box. a_meshSpacing is the cell size in 
  /// each direction, and a_origin is the location in physical coordinates of the 
  /// lower left corner of the problem domain.
  ListBox(const Box& a_box,
          const RealVect& a_meshSpacing,
          const RealVect& a_origin);

  /// Destructor.
  virtual ~ListBox();

  /// Same as the constructor.
  virtual void define(const Box& a_box,
                      const RealVect& a_meshSpacing,
                      const RealVect& a_origin);

  /// Free the memory associated with this ListBox.
  virtual void undefine();

  /// Returns the list of items contained in this ListBox.
  inline List<T>& listItems();

  /// Returns the list of items contained in this ListBox, const version.
  inline const List<T>& listItems() const;

  /// Copy items from this ListBox into a_list only if they are contained in a_valid.
  /// This operation does not remove the items from this ListBox. 
  virtual void getItems(List<T>& a_list, const Box& a_valid) const;

  /// Transfer items from this ListBox into a_list only if they are contained in a_valid.
  /// This operation removes the items from this ListBox. 
  virtual void getItemsDestructive(List<T>& a_list, const Box& a_valid);

  /// Copy items from this ListBox into a_list if they are NOT contained in a_valid.
  /// Does not remove the items from this ListBox. 
  virtual void getInvalid(List<T>& a_list, const Box& a_valid) const;

  /// Transfer items from this ListBox into a_list if they are NOT contained in a_valid.
  /// This removes the items from this ListBox. 
  virtual void getInvalidDestructive(List<T>& a_list, const Box& a_valid);

  /// Add a_item to this ListBox if it is contained within m_box.
  virtual void addItem(const T& a_item);

  /// Add each item in a_list into this ListBox if they are contained within m_box.
  /// Does not remove the items from a_list.
  virtual void addItems(const List<T>& a_list);

  /// Transfer each item in a_list into this ListBox if it is contained within m_box.
  /// This removes the items from a_list.
  virtual void addItemsDestructive(List<T>& a_list);

  /// Transfer each item in a_list into this ListBox if it is contained within a_valid.
  /// This removes the items from a_list.
  virtual void addItemsDestructive(List<T>& a_list, const Box& a_valid);

  /// Performs transfer of the items from src to dest. Notice that as
  /// a result both *this and a_src are modified.
  virtual void transfer(ListBox<T>& a_src,
                        const Box&  a_srcBox,
                        const Box&  a_destBox,
                        const int   a_srcComp =0,
                        const int   a_destComp=0,
                        const int   a_numComp =1);

  /// Remove all the items from this ListBox that are within a_RegionTo, and add
  /// all the particles from the ListBox a_src that are within a_RegionFrom.
  virtual void copy(const Box& a_RegionFrom,
                    const Interval& a_Cdest,
                    const Box& a_RegionTo,
                    const ListBox<T>& a_src,
                    const Interval& a_Csrc);

  /// Return the number of items in this ListBox.
  virtual size_t numItems() const;

  //  Return the number of items in this ListBox that are inside a_box; this one is slow.
  virtual size_t numItems(const Box& a_box, const bool a_in=true) const;

  // Linearization functions

  /// The number of bytes used by \b linearIn()/linearOut().
  /** Returns the size, in number of bytes, of a flat linear
      representation of the items in the specified \p Box.
  */
  virtual int size(const Box&      a_box,
                   const Interval& a_comps = Interval()) const;

  /// Write a serialized (packed) representation into \a a_buf.
  /** Write a linear representation of the items in the specified \p Box.
      Assumes that sufficient memory for the buffer has already been
      allocated by the caller.
  */
  virtual void linearOut(void* a_buf,
                         const Box& a_box,
                         const Interval& a_comps = Interval()) const;

  ///
  virtual void linearOutDestructive(void* buf,
                                    const Box& a_box,
                                    const Interval& a_comps = Interval());

  /// Extract a serialized (packed) representation from \a a_buf.
  /** Given a linear representation of the class data previously made
      using linearIn() in \a a_buf, set the data for this class on the
      specified \p Box.
  */
  virtual void linearIn(void*           a_buf,
                        const Box&      a_box,
                        const Interval& a_comps = Interval());

  /// Retrieve the mesh size.
  virtual const RealVect& meshSpacing() const;

  /// Retrieve the box.
  const Box& box() const;

  /// Delete all items in this ListBox and reset it to the null Box.
  virtual void clear();

  /// Lists are not pre-allocatable. This function is needed so that ListBox 
  /** can be stored inside LevelData<T> containers
  */
  static int preAllocatable() { return 2; }

  /// as in BaseFab; needed for loop-unroll macro
  IntVect size() const
  {
    return this->m_box.size();
  }

protected:

  static Arena* s_Arena;

  static std::string name();

  /// The Box that defines the valid region of this ListBox
  Box m_box;

  /// Pointer to the List that stores our items
  List<T>* m_dptr;

  /// The cell size in physical space
  RealVect m_mesh_spacing;

  /// The center of our coordinate system. Used for binning our items. 
  RealVect m_origin;

  /// Determines whether a_item is enclosed in a_box
  inline virtual bool isEnclosed(const T& a_item, const Box& a_box) const;

private:
  ///
  /// These functions are made private to prevent use of the default
  /// functions provided by the C++ compiler.
  ///
  ListBox<T>& operator= (const ListBox<T>&);
  ListBox(const ListBox<T>&);

};

#include "NamespaceFooter.H"

// Implementation
#include "ListBoxI.H"

#endif
