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

#ifndef _LISTBOXI_H_
#define _LISTBOXI_H_

#include "Arena.H"
#include "BoxIterator.H"
#include <iostream>
#include "parstream.H"
#include "NamespaceHeader.H"

//
// Implementation
//

template <class T>
ListBox<T>::ListBox()
  :
  m_box(Box()),
  m_dptr(0),
  m_mesh_spacing(D_DECL6(0, 0, 0, 0, 0, 0)),
  m_origin(D_DECL6(0, 0, 0, 0, 0, 0))
{}

template <class T>
ListBox<T>::ListBox(const Box&       a_box,
                    const RealVect&  a_meshSpacing,
                    const RealVect&  a_origin)
  :
  m_dptr(0)
{
  define(a_box, a_meshSpacing, a_origin);
}

template <class T> 
Arena* ListBox<T>::s_Arena = NULL;

template <class T>
ListBox<T>::~ListBox()
{
  undefine();
}

template <class T>
void ListBox<T>::define(const Box&  a_box,
                        const RealVect& a_meshSpacing,
                        const RealVect& a_origin)
{

  CH_assert(m_dptr == 0);

#ifdef CH_USE_MEMORY_TRACKING
  if (this->s_Arena == NULL)
    {
      this->s_Arena = new BArena(this->name().c_str());
    }
#else
  if (this->s_Arena == NULL)
    {
      this->s_Arena = new BArena("");
    }
#endif

  if (this->s_Arena == NULL)
    {
      MayDay::Error("Malloc in ListBox failed");
    }

  // allocate our list
  this->m_dptr     =
    static_cast<List<T>*>(this->s_Arena->alloc(sizeof(List<T>)));

#ifdef CH_USE_MEMORY_TRACKING
  this->s_Arena->bytes += sizeof(List<T>) + sizeof(ListBox<T>);
  if (this->s_Arena->bytes > this->s_Arena->peak)
    {
      this->s_Arena->peak = this->s_Arena->bytes;
    }
#endif

  new (m_dptr) List<T>;

  m_box = a_box;
  m_mesh_spacing = a_meshSpacing;
  m_origin = a_origin;
}

template <class T>
void ListBox<T>::undefine()
{

  if (m_dptr == 0)
    {
      return;
    }

    m_dptr->~List<T>();

  s_Arena->free(m_dptr);

#ifdef CH_USE_MEMORY_TRACKING
  s_Arena->bytes -= sizeof(List<T>) + sizeof(ListBox<T>);
  CH_assert(s_Arena->bytes >= 0);
#endif

  m_dptr = 0;
}

// retrieve list of items
template <class T>
List<T>& ListBox<T>::listItems()
{
  CH_assert(!(this->m_dptr == NULL));

  return *(static_cast<List<T>*>(this->m_dptr));
}

template <class T>
const List<T>& ListBox<T>::listItems() const
{
  CH_assert(!(this->m_dptr == NULL));

  return *(static_cast<List<T>*>(this->m_dptr));
}

template <class T>
void ListBox<T>::getItems(List<T>& a_list, const Box& a_valid) const
{
  CH_assert(!(this->m_dptr == NULL));

  if (a_valid == m_box)
    {
      a_list.join( listItems() );
    }
  else
    {
      for (ListIterator<T> li( listItems() ); li; ++li)
        {
          // loop through items
          if ( isEnclosed(li(), a_valid) )
            {
              a_list.add(*li);
            }
        }
    }
}

template <class T>
void ListBox<T>::getItemsDestructive(List<T>& a_list, const Box& a_valid)
{
  CH_assert(!(this->m_dptr == NULL));

  if (a_valid == m_box)
    {
      a_list.catenate( listItems() );
    }
  else
    {
      // loop through items, and add to rlist
      for (ListIterator<T> li( listItems() ); li; )
        {
          // const RealVect& x =  items[li].position();
          if ( isEnclosed(li(),a_valid) )
            {
              a_list.transfer(li);
            }
          else
            {
              // removing a list element has the side-effect of moving
              // the iterator forward, so only increment explicitly if
              // not removed
              ++li;
            }
        }
    }
}

template <class T>
void ListBox<T>::getInvalid(List<T>& a_list, const Box& a_valid) const
{
  CH_assert(!(this->m_dptr == NULL));

  for (ListIterator<T> li( listItems() ); li; ++li)
    {
      // loop through items
      if ( !isEnclosed(li(), a_valid) )
        {
          a_list.add(*li);
        }
    }  
}

template <class T>
void ListBox<T>::getInvalidDestructive(List<T>& a_list, const Box& a_valid)
{
  CH_assert(!(this->m_dptr == NULL));

  // loop through items, and add to rlist
  for (ListIterator<T> li( listItems() ); li; )
    {
      if ( !isEnclosed(li(),a_valid) )
        {
          a_list.transfer(li);
        }
      else
        {
          // removing a list element has the side-effect of moving
          // the iterator forward, so only increment explicitly if
          // not removed
          ++li;
        }
    }    
}

template <class T>
void ListBox<T>::addItem(const T& a_item)
{
  if ( isEnclosed(a_item, m_box) )
    {
      listItems().add(a_item);
    }
}

template <class T>
void ListBox<T>::addItems(const List<T>& a_list)
{
  if (a_list.isNotEmpty())
    {
      List<T>& litems = listItems();

      for (ListIterator<T> li(a_list); li.ok(); ++li)
        {
          if ( isEnclosed(li(), m_box) )
            {
              litems.add(*li);
            }
        }
    }
}

template <class T>
void ListBox<T>::addItemsDestructive(List<T>& a_list)
{
  addItemsDestructive(a_list, m_box);
}

template <class T>
void ListBox<T>::addItemsDestructive(List<T>& a_list, const Box& a_valid)
{
  if (a_list.isNotEmpty())
    {
      List<T>& litems = listItems();

      for (ListIterator<T> li(a_list); li; )
        {
          if ( isEnclosed(li(), a_valid) )
            {
              litems.transfer(li);
            }
          else
            {
              ++li;
            }
        }
    }
}

template <class T>
void ListBox<T>::transfer(ListBox<T>& a_src,
                          const Box&  a_srcBox,
                          const Box&  a_destBox,
                          const int   a_srcComp,
                          const int   a_destComp,
                          const int   a_numComp)
{

  // avoid self-transfering
  if (this->m_dptr != a_src.m_dptr)
    {
      CH_assert(a_destBox.sameSize(a_srcBox));
      CH_assert(this->m_box.contains(a_destBox));

      // when transfering particles don't check destination box because that
      // can't always be right before applying BCs
      a_src.getItemsDestructive(listItems(), a_srcBox);
    }
}

template <class T>
void ListBox<T>::copy(const Box& a_RegionFrom,
                      const Interval& a_Cdest,
                      const Box& a_RegionTo,
                      const ListBox<T>& a_src,
                      const Interval& a_Csrc)
{
  // remove the items from this ListBox that are within a_RegionTo
  getItemsDestructive(listItems(), a_RegionTo);
  List<T> tmpList;
  a_src.getItems(tmpList, a_RegionFrom);
  addItemsDestructive(tmpList, a_RegionTo);
}

template <class T>
size_t ListBox<T>::numItems() const
{
  return listItems().length();
}

template <class T>
size_t ListBox<T>::numItems(const Box& a_box, const bool a_in) const
{

  if (a_in)
    {
      if (a_box.isEmpty())
        {
            return 0;
        }
      else
        {
          size_t nitems = 0;
          for (ListIterator<T> li( listItems() ); li.ok(); ++li)
            {
              if ( isEnclosed(li(), a_box) ) ++nitems;
            }
          return nitems;
        }
    }
  else
    {
      if (a_box.isEmpty())
        {
          return numItems();
        }
      else
        {
          size_t nitems = 0;
          for (ListIterator<T> li( listItems() ); li.ok(); ++li)
            {
              if ( !isEnclosed(li(), a_box) ) ++nitems;
            }
          return nitems;
        }
    }
}

// The number of bytes used by linearIn()/linearOut().
template <class T>
int ListBox<T>::size(const Box&      a_box,
                     const Interval& a_comps) const
{
  // total size = num_items * size_item + integer
  // the integer is to determine how many items to transfer

  const int totalSize = numItems(a_box) * T().size() + sizeof(size_t);

  return totalSize;
}

template <class T> void
ListBox<T>::linearOut(void*           a_buf,
                      const Box&      a_box,
                      const Interval& a_comps) const
{
  // all we do in this function is loop over the list (ignore components),
  // and call linearOut on the individual items...

  // size of template object
  const size_t tsize = T().size();

  // cast buffer to a char*
  char* ch_buf = (char*)a_buf;

  // call functions with box argument
  List<T> items;
  getItems(items,a_box);

  *(size_t*)ch_buf = items.length();
  ch_buf += sizeof(size_t);

  for (ListIterator<T> li(items); li; ++li)
    {
      li().linearOut(ch_buf);
      ch_buf += tsize;
    }
}

template <class T> void
ListBox<T>::linearOutDestructive(void*           a_buf,
                                 const Box&      a_box,
                                 const Interval& a_comps)
{
  // all we do in this function is loop over the box (ignore components),
  // and call linearOut on the individual items...

  // size of template object
  const size_t tsize = T().size();

  // cast buffer to a char*
  char* ch_buf = (char*)a_buf;

  // call functions with box argument
  List<T> items;
  getItemsDestructive(items,a_box);

  *(size_t*)ch_buf = items.length();
  ch_buf += sizeof(size_t);

  for (ListIterator<T> li(items); li; ++li)
    {
      li().linearOut(ch_buf);
      ch_buf += tsize;
    }
}

template <class T>
void ListBox<T>::linearIn(void*           a_buf,
                          const Box&      a_box,
                          const Interval& a_comps)
{
  // size of template object
  const size_t tsize = T().size();

  // should be just the inverse of linearOut
  char* ch_buf = (char*)a_buf;

  List<T>& items = listItems();

  size_t numItems = *((size_t*)ch_buf);
  ch_buf  += sizeof(size_t);

  for (size_t i=0; i<numItems; i++)
    {
      T item;
      item.linearIn( ch_buf );

      // assume care has been taken to ensure the item is inside
      // m_box. Check this below when debugging.
      items.append(item);

      ch_buf += tsize;
    }
}

template <class T>
const RealVect& ListBox<T>::meshSpacing() const
{
  return m_mesh_spacing;
}

template <class T>
const Box& ListBox<T>::box() const
{
  return m_box;
}

template <class T>
void ListBox<T>::clear()
{
  this->undefine();

  m_box = Box();
  m_mesh_spacing = RealVect::Zero;
  m_origin = RealVect::Zero;
}

template <class T> 
inline std::string ListBox<T>::name()
{
  std::string rtn = (typeid(T)).name();

  return rtn;
}

template <class T>
bool ListBox<T>::isEnclosed(const T& a_item, const Box& a_box) const
{
  // corners of valid region
  const IntVect& loEnd = a_box.smallEnd();
  const IntVect& hiEnd = a_box.bigEnd();

  // item position
  const RealVect& x = a_item.position() - m_origin;
  const RealVect& dx= m_mesh_spacing;

  bool is_in = D_TERM6( (x[0]>=loEnd[0]*dx[0] && x[0]<(hiEnd[0]+1)*dx[0]),
			&& (x[1]>=loEnd[1]*dx[1] && x[1]<(hiEnd[1]+1)*dx[1]),
			&& (x[2]>=loEnd[2]*dx[2] && x[2]<(hiEnd[2]+1)*dx[2]),
			&& (x[3]>=loEnd[3]*dx[3] && x[3]<(hiEnd[3]+1)*dx[3]),
			&& (x[4]>=loEnd[4]*dx[4] && x[4]<(hiEnd[4]+1)*dx[4]),
			&& (x[5]>=loEnd[5]*dx[5] && x[5]<(hiEnd[5]+1)*dx[5]));

  return is_in;
}

#include "NamespaceFooter.H"

#endif
