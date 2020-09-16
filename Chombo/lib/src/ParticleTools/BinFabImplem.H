#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _BINFABIMPLEM_H_
#define _BINFABIMPLEM_H_

#include <cmath>

#include "BoxIterator.H"
#include "parstream.H"
#include "NamespaceHeader.H"

//
// Implementation.=====================================
//

template <class T>
BinFab<T>::BinFab()
  : BaseFab<List<T> >(),
    //[NOTE: dont use RealVect::Zero because of unreliable order of static initializations. <dbs>]
    m_origin(D_DECL6(0,0,0,0,0,0)),
    m_mesh_spacing(D_DECL6(0,0,0,0,0,0))
{
}

template <class T>
BinFab<T>::BinFab(const Box&           a_domain,
                  const RealVect&      a_meshSpacing,
                  const RealVect&      a_origin)
{
  define(a_domain, a_meshSpacing, a_origin);
}

#ifdef __IBMCPP__
//[NOTE: this is here only because the IBM xlC compiler tries to compile
//       DefaultDataFactory<BinFab<BinItem>> in particleTest.cpp, even
//       though it isn't used anywhere.  Stupid compiler. <dbs> July2005]
template <class T>
BinFab<T>::BinFab(const Box& a_domain, int a_ncomps)
{
  MayDay::Abort("BinFab(Box,int): Dont use this constructor");
}
#endif

template <class T>
BinFab<T>::~BinFab()
{
  this->undefine();
}

template <class T>
void BinFab<T>::define(const Box&           a_domain,
                       const RealVect&      a_meshSpacing,
                       const RealVect&      a_origin)
{
  // BaseFab vars
  this->m_domain = a_domain;
  this->m_nvar = 1;
  this->m_numpts = this->m_domain.numPts();
  BaseFab<List<T> >::define();
  // BinFab vars
  m_mesh_spacing = a_meshSpacing;
  m_origin = a_origin;
}

// set the mesh spacing
//[NOTE: this is a vector even though Chombo doesn't support isotropic grids yet.]
template <class T>
void BinFab<T>::meshSpacing(const RealVect & a_meshSpacing)
{
  m_mesh_spacing = a_meshSpacing;
}

// retrieve a copy of the mesh spacing
template <class T>
RealVect BinFab<T>::meshSpacing() const
{
  return m_mesh_spacing;
}

// set the coordinate origin;
// i.e. the physical location of the lower corner of cell IntVect::Zero
template <class T>
void BinFab<T>::origin(const RealVect& a_origin)
{
  m_origin = a_origin;
}

// retrieve a copy of the coordinate origin
template <class T>
RealVect BinFab<T>::origin() const
{
  return m_origin;
}

////////////////////////////////////////////////////////////////

/// Insert a single item into this BinFab in the appropriate cell.
template <class T>
void BinFab<T>::addItem(const T& a_item)
{
  int comp = 0;
  IntVect binLoc = locateBin(a_item);
  if (this->m_domain.contains(binLoc))
  {
    this->operator()(binLoc,comp).append(a_item);
  }
}

/// Insert a single item into this BinFab in the specified cell.
template <class T>
void BinFab<T>::addItem(const T& a_item, const IntVect& a_binLoc)
{
  int comp = 0;
  if (this->m_domain.contains(a_binLoc))
  {
    this->operator()(a_binLoc,comp).append(a_item);
  }
}

/// Insert multiple items into this BinFab.
template <class T>
void BinFab<T>::addItems(const List<T>& a_List)
{
  if (a_List.isNotEmpty())
  {
    // loop through items in List, and add to appropriate locations
    ListIterator<T> lit(a_List);
    int comp = 0;
    IntVect binLoc;

    for (lit.rewind(); lit; ++lit)
    {
      const T& thisItem =  a_List[lit];
      binLoc = locateBin(thisItem);

      if (this->m_domain.contains(binLoc))
      {
        this->operator()(binLoc,comp).append(thisItem);
      }
    }
  }
}

////////////////////////////////////////////////////////////////

/// Insert multiple items into BinFab and remove from source list.
//  [NOTE: this is mostly the same as the previous version, except for the inner loop]
template <class T>
void BinFab<T>::addItemsDestructive(List<T>& a_List)
{
  if (a_List.isNotEmpty())
  {
    // loop through items in List, add to appropriate bins,
    // and remove from List
    ListIterator<T> lit(a_List);
    int comp = 0;
    IntVect binLoc;

    for (lit.rewind(); lit; )
    {
      const T& thisItem =  a_List[lit];
      binLoc = locateBin(thisItem);
      if (this->m_domain.contains(binLoc)) //[NOTE: BaseFab domain, not problem domain]
      {
        //this->operator()(binLoc,comp).append(thisItem);
        //a_List.remove(lit);
        this->operator()(binLoc,comp).transfer(lit);
      }
      else
      {
        // removing a list element has the side-effect of moving the
        // iterator forward, so only increment explicitly if not removed
        ++lit ;
      }
    }
  }
}

/// Move or copy items from arg list to *this, depending on in or not-in the given Box.
// [NOTE: this is mostly the same as the previous version, except for the inner loop]
template <class T>
void BinFab<T>::addItemsDestructive(List<T>& a_List, const Box& a_valid ,bool a_in)
{
  CH_assert( a_valid.cellCentered() );
  if (a_List.isNotEmpty())
  {
    // loop through items in List, add to appropriate bins,
    // and remove from List
    ListIterator<T> lit(a_List);
    int comp = 0;
    IntVect binLoc;

        for (lit.rewind(); lit; )
      {
        const T& thisItem =  a_List[lit];
        binLoc = locateBin(thisItem);
        if (this->m_domain.contains(binLoc)) //[NOTE: BaseFab domain, not problem domain]
          {

            // If the a_in flag is true, then delete the particle from the
            // List if it is in the valid Box; if the a_in flag is false,
            // delete the particle if it is NOT in the Box.
            // Transfering a list element has the side-effect of moving the
            // iterator forward, so explicitly increment the iterator only
            // if the particle was not transfered.  Wish C++ have a nice XOR function.
            if ( a_valid.contains(binLoc) )
              {
                if ( a_in )
                  {
                    //this->operator()(binLoc,comp).append(thisItem);
                    //a_List.remove(lit);
                    this->operator()(binLoc,comp).transfer(lit);
                  }
                else
                  {
                    this->operator()(binLoc,comp).append(thisItem);
                    ++lit ;
                  }
              }
            else
              {
                if ( a_in )
                  {
                    this->operator()(binLoc,comp).append(thisItem);
                    ++lit;
                  }
                else
                  {
                    //this->operator()(binLoc,comp).append(thisItem);
                    //a_List.remove(lit);
                    this->operator()(binLoc,comp).transfer(lit);
                  }
              }

          }
        else
          {
            // not added, so just increment the iterator
            ++lit ;
          }
      }
  }
}

////////////////////////////////////////////////////////////////

// Look at every item in the BinFab and move to a different cell
// if it is in the wrong one.  Useful after computing new physical
// locations of particles, e.g. moving in time.
template <class T>
void BinFab<T>::reBin()
{
  // loop through each bin and determine if each binItem is
  // in the correct bin.
  BoxIterator bit(this->m_domain);
  int comp = 0;
  IntVect binLoc;
//  int lost = 0;

  for (BoxIterator bit(this->m_domain); bit.ok(); ++bit)
    {
      const IntVect thisIV = bit();
      List<T>& thisList = this->operator()(thisIV,comp);

      if (thisList.isNotEmpty())
        {
          // now index through the List and figure out which bin each
          // particle should be in
          for (ListIterator<T> thisLit(thisList); thisLit;)
            {
              T& thisItem = thisList[thisLit];
              binLoc = locateBin(thisItem);

              if (binLoc != thisIV)
                {
                  // binItem needs to be moved.
                  // note that if binItem moves outside of domain, it is lost
                  if (this->m_domain.contains(binLoc))
                    {
                      // this->operator()(binLoc,comp).append(thisItem);
                      // thisList.remove(thisLit);
                      this->operator()(binLoc,comp).transfer(thisLit);
                    }
                  else
                    {
//                    ++lost ;
                      thisList.remove(thisLit);
                    }
                }
              else
                {
                  // if this list element is not moved, increment the iterator
                  ++thisLit;
                }
            } // end loop over items in this List of binItems
        } // end if this list has items
    } // end loop over bins
//  return lost;
}

// Sort all the items into the appropriate bins and return items that are in (or not in) the \p Box.
// The returned items are appended to the list, so the caller should clear() first if apprpriate.
// Same as reBin(void) except it keeps track of the lost items.
//[NOTE: items that are moved forward in the box wil be looked at twice. <dbs>]
template <class T>
void BinFab<T>::reBin(List<T>& a_lost, const Box & a_valid, bool a_in)
{
//XXX What to do if a_valid is not contained in m_domain?? <dbs>
  int comp = 0;  //assume only one component
  IntVect binLoc;
  Box valid = a_valid.isEmpty() ? this->m_domain : a_valid ;

  // Loop through _all_ the bins in this BinFab, regardless of the valid Box
  for (BoxIterator bit(this->m_domain); bit.ok(); ++bit)
  {
    const IntVect thisIV = bit();
    List<T>& thisList = this->operator()(thisIV,comp);
    if (thisList.isNotEmpty())
    {
      // Compute the proper cell for each element in this list
      for (ListIterator<T> thisLit(thisList); thisLit;) //[NOTE: no auto-increment operation]
      {
        T& thisItem = thisList[thisLit];
        binLoc = locateBin(thisItem);
        if (binLoc != thisIV)
        {
          // move it or return it, as appropriate
          if (valid.contains(binLoc))
          { // new location is in the arg box
            if (a_in) a_lost.transfer(thisLit) ;//return it
            else      this->operator()(binLoc,comp).transfer(thisLit); //move it
          }
          else if (this->m_domain.contains(binLoc))
          { // new location not in the arg box but in the domain
            this->operator()(binLoc,comp).transfer(thisLit); //move it
          }
          else
          { // new location is not in the domain
            if (!a_in) a_lost.transfer(thisLit); //return it
            else  thisList.remove(thisLit); // and not in the box, so throw it away
          }
          //
        }
        else
        {
          // increment the iterator only if the item stays in place
          ++thisLit;
        }
      } // end loop over items in this List of binItems
    } // end if this list has items
  } // end loop over bins
}

////////////////////////////////////////////////////////////////

// Move items from a_srcBox in the source BinFab to a_destBox in *this.
// This is what ::copy() would do if it allowed a non-const source.
//[NOTE: if the boxes are different, the items will end up in cells with
//       different IntVects from the ones they started in. <dbs>]
template <class T>
void BinFab<T>::transfer(BinFab<T>&  a_src,
                    const Box&  a_srcBox,
                    const Box&  a_destBox)
{
  CH_assert(a_src.box().contains(a_srcBox));
  CH_assert(a_destBox.sameSize(a_srcBox));
  // Ignore components; BinFab assumes only one.
  CH_assert( a_src.nComp() == 1 );
  const int comp = 0, ncomps = 1;
  //XXX stupid lack of operator==
  //XXX CH_assert( a_comps == Interval(0,0) );
  // loop through the cells in the intersection of the two boxes
  ForAllThisBNNXCBN(List<T>,a_destBox,comp,ncomps,a_src,a_srcBox,comp)
    {
      if (a_srcR.isNotEmpty())
      {
        // cast away const-ness that is assumed by the ForAll...
        List<T>& src = (List<T>&)a_srcR ;
//XXX debug
//XXX        pout() << "dbg: transfer: i=" << iR << " has " << a_srcR.length() << endl ;
        thisR.catenate(src);
      }
    } EndForTX
}

////////////////////////////////////////////////////////////////

// Deletes all items from this BinFab and resets it to the undefined state.
template <class T>
void BinFab<T>::clear()
{
  this->undefine();
  this->m_domain = Box();
  this->m_nvar = 0;
  this->m_numpts = 0;
  m_origin = RealVect::Zero;
  m_mesh_spacing = RealVect::Zero;
}

////////////////////////////////////////////////////////////////

// Counts the number of items in the cells of a_box.
template <class T>
int BinFab<T>::numItems(const Box& a_box) const
{
  // Ignore components; BinFab assumes only one.
  const int comp = 0 ;
  // default is an empty box, so use this BinFab's box instead
  // [NOTE: this leaves the caller no way to ask for an empty box;
  //       it must avoid calling this routine in that case. Sorry.]
  Box ibox = a_box;
  if (ibox.isEmpty())
  {
    ibox = this->box();
  }
  // Loop over all bins in the BinFab
  int totalSize = 0;
  for (BoxIterator bi(ibox); bi.ok(); ++bi)
  {
    totalSize += this->operator()(bi(),comp).length();
  }
//XXX -- debug
//XXX  pout() << "dbg: BinFab::numItems is " << totalSize << " for box " << ibox << endl ;
  return totalSize;
}

////////////////////////////////////////////////////////////////

// Computes the linear size (in bytes) of the items in the cells of a_box.
template <class T>
size_t BinFab<T>::size(const Box& a_box, const Interval& a_comps) const
{
  // Ignore components; BinFab assumes only one.
  const int ncomps = 1 ;
  //XXX stupid lack of operator==
  //XXX CH_assert( a_comps == Interval(0,0) );
  // Count the number of items in the specified Box.
  size_t totalSize = numItems(a_box);
  // Compute size as the number of items times the size of the items
  {
    // Compute the size of an item, assuming all are the same size.
    int sizeOfT;
    T tmp;
    sizeOfT = tmp.size();
    // totalSize has total number of items, so now compute total size
    totalSize *= sizeOfT;
  }
  // Will also need to include an integer for each list to specify how
  // many items are in each list.
  //[NOTE: Is this the right thing to do?  An alternate approach would be to
  //       simply unpack one List<T> for all the items and then call addItems()
  //       with that list.  But that probably will be slower. <dfmartin>]
  size_t numBins = a_box.numPts() * ncomps ;
  totalSize += numBins*sizeof(int);
  return totalSize;
}

////////////////////////////////////////////////////////////////

// Linearizes the items in the cells of a_box into a_buf.
template <class T>
void BinFab<T>::linearOut(void* a_buf, const Box& a_box, const Interval& a_comps) const
{
  // Ignore components; BinFab assumes only one.
  const int comp = 0 ;
  //XXX stupid lack of operator==
  //XXX CH_assert( a_comps == Interval(0,0) );
  // make a usable pointer to the buffer
  char* buf_ch = (char*)a_buf;
  // All we do in this function is loop over the Box
  // and call linearOut() on the List elements.
  for (BoxIterator bit(a_box); bit.ok(); ++bit)
  {
    const List<T>& thisList = this->operator()(bit(), comp);
    // write the number of list items first
    int *intBuffer = (int*)buf_ch;
    *intBuffer = thisList.length();
    buf_ch += sizeof(int);
    // now loop over the items in the list and write each one
    for (ListIterator<T> lit(thisList); lit.ok(); ++lit)
    {
      thisList[lit].linearOut(buf_ch);
      buf_ch += thisList[lit].size();
    }
  }
}

////////////////////////////////////////////////////////////////

// Same as linearOut() except it clear()s the List after writing it.
template <class T>
void
BinFab<T>::linearOutDestructive(void* a_buf, const Box& a_box, const Interval& a_comps)
{
  const int comp = 0 ;
  //XXX stupid lack of operator==
  //XXX CH_assert( a_comps == Interval(0,0) );
  char* buf_ch = (char*)a_buf;
  for (BoxIterator bit(a_box); bit.ok(); ++bit)
  {
    List<T>& thisList = this->operator()(bit(), comp);
    int *intBuffer = (int*)buf_ch;
    *intBuffer = thisList.length();
    buf_ch += sizeof(int);
    for (ListIterator<T> lit(thisList); lit.ok(); ++lit)
    {
      thisList[lit].linearOut(buf_ch);
      buf_ch += thisList[lit].size();
    }
    thisList.clear();
  }
}

////////////////////////////////////////////////////////////////

// Fills this BinFab from the linearization in a_buf of the cells in a_box.
template <class T>
void
BinFab<T>::linearIn(void* a_buf, const Box& a_box, const Interval& a_comps)
{
  // Should be just the inverse of linearOut().
  // Ignore components; BinFab assumes only one.
  const int comp = 0 ;
  //XXX stupid lack of operator==
  //XXX CH_assert( a_comps == Interval(0,0) );
  // make a usable pointer to the buffer
  char *buf_ch = (char*)a_buf;
  // All we do in this function is loop over the Box
  // and call linearOut() on the List elements.
  for (BoxIterator bit(a_box); bit.ok(); ++bit)
  {
    List<T>& thisList = this->operator()(bit(), comp);
    // clear the list to prevent accumulation of binItems
    thisList.clear();

    // first extract the length of the list
    int numItems = *((int*)buf_ch);
    buf_ch += sizeof(int);
    // loop over the list items and make the list
    for (int n=0; n<numItems; ++n)
    {
      T thisBinItem;
      thisBinItem.linearIn( buf_ch );
      thisList.append(thisBinItem);
      buf_ch += thisBinItem.size();
    }
  }
}

//X /*******************************************************************************
//X This is an example that I want to keep here to remind that the loops over boxes
//X should be unrolled through the BaseFabMacro.  <fm>
//X *******************************************************************************
//X template <class T>
//X inline
//X void
//X BaseFab<T>::linearIn(void* buf, const Box& R, const Interval& comps)
//X {
//X   T* buffer = (T*)buf;
//X   ForAllThisBNN(T,R,comps.begin(), comps.size())
//X     {
//X       thisR = *buffer;
//X       ++buffer;
//X     } EndFor;
//X }
//X */

////////////////////////////////////////////////////////////////

// Compute the IntVect of the cell containing a_binItem.
template <class T>
inline
IntVect
BinFab<T>::locateBin(const T& a_binItem) const
{
  IntVect binLoc;
  RealVect thisPos = a_binItem.position();
  thisPos -= m_origin;
  thisPos /= m_mesh_spacing;

  for ( int d=0 ; d<SpaceDim ; ++d )
  {
    binLoc[d] = (int)floor(thisPos[d]);
  }
  return binLoc;
}

#include "NamespaceFooter.H"
#endif
