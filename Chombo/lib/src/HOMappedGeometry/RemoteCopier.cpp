#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "RemoteCopier.H"
#include "LayoutIterator.H"
#include "DataIterator.H"
#include "NeighborIterator.H"
#include <iostream>
#include "NamespaceHeader.H"

using namespace std;

//-----------------------------------------------------------------------
RemoteCopier::
RemoteCopier():
  Copier()
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
RemoteCopier::
RemoteCopier(const DisjointBoxLayout& a_level,
             const BoxLayout& a_dest,
             const IntVect& a_ghostVect,
             bool a_exchange):
  Copier()
{
  define(a_level, a_dest, a_ghostVect, a_exchange);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
RemoteCopier::
~RemoteCopier()
{
  clear();
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
RemoteCopier::
define(const DisjointBoxLayout& a_level,
       const BoxLayout& a_dest,
       const IntVect& a_ghostVect,
       bool  a_exchange)
{
  CH_assert(a_level.isClosed());
  CH_assert(a_dest.isClosed());

  clear();
  m_isDefined = true;
  buffersAllocated = false;
  const BoxLayout&         level= a_level;
  const BoxLayout&         dest = a_dest;

  unsigned int myprocID = procID();

#ifdef CH_MPI  // don't need to do this in serial
  // make a vector of boxes (or LayoutIndexes to boxes) from destination layout
  // that are known to reside on this processor.
  vector<DataIndex> vectorDestDI;
  vector<DataIndex> vectorDestOnProcDI;
  for (LayoutIterator to(a_dest.layoutIterator()); to.ok(); ++to)
  {
    vectorDestDI.push_back(DataIndex(to()));
    if (myprocID == dest.procID(to()))
    {
      vectorDestOnProcDI.push_back(DataIndex(to()));
    }
  }

  // make a vector of boxes (or LayoutIndexes to boxes) from "level"/src layout
  // that are known to reside on this processor.
  vector<DataIndex> vectorLevelDI;
  vector<DataIndex> vectorLevelOnProcDI;
  for (LayoutIterator from(a_level.layoutIterator()); from.ok(); ++from)
  {
    vectorLevelDI.push_back(DataIndex(from()));
    if (myprocID == level.procID(from()))
    {
      vectorLevelOnProcDI.push_back(DataIndex(from()));
    }
  }
#else
  // in serial, it's not very interesting as it's all of them.
  vector<DataIndex> vectorDestOnProcDI;
  vector<DataIndex> vectorLevelDI;
  for (LayoutIterator to(a_dest.layoutIterator()); to.ok(); ++to)
  {
    vectorDestOnProcDI.push_back(DataIndex(to()));
  }
  for (LayoutIterator from(a_level.layoutIterator()); from.ok(); ++from)
  {
    vectorLevelDI.push_back(DataIndex(from()));
  }
#endif

  // loop over all dest/to DI's on my processor
  for (vector<DataIndex>::iterator vdi=vectorDestOnProcDI.begin();
      vdi != vectorDestOnProcDI.end(); ++vdi)
  {
    // at this point, i know myprocID == toProcID
    const DataIndex todi(*vdi);
    const Box t = dest[todi];
    Box toBox = t;
    toBox.grow(a_ghostVect); // Grow by ghost vector

    //bool isSorted = (a_level.isSorted() && a_dest.isSorted());
    // then for each level/from DI, see if they intersect
    for (vector<DataIndex>::iterator vli = vectorLevelDI.begin();
        vli != vectorLevelDI.end(); ++vli)
    {

      const DataIndex fromdi(*vli);
      const unsigned int fromProcID = level.procID(fromdi);
//      const Box& fromBox = level[fromdi];
      Box fromBox = level[fromdi];
      fromBox.grow(a_ghostVect);
      MotionItem* item = new (s_motionItemPool.getPtr()) MotionItem(fromdi, todi, fromBox, toBox);
      if (item == NULL)
      {
        MayDay::Error("Out of Memory in copier::define");
      }
      if (fromProcID == myprocID)
      {
        // local move
        if (a_exchange && fromdi == todi)
          s_motionItemPool.returnPtr(item);
        else
          m_localMotionPlan.push_back(item);
      }
      else
      {
        item->procID = fromProcID;
        m_toMotionPlan.push_back(item);
      }
    }
  }

  // Don't need to worry about this in serial as we already
  // took care of the local copy motion items just above.  skip this.
#ifdef CH_MPI
  // loop over all dest/to DI's
  for (vector<DataIndex>::iterator vdi=vectorDestDI.begin();
      vdi != vectorDestDI.end(); ++vdi)
  {

    const DataIndex todi(*vdi);
    const Box t = dest[todi];
    Box toBox = t;
    toBox.grow(a_ghostVect); // Grow by ghost vector.
    const unsigned int toProcID = dest.procID(todi);

    // then for each level/from DI on this processor, see if they intersect
    for (vector<DataIndex>::iterator vli = vectorLevelOnProcDI.begin();
        vli != vectorLevelOnProcDI.end(); ++vli)
    {
      // at this point, i know myprocID == fromProcID

      const DataIndex fromdi(*vli);
      Box fromBox = level[fromdi];
      fromBox.grow(a_ghostVect);
      if (toProcID == myprocID)
      { // local move
        // don't push back here!  or you will get two.
        //     we already did it above...
        //m_localMotionPlan.push_back(item);
      }
      else
      {
        MotionItem* item = new (s_motionItemPool.getPtr())
          MotionItem(fromdi, todi, toBox, fromBox);
        if (item == NULL)
        {
          MayDay::Error("Out of Memory in copier::define");
        }

        item->procID = toProcID;
        m_fromMotionPlan.push_back(item);
      }
    }
  }
#endif

  sort();
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
RemoteCopier::
define(const DisjointBoxLayout& a_from,
       const DisjointBoxLayout& a_to,
       const DisjointBoxLayout& a_sourceBoxes,
       const BoxLayout& a_destBoxes,
       const map<Box, Box>& a_sourceToDestMapping,
       const IntVect& a_ghostVect,
       bool a_exchange)
{
  CH_assert(a_from.isClosed());
  CH_assert(a_to.isClosed());
  CH_assert(a_sourceBoxes.isClosed());
  CH_assert(a_destBoxes.isClosed());

  clear();
  m_isDefined = true;
  buffersAllocated = false;
  //  const BoxLayout& from = a_from;
  //  const BoxLayout& to = a_to;

  unsigned int myprocID = procID();

#ifdef CH_MPI  // don't need to do this in serial

  const BoxLayout& src  = a_sourceBoxes;
  const BoxLayout& dest = a_destBoxes;

  // make a mapping of boxes to their data indices within the layouts.
  map<Box, DataIndex> destDI, localDestDI,
    fromDI, sourceDI, localFromDI, localSourceDI,
    localToDI, toDI;
  for (LayoutIterator it(a_destBoxes.layoutIterator()); it.ok(); ++it)
    {
      Box box = a_destBoxes[it()];
      destDI[box] = DataIndex(it());
      if (myprocID == dest.procID(it()))
        localDestDI[box] = DataIndex(it());
    }
  for (LayoutIterator it(a_from.layoutIterator()); it.ok(); ++it)
    {
      Box box = a_from[it()];
      fromDI[box] = DataIndex(it());
      if (myprocID == a_from.procID(it()))
        localFromDI[box] = DataIndex(it());
    }
  for (LayoutIterator it(a_to.layoutIterator()); it.ok(); ++it)
    {
      Box box = a_to[it()];
      toDI[box] = DataIndex(it());
      if (myprocID == a_to.procID(it()))
        localToDI[box] = DataIndex(it());
    }
  for (LayoutIterator it(a_sourceBoxes.layoutIterator()); it.ok(); ++it)
    {
      Box box = a_sourceBoxes[it()];
      sourceDI[box] = DataIndex(it());
      if (myprocID == src.procID(it()))
        localSourceDI[box] = DataIndex(it());
    }
#else
  // in serial, it's not very interesting as it's all of them.
  map<Box, DataIndex> localDestDI, fromDI, toDI, sourceDI;
  for (LayoutIterator it(a_destBoxes.layoutIterator()); it.ok(); ++it)
    {
      Box box = a_destBoxes[it()];
      localDestDI[box] = DataIndex(it());
    }
  for (LayoutIterator it(a_from.layoutIterator()); it.ok(); ++it)
    {
      Box box = a_from[it()];
      fromDI[box] = DataIndex(it());
    }
  for (LayoutIterator it(a_to.layoutIterator()); it.ok(); ++it)
    {
      Box box = a_to[it()];
      toDI[box] = DataIndex(it());
    }
  for (LayoutIterator it(a_sourceBoxes.layoutIterator()); it.ok(); ++it)
    {
      Box box = a_sourceBoxes[it()];
      sourceDI[box] = DataIndex(it());
    }
#endif

  // Build a mapping between the "from" boxes and the "source" boxes.
  // To each Box in a_sourceBoxes there corresponds a Box in a_from:
  // sourceToFromMapping[sourceBox] = fromBox
  // such that sourceBox is contained in grow(fromBox, 1).
  map<Box, Box> sourceToFromMapping;
  for (LayoutIterator l(a_from.layoutIterator()); l.ok(); ++l)
    {
      Box fromBox = a_from[l()];
      Box fromBox1 = fromBox;
      fromBox1.grow(IntVect::Unit); // Grow ghost cells.
      // Go through a_sourceBoxes to find all Boxes contained in fromBox1.
      for (LayoutIterator s(a_sourceBoxes.layoutIterator()); s.ok(); ++s)
        {
          Box sourceBox = a_sourceBoxes[s()];
          if (fromBox1.contains(sourceBox))
            {
              sourceToFromMapping[sourceBox] = fromBox;
            }
        }
    }

  // Build a mapping between the "to" boxes and the "dest" boxes.
  // To each Box in a_destBoxes there corresponds a Box in a_to:
  // destToToMapping[destBox] = toBox
  // such that destBox is contained in grow(toBox, 1).
  map<Box, Box> destToToMapping;
  for (LayoutIterator l(a_to.layoutIterator()); l.ok(); ++l)
    {
      Box toBox = a_to[l()];
      Box toBox1 = toBox;
      toBox1.grow(IntVect::Unit); // Grow ghost cells.
      // Go through a_destBoxes to find all Boxes contained in toBox1.
      for (LayoutIterator d(a_destBoxes.layoutIterator()); d.ok(); ++d)
        {
          Box destBox = a_destBoxes[d()];
          if (toBox1.contains(destBox))
            {
              destToToMapping[destBox] = toBox;
            }
        }
    }

  // Set up all the motion items for destinations on this processor.
#ifdef CH_MPI
  // changed by petermc, 11 Sep 2012: localSourceDI, not sourceDI.
  // CH_assert(localDestDI.size() == localSourceDI.size());
  CH_assert(destDI.size() == sourceDI.size());
#else
  CH_assert(localDestDI.size() == sourceDI.size());
#endif
  for (map<Box, Box>::const_iterator s2d = a_sourceToDestMapping.begin();
       s2d != a_sourceToDestMapping.end(); ++s2d)
    {
      // Get at all the boxes.
      Box sourceBox = s2d->first;
      Box destBox = s2d->second;
      Box fromBox = sourceToFromMapping[sourceBox];
      Box toBox = destToToMapping[destBox];
      
      // Make sure that the destination is local to this processor.
      // If destIter is NOT in localDestDI, then go to next iteration.
      map<Box, DataIndex>::const_iterator destIter = localDestDI.find(destBox);
      if (destIter == localDestDI.end()) continue;
      // Now check this.
#ifdef CH_MPI
      const DataIndex destdi = destDI[destBox];
      const unsigned int destProcID = dest.procID(destdi);
      CH_assert(destProcID == myprocID);
#endif
      
      // Retrieve the corresponding data indices.
      //    const DataIndex todi = destIter->second;
      const DataIndex todi = toDI[toBox];
      const DataIndex fromdi = fromDI[fromBox];
      // const DataIndex srcdi = sourceDI[sourceBox];

#ifdef CH_MPI
      CH_assert(a_to.procID(todi) == destProcID);
      const DataIndex srcdi = sourceDI[sourceBox];
      CH_assert(a_from.procID(fromdi) == src.procID(srcdi));
#endif
      
      const unsigned int fromProcID = a_from.procID(fromdi);

      // Grow the 'source' and 'dest' boxes by a ghost vector.
      sourceBox.grow(a_ghostVect);
      destBox.grow(a_ghostVect);
      MotionItem* item = new (s_motionItemPool.getPtr())
        MotionItem(fromdi, todi, sourceBox, destBox);
      if (item == NULL)
        {
          MayDay::Error("Out of Memory in copier::define");
        }
      if (fromProcID == myprocID)
        { // In this case, fromProcID == destProcID == myprocID.
          // local move
          if (a_exchange && fromdi == todi)
            s_motionItemPool.returnPtr(item);
          else
            m_localMotionPlan.push_back(item);
        }
      else
        { // non-local move from fromProcID to destProcID == myprocID
          item->procID = fromProcID;
          m_toMotionPlan.push_back(item);
        }
    }
  
  // Set up all the motion items for destinations on other processors.
  // No need to do this in the serial case.
#ifdef CH_MPI
  for (map<Box, Box>::const_iterator s2d = a_sourceToDestMapping.begin();
       s2d != a_sourceToDestMapping.end(); ++s2d)
    {
      // Get at all the boxes.
      Box sourceBox = s2d->first;
      Box destBox = s2d->second;
      // const Box& fromBox = sourceToFromMapping[sourceBox];
      Box fromBox = sourceToFromMapping[sourceBox];
      Box toBox = destToToMapping[destBox];

      // Make sure that the destination is on another processor.
      // If destIter IS in localDestDI, then go to next iteration.
      map<Box, DataIndex>::const_iterator destIter = localDestDI.find(destBox);
      if (destIter != localDestDI.end()) continue;
      // Now check this.
      const DataIndex destdi = destDI[destBox];
      const unsigned int destProcID = dest.procID(destdi);
      CH_assert(destProcID != myprocID);

      // Retrieve the corresponding data indices.
      // const DataIndex todi = toDI[destBox];
      const DataIndex todi = toDI[toBox];
      const DataIndex fromdi = fromDI[fromBox];
      const DataIndex srcdi = sourceDI[sourceBox];

      const unsigned int fromProcID = a_from.procID(fromdi);
      const unsigned int toProcID = a_to.procID(todi);
      CH_assert(toProcID == destProcID);

      const unsigned int srcProcID = src.procID(srcdi);
      CH_assert(fromProcID == srcProcID);

      // petermc, 12 Sep 2012:
      // The loop above dealt with all destinations on this proc.
      // The loop here deals with all sources on this proc other than those
      // with destinations on this proc.
      if (fromProcID == myprocID)
        {
          // If (toProcID == myprocID) then we have already done this.
          if (toProcID != myprocID)
            {
              // Grow the 'source' and 'dest' boxes by a ghost vector.
              sourceBox.grow(a_ghostVect);
              destBox.grow(a_ghostVect);
              MotionItem* item = new (s_motionItemPool.getPtr())
                MotionItem(fromdi, todi, sourceBox, destBox);
              if (item == NULL)
                {
                  MayDay::Error("Out of Memory in copier::define");
                }
              item->procID = toProcID;
              m_fromMotionPlan.push_back(item);
            }
        }
    }
#endif

  sort();
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"
