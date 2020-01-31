#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "DataIterator.H"
#include "IntVect.H"
#include "MultiBlockCopier.H"
#include "MayDay.H"
#include "LayoutIterator.H"
#include "SPMD.H"
#include "CH_Timer.H"

#include <vector>
#include "NamespaceHeader.H"

using std::ostream;


// ---------------------------------------------------------
MultiBlockCopier::MultiBlockCopier(const DisjointBoxLayout& a_level,
                                   const BoxLayout& a_dest,
                                   bool a_exchange)
{
  // don't call define here, since we know that it shouldn't do anything
  // since we don't have any block domains
  MayDay::Warning("Called MultiBlockCopier::define with no block domains!");
  return;
}


// ---------------------------------------------------------
MultiBlockCopier::MultiBlockCopier(const DisjointBoxLayout& a_level,
                                   const BoxLayout& a_dest,
                                   const ProblemDomain& a_domain,
                                   bool a_exchange)
{
  // call define here, since we know that it shouldn't do anything
  // since we don't have any block domains
  define(a_level, a_dest, a_domain, IntVect::Zero, a_exchange);
}


// ---------------------------------------------------------
MultiBlockCopier::MultiBlockCopier(const DisjointBoxLayout& a_level,
                                   const BoxLayout& a_dest,
                                   const IntVect& a_destGhost,
                                   bool a_exchange)
{
  // no block domains, so disallow this function
  MayDay::Error("Can't initialize MultiBlockCopier without block domains");
}


// ---------------------------------------------------------
MultiBlockCopier::MultiBlockCopier(const DisjointBoxLayout& a_level,
                                   const BoxLayout& a_dest,
                                   const ProblemDomain& a_domain,
                                   const IntVect& a_destGhost,
                                   bool a_exchange)
{
  define(a_level, a_dest, a_domain, a_destGhost, a_exchange);
}


// ---------------------------------------------------------
MultiBlockCopier::MultiBlockCopier(const DisjointBoxLayout& a_level,
                                   const BoxLayout& a_dest,
                                   const IntVect& a_ghost,
                                   const Box& a_levelBlockDomain,
                                   const Box& a_destBlockDomain,
                                   bool  a_exchange)
{
  define(a_level, a_dest, a_ghost,
         a_levelBlockDomain, a_destBlockDomain, a_exchange);
}


// ---------------------------------------------------------
MultiBlockCopier::MultiBlockCopier(const DisjointBoxLayout& a_level,
                                   const BoxLayout& a_dest,
                                   const ProblemDomain& a_domain,
                                   const IntVect& a_ghost,
                                   const Box& a_levelBlockDomain,
                                   const Box& a_destBlockDomain,
                                   bool  a_exchange)
{
  define(a_level, a_dest, a_domain, a_ghost,
         a_levelBlockDomain, a_destBlockDomain, a_exchange);
}


// ---------------------------------------------------------
void MultiBlockCopier::clear()
{
  Copier::clear();
}


// ---------------------------------------------------------
void MultiBlockCopier::define(const DisjointBoxLayout& a_level,
                              const BoxLayout& a_dest,
                              bool a_exchange)
{
  // no block domains, so disallow this function
  MayDay::Error("Can't initialize MultiBlockCopier without block domains");
}


// ---------------------------------------------------------
void MultiBlockCopier::define(const DisjointBoxLayout& a_level,
                              const BoxLayout& a_dest,
                              const ProblemDomain& a_domain,
                              bool a_exchange)
{
  // call define, even though it won't do anything
  define(a_level, a_dest, a_domain, IntVect::Zero, a_exchange);
}


// ---------------------------------------------------------
void MultiBlockCopier::define(const DisjointBoxLayout& a_level,
                              const BoxLayout& a_dest,
                              const IntVect& a_ghost,
                              bool a_exchange)
{
  // no block domains, so disallow this function
  MayDay::Error("Can't initialize MultiBlockCopier without block domains");
}


// ---------------------------------------------------------
void MultiBlockCopier::define(const BoxLayout& a_level,
                              const BoxLayout& a_dest,
                              const ProblemDomain& a_domain,
                              const IntVect& a_ghost,
                              bool  a_exchange)
{
  // no block domains, so disallow this function
  MayDay::Error("Can't initialize MultiBlockCopier without block domains");
}


// ---------------------------------------------------------
void MultiBlockCopier::define(const DisjointBoxLayout& a_level,
                              const BoxLayout& a_dest,
                              const IntVect& a_ghost,
                              const Box& a_levelBlockDomain,
                              const Box& a_destBlockDomain,
                              bool  a_exchange)
{
  const ProblemDomain domain = a_level.physDomain();
  define(a_level, a_dest, domain, a_ghost,
         a_levelBlockDomain, a_destBlockDomain, a_exchange);
}


// ---------------------------------------------------------
void MultiBlockCopier::define(const BoxLayout& a_level,
                              const BoxLayout& a_dest,
                              const ProblemDomain& a_domain,
                              const IntVect& a_ghost,
                              const Box& a_levelBlockDomain,
                              const Box& a_destBlockDomain,
                              bool  a_exchange)
{
  CH_TIME("MultiBlockCopier::define");
  CH_assert(a_level.isClosed());
  CH_assert(a_dest.isClosed());
  //  CH_assert(a_level.checkPeriodic(a_domain));

  clear();
  m_levelBlockDomain = a_levelBlockDomain;
  m_destBlockDomain = a_destBlockDomain;
  m_isDefined = true;
  buffersAllocated = false;
  //bool self = a_dest == a_level;
  const BoxLayout&         level= a_level;
  const BoxLayout&         dest = a_dest;

  // set up vector of dataIndexes to keep track of which
  // "to" boxes are not completely contained within the primary
  // domain.  these boxes are then candidates for filling by
  // periodic images of the "from" data.
  Vector<DataIndex> periodicallyFilledToVect;

  // in order to cull which "from" data may be needed to
  // fill the "to" data, keep track of the radius around the
  // primary domain in which all these cells lie.
  // do this by incrementally growing the domain box and
  // keeping track of what this radius is.
  // just to make things simpler, start off with a radius of one
  Box grownDomainCheckBox = a_domain.domainBox();
  grownDomainCheckBox.grow(1);
  int periodicCheckRadius = 1;

  // since valid regions of the "from" DBL may also be outside
  // the primary domain, need to keep track of whether any of these
  // need to be checked separately.
  Vector<DataIndex> periodicFromVect;
  // use same domain trick here as well
  Box grownFromDomainCheckBox = a_domain.domainBox();
  int periodicFromCheckRadius = 1;

  Box domainBox(a_domain.domainBox());
  bool isPeriodic = false;
  if (!domainBox.isEmpty())
    isPeriodic = a_domain.isPeriodic();

  // (dfm -- 9/13/05) as currently written, the Copier won't correctly
  // handle periodic cases where the number of ghost cells is greater
  // than the width of the domain.  We _should_ do multiple wraparounds,
  // but we don't. So, put in this assertion. We can revisit this if it
  // becomes an issue
  if (isPeriodic)
    {
      for (int dir=0; dir<SpaceDim; dir++)
        {
          if (a_domain.isPeriodic(dir))
            {
              CH_assert (a_ghost[dir] <= domainBox.size(dir));
            }
        }
    }

  unsigned int myprocID = procID();

  // The following 4 for loops are the result of a performance optimization.
  // When increasing the size of the problem, we found that the code was
  // looping over every destination box for every source box which was N1*N2
  // loop iterations (essentially an N-squared approach).
  // The following code attempts to simply reduce N1 and N2 by first separating
  // the boxes (or LayoutIndexes to boxes) that reside on the current processor.
  // Then the loop to determine which boxes of the first list intersect with
  // which boxes of the second list can be done in N1' * N2' iterations,
  // where N1' is the reduced N1 and N2' is the reduced N2.
  // We have to break up the assigning of MotionItems into two separate
  // loops and be careful about the local copies.  These 4 loops are
  // significantly faster than the original for loop -- _especially_
  // for large problems.  (ndk)

#ifdef CH_MPI  // don't need to do this in serial
  // make a vector of boxes (or LayoutIndexes to boxes) from destination layout
  // that are known to reside on this processor.
  vector<DataIndex> vectorDestDI;
  vector<DataIndex> vectorDestOnProcDI;
  for (LayoutIterator to(a_dest.layoutIterator()); to.ok(); ++to)
    {
      // special condition for MultiBlockCopier
      if (m_destBlockDomain.contains(a_dest[to]))
        {
          vectorDestDI.push_back(DataIndex(to()));
          if (myprocID == dest.procID(to()))
            {
              vectorDestOnProcDI.push_back(DataIndex(to()));
            }
        }
    }

  // make a vector of boxes (or LayoutIndexes to boxes) from "level"/src layout
  // that are known to reside on this processor.
  vector<DataIndex> vectorLevelDI;
  vector<DataIndex> vectorLevelOnProcDI;
  for (LayoutIterator from(a_level.layoutIterator()); from.ok(); ++from)
    {
      // special condition for MultiBlockCopier
      if (m_levelBlockDomain.contains(a_level[from]))
        {
          vectorLevelDI.push_back(DataIndex(from()));
          if (myprocID == level.procID(from()))
            {
              vectorLevelOnProcDI.push_back(DataIndex(from()));
            }
        }
    }
#else
  // in serial, it's not very interesting as it's all of them.
  vector<DataIndex> vectorDestOnProcDI;
  vector<DataIndex> vectorLevelDI;
  for (LayoutIterator to(a_dest.layoutIterator()); to.ok(); ++to)
    {
      // special condition for MultiBlockCopier
      if (m_destBlockDomain.contains(a_dest[to]))
        {
          vectorDestOnProcDI.push_back(DataIndex(to()));
        }
    }
  for (LayoutIterator from(a_level.layoutIterator()); from.ok(); ++from)
    {
      // special condition for MultiBlockCopier
      if (m_levelBlockDomain.contains(a_level[from]))
        {
          vectorLevelDI.push_back(DataIndex(from()));
        }
    }
#endif

  // loop over all dest/to DI's on my processor
  for (vector<DataIndex>::iterator vdi=vectorDestOnProcDI.begin();
      vdi != vectorDestOnProcDI.end(); ++vdi)
  {

    // at this point, i know myprocID == toProcID
    const DataIndex todi(*vdi);
    Box ghost(dest[todi]);
    ghost.grow(a_ghost);

    // then for each level/from DI, see if they intersect
    for (vector<DataIndex>::iterator vli = vectorLevelDI.begin();
        vli != vectorLevelDI.end(); ++vli)
    {
      const DataIndex fromdi(*vli);
      const unsigned int fromProcID = level.procID(fromdi);
      const Box& fromBox = level[fromdi];
      if (fromBox.bigEnd(0) < ghost.smallEnd(0))
      {
        //can skip rest cuz we haven't gotten to something interesting
        continue;
      }

      if (ghost.intersectsNotEmpty(fromBox))
      {
        Box box(ghost); // ??
        box&=fromBox;   // ??
        MotionItem* item = new (s_motionItemPool.getPtr())
          MotionItem(fromdi, todi, box);
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
      if (fromBox.smallEnd(0) > ghost.bigEnd(0))
      {
        //can break out of loop, since we know that the smallEnd
        // of all the remaining boxes are lexigraphically beyond this ghosted box.
        break;
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
    Box ghost(dest[todi]);
    ghost.grow(a_ghost);

    const unsigned int toProcID = dest.procID(todi);

    // then for each level/from DI on this processor, see if they intersect
    for (vector<DataIndex>::iterator vli = vectorLevelOnProcDI.begin();
        vli != vectorLevelOnProcDI.end(); ++vli)
    {

      // at this point, i know myprocID == fromProcID

      const DataIndex fromdi(*vli);
      const Box& fromBox = level[fromdi];

      if (fromBox.bigEnd(0) < ghost.smallEnd(0))
      {
        //can skip rest cuz we haven't gotten to something interesting
        continue;
      }

      if (ghost.intersectsNotEmpty(fromBox))
      {
        Box box(ghost); // ??
        box&=fromBox;   // ??

        if (toProcID == myprocID)
        {
          // local move
          // don't push back here!  or you will get two.
          //     we already did it above...
          //m_localMotionPlan.push_back(item);
        }
        else
        {
          MotionItem* item = new (s_motionItemPool.getPtr())
            MotionItem(fromdi, todi, box);
          if (item == NULL)
          {
            MayDay::Error("Out of Memory in copier::define");
          }

          item->procID = toProcID;
          m_fromMotionPlan.push_back(item);
        }
      }
      if (fromBox.smallEnd(0) > ghost.bigEnd(0))
      {
        //can break out of loop, since we know that the smallEnd
        // of all the remaining boxes are lexigraphically beyond this ghosted box.
        break;
      }
    }
  }
#endif

  // put periodic intersection checking in here for "to" boxes
  if (isPeriodic)
  {
    for (LayoutIterator to(a_dest.layoutIterator()); to.ok(); ++to)
    {

      Box ghost(dest[to()]);
      ghost.grow(a_ghost);
      //unsigned int toProcID = dest.procID(to());  // unused variable

      // only do this if ghost box hangs over domain edge
      if (!domainBox.contains(ghost))
      {
        // add the dataIndex for this box to the list
        // of boxes which we need to come back to
        periodicallyFilledToVect.push_back(DataIndex(to()));
        // now check to see if we need to grow the
        // periodic check radius
        if (!grownDomainCheckBox.contains(ghost))
        {
          // grow the domainCheckBox until it contains ghost
          while (!grownDomainCheckBox.contains(ghost))
          {
            grownDomainCheckBox.grow(1);
            periodicCheckRadius++;
          }
        } // end if we need to grow radius around domain

      } //end if ghost box is not contained in domain
    } // end if periodic
  }

  // Here ends the so-called N-squared optimizations.  the rest is unchanged. (ndk)

  // now do periodic checking, if necessary
  if (isPeriodic)
    {

      // the only "from" boxes we will need to check
      // will be those within periodicCheckRadius of the
      // domain boundary. so, create a box to screen out
      // those which we will need to check.
      Box shrunkDomainBox = a_domain.domainBox();
      shrunkDomainBox.grow(-periodicCheckRadius);

      ShiftIterator shiftIt = a_domain.shiftIterator();
      IntVect shiftMult(domainBox.size());

      // now loop over "from" boxes
      for (LayoutIterator from(a_level.layoutIterator()); from.ok(); ++from)
        {
          // first check to see whether we need to look at this box
          const Box& fromBox = level[from()];

          // special condition for MultiBlockCopier
          if (m_levelBlockDomain.contains(fromBox))
            {
          if (!shrunkDomainBox.contains(fromBox))
            {
              unsigned int fromProcID = level.procID(from());

              // check to see if fromBox is contained in domain,
              // if not, add it to the list of fromBoxes we need to
              // go back and check separately to see if it will
              // fill one of the "to" boxes
              if (!domainBox.contains(fromBox))
                {
                  periodicFromVect.push_back(DataIndex(from()));

                  if (!grownFromDomainCheckBox.contains(fromBox))
                    {
                      while (!grownFromDomainCheckBox.contains(fromBox))
                        {
                          grownFromDomainCheckBox.grow(1);
                          periodicFromCheckRadius++;
                        }
                    } // end if we need to grow domain check box
                } // end if fromBox is outside domain

              // now loop over those "to" boxes which were not contained
              // in the domain
              for (int toRef=0; toRef<periodicallyFilledToVect.size(); toRef++)
                {
                  DataIndex toIndex = periodicallyFilledToVect[toRef];
                  unsigned int toProcID = dest.procID(toIndex);

                  // don't worry about anything that doesn't involve this proc
                  if (toProcID != myprocID && fromProcID != myprocID)
                    {
                      // do nothing
                    }
                  else
                    {
                      Box ghost(dest[toIndex]);
                      ghost.grow(a_ghost);
                      // now need to loop over shift vectors and look at images
                      for (shiftIt.begin(); shiftIt.ok(); ++shiftIt)
                        {
                          IntVect shiftVect(shiftIt()*shiftMult);
                          ghost.shift(shiftVect);
                          if (ghost.intersectsNotEmpty(fromBox)) // rarely happens
                            {
                              Box intersectBox(ghost);
                              intersectBox &= fromBox;
                              Box toBox(intersectBox);
                              toBox.shift(-shiftVect);
                              MotionItem* item = new (s_motionItemPool.getPtr())
                                MotionItem(DataIndex(from()), DataIndex(toIndex),
                                           intersectBox, toBox);
                              if (item == NULL)
                                {
                                  MayDay::Error("Out of Memory in copier::define");
                                }
                              if (toProcID == fromProcID) // local move
                                m_localMotionPlan.push_back(item);
                              else if (fromProcID == myprocID)
                                {
                                  item->procID = toProcID;
                                  m_fromMotionPlan.push_back(item);
                                }
                              else
                                {
                                  item->procID = fromProcID;
                                  m_toMotionPlan.push_back(item);
                                }

                            } // end if shifted box intersects

                          ghost.shift(-shiftVect);
                        } // end loop over shift vectors
                    } // end if either from box or to box are on this proc
                } // end loop over destination boxes
            } // end if source box is close to domain boundary
            } // end if source block domain contains source box
        } // end loop over destination boxes

      // now go back through the "from" boxes which were outside
      // the domain and see if they intersect any toBoxes
      if (periodicFromVect.size() != 0)
        {
          // the only "to" boxes we will need to check
          // will be those within periodicCheckRadius of the
          // domain boundary. so, create a box to screen out
          // those which we will need to check.
          shrunkDomainBox = a_domain.domainBox();
          shrunkDomainBox.grow(-periodicFromCheckRadius);

          // now loop over the "to" boxes
          for (LayoutIterator to(a_dest.layoutIterator()); to.ok(); ++to)
            {
              // first check to see whether we need to look at this box
              Box ghost(dest[to()]);
              ghost.grow(a_ghost);

              // special condition for MultiBlockCopier
              if (m_destBlockDomain.contains(a_dest[to]))
                {
              if (!shrunkDomainBox.contains(ghost))
                {
                  unsigned int toProcID = a_dest.procID(to());

                  // now loop over those "from" boxes which are not
                  // contained by the domain
                  for (int fromRef = 0; fromRef<periodicFromVect.size(); fromRef++)
                    {
                      DataIndex fromIndex = periodicFromVect[fromRef];
                      const Box& fromBox = level[fromIndex];
                      // special condition for MultiBlockCopier
                      if (m_levelBlockDomain.contains(fromBox))
                        {
                      unsigned int fromProcID = level.procID(fromIndex);

                      // don't worry about anything which doesn't involve
                      // this proc
                      if (toProcID != myprocID && fromProcID != myprocID)
                        {
                          // do nothing
                        }
                      else
                        {
                          // now need to loop over shift vectors and look at images
                          for (shiftIt.begin(); shiftIt.ok(); ++shiftIt)
                            {
                              IntVect shiftVect(shiftIt()*shiftMult);
                              ghost.shift(shiftVect);
                              if (ghost.intersectsNotEmpty(fromBox))
                                {
                                  Box intersectBox(ghost);
                                  intersectBox &= fromBox;
                                  Box toBox(intersectBox);
                                  toBox.shift(-shiftVect);
                                  MotionItem* item = new (s_motionItemPool.getPtr())
                                    MotionItem(DataIndex(fromIndex), DataIndex(to()),
                                               intersectBox, toBox);
                                  if (item == NULL)
                                    {
                                      MayDay::Error("Out of Memory in copier::define");
                                    }
                                  if (toProcID == fromProcID) // local move
                                    m_localMotionPlan.push_back(item);
                                  else if (fromProcID == myprocID)
                                    {
                                      item->procID = toProcID;
                                      m_fromMotionPlan.push_back(item);
                                    }
                                  else
                                    {
                                      item->procID = fromProcID;
                                      m_toMotionPlan.push_back(item);
                                    }

                                } // end if shifted box intersects

                              ghost.shift(-shiftVect);
                            } // end loop over shift vectors
                        } // end if either from box or to box are on this proc
                        } // end if source block contains source box
                    } // end loop over "from" boxes
                } // end if destination box is close to domain boundary
            } // end if dest block domain contains dest box
            } // end loop over destination boxes
        } // end if any of the "From" boxes were outside the domain

    } // end if we need to do anything for periodicity
  sort();
}


// ---------------------------------------------------------
MultiBlockCopier::~MultiBlockCopier()
{
  clear();
}

#include "NamespaceFooter.H"
