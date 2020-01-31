#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cstring>

#include "LevelData.H"
#include "FluxBox.H"
#include "NamespaceHeader.H"

template < > void LevelData<FluxBox>::degenerateLocalOnly( LevelData<FluxBox>& a_to, const SliceSpec& a_sliceSpec ) const

{
  DisjointBoxLayout toDBL;
  m_disjointBoxLayout.degenerate( toDBL, a_sliceSpec, true );
  IntVect toGhost;
  for ( int i=0;i<CH_SPACEDIM;++i )
    {
      if ( i != a_sliceSpec.direction )
        {
          toGhost[i] = m_ghost[i];
        } else
        {
          toGhost[i] = 0;
        }   
    }
  a_to.define( toDBL, this->nComp(), toGhost );
  
  // manage copyTo ourselves to maintain locality
  DataIterator fromDit = this->dataIterator();
  DataIterator toDit = a_to.dataIterator();
  for (toDit.begin(); toDit.ok(); ++toDit)
    {
      fromDit.begin();
      bool done = false;
      while (!done)
        {
          if (m_disjointBoxLayout[fromDit].contains(toDBL[toDit]))
            {
              // boxes intersect, do copy
              // loop over directions
              FluxBox& toFlux = a_to[toDit];
              const FluxBox& fromFlux = this->operator[](fromDit);
              for (int dir=0; dir<SpaceDim; dir++)
                {
                  FArrayBox& toFAB = toFlux[dir];
                  const FArrayBox& fromFAB = fromFlux[dir];
                  // note that we're including ghost cells here
                  Box intersectBox = toFAB.box();
                  intersectBox &= fromFAB.box();
                  // now do copy
                  toFAB.copy(fromFAB, intersectBox);
                } // end loop over directions
              done = true;
            } // end if we found an intersection
          else
            {
              ++fromDit;
              // probably want to check for lexigraphical sorting done-ness
              
              // are we done looping through fromBoxes?
              if (!fromDit.ok()) done = true;
            } // end if we didn't find a match
        } // end while loop over fromBoxes
      
    } // end loop over toBoxes
  
}


#include "NamespaceFooter.H"
