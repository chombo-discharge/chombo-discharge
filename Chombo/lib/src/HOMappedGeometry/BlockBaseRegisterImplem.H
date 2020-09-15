#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _BLOCKBASEREGISTERIMPLEM_H_
#define _BLOCKBASEREGISTERIMPLEM_H_

#include "BlockBaseRegister.H"
#include "BoxIterator.H"
#include "CornerCopier.H"
#include "CH_Timer.H"

#include "NamespaceHeader.H"

inline int getBlockBdryIndex( const int& a_dir, const Side::LoHiSide& a_side );

template <class T>
BlockBaseRegister<T>::BlockBaseRegister( const MultiBlockCoordSys* a_mblock,
                                      /// CELL-centered grids
                                      const DisjointBoxLayout& a_grids,
                                      /// ghost width for state variables
                                      const int& a_ghost )
  : m_verbosity(0),
    m_hasData(false),
    m_hasDefinedDataLayouts(false)
{
   BlockBaseRegister<T>();
   define( a_mblock, a_grids, a_ghost );
}


template <class T>
BlockBaseRegister<T>::~BlockBaseRegister()
{
   //do nothing for now
}


template <class T>
void BlockBaseRegister<T>::define( const MultiBlockCoordSys* a_mblock,
                                /// CELL-centered grids
                                const DisjointBoxLayout& a_grids,
                                /// ghost width for state variables
                                const int& a_ghost )
{
   CH_TIME("BlockBaseRegister<T>::define");

   m_coordsys = (MultiBlockCoordSys*) a_mblock;
   m_grids = a_grids;
   m_ghost = a_ghost;
   for (int dir(0); dir<CH_SPACEDIM; dir++)
   {
      m_ghostVect.setVal( dir, m_ghost );
   }

   const ProblemDomain& inDomain( m_grids.physDomain() );
   Box newBoundingBox( inDomain.domainBox() );
   newBoundingBox.grow( m_ghostVect + IntVect::Unit );
   m_domain.define( newBoundingBox );
   
   m_blockBdryIndexMap.resize( m_coordsys->numBlocks() );

   buildRegisterLayout();
   buildSrcRegisterLayout();
   buildExchangeRegisterLayout();

   m_isDefined = true;
}


template <class T>
void BlockBaseRegister<T>::defineDataLayouts( const int a_ncomp )
{
   for (int j(0); j<2*CH_SPACEDIM; j++)
   {
      IntVect noDirGhost( m_ghost * IntVect::Unit );
      int dir(j%CH_SPACEDIM);  
      noDirGhost[dir] = 0;

      m_register[j].define( m_registerLayout[j], a_ncomp, noDirGhost );
//      m_excRegister[j].define( m_excRegisterLayout[j], a_ncomp );
      m_excRegister[j].define( m_excRegisterLayout[j], a_ncomp, m_ghostVect );
//      m_srcRegister[j].define( m_srcRegisterLayout[j], a_ncomp );
      m_srcRegister[j].define( m_srcRegisterLayout[j], a_ncomp, noDirGhost );
   }
   m_hasDefinedDataLayouts = true;
}


template <class T>
void BlockBaseRegister<T>::clearDataLayouts()
{
   for (int j(0); j<2*CH_SPACEDIM; j++)
   {
      m_register[j].clear();
      m_excRegister[j].clear();
      m_srcRegister[j].clear();
   }
   m_hasDefinedDataLayouts = false;
}


template <class T>
void BlockBaseRegister<T>::store( const T& a_data,
                               const DataIndex& a_dataIndex,
                               int a_dir,
                               Side::LoHiSide a_side )
{
   CH_TIME("BlockBaseRegister<T>::store");

   if ( isDefined() )
   {

      if (!hasStorage())
      {
         defineDataLayouts( a_data.nComp() );
      }

      const int index( getBlockBdryIndex( a_dir, a_side ) );

      LevelData<T >& dataRegister( m_register[index] );

      LevelData<T >& dataExcRegister( m_excRegister[index] );
      if (hasInterface( a_dataIndex, a_dir, a_side ))
      {
         
         const DataIndex& regIndex( m_grdToDstMap[index][a_dataIndex] );
         T& regData( dataRegister[regIndex] );
         
         // We cast away the const-ness in a_flux for the scope of this
         // function so that we can shift without copying.
         T& inData( const_cast<T&>(a_data) ); // Muhahaha.
         inData.shiftHalf( a_dir, -sign(a_side) );
         
         Box srcBox( adjCellBox( m_grids[a_dataIndex], a_dir, a_side, -1 ) );
         if (m_ghost>0)
         {
            for (int gdir(0); gdir<CH_SPACEDIM; ++gdir)
            {
               if (gdir!=a_dir)
               {
                  srcBox.grow( gdir, m_ghost );
               }
               else
               {
                  srcBox.growDir( gdir, flip(a_side), m_ghost );
               }
            }
         }
         
         const Box window( srcBox & inData.box() );
         regData.copy( inData, window & regData.box() );
         
         // Done being evil...for now ;-)
         inData.shiftHalf( a_dir, sign(a_side) );
         
         const DataIndex& excIndex( m_dstToExcMap[index][regIndex] );
         T& excData( dataExcRegister[excIndex] );
         const Box excWindow( regData.box() & excData.box() );
         excData.copy( regData, excWindow );
         
         m_hasData = true;
      }
   }
   else
   {
      MayDay::Error(" BlockBaseRegister<T>::store:  BlockBaseRegister<T> is not defined!");
   }
}


template <class T>
void BlockBaseRegister<T>::increment( const T& a_data,
                                   const DataIndex& a_dataIndex,
                                   int a_dir,
                                   Side::LoHiSide a_side )
{
   CH_TIME("BlockBaseRegister<T>::increment");

   if ( isDefined() )
   {

      if (!hasStorage())
      {
         defineDataLayouts( a_data.nComp() );
      }

      const int index( getBlockBdryIndex( a_dir, a_side ) );

      BoxLayoutData<T >& dataRegister( m_register[index] );
      LevelData<T >& dataExcRegister( m_excRegister[index] );
      if (hasInterface( a_dataIndex, a_dir, a_side ))
      {
         
         const DataIndex& regIndex( m_grdToDstMap[index][a_dataIndex] );
         T& regData( dataRegister[regIndex] );
         
         // We cast away the const-ness in a_flux for the scope of this
         // function so that we can shift without copying.
         T& inData( const_cast<T&>(a_data) ); // Muhahaha.
         inData.shiftHalf( a_dir, -sign(a_side) );
         
         Box srcBox( adjCellBox( m_grids[a_dataIndex], a_dir, a_side, -1 ) );
         if (m_ghost>0)
         {
            for (int gdir(0); gdir<CH_SPACEDIM; ++gdir)
            {
               if (gdir!=a_dir)
               {
                  srcBox.grow( gdir, m_ghost );
               }
               else
               {
                  srcBox.growDir( gdir, flip(a_side), m_ghost );
               }
            }
         }
         
         const Box window( srcBox & inData.box() );
         const int n_comp( a_data.nComp() );
         regData.plus( inData, window & regData.box(), 0, 0, n_comp );
         
         // Done being evil...for now ;-)
         inData.shiftHalf( a_dir, sign(a_side) );
         
         const DataIndex& excIndex( m_dstToExcMap[index][regIndex] );
         T& excData( dataExcRegister[excIndex] );
         excData.copy( regData, regData.box() & excData.box() );
         
         m_hasData = true;
      }
   }
   else
   {
      MayDay::Error(" BlockBaseRegister<T>::increment:  BlockBaseRegister<T> is not defined!");
   }
}


template <class T>
void BlockBaseRegister<T>::exchange()
{
   CH_TIME("BlockBaseRegister<T>::exchange");
   if ( isDefined() )
   {
      if ( hasData() )
      {
         copyToSrcRegister();
         copyToDstRegister();
         if (m_verbosity>0) printRegisters();
      }
   }
   else
   {
      MayDay::Error(" BlockBaseRegister<T>::exchange:  BlockBaseRegister<T> is not defined!");
   }
}


template <class T>
inline
void zeroBoxLayoutData( BoxLayoutData<T >& a_layoutData )
{
   for (DataIterator dit(a_layoutData.dataIterator()); dit.ok(); ++dit)
   {
      a_layoutData[dit].setVal( 0.0 );
   }
}


template <class T>
void BlockBaseRegister<T>::zeroRegister( const int a_n_comp )
{
   CH_TIME("BlockBaseRegister<T>::zeroRegister");
   if ( isDefined() )
   {
      if (!hasStorage())
      {
         if (a_n_comp>0)
         {
            defineDataLayouts( a_n_comp );
         }
         else
         {
            MayDay::Error(" BlockRegister::zeroRegister:  Attempt to zero register of unknown component length!");
         }
      }
      else if ( a_n_comp>0 && a_n_comp != m_register[0].nComp() )
      {
         clearDataLayouts();
         defineDataLayouts( a_n_comp );
      }
      for (int index(0); index<2*CH_SPACEDIM; index++)
      {
         zeroBoxLayoutData( m_register[index] );
         zeroBoxLayoutData( m_excRegister[index] );
         zeroBoxLayoutData( m_srcRegister[index] );
      }
      m_hasData = true;
   }
   else
   {
      MayDay::Error(" BlockRegister::zeroRegister:  Attempt to zero undefined register!");
   }
}


template <class T>
bool BlockBaseRegister<T>::hasInterface( const DataIndex& a_dataIndex,
                                  int a_dir,
                                  Side::LoHiSide a_side ) const
{
   CH_TIME("BlockBaseRegister<T>::hasInterface");
   // JAFH: We assume the a_dataIndex is valid; should we check?

   // construct ghost box grown in a_dir direction on a_side
   const Box& box( m_grids[a_dataIndex] );
   const Box clipBox( adjCellBox( box, a_dir, a_side, 1 ) );

   // construct block box grown in a_dir direction on a_side
   const int block( getBlockID( a_dataIndex ) );
   const Box& blockBox( getBlockBox( block ) );
   const Box blockBdryBox( adjCellBox( blockBox, a_dir, a_side, 1 ) );

   // Has interface if touches interfacial block boundary
   if (blockBdryBox.intersectsNotEmpty( clipBox ))
   {
      if (isInterface( block, a_dir, a_side ))
      {
         return true;
      }
   }
   return false;
}


template <class T>
void BlockBaseRegister<T>::fill( T& a_data,
                          const DataIndex& a_dataIndex,
                          int a_dir,
                          Side::LoHiSide a_side,
                          Side::LoHiSide a_sideData ) const
{
   CH_TIME("BlockBaseRegister<T>::fill");

   if ( isDefined() )
   {
      if ( hasData() )
      {
         if ( hasInterface( a_dataIndex, a_dir, a_side ) )
         {
            const int index( getBlockBdryIndex( a_dir, a_side ) );
            const DataIndex& regIndex( m_grdToDstMap[index][a_dataIndex] );
            const T& regData( m_register[index][regIndex] );

            const int srcSide( sign( a_side ) * sign( a_sideData ) );
            const Box& box( m_grids[a_dataIndex] );
            Box srcBox( adjCellBox( box, a_dir, a_side, srcSide ) );
            const int shift( sign( a_sideData ) );

            if (m_ghost>0)
            {
               for (int gdir(0); gdir<CH_SPACEDIM; ++gdir)
               {
                  if (gdir!=a_dir)
                  {
                     srcBox.grow( m_ghost );
                  }
                  else
                  {
                     Side::LoHiSide side( (srcSide>0) ? Side::Lo : Side::Hi );
                     srcBox.growDir( gdir, side, m_ghost );
                  }
               }
            }

            a_data.shiftHalf( a_dir, shift );
            const Box window( srcBox & a_data.box() & regData.box() );
            a_data.copy( regData, window );
            a_data.shiftHalf( a_dir, -shift );
         }
      }
   }
   else
   {
      MayDay::Error(" BlockBaseRegister<T>::fill:  BlockBaseRegister<T> is not defined!");
   }
}


//////////////////////////////////////////////////////////////////////////////
// Begin Private Methods
//////////////////////////////////////////////////////////////////////////////


template <class T>
inline
int BlockBaseRegister<T>::getBlockID( const DataIndex& a_dataIndex ) const
{
   const Box& box( m_grids[a_dataIndex] );
   return m_coordsys->whichBlock( box );
}


template <class T>
inline
const Box& BlockBaseRegister<T>::getBlockBox( const int& a_block ) const
{
   const Vector<Box>& blocks( m_coordsys->mappingBlocks() );
   return blocks[a_block];
}


template <class T>
inline
const Box& BlockBaseRegister<T>::getNeighborBlockBox( const int& a_block,
                                               const int& a_blockBdryIndex ) const
{
   const Tuple<BlockBoundary, 2*SpaceDim>&
      boundary( m_coordsys->boundaries()[a_block] );
   const int neighborBlock( boundary[a_blockBdryIndex].neighbor() );
   return getBlockBox( neighborBlock );
}


inline
int getBlockBdryIndex( const int& a_dir, const Side::LoHiSide& a_side )
{
   return ( a_side * CH_SPACEDIM + a_dir );
}


template <class T>
inline
IndicesTransformation BlockBaseRegister<T>::getTransformation(
   const int& a_block,
   const int& a_blockBdryIndex ) const
{
   const Tuple<BlockBoundary, 2*SpaceDim>&
      boundary( m_coordsys->boundaries()[a_block] );
   return boundary[a_blockBdryIndex].getTransformation();
}


template <class T>
inline
int BlockBaseRegister<T>::getSrcBlockBdryIndex(
   const int& a_block,
   const int& a_dir,
   const Side::LoHiSide& a_side )
{
   // assumes one neighbor per block boundary
   // assumes connection only at block boundaries
   const int index( getBlockBdryIndex( a_dir, a_side ) );
   BlockBdryIndexMap::iterator it( m_blockBdryIndexMap[a_block].find( index ) );
   if (it == m_blockBdryIndexMap[a_block].end())
   {

      const IndicesTransformation& itrans( getTransformation( a_block, index ) );
      Box mySrcBox( itrans.transformFwd( getBlockBox( a_block ) ) );
      const int neighborDir( itrans.getPermutation()[a_dir] );
      mySrcBox.growLo( neighborDir );
      const Box& neighborBox( getNeighborBlockBox( a_block, index ) );
      const Box intersectBox( mySrcBox & neighborBox );

      Side::LoHiSide neighborSide( Side::Hi );
      if ( intersectBox.isEmpty() )
      {
         neighborSide = flip( neighborSide );
      }
      int src_index( getBlockBdryIndex( neighborDir, neighborSide ) );
      m_blockBdryIndexMap[a_block].insert(
         BlockBdryIndexMap::value_type( index, src_index ) );
   }
   return m_blockBdryIndexMap[a_block][index];
}


template <class T>
inline
bool BlockBaseRegister<T>::isInterface( const int& a_block,
                                 const int& a_dir,
                                 const Side::LoHiSide& a_side ) const

{
   const Tuple<BlockBoundary, 2*SpaceDim>&
      boundary( m_coordsys->boundaries()[a_block] );
   const int index( getBlockBdryIndex( a_dir, a_side ) );
   return boundary[index].isInterface();
}


inline
void buildLocalIndexList( List<int>& a_list,
                          const Vector<Box>& a_boxes,
                          const Vector<int>& a_procs )
{
   int myRank(0);
#ifdef CH_MPI
   MPI_Comm_rank( Chombo_MPI::comm, &myRank );
#endif
   for (int i(0); i<a_boxes.size(); i++)
   {
      if (a_procs[i]==myRank)
      {
         a_list.append(i);
      }
   }
}


template <class T>
inline
void BlockBaseRegister<T>::createDstBoxLayout( DisjointBoxLayout& a_layout,
                                            LayoutData<DataIndex>& a_map,
                                            const Vector<Box>& a_boxes,
                                            const Vector<int>& a_procs,
                                            const Vector<DataIndex>& a_indices )
{
   a_layout.define( a_boxes, a_procs, m_domain );
   a_layout.close();

   List<int> vi;
   buildLocalIndexList( vi, a_boxes, a_procs );

   for (DataIterator dit(a_layout.dataIterator()); dit.ok(); ++dit)
   {
      const Box& box( a_layout[dit] );
      for (ListIterator<int> li( vi.first() ); li.ok(); ++li)
      {
         int i( li() );
         if ( box.eq( a_boxes[i] ) )
         {
            DataIndex foo( dit() );
            a_map[a_indices[i]] = foo;
//             a_map[a_indices[i]] = dit();
//          a_map[a_indices[i]] = dit();
            vi.remove( li() );
            break;
         }
      }
   }
   if (vi.isNotEmpty())
   {
      MayDay::Error(" BlockBaseRegister<T>::createDstBoxLayout:  One or more orphaned boxes!");
   }
}


template <class T>
inline
void BlockBaseRegister<T>::createSrcBoxLayout(
   DisjointBoxLayout& a_srcLayout,
   LayoutData<DataIndex> a_srcMap[2*CH_SPACEDIM],
   const Vector<Box>& a_srcBoxes,
   const Vector<int>& a_srcProcs,
   const Vector<DataIndex>& a_srcDataIndices,
   const Vector<int>& a_dstBlockBdryIndices,
   const LayoutData<DataIndex> a_dstMap[2*CH_SPACEDIM] )
{
   a_srcLayout.define( a_srcBoxes, a_srcProcs, m_domain );
   a_srcLayout.close();

   List<int> vi;
   buildLocalIndexList( vi, a_srcBoxes, a_srcProcs );

   for (DataIterator dit(a_srcLayout.dataIterator()); dit.ok(); ++dit)
   {
      const Box& box( a_srcLayout[dit] );
      for (ListIterator<int> li( vi.first() ); li.ok(); ++li)
      {
         int i( li() );
         if ( box.eq( a_srcBoxes[i] ) )
         {
            const int dstBlockBdryIndex( a_dstBlockBdryIndices[i] );
            const DataIndex& gridDataIndex( a_srcDataIndices[i] );
            const DataIndex& dstDataIndex( a_dstMap[dstBlockBdryIndex][gridDataIndex] );
            a_srcMap[dstBlockBdryIndex][dstDataIndex] = dit();
            vi.remove( li() );
            break;
         }
      }
   }
   if (vi.isNotEmpty())
   {
      MayDay::Error(" BlockBaseRegister<T>::createSrcBoxLayout:  One or more orphaned boxes!");
   }
}


template <class T>
inline
void BlockBaseRegister<T>::createExchangeBoxLayout(
   DisjointBoxLayout& a_exchangeLayout,
   LayoutData<DataIndex>& a_dstToExcMap,
   const Vector<Box>& a_boxes,
   const Vector<int>& a_procs,
   const Vector<DataIndex>& a_dataIndices,
   const LayoutData<DataIndex>& a_grdToDstMap )
{
   a_exchangeLayout.define( a_boxes, a_procs, m_domain );
   a_exchangeLayout.close();

   List<int> vi;
   buildLocalIndexList( vi, a_boxes, a_procs );

   for (DataIterator dit(a_exchangeLayout.dataIterator()); dit.ok(); ++dit)
   {
      const Box& box( a_exchangeLayout[dit] );
      for (ListIterator<int> li( vi.first() ); li.ok(); ++li)
      {
         int i( li() );
         if ( box.eq( a_boxes[i] ) )
         {
            const DataIndex& gridDataIndex( a_dataIndices[i] );
            const DataIndex& dstDataIndex( a_grdToDstMap[gridDataIndex] );
            a_dstToExcMap[dstDataIndex] = dit();
            vi.remove( li() );
            break;
         }
      }
   }
   if (vi.isNotEmpty())
   {
      MayDay::Error(" BlockBaseRegister<T>::createSrcBoxLayout:  One or more orphaned boxes!");
   }
}


template <class T>
void BlockBaseRegister<T>::buildRegisterLayout()
{
   Vector<Box> boxes[2*CH_SPACEDIM];
   Vector<int> boxProcMap[2*CH_SPACEDIM];
   Vector<DataIndex> dataIndex[2*CH_SPACEDIM];
   buildBoxVectors( boxes, boxProcMap, dataIndex, DST );

   // For each block boundary, create a BoxLayout for the destination registers
   // and define a LayoutData that will be used to associate each destination
   // register DataIndex with a source register DataIndex
   for (int dir(0); dir<CH_SPACEDIM; ++dir)
   {
      for (SideIterator side; side.ok(); ++side)
      {
         const int index( getBlockBdryIndex( dir, side() ) );
         m_grdToDstMap[index].define( m_grids );
         createDstBoxLayout( m_registerLayout[index],
                             m_grdToDstMap[index],
                             boxes[index],
                             boxProcMap[index],
                             dataIndex[index] );
         m_dstToSrcMap[index].define( m_registerLayout[index] );
         m_dstToExcMap[index].define( m_registerLayout[index] );
      }
   }
}


template <class T>
void BlockBaseRegister<T>::buildSrcRegisterLayout()
{
   Vector<Box> boxes[2*CH_SPACEDIM];
   Vector<int> boxProcMap[2*CH_SPACEDIM];
   Vector<DataIndex> dataIndex[2*CH_SPACEDIM];
   buildBoxVectors( boxes, boxProcMap, dataIndex, SRC );

   Vector<int> blockBdryIndexMap[2*CH_SPACEDIM];
   buildInverseBlockBdryIndexMap( blockBdryIndexMap );

   // For each block boundary, create a BoxLayout for the source registers.
   // At the same time, we associate the source DataIndexes with destination
   // DataIndexes.
   for (int dir(0); dir<CH_SPACEDIM; ++dir)
   {
      for (SideIterator side; side.ok(); ++side)
      {
         const int index( getBlockBdryIndex( dir, side() ) );
         createSrcBoxLayout( m_srcRegisterLayout[index],
                             m_dstToSrcMap,
                             boxes[index],
                             boxProcMap[index],
                             dataIndex[index],
                             blockBdryIndexMap[index],
                             m_grdToDstMap );
      }
   }
}


template <class T>
void BlockBaseRegister<T>::buildExchangeRegisterLayout()
{
   Vector<Box> boxes[2*CH_SPACEDIM];
   Vector<int> boxProcMap[2*CH_SPACEDIM];
   Vector<DataIndex> dataIndex[2*CH_SPACEDIM];
   buildBoxVectors( boxes, boxProcMap, dataIndex, EXC );

   // For each block boundary, create a BoxLayout for the source registers.
   // At the same time, we associate the source DataIndexes with destination
   // DataIndexes.
   for (int dir(0); dir<CH_SPACEDIM; ++dir)
   {
      for (SideIterator side; side.ok(); ++side)
      {
         const int index( getBlockBdryIndex( dir, side() ) );
         createExchangeBoxLayout( m_excRegisterLayout[index],
                                  m_dstToExcMap[index],
                                  boxes[index],
                                  boxProcMap[index],
                                  dataIndex[index],
                                  m_grdToDstMap[index] );
      }
   }
}


template <class T>
void BlockBaseRegister<T>::buildBoxVectors( Vector<Box> a_boxes[2*CH_SPACEDIM],
                                         Vector<int> a_boxProcMap[2*CH_SPACEDIM],
                                         Vector<DataIndex> a_dataIndex[2*CH_SPACEDIM],
                                         const UseType a_type )
{
   for (LayoutIterator lit(m_grids.layoutIterator()); lit.ok(); ++lit)
   {
      const Box& box( m_grids[lit] );
      const int block( m_coordsys->whichBlock( box ) );
      for (int dir(0); dir<CH_SPACEDIM; ++dir)
      {
         for (SideIterator si; si.ok(); ++si)
         {
            Side::LoHiSide side( si() );
            if (hasInterface( DataIndex( lit() ), dir, side ))
            {
               int index( getBlockBdryIndex( dir, side ) );

               Box registerBox;

               if (a_type==DST)
               {
                  registerBox = adjCellBox( box, dir, side, -1 );
                  registerBox.growDir( dir, side, 1 );
//                  registerBox.grow( m_ghost );
               }

               else if (a_type==SRC)
               {
                  registerBox = adjCellBox( box, dir, side, 1 );
                  if (m_ghost>0)
                  {
#if 0
                     for (int gdir(0); gdir<CH_SPACEDIM; gdir++)
                     {
                        if (gdir!=dir)
                        {
//                           registerBox.grow( gdir, m_ghost );
                        }
                        else
                        {
//                           registerBox.growDir( gdir, side, m_ghost );
                        }
                     }
#endif
                  }
                  const IndicesTransformation& itrans( getTransformation( block, index ) );
                  registerBox = itrans.transformFwd( registerBox );
                  index = getSrcBlockBdryIndex( block, dir, side );
               }

               else if (a_type==EXC)
               {
                  registerBox = adjCellBox( box, dir, side, -(1+m_ghost) );
                  if (m_ghost>0)
                  {
                     for (int gdir(0); gdir<CH_SPACEDIM; ++gdir)
                     {
                        if (gdir!=dir)
                        {
//                           registerBox.grow( gdir, m_ghost );
#if 0
                           if (hasInterface( DataIndex( lit() ), gdir, Side::Lo ))
                           {
                              registerBox.growLo( gdir, m_ghost );
                           }
                           if (hasInterface( DataIndex( lit() ), gdir, Side::Hi ))
                           {
                              registerBox.growHi( gdir, m_ghost );
                           }
#endif
                        }
                     }
                  }
               }

               else
               {
                  MayDay::Error(" BlockBaseRegister<T>::buildBoxVectors:  Unknown use type!");
               }

               a_boxes[index].push_back( registerBox );
               a_boxProcMap[index].push_back( m_grids.procID( lit() ) );
               a_dataIndex[index].push_back( DataIndex( lit() ) );
            }
         }
      }
   }
}


template <class T>
void BlockBaseRegister<T>::buildInverseBlockBdryIndexMap(
   Vector<int> a_map[2*CH_SPACEDIM] )
{
   for (LayoutIterator lit(m_grids.layoutIterator()); lit.ok(); ++lit)
   {
      const Box& box( m_grids[lit] );
      const int block( m_coordsys->whichBlock( box ) );
      for (int dir(0); dir<CH_SPACEDIM; ++dir)
      {
         for (SideIterator si; si.ok(); ++si)
         {
            Side::LoHiSide side( si() );
            if (hasInterface( DataIndex( lit() ), dir, side ))
            {
               const int index( getBlockBdryIndex( dir, side ) );
               const int transformedIndex( getSrcBlockBdryIndex( block, dir, side ) );
               a_map[transformedIndex].push_back( index );
            }
         }
      }
   }
}


template <class T>
inline
void BlockBaseRegister<T>::copyToSrcRegister()
{
   if (m_hasData)
   {
      for (int dir(0); dir<CH_SPACEDIM; dir++)
      {
         for (SideIterator side; side.ok(); ++side)
         {
            const int index( getBlockBdryIndex( dir, side() ) );
            // petermc, 23 May 2013: removed constness in order to do exchange
            LevelData<T >& srcLevelData( m_excRegister[index] );
//            BoxLayoutData<T >& dstBoxLayoutData( m_srcRegister[index] );
            LevelData<T >& dstLevelData( m_srcRegister[index] );
            if (m_ghost==0)
            {
               srcLevelData.copyTo( dstLevelData );
            }
            else
            {
               Copier copier;
               copier.ghostDefine( srcLevelData.disjointBoxLayout(),
                                   dstLevelData.disjointBoxLayout(),
                                   m_domain,
                                   srcLevelData.ghostVect(),
                                   dstLevelData.ghostVect() );
               srcLevelData.exchange(); // added by petermc, 23 May 2013
               srcLevelData.copyTo( dstLevelData, copier );
            }
         }
      }
   }
}


inline
IntVect permute( const IntVect& a_vec, const IntVect& a_permutation )
{
  IntVect tmp_vec;
  for (int dir(0); dir<SpaceDim; ++dir)
  {
     tmp_vec[dir] = a_vec[a_permutation[dir]];
  }
  return tmp_vec;
}

template <class T>
inline
void BlockBaseRegister<T>::copyToDstRegister()
{
   if (m_hasData)
   {
      for (int dir(0); dir<CH_SPACEDIM; dir++)
      {
         for (SideIterator side; side.ok(); ++side)
         {

            const int blockBdryIndex( getBlockBdryIndex( dir, side() ) );

            DataIterator dit( m_grids.dataIterator() );
            for (dit.begin(); dit.ok(); ++dit)
            {

               if (hasInterface( dit(), dir, side() ))
               {
                  const int block( getBlockID( dit() ) );

                  const int srcBlockBdryIndex( getSrcBlockBdryIndex( block, dir, side() ) );
                  const IndicesTransformation& itrans( getTransformation( block, blockBdryIndex ) );

                  const LayoutData<DataIndex>& grdToDst( m_grdToDstMap[blockBdryIndex] );
                  const DataIndex& regDataIndex( grdToDst[dit()] );
                  
                  const LayoutData<DataIndex>& dstToSrc( m_dstToSrcMap[blockBdryIndex] );
                  const DataIndex& srcDataIndex( dstToSrc[regDataIndex] );
                  
                  LevelData<T>& dstRegister( m_register[blockBdryIndex] );
                  const LevelData<T>& srcRegister( m_srcRegister[srcBlockBdryIndex] );
                  
                  T& dst( dstRegister[regDataIndex] );
                  const T& src( srcRegister[srcDataIndex] );
                  
                  const Box& srcWindow( src.box() );
                  const Box& regWindow( itrans.transformBack( srcWindow ) );
                  const IntVect perm( itrans.getPermutation() );
                  const IntVect sign( itrans.getSign() );
                  const IntVect trans( itrans.getTranslation() );

                  for (int n(0); n<dst.nComp(); n++)
                  {
                     BoxIterator bit(regWindow);
                     for (bit.begin(); bit.ok(); ++bit)
                     {
                        const IntVect dst_loc( bit() );
                        IntVect src_loc( permute( dst_loc, perm ) );
                        src_loc *= sign;
                        src_loc += trans;
                        dst(dst_loc,n) = src(src_loc,n);
                     }
                  }
               }
            }
            m_register[blockBdryIndex].exchange();
         }
      }
   }
}


#if 0
template <class T>
inline
void BlockBaseRegister<T>::print( const BoxLayoutData<T > a_register[2*CH_SPACEDIM] ) const
{
   for (int dir(0); dir<CH_SPACEDIM; dir++)
   {
      for (SideIterator side; side.ok(); ++side)
      {
         pout() << "========= dir: " << dir << "   ";
         pout() << "side: " << side() << "   ";
         const int index( getBlockBdryIndex( dir, side() ) );
         pout() << "block boundary index: " << index << std::endl;
         const BoxLayoutData<T >& thisRegister( a_register[index] );
       if (thisRegister.isDefined())
       {
            DataIterator dit( thisRegister.dataIterator() );
            for (dit.begin(); dit.ok(); ++dit)
            {
               const T& data( thisRegister[dit()] );
               const Box& box( data.box() );
               pout() << "--------------------------------------" << std::endl;
               pout() << "Box: " << box << endl;
               for (int ncomp(0); ncomp<data.nComp(); ++ncomp)
               {
                  pout() << "Component: " << ncomp << endl;
                  BoxIterator bit( box );
                  for (bit.begin(); bit.ok(); ++bit)
                  {
                     pout() << bit() << "\t" << data.get( bit(), ncomp ) << std::endl;;
                  }
               }
               pout() << "--------------------------------------" << std::endl;
            }
       }
      }
   }
}
#endif


template <class T>
inline
void BlockBaseRegister<T>::print( const LevelData<T > a_register[2*CH_SPACEDIM] ) const
{
   for (int dir(0); dir<CH_SPACEDIM; dir++)
   {
      for (SideIterator side; side.ok(); ++side)
      {
         pout() << "========= dir: " << dir << "   ";
         pout() << "side: " << side() << "   ";
         const int index( getBlockBdryIndex( dir, side() ) );
         pout() << "block boundary index: " << index << std::endl;
         const LevelData<T >& thisRegister( a_register[index] );
       if (thisRegister.isDefined())
       {
            DataIterator dit( thisRegister.dataIterator() );
            for (dit.begin(); dit.ok(); ++dit)
            {
               const T& data( thisRegister[dit()] );
               const Box& box( data.box() );
               pout() << "--------------------------------------" << std::endl;
               pout() << "Box: " << box << endl;
               for (int ncomp(0); ncomp<data.nComp(); ++ncomp)
               {
                  pout() << "Component: " << ncomp << endl;
                  BoxIterator bit( box );
                  for (bit.begin(); bit.ok(); ++bit)
                  {
                     pout() << bit() << "\t" << data( bit(), ncomp ) << std::endl;;
                  }
               }
               pout() << "--------------------------------------" << std::endl;
            }
       }
      }
   }
}


template <class T>
void BlockBaseRegister<T>::printRegisters() const
{
   pout() << ":::::::::::::::::::::::::::::::::::::::::::::::" << std::endl;
   pout() << "::  Source " << m_name << " Registers" << std::endl;
   pout() << ":::::::::::::::::::::::::::::::::::::::::::::::" << std::endl;
   print( m_srcRegister );

   pout() << ":::::::::::::::::::::::::::::::::::::::::::::::" << std::endl;
   pout() << "::  Destination " << m_name << " Registers" << std::endl;
   pout() << ":::::::::::::::::::::::::::::::::::::::::::::::" << std::endl;
   print( m_register );
}


#include "NamespaceFooter.H"

#endif
