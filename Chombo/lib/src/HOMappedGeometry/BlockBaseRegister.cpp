#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "BlockBaseRegister.H"
#include "BlockBaseRegisterF_F.H"

#include "NamespaceHeader.H"

template < >
inline
void BlockBaseRegister<FArrayBox>::copyToDstRegister()
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
                  
                  LevelData<FArrayBox>& dstRegister( m_register[blockBdryIndex] );
                  const LevelData<FArrayBox>& srcRegister( m_srcRegister[srcBlockBdryIndex] );
                  
                  FArrayBox& dst( dstRegister[regDataIndex] );
                  const FArrayBox& src( srcRegister[srcDataIndex] );
                  
                  const Box& srcWindow( src.box() );
                  const Box& regWindow( itrans.transformBack( srcWindow ) );
                  const int ncomp( dst.nComp() );
                  
                  // copy the transformed box data into the register
                  FORT_TRANSFILLREGISTER(CHF_FRA(dst),
                                         CHF_CONST_FRA(src),
                                         CHF_BOX(regWindow),
                                         CHF_CONST_INTVECT(itrans.getPermutation()),
                                         CHF_CONST_INTVECT(itrans.getSign()),
                                         CHF_CONST_INTVECT(itrans.getTranslation()),
                                         CHF_CONST_INT(ncomp));
               }
            }
            m_register[blockBdryIndex].exchange();
         }
      }
   }
}

#include "NamespaceFooter.H"
