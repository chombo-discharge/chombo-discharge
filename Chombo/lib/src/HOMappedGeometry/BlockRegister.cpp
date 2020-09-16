#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "BlockRegister.H"
#include "BoxIterator.H"
#include "CH_Timer.H"

#include "NamespaceHeader.H"


inline int getBlockBdryIndex( const int& a_dir, const Side::LoHiSide& a_side );


BlockRegister::BlockRegister( const MultiBlockCoordSys* a_mblock,
                              const DisjointBoxLayout& a_grids,
                              const int& a_ghost )
  :
  m_verbosity(0),
  m_isOpen(true)
{
  define( a_mblock, a_grids, a_ghost );
}


BlockRegister::~BlockRegister()
{
   //do nothing for now
}


void BlockRegister::define( const MultiBlockCoordSys* a_mblock,
                            const DisjointBoxLayout& a_grids,
                            const int& a_ghost )
{
   CH_TIME("BlockRegister::define");

   m_fluxRegister.define( a_mblock, a_grids, 0 );
   m_fluxRegister.setName( "Flux" );

   m_auxRegister.define( a_mblock, a_grids, a_ghost );
   m_auxRegister.setName( "Aux" );
}


void BlockRegister::setVerbose( int a_verbose )
{
   m_fluxRegister.setVerbose( a_verbose );
   m_auxRegister.setVerbose( a_verbose );
}


void BlockRegister::storeFlux( const FArrayBox& a_flux,
                               const DataIndex& a_dataIndex,
                               int a_dir,
                               Side::LoHiSide a_side )
{
   CH_TIME("BlockRegister::storeFlux");
   if ( isNotClosed() )
   {
      m_fluxRegister.store( a_flux, a_dataIndex, a_dir, a_side );
   }
   else
   {
      MayDay::Error(" BlockRegister::storeFlux:  BlockRegister is closed!");
   }
}


void BlockRegister::incrementFlux( const FArrayBox& a_flux,
                                   const DataIndex& a_dataIndex,
                                   int a_dir,
                                   Side::LoHiSide a_side )
{
   CH_TIME("BlockRegister::incrementFlux");
   if ( isNotClosed() )
   {
      m_fluxRegister.increment( a_flux, a_dataIndex, a_dir, a_side );
   }
   else
   {
      MayDay::Error(" BlockRegister::incrementFlux:  BlockRegister is closed!");
   }
}


void BlockRegister::storeAux( const FArrayBox& a_data,
                              const DataIndex& a_dataIndex,
                              int a_dir,
                              Side::LoHiSide a_side )
{
   CH_TIME("BlockRegister::storeAux");
   if ( isNotClosed() )
   {
      m_auxRegister.store( a_data, a_dataIndex, a_dir, a_side );
   }
   else
   {
      MayDay::Error(" BlockRegister::storeAux:  BlockRegister is closed!");
   }
}


void BlockRegister::incrementAux( const FArrayBox& a_data,
                                  const DataIndex& a_dataIndex,
                                  int a_dir,
                                  Side::LoHiSide a_side )
{
   CH_TIME("BlockRegister::storeAux");
   if ( isNotClosed() )
   {
      m_auxRegister.increment( a_data, a_dataIndex, a_dir, a_side );
   }
   else
   {
      MayDay::Error(" BlockRegister::storeAux:  BlockRegister is closed!");
   }
}


void BlockRegister::close()
{
   CH_TIME("BlockRegister::close");
   m_fluxRegister.exchange();
   m_auxRegister.exchange();
   m_isOpen = false;
}


void BlockRegister::setToZero( const int a_n_comp )
{
   CH_TIME("BlockRegister::resetToZero");
   m_fluxRegister.zeroRegister( a_n_comp );
   m_auxRegister.zeroRegister( a_n_comp );
   m_isOpen = true;
}


bool BlockRegister::hasInterface( const DataIndex& a_dataIndex,
                                  int a_dir,
                                  Side::LoHiSide a_side ) const
{
   CH_TIME("BlockRegister::hasInterface");
   // JAFH: We assume the a_dataIndex is valid; should we check?
   return m_fluxRegister.hasInterface( a_dataIndex, a_dir, a_side );
}


void BlockRegister::getFlux( FArrayBox& a_flux,
                             const DataIndex& a_dataIndex,
                             int a_dir,
                             Side::LoHiSide a_side,
                             Side::LoHiSide a_sideData ) const
{
   CH_TIME("BlockRegister::getFlux");
   if ( isClosed() )
   {
      m_fluxRegister.fill( a_flux, a_dataIndex, a_dir, a_side, a_sideData );
   }
   else
   {
      MayDay::Error(" BlockRegister::getFlux:  BlockRegister was not closed!");
   }
}


void BlockRegister::getAux( FArrayBox& a_data,
                            const DataIndex& a_dataIndex,
                            int a_dir,
                            Side::LoHiSide a_side,
                            Side::LoHiSide a_sideData ) const
{
   CH_TIME("BlockRegister::getAux");

   if ( isClosed() )
   {
      m_auxRegister.fill( a_data, a_dataIndex, a_dir, a_side, a_sideData );
   }
   else
   {
      MayDay::Error(" BlockRegister::getAux:  BlockRegister was not closed!");
   }
}


#include "NamespaceFooter.H"
