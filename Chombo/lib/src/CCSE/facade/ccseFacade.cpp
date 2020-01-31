//
// Implementation of neutralFacade.H -- the part whose implementation lives in
// CCSE land.
//

#include <iostream>

#include "neutralFacade.H"

#include "IntVect.H"             // CCSE
#include "Array.H"               // CCSE
#include "Box.H"                 // CCSE
#include "BoxArray.H"            // CCSE
#include "MultiFab.H"            // CCSE
#include "ParallelDescriptor.H"  // CCSE
#include "DistributionMapping.H" // CCSE

namespace CCSE
{
  typedef ::DistributionMapping DistributionMapping;
  typedef ::Array<int> Array;
}

//
// BoxLib::Initialize() and Finalize()
//
static bool s_boxlibInitialized;
void BLfacade::initializeBoxLib(int argc, char** argv)
{
  BoxLib::Initialize(argc,argv);
  s_boxlibInitialized = true;
}

void BLfacade::finalizeBoxLib()
{
  BoxLib::Finalize();
}


///////////////////////// IntVect ////////////////////////////

BLfacade::IntVect::IntVect( const CCSE::IntVect& a_ccIV )
{
  m_serialization = new int[CH_SPACEDIM];
  D_EXPR( m_serialization[0] = a_ccIV[0],
          m_serialization[1] = a_ccIV[1],
          m_serialization[2] = a_ccIV[2] );
  initPimpl();
}

/**
  Private method.
  Should be called from last line of constructor.
*/
void
BLfacade::IntVect::initPimpl()
{
  m_pimpl = new CCSE::IntVect( m_serialization );
}


BLfacade::IntVect::~IntVect()
{
  delete[] m_serialization;
  delete m_pimpl;
}


BLfacade::IntVect::operator CCSE::IntVect() const
{
  return *m_pimpl;
}


std::ostream&
BLfacade::IntVect::streamOut( std::ostream& out ) const
{
  out << (*m_pimpl);
  return out;
}


///////////////////////// Box ////////////////////////////

BLfacade::Box::Box( const CCSE::Box& b )
  : m_lo( b.smallEnd() ),
    m_hi( b.bigEnd() )
{
  initPimpl();
}

/**
  Private method.
  Should be called from last line of constructor.
*/
void
BLfacade::Box::initPimpl()
{
  m_pimpl = new CCSE::Box( m_lo, m_hi );
}

BLfacade::Box::~Box()
{
  delete m_pimpl;
}


BLfacade::Box::operator CCSE::Box() const
{
  return *m_pimpl;
}


std::ostream&
BLfacade::Box::streamOut( std::ostream& out ) const
{
  out << (*m_pimpl);
  return out;
}

bool BLfacade::Box::operator<( BLfacade::Box const & that ) const
{
  return m_pimpl->smallEnd() < that.m_pimpl->smallEnd();
}

bool BLfacade::Box::operator==( BLfacade::Box const& that ) const
{
  return *(this->m_pimpl) == *(that.m_pimpl);
}


BLfacade::Box::Box( const BLfacade::Box& that )
  : m_lo( that.m_lo ),
    m_hi( that.m_hi )
{
  m_pimpl = new CCSE::Box( *(that.m_pimpl) );
}


BLfacade::Box& BLfacade::Box::operator=( const BLfacade::Box& that )
{
  if ( this != &that )
  {
    m_lo = that.m_lo;
    m_hi = that.m_hi;
    m_pimpl = new CCSE::Box( *(that.m_pimpl) );
  }
  return *this;
}

///////////////////////// BoxArray ////////////////////////////


/**
  Private method.
  Should be called from last line of constructor.
*/
void
BLfacade::BoxArray::initPimpl()
{
  std::vector<CCSE::Box> boxvector( m_numBoxes );
  for ( int i=0; i<m_numBoxes; ++i )
  {
    boxvector[i] = *(m_boxes[i]);
  }
  m_pimpl = new CCSE::BoxArray( &boxvector[0], m_numBoxes );
}


BLfacade::BoxArray::~BoxArray()
{
  for ( int i=0;i<m_numBoxes;++i ) delete m_boxes[i];
  delete[] m_boxes;
  delete m_pimpl;
}


BLfacade::BoxArray::operator CCSE::BoxArray const&() const
{
  return *m_pimpl;
}


std::ostream&
BLfacade::BoxArray::streamOut( std::ostream& out ) const
{
  out << (*m_pimpl);
  return out;
}


///////////////////////// MultiFab ////////////////////////////

/**
  Private method.
  Should be called from last line of constructor.
  Assign processors as per procIDs, if not NULL.
*/
void
BLfacade::MultiFab::initPimpl( const std::vector<int>* procIDs )
{
  if ( !s_boxlibInitialized )
  {
    std::cerr << "You should put the following line really, really early "
                 "in your main():\n"
                 "BLfacade::initializeBoxLib(argc,argv);\n";
    exit(1);
  }

  m_pimpl = new CCSE::MultiFab;

  if ( procIDs ) // Assign 'em.  (Got this code from Mike Lijewski.)
  {
    CCSE::Array pmap( m_boxarray.numBoxes()+1 );
    for ( int b=0; b<m_boxarray.numBoxes(); ++b )
    {
      pmap[b] = (*procIDs)[b];
    }   
    pmap[m_boxarray.numBoxes()] = ParallelDescriptor::MyProc(); //sentinel value
    CCSE::DistributionMapping distMapping( pmap );

    m_pimpl->define( m_boxarray, m_ncomps, m_nghosts, distMapping, Fab_allocate );
  } else
  {
    m_pimpl->define( m_boxarray, m_ncomps, m_nghosts, Fab_allocate );
  }
}


std::vector<int>
BLfacade::MultiFab::procIDs() const
{
  CCSE::Array distMap( m_pimpl->DistributionMap().ProcessorMap() );
  std::vector<int> result( m_boxarray.numBoxes() );
  for ( int b=0; b<distMap.size()-1; ++b ) // -1 cuz last value is a sentinel
  {
    result[b] = distMap[b];
  }
  return result;
}


BLfacade::MultiFab::~MultiFab()
{
  delete m_pimpl;
  chomboDelete();
}


std::vector<BLfacade::BoxNDptr>
BLfacade::MultiFab::getBoxNDptrs() const
{
  std::vector<BLfacade::BoxNDptr> result;
  // FIXME: reserve()!

  for ( CCSE::MFIter mfi(*m_pimpl); mfi.isValid(); ++mfi )
  {
    CCSE::Box boxNoGhosts( (*m_pimpl)[mfi].box() );
    boxNoGhosts.grow(-m_nghosts);
    result.push_back( BLfacade::BoxNDptr(boxNoGhosts, (*m_pimpl)[mfi].dataPtr()) );
    // We store boxes net of ghost cells; that's the Chombo convention and it'll
    // be the convention here too.  (CCSE::MultiFab stores its boxes gross of
    // ghost cells.)
  }

  return result;
}


BLfacade::MultiFab::operator CCSE::MultiFab&() const
{
  return *m_pimpl;
}
