//
// Implementation of neutralFacade.H -- the part whose implementation lives in
// Chombo land.
//

#include "neutralFacade.H"
#include "IntVect.H"   // Chombo
#include "Box.H"       //  "
#include "BoxLayout.H" //  "
#include "FArrayBox.H" //  "
#include "DisjointBoxLayout.H" //  "
#include "LevelData.H" //  "
#include "LoadBalance.H" //  "

///////////////////////// IntVect ////////////////////////////

BLfacade::IntVect::IntVect( const Chombo::IntVect& a_chIV )
  : m_serialization(0),
    m_pimpl(0)
{
  m_serialization = new int[CH_SPACEDIM];
  D_EXPR( m_serialization[0] = a_chIV.vect[0],
          m_serialization[1] = a_chIV.vect[1],
          m_serialization[2] = a_chIV.vect[2] );
  initPimpl();
  CH_assert( m_pimpl );
}


BLfacade::IntVect::operator Chombo::IntVect() const
{
  return Chombo::IntVect( m_serialization );
}

///////////////////////// Box ////////////////////////////

BLfacade::Box::Box( const Chombo::Box& a_chBox )
  : m_lo( a_chBox.smallEnd() ),
    m_hi( a_chBox.bigEnd() ),
    m_pimpl(0)
{
  initPimpl();
  CH_assert( m_pimpl );
}

BLfacade::Box::operator Chombo::Box() const
{
  return Chombo::Box( Chombo::IntVect(m_lo),
                      Chombo::IntVect(m_hi) );
}


///////////////////////// BoxArray ////////////////////////////

/**
 *  A CCSE::BoxArray has three data members, but two of them -- m_crsn and
 *  m_hash -- look like they're generated on the fly.  So we need only provide
 *  the array of actual boxes.
*/
BLfacade::BoxArray::BoxArray( const Chombo::BoxLayout& a_chBL )
  : m_numBoxes( a_chBL.size() ),
    m_boxes(0),
    m_pimpl(0),
    m_dbl(0)
{
  Chombo::Vector<Chombo::Box> chBoxVector(a_chBL.boxArray());
  typedef Box* pBox;
  m_boxes = new pBox[m_numBoxes];  // Cuz we don't have a Box copy ctor.
  for ( int i=0; i<m_numBoxes; ++i )
  {
    m_boxes[i] = new Box( chBoxVector[i] );
  }

  initPimpl();
  CH_assert( m_pimpl );
}


Chombo::Box
BLfacade::BoxArray::operator[](int b) const
{
  return *((m_boxes)[b]);
}


///////////////////////// MultiFab ////////////////////////////

/**
 The BoxArray arg must not be destroyed during the lifetime of *this.  It's not
 const because we assign it a DisjointBoxLayout pointer in this->pLevelData().
 (That's to ensure that LevelData's constructed by MultiFab's built from the
 same BoxArray, hold equivalent DisjointBoxLayouts.)

 If arg procIDs is not null, then if and when we call this->pLevelData(),
 the DisjointBoxLayout under that LevelData (and to which m_boxarray will hold
 a pointer) gets assigned those procIDs.

 It's an error to construct two MultiFabs from the same BoxArray, but assign
 procIDs in one but not the other constructor, or assign procIDs each time but
 not the same procIDs.
*/
BLfacade::MultiFab::MultiFab( BLfacade::BoxArray& ba,
                              int ncomps,
                              int nghosts,
                              const std::vector<int>* procIDs /*=0*/ )
  : m_boxarray( ba ),
    m_ncomps( ncomps ),
    m_nghosts( nghosts ),
    m_pLevelData( 0 )
{
  initPimpl( procIDs );
  CH_assert( m_pimpl );
}
 

/** Called from destructor to do things that require Chombo headers. */
void
BLfacade::MultiFab::chomboDelete()
{
  delete m_pLevelData;
}

namespace BLfacade
{
struct BoxNDit
{
  BoxNDit( Box const& a_box, Chombo::DataIterator const& a_dit )
   : box(a_box), dit(a_dit) {}
  Box box;
  Chombo::DataIterator dit;
  bool operator<( BoxNDit const& ) const;
};


bool BoxNDit::operator<( BoxNDit const& that ) const
{
  return this->box < that.box;
}

}


/** Creates a LevelData that corresponds to this MultiFab.
    
    MultiFab instance owns the return value.
    
    It's an error to call this function twice with different arguments.
 */
Chombo::LevelData<Chombo::FArrayBox>*
BLfacade::MultiFab::pLevelData( const Chombo::Box&     a_domain,
                                const Chombo::IntVect& a_periodicDirs )
{
  if ( !m_pLevelData )
  {
    Chombo::Vector<Chombo::Box> boxes( m_boxarray.numBoxes() );
    for ( int b=0; b<m_boxarray.numBoxes(); ++b )
    {
       boxes[b] = m_boxarray[b];
    }

    // All MultiFabs constructed from the same BLfacade::BoxArray must return
    // LevelData's whose DisjointBoxLayouts are equivalent.  Thus, we let
    // m_boxarray own a DisjointBoxLayout, which we construct the first time
    // it's needed, and keep cached.  Of course, such MultiFabs need to have
    // the same processor assignments or none of this makes sense.  According
    // to Mike Lijewski, two CCSE::MultiFabs constructed from the same
    // CCSE::BoxArray are guaranteed to have the same procIDs.  That means the
    // only way to obtain two BLfacade::MultiFabs, constructed from the same
    // BLfacade::BoxArray, with different procIDs is if we assign procIDs in
    // the BLfacade::MultiFabs' constructors, and either assign different
    // procIDs each time, or assign them to one BLfacade::MultiFab but not the
    // other (in which case the other gets what CCSE calculates, and that might
    // not be what we assigned to the first BLfacade::MultiFab).  Errors like
    // these will be caught by the assertion 19 lines down from here (the one
    // that mentions "procIDs()").
    if ( !m_boxarray.dbl() )
    {
      // Need to assign processor assignments as in the MultiFab, since we're
      // about to set pointers to the FABs, and those pointers are only valid
      // on the processors that actually have those FABs.
      Chombo::Vector<int> procIDs( this->procIDs() );
      // Processor info is in FabArray::DistributionMap().  Of course you won't
      // be able to work with that in this file.

      bool bPeriodic[CH_SPACEDIM];
      for ( int d=0;d<CH_SPACEDIM;++d ) bPeriodic[d] = a_periodicDirs[d];
      Chombo::ProblemDomain dom( a_domain, bPeriodic );

      m_boxarray.dbl( new Chombo::DisjointBoxLayout( boxes, procIDs, dom ) );
      // Now assign the these procIDs to m_boxarray?
    }

    // It's an error to construct two MultiFabs off the same BoxArray, but let
    // (or force) those MultiFabs to have different processor assignments:
    CH_assert( m_boxarray.dbl()->procIDs().stdVector() == this->procIDs() );

    // Make sure args are as they were the first time this function was called.
    CH_assert( a_domain == m_boxarray.dbl()->physDomain().domainBox() );
#ifdef CH_DEBUG
    {
    Chombo::IntVect periodicDirsCheck;
    for ( int d=0;d<CH_SPACEDIM;++d )
    {
        periodicDirsCheck[d] = m_boxarray.dbl()->physDomain.isPeriodic(d);
    }
    CH_assert( a_periodicDirs == periodicDirsCheck );
    }
#endif

    //
    // Match up our boxes to their boxes.
    //
    struct ChomboPtrLayout
    {
      Chombo::DataIterator dit;
      Box          box;
      bool operator<( ChomboPtrLayout const& that )
        { return this->box < that.box; }
    };      

    // We'll build the LevelData using a FABAliasDataFactory.  But to do that,
    // first we need to construct a LayoutData<Real*> with pointers set to the
    // addresses of the CCSE::MultiFab's data.
    // The box-to-processor assignment is the same between
    // the MultiFab and the LevelData, but the order of boxes on each processor
    // is probably different.  Make two arrays of FAB data pointers -- one for
    // the MultiFab, one for the LevelData -- and sort them by box low corner.
    //
    // (Example of how FABAliasDataFactory is used: see AMR-HDF5.cpp in the
    // AMR-CCA CVS module.)
    //

    Chombo::LayoutData<Real*> lyd( *m_boxarray.dbl() );
    std::vector<BoxNDptr>  ccsePtrLayout( getBoxNDptrs() );
    std::vector<BoxNDit> chomboPtrLayout( getBoxNDits(lyd) );
    std::sort( ccsePtrLayout.begin(), ccsePtrLayout.end() );
    std::sort( chomboPtrLayout.begin(), chomboPtrLayout.end() );
    CH_assert( ccsePtrLayout.size() == chomboPtrLayout.size() );

    for ( unsigned b=0; b<chomboPtrLayout.size(); ++b )
    {
      CH_assert( chomboPtrLayout[b].box == ccsePtrLayout[b].box );
      lyd[ chomboPtrLayout[b].dit() ] = ccsePtrLayout[b].dptr;
    }
    Chombo::FABAliasDataFactory factory( lyd );
    m_pLevelData = new Chombo::LevelData<Chombo::FArrayBox>(
      *m_boxarray.dbl(), m_ncomps, Chombo::IntVect::Unit*m_nghosts, factory );
  }

  return m_pLevelData;
}


/**
 * Private member that helps with creation of LevelData.
*/
std::vector<BLfacade::BoxNDit>
BLfacade::MultiFab::getBoxNDits( const Chombo::LayoutData<Real*>& lyd ) const
{
  std::vector<BLfacade::BoxNDit> result;
  // FIXME: reserve()!
  for ( Chombo::DataIterator dit(lyd.dataIterator()); dit.ok(); ++dit )
  {
    result.push_back( BLfacade::BoxNDit( lyd.box(dit()), dit ) );
  }

  return result;
}
