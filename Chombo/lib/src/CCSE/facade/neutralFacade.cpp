//
// Implementation of *some* of what's declared in neutralFacade.H.  Other
// parts of it are implemented in ccseFacade.cpp and chomboFacade.cpp.
//

#include <iostream>
#include <cassert>
#include "neutralFacade.H"

std::ostream& BLfacade::Streamable::streamOut( std::ostream& out ) const
{
  out << "Unimplemented.";
  return out;
}

std::ostream& operator<<(std::ostream& out, const BLfacade::Streamable& obj)
{
  return obj.streamOut( out );
}


bool BLfacade::BoxNDptr::operator<( BLfacade::BoxNDptr const& that ) const
{
  return this->box < that.box;
}


BLfacade::IntVect::IntVect( const BLfacade::IntVect& that )
{
  m_serialization = new int[CH_SPACEDIM];
  memcpy( m_serialization, that.m_serialization, CH_SPACEDIM*sizeof(int) );
  initPimpl();
  assert( m_pimpl );
}

BLfacade::IntVect&
BLfacade::IntVect::operator=( const BLfacade::IntVect& that )
{
  if ( this != &that )
  {
    m_serialization = new int[CH_SPACEDIM];
    memcpy( m_serialization, that.m_serialization, CH_SPACEDIM*sizeof(int) );
    initPimpl();
    assert( m_pimpl );
  }
  return *this;
}
