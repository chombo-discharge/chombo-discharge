//
// $Id: BArena.cpp,v 1.1 2007-05-24 18:12:02 tdsternberg Exp $
//
#include <BArena.H>
#include <BoxLib.H>

void*
BArena::alloc (std::size_t _sz)
{
    return ::operator new(_sz);
}

void
BArena::free (void* pt)
{
    ::operator delete(pt);
}
