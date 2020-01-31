//
// $Id: Arena.cpp,v 1.1 2007-05-24 18:12:02 tdsternberg Exp $
//

#include <Arena.H>
#include <BoxLib.H>

Arena::~Arena () {}

std::size_t
Arena::align (std::size_t s)
{
    std::size_t x = s + sizeof(Word) - 1;
    x -= x%sizeof(Word);
    return x;
}
