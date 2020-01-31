
#ifndef _MACBNDRYDATA_H_
#define _MACBNDRYDATA_H_

//
// $Id: MacBndry.H,v 1.2 2007-05-24 23:54:48 tdsternberg Exp $
//

#include <InterpBndryData.H>

class MacBndry
    :
    public InterpBndryData
{
public:

    MacBndry ();

    MacBndry (const BoxArray& _grids,
              int             _ncomp,
              const Geometry& _geom);

    virtual void setBndryConds (const BCRec& phys_bc,
                                int          ratio);

    virtual void setBndryConds (const BCRec& phys_bc,
                                IntVect&     ratio,
                                int          comp=0);
};

#endif
