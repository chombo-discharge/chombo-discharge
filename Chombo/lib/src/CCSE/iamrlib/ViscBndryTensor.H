#ifndef _VISCBNDRYTENSOR_H_
#define _VISCBNDRYTENSOR_H_

//
// $Id: ViscBndryTensor.H,v 1.2 2007-05-24 23:54:48 tdsternberg Exp $
//

#include <MCInterpBndryData.H>
#include <BC_TYPES.H>
#include <MCLinOp.H>

class ViscBndryTensor
    :
    public MCInterpBndryData
{
public:
    ViscBndryTensor () {}

    ViscBndryTensor (const BoxArray& grids,
                     const Geometry& geom)
        :
        MCInterpBndryData(grids,MCLinOp::bcComponentsNeeded(),geom) {}

    virtual void setBndryConds (const BCRec& phys_bc,
                                int          ratio,
                                int          comp=0);

    void setHomogValues (const Array<BCRec>& phys_bc,
                         int                 ratio);

};

#endif /*_VISCBNDRYTENSOR_H_*/
