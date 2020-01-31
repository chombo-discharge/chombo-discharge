#ifndef _CELLTOEDGEF_F_H_
#define _CELLTOEDGEF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDCELLTOEDGE 
#define GUARDCELLTOEDGE 
// Prototype for Fortran procedure CELLTOEDGE ...
//
void FORTRAN_NAME( CELLTOEDGE ,celltoedge )(
      CHFp_CONST_FRA1(cellData)
      ,CHFp_FRA1(edgeData)
      ,CHFp_BOX(edgeBox)
      ,CHFp_CONST_INT(dir) );

#define FORT_CELLTOEDGE FORTRAN_NAME( inlineCELLTOEDGE, inlineCELLTOEDGE)
#define FORTNT_CELLTOEDGE FORTRAN_NAME( CELLTOEDGE, celltoedge)

inline void FORTRAN_NAME(inlineCELLTOEDGE, inlineCELLTOEDGE)(
      CHFp_CONST_FRA1(cellData)
      ,CHFp_FRA1(edgeData)
      ,CHFp_BOX(edgeBox)
      ,CHFp_CONST_INT(dir) )
{
 CH_TIMELEAF("FORT_CELLTOEDGE");
 FORTRAN_NAME( CELLTOEDGE ,celltoedge )(
      CHFt_CONST_FRA1(cellData)
      ,CHFt_FRA1(edgeData)
      ,CHFt_BOX(edgeBox)
      ,CHFt_CONST_INT(dir) );
}
#endif  // GUARDCELLTOEDGE 

}

#endif
