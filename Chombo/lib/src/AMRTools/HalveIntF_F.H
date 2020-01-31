#ifndef _HALVEINTF_F_H_
#define _HALVEINTF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDHALVEINT 
#define GUARDHALVEINT 
// Prototype for Fortran procedure halveint ...
//
void FORTRAN_NAME( HALVEINT ,halveint )(
      CHFp_FIA(arr)
      ,CHFp_BOX(bx) );

#define FORT_HALVEINT FORTRAN_NAME( inlineHALVEINT, inlineHALVEINT)
#define FORTNT_HALVEINT FORTRAN_NAME( HALVEINT, halveint)

inline void FORTRAN_NAME(inlineHALVEINT, inlineHALVEINT)(
      CHFp_FIA(arr)
      ,CHFp_BOX(bx) )
{
 CH_TIMELEAF("FORT_HALVEINT");
 FORTRAN_NAME( HALVEINT ,halveint )(
      CHFt_FIA(arr)
      ,CHFt_BOX(bx) );
}
#endif  // GUARDHALVEINT 

}

#endif
