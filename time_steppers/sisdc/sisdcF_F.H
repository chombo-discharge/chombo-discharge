#ifndef _SISDCF_F_H_
#define _SISDCF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C" {

#ifndef GUARDUPWIND_SOURCE 
#define GUARDUPWIND_SOURCE 
// Prototype for Fortran procedure upwind_source ...
//
void FORTRAN_NAME( UPWIND_SOURCE ,upwind_source )(
      CHFp_FRA1(dst)
      ,CHFp_CONST_FRA1(src)
      ,CHFp_CONST_FRA1(vel)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(dcalc) );

#define FORT_UPWIND_SOURCE FORTRAN_NAME( inlineUPWIND_SOURCE, inlineUPWIND_SOURCE)
#define FORTNT_UPWIND_SOURCE FORTRAN_NAME( UPWIND_SOURCE, upwind_source)

inline void FORTRAN_NAME(inlineUPWIND_SOURCE, inlineUPWIND_SOURCE)(
      CHFp_FRA1(dst)
      ,CHFp_CONST_FRA1(src)
      ,CHFp_CONST_FRA1(vel)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(dcalc) )
{
 CH_TIMELEAF("FORT_UPWIND_SOURCE");
 FORTRAN_NAME( UPWIND_SOURCE ,upwind_source )(
      CHFt_FRA1(dst)
      ,CHFt_CONST_FRA1(src)
      ,CHFt_CONST_FRA1(vel)
      ,CHFt_CONST_INT(dir)
      ,CHFt_BOX(dcalc) );
}
#endif  // GUARDUPWIND_SOURCE 

}

#endif
