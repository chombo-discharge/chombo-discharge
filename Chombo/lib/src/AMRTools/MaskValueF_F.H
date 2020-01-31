#ifndef _MASKVALUEF_F_H_
#define _MASKVALUEF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDMASKVALUE 
#define GUARDMASKVALUE 
// Prototype for Fortran procedure maskvalue ...
//
void FORTRAN_NAME( MASKVALUE ,maskvalue )(
      CHFp_FRA1(mask)
      ,CHFp_CONST_FIA1(test)
      ,CHFp_BOX(bx)
      ,CHFp_CONST_INT(val)
      ,CHFp_CONST_INT(onoff) );

#define FORT_MASKVALUE FORTRAN_NAME( inlineMASKVALUE, inlineMASKVALUE)
#define FORTNT_MASKVALUE FORTRAN_NAME( MASKVALUE, maskvalue)

inline void FORTRAN_NAME(inlineMASKVALUE, inlineMASKVALUE)(
      CHFp_FRA1(mask)
      ,CHFp_CONST_FIA1(test)
      ,CHFp_BOX(bx)
      ,CHFp_CONST_INT(val)
      ,CHFp_CONST_INT(onoff) )
{
 CH_TIMELEAF("FORT_MASKVALUE");
 FORTRAN_NAME( MASKVALUE ,maskvalue )(
      CHFt_FRA1(mask)
      ,CHFt_CONST_FIA1(test)
      ,CHFt_BOX(bx)
      ,CHFt_CONST_INT(val)
      ,CHFt_CONST_INT(onoff) );
}
#endif  // GUARDMASKVALUE 

}

#endif
