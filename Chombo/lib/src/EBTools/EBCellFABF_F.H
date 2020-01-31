#ifndef _EBCELLFABF_F_H_
#define _EBCELLFABF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSETCOVERED 
#define GUARDSETCOVERED 
// Prototype for Fortran procedure setcovered ...
//
void FORTRAN_NAME( SETCOVERED ,setcovered )(
      CHFp_CONST_FRA1(vel)
      ,CHFp_CONST_REAL(m)
      ,CHFp_CONST_FBA1(mask)
      ,CHFp_BOX(box) );

#define FORT_SETCOVERED FORTRAN_NAME( inlineSETCOVERED, inlineSETCOVERED)
#define FORTNT_SETCOVERED FORTRAN_NAME( SETCOVERED, setcovered)

inline void FORTRAN_NAME(inlineSETCOVERED, inlineSETCOVERED)(
      CHFp_CONST_FRA1(vel)
      ,CHFp_CONST_REAL(m)
      ,CHFp_CONST_FBA1(mask)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SETCOVERED");
 FORTRAN_NAME( SETCOVERED ,setcovered )(
      CHFt_CONST_FRA1(vel)
      ,CHFt_CONST_REAL(m)
      ,CHFt_CONST_FBA1(mask)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSETCOVERED 

}

#endif
