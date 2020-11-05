#ifndef _EBCONSTANTCFINTERPF_F_H_
#define _EBCONSTANTCFINTERPF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDREGCONSTANTINTERP 
#define GUARDREGCONSTANTINTERP 
// Prototype for Fortran procedure regconstantinterp ...
//
void FORTRAN_NAME( REGCONSTANTINTERP ,regconstantinterp )(
      CHFp_FRA1(phi)
      ,CHFp_BOX(validbox)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(nghost)
      ,CHFp_CONST_INT(side) );

#define FORT_REGCONSTANTINTERP FORTRAN_NAME( inlineREGCONSTANTINTERP, inlineREGCONSTANTINTERP)
#define FORTNT_REGCONSTANTINTERP FORTRAN_NAME( REGCONSTANTINTERP, regconstantinterp)

inline void FORTRAN_NAME(inlineREGCONSTANTINTERP, inlineREGCONSTANTINTERP)(
      CHFp_FRA1(phi)
      ,CHFp_BOX(validbox)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(nghost)
      ,CHFp_CONST_INT(side) )
{
 CH_TIMELEAF("FORT_REGCONSTANTINTERP");
 FORTRAN_NAME( REGCONSTANTINTERP ,regconstantinterp )(
      CHFt_FRA1(phi)
      ,CHFt_BOX(validbox)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(nghost)
      ,CHFt_CONST_INT(side) );
}
#endif  // GUARDREGCONSTANTINTERP 

}

#endif
