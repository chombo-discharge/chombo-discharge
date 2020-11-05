#ifndef _EBMGINTERPF_F_H_
#define _EBMGINTERPF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDREGPROLONG 
#define GUARDREGPROLONG 
// Prototype for Fortran procedure regprolong ...
//
void FORTRAN_NAME( REGPROLONG ,regprolong )(
      CHFp_FRA1(phithislevel)
      ,CHFp_CONST_FRA1(correctcoarse)
      ,CHFp_BOX(coarbox)
      ,CHFp_BOX(refbox)
      ,CHFp_CONST_INT(reftocoar) );

#define FORT_REGPROLONG FORTRAN_NAME( inlineREGPROLONG, inlineREGPROLONG)
#define FORTNT_REGPROLONG FORTRAN_NAME( REGPROLONG, regprolong)

inline void FORTRAN_NAME(inlineREGPROLONG, inlineREGPROLONG)(
      CHFp_FRA1(phithislevel)
      ,CHFp_CONST_FRA1(correctcoarse)
      ,CHFp_BOX(coarbox)
      ,CHFp_BOX(refbox)
      ,CHFp_CONST_INT(reftocoar) )
{
 CH_TIMELEAF("FORT_REGPROLONG");
 FORTRAN_NAME( REGPROLONG ,regprolong )(
      CHFt_FRA1(phithislevel)
      ,CHFt_CONST_FRA1(correctcoarse)
      ,CHFt_BOX(coarbox)
      ,CHFt_BOX(refbox)
      ,CHFt_CONST_INT(reftocoar) );
}
#endif  // GUARDREGPROLONG 

#ifndef GUARDPROLONGADDSLOPE 
#define GUARDPROLONGADDSLOPE 
// Prototype for Fortran procedure prolongaddslope ...
//
void FORTRAN_NAME( PROLONGADDSLOPE ,prolongaddslope )(
      CHFp_FRA1(phithislevel)
      ,CHFp_CONST_FRA1(correctcoarse)
      ,CHFp_BOX(coarbox)
      ,CHFp_BOX(refbox)
      ,CHFp_INT(idir)
      ,CHFp_REAL(dxf)
      ,CHFp_REAL(dxc)
      ,CHFp_CONST_INT(reftocoar) );

#define FORT_PROLONGADDSLOPE FORTRAN_NAME( inlinePROLONGADDSLOPE, inlinePROLONGADDSLOPE)
#define FORTNT_PROLONGADDSLOPE FORTRAN_NAME( PROLONGADDSLOPE, prolongaddslope)

inline void FORTRAN_NAME(inlinePROLONGADDSLOPE, inlinePROLONGADDSLOPE)(
      CHFp_FRA1(phithislevel)
      ,CHFp_CONST_FRA1(correctcoarse)
      ,CHFp_BOX(coarbox)
      ,CHFp_BOX(refbox)
      ,CHFp_INT(idir)
      ,CHFp_REAL(dxf)
      ,CHFp_REAL(dxc)
      ,CHFp_CONST_INT(reftocoar) )
{
 CH_TIMELEAF("FORT_PROLONGADDSLOPE");
 FORTRAN_NAME( PROLONGADDSLOPE ,prolongaddslope )(
      CHFt_FRA1(phithislevel)
      ,CHFt_CONST_FRA1(correctcoarse)
      ,CHFt_BOX(coarbox)
      ,CHFt_BOX(refbox)
      ,CHFt_INT(idir)
      ,CHFt_REAL(dxf)
      ,CHFt_REAL(dxc)
      ,CHFt_CONST_INT(reftocoar) );
}
#endif  // GUARDPROLONGADDSLOPE 

}

#endif
