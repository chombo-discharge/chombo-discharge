#ifndef _DIVFLOWIBCF_F_H_
#define _DIVFLOWIBCF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDVECLATLONDIVFLOW 
#define GUARDVECLATLONDIVFLOW 
// Prototype for Fortran procedure veclatlondivflow ...
//
void FORTRAN_NAME( VECLATLONDIVFLOW ,veclatlondivflow )(
      CHFp_FRA(vecRLL)
      ,CHFp_CONST_FRA(ptRLL)
      ,CHFp_CONST_REAL(k)
      ,CHFp_CONST_REAL(time)
      ,CHFp_CONST_REAL(period) );

#define FORT_VECLATLONDIVFLOW FORTRAN_NAME( inlineVECLATLONDIVFLOW, inlineVECLATLONDIVFLOW)
#define FORTNT_VECLATLONDIVFLOW FORTRAN_NAME( VECLATLONDIVFLOW, veclatlondivflow)

inline void FORTRAN_NAME(inlineVECLATLONDIVFLOW, inlineVECLATLONDIVFLOW)(
      CHFp_FRA(vecRLL)
      ,CHFp_CONST_FRA(ptRLL)
      ,CHFp_CONST_REAL(k)
      ,CHFp_CONST_REAL(time)
      ,CHFp_CONST_REAL(period) )
{
 CH_TIMELEAF("FORT_VECLATLONDIVFLOW");
 FORTRAN_NAME( VECLATLONDIVFLOW ,veclatlondivflow )(
      CHFt_FRA(vecRLL)
      ,CHFt_CONST_FRA(ptRLL)
      ,CHFt_CONST_REAL(k)
      ,CHFt_CONST_REAL(time)
      ,CHFt_CONST_REAL(period) );
}
#endif  // GUARDVECLATLONDIVFLOW 

}

#endif
