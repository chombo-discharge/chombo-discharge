#ifndef _NODEAVERAGEF_F_H_
#define _NODEAVERAGEF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDNODEAVERAGE 
#define GUARDNODEAVERAGE 
// Prototype for Fortran procedure NODEAVERAGE ...
//
void FORTRAN_NAME( NODEAVERAGE ,nodeaverage )(
      CHFp_FRA(coarse)
      ,CHFp_CONST_FRA(fine)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(ref_ratio)
      ,CHFp_CONST_FRA(weight) );

#define FORT_NODEAVERAGE FORTRAN_NAME( inlineNODEAVERAGE, inlineNODEAVERAGE)
#define FORTNT_NODEAVERAGE FORTRAN_NAME( NODEAVERAGE, nodeaverage)

inline void FORTRAN_NAME(inlineNODEAVERAGE, inlineNODEAVERAGE)(
      CHFp_FRA(coarse)
      ,CHFp_CONST_FRA(fine)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(ref_ratio)
      ,CHFp_CONST_FRA(weight) )
{
 CH_TIMELEAF("FORT_NODEAVERAGE");
 FORTRAN_NAME( NODEAVERAGE ,nodeaverage )(
      CHFt_FRA(coarse)
      ,CHFt_CONST_FRA(fine)
      ,CHFt_BOX(b)
      ,CHFt_CONST_INT(ref_ratio)
      ,CHFt_CONST_FRA(weight) );
}
#endif  // GUARDNODEAVERAGE 

#ifndef GUARDNODEAVERAGEPOINT 
#define GUARDNODEAVERAGEPOINT 
// Prototype for Fortran procedure NODEAVERAGEPOINT ...
//
void FORTRAN_NAME( NODEAVERAGEPOINT ,nodeaveragepoint )(
      CHFp_FRA(coarse)
      ,CHFp_CONST_FRA(fine)
      ,CHFp_CONST_INTVECT(pcrse)
      ,CHFp_CONST_INT(ref_ratio)
      ,CHFp_CONST_FRA(weight) );

#define FORT_NODEAVERAGEPOINT FORTRAN_NAME( inlineNODEAVERAGEPOINT, inlineNODEAVERAGEPOINT)
#define FORTNT_NODEAVERAGEPOINT FORTRAN_NAME( NODEAVERAGEPOINT, nodeaveragepoint)

inline void FORTRAN_NAME(inlineNODEAVERAGEPOINT, inlineNODEAVERAGEPOINT)(
      CHFp_FRA(coarse)
      ,CHFp_CONST_FRA(fine)
      ,CHFp_CONST_INTVECT(pcrse)
      ,CHFp_CONST_INT(ref_ratio)
      ,CHFp_CONST_FRA(weight) )
{
 CH_TIMELEAF("FORT_NODEAVERAGEPOINT");
 FORTRAN_NAME( NODEAVERAGEPOINT ,nodeaveragepoint )(
      CHFt_FRA(coarse)
      ,CHFt_CONST_FRA(fine)
      ,CHFt_CONST_INTVECT(pcrse)
      ,CHFt_CONST_INT(ref_ratio)
      ,CHFt_CONST_FRA(weight) );
}
#endif  // GUARDNODEAVERAGEPOINT 

#ifndef GUARDNODEAVERAGE_GETWEIGHTS 
#define GUARDNODEAVERAGE_GETWEIGHTS 
// Prototype for Fortran procedure NODEAVERAGE_GETWEIGHTS ...
//
void FORTRAN_NAME( NODEAVERAGE_GETWEIGHTS ,nodeaverage_getweights )(
      CHFp_FRA(weight)
      ,CHFp_CONST_INT(ref_ratio) );

#define FORT_NODEAVERAGE_GETWEIGHTS FORTRAN_NAME( inlineNODEAVERAGE_GETWEIGHTS, inlineNODEAVERAGE_GETWEIGHTS)
#define FORTNT_NODEAVERAGE_GETWEIGHTS FORTRAN_NAME( NODEAVERAGE_GETWEIGHTS, nodeaverage_getweights)

inline void FORTRAN_NAME(inlineNODEAVERAGE_GETWEIGHTS, inlineNODEAVERAGE_GETWEIGHTS)(
      CHFp_FRA(weight)
      ,CHFp_CONST_INT(ref_ratio) )
{
 CH_TIMELEAF("FORT_NODEAVERAGE_GETWEIGHTS");
 FORTRAN_NAME( NODEAVERAGE_GETWEIGHTS ,nodeaverage_getweights )(
      CHFt_FRA(weight)
      ,CHFt_CONST_INT(ref_ratio) );
}
#endif  // GUARDNODEAVERAGE_GETWEIGHTS 

}

#endif
