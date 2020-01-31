#ifndef _LEVELFLUXREGISTERF_F_H_
#define _LEVELFLUXREGISTERF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDINCREMENTFINE 
#define GUARDINCREMENTFINE 
// Prototype for Fortran procedure INCREMENTFINE ...
//
void FORTRAN_NAME( INCREMENTFINE ,incrementfine )(
      CHFp_CONST_FRA(fine)
      ,CHFp_FRA(cFine)
      ,CHFp_BOX(fineBox)
      ,CHFp_CONST_INTVECT(nRef)
      ,CHFp_CONST_REAL(scale)
      ,CHFp_CONST_INT(srcStart)
      ,CHFp_CONST_INT(destStart)
      ,CHFp_CONST_INT(ncomp) );

#define FORT_INCREMENTFINE FORTRAN_NAME( inlineINCREMENTFINE, inlineINCREMENTFINE)
#define FORTNT_INCREMENTFINE FORTRAN_NAME( INCREMENTFINE, incrementfine)

inline void FORTRAN_NAME(inlineINCREMENTFINE, inlineINCREMENTFINE)(
      CHFp_CONST_FRA(fine)
      ,CHFp_FRA(cFine)
      ,CHFp_BOX(fineBox)
      ,CHFp_CONST_INTVECT(nRef)
      ,CHFp_CONST_REAL(scale)
      ,CHFp_CONST_INT(srcStart)
      ,CHFp_CONST_INT(destStart)
      ,CHFp_CONST_INT(ncomp) )
{
 CH_TIMELEAF("FORT_INCREMENTFINE");
 FORTRAN_NAME( INCREMENTFINE ,incrementfine )(
      CHFt_CONST_FRA(fine)
      ,CHFt_FRA(cFine)
      ,CHFt_BOX(fineBox)
      ,CHFt_CONST_INTVECT(nRef)
      ,CHFt_CONST_REAL(scale)
      ,CHFt_CONST_INT(srcStart)
      ,CHFt_CONST_INT(destStart)
      ,CHFt_CONST_INT(ncomp) );
}
#endif  // GUARDINCREMENTFINE 

}

#endif
