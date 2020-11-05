#ifndef _EBMGAVERAGEF_F_H_
#define _EBMGAVERAGEF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDREGAVERAGE 
#define GUARDREGAVERAGE 
// Prototype for Fortran procedure regaverage ...
//
void FORTRAN_NAME( REGAVERAGE ,regaverage )(
      CHFp_FRA1(coarse)
      ,CHFp_CONST_FRA1(fine)
      ,CHFp_BOX(coarbox)
      ,CHFp_BOX(refbox)
      ,CHFp_CONST_INT(numfinepercoar)
      ,CHFp_CONST_INT(reftocoar) );

#define FORT_REGAVERAGE FORTRAN_NAME( inlineREGAVERAGE, inlineREGAVERAGE)
#define FORTNT_REGAVERAGE FORTRAN_NAME( REGAVERAGE, regaverage)

inline void FORTRAN_NAME(inlineREGAVERAGE, inlineREGAVERAGE)(
      CHFp_FRA1(coarse)
      ,CHFp_CONST_FRA1(fine)
      ,CHFp_BOX(coarbox)
      ,CHFp_BOX(refbox)
      ,CHFp_CONST_INT(numfinepercoar)
      ,CHFp_CONST_INT(reftocoar) )
{
 CH_TIMELEAF("FORT_REGAVERAGE");
 FORTRAN_NAME( REGAVERAGE ,regaverage )(
      CHFt_FRA1(coarse)
      ,CHFt_CONST_FRA1(fine)
      ,CHFt_BOX(coarbox)
      ,CHFt_BOX(refbox)
      ,CHFt_CONST_INT(numfinepercoar)
      ,CHFt_CONST_INT(reftocoar) );
}
#endif  // GUARDREGAVERAGE 

}

#endif
