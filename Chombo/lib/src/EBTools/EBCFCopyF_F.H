#ifndef _EBCFCOPYF_F_H_
#define _EBCFCOPYF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDCOPYCFFAB 
#define GUARDCOPYCFFAB 
// Prototype for Fortran procedure copycffab ...
//
void FORTRAN_NAME( COPYCFFAB ,copycffab )(
      CHFp_FRA1(finelevel)
      ,CHFp_CONST_FRA1(coarlevel)
      ,CHFp_BOX(coarbox)
      ,CHFp_BOX(refbox)
      ,CHFp_CONST_INT(reftocoar) );

#define FORT_COPYCFFAB FORTRAN_NAME( inlineCOPYCFFAB, inlineCOPYCFFAB)
#define FORTNT_COPYCFFAB FORTRAN_NAME( COPYCFFAB, copycffab)

inline void FORTRAN_NAME(inlineCOPYCFFAB, inlineCOPYCFFAB)(
      CHFp_FRA1(finelevel)
      ,CHFp_CONST_FRA1(coarlevel)
      ,CHFp_BOX(coarbox)
      ,CHFp_BOX(refbox)
      ,CHFp_CONST_INT(reftocoar) )
{
 CH_TIMELEAF("FORT_COPYCFFAB");
 FORTRAN_NAME( COPYCFFAB ,copycffab )(
      CHFt_FRA1(finelevel)
      ,CHFt_CONST_FRA1(coarlevel)
      ,CHFt_BOX(coarbox)
      ,CHFt_BOX(refbox)
      ,CHFt_CONST_INT(reftocoar) );
}
#endif  // GUARDCOPYCFFAB 

}

#endif
