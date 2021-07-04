#ifndef _CD_EBMULTIGRIDRESTRICTF_F_H_
#define _CD_EBMULTIGRIDRESTRICTF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDRESTRICTREGULAR 
#define GUARDRESTRICTREGULAR 
// Prototype for Fortran procedure restrictregular ...
//
void FORTRAN_NAME( RESTRICTREGULAR ,restrictregular )(
      CHFp_FRA1(coarse)
      ,CHFp_CONST_FRA1(fine)
      ,CHFp_BOX(coarbox)
      ,CHFp_BOX(refbox)
      ,CHFp_CONST_INT(numfinepercoar)
      ,CHFp_CONST_INT(reftocoar) );

#define FORT_RESTRICTREGULAR FORTRAN_NAME( inlineRESTRICTREGULAR, inlineRESTRICTREGULAR)
#define FORTNT_RESTRICTREGULAR FORTRAN_NAME( RESTRICTREGULAR, restrictregular)

inline void FORTRAN_NAME(inlineRESTRICTREGULAR, inlineRESTRICTREGULAR)(
      CHFp_FRA1(coarse)
      ,CHFp_CONST_FRA1(fine)
      ,CHFp_BOX(coarbox)
      ,CHFp_BOX(refbox)
      ,CHFp_CONST_INT(numfinepercoar)
      ,CHFp_CONST_INT(reftocoar) )
{
 CH_TIMELEAF("FORT_RESTRICTREGULAR");
 FORTRAN_NAME( RESTRICTREGULAR ,restrictregular )(
      CHFt_FRA1(coarse)
      ,CHFt_CONST_FRA1(fine)
      ,CHFt_BOX(coarbox)
      ,CHFt_BOX(refbox)
      ,CHFt_CONST_INT(numfinepercoar)
      ,CHFt_CONST_INT(reftocoar) );
}
#endif  // GUARDRESTRICTREGULAR 

}

#endif
