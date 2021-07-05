#ifndef _CD_EBMULTIGRIDDANCEF_F_H_
#define _CD_EBMULTIGRIDDANCEF_F_H_

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

#ifndef GUARDPROLONGREGULAR 
#define GUARDPROLONGREGULAR 
// Prototype for Fortran procedure prolongregular ...
//
void FORTRAN_NAME( PROLONGREGULAR ,prolongregular )(
      CHFp_FRA1(fineData)
      ,CHFp_CONST_FRA1(coarData)
      ,CHFp_BOX(coarBox)
      ,CHFp_BOX(refBox)
      ,CHFp_CONST_INT(refToCoar) );

#define FORT_PROLONGREGULAR FORTRAN_NAME( inlinePROLONGREGULAR, inlinePROLONGREGULAR)
#define FORTNT_PROLONGREGULAR FORTRAN_NAME( PROLONGREGULAR, prolongregular)

inline void FORTRAN_NAME(inlinePROLONGREGULAR, inlinePROLONGREGULAR)(
      CHFp_FRA1(fineData)
      ,CHFp_CONST_FRA1(coarData)
      ,CHFp_BOX(coarBox)
      ,CHFp_BOX(refBox)
      ,CHFp_CONST_INT(refToCoar) )
{
 CH_TIMELEAF("FORT_PROLONGREGULAR");
 FORTRAN_NAME( PROLONGREGULAR ,prolongregular )(
      CHFt_FRA1(fineData)
      ,CHFt_CONST_FRA1(coarData)
      ,CHFt_BOX(coarBox)
      ,CHFt_BOX(refBox)
      ,CHFt_CONST_INT(refToCoar) );
}
#endif  // GUARDPROLONGREGULAR 

}

#endif
