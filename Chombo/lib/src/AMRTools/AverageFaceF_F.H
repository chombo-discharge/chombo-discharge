#ifndef _AVERAGEFACEF_F_H_
#define _AVERAGEFACEF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDAVERAGEFACE 
#define GUARDAVERAGEFACE 
// Prototype for Fortran procedure AVERAGEFACE ...
//
void FORTRAN_NAME( AVERAGEFACE ,averageface )(
      CHFp_FRA(coarse)
      ,CHFp_CONST_FRA(fine)
      ,CHFp_BOX(crseBox)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(nRef)
      ,CHFp_CONST_INT(refFactor)
      ,CHFp_BOX(refBox) );

#define FORT_AVERAGEFACE FORTRAN_NAME( inlineAVERAGEFACE, inlineAVERAGEFACE)
#define FORTNT_AVERAGEFACE FORTRAN_NAME( AVERAGEFACE, averageface)

inline void FORTRAN_NAME(inlineAVERAGEFACE, inlineAVERAGEFACE)(
      CHFp_FRA(coarse)
      ,CHFp_CONST_FRA(fine)
      ,CHFp_BOX(crseBox)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(nRef)
      ,CHFp_CONST_INT(refFactor)
      ,CHFp_BOX(refBox) )
{
 CH_TIMELEAF("FORT_AVERAGEFACE");
 FORTRAN_NAME( AVERAGEFACE ,averageface )(
      CHFt_FRA(coarse)
      ,CHFt_CONST_FRA(fine)
      ,CHFt_BOX(crseBox)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(nRef)
      ,CHFt_CONST_INT(refFactor)
      ,CHFt_BOX(refBox) );
}
#endif  // GUARDAVERAGEFACE 

#ifndef GUARDAVERAGEFACEHARMONIC 
#define GUARDAVERAGEFACEHARMONIC 
// Prototype for Fortran procedure AVERAGEFACEHARMONIC ...
//
void FORTRAN_NAME( AVERAGEFACEHARMONIC ,averagefaceharmonic )(
      CHFp_FRA(coarse)
      ,CHFp_CONST_FRA(fine)
      ,CHFp_BOX(crseBox)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(nRef)
      ,CHFp_CONST_INT(refFactor)
      ,CHFp_BOX(refBox) );

#define FORT_AVERAGEFACEHARMONIC FORTRAN_NAME( inlineAVERAGEFACEHARMONIC, inlineAVERAGEFACEHARMONIC)
#define FORTNT_AVERAGEFACEHARMONIC FORTRAN_NAME( AVERAGEFACEHARMONIC, averagefaceharmonic)

inline void FORTRAN_NAME(inlineAVERAGEFACEHARMONIC, inlineAVERAGEFACEHARMONIC)(
      CHFp_FRA(coarse)
      ,CHFp_CONST_FRA(fine)
      ,CHFp_BOX(crseBox)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(nRef)
      ,CHFp_CONST_INT(refFactor)
      ,CHFp_BOX(refBox) )
{
 CH_TIMELEAF("FORT_AVERAGEFACEHARMONIC");
 FORTRAN_NAME( AVERAGEFACEHARMONIC ,averagefaceharmonic )(
      CHFt_FRA(coarse)
      ,CHFt_CONST_FRA(fine)
      ,CHFt_BOX(crseBox)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(nRef)
      ,CHFt_CONST_INT(refFactor)
      ,CHFt_BOX(refBox) );
}
#endif  // GUARDAVERAGEFACEHARMONIC 

}

#endif
