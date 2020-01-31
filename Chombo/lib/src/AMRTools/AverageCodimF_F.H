#ifndef _AVERAGECODIMF_F_H_
#define _AVERAGECODIMF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDAVERAGECODIM 
#define GUARDAVERAGECODIM 
// Prototype for Fortran procedure averagecodim ...
//
void FORTRAN_NAME( AVERAGECODIM ,averagecodim )(
      CHFp_FRA(coarse)
      ,CHFp_CONST_FRA(fine)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(refRatio)
      ,CHFp_CONST_INT(refFactor)
      ,CHFp_CONST_INT(refDim)
      ,CHFp_BOX(bref) );

#define FORT_AVERAGECODIM FORTRAN_NAME( inlineAVERAGECODIM, inlineAVERAGECODIM)
#define FORTNT_AVERAGECODIM FORTRAN_NAME( AVERAGECODIM, averagecodim)

inline void FORTRAN_NAME(inlineAVERAGECODIM, inlineAVERAGECODIM)(
      CHFp_FRA(coarse)
      ,CHFp_CONST_FRA(fine)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(refRatio)
      ,CHFp_CONST_INT(refFactor)
      ,CHFp_CONST_INT(refDim)
      ,CHFp_BOX(bref) )
{
 CH_TIMELEAF("FORT_AVERAGECODIM");
 FORTRAN_NAME( AVERAGECODIM ,averagecodim )(
      CHFt_FRA(coarse)
      ,CHFt_CONST_FRA(fine)
      ,CHFt_BOX(box)
      ,CHFt_CONST_INT(refRatio)
      ,CHFt_CONST_INT(refFactor)
      ,CHFt_CONST_INT(refDim)
      ,CHFt_BOX(bref) );
}
#endif  // GUARDAVERAGECODIM 

#ifndef GUARDAVERAGECODIMHARMONIC 
#define GUARDAVERAGECODIMHARMONIC 
// Prototype for Fortran procedure averagecodimharmonic ...
//
void FORTRAN_NAME( AVERAGECODIMHARMONIC ,averagecodimharmonic )(
      CHFp_FRA(coarse)
      ,CHFp_CONST_FRA(fine)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(refRatio)
      ,CHFp_CONST_INT(refFactor)
      ,CHFp_CONST_INT(refDim)
      ,CHFp_BOX(bref) );

#define FORT_AVERAGECODIMHARMONIC FORTRAN_NAME( inlineAVERAGECODIMHARMONIC, inlineAVERAGECODIMHARMONIC)
#define FORTNT_AVERAGECODIMHARMONIC FORTRAN_NAME( AVERAGECODIMHARMONIC, averagecodimharmonic)

inline void FORTRAN_NAME(inlineAVERAGECODIMHARMONIC, inlineAVERAGECODIMHARMONIC)(
      CHFp_FRA(coarse)
      ,CHFp_CONST_FRA(fine)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(refRatio)
      ,CHFp_CONST_INT(refFactor)
      ,CHFp_CONST_INT(refDim)
      ,CHFp_BOX(bref) )
{
 CH_TIMELEAF("FORT_AVERAGECODIMHARMONIC");
 FORTRAN_NAME( AVERAGECODIMHARMONIC ,averagecodimharmonic )(
      CHFt_FRA(coarse)
      ,CHFt_CONST_FRA(fine)
      ,CHFt_BOX(box)
      ,CHFt_CONST_INT(refRatio)
      ,CHFt_CONST_INT(refFactor)
      ,CHFt_CONST_INT(refDim)
      ,CHFt_BOX(bref) );
}
#endif  // GUARDAVERAGECODIMHARMONIC 

}

#endif
