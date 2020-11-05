#ifndef _EDGETOCELLF_F_H_
#define _EDGETOCELLF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDEDGETOCELL 
#define GUARDEDGETOCELL 
// Prototype for Fortran procedure EDGETOCELL ...
//
void FORTRAN_NAME( EDGETOCELL ,edgetocell )(
      CHFp_CONST_FRA1(edgeData)
      ,CHFp_FRA1(cellData)
      ,CHFp_BOX(cellBox)
      ,CHFp_CONST_INT(dir) );

#define FORT_EDGETOCELL FORTRAN_NAME( inlineEDGETOCELL, inlineEDGETOCELL)
#define FORTNT_EDGETOCELL FORTRAN_NAME( EDGETOCELL, edgetocell)

inline void FORTRAN_NAME(inlineEDGETOCELL, inlineEDGETOCELL)(
      CHFp_CONST_FRA1(edgeData)
      ,CHFp_FRA1(cellData)
      ,CHFp_BOX(cellBox)
      ,CHFp_CONST_INT(dir) )
{
 CH_TIMELEAF("FORT_EDGETOCELL");
 FORTRAN_NAME( EDGETOCELL ,edgetocell )(
      CHFt_CONST_FRA1(edgeData)
      ,CHFt_FRA1(cellData)
      ,CHFt_BOX(cellBox)
      ,CHFt_CONST_INT(dir) );
}
#endif  // GUARDEDGETOCELL 

#ifndef GUARDEDGETOINCREMENTCELL 
#define GUARDEDGETOINCREMENTCELL 
// Prototype for Fortran procedure EDGETOINCREMENTCELL ...
//
void FORTRAN_NAME( EDGETOINCREMENTCELL ,edgetoincrementcell )(
      CHFp_CONST_FRA1(edgeData)
      ,CHFp_FRA1(cellData)
      ,CHFp_BOX(cellBox)
      ,CHFp_CONST_INT(dir) );

#define FORT_EDGETOINCREMENTCELL FORTRAN_NAME( inlineEDGETOINCREMENTCELL, inlineEDGETOINCREMENTCELL)
#define FORTNT_EDGETOINCREMENTCELL FORTRAN_NAME( EDGETOINCREMENTCELL, edgetoincrementcell)

inline void FORTRAN_NAME(inlineEDGETOINCREMENTCELL, inlineEDGETOINCREMENTCELL)(
      CHFp_CONST_FRA1(edgeData)
      ,CHFp_FRA1(cellData)
      ,CHFp_BOX(cellBox)
      ,CHFp_CONST_INT(dir) )
{
 CH_TIMELEAF("FORT_EDGETOINCREMENTCELL");
 FORTRAN_NAME( EDGETOINCREMENTCELL ,edgetoincrementcell )(
      CHFt_CONST_FRA1(edgeData)
      ,CHFt_FRA1(cellData)
      ,CHFt_BOX(cellBox)
      ,CHFt_CONST_INT(dir) );
}
#endif  // GUARDEDGETOINCREMENTCELL 

#ifndef GUARDEDGETOCELLMAX 
#define GUARDEDGETOCELLMAX 
// Prototype for Fortran procedure EDGETOCELLMAX ...
//
void FORTRAN_NAME( EDGETOCELLMAX ,edgetocellmax )(
      CHFp_CONST_FRA1(edgeData)
      ,CHFp_FRA1(cellData)
      ,CHFp_BOX(cellBox)
      ,CHFp_CONST_INT(dir) );

#define FORT_EDGETOCELLMAX FORTRAN_NAME( inlineEDGETOCELLMAX, inlineEDGETOCELLMAX)
#define FORTNT_EDGETOCELLMAX FORTRAN_NAME( EDGETOCELLMAX, edgetocellmax)

inline void FORTRAN_NAME(inlineEDGETOCELLMAX, inlineEDGETOCELLMAX)(
      CHFp_CONST_FRA1(edgeData)
      ,CHFp_FRA1(cellData)
      ,CHFp_BOX(cellBox)
      ,CHFp_CONST_INT(dir) )
{
 CH_TIMELEAF("FORT_EDGETOCELLMAX");
 FORTRAN_NAME( EDGETOCELLMAX ,edgetocellmax )(
      CHFt_CONST_FRA1(edgeData)
      ,CHFt_FRA1(cellData)
      ,CHFt_BOX(cellBox)
      ,CHFt_CONST_INT(dir) );
}
#endif  // GUARDEDGETOCELLMAX 

}

#endif
