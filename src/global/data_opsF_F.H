#ifndef _DATA_OPSF_F_H_
#define _DATA_OPSF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDVECTOR_LENGTH 
#define GUARDVECTOR_LENGTH 
// Prototype for Fortran procedure vector_length ...
//
void FORTRAN_NAME( VECTOR_LENGTH ,vector_length )(
      CHFp_FRA1(dst)
      ,CHFp_CONST_FRA(src)
      ,CHFp_BOX(region) );

#define FORT_VECTOR_LENGTH FORTRAN_NAME( inlineVECTOR_LENGTH, inlineVECTOR_LENGTH)
#define FORTNT_VECTOR_LENGTH FORTRAN_NAME( VECTOR_LENGTH, vector_length)

inline void FORTRAN_NAME(inlineVECTOR_LENGTH, inlineVECTOR_LENGTH)(
      CHFp_FRA1(dst)
      ,CHFp_CONST_FRA(src)
      ,CHFp_BOX(region) )
{
 CH_TIMELEAF("FORT_VECTOR_LENGTH");
 FORTRAN_NAME( VECTOR_LENGTH ,vector_length )(
      CHFt_FRA1(dst)
      ,CHFt_CONST_FRA(src)
      ,CHFt_BOX(region) );
}
#endif  // GUARDVECTOR_LENGTH 

#ifndef GUARDFLOOR 
#define GUARDFLOOR 
// Prototype for Fortran procedure floor ...
//
void FORTRAN_NAME( FLOOR ,floor )(
      CHFp_FRA(src)
      ,CHFp_CONST_INT(ncomps)
      ,CHFp_CONST_REAL(value)
      ,CHFp_BOX(region) );

#define FORT_FLOOR FORTRAN_NAME( inlineFLOOR, inlineFLOOR)
#define FORTNT_FLOOR FORTRAN_NAME( FLOOR, floor)

inline void FORTRAN_NAME(inlineFLOOR, inlineFLOOR)(
      CHFp_FRA(src)
      ,CHFp_CONST_INT(ncomps)
      ,CHFp_CONST_REAL(value)
      ,CHFp_BOX(region) )
{
 CH_TIMELEAF("FORT_FLOOR");
 FORTRAN_NAME( FLOOR ,floor )(
      CHFt_FRA(src)
      ,CHFt_CONST_INT(ncomps)
      ,CHFt_CONST_REAL(value)
      ,CHFt_BOX(region) );
}
#endif  // GUARDFLOOR 

#ifndef GUARDMAX_MIN_NORM 
#define GUARDMAX_MIN_NORM 
// Prototype for Fortran procedure max_min_norm ...
//
void FORTRAN_NAME( MAX_MIN_NORM ,max_min_norm )(
      CHFp_REAL(max_val)
      ,CHFp_REAL(min_val)
      ,CHFp_CONST_FRA(src)
      ,CHFp_CONST_INT(ncomps)
      ,CHFp_BOX(box) );

#define FORT_MAX_MIN_NORM FORTRAN_NAME( inlineMAX_MIN_NORM, inlineMAX_MIN_NORM)
#define FORTNT_MAX_MIN_NORM FORTRAN_NAME( MAX_MIN_NORM, max_min_norm)

inline void FORTRAN_NAME(inlineMAX_MIN_NORM, inlineMAX_MIN_NORM)(
      CHFp_REAL(max_val)
      ,CHFp_REAL(min_val)
      ,CHFp_CONST_FRA(src)
      ,CHFp_CONST_INT(ncomps)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_MAX_MIN_NORM");
 FORTRAN_NAME( MAX_MIN_NORM ,max_min_norm )(
      CHFt_REAL(max_val)
      ,CHFt_REAL(min_val)
      ,CHFt_CONST_FRA(src)
      ,CHFt_CONST_INT(ncomps)
      ,CHFt_BOX(box) );
}
#endif  // GUARDMAX_MIN_NORM 

}

#endif
