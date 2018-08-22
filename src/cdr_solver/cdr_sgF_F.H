#ifndef _CDR_SGF_F_H_
#define _CDR_SGF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSGFLUX 
#define GUARDSGFLUX 
// Prototype for Fortran procedure sgflux ...
//
void FORTRAN_NAME( SGFLUX ,sgflux )(
      CHFp_FRA1(flux)
      ,CHFp_CONST_FRA1(state)
      ,CHFp_CONST_FRA1(velo)
      ,CHFp_CONST_FRA1(diffco)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(face_box) );

#define FORT_SGFLUX FORTRAN_NAME( inlineSGFLUX, inlineSGFLUX)
#define FORTNT_SGFLUX FORTRAN_NAME( SGFLUX, sgflux)

inline void FORTRAN_NAME(inlineSGFLUX, inlineSGFLUX)(
      CHFp_FRA1(flux)
      ,CHFp_CONST_FRA1(state)
      ,CHFp_CONST_FRA1(velo)
      ,CHFp_CONST_FRA1(diffco)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(face_box) )
{
 CH_TIMELEAF("FORT_SGFLUX");
 FORTRAN_NAME( SGFLUX ,sgflux )(
      CHFt_FRA1(flux)
      ,CHFt_CONST_FRA1(state)
      ,CHFt_CONST_FRA1(velo)
      ,CHFt_CONST_FRA1(diffco)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(face_box) );
}
#endif  // GUARDSGFLUX 

#ifndef GUARDSGUPWIND 
#define GUARDSGUPWIND 
// Prototype for Fortran procedure sgupwind ...
//
void FORTRAN_NAME( SGUPWIND ,sgupwind )(
      CHFp_FRA1(flux)
      ,CHFp_CONST_FRA1(state)
      ,CHFp_CONST_FRA1(velo)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(face_box) );

#define FORT_SGUPWIND FORTRAN_NAME( inlineSGUPWIND, inlineSGUPWIND)
#define FORTNT_SGUPWIND FORTRAN_NAME( SGUPWIND, sgupwind)

inline void FORTRAN_NAME(inlineSGUPWIND, inlineSGUPWIND)(
      CHFp_FRA1(flux)
      ,CHFp_CONST_FRA1(state)
      ,CHFp_CONST_FRA1(velo)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(face_box) )
{
 CH_TIMELEAF("FORT_SGUPWIND");
 FORTRAN_NAME( SGUPWIND ,sgupwind )(
      CHFt_FRA1(flux)
      ,CHFt_CONST_FRA1(state)
      ,CHFt_CONST_FRA1(velo)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(face_box) );
}
#endif  // GUARDSGUPWIND 

}

#endif
