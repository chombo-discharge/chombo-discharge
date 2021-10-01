#ifndef _CD_EBGRADIENTF_F_H_
#define _CD_EBGRADIENTF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDGRADIENT 
#define GUARDGRADIENT 
// Prototype for Fortran procedure gradient ...
//
void FORTRAN_NAME( GRADIENT ,gradient )(
      CHFp_FRA(grad)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(region) );

#define FORT_GRADIENT FORTRAN_NAME( inlineGRADIENT, inlineGRADIENT)
#define FORTNT_GRADIENT FORTRAN_NAME( GRADIENT, gradient)

inline void FORTRAN_NAME(inlineGRADIENT, inlineGRADIENT)(
      CHFp_FRA(grad)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(region) )
{
 CH_TIMELEAF("FORT_GRADIENT");
 FORTRAN_NAME( GRADIENT ,gradient )(
      CHFt_FRA(grad)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(region) );
}
#endif  // GUARDGRADIENT 

#ifndef GUARDFACEGRADIENT 
#define GUARDFACEGRADIENT 
// Prototype for Fortran procedure facegradient ...
//
void FORTRAN_NAME( FACEGRADIENT ,facegradient )(
      CHFp_FRA(grad)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(facebox) );

#define FORT_FACEGRADIENT FORTRAN_NAME( inlineFACEGRADIENT, inlineFACEGRADIENT)
#define FORTNT_FACEGRADIENT FORTRAN_NAME( FACEGRADIENT, facegradient)

inline void FORTRAN_NAME(inlineFACEGRADIENT, inlineFACEGRADIENT)(
      CHFp_FRA(grad)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(facebox) )
{
 CH_TIMELEAF("FORT_FACEGRADIENT");
 FORTRAN_NAME( FACEGRADIENT ,facegradient )(
      CHFt_FRA(grad)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(facebox) );
}
#endif  // GUARDFACEGRADIENT 

}

#endif
