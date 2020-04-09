#ifndef _EBMESHINTERPF_F_H_
#define _EBMESHINTERPF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDNGP_INTERPOLATE_SCALAR 
#define GUARDNGP_INTERPOLATE_SCALAR 
// Prototype for Fortran procedure ngp_interpolate_scalar ...
//
void FORTRAN_NAME( NGP_INTERPOLATE_SCALAR ,ngp_interpolate_scalar )(
      CHFp_REAL(particle_field)
      ,CHFp_CONST_FRA1(field)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x) );

#define FORT_NGP_INTERPOLATE_SCALAR FORTRAN_NAME( inlineNGP_INTERPOLATE_SCALAR, inlineNGP_INTERPOLATE_SCALAR)
#define FORTNT_NGP_INTERPOLATE_SCALAR FORTRAN_NAME( NGP_INTERPOLATE_SCALAR, ngp_interpolate_scalar)

inline void FORTRAN_NAME(inlineNGP_INTERPOLATE_SCALAR, inlineNGP_INTERPOLATE_SCALAR)(
      CHFp_REAL(particle_field)
      ,CHFp_CONST_FRA1(field)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x) )
{
 CH_TIMELEAF("FORT_NGP_INTERPOLATE_SCALAR");
 FORTRAN_NAME( NGP_INTERPOLATE_SCALAR ,ngp_interpolate_scalar )(
      CHFt_REAL(particle_field)
      ,CHFt_CONST_FRA1(field)
      ,CHFt_CONST_REALVECT(left_edge)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(x) );
}
#endif  // GUARDNGP_INTERPOLATE_SCALAR 

#ifndef GUARDCIC_INTERPOLATE_SCALAR 
#define GUARDCIC_INTERPOLATE_SCALAR 
// Prototype for Fortran procedure cic_interpolate_scalar ...
//
void FORTRAN_NAME( CIC_INTERPOLATE_SCALAR ,cic_interpolate_scalar )(
      CHFp_REAL(particle_field)
      ,CHFp_CONST_FRA1(field)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x) );

#define FORT_CIC_INTERPOLATE_SCALAR FORTRAN_NAME( inlineCIC_INTERPOLATE_SCALAR, inlineCIC_INTERPOLATE_SCALAR)
#define FORTNT_CIC_INTERPOLATE_SCALAR FORTRAN_NAME( CIC_INTERPOLATE_SCALAR, cic_interpolate_scalar)

inline void FORTRAN_NAME(inlineCIC_INTERPOLATE_SCALAR, inlineCIC_INTERPOLATE_SCALAR)(
      CHFp_REAL(particle_field)
      ,CHFp_CONST_FRA1(field)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x) )
{
 CH_TIMELEAF("FORT_CIC_INTERPOLATE_SCALAR");
 FORTRAN_NAME( CIC_INTERPOLATE_SCALAR ,cic_interpolate_scalar )(
      CHFt_REAL(particle_field)
      ,CHFt_CONST_FRA1(field)
      ,CHFt_CONST_REALVECT(left_edge)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(x) );
}
#endif  // GUARDCIC_INTERPOLATE_SCALAR 

#ifndef GUARDTSC_INTERPOLATE_SCALAR 
#define GUARDTSC_INTERPOLATE_SCALAR 
// Prototype for Fortran procedure tsc_interpolate_scalar ...
//
void FORTRAN_NAME( TSC_INTERPOLATE_SCALAR ,tsc_interpolate_scalar )(
      CHFp_REAL(particle_field)
      ,CHFp_CONST_FRA1(field)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x) );

#define FORT_TSC_INTERPOLATE_SCALAR FORTRAN_NAME( inlineTSC_INTERPOLATE_SCALAR, inlineTSC_INTERPOLATE_SCALAR)
#define FORTNT_TSC_INTERPOLATE_SCALAR FORTRAN_NAME( TSC_INTERPOLATE_SCALAR, tsc_interpolate_scalar)

inline void FORTRAN_NAME(inlineTSC_INTERPOLATE_SCALAR, inlineTSC_INTERPOLATE_SCALAR)(
      CHFp_REAL(particle_field)
      ,CHFp_CONST_FRA1(field)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x) )
{
 CH_TIMELEAF("FORT_TSC_INTERPOLATE_SCALAR");
 FORTRAN_NAME( TSC_INTERPOLATE_SCALAR ,tsc_interpolate_scalar )(
      CHFt_REAL(particle_field)
      ,CHFt_CONST_FRA1(field)
      ,CHFt_CONST_REALVECT(left_edge)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(x) );
}
#endif  // GUARDTSC_INTERPOLATE_SCALAR 

#ifndef GUARDW4_INTERPOLATE_SCALAR 
#define GUARDW4_INTERPOLATE_SCALAR 
// Prototype for Fortran procedure w4_interpolate_scalar ...
//
void FORTRAN_NAME( W4_INTERPOLATE_SCALAR ,w4_interpolate_scalar )(
      CHFp_REAL(particle_field)
      ,CHFp_CONST_FRA1(field)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x) );

#define FORT_W4_INTERPOLATE_SCALAR FORTRAN_NAME( inlineW4_INTERPOLATE_SCALAR, inlineW4_INTERPOLATE_SCALAR)
#define FORTNT_W4_INTERPOLATE_SCALAR FORTRAN_NAME( W4_INTERPOLATE_SCALAR, w4_interpolate_scalar)

inline void FORTRAN_NAME(inlineW4_INTERPOLATE_SCALAR, inlineW4_INTERPOLATE_SCALAR)(
      CHFp_REAL(particle_field)
      ,CHFp_CONST_FRA1(field)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x) )
{
 CH_TIMELEAF("FORT_W4_INTERPOLATE_SCALAR");
 FORTRAN_NAME( W4_INTERPOLATE_SCALAR ,w4_interpolate_scalar )(
      CHFt_REAL(particle_field)
      ,CHFt_CONST_FRA1(field)
      ,CHFt_CONST_REALVECT(left_edge)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(x) );
}
#endif  // GUARDW4_INTERPOLATE_SCALAR 

}

#endif
