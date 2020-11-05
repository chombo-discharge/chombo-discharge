#ifndef _MESHINTERPF_F_H_
#define _MESHINTERPF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDNGP_DEPOSIT 
#define GUARDNGP_DEPOSIT 
// Prototype for Fortran procedure ngp_deposit ...
//
void FORTRAN_NAME( NGP_DEPOSIT ,ngp_deposit )(
      CHFp_FRA1(rho)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(q) );

#define FORT_NGP_DEPOSIT FORTRAN_NAME( inlineNGP_DEPOSIT, inlineNGP_DEPOSIT)
#define FORTNT_NGP_DEPOSIT FORTRAN_NAME( NGP_DEPOSIT, ngp_deposit)

inline void FORTRAN_NAME(inlineNGP_DEPOSIT, inlineNGP_DEPOSIT)(
      CHFp_FRA1(rho)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(q) )
{
 CH_TIMELEAF("FORT_NGP_DEPOSIT");
 FORTRAN_NAME( NGP_DEPOSIT ,ngp_deposit )(
      CHFt_FRA1(rho)
      ,CHFt_CONST_REALVECT(left_edge)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(x)
      ,CHFt_CONST_REAL(q) );
}
#endif  // GUARDNGP_DEPOSIT 

#ifndef GUARDNGP_INTERPOLATE 
#define GUARDNGP_INTERPOLATE 
// Prototype for Fortran procedure ngp_interpolate ...
//
void FORTRAN_NAME( NGP_INTERPOLATE ,ngp_interpolate )(
      CHFp_REALVECT(particle_field)
      ,CHFp_CONST_FRA(field)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x) );

#define FORT_NGP_INTERPOLATE FORTRAN_NAME( inlineNGP_INTERPOLATE, inlineNGP_INTERPOLATE)
#define FORTNT_NGP_INTERPOLATE FORTRAN_NAME( NGP_INTERPOLATE, ngp_interpolate)

inline void FORTRAN_NAME(inlineNGP_INTERPOLATE, inlineNGP_INTERPOLATE)(
      CHFp_REALVECT(particle_field)
      ,CHFp_CONST_FRA(field)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x) )
{
 CH_TIMELEAF("FORT_NGP_INTERPOLATE");
 FORTRAN_NAME( NGP_INTERPOLATE ,ngp_interpolate )(
      CHFt_REALVECT(particle_field)
      ,CHFt_CONST_FRA(field)
      ,CHFt_CONST_REALVECT(left_edge)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(x) );
}
#endif  // GUARDNGP_INTERPOLATE 

#ifndef GUARDCIC_DEPOSIT 
#define GUARDCIC_DEPOSIT 
// Prototype for Fortran procedure cic_deposit ...
//
void FORTRAN_NAME( CIC_DEPOSIT ,cic_deposit )(
      CHFp_FRA1(rho)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(q) );

#define FORT_CIC_DEPOSIT FORTRAN_NAME( inlineCIC_DEPOSIT, inlineCIC_DEPOSIT)
#define FORTNT_CIC_DEPOSIT FORTRAN_NAME( CIC_DEPOSIT, cic_deposit)

inline void FORTRAN_NAME(inlineCIC_DEPOSIT, inlineCIC_DEPOSIT)(
      CHFp_FRA1(rho)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(q) )
{
 CH_TIMELEAF("FORT_CIC_DEPOSIT");
 FORTRAN_NAME( CIC_DEPOSIT ,cic_deposit )(
      CHFt_FRA1(rho)
      ,CHFt_CONST_REALVECT(left_edge)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(x)
      ,CHFt_CONST_REAL(q) );
}
#endif  // GUARDCIC_DEPOSIT 

#ifndef GUARDCIC_INTERPOLATE 
#define GUARDCIC_INTERPOLATE 
// Prototype for Fortran procedure cic_interpolate ...
//
void FORTRAN_NAME( CIC_INTERPOLATE ,cic_interpolate )(
      CHFp_REALVECT(particle_field)
      ,CHFp_CONST_FRA(field)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x) );

#define FORT_CIC_INTERPOLATE FORTRAN_NAME( inlineCIC_INTERPOLATE, inlineCIC_INTERPOLATE)
#define FORTNT_CIC_INTERPOLATE FORTRAN_NAME( CIC_INTERPOLATE, cic_interpolate)

inline void FORTRAN_NAME(inlineCIC_INTERPOLATE, inlineCIC_INTERPOLATE)(
      CHFp_REALVECT(particle_field)
      ,CHFp_CONST_FRA(field)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x) )
{
 CH_TIMELEAF("FORT_CIC_INTERPOLATE");
 FORTRAN_NAME( CIC_INTERPOLATE ,cic_interpolate )(
      CHFt_REALVECT(particle_field)
      ,CHFt_CONST_FRA(field)
      ,CHFt_CONST_REALVECT(left_edge)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(x) );
}
#endif  // GUARDCIC_INTERPOLATE 

#ifndef GUARDTSC_DEPOSIT 
#define GUARDTSC_DEPOSIT 
// Prototype for Fortran procedure tsc_deposit ...
//
void FORTRAN_NAME( TSC_DEPOSIT ,tsc_deposit )(
      CHFp_FRA1(rho)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(q) );

#define FORT_TSC_DEPOSIT FORTRAN_NAME( inlineTSC_DEPOSIT, inlineTSC_DEPOSIT)
#define FORTNT_TSC_DEPOSIT FORTRAN_NAME( TSC_DEPOSIT, tsc_deposit)

inline void FORTRAN_NAME(inlineTSC_DEPOSIT, inlineTSC_DEPOSIT)(
      CHFp_FRA1(rho)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(q) )
{
 CH_TIMELEAF("FORT_TSC_DEPOSIT");
 FORTRAN_NAME( TSC_DEPOSIT ,tsc_deposit )(
      CHFt_FRA1(rho)
      ,CHFt_CONST_REALVECT(left_edge)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(x)
      ,CHFt_CONST_REAL(q) );
}
#endif  // GUARDTSC_DEPOSIT 

#ifndef GUARDTSC_INTERPOLATE 
#define GUARDTSC_INTERPOLATE 
// Prototype for Fortran procedure tsc_interpolate ...
//
void FORTRAN_NAME( TSC_INTERPOLATE ,tsc_interpolate )(
      CHFp_REALVECT(particle_field)
      ,CHFp_CONST_FRA(field)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x) );

#define FORT_TSC_INTERPOLATE FORTRAN_NAME( inlineTSC_INTERPOLATE, inlineTSC_INTERPOLATE)
#define FORTNT_TSC_INTERPOLATE FORTRAN_NAME( TSC_INTERPOLATE, tsc_interpolate)

inline void FORTRAN_NAME(inlineTSC_INTERPOLATE, inlineTSC_INTERPOLATE)(
      CHFp_REALVECT(particle_field)
      ,CHFp_CONST_FRA(field)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x) )
{
 CH_TIMELEAF("FORT_TSC_INTERPOLATE");
 FORTRAN_NAME( TSC_INTERPOLATE ,tsc_interpolate )(
      CHFt_REALVECT(particle_field)
      ,CHFt_CONST_FRA(field)
      ,CHFt_CONST_REALVECT(left_edge)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(x) );
}
#endif  // GUARDTSC_INTERPOLATE 

#ifndef GUARDW4_DEPOSIT 
#define GUARDW4_DEPOSIT 
// Prototype for Fortran procedure w4_deposit ...
//
void FORTRAN_NAME( W4_DEPOSIT ,w4_deposit )(
      CHFp_FRA1(rho)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(q) );

#define FORT_W4_DEPOSIT FORTRAN_NAME( inlineW4_DEPOSIT, inlineW4_DEPOSIT)
#define FORTNT_W4_DEPOSIT FORTRAN_NAME( W4_DEPOSIT, w4_deposit)

inline void FORTRAN_NAME(inlineW4_DEPOSIT, inlineW4_DEPOSIT)(
      CHFp_FRA1(rho)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(q) )
{
 CH_TIMELEAF("FORT_W4_DEPOSIT");
 FORTRAN_NAME( W4_DEPOSIT ,w4_deposit )(
      CHFt_FRA1(rho)
      ,CHFt_CONST_REALVECT(left_edge)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(x)
      ,CHFt_CONST_REAL(q) );
}
#endif  // GUARDW4_DEPOSIT 

#ifndef GUARDW4_INTERPOLATE 
#define GUARDW4_INTERPOLATE 
// Prototype for Fortran procedure w4_interpolate ...
//
void FORTRAN_NAME( W4_INTERPOLATE ,w4_interpolate )(
      CHFp_REALVECT(particle_field)
      ,CHFp_CONST_FRA(field)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x) );

#define FORT_W4_INTERPOLATE FORTRAN_NAME( inlineW4_INTERPOLATE, inlineW4_INTERPOLATE)
#define FORTNT_W4_INTERPOLATE FORTRAN_NAME( W4_INTERPOLATE, w4_interpolate)

inline void FORTRAN_NAME(inlineW4_INTERPOLATE, inlineW4_INTERPOLATE)(
      CHFp_REALVECT(particle_field)
      ,CHFp_CONST_FRA(field)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x) )
{
 CH_TIMELEAF("FORT_W4_INTERPOLATE");
 FORTRAN_NAME( W4_INTERPOLATE ,w4_interpolate )(
      CHFt_REALVECT(particle_field)
      ,CHFt_CONST_FRA(field)
      ,CHFt_CONST_REALVECT(left_edge)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(x) );
}
#endif  // GUARDW4_INTERPOLATE 

}

#endif
