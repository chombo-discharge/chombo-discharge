#ifndef _EBPARTICLEINTERPF_F_H_
#define _EBPARTICLEINTERPF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDNGP_DEPOSIT_SCALAR 
#define GUARDNGP_DEPOSIT_SCALAR 
// Prototype for Fortran procedure ngp_deposit_scalar ...
//
void FORTRAN_NAME( NGP_DEPOSIT_SCALAR ,ngp_deposit_scalar )(
      CHFp_FRA1(rho)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(q) );

#define FORT_NGP_DEPOSIT_SCALAR FORTRAN_NAME( inlineNGP_DEPOSIT_SCALAR, inlineNGP_DEPOSIT_SCALAR)
#define FORTNT_NGP_DEPOSIT_SCALAR FORTRAN_NAME( NGP_DEPOSIT_SCALAR, ngp_deposit_scalar)

inline void FORTRAN_NAME(inlineNGP_DEPOSIT_SCALAR, inlineNGP_DEPOSIT_SCALAR)(
      CHFp_FRA1(rho)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(q) )
{
 CH_TIMELEAF("FORT_NGP_DEPOSIT_SCALAR");
 FORTRAN_NAME( NGP_DEPOSIT_SCALAR ,ngp_deposit_scalar )(
      CHFt_FRA1(rho)
      ,CHFt_CONST_REALVECT(left_edge)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(x)
      ,CHFt_CONST_REAL(q) );
}
#endif  // GUARDNGP_DEPOSIT_SCALAR 

#ifndef GUARDCIC_DEPOSIT_SCALAR 
#define GUARDCIC_DEPOSIT_SCALAR 
// Prototype for Fortran procedure cic_deposit_scalar ...
//
void FORTRAN_NAME( CIC_DEPOSIT_SCALAR ,cic_deposit_scalar )(
      CHFp_FRA1(rho)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(q) );

#define FORT_CIC_DEPOSIT_SCALAR FORTRAN_NAME( inlineCIC_DEPOSIT_SCALAR, inlineCIC_DEPOSIT_SCALAR)
#define FORTNT_CIC_DEPOSIT_SCALAR FORTRAN_NAME( CIC_DEPOSIT_SCALAR, cic_deposit_scalar)

inline void FORTRAN_NAME(inlineCIC_DEPOSIT_SCALAR, inlineCIC_DEPOSIT_SCALAR)(
      CHFp_FRA1(rho)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(q) )
{
 CH_TIMELEAF("FORT_CIC_DEPOSIT_SCALAR");
 FORTRAN_NAME( CIC_DEPOSIT_SCALAR ,cic_deposit_scalar )(
      CHFt_FRA1(rho)
      ,CHFt_CONST_REALVECT(left_edge)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(x)
      ,CHFt_CONST_REAL(q) );
}
#endif  // GUARDCIC_DEPOSIT_SCALAR 

#ifndef GUARDTSC_DEPOSIT_SCALAR 
#define GUARDTSC_DEPOSIT_SCALAR 
// Prototype for Fortran procedure tsc_deposit_scalar ...
//
void FORTRAN_NAME( TSC_DEPOSIT_SCALAR ,tsc_deposit_scalar )(
      CHFp_FRA1(rho)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(q) );

#define FORT_TSC_DEPOSIT_SCALAR FORTRAN_NAME( inlineTSC_DEPOSIT_SCALAR, inlineTSC_DEPOSIT_SCALAR)
#define FORTNT_TSC_DEPOSIT_SCALAR FORTRAN_NAME( TSC_DEPOSIT_SCALAR, tsc_deposit_scalar)

inline void FORTRAN_NAME(inlineTSC_DEPOSIT_SCALAR, inlineTSC_DEPOSIT_SCALAR)(
      CHFp_FRA1(rho)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(q) )
{
 CH_TIMELEAF("FORT_TSC_DEPOSIT_SCALAR");
 FORTRAN_NAME( TSC_DEPOSIT_SCALAR ,tsc_deposit_scalar )(
      CHFt_FRA1(rho)
      ,CHFt_CONST_REALVECT(left_edge)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(x)
      ,CHFt_CONST_REAL(q) );
}
#endif  // GUARDTSC_DEPOSIT_SCALAR 

#ifndef GUARDW4_DEPOSIT_SCALAR 
#define GUARDW4_DEPOSIT_SCALAR 
// Prototype for Fortran procedure w4_deposit_scalar ...
//
void FORTRAN_NAME( W4_DEPOSIT_SCALAR ,w4_deposit_scalar )(
      CHFp_FRA1(rho)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(q) );

#define FORT_W4_DEPOSIT_SCALAR FORTRAN_NAME( inlineW4_DEPOSIT_SCALAR, inlineW4_DEPOSIT_SCALAR)
#define FORTNT_W4_DEPOSIT_SCALAR FORTRAN_NAME( W4_DEPOSIT_SCALAR, w4_deposit_scalar)

inline void FORTRAN_NAME(inlineW4_DEPOSIT_SCALAR, inlineW4_DEPOSIT_SCALAR)(
      CHFp_FRA1(rho)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(q) )
{
 CH_TIMELEAF("FORT_W4_DEPOSIT_SCALAR");
 FORTRAN_NAME( W4_DEPOSIT_SCALAR ,w4_deposit_scalar )(
      CHFt_FRA1(rho)
      ,CHFt_CONST_REALVECT(left_edge)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(x)
      ,CHFt_CONST_REAL(q) );
}
#endif  // GUARDW4_DEPOSIT_SCALAR 

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

#ifndef GUARDNGP_INTERPOLATE_VECTOR 
#define GUARDNGP_INTERPOLATE_VECTOR 
// Prototype for Fortran procedure ngp_interpolate_vector ...
//
void FORTRAN_NAME( NGP_INTERPOLATE_VECTOR ,ngp_interpolate_vector )(
      CHFp_REALVECT(particle_field)
      ,CHFp_CONST_FRA(field)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x) );

#define FORT_NGP_INTERPOLATE_VECTOR FORTRAN_NAME( inlineNGP_INTERPOLATE_VECTOR, inlineNGP_INTERPOLATE_VECTOR)
#define FORTNT_NGP_INTERPOLATE_VECTOR FORTRAN_NAME( NGP_INTERPOLATE_VECTOR, ngp_interpolate_vector)

inline void FORTRAN_NAME(inlineNGP_INTERPOLATE_VECTOR, inlineNGP_INTERPOLATE_VECTOR)(
      CHFp_REALVECT(particle_field)
      ,CHFp_CONST_FRA(field)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x) )
{
 CH_TIMELEAF("FORT_NGP_INTERPOLATE_VECTOR");
 FORTRAN_NAME( NGP_INTERPOLATE_VECTOR ,ngp_interpolate_vector )(
      CHFt_REALVECT(particle_field)
      ,CHFt_CONST_FRA(field)
      ,CHFt_CONST_REALVECT(left_edge)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(x) );
}
#endif  // GUARDNGP_INTERPOLATE_VECTOR 

#ifndef GUARDCIC_INTERPOLATE_VECTOR 
#define GUARDCIC_INTERPOLATE_VECTOR 
// Prototype for Fortran procedure cic_interpolate_vector ...
//
void FORTRAN_NAME( CIC_INTERPOLATE_VECTOR ,cic_interpolate_vector )(
      CHFp_REALVECT(particle_field)
      ,CHFp_CONST_FRA(field)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x) );

#define FORT_CIC_INTERPOLATE_VECTOR FORTRAN_NAME( inlineCIC_INTERPOLATE_VECTOR, inlineCIC_INTERPOLATE_VECTOR)
#define FORTNT_CIC_INTERPOLATE_VECTOR FORTRAN_NAME( CIC_INTERPOLATE_VECTOR, cic_interpolate_vector)

inline void FORTRAN_NAME(inlineCIC_INTERPOLATE_VECTOR, inlineCIC_INTERPOLATE_VECTOR)(
      CHFp_REALVECT(particle_field)
      ,CHFp_CONST_FRA(field)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x) )
{
 CH_TIMELEAF("FORT_CIC_INTERPOLATE_VECTOR");
 FORTRAN_NAME( CIC_INTERPOLATE_VECTOR ,cic_interpolate_vector )(
      CHFt_REALVECT(particle_field)
      ,CHFt_CONST_FRA(field)
      ,CHFt_CONST_REALVECT(left_edge)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(x) );
}
#endif  // GUARDCIC_INTERPOLATE_VECTOR 

#ifndef GUARDTSC_INTERPOLATE_VECTOR 
#define GUARDTSC_INTERPOLATE_VECTOR 
// Prototype for Fortran procedure tsc_interpolate_vector ...
//
void FORTRAN_NAME( TSC_INTERPOLATE_VECTOR ,tsc_interpolate_vector )(
      CHFp_REALVECT(particle_field)
      ,CHFp_CONST_FRA(field)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x) );

#define FORT_TSC_INTERPOLATE_VECTOR FORTRAN_NAME( inlineTSC_INTERPOLATE_VECTOR, inlineTSC_INTERPOLATE_VECTOR)
#define FORTNT_TSC_INTERPOLATE_VECTOR FORTRAN_NAME( TSC_INTERPOLATE_VECTOR, tsc_interpolate_vector)

inline void FORTRAN_NAME(inlineTSC_INTERPOLATE_VECTOR, inlineTSC_INTERPOLATE_VECTOR)(
      CHFp_REALVECT(particle_field)
      ,CHFp_CONST_FRA(field)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x) )
{
 CH_TIMELEAF("FORT_TSC_INTERPOLATE_VECTOR");
 FORTRAN_NAME( TSC_INTERPOLATE_VECTOR ,tsc_interpolate_vector )(
      CHFt_REALVECT(particle_field)
      ,CHFt_CONST_FRA(field)
      ,CHFt_CONST_REALVECT(left_edge)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(x) );
}
#endif  // GUARDTSC_INTERPOLATE_VECTOR 

#ifndef GUARDW4_INTERPOLATE_VECTOR 
#define GUARDW4_INTERPOLATE_VECTOR 
// Prototype for Fortran procedure w4_interpolate_vector ...
//
void FORTRAN_NAME( W4_INTERPOLATE_VECTOR ,w4_interpolate_vector )(
      CHFp_REALVECT(particle_field)
      ,CHFp_CONST_FRA(field)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x) );

#define FORT_W4_INTERPOLATE_VECTOR FORTRAN_NAME( inlineW4_INTERPOLATE_VECTOR, inlineW4_INTERPOLATE_VECTOR)
#define FORTNT_W4_INTERPOLATE_VECTOR FORTRAN_NAME( W4_INTERPOLATE_VECTOR, w4_interpolate_vector)

inline void FORTRAN_NAME(inlineW4_INTERPOLATE_VECTOR, inlineW4_INTERPOLATE_VECTOR)(
      CHFp_REALVECT(particle_field)
      ,CHFp_CONST_FRA(field)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x) )
{
 CH_TIMELEAF("FORT_W4_INTERPOLATE_VECTOR");
 FORTRAN_NAME( W4_INTERPOLATE_VECTOR ,w4_interpolate_vector )(
      CHFt_REALVECT(particle_field)
      ,CHFt_CONST_FRA(field)
      ,CHFt_CONST_REALVECT(left_edge)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(x) );
}
#endif  // GUARDW4_INTERPOLATE_VECTOR 

}

#endif
