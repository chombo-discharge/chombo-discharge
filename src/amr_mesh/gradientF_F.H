#ifndef _GRADIENTF_F_H_
#define _GRADIENTF_F_H_

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

}

#endif
