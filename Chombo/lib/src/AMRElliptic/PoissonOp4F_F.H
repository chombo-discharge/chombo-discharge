#ifndef _POISSONOP4F_F_H_
#define _POISSONOP4F_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDOPERATORLAP4 
#define GUARDOPERATORLAP4 
// Prototype for Fortran procedure OPERATORLAP4 ...
//
void FORTRAN_NAME( OPERATORLAP4 ,operatorlap4 )(
      CHFp_FRA(lofphi)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta) );

#define FORT_OPERATORLAP4 FORTRAN_NAME( inlineOPERATORLAP4, inlineOPERATORLAP4)
#define FORTNT_OPERATORLAP4 FORTRAN_NAME( OPERATORLAP4, operatorlap4)

inline void FORTRAN_NAME(inlineOPERATORLAP4, inlineOPERATORLAP4)(
      CHFp_FRA(lofphi)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta) )
{
 CH_TIMELEAF("FORT_OPERATORLAP4");
 FORTRAN_NAME( OPERATORLAP4 ,operatorlap4 )(
      CHFt_FRA(lofphi)
      ,CHFt_CONST_FRA(phi)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta) );
}
#endif  // GUARDOPERATORLAP4 

#ifndef GUARDRESTRICTRES4 
#define GUARDRESTRICTRES4 
// Prototype for Fortran procedure RESTRICTRES4 ...
//
void FORTRAN_NAME( RESTRICTRES4 ,restrictres4 )(
      CHFp_FRA(res)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) );

#define FORT_RESTRICTRES4 FORTRAN_NAME( inlineRESTRICTRES4, inlineRESTRICTRES4)
#define FORTNT_RESTRICTRES4 FORTRAN_NAME( RESTRICTRES4, restrictres4)

inline void FORTRAN_NAME(inlineRESTRICTRES4, inlineRESTRICTRES4)(
      CHFp_FRA(res)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_RESTRICTRES4");
 FORTRAN_NAME( RESTRICTRES4 ,restrictres4 )(
      CHFt_FRA(res)
      ,CHFt_CONST_FRA(phi)
      ,CHFt_CONST_FRA(rhs)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDRESTRICTRES4 

#ifndef GUARDGSRBLAPLACIAN4 
#define GUARDGSRBLAPLACIAN4 
// Prototype for Fortran procedure GSRBLAPLACIAN4 ...
//
void FORTRAN_NAME( GSRBLAPLACIAN4 ,gsrblaplacian4 )(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_FRA(tmp)
      ,CHFp_CONST_INT(redBlack) );

#define FORT_GSRBLAPLACIAN4 FORTRAN_NAME( inlineGSRBLAPLACIAN4, inlineGSRBLAPLACIAN4)
#define FORTNT_GSRBLAPLACIAN4 FORTRAN_NAME( GSRBLAPLACIAN4, gsrblaplacian4)

inline void FORTRAN_NAME(inlineGSRBLAPLACIAN4, inlineGSRBLAPLACIAN4)(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_FRA(tmp)
      ,CHFp_CONST_INT(redBlack) )
{
 CH_TIMELEAF("FORT_GSRBLAPLACIAN4");
 FORTRAN_NAME( GSRBLAPLACIAN4 ,gsrblaplacian4 )(
      CHFt_FRA(phi)
      ,CHFt_CONST_FRA(rhs)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_FRA(tmp)
      ,CHFt_CONST_INT(redBlack) );
}
#endif  // GUARDGSRBLAPLACIAN4 

}

#endif
