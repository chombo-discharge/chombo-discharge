#ifndef _DOTPRODF_F_H_
#define _DOTPRODF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDDOTPRODUCT 
#define GUARDDOTPRODUCT 
// Prototype for Fortran procedure DOTPRODUCT ...
//
void FORTRAN_NAME( DOTPRODUCT ,dotproduct )(
      CHFp_REAL(dotprodout)
      ,CHFp_CONST_FRA(afab)
      ,CHFp_CONST_FRA(bfab)
      ,CHFp_BOX(region)
      ,CHFp_CONST_INT(startcomp)
      ,CHFp_CONST_INT(endcomp) );

#define FORT_DOTPRODUCT FORTRAN_NAME( inlineDOTPRODUCT, inlineDOTPRODUCT)
#define FORTNT_DOTPRODUCT FORTRAN_NAME( DOTPRODUCT, dotproduct)

inline void FORTRAN_NAME(inlineDOTPRODUCT, inlineDOTPRODUCT)(
      CHFp_REAL(dotprodout)
      ,CHFp_CONST_FRA(afab)
      ,CHFp_CONST_FRA(bfab)
      ,CHFp_BOX(region)
      ,CHFp_CONST_INT(startcomp)
      ,CHFp_CONST_INT(endcomp) )
{
 CH_TIMELEAF("FORT_DOTPRODUCT");
 FORTRAN_NAME( DOTPRODUCT ,dotproduct )(
      CHFt_REAL(dotprodout)
      ,CHFt_CONST_FRA(afab)
      ,CHFt_CONST_FRA(bfab)
      ,CHFt_BOX(region)
      ,CHFt_CONST_INT(startcomp)
      ,CHFt_CONST_INT(endcomp) );
}
#endif  // GUARDDOTPRODUCT 

}

#endif
