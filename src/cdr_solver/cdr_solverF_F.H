#ifndef _CDR_SOLVERF_F_H_
#define _CDR_SOLVERF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDADVECTIVEDERIV 
#define GUARDADVECTIVEDERIV 
// Prototype for Fortran procedure advectivederiv ...
//
void FORTRAN_NAME( ADVECTIVEDERIV ,advectivederiv )(
      CHFp_FRA1(udelrho)
      ,CHFp_CONST_FRA1(facerho)
      ,CHFp_CONST_FRA1(facevel)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(ncomp)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(dcalc) );

#define FORT_ADVECTIVEDERIV FORTRAN_NAME( inlineADVECTIVEDERIV, inlineADVECTIVEDERIV)
#define FORTNT_ADVECTIVEDERIV FORTRAN_NAME( ADVECTIVEDERIV, advectivederiv)

inline void FORTRAN_NAME(inlineADVECTIVEDERIV, inlineADVECTIVEDERIV)(
      CHFp_FRA1(udelrho)
      ,CHFp_CONST_FRA1(facerho)
      ,CHFp_CONST_FRA1(facevel)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(ncomp)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(dcalc) )
{
 CH_TIMELEAF("FORT_ADVECTIVEDERIV");
 FORTRAN_NAME( ADVECTIVEDERIV ,advectivederiv )(
      CHFt_FRA1(udelrho)
      ,CHFt_CONST_FRA1(facerho)
      ,CHFt_CONST_FRA1(facevel)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_CONST_INT(ncomp)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(dcalc) );
}
#endif  // GUARDADVECTIVEDERIV 

#ifndef GUARDCONSDIV_REG 
#define GUARDCONSDIV_REG 
// Prototype for Fortran procedure consdiv_reg ...
//
void FORTRAN_NAME( CONSDIV_REG ,consdiv_reg )(
      CHFp_FRA1(udelrho)
      ,CHFp_CONST_FRA1(faceflux)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(ncomp)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(dcalc) );

#define FORT_CONSDIV_REG FORTRAN_NAME( inlineCONSDIV_REG, inlineCONSDIV_REG)
#define FORTNT_CONSDIV_REG FORTRAN_NAME( CONSDIV_REG, consdiv_reg)

inline void FORTRAN_NAME(inlineCONSDIV_REG, inlineCONSDIV_REG)(
      CHFp_FRA1(udelrho)
      ,CHFp_CONST_FRA1(faceflux)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(ncomp)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(dcalc) )
{
 CH_TIMELEAF("FORT_CONSDIV_REG");
 FORTRAN_NAME( CONSDIV_REG ,consdiv_reg )(
      CHFt_FRA1(udelrho)
      ,CHFt_CONST_FRA1(faceflux)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_CONST_INT(ncomp)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(dcalc) );
}
#endif  // GUARDCONSDIV_REG 

#ifndef GUARDADVECTIVE_CFL 
#define GUARDADVECTIVE_CFL 
// Prototype for Fortran procedure advective_cfl ...
//
void FORTRAN_NAME( ADVECTIVE_CFL ,advective_cfl )(
      CHFp_CONST_FRA(velo)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_REAL(regionchf_realmin_dt) );

#define FORT_ADVECTIVE_CFL FORTRAN_NAME( inlineADVECTIVE_CFL, inlineADVECTIVE_CFL)
#define FORTNT_ADVECTIVE_CFL FORTRAN_NAME( ADVECTIVE_CFL, advective_cfl)

inline void FORTRAN_NAME(inlineADVECTIVE_CFL, inlineADVECTIVE_CFL)(
      CHFp_CONST_FRA(velo)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_REAL(regionchf_realmin_dt) )
{
 CH_TIMELEAF("FORT_ADVECTIVE_CFL");
 FORTRAN_NAME( ADVECTIVE_CFL ,advective_cfl )(
      CHFt_CONST_FRA(velo)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_REAL(regionchf_realmin_dt) );
}
#endif  // GUARDADVECTIVE_CFL 

}

#endif
