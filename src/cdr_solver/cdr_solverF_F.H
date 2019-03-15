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
      ,CHFp_BOX(dcalc)
      ,CHFp_REAL(min_dt) );

#define FORT_ADVECTIVE_CFL FORTRAN_NAME( inlineADVECTIVE_CFL, inlineADVECTIVE_CFL)
#define FORTNT_ADVECTIVE_CFL FORTRAN_NAME( ADVECTIVE_CFL, advective_cfl)

inline void FORTRAN_NAME(inlineADVECTIVE_CFL, inlineADVECTIVE_CFL)(
      CHFp_CONST_FRA(velo)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(dcalc)
      ,CHFp_REAL(min_dt) )
{
 CH_TIMELEAF("FORT_ADVECTIVE_CFL");
 FORTRAN_NAME( ADVECTIVE_CFL ,advective_cfl )(
      CHFt_CONST_FRA(velo)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(dcalc)
      ,CHFt_REAL(min_dt) );
}
#endif  // GUARDADVECTIVE_CFL 

#ifndef GUARDSOURCE_DT 
#define GUARDSOURCE_DT 
// Prototype for Fortran procedure source_dt ...
//
void FORTRAN_NAME( SOURCE_DT ,source_dt )(
      CHFp_REAL(min_dt)
      ,CHFp_CONST_FRA1(state)
      ,CHFp_CONST_FRA1(source)
      ,CHFp_CONST_REAL(tolerance)
      ,CHFp_CONST_REAL(max_state)
      ,CHFp_BOX(region) );

#define FORT_SOURCE_DT FORTRAN_NAME( inlineSOURCE_DT, inlineSOURCE_DT)
#define FORTNT_SOURCE_DT FORTRAN_NAME( SOURCE_DT, source_dt)

inline void FORTRAN_NAME(inlineSOURCE_DT, inlineSOURCE_DT)(
      CHFp_REAL(min_dt)
      ,CHFp_CONST_FRA1(state)
      ,CHFp_CONST_FRA1(source)
      ,CHFp_CONST_REAL(tolerance)
      ,CHFp_CONST_REAL(max_state)
      ,CHFp_BOX(region) )
{
 CH_TIMELEAF("FORT_SOURCE_DT");
 FORTRAN_NAME( SOURCE_DT ,source_dt )(
      CHFt_REAL(min_dt)
      ,CHFt_CONST_FRA1(state)
      ,CHFt_CONST_FRA1(source)
      ,CHFt_CONST_REAL(tolerance)
      ,CHFt_CONST_REAL(max_state)
      ,CHFt_BOX(region) );
}
#endif  // GUARDSOURCE_DT 

#ifndef GUARDGET_MAX_VEL 
#define GUARDGET_MAX_VEL 
// Prototype for Fortran procedure get_max_vel ...
//
void FORTRAN_NAME( GET_MAX_VEL ,get_max_vel )(
      CHFp_REAL(max_vel)
      ,CHFp_CONST_FRA1(velo)
      ,CHFp_BOX(box)
      ,CHFp_CONST_FBA1(mask) );

#define FORT_GET_MAX_VEL FORTRAN_NAME( inlineGET_MAX_VEL, inlineGET_MAX_VEL)
#define FORTNT_GET_MAX_VEL FORTRAN_NAME( GET_MAX_VEL, get_max_vel)

inline void FORTRAN_NAME(inlineGET_MAX_VEL, inlineGET_MAX_VEL)(
      CHFp_REAL(max_vel)
      ,CHFp_CONST_FRA1(velo)
      ,CHFp_BOX(box)
      ,CHFp_CONST_FBA1(mask) )
{
 CH_TIMELEAF("FORT_GET_MAX_VEL");
 FORTRAN_NAME( GET_MAX_VEL ,get_max_vel )(
      CHFt_REAL(max_vel)
      ,CHFt_CONST_FRA1(velo)
      ,CHFt_BOX(box)
      ,CHFt_CONST_FBA1(mask) );
}
#endif  // GUARDGET_MAX_VEL 

#ifndef GUARDGET_MAXNORM 
#define GUARDGET_MAXNORM 
// Prototype for Fortran procedure get_maxnorm ...
//
void FORTRAN_NAME( GET_MAXNORM ,get_maxnorm )(
      CHFp_REAL(max_vel)
      ,CHFp_CONST_FRA1(velo)
      ,CHFp_BOX(box) );

#define FORT_GET_MAXNORM FORTRAN_NAME( inlineGET_MAXNORM, inlineGET_MAXNORM)
#define FORTNT_GET_MAXNORM FORTRAN_NAME( GET_MAXNORM, get_maxnorm)

inline void FORTRAN_NAME(inlineGET_MAXNORM, inlineGET_MAXNORM)(
      CHFp_REAL(max_vel)
      ,CHFp_CONST_FRA1(velo)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_GET_MAXNORM");
 FORTRAN_NAME( GET_MAXNORM ,get_maxnorm )(
      CHFt_REAL(max_vel)
      ,CHFt_CONST_FRA1(velo)
      ,CHFt_BOX(box) );
}
#endif  // GUARDGET_MAXNORM 

#ifndef GUARDADVECTIVE_CFL_DT 
#define GUARDADVECTIVE_CFL_DT 
// Prototype for Fortran procedure advective_cfl_dt ...
//
void FORTRAN_NAME( ADVECTIVE_CFL_DT ,advective_cfl_dt )(
      CHFp_CONST_FRA(velo)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(dcalc)
      ,CHFp_REAL(min_dt) );

#define FORT_ADVECTIVE_CFL_DT FORTRAN_NAME( inlineADVECTIVE_CFL_DT, inlineADVECTIVE_CFL_DT)
#define FORTNT_ADVECTIVE_CFL_DT FORTRAN_NAME( ADVECTIVE_CFL_DT, advective_cfl_dt)

inline void FORTRAN_NAME(inlineADVECTIVE_CFL_DT, inlineADVECTIVE_CFL_DT)(
      CHFp_CONST_FRA(velo)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(dcalc)
      ,CHFp_REAL(min_dt) )
{
 CH_TIMELEAF("FORT_ADVECTIVE_CFL_DT");
 FORTRAN_NAME( ADVECTIVE_CFL_DT ,advective_cfl_dt )(
      CHFt_CONST_FRA(velo)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(dcalc)
      ,CHFt_REAL(min_dt) );
}
#endif  // GUARDADVECTIVE_CFL_DT 

#ifndef GUARDUPWIND_SOURCE 
#define GUARDUPWIND_SOURCE 
// Prototype for Fortran procedure upwind_source ...
//
void FORTRAN_NAME( UPWIND_SOURCE ,upwind_source )(
      CHFp_FRA1(dst)
      ,CHFp_CONST_FRA1(src)
      ,CHFp_CONST_FRA1(vel)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(dcalc) );

#define FORT_UPWIND_SOURCE FORTRAN_NAME( inlineUPWIND_SOURCE, inlineUPWIND_SOURCE)
#define FORTNT_UPWIND_SOURCE FORTRAN_NAME( UPWIND_SOURCE, upwind_source)

inline void FORTRAN_NAME(inlineUPWIND_SOURCE, inlineUPWIND_SOURCE)(
      CHFp_FRA1(dst)
      ,CHFp_CONST_FRA1(src)
      ,CHFp_CONST_FRA1(vel)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(dcalc) )
{
 CH_TIMELEAF("FORT_UPWIND_SOURCE");
 FORTRAN_NAME( UPWIND_SOURCE ,upwind_source )(
      CHFt_FRA1(dst)
      ,CHFt_CONST_FRA1(src)
      ,CHFt_CONST_FRA1(vel)
      ,CHFt_CONST_INT(dir)
      ,CHFt_BOX(dcalc) );
}
#endif  // GUARDUPWIND_SOURCE 

}

#endif
