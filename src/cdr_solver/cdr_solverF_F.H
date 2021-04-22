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
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(dcalc) );

#define FORT_CONSDIV_REG FORTRAN_NAME( inlineCONSDIV_REG, inlineCONSDIV_REG)
#define FORTNT_CONSDIV_REG FORTRAN_NAME( CONSDIV_REG, consdiv_reg)

inline void FORTRAN_NAME(inlineCONSDIV_REG, inlineCONSDIV_REG)(
      CHFp_FRA1(udelrho)
      ,CHFp_CONST_FRA1(faceflux)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(dcalc) )
{
 CH_TIMELEAF("FORT_CONSDIV_REG");
 FORTRAN_NAME( CONSDIV_REG ,consdiv_reg )(
      CHFt_FRA1(udelrho)
      ,CHFt_CONST_FRA1(faceflux)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(dcalc) );
}
#endif  // GUARDCONSDIV_REG 

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

#ifndef GUARDADVECTION_DT 
#define GUARDADVECTION_DT 
// Prototype for Fortran procedure advection_dt ...
//
void FORTRAN_NAME( ADVECTION_DT ,advection_dt )(
      CHFp_FRA1(min_dt)
      ,CHFp_CONST_FRA(velo)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(dcalc) );

#define FORT_ADVECTION_DT FORTRAN_NAME( inlineADVECTION_DT, inlineADVECTION_DT)
#define FORTNT_ADVECTION_DT FORTRAN_NAME( ADVECTION_DT, advection_dt)

inline void FORTRAN_NAME(inlineADVECTION_DT, inlineADVECTION_DT)(
      CHFp_FRA1(min_dt)
      ,CHFp_CONST_FRA(velo)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(dcalc) )
{
 CH_TIMELEAF("FORT_ADVECTION_DT");
 FORTRAN_NAME( ADVECTION_DT ,advection_dt )(
      CHFt_FRA1(min_dt)
      ,CHFt_CONST_FRA(velo)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(dcalc) );
}
#endif  // GUARDADVECTION_DT 

#ifndef GUARDDIFFUSION_DT 
#define GUARDDIFFUSION_DT 
// Prototype for Fortran procedure diffusion_dt ...
//
void FORTRAN_NAME( DIFFUSION_DT ,diffusion_dt )(
      CHFp_CONST_FRA1(diffco)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(facebox)
      ,CHFp_REAL(min_dt) );

#define FORT_DIFFUSION_DT FORTRAN_NAME( inlineDIFFUSION_DT, inlineDIFFUSION_DT)
#define FORTNT_DIFFUSION_DT FORTRAN_NAME( DIFFUSION_DT, diffusion_dt)

inline void FORTRAN_NAME(inlineDIFFUSION_DT, inlineDIFFUSION_DT)(
      CHFp_CONST_FRA1(diffco)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(facebox)
      ,CHFp_REAL(min_dt) )
{
 CH_TIMELEAF("FORT_DIFFUSION_DT");
 FORTRAN_NAME( DIFFUSION_DT ,diffusion_dt )(
      CHFt_CONST_FRA1(diffco)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(facebox)
      ,CHFt_REAL(min_dt) );
}
#endif  // GUARDDIFFUSION_DT 

#ifndef GUARDCONSDIV_BLOCK 
#define GUARDCONSDIV_BLOCK 
// Prototype for Fortran procedure consdiv_block ...
//
void FORTRAN_NAME( CONSDIV_BLOCK ,consdiv_block )(
      CHFp_FRA1(udelrho)
      ,CHFp_CONST_FRA1(faceflux)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(tilesize)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(dcalc) );

#define FORT_CONSDIV_BLOCK FORTRAN_NAME( inlineCONSDIV_BLOCK, inlineCONSDIV_BLOCK)
#define FORTNT_CONSDIV_BLOCK FORTRAN_NAME( CONSDIV_BLOCK, consdiv_block)

inline void FORTRAN_NAME(inlineCONSDIV_BLOCK, inlineCONSDIV_BLOCK)(
      CHFp_FRA1(udelrho)
      ,CHFp_CONST_FRA1(faceflux)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(tilesize)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(dcalc) )
{
 CH_TIMELEAF("FORT_CONSDIV_BLOCK");
 FORTRAN_NAME( CONSDIV_BLOCK ,consdiv_block )(
      CHFt_FRA1(udelrho)
      ,CHFt_CONST_FRA1(faceflux)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_CONST_INT(tilesize)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(dcalc) );
}
#endif  // GUARDCONSDIV_BLOCK 

#ifndef GUARDDFLUX_REG 
#define GUARDDFLUX_REG 
// Prototype for Fortran procedure dflux_reg ...
//
void FORTRAN_NAME( DFLUX_REG ,dflux_reg )(
      CHFp_FRA1(divD)
      ,CHFp_CONST_FRA1(state)
      ,CHFp_CONST_FRA1(dco)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(facebox) );

#define FORT_DFLUX_REG FORTRAN_NAME( inlineDFLUX_REG, inlineDFLUX_REG)
#define FORTNT_DFLUX_REG FORTRAN_NAME( DFLUX_REG, dflux_reg)

inline void FORTRAN_NAME(inlineDFLUX_REG, inlineDFLUX_REG)(
      CHFp_FRA1(divD)
      ,CHFp_CONST_FRA1(state)
      ,CHFp_CONST_FRA1(dco)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(facebox) )
{
 CH_TIMELEAF("FORT_DFLUX_REG");
 FORTRAN_NAME( DFLUX_REG ,dflux_reg )(
      CHFt_FRA1(divD)
      ,CHFt_CONST_FRA1(state)
      ,CHFt_CONST_FRA1(dco)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(facebox) );
}
#endif  // GUARDDFLUX_REG 

}

#endif
