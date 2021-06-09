#ifndef _CD_EBHELMHOLTZOPF_F_H_
#define _CD_EBHELMHOLTZOPF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDADDBCOTERMTOINVRELCOEF 
#define GUARDADDBCOTERMTOINVRELCOEF 
// Prototype for Fortran procedure addbcotermtoinvrelcoef ...
//
void FORTRAN_NAME( ADDBCOTERMTOINVRELCOEF ,addbcotermtoinvrelcoef )(
      CHFp_FRA1(relcoef)
      ,CHFp_CONST_FRA1(bcoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(cellbox) );

#define FORT_ADDBCOTERMTOINVRELCOEF FORTRAN_NAME( inlineADDBCOTERMTOINVRELCOEF, inlineADDBCOTERMTOINVRELCOEF)
#define FORTNT_ADDBCOTERMTOINVRELCOEF FORTRAN_NAME( ADDBCOTERMTOINVRELCOEF, addbcotermtoinvrelcoef)

inline void FORTRAN_NAME(inlineADDBCOTERMTOINVRELCOEF, inlineADDBCOTERMTOINVRELCOEF)(
      CHFp_FRA1(relcoef)
      ,CHFp_CONST_FRA1(bcoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(cellbox) )
{
 CH_TIMELEAF("FORT_ADDBCOTERMTOINVRELCOEF");
 FORTRAN_NAME( ADDBCOTERMTOINVRELCOEF ,addbcotermtoinvrelcoef )(
      CHFt_FRA1(relcoef)
      ,CHFt_CONST_FRA1(bcoef)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(cellbox) );
}
#endif  // GUARDADDBCOTERMTOINVRELCOEF 

#ifndef GUARDINVERTRELAXATIONCOEFFICIENT 
#define GUARDINVERTRELAXATIONCOEFFICIENT 
// Prototype for Fortran procedure invertrelaxationcoefficient ...
//
void FORTRAN_NAME( INVERTRELAXATIONCOEFFICIENT ,invertrelaxationcoefficient )(
      CHFp_FRA1(relcoef)
      ,CHFp_BOX(cellbox) );

#define FORT_INVERTRELAXATIONCOEFFICIENT FORTRAN_NAME( inlineINVERTRELAXATIONCOEFFICIENT, inlineINVERTRELAXATIONCOEFFICIENT)
#define FORTNT_INVERTRELAXATIONCOEFFICIENT FORTRAN_NAME( INVERTRELAXATIONCOEFFICIENT, invertrelaxationcoefficient)

inline void FORTRAN_NAME(inlineINVERTRELAXATIONCOEFFICIENT, inlineINVERTRELAXATIONCOEFFICIENT)(
      CHFp_FRA1(relcoef)
      ,CHFp_BOX(cellbox) )
{
 CH_TIMELEAF("FORT_INVERTRELAXATIONCOEFFICIENT");
 FORTRAN_NAME( INVERTRELAXATIONCOEFFICIENT ,invertrelaxationcoefficient )(
      CHFt_FRA1(relcoef)
      ,CHFt_BOX(cellbox) );
}
#endif  // GUARDINVERTRELAXATIONCOEFFICIENT 

#ifndef GUARDGETINTERIORREGFLUX 
#define GUARDGETINTERIORREGFLUX 
// Prototype for Fortran procedure getinteriorregflux ...
//
void FORTRAN_NAME( GETINTERIORREGFLUX ,getinteriorregflux )(
      CHFp_FRA1(regflux)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(bco)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(facebox) );

#define FORT_GETINTERIORREGFLUX FORTRAN_NAME( inlineGETINTERIORREGFLUX, inlineGETINTERIORREGFLUX)
#define FORTNT_GETINTERIORREGFLUX FORTRAN_NAME( GETINTERIORREGFLUX, getinteriorregflux)

inline void FORTRAN_NAME(inlineGETINTERIORREGFLUX, inlineGETINTERIORREGFLUX)(
      CHFp_FRA1(regflux)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(bco)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(facebox) )
{
 CH_TIMELEAF("FORT_GETINTERIORREGFLUX");
 FORTRAN_NAME( GETINTERIORREGFLUX ,getinteriorregflux )(
      CHFt_FRA1(regflux)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_FRA1(bco)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(facebox) );
}
#endif  // GUARDGETINTERIORREGFLUX 

#ifndef GUARDLAPLACIANINPLACE 
#define GUARDLAPLACIANINPLACE 
// Prototype for Fortran procedure laplacianinplace ...
//
void FORTRAN_NAME( LAPLACIANINPLACE ,laplacianinplace )(
      CHFp_FRA1(Lphi)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(bcoX)
      ,CHFp_CONST_FRA1(bcoY)
      ,CHFp_CONST_FRA1(bcoZ)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_BOX(cellBox) );

#define FORT_LAPLACIANINPLACE FORTRAN_NAME( inlineLAPLACIANINPLACE, inlineLAPLACIANINPLACE)
#define FORTNT_LAPLACIANINPLACE FORTRAN_NAME( LAPLACIANINPLACE, laplacianinplace)

inline void FORTRAN_NAME(inlineLAPLACIANINPLACE, inlineLAPLACIANINPLACE)(
      CHFp_FRA1(Lphi)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(bcoX)
      ,CHFp_CONST_FRA1(bcoY)
      ,CHFp_CONST_FRA1(bcoZ)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_BOX(cellBox) )
{
 CH_TIMELEAF("FORT_LAPLACIANINPLACE");
 FORTRAN_NAME( LAPLACIANINPLACE ,laplacianinplace )(
      CHFt_FRA1(Lphi)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_FRA1(bcoX)
      ,CHFt_CONST_FRA1(bcoY)
      ,CHFt_CONST_FRA1(bcoZ)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_BOX(cellBox) );
}
#endif  // GUARDLAPLACIANINPLACE 

#ifndef GUARDHELMHOLTZAPPLYDOMAINFLUX 
#define GUARDHELMHOLTZAPPLYDOMAINFLUX 
// Prototype for Fortran procedure helmholtzapplydomainflux ...
//
void FORTRAN_NAME( HELMHOLTZAPPLYDOMAINFLUX ,helmholtzapplydomainflux )(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(flux)
      ,CHFp_CONST_FRA1(bco)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(side)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(cellbox) );

#define FORT_HELMHOLTZAPPLYDOMAINFLUX FORTRAN_NAME( inlineHELMHOLTZAPPLYDOMAINFLUX, inlineHELMHOLTZAPPLYDOMAINFLUX)
#define FORTNT_HELMHOLTZAPPLYDOMAINFLUX FORTRAN_NAME( HELMHOLTZAPPLYDOMAINFLUX, helmholtzapplydomainflux)

inline void FORTRAN_NAME(inlineHELMHOLTZAPPLYDOMAINFLUX, inlineHELMHOLTZAPPLYDOMAINFLUX)(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(flux)
      ,CHFp_CONST_FRA1(bco)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(side)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(cellbox) )
{
 CH_TIMELEAF("FORT_HELMHOLTZAPPLYDOMAINFLUX");
 FORTRAN_NAME( HELMHOLTZAPPLYDOMAINFLUX ,helmholtzapplydomainflux )(
      CHFt_FRA1(phi)
      ,CHFt_CONST_FRA1(flux)
      ,CHFt_CONST_FRA1(bco)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(side)
      ,CHFt_CONST_INT(dir)
      ,CHFt_BOX(cellbox) );
}
#endif  // GUARDHELMHOLTZAPPLYDOMAINFLUX 

#ifndef GUARDHELMHOLTZGSRBCOLOR 
#define GUARDHELMHOLTZGSRBCOLOR 
// Prototype for Fortran procedure helmholtzgsrbcolor ...
//
void FORTRAN_NAME( HELMHOLTZGSRBCOLOR ,helmholtzgsrbcolor )(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(lphi)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_CONST_FRA1(relCoef)
      ,CHFp_BOX(coloredBox) );

#define FORT_HELMHOLTZGSRBCOLOR FORTRAN_NAME( inlineHELMHOLTZGSRBCOLOR, inlineHELMHOLTZGSRBCOLOR)
#define FORTNT_HELMHOLTZGSRBCOLOR FORTRAN_NAME( HELMHOLTZGSRBCOLOR, helmholtzgsrbcolor)

inline void FORTRAN_NAME(inlineHELMHOLTZGSRBCOLOR, inlineHELMHOLTZGSRBCOLOR)(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(lphi)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_CONST_FRA1(relCoef)
      ,CHFp_BOX(coloredBox) )
{
 CH_TIMELEAF("FORT_HELMHOLTZGSRBCOLOR");
 FORTRAN_NAME( HELMHOLTZGSRBCOLOR ,helmholtzgsrbcolor )(
      CHFt_FRA1(phi)
      ,CHFt_CONST_FRA1(lphi)
      ,CHFt_CONST_FRA1(rhs)
      ,CHFt_CONST_FRA1(relCoef)
      ,CHFt_BOX(coloredBox) );
}
#endif  // GUARDHELMHOLTZGSRBCOLOR 

}

#endif
