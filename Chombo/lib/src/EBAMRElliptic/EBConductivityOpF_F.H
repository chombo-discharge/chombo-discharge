#ifndef _EBCONDUCTIVITYOPF_F_H_
#define _EBCONDUCTIVITYOPF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDCONDUCTIVITYGSRB 
#define GUARDCONDUCTIVITYGSRB 
// Prototype for Fortran procedure conductivitygsrb ...
//
void FORTRAN_NAME( CONDUCTIVITYGSRB ,conductivitygsrb )(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_CONST_FRA1(relcoef)
      ,CHFp_CONST_FRA1(acoef)
      ,CHFp_CONST_FRA1(b0)
      ,CHFp_CONST_FRA1(b1)
      ,CHFp_CONST_FRA1(b2)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(region)
      ,CHFp_CONST_INT(redBlack) );

#define FORT_CONDUCTIVITYGSRB FORTRAN_NAME( inlineCONDUCTIVITYGSRB, inlineCONDUCTIVITYGSRB)
#define FORTNT_CONDUCTIVITYGSRB FORTRAN_NAME( CONDUCTIVITYGSRB, conductivitygsrb)

inline void FORTRAN_NAME(inlineCONDUCTIVITYGSRB, inlineCONDUCTIVITYGSRB)(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_CONST_FRA1(relcoef)
      ,CHFp_CONST_FRA1(acoef)
      ,CHFp_CONST_FRA1(b0)
      ,CHFp_CONST_FRA1(b1)
      ,CHFp_CONST_FRA1(b2)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(region)
      ,CHFp_CONST_INT(redBlack) )
{
 CH_TIMELEAF("FORT_CONDUCTIVITYGSRB");
 FORTRAN_NAME( CONDUCTIVITYGSRB ,conductivitygsrb )(
      CHFt_FRA1(phi)
      ,CHFt_CONST_FRA1(rhs)
      ,CHFt_CONST_FRA1(relcoef)
      ,CHFt_CONST_FRA1(acoef)
      ,CHFt_CONST_FRA1(b0)
      ,CHFt_CONST_FRA1(b1)
      ,CHFt_CONST_FRA1(b2)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(region)
      ,CHFt_CONST_INT(redBlack) );
}
#endif  // GUARDCONDUCTIVITYGSRB 

#ifndef GUARDEBCOREGAPPLYDOMAINFLUX 
#define GUARDEBCOREGAPPLYDOMAINFLUX 
// Prototype for Fortran procedure ebcoregapplydomainflux ...
//
void FORTRAN_NAME( EBCOREGAPPLYDOMAINFLUX ,ebcoregapplydomainflux )(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(faceflux)
      ,CHFp_CONST_FRA1(bc)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(side)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_EBCOREGAPPLYDOMAINFLUX FORTRAN_NAME( inlineEBCOREGAPPLYDOMAINFLUX, inlineEBCOREGAPPLYDOMAINFLUX)
#define FORTNT_EBCOREGAPPLYDOMAINFLUX FORTRAN_NAME( EBCOREGAPPLYDOMAINFLUX, ebcoregapplydomainflux)

inline void FORTRAN_NAME(inlineEBCOREGAPPLYDOMAINFLUX, inlineEBCOREGAPPLYDOMAINFLUX)(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(faceflux)
      ,CHFp_CONST_FRA1(bc)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(side)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_EBCOREGAPPLYDOMAINFLUX");
 FORTRAN_NAME( EBCOREGAPPLYDOMAINFLUX ,ebcoregapplydomainflux )(
      CHFt_FRA1(phi)
      ,CHFt_CONST_FRA1(faceflux)
      ,CHFt_CONST_FRA1(bc)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(side)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDEBCOREGAPPLYDOMAINFLUX 

#ifndef GUARDAPPLYOPEBCONDNOBCS 
#define GUARDAPPLYOPEBCONDNOBCS 
// Prototype for Fortran procedure applyopebcondnobcs ...
//
void FORTRAN_NAME( APPLYOPEBCONDNOBCS ,applyopebcondnobcs )(
      CHFp_FRA1(opphidir)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(aco)
      ,CHFp_CONST_FRA1(b0)
      ,CHFp_CONST_FRA1(b1)
      ,CHFp_CONST_FRA1(b2)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_BOX(box) );

#define FORT_APPLYOPEBCONDNOBCS FORTRAN_NAME( inlineAPPLYOPEBCONDNOBCS, inlineAPPLYOPEBCONDNOBCS)
#define FORTNT_APPLYOPEBCONDNOBCS FORTRAN_NAME( APPLYOPEBCONDNOBCS, applyopebcondnobcs)

inline void FORTRAN_NAME(inlineAPPLYOPEBCONDNOBCS, inlineAPPLYOPEBCONDNOBCS)(
      CHFp_FRA1(opphidir)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(aco)
      ,CHFp_CONST_FRA1(b0)
      ,CHFp_CONST_FRA1(b1)
      ,CHFp_CONST_FRA1(b2)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_APPLYOPEBCONDNOBCS");
 FORTRAN_NAME( APPLYOPEBCONDNOBCS ,applyopebcondnobcs )(
      CHFt_FRA1(opphidir)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_FRA1(aco)
      ,CHFt_CONST_FRA1(b0)
      ,CHFt_CONST_FRA1(b1)
      ,CHFt_CONST_FRA1(b2)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_BOX(box) );
}
#endif  // GUARDAPPLYOPEBCONDNOBCS 

#ifndef GUARDGSCOLOREBCONDOP 
#define GUARDGSCOLOREBCONDOP 
// Prototype for Fortran procedure gscolorebcondop ...
//
void FORTRAN_NAME( GSCOLOREBCONDOP ,gscolorebcondop )(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_CONST_FRA1(relco)
      ,CHFp_CONST_FRA1(aco)
      ,CHFp_CONST_FRA1(b0)
      ,CHFp_CONST_FRA1(b1)
      ,CHFp_CONST_FRA1(b2)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_BOX(box) );

#define FORT_GSCOLOREBCONDOP FORTRAN_NAME( inlineGSCOLOREBCONDOP, inlineGSCOLOREBCONDOP)
#define FORTNT_GSCOLOREBCONDOP FORTRAN_NAME( GSCOLOREBCONDOP, gscolorebcondop)

inline void FORTRAN_NAME(inlineGSCOLOREBCONDOP, inlineGSCOLOREBCONDOP)(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_CONST_FRA1(relco)
      ,CHFp_CONST_FRA1(aco)
      ,CHFp_CONST_FRA1(b0)
      ,CHFp_CONST_FRA1(b1)
      ,CHFp_CONST_FRA1(b2)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_GSCOLOREBCONDOP");
 FORTRAN_NAME( GSCOLOREBCONDOP ,gscolorebcondop )(
      CHFt_FRA1(phi)
      ,CHFt_CONST_FRA1(rhs)
      ,CHFt_CONST_FRA1(relco)
      ,CHFt_CONST_FRA1(aco)
      ,CHFt_CONST_FRA1(b0)
      ,CHFt_CONST_FRA1(b1)
      ,CHFt_CONST_FRA1(b2)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_BOX(box) );
}
#endif  // GUARDGSCOLOREBCONDOP 

#ifndef GUARDCONDUCTIVITYINPLACE 
#define GUARDCONDUCTIVITYINPLACE 
// Prototype for Fortran procedure conductivityinplace ...
//
void FORTRAN_NAME( CONDUCTIVITYINPLACE ,conductivityinplace )(
      CHFp_FRA1(opphidir)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(b0)
      ,CHFp_CONST_FRA1(b1)
      ,CHFp_CONST_FRA1(b2)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_CONDUCTIVITYINPLACE FORTRAN_NAME( inlineCONDUCTIVITYINPLACE, inlineCONDUCTIVITYINPLACE)
#define FORTNT_CONDUCTIVITYINPLACE FORTRAN_NAME( CONDUCTIVITYINPLACE, conductivityinplace)

inline void FORTRAN_NAME(inlineCONDUCTIVITYINPLACE, inlineCONDUCTIVITYINPLACE)(
      CHFp_FRA1(opphidir)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(b0)
      ,CHFp_CONST_FRA1(b1)
      ,CHFp_CONST_FRA1(b2)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_CONDUCTIVITYINPLACE");
 FORTRAN_NAME( CONDUCTIVITYINPLACE ,conductivityinplace )(
      CHFt_FRA1(opphidir)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_FRA1(b0)
      ,CHFt_CONST_FRA1(b1)
      ,CHFt_CONST_FRA1(b2)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDCONDUCTIVITYINPLACE 

#ifndef GUARDINCRAPPLYEBCO 
#define GUARDINCRAPPLYEBCO 
// Prototype for Fortran procedure incrapplyebco ...
//
void FORTRAN_NAME( INCRAPPLYEBCO ,incrapplyebco )(
      CHFp_FRA1(lhs)
      ,CHFp_CONST_FRA1(interiorflux)
      ,CHFp_CONST_FRA1(domainfluxlo)
      ,CHFp_CONST_FRA1(domainfluxhi)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(lobox)
      ,CHFp_BOX(hibox)
      ,CHFp_BOX(centerbox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_CONST_INT(facedir) );

#define FORT_INCRAPPLYEBCO FORTRAN_NAME( inlineINCRAPPLYEBCO, inlineINCRAPPLYEBCO)
#define FORTNT_INCRAPPLYEBCO FORTRAN_NAME( INCRAPPLYEBCO, incrapplyebco)

inline void FORTRAN_NAME(inlineINCRAPPLYEBCO, inlineINCRAPPLYEBCO)(
      CHFp_FRA1(lhs)
      ,CHFp_CONST_FRA1(interiorflux)
      ,CHFp_CONST_FRA1(domainfluxlo)
      ,CHFp_CONST_FRA1(domainfluxhi)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(lobox)
      ,CHFp_BOX(hibox)
      ,CHFp_BOX(centerbox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_CONST_INT(facedir) )
{
 CH_TIMELEAF("FORT_INCRAPPLYEBCO");
 FORTRAN_NAME( INCRAPPLYEBCO ,incrapplyebco )(
      CHFt_FRA1(lhs)
      ,CHFt_CONST_FRA1(interiorflux)
      ,CHFt_CONST_FRA1(domainfluxlo)
      ,CHFt_CONST_FRA1(domainfluxhi)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(lobox)
      ,CHFt_BOX(hibox)
      ,CHFt_BOX(centerbox)
      ,CHFt_CONST_INT(haslo)
      ,CHFt_CONST_INT(hashi)
      ,CHFt_CONST_INT(facedir) );
}
#endif  // GUARDINCRAPPLYEBCO 

#ifndef GUARDDECRINVRELCOEFEBCO 
#define GUARDDECRINVRELCOEFEBCO 
// Prototype for Fortran procedure decrinvrelcoefebco ...
//
void FORTRAN_NAME( DECRINVRELCOEFEBCO ,decrinvrelcoefebco )(
      CHFp_FRA1(relcoef)
      ,CHFp_FRA1(bcoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_BOX(box)
      ,CHFp_REAL(dx)
      ,CHFp_INT(idir) );

#define FORT_DECRINVRELCOEFEBCO FORTRAN_NAME( inlineDECRINVRELCOEFEBCO, inlineDECRINVRELCOEFEBCO)
#define FORTNT_DECRINVRELCOEFEBCO FORTRAN_NAME( DECRINVRELCOEFEBCO, decrinvrelcoefebco)

inline void FORTRAN_NAME(inlineDECRINVRELCOEFEBCO, inlineDECRINVRELCOEFEBCO)(
      CHFp_FRA1(relcoef)
      ,CHFp_FRA1(bcoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_BOX(box)
      ,CHFp_REAL(dx)
      ,CHFp_INT(idir) )
{
 CH_TIMELEAF("FORT_DECRINVRELCOEFEBCO");
 FORTRAN_NAME( DECRINVRELCOEFEBCO ,decrinvrelcoefebco )(
      CHFt_FRA1(relcoef)
      ,CHFt_FRA1(bcoef)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_BOX(box)
      ,CHFt_REAL(dx)
      ,CHFt_INT(idir) );
}
#endif  // GUARDDECRINVRELCOEFEBCO 

#ifndef GUARDINVERTLAMBDAEBCO 
#define GUARDINVERTLAMBDAEBCO 
// Prototype for Fortran procedure invertlambdaebco ...
//
void FORTRAN_NAME( INVERTLAMBDAEBCO ,invertlambdaebco )(
      CHFp_FRA1(lambda)
      ,CHFp_REAL(safety)
      ,CHFp_BOX(box) );

#define FORT_INVERTLAMBDAEBCO FORTRAN_NAME( inlineINVERTLAMBDAEBCO, inlineINVERTLAMBDAEBCO)
#define FORTNT_INVERTLAMBDAEBCO FORTRAN_NAME( INVERTLAMBDAEBCO, invertlambdaebco)

inline void FORTRAN_NAME(inlineINVERTLAMBDAEBCO, inlineINVERTLAMBDAEBCO)(
      CHFp_FRA1(lambda)
      ,CHFp_REAL(safety)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_INVERTLAMBDAEBCO");
 FORTRAN_NAME( INVERTLAMBDAEBCO ,invertlambdaebco )(
      CHFt_FRA1(lambda)
      ,CHFt_REAL(safety)
      ,CHFt_BOX(box) );
}
#endif  // GUARDINVERTLAMBDAEBCO 

#ifndef GUARDGETFLUXEBCO 
#define GUARDGETFLUXEBCO 
// Prototype for Fortran procedure getfluxebco ...
//
void FORTRAN_NAME( GETFLUXEBCO ,getfluxebco )(
      CHFp_FRA1(flux)
      ,CHFp_CONST_FRA1(bcoef)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_BOX(opphibox)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(idir) );

#define FORT_GETFLUXEBCO FORTRAN_NAME( inlineGETFLUXEBCO, inlineGETFLUXEBCO)
#define FORTNT_GETFLUXEBCO FORTRAN_NAME( GETFLUXEBCO, getfluxebco)

inline void FORTRAN_NAME(inlineGETFLUXEBCO, inlineGETFLUXEBCO)(
      CHFp_FRA1(flux)
      ,CHFp_CONST_FRA1(bcoef)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_BOX(opphibox)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(idir) )
{
 CH_TIMELEAF("FORT_GETFLUXEBCO");
 FORTRAN_NAME( GETFLUXEBCO ,getfluxebco )(
      CHFt_FRA1(flux)
      ,CHFt_CONST_FRA1(bcoef)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_BOX(opphibox)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(idir) );
}
#endif  // GUARDGETFLUXEBCO 

#ifndef GUARDGSRBEBCO 
#define GUARDGSRBEBCO 
// Prototype for Fortran procedure gsrbebco ...
//
void FORTRAN_NAME( GSRBEBCO ,gsrbebco )(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(lph)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_CONST_FRA1(lam)
      ,CHFp_BOX(coloredbox) );

#define FORT_GSRBEBCO FORTRAN_NAME( inlineGSRBEBCO, inlineGSRBEBCO)
#define FORTNT_GSRBEBCO FORTRAN_NAME( GSRBEBCO, gsrbebco)

inline void FORTRAN_NAME(inlineGSRBEBCO, inlineGSRBEBCO)(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(lph)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_CONST_FRA1(lam)
      ,CHFp_BOX(coloredbox) )
{
 CH_TIMELEAF("FORT_GSRBEBCO");
 FORTRAN_NAME( GSRBEBCO ,gsrbebco )(
      CHFt_FRA1(phi)
      ,CHFt_CONST_FRA1(lph)
      ,CHFt_CONST_FRA1(rhs)
      ,CHFt_CONST_FRA1(lam)
      ,CHFt_BOX(coloredbox) );
}
#endif  // GUARDGSRBEBCO 

}

#endif
