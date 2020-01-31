#ifndef _EBVISCOUSTENSOROPF_F_H_
#define _EBVISCOUSTENSOROPF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDNWOEBVTOPOINTLPH 
#define GUARDNWOEBVTOPOINTLPH 
// Prototype for Fortran procedure nwoebvtopointlph ...
//
void FORTRAN_NAME( NWOEBVTOPOINTLPH ,nwoebvtopointlph )(
      CHFp_REALVECT(lphi)
      ,CHFp_CONST_INTVECT(iv)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_FRA1(acofab)
      ,CHFp_CONST_FRA1(eta0fab)
      ,CHFp_CONST_FRA1(eta1fab)
      ,CHFp_CONST_FRA1(eta2fab)
      ,CHFp_CONST_FRA1(lam0fab)
      ,CHFp_CONST_FRA1(lam1fab)
      ,CHFp_CONST_FRA1(lam2fab)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta) );

#define FORT_NWOEBVTOPOINTLPH FORTRAN_NAME( inlineNWOEBVTOPOINTLPH, inlineNWOEBVTOPOINTLPH)
#define FORTNT_NWOEBVTOPOINTLPH FORTRAN_NAME( NWOEBVTOPOINTLPH, nwoebvtopointlph)

inline void FORTRAN_NAME(inlineNWOEBVTOPOINTLPH, inlineNWOEBVTOPOINTLPH)(
      CHFp_REALVECT(lphi)
      ,CHFp_CONST_INTVECT(iv)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_FRA1(acofab)
      ,CHFp_CONST_FRA1(eta0fab)
      ,CHFp_CONST_FRA1(eta1fab)
      ,CHFp_CONST_FRA1(eta2fab)
      ,CHFp_CONST_FRA1(lam0fab)
      ,CHFp_CONST_FRA1(lam1fab)
      ,CHFp_CONST_FRA1(lam2fab)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta) )
{
 CH_TIMELEAF("FORT_NWOEBVTOPOINTLPH");
 FORTRAN_NAME( NWOEBVTOPOINTLPH ,nwoebvtopointlph )(
      CHFt_REALVECT(lphi)
      ,CHFt_CONST_INTVECT(iv)
      ,CHFt_CONST_FRA(phi)
      ,CHFt_CONST_FRA1(acofab)
      ,CHFt_CONST_FRA1(eta0fab)
      ,CHFt_CONST_FRA1(eta1fab)
      ,CHFt_CONST_FRA1(eta2fab)
      ,CHFt_CONST_FRA1(lam0fab)
      ,CHFt_CONST_FRA1(lam1fab)
      ,CHFt_CONST_FRA1(lam2fab)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta) );
}
#endif  // GUARDNWOEBVTOPOINTLPH 

#ifndef GUARDGSRBNWOEBVTOP 
#define GUARDGSRBNWOEBVTOP 
// Prototype for Fortran procedure gsrbnwoebvtop ...
//
void FORTRAN_NAME( GSRBNWOEBVTOP ,gsrbnwoebvtop )(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_FRA(lambda)
      ,CHFp_CONST_FRA1(acofab)
      ,CHFp_CONST_FRA1(eta0fab)
      ,CHFp_CONST_FRA1(eta1fab)
      ,CHFp_CONST_FRA1(eta2fab)
      ,CHFp_CONST_FRA1(lam0fab)
      ,CHFp_CONST_FRA1(lam1fab)
      ,CHFp_CONST_FRA1(lam2fab)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_BOX(coloredbox) );

#define FORT_GSRBNWOEBVTOP FORTRAN_NAME( inlineGSRBNWOEBVTOP, inlineGSRBNWOEBVTOP)
#define FORTNT_GSRBNWOEBVTOP FORTRAN_NAME( GSRBNWOEBVTOP, gsrbnwoebvtop)

inline void FORTRAN_NAME(inlineGSRBNWOEBVTOP, inlineGSRBNWOEBVTOP)(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_FRA(lambda)
      ,CHFp_CONST_FRA1(acofab)
      ,CHFp_CONST_FRA1(eta0fab)
      ,CHFp_CONST_FRA1(eta1fab)
      ,CHFp_CONST_FRA1(eta2fab)
      ,CHFp_CONST_FRA1(lam0fab)
      ,CHFp_CONST_FRA1(lam1fab)
      ,CHFp_CONST_FRA1(lam2fab)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_BOX(coloredbox) )
{
 CH_TIMELEAF("FORT_GSRBNWOEBVTOP");
 FORTRAN_NAME( GSRBNWOEBVTOP ,gsrbnwoebvtop )(
      CHFt_FRA(phi)
      ,CHFt_CONST_FRA(rhs)
      ,CHFt_CONST_FRA(lambda)
      ,CHFt_CONST_FRA1(acofab)
      ,CHFt_CONST_FRA1(eta0fab)
      ,CHFt_CONST_FRA1(eta1fab)
      ,CHFt_CONST_FRA1(eta2fab)
      ,CHFt_CONST_FRA1(lam0fab)
      ,CHFt_CONST_FRA1(lam1fab)
      ,CHFt_CONST_FRA1(lam2fab)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_BOX(coloredbox) );
}
#endif  // GUARDGSRBNWOEBVTOP 

#ifndef GUARDCELLGRADEBVTOP 
#define GUARDCELLGRADEBVTOP 
// Prototype for Fortran procedure cellgradebvtop ...
//
void FORTRAN_NAME( CELLGRADEBVTOP ,cellgradebvtop )(
      CHFp_FRA1(grad)
      ,CHFp_CONST_FRA1(vel)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(lobox)
      ,CHFp_BOX(hibox)
      ,CHFp_BOX(centerbox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_CONST_INT(divdir) );

#define FORT_CELLGRADEBVTOP FORTRAN_NAME( inlineCELLGRADEBVTOP, inlineCELLGRADEBVTOP)
#define FORTNT_CELLGRADEBVTOP FORTRAN_NAME( CELLGRADEBVTOP, cellgradebvtop)

inline void FORTRAN_NAME(inlineCELLGRADEBVTOP, inlineCELLGRADEBVTOP)(
      CHFp_FRA1(grad)
      ,CHFp_CONST_FRA1(vel)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(lobox)
      ,CHFp_BOX(hibox)
      ,CHFp_BOX(centerbox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_CONST_INT(divdir) )
{
 CH_TIMELEAF("FORT_CELLGRADEBVTOP");
 FORTRAN_NAME( CELLGRADEBVTOP ,cellgradebvtop )(
      CHFt_FRA1(grad)
      ,CHFt_CONST_FRA1(vel)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(lobox)
      ,CHFt_BOX(hibox)
      ,CHFt_BOX(centerbox)
      ,CHFt_CONST_INT(haslo)
      ,CHFt_CONST_INT(hashi)
      ,CHFt_CONST_INT(divdir) );
}
#endif  // GUARDCELLGRADEBVTOP 

#ifndef GUARDINCRAPPLYEBVTOP 
#define GUARDINCRAPPLYEBVTOP 
// Prototype for Fortran procedure incrapplyebvtop ...
//
void FORTRAN_NAME( INCRAPPLYEBVTOP ,incrapplyebvtop )(
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

#define FORT_INCRAPPLYEBVTOP FORTRAN_NAME( inlineINCRAPPLYEBVTOP, inlineINCRAPPLYEBVTOP)
#define FORTNT_INCRAPPLYEBVTOP FORTRAN_NAME( INCRAPPLYEBVTOP, incrapplyebvtop)

inline void FORTRAN_NAME(inlineINCRAPPLYEBVTOP, inlineINCRAPPLYEBVTOP)(
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
 CH_TIMELEAF("FORT_INCRAPPLYEBVTOP");
 FORTRAN_NAME( INCRAPPLYEBVTOP ,incrapplyebvtop )(
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
#endif  // GUARDINCRAPPLYEBVTOP 

#ifndef GUARDNORMALGRADVISCDIRCH 
#define GUARDNORMALGRADVISCDIRCH 
// Prototype for Fortran procedure normalgradviscdirch ...
//
void FORTRAN_NAME( NORMALGRADVISCDIRCH ,normalgradviscdirch )(
      CHFp_FRA1(gradvelface)
      ,CHFp_CONST_FRA1(velcomp)
      ,CHFp_CONST_FRA1(velvalu)
      ,CHFp_BOX(facebox)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(iside)
      ,CHFp_CONST_INT(divdir) );

#define FORT_NORMALGRADVISCDIRCH FORTRAN_NAME( inlineNORMALGRADVISCDIRCH, inlineNORMALGRADVISCDIRCH)
#define FORTNT_NORMALGRADVISCDIRCH FORTRAN_NAME( NORMALGRADVISCDIRCH, normalgradviscdirch)

inline void FORTRAN_NAME(inlineNORMALGRADVISCDIRCH, inlineNORMALGRADVISCDIRCH)(
      CHFp_FRA1(gradvelface)
      ,CHFp_CONST_FRA1(velcomp)
      ,CHFp_CONST_FRA1(velvalu)
      ,CHFp_BOX(facebox)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(iside)
      ,CHFp_CONST_INT(divdir) )
{
 CH_TIMELEAF("FORT_NORMALGRADVISCDIRCH");
 FORTRAN_NAME( NORMALGRADVISCDIRCH ,normalgradviscdirch )(
      CHFt_FRA1(gradvelface)
      ,CHFt_CONST_FRA1(velcomp)
      ,CHFt_CONST_FRA1(velvalu)
      ,CHFt_BOX(facebox)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(iside)
      ,CHFt_CONST_INT(divdir) );
}
#endif  // GUARDNORMALGRADVISCDIRCH 

#ifndef GUARDSLIPWALLGRAD 
#define GUARDSLIPWALLGRAD 
// Prototype for Fortran procedure slipwallgrad ...
//
void FORTRAN_NAME( SLIPWALLGRAD ,slipwallgrad )(
      CHFp_FRA1(gradvelface)
      ,CHFp_CONST_FRA1(velcomp)
      ,CHFp_BOX(facebox)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(veldir)
      ,CHFp_CONST_INT(divdir)
      ,CHFp_CONST_INT(iside)
      ,CHFp_BOX(lobox)
      ,CHFp_BOX(hibox)
      ,CHFp_BOX(centerbox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_CONST_INT(facedir) );

#define FORT_SLIPWALLGRAD FORTRAN_NAME( inlineSLIPWALLGRAD, inlineSLIPWALLGRAD)
#define FORTNT_SLIPWALLGRAD FORTRAN_NAME( SLIPWALLGRAD, slipwallgrad)

inline void FORTRAN_NAME(inlineSLIPWALLGRAD, inlineSLIPWALLGRAD)(
      CHFp_FRA1(gradvelface)
      ,CHFp_CONST_FRA1(velcomp)
      ,CHFp_BOX(facebox)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(veldir)
      ,CHFp_CONST_INT(divdir)
      ,CHFp_CONST_INT(iside)
      ,CHFp_BOX(lobox)
      ,CHFp_BOX(hibox)
      ,CHFp_BOX(centerbox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_CONST_INT(facedir) )
{
 CH_TIMELEAF("FORT_SLIPWALLGRAD");
 FORTRAN_NAME( SLIPWALLGRAD ,slipwallgrad )(
      CHFt_FRA1(gradvelface)
      ,CHFt_CONST_FRA1(velcomp)
      ,CHFt_BOX(facebox)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(veldir)
      ,CHFt_CONST_INT(divdir)
      ,CHFt_CONST_INT(iside)
      ,CHFt_BOX(lobox)
      ,CHFt_BOX(hibox)
      ,CHFt_BOX(centerbox)
      ,CHFt_CONST_INT(haslo)
      ,CHFt_CONST_INT(hashi)
      ,CHFt_CONST_INT(facedir) );
}
#endif  // GUARDSLIPWALLGRAD 

#ifndef GUARDVELDOTSIGMA 
#define GUARDVELDOTSIGMA 
// Prototype for Fortran procedure veldotsigma ...
//
void FORTRAN_NAME( VELDOTSIGMA ,veldotsigma )(
      CHFp_FRA(veldotsig)
      ,CHFp_CONST_FRA(vel)
      ,CHFp_CONST_FRA(sig)
      ,CHFp_BOX(facebox) );

#define FORT_VELDOTSIGMA FORTRAN_NAME( inlineVELDOTSIGMA, inlineVELDOTSIGMA)
#define FORTNT_VELDOTSIGMA FORTRAN_NAME( VELDOTSIGMA, veldotsigma)

inline void FORTRAN_NAME(inlineVELDOTSIGMA, inlineVELDOTSIGMA)(
      CHFp_FRA(veldotsig)
      ,CHFp_CONST_FRA(vel)
      ,CHFp_CONST_FRA(sig)
      ,CHFp_BOX(facebox) )
{
 CH_TIMELEAF("FORT_VELDOTSIGMA");
 FORTRAN_NAME( VELDOTSIGMA ,veldotsigma )(
      CHFt_FRA(veldotsig)
      ,CHFt_CONST_FRA(vel)
      ,CHFt_CONST_FRA(sig)
      ,CHFt_BOX(facebox) );
}
#endif  // GUARDVELDOTSIGMA 

#ifndef GUARDINCRPOINTDOTPROD 
#define GUARDINCRPOINTDOTPROD 
// Prototype for Fortran procedure incrpointdotprod ...
//
void FORTRAN_NAME( INCRPOINTDOTPROD ,incrpointdotprod )(
      CHFp_FRA1(sigmadotgradu)
      ,CHFp_CONST_FRA1(gradu)
      ,CHFp_CONST_FRA1(sigma)
      ,CHFp_BOX(box) );

#define FORT_INCRPOINTDOTPROD FORTRAN_NAME( inlineINCRPOINTDOTPROD, inlineINCRPOINTDOTPROD)
#define FORTNT_INCRPOINTDOTPROD FORTRAN_NAME( INCRPOINTDOTPROD, incrpointdotprod)

inline void FORTRAN_NAME(inlineINCRPOINTDOTPROD, inlineINCRPOINTDOTPROD)(
      CHFp_FRA1(sigmadotgradu)
      ,CHFp_CONST_FRA1(gradu)
      ,CHFp_CONST_FRA1(sigma)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_INCRPOINTDOTPROD");
 FORTRAN_NAME( INCRPOINTDOTPROD ,incrpointdotprod )(
      CHFt_FRA1(sigmadotgradu)
      ,CHFt_CONST_FRA1(gradu)
      ,CHFt_CONST_FRA1(sigma)
      ,CHFt_BOX(box) );
}
#endif  // GUARDINCRPOINTDOTPROD 

#ifndef GUARDINCRCCSIGMALAMBDA 
#define GUARDINCRCCSIGMALAMBDA 
// Prototype for Fortran procedure incrccsigmalambda ...
//
void FORTRAN_NAME( INCRCCSIGMALAMBDA ,incrccsigmalambda )(
      CHFp_FRA1(sigma)
      ,CHFp_CONST_FRA(gradu)
      ,CHFp_CONST_FRA1(lambda)
      ,CHFp_BOX(grid)
      ,CHFp_CONST_INT(diagcomp) );

#define FORT_INCRCCSIGMALAMBDA FORTRAN_NAME( inlineINCRCCSIGMALAMBDA, inlineINCRCCSIGMALAMBDA)
#define FORTNT_INCRCCSIGMALAMBDA FORTRAN_NAME( INCRCCSIGMALAMBDA, incrccsigmalambda)

inline void FORTRAN_NAME(inlineINCRCCSIGMALAMBDA, inlineINCRCCSIGMALAMBDA)(
      CHFp_FRA1(sigma)
      ,CHFp_CONST_FRA(gradu)
      ,CHFp_CONST_FRA1(lambda)
      ,CHFp_BOX(grid)
      ,CHFp_CONST_INT(diagcomp) )
{
 CH_TIMELEAF("FORT_INCRCCSIGMALAMBDA");
 FORTRAN_NAME( INCRCCSIGMALAMBDA ,incrccsigmalambda )(
      CHFt_FRA1(sigma)
      ,CHFt_CONST_FRA(gradu)
      ,CHFt_CONST_FRA1(lambda)
      ,CHFt_BOX(grid)
      ,CHFt_CONST_INT(diagcomp) );
}
#endif  // GUARDINCRCCSIGMALAMBDA 

#ifndef GUARDINCRCCSIGMAETA 
#define GUARDINCRCCSIGMAETA 
// Prototype for Fortran procedure incrccsigmaeta ...
//
void FORTRAN_NAME( INCRCCSIGMAETA ,incrccsigmaeta )(
      CHFp_FRA1(sigma)
      ,CHFp_CONST_FRA(gradu)
      ,CHFp_CONST_FRA1(eta)
      ,CHFp_BOX(grid)
      ,CHFp_CONST_INT(gradcomp)
      ,CHFp_CONST_INT(trancomp) );

#define FORT_INCRCCSIGMAETA FORTRAN_NAME( inlineINCRCCSIGMAETA, inlineINCRCCSIGMAETA)
#define FORTNT_INCRCCSIGMAETA FORTRAN_NAME( INCRCCSIGMAETA, incrccsigmaeta)

inline void FORTRAN_NAME(inlineINCRCCSIGMAETA, inlineINCRCCSIGMAETA)(
      CHFp_FRA1(sigma)
      ,CHFp_CONST_FRA(gradu)
      ,CHFp_CONST_FRA1(eta)
      ,CHFp_BOX(grid)
      ,CHFp_CONST_INT(gradcomp)
      ,CHFp_CONST_INT(trancomp) )
{
 CH_TIMELEAF("FORT_INCRCCSIGMAETA");
 FORTRAN_NAME( INCRCCSIGMAETA ,incrccsigmaeta )(
      CHFt_FRA1(sigma)
      ,CHFt_CONST_FRA(gradu)
      ,CHFt_CONST_FRA1(eta)
      ,CHFt_BOX(grid)
      ,CHFt_CONST_INT(gradcomp)
      ,CHFt_CONST_INT(trancomp) );
}
#endif  // GUARDINCRCCSIGMAETA 

#ifndef GUARDINCRFACETOCELLAVERAGE 
#define GUARDINCRFACETOCELLAVERAGE 
// Prototype for Fortran procedure incrfacetocellaverage ...
//
void FORTRAN_NAME( INCRFACETOCELLAVERAGE ,incrfacetocellaverage )(
      CHFp_FRA1(cellcoef)
      ,CHFp_CONST_FRA1(facecoef)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(grid)
      ,CHFp_REAL(nfacesper) );

#define FORT_INCRFACETOCELLAVERAGE FORTRAN_NAME( inlineINCRFACETOCELLAVERAGE, inlineINCRFACETOCELLAVERAGE)
#define FORTNT_INCRFACETOCELLAVERAGE FORTRAN_NAME( INCRFACETOCELLAVERAGE, incrfacetocellaverage)

inline void FORTRAN_NAME(inlineINCRFACETOCELLAVERAGE, inlineINCRFACETOCELLAVERAGE)(
      CHFp_FRA1(cellcoef)
      ,CHFp_CONST_FRA1(facecoef)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(grid)
      ,CHFp_REAL(nfacesper) )
{
 CH_TIMELEAF("FORT_INCRFACETOCELLAVERAGE");
 FORTRAN_NAME( INCRFACETOCELLAVERAGE ,incrfacetocellaverage )(
      CHFt_FRA1(cellcoef)
      ,CHFt_CONST_FRA1(facecoef)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_BOX(grid)
      ,CHFt_REAL(nfacesper) );
}
#endif  // GUARDINCRFACETOCELLAVERAGE 

}

#endif
