#ifndef _EBAMRPOISSONOPF_F_H_
#define _EBAMRPOISSONOPF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDEBAMRPDOTPROD 
#define GUARDEBAMRPDOTPROD 
// Prototype for Fortran procedure ebamrpdotprod ...
//
void FORTRAN_NAME( EBAMRPDOTPROD ,ebamrpdotprod )(
      CHFp_REAL(sum)
      ,CHFp_CONST_FRA1(aone)
      ,CHFp_CONST_FRA1(atwo)
      ,CHFp_BOX(box) );

#define FORT_EBAMRPDOTPROD FORTRAN_NAME( inlineEBAMRPDOTPROD, inlineEBAMRPDOTPROD)
#define FORTNT_EBAMRPDOTPROD FORTRAN_NAME( EBAMRPDOTPROD, ebamrpdotprod)

inline void FORTRAN_NAME(inlineEBAMRPDOTPROD, inlineEBAMRPDOTPROD)(
      CHFp_REAL(sum)
      ,CHFp_CONST_FRA1(aone)
      ,CHFp_CONST_FRA1(atwo)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_EBAMRPDOTPROD");
 FORTRAN_NAME( EBAMRPDOTPROD ,ebamrpdotprod )(
      CHFt_REAL(sum)
      ,CHFt_CONST_FRA1(aone)
      ,CHFt_CONST_FRA1(atwo)
      ,CHFt_BOX(box) );
}
#endif  // GUARDEBAMRPDOTPROD 

#ifndef GUARDGETINVDIAGRHS 
#define GUARDGETINVDIAGRHS 
// Prototype for Fortran procedure getinvdiagrhs ...
//
void FORTRAN_NAME( GETINVDIAGRHS ,getinvdiagrhs )(
      CHFp_FRA(lhs)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_REAL(scale)
      ,CHFp_CONST_INT(ncomp)
      ,CHFp_BOX(box) );

#define FORT_GETINVDIAGRHS FORTRAN_NAME( inlineGETINVDIAGRHS, inlineGETINVDIAGRHS)
#define FORTNT_GETINVDIAGRHS FORTRAN_NAME( GETINVDIAGRHS, getinvdiagrhs)

inline void FORTRAN_NAME(inlineGETINVDIAGRHS, inlineGETINVDIAGRHS)(
      CHFp_FRA(lhs)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_REAL(scale)
      ,CHFp_CONST_INT(ncomp)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_GETINVDIAGRHS");
 FORTRAN_NAME( GETINVDIAGRHS ,getinvdiagrhs )(
      CHFt_FRA(lhs)
      ,CHFt_CONST_FRA(rhs)
      ,CHFt_REAL(scale)
      ,CHFt_CONST_INT(ncomp)
      ,CHFt_BOX(box) );
}
#endif  // GUARDGETINVDIAGRHS 

#ifndef GUARDMAXNORM 
#define GUARDMAXNORM 
// Prototype for Fortran procedure maxnorm ...
//
void FORTRAN_NAME( MAXNORM ,maxnorm )(
      CHFp_REAL(m)
      ,CHFp_CONST_FRA1(vel)
      ,CHFp_BOX(box) );

#define FORT_MAXNORM FORTRAN_NAME( inlineMAXNORM, inlineMAXNORM)
#define FORTNT_MAXNORM FORTRAN_NAME( MAXNORM, maxnorm)

inline void FORTRAN_NAME(inlineMAXNORM, inlineMAXNORM)(
      CHFp_REAL(m)
      ,CHFp_CONST_FRA1(vel)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_MAXNORM");
 FORTRAN_NAME( MAXNORM ,maxnorm )(
      CHFt_REAL(m)
      ,CHFt_CONST_FRA1(vel)
      ,CHFt_BOX(box) );
}
#endif  // GUARDMAXNORM 

#ifndef GUARDMAXNORMMASK 
#define GUARDMAXNORMMASK 
// Prototype for Fortran procedure maxnormmask ...
//
void FORTRAN_NAME( MAXNORMMASK ,maxnormmask )(
      CHFp_REAL(m)
      ,CHFp_CONST_FRA1(vel)
      ,CHFp_BOX(box)
      ,CHFp_CONST_FBA1(mask) );

#define FORT_MAXNORMMASK FORTRAN_NAME( inlineMAXNORMMASK, inlineMAXNORMMASK)
#define FORTNT_MAXNORMMASK FORTRAN_NAME( MAXNORMMASK, maxnormmask)

inline void FORTRAN_NAME(inlineMAXNORMMASK, inlineMAXNORMMASK)(
      CHFp_REAL(m)
      ,CHFp_CONST_FRA1(vel)
      ,CHFp_BOX(box)
      ,CHFp_CONST_FBA1(mask) )
{
 CH_TIMELEAF("FORT_MAXNORMMASK");
 FORTRAN_NAME( MAXNORMMASK ,maxnormmask )(
      CHFt_REAL(m)
      ,CHFt_CONST_FRA1(vel)
      ,CHFt_BOX(box)
      ,CHFt_CONST_FBA1(mask) );
}
#endif  // GUARDMAXNORMMASK 

#ifndef GUARDAMRPZEROSUB 
#define GUARDAMRPZEROSUB 
// Prototype for Fortran procedure amrpzerosub ...
//
void FORTRAN_NAME( AMRPZEROSUB ,amrpzerosub )(
      CHFp_FRA(phi)
      ,CHFp_BOX(overlaybox)
      ,CHFp_INT(ncomp) );

#define FORT_AMRPZEROSUB FORTRAN_NAME( inlineAMRPZEROSUB, inlineAMRPZEROSUB)
#define FORTNT_AMRPZEROSUB FORTRAN_NAME( AMRPZEROSUB, amrpzerosub)

inline void FORTRAN_NAME(inlineAMRPZEROSUB, inlineAMRPZEROSUB)(
      CHFp_FRA(phi)
      ,CHFp_BOX(overlaybox)
      ,CHFp_INT(ncomp) )
{
 CH_TIMELEAF("FORT_AMRPZEROSUB");
 FORTRAN_NAME( AMRPZEROSUB ,amrpzerosub )(
      CHFt_FRA(phi)
      ,CHFt_BOX(overlaybox)
      ,CHFt_INT(ncomp) );
}
#endif  // GUARDAMRPZEROSUB 

#ifndef GUARDDOALLREGULARMULTICOLOR 
#define GUARDDOALLREGULARMULTICOLOR 
// Prototype for Fortran procedure doallregularmulticolor ...
//
void FORTRAN_NAME( DOALLREGULARMULTICOLOR ,doallregularmulticolor )(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_CONST_REAL(weight)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_BOX(coloredbox) );

#define FORT_DOALLREGULARMULTICOLOR FORTRAN_NAME( inlineDOALLREGULARMULTICOLOR, inlineDOALLREGULARMULTICOLOR)
#define FORTNT_DOALLREGULARMULTICOLOR FORTRAN_NAME( DOALLREGULARMULTICOLOR, doallregularmulticolor)

inline void FORTRAN_NAME(inlineDOALLREGULARMULTICOLOR, inlineDOALLREGULARMULTICOLOR)(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_CONST_REAL(weight)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_BOX(coloredbox) )
{
 CH_TIMELEAF("FORT_DOALLREGULARMULTICOLOR");
 FORTRAN_NAME( DOALLREGULARMULTICOLOR ,doallregularmulticolor )(
      CHFt_FRA1(phi)
      ,CHFt_CONST_FRA1(rhs)
      ,CHFt_CONST_REAL(weight)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_BOX(coloredbox) );
}
#endif  // GUARDDOALLREGULARMULTICOLOR 

#ifndef GUARDDOALLREGULARUPDATE 
#define GUARDDOALLREGULARUPDATE 
// Prototype for Fortran procedure doallregularupdate ...
//
void FORTRAN_NAME( DOALLREGULARUPDATE ,doallregularupdate )(
      CHFp_FRA1(phinew)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_CONST_REAL(weight)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_BOX(coloredbox) );

#define FORT_DOALLREGULARUPDATE FORTRAN_NAME( inlineDOALLREGULARUPDATE, inlineDOALLREGULARUPDATE)
#define FORTNT_DOALLREGULARUPDATE FORTRAN_NAME( DOALLREGULARUPDATE, doallregularupdate)

inline void FORTRAN_NAME(inlineDOALLREGULARUPDATE, inlineDOALLREGULARUPDATE)(
      CHFp_FRA1(phinew)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_CONST_REAL(weight)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_BOX(coloredbox) )
{
 CH_TIMELEAF("FORT_DOALLREGULARUPDATE");
 FORTRAN_NAME( DOALLREGULARUPDATE ,doallregularupdate )(
      CHFt_FRA1(phinew)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_FRA1(rhs)
      ,CHFt_CONST_REAL(weight)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_BOX(coloredbox) );
}
#endif  // GUARDDOALLREGULARUPDATE 

#ifndef GUARDDOALLREGULARGSRB 
#define GUARDDOALLREGULARGSRB 
// Prototype for Fortran procedure doallregulargsrb ...
//
void FORTRAN_NAME( DOALLREGULARGSRB ,doallregulargsrb )(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_CONST_REAL(weight)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_BOX(region)
      ,CHFp_CONST_INT(redBlack) );

#define FORT_DOALLREGULARGSRB FORTRAN_NAME( inlineDOALLREGULARGSRB, inlineDOALLREGULARGSRB)
#define FORTNT_DOALLREGULARGSRB FORTRAN_NAME( DOALLREGULARGSRB, doallregulargsrb)

inline void FORTRAN_NAME(inlineDOALLREGULARGSRB, inlineDOALLREGULARGSRB)(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_CONST_REAL(weight)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_BOX(region)
      ,CHFp_CONST_INT(redBlack) )
{
 CH_TIMELEAF("FORT_DOALLREGULARGSRB");
 FORTRAN_NAME( DOALLREGULARGSRB ,doallregulargsrb )(
      CHFt_FRA1(phi)
      ,CHFt_CONST_FRA1(rhs)
      ,CHFt_CONST_REAL(weight)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_BOX(region)
      ,CHFt_CONST_INT(redBlack) );
}
#endif  // GUARDDOALLREGULARGSRB 

#ifndef GUARDSLOWGSRBEBAMRPO 
#define GUARDSLOWGSRBEBAMRPO 
// Prototype for Fortran procedure slowgsrbebamrpo ...
//
void FORTRAN_NAME( SLOWGSRBEBAMRPO ,slowgsrbebamrpo )(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(lph)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_CONST_FRA1(lam)
      ,CHFp_BOX(coloredbox) );

#define FORT_SLOWGSRBEBAMRPO FORTRAN_NAME( inlineSLOWGSRBEBAMRPO, inlineSLOWGSRBEBAMRPO)
#define FORTNT_SLOWGSRBEBAMRPO FORTRAN_NAME( SLOWGSRBEBAMRPO, slowgsrbebamrpo)

inline void FORTRAN_NAME(inlineSLOWGSRBEBAMRPO, inlineSLOWGSRBEBAMRPO)(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(lph)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_CONST_FRA1(lam)
      ,CHFp_BOX(coloredbox) )
{
 CH_TIMELEAF("FORT_SLOWGSRBEBAMRPO");
 FORTRAN_NAME( SLOWGSRBEBAMRPO ,slowgsrbebamrpo )(
      CHFt_FRA1(phi)
      ,CHFt_CONST_FRA1(lph)
      ,CHFt_CONST_FRA1(rhs)
      ,CHFt_CONST_FRA1(lam)
      ,CHFt_BOX(coloredbox) );
}
#endif  // GUARDSLOWGSRBEBAMRPO 

#ifndef GUARDDOALLREGULARJACOBI 
#define GUARDDOALLREGULARJACOBI 
// Prototype for Fortran procedure doallregularjacobi ...
//
void FORTRAN_NAME( DOALLREGULARJACOBI ,doallregularjacobi )(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_CONST_REAL(weight)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_BOX(box) );

#define FORT_DOALLREGULARJACOBI FORTRAN_NAME( inlineDOALLREGULARJACOBI, inlineDOALLREGULARJACOBI)
#define FORTNT_DOALLREGULARJACOBI FORTRAN_NAME( DOALLREGULARJACOBI, doallregularjacobi)

inline void FORTRAN_NAME(inlineDOALLREGULARJACOBI, inlineDOALLREGULARJACOBI)(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_CONST_REAL(weight)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_DOALLREGULARJACOBI");
 FORTRAN_NAME( DOALLREGULARJACOBI ,doallregularjacobi )(
      CHFt_FRA1(phi)
      ,CHFt_CONST_FRA1(rhs)
      ,CHFt_CONST_REAL(weight)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDDOALLREGULARJACOBI 

#ifndef GUARDUNDOREGULARGS 
#define GUARDUNDOREGULARGS 
// Prototype for Fortran procedure undoregulargs ...
//
void FORTRAN_NAME( UNDOREGULARGS ,undoregulargs )(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_CONST_REAL(weight)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_INTVECT(iv) );

#define FORT_UNDOREGULARGS FORTRAN_NAME( inlineUNDOREGULARGS, inlineUNDOREGULARGS)
#define FORTNT_UNDOREGULARGS FORTRAN_NAME( UNDOREGULARGS, undoregulargs)

inline void FORTRAN_NAME(inlineUNDOREGULARGS, inlineUNDOREGULARGS)(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_CONST_REAL(weight)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_INTVECT(iv) )
{
 CH_TIMELEAF("FORT_UNDOREGULARGS");
 FORTRAN_NAME( UNDOREGULARGS ,undoregulargs )(
      CHFt_FRA1(phi)
      ,CHFt_CONST_FRA1(rhs)
      ,CHFt_CONST_REAL(weight)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_INTVECT(iv) );
}
#endif  // GUARDUNDOREGULARGS 

#ifndef GUARDREGAPPLYDOMAINFLUX_INPLACE 
#define GUARDREGAPPLYDOMAINFLUX_INPLACE 
// Prototype for Fortran procedure regapplydomainflux_inplace ...
//
void FORTRAN_NAME( REGAPPLYDOMAINFLUX_INPLACE ,regapplydomainflux_inplace )(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(faceflux)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(side)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_REGAPPLYDOMAINFLUX_INPLACE FORTRAN_NAME( inlineREGAPPLYDOMAINFLUX_INPLACE, inlineREGAPPLYDOMAINFLUX_INPLACE)
#define FORTNT_REGAPPLYDOMAINFLUX_INPLACE FORTRAN_NAME( REGAPPLYDOMAINFLUX_INPLACE, regapplydomainflux_inplace)

inline void FORTRAN_NAME(inlineREGAPPLYDOMAINFLUX_INPLACE, inlineREGAPPLYDOMAINFLUX_INPLACE)(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(faceflux)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(side)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_REGAPPLYDOMAINFLUX_INPLACE");
 FORTRAN_NAME( REGAPPLYDOMAINFLUX_INPLACE ,regapplydomainflux_inplace )(
      CHFt_FRA1(phi)
      ,CHFt_CONST_FRA1(faceflux)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(side)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDREGAPPLYDOMAINFLUX_INPLACE 

#ifndef GUARDMVOPERATORLAP 
#define GUARDMVOPERATORLAP 
// Prototype for Fortran procedure mvoperatorlap ...
//
void FORTRAN_NAME( MVOPERATORLAP ,mvoperatorlap )(
      CHFp_FRA(lph)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_INT(ncomps)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_BOX(box) );

#define FORT_MVOPERATORLAP FORTRAN_NAME( inlineMVOPERATORLAP, inlineMVOPERATORLAP)
#define FORTNT_MVOPERATORLAP FORTRAN_NAME( MVOPERATORLAP, mvoperatorlap)

inline void FORTRAN_NAME(inlineMVOPERATORLAP, inlineMVOPERATORLAP)(
      CHFp_FRA(lph)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_INT(ncomps)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_MVOPERATORLAP");
 FORTRAN_NAME( MVOPERATORLAP ,mvoperatorlap )(
      CHFt_FRA(lph)
      ,CHFt_CONST_FRA(phi)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_INT(ncomps)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDMVOPERATORLAP 

#ifndef GUARDREGINCRLAPLACIAN_INPLACE 
#define GUARDREGINCRLAPLACIAN_INPLACE 
// Prototype for Fortran procedure regincrlaplacian_inplace ...
//
void FORTRAN_NAME( REGINCRLAPLACIAN_INPLACE ,regincrlaplacian_inplace )(
      CHFp_FRA1(opphidir)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_BOX(box) );

#define FORT_REGINCRLAPLACIAN_INPLACE FORTRAN_NAME( inlineREGINCRLAPLACIAN_INPLACE, inlineREGINCRLAPLACIAN_INPLACE)
#define FORTNT_REGINCRLAPLACIAN_INPLACE FORTRAN_NAME( REGINCRLAPLACIAN_INPLACE, regincrlaplacian_inplace)

inline void FORTRAN_NAME(inlineREGINCRLAPLACIAN_INPLACE, inlineREGINCRLAPLACIAN_INPLACE)(
      CHFp_FRA1(opphidir)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_REGINCRLAPLACIAN_INPLACE");
 FORTRAN_NAME( REGINCRLAPLACIAN_INPLACE ,regincrlaplacian_inplace )(
      CHFt_FRA1(opphidir)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDREGINCRLAPLACIAN_INPLACE 

#ifndef GUARDREGGET1DLAPLACIAN 
#define GUARDREGGET1DLAPLACIAN 
// Prototype for Fortran procedure regget1dlaplacian ...
//
void FORTRAN_NAME( REGGET1DLAPLACIAN ,regget1dlaplacian )(
      CHFp_FRA1(opphidir)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(lobox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_BOX(hibox)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_BOX(centerbox) );

#define FORT_REGGET1DLAPLACIAN FORTRAN_NAME( inlineREGGET1DLAPLACIAN, inlineREGGET1DLAPLACIAN)
#define FORTNT_REGGET1DLAPLACIAN FORTRAN_NAME( REGGET1DLAPLACIAN, regget1dlaplacian)

inline void FORTRAN_NAME(inlineREGGET1DLAPLACIAN, inlineREGGET1DLAPLACIAN)(
      CHFp_FRA1(opphidir)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(lobox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_BOX(hibox)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_BOX(centerbox) )
{
 CH_TIMELEAF("FORT_REGGET1DLAPLACIAN");
 FORTRAN_NAME( REGGET1DLAPLACIAN ,regget1dlaplacian )(
      CHFt_FRA1(opphidir)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(lobox)
      ,CHFt_CONST_INT(haslo)
      ,CHFt_BOX(hibox)
      ,CHFt_CONST_INT(hashi)
      ,CHFt_BOX(centerbox) );
}
#endif  // GUARDREGGET1DLAPLACIAN 

#ifndef GUARDREGAPPLYDOMAINFLUX 
#define GUARDREGAPPLYDOMAINFLUX 
// Prototype for Fortran procedure regapplydomainflux ...
//
void FORTRAN_NAME( REGAPPLYDOMAINFLUX ,regapplydomainflux )(
      CHFp_FRA1(opphidir)
      ,CHFp_CONST_FRA1(faceflux)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(side)
      ,CHFp_BOX(lobox) );

#define FORT_REGAPPLYDOMAINFLUX FORTRAN_NAME( inlineREGAPPLYDOMAINFLUX, inlineREGAPPLYDOMAINFLUX)
#define FORTNT_REGAPPLYDOMAINFLUX FORTRAN_NAME( REGAPPLYDOMAINFLUX, regapplydomainflux)

inline void FORTRAN_NAME(inlineREGAPPLYDOMAINFLUX, inlineREGAPPLYDOMAINFLUX)(
      CHFp_FRA1(opphidir)
      ,CHFp_CONST_FRA1(faceflux)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(side)
      ,CHFp_BOX(lobox) )
{
 CH_TIMELEAF("FORT_REGAPPLYDOMAINFLUX");
 FORTRAN_NAME( REGAPPLYDOMAINFLUX ,regapplydomainflux )(
      CHFt_FRA1(opphidir)
      ,CHFt_CONST_FRA1(faceflux)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(side)
      ,CHFt_BOX(lobox) );
}
#endif  // GUARDREGAPPLYDOMAINFLUX 

#ifndef GUARDREGSUMLAPLACIAN 
#define GUARDREGSUMLAPLACIAN 
// Prototype for Fortran procedure regsumlaplacian ...
//
void FORTRAN_NAME( REGSUMLAPLACIAN ,regsumlaplacian )(
      CHFp_FRA1(opphi)
      ,CHFp_CONST_FRA1(opphidir)
      ,CHFp_BOX(opphibox) );

#define FORT_REGSUMLAPLACIAN FORTRAN_NAME( inlineREGSUMLAPLACIAN, inlineREGSUMLAPLACIAN)
#define FORTNT_REGSUMLAPLACIAN FORTRAN_NAME( REGSUMLAPLACIAN, regsumlaplacian)

inline void FORTRAN_NAME(inlineREGSUMLAPLACIAN, inlineREGSUMLAPLACIAN)(
      CHFp_FRA1(opphi)
      ,CHFp_CONST_FRA1(opphidir)
      ,CHFp_BOX(opphibox) )
{
 CH_TIMELEAF("FORT_REGSUMLAPLACIAN");
 FORTRAN_NAME( REGSUMLAPLACIAN ,regsumlaplacian )(
      CHFt_FRA1(opphi)
      ,CHFt_CONST_FRA1(opphidir)
      ,CHFt_BOX(opphibox) );
}
#endif  // GUARDREGSUMLAPLACIAN 

#ifndef GUARDREGMULTICOLORGS 
#define GUARDREGMULTICOLORGS 
// Prototype for Fortran procedure regmulticolorgs ...
//
void FORTRAN_NAME( REGMULTICOLORGS ,regmulticolorgs )(
      CHFp_FRA1(newphi)
      ,CHFp_CONST_REAL(weight)
      ,CHFp_CONST_FRA1(resid)
      ,CHFp_BOX(coloredbox) );

#define FORT_REGMULTICOLORGS FORTRAN_NAME( inlineREGMULTICOLORGS, inlineREGMULTICOLORGS)
#define FORTNT_REGMULTICOLORGS FORTRAN_NAME( REGMULTICOLORGS, regmulticolorgs)

inline void FORTRAN_NAME(inlineREGMULTICOLORGS, inlineREGMULTICOLORGS)(
      CHFp_FRA1(newphi)
      ,CHFp_CONST_REAL(weight)
      ,CHFp_CONST_FRA1(resid)
      ,CHFp_BOX(coloredbox) )
{
 CH_TIMELEAF("FORT_REGMULTICOLORGS");
 FORTRAN_NAME( REGMULTICOLORGS ,regmulticolorgs )(
      CHFt_FRA1(newphi)
      ,CHFt_CONST_REAL(weight)
      ,CHFt_CONST_FRA1(resid)
      ,CHFt_BOX(coloredbox) );
}
#endif  // GUARDREGMULTICOLORGS 

#ifndef GUARDREGGSRB 
#define GUARDREGGSRB 
// Prototype for Fortran procedure reggsrb ...
//
void FORTRAN_NAME( REGGSRB ,reggsrb )(
      CHFp_FRA1(newphi)
      ,CHFp_CONST_FRA1(resid)
      ,CHFp_CONST_REAL(weight)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(color) );

#define FORT_REGGSRB FORTRAN_NAME( inlineREGGSRB, inlineREGGSRB)
#define FORTNT_REGGSRB FORTRAN_NAME( REGGSRB, reggsrb)

inline void FORTRAN_NAME(inlineREGGSRB, inlineREGGSRB)(
      CHFp_FRA1(newphi)
      ,CHFp_CONST_FRA1(resid)
      ,CHFp_CONST_REAL(weight)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(color) )
{
 CH_TIMELEAF("FORT_REGGSRB");
 FORTRAN_NAME( REGGSRB ,reggsrb )(
      CHFt_FRA1(newphi)
      ,CHFt_CONST_FRA1(resid)
      ,CHFt_CONST_REAL(weight)
      ,CHFt_BOX(box)
      ,CHFt_CONST_INT(color) );
}
#endif  // GUARDREGGSRB 

#ifndef GUARDREGGETFLUX 
#define GUARDREGGETFLUX 
// Prototype for Fortran procedure reggetflux ...
//
void FORTRAN_NAME( REGGETFLUX ,reggetflux )(
      CHFp_FRA(flux)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_BOX(opphibox)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(ncomp) );

#define FORT_REGGETFLUX FORTRAN_NAME( inlineREGGETFLUX, inlineREGGETFLUX)
#define FORTNT_REGGETFLUX FORTRAN_NAME( REGGETFLUX, reggetflux)

inline void FORTRAN_NAME(inlineREGGETFLUX, inlineREGGETFLUX)(
      CHFp_FRA(flux)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_BOX(opphibox)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(ncomp) )
{
 CH_TIMELEAF("FORT_REGGETFLUX");
 FORTRAN_NAME( REGGETFLUX ,reggetflux )(
      CHFt_FRA(flux)
      ,CHFt_CONST_FRA(phi)
      ,CHFt_BOX(opphibox)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(ncomp) );
}
#endif  // GUARDREGGETFLUX 

}

#endif
