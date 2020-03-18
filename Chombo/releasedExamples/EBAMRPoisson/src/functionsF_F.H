#ifndef _FUNCTIONSF_F_H_
#define _FUNCTIONSF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDGETETARESIST 
#define GUARDGETETARESIST 
// Prototype for Fortran procedure getetaresist ...
//
void FORTRAN_NAME( GETETARESIST ,getetaresist )(
      CHFp_FRA1(eta)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_BOX(box)
      ,CHFp_INT(idir)
      ,CHFp_REAL(etaval)
      ,CHFp_REAL(eps)
      ,CHFp_INT(whicheta) );

#define FORT_GETETARESIST FORTRAN_NAME( inlineGETETARESIST, inlineGETETARESIST)
#define FORTNT_GETETARESIST FORTRAN_NAME( GETETARESIST, getetaresist)

inline void FORTRAN_NAME(inlineGETETARESIST, inlineGETETARESIST)(
      CHFp_FRA1(eta)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_BOX(box)
      ,CHFp_INT(idir)
      ,CHFp_REAL(etaval)
      ,CHFp_REAL(eps)
      ,CHFp_INT(whicheta) )
{
 CH_TIMELEAF("FORT_GETETARESIST");
 FORTRAN_NAME( GETETARESIST ,getetaresist )(
      CHFt_FRA1(eta)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(problo)
      ,CHFt_BOX(box)
      ,CHFt_INT(idir)
      ,CHFt_REAL(etaval)
      ,CHFt_REAL(eps)
      ,CHFt_INT(whicheta) );
}
#endif  // GUARDGETETARESIST 

#ifndef GUARDGETETAPOINTRESIST 
#define GUARDGETETAPOINTRESIST 
// Prototype for Fortran procedure getetapointresist ...
//
void FORTRAN_NAME( GETETAPOINTRESIST ,getetapointresist )(
      CHFp_REAL(eta)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(xval)
      ,CHFp_REAL(etaval)
      ,CHFp_REAL(eps)
      ,CHFp_INT(whicheta) );

#define FORT_GETETAPOINTRESIST FORTRAN_NAME( inlineGETETAPOINTRESIST, inlineGETETAPOINTRESIST)
#define FORTNT_GETETAPOINTRESIST FORTRAN_NAME( GETETAPOINTRESIST, getetapointresist)

inline void FORTRAN_NAME(inlineGETETAPOINTRESIST, inlineGETETAPOINTRESIST)(
      CHFp_REAL(eta)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(xval)
      ,CHFp_REAL(etaval)
      ,CHFp_REAL(eps)
      ,CHFp_INT(whicheta) )
{
 CH_TIMELEAF("FORT_GETETAPOINTRESIST");
 FORTRAN_NAME( GETETAPOINTRESIST ,getetapointresist )(
      CHFt_REAL(eta)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(xval)
      ,CHFt_REAL(etaval)
      ,CHFt_REAL(eps)
      ,CHFt_INT(whicheta) );
}
#endif  // GUARDGETETAPOINTRESIST 

#ifndef GUARDGETKLBRESIST 
#define GUARDGETKLBRESIST 
// Prototype for Fortran procedure getklbresist ...
//
void FORTRAN_NAME( GETKLBRESIST ,getklbresist )(
      CHFp_FRA1(klb)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_BOX(box)
      ,CHFp_INT(icomp)
      ,CHFp_REAL(eps)
      ,CHFp_INT(whichmag)
      ,CHFp_INT(whicheta) );

#define FORT_GETKLBRESIST FORTRAN_NAME( inlineGETKLBRESIST, inlineGETKLBRESIST)
#define FORTNT_GETKLBRESIST FORTRAN_NAME( GETKLBRESIST, getklbresist)

inline void FORTRAN_NAME(inlineGETKLBRESIST, inlineGETKLBRESIST)(
      CHFp_FRA1(klb)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_BOX(box)
      ,CHFp_INT(icomp)
      ,CHFp_REAL(eps)
      ,CHFp_INT(whichmag)
      ,CHFp_INT(whicheta) )
{
 CH_TIMELEAF("FORT_GETKLBRESIST");
 FORTRAN_NAME( GETKLBRESIST ,getklbresist )(
      CHFt_FRA1(klb)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(problo)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_BOX(box)
      ,CHFt_INT(icomp)
      ,CHFt_REAL(eps)
      ,CHFt_INT(whichmag)
      ,CHFt_INT(whicheta) );
}
#endif  // GUARDGETKLBRESIST 

#ifndef GUARDGETMAGRESIST 
#define GUARDGETMAGRESIST 
// Prototype for Fortran procedure getmagresist ...
//
void FORTRAN_NAME( GETMAGRESIST ,getmagresist )(
      CHFp_FRA1(mag)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_BOX(box)
      ,CHFp_INT(icomp)
      ,CHFp_INT(whichmag) );

#define FORT_GETMAGRESIST FORTRAN_NAME( inlineGETMAGRESIST, inlineGETMAGRESIST)
#define FORTNT_GETMAGRESIST FORTRAN_NAME( GETMAGRESIST, getmagresist)

inline void FORTRAN_NAME(inlineGETMAGRESIST, inlineGETMAGRESIST)(
      CHFp_FRA1(mag)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_BOX(box)
      ,CHFp_INT(icomp)
      ,CHFp_INT(whichmag) )
{
 CH_TIMELEAF("FORT_GETMAGRESIST");
 FORTRAN_NAME( GETMAGRESIST ,getmagresist )(
      CHFt_FRA1(mag)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(problo)
      ,CHFt_BOX(box)
      ,CHFt_INT(icomp)
      ,CHFt_INT(whichmag) );
}
#endif  // GUARDGETMAGRESIST 

#ifndef GUARDGETMAGPOINTRESIST 
#define GUARDGETMAGPOINTRESIST 
// Prototype for Fortran procedure getmagpointresist ...
//
void FORTRAN_NAME( GETMAGPOINTRESIST ,getmagpointresist )(
      CHFp_REAL(mag)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(xval)
      ,CHFp_CONST_INT(icomp)
      ,CHFp_CONST_INT(whichmag) );

#define FORT_GETMAGPOINTRESIST FORTRAN_NAME( inlineGETMAGPOINTRESIST, inlineGETMAGPOINTRESIST)
#define FORTNT_GETMAGPOINTRESIST FORTRAN_NAME( GETMAGPOINTRESIST, getmagpointresist)

inline void FORTRAN_NAME(inlineGETMAGPOINTRESIST, inlineGETMAGPOINTRESIST)(
      CHFp_REAL(mag)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(xval)
      ,CHFp_CONST_INT(icomp)
      ,CHFp_CONST_INT(whichmag) )
{
 CH_TIMELEAF("FORT_GETMAGPOINTRESIST");
 FORTRAN_NAME( GETMAGPOINTRESIST ,getmagpointresist )(
      CHFt_REAL(mag)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(xval)
      ,CHFt_CONST_INT(icomp)
      ,CHFt_CONST_INT(whichmag) );
}
#endif  // GUARDGETMAGPOINTRESIST 

#ifndef GUARDGETDVDXPOINTRESIST 
#define GUARDGETDVDXPOINTRESIST 
// Prototype for Fortran procedure getdvdxpointresist ...
//
void FORTRAN_NAME( GETDVDXPOINTRESIST ,getdvdxpointresist )(
      CHFp_REAL(mag)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(xval)
      ,CHFp_CONST_INT(icomp)
      ,CHFp_CONST_INT(ideriv)
      ,CHFp_CONST_INT(whichmag) );

#define FORT_GETDVDXPOINTRESIST FORTRAN_NAME( inlineGETDVDXPOINTRESIST, inlineGETDVDXPOINTRESIST)
#define FORTNT_GETDVDXPOINTRESIST FORTRAN_NAME( GETDVDXPOINTRESIST, getdvdxpointresist)

inline void FORTRAN_NAME(inlineGETDVDXPOINTRESIST, inlineGETDVDXPOINTRESIST)(
      CHFp_REAL(mag)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(xval)
      ,CHFp_CONST_INT(icomp)
      ,CHFp_CONST_INT(ideriv)
      ,CHFp_CONST_INT(whichmag) )
{
 CH_TIMELEAF("FORT_GETDVDXPOINTRESIST");
 FORTRAN_NAME( GETDVDXPOINTRESIST ,getdvdxpointresist )(
      CHFt_REAL(mag)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(xval)
      ,CHFt_CONST_INT(icomp)
      ,CHFt_CONST_INT(ideriv)
      ,CHFt_CONST_INT(whichmag) );
}
#endif  // GUARDGETDVDXPOINTRESIST 

#ifndef GUARDGETKLBPOINTRESIST 
#define GUARDGETKLBPOINTRESIST 
// Prototype for Fortran procedure getklbpointresist ...
//
void FORTRAN_NAME( GETKLBPOINTRESIST ,getklbpointresist )(
      CHFp_REAL(klb)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(xvec)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_INT(icomp)
      ,CHFp_REAL(eps)
      ,CHFp_INT(whichmag)
      ,CHFp_INT(whicheta) );

#define FORT_GETKLBPOINTRESIST FORTRAN_NAME( inlineGETKLBPOINTRESIST, inlineGETKLBPOINTRESIST)
#define FORTNT_GETKLBPOINTRESIST FORTRAN_NAME( GETKLBPOINTRESIST, getklbpointresist)

inline void FORTRAN_NAME(inlineGETKLBPOINTRESIST, inlineGETKLBPOINTRESIST)(
      CHFp_REAL(klb)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(xvec)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_INT(icomp)
      ,CHFp_REAL(eps)
      ,CHFp_INT(whichmag)
      ,CHFp_INT(whicheta) )
{
 CH_TIMELEAF("FORT_GETKLBPOINTRESIST");
 FORTRAN_NAME( GETKLBPOINTRESIST ,getklbpointresist )(
      CHFt_REAL(klb)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(xvec)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_INT(icomp)
      ,CHFt_REAL(eps)
      ,CHFt_INT(whichmag)
      ,CHFt_INT(whicheta) );
}
#endif  // GUARDGETKLBPOINTRESIST 

#ifndef GUARDGETKLVVISCOUS 
#define GUARDGETKLVVISCOUS 
// Prototype for Fortran procedure getklvviscous ...
//
void FORTRAN_NAME( GETKLVVISCOUS ,getklvviscous )(
      CHFp_FRA1(klb)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_BOX(box)
      ,CHFp_INT(icomp)
      ,CHFp_REAL(eps)
      ,CHFp_INT(whichvel)
      ,CHFp_INT(whicheta)
      ,CHFp_INT(whichlambda)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_REAL(lambdafactor) );

#define FORT_GETKLVVISCOUS FORTRAN_NAME( inlineGETKLVVISCOUS, inlineGETKLVVISCOUS)
#define FORTNT_GETKLVVISCOUS FORTRAN_NAME( GETKLVVISCOUS, getklvviscous)

inline void FORTRAN_NAME(inlineGETKLVVISCOUS, inlineGETKLVVISCOUS)(
      CHFp_FRA1(klb)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_BOX(box)
      ,CHFp_INT(icomp)
      ,CHFp_REAL(eps)
      ,CHFp_INT(whichvel)
      ,CHFp_INT(whicheta)
      ,CHFp_INT(whichlambda)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_REAL(lambdafactor) )
{
 CH_TIMELEAF("FORT_GETKLVVISCOUS");
 FORTRAN_NAME( GETKLVVISCOUS ,getklvviscous )(
      CHFt_FRA1(klb)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(problo)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_BOX(box)
      ,CHFt_INT(icomp)
      ,CHFt_REAL(eps)
      ,CHFt_INT(whichvel)
      ,CHFt_INT(whicheta)
      ,CHFt_INT(whichlambda)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_REAL(lambdafactor) );
}
#endif  // GUARDGETKLVVISCOUS 

#ifndef GUARDGETKLVPOINTVISCOUS 
#define GUARDGETKLVPOINTVISCOUS 
// Prototype for Fortran procedure getklvpointviscous ...
//
void FORTRAN_NAME( GETKLVPOINTVISCOUS ,getklvpointviscous )(
      CHFp_REAL(klv)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(xvec)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_INT(icomp)
      ,CHFp_REAL(eps)
      ,CHFp_INT(whichvel)
      ,CHFp_INT(whicheta)
      ,CHFp_INT(whichlambda)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_REAL(lambdafactor) );

#define FORT_GETKLVPOINTVISCOUS FORTRAN_NAME( inlineGETKLVPOINTVISCOUS, inlineGETKLVPOINTVISCOUS)
#define FORTNT_GETKLVPOINTVISCOUS FORTRAN_NAME( GETKLVPOINTVISCOUS, getklvpointviscous)

inline void FORTRAN_NAME(inlineGETKLVPOINTVISCOUS, inlineGETKLVPOINTVISCOUS)(
      CHFp_REAL(klv)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(xvec)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_INT(icomp)
      ,CHFp_REAL(eps)
      ,CHFp_INT(whichvel)
      ,CHFp_INT(whicheta)
      ,CHFp_INT(whichlambda)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_REAL(lambdafactor) )
{
 CH_TIMELEAF("FORT_GETKLVPOINTVISCOUS");
 FORTRAN_NAME( GETKLVPOINTVISCOUS ,getklvpointviscous )(
      CHFt_REAL(klv)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(xvec)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_INT(icomp)
      ,CHFt_REAL(eps)
      ,CHFt_INT(whichvel)
      ,CHFt_INT(whicheta)
      ,CHFt_INT(whichlambda)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_REAL(lambdafactor) );
}
#endif  // GUARDGETKLVPOINTVISCOUS 

#ifndef GUARDGETPHI 
#define GUARDGETPHI 
// Prototype for Fortran procedure getphi ...
//
void FORTRAN_NAME( GETPHI ,getphi )(
      CHFp_FRA1(phi)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(time)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_CONST_REALVECT(probhi)
      ,CHFp_BOX(box) );

#define FORT_GETPHI FORTRAN_NAME( inlineGETPHI, inlineGETPHI)
#define FORTNT_GETPHI FORTRAN_NAME( GETPHI, getphi)

inline void FORTRAN_NAME(inlineGETPHI, inlineGETPHI)(
      CHFp_FRA1(phi)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(time)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_CONST_REALVECT(probhi)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_GETPHI");
 FORTRAN_NAME( GETPHI ,getphi )(
      CHFt_FRA1(phi)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REAL(time)
      ,CHFt_CONST_REALVECT(problo)
      ,CHFt_CONST_REALVECT(probhi)
      ,CHFt_BOX(box) );
}
#endif  // GUARDGETPHI 

#ifndef GUARDGETPHIPOINT 
#define GUARDGETPHIPOINT 
// Prototype for Fortran procedure getphipoint ...
//
void FORTRAN_NAME( GETPHIPOINT ,getphipoint )(
      CHFp_REAL(phi)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(time) );

#define FORT_GETPHIPOINT FORTRAN_NAME( inlineGETPHIPOINT, inlineGETPHIPOINT)
#define FORTNT_GETPHIPOINT FORTRAN_NAME( GETPHIPOINT, getphipoint)

inline void FORTRAN_NAME(inlineGETPHIPOINT, inlineGETPHIPOINT)(
      CHFp_REAL(phi)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(time) )
{
 CH_TIMELEAF("FORT_GETPHIPOINT");
 FORTRAN_NAME( GETPHIPOINT ,getphipoint )(
      CHFt_REAL(phi)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(x)
      ,CHFt_CONST_REAL(time) );
}
#endif  // GUARDGETPHIPOINT 

#ifndef GUARDGETSHPHI 
#define GUARDGETSHPHI 
// Prototype for Fortran procedure getshphi ...
//
void FORTRAN_NAME( GETSHPHI ,getshphi )(
      CHFp_FRA1(phi)
      ,CHFp_CONST_REALVECT(lmp)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(time)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_CONST_REALVECT(probhi)
      ,CHFp_BOX(box) );

#define FORT_GETSHPHI FORTRAN_NAME( inlineGETSHPHI, inlineGETSHPHI)
#define FORTNT_GETSHPHI FORTRAN_NAME( GETSHPHI, getshphi)

inline void FORTRAN_NAME(inlineGETSHPHI, inlineGETSHPHI)(
      CHFp_FRA1(phi)
      ,CHFp_CONST_REALVECT(lmp)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(time)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_CONST_REALVECT(probhi)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_GETSHPHI");
 FORTRAN_NAME( GETSHPHI ,getshphi )(
      CHFt_FRA1(phi)
      ,CHFt_CONST_REALVECT(lmp)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REAL(time)
      ,CHFt_CONST_REALVECT(problo)
      ,CHFt_CONST_REALVECT(probhi)
      ,CHFt_BOX(box) );
}
#endif  // GUARDGETSHPHI 

#ifndef GUARDGETSHPHIPOINT 
#define GUARDGETSHPHIPOINT 
// Prototype for Fortran procedure getshphipoint ...
//
void FORTRAN_NAME( GETSHPHIPOINT ,getshphipoint )(
      CHFp_REAL(phi)
      ,CHFp_CONST_REALVECT(lmp)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(time) );

#define FORT_GETSHPHIPOINT FORTRAN_NAME( inlineGETSHPHIPOINT, inlineGETSHPHIPOINT)
#define FORTNT_GETSHPHIPOINT FORTRAN_NAME( GETSHPHIPOINT, getshphipoint)

inline void FORTRAN_NAME(inlineGETSHPHIPOINT, inlineGETSHPHIPOINT)(
      CHFp_REAL(phi)
      ,CHFp_CONST_REALVECT(lmp)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(time) )
{
 CH_TIMELEAF("FORT_GETSHPHIPOINT");
 FORTRAN_NAME( GETSHPHIPOINT ,getshphipoint )(
      CHFt_REAL(phi)
      ,CHFt_CONST_REALVECT(lmp)
      ,CHFt_CONST_REALVECT(x)
      ,CHFt_CONST_REAL(time) );
}
#endif  // GUARDGETSHPHIPOINT 

#ifndef GUARDGETLOFPHIRZPOLY 
#define GUARDGETLOFPHIRZPOLY 
// Prototype for Fortran procedure getlofphirzpoly ...
//
void FORTRAN_NAME( GETLOFPHIRZPOLY ,getlofphirzpoly )(
      CHFp_REAL(lofphi)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta) );

#define FORT_GETLOFPHIRZPOLY FORTRAN_NAME( inlineGETLOFPHIRZPOLY, inlineGETLOFPHIRZPOLY)
#define FORTNT_GETLOFPHIRZPOLY FORTRAN_NAME( GETLOFPHIRZPOLY, getlofphirzpoly)

inline void FORTRAN_NAME(inlineGETLOFPHIRZPOLY, inlineGETLOFPHIRZPOLY)(
      CHFp_REAL(lofphi)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta) )
{
 CH_TIMELEAF("FORT_GETLOFPHIRZPOLY");
 FORTRAN_NAME( GETLOFPHIRZPOLY ,getlofphirzpoly )(
      CHFt_REAL(lofphi)
      ,CHFt_CONST_REALVECT(x)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta) );
}
#endif  // GUARDGETLOFPHIRZPOLY 

#ifndef GUARDGETPHIRZPOLY 
#define GUARDGETPHIRZPOLY 
// Prototype for Fortran procedure getphirzpoly ...
//
void FORTRAN_NAME( GETPHIRZPOLY ,getphirzpoly )(
      CHFp_REAL(phi)
      ,CHFp_CONST_REALVECT(x) );

#define FORT_GETPHIRZPOLY FORTRAN_NAME( inlineGETPHIRZPOLY, inlineGETPHIRZPOLY)
#define FORTNT_GETPHIRZPOLY FORTRAN_NAME( GETPHIRZPOLY, getphirzpoly)

inline void FORTRAN_NAME(inlineGETPHIRZPOLY, inlineGETPHIRZPOLY)(
      CHFp_REAL(phi)
      ,CHFp_CONST_REALVECT(x) )
{
 CH_TIMELEAF("FORT_GETPHIRZPOLY");
 FORTRAN_NAME( GETPHIRZPOLY ,getphirzpoly )(
      CHFt_REAL(phi)
      ,CHFt_CONST_REALVECT(x) );
}
#endif  // GUARDGETPHIRZPOLY 

#ifndef GUARDGETGRADPHIRZPOLY 
#define GUARDGETGRADPHIRZPOLY 
// Prototype for Fortran procedure getgradphirzpoly ...
//
void FORTRAN_NAME( GETGRADPHIRZPOLY ,getgradphirzpoly )(
      CHFp_REALVECT(gradphi)
      ,CHFp_CONST_REALVECT(x) );

#define FORT_GETGRADPHIRZPOLY FORTRAN_NAME( inlineGETGRADPHIRZPOLY, inlineGETGRADPHIRZPOLY)
#define FORTNT_GETGRADPHIRZPOLY FORTRAN_NAME( GETGRADPHIRZPOLY, getgradphirzpoly)

inline void FORTRAN_NAME(inlineGETGRADPHIRZPOLY, inlineGETGRADPHIRZPOLY)(
      CHFp_REALVECT(gradphi)
      ,CHFp_CONST_REALVECT(x) )
{
 CH_TIMELEAF("FORT_GETGRADPHIRZPOLY");
 FORTRAN_NAME( GETGRADPHIRZPOLY ,getgradphirzpoly )(
      CHFt_REALVECT(gradphi)
      ,CHFt_CONST_REALVECT(x) );
}
#endif  // GUARDGETGRADPHIRZPOLY 

#ifndef GUARDGETGRADPHIPOINT 
#define GUARDGETGRADPHIPOINT 
// Prototype for Fortran procedure getgradphipoint ...
//
void FORTRAN_NAME( GETGRADPHIPOINT ,getgradphipoint )(
      CHFp_REALVECT(gradphi)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(time) );

#define FORT_GETGRADPHIPOINT FORTRAN_NAME( inlineGETGRADPHIPOINT, inlineGETGRADPHIPOINT)
#define FORTNT_GETGRADPHIPOINT FORTRAN_NAME( GETGRADPHIPOINT, getgradphipoint)

inline void FORTRAN_NAME(inlineGETGRADPHIPOINT, inlineGETGRADPHIPOINT)(
      CHFp_REALVECT(gradphi)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(time) )
{
 CH_TIMELEAF("FORT_GETGRADPHIPOINT");
 FORTRAN_NAME( GETGRADPHIPOINT ,getgradphipoint )(
      CHFt_REALVECT(gradphi)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(x)
      ,CHFt_CONST_REAL(time) );
}
#endif  // GUARDGETGRADPHIPOINT 

#ifndef GUARDGETGRADSHPHIPOINT 
#define GUARDGETGRADSHPHIPOINT 
// Prototype for Fortran procedure getgradshphipoint ...
//
void FORTRAN_NAME( GETGRADSHPHIPOINT ,getgradshphipoint )(
      CHFp_REALVECT(gradphi)
      ,CHFp_CONST_REALVECT(lmp)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(time) );

#define FORT_GETGRADSHPHIPOINT FORTRAN_NAME( inlineGETGRADSHPHIPOINT, inlineGETGRADSHPHIPOINT)
#define FORTNT_GETGRADSHPHIPOINT FORTRAN_NAME( GETGRADSHPHIPOINT, getgradshphipoint)

inline void FORTRAN_NAME(inlineGETGRADSHPHIPOINT, inlineGETGRADSHPHIPOINT)(
      CHFp_REALVECT(gradphi)
      ,CHFp_CONST_REALVECT(lmp)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(time) )
{
 CH_TIMELEAF("FORT_GETGRADSHPHIPOINT");
 FORTRAN_NAME( GETGRADSHPHIPOINT ,getgradshphipoint )(
      CHFt_REALVECT(gradphi)
      ,CHFt_CONST_REALVECT(lmp)
      ,CHFt_CONST_REALVECT(x)
      ,CHFt_CONST_REAL(time) );
}
#endif  // GUARDGETGRADSHPHIPOINT 

#ifndef GUARDGETMARSHAGRADPHIPOINT 
#define GUARDGETMARSHAGRADPHIPOINT 
// Prototype for Fortran procedure getmarshagradphipoint ...
//
void FORTRAN_NAME( GETMARSHAGRADPHIPOINT ,getmarshagradphipoint )(
      CHFp_REALVECT(gradphi)
      ,CHFp_CONST_REALVECT(x) );

#define FORT_GETMARSHAGRADPHIPOINT FORTRAN_NAME( inlineGETMARSHAGRADPHIPOINT, inlineGETMARSHAGRADPHIPOINT)
#define FORTNT_GETMARSHAGRADPHIPOINT FORTRAN_NAME( GETMARSHAGRADPHIPOINT, getmarshagradphipoint)

inline void FORTRAN_NAME(inlineGETMARSHAGRADPHIPOINT, inlineGETMARSHAGRADPHIPOINT)(
      CHFp_REALVECT(gradphi)
      ,CHFp_CONST_REALVECT(x) )
{
 CH_TIMELEAF("FORT_GETMARSHAGRADPHIPOINT");
 FORTRAN_NAME( GETMARSHAGRADPHIPOINT ,getmarshagradphipoint )(
      CHFt_REALVECT(gradphi)
      ,CHFt_CONST_REALVECT(x) );
}
#endif  // GUARDGETMARSHAGRADPHIPOINT 

#ifndef GUARDGETLOFPHI 
#define GUARDGETLOFPHI 
// Prototype for Fortran procedure getlofphi ...
//
void FORTRAN_NAME( GETLOFPHI ,getlofphi )(
      CHFp_FRA1(lofphi)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_CONST_REALVECT(probhi)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REAL(time)
      ,CHFp_BOX(box) );

#define FORT_GETLOFPHI FORTRAN_NAME( inlineGETLOFPHI, inlineGETLOFPHI)
#define FORTNT_GETLOFPHI FORTRAN_NAME( GETLOFPHI, getlofphi)

inline void FORTRAN_NAME(inlineGETLOFPHI, inlineGETLOFPHI)(
      CHFp_FRA1(lofphi)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_CONST_REALVECT(probhi)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REAL(time)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_GETLOFPHI");
 FORTRAN_NAME( GETLOFPHI ,getlofphi )(
      CHFt_FRA1(lofphi)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(problo)
      ,CHFt_CONST_REALVECT(probhi)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_REAL(time)
      ,CHFt_BOX(box) );
}
#endif  // GUARDGETLOFPHI 

#ifndef GUARDGETLOFSHPHI 
#define GUARDGETLOFSHPHI 
// Prototype for Fortran procedure getlofshphi ...
//
void FORTRAN_NAME( GETLOFSHPHI ,getlofshphi )(
      CHFp_FRA1(lofphi)
      ,CHFp_CONST_REALVECT(lmp)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_CONST_REALVECT(probhi)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REAL(time)
      ,CHFp_BOX(box) );

#define FORT_GETLOFSHPHI FORTRAN_NAME( inlineGETLOFSHPHI, inlineGETLOFSHPHI)
#define FORTNT_GETLOFSHPHI FORTRAN_NAME( GETLOFSHPHI, getlofshphi)

inline void FORTRAN_NAME(inlineGETLOFSHPHI, inlineGETLOFSHPHI)(
      CHFp_FRA1(lofphi)
      ,CHFp_CONST_REALVECT(lmp)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_CONST_REALVECT(probhi)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REAL(time)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_GETLOFSHPHI");
 FORTRAN_NAME( GETLOFSHPHI ,getlofshphi )(
      CHFt_FRA1(lofphi)
      ,CHFt_CONST_REALVECT(lmp)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(problo)
      ,CHFt_CONST_REALVECT(probhi)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_REAL(time)
      ,CHFt_BOX(box) );
}
#endif  // GUARDGETLOFSHPHI 

#ifndef GUARDGETMARSHALOFPHI 
#define GUARDGETMARSHALOFPHI 
// Prototype for Fortran procedure getmarshalofphi ...
//
void FORTRAN_NAME( GETMARSHALOFPHI ,getmarshalofphi )(
      CHFp_FRA1(lofphi)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_CONST_REALVECT(probhi)
      ,CHFp_BOX(box) );

#define FORT_GETMARSHALOFPHI FORTRAN_NAME( inlineGETMARSHALOFPHI, inlineGETMARSHALOFPHI)
#define FORTNT_GETMARSHALOFPHI FORTRAN_NAME( GETMARSHALOFPHI, getmarshalofphi)

inline void FORTRAN_NAME(inlineGETMARSHALOFPHI, inlineGETMARSHALOFPHI)(
      CHFp_FRA1(lofphi)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_CONST_REALVECT(probhi)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_GETMARSHALOFPHI");
 FORTRAN_NAME( GETMARSHALOFPHI ,getmarshalofphi )(
      CHFt_FRA1(lofphi)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(problo)
      ,CHFt_CONST_REALVECT(probhi)
      ,CHFt_BOX(box) );
}
#endif  // GUARDGETMARSHALOFPHI 

#ifndef GUARDGETMARSHALOFPHIPOINT 
#define GUARDGETMARSHALOFPHIPOINT 
// Prototype for Fortran procedure getmarshalofphipoint ...
//
void FORTRAN_NAME( GETMARSHALOFPHIPOINT ,getmarshalofphipoint )(
      CHFp_REAL(lofphi)
      ,CHFp_CONST_REALVECT(x) );

#define FORT_GETMARSHALOFPHIPOINT FORTRAN_NAME( inlineGETMARSHALOFPHIPOINT, inlineGETMARSHALOFPHIPOINT)
#define FORTNT_GETMARSHALOFPHIPOINT FORTRAN_NAME( GETMARSHALOFPHIPOINT, getmarshalofphipoint)

inline void FORTRAN_NAME(inlineGETMARSHALOFPHIPOINT, inlineGETMARSHALOFPHIPOINT)(
      CHFp_REAL(lofphi)
      ,CHFp_CONST_REALVECT(x) )
{
 CH_TIMELEAF("FORT_GETMARSHALOFPHIPOINT");
 FORTRAN_NAME( GETMARSHALOFPHIPOINT ,getmarshalofphipoint )(
      CHFt_REAL(lofphi)
      ,CHFt_CONST_REALVECT(x) );
}
#endif  // GUARDGETMARSHALOFPHIPOINT 

#ifndef GUARDGETMARSHAPHI 
#define GUARDGETMARSHAPHI 
// Prototype for Fortran procedure getmarshaphi ...
//
void FORTRAN_NAME( GETMARSHAPHI ,getmarshaphi )(
      CHFp_FRA1(phi)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_CONST_REALVECT(probhi)
      ,CHFp_BOX(box) );

#define FORT_GETMARSHAPHI FORTRAN_NAME( inlineGETMARSHAPHI, inlineGETMARSHAPHI)
#define FORTNT_GETMARSHAPHI FORTRAN_NAME( GETMARSHAPHI, getmarshaphi)

inline void FORTRAN_NAME(inlineGETMARSHAPHI, inlineGETMARSHAPHI)(
      CHFp_FRA1(phi)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_CONST_REALVECT(probhi)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_GETMARSHAPHI");
 FORTRAN_NAME( GETMARSHAPHI ,getmarshaphi )(
      CHFt_FRA1(phi)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(problo)
      ,CHFt_CONST_REALVECT(probhi)
      ,CHFt_BOX(box) );
}
#endif  // GUARDGETMARSHAPHI 

#ifndef GUARDGETMARSHAPHIPOINT 
#define GUARDGETMARSHAPHIPOINT 
// Prototype for Fortran procedure getmarshaphipoint ...
//
void FORTRAN_NAME( GETMARSHAPHIPOINT ,getmarshaphipoint )(
      CHFp_REAL(phi)
      ,CHFp_CONST_REALVECT(x) );

#define FORT_GETMARSHAPHIPOINT FORTRAN_NAME( inlineGETMARSHAPHIPOINT, inlineGETMARSHAPHIPOINT)
#define FORTNT_GETMARSHAPHIPOINT FORTRAN_NAME( GETMARSHAPHIPOINT, getmarshaphipoint)

inline void FORTRAN_NAME(inlineGETMARSHAPHIPOINT, inlineGETMARSHAPHIPOINT)(
      CHFp_REAL(phi)
      ,CHFp_CONST_REALVECT(x) )
{
 CH_TIMELEAF("FORT_GETMARSHAPHIPOINT");
 FORTRAN_NAME( GETMARSHAPHIPOINT ,getmarshaphipoint )(
      CHFt_REAL(phi)
      ,CHFt_CONST_REALVECT(x) );
}
#endif  // GUARDGETMARSHAPHIPOINT 

#ifndef GUARDGETLOFPHIPOINT 
#define GUARDGETLOFPHIPOINT 
// Prototype for Fortran procedure getlofphipoint ...
//
void FORTRAN_NAME( GETLOFPHIPOINT ,getlofphipoint )(
      CHFp_REAL(lofphi)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REAL(time) );

#define FORT_GETLOFPHIPOINT FORTRAN_NAME( inlineGETLOFPHIPOINT, inlineGETLOFPHIPOINT)
#define FORTNT_GETLOFPHIPOINT FORTRAN_NAME( GETLOFPHIPOINT, getlofphipoint)

inline void FORTRAN_NAME(inlineGETLOFPHIPOINT, inlineGETLOFPHIPOINT)(
      CHFp_REAL(lofphi)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REAL(time) )
{
 CH_TIMELEAF("FORT_GETLOFPHIPOINT");
 FORTRAN_NAME( GETLOFPHIPOINT ,getlofphipoint )(
      CHFt_REAL(lofphi)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(x)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_REAL(time) );
}
#endif  // GUARDGETLOFPHIPOINT 

#ifndef GUARDGETLOFSHPHIPOINT 
#define GUARDGETLOFSHPHIPOINT 
// Prototype for Fortran procedure getlofshphipoint ...
//
void FORTRAN_NAME( GETLOFSHPHIPOINT ,getlofshphipoint )(
      CHFp_REAL(lofphi)
      ,CHFp_CONST_REALVECT(lmp)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REAL(time) );

#define FORT_GETLOFSHPHIPOINT FORTRAN_NAME( inlineGETLOFSHPHIPOINT, inlineGETLOFSHPHIPOINT)
#define FORTNT_GETLOFSHPHIPOINT FORTRAN_NAME( GETLOFSHPHIPOINT, getlofshphipoint)

inline void FORTRAN_NAME(inlineGETLOFSHPHIPOINT, inlineGETLOFSHPHIPOINT)(
      CHFp_REAL(lofphi)
      ,CHFp_CONST_REALVECT(lmp)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REAL(time) )
{
 CH_TIMELEAF("FORT_GETLOFSHPHIPOINT");
 FORTRAN_NAME( GETLOFSHPHIPOINT ,getlofshphipoint )(
      CHFt_REAL(lofphi)
      ,CHFt_CONST_REALVECT(lmp)
      ,CHFt_CONST_REALVECT(x)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_REAL(time) );
}
#endif  // GUARDGETLOFSHPHIPOINT 

#ifndef GUARDGETDBGPHI 
#define GUARDGETDBGPHI 
// Prototype for Fortran procedure getdbgphi ...
//
void FORTRAN_NAME( GETDBGPHI ,getdbgphi )(
      CHFp_FRA1(dbgphi)
      ,CHFp_FRA1(beta)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_CONST_REALVECT(probhi)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(time)
      ,CHFp_BOX(box) );

#define FORT_GETDBGPHI FORTRAN_NAME( inlineGETDBGPHI, inlineGETDBGPHI)
#define FORTNT_GETDBGPHI FORTRAN_NAME( GETDBGPHI, getdbgphi)

inline void FORTRAN_NAME(inlineGETDBGPHI, inlineGETDBGPHI)(
      CHFp_FRA1(dbgphi)
      ,CHFp_FRA1(beta)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_CONST_REALVECT(probhi)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(time)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_GETDBGPHI");
 FORTRAN_NAME( GETDBGPHI ,getdbgphi )(
      CHFt_FRA1(dbgphi)
      ,CHFt_FRA1(beta)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(problo)
      ,CHFt_CONST_REALVECT(probhi)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(time)
      ,CHFt_BOX(box) );
}
#endif  // GUARDGETDBGPHI 

#ifndef GUARDGETDBGPHIPOINT 
#define GUARDGETDBGPHIPOINT 
// Prototype for Fortran procedure getdbgphipoint ...
//
void FORTRAN_NAME( GETDBGPHIPOINT ,getdbgphipoint )(
      CHFp_REAL(dbgphi)
      ,CHFp_REAL(beta)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(time) );

#define FORT_GETDBGPHIPOINT FORTRAN_NAME( inlineGETDBGPHIPOINT, inlineGETDBGPHIPOINT)
#define FORTNT_GETDBGPHIPOINT FORTRAN_NAME( GETDBGPHIPOINT, getdbgphipoint)

inline void FORTRAN_NAME(inlineGETDBGPHIPOINT, inlineGETDBGPHIPOINT)(
      CHFp_REAL(dbgphi)
      ,CHFp_REAL(beta)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(time) )
{
 CH_TIMELEAF("FORT_GETDBGPHIPOINT");
 FORTRAN_NAME( GETDBGPHIPOINT ,getdbgphipoint )(
      CHFt_REAL(dbgphi)
      ,CHFt_REAL(beta)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(x)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(time) );
}
#endif  // GUARDGETDBGPHIPOINT 

#ifndef GUARDGETBETAPOINT 
#define GUARDGETBETAPOINT 
// Prototype for Fortran procedure getbetapoint ...
//
void FORTRAN_NAME( GETBETAPOINT ,getbetapoint )(
      CHFp_REAL(beta)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(time) );

#define FORT_GETBETAPOINT FORTRAN_NAME( inlineGETBETAPOINT, inlineGETBETAPOINT)
#define FORTNT_GETBETAPOINT FORTRAN_NAME( GETBETAPOINT, getbetapoint)

inline void FORTRAN_NAME(inlineGETBETAPOINT, inlineGETBETAPOINT)(
      CHFp_REAL(beta)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(time) )
{
 CH_TIMELEAF("FORT_GETBETAPOINT");
 FORTRAN_NAME( GETBETAPOINT ,getbetapoint )(
      CHFt_REAL(beta)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(x)
      ,CHFt_CONST_REAL(time) );
}
#endif  // GUARDGETBETAPOINT 

#ifndef GUARDGETGRADBETAPOINT 
#define GUARDGETGRADBETAPOINT 
// Prototype for Fortran procedure getgradbetapoint ...
//
void FORTRAN_NAME( GETGRADBETAPOINT ,getgradbetapoint )(
      CHFp_REALVECT(gradbeta)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(time) );

#define FORT_GETGRADBETAPOINT FORTRAN_NAME( inlineGETGRADBETAPOINT, inlineGETGRADBETAPOINT)
#define FORTNT_GETGRADBETAPOINT FORTRAN_NAME( GETGRADBETAPOINT, getgradbetapoint)

inline void FORTRAN_NAME(inlineGETGRADBETAPOINT, inlineGETGRADBETAPOINT)(
      CHFp_REALVECT(gradbeta)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(time) )
{
 CH_TIMELEAF("FORT_GETGRADBETAPOINT");
 FORTRAN_NAME( GETGRADBETAPOINT ,getgradbetapoint )(
      CHFt_REALVECT(gradbeta)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(x)
      ,CHFt_CONST_REAL(time) );
}
#endif  // GUARDGETGRADBETAPOINT 

#ifndef GUARDGETBETAGRADPHIPOINT 
#define GUARDGETBETAGRADPHIPOINT 
// Prototype for Fortran procedure getbetagradphipoint ...
//
void FORTRAN_NAME( GETBETAGRADPHIPOINT ,getbetagradphipoint )(
      CHFp_REALVECT(gradphi)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(time) );

#define FORT_GETBETAGRADPHIPOINT FORTRAN_NAME( inlineGETBETAGRADPHIPOINT, inlineGETBETAGRADPHIPOINT)
#define FORTNT_GETBETAGRADPHIPOINT FORTRAN_NAME( GETBETAGRADPHIPOINT, getbetagradphipoint)

inline void FORTRAN_NAME(inlineGETBETAGRADPHIPOINT, inlineGETBETAGRADPHIPOINT)(
      CHFp_REALVECT(gradphi)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(time) )
{
 CH_TIMELEAF("FORT_GETBETAGRADPHIPOINT");
 FORTRAN_NAME( GETBETAGRADPHIPOINT ,getbetagradphipoint )(
      CHFt_REALVECT(gradphi)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(x)
      ,CHFt_CONST_REAL(time) );
}
#endif  // GUARDGETBETAGRADPHIPOINT 

#ifndef GUARDGETBETAGRADSHPHIPOINT 
#define GUARDGETBETAGRADSHPHIPOINT 
// Prototype for Fortran procedure getbetagradshphipoint ...
//
void FORTRAN_NAME( GETBETAGRADSHPHIPOINT ,getbetagradshphipoint )(
      CHFp_REALVECT(gradphi)
      ,CHFp_CONST_REALVECT(lmp)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(time) );

#define FORT_GETBETAGRADSHPHIPOINT FORTRAN_NAME( inlineGETBETAGRADSHPHIPOINT, inlineGETBETAGRADSHPHIPOINT)
#define FORTNT_GETBETAGRADSHPHIPOINT FORTRAN_NAME( GETBETAGRADSHPHIPOINT, getbetagradshphipoint)

inline void FORTRAN_NAME(inlineGETBETAGRADSHPHIPOINT, inlineGETBETAGRADSHPHIPOINT)(
      CHFp_REALVECT(gradphi)
      ,CHFp_CONST_REALVECT(lmp)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(time) )
{
 CH_TIMELEAF("FORT_GETBETAGRADSHPHIPOINT");
 FORTRAN_NAME( GETBETAGRADSHPHIPOINT ,getbetagradshphipoint )(
      CHFt_REALVECT(gradphi)
      ,CHFt_CONST_REALVECT(lmp)
      ,CHFt_CONST_REALVECT(x)
      ,CHFt_CONST_REAL(time) );
}
#endif  // GUARDGETBETAGRADSHPHIPOINT 

#ifndef GUARDGETSRC 
#define GUARDGETSRC 
// Prototype for Fortran procedure getsrc ...
//
void FORTRAN_NAME( GETSRC ,getsrc )(
      CHFp_FRA1(src)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(time)
      ,CHFp_CONST_REAL(diffconst)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_CONST_REALVECT(probhi)
      ,CHFp_BOX(box) );

#define FORT_GETSRC FORTRAN_NAME( inlineGETSRC, inlineGETSRC)
#define FORTNT_GETSRC FORTRAN_NAME( GETSRC, getsrc)

inline void FORTRAN_NAME(inlineGETSRC, inlineGETSRC)(
      CHFp_FRA1(src)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(time)
      ,CHFp_CONST_REAL(diffconst)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_CONST_REALVECT(probhi)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_GETSRC");
 FORTRAN_NAME( GETSRC ,getsrc )(
      CHFt_FRA1(src)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REAL(time)
      ,CHFt_CONST_REAL(diffconst)
      ,CHFt_CONST_REALVECT(problo)
      ,CHFt_CONST_REALVECT(probhi)
      ,CHFt_BOX(box) );
}
#endif  // GUARDGETSRC 

#ifndef GUARDGETSRCPOINT 
#define GUARDGETSRCPOINT 
// Prototype for Fortran procedure getsrcpoint ...
//
void FORTRAN_NAME( GETSRCPOINT ,getsrcpoint )(
      CHFp_REAL(src)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(time)
      ,CHFp_CONST_REAL(diffconst) );

#define FORT_GETSRCPOINT FORTRAN_NAME( inlineGETSRCPOINT, inlineGETSRCPOINT)
#define FORTNT_GETSRCPOINT FORTRAN_NAME( GETSRCPOINT, getsrcpoint)

inline void FORTRAN_NAME(inlineGETSRCPOINT, inlineGETSRCPOINT)(
      CHFp_REAL(src)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(time)
      ,CHFp_CONST_REAL(diffconst) )
{
 CH_TIMELEAF("FORT_GETSRCPOINT");
 FORTRAN_NAME( GETSRCPOINT ,getsrcpoint )(
      CHFt_REAL(src)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(x)
      ,CHFt_CONST_REAL(time)
      ,CHFt_CONST_REAL(diffconst) );
}
#endif  // GUARDGETSRCPOINT 

#ifndef GUARDGETSHSRC 
#define GUARDGETSHSRC 
// Prototype for Fortran procedure getshsrc ...
//
void FORTRAN_NAME( GETSHSRC ,getshsrc )(
      CHFp_FRA1(src)
      ,CHFp_CONST_REALVECT(lmp)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(time)
      ,CHFp_CONST_REAL(diffconst)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_CONST_REALVECT(probhi)
      ,CHFp_BOX(box) );

#define FORT_GETSHSRC FORTRAN_NAME( inlineGETSHSRC, inlineGETSHSRC)
#define FORTNT_GETSHSRC FORTRAN_NAME( GETSHSRC, getshsrc)

inline void FORTRAN_NAME(inlineGETSHSRC, inlineGETSHSRC)(
      CHFp_FRA1(src)
      ,CHFp_CONST_REALVECT(lmp)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(time)
      ,CHFp_CONST_REAL(diffconst)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_CONST_REALVECT(probhi)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_GETSHSRC");
 FORTRAN_NAME( GETSHSRC ,getshsrc )(
      CHFt_FRA1(src)
      ,CHFt_CONST_REALVECT(lmp)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REAL(time)
      ,CHFt_CONST_REAL(diffconst)
      ,CHFt_CONST_REALVECT(problo)
      ,CHFt_CONST_REALVECT(probhi)
      ,CHFt_BOX(box) );
}
#endif  // GUARDGETSHSRC 

#ifndef GUARDGETSHSRCPOINT 
#define GUARDGETSHSRCPOINT 
// Prototype for Fortran procedure getshsrcpoint ...
//
void FORTRAN_NAME( GETSHSRCPOINT ,getshsrcpoint )(
      CHFp_REAL(src)
      ,CHFp_CONST_REALVECT(lmp)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(time)
      ,CHFp_CONST_REAL(diffconst) );

#define FORT_GETSHSRCPOINT FORTRAN_NAME( inlineGETSHSRCPOINT, inlineGETSHSRCPOINT)
#define FORTNT_GETSHSRCPOINT FORTRAN_NAME( GETSHSRCPOINT, getshsrcpoint)

inline void FORTRAN_NAME(inlineGETSHSRCPOINT, inlineGETSHSRCPOINT)(
      CHFp_REAL(src)
      ,CHFp_CONST_REALVECT(lmp)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(time)
      ,CHFp_CONST_REAL(diffconst) )
{
 CH_TIMELEAF("FORT_GETSHSRCPOINT");
 FORTRAN_NAME( GETSHSRCPOINT ,getshsrcpoint )(
      CHFt_REAL(src)
      ,CHFt_CONST_REALVECT(lmp)
      ,CHFt_CONST_REALVECT(x)
      ,CHFt_CONST_REAL(time)
      ,CHFt_CONST_REAL(diffconst) );
}
#endif  // GUARDGETSHSRCPOINT 

#ifndef GUARDCONVERT_SPHERICAL 
#define GUARDCONVERT_SPHERICAL 
// Prototype for Fortran procedure convert_spherical ...
//
void FORTRAN_NAME( CONVERT_SPHERICAL ,convert_spherical )(
      CHFp_CONST_REALVECT(x)
      ,CHFp_REALVECT(r) );

#define FORT_CONVERT_SPHERICAL FORTRAN_NAME( inlineCONVERT_SPHERICAL, inlineCONVERT_SPHERICAL)
#define FORTNT_CONVERT_SPHERICAL FORTRAN_NAME( CONVERT_SPHERICAL, convert_spherical)

inline void FORTRAN_NAME(inlineCONVERT_SPHERICAL, inlineCONVERT_SPHERICAL)(
      CHFp_CONST_REALVECT(x)
      ,CHFp_REALVECT(r) )
{
 CH_TIMELEAF("FORT_CONVERT_SPHERICAL");
 FORTRAN_NAME( CONVERT_SPHERICAL ,convert_spherical )(
      CHFt_CONST_REALVECT(x)
      ,CHFt_REALVECT(r) );
}
#endif  // GUARDCONVERT_SPHERICAL 

}

#endif
