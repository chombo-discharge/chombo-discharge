#ifndef _EBPATCHGODUNOVF_F_H_
#define _EBPATCHGODUNOVF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDMINFLAT 
#define GUARDMINFLAT 
// Prototype for Fortran procedure minflat ...
//
void FORTRAN_NAME( MINFLAT ,minflat )(
      CHFp_FRA1(flattening)
      ,CHFp_CONST_FRA(zetadir)
      ,CHFp_CONST_FRA(du)
      ,CHFp_BOX(box) );

#define FORT_MINFLAT FORTRAN_NAME( inlineMINFLAT, inlineMINFLAT)
#define FORTNT_MINFLAT FORTRAN_NAME( MINFLAT, minflat)

inline void FORTRAN_NAME(inlineMINFLAT, inlineMINFLAT)(
      CHFp_FRA1(flattening)
      ,CHFp_CONST_FRA(zetadir)
      ,CHFp_CONST_FRA(du)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_MINFLAT");
 FORTRAN_NAME( MINFLAT ,minflat )(
      CHFt_FRA1(flattening)
      ,CHFt_CONST_FRA(zetadir)
      ,CHFt_CONST_FRA(du)
      ,CHFt_BOX(box) );
}
#endif  // GUARDMINFLAT 

#ifndef GUARDGETDPTWO 
#define GUARDGETDPTWO 
// Prototype for Fortran procedure getdptwo ...
//
void FORTRAN_NAME( GETDPTWO ,getdptwo )(
      CHFp_FRA1(delta2p)
      ,CHFp_CONST_FRA1(delta1p)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(lobox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_BOX(hibox)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_BOX(centerbox) );

#define FORT_GETDPTWO FORTRAN_NAME( inlineGETDPTWO, inlineGETDPTWO)
#define FORTNT_GETDPTWO FORTRAN_NAME( GETDPTWO, getdptwo)

inline void FORTRAN_NAME(inlineGETDPTWO, inlineGETDPTWO)(
      CHFp_FRA1(delta2p)
      ,CHFp_CONST_FRA1(delta1p)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(lobox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_BOX(hibox)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_BOX(centerbox) )
{
 CH_TIMELEAF("FORT_GETDPTWO");
 FORTRAN_NAME( GETDPTWO ,getdptwo )(
      CHFt_FRA1(delta2p)
      ,CHFt_CONST_FRA1(delta1p)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(lobox)
      ,CHFt_CONST_INT(haslo)
      ,CHFt_BOX(hibox)
      ,CHFt_CONST_INT(hashi)
      ,CHFt_BOX(centerbox) );
}
#endif  // GUARDGETDPTWO 

#ifndef GUARDGETFLAT 
#define GUARDGETFLAT 
// Prototype for Fortran procedure getflat ...
//
void FORTRAN_NAME( GETFLAT ,getflat )(
      CHFp_FRA1(zetatwiddle)
      ,CHFp_CONST_FRA1(delta1p)
      ,CHFp_CONST_FRA1(delta2p)
      ,CHFp_CONST_FRA1(bulkmin)
      ,CHFp_BOX(box) );

#define FORT_GETFLAT FORTRAN_NAME( inlineGETFLAT, inlineGETFLAT)
#define FORTNT_GETFLAT FORTRAN_NAME( GETFLAT, getflat)

inline void FORTRAN_NAME(inlineGETFLAT, inlineGETFLAT)(
      CHFp_FRA1(zetatwiddle)
      ,CHFp_CONST_FRA1(delta1p)
      ,CHFp_CONST_FRA1(delta2p)
      ,CHFp_CONST_FRA1(bulkmin)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_GETFLAT");
 FORTRAN_NAME( GETFLAT ,getflat )(
      CHFt_FRA1(zetatwiddle)
      ,CHFt_CONST_FRA1(delta1p)
      ,CHFt_CONST_FRA1(delta2p)
      ,CHFt_CONST_FRA1(bulkmin)
      ,CHFt_BOX(box) );
}
#endif  // GUARDGETFLAT 

#ifndef GUARDGETGRAD 
#define GUARDGETGRAD 
// Prototype for Fortran procedure getgrad ...
//
void FORTRAN_NAME( GETGRAD ,getgrad )(
      CHFp_FRA1(du)
      ,CHFp_CONST_FRA1(u)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(lobox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_BOX(hibox)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_BOX(centerbox) );

#define FORT_GETGRAD FORTRAN_NAME( inlineGETGRAD, inlineGETGRAD)
#define FORTNT_GETGRAD FORTRAN_NAME( GETGRAD, getgrad)

inline void FORTRAN_NAME(inlineGETGRAD, inlineGETGRAD)(
      CHFp_FRA1(du)
      ,CHFp_CONST_FRA1(u)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(lobox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_BOX(hibox)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_BOX(centerbox) )
{
 CH_TIMELEAF("FORT_GETGRAD");
 FORTRAN_NAME( GETGRAD ,getgrad )(
      CHFt_FRA1(du)
      ,CHFt_CONST_FRA1(u)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(lobox)
      ,CHFt_CONST_INT(haslo)
      ,CHFt_BOX(hibox)
      ,CHFt_CONST_INT(hashi)
      ,CHFt_BOX(centerbox) );
}
#endif  // GUARDGETGRAD 

#ifndef GUARDGETRELATIVEGRAD 
#define GUARDGETRELATIVEGRAD 
// Prototype for Fortran procedure getrelativegrad ...
//
void FORTRAN_NAME( GETRELATIVEGRAD ,getrelativegrad )(
      CHFp_FRA1(du)
      ,CHFp_CONST_FRA1(u)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(lobox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_BOX(hibox)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_BOX(centerbox) );

#define FORT_GETRELATIVEGRAD FORTRAN_NAME( inlineGETRELATIVEGRAD, inlineGETRELATIVEGRAD)
#define FORTNT_GETRELATIVEGRAD FORTRAN_NAME( GETRELATIVEGRAD, getrelativegrad)

inline void FORTRAN_NAME(inlineGETRELATIVEGRAD, inlineGETRELATIVEGRAD)(
      CHFp_FRA1(du)
      ,CHFp_CONST_FRA1(u)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(lobox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_BOX(hibox)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_BOX(centerbox) )
{
 CH_TIMELEAF("FORT_GETRELATIVEGRAD");
 FORTRAN_NAME( GETRELATIVEGRAD ,getrelativegrad )(
      CHFt_FRA1(du)
      ,CHFt_CONST_FRA1(u)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(lobox)
      ,CHFt_CONST_INT(haslo)
      ,CHFt_BOX(hibox)
      ,CHFt_CONST_INT(hashi)
      ,CHFt_BOX(centerbox) );
}
#endif  // GUARDGETRELATIVEGRAD 

#ifndef GUARDMAGNITUDE 
#define GUARDMAGNITUDE 
// Prototype for Fortran procedure magnitude ...
//
void FORTRAN_NAME( MAGNITUDE ,magnitude )(
      CHFp_FRA1(magdata)
      ,CHFp_CONST_FRA(data)
      ,CHFp_BOX(box) );

#define FORT_MAGNITUDE FORTRAN_NAME( inlineMAGNITUDE, inlineMAGNITUDE)
#define FORTNT_MAGNITUDE FORTRAN_NAME( MAGNITUDE, magnitude)

inline void FORTRAN_NAME(inlineMAGNITUDE, inlineMAGNITUDE)(
      CHFp_FRA1(magdata)
      ,CHFp_CONST_FRA(data)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_MAGNITUDE");
 FORTRAN_NAME( MAGNITUDE ,magnitude )(
      CHFt_FRA1(magdata)
      ,CHFt_CONST_FRA(data)
      ,CHFt_BOX(box) );
}
#endif  // GUARDMAGNITUDE 

#ifndef GUARDMIN3PTS 
#define GUARDMIN3PTS 
// Prototype for Fortran procedure min3pts ...
//
void FORTRAN_NAME( MIN3PTS ,min3pts )(
      CHFp_FRA1(mindata)
      ,CHFp_CONST_FRA1(data)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(lobox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_BOX(hibox)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_BOX(centerbox) );

#define FORT_MIN3PTS FORTRAN_NAME( inlineMIN3PTS, inlineMIN3PTS)
#define FORTNT_MIN3PTS FORTRAN_NAME( MIN3PTS, min3pts)

inline void FORTRAN_NAME(inlineMIN3PTS, inlineMIN3PTS)(
      CHFp_FRA1(mindata)
      ,CHFp_CONST_FRA1(data)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(lobox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_BOX(hibox)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_BOX(centerbox) )
{
 CH_TIMELEAF("FORT_MIN3PTS");
 FORTRAN_NAME( MIN3PTS ,min3pts )(
      CHFt_FRA1(mindata)
      ,CHFt_CONST_FRA1(data)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(lobox)
      ,CHFt_CONST_INT(haslo)
      ,CHFt_BOX(hibox)
      ,CHFt_CONST_INT(hashi)
      ,CHFt_BOX(centerbox) );
}
#endif  // GUARDMIN3PTS 

#ifndef GUARDSECONDSLOPEDIFFS 
#define GUARDSECONDSLOPEDIFFS 
// Prototype for Fortran procedure secondslopediffs ...
//
void FORTRAN_NAME( SECONDSLOPEDIFFS ,secondslopediffs )(
      CHFp_FRA(deltawc)
      ,CHFp_FRA(deltawl)
      ,CHFp_FRA(deltawr)
      ,CHFp_CONST_FRA(w)
      ,CHFp_CONST_INT(numslopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(lobox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_BOX(hibox)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_BOX(centerbox) );

#define FORT_SECONDSLOPEDIFFS FORTRAN_NAME( inlineSECONDSLOPEDIFFS, inlineSECONDSLOPEDIFFS)
#define FORTNT_SECONDSLOPEDIFFS FORTRAN_NAME( SECONDSLOPEDIFFS, secondslopediffs)

inline void FORTRAN_NAME(inlineSECONDSLOPEDIFFS, inlineSECONDSLOPEDIFFS)(
      CHFp_FRA(deltawc)
      ,CHFp_FRA(deltawl)
      ,CHFp_FRA(deltawr)
      ,CHFp_CONST_FRA(w)
      ,CHFp_CONST_INT(numslopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(lobox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_BOX(hibox)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_BOX(centerbox) )
{
 CH_TIMELEAF("FORT_SECONDSLOPEDIFFS");
 FORTRAN_NAME( SECONDSLOPEDIFFS ,secondslopediffs )(
      CHFt_FRA(deltawc)
      ,CHFt_FRA(deltawl)
      ,CHFt_FRA(deltawr)
      ,CHFt_CONST_FRA(w)
      ,CHFt_CONST_INT(numslopes)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(lobox)
      ,CHFt_CONST_INT(haslo)
      ,CHFt_BOX(hibox)
      ,CHFt_CONST_INT(hashi)
      ,CHFt_BOX(centerbox) );
}
#endif  // GUARDSECONDSLOPEDIFFS 

#ifndef GUARDFORTHSLOPEDIFFS 
#define GUARDFORTHSLOPEDIFFS 
// Prototype for Fortran procedure forthslopediffs ...
//
void FORTRAN_NAME( FORTHSLOPEDIFFS ,forthslopediffs )(
      CHFp_FRA(delta4wc)
      ,CHFp_CONST_FRA(w)
      ,CHFp_CONST_FRA(delta2w)
      ,CHFp_CONST_INT(numslopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(lobox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_BOX(hibox)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_BOX(centerbox) );

#define FORT_FORTHSLOPEDIFFS FORTRAN_NAME( inlineFORTHSLOPEDIFFS, inlineFORTHSLOPEDIFFS)
#define FORTNT_FORTHSLOPEDIFFS FORTRAN_NAME( FORTHSLOPEDIFFS, forthslopediffs)

inline void FORTRAN_NAME(inlineFORTHSLOPEDIFFS, inlineFORTHSLOPEDIFFS)(
      CHFp_FRA(delta4wc)
      ,CHFp_CONST_FRA(w)
      ,CHFp_CONST_FRA(delta2w)
      ,CHFp_CONST_INT(numslopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(lobox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_BOX(hibox)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_BOX(centerbox) )
{
 CH_TIMELEAF("FORT_FORTHSLOPEDIFFS");
 FORTRAN_NAME( FORTHSLOPEDIFFS ,forthslopediffs )(
      CHFt_FRA(delta4wc)
      ,CHFt_CONST_FRA(w)
      ,CHFt_CONST_FRA(delta2w)
      ,CHFt_CONST_INT(numslopes)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(lobox)
      ,CHFt_CONST_INT(haslo)
      ,CHFt_BOX(hibox)
      ,CHFt_CONST_INT(hashi)
      ,CHFt_BOX(centerbox) );
}
#endif  // GUARDFORTHSLOPEDIFFS 

#ifndef GUARDAPPLYFLAT 
#define GUARDAPPLYFLAT 
// Prototype for Fortran procedure applyflat ...
//
void FORTRAN_NAME( APPLYFLAT ,applyflat )(
      CHFp_FRA(dw)
      ,CHFp_CONST_FRA1(flattening)
      ,CHFp_CONST_INT(numslopes)
      ,CHFp_BOX(box) );

#define FORT_APPLYFLAT FORTRAN_NAME( inlineAPPLYFLAT, inlineAPPLYFLAT)
#define FORTNT_APPLYFLAT FORTRAN_NAME( APPLYFLAT, applyflat)

inline void FORTRAN_NAME(inlineAPPLYFLAT, inlineAPPLYFLAT)(
      CHFp_FRA(dw)
      ,CHFp_CONST_FRA1(flattening)
      ,CHFp_CONST_INT(numslopes)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_APPLYFLAT");
 FORTRAN_NAME( APPLYFLAT ,applyflat )(
      CHFt_FRA(dw)
      ,CHFt_CONST_FRA1(flattening)
      ,CHFt_CONST_INT(numslopes)
      ,CHFt_BOX(box) );
}
#endif  // GUARDAPPLYFLAT 

#ifndef GUARDINCSOURCE 
#define GUARDINCSOURCE 
// Prototype for Fortran procedure incsource ...
//
void FORTRAN_NAME( INCSOURCE ,incsource )(
      CHFp_FRA(prim)
      ,CHFp_CONST_FRA(source)
      ,CHFp_CONST_REAL(scale)
      ,CHFp_BOX(dcalc) );

#define FORT_INCSOURCE FORTRAN_NAME( inlineINCSOURCE, inlineINCSOURCE)
#define FORTNT_INCSOURCE FORTRAN_NAME( INCSOURCE, incsource)

inline void FORTRAN_NAME(inlineINCSOURCE, inlineINCSOURCE)(
      CHFp_FRA(prim)
      ,CHFp_CONST_FRA(source)
      ,CHFp_CONST_REAL(scale)
      ,CHFp_BOX(dcalc) )
{
 CH_TIMELEAF("FORT_INCSOURCE");
 FORTRAN_NAME( INCSOURCE ,incsource )(
      CHFt_FRA(prim)
      ,CHFt_CONST_FRA(source)
      ,CHFt_CONST_REAL(scale)
      ,CHFt_BOX(dcalc) );
}
#endif  // GUARDINCSOURCE 

#ifndef GUARDVLLIMITER 
#define GUARDVLLIMITER 
// Prototype for Fortran procedure vllimiter ...
//
void FORTRAN_NAME( VLLIMITER ,vllimiter )(
      CHFp_FRA(slopeprim)
      ,CHFp_CONST_FRA(slopeleft)
      ,CHFp_CONST_FRA(sloperigh)
      ,CHFp_BOX(dcalc) );

#define FORT_VLLIMITER FORTRAN_NAME( inlineVLLIMITER, inlineVLLIMITER)
#define FORTNT_VLLIMITER FORTRAN_NAME( VLLIMITER, vllimiter)

inline void FORTRAN_NAME(inlineVLLIMITER, inlineVLLIMITER)(
      CHFp_FRA(slopeprim)
      ,CHFp_CONST_FRA(slopeleft)
      ,CHFp_CONST_FRA(sloperigh)
      ,CHFp_BOX(dcalc) )
{
 CH_TIMELEAF("FORT_VLLIMITER");
 FORTRAN_NAME( VLLIMITER ,vllimiter )(
      CHFt_FRA(slopeprim)
      ,CHFt_CONST_FRA(slopeleft)
      ,CHFt_CONST_FRA(sloperigh)
      ,CHFt_BOX(dcalc) );
}
#endif  // GUARDVLLIMITER 

#ifndef GUARDPOINTVLLIMITER 
#define GUARDPOINTVLLIMITER 
// Prototype for Fortran procedure pointvllimiter ...
//
void FORTRAN_NAME( POINTVLLIMITER ,pointvllimiter )(
      CHFp_REAL(dqlim)
      ,CHFp_CONST_REAL(dql)
      ,CHFp_CONST_REAL(dqr) );

#define FORT_POINTVLLIMITER FORTRAN_NAME( inlinePOINTVLLIMITER, inlinePOINTVLLIMITER)
#define FORTNT_POINTVLLIMITER FORTRAN_NAME( POINTVLLIMITER, pointvllimiter)

inline void FORTRAN_NAME(inlinePOINTVLLIMITER, inlinePOINTVLLIMITER)(
      CHFp_REAL(dqlim)
      ,CHFp_CONST_REAL(dql)
      ,CHFp_CONST_REAL(dqr) )
{
 CH_TIMELEAF("FORT_POINTVLLIMITER");
 FORTRAN_NAME( POINTVLLIMITER ,pointvllimiter )(
      CHFt_REAL(dqlim)
      ,CHFt_CONST_REAL(dql)
      ,CHFt_CONST_REAL(dqr) );
}
#endif  // GUARDPOINTVLLIMITER 

#ifndef GUARDDIVUEDGE 
#define GUARDDIVUEDGE 
// Prototype for Fortran procedure divuedge ...
//
void FORTRAN_NAME( DIVUEDGE ,divuedge )(
      CHFp_FRA1(divu)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(lobox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_BOX(hibox)
      ,CHFp_CONST_INT(hashi) );

#define FORT_DIVUEDGE FORTRAN_NAME( inlineDIVUEDGE, inlineDIVUEDGE)
#define FORTNT_DIVUEDGE FORTRAN_NAME( DIVUEDGE, divuedge)

inline void FORTRAN_NAME(inlineDIVUEDGE, inlineDIVUEDGE)(
      CHFp_FRA1(divu)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(lobox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_BOX(hibox)
      ,CHFp_CONST_INT(hashi) )
{
 CH_TIMELEAF("FORT_DIVUEDGE");
 FORTRAN_NAME( DIVUEDGE ,divuedge )(
      CHFt_FRA1(divu)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_BOX(lobox)
      ,CHFt_CONST_INT(haslo)
      ,CHFt_BOX(hibox)
      ,CHFt_CONST_INT(hashi) );
}
#endif  // GUARDDIVUEDGE 

#ifndef GUARDAVEFLUXTOFACE 
#define GUARDAVEFLUXTOFACE 
// Prototype for Fortran procedure avefluxtoface ...
//
void FORTRAN_NAME( AVEFLUXTOFACE ,avefluxtoface )(
      CHFp_FRA1(faceflux)
      ,CHFp_CONST_FRA1(ccflux)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(facebox) );

#define FORT_AVEFLUXTOFACE FORTRAN_NAME( inlineAVEFLUXTOFACE, inlineAVEFLUXTOFACE)
#define FORTNT_AVEFLUXTOFACE FORTRAN_NAME( AVEFLUXTOFACE, avefluxtoface)

inline void FORTRAN_NAME(inlineAVEFLUXTOFACE, inlineAVEFLUXTOFACE)(
      CHFp_FRA1(faceflux)
      ,CHFp_CONST_FRA1(ccflux)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(facebox) )
{
 CH_TIMELEAF("FORT_AVEFLUXTOFACE");
 FORTRAN_NAME( AVEFLUXTOFACE ,avefluxtoface )(
      CHFt_FRA1(faceflux)
      ,CHFt_CONST_FRA1(ccflux)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_BOX(facebox) );
}
#endif  // GUARDAVEFLUXTOFACE 

#ifndef GUARDDIVUONED 
#define GUARDDIVUONED 
// Prototype for Fortran procedure divuoned ...
//
void FORTRAN_NAME( DIVUONED ,divuoned )(
      CHFp_FRA1(divu)
      ,CHFp_CONST_FRA1(velnorm)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(centerbox) );

#define FORT_DIVUONED FORTRAN_NAME( inlineDIVUONED, inlineDIVUONED)
#define FORTNT_DIVUONED FORTRAN_NAME( DIVUONED, divuoned)

inline void FORTRAN_NAME(inlineDIVUONED, inlineDIVUONED)(
      CHFp_FRA1(divu)
      ,CHFp_CONST_FRA1(velnorm)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(centerbox) )
{
 CH_TIMELEAF("FORT_DIVUONED");
 FORTRAN_NAME( DIVUONED ,divuoned )(
      CHFt_FRA1(divu)
      ,CHFt_CONST_FRA1(velnorm)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_BOX(centerbox) );
}
#endif  // GUARDDIVUONED 

#ifndef GUARDDIVUTRAN 
#define GUARDDIVUTRAN 
// Prototype for Fortran procedure divutran ...
//
void FORTRAN_NAME( DIVUTRAN ,divutran )(
      CHFp_FRA1(divu)
      ,CHFp_CONST_FRA1(slopevel)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(centerbox) );

#define FORT_DIVUTRAN FORTRAN_NAME( inlineDIVUTRAN, inlineDIVUTRAN)
#define FORTNT_DIVUTRAN FORTRAN_NAME( DIVUTRAN, divutran)

inline void FORTRAN_NAME(inlineDIVUTRAN, inlineDIVUTRAN)(
      CHFp_FRA1(divu)
      ,CHFp_CONST_FRA1(slopevel)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(centerbox) )
{
 CH_TIMELEAF("FORT_DIVUTRAN");
 FORTRAN_NAME( DIVUTRAN ,divutran )(
      CHFt_FRA1(divu)
      ,CHFt_CONST_FRA1(slopevel)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_BOX(centerbox) );
}
#endif  // GUARDDIVUTRAN 

#ifndef GUARDARTVISC 
#define GUARDARTVISC 
// Prototype for Fortran procedure artvisc ...
//
void FORTRAN_NAME( ARTVISC ,artvisc )(
      CHFp_FRA(f)
      ,CHFp_CONST_FRA(u)
      ,CHFp_CONST_FRA1(divu)
      ,CHFp_CONST_REAL(coeff)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(numcons)
      ,CHFp_CONST_REAL(dx) );

#define FORT_ARTVISC FORTRAN_NAME( inlineARTVISC, inlineARTVISC)
#define FORTNT_ARTVISC FORTRAN_NAME( ARTVISC, artvisc)

inline void FORTRAN_NAME(inlineARTVISC, inlineARTVISC)(
      CHFp_FRA(f)
      ,CHFp_CONST_FRA(u)
      ,CHFp_CONST_FRA1(divu)
      ,CHFp_CONST_REAL(coeff)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(numcons)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_ARTVISC");
 FORTRAN_NAME( ARTVISC ,artvisc )(
      CHFt_FRA(f)
      ,CHFt_CONST_FRA(u)
      ,CHFt_CONST_FRA1(divu)
      ,CHFt_CONST_REAL(coeff)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box)
      ,CHFt_CONST_INT(numcons)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDARTVISC 

#ifndef GUARDUPDATE 
#define GUARDUPDATE 
// Prototype for Fortran procedure update ...
//
void FORTRAN_NAME( UPDATE ,update )(
      CHFp_BOX(dcalc)
      ,CHFp_FRA(state)
      ,CHFp_CONST_FRA(flux)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(nconserved)
      ,CHFp_CONST_REAL(dtbydx) );

#define FORT_UPDATE FORTRAN_NAME( inlineUPDATE, inlineUPDATE)
#define FORTNT_UPDATE FORTRAN_NAME( UPDATE, update)

inline void FORTRAN_NAME(inlineUPDATE, inlineUPDATE)(
      CHFp_BOX(dcalc)
      ,CHFp_FRA(state)
      ,CHFp_CONST_FRA(flux)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(nconserved)
      ,CHFp_CONST_REAL(dtbydx) )
{
 CH_TIMELEAF("FORT_UPDATE");
 FORTRAN_NAME( UPDATE ,update )(
      CHFt_BOX(dcalc)
      ,CHFt_FRA(state)
      ,CHFt_CONST_FRA(flux)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_CONST_INT(nconserved)
      ,CHFt_CONST_REAL(dtbydx) );
}
#endif  // GUARDUPDATE 

#ifndef GUARDDIVERGEF 
#define GUARDDIVERGEF 
// Prototype for Fortran procedure divergef ...
//
void FORTRAN_NAME( DIVERGEF ,divergef )(
      CHFp_BOX(dcalc)
      ,CHFp_FRA(divf)
      ,CHFp_CONST_FRA(flux)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(nconserved)
      ,CHFp_CONST_REAL(dx) );

#define FORT_DIVERGEF FORTRAN_NAME( inlineDIVERGEF, inlineDIVERGEF)
#define FORTNT_DIVERGEF FORTRAN_NAME( DIVERGEF, divergef)

inline void FORTRAN_NAME(inlineDIVERGEF, inlineDIVERGEF)(
      CHFp_BOX(dcalc)
      ,CHFp_FRA(divf)
      ,CHFp_CONST_FRA(flux)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(nconserved)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_DIVERGEF");
 FORTRAN_NAME( DIVERGEF ,divergef )(
      CHFt_BOX(dcalc)
      ,CHFt_FRA(divf)
      ,CHFt_CONST_FRA(flux)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_CONST_INT(nconserved)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDDIVERGEF 

#ifndef GUARDREGUPDATE 
#define GUARDREGUPDATE 
// Prototype for Fortran procedure regupdate ...
//
void FORTRAN_NAME( REGUPDATE ,regupdate )(
      CHFp_BOX(dcalc)
      ,CHFp_FRA(consstate)
      ,CHFp_CONST_FRA(divf)
      ,CHFp_CONST_INT(nconserved)
      ,CHFp_CONST_REAL(dt) );

#define FORT_REGUPDATE FORTRAN_NAME( inlineREGUPDATE, inlineREGUPDATE)
#define FORTNT_REGUPDATE FORTRAN_NAME( REGUPDATE, regupdate)

inline void FORTRAN_NAME(inlineREGUPDATE, inlineREGUPDATE)(
      CHFp_BOX(dcalc)
      ,CHFp_FRA(consstate)
      ,CHFp_CONST_FRA(divf)
      ,CHFp_CONST_INT(nconserved)
      ,CHFp_CONST_REAL(dt) )
{
 CH_TIMELEAF("FORT_REGUPDATE");
 FORTRAN_NAME( REGUPDATE ,regupdate )(
      CHFt_BOX(dcalc)
      ,CHFt_FRA(consstate)
      ,CHFt_CONST_FRA(divf)
      ,CHFt_CONST_INT(nconserved)
      ,CHFt_CONST_REAL(dt) );
}
#endif  // GUARDREGUPDATE 

}

#endif
